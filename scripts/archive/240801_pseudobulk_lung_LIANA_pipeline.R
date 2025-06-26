# # Install packages ####
# install.packages("renv")
renv::init() #Select 1
install.packages("stringi")
install.packages("./packages/Matrix_1.6-4.tar.gz", repos = NULL, type = "source")
install.packages("./packages/SeuratObject_5.0.1.tar.gz", repos = NULL, type = "source")
install.packages("./packages/Seurat_5.0.3.tar.gz", repos = NULL, type = "source")
install.packages("./packages/OmnipathR-master_april_2024.tar.gz", repos = NULL, type = "source")
install.packages("./packages/liana-master_sep_2024.tar.gz", repos = NULL, type = "source", verbose = T) #might take a while
install.packages("./packages/SCpubr-1.1.2-dev-stable.tar.gz", repos = NULL, type = "source", verbose = T) #might take a while
renv::snapshot(type = "all")# exclude = c("liana", "OmnipathR"))

get_freq <- function(liana_trunc) {
  all_cell_types <- unique(c(liana_trunc$source, liana_trunc$target))
  freq_table <- table(liana_trunc$source, liana_trunc$target)
  freq_matrix <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                        dimnames = list(all_cell_types, all_cell_types))
  freq_matrix[rownames(freq_table), colnames(freq_table)] <- freq_table
  return(freq_matrix)
}
heat_freq <- function(liana_trunc, ...) {
  freqs <- get_freq(liana_trunc)
  liana_heatmap(mat = freqs, ...)
}


# Load packages ####
library(tidyverse)
library(Seurat)
library(qs) #archive for saving objects
library(ProfilerAPI2)
library(liana)
suppressPackageStartupMessages(library('SCpubr'))
conflicts_prefer(dplyr::filter)


# Add Merck Colors ####
merck_colors <- c(violet = "#503291", blue = "#0F69AF", green = "#149B5F", red = "#E61E50", 
                  pink = "#EB3C96", lightblue = "#2DBECD", lightgreen = "#A5CD50",
                  yellow = "#FFC832", palepink = "#E1C3CD", paleblue = "#96D7D2",
                  palegreen = "#B4DC96", paleyellow = "#FFDCB9", grey = "#999999")
# 
# 
# ##### Get spatial data from datalake ####
# api = ProfilerAPI2::profiler_api(profile = "default")   ## needed for Rmarkdown connection to ProfilerAPI2
# 
# Views_spatial <- api$studies$list() %>%
#   filter(grepl("Publication - Spatial", name)) %>%
#   pull("id")
# 
# profiler_folder_id_spatial <- api$data_lake$object_id(Views_spatial, '')
# 
# Datasets_spatial <- api$data_lake$list(profiler_folder_id_spatial, recursive = TRUE) %>%
#   dplyr::filter(grepl(".qs", name)) %>%
#   as.data.frame()
# 
# api$data_lake$download_file(Datasets_spatial$id[1], path = "./seurat_files_spatial/", overwrite = TRUE, show_progress = T)
# 
# 
##### Get scRNAseq data from datalake ####
api = ProfilerAPI2::profiler_api(profile = "default")   ## needed for Rmarkdown connection to ProfilerAPI2

Views_singlecell <- api$studies$list() %>%
  filter(grepl("Publication - Single Cell", name)) %>%
  pull("id")

profiler_folder_id_sc <- api$data_lake$object_id(Views_singlecell, '/App_Development_SPOT')

Datasets_singlecell <- api$data_lake$list(profiler_folder_id_sc, recursive = TRUE) %>%
  dplyr::filter(grepl(".qs", name)) %>%
  as.data.frame()

if (!dir.exists("seurat_files")) {dir.create("seurat_files")}
api$data_lake$download_file("file-5f5d34df-f00e-41df-b935-4d8697f99448", path = "./seurat_files/", overwrite = TRUE, show_progress = T)


######## Analysis ########
impactQ3 <- qread("./seurat_files/Seurat_IMPACT_Atlas_Q3_2024.qs")

# Subset to only look at cancer data and indication type in particular ####
impactQ3 <- impactQ3[,which(impactQ3$sample_status == 'Cancer')]
impactQ3 <- impactQ3[,which(!is.na(impactQ3$cell_type_basic))]
impactQ3 <- impactQ3[,which(impactQ3$sample_disease_oncotree %in% c('SCLC', 'NSCLC', 'LUSC', 'LUAD'))]
qsave(impactQ3, "./seurat_files/pseudobulk_lung_Q3.qs")
impactQ3 <- qread("./seurat_files/pseudobulk_lung_Q3.qs")

# Rename active.ident
current_cell_types <- levels(impactQ3@active.ident)

new_cell_types <- setNames(c("B Cell", "Other", "Other", "Endothelial Cell", "Epithelial Cell", 
                             "Fibroblast", "Other", "Other", "Malignant Cell", "Monocyte/Macrophage", 
                             "Specialized Immune Cell", "Other", "Other", "Other", 
                             "Specialized Immune Cell", "T Cell"), current_cell_types)

pseudo_lung <- RenameIdents(impactQ3, new_cell_types)
qsave(pseudo_lung, "./seurat_files/renamed_pseudobulk_lung_Q3.qs")

# Load pre-subsetted seurat file ####
pseudo_lung <- qread("./seurat_files/renamed_pseudobulk_lung_Q3.qs")

######## LIANA: ALL METHODS ########
# This is a very long step! Especially if you have many cells!  ####
liana_data <- liana_wrap(pseudo_lung,
                         permutation.params = list(nperms=100,
                                                   parallelize=FALSE,
                                                   workers=8),
                         return_all = T,
                         min_cells = 0)

# Aggregate and Obtain Consensus Ranks####
liana_data_aggregated <- liana_data %>%
  liana_aggregate() %>%
  mutate(original_rank = row_number())

# Save and Load LIANA pre-saved results ####
qsave(liana_data,"./liana_files/240902_pseudo_lung_Q3_liana_data.qs")
liana_data <- qread("./liana_files/240902_pseudo_lung_Q3_liana_data.qs")

qsave(liana_data_aggregated,"./liana_files/240902_pseudo_lung_Q3_liana_data_aggregated.qs")
liana_data_aggregated <- qread("./liana_files/240902_pseudo_lung_Q3_liana_data_aggregated.qs")


######## Visualisation ########

# all LIANA methods, unfiltered  ####
p1 <- liana_data_aggregated %>%
  #  filter(pvalue <= 0.05) %>%
  liana_dotplot(source_groups = "Malignant Cell",
                target_groups = c("B Cell", "Other", "Endothelial Cell", "Epithelial Cell", "Fibroblast", "Monocyte/Macrophage", "Specialized Immune Cell", "T Cell"),
                ntop = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 12),
        plot.caption.position = "plot",
        plot.subtitle = element_text(hjust = 0.5, size = 23),
        plot.title = element_text(face = "bold", size = 30, color = merck_colors[1]))  +
  labs(caption = "Dataset: IMPACT_lung_cancers", title = "Top 20 unfiltered interactions", subtitle = "Source")

ggsave("01_unfiltered_interactions.png", plot = p1, create.dir = T, path = "./res/", height = 13, width = 12)


p2 <- liana_data_aggregated %>%
  liana_dotplot(target_groups = "Malignant Cell",
                source_groups = c("B Cell", "Other", "Endothelial Cell", "Epithelial Cell", "Fibroblast", "Monocyte/Macrophage", "Specialized Immune Cell", "T Cell"),
                ntop = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 12),
        plot.caption.position = "plot",
        plot.subtitle = element_text(hjust = 0.5, size = 23),
        plot.title = element_text(face = "bold", size = 30, color = merck_colors[1])) +
  labs(caption = "Dataset: IMPACT_lung_cancers", title = "Top 20 unfiltered interactions", subtitle = "Source")

ggsave("02_unfiltered_interactions.png", plot = p2, create.dir = T, path = "./res/", height = 13, width = 12)

# all LIANA methods, FILTERED  ####
p1 <- liana_data_aggregated %>%
  filter(aggregate_rank <= 0.001) %>%
  liana_dotplot(source_groups = "Malignant Cell",
                target_groups = c("B Cell", "Other", "Endothelial Cell", "Epithelial Cell", "Fibroblast", "Monocyte/Macrophage", "Specialized Immune Cell", "T Cell"),
                ntop = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 12),
        plot.caption.position = "plot",
        plot.subtitle = element_text(hjust = 0.5, size = 23),
        plot.title = element_text(face = "bold", size = 30, color = merck_colors[1]))  +
  labs(caption = "Dataset: IMPACT_lung_cancers", title = "Top 20 interactions, p ≤ 0.001", subtitle = "Source")

ggsave("03_filtered_interactions.png", plot = p1, create.dir = T, path = "./res/", height = 13, width = 12)


p2 <- liana_data_aggregated %>%
  filter(aggregate_rank <= 0.001) %>%
  liana_dotplot(target_groups = "Malignant Cell",
                source_groups = c("B Cell", "Other", "Endothelial Cell", "Epithelial Cell", "Fibroblast", "Monocyte/Macrophage", "Specialized Immune Cell", "T Cell"),
                ntop = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 12),
        plot.caption.position = "plot",
        plot.subtitle = element_text(hjust = 0.5, size = 23),
        plot.title = element_text(face = "bold", size = 30, color = merck_colors[1])) +
  labs(caption = "Dataset: IMPACT_lung_cancers", title = "Top 20 interactions, p ≤ 0.001", subtitle = "Source")

ggsave("04_filtered_interactions.png", plot = p2, create.dir = T, path = "./res/", height = 13, width = 12)




#############################

liana_trunc <- liana_data_aggregated %>% # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.05) # this can be FDR-corr if n is too high

# Frequency Heatmap  ####
png(file="./res/03_filtered_freq_heatmap.png",
    width=750, height=450, res = 100, units = "px")
heat_freq(liana_trunc)
dev.off()

# SCpubr plots, arranged by specificity and expression ####
out_both <- SCpubr::do_LigandReceptorPlot(liana_output = liana_data,
                                     top_interactions = 25,
                                     compute_ChordDiagrams = TRUE,
                                     arrange_interactions_by = "both")

png("./res/04_SCpubr_chord_total_interactions.png", height = 13, width = 12, units = "in", res = 300)
  print(out_both$chord_total_interactions)  # or recreate the plot here
  dev.off()
png("./res/05_SCpubr_chord_LR_interactions.png", height = 13, width = 12, units = "in", res = 300)
  print(out_both$chord_ligand_receptor)  # or recreate the plot here
  dev.off()
ggsave("06_SCpubr_dotplot.png", plot = out_both$dotplot, create.dir = T, path = "./res/", height = 13, width = 20)


# SCpubr plots, arranged by aggregate_rank (p-value) ####
out <- SCpubr::do_LigandReceptorPlot(liana_output = liana_data,
                                          top_interactions = 25,
                                          compute_ChordDiagrams = TRUE,
                                          arrange_interactions_by = "aggregate_rank")

png("./res/07_SCpubr_chord_total_interactions.png", height = 13, width = 12, units = "in", res = 300)
  print(out$chord_total_interactions)  # or recreate the plot here
  dev.off()
png("./res/08_SCpubr_chord_LR_interactions.png", height = 13, width = 12, units = "in", res = 300)
  print(out$chord_ligand_receptor)  # or recreate the plot here
  dev.off()
ggsave("09_SCpubr_dotplot.png", plot = out$dotplot, create.dir = T, path = "./res/", height = 13, width = 20)


###### CellPhoneDB only ######
# identify interactions of interest ####
liana_interactions <- liana_data$cellphonedb %>%
  # only keep interactions with p-val <= 0.05
  filter(pvalue <= 0.05) %>% # this reflects interactions `specificity`
  # then rank according to `magnitude` or 'specificity'
  rank_method(method_name = "cellphonedb",
              mode = "magnitude") %>%
  # keep top 20 interactions (regardless of cell type)
  distinct_at(c("ligand.complex", "receptor.complex")) %>%
  head(20)

p3 <- liana_interactions %>%
  # keep only the interactions of interest
  inner_join(liana_data$cellphonedb, 
             by = c("ligand.complex", "receptor.complex")) %>%
  # invert size (low p-value/high specificity = larger dot size)
  # + add a small value to avoid Infinity for 0s
  mutate(pvalue = -log10(pvalue + 1e-10)) %>% 
  liana_dotplot(source_groups = "Malignant Cell",
                specificity = "pvalue",
                magnitude = "lr.mean",
                show_complex = TRUE,
                size.label = "-log10(p-value)",
                ntop =10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 12),
        plot.caption.position = "plot",
        plot.subtitle = element_text(hjust = 0.5, size = 23),
        plot.title = element_text(face = "bold", size = 30, color = merck_colors[1])) +
  labs(caption = "Dataset: IMPACT_lung_cancers", title = "Top 10 CellPhoneDB interactions (pvalue <= 0.05)", subtitle = "Source: Malignant")

ggsave("10_CellPhoneDB_filtered_interactions.png", plot = p3, create.dir = T, path = "./res/", height = 13, width = 12)



p4 <- liana_interactions %>%
  # keep only the interactions of interest
  inner_join(liana_data$cellphonedb, 
             by = c("ligand.complex", "receptor.complex")) %>%
  # invert size (low p-value/high specificity = larger dot size)
  # + add a small value to avoid Infinity for 0s
  mutate(pvalue = -log10(pvalue + 1e-10)) %>% 
  liana_dotplot(target_groups = "Malignant Cell",
                source_groups = c("B Cell", "Other", "Endothelial Cell", "Epithelial Cell", "Fibroblast", "Monocyte/Macrophage", "Specialized Immune Cell", "T Cell"),
                specificity = "pvalue",
                magnitude = "lr.mean",
                show_complex = TRUE,
                size.label = "-log10(p-value)",
                ntop =10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 12),
        plot.caption.position = "plot",
        plot.subtitle = element_text(hjust = 0.5, size = 23),
        plot.title = element_text(face = "bold", size = 30, color = merck_colors[1])) +
  labs(caption = "Dataset: IMPACT_lung_cancers", title = "Top 10 CellPhoneDB interactions (pvalue <= 0.05)" , x = "Target: Malignant", subtitle = "Source")

ggsave("11_CellPhoneDB_filtered_interactions.png", plot = p4, scale = 2, create.dir = T, path = "./res/", height = 8, width = 8)
