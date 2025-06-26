# # Install packages ####
# install.packages("renv")
renv::init() #Select 1
install.packages("./packages/OmnipathR-master_april_2024.tar.gz", repos = NULL, type = "source")
install.packages("./packages/liana-master_april_2024.tar.gz", repos = NULL, type = "source")
renv::snapshot(type = "all", exclude = c("liana", "OmnipathR"))
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

# Add Merck Colors ####
merck_colors <- c(violet = "#503291", blue = "#0F69AF", green = "#149B5F", red = "#E61E50", 
                  pink = "#EB3C96", lightblue = "#2DBECD", lightgreen = "#A5CD50",
                  yellow = "#FFC832", palepink = "#E1C3CD", paleblue = "#96D7D2",
                  palegreen = "#B4DC96", paleyellow = "#FFDCB9", grey = "#999999")


LUAD <- qread("./seurat_files/2022_Zhu_GSE189487_LUAD_ST1_deconv.qs")

# cluster0 <- LUAD@meta.data[LUAD@meta.data$seurat_clusters == 0, ]
# mean(cluster0$Malignant)

# Calculate mean x_percentage for each cluster
# clusters <- LUAD@meta.data %>%
#   group_by(seurat_clusters) %>%
#   summarise(Malignant = mean(Malignant))


######## LIANA: ALL METHODS ########
# This is a very long step! Especially if you have many cells!  ####
liana_data <- liana_wrap(LUAD,
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
if (!dir.exists('./liana_files/')) {dir.create("./liana_files/")}
qsave(liana_data,"./liana_files/240829_LUAD_ST1_liana_data.qs")
liana_data <- qread("./liana_files/240829_LUAD_ST1_liana_data.qs")
# 
# qsave(liana_data_aggregated,"./seurat_files/240805_pseudo_lung_liana_data_aggregated.qs")
qsave(liana_data_aggregated, "./liana_files/240829_LUAD_ST1_liana_data_aggregated.qs")
liana_data_aggregated <- qread("./liana_files/240829_LUAD_ST1_liana_data_aggregated.qs")

######## Visualisation ########
# Simple Dot Plot  ####
p1 <- liana_data_aggregated %>%
  #  filter(pvalue <= 0.05) %>%
  liana_dotplot(source_groups = "1",
                #target_groups = cluster_ids,
                ntop = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 12),
        plot.caption.position = "plot",
        plot.subtitle = element_text(hjust = 0.5, size = 23),
        plot.title = element_text(face = "bold", size = 30, color = merck_colors[1]))  +
  labs(caption = "Dataset: LUAD_ST1", title = "Top 20 unfiltered interactions", subtitle = "Source: Malignant")

ggsave("p1.pdf", plot = p1, create.dir = T, path = "./res/", height = 13, width = 12)


# Simple Dot Plot  ####
p2 <- liana_data_aggregated %>%
  liana_dotplot(target_groups = "1",
                #target_groups = cluster_ids,
                ntop = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text = element_text(size = 12),
        plot.caption.position = "plot",
        plot.subtitle = element_text(hjust = 0.5, size = 23),
        plot.title = element_text(face = "bold", size = 30, color = merck_colors[1])) +
  labs(caption = "Dataset: LUAD_ST1", title = "Top 20 unfiltered interactions", x = "Target: Malignant", subtitle = "Source")

ggsave("p2.png", plot = p2, create.dir = T, path = "./res/", height = 13, width = 12)


# Frequency Heatmap  ####
liana_trunc <- liana_data_aggregated %>% # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.05) # this can be FDR-corr if n is too high

heat_freq(liana_trunc) +
  
  # Frequency chord plot  ####
chord_freq(liana_trunc,
           #source_groups = c("Endothelial cell of vascular tree"),
           #target_groups = cluster_ids
)



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
  liana_dotplot(source_groups = "1",
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
  labs(caption = "Dataset: LUAD_ST1", title = "Top 10 CellPhoneDB interactions (pvalue <= 0.05)", subtitle = "Source: Malignant")

ggsave("p3.pdf", plot = p3, create.dir = T, path = "./res/", height = 13, width = 12)



p4 <- liana_interactions %>%
  # keep only the interactions of interest
  inner_join(liana_data$cellphonedb, 
             by = c("ligand.complex", "receptor.complex")) %>%
  # invert size (low p-value/high specificity = larger dot size)
  # + add a small value to avoid Infinity for 0s
  mutate(pvalue = -log10(pvalue + 1e-10)) %>% 
  liana_dotplot(target_groups = "1",
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
  labs(caption = "Dataset: LUAD_ST1", title = "Top 10 CellPhoneDB interactions (pvalue <= 0.05)" , x = "Target: Malignant", subtitle = "Source")

ggsave("p4.png", plot = p4, scale = 2, create.dir = T, path = "./res/", height = 13, width = 12)

################################################################
# 
# ######## Visualisation ########
# 
# # Simple Dot Plot  ####
# liana_data_aggregated %>%
#   #  filter(pvalue <= 0.05) %>%
#   liana_dotplot(source_groups = "Malignant cell",
#                 #target_groups = cluster_ids,
#                 ntop = 20) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         strip.text = element_text(size = 12)) #+
# #scale_x_discrete(labels = c("Endothelial cell of vascular tree" = "Endothelial cell \nof vascular tree", "Effector memory CD8+ alpha-beta T cell" = "Effector memory CD8+ \nalpha-beta T cell"))
# 
# # Simple Dot Plot  ####
# liana_data_aggregated %>%
#   liana_dotplot(target_groups = "Malignant cell",
#                 #target_groups = cluster_ids,
#                 ntop = 20) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         strip.text = element_text(size = 12)) #+
# #scale_x_discrete(labels = c("Endothelial cell of vascular tree" = "Endothelial cell \nof vascular tree", "Effector memory CD8+ alpha-beta T cell" = "Effector memory CD8+ \nalpha-beta T cell"))
# 
# 
# # Frequency Heatmap  ####
# liana_trunc <- liana_data_aggregated %>% # only keep interactions concordant between methods
#   filter(aggregate_rank <= 0.05) # this can be FDR-corr if n is too high
# 
# heat_freq(liana_trunc)
# 
# # Frequency chord plot  ####
# chord_freq(liana_trunc,
#            #source_groups = c("Endothelial cell of vascular tree"),
#            #target_groups = cluster_ids
# )
# 
# 
# 
# ###### CellPhoneDB only ######
# # identify interactions of interest ####
# liana_interactions <- liana_data$cellphonedb %>%
#   # only keep interactions with p-val <= 0.05
#   filter(pvalue <= 0.05) %>% # this reflects interactions `specificity`
#   # then rank according to `magnitude` or 'specificity'
#   rank_method(method_name = "cellphonedb",
#               mode = "magnitude") %>%
#   # keep top 20 interactions (regardless of cell type)
#   distinct_at(c("ligand.complex", "receptor.complex")) %>%
#   head(20)
# 
# liana_interactions %>%
#   # keep only the interactions of interest
#   inner_join(liana_data$cellphonedb, 
#              by = c("ligand.complex", "receptor.complex")) %>%
#   # invert size (low p-value/high specificity = larger dot size)
#   # + add a small value to avoid Infinity for 0s
#   mutate(pvalue = -log10(pvalue + 1e-10)) %>% 
#   liana_dotplot(source_groups = "1",
#                 specificity = "pvalue",
#                 magnitude = "lr.mean",
#                 show_complex = TRUE,
#                 size.label = "-log10(p-value)",
#                 ntop =10) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#         strip.text = element_text(size = 12),
#         plot.caption.position = "plot",
#         plot.subtitle = element_text(hjust = 0.5, size = 23),
#         plot.title = element_text(face = "bold", size = 30, color = merck_colors[1])) +
#   labs(caption = "Dataset: LUAD_ST1", title = "Top 10 CellPhoneDB interactions (pvalue <= 0.05)", subtitle = "Source: Malignant")
# 
# liana_interactions %>%
#   # keep only the interactions of interest
#   inner_join(liana_data$cellphonedb, 
#              by = c("ligand.complex", "receptor.complex")) %>%
#   # invert size (low p-value/high specificity = larger dot size)
#   # + add a small value to avoid Infinity for 0s
#   mutate(pvalue = -log10(pvalue + 1e-10)) %>% 
#   liana_dotplot(target_groups = "Malignant cell",
#                 specificity = "pvalue",
#                 magnitude = "lr.mean",
#                 show_complex = TRUE,
#                 size.label = "-log10(p-value)",
#                 ntop =10)
# 
