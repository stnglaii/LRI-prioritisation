library(readr)
library(tidyverse)
library(limma)

sp_healthy_samples <- setdiff(list.files(path = "~/lri-prioritization/liana_results/spatial/normal/"), list.dirs(path = "~/lri-prioritization/liana_results/spatial/normal", recursive = FALSE, full.names = FALSE))
sp_cancer_samples <- setdiff(list.files(path = "~/lri-prioritization/liana_results/spatial/cancer/"), list.dirs(path= "~/lri-prioritization/liana_results/spatial/cancer", recursive = FALSE, full.names = FALSE))

sp_all_lri <- data.frame()
sp_cancer_results <- data.frame()
sp_healthy_results <- data.frame()

setwd("~/lri-prioritization/liana_results/spatial/normal")
for (file in sp_healthy_samples) {
  file_name <- tools::file_path_sans_ext(basename(file))
  temp <- read.csv(file)
  temp <- mutate(temp, sample = file_name, condition = "Healthy")
  assign(file_name, temp)
  sp_healthy_results <- rbind(sp_healthy_results, temp)
}

setwd("~/lri-prioritization/liana_results/spatial/cancer")
for (file in sp_cancer_samples) {
  file_name <- tools::file_path_sans_ext(basename(file))
  temp <- read.csv(file)
  temp <- mutate(temp, sample = file_name, condition = "Cancer")
  assign(file_name, temp)
  sp_cancer_results <- rbind(sp_cancer_results, temp)
}

sp_all_lri <- rbind(sp_cancer_results, sp_healthy_results)

sp_healthy_samples <- tools::file_path_sans_ext(basename(sp_healthy_samples))
sp_cancer_samples <- tools::file_path_sans_ext(basename(sp_cancer_samples))

sp_lri_limma <- sp_all_lri %>%
  filter(morans > 0,
         #morans_pvals < 0.05,
         mean > 0) %>% 
  #group_by(interaction) %>%
  #filter(n_distinct(condition) == 2) %>%
  #ungroup() %>%
  select(sample, interaction, morans) %>%
  spread(sample, morans) %>%
  # rowwise() %>%
  # filter(
  #   sum(!is.na(c_across(all_of(sp_cancer_samples)))) > 4,
  #   sum(!is.na(c_across(all_of(sp_healthy_samples)))) > 4) %>%
  # mutate(
  #   cancer_median = median(c_across(all_of(tools::file_path_sans_ext(basename(sp_cancer_samples)))), na.rm = TRUE),
  #   healthy_median = median(c_across(all_of(tools::file_path_sans_ext(basename(sp_healthy_samples)))), na.rm = TRUE),
  #   # Overwrite cancer columns
  #   across(
  #     all_of(tools::file_path_sans_ext(basename(sp_cancer_samples))),
  #     ~ ifelse(is.na(.x), cancer_median, .x)
#   ),
#   # Overwrite healthy columns
#   across(
#     all_of(tools::file_path_sans_ext(basename(sp_healthy_samples))),
#     ~ ifelse(is.na(.x), healthy_median, .x)
#   )) %>%
# ungroup() %>%
# select(-cancer_median, -healthy_median) %>%
as.data.frame()

rownames(sp_lri_limma) <- sp_lri_limma$interaction
sp_lri_limma <- sp_lri_limma[,-1]

sample_info <- unique(sp_all_lri[, c("sample", "condition")])
condition <- factor(sample_info$condition)
design <- model.matrix(~0 + condition)  # Use model without intercept
colnames(design) <- levels(condition)

design

fit <- lmFit(sp_lri_limma, design)

contrast_matrix <- makeContrasts(Cancer_vs_Healthy = Cancer - Healthy, levels=design)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Examine top differentially expressed interactions
sp_limma_results <- topTable(fit2, coef = "Cancer_vs_Healthy", number = Inf, adjust.method = "BH")
sp_limma_results$interaction <- row.names(sp_limma_results)
sp_limma_results <- sp_limma_results %>% 
  arrange(desc(logFC))
View(sp_limma_results)

sp_cancer_ranking <- sp_limma_results %>% 
  arrange(desc(logFC)) %>% 
  filter(logFC > 0) %>% 
  mutate(rank = rank(-logFC, ties.method = "min"))

sp_cancer_standalone <- sp_lri_limma %>%
  filter(rowSums(is.na(select(., contains("normal")))) == ncol(select(., contains("normal"))))

sp_healthy_ranking <- sp_limma_results %>% 
  arrange(logFC) %>% 
  filter(logFC < 0) %>% 
  mutate(rank = rank(logFC, ties.method = "min"))

sp_healthy_standalone <- sp_lri_limma %>%
  filter(rowSums(is.na(select(., contains("cancer")))) == ncol(select(., contains("cancer"))))


######### Visualisations ##########
# Scatterplot
cancer_cols <- sp_lri_limma %>%
  select(contains("cancer"))
cancer_means <- apply(cancer_cols, 1, mean, na.rm = TRUE)
cancer_medians <- apply(cancer_cols, 1, median, na.rm = TRUE)
sp_cancer <- cbind(cancer_cols, cancer_means = cancer_means, cancer_medians = cancer_medians)

healthy_cols <- sp_lri_limma %>%
  select(contains("normal"))
healthy_means <- apply(healthy_cols, 1, mean, na.rm = TRUE)
healthy_medians <- apply(healthy_cols, 1, median, na.rm = TRUE)
sp_healthy <- cbind(healthy_cols, healthy_means = healthy_means, healthy_medians = healthy_medians)

sp_merged <- merge(sp_cancer, sp_healthy, by = 0)
rownames(sp_merged) <- sp_merged$Row.names
names(sp_merged)[names(sp_merged) == "Row.names"] <- "interaction"

sp_merged <- sp_merged %>%
  mutate(log2FC_means = log2(cancer_means/healthy_means),
         log2FC_medians = log2(cancer_medians/healthy_medians)) %>%
  arrange(desc(log2FC_means))

scatterplot_spatial <- ggplot(sp_merged, aes(x = healthy_means, y = cancer_means, label = interaction)) +
  geom_point(data = sp_merged[which((sp_merged$log2FC_means)<=(-1)),],color = 'blue', size = 0.7)+
  geom_point(data = sp_merged[which((sp_merged$log2FC_means)>(-1) & (sp_merged$log2FC_means)<1),],color = 'gray', size = .7)+
  geom_point(data = sp_merged[which((sp_merged$log2FC_means)>1),],color = 'red', size = 0.7)+
  theme_classic()+
  geom_label_repel(data = sp_merged[which((sp_merged$log2FC_means)>(5.5)),],  size = 3, min.segment.length = unit(0, 'lines'))+
  ggtitle('Ligand-Receptor analysis: Lung cancer cohort from Spatial Transcriptomics data',
          subtitle = 'Ranked by global Moran\'s R\nMalignant/Epithelial cells <-> Immune cells \nred: Ligand-Receptor pairs specific to tumors')

scatterplot_spatial



barplot_pseudobulk_cancer <- ggplot(sp_cancer_ranking[1:50,], aes(y = factor(sp_cancer_ranking$interaction[1:50], levels = sp_cancer_ranking$interaction[50:1]), x = logFC)) +
  geom_bar(stat = "identity") +
  labs(title = "Top 50 LRI in cancerous tissue",
       subtitle = "Ranked by global Moran's R",
       y = "Interaction Pairs", 
       x = "logFC") +
  #guides(y = guide_axis(angle = 45)) +
  xlim(c(0,0.7))+
  theme_classic() +                                                                                        
  geom_text(aes(label = format(logFC, digits = 3)), size = 3, hjust = -0.5) +
  theme(plot.margin = unit(c(0,1,0,0), "cm"))

barplot_pseudobulk_healthy <- ggplot(sp_healthy_ranking[1:50,], aes(y = factor(sp_healthy_ranking$interaction[1:50], levels = sp_healthy_ranking$interaction[50:1]), x = logFC)) +
  geom_bar(stat = "identity") +
  labs(title = "Top 50 LRI in healthy tissue",
       subtitle = "Ranked by global Moran's R",
       y = "Interaction Pairs", 
       x = "logFC") +
  #guides(y = guide_axis(angle = 45)) +
  theme_classic() +   
  xlim(c(-0.4,0))+
  geom_text(aes(label = format(logFC, digits = 3)), size = 3, hjust = 1.5) +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))

barplot_pseudobulk_cancer+barplot_pseudobulk_healthy


# 
# #### manual calculation ####
# cancer_cols <- sp_lri_limma %>%
#   select(contains("cancer"))
# cancer_means <- apply(cancer_cols, 1, mean, na.rm = TRUE)
# cancer_medians <- apply(cancer_cols, 1, median, na.rm = TRUE)
# sp_cancer <- cbind(cancer_cols, cancer_means = cancer_means, cancer_medians = cancer_medians)
# 
# healthy_cols <- sp_lri_limma %>%
#   select(contains("normal"))
# healthy_means <- apply(healthy_cols, 1, mean, na.rm = TRUE)
# healthy_medians <- apply(healthy_cols, 1, median, na.rm = TRUE)
# sp_healthy <- cbind(healthy_cols, healthy_means = healthy_means, healthy_medians = healthy_medians)
# 
# sp_merged <- merge(sp_cancer, sp_healthy, by = 0)
# rownames(sp_merged) <- sp_merged$Row.names
# names(sp_merged)[names(sp_merged) == "Row.names"] <- "interaction"
# 
# sp_merged <- sp_merged %>%
#   mutate(log2FC_means = log2(cancer_means/healthy_means),
#          log2FC_medians = log2(cancer_medians/healthy_medians)) %>%
#   arrange(desc(log2FC_means)) #%>%
#   #select(cancer_means, healthy_means, log2FC_means)
# 
# p_value <- sp_all_lri %>% 
#   filter(morans > 0,
#          morans_pvals < 0.05,
#          mean > 0) %>% 
#   group_by(interaction) %>% #  necessary because we'll perform statistical tests on each gene separately
#   filter(n_distinct(condition) == 2) %>%
#   summarise(pvalue = suppressWarnings(wilcox.test(mean ~ condition)$p.value), .groups = "drop") %>%
#   #ungroup() %>%
#   mutate(padjust = p.adjust(pvalue, method = "fdr")) %>% 
#   as.data.frame()
# 
# rownames(p_value) <- p_value$interaction
# p_value <- p_value[,-1]
# 
# sp_merged <- merge(sp_merged, p_value, by = 0, all.x = T)
# 
# 

