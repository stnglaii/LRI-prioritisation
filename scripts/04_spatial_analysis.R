if (run_spatial_analysis == T) {
  setwd("~/lri-prioritization")
  
  # list generated liana result files
  sp_healthy_samples <- list.files(path = "./liana_results/spatial/healthy/")
  sp_cancer_samples <- list.files(path = "./liana_results/spatial/cancer/")
  
  # create empty dataframes for data storage
  sp_cancer_results <- data.frame()
  sp_healthy_results <- data.frame()
  
  # read liana results, add file name and condition, concatenate tables
  setwd("~/lri-prioritization/liana_results/spatial/healthy")
  for (file in sp_healthy_samples) {
    file_name <- tools::file_path_sans_ext(basename(file))
    temp <- read.csv(file)
    temp <- mutate(temp, sample = file_name, condition = "healthy")
    assign(file_name, temp)
    sp_healthy_results <- rbind(sp_healthy_results, temp)
  }
  
  setwd("~/lri-prioritization/liana_results/spatial/cancer")
  for (file in sp_cancer_samples) {
    file_name <- tools::file_path_sans_ext(basename(file))
    temp <- read.csv(file)
    temp <- mutate(temp, sample = file_name, condition = "cancer")
    assign(file_name, temp)
    sp_cancer_results <- rbind(sp_cancer_results, temp)
  }
  
  # combine both conditions into one table
  sp_all_lri <- rbind(sp_cancer_results, sp_healthy_results)
  
  # delete file extension from file names
  setwd("~/lri-prioritization")
  sp_healthy_samples <- tools::file_path_sans_ext(basename(sp_healthy_samples))
  sp_cancer_samples <- tools::file_path_sans_ext(basename(sp_cancer_samples))
  
  # calculate scores and process table for limma
  sp_lri_limma <- sp_all_lri %>%
    filter(morans > 0, mean > 0) %>%         # filter for positively expressed & co-localized LRIs
    group_by(interaction, condition) %>%  
    filter(n() >= 2) %>%                     # filter for LRIs that occur at least twice per condition
    ungroup() %>% 
    group_by(interaction) %>%
    filter(all(c("cancer", "healthy") %in% condition)) %>%  # filter for LRIs that occur in both conditions
    ungroup() %>%
    mutate(norm_mean = rescale(mean),        # calculate score
           norm_morans = rescale(morans),
           CV = std/mean,
           penalty = rescale(CV),
           score = rescale(norm_morans*norm_mean*penalty)) %>% 
    select(interaction, sample, score) %>%   # choose relevant columns
    spread(sample, score) %>%                # convert table to wide format
    as.data.frame()                          # convert from tibble to dataframe
  
  rownames(sp_lri_limma) <- sp_lri_limma$interaction        # set interaction as row names
  sp_lri_limma <- sp_lri_limma[,-1]                         # delete interaction column (here: first column)
  
  
  # create design matrix for limma
  sample_info <- unique(sp_all_lri[, c("sample", "condition")])
  condition <- factor(sample_info$condition)
  design <- model.matrix(~0 + condition)  # Use model without intercept
  colnames(design) <- levels(condition)
  
  # fit linear model to processed liana results with limma
  fit <- lmFit(sp_lri_limma, design)
  
  # create contrast matrix for limma
  contrast_matrix <- makeContrasts(cancer_vs_healthy = cancer - healthy, levels=design)
  
  # fit contrasts
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # examine top differentially expressed interactions
  sp_limma_results <- topTable(fit2, coef = "cancer_vs_healthy", number = Inf, adjust.method = "BH")
  sp_limma_results$interaction <- row.names(sp_limma_results)
  sp_limma_results <- sp_limma_results %>% 
    arrange(desc(logFC))                          # list all LRIs with resulting log2FCs in descending order
  
  # filter and rank those LRIs that have a positive log2FC for cancer-specific LRIs
  sp_cancer_ranking <- sp_limma_results %>% 
    arrange(desc(logFC)) %>% 
    filter(logFC > 0) %>%      
    mutate(rank = rank(-logFC, ties.method = "min"))
  
  # filter and rank those LRIs that have a positive log2FC for healthy-specific LRIs
  sp_healthy_ranking <- sp_limma_results %>% 
    arrange(logFC) %>% 
    filter(logFC < 0) %>% 
    mutate(rank = rank(logFC, ties.method = "min"))
  
  # save spatial ranking results
  write.csv(sp_lri_limma, paste0("./out/", tissue, "/spatial_limma_table.csv"), row.names=F, quote=F)
  write.csv(sp_limma_results, paste0("./out/", tissue, "/spatial_limma_results.csv"), row.names=F, quote=F)
  write.csv(sp_cancer_ranking, paste0("./out/", tissue, "/spatial_cancer_ranking.csv"), row.names=F, quote=F)
} else {
  sp_lri_limma <- read_csv(paste0("./out/", tissue, "/spatial_limma_table.csv"))
  sp_limma_results <- read_csv(paste0("./out/", tissue, "/sp_limma_results.csv"))
  sp_cancer_ranking <- read_csv(paste0("./out/", tissue, "/sp_cancer_ranking.csv"))
}



#### Plots ####
if (create_sp_plots == T) {
  # bar plot of top 50 interactions
  barplot_spatial_cancer <- ggplot(sp_limma_results[50:1,], aes(y = factor(sp_limma_results$interaction[50:1], levels = sp_limma_results$interaction[50:1]), x = logFC)) +
    geom_bar(stat = "identity", fill = "#503291") +
    labs(y = "Interaction", 
         x = "log2FC") +
    xlim(c(0,0.3))+
    #guides(y = guide_axis(angle = 45)) +
    theme_classic()+                                                                                        
    geom_text(aes(label = format(logFC, digits = 2)), size = 3, hjust = -0.5) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  
  barplot_spatial_healthy <- ggplot(sp_healthy_ranking[1:50,], aes(y = factor(sp_healthy_ranking$interaction[1:50], levels = sp_healthy_ranking$interaction[50:1]), x = logFC)) +
    geom_bar(stat = "identity", fill = "#503291") +
    labs(y = "Interaction", 
         x = "log2FC") +
    theme_classic() +   
    xlim(c(-0.15,0))+
    geom_text(aes(label = format(logFC, digits = 2)), size = 3, hjust = 1.5) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  
  # scatterplot of all interactions
  cancer_cols <- sp_lri_limma %>%
    select(contains("cancer"))
  cancer_means <- apply(cancer_cols, 1, mean, na.rm = TRUE)
  cancer_medians <- apply(cancer_cols, 1, median, na.rm = TRUE)
  sp_cancer <- cbind(cancer_cols, cancer_means = cancer_means, cancer_medians = cancer_medians)
  
  healthy_cols <- sp_lri_limma %>%
    select(contains("healthy"))
  healthy_means <- apply(healthy_cols, 1, mean, na.rm = TRUE)
  healthy_medians <- apply(healthy_cols, 1, median, na.rm = TRUE)
  sp_healthy <- cbind(healthy_cols, healthy_means = healthy_means, healthy_medians = healthy_medians)
  
  sp_merged <- merge(sp_cancer, sp_healthy, by = 0)
  rownames(sp_merged) <- sp_merged$Row.names
  names(sp_merged)[names(sp_merged) == "Row.names"] <- "interaction"
  
  sp_merged <- sp_merged %>%
    mutate(mean_score_cancer = rescale((cancer_means)),
           mean_score_healthy = rescale((healthy_means)),
           log2FC_means = log2(mean_score_cancer/mean_score_healthy),
           log2FC_medians = log2(mean_score_cancer/mean_score_healthy)) %>%
    arrange(desc(log2FC_means)) %>% 
    select(mean_score_healthy, mean_score_cancer, log2FC_means, interaction)
  
  sp_merged$color_category <- cut(sp_merged$log2FC_means,
                                  breaks = c(-Inf, -0.2, 0.2, Inf),
                                  labels = c("healthy-specific", "condition-unspecific", "cancer-specific"))
  
  
  scatterplot_spatial <- ggplot(sp_merged, aes(x = mean_score_healthy, y = mean_score_cancer, label = interaction, color = color_category)) +
    geom_point(size = 0.7) +
    scale_color_manual(values = c("healthy-specific" = "#0F69AF", "condition-unspecific" = "gray", "cancer-specific" = "#E61E50")) +  # Specify colors
    theme_classic()+
    coord_fixed(clip = "off") +
    geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed", color = 'gray') +
    geom_label_repel(data = sp_merged[which((sp_merged$log2FC_means)>(1) & (sp_merged$mean_score_cancer)>(0.6)),],  size = 3, min.segment.length = unit(0, 'lines'), colour = "#E61E50", max.overlaps = 6, xlim = c(-0.1,0.5)) +
    labs(color = "Tissue Specificity") 
  
}

if (display_sp_plots == T) {
  print(barplot_spatial_cancer)
  print(barplot_spatial_healthy)
  print(scatterplot_spatial)
}

if (save_sp_plots == T) {
  ggsave(paste0("./out/", tissue, "/barplot_spatial_cancer.png"), barplot_spatial_cancer)
  ggsave(paste0("./out/", tissue, "/barplot_spatial_healthy.png"), barplot_spatial_healthy)
  ggsave(paste0("./out/", tissue, "/scatterplot_spatial.png"), scatterplot_spatial)
}