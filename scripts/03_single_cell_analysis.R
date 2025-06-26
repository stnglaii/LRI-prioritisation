if (run_download_IMPACT == T) {
  api$data_lake$download_file("file-d7bb06f9-ae20-41ec-86cc-1526feb0f3e4", path = "./seurat_files/", overwrite = TRUE, show_progress = T)
  
  impact <- qread("./seurat_files/Seurat_IMPACT_Atlas_Q1_2025.qs") #28372 cells
  impact <- impact[,which(impact$sample_anatomical_entity_efo_name == tissue)] # 3822 cells left
  cancer <- impact[,which(impact$sample_status == 'tumor' | impact$sample_status == 'metastasis')] #2172 cells
  cancer <- SetIdent(cancer, value = cancer@meta.data[["cell_type_general"]])
  current_cell_types <- levels(cancer@active.ident)
  new_cell_types <- setNames(c("Immune Cell", "Immune Cell", "Immune Cell", "Immune Cell",
                               "Immune Cell", "Endothelial Cell", "Epithelial Cell", "Malignant Cell",
                               "Smooth Muscle Cell", "Fibroblast", "Progenitor Cell", "Blood Cell",
                               "Neuron", "Glial Cell"), current_cell_types)
  cancer <- RenameIdents(cancer, new_cell_types)
  healthy <- impact[,which(impact$sample_status == 'healthy')]

  if (!dir.exists(paste0("./seurat_files/", tissue))) {dir.create(paste0("./seurat_files/", tissue))}
  
  qsave(cancer, paste0("./seurat_files/", tissue, "/cancers_Q1_2025.qs"))
  qsave(healthy, paste0("./seurat_files/", tissue, "/healthy_Q1_2025.qs"))
  rm(impact)
} else {
  cancer <- qread(paste0("./seurat_files/", tissue, "/cancers_Q1_2025.qs"))
  healthy <- qread(paste0("./seurat_files/", tissue, "/healthy_Q1_2025.qs"))
}

#### run analysis for single-cell LIANA + results?
run_sc_analysis <- F

if (run_sc_analysis == T) {
  # read liana results
  sc_cancer_results <- read_csv("./liana_results/single_cell/lung_cancers_cell_type_general_Q1_2025.csv")
  sc_healthy_results <- read_csv("./liana_results/single_cell/lung_healthy_cell_type_general_Q1_2025.csv")
  
  # select relevant columns, add interaction column, add integer ranks to specificity and magnitude scores
  sc_cancer_results <- sc_cancer_results %>%
    select(source, target, ligand_complex, receptor_complex, specificity_rank, magnitude_rank) %>%
    #filter(specificity_rank <0.7,
    #       magnitude_rank <0.7) %>% 
    mutate(interaction = paste0(source, "^", target, "_", ligand_complex, "^", receptor_complex ),
           specificity_score = rescale(1-specificity_rank),
           magnitude_score = rescale(1-magnitude_rank),
           #aggregate_score = rescale(sqrt(specificity_score*magnitude_score)),
           #aggregate_rank = as.integer(rank(-aggregate_score))
           ) %>%
    relocate(interaction, .after = receptor_complex)
  
  sc_healthy_results <- sc_healthy_results %>%
    #select(source, target, ligand_complex, receptor_complex, specificity_rank, magnitude_rank) %>%
    #filter(specificity_rank <0.7,
    #       magnitude_rank <0.7) %>% 
    mutate(interaction = paste0(source, "^", target, "_", ligand_complex, "^", receptor_complex ),
           specificity_score = rescale(1-specificity_rank),
           magnitude_score = rescale(1-magnitude_rank),
           #aggregate_score = rescale(sqrt(specificity_score*magnitude_score)),
           #aggregate_rank = as.integer(rank(-aggregate_score))
           ) %>%
    relocate(interaction, .after = receptor_complex)
  
  # rename cell types for analysis
  sc_cancer_results <- sc_cancer_results %>%
    mutate(
      source = recode(source,
                      "B cell" = "Immune",
                      "CD4+ T cell" = "Immune",
                      "CD8+ T cell" = "Immune",
                      "Myeloid cell" = "Immune",
                      "T cell" = "Immune",
                      "NK cell" = "Immune",
                      'Malignant cell' = 'Mal/Epi'),
      target = recode(target,
                      "B cell" = "Immune",
                      "CD4+ T cell" = "Immune",
                      "CD8+ T cell" = "Immune",
                      "Myeloid cell" = "Immune",
                      "T cell" = "Immune",
                      "NK cell" = "Immune",
                      'Malignant cell' = 'Mal/Epi'),
      interaction = paste0(source, "^", target, "_", ligand_complex, "^", receptor_complex )
    )
  
  sc_healthy_results <- sc_healthy_results %>%
    mutate(
      source = recode(source,
                      "B cell" = "Immune",
                      "CD4+ T cell" = "Immune",
                      "CD8+ T cell" = "Immune",
                      "Myeloid cell" = "Immune",
                      "T cell" = "Immune",
                      "NK cell" = "Immune",
                      'Epithelial cell' = 'Mal/Epi'),
      target = recode(target,
                      "B cell" = "Immune",
                      "CD4+ T cell" = "Immune",
                      "CD8+ T cell" = "Immune",
                      "Myeloid cell" = "Immune",
                      "T cell" = "Immune",
                      "NK cell" = "Immune",
                      'Epithelial cell' = 'Mal/Epi'),
      interaction = paste0(source, "^", target, "_", ligand_complex, "^", receptor_complex )
    )
  
  # preprocess dataframe:
    # calculate mean specificity and magnitude scores by interaction
    # create necessary columns for subsequent steps
    # filter out any interactions that aren't between malignant and immune cells
  sc_cancer_results <- sc_cancer_results %>% 
    group_by(interaction) %>% 
    summarise(specificity_score = mean(specificity_score),
              magnitude_score = mean(magnitude_score)) %>% 
    mutate(aggregate_score = rescale(sqrt(specificity_score*magnitude_score)),
           source = str_extract(interaction, "^[^\\^]+"),
           target = str_extract(interaction, "(?<=\\^)[^_]+"),
           ligand = str_extract(interaction, "(?<=_)[^\\^]+"),
           ligand_expr = NA,
           ligand_log2FC_min = NA,
           ligand_log2FC_max = NA,
           receptor = str_extract(interaction, "(?<=\\^)[^\\^]+$"),
           receptor_expr = NA,
           receptor_log2FC_min = NA,
           receptor_log2FC_max = NA) %>% 
    relocate(interaction, .after = receptor_log2FC_max) %>%
    relocate(c(specificity_score, magnitude_score, aggregate_score), .after = interaction,) %>% 
    filter(!grepl("^[0-9]", interaction), # Filter out numeric rows
           aggregate_score != 0,
           source %in% c("Mal/Epi", "Immune"),
           target %in% c("Mal/Epi", "Immune"),
           source != target) %>% 
    arrange(desc(aggregate_score))
  rownames(sc_cancer_results) <- NULL
  
  sc_healthy_results <- sc_healthy_results %>% 
    group_by(interaction) %>% 
    summarise(specificity_score = mean(specificity_score),
              magnitude_score = mean(magnitude_score)) %>% 
    mutate(aggregate_score = rescale(sqrt(specificity_score*magnitude_score)),
           source = str_extract(interaction, "^[^\\^]+"),
           target = str_extract(interaction, "(?<=\\^)[^_]+"),
           ligand = str_extract(interaction, "(?<=_)[^\\^]+"),
           receptor = str_extract(interaction, "(?<=\\^)[^\\^]+$")) %>% 
    relocate(interaction, .after = receptor) %>%
    relocate(c(specificity_score, magnitude_score, aggregate_score), .after = interaction,) %>% 
    filter(!grepl("^[0-9]", interaction), # Filter out numeric rows
           aggregate_score != 0,
           source %in% c("Mal/Epi", "Immune"),
           target %in% c("Mal/Epi", "Immune"),
           source != target)
  rownames(sc_healthy_results) <- NULL
  
  # fetch expression data from seurat object
  for (row in 1:nrow(sc_cancer_results)) {
    ligand <- sc_cancer_results$ligand[row]
    receptor <- sc_cancer_results$receptor[row]
    if (grepl('_', ligand, fixed=TRUE) & grepl('_', receptor, fixed=TRUE)) {
      split_strings <- strsplit(ligand, "_")[[1]]
      ligand1 <- split_strings[1]
      ligand2 <- split_strings[2]
      split_strings <- strsplit(receptor, "_")[[1]]
      receptor1 <- split_strings[1]
      receptor2 <- split_strings[2]
      result <- FetchData(cancer, vars = c("ident", ligand1, ligand2, receptor1, receptor2)) %>% 
        group_by(ident) %>% 
        summarise(
          mean_ligand1 = mean(get(ligand1)),
          mean_ligand2 = mean(get(ligand2)),
          mean_receptor1 = mean(get(receptor1)),
          mean_receptor2 = mean(get(receptor2))
        ) %>%
        mutate(
          mean_ligand = pmin(mean_ligand1, mean_ligand2),
          comp1 = ifelse(mean_ligand1 <= mean_ligand2, ligand1, ligand2),
          mean_receptor = pmin(mean_receptor1, mean_receptor2),
          comp2 = ifelse(mean_receptor1 <= mean_receptor2, receptor1, receptor2)
        ) %>%
        select(-mean_receptor1, -mean_receptor2)
    } else  if (!grepl('_', ligand, fixed=TRUE) & !grepl('_', receptor, fixed=TRUE)) {
      result <- FetchData(cancer, vars = c("ident", ligand, receptor)) %>% 
        group_by(ident) %>% 
        summarise(mean_ligand = mean(get(ligand)),
                  mean_receptor = mean(get(receptor)))
    } else if (grepl('_', ligand, fixed=TRUE)) {
      split_strings <- strsplit(ligand, "_")[[1]]
      ligand1 <- split_strings[1]
      ligand2 <- split_strings[2]
      result <- FetchData(cancer, vars = c("ident", ligand1, ligand2, receptor)) %>% 
        group_by(ident) %>% 
        summarise(
          mean_ligand1 = mean(get(ligand1)),
          mean_ligand2 = mean(get(ligand2)),
          mean_receptor = mean(get(receptor))
        ) %>%
        mutate(
          mean_ligand = pmin(mean_ligand1, mean_ligand2),
          comp = ifelse(mean_ligand1 <= mean_ligand2, ligand1, ligand2)
        ) %>%
        select(-mean_ligand1, -mean_ligand2, -comp)
    } else if (grepl('_', receptor, fixed=TRUE)) {
      split_strings <- strsplit(receptor, "_")[[1]]
      receptor1 <- split_strings[1]
      receptor2 <- split_strings[2]
      result <- FetchData(cancer, vars = c("ident", ligand, receptor1, receptor2)) %>% 
        group_by(ident) %>% 
        summarise(
          mean_ligand = mean(get(ligand)),
          mean_receptor1 = mean(get(receptor1)),
          mean_receptor2 = mean(get(receptor2))
        ) %>%
        mutate(
          mean_receptor = pmin(mean_receptor1, mean_receptor2),
          comp = ifelse(mean_receptor1 <= mean_receptor2, receptor1, receptor2)
        ) %>%
        select(-mean_receptor1, -mean_receptor2, -comp)
    }
    
    if (sc_cancer_results$source[row] == "Mal/Epi") {
      result <- result %>%
        mutate(ligand_log2FC = log2(as.numeric(result$mean_ligand[result$ident == "Malignant cell"])/result$mean_ligand),
               receptor_log2FC = log2(as.numeric(result$mean_receptor[result$ident == "Immune cell"])/result$mean_receptor))
      #receptor_log2FC = log2(as.numeric(result$mean_receptor[c("B cell","CD4+ cell", "CD8+ cell", "Myeloid cell", "T cell", "NK cell") %in% result$ident])/result$mean_receptor))
      sc_cancer_results$ligand_expr[row] <- as.numeric(result$mean_ligand[result$ident == "Malignant cell"])
      sc_cancer_results$receptor_expr[row] <- as.numeric(result$mean_receptor[result$ident == "Immune cell"])
    } else {
      result <- result %>%
        mutate(ligand_log2FC = log2(as.numeric(result$mean_ligand[result$ident == "Immune cell"])/result$mean_ligand),
               #ligand_log2FC = log2(as.numeric(result$mean_ligand[c("B cell","CD4+ cell", "CD8+ cell", "Myeloid cell", "T cell", "NK cell") %in% result$ident])/result$mean_ligand),
               receptor_log2FC = log2(as.numeric(result$mean_receptor[result$ident == "Malignant cell"])/result$mean_receptor))
      sc_cancer_results$ligand_expr[row] <- as.numeric(result$mean_ligand[result$ident == "Immune cell"])
      sc_cancer_results$receptor_expr[row] <- as.numeric(result$mean_receptor[result$ident == "Malignant cell"])
    }
    sc_cancer_results$ligand_log2FC_min[row] <- min(result$ligand_log2FC)
    sc_cancer_results$ligand_log2FC_max[row] <- max(result$ligand_log2FC)
    sc_cancer_results$receptor_log2FC_min[row] <- min(result$receptor_log2FC)
    sc_cancer_results$receptor_log2FC_max[row] <- max(result$receptor_log2FC)
    cat("done with row", row, "\n")
  }
  
  # save analysis results
  write.csv(sc_cancer_results, paste0("./out/", tissue, "/sc_cancer_results_processed.csv"), row.names=F, quote=F)
  write.csv(sc_healthy_results, paste0("./out/", tissue, "/sc_healthy_results_processed.csv"), row.names=F, quote=F)
  
  # merge filtered cancer and healthy data, filter for most relevant interactions, calculate score ratios and log2FCs
  sc_common_lri <- merge(sc_cancer_results, sc_healthy_results,  by = c("interaction", "source", "target", "ligand", "receptor"), suffixes = c("_cancer", "_healthy"), all = T) %>% 
    select(-specificity_score_cancer, -specificity_score_healthy, -magnitude_score_cancer, -magnitude_score_healthy) %>% 
    mutate(score_log2FC = log2(aggregate_score_cancer/aggregate_score_healthy)) %>% 
    arrange(score_log2FC) %>% 
    na.omit()
  
  # categorize log2FCs
  sc_common_lri$color_category <- cut(sc_common_lri$score_log2FC,
                                      breaks = c(-Inf, -0.2, 0.2, Inf),
                                      labels = c("healthy-specific", "condition-unspecific", "cancer-specific"))
  # save standalone interactions
  sc_common_lri$score_log2FC[is.na(sc_common_lri$aggregate_score_healthy)] <- Inf
  sc_cancer_standalone <- sc_common_lri[which(sc_common_lri$score_log2FC == Inf),]
  sc_cancer_standalone <- sc_cancer_standalone %>% 
    arrange(desc(aggregate_score_cancer))
  sc_common_lri$score_log2FC[(sc_common_lri$score_log2FC == Inf)] <- NA
  
  sc_common_lri$score_log2FC[is.na(sc_common_lri$aggregate_score_cancer)] <- -Inf
  sc_healthy_standalone <- sc_common_lri[which(sc_common_lri$score_log2FC == -Inf),]
  sc_healthy_standalone <- sc_healthy_standalone %>% 
    arrange(desc(aggregate_score_healthy))
  sc_common_lri$score_log2FC[(sc_common_lri$score_log2FC == -Inf)] <- NA
  
  
  # filter for interactions relevant for cancer (log2FC > 0) 
  sc_cancer_all_lri <- sc_common_lri %>%
    filter(score_log2FC > 0) %>%
    filter(ligand_expr >0.5,
           receptor_expr>0.5) %>%           
    filter(ligand_log2FC_min == 0) %>%      
    filter(receptor_log2FC_min == 0) %>%    
    filter(ligand_log2FC_max >= 0.5) %>%    
    filter(receptor_log2FC_max >= 0.5) %>%  
    arrange(desc(score_log2FC)) %>% 
    mutate(rank = rank(-score_log2FC, ties.method = "min"),
           interaction = paste0(ligand, "^", receptor)) %>% 
    na.omit()

  # filter for interactions relevant for cancer (logFC > 0) with malignant cells as a source
  sc_cancer_mal_source <- sc_cancer_all_lri[sc_cancer_all_lri$source == "Mal/Epi",] %>% 
    mutate(interaction = paste0(ligand, "^", receptor)) %>% # shortens interaction name since directionality is clear
    select(source, target, interaction, aggregate_score_healthy, aggregate_score_cancer, score_log2FC) %>% 
    mutate(rank = rank(-score_log2FC, ties.method = "min")) %>% 
    arrange(desc(score_log2FC)) %>% 
    na.omit()
  # Convert interactions to factor levels for easier visualisation in bar plots
  sc_cancer_mal_source$interaction <- factor(sc_cancer_mal_source$interaction, levels = sc_cancer_mal_source$interaction)
  
  # Filter for interactions relevant for cancer (logFC > 0) with malignant cells as a source
  sc_cancer_mal_target <- sc_cancer_all_lri[sc_cancer_all_lri$target == "Mal/Epi",] %>% 
    mutate(interaction = paste0(ligand, "^", receptor)) %>% # shortens interaction name since directionality is clear
    select(source, target, interaction, aggregate_score_healthy, aggregate_score_cancer, score_log2FC) %>% 
    mutate(rank = rank(-score_log2FC, ties.method = "min")) %>% 
    arrange(desc(score_log2FC)) %>% 
    na.omit() #%>% 
  # Convert interactions to factor levels for easier visualisation in bar plots
  sc_cancer_mal_target$interaction <- factor(sc_cancer_mal_target$interaction, levels = sc_cancer_mal_target$interaction)
  
  # save single-cell ranking results
  write.csv(sc_cancer_all_lri, paste0("./out/", tissue, "/sc_cancer_full_ranking.csv"), row.names=F, quote=F)
  write.csv(sc_cancer_mal_source, paste0("./out/", tissue, "/sc_cancer_mal_source_ranking.csv"), row.names=F, quote=F)
  write.csv(sc_cancer_mal_target, paste0("./out/", tissue, "/sc_cancer_mal_target_ranking.csv"), row.names=F, quote=F)
  write.csv(sc_common_lri, paste0("./out/", tissue, "/sc_common_lri.csv"), row.names=F, quote=F)
  
  } else {
  sc_cancer_results <- read_csv(paste0("./out/", tissue, "/sc_cancer_results_processed.csv"))
  sc_healthy_results <- read_csv(paste0("./out/", tissue, "/sc_healthy_results_processed.csv"))
  sc_common_lri <- read_csv(paste0("./out/", tissue, "/sc_common_lri.csv"))
  sc_cancer_all_lri <- read_csv(paste0("./out/", tissue, "/sc_cancer_full_ranking.csv"))
  sc_cancer_mal_source <- read_csv(paste0("./out/", tissue, "/sc_cancer_mal_source_ranking.csv"))
  sc_cancer_mal_target <- read_csv(paste0("./out/", tissue, "/sc_cancer_mal_target_ranking.csv"))
  }


#### Plots ####
if (create_sc_plots == T) {
  # scatterplot of all interactions
  scatterplot_pseudobulk <- ggplot(sc_common_lri, aes(x = aggregate_score_healthy, y = aggregate_score_cancer, color = color_category, label = interaction)) +
    geom_point(size = 1) +
    scale_color_manual(values = c("healthy-specific" = "#0F69AF", "condition-unspecific" = "gray", "cancer-specific" = "#E61E50")) +  # Specify colors
    theme_classic()+
    geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed", color = 'gray', ) +
    geom_label_repel(data = sc_common_lri[which((sc_common_lri$score_log2FC)>(1) & sc_common_lri$aggregate_score_cancer>0.9),], 
                     size = 3,
                     min.segment.length = unit(0, 'lines'),
                     show.legend = FALSE,
                     xlim = c(-0.3,0.5),
                     ylim = c(0.5,1.5),
                     max.overlaps = 10,
                     fill = "white",
                     force = 2,
                     max.time = 20,
                     max.iter = 1000000)+
    labs(color = "Tissue Specificity") +
    coord_fixed(clip = "off") +
    theme(plot.margin = unit(c(0.2, 0, 0, 1), "cm"))  # Increases the margin around the plot

  
  # scatterplot of interactions after filtering
  sc_common_lri_categorized <- sc_common_lri
  sc_common_lri_categorized$color_category_filtered <- NA
  sc_common_lri_categorized$color_category_filtered[sc_common_lri_categorized$ligand_log2FC_min != 0 | sc_common_lri_categorized$receptor_log2FC_min != 0 | sc_common_lri_categorized$ligand_log2FC_max < 0.5 | sc_common_lri_categorized$receptor_log2FC_max < 0.5] <- "cell type-unspecific"
  sc_common_lri_categorized$color_category_filtered[sc_common_lri_categorized$ligand_expr<0.5 | sc_common_lri_categorized$receptor_expr<0.5] <- "low expression"
  sc_common_lri_categorized$color_category_filtered[is.na(sc_common_lri_categorized$color_category_filtered)] <- as.character(sc_common_lri_categorized$color_category[is.na(sc_common_lri_categorized$color_category_filtered)])
  sc_common_lri_categorized$color_category_filtered[sc_common_lri_categorized$color_category == "healthy-specific"] <-  as.character(sc_common_lri_categorized$color_category[sc_common_lri_categorized$color_category == "healthy-specific"])
  sc_common_lri_categorized$color_category_filtered[sc_common_lri_categorized$color_category == "condition-unspecific"] <-  as.character(sc_common_lri_categorized$color_category[sc_common_lri_categorized$color_category == "condition-unspecific"])
  
  scatterplot_pseudobulk_filtered <- ggplot(sc_common_lri_categorized, aes(x = aggregate_score_healthy, y = aggregate_score_cancer, color = color_category_filtered, label = interaction)) +
    geom_point(size = 0.7) +
    scale_color_manual(values = c("healthy-specific" = "#0F69AF", "condition-unspecific" = "gray", "cancer-specific" = "#E61E50", "low expression" = "#E1C3CD", "cell type-unspecific" = "#FFDCB9")) +  
    theme_classic()+
    geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed", color = 'gray') +
    geom_label_repel(data = sc_common_lri_categorized[sc_common_lri_categorized$color_category_filtered == "cancer-specific" & sc_common_lri_categorized$score_log2FC > 0.6,],  size = 3, min.segment.length = unit(0, 'lines'), show.legend = FALSE, ylim = c(0.75, 1.2))+
    #ggtitle('Ligand-Receptor analysis: Lung cancer cohort from the IMPACT atlas',
    #        subtitle = 'Malignant/Epithelial cells <-> Immune cells \nred dots are Ligand-Receptor pairs, very specific to tumors') +
    labs(color = "Tissue Specificity") +
    coord_fixed(clip = "off")
  
  # ranked interactions bar plot
  n_interactions <- length(sc_cancer_all_lri$interaction)
  if (n_interactions > 50) {
    barplot_pseudobulk_top_LRIs <- ggplot(sc_cancer_all_lri[50:1,], aes(y = factor(sc_cancer_all_lri$interaction[50:1], levels = sc_cancer_all_lri$interaction[50:1]), x = score_log2FC, fill = source)) +
      geom_bar(stat = "identity") +
      labs(y = "Interaction Pairs", 
           x = "score_log2FC") +
      theme_classic() +                                                                                        
      geom_text(aes(label = format_significant_digits(score_log2FC, 2)), size = 3, hjust = -0.5) +
      xlim(c(0, 2)) +
      scale_fill_manual(values = c("#0F69AF", "#E61E50"))
    print("Bar plot is created with top 50 interactions.")
  } else {
    barplot_pseudobulk_top_LRIs <- ggplot(sc_cancer_all_lri[n_interactions:1,], aes(y = factor(sc_cancer_all_lri$interaction[n_interactions:1], levels = sc_cancer_all_lri$interaction[n_interactions:1]), x = score_log2FC, fill = source)) +
      geom_bar(stat = "identity") +
      labs(y = "Interaction Pairs", 
           x = "score_log2FC") +
      theme_classic() +                                                                                        
      geom_text(aes(label = format_significant_digits(score_log2FC, 2)), size = 3, hjust = -0.5) +
      xlim(c(0, 2)) +
      scale_fill_manual(values = c("#0F69AF", "#E61E50"))
  }

  
  # bar plot with malignant cells as source
  n_interactions_mal_source <- length(sc_cancer_mal_source$interaction)
  barplot_pseudobulk_mal_source <- ggplot(sc_cancer_mal_source[n_interactions_mal_source:1,], aes(y = factor(sc_cancer_mal_source$interaction[n_interactions_mal_source:1], levels = sc_cancer_mal_source$interaction[n_interactions:1]), x = score_log2FC)) +
    geom_bar(stat = "identity", fill = "#E61E50") +
    labs(subtitle = "Malignant Cells as sources of ligands",
         y = "Interaction", 
         x = "log2FC") +
    #guides(y = guide_axis(angle = 45)) +
    xlim(c(0,2))+
    theme_classic() +                                                                                        
    geom_text(aes(label = format_significant_digits(score_log2FC, digits = 2)), size = 3, hjust = -0.5) +
    theme(plot.margin = unit(c(0,1,0,0), "cm"))
  
  # bar plot with malignant cells as target
  n_interactions_mal_target <- length(sc_cancer_mal_target$interaction)
  barplot_pseudobulk_mal_target <- ggplot(sc_cancer_mal_target[n_interactions_mal_target:1,], aes(y = factor(sc_cancer_mal_target$interaction[n_interactions_mal_target:1], levels = sc_cancer_mal_target$interaction[n_interactions_mal_target:1]), x = score_log2FC)) +
    geom_bar(stat = "identity", fill = "#0F69AF") +
    labs(subtitle = "Immune Cells as sources of ligands",
         y = "Interaction", 
         x = "log2FC") +
    xlim(c(0,2))+
    #guides(y = guide_axis(angle = 45)) +
    theme_classic()+                                                                                        
    geom_text(aes(label = format_significant_digits(score_log2FC, digits = 2)), size = 3, hjust = -0.5) +
    theme(plot.margin = unit(c(0,0,0,1), "cm"))
}

if (display_sc_plots == T) {
  print(scatterplot_pseudobulk)
  print(scatterplot_pseudobulk_filtered)
  print(barplot_pseudobulk_top_LRIs)
  print(barplot_pseudobulk_mal_source)
  print(barplot_pseudobulk_mal_target)
}

if (save_sc_plots == T) {
  ggsave(paste0("./out/", tissue, "/scatterplot_pseudobulk.png"), scatterplot_pseudobulk)
  ggsave(paste0("./out/", tissue, "/scatterplot_pseudobulk_filtered.png"), scatterplot_pseudobulk_filtered)
  ggsave(paste0("./out/", tissue, "/barplot_pseudobulk_top50.png"), barplot_pseudobulk_top50)
  
}

