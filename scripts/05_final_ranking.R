if (run_final_analysis == T) {
  #### all LRIs ####
  final_ranking_all_lri <- merge(sp_cancer_ranking, sc_cancer_all_lri, by = "interaction", suffixes = c("_spatial", "_singlecell"), all = F) %>% 
    select(source, target, interaction, score_log2FC, rank_singlecell, logFC, rank_spatial) %>% 
    mutate(#interaction_long = paste0(source, "^", target, "_", interaction),
      score_log2FC = as.numeric(format_significant_digits(score_log2FC, 2)),
      logFC = as.numeric(format_significant_digits(logFC, 2)),
      score_singlecell = rescale(score_log2FC, to = c(min(score_log2FC), 1)),
      score_singlecell = as.numeric(format_significant_digits(score_singlecell, 2)),
      score_spatial = rescale(logFC, to = c(min(logFC), 1)),
      score_spatial = as.numeric(format_significant_digits(score_spatial, 2)),
      combined_score = sqrt(score_singlecell*score_spatial),
      combined_score = as.numeric(format_significant_digits(combined_score, 2)),
      final_score = rescale(combined_score, to = c(min(combined_score), 1)),
      final_score = as.numeric(format_significant_digits(final_score, 2))
    ) %>% 
    arrange(desc(final_score)) %>% 
    na.omit()
  
  #### Malignant Cells as Sources ####
  final_ranking_mal_source <- merge(sp_cancer_ranking, sc_cancer_mal_source, by = "interaction", suffixes = c("_spatial", "_singlecell"), all = T) %>% 
    select(interaction, score_log2FC, rank_singlecell, logFC, rank_spatial) %>% 
    mutate(score_singlecell = rescale(score_log2FC),
           score_spatial = rescale(logFC),
           final_score = rescale(sqrt(score_singlecell*score_spatial))) %>% 
    arrange(desc(final_score)) %>% 
    arrange(desc(final_score)) %>% 
    na.omit()
  final_ranking_mal_source$interaction <- factor(final_ranking_mal_source$interaction, levels = final_ranking_mal_source$interaction)
  
  #### Malignant Cells as Targets ####
  final_ranking_mal_target <- merge(sp_cancer_ranking, sc_cancer_mal_target, by = "interaction", suffixes = c("_spatial", "_singlecell"), all = T) %>% 
    select(interaction, score_log2FC, rank_singlecell, logFC, rank_spatial) %>% 
    mutate(score_singlecell = rescale(score_log2FC),
           score_spatial = rescale(logFC),
           final_score = rescale(sqrt(score_singlecell*score_spatial))) %>% 
    arrange(desc(final_score)) %>% 
    na.omit()
  final_ranking_mal_target$interaction <- factor(final_ranking_mal_target$interaction, levels = final_ranking_mal_target$interaction)

  write.csv(final_ranking_all_lri, paste0("./out/", tissue, "/final_ranking_all_lri.csv"), row.names=F, quote=F)
  write.csv(final_ranking_mal_source, paste0("./out/", tissue, "/final_ranking_mal_source.csv"), row.names=F, quote=F)
  write.csv(final_ranking_mal_target, paste0("./out/", tissue, "/final_ranking_mal_target.csv"), row.names=F, quote=F)
} else {
  final_ranking_all_lri <- read_csv(paste0("./out/", tissue, "/final_ranking_all_lri.csv"))
  final_ranking_mal_source <- read_csv(paste0("./out/", tissue, "/final_ranking_mal_source.csv"))
  final_ranking_mal_target <- read_csv(paste0("./out/", tissue, "/final_ranking_mal_target.csv"))
}


#### Plots ####
# Barplot
n_top_LRIs <- length(final_ranking_all_lri$interaction)
top_LRIs <- ggplot(final_ranking_all_lri[n_top_LRIs:1,], aes(y = factor(final_ranking_all_lri$interaction[n_top_LRIs:1], levels= final_ranking_all_lri$interaction[n_top_LRIs:1]), x = final_score, fill = source)) +
  geom_bar(stat = "identity") +
  labs(y = "Interaction Pairs", 
       x = "log2FC") +
  theme(panel.background = element_blank(), legend.position=c(0.7, 0.3), panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", 
                                                                     linewidth = rel(1)), complete = TRUE)+                                                                                        
  geom_text(aes(label = format_significant_digits(final_score, 2)), size = 3, hjust = -0.5) +
  scale_fill_manual(values = c("Mal/Epi" = "#E61E34", "Immune" = "#0F69AF"), name = "sender cell") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(c(10,10,10,10)))

print(top_LRIs)

if (save_final_plots == T) {
  ggsave(paste0("./out/", tissue, "/barplot_top_LRIs.png"), top_LRIs)
}
