if (run_ORA == T) {
  genes_healthy <- rownames(healthy)
  
  genes_cancer <- final_ranking_all_lri %>% 
    arrange(desc(final_score))
  genes_cancer$interaction <- sub("^[^_]*_", "", genes_cancer$interaction)
  genes_cancer <- unlist(strsplit(as.character(genes_cancer$interaction), "\\^"))
  
  ego_cancer <- enrichGO(gene = genes_cancer, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", universe = genes_healthy)
  ego_cancer <- pairwise_termsim(ego_cancer)
  qsave(ego_cancer, paste0("./out/", tissue, "/ORA_results.qs"))
} else {
  ego_cancer <- qread(paste0("./out/", tissue, "/ORA_results.qs"))
}

if (display_ORA_plots == T) {
  ORA_barplot_FoldEnrichment <- ego_cancer %>% 
    arrange(desc(FoldEnrichment)) %>% 
    # mutate(Description = case_when(
    #   Description == "T cell activation via T cell receptor contact with antigen bound to MHC molecule on antigen presenting cell" ~ "New Label 1",
    #   #Description == "Old Label 2" ~ "New Label 2",
    #   TRUE ~ Description)) %>%  # Keep other labels unchanged
    barplot(x = "FoldEnrichment", 
            showCategory = 30, 
            #title = "Top GO Biological Processes in Cancer", 
            font.size = 10,
            label_format = 70) +
    scale_fill_continuous(low="#E61E50", high="#0F69AF", name = 'p.adjust',
                          guide=guide_colorbar(reverse=TRUE)) +
    theme(legend.position=c(0.8, 0.3), panel.border = element_blank(),  panel.background = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", 
                                                                       linewidth = rel(1)), strip.background = element_rect(fill = "white", colour = "black", linewidth = rel(2)), complete = TRUE)
  
  ORA_barplot_GeneRatio <- ego_cancer %>% 
    arrange(desc(GeneRatio)) %>% 
    barplot(x = "GeneRatio", showCategory=30, title="Top GO Biological Processes in Cancer", font.size = 10, label_format = 70)
  
  ORA_dotplot_FoldEnrichment <- dotplot(ego_cancer, x = "FoldEnrichment", showCategory=30, title="GO Enrichment in Cancer")
  
  ORA_emapplot <- emapplot(ego_cancer, showCategory=30, title="Enrichment Map: Cancer")
  
  ORA_cnetplot <- cnetplot(ego_cancer, showCategory=10, title="Gene-Concept Network: Cancer")
  
  print(ORA_barplot_FoldEnrichment)
  print(ORA_barplot_GeneRatio)
  print(ORA_dotplot_FoldEnrichment)
  print(ORA_emapplot)
  print(ORA_cnetplot)
}