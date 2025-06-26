pseudobulk = api$conn %>%
  tbl("literature_singlecelltypepseudobulk_rna_wts_gexp_20240828_prod_mura_0") %>%
  dplyr::filter(id_gene_symbol == "EGFR") %>%
  collect()
