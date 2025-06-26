# Define the specific ligand-receptor pairs you are interested in
ligand_receptor_pairs <- list(
  c("EGF", "EGFR"),
  c("VEGFA", "FLT1"),
  c("VEGFA", "KDR"),
  c("HGF", "MET"),
  c("CD274", "PDCD1"),
  c("TGFB1", "TGFBR1"),
  c("TGFB1", "TGFBR2"),
  c("IGF1", "IGF1R"),
  c("KITLG", "KIT"),
  c("CXCL12", "CXCR4"),
  c("FGF2", "FGFR1"),
  c("TGFB1", "TGFBR1"),
  c("BCL2", "BAX"),
  c("NTS", "NTSR1"),
  c("GRP", "CRPR"),
  c("DLL3", "NOTCH1"),
  c("DLL3", "NOTCH2"),
  c("DLL3", "NOTCH3"),
  c("DLL3", "NOTCH4"),
  c("NPY", "NPY1R")
)



# Initialize an empty dataframe to store the results
pos_controls <- data.frame()

# Loop through each pair and query the data
for (pair in ligand_receptor_pairs) {
  GOI1 <- pair[1]
  GOI2 <- pair[2]
  
  # Query for the specific ligand-receptor pair and push them into the dataframe
  LRI_df <- liana_data_aggregated %>%
    filter((str_detect(ligand.complex, GOI1) & str_detect(receptor.complex, GOI2)) |
             (str_detect(ligand.complex, GOI2) & str_detect(receptor.complex, GOI1))) %>%
    filter(source == "Malignant Cell" | target == "Malignant Cell") #%>%
    distinct_at(c("ligand.complex", "receptor.complex"), .keep_all = TRUE) #keeps other columns as well
  
  # Append the results to the pos_controls dataframe
  pos_controls <- bind_rows(pos_controls, LRI_df)
  pos_controls <- arrange(pos_controls, original_rank)
}

# View the resulting dataframe
print(pos_controls)

liana_dotplot(pos_controls, ntop = 20) +   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
                                                 strip.text = element_text(size = 12)) #+

liana_data_aggregated_filtered <- liana_data_aggregated %>%
  mutate(labels = paste0(ligand.complex, " | ", receptor.complex)) %>%
  mutate(known_in_literature = (ligand.complex %in% pos_controls$ligand.complex & receptor.complex %in% pos_controls$receptor.complex)) %>% 
  filter(source == "Malignant Cell" & target != "Malignant Cell" |
           source != "Malignant Cell" & target == "Malignant Cell") %>%
  arrange(aggregate_rank)


p5 <- ggplot(liana_data_aggregated_filtered[1:50,], aes(x = reorder(original_rank, aggregate_rank), y = -log10(aggregate_rank), fill = literature_LC)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("1" = "#B4DC96", "2" = "#FFC832", "3" = "#EB3C96")) +
  labs(title = "Top 50 LRI between malignant and other cells", 
       x = "Interaction Pairs", 
       y = "-log10(aggregate rank)") +
  scale_x_discrete(labels = liana_data_aggregated_filtered$labels[1:50]) +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() + 
  geom_text_repel(
    aes(label = format(aggregate_rank, scientific = TRUE, digits = 2)),
    force_pull   = 0, 
    nudge_y      = 0.05,
    angle        = 90,
    hjust        = 0.5,
    segment.size = 0.1,
    max.iter = 1e4, max.time = 1
  )   +
  ylim(0, 9.5) 
p5

ggsave("12_top50_pos_controls.png", plot = p5, create.dir = T, path = "./res/", height = 9, width = 18)



liana_data_aggregated_filtered$literature_LC <- NA
#1 = known in cancers in general
#2 = known for lung cancer according to database or literature
#3 = known for cancers, known in cancers acc. to DB and lit.
liana_data_aggregated_filtered$literature_LC[1:50] <- c(
  3,  # "COL1A2 | SDC1"
  2,  # "COL1A1 | DDR1"
  3,  # "COL1A1 | SDC1"
  3,  # "COL3A1 | DDR1"
  2,  # "COL6A3 | SDC1"
  1,  # "COL6A1 | SDC1"
  1,  # "COL6A2 | SDC1"
  3,  # "DCN | MET"
  3,  # "HSPG2 | SDC1"
  1,  # "DCN | EGFR"
  1,  # "SFTPA2 | CD93"
  2,  # "CALCA | CALCRL_RAMP3"
  2,  # "CALCA | CALCRL_RAMP2"
  3,  # "COL4A1 | SDC1"
  2,  # "COL4A2 | SDC1"
  1,  # "CXCL2 | ACKR1"
  1,  # "MDK | TSPAN1"
  2,  # "CXCL8 | ACKR1"
  2,  # "COL1A2 | ITGA3_ITGB1"
  1,  # "COL4A1 | SDC1"
  1,  # "AVP | RAMP2"
  2,  # "PTH2 | RAMP2"
  2,  # "TIMP3 | MET"
  2,  # "CALCB | RAMP2"
  2,  # "TIMP3 | MET"
  1,  # "GHRH | RAMP2"
  1,  # "ADM2 | RAMP2"
  2,  # "GCG | RAMP2"
  2,  # "CXCL1 | ACKR1"
  1,  # "TIMP3 | DDR1"
  1,  # "COL1A1 | ITGA3_ITGB1"
  2,  # "COL4A2 | SDC1"
  3,  # "TIMP3 | DDR1"
  3,  # "THBS2 | SDC1"
  2,  # "THBS2 | ITGA3_ITGB1"
  1,  # "COL1A2 | ITGAV_ITGB8"
  1,  # "COL5A2 | DDR1"
  3,  # "GRP | FAP"
  2,  # "THBS2 | CD47"
  1,  # "NPS | RAMP2"
  1,  # "LHB | RAMP2"
  2,  # "CRH | RAMP2"
  1,  # "GIP | RAMP2"
  2,  # "VEGFA | FLT1_KDR"
  1,  # "CD24 | SELP"
  2,  # "SFRP2 | FZD5"
  1,  # "PDAP1 | PDGFRB"
  2,  # "VWF | ITGB1"
  1,  # "ADM | RAMP2"
  2   # "VCAN | EGFR"
)
liana_data_aggregated_filtered$literature_LC <- as.factor(liana_data_aggregated_filtered$literature_LC)
