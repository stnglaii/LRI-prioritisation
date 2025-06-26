# Interactions in each dataset
interactions_cancer <- unique(cancer_2$interaction)
interactions_healthy <- unique(healthy$interaction)

# Common interactions
common_interactions <- intersect(interactions_cancer, interactions_healthy)

# Unique interactions
unique_cancer <- setdiff(interactions_cancer, interactions_healthy)
unique_healthy <- setdiff(interactions_healthy, interactions_cancer)

# For common interactions
data_common <- merge(
  cancer_2[cancer_2$interaction %in% common_interactions, ],
  healthy[healthy$interaction %in% common_interactions, ],
  by = "interaction",
  suffixes = c("_cancer", "_healthy")
)

# For unique interactions
data_unique_cancer <- cancer_2[cancer_2$interaction %in% unique_cancer, ]
data_unique_healthy <- healthy[healthy$interaction %in% unique_healthy, ]


data_common <- data_common %>%
  mutate(rank_difference = Score_cancer - Score_healthy,
         rank_ratio = Score_cancer / Score_healthy) %>% 
  filter(rank_ratio != 1) %>% 
  select(interaction, Score_cancer, Score_healthy)

combined_data <- merge(cancer_2, healthy, by = "interaction", suffixes = c("_cancer", "_healthy"), all = TRUE)

combined_data <- combined_data %>%
  mutate(rank_ratio = Score_cancer / Score_healthy
  )












# Assuming 'cancer' and 'healthy' are data frames with interactions and ranks
# Merge datasets on interaction identifier
# Perform Wilcoxon signed-rank test
test_result <- wilcox.test(data_common$Score_cancer, data_common$Score_healthy, paired = TRUE)

library(pheatmap)

# Prepare data matrix
rank_matrix <- data_common %>%
  select(Score_cancer, Score_healthy)

pheatmap(rank_matrix, scale = "row", clustering_distance_rows = "euclidean",
         clustering_method = "complete", main = "Heatmap of Interaction Ranks")
