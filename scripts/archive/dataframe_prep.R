#install.packages("RobustRankAggreg")
library(readr)
library(tidyverse)
library(ggrepel)
library(RobustRankAggreg)
conflicts_prefer(dplyr::filter)

cancer <- read_csv("results/pseudobulk_cancer_results.csv")

# #### Problem: we can't just multiply them or easily take a mean because the distributions are really different
# ggplot(cancer, aes(x=specificity_rank)) +
#   geom_histogram(binwidth=.001) +
#   xlim(0, 0.05) +
#   labs(title = "Distribution of specificity rank", subtitle = "pseudobulk_cancer" )
# 
# ggplot(cancer, aes(x=magnitude_rank)) +
#   geom_histogram(binwidth=.001) +
#   xlim(0, 0.05) +
#   labs(title = "Distribution of magnitude rank", subtitle = "pseudobulk_cancer" )
  
cancer <- cancer %>%
  select(source, target, ligand_complex, receptor_complex, specificity_rank, magnitude_rank) %>%
  mutate(interaction = paste0(ligand_complex, "^", receptor_complex )) %>%
  #filter(specificity_rank<0.05) %>% 
  #filter(magnitude_rank<0.05) %>% 
  arrange((specificity_rank)) %>%
  mutate(rank_specificity = row_number()) %>%
  arrange((magnitude_rank)) %>%
  mutate(rank_magnitude = row_number())

ranked_lists <- list(
  specificity = cancer$interaction[order(cancer$rank_specificity)],
  magnitude = cancer$interaction[order(cancer$rank_magnitude)]
)

rra_results <- RobustRankAggreg::aggregateRanks(ranked_lists)
colnames(rra_results) <- c("interaction", "Score")
cancer <- left_join(cancer, rra_results, by = "interaction")

