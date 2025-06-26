combination <- cancer_2 %>% 
  mutate(interaction = paste0(ligand_complex, "^", receptor_complex )) %>%
  select(interaction, Score)

combi <- healthy %>% 
  mutate(interaction = paste0(ligand_complex, "^", receptor_complex )) %>%
  select(interaction, Score)

test <- left_join(combination, combi, by = 'interaction')

test <- test %>% c
arrange((specificity_rank)) %>%
  mutate(rank_specificity = row_number()) %>%
  arrange((magnitude_rank)) %>%
  mutate(rank_magnitude = row_number())
