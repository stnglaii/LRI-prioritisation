Index <- read_table("Index_LR.txt")
Index_lung <-  Index %>% 
  filter(Tissue == "Lung") %>%
  filter(Ligand == "TF") #%>% 
  #filter(Receptor == "ROBO4")
cat(unique(Index_lung$DiseaseStage),sep = "\n")
  
  
  
  
  
Index_lung <- Index[which(Index$Tissue == "Lung"), ]
test <- Index_lung[which(Index_lung$LRpair == "Lung"), ]