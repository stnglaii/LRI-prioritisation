# # Install packages ####
install.packages("./packages/Matrix_1.6-4.tar.gz", repos = NULL, type = "source")
install.packages("./packages/SeuratObject_5.0.1.tar.gz", repos = NULL, type = "source")
install.packages("./packages/Seurat_5.0.3.tar.gz", repos = NULL, type = "source")
renv::snapshot(type = "all")# exclude = c("liana", "OmnipathR"))

