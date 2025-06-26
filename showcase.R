#### project startup and installation ####
# select 1 on prompt
renv::init()

source("scripts/01_project_setup_installation.R")

#### load libraries, modified functions, settings ####
source("scripts/02_libraries_functions_settings.R")

#### set working directory ####
setwd("~/spatially-informed-interactions")

#### check current directory structure ####
if (!dir.exists("seurat_files")) {dir.create("seurat_files")}
if (!dir.exists("out")) {dir.create("out")}
if (!dir.exists("liana_results")) {
  dir.create("liana_results")
  stop("Please upload LIANA+ results to /liana_results.")
}

###########################################
########### PSEUDOBULK ANALYSIS ########### 
###########################################

#### Which tissue to inspect? ####
tissue <- "lung"
if (!dir.exists(paste0("./out/", tissue))) {dir.create(paste0("./out/", tissue))}

#### download IMPACT atlas and run analysis? (T) Or load previous results? (F) ####
run_download_IMPACT <- T

#### Plots ####
# create plots?
create_sc_plots <- T

# display plots?
display_sc_plots <- T

# save_plots?
save_sc_plots <- F

source("scripts/03_single_cell_analysis.R")


###########################################
############# SPATIAL ANALYSIS ############
###########################################

# run spatial analysis (T) or load previous spatial results (F)?
run_spatial_analysis <- T

#### Plots ####
# create plots?
create_sp_plots <- T

# display plots?
display_sp_plots <- T

# save_plots?
save_sp_plots <- F

source("scripts/04_spatial_analysis.R")


###########################################
############## FINAL RANKING ##############
###########################################

# run final analysis (T) or load previous results (F)?
run_final_analysis <- T

# save_plots?
save_final_plots <- F

source("scripts/05_final_ranking.R")


###########################################
###### OVER-REPRESENTATION ANALYSIS #######
###########################################

# run ORA (T) or load ORA results (F)?
run_ORA = T

# display plots?
display_ORA_plots = T

# save_plots?
save_ORA_plots = F

source("scripts/06_ORA_analysis.R")


