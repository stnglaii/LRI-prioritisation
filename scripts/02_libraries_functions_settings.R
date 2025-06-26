# Load packages ####
library(readr)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(scales)
library(Seurat)
library(qs)
library(ProfilerAPI2)
library(limma)
library(clusterProfiler)
library(enrichplot)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
setwd("~/lri-prioritization")


# Add Merck Colors ####
merck_colors <- c(violet = "#503291", blue = "#0F69AF", green = "#149B5F", red = "#E61E50", 
                  pink = "#EB3C96", lightblue = "#2DBECD", lightgreen = "#A5CD50",
                  yellow = "#FFC832", palepink = "#E1C3CD", paleblue = "#96D7D2",
                  palegreen = "#B4DC96", paleyellow = "#FFDCB9", grey = "#999999")


# Connection to ProfilerAPI2 ####
api = ProfilerAPI2::profiler_api(profile = "default")

# Custom function to format numbers to a specified number of significant digits
format_significant_digits <- function(x, digits) {
  formatC(x, format = "fg", digits = digits, flag = "#")  # Use flag = "#" to keep trailing zeros
}
