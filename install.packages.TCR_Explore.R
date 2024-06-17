## volcano plots
install.packages(c("tidyverse","ggplot2" ,"ggrepel","shiny","shinyBS","gridExtra","DT","plyr","dplyr","reshape2","treemapify","circlize","scales", "readxl","vegan","VLF","randomcoloR","colourpicker","devtools","muscle","markdown","umap"))

install.packages("circlize")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("motifStack")
BiocManager::install("flowCore")
BiocManager::install("muscle")
BiocManager::install("sangerseqR")
BiocManager::install("muscle")
BiocManager::install("fossil")
BiocManager::install("DiffLogo")
BiocManager::install("ComplexHeatmap")

library("devtools")
# install_github("jokergoo/ComplexHeatmap")
# install_github("mgledi/DiffLogo")

