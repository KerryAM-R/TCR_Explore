## volcano plots
install.packages(c("tidyverse","ggplot2" ,"ggrepel","shiny","shinyBS","gridExtra","DT","plyr","dplyr","reshape2","treemapify","circlize","scales", "readxl","vegan","VLF","randomcoloR","colourpicker","devtools","muscle","markdown"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("motifStack")
BiocManager::install("flowCore")
BiocManager::install("muscle")

library("devtools")
install_github("jokergoo/ComplexHeatmap")
install_github("mgledi/DiffLogo")
