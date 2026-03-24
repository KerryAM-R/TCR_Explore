# CRAN packages
cran_packages <- c(
  "markdown","rmarkdown","tidyverse","ggplot2","ggrepel",
  "shiny","shinyBS","gridExtra","DT","plyr","dplyr","reshape2",
  "treemapify","circlize","scales","readxl","RColorBrewer",
  "randomcoloR","colourpicker","ComplexHeatmap","muscle",
  "vegan","VLF","shinyWidgets","showtext","ggseqlogo",
  "umap","fpc","fossil","shinybusy","ggridges"
)

# Bioconductor packages
bioc_packages <- c(
  "motifStack","DiffLogo","sangerseqR"
)


# Bioconductor packages
bioc_packages <- c("motifStack","DiffLogo","sangerseqR")

# Install BiocManager if not already installed
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org/")

# Install missing Bioconductor packages
to_install_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(to_install_bioc) > 0){
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
} else {
  message("All Bioconductor packages already installed")
}


# Install missing CRAN packages
to_install_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]

if(length(to_install_cran) > 0){
  install.packages(to_install_cran, dependencies = TRUE, repos="https://cloud.r-project.org/")
} else {
  message("All CRAN packages already installed")
}

# Install Bioconductor packages if missing
if (!requireNamespace(
  
  
  
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

