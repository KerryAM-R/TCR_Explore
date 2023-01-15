install.packages("immunarch")           # Install the package
library(immunarch); data(immdata)       # Load the package and the test dataset
repOverlap(immdata$data) %>% vis()      # Compute and visualise the most important statistics:
geneUsage(immdata$data[[1]]) %>% vis()  #     public clonotypes, gene usage, sample diversity
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta) 
repExplore(immdata$data, "lens") %>% vis()  # Visualise the length distribution of CDR3
repClonality(immdata$data, "homeo") %>% vis()  
repOverlap(immdata$data) %>% vis()  # Build the heatmap of public clonotypes shared between repertoires
geneUsage(immdata$data[[1]]) %>% vis()  # Visualise the V-gene distribution for the first repertoire
repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta)  # Visualise the Chao1 diversity of repertoires, grouped by the patient status

install.packages("shiny")
install.packages("DT")
install.packages("RColorBrewer")
install.packages("reticulate")
