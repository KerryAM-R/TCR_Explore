<img src="IMAGES/Logo.png" width="400">

### TCR_Explore Shiny R application

TCR_explore was designed to aid in the processing and analysis of TCR repertiore for both alpha beta and gamma-delta chains

there are three sections to the application:
- Quality control steps (IMGT and IMGT+MIXCR)
 + this creates the file needed to graph the data
 + Please see the Workflow -> QC tab for details on the file name

Graphical sections: 
- The scTCR section 
    + The user needs to upload the unsummarised dataset 
    + This contains several graphs incluyding: Treemap, circular plot, length distributions, motif plots, pie graphs, heatmap and upset plot
    + for the inverse simpson cancluation a differnet file will need to be uploaded. Please see Workflow -> scTCR plots for details
- FACS index data
    + This section contains combining the paired TCR data from the QC section with the index sort data (FAC file)
    + Please see Workflow -> FACS index QC and plots for details 
