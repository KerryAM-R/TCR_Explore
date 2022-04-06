
<img src="www/Logo.png" width="200">

### TCR_Explore Shiny R application

TCR_Explore was designed as an open-access web server that analyses and visualises TCR repertoire data without the need for coding expertise. TCR_Explore introduces multiple pipelines using an automated process that includes pairing of αβ or γδ chains, as well as facilitating interrogation of linked flow cytometric index data for immunophenotyping analyses. Additionally, automated summarisation process from a single input file enables the creation of a variety of publication-ready analytical plots. 

There are three main sections:
- Quality control (QC) processes
 + Uses output files generated from IMGT^1
 + Workflow → QC tab
 + Creates a universal input file for TCR repertoire data analysis
 + Tutorial video available
- TCR analysis 
    + User uploads the paired file generated from TCR_Explore QC process
    + Alternatively the user can  upload a file from other outputs (e.g. iRepertiore), which need to include the following column names: cloneCount, Indiv, group, Indiv.group
    + Several analytical graph features available including Treemap, Chord diagram, Pie chart, Motif analysis, Diversity and chain usage, and Overlap for comparison of multiple datasets (Heatmap and Upset plots)
    + For more information on the functions, see the TCR analysis information tab
- Paired TCR with Index data 
    + User uploads the paired file generated from TCR_Explore QC process and a corresponding .fcs (FACS index data) file
    + The merged file undergoes further QC process in the 'data cleaning steps'
        1. Changes the flow cytometric values from negative to small positive 
        2. User can filter using the clone count for coloring purposes (0=all values included)
    + This clean file is then used to create the dotplot, which has over 20 cusomisable features
    + For more information on the functions, see Paired TCR with Index data information tab

Please contact: Kerry.Mullan@monash.edu or Nicole.Mifsud@monash.edu to report errors.

Biomedicine Discovery Institute and Department of Biochemistry and Molecular Biology, Monash University, Melbourne, VIC 3800, Australia

<img src="www/Monash-BDI-logo-2016-1.png" width="600">

##### References:
1.	Lefranc MP, Giudicelli V, Duroux P, Jabado-Michaloud J, Folch G, Aouinti S, et al. IMGT(R), the international ImMunoGeneTics information system(R) 25 years on. Nucleic Acids Res. 2015;43(Database issue):D413-22.
