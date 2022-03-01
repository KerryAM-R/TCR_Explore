
## volcano plots
require("tidyverse")
require("ggplot2") #Best plots
require("ggrepel") #Avoid overlapping labels
require("shiny")
require("shinyBS")
require("gridExtra")
require("DT")
require("plyr")
require("dplyr")
require("reshape2")
require("treemapify") # tremap plot
require("circlize")
require("motifStack")
require("scales") # to access break formatting functions
require("flowCore")
require("readxl")
require("RColorBrewer")
require("randomcoloR") 
require("colourpicker") # selectively colour
require("ComplexHeatmap")
require("muscle") # aligning sequences
require("DiffLogo") # comparing motif plots
require("vegan") # diversity statistic
require("VLF") ## aa.count.function

test_fun <- function()
{
  for (i in 1:15) {
    incProgress(1/15)
    sum(runif(1000000,0,1))
  }
}

index <- c("shannon", "simpson","invsimpson")
ASN$cols <- colorset(alphabet="AA",
                     colorScheme="chemistry")
draw_colnames_rotate <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, hjust = .5, rot = 0, gp = gpar(...)) # rot = rotation for # degrees
  return(res)}
gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Nucleotide <- function (Nucleotide, seqlength) {
  nt <- c("A", "C", "G", "T")
  spec.no <- nrow(Nucleotide)
  count <- mat.or.vec(nr = 4, nc = seqlength)
  for (i in 1:seqlength) {
    count[1, i] <- length(which(Nucleotide[, i + 2] == nt[1]))
    count[2, i] <- length(which(Nucleotide[, i + 2] == nt[2]))
    count[3, i] <- length(which(Nucleotide[, i + 2] == nt[3]))
    count[4, i] <- length(which(Nucleotide[, i + 2] == nt[4]))
    
  }
  rownames(count) <- nt
  return(count)
}

options(shiny.maxRequestSize=10*1024^2)
credentials <- data.frame(
  user = c("shiny", "shinymanager"),
  password = c("azerty", "12345"),
  stringsAsFactors = FALSE
)

graph_type <- c("histogram","density")
axis_density_group <- c("x-axis","y-axis")

angle <- c(0,90,180,270)

error_message_val1 <- "No data found"
error_message_val2 <- "Uploading file"
error_message_val3 <- "Upload clone file"
error_message_val4 <- "no own list found\n \nSuggest uploading file\nheaders=ID"

simp.index.names <- c("inv.simpson.index","total # clones","unique # clones","V1","V2","Indiv_group")

# user interface  ----
ui <- navbarPage(title = tags$img(src = "Logo.png", height = 70, width = 120,style = "margin:-25px 10px"), position = "fixed-top",collapsible = TRUE,
                 tags$head(
                   tags$style(HTML(' .navbar {
                          height: 80px;
                          min-height:80px !important;
                        }
                      .navbar-nav > li > a, .navbar-brand {
                            padding-top:30px !important; 
                            padding-bottom:30px !important;
                            height: 20px;
                            }'))),
                 
                 navbarMenu("TCR_Explore workflow",
                            tabPanel("Overview",
                                     fluidRow(includeMarkdown("README.md")),
                                     # tags$video(id="video2", type = "video/mp4",src = "test.mp4", controls = "controls", height="720px")
                            ),     
                            tabPanel("QC",
                                     h3("Tutorial video of Quality control processes"),
                                     uiOutput("video"),
                                     fluidRow(includeMarkdown("READMEQC.md")),
                                     
                                     # tags$video(id="video2", type = "video/mp4",src = "test.mp4", controls = "controls", height="720px")
                            ),     
                            tabPanel("scTCR plots",
                                     fluidRow(includeMarkdown("README.scTCR.md"))),
                            tabPanel("FACS index QC and plots",
                                     fluidRow(includeMarkdown("README.FACS.md"))),
                            
                            tabPanel("Session info", 
                                     tabPanel("Session info", verbatimTextOutput("sessionInfo"))
                            )
                            
                            #          fluidRow(includeMarkdown("README.FACS.md"))
                            # )
                            
                 ),
                 # QC ----
                 navbarMenu("QC",
                            tags$head(
                              tags$style(type = 'text/css', 
                                         HTML('.navbar { background-color: white;}
                          .navbar-default .navbar-brand{color: white;}
                          .tab-panel{ background-color: red; color: white}
                          .navbar-default .navbar-nav > .active > a, 
                           .navbar-default .navbar-nav > .active > a:focus, 
                           .navbar-default .navbar-nav > .active > a:hover {
                                color: #555;
                                background-color: darkblue;
                                color:white
                                
                            }')
                              )
                            ),
                            
                            # tabPanel("Making fasta files",
                            #          directoryInput('directory', label = 'select a directory'),
                            #          verbatimTextOutput("dir", placeholder = TRUE),  
                            #          actionButton("do", "Click Me to make fasta file (50 per file)"),
                            #          h4('Merging statistics'), 
                            #          uiOutput('textWithHTML') # ui output as a list of HTML p() tags
                            #          ),
                            # 
                            # UI IMGT only ----
                            tabPanel("IMGT",
                                     sidebarLayout(
                                       sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                                    
                                                    selectInput("dataset_IMGT3", "Choose a dataset:", choices = c("ab-test-data1", "own_data")),
                                                    fileInput('file_IMGT3', 'Select file for IMGT datafile',
                                                              accept=c('xls/xlsx', '.xls')),
                                                    
                                                    selectInput("dataset_IMGT_afterQC", "Choose a dataset:", choices = c("ab-test-data1", "own1")),
                                                    
                                                    fileInput('file_IMGT_afterQC', 'Chromatogram checked file (.csv)',
                                                              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                                    h5("option for paired and TCRdist outputs"),
                                                    selectInput("IMGT_chain2","Alpha-beta or gamma-delta",choices = c("ab","gd")),
                                                    selectInput("sheet2","Sheets included", choices = c("Summary+JUNCTION","Summary"))
                                                    
                                       ),
                                       mainPanel(
                                         tabsetPanel(
                                           tabPanel("IMGT create QC file",
                                                    h4("Fill in the 'clone_quality' column: pass or fail"), 
                                                    h4("Add comments if desired"),
                                                    fluidRow(column(4, selectInput("sheet","Sheets included", choices = c("Summary+JUNCTION","Summary"))),
                                                             column(8, selectInput("include.origin","Include VDJ (n/p) origins (Summary+JUNCTION only)",choices = c("no",'yes'), width = "800px")),
                                                    ),
                                                    tags$head(tags$style("#IMGT2_out  {white-space: nowrap;  }")),
                                                    div(DT::dataTableOutput("IMGT2_out")),
                                                    downloadButton('downloadTABLE_IMGTonly','Download table')
                                           ),
                                           tabPanel("Paired chain file",
                                                    
                                                    
                                                    
                                                    tags$head(tags$style("#chain_table_IMGT.QC1  {white-space: nowrap;  }")),
                                                    div(DT::dataTableOutput("chain_table_IMGT.QC1")),
                                                    p(" "),
                                                    downloadButton('downloadTABLE.QC1','Download paired chain file'),
                                                    
                                                    
                                           ),
                                           tabPanel("TCRdist output file",
                                                    textInput("tcr_lab","ID for TCRdist","human_tcr"),
                                                    
                                                    tags$head(tags$style("#chain_table_IMGT.QC1  {white-space: nowrap;  }")),
                                                    div(DT::dataTableOutput("chain_table_IMGT.tcrdist")),
                                                    
                                                    
                                                    h5(" "),
                                                    downloadButton('downloadTABLE.TSV','Download tsv file for TCRdist')
                                           )
                                         )
                                       )
                                       
                                     )
                            ),
                            
                            # do MIXCR+IMGT ------           
                            # tabPanel("IMGT+MiXCR",
                            #          sidebarLayout(
                            #            sidebarPanel(
                            # # IMGT+MIXCR file -----
                            #              id = "tPanel2",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                            #              selectInput("chain","alpha-beta or gamma-delta",choices = c("ab","gd")),
                            #              
                            #              selectInput("dataset_IMGT", "Choose a dataset:", choices = c("ab-test-data","gd-test-data", "own")),
                            #              fileInput('file_IMGT1', 'Select file for IMGT datafile (sheet 1)',
                            #                        accept=c('xls/xlsx', '.xls')),
                            #              fileInput('file_IMGT2', 'Select MiXCR txt file',
                            #                        accept=c('text/csv', 'text/comma-separated-values,text/plain', '.txt')),
                            #              fluidRow(
                            #                column(4,radioButtons('sep_IMGT2', 'Separator', c( Tab='\t', Comma=','), '\t')),
                            #                column(4,radioButtons('quote_IMGT2', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), '"')),
                            #                column(4, radioButtons("header_IMGT2",'header',c(yes=TRUE,no=FALSE)))
                            #              ),
                            #              
                            #              # merged MIXCR/IMGT or IMGT file with QC 
                            #              fileInput('file_IMGT.MiXCR', 'Select chromatogram checked file (.csv)',
                            #                        accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                            #            ),
                            #            mainPanel(
                            #              tabsetPanel(
                            #                tabPanel("IMGT+MiXCR",
                            #                         h5("Pre-processed file"),
                            #                         tags$head(tags$style("#unmerged_chain1  {white-space: nowrap;  }")),
                            #                         div(DT::dataTableOutput("unmerged_chain1")),
                            #                         tags$head(tags$style("#unmerged_chain2  {white-space: nowrap;  }")),
                            #                         div(DT::dataTableOutput("unmerged_chain2")),
                            #                         tags$head(tags$style("#merged_IMGT  {white-space: nowrap;  }")),
                            #                         div(DT::dataTableOutput("merged_IMGT")),
                            #                         downloadButton('downloadTABLE2','Download merged table')
                            #                         
                            #                ),
                            #                tabPanel("Unsummarised merged A/B or G/D",
                            #                         tags$head(tags$style("#chain_table  {white-space: nowrap;  }")),
                            #                         div(DT::dataTableOutput("chain_table")),
                            #                         downloadButton('downloadTABLE3','Download table'),
                            #                         
                            #                ),
                            #                tabPanel("summarised file for records",
                            #                         verbatimTextOutput("names.in.file2"),
                            #                         fluidRow(column(12, textInput("string.data2","column names for summary","Indiv.group, V.GENE.and.allele_A, J.GENE.and.allele_A, JUNCTION..AA._A, V.GENE.and.allele_B, J.GENE.and.allele_B, D.GENE.and.allele_B, JUNCTION..AA._B", width = "1200px") )),
                            #                         tags$head(tags$style("#chain_table2  {white-space: nowrap;  }")),
                            #                         div(DT::dataTableOutput("chain_table2")),
                            #                         downloadButton('downloadTABLE4','Download table'),
                            #                )
                            #                
                            #              )
                            #              
                            #            )
                            #          )
                            #          
                            # )     
                            #            
                 ),
                 
                 
                 
                 
                 # UI TCR plots ----
                 
                 tabPanel("TCR analysis",
                          
                          tags$style(HTML("
    .tabbable > .nav > li > a                  {background-color: white;  color:black}
    .tabbable > .nav > li[class=active]    > a {background-color: darkred; color:white}
  ")),
                          
                          sidebarLayout(
                            sidebarPanel(id = "tPanel",style = "overflow-y:scroll; max-height: 700px; position:relative;", width=3,
                                         tags$style(type="text/css", "body {padding-top: 80px; padding-left: 10px;}"),
                                         #textInput(inputId = "lab1", label = "Group label of file 1",value = "Ex.vivo"),
                                         tags$head(tags$style(HTML(".shiny-notification {position:fixed;top: 50%;left: 30%;right: 30%;}"))),
                                         tags$head(tags$style(HTML('.progress-bar {background-color: purple;}'))),
                                         selectInput("dataset", "Choose a dataset:", choices = c("test-data", "own")),
                                         fileInput('file2', 'Select file for single samples',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                         
                                         fluidRow(
                                           column(6,radioButtons('sep', 'Separator', c( Tab='\t', Comma=','), ',')),
                                           column(6,radioButtons('quote', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), '"'))
                                         ),
                                         
                                         selectInput("group_column",label = h5("Column of group"), ""),
                                         selectInput("type.tree",label = h5("Type of input"), choices =  c("scTCR","bulk")),
                                         tags$hr(),
                                         selectInput("shannon.index", "Choose a dataset:", choices = c("ab.test.index", "own.index")),
                                         fileInput('file_diversity.index', 'Diversity statistics file',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                            ),
                            
                            mainPanel(tabsetPanel(
                              tabPanel("Overview of TCR pairing",tabsetPanel(
                                # UI Summary table -----
                                tabPanel("Summary table",
                                         verbatimTextOutput("names.in.file3"),
                                         selectInput("type.chain","Alpha-beta or gamma-delta",choices = c("ab","gd")),
                                         fluidRow(column(12, selectInput("string.data3","column names for summary","",multiple = T, width = "1200px") )),
                                         tags$head(tags$style("#chain_table_IMGT.QC3  {white-space: nowrap;  }")),
                                         div(DT::dataTableOutput("chain_table_IMGT.QC3")),
                                         downloadButton('downloadTABLE.QC3','Download table')
                                         
                                ),
                                # UI Treemap -----
                                tabPanel("Treemap",
                                         fluidRow(
                                           
                                           column(3,  numericInput("nrow.tree",label = h5("Rows"), value = 1))
                                         ),
                                         
                                         
                                         fluidRow( 
                                           column(3, selectInput("tree_colour.choise",label = h5("Colour"), choices =  c("default","random","grey"))),
                                           column(3, selectInput("fill2",label = h5("Colour treemap by"),"" )),
                                           column(3, selectInput("sub_group2",label = h5("Separate panels by"),"" )),
                                           column(3,selectInput( "wrap",label = h5("Group"),"sub" )),
                                           column(3,selectInput( "count2",label = h5("Count column"),"")),
                                           column(3, selectInput("tree.lab",label = h5 ("Add label"),choices = c("yes","no")))
                                           
                                         ),
                                         fluidRow(column(3,
                                                         wellPanel(id = "tPanel21",style = "overflow-y:scroll; max-height: 600px",
                                                                   uiOutput('myPanel'))),
                                                  column(9,plotOutput("Treemap2", height="600px"))),
                                         fluidRow(
                                           column(3,numericInput("width_tree", "Width of PDF", value=10)),
                                           column(3,numericInput("height_tree", "Height of PDF", value=8)),
                                           column(3),
                                           column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_scTREE','Download PDF'))
                                         ),
                                         
                                         fluidRow(
                                           column(3,numericInput("width_png_tree","Width of PNG", value = 1600)),
                                           column(3,numericInput("height_png_tree","Height of PNG", value = 1200)),
                                           column(3,numericInput("resolution_PNG_tree","Resolution of PNG", value = 144)),
                                           column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_scTREE','Download PNG'))
                                         ),
                                ),
                                # UI circular plot -----
                                tabPanel("Chord diagram", 
                                         fluidRow(
                                           
                                           column(2,selectInput( "group_selected2",label = h5("Group"),"" )),
                                           column(2,selectInput( "chain1",label = h5("Chain one"),"" )),
                                           column(2,selectInput( "chain2",label = h5("Chain two"),"" )),
                                           column(2,style = "margin-top: 15px;", numericInput("chord.transparancy","Transparancy",value = 0.5, min=0,max=0.9)),
                                           column(2,tableOutput("table_display")),
                                           
                                         ),
                                         fluidRow(
                                           column(3, selectInput("circ.lab",label = h5("Add label"),choices = c("yes","no"))),
                                           column(3,selectInput( "colour_cir",label = h5("Colour"),choices = c("rainbow","random","grey"))),  
                                           column(3,style = "margin-top: 15px;", numericInput("seed.numb.chord","Random colour generator",value = 123)),
                                         ),
                                         fluidRow(column(3,
                                                         wellPanel(id = "tPanel22",style = "overflow-y:scroll; max-height: 600px",
                                                                   uiOutput('myPanel_circ'))),
                                                  column(9,plotOutput("Circular",height="600px"))),
                                         h4("Exporting the Circular plot"),
                                         fluidRow(
                                           
                                           column(3,numericInput("width_circ", "Width of PDF", value=10)),
                                           column(3,numericInput("height_circ", "Height of PDF", value=8)),
                                           column(3),
                                           column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_circ','Download PDF'))
                                         ),
                                         
                                         fluidRow(
                                           column(3,numericInput("width_png_circ","Width of PNG", value = 1600)),
                                           column(3,numericInput("height_png_circ","Height of PNG", value = 1200)),
                                           column(3,numericInput("resolution_PNG_circ","Resolution of PNG", value = 144)),
                                           column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_circ','Download PNG'))
                                         ),
                                         # tableOutput("out.col.table1")
                                         
                                ),
                                
                                
                                # UI Pie ----
                                tabPanel("Pie chart",
                                         fluidRow(column(2,selectInput("pie_chain",label = h5("Colour by this group"),"")),
                                                  column(2,selectInput("pie_colour.choise",label = h5("Colour"), choices =  c("default","random","grey"))),
                                                  column(2, selectInput("cir.legend",label=h5("Legend location"),choices = c("top","bottom","left","right","none"),selected = "bottom")),
                                                  column(2,  numericInput("nrow.pie",label = h5("Rows"), value = 3)),
                                                  column(2,  numericInput("size.circ",label = h5("Size of legend text"), value = 6))
                                                  
                                         ),
                                         
                                         
                                         
                                         fluidRow(column(3,
                                                         wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                   uiOutput('myPanel_pie'))),
                                                  column(9, plotOutput("pie_out",height="800px"))),
                                         fluidRow(
                                           column(3,numericInput("width_pie", "Width of PDF", value=10)),
                                           column(3,numericInput("height_pie", "Height of PDF", value=8)),
                                           column(3),
                                           column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_pie','Download PDF'))
                                         ),
                                         fluidRow(
                                           column(3,numericInput("width_png_pie","Width of PNG", value = 1600)),
                                           column(3,numericInput("height_png_pie","Height of PNG", value = 1200)),
                                           column(3,numericInput("resolution_PNG_pie","Resolution of PNG", value = 144)),
                                           column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_pie','Download PNG'))
                                         ),
                                ),
                                
                                
                              )),
                              
                              
                              
                              
                              
                              
                              tabPanel("Motif analysis",
                                       p("This section contains 4 tabs for motif analysis"),
                                       tabsetPanel(
                                         # UI CDR3 length distribution graphs ----- 
                                         tabPanel("CDR3 length distribution",
                                                  h5("length distribution plot requires an unsummarised dataset"),
                                                  fluidRow(
                                                    column(3,selectInput( "aa.or.nt",label = h5("Amino acid or nucleotide sequence column"),"" )),
                                                    column(3,selectInput( "selected_group_len",label = h5("Group"),"" ))
                                                  ),
                                                  fluidRow(
                                                    column(3, numericInput("bin","Bins of histogram",value=30)),
                                                    column(3, colourInput("hist_col","Colour of histogram","red")),
                                                    column(3, selectInput('graph_type', 'Type of graph', graph_type))
                                                  ),
                                                  plotOutput("Chain1_length"),
                                                  downloadButton("table_length","Download length table"),
                                                  fluidRow(
                                                    column(3,numericInput("width_length", "Width of PDF", value=6)),
                                                    column(3,numericInput("height_length", "Height of PDF", value=4)),
                                                    column(3),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_length','Download PDF'))
                                                  ),
                                                  fluidRow(
                                                    column(3,numericInput("width_png_length","Width of PNG", value = 960)),
                                                    column(3,numericInput("height_png_length","Height of PNG", value = 600)),
                                                    column(3,numericInput("resolution_PNG_length","Resolution of PNG", value = 144)),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_length','Download PNG'))
                                                  ),
                                         ),
                                         # UI motif -----
                                         tabPanel("Motif (amino acid)",
                                                  h5("Select amino acid column and CDR3 length"),
                                                  verbatimTextOutput("length"),
                                                  fluidRow(
                                                    column(3,selectInput( "aa.or.nt2",label = h5("Amino acid CDR3 column"),"" )),
                                                    column(3, numericInput("len","CDR3 amino acid length", value = 15)),                               
                                                    column(3,selectInput( "group_selected_motif",label = h5("Group"),"" ))
                                                  ),
                                                  fluidRow(
                                                    column(6,div(DT::dataTableOutput("length.table"))),
                                                    column(6,div(DT::dataTableOutput("Motif"))),
                                                  ),
                                                  plotOutput("Motif_plot"),
                                                  h4("Exporting amino acid plot"),
                                                  fluidRow(
                                                    column(3,numericInput("width_motif", "Width of PDF", value=10)),
                                                    column(3,numericInput("height_motif", "Height of PDF", value=3.5)),
                                                    column(3),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_motif','Download PDF'))),
                                                  
                                                  fluidRow(
                                                    column(3,numericInput("width_png_motif","Width of PNG", value = 1600)),
                                                    column(3,numericInput("height_png_motif","Height of PNG", value = 600)),
                                                    column(3,numericInput("resolution_PNG_motif","Resolution of PNG", value = 144)),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_motif','Download PNG'))
                                                  )),
                                         # UI motif NT -----
                                         tabPanel("Motif (nucleotide sequence)",
                                                  h5("Select nucleotide column and CDR3 length"),
                                                  verbatimTextOutput("length_nt"),
                                                  fluidRow(
                                                    column(3,selectInput( "aa.or.nt3",label = h5("Nucleotide CDR3 column"),"")),
                                                    column(3,selectInput( "group_selected",label = h5("Group"),"" )),
                                                    column(3, numericInput("len_nt","CDR3 nucleotide length", value = 30))
                                                  ),
                                                  fluidRow(
                                                    column(6,div(DT::dataTableOutput("length.table_nt"))),
                                                    column(6,div(DT::dataTableOutput("Motif_nt"))),
                                                  ),
                                                  plotOutput("Motif_plot_nt"),
                                                  h4("Exporting plot"),
                                                  fluidRow(
                                                    column(3,numericInput("width_motif_nt", "Width of PDF", value=10)),
                                                    column(3,numericInput("height_motif_nt", "Height of PDF", value=3.5)),
                                                    column(3),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_motif_nt','Download PDF'))),
                                                  
                                                  fluidRow(
                                                    column(3,numericInput("width_png_motif_nt","Width of PNG", value = 1600)),
                                                    column(3,numericInput("height_png_motif_nt","Height of PNG", value = 600)),
                                                    column(3,numericInput("resolution_PNG_motif_nt","Resolution of PNG", value = 144)),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_motif_nt','Download PNG'))
                                                  )
                                         ),
                                         # motif align with muscle -----
                                         tabPanel("Motif (AA or NT alignment)",
                                                  fluidRow(
                                                    column(3,selectInput("aa.or.nt4",label = h5("CDR3 column"),"")),
                                                    column(3,selectInput("group_selected_one",label = h5("First group (bottom of plot)"),"" )),
                                                    column(3,selectInput("group_selected_two",label = h5("Second group (top of plot)"),"" )),
                                                    
                                                  ),
                                                  p("ASN = amino acid data and DNA = DNA data"),
                                                  fluidRow(
                                                    column(3,selectInput("aa.nt.col",label=h5("Type of data:"),choices =c("ASN","DNA"))),
                                                    column(3,selectInput("diff",label=h5("Type of plot"),choices =c("compare","plot_one","plot_two"))),
                                                  ),
                                                  fluidRow(
                                                    column(12,div(DT::dataTableOutput("Motif_align"))),
                                                  ),
                                                  plotOutput("Motif_plot_align",height="600px"),
                                                  h4("Exporting plot"),
                                                  fluidRow(
                                                    column(3,numericInput("width_motif_align", "Width of PDF", value=10)),
                                                    column(3,numericInput("height_motif_align", "Height of PDF", value=7)),
                                                    column(3),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_motif_align','Download PDF'))),
                                                  fluidRow(
                                                    column(3,numericInput("width_png_motif_align","Width of PNG", value = 1600)),
                                                    column(3,numericInput("height_png_motif_align","Height of PNG", value = 600)),
                                                    column(3,numericInput("resolution_PNG_motif_align","Resolution of PNG", value = 144)),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_motif_align','Download PNG'))
                                                  )
                                         ),
                                       )
                              ),

                              tabPanel("Diversity and chain usage",
                                       tabsetPanel(
                                         # UI inverse simpson index -----
                                         tabPanel("Inverse Simpson Index",
                                                  p("Inverse simpson index; ∞=infinite diversity and 1=limited diversity"),
                                                  fluidRow(
                                                    column(3,selectInput("group_column_simp",label = h5("Group colum"),
                                                                                choices = c("group","Indiv","Indiv.group"),
                                                                                selected = "Indiv.group")),
                                                           column(3,selectInput("group_column_simp2",label = h5("Unique clone column"),
                                                                                         "")),                            
                                         
                                         
                                         ),
                                                  fluidRow(column(12, div(DT::dataTableOutput("table_display.diversity")))),
                                                  downloadButton('downloadTABLE_simpson.inv','Download table'),
                                                  fluidRow(
                                                    column(3,selectInput("inv.simp_colour.choise",label = h5("Colour"), choices =  c("default","random","grey"))),
                                                    column(3,selectInput("group.index",label = h5("Select x-axis"),
                                                                         choices = simp.index.names,
                                                                         selected = "V2")),
                                                    column(3,selectInput("group2.index",label = h5("Colour by this group"),
                                                                         choices = simp.index.names,
                                                                         selected = "V1")),
                                                    column(3,selectInput("x.axis.index",label = h5("Select x-axis (total or unique clones"),
                                                                         choices = simp.index.names,
                                                                         selected = "total # clones"
                                                                         ))),
                                                  fluidRow(
                                                    column(3,
                                                           wellPanel(id = "tPanel22",style = "overflow-y:scroll; max-height: 400px",
                                                                     uiOutput('myPanel.inv.simp'))),
                                                    column(4,plotOutput("simpson.index1", height="400px")),
                                                    column(4,plotOutput("simpson.index2", height="400px"))),
                                                  fluidRow(
                                                    
                                                    
                                                    column(2, numericInput("conf","confidence of T test", value =0.95, max = 0.99)),
                                                    column(2,selectInput("group1_column",label = h5("Column of group"), 
                                                                         choices = simp.index.names,
                                                                         selected = "V2")),
                                                    column(2,selectInput("group3_column",label = h5("Group selected"),
                                                                         choices=c("group","Indiv","Indiv.group"),selected = "group")),
                                                    column(2,selectInput( "group1_selected",label = h5("Group1"),"" )),
                                                    column(2,selectInput( "group2_selected",label = h5("Group2"),"" )),
                                                    column(2,  selectInput("tail",
                                                                           label = "Please Select a relationship you want to test:",
                                                                           choices = c("Two.tailed" = "two.sided", 
                                                                                       "one.tailed(Less)" = "less",
                                                                                       "one.tailed(Greater)" = "greater")))
                                                  ),
                                                  
                                                  fluidRow(
                                                    column(2,   radioButtons("varequal",
                                                                             "Assume equal variance:",
                                                                             choices = c("Yes" = "y",
                                                                                         "No" = "n"))),
                                                    column(2,radioButtons("paired",
                                                                          "Paired?",
                                                                          choices = c("Yes" = "y",
                                                                                      "No" = "n"))),
                                                    
                                                  ),
                                                  p("The observed t test statistic :"),
                                                  textOutput('tvalue'),
                                                  p("The p-value is"),
                                                  textOutput('pvalue'),
                                                  p("The confidence intervalue is:"),
                                                  textOutput("confidence.int"),
                                                  
                                                  fluidRow(
                                                    column(3,numericInput("width_simpson.inv", "Width of PDF", value=10)),
                                                    column(3,numericInput("height_simpson.inv", "Height of PDF", value=6)),
                                                    column(3),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_simpson.inv','Download PDF'))
                                                  ),
                                                  fluidRow(
                                                    column(3,numericInput("width_png_simpson.inv","Width of PNG", value = 1600)),
                                                    column(3,numericInput("height_png_simpson.inv","Height of PNG", value = 900)),
                                                    column(3,numericInput("resolution_PNG_simpson.inv","Resolution of PNG", value = 144)),
                                                    column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_simpson.inv','Download PNG'))
                                                  )
                                                  
                                                  # gini index is created from Lorentz Surface Calculation pone.0125373.s004.xlsx
                                                  
                                         )
                                         #####   
                                         
                                       )
                                       
                                       
                              )
                              
                            )
                            )
                          )
                 ),
                 # UI Index data graphs -----
                 tabPanel("FACS Index data",
                          sidebarLayout(
                            sidebarPanel(id = "tPanel3",style = "overflow-y:scroll; max-height: 1000px; position:relative;", width=3,
                                         tags$style(type="text/css", "body {padding-top: 80px; padding-left: 10px;}"),
                                         selectInput("dataset3", "FACS file:", choices = c("test-FACS", "own_FACS")),
                                         fileInput('file_FACS', 'Raw index FACS file',
                                                   accept=c('FACS files', '.fcs')),
                                         
                                         selectInput("data_clone.index", "Unsummarised clone file:", choices = c("gd.test.clone" ,"own.clone.file")),
                                         fileInput('file_diversity.index.2', 'Upload unsummarised clone file',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                         selectInput("dataset7", "Merged FACS and clone file for colouring", choices = c("test-csv" ,"own_csv")),
                                         fileInput('file_FACS.csv1', 'FACS+clone file',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                         selectInput("dataset_index.2", "Choose a dataset for complex plot:", choices = c("test-csv" ,"own_csv_file")),
                                         fileInput('file_FACS.csv2', 'File for dot plot',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                         fluidRow(
                                           column(4, numericInput("yintercept",label = h5("y-intercept line"),value = 1000 )),
                                           column(4, numericInput("xintercept",label = h5("x-intercept line"),value = 1000 )),
                                           column(4, selectInput("int.type" ,label = h5("Line type"), choices = c("solid","dotted","dashed")))
                                         ),
                                         fluidRow(
                                           column(4, colourInput("intercept.col",label = h5("Line colour"),value = "grey" )),
                                           column(4, numericInput("min.y",label = h5("Min range (y-axis)"),value = 1 )),
                                           column(4, numericInput("min.x",label = h5("min range (x-axis)"),value = 1 ))),                                        
                                         
                                         
                                         fluidRow(
                                           
                                           column(4, numericInput("max.y",label = h5("Max range (y-axis)"),value = 5 )),
                                           column(4, numericInput("max.x",label = h5("Max range (x-axis)"),value = 5 )),  
                                           column(4, numericInput("leg.dot.size",label = h5("Legend dot size"),value = 5 ))),
                                         
                                         fluidRow(
                                           column(4, numericInput("axis.numeric.size",label = h5("Numeric text size"),value = 28 )),
                                           column(4, numericInput("axis.title.size",label = h5("Label text size"),value = 40 )),
                                           column(4, numericInput("dot.alpha",label = h5("Transparancy of point"),value = 1 )),
                                         ),
                                         
                                         fluidRow(
                                           
                                           
                                           column(4, selectInput("legend.dot",label=h5("Legend location"),choices = c("top","bottom","left","right","none"),selected = "right")),
                                           column(4,numericInput("legend.size.cd","Legend text size",value=12)),
                                           column(4,numericInput("legend.column", "# of legend columns", value=1)),
                                         ),
                                         
                                         
                            ),
                            mainPanel(tabsetPanel(
                              
                              # merging FACS file with clone file -----
                              tabPanel("Combined scTCR with Index data",
                                       div(DT::dataTableOutput("FACS.CSV")),
                                       # div(DT::dataTableOutput("merged.clone")),
                                       div(DT::dataTableOutput("merged.index.clone")),
                                       
                                       textInput("name.colour2","Prefix of file name","ID.780_plate1.section1."),
                                       downloadButton('downloadTABLE_FACS','Download table')
                              ),
                              
                              # UI complex dotplot add columns if needed -----
                              tabPanel("Adding columns for colouring",
                                       
                                       selectInput("string.data","Recommended selecting for ab TCR data: Indiv (or group),TRBV,CDR3b.Sequence, TRBJ, TRAV, CDR3a.Sequence, TRAJ, AJ, BJ and AJBJ. \nDo not select flurochrome columns, clone, cloneCount, LocX or LocY, row or column","",multiple = T, width = "1200px"),
                                       
                                       fluidRow(
                                         column(2,selectInput("group.col.dot",label = h5("Group"),"")),
                                         column(2,selectInput("V.gene.1",label = h5("V Gene 1"),"")),
                                         column(2,selectInput("CDR3.1",label = h5("CDR3 1"),"")),
                                         column(2,selectInput('V.gene.2', label = h5("V Gene 2"), "")),
                                         column(2,selectInput("CDR3.2",label = h5("CDR3 1"),"")),
                                         
                                       ),
                                       fluidRow(
                                         column(6,numericInput("numeric.cloneCount","Filter based on number of times a clone was observed: select 0 for all",value=1))
                                       ),
                                       div(DT::dataTableOutput("table.index.1")),
                                       verbatimTextOutput("NAMES.df"),
                                       
                                       textInput("name.colour","Prefix of file name","ID.780_"),
                                       downloadButton('downloadTABLE_cleaning','Download table')),
                              # UI complex dotplot -----
                              tabPanel("Complex dotplot",
                                       fluidRow(
                                         column(2,selectInput("x.axis2",label = h5("Select x-axis"),"")),
                                         column(2,selectInput("y.axis2",label = h5("Select y-axis"),"")),
                                         column(2,selectInput("density_dotplot",label = h5("Add histogram"), choices = c("no","yes"))),
                                         column(2, selectInput("grid.lines.dot", label = h5("Add gridlines?"), choices = c("no","yes"))),
                                         column(2,selectInput("group_complex_dot",label = h5("Colour by:"),"")),
                                         column(2,selectInput( "FACS.index_colour.choise",label = h5("Colour"),choices = c("default","random","grey"), selected = "random")),
                                         
                                         
                                       ),
                                       
                                       fluidRow(column(3,
                                                       wellPanel(id = "tPanel222",style = "overflow-y:scroll; max-height: 350px",
                                                                 h4("Colour"),
                                                                 uiOutput('myPanel.FACS.index'),
                                                       )),
                                                
                                                column(3,
                                                       wellPanel(id = "tPanel222",style = "overflow-y:scroll; max-height: 350px",
                                                                 h4("Shape"),
                                                                 uiOutput('myPanel.FACS.index.shape')
                                                                 
                                                       )),
                                                
                                                column(3,
                                                       wellPanel(id = "tPanel222",style = "overflow-y:scroll; max-height: 350px",
                                                                 h4("Size"),
                                                                 uiOutput('myPanel.FACS.index.size')
                                                                 
                                                       )),
                                                column(3, tags$img(src = "shape.png", height = "300px"))
                                                
                                                
                                       ),
                                       fluidRow(column(12, plotOutput("dot_plot.complex2",height = "600px"))),
                                       
                                       textInput("name.colour3","Prefix of file name","ID.780_"),
                                       
                                       fluidRow(
                                         column(3,numericInput("width_complex.dotplot", "Width of PDF", value=10)),
                                         column(3,numericInput("height_complex.dotplot", "Height of PDF", value=8)),
                                         column(3),
                                         column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_complex.dotplot','Download PDF'))
                                       ),
                                       fluidRow(
                                         column(3,numericInput("width_png_complex.dotplot","Width of PNG", value = 1600)),
                                         column(3,numericInput("height_png_complex.dotplot","Height of PNG", value = 1200)),
                                         column(3,numericInput("resolution_PNG_complex.dotplot","Resolution of PNG", value = 144)),
                                         column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_complex.dotplot','Download PNG'))
                                       ),
                              )
                            )
                            )
                          )
                 )
)
##### 

server  <- function(input, output, session) {
  # reactive variables -----
  vals <- reactiveValues(Treemap=NULL)
  vals2 <- reactiveValues(Treemap2=NULL)
  vals3 <- reactiveValues(Circular2=NULL)
  vals4 <- reactiveValues(bar.len=NULL)
  vals5 <- reactiveValues(bar.usage=NULL)
  vals6 <- reactiveValues(bar.aa_percent=NULL)
  vals7 <- reactiveValues(basic=NULL)
  vals8 <- reactiveValues(bar.density=NULL)
  vals9 <- reactiveValues(pie=NULL)
  vals10 <- reactiveValues(heatmap_clonal=NULL)
  vals22 <- reactiveValues(Treemap22=NULL)
  options(shiny.sanitize.errors = F)
  output$sessionInfo <- renderPrint({
    print(sessionInfo())
  })
  
  output$video <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/mMkHpiLt_Hg", width = 1080, height = 720)
  })
  
  
  # Tree map ------
  input.data2 <- reactive({switch(input$dataset,"test-data" = test.data2(),"own" = own.data2())})
  test.data2 <- reactive({
    dataframe = read.csv("test-data/Group/paired_unsummarised2021.09.22.csv",header=T) 
  })
  own.data2 <- reactive({
    inFile2 <- input$file2 
    if (is.null(inFile2)) return(NULL)
    
    else {
      dataframe <- read.csv(
        inFile2$datapath,
        header=TRUE,
        sep=input$sep,
        quote=input$quote)}
    
  })
  # simpson calc -----
  
  input.data.diversity.index <- reactive({switch(input$shannon.index,"ab.test.index" = test.data.gd.index.csv(),"own.index" = own.data.index.csv())})
  
  test.data.gd.index.csv <- reactive({
    dataframe = read.csv("test-data/Group/diversity.csv",header=T)   
  })
  
  own.data.index.csv <- reactive({
    inFile7 <- input$file_diversity.index
    if (is.null(inFile7)) return(NULL)
    
    else {
      dataframe <- read.csv(
        inFile7$datapath)}
    
  })
  
  vals11 <- reactiveValues(Simp1=NULL)
  vals12 <- reactiveValues(Simp2=NULL)
  
  #  
  
  observe({
    updateSelectInput(
      session,
      "group_column_simp2",
      choices=names(input.data2()),
      selected = "AJBJ_aCDR3_BJ_bCDR3")

  })
  
  inv.simpson.index <- function() {
    dataframe = input.data2();
    head(dataframe)
    
    df.names <-  dataframe[ , -which(names(dataframe) %in% c("cloneCount","clone"))]
    df1 <- ddply(dataframe,names(df.names) ,numcolwise(sum))
    df1 <- df1[order(df1$cloneCount, decreasing = T),]
    
    names(df1)
    
    df.group <- unique(df1[names(df1) %in% input$group_column_simp])
    names(df.group) <- "V1"
    
    column.length <- length(df.group$V1)
    column.length
    
    df.group2 <- unique(df1[names(df1) %in% input$group_column_simp2])
    names(df.group2) <- "V1"
    
    row.length <- length(df.group2$V1)
    row.length
    
    m = matrix(NA,ncol=column.length, nrow=row.length)
    samps <- df.group$V1
    
    for (j in 1:column.length){
      df2 <- subset(df1,get(input$group_column_simp)==samps[j])
      m[,j] <- c(df2$cloneCount, rep(NA, row.length - length(df2$cloneCount)))
    }
    m <- as.data.frame(m)
    names(m) <- samps
    head(m)
    m
    
    a <- matrix(nrow=1,ncol=dim(m)[2])
    b <- matrix(nrow=1,ncol=dim(m)[2])
    d <- matrix(nrow=1,ncol=dim(m)[2])
    
    for( i in 1:dim(m)[2]) {
      
      samp <- m[,i]
      samp <- na.omit(samp)
      a[,i] <- diversity(samp,"invsimpson")
      b[,i] <- sum(samp)
      d[,i] <- nrow(as.data.frame(samp))
    }
    
    a1 <- rbind(a,b,d)  
    a1 <- as.data.frame(a1)
    names(a1) <- names(m)
    a1 <- rbind(a,b,d)  
    a1 <- as.data.frame(a1)
    names(a1) <- names(m)
    a1
    df_name <- as.data.frame(do.call(rbind, strsplit(as.character(names(m)), "\\.")))
    head(df_name) 
    
    a2 <- as.data.frame(t(a1))
    names(a2) <- c("inv.simpson.index","total # clones","unique # clones")
    a2
    
    both <- cbind(a2,df_name)
    both$Indiv_group <- paste(both$V1,both$V2,sep = "_")
    as.data.frame(both)
    
    
  }
  
  output$table_display.diversity <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    dat <- inv.simpson.index()
    dat <- as.data.frame(dat)
    dat
  })
  
  # other ------
  cols_simp.index <- reactive({
    dat <- inv.simpson.index();
    dat <- as.data.frame(dat)

    selected.col <- dat[names(dat) %in% input$group2.index]
    names(selected.col) <- "V1"
    dat[names(dat) %in% input$group2.index] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))


    num <- unique(dat[names(dat) %in% input$group2.index])
    col.gg <- gg_fill_hue(dim(num)[1])


    if (input$inv.simp_colour.choise == "default") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.inv.simpson", i, sep="_"), paste(num[i,]), col.gg[i])
      })
    }
    else if (input$inv.simp_colour.choise == "random") {
      palette1 <- distinctColorPalette(dim(num)[1])
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.inv.simpson", i, sep="_"), paste(num[i,]), palette1[i])
      })

    }

    else {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.inv.simpson", i, sep="_"), paste(num[i,]), "grey")
      })


    }

  })

  output$myPanel.inv.simp <- renderUI({cols_simp.index()})

  colors_inv.simp <- reactive({
    dat <- inv.simpson.index()
    dat <- as.data.frame(dat)
    selected.col <- dat[names(dat) %in% input$group2.index]
    names(selected.col) <- "V1"
    dat[names(dat) %in% input$group2.index] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))


    num <- unique(dat[names(dat) %in% input$group2.index])
    lapply(1:dim(num)[1], function(i) {
      input[[paste("col.inv.simpson", i, sep="_")]]
    })
  })
  group.diversity1 <- function() {
    both <- inv.simpson.index()
    both <- as.data.frame(both)

    cols <- unlist(colors_inv.simp())

    selected.col <- both[names(both) %in% input$group2.index]
    names(selected.col) <- "V1"
    both[names(both) %in% input$group2.index] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))

    unique.col <- as.data.frame(unique(both[names(both) %in% input$group2.index]))
    names(unique.col) <- "V1"
    unique.col$simp.inv_palette <- cols
    #    df3 <- as.data.frame(merge(both,unique.col,by.x="V1",by.y = "V1"))


    vals11$Simp1 <- ggplot(both,aes(x=get(input$group.index),y=inv.simpson.index))+
      geom_boxplot(show.legend = F)+
      geom_dotplot(aes(fill=get(input$group2.index)),binaxis = 'y',
                   dotsize = 1,
                   stackdir = 'center',show.legend = T) +
      theme_classic() +
      scale_fill_manual(values = c(unique.col$simp.inv_palette)) +
      theme(text=element_text(size=20,family="serif"),
            axis.title = element_text(colour="black", size=20,family="serif"),
            axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
            axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
            axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
            axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
            legend.title  =element_blank(),
            legend.position = "bottom") +

      xlab("")+
      ylab("Inverse simpson index")

    vals11$Simp1

  }
  group.diversity2 <- function() {
    both <- inv.simpson.index()
    cols <- unlist(colors_inv.simp())
    validate(
      need(nrow(both)>0,
           "select correct chain")
    )

    both <- as.data.frame(both)

    selected.col <- both[names(both) %in% input$group2.index]
    names(selected.col) <- "V1"
    both[names(both) %in% input$group2.index] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))

    unique.col <- as.data.frame(unique(both[names(both) %in% input$group2.index]))
    names(unique.col) <- "V1"
    unique.col$simp.inv_palette <- cols


    vals12$Simp2 <- ggplot(both,aes(x=get(input$x.axis.index), y=inv.simpson.index,color=get(input$group2.index)))+
      geom_point(size =3, alpha =1, show.legend =T)+
      # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
      #               limits = c(1,10^6),
      #               labels = trans_format("log10", math_format(10^.x))) +
      theme_bw() +
      scale_color_manual(values=unique.col$simp.inv_palette) +
      #annotation_logticks()  +
      theme(text=element_text(size=20,family="serif"),
            axis.title = element_text(colour="black", size=20,family="serif"),
            axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
            axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
            axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
            axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
            legend.title  =element_blank(),
            legend.position = "bottom") +
      scale_alpha(guide = 'none') +
      labs(x="number of clones",
           y="inverse simpson index")

    vals12$Simp2

  }
  output$simpson.index1 <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    group.diversity1()
  })
  output$simpson.index2 <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    group.diversity2()
  })
  
  table.inv.simpson <- function () {
    dat <- inv.simpson.index()
    dat <- as.data.frame(dat)
    dat
    
  }
  # 
  
  select_group2 <- function () {
    df <- input.data2();
    
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    
    df2 <- as.data.frame(unique(df[names(df) %in% input$group3_column]))
    df2 <- as.data.frame(df2)
    #names(df2) <- "V1"
    df2
  }
  
  observe({
    updateSelectInput(
      session,
      "group1_selected",
      choices=select_group2(),
      selected = "CD8")

  }) # group
  
  
  
  observe({
    updateSelectInput(
      session,
      "group2_selected",
      choices=select_group2(),
      selected = "IFN")

  }) # group
  # 
  ttestout <- reactive({
    dat <- table.inv.simpson()
    conf <- input$conf
    dat <- dat[order(dat$V1),]
    ve <- ifelse(input$varequal == 'y', TRUE, FALSE)
    pair_samp <- ifelse(input$paired == 'y', TRUE, FALSE)
    group1 <- subset(dat, get(input$group1_column)==input$group1_selected) # group 1
    group2 <- subset(dat, get(input$group1_column)==input$group2_selected) # group 2
    t.test(group1$inv.simpson.index, group2$inv.simpson.index, paired = pair_samp, var.equal = ve, alternative = input$tail,conf.level = conf)
  })

  output$tvalue <- renderPrint({
    vals <- ttestout()
    if (is.null(vals)){return(NULL)}
    vals$statistic
  })

  # Output of p value
  output$pvalue <- renderPrint({
    vals <- ttestout()
    if (is.null(vals)){return(NULL)}
    vals$p.value
  })
  output$confidence.int <- renderPrint({
    vals <- ttestout()
    if (is.null(vals)){return(NULL)}
    vals$conf.int
  })
  

  

  
 
  
 
  
  
  
}



shinyApp(ui, server)



