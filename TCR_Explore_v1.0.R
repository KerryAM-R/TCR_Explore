
## volcano plots
require("markdown")
require("rmarkdown")
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
require("motifStack") # pcm2pfm function
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
library("shinyWidgets")
library("showtext")
library("ggseqlogo") # motif plot
library("sangerseqR")
require("scales")

require("umap")
require("fpc")
require("fossil")

library(shinybusy)
library(ggridges)

font_add_google("Gochi Hand", "gochi")
font_add_google("Schoolbell", "bell")
font_add_google("Press Start 2P", "Game")

showtext_auto() 

credentials <- data.frame(
  user = c("shiny"),
  password = c("shiny"),
  stringsAsFactors = FALSE
)

font <- as.data.frame(font_families())
font
names(font) <- "Fonts"
test_fun <- function() {
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

options(shiny.maxRequestSize=200*1024^2)

credentials <- data.frame(
  user = c("shiny", "shinymanager"),
  password = c("azerty", "12345"),
  stringsAsFactors = FALSE
)

# 95% confidence interval
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
middle <- function(x) { r <- quantile(x, probs=c(0.25, 0.25, 0.5, 0.75, 0.75))
names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
r
}
o <- function(x) {
  subset(x, x < quantiles_95(x)[1] | quantiles_95(x)[5] < x)
}

graph_type <- c("histogram","density")
axis_density_group <- c("x-axis","y-axis")

angle <- c(0,90,180,270)

error_message_val1 <- "No data found"
error_message_val2 <- "Uploading file"
error_message_val3 <- "Upload clone file"
error_message_val4 <- "no own list found\n \nSuggest uploading file\nheaders=ID"

simp.index.names <- c("total # clones","unique # clones")
# user interface  ----
ui <- navbarPage(title = tags$img(src = "Logo.png",window_title="TCR_Explore", height = 90, width = 140,
                                  
                                  style = "margin:-35px 10px"
                                  
),

position = "static-top",collapsible = F, 
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
# tutorials -----
tabPanel("Tutorials",
         navlistPanel(id = "Markdown_panel",widths = c(2, 10),
                      tabPanel("Overview",
                               includeMarkdown("README.md"),
                               # tags$video(id="video2", type = "video/mp4",src = "test.mp4", controls = "controls", height="720px")
                      ),
                      tabPanel("Quality control information (includes video tutorial)",
                               h3("Tutorial video of Quality control processes"),
                               uiOutput("video"),
                               fluidRow(includeMarkdown("READMEQC.md")),
                               
                               # tags$video(id="video2", type = "video/mp4",src = "test.mp4", controls = "controls", height="720px")
                      ),     
                      tabPanel("TCR analysis information",
                               includeMarkdown("README.scTCR.md")),
                      
                      tabPanel("Paired TCR with Index data information",
                               includeMarkdown("README.FACS.md")),
                      
                      
                      tabPanel("Video examples",
                               tabsetPanel(
                                 tabPanel("Overview of Pairing",
                                          uiOutput("video2"),
                                 ),
                                 tabPanel("Motif analysis",
                                          uiOutput("video3"),
                                 ),
                                 
                                 tabPanel("Diversity and chain usage",
                                          uiOutput("video4"),
                                 ),
                                 
                                 tabPanel("Overlap",
                                          uiOutput("video5"),
                                 ),
                                 
                                 tabPanel("Paired TCR with Index data",
                                          uiOutput("video6"),
                                 ),
                               ),
                      ),
                      
                      tabPanel("Session info", 
                               tabPanel("Session info", 
                                        div(style="width:800px",
                                            verbatimTextOutput("sessionInfo")),
                                        tags$head(includeHTML(("google-analytics.html")))
                               )
                      )
         )
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
           
# seq to fasta file merger -----
           tabPanel("SEQ to FASTA file merger",
                    fluidPage(
                      sidebarPanel(
                        fileInput("file1_seq.file",
                                  "Choose .seq files from directory",
                                  multiple = TRUE,
                                  accept=c('.seq')),
                        h5("Add Indiv and group/chain name"),
                        h6("IndividualID.groupChain-initialwell"),
                        
                        h4("Select range of 50"),
                        fluidRow(
                          column(6,numericInput("lower.seq","First file",value = 1, step = 50)),
                          column(6,numericInput("upper.seq","Last file",value = 50, step = 50)),
                          
                        ),
                        fluidRow(
                          column(4,selectInput("indiv_miss","Add Indiv label", choices = c("No","Yes"))),
                          column(8,selectInput("group_miss","Add GroupChain-MultipleWells label", choices = c("No","Yes"))),
                          
                        ),
                        
                        
                        fluidRow(column(4,
                                        conditionalPanel(
                                          condition = "input.indiv_miss == 'Yes'",
                                          textInput("indiv.miss.name","Individual ID","Other"),
                                          
                                        )),
                                 column(8,
                                        conditionalPanel(
                                          condition = "input.group_miss == 'Yes'",
                                          textInput("group.miss.name","Group ID","LabA-1"),
                                          
                                        ),
                                        
                                 )
                        ),
                        
                        textInput("seq.name","name of file","test-data"),
                        downloadButton('downloadData_fasta.files', 'Download')
                      ),
                      
                      mainPanel(
                        tableOutput('contents')
                      )
                      
                    )
                    
                    
           ),
# automating the .ab1 QC procecss -----
tabPanel("Automated .ab1 QC", 
         fluidPage(
           sidebarPanel(
            
                              fileInput("file1_ab1.file",
                                        "Choose .ab1 files from directory",
                                        multiple = TRUE,
                                        accept=c('.ab1')),
                              actionButton("go_ab1","Upload or update file"),
                              fluidRow(
                                column(4,selectInput("indiv_miss_ab1","Add Indiv label", choices = c("No","Yes"))),
                                column(8,selectInput("group_miss_ab1","Add GroupChain-MultipleWells label", choices = c("No","Yes"))),
                                
                              ),
                              
                              
                              fluidRow(column(4,
                                              conditionalPanel(
                                                condition = "input.indiv_miss_ab1 == 'Yes'",
                                                textInput("indiv.miss.name_ab1","Individual ID","Other"),
                                                
                                              )),
                                       column(8,
                                              conditionalPanel(
                                                condition = "input.group_miss_ab1 == 'Yes'",
                                                textInput("group.miss.name_ab1","Group ID","LabA-1"),
                                                
                                              ),
                                              
                                       )
                              ),
                              textInput("ab1.name","name of file","test-data"),
                              downloadButton('downloadData_ab1_QC_file','Download single chain file')
            
           ),
         
           mainPanel(
             add_busy_spinner(spin = "fading-circle"),
           tableOutput('contents_ab1')
           
           ),
 
  
         )
  
),
# UI IMGT only ----
           tabPanel("IMGT (Sanger Sequencing)",
                    sidebarLayout(
                      sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                   conditionalPanel(condition="input.QC_panel==1",
                                                    selectInput("dataset_IMGT3", "Choose a dataset:", choices = c("ab-test-data1", "own_data")),
                                                    fileInput('file_IMGT3', 'Select file for IMGT datafile',
                                                              accept=c('xls/xlsx', '.xls')),
                                                    
                                                    fileInput('file_ab1_QC', 'Select file QC ab1 file',
                                                              accept=c('csv', '.csv')),
                                                    conditionalPanel(condition="input.auto_QC==11",
                                                                     downloadButton('downloadTABLE_IMGTonly','Download man QC table'),
                                                    ),
                                                    conditionalPanel(condition="input.auto_QC==12",
                                                                     downloadButton('downloadTABLE_IMGTonly_auto','Download auto QC table')
                                                    ),
                                                    
                                                    
                                   ),
                                    
                                   conditionalPanel(condition="input.QC_panel==2 || input.QC_panel==3 ||input.QC_panel==4",
                                                    
                                                    selectInput("dataset_IMGT_afterQC", "Choose a dataset:", choices = c("ab-test-data1", "own_data1")),
                                                    
                                                    fileInput('file_IMGT_afterQC', 'Completed QC file (.csv)',
                                                              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                                    h5("option for paired and TCRdist outputs"),
                                                    selectInput("IMGT_chain2","Alpha-beta or gamma-delta",choices = c("ab","gd")),
                                                    selectInput("sheet2","Information included", choices = c("Summary+JUNCTION","Summary"))),
                                   conditionalPanel(condition="input.QC_panel==2",
                                                    downloadButton('downloadTABLE.QC1','Download paired chain file')
                                   ),
                                   
                                   conditionalPanel(condition="input.QC_panel==3",
                                                    # textInput("tcr_lab","ID for TCRdist","human_tcr"),
                                                    downloadButton('downloadTABLE.TSV','Download tsv file for TCRdist')
                                   ), 
                                   conditionalPanel(condition = "input.QC_panel==4",
                                                    downloadButton('downloadTABLE.QC1.single.chain','Download single chain file')
                                                    ),
                                   
                                   textInput("IMGT_name_df","File name",""),
                      ),
                      
# main Panel TCR QC -----
                      mainPanel(
                        tabsetPanel(id = "QC_panel",
                                    tabPanel("IMGT create QC file",value = 1,
                                             h4("Fill in the 'clone_quality' column with lowercase: pass or fail"), 
                                             # h4("Add comments if desired"),
                                             fluidRow(column(4, selectInput("sheet","Information included", choices = c("Summary+JUNCTION","Summary"))),
                                                      # column(8, selectInput("include.origin","Include VDJ (n/p) origins (Summary+JUNCTION only)",choices = c("no",'yes'), width = "800px")),
                                             ),
                                             tabsetPanel(id= "auto_QC",
                                               tabPanel("Check uploaded files",value =10,
                                                        tags$head(tags$style("#IMGT2_out  {white-space: nowrap;  }")),
                                                        div(DT::dataTableOutput("IMGT2_out")),
                                                        div(DT::dataTableOutput("IMGT2_out_ab1"))
                                                        
                                                        ),
                                               tabPanel("QC file (manual interrogation)", value = 11,
                                                        div(DT::dataTableOutput("IMGT2_original_out")),
                                               ),
                                               
                                               
                                               tabPanel("QC file with .ab1 scoring", value = 12,
                                                        div(DT::dataTableOutput("IMGT2_out2")),
                                                        )
                                             ),
                                             
                                            
                                            
                                             
                                    ),
# pairing the file -----
                                                  tabPanel("Paired chain file",value = 2,
                                                           tags$head(tags$style("#chain_table_IMGT.QC1  {white-space: nowrap;  }")),
                                                           div(DT::dataTableOutput("Pass.Fail.NA_table")),
                                                           
                                                           
                                                           div(DT::dataTableOutput("chain_table_IMGT.QC1"))
                                                  ),
                                                  
                                                  tabPanel("Single chain file", value = 4,
                                                           div(DT::dataTableOutput("single.chain_table_IMGT.QC1"))
                                                  ),

 
                                                  
                                                  tabPanel("TCRdist output file",value = 3,
                                                           tags$head(tags$style("#chain_table_IMGT.tcrdist  {white-space: nowrap;  }")),
                                                           div(DT::dataTableOutput("chain_table_IMGT.tcrdist")),
                                                           
                                                  ),

# merging the file -----
              
                                                  
                                      )
                                    )
                                    
                                  )
                         ),

# .ab1 chromatogram file -----
           tabPanel("Check .ab1 files",
                    sidebarLayout(
                      sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                   selectInput("dataset_.ab1", "Choose a dataset:", choices = c(".ab1-test-data", ".ab1-own_data")),
                                   fileInput('file_.ab1', 'Chromatogram .ab1 file',
                                             accept=c('.ab1')),
                                   
                      ),
                      mainPanel(
                        tabsetPanel(
                          
                          tabPanel("Chromatogram",
                                   p("A high quality chormatogram will display few mismatches (Blue) between the primary and secondary sequences."),
                                   p("Showcasing the heterogeneous sequences in the .ab1 file"),
                                   fluidRow(
                                     column(3,numericInput("Number.seq.line","Number of sequences per row",value = 100)),
                                     column(3,numericInput("trim5.seq","Trim 5` sequences",value = 0)),
                                     column(3,numericInput("trim3.seq","Trim 3` sequences",value = 0)),
                                     
                                     
                                   ),
                                   plotOutput("chromatogram.seq", height="600px"),
                                   
                                   
                                   fluidRow(
                                     column(3,numericInput("width_chromatogram.seq", "Width of PDF", value=10)),
                                     column(3,numericInput("height_chromatogram.seq", "Height of PDF", value=8)),
                                     column(3),
                                     column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_chromatogram.seq','Download PDF'))
                                   ),
                                   
                                   fluidRow(
                                     column(3,numericInput("width_png_chromatogram.seq","Width of PNG", value = 1600)),
                                     column(3,numericInput("height_png_chromatogram.seq","Height of PNG", value = 1200)),
                                     column(3,numericInput("resolution_PNG_chromatogram.seq","Resolution of PNG", value = 144)),
                                     column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_chromatogram.seq','Download PNG'))
                                   ),
                                   
                                   
                                   
                          ),
                          tabPanel("Primary and secondary sequence aligment",
                                   p("Overlap of the primary and secondary sequence."),
                                   
                                   verbatimTextOutput("alignment"),
                                   
                          ),
                          
                          tabPanel("Chromatogram sequence",
                                   p("Checking for heterozygosity in sequence overlap."),
                                   
                                   p("Check heterogenatiy of sequences with a 0.33 ratio cut-off as per the 'sangerseqR' package recommendation"),
                                   verbatimTextOutput("hetsangerseq"),
                          ),
                        )
                      )
                    ),
                    
           ),
# ImmunoSEQ ====

tabPanel("Convert to TCR_Explore file format",
         sidebarLayout(
           sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                        selectInput("dataset_TSV.Immunoseq", "Choose a dataset:", choices = c("Demo.gd.test-data", "Immunoseq.own-data")),
                        fileInput('file_TSV.Immunoseq', 'File to upload',
                                  accept=c('.tsv',".csv",".txt")),
                        
                        radioButtons("sep.imm", "Separator",
                                     choices = c(Comma = ",",
                                                 Semicolon = ";",
                                                 Tab = "\t"),
                                     selected = "\t"),
                        
                        
                        radioButtons("quote.imm", "Quote",
                                     choices = c(None = "",
                                                 "Double Quote" = '"',
                                                 "Single Quote" = "'"),
                                     selected = '"'),
                        textInput("group.imm","Group ID","Group"),
                        textInput("indiv.imm","Individual ID","ID"),

                        downloadButton('downloadTABLE.Immunoseq','Download filtered table')
                        
                        
           ),
           mainPanel(
             tabsetPanel(
               tabPanel("Converting to TCR_Explore",
               h5("Upload eiter .tsv, .csv or .txt files to convert to TCR_Explore format"),
               p("Rows with missing sequences are removed from V and J gene columns"),
               p("If using ImmunoSEQ data, there is an additional filtering step to only keep in-frame sequences"),
               p(" "),
               
               fluidRow(
                 column(3, selectInput("datasource","Input type",choices = c("ImmunoSEQ","MiXCR","10x_scSeq","Other"))),
                        
               ),
               
               fluidRow(
                 column(4, selectInput("countcolumn","Count column",choices = "")),
                 column(4, selectInput("inframe_immseq","In-frame column (e.g. sequenceStatus)","")),
                 column(4, selectInput("CDR3.gene.clean","CDR3 amino acid column",choices = "")),
               ),
               

               
               fluidRow(
                    column(3, selectInput("D_chain_present","D chain present?",choices = c("yes","no"))),
                    column(3, selectInput("V.GENE.clean","Variable gene column",choices = "")),
                    conditionalPanel("input.D_chain_present == 'yes'",
                                     column(3, selectInput("D.GENE.clean","Diversity gene column",choices = ""))
                        ),
                    column(3, selectInput("J.GENE.clean","Junction gene column",choices = "")),
                        
                        
                        # conditionalPanel(condition = "input.datasource == '10x_scSeq'",
                        #                  column(3,selectInput("V.GENE.clean2","Variable gene column (2)",choices = "")),
                        #                  column(3,selectInput("J.GENE.clean2","Junction gene column (2)",choices = "")),
                        #                  column(3, selectInput("CDR3.gene.clean2","CDR3 amino acid column",choices = ""))
                        #                  )
                        
                        ),
               fluidRow(column(12, selectInput("col.to.remove","Columns to remove","",multiple = T, width = "1200px") )),
               div(DT::dataTableOutput("ImmunoSeq.table"))
             ),
             tabPanel("Video of conversion process",
                      uiOutput("video7"),
             ),
             
             ),
             
           ),
         ),
),


           
           
),
#### TCRex merge files ----
tabPanel("TCR_Explore merge",
         sidebarLayout(
           sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                        # selectInput("dataset2", "Choose a dataset:", choices = c("test_data_clusTCR2","own_data_clusTCR2")),
                        fileInput('file2_TCRexMerge', 'Select files to merge',
                                  multiple = TRUE,
                                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                        downloadButton('downloaddf_TCRexFiltered','Download table')
           ),
           mainPanel(
             tabsetPanel(id = "TCRex_tabs",
                         tabPanel("Merge Multiple Files",value = "merge",
                                  div(DT::dataTableOutput("DEx_TCRexFiltered")),
                                  # downloadButton('downloaddf_multiple_ClusTCR2','Download table')
                         ),
                         
                         
             )
           )
         )
),

# UI TCR plots ----
tabPanel("TCR analysis",
         
         tags$style(HTML("
    .tabbable > .nav > li > a                  {background-color: white;  color:black}
    .tabbable > .nav > li[class=active]    > a {background-color: darkred; color:white}
  ")),
         
         sidebarLayout(
           sidebarPanel(id = "tPanel",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                        # tags$style(type="text/css", "body {padding-top: 80px; padding-left: 10px;}"),
                        #textInput(inputId = "lab1", label = "Group label of file 1",value = "Ex.vivo"),
                        tags$head(tags$style(HTML(".shiny-notification {position:fixed;top: 50%;left: 30%;right: 30%;}"))),
                        tags$head(tags$style(HTML('.progress-bar {background-color: purple;}'))),
                        selectInput("dataset", "Choose a dataset:", choices = c("ab-test-data2",
                                                                                # "ImmunoSEQ-test-data", 
                                                                                "own_data2")),
                        fileInput('file2', 'Select file for single samples',
                                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                        
                        fluidRow(
                          column(6,radioButtons('sep', 'Separator', c( Tab='\t', Comma=','), ',')),
                          column(6,radioButtons('quote', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), '"'))
                        ),
                        selectInput("group_column",label = h4("Column of group"), ""),
                        selectInput("Filter_NGS","Filtering",choices = c("no","yes")),
                        fluidRow(
                         
                          column(6, selectInput("clonotypes_column","Clonotypes (e.g. vdj_CDR3_AG_BD)","")),
                          column(6, selectInput("group_column_first","Summarise by (e.g. Indiv)","")),
                          
                        ),
                        fluidRow(
                          column(6, selectInput("filter_by_cumsum_or_order","Filter by either:",c("cummulative_freq","Top_clontypes"))),
                          conditionalPanel("input.filter_by_cumsum_or_order == 'cummulative_freq'",
                                           column(6,numericInput("filter_below_cumsum","Display clonotypes below frequency:",value = 0.5, min = 0.1,max=1,step = 0.01)),
                                           ),
                          conditionalPanel("input.filter_by_cumsum_or_order == 'Top_clontypes'",
                                           column(6,numericInput("top_clones_to_include","Display clonotypes below frequency:",value = 100, min = 1)),
                          ),
                          
                        ),
                        
                        colourInput("one.colour.default","One colour","grey"),
                        
                        selectInput("type.tree",label = h4("Type of input"), choices =  c("raw data","Summarised data")),
                        
                        selectInput("font_type",label = h4("Type of font"),choices = font,selected = "serif"),
                        downloadButton("table_length","Download summarised table with length"),
                        
                        tags$hr()
           ),
           
           mainPanel(tabsetPanel(
             tabPanel("Overview of TCR pairing",tabsetPanel(
# UI Summary table -----
               tabPanel("Summary table",
                        # verbatimTextOutput("names.in.file3"),
                        fluidRow(
                          column(3,selectInput("type.chain","Alpha-beta or gamma-delta",choices = c("ab","gd","NGS_immunoseq"))),
                          column(3,selectInput("type.of.graph", "Summary table output",choices = c("general summary","TCRdist3")))
                        ),
                        fluidRow(column(12, selectInput("string.data3","column names for summary","",multiple = T, width = "1200px") )),
                        tags$head(tags$style("#chain_table_IMGT.QC3  {white-space: nowrap;  }")),
                        div(DT::dataTableOutput("chain_table_IMGT.QC3")),
                        downloadButton('downloadTABLE.QC3','Download table')
                        
               ),
tabPanel("Filtered table",
         div(DT::dataTableOutput("Test_table")),
         ),
# UI Treemap -----
               tabPanel("Treemap",
                        fluidRow(
                          
                          column(3,  numericInput("nrow.tree",label = h5("Rows"), value = 1))
                        ),
                        
                        fluidRow(
                          column(6, selectInput("string.data.tree.order","Order of group in graph",choices = "",multiple = T, width = "600px")),
                          column(2, colourInput("strip.colour.tree","Strip colour",value = "lightgrey")),
                          column(2,  colourInput("strip.text.colour.tree","Text strip colour",value = "black")),
                          column(2, numericInput("panel.text.size.tree","Size of panel text", value = 20))
                        ),
                        
                        fluidRow( 
                          column(2, selectInput("tree_colour.choise",label = h5("Colour"), choices =  c("default","rainbow","random","one colour"))),
                          column(2, selectInput("fill2",label = h5("Colour treemap by"),"" )),
                          column(2, selectInput("sub_group2",label = h5("Separate panels by"),"" )),
                          # column(3,selectInput( "wrap",label = h5("Group"),"" )),
                          column(2,selectInput( "count2",label = h5("Count column"),"")),
                          column(2, selectInput("tree.lab",label = h5 ("Add label"),choices = c("yes","no"))),
                          column(2, colourInput("Treemap.text.colour","Text colour",value = "black"))
                          
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
                        h5("If you see this error: 'not enough space for cells at track index '1'. 
                                           Adjust Text size (cex)"),
                        p(" "),
                        fluidRow(
                          
                          column(2,selectInput( "group_selected2",label = h5("Group"),"" )),
                          column(2,selectInput( "chain1",label = h5("Chain one"),"" )),
                          column(2,selectInput( "chain2",label = h5("Chain two"),"" )),
                          column(2,style = "margin-top: 15px;", sliderInput("chord.transparancy","Transparancy",value = 0.5,step = 0.05, min=0,max=1)),
                          column(2,style = "margin-top: 15px;", numericInput("CHORD.cex","Text size (cex)",value = 1, min=0,step = 0.05)),
                          # column(2,tableOutput("table_display")),
                          
                        ),
                        fluidRow(
                          column(2, selectInput("circ_lab",
                                                label = h5("Type of label"),
                                                choices = c("Label","colour selected clone/s (label)","colour selected clone/s (no label)","no labels"))),
                          column(2,selectInput( "colour_cir",label = h5("Colour"),choices = c("default","rainbow","random","one colour"))),  
                          column(2,style = "margin-top: 15px;", numericInput("seed.numb.chord","Random colour generator",value = 123)),
                        ),
                        
                        conditionalPanel(
                          condition = "input.circ_lab == 'colour selected clone/s (label)' || input.circ_lab == 'colour selected clone/s (no label)'",
                          selectInput("string.data.circ.order","Chains to highlight",choices = "",multiple = T,width = "800"),
                          
                          fluidRow(
                            column(2, colourInput("colour.chord.line","Line colour","black")),
                            column(2, sliderInput("line.chord.type","Line type (0 = no line)",min=0,max=6,value=1)),
                            column(2,numericInput("thickness.chord.line","Thickness of line", value = 2)),
                            column(2, sliderInput("unselected.chord.transparacy","Transparancy unselected",min=0,max=1,value=0.1,step = 0.05)),
                            column(2, sliderInput("selected.chord.transparacy","Transparancy selected",min=0,max=1,value=1,step = 0.05)),
                          ),
                        ),
                        
                        # plotOutput("colour.trans.test"),
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
                          column(3,numericInput("width_png_circ","Width of PNG", value = 1200)),
                          column(3,numericInput("height_png_circ","Height of PNG", value = 1200)),
                          column(3,numericInput("resolution_PNG_circ","Resolution of PNG", value = 144)),
                          column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_circ','Download PNG'))
                        ),
                        # tableOutput("out.col.table1")
                        
               ),
# UI Pie ----
               tabPanel("Pie chart",
                        fluidRow(column(2,selectInput("pie_chain",label = h5("Colour by this chain"),"")),
                                 column(2,selectInput("pie_colour.choise",label = h5("Colour"), choices =  c("default",'rainbow',"random","one colour"), selected = "random")),
                                 column(2, selectInput("cir.legend",label=h5("Legend location"),choices = c("top","bottom","left","right","none"),selected = "none")),
                                 column(2,  numericInput("nrow.pie",label = h5("Rows"), value = 1)),
                                 column(2,  numericInput("size.circ",label = h5("Size of legend text"), value = 6))
                                 
                        ),
                        fluidRow(
                          column(6, selectInput("string.data.pie.order","Order of group in graph",choices = "",multiple = T, width = "600px")),
                          column(2, colourInput("strip.colour.pie","Strip colour",value = "lightgrey")),
                          column(2,  colourInput("strip.text.colour.pie","Text strip colour",value = "black")),
                          column(2, numericInput("panel.text.size.pie","Size of panel text", value = 20))
                        ),
                        fluidRow(column(3,
                                        wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                  uiOutput('myPanel_pie'))),
                                 column(9, plotOutput("pie_out",height="600px"))),
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
                      # p("This section contains 4 tabs for motif analysis"),
                      tabsetPanel(
# UI CDR3 length distribution graphs ----- 
                        tabPanel("CDR3 length distribution",
                                 p(" "),
                                 h6("The amino acid CDR3  columns are callled: AA.JUNCTION, JUNCTION..AA. or CDR3_IMGT."),
                                 h6("The _A (alpha), _B (beta), _G (gamma), _D (delta)"),
                                 
                                 fluidRow(
                                   column(2,selectInput('graph_type', 'Type of graph', graph_type)),
                                   column(2,selectInput( "aa.or.nt","CDR3 length column","" )),
                                   
                                   conditionalPanel(
                                     condition = "input.graph_type == 'histogram'",
                                     column(2,selectInput( "selected_group_len","Group","" )),
                                     column(2, selectInput("sep_hist_group","Wrap by group",choices = c("yes","no"), selected = "no")),
                                     column(2,selectInput("chain.hist.col","Colour by:",""))),
                                   
                                   conditionalPanel(
                                     condition = "input.graph_type == 'density'",
                                     column(3,sliderInput("alpha.density","Transparency",min=0, max = 1,value = 0.25, step = 0.05)),
                                     column(2, selectInput("facet_by_group","Wrap by group",choices = c("yes","no"), selected = "no")),
                                     
                                     
                                   ),
                                   
                                 ),
                                 fluidRow(
                                   column(2, selectInput("hist.density.legend","Legend location",choices = c("top","bottom","left","right","none"),selected = "right")),
                                   column(2, numericInput("col.num.CDR3len","# of Legend columns",value = 3)),
                                   column(2,numericInput("legend.text.hist","Legend text size",value = 6)),
                                   column(2,selectInput("hist_colour.choise","Colour", choices =  c("default","rainbow","random","grey")))
                                 ),
                                 
                                 fluidRow(
                                   column(2,numericInput("hist.text.sizer","Size of #",value=16)),
                                   column(2,numericInput("hist.text.sizer2","Axis text size",value=30)),
                                   column(2, numericInput("xlow","x-axis (min)",value=0)),
                                   column(2, numericInput("xhigh","x-axis (max)",value=30)),
                                   column(2, numericInput("xbreaks","x-axis tick marks",value=5)),
                                   column(2, numericInput("ybreaks","y-axis tick marks",value=2)),
                                 ),
                                 
                                 
                                 fluidRow(
                                   conditionalPanel(
                                     condition = "input.graph_type == 'histogram'",
                                     
                                     column(3,
                                            wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                      uiOutput('myPanel.hist')))),
                                   
                                   conditionalPanel(
                                     condition = "input.graph_type == 'density'",
                                     
                                     column(3,
                                            wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                      uiOutput('myPanel.hist2')))),
                                   
                                   
                                   column(9, plotOutput("Chain1_length",height="600px"))),
                                 
                                 conditionalPanel(
                                   condition = "input.graph_type == 'histogram'",
                                   div(DT::dataTableOutput("hist.table")),
                                   
                                 ),
                                 
                                 fluidRow(
                                   column(3,numericInput("width_length", "Width of PDF", value=10)),
                                   column(3,numericInput("height_length", "Height of PDF", value=4)),
                                   column(3),
                                   column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_length','Download PDF'))
                                 ),
                                 fluidRow(
                                   column(3,numericInput("width_png_length","Width of PNG", value = 1600)),
                                   column(3,numericInput("height_png_length","Height of PNG", value = 600)),
                                   column(3,numericInput("resolution_PNG_length","Resolution of PNG", value = 144)),
                                   column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_length','Download PNG'))
                                 ),
                        ),
# UI motif -----
                        tabPanel("Motif (amino acid)",
                                 p(" "),
                                 h6("The amino acid CDR3 columns are callled: AA.JUNCTION, JUNCTION..AA. or CDR3_IMGT."),
                                 h6("The _A (alpha), _B (beta), _G (gamma), _D (delta)"),
                                 h5("Select amino acid column and CDR3 length"),
                                 verbatimTextOutput("length"),
                                 fluidRow(
                                   column(6,div(DT::dataTableOutput("length.table"))),
                                   column(6,div(DT::dataTableOutput("Motif"))),
                                 ),
                                 
                                 fluidRow(
                                   column(2, checkboxInput("compar.lab.motif.aa.single","Add Label",value = T)),
                                   column(2,selectInput( "aa.or.nt2",label = h5("Amino acid CDR3 column"),"" )),
                                   column(2,style = "margin-top: 15px;",numericInput("len","CDR3 amino acid length", value = 15)),                               
                                   column(2,selectInput( "group_selected_motif",label = h5("Group 1 (top)"),"" )),
                                   column(2,selectInput( "group_selected_motif2",label = h5("Group 2 (bottom)"),"" )),
                                   column(2, selectInput("comarpison.aa.motif",label = h5("Type of comparison"), choices= c("single.group1","compare two groups"))),
                                   column(2,style = "margin-top: 15px;",numericInput("font_size_motif","Font size", value = 6)),  
                                   column(2,style = "margin-top: 15px;",numericInput("font_legend_motif","Legend size", value = 12)),  
                                   
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
                                 p(" "),
                                 h6("The amino acid CDR3  columns are callled: AA.JUNCTION, JUNCTION..AA. or CDR3_IMGT."),
                                 h6("The _A (alpha), _B (beta), _G (gamma), _D (delta)"),
                                 fluidRow(
                                   column(3,selectInput("restricted.length.range",label = h5("Restrict range"), choices = c("no","yes"))),
                                   column(3,selectInput("group_selected_three",label = h5("select length"),"", multiple = T)),
                                 ),
                                 fluidRow(
                                   column(3,selectInput("aa.or.nt4",label = h5("CDR3 column"),"")),
                                   column(3,selectInput("group_selected_one",label = h5("First group (top of plot)"),"" )),
                                   column(3,selectInput("group_selected_two",label = h5("Second group (bottom of plot)"),"" )),
                                   
                                   # column(2,checkboxInput("",select = T))
                                 ),
                                 p("ASN = amino acid data and DNA = DNA data"),
                                 fluidRow(
                                   column(3,selectInput("aa.nt.col",label=h5("Type of data:"),choices =c("ASN","DNA"))),
                                   column(3,selectInput("diff",label=h5("Type of plot"),choices =c("compare","plot_one","plot_two"))),
                                   column(3, checkboxInput("compar.lab.motif.all","Add Label (compare only)",value = T)),
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
             
# diversity and chain usage -----
             
             tabPanel("Diversity and chain usage",
                      tabsetPanel(
# UI bar graphs ----- 
                        tabPanel("Chain bar graph",
                                 fluidRow(
                                   column(3,selectInput("stat",label = h5("Plot output"),choices=c("chains","frequency","stacked")),),
                                   column(3,selectInput( "selected_group_chain",label = h5("Group"),"" )),
                                   column(3,selectInput( "variable_chain",label = h5("Select y-axis"),"" )),
                                 ),
                                 
                                 
                                 
                                 
# chain usage -----
                                 conditionalPanel(
                                   condition = "input.stat == 'chains'",
                                   h5("Individual chains"),
                                   fluidRow(
                                     
                                     column(3,selectInput( "graph_bar_type",label = h5("Select x-axis"),choices = c("count","percentage"))),
                                     column(3,style = "margin-top: 10px;", numericInput("bar.numeric.size","Size of axis label", value = 12)),
                                     column(3,style = "margin-top: 10px;", colourInput("colour_bar.usage","Colour of bars", value = "black"))),
                                   
                                 ),
# cummulative freq  -----
                                 conditionalPanel(
                                   condition = "input.stat == 'frequency'",
                                   h5("Cummulative frequency"),
                                   fluidRow(
                                     column(3,numericInput("numeric.adjust","Adjust # clones label",value=-1)),
                                     column(3, colourInput("colour.numeric.bar","Colour numeric", value = "black")),
                                     column(3, numericInput("label.size","Size of numeric label", value = 6)),
                                     column(3, numericInput("label.size.axis","Size of axis label", value = 20)),
                                   ),
                                   
                                 ),
# stacked bar graph  -----
                                 conditionalPanel(
                                   
                                   condition = "input.stat == 'stacked'",
                                   h5("Stacked bar plot"),
                                   selectInput("string.data2","Order of group in graph",choices = "",multiple = T, width = "1200px"),
                                   fluidRow(
                                     column(3, numericInput("label.size.axis2","Size of axis label", value = 20)),
                                     column(3,selectInput( "bar.stacked_colour.choise",label = h5("Colour"),choices = c(
                                       "default","rainbow","random","grey"))),
                                     fluidRow(
                                       column(2,numericInput("bar.stack.angle","Angle of text",value = 90)),
                                       column(2,numericInput("hight.bar.stack.adj","Position of text",value = 0)),
                                       column(2,selectInput("lines.bar.graph","Display black lines?",
                                                            choices = c("yes","no"),
                                                            selected = "no"))
                                     ),
                                     fluidRow(column(3,numericInput("stacked.no.legend","legend columns",value = 3)),
                                              column(3, selectInput("stacked.legend",label=h5("Legend location"),choices = c("top","bottom","left","right","none"),selected = "right")),
                                              column(3,numericInput("stacked.legend.size","Legend text size",value = 12)),
                                     ),
                                     
                                     
                                   ),
                                 ),
                                 
                                 fluidRow(
                                   
                                   conditionalPanel(
                                     condition = "input.stat == 'stacked'",
                                     column(3,
                                            wellPanel(id = "tPanel22",style = "overflow-y:scroll; max-height: 400px",
                                                      uiOutput('myPanel_cols_stacked_bar'))),
                                   ), 
                                   column(9,plotOutput("Chain1_usage",height="400px")),
                                 ),
                                 
                                 
                                 
                                 
                                 
                                 fluidRow(
                                   column(3,numericInput("width_chain.usage", "Width of PDF", value=10)),
                                   column(3,numericInput("height_chain.usage", "Height of PDF", value=8)),
                                   column(3),
                                   column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_chain.usage','Download PDF'))
                                 ),
                                 fluidRow(
                                   column(3,numericInput("width_png_chain.usage","Width of PNG", value = 1600)),
                                   column(3,numericInput("height_png_chain.usage","Height of PNG", value = 1200)),
                                   column(3,numericInput("resolution_PNG_chain.usage","Resolution of PNG", value = 144)),
                                   column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_chain.usage','Download PNG'))
                                 ),
                                 
                                 
                        ),
                        
                        
# UI inverse simpson index -----
                        tabPanel("Diversity Index",
                                 p("Inverse Simpson Diversity Index: ∞=infinite diversity and 1=limited diversity"),
                                 conditionalPanel(
                                       condition = "input.QC_panel_Simp == 2 || input.QC_panel_Simp == 3",
                                       fluidRow(column(3,selectInput("index.type",label = h5("Type of Diversity"), 
                                                                     choices =  c("shannon_div", "Pielou_even", "inverse_simp", "inverse_simp_corrected","chao1_div"), selected = "Pielou_even"),),
                                                column(3, selectInput("group1_column",label = h5("Column of group"), 
                                                                      ""))),
                                       
                                       
                                      
                                       ),
                                   
                                
                                 tabsetPanel(id = "QC_panel_Simp",
                                   tabPanel("Table and selecting groups",value =1,
                                   fluidRow(
                                     column(3,selectInput("group_column_simp",label = h5("First ID column"),
                                                          "")),
                                     column(3,selectInput("group_column_simp3",label = h5("Second ID column"),
                                                          "")),
                                     column(3,selectInput("group_column_simp2",label = h5("Unique clone column"),
                                                          "")),                            
                                   ),
                                   fluidRow(column(12, div(DT::dataTableOutput("table_display.diversity")))),
                                   downloadButton('downloadTABLE_simpson.inv','Download table'),
                                 ),
                                 
                                 tabPanel("Graph",value =2,
                                          fluidRow(
                                            column(3,selectInput("inv.simp_colour.choise",label = h5("Colour"), choices =  c("default","random","grey"))),
                                            
                                            column(2,selectInput("group.index",label = h5("x-axis group"),
                                                                 "")),
                                            
                                            column(2,selectInput("group2.index",label = h5("Colour by this group"),
                                                                 "")),
                                            
                                            
                                            column(2,selectInput("x.axis.index",label = h5("Select x-axis (total or unique clones"),
                                                                 choices = simp.index.names,
                                                                 selected = "total # clones"
                                            )),
                                            column(3,checkboxInput("scale_x_continuous_x","Change to scientific values",value = F)),
                                            
                                            ),
                                          # fluidRow(
                                          # column(3,selectInput("IQR_95","Select type of error bars",choices = c("IQR","95% CI"))),
                                          #   conditionalPanel(
                                          #     condition = "IQR_95 == '95% CI'",
                                          #   column(3,selectInput("show_all_dots","Type of 95%CI", choices= c("Show outliers only","Colour all dots"))))
                                          # ),
                                          fluidRow(
                                            column(3,numericInput("col.num.simp",label = h5("Legend columns"),value = 1)),
                                            column(3, selectInput("legend.placement.simp",label=h5("Legend location"),choices = c("top","bottom","left","right","none"),selected = "right")),
                                            column(3,numericInput("legend.text.simp",label = h5("Legend text size"),value = 12)),
                                          ),
                                          fluidRow(
                                            column(3,
                                                   wellPanel(id = "tPanel22",style = "overflow-y:scroll; max-height: 400px",
                                                             uiOutput('myPanel.inv.simp'))),
                                            column(4,plotOutput("simpson.index1", height="400px")),
                                            column(4,plotOutput("simpson.index2", height="400px"))
                                          ),
                                          
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
                                          
                                          ),
                                 tabPanel("Statistics",value =3,
                                         
                                          tabsetPanel(
                                            tabPanel("T-test",
                                                     fluidRow(

                                                       column(2, numericInput("conf","confidence of T test", value =0.95, max = 0.99)),
                                                       
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
                                                       column(3, selectInput( "test_ttest","Type of test",choices= c("parametric","non-parametric")))
                                                       
                                                     ),
                                                     
                                                     verbatimTextOutput('tvalue'),
                                                     
                                                     ),
                                            tabPanel("ANOVA",
                                                     verbatimTextOutput('aov_test'),
                                                     verbatimTextOutput('postaov_test')
                                                     )
                                          ),
                                          

                                          
                                          
                                          
                                          )
                                 ),



      
                                
                                 
                                 # gini index is created from Lorentz Surface Calculation pone.0125373.s004.xlsx
                                 
                        )
                        
                        
                      )
                      
                      
             ),
             
# Overlap ----          
             tabPanel("Overlap",
                      tabsetPanel(
# UI heatmap -----
                        tabPanel("Heatmap",
                                 selectInput("group_hm", "Select specific groups", choices = c("yes", "no")),
                                 fluidRow(
                                   column(3,selectInput("group_selected3",label = h5("Select group"),"" )),
                                   column(3,selectInput( "heatmap_2",label = h5("Select x-axis"),"" )),
                                   column(3,selectInput("group.heatmap",label = h5("Select y-axis"),"" )),
                                   column(3, colourInput("col.heatmap",label = h5("Colour"), value = "red"))
                                 ),
                                 fluidRow(
                                   column(3,numericInput("heat.font.size.row","Font size (row)",value=8)),
                                   column(3,numericInput("heat.font.size.col","Font size (col)",value=8)),
                                   
                                 ),
                                 
                                 plotOutput("heatmap_out2",height="800px"),
                                 fluidRow(
                                   
                                   column(3,numericInput("width_heatmap", "Width of PDF", value=10)),
                                   column(3,numericInput("height_heatmap", "Height of PDF", value=8)),
                                   column(3),
                                   column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_heatmap','Download PDF'))
                                 ),
                                 fluidRow(
                                   column(3,numericInput("width_png_heatmap","Width of PNG", value = 1600)),
                                   column(3,numericInput("height_png_heatmap","Height of PNG", value = 1200)),
                                   column(3,numericInput("resolution_PNG_heatmap","Resolution of PNG", value = 144)),
                                   column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_heatmap','Download PNG'))
                                 ),
                        ),
# upset plot -----
                        tabPanel("Upset plot",
                                 fluidRow(
                                   column(3,selectInput("upset.select",label = h5("Select chain"), choices = "", selected = "")),
                                   column(3,selectInput("upset.group.select",label = h5("Group column (max 31 groups)"), choices = "",selected= "")),
                                   column(3, selectInput("upset_anno","Colouring options",choices = c("default", "Colour by degree"))),
                                 ),
                                 
                                 selectInput("order.of.group",label = h5("Group column (max 31 groups)"), choices = "",selected= "", multiple = T, width = "1200px"),
                                 
                                fluidRow(
                                  
                                  column(3, colourInput("right_annotation_colour","Colour of bar (right)",value = "black")),
                                  column(3, colourInput("top_annotation_colour","Colour of bar (top)",value = "black")),
                                  conditionalPanel(condition = "input.upset_anno == 'Colour by degree'",
                                                   column(6,textInput("upset_colours_list1","Colour of dots by degree", value="darkgreen,purple")),
                                                   
                                  ),
                                ),
                                 
                                 fluidRow(
                                   column(2, numericInput("upset.point.size","Size of point",value = 5)),
                                   column(2, numericInput("upset.lwd","Line width",value = 2)),
                                   column(3,numericInput("upset.text.size","Size of text",value = 20)),
                                   column(3,numericInput("upset.font.size","Size of number (top)",value = 12)),
                                   column(3,numericInput("font.size.anno.upset","Size of number (right)",value = 12)),
                                 ),
                                 
                                 plotOutput("UpSet.plot", height = "600px"),
                                 fluidRow(
                                   column(3,numericInput("width_upset", "Width of PDF", value=10)),
                                   column(3,numericInput("height_upset", "Height of PDF", value=8)),
                                   column(3),
                                   column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_upset','Download PDF'))
                                 ),
                                 
                                 fluidRow(
                                   column(3,numericInput("width_png_upset","Width of PNG", value = 1600)),
                                   column(3,numericInput("height_png_upset","Height of PNG", value = 1200)),
                                   column(3,numericInput("resolution_PNG_upset","Resolution of PNG", value = 144)),
                                   column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_upset','Download PNG'))
                                 ),
                                 tags$head(tags$style("#upset.datatable  {white-space: nowrap;  }")),
                                 div(DT::dataTableOutput("upset.datatable")),
                        )
                      ))
             
           )
           )
         )
),
# UI Index data graphs -----
tabPanel("Paired TCR with Index data",
         sidebarLayout(
           sidebarPanel(id = "tPanel3",style = "overflow-y:scroll; max-height: 850px; position:relative;", width=3,
                        # tags$style(type="text/css", "body {padding-top: 80px; padding-left: 10px;}"),
                        
                        conditionalPanel(condition="input.tabselected==1",
                                         selectInput("dataset3", "FACS file:", choices = c("test-FACS", "own_FACS")),
                                         fileInput('file_FACS', 'Raw index FACS file',
                                                   accept=c('FACS files', '.fcs')),
                                         
                                         selectInput("data_clone.index", "Unsummarised clone file:", choices = c("ab.test.clone3" ,"own.clone.file")),
                                         fileInput('file_diversity.index.2', 'Upload unsummarised clone file',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                         textInput("name.colour2","Prefix of file name","ID.780_plate1.section1."),
                                         downloadButton('downloadTABLE_FACS','Download table')
                        ),
# tab panel 2 (Other QC steps) ------
                        conditionalPanel(condition="input.tabselected==2",
                                         selectInput("dataset7", "Merged FACS and clone file for colouring", choices = c("test-csv" ,"own_csv")),
                                         fileInput('file_FACS.csv1', 'FACS+clone file',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                         
                                         fluidRow(
                                           
                                           column(6,selectInput("V.gene.1",label = h5("V Gene 1"),"")),
                                           column(6,selectInput("CDR3.1",label = h5("CDR3 1"),"")),
                                           column(6,selectInput('V.gene.2', label = h5("V Gene 2"), "")),
                                           column(6,selectInput("CDR3.2",label = h5("CDR3 2"),"")),
                                           
                                         ),
                                         fluidRow(
                                           column(6,selectInput("group.col.dot",label = h5("Group"),"")),
                                           column(6,numericInput("numeric.cloneCount","Filter based on number of times a clone was observed: select 0 for all",value=1))
                                         ),
                                         verbatimTextOutput("NAMES.df"),
                                         textInput("name.colour","Prefix of file name","ID.780_"),
                                         downloadButton('downloadTABLE_cleaning','Download table')
                        ),
# conditional panel 3 -----
                        conditionalPanel(condition="input.tabselected==3",
                                         selectInput("dataset_index.2", "Choose a dataset for complex plot:", choices = c("test-csv" ,"own_csv_file")),
                                         fileInput('file_FACS.csv2', 'File for dot plot',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                         selectInput("font_type2","Type of font",choices = font,selected = "Times"),
                                         fluidRow(
                                           column(6,selectInput("x.axis2",label = h5("Select x-axis"),"")),
                                           column(6,selectInput("y.axis2",label = h5("Select y-axis"),"")),
                                           
                                         ),
                                         fluidRow(
                                           column(4, selectInput("grid.lines.dot", label = h5("Add gridlines?"), choices = c("no","yes"))),
                                           column(4,selectInput("group_complex_dot",label = h5("Colour by:"),"")),
                                           column(4,selectInput( "FACS.index_colour.choise",label = h5("Colour"),choices = c("default","random","grey"), selected = "random")),
                                           
                                         ),
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
                        
           ),
           mainPanel(tabsetPanel(id = "tabselected",
                                 
# merging FACS file with clone file -----
                                 tabPanel("Merging paired TCR with Index data",value = 1,
                                          fluidRow(column(4, selectInput("group_FACS","Group of data","")),
                                                   column(4, selectInput("indiv_FACS","Individual of data","780")),
                                                   column(2, checkboxInput("multiple_plates","Multiple plates",value = F)),
                                                   column(2, numericInput("Plate_FACS","Plate #","1")),
                                          ),
                                          
                                          
                                          div(DT::dataTableOutput("FACS.CSV")),
                                          # div(DT::dataTableOutput("merged.clone")),
                                          div(DT::dataTableOutput("merged.index.clone")),
                                          
                                          
                                 ),
                                 
# UI complex dotplot add columns if needed -----
                                 tabPanel("Data cleaning steps",value = 2,
                                          tabsetPanel(
                                            tabPanel("Conversion to bioexponetional",
                                                     selectInput("string.data","Recommended selecting for ab TCR data: Indiv, group,TRBV,CDR3b.Sequence, TRBJ, TRAV, CDR3a.Sequence, TRAJ, AJ, BJ and AJBJ. Do not select flurochrome columns, or cloneCount","",multiple = T, width = "1200px"),
                                                     div(DT::dataTableOutput("table.index.1")),
                                                     
                                                     ),
                                            tabPanel("UMAP reduction",
                                                     selectInput("string.data_UMAP","Select flurochrome columns for dimension reduction (UMAP)","",multiple = T, width = "1200px"),
                                                     fluidRow(
                                                       column(4,numericInput("lower_range","Lower number of clusters",value=2)),
                                                       column(4,numericInput("upper_range","Upper number of clusters",value=10)),
                                                     ),
                                                     div(DT::dataTableOutput("table.index.UMAP")),
                                                     h4("Average of each flurochrome per cluster (log10 transformed)"),
                                                     div(DT::dataTableOutput("table.index.UMAP2")),
                                                     )  
                                                      ),
                                 ),
# UI complex dotplot -----
                                 tabPanel("TCR with Index data plot",value = 3,
                                          fluidRow(
                                            column(4,selectInput( "plot_type_umap","type of plot",choices = c("overlaid dot plot","Normal","Ridge plot"))),
                                            
                                            conditionalPanel(
                                              condition = "input.plot_type_umap == 'overlaid dot plot'",
                                              column(4,selectInput("density_dotplot",label = h5("Add histogram"), choices = c("no","yes"))),
                                              
                                            ),
                                            
                                              
                                              column(4,selectInput( "Add_ellipse","Add ellipse to clusters",choices = c("no","yes"))),
                                            
                                          ),
                                          # div(DT::dataTableOutput("Add_size")),
                                          
                                          
                                          fluidRow(column(3,
                                                          wellPanel(id = "tPanel222",style = "overflow-y:scroll; max-height: 250px",
                                                                    h4("Colour"),
                                                                    uiOutput('myPanel.FACS.index'),
                                                          )),
                                                   
                                                   column(3,
                                                          wellPanel(id = "tPanel222",style = "overflow-y:scroll; max-height: 250px",
                                                                    h4("Shape"),
                                                                    uiOutput('myPanel.FACS.index.shape')
                                                                    
                                                          )),
                                                   
                                                   column(3,
                                                          wellPanel(id = "tPanel222",style = "overflow-y:scroll; max-height: 250px",
                                                                    h4("Size"),
                                                                    uiOutput('myPanel.FACS.index.size')
                                                                    
                                                          )),
                                                   column(3, tags$img(src = "shape.png", height = "250px"))
                                                   
                                                   
                                          ),
                                          conditionalPanel(
                                            condition = "input.plot_type_umap == 'overlaid dot plot' || input.plot_type_umap == 'Normal'",
                                            fluidRow(column(12, plotOutput("dot_plot.complex2",height = "600px"))),
                                            
                                          ),
                                          
                                          conditionalPanel(
                                            condition = "input.plot_type_umap == 'Ridge plot'", 
                                            fluidRow(column(12, plotOutput("dot_plot_ridge",height = "400px")),
                                                     column(12, div(DT::dataTableOutput("dot_plot_ridge_tab")))),
                                            
                                          ),

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
                                 ),


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
  vals33 <- reactiveValues(geom_comp=NULL)
  vals44 <- reactiveValues(plot.ggseq.2=NULL)
  options(shiny.sanitize.errors = F)
  output$sessionInfo <- renderPrint({
    print(sessionInfo())
  })
  
  # video outputs -----
  output$video <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/IGNzlz3hQAc", width = 1000, height = 666.6666)
  })
  
  output$video2 <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/fjH06AbMjYE",  width = 1000, height = 666.6666)
  })
  
  output$video3 <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/d5IOs72kl5A", width = 1000, height = 666.6666)
  })
  
  output$video4 <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/wjCGSgmKrZA",  width = 1000, height = 666.6666)
  })

  output$video5 <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/1qVL4m8Bf4o",  width = 1000, height = 666.6666)
  })
  
  output$video6 <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/H_vzlJ3mOrQ",  width = 1000, height = 666.6666)
  })
  
  output$video7 <- renderUI({
    tags$iframe(src = "https://www.youtube.com/embed/g5xJHdMj4uA",  width = 1000, height = 666.6666)
  })
  
  
  # read in multiple .ab1 fles -----
  getData_ab1 <- eventReactive(input$go_ab1, {
    inFile.ab1 <- input$file1_ab1.file
    if (is.null(inFile.ab1)){
      return(NULL)
    }
    
    else {
      
      numfiles = nrow(inFile.ab1) 
      df_total = data.frame()
      
      for (i in 1:numfiles) { 
        tryCatch({
          hetsangerseq<- readsangerseq(input$file1_ab1.file[[i, 'datapath']]) 
          hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
          Primaryseq <- primarySeq(hetcalls, string = TRUE)
          secondary_seq <- secondarySeq(hetcalls, string = TRUE)
          
          pa <- pairwiseAlignment(primarySeq(hetcalls), secondarySeq(hetcalls))
          pa_score <- pa@score
          len_pa <- length(hetcalls@secondarySeq)[1]
          pa_score_base <- pa_score/len_pa
          
          
          if (input$indiv_miss_ab1 == "Yes" && input$group_miss_ab1 == "No") {
            name_temp <- paste(input$indiv.miss.name_ab1,".",inFile.ab1[i,1],sep = "")
          }
          
          else if (input$indiv_miss_ab1 == "No" && input$group_miss_ab1 == "Yes") {
            name_temp <- paste(input$group.miss.name_ab1,inFile.ab1[i,1],sep = "")
          }
          
          else if (input$indiv_miss_ab1 == "Yes" && input$group_miss_ab1 == "Yes") {
            name_temp <- paste(input$indiv.miss.name_ab1,".", input$group.miss.name_ab1,inFile.ab1[i,1],sep = "")
          }
          else {
            name_temp <- paste(inFile.ab1[i,1],sep = "") 
          }
          
          name_temp <- gsub(".ab1","", name_temp)
          
          df <- cbind(name_temp,pa_score,len_pa,pa_score_base)
          df_total <- rbind(df_total,df)
          
        }, error=function(e){}) 
      }
      df_total 
      
      }

  })
  
  output$contents_ab1 <- renderTable( 
    getData_ab1() 
  )
  
  output$downloadData_ab1_QC_file <- downloadHandler(
    filename = function() { 
      paste(input$ab1.name, " .ab1 QC ", Sys.Date(), ".csv", sep="")
    },
    content = function(file) { 
      write.csv(getData_ab1(), file, row.names=F)   
    })
  
  # .ab1 files for checking heitogenity ----
  input.data_ab1 <- reactive({switch(input$dataset_.ab1,".ab1-test-data" = test.data_ab.ab1(), ".ab1-own_data" = own.data.ab1_2())})
  test.data_ab.ab1 <- reactive({
    hetsangerseq <- readsangerseq("test-data/QC/SJS.TEN/E10630/Micromon/IFNg/IFNA-A10_C07.ab1") 
  })
  own.data.ab1_2 <- reactive({
    inFile_.ab1 <- input$file_.ab1
    if (is.null(inFile_.ab1)) return(NULL)
    
    else {
      hetsangerseq <- readsangerseq(
        inFile_.ab1$datapath
        
      )}
    
  })
  
  output$hetsangerseq <- renderPrint({
    
    hetsangerseq <- input.data_ab1()
    hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
    
    print(hetcalls)
  })
  
  chromatogram.seq1 <- reactive({
    
    hetsangerseq <- input.data_ab1()
    hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
    
    chromatogram(hetcalls, width = input$Number.seq.line, height = 4, cex.mtext = 1, cex.base = 10, showcalls = "both", trim5 =input$trim5.seq, trim3 = input$trim3.seq)
    
  })
  
  output$chromatogram.seq <- renderPlot({ 
    
    chromatogram.seq1()
    
  })
  
  
  output$downloadPlot_chromatogram.seq <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_chromatogram_",Sys.Date(), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_chromatogram.seq,height=input$height_chromatogram.seq, onefile = FALSE) # open the pdf device
      hetsangerseq <- input.data_ab1()
      
      validate(
        need(nrow(hetsangerseq@traceMatrix)>0,
             error_message_val1)
      )
      
      hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
      
      chromatogram(hetcalls, width = input$Number.seq.line, height = 4, cex.mtext = 1, cex.base = 3, showcalls = "both", trim5 =input$trim5.seq, trim3 = input$trim3.seq)
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_chromatogram.seq <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_chromatogram.seq_",Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_chromatogram.seq, height = input$height_png_chromatogram.seq, res = input$resolution_PNG_chromatogram.seq)
      
      
      hetsangerseq <- input.data_ab1()
      hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
      chromatogram(hetcalls, width = input$Number.seq.line, height = 4, cex.mtext = 1, cex.base = 3, showcalls = "both", trim5 =input$trim5.seq, trim3 = input$trim3.seq)
      
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  output$alignment <- renderPrint({
    
    hetsangerseq <- input.data_ab1()
    hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
    hetcalls
    chromatogram(hetcalls, width = 100, height = 2, showcalls = "both", trim5 =20, trim3 = 20)
    Primaryseq <- primarySeq(hetcalls, string = TRUE)
    secondary_seq <- secondarySeq(hetcalls, string = TRUE)
    
    pa <- pairwiseAlignment(primarySeq(hetcalls), secondarySeq(hetcalls))
    print(writePairwiseAlignments(pa))
  })
  
  
  
  # create and merge .seq to fasta files -----
  getData <- reactive({
    inFile.seq <- input$file1_seq.file
    if (is.null(inFile.seq)){
      return(NULL)
    }
    
    else {
      
      numfiles = nrow(inFile.seq) 
      
      
      df_total = data.frame()
      
      
      for (i in input$lower.seq:input$upper.seq) { 
        tryCatch({
          temp<- read.table(input$file1_seq.file[[i, 'datapath']],header = F) 
          temp
          
          
          if (input$indiv_miss == "Yes" && input$group_miss == "No") {
            
            name_temp <- paste("> ", input$indiv.miss.name,".",inFile.seq[i,1],"#",i,sep = "")
            
          }
          
          
          else if (input$indiv_miss == "No" && input$group_miss == "Yes") {
            
            name_temp <- paste("> ", input$group.miss.name,inFile.seq[i,1],"#",i,sep = "")
            
          }
          
          else if (input$indiv_miss == "Yes" && input$group_miss == "Yes") {
            
            name_temp <- paste("> ", input$indiv.miss.name,".", input$group.miss.name,inFile.seq[i,1],"#",i,sep = "")
            
          }
          
          
          
          else {
            name_temp <- paste("> ", inFile.seq[i,1],"#",i,sep = "") 
          }
          
          
          df <- rbind(name_temp,temp)
          df_total <- rbind(df_total,df)
          
        }, error=function(e){}) 
      }
      df_total }
    
  })
  
  output$contents <- renderTable( 
    getData() 
  )
  output$downloadData_fasta.files <- downloadHandler(
    filename = function() { 
      paste(input$seq.name, Sys.Date(), ".fasta", sep="")
    },
    content = function(file) { 
      write.table(getData(), file, row.names=F, col.names = F,quote = F)   
    })
  
  # IMGT only  -----
  input.data_IMGT.xls3 <- reactive({switch(input$dataset_IMGT3,"ab-test-data1" = test.data_ab.xls3(), "own_data" = own.data.IMGT3())})
  test.data_ab.xls3 <- reactive({
    dataframe = read_excel("test-data/QC/Vquest_data/CD8_E10630_A.xls") 
  })
  own.data.IMGT3 <- reactive({
    inFile_IMGT3 <- input$file_IMGT3
    if (is.null(inFile_IMGT3)) return(NULL)
    
    else {
      dataframe <- read_excel(
        inFile_IMGT3$datapath
        
      )}
    
  })
  input.data_IMGT.xls4 <- reactive({switch(input$dataset_IMGT3,"ab-test-data1" = test.data_ab.xls4(), "own_data" = own.data.IMGT4())})
  test.data_ab.xls4 <- reactive({
    dataframe = read_xls("test-data/QC/Vquest_data/CD8_E10630_A.xls",sheet = 2) 
  })
  
  own.data.IMGT4 <- reactive({
    inFile_IMGT4 <- input$file_IMGT3
    if (is.null(inFile_IMGT4)) return(NULL)
    
    else {
      dataframe <- read_excel(
        inFile_IMGT4$datapath, sheet = 2
        
      )}
    
  })
  IMGT2 <- reactive({
    df1 <- input.data_IMGT.xls3();
    
    validate(
      need(nrow(df1)>0,
           "Upload file")
    )
    df1
    if (input$sheet == "Summary+JUNCTION") {
      df2 <- input.data_IMGT.xls4();

      df3 <- df1[names(df1) %in% c("Sequence number","Sequence ID","V-DOMAIN Functionality", "V-GENE and allele","V-REGION identity %","J-GENE and allele","J-REGION identity %","D-GENE and allele","JUNCTION frame","JUNCTION (with frameshift)","CDR3-IMGT (with frameshift)","Sequence")]
      
      df4 <- df2[names(df2) %in% c("Sequence number","Sequence ID","JUNCTION","JUNCTION (AA)","JUNCTION (with frameshift)","JUNCTION (AA) (with frameshift)","CDR3-IMGT","CDR3-IMGT (AA)","V-REGION")]
      df_chain1 <- merge(df3,df4,by=c("Sequence number","Sequence ID"))
      df_chain1
      df_chain1 <- as.data.frame(df_chain1)
      df_chain1$`J-GENE and allele` <- gsub('Homsap ','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('Homsap ','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('Homsap ','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('Musmus ','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('Musmus ','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('Musmus ','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('see comment','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('see comment','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('see comment','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('[(]','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('[)]','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub(' F','',df_chain1$`D-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[(]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[)]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' F,','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' F','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[[]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[]]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('F','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' or ',', ',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' ','',df_chain1$`V-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('TR','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' or ',', ',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('TR','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('TR','',df_chain1$`D-GENE and allele`)
      df_chain1$`V-DOMAIN Functionality` <- gsub(' [(]see comment','',df_chain1$`V-DOMAIN Functionality`)
      df_chain1$`V-DOMAIN Functionality` <- gsub('[)]','',df_chain1$`V-DOMAIN Functionality`)

      df_chain1$JUNCTION <- toupper(df_chain1$JUNCTION)
      df_chain1
      name_temp2  <- as.data.frame(do.call(rbind, strsplit(as.character(df_chain1$`Sequence ID`), ".seq")))
      
      df_chain1$name_temp <- name_temp2$V1
      df_chain1
    }

    else  {
      df3 <- df1[names(df1) %in% c("Sequence number","Sequence ID","V-DOMAIN Functionality", "V-GENE and allele","V-REGION identity %","J-GENE and allele","J-REGION identity %","D-GENE and allele","JUNCTION frame","Sequence", "AA JUNCTION" )]

      df_chain1 <- df3
      df_chain1 <- as.data.frame(df_chain1)
      df_chain1$`J-GENE and allele` <- gsub('Homsap ','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('Homsap ','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('Homsap ','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('Musmus ','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('Musmus ','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('Musmus ','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('see comment','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('see comment','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('see comment','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('[(]','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('[)]','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub(' F','',df_chain1$`D-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[(]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[)]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' F,','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' F','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[[]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[]]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('F','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' or ',', ',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' ','',df_chain1$`V-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' or ',', ',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('TR','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('TR','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('TR','',df_chain1$`D-GENE and allele`)
      df_chain1$`V-DOMAIN Functionality` <- gsub(' [(]see comment','',df_chain1$`V-DOMAIN Functionality`)
      df_chain1$`V-DOMAIN Functionality` <- gsub('[)]','',df_chain1$`V-DOMAIN Functionality`)

      df_chain1$JUNCTION <- toupper(df_chain1$JUNCTION)
      name_temp2  <- as.data.frame(do.call(rbind, strsplit(as.character(df_chain1$`Sequence ID`), ".seq")))
      df_chain1$name_temp <- name_temp2$V1

      df_chain1

    }
    
    
  })
  
  IMGT2_original <- reactive({
    df1 <- input.data_IMGT.xls3();
    
    validate(
      need(nrow(df1)>0,
           "Upload file")
    )
    df1
    if (input$sheet == "Summary+JUNCTION") {
      df2 <- input.data_IMGT.xls4();
      
      df3 <- df1[names(df1) %in% c("Sequence number","Sequence ID","V-DOMAIN Functionality", "V-GENE and allele","V-REGION identity %","J-GENE and allele","J-REGION identity %","D-GENE and allele","JUNCTION frame","JUNCTION (with frameshift)","CDR3-IMGT (with frameshift)","Sequence")]
      
      df4 <- df2[names(df2) %in% c("Sequence number","Sequence ID","JUNCTION","JUNCTION (AA)","JUNCTION (with frameshift)","JUNCTION (AA) (with frameshift)","CDR3-IMGT","CDR3-IMGT (AA)","V-REGION")]
      df_chain1 <- merge(df3,df4,by=c("Sequence number","Sequence ID"))
      df_chain1
      df_chain1 <- as.data.frame(df_chain1)
      df_chain1$`J-GENE and allele` <- gsub('Homsap ','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('Homsap ','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('Homsap ','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('Musmus ','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('Musmus ','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('Musmus ','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('see comment','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('see comment','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('see comment','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('[(]','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('[)]','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub(' F','',df_chain1$`D-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[(]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[)]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' F,','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' F','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[[]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[]]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('F','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' or ',', ',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' ','',df_chain1$`V-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('TR','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' or ',', ',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('TR','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('TR','',df_chain1$`D-GENE and allele`)
      df_chain1$`V-DOMAIN Functionality` <- gsub(' [(]see comment','',df_chain1$`V-DOMAIN Functionality`)
      df_chain1$`V-DOMAIN Functionality` <- gsub('[)]','',df_chain1$`V-DOMAIN Functionality`)
      
      df_chain1$JUNCTION <- toupper(df_chain1$JUNCTION)
      df_chain1
      name_temp2  <- as.data.frame(do.call(rbind, strsplit(as.character(df_chain1$`Sequence ID`), ".seq")))
      
      df_chain1$name_temp <- name_temp2$V1
      df_chain1
    }
    
    else  {
      df3 <- df1[names(df1) %in% c("Sequence number","Sequence ID","V-DOMAIN Functionality", "V-GENE and allele","V-REGION identity %","J-GENE and allele","J-REGION identity %","D-GENE and allele","JUNCTION frame","Sequence", "AA JUNCTION" )]
      
      df_chain1 <- df3
      df_chain1 <- as.data.frame(df_chain1)
      df_chain1$`J-GENE and allele` <- gsub('Homsap ','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('Homsap ','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('Homsap ','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('Musmus ','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('Musmus ','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('Musmus ','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('see comment','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('see comment','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('see comment','',df_chain1$`D-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('[(]','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('[)]','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub(' F','',df_chain1$`D-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[(]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[)]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' F,','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' F','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[[]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('[]]','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('F','',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' or ',', ',df_chain1$`V-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub(' ','',df_chain1$`V-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub(' or ',', ',df_chain1$`J-GENE and allele`)
      df_chain1$`J-GENE and allele` <- gsub('TR','',df_chain1$`J-GENE and allele`)
      df_chain1$`V-GENE and allele` <- gsub('TR','',df_chain1$`V-GENE and allele`)
      df_chain1$`D-GENE and allele` <- gsub('TR','',df_chain1$`D-GENE and allele`)
      df_chain1$`V-DOMAIN Functionality` <- gsub(' [(]see comment','',df_chain1$`V-DOMAIN Functionality`)
      df_chain1$`V-DOMAIN Functionality` <- gsub('[)]','',df_chain1$`V-DOMAIN Functionality`)
      
      df_chain1$JUNCTION <- toupper(df_chain1$JUNCTION)
      # name_temp2  <- as.data.frame(do.call(rbind, strsplit(as.character(df_chain1$`Sequence ID`), ".seq")))
      # df_chain1$name_temp <- name_temp2$V1
      
      df_chain1
      
    }
    df_chain1$V.sequence.quality.check <- ifelse(df_chain1$`V-DOMAIN Functionality`== "No rearrangement found","No arrangement",
                                                 ifelse(df_chain1$`V-DOMAIN Functionality`=="unproductive", "Unproductive issue",
                                                        ifelse(df_chain1$`V-DOMAIN Functionality`=="rearranged sequence (but no junction found","Junction not found",
                                                               ifelse(df_chain1$`V-DOMAIN Functionality`=="No results", "No alignment",
                                                                      ifelse(as.numeric(df_chain1$`V-REGION identity %`<=90),"V Identity issue",
                                                                             ifelse(as.numeric(df_chain1$`J-REGION identity %`)<=80,"J Identity issue",
                                                                                    "No issue flagged by IMGT"))))))
    
    
    df_chain1$chromatogram_check <- NA
    df_chain1$clone_quality <- ifelse(df_chain1$V.sequence.quality.check=="No issue flagged by IMGT","pass",NA)
    df_chain1$comments <- NA
    df_chain1
    
  })
  
  output$IMGT2_original_out <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
    IMGT2_original()
  })
  
  output$downloadTABLE_IMGTonly <- downloadHandler(
    filename = function(){
      paste(input$IMGT_name_df," IMGT_manual ",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- IMGT2_original()
      write.csv(df,file, row.names = FALSE)
    } )
  
  output$IMGT2_out <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
    IMGT2()
  })
  
  input.data_IMGT.ab1 <-  reactive({switch(input$dataset_IMGT3,"ab-test-data1" = test.data_ab1(), "own_data" = own.data.ab1())})
  test.data_ab1 <- reactive({
    dataframe = read.csv("test-data/QC/QC .ab1 files/E0630 .ab1 QC 2022-12-30.csv") 
  })
  own.data.ab1 <- reactive({
    inFile_ab1 <- input$file_ab1_QC
    if (is.null(inFile_ab1)) return(NULL)
    
    else {
      dataframe <- read.csv(
        inFile_ab1$datapath
        
      )}
    
  })
  
  output$IMGT2_out_ab1 <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
    df1 <- as.data.frame(input.data_IMGT.ab1())
    
    validate(
      need(nrow(df1)>0,
           "Upload file")
    )
    
     df1
  })

  IMGT_ab1_QC_merge <- reactive({
    df_chain1 <- as.data.frame(IMGT2())
    ab1_File <- as.data.frame(input.data_IMGT.ab1())
    validate(
      need(nrow(df_chain1)>0 & nrow(ab1_File)>0,
           "Upload file")
    )
    
    df_chain1 <- merge(df_chain1,ab1_File,by = "name_temp")
    df_chain1 <- df_chain1[ , !(names(df_chain1) %in%  c("name_temp"))]

    df_chain1$V.sequence.quality.check <- ifelse(df_chain1$`V-DOMAIN Functionality`== "No rearrangement found","No arrangement",
                                                 ifelse(df_chain1$`V-DOMAIN Functionality`=="unproductive", "Unproductive issue",
                                                        ifelse(df_chain1$`V-DOMAIN Functionality`=="rearranged sequence (but no junction found","Junction not found",
                                                                ifelse(df_chain1$`V-DOMAIN Functionality`=="No results", "No alignment",
                                                                       ifelse(as.numeric(df_chain1$`V-REGION identity %`<=90),"V Identity issue",
                                                                              ifelse(as.numeric(df_chain1$`J-REGION identity %`)<=80,"J Identity issue",
                                                                                     "No issue flagged by IMGT"))))))


    df_chain1$chromatogram_check <- ifelse(df_chain1$pa_score_base < 0.2,"Poor",
                                           ifelse(df_chain1$pa_score_base > 1,"Very high",
                                                  ifelse(df_chain1$pa_score_base > 0.9, "High",
                                                         ifelse(df_chain1$pa_score_base >0.7,"Moderate","Low"))))

    df_chain1$clone_quality <- ifelse(df_chain1$V.sequence.quality.check=="No issue flagged by IMGT" & df_chain1$chromatogram_check=="Very high",'pass',
                                      ifelse(df_chain1$V.sequence.quality.check=="No issue flagged by IMGT" & df_chain1$chromatogram_check=="High",'pass',
                                      ifelse(df_chain1$V.sequence.quality.check=="No issue flagged by IMGT" & df_chain1$chromatogram_check=="Moderate",'pass',
                                      ifelse(df_chain1$V.sequence.quality.check=="No alignment" |
                                                      df_chain1$chromatogram_check=="Poor" |
                                                      df_chain1$V.sequence.quality.check=="No arrangement" |
                                                      df_chain1$chromatogram_check=="Low" ,'fail',
                                                    ifelse(df_chain1$V.sequence.quality.check=="unproductive" &
                                                             df_chain1$chromatogram_check=="Very high","fail",

                                                    ifelse( df_chain1$chromatogram_check=="Very high" & df_chain1$V.sequence.quality.check=="Junction not found", "fail",

                                                    ifelse(df_chain1$chromatogram_check=="Very high" & df_chain1$V.sequence.quality.check=="J Identity issue","pass",
                                                    ifelse(df_chain1$chromatogram_check=="Very high" & df_chain1$V.sequence.quality.check=="V Identity issue","pass",

                                                     ifelse(df_chain1$chromatogram_check=="High" & df_chain1$V.sequence.quality.check=="J Identity issue","pass",
                                                     ifelse(df_chain1$chromatogram_check=="High" & df_chain1$V.sequence.quality.check=="V Identity issue","pass",
                                                     
                                                           	
                                                    'fail'))))))))))
    #

    df_chain1$comment <- ifelse(df_chain1$V.sequence.quality.check=="No issue flagged by IMGT" & df_chain1$chromatogram_check=="Very high",NA,
                                      ifelse(df_chain1$V.sequence.quality.check=="No issue flagged by IMGT" & df_chain1$chromatogram_check=="High",NA,
                                             ifelse(df_chain1$chromatogram_check=="Moderate" & df_chain1$V.sequence.quality.check=="No issue flagged by IMGT","possible need to check moderate quality",
                                             ifelse(df_chain1$V.sequence.quality.check=="No alignment" |
                                                      df_chain1$chromatogram_check=="Poor" |
                                                      df_chain1$V.sequence.quality.check=="No arrangement" |
                                                      df_chain1$chromatogram_check=="Low" ,'Either poor sequence quality or no sequence called',
                                                    ifelse(df_chain1$V.sequence.quality.check=="Unproductive issue" &
                                                             df_chain1$chromatogram_check=="Very high","Unproductive high quality sequence, possibly resolvable",

                                                           ifelse( df_chain1$chromatogram_check=="Very high" & df_chain1$V.sequence.quality.check=="Junction not found", "No Junction found",

                                                      ifelse(df_chain1$chromatogram_check=="Very high" & df_chain1$V.sequence.quality.check=="J Identity issue","Possible J identity issue",
                                                      ifelse(df_chain1$chromatogram_check=="Very high" & df_chain1$V.sequence.quality.check=="V Identity issue","Possible V identity issue",

                                                      ifelse(df_chain1$chromatogram_check=="High" & df_chain1$V.sequence.quality.check=="J Identity issue","Possible J identity issue",
                                                      ifelse(df_chain1$chromatogram_check=="High" & df_chain1$V.sequence.quality.check=="V Identity issue","Possible V identity issue",
                                                                                     'Other possible issue'))))))))))

    df_chain1
    # ab1_File
  })
  
  output$IMGT2_out2 <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    IMGT_ab1_QC_merge()
  })
  
  output$downloadTABLE_IMGTonly_auto <- downloadHandler(
    filename = function(){
      paste(input$IMGT_name_df," IMGT_auto ",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- IMGT_ab1_QC_merge()
      write.csv(df,file, row.names = FALSE)
    } )
  
  
  # table of IMGT for pairing -----
  input.data.IMGT_afterQC <- reactive({switch(input$dataset_IMGT_afterQC,"ab-test-data1" = test.data.ab.csv3(), "own_data1" = own.data.csv3())})
  test.data.ab.csv3 <- reactive({
    dataframe = read.csv("test-data/QC/QC.csv_files/SJS.TEN.three.samps.csv",header=T) 
  })
  own.data.csv3 <- reactive({
    inFile12 <- input$file_IMGT_afterQC
    if (is.null(inFile12)) return(NULL)
    
    else {
      dataframe <- read.csv(
        inFile12$datapath)}
    
  })

  Pass.Fail.NA <- reactive({
    df1 <- input.data.IMGT_afterQC();
    
    validate(
      need(nrow(df1)>0,
           "Upload file")
    )
    
    df1$clone_quality <- gsub("pass","pass",df1$clone_quality,ignore.case = T)
    df1$clone_quality <- gsub("fail","fail",df1$clone_quality,ignore.case = T)
    df1$cloneCount <- 1
    df2 <- df1[,c("cloneCount","clone_quality","V.sequence.quality.check","chromatogram_check")] 
    
    as.data.frame(ddply(df2,c("clone_quality","V.sequence.quality.check","chromatogram_check"),numcolwise(sum)))
    

    
  })
  
  output$Pass.Fail.NA_table<- DT::renderDataTable( {
    datatable(Pass.Fail.NA(), extensions = "Buttons", options = list(searching = TRUE,
                                                                     ordering = TRUE,
                                                                     buttons = c('copy','csv', 'excel'),
                                                                     dom = 'Bfrtip',
                                                                     pageLength=10, 
                                                                     lengthMenu=c(2,5,10,20,50,100), 
                                                                     scrollX = TRUE
    ))
  }, server = FALSE)
  
  chain_merge_IMGTonly <- reactive({
    df1 <- input.data.IMGT_afterQC();
    
    validate(
      need(nrow(df1)>0,
           "Upload file")
    )
    
    df1 <- as.data.frame(df1)
    df1$clone_quality <- gsub("pass","pass",df1$clone_quality,ignore.case = T)
    df1$clone_quality <- gsub("fail","fail",df1$clone_quality,ignore.case = T)
    df <- subset(df1,df1$clone_quality=="pass")
    df <- as.data.frame(df)
    df2 <- df[!names(df) %in% c("V.sequence.quality.check","clone_quality","comments","JUNCTION..with.frameshift.","CDR3.IMGT..with.frameshift.","JUNCTION..AA...with.frameshift.","Sequence.number","V.REGION.identity..","J.REGION.identity..","comment","pa_score","len_pa","pa_score_base","chromatogram_check")]
    
    df.Vgene <- as.data.frame(do.call(rbind, strsplit(as.character(df2$V.GENE.and.allele), ",")))
    df2$V.GENE <- df.Vgene$V1
    y = dim(df2)[2]
    y
    df2$V.GENE <- gsub(" ","",df2$V.GENE)
    df2$cloneCount <- 1
    if (input$IMGT_chain2 =="ab" & input$sheet2 == "Summary+JUNCTION") {
      
      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), ".-")))
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      
      df2$ID <- df_name2$V1
      head(df2)
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$well <- df_name3$V2
      
      chain1 <- df2[grep("A",df2$V.GENE),]
      chain2 <- df2[grep("B",df2$V.GENE),]
      # paste chain into name
      chain1$ID <- gsub("A-","-",chain1$ID)
      names(chain1)[1:y] <- paste(names(chain1)[1:y],"A",sep="_")
      head(chain1)
      chain2$ID <- gsub("B-","-",chain2$ID)
      names(chain2)[1:y] <- paste(names(chain2)[1:y],"B",sep="_")
      z = y+1
      x <- names(chain2)[z:dim(chain2)[2]]
      
      merged_chain <- merge(chain1,chain2,by =x)
      head(merged_chain)
      merged_chain2 <- merged_chain[ , -which(names(merged_chain) %in% c("ID","Sequence.ID_A","Sequence.ID_B","V.DOMAIN.Functionality_A","V.DOMAIN.Functionality_B","D.GENE.and.allele_A","JUNCTION.frame_A","JUNCTION.frame_B","comment_A","comment_B"))]
      names(merged_chain2)
      dat <- merged_chain2
      dat$AV <- paste(dat$V.GENE_A)
      dat$AJ <- paste(dat$J.GENE.and.allele_A,sep="")
      dat$AJ <- gsub("[*]0.","",dat$AJ)
      dat$AJ <- gsub(",, AJ..","",dat$AJ)
      dat$AJ <- gsub(",, AJ.","",dat$AJ)
      
      dat$AVJ <- paste(dat$AV,".",dat$AJ,sep="")
      dat$AV <- gsub("[*]0.","",dat$AV)
      dat$AVJ <- gsub("[*]0.","",dat$AVJ)
      dat$BV <- paste(dat$V.GENE_B)
      dat$BJ <- paste(dat$J.GENE.and.allele_B)
      dat$BD <- paste(dat$D.GENE.and.allele_B)
      dat$BVJ <- paste(dat$BV,".",dat$BJ,sep="")
      dat$BVDJ <- paste(dat$BV,".",dat$BD,".",dat$BJ,sep="")
      
      dat$BV <- gsub("[*]0.","",dat$BV)
      dat$BJ <- gsub("[*]0.","",dat$BJ)
      dat$BD <- gsub("[*]0.","",dat$BD)
      dat$BVJ <- gsub("[*]0.","",dat$BVJ)
      dat$BVDJ <- gsub("[*]0.","",dat$BVDJ)
      dat$BVDJ <- gsub(".NA.",".",dat$BVDJ)
      
      dat$AJ <- gsub("TR","",dat$AJ)
      dat$AVJ <- gsub("TR","",dat$AVJ)
      dat$AVJ <- gsub("AJ","J",dat$AVJ)
      dat$AVJ.BVJ <- paste(dat$AVJ,"_",dat$BVJ,sep="")
      dat$AVJ.BVDJ <- paste(dat$AVJ,"_",dat$BVDJ,sep="")
      dat$AVJ_aCDR3 <- paste(dat$AVJ,dat$JUNCTION..AA._A,sep="_")
      dat$BVJ_bCDR3 <- paste(dat$BVJ,dat$JUNCTION..AA._B,sep="_")
      dat$BVDJ_bCDR3 <- paste(dat$BVDJ,dat$JUNCTION..AA._B,sep="_")
      
      dat$AVJ_aCDR3_BVJ_bCDR3 <- paste(dat$AVJ_aCDR3,dat$BVJ_bCDR3,sep=" & ")
      dat$AVJ_aCDR3_BVDJ_bCDR3 <- paste(dat$AVJ_aCDR3,dat$BVDJ_bCDR3,sep=" & ")
      dat$BD <- gsub("NA","-",dat$BD)
      head(dat)
      
      dat
      
    }
    else if (input$IMGT_chain2 =="ab" & input$sheet2 == "Summary") {
      
      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), ".-")))
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      
      df2$ID <- df_name2$V1
      head(df2)
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$well <- df_name3$V2
      
      chain1 <- df2[grep("A",df2$V.GENE),]
      chain2 <- df2[grep("B",df2$V.GENE),]
      # paste chain into name
      chain1$ID <- gsub("A-","-",chain1$ID)
      names(chain1)[1:y] <- paste(names(chain1)[1:y],"A",sep="_")
      head(chain1)
      chain2$ID <- gsub("B-","-",chain2$ID)
      names(chain2)[1:y] <- paste(names(chain2)[1:y],"B",sep="_")
      z = y+1
      x <- names(chain2)[z:dim(chain2)[2]]
      
      merged_chain <- merge(chain1,chain2,by =x)
      head(merged_chain)
      merged_chain2 <- merged_chain[ , -which(names(merged_chain) %in% c("ID","Sequence.ID_A","Sequence.ID_B","V.DOMAIN.Functionality_A","V.DOMAIN.Functionality_B","D.GENE.and.allele_A","JUNCTION.frame_A","JUNCTION.frame_B","JUNCTION_A","JUNCTION_B"))]
      names(merged_chain2)
      dat <- merged_chain2
      dat$AV <- paste(dat$V.GENE_A)
      dat$AJ <- paste(dat$J.GENE.and.allele_A,sep="")
      dat$AVJ <- paste(dat$AV,".",dat$AJ,sep="")
      dat$AV <- gsub("[*]0.","",dat$AV)
      dat$AJ <- gsub("[*]0.","",dat$AJ)
      dat$AVJ <- gsub("[*]0.","",dat$AVJ)
      
      dat$BV <- paste(dat$V.GENE_B)
      dat$BJ <- paste(dat$J.GENE.and.allele_B)
      dat$BD <- paste(dat$D.GENE.and.allele_B)
      dat$BVJ <- paste(dat$BV,".",dat$BJ,sep="")
      dat$BVDJ <- paste(dat$BV,".",dat$BD,".",dat$BJ,sep="")
      dat$BV <- gsub("[*]0.","",dat$BV)
      dat$BJ <- gsub("[*]0.","",dat$BJ)
      dat$BD <- gsub("[*]0.","",dat$BD)
      dat$BVJ <- gsub("[*]0.","",dat$BVJ)
      dat$BVDJ <- gsub("[*]0.","",dat$BVDJ)
      dat$BVDJ <- gsub(".NA.",".",dat$BVDJ)
      
      
      
      dat$AJ <- gsub("TR","",dat$AJ)
      dat$AVJ <- gsub("TR","",dat$AVJ)
      dat$AVJ <- gsub("AJ","J",dat$AVJ)
      dat$AVJ.BVJ <- paste(dat$AVJ,"_",dat$BVJ,sep="")
      dat$AVJ.BVDJ <- paste(dat$AVJ,"_",dat$BVDJ,sep="")
      dat$AVJ_aCDR3 <- paste(dat$AVJ,dat$AA.JUNCTION_A,sep="_")
      dat$BVJ_bCDR3 <- paste(dat$BVJ,dat$AA.JUNCTION_B,sep="_")
      dat$BVDJ_bCDR3 <- paste(dat$BVDJ,dat$AA.JUNCTION_B,sep="_")
      
      dat$AVJ_aCDR3_BVJ_bCDR3 <- paste(dat$AVJ_aCDR3,dat$BVJ_bCDR3,sep=" & ")
      dat$AVJ_aCDR3_BVDJ_bCDR3 <- paste(dat$AVJ_aCDR3,dat$BVDJ_bCDR3,sep=" & ")
      dat$BD <- gsub("NA","-",dat$BD)
      head(dat)
      
      dat
      
    }
    else if (input$IMGT_chain2 =="gd" & input$sheet2 == "Summary+JUNCTION") {
      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), ".-")))
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      
      df2$ID <- df_name2$V1
      head(df2)
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$well <- df_name3$V2
      
      chain1 <- df2[grep("G",df2$V.GENE.and.allele),]
      chain2 <- df2[grep("D",df2$V.GENE.and.allele),]
      chain1$ID <- gsub("G-","-",chain1$ID)
      names(chain1)[1:y] <- paste(names(chain1)[1:y],"G",sep="_")
      
      chain2$ID <- gsub("D-","-",chain2$ID)
      names(chain2)[1:y] <- paste(names(chain2)[1:y],"D",sep="_")
      head(chain2)
      
      z = y+1
      x <- names(chain2)[z:dim(chain2)[2]]
      merged_chain <- merge(chain1,chain2,by =x)
      head(merged_chain)
      merged_chain2 <- merged_chain[ , -which(names(merged_chain) %in% c("ID","Sequence.ID_G","Sequence.ID_D","V.DOMAIN.Functionality_G","V.DOMAIN.Functionality_D","D.GENE.and.allele_G","JUNCTION.frame_G","JUNCTION.frame_D"))]
      dat <- merged_chain2
      dat$GV <- paste(dat$V.GENE_G)
      dat$GJ <- paste(dat$J.GENE.and.allele_G,sep="")
      
      dat$GV <- gsub("[*]0.","",dat$GV)
      dat$GJ <- gsub("[*]0.","",dat$GJ)
      dat$GJ <- gsub("TR","",dat$GJ)
      dat$GJ <- gsub("GJ","J",dat$GJ)
      dat$GJ <- gsub(", or J.","",dat$GJ)
      dat$GVJ <- paste(dat$GV,".",dat$GJ,sep="")
      
      
      dat$DV <- paste(dat$V.GENE_D)
      dat$DJ <- paste(dat$J.GENE.and.allele_D)
      dat$DD <- paste(dat$D.GENE.and.allele_D)
      
      dat$DV <- gsub("[*]0.","",dat$DV)
      dat$DJ <- gsub("[*]0.","",dat$DJ)
      dat$DD <- gsub("[*]0.","",dat$DD)
      dat$DD <- gsub(" and ",".",dat$DD)
      
      dat$DVJ <- paste(dat$DV,".",dat$DJ,sep="")
      dat$DVDJ <- paste(dat$DV,".",dat$DD,".",dat$DJ,sep="")
      
      dat$DVJ <- gsub("[*]0.","",dat$DVJ)
      dat$DVDJ <- gsub("[*]0.","",dat$DVDJ)
      dat$DVDJ <- gsub(".NA.",".",dat$DVDJ)
      
      dat$GJ <- gsub("TR","",dat$GJ)
      dat$GVJ <- gsub("TR","",dat$GVJ)
      dat$GVJ.DVJ <- paste(dat$GVJ,"_",dat$DVJ,sep="")
      dat$GVJ.DVDJ <- paste(dat$GVJ,"_",dat$DVDJ,sep="")
      dat$GVJ_gCDR3 <- paste(dat$GVJ,dat$JUNCTION..AA._G,sep="_")
      dat$DVJ_dCDR3 <- paste(dat$DVJ,dat$JUNCTION..AA._D,sep="_")
      dat$DVDJ_dCDR3 <- paste(dat$DVDJ,dat$JUNCTION..AA._D,sep="_")
      
      dat$GVJ_gCDR3_DVJ_dCDR3 <- paste(dat$GVJ_gCDR3,dat$DVJ_dCDR3,sep=" & ")
      dat$GVJ_gCDR3_DVDJ_dCDR3 <- paste(dat$GVJ_gCDR3,dat$DVDJ_dCDR3,sep=" & ")
      dat$DD <- gsub("NA","-",dat$DD)
      
      dat
    }
    else  {
      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), ".-")))
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      
      df2$ID <- df_name2$V1
      head(df2)
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$well <- df_name3$V2
      
      chain1 <- df2[grep("G",df2$V.GENE.and.allele),]
      chain2 <- df2[grep("D",df2$V.GENE.and.allele),]
      chain1$ID <- gsub("G-","-",chain1$ID)
      names(chain1)[1:y] <- paste(names(chain1)[1:y],"G",sep="_")
      
      chain2$ID <- gsub("D-","-",chain2$ID)
      names(chain2)[1:y] <- paste(names(chain2)[1:y],"D",sep="_")
      head(chain2)
      
      z = y+1
      x <- names(chain2)[z:dim(chain2)[2]]
      merged_chain <- merge(chain1,chain2,by =x)
      head(merged_chain)
      merged_chain2 <- merged_chain[ , -which(names(merged_chain) %in% c("ID","Sequence.ID_G","Sequence.ID_D","V.DOMAIN.Functionality_G","V.DOMAIN.Functionality_D","D.GENE.and.allele_G","JUNCTION.frame_G","JUNCTION.frame_D","JUNCTION_G","JUNCTION_D"))]
      
      dat <- merged_chain2
      
      dat$GV <- paste(dat$V.GENE_G)
      dat$GJ <- paste(dat$J.GENE.and.allele_G,sep="")
      
      dat$GV <- gsub("[*]0.","",dat$GV)
      dat$GJ <- gsub("[*]0.","",dat$GJ)
      dat$GJ <- gsub("TR","",dat$GJ)
      dat$GJ <- gsub("GJ","J",dat$GJ)
      dat$GJ <- gsub(", or J.","",dat$GJ)
      dat$GVJ <- paste(dat$GV,".",dat$GJ,sep="")
      
      dat$DV <- paste(dat$V.GENE_D)
      dat$DJ <- paste(dat$J.GENE.and.allele_D)
      dat$DD <- paste(dat$D.GENE.and.allele_D)
      
      dat$DV <- gsub("[*]0.","",dat$DV)
      dat$DJ <- gsub("[*]0.","",dat$DJ)
      dat$DD <- gsub("[*]0.","",dat$DD)
      dat$DD <- gsub(" and ",".",dat$DD)
      
      dat$DVJ <- paste(dat$DV,".",dat$DJ,sep="")
      dat$DVDJ <- paste(dat$DV,".",dat$DD,".",dat$DJ,sep="")
      
      dat$DVJ <- gsub("[*]0.","",dat$DVJ)
      dat$DVDJ <- gsub("[*]0.","",dat$DVDJ)
      dat$DVDJ <- gsub(".NA.",".",dat$DVDJ)
      
      dat$DD <- gsub("NA","no DD",dat$DD)
      
      dat$GVJ <- gsub("TR","",dat$GVJ)
      dat$GVJ.DVJ <- paste(dat$GVJ,"_",dat$DVJ,sep="")
      dat$GVJ.DVDJ <- paste(dat$GVJ,"_",dat$DVDJ,sep="")
      dat$GVJ_aCDR3 <- paste(dat$GVJ,dat$AA.JUNCTION_G,sep="_")
      dat$DVJ_dCDR3 <- paste(dat$DVJ,dat$AA.JUNCTION_D,sep="_")
      dat$DVDJ_dCDR3 <- paste(dat$DVDJ,dat$AA.JUNCTION_D,sep="_")
      
      dat$GVJ_gCDR3_DVJ_dCDR3 <- paste(dat$GVJ_gCDR3,dat$DVJ_dCDR3,sep=" & ")
      dat$GVJ_gCDR3_DVDJ_dCDR3 <- paste(dat$GVJ_gCDR3,dat$DVDJ_dCDR3,sep=" & ")
      dat
    }
    
  })
  
  chain_single.chain_IMGT <- reactive({
    df1 <- input.data.IMGT_afterQC();
    
    validate(
      need(nrow(df1)>0,
           "Upload file")
    )
    
    df1 <- as.data.frame(df1)
    df1$clone_quality <- gsub("pass","pass",df1$clone_quality,ignore.case = T)
    df1$clone_quality <- gsub("fail","fail",df1$clone_quality,ignore.case = T)
    df <- subset(df1,df1$clone_quality=="pass")
    df <- as.data.frame(df)
    df2 <- df[!names(df) %in% c("V.sequence.quality.check","clone_quality","comments","JUNCTION..with.frameshift.","CDR3.IMGT..with.frameshift.","JUNCTION..AA...with.frameshift.","Sequence.number","V.REGION.identity..","J.REGION.identity..")]
    
    df.Vgene <- as.data.frame(do.call(rbind, strsplit(as.character(df2$V.GENE.and.allele), ",")))
    df2$V.GENE <- df.Vgene$V1
    y = dim(df2)[2]
    y
    df2$V.GENE <- gsub(" ","",df2$V.GENE)
    df2$cloneCount <- 1
    
    if (input$sheet2 == "Summary+JUNCTION") {
      
      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), ".-")))
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      
      df2$ID <- df_name2$V1
      head(df2)
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$well <- df_name3$V2
      
      df2 <- df2[!names(df2) %in% c("ID","Sequence.number","Sequence.ID","V.DOMAIN.Functionality","JUNCTION.frame","JUNCTION")]
      df2 <- df2[,c(9:13,1:8)]
      df2

      dat <- as.data.frame(df2)
      dat$TRV <- paste(dat$V.GENE)
      dat$TRJ <- paste(dat$J.GENE.and.allele,sep="")
      dat$TRD <- paste(dat$D.GENE.and.allele,sep="")
      
      dat$TRV <- gsub("[*]0.","",dat$TRV)
      dat$TRJ <- gsub("[*]0.","",dat$TRJ)
      dat$TRD <- gsub("[*]0.","",dat$TRD)
      
      dat$TRVJ <- paste(dat$TRV,dat$TRJ,sep=".")
      dat$TRVDJ <- paste(dat$TRV,dat$TRD,dat$TRJ,sep=".")
    

       dat$TRVDJ <- gsub(".NA.",".",dat$TRVDJ)


      dat$TRVJ_CDR3 <- paste(dat$TRVJ,dat$JUNCTION..AA.,sep="_")
      dat$TRVDJ_CDR3 <- paste(dat$TRVDJ,dat$JUNCTION..AA.,sep="_")

      dat$TRD <- gsub("NA","-",dat$TRD)
      head(dat)
      # 
      dat
      
    }
    else {
      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), ".-")))
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      
      df2$ID <- df_name2$V1
      head(df2)
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$well <- df_name3$V2
      
      df2 <- df2[!names(df2) %in% c("ID","Sequence.number","Sequence.ID","V.DOMAIN.Functionality","JUNCTION.frame","JUNCTION")]
      df2 <- df2[,c(9:13,1:8)]
      
      dat <- as.data.frame(df2)
      dat$TRV <- paste(dat$V.GENE)
      dat$TRJ <- paste(dat$J.GENE.and.allele,sep="")
      dat$TRD <- paste(dat$D.GENE.and.allele,sep="")
      
      dat$TRV <- gsub("[*]0.","",dat$TRV)
      dat$TRJ <- gsub("[*]0.","",dat$TRJ)
      dat$TRD <- gsub("[*]0.","",dat$TRD)
      
      dat$TRVJ <- paste(dat$TRV,dat$TRJ,sep=".")
      dat$TRVDJ <- paste(dat$TRV,dat$TRD,dat$TRJ,sep=".")
      dat$TRVDJ <- gsub(".NA.",".",dat$TRVDJ)
      dat$TRVJ_CDR3 <- paste(dat$TRVJ,dat$AA.JUNCTION,sep="_")
      dat$TRVDJ_CDR3 <- paste(dat$TRVDJ,dat$AA.JUNCTION,sep="_")
      
      dat$TRD <- gsub("NA","-",dat$TRD)
      head(dat)
      # 
      dat
      
    }

  })
  
  
  # TSV output
  TSV.file.chain <- reactive({
    dat <- chain_merge_IMGTonly()
    dat$id <- paste0(dat$Indiv,".",dat$group,"-",dat$well) 
    
    
    
    if (input$IMGT_chain2 == "ab") {
      dat <-  dat[,c("id","group","Indiv","Sequence_A","Sequence_B")]
      names(dat) <- c("id","epitope","subject","a_nucseq","b_nucseq")
      dat$a_quals <- ""
      dat$b_quals <- ""
      dat 
      
    }
    
    else {
      dat <-  dat[,c("id","group","Indiv","Sequence_G","Sequence_D")]
      names(dat) <- c("id","epitope","subject","g_nucseq","d_nucseq")
      dat$g_quals <- ""
      dat$d_quals <- ""
      dat
    }
    
  })
  
  output$chain_table_IMGT.tcrdist <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    TSV.file.chain()
  })
  output$downloadTABLE.TSV <- downloadHandler(
    filename = function(){
      paste(input$IMGT_name_df," TCRdist.tsv", sep = "")
    },
    content = function(file){
      df <- TSV.file.chain()
      df <- as.data.frame(df)
      
      write.table(df, file, quote=FALSE, sep='\t', row.names = F)
      
    } )
  
  # paired chain output
  output$chain_table_IMGT.QC1 <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    df1 <- input.data.IMGT_afterQC();
    df1 <- as.data.frame(df1)
    a <- subset(df1 ,is.na(df1$clone_quality)==TRUE)
    if (dim(a)[1]>0) {
      df <- as.data.frame("please complete QC analysis")
      names(df) <- " "
      df
    }
    else {
      df <- chain_merge_IMGTonly()
      df <- as.data.frame(df)
      df
      
    }
  })
  output$downloadTABLE.QC1 <- downloadHandler(
    filename = function(){
      paste(input$IMGT_name_df," paired_TCR_file ",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- chain_merge_IMGTonly()
      df <- as.data.frame(df)
      write.csv(df,file, row.names = FALSE)
    } )
  

  
  # single chain oytput 
  output$single.chain_table_IMGT.QC1 <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    df1 <- input.data.IMGT_afterQC();
    df1 <- as.data.frame(df1)
    a <- subset(df1 ,is.na(df1$clone_quality)==TRUE)
    if (dim(a)[1]>0) {
      df <- as.data.frame("please complete QC analysis")
      names(df) <- " "
      df
    }
    else {
      df <- chain_single.chain_IMGT()
      df <- as.data.frame(df)
      df
      
    }
  })
  output$downloadTABLE.QC1.single.chain <- downloadHandler(
    filename = function(){
      paste(input$IMGT_name_df," single.chain_TCR_file ",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- chain_single.chain_IMGT()
      df <- as.data.frame(df)
      write.csv(df,file, row.names = FALSE)
    } )
  
  
  
  # Immunoseq QC -----
  input.data.Immunoseq <- reactive({switch(input$dataset_TSV.Immunoseq,"Demo.gd.test-data" = test.data.ImmunoSeq(), "Immunoseq.own-data" = own.data.immmunoseq())})
    test.data.ImmunoSeq <- reactive({
    dataframe = read.table("test-data/QC/ImmunoSEQ/ES8_TSNLQEQIGW_3.tsv",sep="\t",header=T)
  })
    own.data.immmunoseq <- reactive({
    inFile_immunoseq <- input$file_TSV.Immunoseq
    if (is.null(inFile_immunoseq)) return(NULL)
    
    else {
      dataframe <- read.table(
        inFile_immunoseq$datapath,
        sep = input$sep.imm,
        quote = input$quote.imm,
        header = T)}
    
  })
    
    # for the observe event
    TSV.col.names <- reactive({
      x <- as.data.frame(input.data.Immunoseq())
      x2 <- x %>%
        select_if(~ !any(is.na(.)))
      
      x2
    })
    
 # count column ----
    observe({
      if (input$datasource == "ImmunoSEQ") {
        updateSelectInput(
          session,
          "countcolumn",
          choices=names(TSV.col.names()),
          selected = c("count..templates.reads."))
        
      }
      else if (input$datasource == "MiXCR") {
      
        updateSelectInput(
          session,
          "countcolumn",
          choices=names(TSV.col.names()),
          selected = c("cloneCount"))
        
      }
      
      else if (input$datasource == "10x_scSeq") {
        
        updateSelectInput(
          session,
          "countcolumn",
          choices=names(TSV.col.names()),
          selected = c("number_clonotypes"))
        
      }  
      
      
      else {
        updateSelectInput(
          session,
          "countcolumn",
          choices=names(TSV.col.names()),
          selected = c("count"))
      }
      
    }) 
    # J gene
    observe({
      
      if (input$datasource == "ImmunoSEQ") {
        updateSelectInput(
          session,
          "J.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("jFamilyName"))
        
      }
      else if (input$datasource == "MiXCR") {
        
        updateSelectInput(
          session,
          "J.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("bestJHit"))
        
      }
      
      else if (input$datasource == "10x_scSeq") {
        
        updateSelectInput(
          session,
          "J.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("j_gene_B"))
        
      }  
      
      
      else {
        updateSelectInput(
          session,
          "J.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("j_gene"))
      }

    })
    # V gene
    observe({
      if (input$datasource == "ImmunoSEQ") {
        updateSelectInput(
          session,
          "V.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("vFamilyName"))
        
      }
      else if (input$datasource == "MiXCR") {
        
        updateSelectInput(
          session,
          "V.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("bestVHit"))
        
      }
      
      else if (input$datasource == "10x_scSeq") {
        
        updateSelectInput(
          session,
          "V.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("v_gene_B"))
        
      }
      
      
      else {
        updateSelectInput(
          session,
          "V.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("v_gene"))
      }
    })
    # D gene
    observe({
      if (input$datasource == "ImmunoSEQ") {
        updateSelectInput(
          session,
          "D.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("dFamilyName"))
        
      }
      else if (input$datasource == "MiXCR") {
        
        updateSelectInput(
          session,
          "D.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("bestDHit"))
        
      }
      
      else if (input$datasource == "10x_scSeq") {
        
        updateSelectInput(
          session,
          "D.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("d_gene_B"))
        
      }
      
      else {
        updateSelectInput(
          session,
          "D.GENE.clean",
          choices=names(TSV.col.names()),
          selected = c("d_gene"))
      }

    })
    # amino acid 
    observe({
      if (input$datasource == "ImmunoSEQ") {
      updateSelectInput(
        session,
        "CDR3.gene.clean",
        choices=names(TSV.col.names()),
        
        selected = c("aminoAcid"))
        
        }
        else if (input$datasource == "MiXCR") {
          updateSelectInput(
            session,
            "CDR3.gene.clean",
            choices=names(TSV.col.names()),
          selected = c("aaSeqCDR3"))
        }
      
      
      else if (input$datasource == "10x_scSeq") {
        updateSelectInput(
          session,
          "CDR3.gene.clean",
          choices=names(TSV.col.names()),
          selected = c("cdr3_amino_acid_B"))
      }
      else {
        updateSelectInput(
          session,
          "CDR3.gene.clean",
          choices=names(TSV.col.names()),
          selected = c("JUNCTION"))
      }
    })
    # columns to remove
    observe({
      
      if (input$datasource == "ImmunoSEQ") {
      
      updateSelectInput(
        session,
        "col.to.remove",
        choices=names(TSV.col.names()),
        
        
          selected = c("count..templates.reads.","frequencyCount...."
          # c("product_subtype","frame_type","total_dj_reads",
          #              "productive_entropy","rearrangement_type",
          #              "order_name","release_date",
          #              "upload_date","primer_set","cdr3_length","frequency",
          #              "total_outofframe_reads","sample_catalog_tags","sample_rich_tags_json",
          #              "fraction_productive","sample_tags","sku","total_templates",
          #              "sequence_result_status","productive_clonality","stop_rearrangements",
          #              "outofframe_rearrangements","total_rearrangements","total_reads","sample_cells","fraction_productive_of_cells_mass_estimate", "sample_cells_mass_estimate","sample_amount_ng",
          #              "productive_rearrangements","counting_method","v_allele_ties","v_gene_ties","antibody",
          #              "sample_clonality", "max_productive_frequency","sample_entropy","sample_simpson_clonality",
          #              "max_frequency","productive_simpson_clonality","total_stop_reads","total_productive_reads",
          #              "v_deletions",	"d5_deletions",	"d3_deletions",	"j_deletions",	"n2_insertions",
          #              "n1_insertions",	"v_index",	"n1_index",	"n2_index",	"d_index",	"j_index",	"v_family_ties",	"d_family_ties",	"d_gene_ties",	"d_allele_ties",	"j_gene_ties"
          )
        )
        }

      else {
        
        updateSelectInput(
          session,
          "col.to.remove",
          choices=names(TSV.col.names()),
          selected = c()
        )
        
      }
      
    })
    observe({
      
      updateSelectInput(
        session,
        "V.GENE.clean2",
        choices=names(TSV.col.names()),
        selected = c("v_gene_A")
      )
      
    })
    observe({
      
      updateSelectInput(
        session,
        "J.GENE.clean2",
        choices=names(TSV.col.names()),
        selected = c("j_gene_A")
      )
      
    })
    observe({
      
      updateSelectInput(
        session,
        "CDR3.gene.clean2",
        choices=names(TSV.col.names()),
        selected = c("cdr3_amino_acid_A")
      )
    # in frame for immunoseq  
    })
    observe({
      updateSelectInput(
        session,
        "inframe_immseq",
        choices=names(TSV.col.names()),
        selected = c("sequenceStatus")
      )
      
    })

  TSV.file.Immunoseq <- reactive({
    x <- as.data.frame(input.data.Immunoseq())
    validate(
      need(nrow(x)>0,
           "upload file")
    )
    x2 <- x %>%
      select_if(~ !any(is.na(.)))
    #ImmunoSEQ 
    if (input$datasource == "ImmunoSEQ") {
      x2 <- subset(x2, x2[names(x2) %in% input$inframe_immseq]=="In")
      x2 <- data.frame(cloneCount = x2[,names(x2) %in% input$countcolumn], x2)
      x2$group <- input$group.imm
      x2$Indiv <- input$indiv.imm
      x2$Indiv.group <- paste(x2$group,x2$Indiv,sep=".")
      names(x2)[1] <- "cloneCount"
      head(x2)
      x3 <- x2

      # x3 <- x3[!is.na(x3[names(x3) %in% c("v_gene","j_gene",input$V.GENE.clean,input$J.GENE.clean)]),]
      
      x3$TRV <- x3[,names(x3) %in% input$V.GENE.clean]
      x3$TRV <- gsub("^TCR","",x3$TRV)
      
      x3$TRV[x3$TRV == ''] <- "TRV"
      
      # x3$TRV <- gsub("","TRA",x3$TRV)
      
      x3$TRJ <- x3[,names(x3) %in% input$J.GENE.clean]
      x3$TRJ <- gsub("^TCR","",x3$TRJ)
      x3$TRJ[x3$TRJ == ''] <- "TRJ"
      # x3$TRJ <- gsub("","TRJ",x3$TRJ)
      if (input$D_chain_present == "yes") {
        x3$TRD <- x3[,names(x3) %in% input$D.GENE.clean]
        # x3[is.na(x3$TRD)] <- "-"
        x3$TRD <- gsub("^TCR","",x3$TRD)
        x3$TRD[x3$TRD == ''] <- "-"
      }
     
      
      x3 <- x3[!is.na(x3[names(x3) %in% c("TRV")]),]
      x3 <- x3[!is.na(x3[names(x3) %in% c("TRJ")]),]
      
      x3$TRVJ <- paste(x3$TRV,x3$TRJ,sep=".")
      if (input$D_chain_present == "yes") {
      x3$TRVDJ <- paste(x3$TRV,x3$TRD,x3$TRJ,sep=".")
      }
      # x3$TRVDJ <- gsub(".-.",".",x3$TRVDJ)
    
      x3$TRVJ_CDR3 <- paste(x3$TRVJ, x3[,names(x3) %in% input$CDR3.gene.clean],sep="_")
      if (input$D_chain_present == "yes") {
      x3$TRVDJ_CDR3 <- paste(x3$TRVDJ, x3[,names(x3) %in% input$CDR3.gene.clean],sep="_")
      }
      x3 <- x3[!names(x3) %in% input$col.to.remove]
      x3[x3 == ''] <- "Missing"
      # x3[is.na(x3)] <- "Missing"
      x3
    }
    
    # mixcr 
    else if (input$datasource == "MiXCR") {
      x2 <- data.frame(cloneCount = x2[,names(x2) %in% input$countcolumn], x2)
      names(x2)[1] <- "cloneCount"
      
      x3 <- x2
      
      x3 <- x3[!is.na(x3[names(x2) %in% c(input$V.GENE.clean,input$J.GENE.clean)]),]
      
      x3$TRV <- str_remove(x3[,names(x3) %in% input$V.GENE.clean], "\\*00")
      x3$TRV <- gsub("^TR","",x3$TRV)
      
      if (input$D_chain_present == "yes") {
        x3$TRD <- str_remove(x3[,names(x3) %in% input$D.GENE.clean], "\\*00")
        x3$TRD <- gsub("^TR","",x3$TRD)
      }
      
    
      
      x3$TRJ <- str_remove(x3[,names(x3) %in% input$J.GENE.clean], "\\*00")
      x3$TRJ <- gsub("^TR","",x3$TRJ)
      
      
      x3 <- x3[-c(grep("\\_",x3[,names(x3) %in% input$CDR3.gene.clean])),]
      x3 <- x3[-c(grep("\\*",x3[,names(x3) %in% input$CDR3.gene.clean])),]
      
      x3$TRVJ <- paste(x3$TRV,x3$TRJ,sep=".")
      if (input$D_chain_present == "yes") {
      x3$TRVDJ <- paste(x3$TRV,x3$TRD,x3$TRJ,sep=".")
      x3$TRVDJ <- gsub(".NA.",".",x3$TRVDJ)
      x3$TRD <- gsub("NA","-",x3$TRD)
      }
      
      x3$TRVJ_CDR3 <- paste(x3$TRVJ, x3[,names(x3) %in% input$CDR3.gene.clean],sep="_")
      if (input$D_chain_present == "yes") {
      x3$TRVDJ_CDR3 <- paste(x3$TRVDJ, x3[,names(x3) %in% input$CDR3.gene.clean],sep="_")
      }
      
      x3 <- x3[!names(x3) %in% c(input$col.to.remove,"cloneCount.1")]
      x3

    }
    # paired 10x scSeq
    
    else if (input$datasource == "10x_scSeq") {

      x2$cloneCount <- 1
      contigs <- x2
      contigs <- subset(contigs,contigs$productive==T)
      contigs_lim <- contigs[!names(contigs) %in% c("is_cell","contig_id","high_confidence","raw_consensus_id","exact_subclonotype_id","umis","reads","length","cdr3_nt",names(contigs[grep("fwr",names(contigs))]),names(contigs[grep("cdr1",names(contigs))]),names(contigs[grep("cdr2",names(contigs))])
      )]
      contigs_lim
      contig_AG <- subset(contigs_lim,contigs_lim$chain=="TRA" | contigs_lim$chain=="TRG")
      name.list <- names(contig_AG[c(names(contig_AG[grep("gene",names(contig_AG))]),
                                     names(contig_AG[grep("cdr3",names(contig_AG))]),
                                     "chain")])
      contig_AG <- contig_AG %>%
        select(all_of(name.list), everything())
      names(contig_AG)[1:summary(name.list)[1]] <-paste(names(contig_AG[names(contig_AG) %in% name.list]),"_AG",sep="")
      contig_AG
      
      contig_BD <- subset(contigs_lim,contigs_lim$chain=="TRB" | contigs_lim$chain=="TRD")
      name.list <- names(contig_BD[c(names(contig_BD[grep("gene",names(contig_BD))]),
                                     names(contig_BD[grep("cdr3",names(contig_BD))]),
                                     "chain")])
      contig_BD <- contig_BD %>%
        select(all_of(name.list), everything())
      
      
      names(contig_BD)[1:summary(name.list)[1]] <-paste(names(contig_BD[names(contig_BD) %in% name.list]),"_BD",sep="")
      contig_BD
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      
      contig_paired$pairing <- ifelse(contig_paired$chain_BD=="TRB" & contig_paired$chain_AG=="TRA","abTCR Paired",
                                      ifelse(contig_paired$chain_BD=="TRD" & contig_paired$chain_AG=="TRG","gdTCR Paired"
                                      ))
      contig_paired
      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
      contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_AG")]
      contig_paired_only <- contig_paired
      contig_paired_only$d_gene_BD <- sub("^$","NA", contig_paired_only$d_gene_BD)
      # 
      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG,contig_paired_only$j_gene_AG,sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("NA.NA","",contig_paired_only$vj_gene_AG)
      # 
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.",".",contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("NA.NA","",contig_paired_only$vj_gene_BD)
      # 
      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$d_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub(".NA.",".",contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("NA.NA","",contig_paired_only$vdj_gene_BD)
      # 
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$cdr3_AG,sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_AG)
      # 
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_BD)
      # 
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vdj_gene_cdr3_BD)
      # 
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vdj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_AG_BD)
      # 
      # #updating names to be consistant.... 
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vj_gene_cdr3_AG_BD)
      
      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_cdr3_AG_BD)
      # contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      names(contig_paired_only)[names(contig_paired_only) %in% "barcode"] <- "Cell_Index"
      contig_paired_only <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index),] # remove duplicates
      # contig_paired_only <- contig_paired_only %>%
      #   select(all_of(c("Cell_Index","Indiv.group")), everything())
      
      x3 <- contig_paired_only
      x3 <- x3[!names(x3) %in% input$col.to.remove]
      x3$group <- input$group.imm
      x3$indiv <- input$indiv.imm
      x3$Indiv.group <- paste(x3$group,x3$indiv,sep=".")
      x3
      # x2 <- data.frame(cloneCount = x2[,names(x2) %in% input$countcolumn], x2)
      # names(x2)[1] <- "cloneCount"
      # 
      # for (i in 1:dim(x2)[2]) {
      #   x2[,i]   <- gsub(";",",",x2[,i])
      # }
      # 
      # x2$group <- input$group.imm
      # x2$indiv <- input$indiv.imm
      # x2$Indiv.group <- paste(x2$group,x2$indiv,sep=".")
      # names(x2)[1] <- "cloneCount"
      # 
      # x3 <- x2
      # 
      # df_name <- as.data.frame(do.call(rbind, strsplit(as.character(x3[,names(x3) %in% input$V.GENE.clean2]), ",")))
      # df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(x3[,names(x3) %in% input$J.GENE.clean2]), ",")))
      # 
      # x3$AV <- df_name[,1]
      # x3$AJ <- df_name2[,1]
      # x3$AJ <- gsub("[*]0.","",x3$AJ)
      # x3$AJ <- gsub(",, AJ..","",x3$AJ)
      # x3$AJ <- gsub(",, AJ.","",x3$AJ)
      # 
      # x3$AVJ <- paste(x3$AV,".",x3$AJ,sep="")
      # x3$AV <- gsub("[*]0.","",x3$AV)
      # x3$AVJ <- gsub("[*]0.","",x3$AVJ)
      # 
      # df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(x3[,names(x3) %in% input$V.GENE.clean]), ",")))
      # df_name4 <- as.data.frame(do.call(rbind, strsplit(as.character(x3[,names(x3) %in% input$J.GENE.clean]), ",")))
      # df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(x3[,names(x3) %in% input$D.GENE.clean]), ",")))
      # 
      # x3$BV <- df_name3[,1]
      # x3$BJ <- df_name4[,1]
      # x3$BD <- df_name5[,1]
      # 
      # 
      # x3$BVJ <- paste(x3$BV,".",x3$BJ,sep="")
      # x3$BVDJ <- paste(x3$BV,".",x3$BD,".",x3$BJ,sep="")
      # 
      # x3$BV <- gsub("[*]0.","",x3$BV)
      # x3$BJ <- gsub("[*]0.","",x3$BJ)
      # x3$BD <- gsub("[*]0.","",x3$BD)
      # x3$BD <- gsub(" ","",x3$BD)
      # x3$BVJ <- gsub("[*]0.","",x3$BVJ)
      # x3$BVDJ <- gsub("[*]0.","",x3$BVDJ)
      # x3$BVDJ <- gsub(".NA.",".",x3$BVDJ)
      # 
      # x3$AJ <- gsub("TR","",x3$AJ)
      # x3$AVJ <- gsub("TR","",x3$AVJ)
      # x3$AVJ <- gsub("AJ","J",x3$AVJ)
      # x3$AVJ.BVJ <- paste(x3$AVJ,"_",x3$BVJ,sep="")
      # x3$AVJ.BVDJ <- paste(x3$AVJ,"_",x3$BVDJ,sep="")
      # 
      # x3$AVJ_aCDR3 <- paste(x3$AVJ,x3[,names(x3) %in% input$CDR3.gene.clean2],sep="_")
      # x3$BVJ_bCDR3 <- paste(x3$BVJ,x3[,names(x3) %in% input$CDR3.gene.clean],sep="_")
      # x3$BVDJ_bCDR3 <- paste(x3$BVDJ,x3[,names(x3) %in% input$CDR3.gene.clean],sep="_")
      # 
      # x3$AVJ_aCDR3_BVJ_bCDR3 <- paste(x3$AVJ_aCDR3,x3$BVJ_bCDR3,sep=" & ")
      # x3$AVJ_aCDR3_BVDJ_bCDR3 <- paste(x3$AVJ_aCDR3,x3$BVDJ_bCDR3,sep=" & ")
      # x3$BD <- gsub("NA","-",x3$BD)
      # 
      # x3 <- x3[!names(x3) %in% input$col.to.remove]
      # x3[is.na(x3)] <- "Missing"
      # 
      # x3 <- subset(x3, !x3$AJ=="Missing")
      # x3 <- subset(x3, !x3$BJ=="Missing")
      # x3[is.na(x3)] <- " "
      
    }
    
    # other data 
    else {
      
      x2$InFrame <- tolower(x2[,names(x2) %in% input$inframe_immseq])
      # x2 <- subset(x2, x2[names(x2) %in% input$inframe_immseq]==TRUE)
      
      if (nrow(x2[-c(grep("false",x2$InFrame)),]>0)) {
        x2 <- x2[-c(grep("false",x2$InFrame)),]
      }
      
      if (nrow(x2[-c(grep("[*]",x2$junction_aa)),]>0)) {
        x2 <- x2[-c(grep("[*]",x2$junction_aa)),]
      }
      
      if (nrow(x2[-c(grep("",x2$junction_aa)),]>0)) {
        x2 <- x2[-c(grep("",x2$junction_aa)),]
      }
      
      x2 <- data.frame(cloneCount = x2[,names(x2) %in% input$countcolumn], x2)
      names(x2)[1] <- "cloneCount"

      x3 <- x2
      x3
      x3 <- x3[!names(x3) %in% input$col.to.remove]
      # x3 <- x3[!is.na(x3[names(x3) %in% c(input$V.GENE.clean,input$J.GENE.clean)]),]
      x3$TRJ <- x3[,names(x3) %in% input$J.GENE.clean]
      x3$TRV <- x3[,names(x3) %in% input$V.GENE.clean]
      
      if (input$D_chain_present == "yes") {
        x3$TRD <- x3[,names(x3) %in% input$D.GENE.clean]
      }

      x3$TRVJ <- paste(x3$TRV,x3$TRJ,sep=".")
      
      if (input$D_chain_present == "yes") {
      x3$TRVDJ <- paste(x3$TRV,x3$TRD,x3$TRJ,sep=".")
      x3$TRVDJ <- gsub(".NA.",".",x3$TRVDJ)
      x3$TRD <- gsub("NA","-",x3$TRD)
}
      x3$TRVJ_CDR3 <- paste(x3$TRVJ, x3[,names(x3) %in% input$CDR3.gene.clean],sep="_")
      if (input$D_chain_present == "yes") {
      x3$TRVDJ_CDR3 <- paste(x3$TRVDJ, x3[,names(x3) %in% input$CDR3.gene.clean],sep="_")
}
      x3 <- x3[!names(x3) %in% c(input$col.to.remove,"cloneCount.1")]
      x3 <- subset(x3,x3$TRV!="None")
      # x3 <- x3[-c(grep("\\_",x3[,names(x3) %in% input$CDR3.gene.clean])),]
      # x3 <- x3[-c(grep("\\*",x3[,names(x3) %in% input$CDR3.gene.clean])),]
    }

  x3
  
  })

  output$ImmunoSeq.table <- DT::renderDataTable(escape = FALSE, filter = "top", options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    TSV.file.Immunoseq()
  })
  
  output$downloadTABLE.Immunoseq <- downloadHandler(
    filename = function(){
      paste(input$IMGT_name_df," TCR_Explore.analysis.file-",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- TSV.file.Immunoseq()
      df <- as.data.frame(df)
      write.csv(df,file, row.names = FALSE)
    } )
  
  # summarised table -----
  output$names.in.file3 <- renderPrint( {
    df <- input.data2()
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    names(df)
    
  })
  
  observe({
    
    if (input$type.chain == 'ab') {
      updateSelectInput(
        session,
        "string.data3",
        choices=names(input.data_old2()),
        selected = c("group","JUNCTION_A"))
    }
    
    else if (input$type.chain == 'NGS_immunoseq') {
      updateSelectInput(
        session,
        "string.data3",
        choices=names(input.data_old2()),
        selected = c("group","aminoAcid"))
    }
    
    else {
      updateSelectInput(
        session,
        "string.data3",
        choices=names(input.data_old2()),
        selected = c("group","JUNCTION_G")) 
    }
  }) 
  
  chain_table_summary <- reactive({
    df <- input.data_old2()
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df2 <- df[,c("cloneCount",input$string.data3)] 
    df2
    df3 <- as.data.frame(ddply(df2,input$string.data3,numcolwise(sum)))
    df3
  })
  
  chain_table_summary.TCRdist3.ab <- reactive({
    df <- input.data_old2()
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    
    df.Vgene_A <- as.data.frame(do.call(rbind, strsplit(as.character(df$V.GENE.and.allele_A), ",")))
    df$V.GENE.and.allele_A <- df.Vgene_A$V1
    df.Vgene_B <- as.data.frame(do.call(rbind, strsplit(as.character(df$V.GENE.and.allele_B), ",")))
    df$V.GENE.and.allele_B <- df.Vgene_B$V1
    df.Jgene_A <- as.data.frame(do.call(rbind, strsplit(as.character(df$J.GENE.and.allele_A), ",")))
    df$J.GENE.and.allele_A <- df.Jgene_A$V1
    
    
    df$V.GENE.and.allele_A <- gsub("AV","TRAV", df$V.GENE.and.allele_A)
    df$V.GENE.and.allele_A <- gsub("DV","TRDV", df$V.GENE.and.allele_A)
    df$V.GENE.and.allele_B <- gsub("BV","TRBV", df$V.GENE.and.allele_B)
    
    df$J.GENE.and.allele_A <- gsub("AJ","TRAJ", df$J.GENE.and.allele_A)
    df$J.GENE.and.allele_B <- gsub("BJ","TRBJ", df$J.GENE.and.allele_B)
    
    df2 <- df[,c("Indiv","group","cloneCount","V.GENE.and.allele_A","J.GENE.and.allele_A","JUNCTION..AA._A","JUNCTION_A","V.GENE.and.allele_B","J.GENE.and.allele_B","JUNCTION..AA._B","JUNCTION_B")] 
    
    names(df2) <- c("subject","epitope","count","v_a_gene","j_a_gene","cdr3_a_aa","cdr3_a_nucseq","v_b_gene","j_b_gene","cdr3_b_aa","cdr3_b_nucseq")
    
    df2$well_id <- paste("well")
  
    
    df3 <- as.data.frame(ddply(df2,c("subject","epitope","v_a_gene","j_a_gene","cdr3_a_aa","cdr3_a_nucseq","v_b_gene","j_b_gene","cdr3_b_aa","cdr3_b_nucseq","well_id"),numcolwise(sum)))
    
    df3$well_id <- paste("Clone_",rownames(df3),".",df3$subject,"-",df3$epitope,sep="")
    
    df3 <- df3[,c(1,2,12,3:11)]
    
    df3
    
  })
  
  chain_table_summary.TCRdist3.gd <- reactive({
    df <- input.data.IMGT_afterQC()
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    
    df.Vgene_G <- as.data.frame(do.call(rbind, strsplit(as.character(df$V.GENE.and.allele_G), ",")))
    df$V.GENE.and.allele_G <- df.Vgene_G$V1
    df.Vgene_D <- as.data.frame(do.call(rbind, strsplit(as.character(df$V.GENE.and.allele_D), ",")))
    df$V.GENE.and.allele_D <- df.Vgene_D$V1
    df.Jgene_D <- as.data.frame(do.call(rbind, strsplit(as.character(df$J.GENE.and.allele_D), ",")))
    df$J.GENE.and.allele_D <- df.Jgene_D$V1
    
    
    df$V.GENE.and.allele_D <- gsub("DV","TRDV", df$V.GENE.and.allele_D)
    df$V.GENE.and.allele_G <- gsub("GV","TRGV", df$V.GENE.and.allele_G)
    
    df$J.GENE.and.allele_D <- gsub("DJ","TRDJ", df$J.GENE.and.allele_D)
    df$J.GENE.and.allele_G <- gsub("GJ","TRGJ", df$J.GENE.and.allele_G)
    
    
    df2 <- df[,c("Indiv","group","cloneCount","V.GENE.and.allele_G","J.GENE.and.allele_G","JUNCTION..AA._G","JUNCTION_G","V.GENE.and.allele_D","J.GENE.and.allele_D","JUNCTION..AA._D","JUNCTION_D")] 
    
    names(df2) <- c("subject","epitope","count","v_g_gene","j_g_gene","cdr3_g_aa","cdr3_g_nucseq","v_d_gene","j_d_gene","cdr3_d_aa","cdr3_d_nucseq")
    
    df2$well_id <- paste("well")
    df3 <- as.data.frame(ddply(df2,c("subject","epitope","v_g_gene","j_g_gene","cdr3_g_aa","cdr3_g_nucseq","v_d_gene","j_d_gene","cdr3_d_aa","cdr3_d_nucseq","well_id"),numcolwise(sum)))
    
    df3$well_id <- paste("Clone_",rownames(df3),".",df3$subject,"-",df3$epitope,sep="")
    
    df3 <- df3[,c(1,2,12,3:11)]
    df3
    
  })
  
  output$chain_table_IMGT.QC3 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    df1 <- input.data.IMGT_afterQC();
    df1 <- as.data.frame(df1)
    
    if (input$type.of.graph == "general summary") {
      chain_table_summary()
    }
    
    else if (input$type.of.graph == "TCRdist3" && input$type.chain == "ab") {
      chain_table_summary.TCRdist3.ab()
    }
    
    
    else if (input$type.of.graph == "TCRdist3" && input$type.chain == "gd") {
      chain_table_summary.TCRdist3.gd()
    }
    
    
    else {
      
      chain_table_summary.TCRdist3()
      
    }
    
  })
  
  # summary table download file -----
  output$downloadTABLE.QC3 <- downloadHandler(
    filename = function(){
      if (input$type.of.graph == "general summary") {
        paste("paired_chain_CDR3",gsub("-", ".", Sys.Date()),".csv", sep = "")
      }
      
      else if (input$type.of.graph == "TCRdist3" && input$type.chain == "ab") {
        paste("paired_TCRdist.ab",gsub("-", ".", Sys.Date()),".csv", sep = "")
      }
      
      else if (input$type.of.graph == "TCRdist3" && input$type.chain == "gd") {
        paste("paired_TCRdist.gd",gsub("-", ".", Sys.Date()),".csv", sep = "")
      }
      
      
      else {
        paste("paired_TCRdist",gsub("-", ".", Sys.Date()),".csv", sep = "")
      }
      
    },
    content = function(file){
      
      if (input$type.of.graph == "general summary") {
        df <- chain_table_summary()
        write.csv(df,file, row.names = FALSE)
      }
      
      else if (input$type.of.graph == "TCRdist3" && input$type.chain == "ab") {
        df <- chain_table_summary.TCRdist3.ab()
        write.csv(df,file, row.names = FALSE)
      }
      
      else if (input$type.of.graph == "TCRdist3" && input$type.chain == "gd") {
        df <- chain_table_summary.TCRdist3.gd()
        write.csv(df,file, row.names = FALSE)
      }
      
      
      else {
        df <- chain_table_summary()
        write.csv(df,file, row.names = FALSE)
      }
      
      
    })
  
  
  # merging Multiple TCR_Explore files
  # TCRex Merge  ------
  input.data_TCRexMerge <- reactive({
    inFile2_TCRexMerge <- input$file2_TCRexMerge
    num <- dim(inFile2_TCRexMerge)[1]
    samples_list <- vector("list", length = num)
    samples_list
    for (i in 1:num) {
      sc <- read.csv(input$file2_TCRexMerge[[i, 'datapath']],header = T)
      samples_list[[i]] <- sc
    }
    samples_list
  })
  
  merged_TCRexFiltered <- reactive({
    inFile2_TCRexMerge <- input$file2_TCRexMerge
    validate(
      need(nrow(inFile2_TCRexMerge)>0,
           "Upload 2 or more files")
    )
    
    myfiles <- input.data_TCRexMerge()
    df <- rbind(myfiles[[1]])
    for (i in 2:length(myfiles)) {
      df <- rbind(df,myfiles[[i]])
    }
    
    df
  })
  
  output$DEx_TCRexFiltered <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    df <- merged_TCRexFiltered()
    df
  })
  
  
  output$downloaddf_TCRexFiltered <- downloadHandler(
    filename = function(){
      paste("TCR_Explore_merged",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- merged_TCRexFiltered()
      write.csv(df,file, row.names = F)
    })
  
  
  # file for analytical plots -----
  input.data_old2 <- reactive({switch(input$dataset,"ab-test-data2" = test.data2_TCR.Explore(), "ImmunoSEQ-test-data" = test.data2_ImmunoSEQ(),"own_data2" = own.data2())})
  test.data2_TCR.Explore <- reactive({
    # dataframe = read.csv("test-data/Group/paired_unsummarised2021.09.22.csv",header=T) 
    dataframe = read.csv("test-data/Group/paired_TCR_file2022.05.24.csv",header=T) 
  })
  
  test.data2_ImmunoSEQ <- reactive({
    # dataframe = read.csv("test-data/Group/paired_unsummarised2021.09.22.csv",header=T) 
    dataframe = read.csv("test-data/Group/ImmunoSEQ.test.TCR_Explore.analysis.file-2022.08.04.csv",header=T) 
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
  # filtering out sequences 
  observe({
    
    
    
    updateSelectInput(
      session,
      "clonotypes_column",
      choices=names(input.data_old2()),
      selected = names(input.data_old2())[dim(input.data_old2())[2]] 
      # selected = "vdj_gene_AG_BD"
      )
  })
  observe({
    updateSelectInput(
      session,
      "group_column_first",
      choices=names(input.data_old2()),
      selected = "Indiv" )
  })
  
  filtering <- reactive({
    
    dataframe1 <- input.data_old2()
    validate(
      need(nrow(dataframe1)>0,
           "Upload file")
    )
    
    totals <- dataframe1[,names(dataframe1) %in% c("cloneCount",input$clonotypes_column,input$group_column_first)]
    total.condition <- as.data.frame(ddply(totals,c(input$clonotypes_column,input$group_column_first),numcolwise(sum)))
    totals <- dataframe1[,names(dataframe1) %in% c("cloneCount",input$group_column_first)]
    indiv.sum <-  as.data.frame(ddply(totals,c(input$group_column_first),numcolwise(sum)))
    names(indiv.sum)[names(indiv.sum) %in% input$group_column_first] <- "Var1"
    names(total.condition)[names(total.condition) %in% input$group_column_first] <- "Var1"
    emtpy <- matrix(nrow =dim(total.condition)[1],ncol=dim(indiv.sum)[1])
    for (j in 1:dim(total.condition)[1]){
      for (i in 1:dim(indiv.sum)[1]){
        emtpy[j,i] <- ifelse(indiv.sum$Var1[i]==total.condition$Var1[j],indiv.sum$cloneCount[i],F)
      }
    }
    total.condition$n <- total.condition$cloneCount/rowSums(emtpy)
    total.condition <- total.condition[order(total.condition$n, decreasing = T),]
    total.condition <- total.condition[order(total.condition$Var1, decreasing = T),]
    p <- list()
    kl <- list()
    indiv_uniq <- unique(total.condition$Var1)
    indiv_sub <- subset(total.condition,total.condition$Var1==indiv_uniq[3])
    for (i in 1:dim(indiv.sum)[1]) {
      indiv_sub <- subset(total.condition,total.condition$Var1==indiv_uniq[i])
      indiv_sub <- indiv_sub[order(indiv_sub$n, decreasing = T),]
      p[[i]]  <- cumsum(indiv_sub$n)
      kl[[i]] <- 1:dim(indiv_sub)[1]
    }
    total.condition$cumsum <- unlist(p)
    total.condition$order <- unlist(kl)
    if (input$filter_by_cumsum_or_order =="cummulative_freq") {
      total.condition_sub <- subset(total.condition,total.condition$cumsum<input$filter_below_cumsum)
    }
    else {
      total.condition_sub <- subset(total.condition,total.condition$order <= input$top_clones_to_include)
    }
    total.condition_sub
   
    
    
  })
  
  input.data2 <- reactive({
    dataframe1 <- input.data_old2()
    validate(
      need(nrow(dataframe1)>0,
           "Upload file")
    )
    if(input$Filter_NGS == "yes") {
      total.condition_sub <- filtering()
      
      names(total.condition_sub)[names(total.condition_sub) %in% "cloneCount"] <- "total_count"
      names(total.condition_sub)[names(total.condition_sub) %in% "Var1"] <- input$group_column_first
      total.condition_sub <- total.condition_sub[,names(total.condition_sub) %in% c(input$clonotypes_column,input$group_column_first)]
      subsetted_file <- merge(dataframe1,total.condition_sub,by = c(input$clonotypes_column,input$group_column_first))
      subsetted_file
    }
    
    else {
      dataframe1
    }
    
    
  })
  
  output$Test_table <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE), {
    filtering()
  })
  
  # Treemap ------
  
  
  observe({
    updateSelectInput(
      session,
      "count2",
      choices=names(input.data_old2()),
      selected = "cloneCount")
  })
  
  observe({
    updateSelectInput(
      session,
      "fill2",
      choices=names(input.data_old2()),
      selected = "AVJ" )
  })
  observe({
    updateSelectInput(
      session,
      "sub_group2",
      choices=names(input.data_old2()),
      selected = "AVJ.BVJ")
    
  })
  observe({
    updateSelectInput(
      session,
      "string.data.tree.order",
      choices=select_group(),
      selected = select_group()) 
  }) 
  
  cols <- reactive({
    dat <- input.data2()
    dat <- as.data.frame(dat)
    num <- unique(dat[names(dat) %in% input$fill2])
    col.gg <- gg_fill_hue(dim(num)[1])
    unique.col <- as.data.frame(unique(dat[grep(input$fill2,names(dat))]))
    
    palette_rainbow <- rev(rainbow(dim(num)[1]))
    
    if (input$tree_colour.choise == "rainbow") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col", i, sep="_"), paste(num[i,]), palette_rainbow[i])        
      }) }
    
    else if (input$tree_colour.choise == "default") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col", i, sep="_"), paste(num[i,]), col.gg[i])        
      })
    }
    else if (input$tree_colour.choise == "random") {
      palette1 <- distinctColorPalette(dim(unique.col)[1])
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col", i, sep="_"), paste(num[i,]), palette1[i])        
      })
      
    }
    
    else {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col", i, sep="_"), paste(num[i,]), input$one.colour.default)        
      })
      
      
    }
    
  })
  output$myPanel <- renderUI({cols()})
  colors <- reactive({
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    num <- unique(dat[names(dat) %in% input$fill2])
    lapply(1:dim(num)[1], function(i) {
      input[[paste("col", i, sep="_")]]
    })
  })
  tree_plot_dynamic <- reactive({
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    
    if (is.null(input$col_1)) {
      cols <- rep("#000000", ncol(dat))
    } 
    else {
      cols <- unlist(colors())
    }
    
    if (input$tree.lab == "yes" & input$type.tree == "raw data") {
      df1 <- dat[names(dat) %in% c(input$count2,input$fill2,input$sub_group2,input$group_column)]
      df2 <- as.data.frame(ddply(dat,names(df1)[-c(1)],numcolwise(sum)))
      unique.col <- as.data.frame(unique(dat[names(dat) %in% input$fill2]))
      names(unique.col) <- "V1"
      unique.col$tree_palette <- cols
      
      
      
      df3 <- as.data.frame(merge(df2,unique.col,by.x=input$fill2,by.y = "V1"))
      
      
      df3$ID.names <- factor(df3[,names(df3) %in% input$group_column],levels = input$string.data.tree.order)
      
      vals22$Treemap22 <- ggplot(df3, aes(area = get(input$count2),
                                          fill = get(input$fill2),
                                          subgroup = get(input$sub_group2))) +
        geom_treemap(colour="white",show.legend = F, fill = df3$tree_palette) +
        geom_treemap_subgroup_border(colour = "white", show.legend = F,size=12) +
        geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 1, family = input$font_type,
                                   colour = input$Treemap.text.colour, fontface = "italic", min.size = 0,show.legend = F) +
        facet_wrap(~df3$ID.names,nrow = input$nrow.tree) +
        theme(strip.text = element_text(size = input$panel.text.size.tree, family = input$font_type))+
        theme(strip.background =element_rect(fill=input$strip.colour.tree))+
        theme(strip.text = element_text(colour = input$strip.text.colour.tree))
      
      vals22$Treemap22
      
    }
    else if (input$tree.lab == "no" & input$type.tree == "raw data") {
      df1 <- dat[names(dat) %in% c(input$count2,input$fill2,input$sub_group2,input$group_column)]
      df2 <- as.data.frame(ddply(dat,names(df1)[-c(1)],numcolwise(sum)))
      unique.col <- as.data.frame(unique(dat[names(dat) %in% input$fill2]))
      names(unique.col) <- "V1"
      unique.col$tree_palette <- cols
      df3 <- as.data.frame(merge(df2,unique.col,by.x=input$fill2,by.y = "V1"))
      
      df3$ID.names <- factor(df3[,names(df3) %in% input$group_column],levels = input$string.data.tree.order)
      
      vals22$Treemap22 <- ggplot(df3, aes(area = get(input$count2),
                                          fill = get(input$fill2),
                                          subgroup = get(input$sub_group2))) +
        geom_treemap(aes(alpha = 1),colour="white",show.legend = F, fill = df3$tree_palette) +
        geom_treemap_subgroup_border(colour = "white", show.legend = F,size=12) +
        facet_wrap(~df3$ID.names,nrow = input$nrow.tree) +
        theme(strip.text = element_text(size = input$panel.text.size.tree, family = input$font_type))+
        theme(strip.background =element_rect(fill=input$strip.colour.tree))+
        theme(strip.text = element_text(colour = input$strip.text.colour.tree))
      vals22$Treemap22
      
      
    }
    else if (input$tree.lab == "yes" & input$type.tree == "Summarised data") {
      df1 <- dat[names(dat) %in% c(input$count2,input$fill2,input$sub_group2,input$group_column)]
      unique.col <- as.data.frame(unique(dat[names(dat) %in% input$fill2]))
      names(unique.col) <- "V1"
      unique.col$tree_palette <- cols
      df3 <- as.data.frame(merge(df1,unique.col,by.x=input$fill2,by.y = "V1"), replace = FALSE)
      df3$ID.names <- factor(df3[,names(df3) %in% input$group_column],levels = input$string.data.tree.order)
      vals22$Treemap22 <- ggplot(df3, aes(area = get(input$count2),
                                          fill = get(input$fill2),
                                          subgroup = get(input$sub_group2))) +
        geom_treemap(aes(alpha = 1),colour="white",show.legend = F, fill = df3$tree_palette) +
        geom_treemap_subgroup_border(colour = "white", show.legend = F,size=12) +
        geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 1, family = input$font_type,
                                   colour = "black", fontface = "italic", min.size = 0,show.legend = F) +
        facet_wrap(~df3$ID.names,nrow = input$nrow.tree) +
        theme(strip.text = element_text(size = input$panel.text.size.tree, family = input$font_type))+
        theme(strip.background =element_rect(fill=input$strip.colour.tree))+
        theme(strip.text = element_text(colour = input$strip.text.colour.tree))
      vals22$Treemap22
      
    }
    else {
      df1 <- dat[names(dat) %in% c(input$count2,input$fill2,input$sub_group2,input$group_column)]
      unique.col <- as.data.frame(unique(dat[names(dat) %in% input$fill2]))
      names(unique.col) <- "V1"
      unique.col$tree_palette <- cols
      df3 <- as.data.frame(merge(df1,unique.col,by.x=input$fill2,by.y = "V1"), replace = FALSE)
      df3$ID.names <- factor(df3[,names(df3) %in% input$group_column],levels = input$string.data.tree.order)
      vals22$Treemap22 <- ggplot(df3, aes(area = get(input$count2),
                                          fill = get(input$fill2),
                                          subgroup = get(input$sub_group2))) +
        geom_treemap(aes(alpha = 1),colour="white",show.legend = F, fill = df3$tree_palette) +
        geom_treemap_subgroup_border(colour = "white", show.legend = F,size=12) +
        facet_wrap(~df3$ID.names,nrow = input$nrow.tree) +
        theme(strip.text = element_text(size = input$panel.text.size.tree, family = input$font_type))+
        theme(strip.background =element_rect(fill=input$strip.colour.tree))+
        theme(strip.text = element_text(colour = input$strip.text.colour.tree))
      vals22$Treemap22
      
    }
    
    
  })
  output$Treemap2 <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    tree_plot_dynamic()
  })
  # download Treeplots -----
  output$downloadPlot_scTREE <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_treemap_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_tree,height=input$height_tree, onefile = FALSE) # open the pdf device
      print(tree_plot_dynamic())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_scTREE <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_treemap_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_tree, height = input$height_png_tree, res = input$resolution_PNG_tree)
      print(tree_plot_dynamic())
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  # chord plot =====
  observe({
    updateSelectInput(
      session,
      "group_column",
      choices=names(input.data_old2()),
      selected = "group")
    
  }) # group 
  select_group <- reactive({
    df <- input.data_old2();
    
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    
    df2 <- as.data.frame(unique(df[names(df) %in% input$group_column]))
    df2 <- as.data.frame(df2)
    #names(df2) <- "V1"
    as.list(df2)
  })
  
  selected_chain_1 <- function () {
    df <- input.data_old2();
    
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    
    df2 <- as.data.frame(unique(df[names(df) %in% input$chain1]))
    df2 <- as.data.frame(df2)
    #names(df2) <- "V1"
    df2
  }
  
  
  observe({
    updateSelectInput(
      session,
      "group_selected2",
      choices= select_group(),
      selected = select_group()[1]
      )
    
  }) # group 
  observe({
    updateSelectInput(
      session,
      "chain1",
      choices=names(input.data_old2()),
      selected = "AV")
    
  }) # chain 1
  observe({
    updateSelectInput(
      session,
      "chain2",
      choices=names(input.data_old2()),
      selected = "AJ")
    
  }) # chain 2
  output$table_display <- renderTable({
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    dat <- subset(dat, get(input$group_column)==input$group_selected2)
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    hierarchy <- hierarchy %>%
      select(input$chain1, everything())
    
    head(hierarchy, n=6)
  })
  cols_circ <- reactive({
    set.seed(input$seed.numb.chord)
    dat <- input.data2()
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    # dat <- subset(dat, get(input$group_column)==input$group_selected2)
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    hierarchy <- hierarchy[,c(input$chain1,input$chain2)]
    hierarchy <- hierarchy %>%
      select(input$chain1, everything())
    df.col1 <- as.data.frame(unique(hierarchy[,1]))
    names(df.col1) <- "V1"
    
    df.col.j <- as.data.frame(unique(hierarchy[,2]))
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    palette_rainbow <- rev(rainbow(length(t(df.col.2))))
    if (input$colour_cir == "rainbow") {
      lapply(1:dim(df.col.2)[1], function(i) {
        colourInput(paste("col.cir", i, sep="_"), paste(df.col.2[i,]), palette_rainbow[i])        
      }) }
    
    else if (input$colour_cir == "default") {
      lapply(1:dim(df.col.2)[1], function(i) {
        col.gg <- gg_fill_hue(dim(df.col.2)[1])
        colourInput(paste("col.cir", i, sep="_"), paste(df.col.2[i,]), col.gg[i])        
      })
    }
    
    
    else if (input$colour_cir == "random") {
      
      lapply(1:dim(df.col.2)[1], function(i) {
        palette2 <- distinctColorPalette(dim(df.col.2)[1])
        colourInput(paste("col.cir", i, sep="_"), paste(df.col.2[i,]), palette2[i])        
      }) }
    
    else  {
      lapply(1:dim(df.col.2)[1], function(i) {
        colourInput(paste("col.cir", i, sep="_"), paste(df.col.2[i,]), input$one.colour.default)        
      }) }
    
  })
  output$myPanel_circ <- renderUI({cols_circ()})
  
  
  observe({
    updateSelectInput(
      session,
      "string.data.circ.order",
      choices=selected_chain_1(),
      selected = c("AV4","AV22","AV19")) 
  }) 
  
  
  colors_cir <- reactive({
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    # dat <- subset(dat, get(input$group_column)==input$group_selected2)
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    hierarchy <- hierarchy[,c(input$chain1,input$chain2)]
    hierarchy <- hierarchy %>%
      select(input$chain1, everything())
    df.col1 <- as.data.frame(unique(hierarchy[,1]))
    names(df.col1) <- "V1"
    
    df.col.j <- as.data.frame(unique(hierarchy[,2]))
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    df.col.2
    
    lapply(1:dim(df.col.2)[1], function(i) {
      input[[paste("col.cir", i, sep="_")]]
    })
  })
  col.table1 <- reactive({
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    hierarchy <- hierarchy[,c(input$chain1,input$chain2)]
    hierarchy <- hierarchy %>%
      select(input$chain1, everything())
    df.col1 <- as.data.frame(unique(hierarchy[,1]))
    names(df.col1) <- "V1"
    df.col.j <- as.data.frame(unique(hierarchy[,2]))
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    length(t(df.col.2))
    col2 <- unlist(colors_cir())
    
    if (input$circ_lab =="colour selected clone/s (label)"|input$circ_lab =="colour selected clone/s (no label)") {
      my_col_alpha_all <- col2
      
      for(i in 1:length(col2)) {
        my_col_alpha_all[i] <- ifelse(df.col.2[i,1] %in% c(input$string.data.circ.order),
                                      adjustcolor(col2[i], alpha.f =input$selected.chord.transparacy),
                                      adjustcolor(col2[i], alpha.f =input$unselected.chord.transparacy))
      }
      
      df.col.2$colour <- my_col_alpha_all
      grid.col <- as.data.frame(as.matrix(t(as.data.frame(df.col.2$colour))))
      names(grid.col) <- df.col.2$V1
      grid.col <- as.data.frame(grid.col)
      grid.col
      
    }
    
    else {
      df.col.2$colour <- col2
      grid.col <- as.data.frame(as.matrix(t(as.data.frame(df.col.2$colour))))
      names(grid.col) <- df.col.2$V1
      grid.col <- as.data.frame(grid.col)
      grid.col
      
    }
    

    
  })
  
  
  output$colour.trans.test <- renderPlot({
    dat <- input.data2()
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    hierarchy <- hierarchy[,c(input$chain1,input$chain2)]
    hierarchy <- hierarchy %>%
      select(input$chain1, everything())
    
    df.col1 <- as.data.frame(unique(hierarchy[,1]))
    names(df.col1) <- "V1"
    df.col.j <- as.data.frame(unique(hierarchy[,2]))
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    length(t(df.col.2))
    col2 <- unlist(colors_cir())
    
    if (input$circ_lab =="colour selected clone/s (label)" |input$circ_lab =="colour selected clone/s (no label)" ) {
      my_col_alpha_all <- col2

      for(i in 1:length(col2)) {
        my_col_alpha_all[i] <- ifelse(df.col.2[i,1] %in% c(input$string.data.circ.order),
                                      adjustcolor(col2[i], alpha.f =input$selected.chord.transparacy),
                                      adjustcolor(col2[i], alpha.f =input$unselected.chord.transparacy))
      }


      show_col(my_col_alpha_all)

    }

    else {
    show_col(col2)
    }
  })

  output$out.col.table1 <- renderTable({
    dat <- input.data2()
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    hierarchy <- hierarchy[,c(input$chain1,input$chain2)]
    hierarchy <- hierarchy %>%
      select(input$chain1, everything())
    df.col1 <- as.data.frame(unique(hierarchy[,1]))
    names(df.col1) <- "V1"
    df.col.j <- as.data.frame(unique(hierarchy[,2]))
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    length(t(df.col.2))
    col2 <- unlist(colors_cir())
    df.col.2$colour <- col2
    grid.col <- as.data.frame(as.matrix(t(as.data.frame(df.col.2$colour))))
    names(grid.col) <- df.col.2$V1
    grid.col <- as.data.frame(grid.col)
    grid.col
    dat2 <- subset(dat, get(input$group_column)==input$group_selected2)
    hierarchy2 <- dat2[names(dat2) %in% c(input$chain1,input$chain2)]
    hierarchy2 <- hierarchy2 %>%
      select(input$chain1, everything())
    
    df.col1 <- as.data.frame(unique(hierarchy2[,1]))
    names(df.col1) <- "V1"
    df.col.j <- as.data.frame(unique(hierarchy2[,2]))
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    df.col.2
    #df.col.2$V1 <- gsub("[.]","-",df.col.2$V1)
    df.col.2
    grid.col1 <- as.data.frame(col.table1())
    grid.col1
    df2 <- grid.col1[,names(grid.col1) %in% df.col.2$V1]
    df2[,order(names(df2))]
    # grid.col3 <- as.matrix(grid.col2)
    # grid.col3
    
    
  })
  
  Circular_plot2 <- function () {
    
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    dat <- subset(dat, get(input$group_column)==input$group_selected2)
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    hierarchy <- hierarchy[,c(input$chain1,input$chain2)]
    hierarchy <- hierarchy %>%
      select(input$chain1, everything())
    hierarchy$cloneCount <- 1
    chain1 <- as.data.frame(ddply(hierarchy,names(hierarchy)[-c(2,3)],numcolwise(sum)))
    chain1 <- chain1[order(chain1$cloneCount, decreasing = T),]
    
    chain2 <- as.data.frame(ddply(hierarchy,names(hierarchy)[-c(1,3)],numcolwise(sum)))
    chain2 <- chain2[order(chain2$cloneCount, decreasing = T),]
    
    df.col1 <- as.data.frame(chain1[,1])
    names(df.col1) <- "V1"
    df.col.j <- as.data.frame(chain2[,1])
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    
    grid.col1 <- as.data.frame(col.table1())
    grid.col2 <- grid.col1[,names(grid.col1) %in% df.col.2$V1]
    grid.col2 <-grid.col2[,order(names(grid.col2))]
    grid.col3 <- as.matrix(grid.col2)
    names(grid.col3) <- names(grid.col2)
    
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    hierarchy <- hierarchy[,c(input$chain1,input$chain2)]
    hierarchy <- hierarchy %>%
      select(input$chain1, everything())
    hierarchy <- as.matrix(table(hierarchy[,1], hierarchy[,2]))
    
    par(mar = rep(0, 4), cex=input$CHORD.cex, family = input$font_type)
    
    if (input$circ_lab=="Label") {
      circos.clear()
      #par(new = TRUE) # <- magic
      circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
      chordDiagram(hierarchy, annotationTrack = "grid", grid.col = grid.col3,
                   order = df.col.2$V1,
                   transparency = input$chord.transparancy,
                   # transparency = 0.5,
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(hierarchy))))))
      # we go back to the first track and customize sector labels
      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.par(track.margin=c(0,0)) 
        xlim = get.cell.meta.data("xlim")
        sector.index = get.cell.meta.data("sector.index")
        #text direction (dd) and adjusmtents (aa)
        theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
        dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
        aa = c(1, 0.5)
        if(theta < 90 || theta > 270)  aa =c(0, 0.5)
        circos.text(x = mean(xlim), y = 0.1, labels = sector.index, facing = dd, adj = aa)
        
      }, bg.border = NA)
      
    }
    # 'colour selected clone/s (label)' || 'colour selected clone/s (no label)''
    else if (input$circ_lab =="colour selected clone/s (label)") {
      lwd_mat = hierarchy
      
      # line thickness
      lwd_mat = hierarchy
      lwd_mat[lwd_mat>0] <- "x"
      lwd_mat[rownames(lwd_mat) %in% input$string.data.circ.order & lwd_mat=="x"] <- input$thickness.chord.line
      lwd_mat[!rownames(lwd_mat) %in% input$string.data.circ.order & lwd_mat=="x"] <- 0
      lwd_mat[lwd_mat==0] <- 1
      
      
      # boarder colour
      border_mat <- hierarchy
      border_mat[border_mat>0] <- 1
      border_mat[rownames(border_mat) %in% input$string.data.circ.order & border_mat==1] <- input$colour.chord.line
      border_mat[!rownames(border_mat) %in% input$string.data.circ.order & border_mat==1] <- 0
      border_mat[border_mat==0] <- NA
      border_mat
      
      # line type 
      lty_mat = hierarchy
      lty_mat[lty_mat>0] <- input$line.chord.type
      
      # transparancy 
      alpha_mat <- hierarchy
      alpha_mat[alpha_mat>0] <- 1
      alpha_mat[rownames(alpha_mat) %in% input$string.data.circ.order & alpha_mat==1] <- input$selected.chord.transparacy
      alpha_mat[!rownames(alpha_mat) %in% input$string.data.circ.order & alpha_mat==1] <- input$unselected.chord.transparacy
      alpha_mat

      
      
      
      circos.clear()
      #par(new = TRUE) # <- magic
      circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
      chordDiagram(hierarchy, annotationTrack = "grid", grid.col = grid.col3,
                   order = df.col.2$V1,
                   link.lty = lty_mat,
                   link.lwd = lwd_mat,
                   link.border = border_mat,
                   # transparency = alpha_mat,
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(hierarchy))))))
      # we go back to the first track and customize sector labels
      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.par(track.margin=c(0,0)) 
        xlim = get.cell.meta.data("xlim")
        sector.index = get.cell.meta.data("sector.index")
        #text direction (dd) and adjusmtents (aa)
        theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
        dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
        aa = c(1, 0.5)
        if(theta < 90 || theta > 270)  aa =c(0, 0.5)
        circos.text(x = mean(xlim), y = 0.1, labels = sector.index, facing = dd, adj = aa)
      }, bg.border = NA)
      
    }
    
    else if (input$circ_lab =="colour selected clone/s (no label)") {
      lwd_mat = hierarchy
      
      # line thickness
      lwd_mat = hierarchy
      lwd_mat[lwd_mat>0] <- "x"
      lwd_mat[rownames(lwd_mat) %in% input$string.data.circ.order & lwd_mat=="x"] <- input$thickness.chord.line
      lwd_mat[!rownames(lwd_mat) %in% input$string.data.circ.order & lwd_mat=="x"] <- 0
      lwd_mat[lwd_mat==0] <- 1
      
      
      # boarder colour
      border_mat <- hierarchy
      border_mat[border_mat>0] <- 1
      border_mat[rownames(border_mat) %in% input$string.data.circ.order & border_mat==1] <- input$colour.chord.line
      border_mat[!rownames(border_mat) %in% input$string.data.circ.order & border_mat==1] <- 0
      border_mat[border_mat==0] <- NA
      border_mat
      
      # line type 
      lty_mat = hierarchy
      lty_mat[lty_mat>0] <- input$line.chord.type
      
      # transparancy 
      alpha_mat <- hierarchy
      alpha_mat[alpha_mat>0] <- 1
      alpha_mat[rownames(alpha_mat) %in% input$string.data.circ.order & alpha_mat==1] <- input$selected.chord.transparacy
      alpha_mat[!rownames(alpha_mat) %in% input$string.data.circ.order & alpha_mat==1] <- input$unselected.chord.transparacy
      alpha_mat
      
      circos.clear()
      #par(new = TRUE) # <- magic
      circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
      chordDiagram(hierarchy, annotationTrack = "grid", grid.col = grid.col3,
                   order = df.col.2$V1,
                   link.lty = lty_mat,
                   link.lwd = lwd_mat,
                   link.border = border_mat,
                   # transparency = alpha_mat,
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(hierarchy))))))
      # we go back to the first track and customize sector labels
      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.par(track.margin=c(0,0))
        xlim = get.cell.meta.data("xlim")
        sector.index = get.cell.meta.data("sector.index")
        #text direction (dd) and adjusmtents (aa)
        theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
        dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
        aa = c(1, 0.5)
        if(theta < 90 || theta > 270)  aa =c(0, 0.5)
        
      }, bg.border = NA)
      
      
    }
    
    else {
      
      
      chordDiagram(hierarchy, annotationTrack = "grid", grid.col = grid.col3,
                   order = df.col.2$V1,
                   
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(hierarchy))))))
      
    }
    
  }
  output$Circular <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    Circular_plot2()
  })
  output$downloadPlot_circ <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_circular_plot_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_circ,height=input$height_circ, onefile = FALSE) # open the pdf device
      print(Circular_plot2())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_circ <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_circular_plot_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_circ, height = input$height_png_circ, res = input$resolution_PNG_circ)
      print(Circular_plot2())
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  # CDR3 distribution =====
  observe({
    updateSelectInput(
      session,
      "aa.or.nt",
      choices=names(input.data_old2()),
      selected = "JUNCTION..AA._A") 
  }) # amino acid or nucleotides column
  observe({
    updateSelectInput( 
      session,
      "selected_group_len",
      choices=select_group()) }) # group
  observe({
    updateSelectInput(
      session,
      "chain.hist.col",
      choices=names(input.data_old2()),
      selected = "AVJ")
  })
  
  cols.hist <- reactive({
    df <- input.data_old2(); 
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    # 
    # df.names <-  df[ , -which(names(df) %in% c("cloneCount","clone"))]
    df1 <- df
    df1$len1 <- nchar(df1[,grep(input$aa.or.nt,names(df1))])
    
    df1$chain <- df1[,names(df1) %in% input$chain.hist.col]
    df1 <- df1[order(df1$chain, decreasing = F),]
    df1$chain <- factor(df1$chain,levels = unique(df1$chain))
    
    num <- as.data.frame(unique(df1$chain))
    
    col.gg <- gg_fill_hue(dim(num)[1])
    
    unique.col <- as.data.frame(unique(df$chain))
    
    palette_rainbow <- rev(rainbow(dim(num)[1]))
    
    if (input$hist_colour.choise == "rainbow") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.hist", i, sep="_"), paste(num[i,]), palette_rainbow[i])        
      }) }
    else if (input$hist_colour.choise == "default") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.hist", i, sep="_"), paste(num[i,]), col.gg[i])        
      })
    }
    else if (input$hist_colour.choise == "random") {
      palette1 <- distinctColorPalette(dim(num)[1])
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.hist", i, sep="_"), paste(num[i,]), palette1[i])        
      })
      
    }
    else {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.hist", i, sep="_"), paste(num[i,]), input$one.colour.default)        
      })
      
      
    }
    
  })
  output$myPanel.hist <- renderUI({cols.hist()})
  colors.hist <- reactive({
    df <- input.data_old2(); 
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df.names <-  df[ , -which(names(df) %in% c("cloneCount","well"))]
    df1 <- ddply(df,names(df.names) ,numcolwise(sum))
    df1$len1 <- nchar(df1[,grep(input$aa.or.nt,names(df1))])
    
    df1$chain <- df1[,names(df1) %in% input$chain.hist.col]
    df1 <- df1[order(df1$chain,decreasing = F),]
    df1$chain <- factor(df1$chain,levels = unique(df1$chain))
    
    num <- as.data.frame(unique(df1$chain))
    
    lapply(1:dim(num)[1], function(i) {
      input[[paste("col.hist", i, sep="_")]]
    })
  })
  
  cols.hist2 <- reactive({
    df <- input.data_old2(); 
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df1 <- df
    
    num <- as.data.frame(unique(df1[names(df1) %in% input$group_column]))
    
    col.gg <- gg_fill_hue(dim(num)[1])
    unique.col <- as.data.frame(unique(df1[names(df1) %in% input$group_column]))
    
    palette_rainbow <- rev(rainbow(dim(num)[1]))
    
    if (input$hist_colour.choise == "rainbow") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.hist2", i, sep="_"), paste(num[i,]), palette_rainbow[i])        
      }) }
    else if (input$hist_colour.choise == "default") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.hist2", i, sep="_"), paste(num[i,]), col.gg[i])        
      })
    }
    else if (input$hist_colour.choise == "random") {
      palette1 <- distinctColorPalette(dim(num)[1])
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.hist2", i, sep="_"), paste(num[i,]), palette1[i])        
      })
      
    }
    else {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.hist2", i, sep="_"), paste(num[i,]), input$one.colour.default)        
      })
      
      
    }
    
  })
  output$myPanel.hist2 <- renderUI({cols.hist2()})
  colors.hist2 <- reactive({
    df <- input.data_old2(); 
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df1 <- df
    
    num <- as.data.frame(unique(df1[names(df1) %in% input$group_column]))
    lapply(1:dim(num)[1], function(i) {
      input[[paste("col.hist2", i, sep="_")]]
    })
  })
  
  output$hist.table <- DT::renderDataTable( {
    df <- input.data_old2(); 
    df <- as.data.frame(df)
    df <- df[names(df) %in% c("cloneCount",input$group_column,input$aa.or.nt,input$chain.hist.col)]
    df.names <-  df[ , -which(names(df) %in% c("cloneCount"))]
    df1 <- ddply(df,names(df.names) ,numcolwise(sum))
    df1
    df1$len1 <- nchar(df1[,grep(input$aa.or.nt,names(df1))])
    df1$chain <- df1[,names(df1) %in% input$chain.hist.col]
    df1
    
    datatable(df1, extensions = "Buttons",filter = "top", options = list(searching = TRUE,
                                                          ordering = TRUE,
                                                          buttons = c('copy','csv', 'excel'),
                                                          dom = 'Bfrtip',
                                                          pageLength=5,
                                                          lengthMenu=c(2,5,10,20,50,100),
                                                          scrollX = TRUE
    ))
  }, server = FALSE) 
  
  hist.col.table <- function () {
    df <- input.data_old2();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df1 <- df
    df1$len1 <- nchar(df1[,grep(input$aa.or.nt,names(df1))])
    
    df1$chain <- df1[,names(df1) %in% input$chain.hist.col]
    df1 <- df1[order(df1$chain),]
    df1$chain <- factor(df1$chain,levels = unique(df1$chain))
    
    df.col.2 <- as.data.frame(unique(df1$chain))
    names(df.col.2) <- "V1"
    col2 <- unlist(colors.hist())
    as.data.frame(col2)
    df.col.2$col <- col2
    df.col.2
  }
  
  
  # CHAIN LENGTH HISTOGRAM ----
  Chain1_length <- function () {
    df <- input.data_old2(); 
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    
    if (input$graph_type == "histogram" & input$type.tree == "raw data") {
      df <- as.data.frame(df)
      df <- df[names(df) %in% c("cloneCount",input$group_column,input$aa.or.nt,input$chain.hist.col)]
      df.names <-  df[ , -which(names(df) %in% c("cloneCount"))]
      df1 <- ddply(df,names(df.names) ,numcolwise(sum))
      
      df1$len1 <- nchar(df1[,grep(input$aa.or.nt,names(df1))])
      df1$chain <- df1[,names(df1) %in% input$chain.hist.col]
      df1 <- df1[order(df1$chain, decreasing = F),]
      df1$chain <- factor(df1$chain,levels = unique(df1$chain))
      df.col.2 <- as.data.frame(hist.col.table())
      names(df.col.2) <- c("V1","col")
      df2 <- merge(df.col.2,df1,by.x="V1",by.y=input$chain.hist.col,sort = F)
      df1 <- subset(df2, get(input$group_column)==input$selected_group_len)
      
      df.col.hist <- df.col.2[df.col.2$V1 %in% unique(df1$chain),]
      df1$unique <- 1
      max.1 <- ddply(df1, c(input$group_column,"len1"),numcolwise(sum))
      
      if (input$sep_hist_group == "yes") {
        df <- input.data2(); 
        df <- as.data.frame(df)
        df <- df[names(df) %in% c("cloneCount",input$group_column,input$aa.or.nt,input$chain.hist.col)]
        df.names <-  df[ , -which(names(df) %in% c("cloneCount"))]
        df1 <- ddply(df,names(df.names) ,numcolwise(sum))
        
        df1$len1 <- nchar(df1[,grep(input$aa.or.nt,names(df1))])
        df1$chain <- df1[,names(df1) %in% input$chain.hist.col]
        df1 <- df1[order(df1$chain, decreasing = F),]
        df1$chain <- factor(df1$chain,levels = unique(df1$chain))
        df.col.2 <- as.data.frame(hist.col.table())
        names(df.col.2) <- c("V1","col")
        df2 <- merge(df.col.2,df1,by.x="V1",by.y=input$chain.hist.col,sort = F)
        df1 <- df2
        
        df.col.hist <- df.col.2[df.col.2$V1 %in% unique(df1$chain),]
        df1$unique <- 1
        max.1 <- ddply(df1, c(input$group_column,"len1"),numcolwise(sum))
        
        max.2 <- max.1
        max.hist <- max(max.2$unique)+1
        
        vals4$bar.len <- ggplot(df1,aes(x=len1,fill = chain)) +
          geom_bar() + 
          scale_fill_manual(values = df.col.hist$col) +
          theme_bw()  +
          theme(legend.title = element_blank(),
                legend.position = input$hist.density.legend) +
          labs(y="Count",
               x="CDR3 length distribution",
               title="") +
          theme(
            axis.title.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer2),
            axis.text.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer),
            axis.text.x = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer,angle=0),
            axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$hist.text.sizer2),
            legend.text = element_text(colour="black", size=input$legend.text.hist,family=input$font_type) 
          ) +
          scale_y_continuous(limits = c(0, max.hist), breaks = seq(0,max.hist,by = input$ybreaks), expand = c(0, 0))+
          scale_x_continuous(limits = c(input$xlow, input$xhigh), breaks = seq(input$xlow, input$xhigh, by = input$xbreaks),expand = c(0, 0)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          guides(fill=guide_legend(ncol=input$col.num.CDR3len)) +
          facet_grid(get(input$group_column)~.)
        vals4$bar.len
        
        
      }
      
      else {
        max.2 <- subset(max.1, get(input$group_column)==input$selected_group_len)
        max.hist <- max(max.2$unique)+1
        
        vals4$bar.len <- ggplot(df1,aes(x=len1,fill = chain)) +
          geom_bar() + 
          scale_fill_manual(values = df.col.hist$col) +
          theme_bw()  +
          theme(legend.title = element_blank(),
                legend.position = input$hist.density.legend) +
          labs(y="Count",
               x="CDR3 length distribution",
               title="") +
          theme(
            axis.title.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer2),
            axis.text.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer),
            axis.text.x = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer,angle=0),
            axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$hist.text.sizer2),
            legend.text = element_text(colour="black", size=input$legend.text.hist,family=input$font_type) 
          ) +
          scale_y_continuous(limits = c(0, max.hist), breaks = seq(0,max.hist,by = input$ybreaks), expand = c(0, 0))+
          scale_x_continuous(limits = c(input$xlow, input$xhigh), breaks = seq(input$xlow, input$xhigh, by = input$xbreaks),expand = c(0, 0)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          guides(fill=guide_legend(ncol=input$col.num.CDR3len)) 
        vals4$bar.len
        
      }
      
    }
    else if (input$graph_type == "density" & input$type.tree == "raw data") {
      df <- as.data.frame(df)
      head(df)
      df <- df[names(df) %in% c("cloneCount",input$group_column,input$aa.or.nt,input$chain.hist.col)]
      df.names <-  df[ , -which(names(df) %in% c("cloneCount"))]
      df1 <- ddply(df,names(df.names) ,numcolwise(sum))
      names(df1)
      
      df1$len1 <- nchar(df1[, which(names(df1) %in% c(input$aa.or.nt))])
      
      df1$unique <- 1
      max.1 <- ddply(df1, c(input$group_column,"len1"),numcolwise(sum))
      max.2 <- subset(max.1, get(input$group_column)==input$selected_group_len)
      max.2$feq <- max.2$unique/sum(max.2$unique)
      max.hist <- max(max.2$feq)+0.05
      
      col2 <- unlist(colors.hist2())
      num <- as.data.frame(unique(df1[names(df1) %in% input$group_column]))
      num$col2 <- col2
      
     
        
      
      if (input$facet_by_group =="yes") {
        vals4$bar.len <-  vals4$bar.len <- ggplot(df1,aes(x=len1,colour = get(input$group_column),fill = get(input$group_column))) +
          geom_density(alpha = input$alpha.density) +
          scale_alpha(guide = 'none') + 
          scale_color_manual(values = num$col2) +
          scale_fill_manual(values = num$col2) +
          theme_bw()  +
          theme(legend.title = element_blank(),
                legend.position = input$hist.density.legend) +
          labs(y="Frequency",
               x="CDR3 length distribution",
               title="") +
          theme(
            axis.title.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer2),
            axis.text.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer),
            axis.text.x = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer,angle=0),
            axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type,size = input$hist.text.sizer2),
            legend.text = element_text(colour="black", size=input$legend.text.hist,family=input$font_type) 
          ) +
          scale_y_continuous(expand = c(0,0), limits=c(0,max.hist)) +
          scale_x_continuous(limits = c(input$xlow, input$xhigh), breaks = seq(input$xlow, input$xhigh, by = input$xbreaks),expand = c(0,0)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          guides(fill=guide_legend(ncol=input$col.num.CDR3len)) +
          facet_grid(get(input$group_column)~.)
        
      }
      
      else {
        vals4$bar.len <- ggplot(df1,aes(x=len1,colour = get(input$group_column),fill = get(input$group_column))) +
          geom_density(alpha = input$alpha.density) +
          scale_alpha(guide = 'none') + 
          scale_color_manual(values = num$col2) +
          scale_fill_manual(values = num$col2) +
          theme_bw()  +
          theme(legend.title = element_blank(),
                legend.position = input$hist.density.legend) +
          labs(y="Frequency",
               x="CDR3 length distribution",
               title="") +
          theme(
            axis.title.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer2),
            axis.text.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer),
            axis.text.x = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer,angle=0),
            axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type,size = input$hist.text.sizer2),
            legend.text = element_text(colour="black", size=input$legend.text.hist,family=input$font_type) 
          ) +
          scale_y_continuous(expand = c(0,0), limits=c(0,max.hist)) +
          scale_x_continuous(limits = c(input$xlow, input$xhigh), breaks = seq(input$xlow, input$xhigh, by = input$xbreaks),expand = c(0,0)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
          guides(fill=guide_legend(ncol=input$col.num.CDR3len)) 
      }
      
      vals4$bar.len
    }
    else if (input$graph_type == "histogram" & input$type.tree == "Summarised data") {
      
      df <- as.data.frame(df)
      df1 <- df
      
      df1$len1 <- nchar(df1[,grep(input$aa.or.nt,names(df1))])
      df1$chain <- df1[,names(df1) %in% input$chain.hist.col]
      df1 <- df1[order(df1$chain, decreasing = F),]
      df1$chain <- factor(df1$chain,levels = unique(df1$chain))
      df.col.2 <- as.data.frame(hist.col.table())
      names(df.col.2) <- c("V1","col")
      df2 <- merge(df.col.2,df1,by.x="V1",by.y=input$chain.hist.col,sort = F)
      df1 <- subset(df2, get(input$group_column)==input$selected_group_len)
      
      df.col.hist <- df.col.2[df.col.2$V1 %in% unique(df1$chain),]
      
      df1$unique <- 1
      max.1 <- ddply(df1, c(input$group_column,"len1"),numcolwise(sum))
      max.2 <- subset(max.1, get(input$group_column)==input$selected_group_len)
      max.hist <- max(max.2$unique)+1
      
      vals4$bar.len <- ggplot(df1,aes(x=len1,fill = chain)) +
        geom_bar() + 
        theme_bw()  +
        theme(legend.title = element_blank(),
              legend.position = input$hist.density.legend) +
        scale_fill_manual(values = df.col.hist$col) +
        labs(y="Count",
             x="CDR3 length distribution",
             title="") +
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer,angle=0),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type,size = input$hist.text.sizer2),
          legend.text = element_text(colour="black", size=input$legend.text.hist,family=input$font_type) 
        ) +
        scale_y_continuous(limits = c(0, max.hist), breaks = seq(0,max.hist,by = input$ybreaks), expand = c(0, 0))+
        scale_x_continuous(limits = c(input$xlow, input$xhigh), breaks = seq(input$xlow, input$xhigh, by = input$xbreaks),expand = c(0, 0)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        guides(fill=guide_legend(ncol=input$col.num.CDR3len)) 
      vals4$bar.len
      
    }
    else {
      df1 <- df
      df1$len1 <- nchar(df1[, which(names(df1) %in% c(input$aa.or.nt))])
      
      df1$unique <- 1
      max.1 <- ddply(df1, c(input$group_column,"len1"),numcolwise(sum))
      max.2 <- subset(max.1, get(input$group_column)==input$selected_group_len)
      max.2$feq <- max.2$unique/sum(max.2$unique)
      max.hist <- max(max.2$feq)+0.05
      
      col2 <- unlist(colors.hist2())
      num <- as.data.frame(unique(df1[names(df1) %in% input$group_column]))
      num$col2 <- col2
      
      vals4$bar.len <- ggplot(df1,aes(x=len1,colour = get(input$group_column),fill = get(input$group_column))) +
        geom_density(alpha = 0.25) +
        scale_alpha(guide = 'none') + 
        scale_color_manual(values = num$col2) +
        scale_fill_manual(values = num$col2) +
        theme_bw()  +
        theme(legend.title = element_blank(),
              legend.position = input$hist.density.legend) +
        labs(y="Frequency",
             x="CDR3 length distribution",
             title="") +
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$hist.text.sizer,angle=0),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type,size = input$hist.text.sizer2),
          legend.text = element_text(colour="black", size=input$legend.text.hist,family=input$font_type) 
        ) +
        scale_y_continuous(expand = c(0,0), limits=c(0,max.hist)) +
        scale_x_continuous(limits = c(input$xlow, input$xhigh), breaks = seq(input$xlow, input$xhigh, by = input$xbreaks),expand = c(0,0)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        guides(fill=guide_legend(ncol=input$col.num.CDR3len)) 
      vals4$bar.len
    }
    
  }
  
  output$Chain1_length <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    Chain1_length()
  })
  
  output$downloadPlot_length <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_length_plot_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_length,height=input$height_length, onefile = FALSE) # open the pdf device
      print(Chain1_length())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_length <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_length_plot_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_length, height = input$height_png_length, res = input$resolution_PNG_length)
      print(Chain1_length())
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  # download summary table
  table.len.download <- reactive( {
    df <- input.data_old2()
    df <- as.data.frame(df)

    if (input$type.tree == "raw data") {
      df1 <-  df[names(df) %in% c("cloneCount",
                                  "Indiv",
                                  "group",
                                  "GVJ",
                                  "DVDJ",
                                  "AVJ",
                                  "BVDJ",
                                  names(df[grep("JUNCTION",names(df))]),
                                  names(df[grep("IMGT",names(df))])
      )]   
      
      df2 <- ddply(df1,names(df1[-c(1)]),numcolwise(sum))
      df2
      
      df3 <-  df2[names(df2) %in% c(names(df2[grep("JUNCTION",names(df2))]),
                                    names(df2[grep("IMGT",names(df2))])
                 )]
      
      for (i in 1:dim(df3)[1]) {
        df3[i,] <- nchar(df3[i,])
        df3
      }
      names(df3) <- paste(names(df3),"length", sep="_")
      df4 <- cbind(df2,df3)
      df4
    }
    
    else {
      df1 <- df
      for (i in 1:dim(df1)[1]) {
        df1[i,] <- nchar(df1[i,])
        df1
      }
      names(df1) <- paste(names(df1),"length", sep="_")
      df4 <- cbind(df,df1)
      df4
    }
    
    
  })
  
  output$table_length <- downloadHandler(
    filename = function(){
      paste("TCR_Explore_length_plot_",".csv", sep = "")
    },
    content = function(file){
      write.csv(table.len.download(),file, row.names = FALSE)
    }
  )
  
  
  # bar plot -----
  observe({
    updateSelectInput(
      session,
      "variable_chain",
      choices=names(input.data_old2()),
      selected = "AVJ")
    
  })
  observe({
    updateSelectInput( 
      session,
      "selected_group_chain",
      choices=select_group()) }) # group
  
  Chain_usage <- function () {
    df <- input.data2(); 
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df <- subset(df, get(input$group_column)==input$selected_group_chain)
    df2 <- as.data.frame(ddply(df,c(input$variable_chain),numcolwise(sum)))[1:2]
    names(df2) <- c("chain","cloneCount")
    
    df2 <- df2[order(df2$cloneCount),]
    df2$chain <- factor(df2$chain, levels = unique(df2$chain),labels = df2$chain)
    
    if (input$graph_bar_type == "count") {
      
      
      vals5$bar.usage <- ggplot(df2,aes(x=chain,y=cloneCount,fill=input$colour_bar.usage)) +
        geom_bar(stat="identity", position = "dodge") +
        theme_bw()  +
        theme(legend.title = element_blank(),
              legend.position = "none") +
        scale_fill_manual(values=input$colour_bar.usage) +
        labs(y="count",
             x="",
             title="") +
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$bar.numeric.size),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$bar.numeric.size),
          axis.text.x = element_text(colour="black",family=input$font_type,angle=0,size = input$bar.numeric.size),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type,size = input$bar.numeric.size),
          legend.text = element_text(colour="black",family=input$font_type)
        ) +
        coord_flip()
      vals5$bar.usage
      
      
      
    }
    
    else {
      df2$percent <- round(df2$cloneCount/sum(df2$cloneCount)*100,2)
      vals5$bar.usage <- ggplot(df2,aes(x=chain,y=percent)) +
        geom_bar(stat="identity", position = "dodge",fill=input$colour_bar.usage) +
        theme_bw()  +
        theme(legend.title = element_blank(),
              legend.position = "bottom") +
        #scale_fill_manual(values=c("red","blue","darkgreen")) +
        labs(y="Percentage",
             x="",
             title="") +
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = 12),
          axis.text.y = element_text(colour="black",family=input$font_type,size = 12),
          axis.text.x = element_text(colour="black",family=input$font_type,angle=0,size = 12),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type,size = 12),
          legend.text = element_text(colour="black",family=input$font_type)
        ) +
        coord_flip()
      vals5$bar.usage
      
    }
    
  }
  
  vals30 <- reactiveValues(bar.usage2=NULL)
  
  Chain2_usage <- function () {
    df <- input.data2(); 
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df <- subset(df, get(input$group_column)==input$selected_group_chain)
    
    
    df2 <- df[names(df) %in% c("cloneCount",input$variable_chain)]
    df3 <- as.data.frame(ddply(df2,names(df2)[2],numcolwise(sum)))
    df3$count2 <- 1
    df3$percent <- df3$cloneCount/sum(df3$cloneCount)
    df4 <- as.data.frame(ddply(df3,names(df3)[2],numcolwise(sum)))
    df5 <- df4 %>% group_by(percent) %>% mutate(csum = cumsum(percent))
    
    vals30$bar.usage2 <- ggplot()+
      geom_bar(aes(x=df3$cloneCount,y=df3$percent),stat="identity",fill=input$colour_bar.usage)+
      geom_line(aes(x=df5$cloneCount,y=df5$csum))+
      geom_point(aes(x=df5$cloneCount,y=df5$csum)) +
      geom_text(aes(x=df5$cloneCount,y=df5$percent,label=df5$count2),vjust=input$numeric.adjust, family=input$font_type,colour=input$colour.numeric.bar, size=input$label.size)+
      xlab("Distinct CDR3")+
      ylab("Frequency in repertoire")+
      theme_bw()+
      theme(
        axis.title.y = element_text(colour="black",family=input$font_type,size = input$label.size.axis),
        axis.text.y = element_text(colour="black",family=input$font_type,size = input$label.size.axis),
        axis.text.x = element_text(colour="black",family=input$font_type,angle=0,size = input$label.size.axis),
        axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type,size = input$label.size.axis),
        legend.text = element_text(colour="black",family=input$font_type)
      ) 
    
    vals30$bar.usage2
    
    
  }
  
  vals31 <- reactiveValues(bar.usage3=NULL)
  
  observe({
    updateSelectInput(
      session,
      "string.data2",
      choices=select_group(),
      selected = c("IFN","CD8")) 
  }) 
  
  cols_stacked_bar <- reactive({
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    df <- as.data.frame(ddply(dat,(c(input$group_column,input$variable_chain)),numcolwise(sum)))
    names(df) <- c("group","chain","cloneCount")
    
    df <-df[order(df$chain),]
    
    num <- unique(df$chain)
    col.gg <- gg_fill_hue(length(num))
    
    palette_rainbow <- rev(rainbow(length(t(num))))
    
    if (input$bar.stacked_colour.choise == "default") {
      lapply(1:length(num), function(i) {
        colourInput(paste("cols_stacked_bar", i, sep="_"), paste(num[i]), col.gg[i])        
      })}
    
    else  if (input$bar.stacked_colour.choise == "rainbow") {
      lapply(1:length(num), function(i) {
        colourInput(paste("cols_stacked_bar", i, sep="_"), paste(num[i]), palette_rainbow[i])        
      }) }
    
    else if (input$bar.stacked_colour.choise == "random") {
      
      lapply(1:length(num), function(i) {
        palette1 <- distinctColorPalette(length(num))
        colourInput(paste("cols_stacked_bar", i, sep="_"), paste(num[i]), palette1[i])        
      }) }
    else  {
      lapply(1:length(num), function(i) {
        colourInput(paste("cols_stacked_bar", i, sep="_"), paste(num[i]), input$one.colour.default)        
      }) }
    
  })
  
  output$myPanel_cols_stacked_bar <- renderUI({cols_stacked_bar()})
  
  colors_bar.stacked <- reactive({
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    df <- as.data.frame(ddply(dat,(c(input$group_column,input$variable_chain)),numcolwise(sum)))
    names(df) <- c("group","chain","cloneCount")
    
    df <-df[order(df$chain),]
    
    num <- unique(df$chain)
    lapply(1:length(num), function(i) {
      input[[paste("cols_stacked_bar", i, sep="_")]]
    })
  })
  
  Chain3_usage <- function () {
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    df <- as.data.frame(ddply(dat,(c(input$group_column,input$variable_chain)),numcolwise(sum)))
    names(df) <- c("group","chain","cloneCount")
    
    df <-df[order(df$chain),]
    df <-df[df$group %in% input$string.data2,]
    df$group <- factor(df$group, levels = input$string.data2)
    cols <- unlist(colors_bar.stacked())
    palette <- cols
    
    if (input$lines.bar.graph == 'yes') {
      vals31$bar.usage3 <-  ggplot(df, aes(fill=chain, y=cloneCount, x=group,colour="black")) +
        geom_bar(position="fill", stat="identity") +
        xlab("")+
        ylab("Frequency")+
        guides(fill=guide_legend(ncol=input$stacked.no.legend))+
        theme_bw()+
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$label.size.axis2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$label.size.axis2),
          axis.text.x = element_text(colour="black",family=input$font_type,angle=input$bar.stack.angle,size = input$label.size.axis2, hjust=input$hight.bar.stack.adj),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type,size = input$label.size.axis2),
          legend.text = element_text(colour="black",family=input$font_type,size = input$stacked.legend.size),
          legend.title = element_blank()
        )+ 
        scale_fill_manual(values=palette) +
        scale_color_manual(values="black", guide = "none")
      vals31$bar.usage3
      
      
    }
    
    else {
      vals31$bar.usage3 <-  ggplot(df, aes(fill=chain, y=cloneCount, x=group)) +
        geom_bar(position="fill", stat="identity") +
        xlab("")+
        ylab("Frequency")+
        guides(fill=guide_legend(ncol=input$stacked.no.legend))+
        theme_bw()+
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$label.size.axis2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$label.size.axis2),
          axis.text.x = element_text(colour="black",family=input$font_type,angle=input$bar.stack.angle,size = input$label.size.axis2, hjust=input$hight.bar.stack.adj),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type,size = input$label.size.axis2),
          legend.text = element_text(colour="black",family=input$font_type,size = input$stacked.legend.size),
          legend.title = element_blank(),
          legend.position = input$stacked.legend,
          
        ) +
        scale_fill_manual(values=palette) 
      
      vals31$bar.usage3
      
    }
    
  }
  
  output$Chain1_usage <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    if (input$stat=="chains") {
      print(Chain_usage())
    }
    
    
    else if (input$stat=="frequency") {
      print(Chain2_usage())
    }
    
    else {
      print(Chain3_usage())
      
    }
    
  }) 
  
  # downloading bar plots =====
  output$downloadPlot_chain.usage <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("bar.plot_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_chain.usage,height=input$height_chain.usage, onefile = FALSE) # open the pdf device
      if (input$stat=="chains") {
        print(Chain_usage())
      }
      
      
      else if (input$stat=="frequency") {
        print(Chain2_usage())
      }
      
      else {
        print(Chain3_usage())
        
      }
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_chain.usage <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("bar.plot_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_chain.usage, height = input$height_png_chain.usage, res = input$resolution_PNG_chain.usage)
      if (input$stat=="chains") {
        print(Chain_usage())
      }
      
      
      else if (input$stat=="frequency") {
        print(Chain2_usage())
      }
      
      else {
        print(Chain3_usage())
        
      }
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  # skewness of data (to be possibly added as stat) ------
  skewness.data <- function () {
    df <- input.data2(); 
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df <- subset(df, get(input$group_column)==input$selected_group_chain)
    df2 <- as.data.frame(ddply(df,c(input$variable_chain),numcolwise(sum)))[1:2]
    names(df2) <- c("chain","cloneCount")
    
    df2 <- df2[order(df2$cloneCount),]
    df2$chain <- factor(df2$chain, levels = unique(df2$chain),labels = df2$chain)
    
    
  }
  # motif amino acid-----
  observe({
    updateSelectInput(
      session,
      "aa.or.nt2",
      choices=names(input.data_old2()),
      selected = "JUNCTION..AA._A")})
  observe({
    updateSelectInput(
      session,
      "group_selected_motif",
      choices=select_group(),
      selected = "IFN") })
  observe({
    updateSelectInput(
      session,
      "group_selected_motif2",
      choices=select_group(),
      selected = "CD8") })
  
  output$Motif <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE), {
    df <- input.data_old2();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt2)),numcolwise(sum)))
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt2])
    df_subset <- subset(df_unique,df_unique$len1==input$len)
    df_subset <- subset(df_subset,get(input$group_column)==input$group_selected_motif)
    
    motif <- as.data.frame(t(as.data.frame(strsplit(df_subset[,grep(input$aa.or.nt2,names(df_subset))], ""))))
    rownames(motif) <- 1:dim(motif)[1]
    
    
    datatable(motif, extensions = "Buttons", options = list(searching = TRUE,
                                                                        ordering = TRUE,
                                                                        buttons = c('copy','csv', 'excel'),
                                                                        dom = 'Bfrtip',
                                                                        pageLength=2, 
                                                                        lengthMenu=c(2,5,10,20,50,100), 
                                                                        scrollX = TRUE
    ))
    
    
  })
  
  output$length <- renderPrint({
    df <- input.data_old2();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df$len1 <- nchar(df[,grep(input$aa.or.nt2,names(df))])
    df <- df[order(df$len1),]
    df <- subset(df,get(input$group_column)==input$group_selected_motif)
    df_len <- unique(df$len1)
    cat(input$group_selected_motif,"dataset contains CDR3 lengths of:",  df_len)
  })
  
  output$length.table <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE), {
    df <- input.data_old2();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt2)),numcolwise(sum)))
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt2])
    df_unique$Unique_clones <- 1
    motif <- as.data.frame(ddply(df_unique,(c(input$group_column,"len1")),numcolwise(sum)))
    
    datatable(motif, extensions = "Buttons", filter = "top", options = list(searching = TRUE,
                                                            ordering = TRUE,
                                                            buttons = c('copy','csv', 'excel'),
                                                            dom = 'Bfrtip',
                                                            pageLength=2, 
                                                            lengthMenu=c(2,5,10,20,50,100), 
                                                            scrollX = TRUE
    ))
  })
  
  Motif_plot2 <- reactive({
    df <- input.data_old2();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt2)),numcolwise(sum)))
    
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt2])
    df_subset <- subset(df_unique,df_unique$len1==input$len)
    df_subset <- subset(df_subset,get(input$group_column)==input$group_selected_motif)
    
    motif <- as.data.frame(t(as.data.frame(strsplit(df_subset[,grep(input$aa.or.nt2,names(df_subset))], ""))))
    cbind(x=1,y=2,motif)
    
    
    motif_count <- aa.count.function(cbind(x=1,y=2,motif), input$len)
    motif_count<-pcm2pfm(motif_count)
    motif_count
    
    ggseqlogo(motif_count, seq_type='aa', method='p') + 
      ylab('bits')+ 
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      theme(
        axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
        axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
        legend.title  =element_blank(),
        legend.position = "right",
        legend.text = element_text(colour="black", size=input$font_legend_motif,family=input$font_type)
        )
  })
  
  Motif_compare_aa_group1 <- reactive({
    df <- input.data_old2()
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt2)),numcolwise(sum)))
    
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt2])
    df_subset <- subset(df_unique,df_unique$len1==input$len)
    df_subset <- subset(df_subset,get(input$group_column)==input$group_selected_motif)
    
    motif <- as.data.frame(t(as.data.frame(strsplit(df_subset[,grep(input$aa.or.nt2,names(df_subset))], ""))))
    motif_count <- aa.count.function(cbind(x=1,y=2,motif), input$len)
    motif_count1_aa<-pcm2pfm(motif_count)
    as.data.frame(motif_count1_aa)
  })
  
  Motif_compare_aa_group2 <- reactive({
    df <- input.data_old2()
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt2)),numcolwise(sum)))
    
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt2])
    df_subset <- subset(df_unique,df_unique$len1==input$len)
    df_subset <- subset(df_subset,get(input$group_column)==input$group_selected_motif2)
    
    motif <- as.data.frame(t(as.data.frame(strsplit(df_subset[,grep(input$aa.or.nt2,names(df_subset))], ""))))
    motif_count <- aa.count.function(cbind(x=1,y=2,motif), input$len)
    motif_count2_aa<-pcm2pfm(motif_count)
    as.data.frame(motif_count2_aa)
  })
  
  motif.compar.plot <- reactive({
    
    motif_count1_aa <- Motif_compare_aa_group1()
    motif_count2_aa <- Motif_compare_aa_group2()
    diffLogoObj = createDiffLogoObject(pwm1 = as.data.frame(motif_count1_aa), 
                                       pwm2 = as.data.frame(motif_count2_aa), 
                                       alphabet = ASN
                                       
    )
    
    mat <- (diffLogoObj$pwm1 - diffLogoObj$pwm2)
    names(mat) <- 1:dim(mat)[2]
    
    if (input$compar.lab.motif.aa.single == T) {
      
      vals33$geom_comp <- ggseqlogo(mat, method='custom', seq_type='aa') + 
        ylab('JS divergence') + 
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        annotate(geom="text",x=1,y=Inf,vjust=2,label=input$group_selected_motif,size=input$font_size_motif,face="plain",family=input$font_type)+
        annotate(geom="text",x=1,y=-Inf,vjust=-2,label=input$group_selected_motif2,size=input$font_size_motif,face="plain",family=input$font_type)+
        theme(
          axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
          axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
          axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
          axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
          legend.title  =element_blank(),
          legend.position = "right",
          legend.text = element_text(colour="black", size=input$font_legend_motif,family=input$font_type)) 
      
      vals33$geom_comp 
      
    }
    
    else {
      vals33$geom_comp <- ggseqlogo(mat, method='custom', seq_type='aa') + 
        ylab('JS divergence') + 
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        theme(
          axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
          axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
          axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
          axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
          legend.title  =element_blank(),
          legend.position = "right",
          legend.text = element_text(colour="black", size=input$font_legend_motif,family=input$font_type)) 
      
      vals33$geom_comp 
      
    }
    
    
    
  })
  
  output$Motif_plot <- renderPlot( {
    
    if (input$comarpison.aa.motif == "single.group1") {
      withProgress(message = 'Figure is being generated...',
                   detail = '', value = 0, {
                     test_fun()
                   })
      Motif_plot2()
      
    }
    
    else {
      
      motif.compar.plot()
      
    }
    
    
  })
  
  output$downloadPlot_motif <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_motif_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_motif,height=input$height_motif, onefile = FALSE) # open the pdf device
      if (input$comarpison.aa.motif == "single.group1") {
        withProgress(message = 'Figure is being generated...',
                     detail = '', value = 0, {
                       test_fun()
                     })
        plot(Motif_plot2())
        
      }
      
      else {
        plot(motif.compar.plot())
        
      }
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_motif <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_motif_", gsub("/", "-", x), ".png", sep = "")
    },    
    content = function(file) {
      png(file, width = input$width_png_motif, 
          height = input$height_png_motif, 
          res = input$resolution_PNG_motif)
      if (input$comarpison.aa.motif == "single.group1") {
        motif <- Motif_plot2()
        withProgress(message = 'Figure is being generated...',
                     detail = '', value = 0, {
                       test_fun()
                     })
        plot(motif)
        
      }
      
      else {
        plot(motif.compar.plot())
        
      }
      dev.off()},   contentType = "application/png" # MIME type of the image
  )
  
  # motif NT ------
  observe({
    updateSelectInput(
      session,
      "aa.or.nt3",
      choices=names(input.data_old2()),
      selected = "JUNCTION_A")
  })
  
  observe({
    updateSelectInput(
      session,
      "group_selected",
      choices=select_group())
    
  })
  
  output$Motif_nt <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5), {
    df <- input.data_old2();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt3)),numcolwise(sum)))
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt3])
    df_subset <- subset(df_unique,df_unique$len1==input$len_nt)
    df_subset <- subset(df_subset,get(input$group_column)==input$group_selected)
    motif <- as.data.frame(t(as.data.frame(strsplit(df_subset[,grep(input$aa.or.nt3,names(df_subset))], ""))))
    rownames(motif) <- 1:dim(motif)[1]
    motif
  })
  output$length_nt <- renderPrint( {
    df <- input.data_old2();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df$len1 <- nchar(df[,grep(input$aa.or.nt3,names(df))])
    df <- subset(df,get(input$group_column)==input$group_selected)
    df <- df[order(df$len1),]
    df_len <- unique(df$len1)
    cat("The dataset contains nucleotide CDR3 lengths of:",  df_len)
  })
  output$length.table_nt <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE), {
    df <- input.data_old2();
    df <- as.data.frame(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt3)),numcolwise(sum)))
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt3])
    df_unique
  })
  Motif_plot2_nt <- reactive( {
    df <- input.data_old2();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt3)),numcolwise(sum)))
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt3])
    df_subset <- subset(df_unique,df_unique$len1==input$len_nt)
    df_subset <- subset(df_subset,get(input$group_column)==input$group_selected)
    motif <- as.data.frame(t(as.data.frame(strsplit(df_subset[,grep(input$aa.or.nt3,names(df_subset))], ""))))
    motif
  })
  output$Motif_plot_nt <- renderPlot( {
    motif <- Motif_plot2_nt()
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    motif_count <- Nucleotide(cbind(x=1,y=2,motif), input$len_nt)
    motif<-new("pcm", mat=motif_count, name="")
    plot(motif)
  })
  output$downloadPlot_motif_nt <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_motif_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_motif_nt,height=input$height_motif_nt, onefile = FALSE) # open the pdf device
      motif <- Motif_plot2_nt()
      motif_count <- Nucleotide(cbind(x=1,y=2,motif), input$len_nt)
      motif<-new("pcm", mat=motif_count, name="")
      plot(motif)
      dev.off()},
    contentType = "application/pdf"
  )
  
  output$downloadPlotPNG_motif_nt <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_motif_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_motif_nt, 
          height = input$height_png_motif_nt, 
          res = input$resolution_PNG_motif_nt)
      motif <- Motif_plot2_nt()
      motif_count <- Nucleotide(cbind(x=1,y=2,motif), input$len_nt)
      motif<-new("pcm", mat=motif_count, name="")
      plot(motif)
      dev.off()},
    
    contentType = "application/png" # MIME type of the image
    
  )
  
  # motif aligned ------
  observe({
    updateSelectInput(
      session,
      "aa.or.nt4",
      choices=names(input.data_old2()),
      selected = "JUNCTION..AA._A")
  })
  observe({
    updateSelectInput(
      session,
      "group_selected_one",
      choices=select_group(),
      selected = "IFN")
  })
  observe({
    updateSelectInput(
      session,
      "group_selected_two",
      choices=select_group(),
      selected = "CD8")
    
  })
  
  select_group_muscle <- function () {
    df <- input.data_old2();
    df <- as.data.frame(df)
    names(df)
    df_unique2 <- as.data.frame(unique (nchar(df[,names(df) %in% input$aa.or.nt4])))
    df_unique2
    names(df_unique2) <- "len"
    
    df_unique2$len<- df_unique2[order(df_unique2$len),]
    df_unique2
  }
  
  observe({
    updateSelectInput(
      session,
      "group_selected_three",
      choices=select_group_muscle(),
      selected = c(11,15))
  })
  
  chain_muscle <- reactive({
    df <- input.data_old2();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    names(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt4)),numcolwise(sum)))
    names(df_unique) <- c("group","chain","cloneCount")
    df_unique <- df_unique[1:3]
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% "chain"])
    
    if (input$restricted.length.range == "yes") {
      df_unique <- subset(df_unique,df_unique$len1==input$group_selected_three)
      x <- AAStringSet(df_unique$chain)
    }
    
    else {
      x <- AAStringSet(df_unique$chain)
      
    }
    
    
    
    
    x <- AAStringSet(df_unique$chain)
    
    if (dim(df_unique)[1] < 501) {
      aln <- muscle(x)
      df1 <- as.data.frame(aln@unmasked)
      df_unique$chain1 <- df1$x
      df_unique
      
    }
    else {
      x <- as.data.frame(c(">500 sequences",
                           "Online tool cannot align >500 sequences",
                           "Use local application as more sequences"
      ))
      names(x) <- "error message"
      x
    } 
    
    
  })
  
  output$Motif_align <- DT::renderDataTable( {
    
    datatable(chain_muscle(), extensions = "Buttons", options = list(searching = TRUE,
                                                                     ordering = TRUE,
                                                                     buttons = c('copy','csv', 'excel'),
                                                                     dom = 'Bfrtip',
                                                                     pageLength=5, 
                                                                     lengthMenu=c(2,5,10,20,50,100), 
                                                                     scrollX = TRUE
    ))
  }, server = FALSE)
  
  chain1_align_aa <- reactive({
    df_unique <- chain_muscle()
    validate(
      need(nrow(df_unique)>0,
           error_message_val1)
    )
    
    motif <- as.data.frame(t(as.data.frame(strsplit(df_unique[,grep("chain1",names(df_unique))], ""))))
    z=dim(motif)[2]
    z
    df_unique1 <- subset(df_unique,df_unique$group==input$group_selected_one)
    motif1 <- as.data.frame(t(as.data.frame(strsplit(df_unique1[,grep("chain1",names(df_unique1))], ""))))
    motif_count1 <- aa.count.function(cbind(x=1,y=2,motif1), z)
    motif_count1<-pcm2pfm(motif_count1)
    as.data.frame(motif_count1)
  })
  chain2_align_aa <- reactive({
    df_unique <- chain_muscle()
    validate(
      need(nrow(df_unique)>0,
           error_message_val1)
    )
    motif <- as.data.frame(t(as.data.frame(strsplit(df_unique[,grep("chain1",names(df_unique))], ""))))
    z=dim(motif)[2]
    z
    df_unique1 <- subset(df_unique,df_unique$group==input$group_selected_two)
    motif1 <- as.data.frame(t(as.data.frame(strsplit(df_unique1[,grep("chain1",names(df_unique1))], ""))))
    motif_count1 <- aa.count.function(cbind(x=1,y=2,motif1), z)
    motif_count1<-pcm2pfm(motif_count1)
    as.data.frame(motif_count1)
    
  })
  chain1_align_nt <- reactive({
    df_unique <- chain_muscle()
    validate(
      need(nrow(df_unique)>0,
           error_message_val1)
    )
    motif <- as.data.frame(t(as.data.frame(strsplit(df_unique[,grep("chain1",names(df_unique))], ""))))
    z=dim(motif)[2]
    z
    df_unique1 <- subset(df_unique,df_unique$group==input$group_selected_one)
    motif1 <- as.data.frame(t(as.data.frame(strsplit(df_unique1[,grep("chain1",names(df_unique1))], ""))))
    motif_count1 <- Nucleotide(cbind(x=1,y=2,motif1), z)
    motif_count1<-pcm2pfm(motif_count1)
    as.data.frame(motif_count1)
  })
  chain2_align_nt <- reactive({
    df_unique <- chain_muscle()
    validate(
      need(nrow(df_unique)>0,
           error_message_val1)
    )
    motif <- as.data.frame(t(as.data.frame(strsplit(df_unique[,grep("chain1",names(df_unique))], ""))))
    z=dim(motif)[2]
    z
    df_unique1 <- subset(df_unique,df_unique$group==input$group_selected_two)
    motif1 <- as.data.frame(t(as.data.frame(strsplit(df_unique1[,grep("chain1",names(df_unique1))], ""))))
    motif_count1 <- Nucleotide(cbind(x=1,y=2,motif1), z)
    motif_count1<-pcm2pfm(motif_count1)
    as.data.frame(motif_count1)
    
  })
  
  Motif_plot_align1 <- reactive({
    motif_count1_aa <- chain1_align_aa()
    motif_count2_aa <- chain2_align_aa()
    motif_count1_nt <- chain1_align_nt()
    motif_count2_nt <- chain2_align_nt()
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    
    
    if (input$compar.lab.motif.all == T) {
      if (input$diff == "compare" && input$aa.nt.col=="ASN") {
        diffLogoObj = createDiffLogoObject(pwm1 = as.data.frame(motif_count1_aa), pwm2 = as.data.frame(motif_count2_aa), alphabet = ASN)
        mat <- (diffLogoObj$pwm1 - diffLogoObj$pwm2)
        names(mat) <- 1:dim(mat)[2]
        
        vals44$plot.ggseq.2 <- ggseqlogo(mat, method='custom', seq_type='aa') + 
          ylab('JS divergence') + 
          geom_hline(yintercept=0) +
          geom_vline(xintercept=0) +
          annotate(geom="text",x=1,y=Inf,vjust=2,label=input$group_selected_one,size=10,face="plain",family=input$font_type)+
          annotate(geom="text",x=1,y=-Inf,vjust=-2,label=input$group_selected_two,size=10,face="plain",family=input$font_type)+
          theme(
            axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
            axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            legend.title  =element_blank(),
            legend.position = "right",
            legend.text = element_text(colour="black", size=12,family=input$font_type)) 
        vals44$plot.ggseq.2
        
      }
      else if (input$diff == "compare" && input$aa.nt.col=="DNA") {
        diffLogoObj = createDiffLogoObject(pwm1 = as.data.frame(motif_count1_nt), pwm2 = as.data.frame(motif_count2_nt), alphabet = DNA)
        mat <- (diffLogoObj$pwm1 - diffLogoObj$pwm2)
        names(mat) <- 1:dim(mat)[2]
        
        vals44$plot.ggseq.2 <- ggseqlogo(mat, method='custom', seq_type='dna') + 
          ylab('JS divergence') + 
          geom_hline(yintercept=0) +
          geom_vline(xintercept=0) +
          theme(
            axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
            axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            legend.title  =element_blank(),
            legend.position = "bottom",
            legend.text = element_text(colour="black", size=12,family=input$font_type)) 
        vals44$plot.ggseq.2
        
      }
      else if (input$diff == "plot_one" && input$aa.nt.col=="ASN") {
        
        
        df_unique <- chain_muscle()
        validate(
          need(nrow(df_unique)>0,
               error_message_val1)
        )
        motif <- as.data.frame(t(as.data.frame(strsplit(df_unique[,grep("chain1",names(df_unique))], ""))))
        z=dim(motif)[2]
        z
        df_unique1 <- subset(df_unique,df_unique$group==input$group_selected_two)
        motif1 <- as.data.frame(t(as.data.frame(strsplit(df_unique1[,grep("chain1",names(df_unique1))], ""))))
        motif_count1 <- aa.count.function(cbind(x=1,y=2,motif1), z)
        motif_count1<-pcm2pfm(motif_count1)
        
        ggseqlogo(motif_count1, seq_type='aa', method='p') + 
          ylab('bits')+ 
          geom_hline(yintercept=0) +
          geom_vline(xintercept=0) +
          theme(
            axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
            axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            legend.title  =element_blank(),
            legend.position = "right",
            legend.text = element_text(colour="black", size=12,family=input$font_type))
        
        
        # motif <- new("pcm", mat=as.matrix(motif_count1_aa), name="")
        # 
        # vals44$plot.ggseq.2 <- seqLogo(as.data.frame(motif_count1_aa), sparse = FALSE, drawLines = 1,
        #                               baseDistribution = probabilities,
        #                               alphabet = ASN, main = NULL)
        # vals44$plot.ggseq.2
      }
      else if (input$diff == "plot_one" && input$aa.nt.col=="DNA") {
        vals44$plot.ggseq.2 <-seqLogo(as.data.frame(motif_count1_nt), sparse = FALSE, drawLines = 1,
                                      baseDistribution = probabilities,
                                      alphabet = DNA, main = NULL)
        vals44$plot.ggseq.2
      }
      else if (input$diff == "plot_two" && input$aa.nt.col=="ASN") {
        vals44$plot.ggseq.2 <- seqLogo(as.data.frame(motif_count2_aa), sparse = FALSE, drawLines = 1,
                                       baseDistribution = probabilities,
                                       alphabet = ASN, main = NULL)
        vals44$plot.ggseq.2
      }
      else {
        vals44$plot.ggseq.2 <-seqLogo(as.data.frame(motif_count2_nt), sparse = FALSE, drawLines = 1,
                                      baseDistribution = probabilities,
                                      alphabet = DNA, main = NULL)
        vals44$plot.ggseq.2
      }
    }
    
    else {
      if (input$diff == "compare" && input$aa.nt.col=="ASN") {
        diffLogoObj = createDiffLogoObject(pwm1 = as.data.frame(motif_count1_aa), pwm2 = as.data.frame(motif_count2_aa), alphabet = ASN)
        mat <- (diffLogoObj$pwm1 - diffLogoObj$pwm2)
        names(mat) <- 1:dim(mat)[2]
        
        vals44$plot.ggseq.2 <- ggseqlogo(mat, method='custom', seq_type='aa') + 
          ylab('JS divergence') + 
          geom_hline(yintercept=0) +
          geom_vline(xintercept=0) +
          theme(
            axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
            axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            legend.title  =element_blank(),
            legend.position = "right",
            legend.text = element_text(colour="black", size=12,family=input$font_type)) 
        vals44$plot.ggseq.2
        
      }
      else if (input$diff == "compare" && input$aa.nt.col=="DNA") {
        diffLogoObj = createDiffLogoObject(pwm1 = as.data.frame(motif_count1_nt), pwm2 = as.data.frame(motif_count2_nt), alphabet = DNA)
        mat <- (diffLogoObj$pwm1 - diffLogoObj$pwm2)
        names(mat) <- 1:dim(mat)[2]
        
        vals44$plot.ggseq.2 <- ggseqlogo(mat, method='custom', seq_type='dna') + 
          ylab('JS divergence') + 
          geom_hline(yintercept=0) +
          geom_vline(xintercept=0) +
          theme(
            axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
            axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            legend.title  =element_blank(),
            legend.position = "bottom",
            legend.text = element_text(colour="black", size=12,family=input$font_type)) 
        vals44$plot.ggseq.2
        
      }
      else if (input$diff == "plot_one" && input$aa.nt.col=="ASN") {
        vals44$plot.ggseq.2 <-seqLogo(as.data.frame(motif_count1_aa), sparse = FALSE, drawLines = 1,
                                      baseDistribution = probabilities,
                                      alphabet = ASN, main = NULL)
        vals44$plot.ggseq.2
      }
      else if (input$diff == "plot_one" && input$aa.nt.col=="DNA") {
        vals44$plot.ggseq.2 <-seqLogo(as.data.frame(motif_count1_nt), sparse = FALSE, drawLines = 1,
                                      baseDistribution = probabilities,
                                      alphabet = DNA, main = NULL)
        vals44$plot.ggseq.2
      }
      else if (input$diff == "plot_two" && input$aa.nt.col=="ASN") {
        vals44$plot.ggseq.2 <- seqLogo(as.data.frame(motif_count2_aa), sparse = FALSE, drawLines = 1,
                                       baseDistribution = probabilities,
                                       alphabet = ASN, main = NULL)
        vals44$plot.ggseq.2
      }
      else {
        vals44$plot.ggseq.2 <-seqLogo(as.data.frame(motif_count2_nt), sparse = FALSE, drawLines = 1,
                                      baseDistribution = probabilities,
                                      alphabet = DNA, main = NULL)
        vals44$plot.ggseq.2
      }
      
    }
    
    
  })
  output$Motif_plot_align <- renderPlot( {
    Motif_plot_align1()
    
  })
  
  output$downloadPlotPNG_motif_align <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_aligned_motif_", gsub("/", "-", x), ".png", sep = "")
    },    
    content = function(file) {
      png(file, width = input$width_png_motif_align, 
          height = input$height_png_motif_align, 
          res = input$resolution_PNG_motif_align)
      
      motif <- Motif_plot_align1()
      withProgress(message = 'Figure is being generated...',
                   detail = '', value = 0, {
                     test_fun()
                   })
      plot(motif)
      
      
      
      
      dev.off()},   contentType = "application/png" # MIME type of the image
  )
  
  
  output$downloadPlot_motif_align <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_aligned_motif_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_motif_align,height=input$height_motif_align, onefile = FALSE) # open the pdf device
      plot(Motif_plot_align1())
      dev.off()
    },
    contentType = "application/pdf"
    
  )
  
  # output$downloadPlotPNG_motif_align <- downloadHandler(
  #   filename = function() {
  #     x <- gsub(":", ".", Sys.time())
  #     paste("TCR_Explore_aligned_motif_", gsub("/", "-", x), ".png", sep = "")
  #   },
  #   content = function(file) {
  #     
  #     png(file, width = input$width_png_motif_align,
  #         height = input$height_png_motif_align,
  #         res = input$resolution_PNG_motif_align)
  #     Motif_plot_align1()
  #     dev.off()
  #     },
  
  #   contentType = "application/png" # MIME type of the image
  #   
  # )
  #
  # pie graph -----
  observe({
    updateSelectInput(
      session,
      "pie_chain",
      choices=names(input.data_old2()),
      selected = "AV")
    
  })
  
  observe({
    updateSelectInput(
      session,
      "string.data.pie.order",
      choices=select_group(),
      selected = c("CD8","IFN")) 
  }) 
  
  cols_pie <- reactive({
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    num <- unique(dat[names(dat) %in% input$pie_chain])
    
    col.gg <- gg_fill_hue(dim(num)[1])
    
    palette_rainbow <- rev(rainbow(dim(num)[1]))
    
   if (input$pie_colour.choise == "default") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), col.gg[i])        
      })
   }
    
    
    else if (input$pie_colour.choise == "rainbow") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), palette_rainbow[i])        
      }) }
    
    
    else if (input$pie_colour.choise == "random") {
      palette1 <- distinctColorPalette(dim(num)[1])
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), palette1[i])        
      })
      
    }
    
    else {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), input$one.colour.default)        
      })
      
      
    }
    
  })
  output$myPanel_pie <- renderUI({cols_pie()})
  
  colors_pie <- reactive({
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    num <- unique(dat[names(dat) %in% input$pie_chain])
    col.gg <- gg_fill_hue(dim(num)[1])
    palette_rainbow <- rev(rainbow(dim(num)[1]))

    lapply(1:dim(num)[1], function(i) {
      input[[paste("col.pie", i, sep="_")]]
    })
  })
  
  # df1 <- dat[names(dat) %in% c(input$count2,input$fill2,input$sub_group2,input$group_column)]
  # df2 <- as.data.frame(ddply(dat,names(df1)[-c(1)],numcolwise(sum)))
  # unique.col <- as.data.frame(unique(dat[names(dat) %in% input$fill2]))
  # names(unique.col) <- "V1"
  # unique.col$tree_palette <- cols
  # 
  # 
  # 
  # df3 <- as.data.frame(merge(df2,unique.col,by.x=input$fill2,by.y = "V1"))
  # 
  # 
  # df3$ID.names <- factor(df3[,names(df3) %in% input$group_column],levels = input$string.data.tree.order)
  
  
  pie_chart <- function() {
    set.seed(123)
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    
    df1 <- dat[names(dat) %in% c("cloneCount",input$pie_chain,input$group_column)]
    df2 <- as.data.frame(ddply(dat,names(df1)[-c(1)],numcolwise(sum)))
    unique.col <- unique(dat[names(dat) %in% input$pie_chain])
    names(unique.col) <- "V1"
    
    cols <- unlist(colors_pie())

    unique.col$palette <- cols
    
    df3 <- as.data.frame(merge(df2,unique.col,by.x=input$pie_chain,by.y = "V1"))
    
    df <- transform(df3, percent = ave(cloneCount, get(input$group_column), FUN = prop.table))
    
    df$group2 <- factor(df3[,names(df3) %in% input$group_column],levels = input$string.data.pie.order)
    
    vals9$pie <- ggplot(df, aes(x="", y=percent, fill=get(input$pie_chain))) +
      geom_bar(width = 1, stat = "identity",aes(colour = "black")) +
      scale_fill_manual(values=unique(df$palette)) +
      scale_color_manual(values = "black") +
      coord_polar("y", start=0) + 
      facet_wrap(~df$group2,nrow = input$nrow.pie) +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title.y= element_blank(),
        legend.position = input$cir.legend,
        legend.text = element_text(size = input$size.circ, family = input$font_type),
        legend.title = element_blank(),
        axis.title = element_blank()) +
      guides(color = "none", size = "none")+
      theme(strip.text = element_text(size = input$panel.text.size.pie, family = input$font_type))+
      theme(panel.background = element_blank()) +
      theme(strip.background =element_rect(fill=input$strip.colour.pie))+
      theme(strip.text = element_text(colour = input$strip.text.colour.pie))
    vals9$pie
    
  }
  
  output$pie_out <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    pie_chart()
  })
  
  output$downloadPlot_pie <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_pie_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_pie,height=input$height_pie, onefile = FALSE) # open the pdf device
      plot(pie_chart())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_pie <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_pie_", gsub("/", "-", x), ".png", sep = "")
    },    
    content = function(file) {
      png(file, width = input$width_png_pie, 
          height = input$height_png_pie, 
          res = input$resolution_PNG_pie)
      plot(pie_chart())
      dev.off()},   contentType = "application/png" # MIME type of the image
  )
  
  
  
  # heatmap ----
  
  # y.axis heatmap
  observe({
    updateSelectInput(
      session,
      "group.heatmap",
      choices=names(input.data_old2()),
      selected = "AVJ")
    
  })
  # x-axis heatmap
  observe({
    updateSelectInput(
      session,
      "heatmap_2",
      choices=names(input.data_old2()),
      selected = "BVJ")
    
  })
  # group 
  observe({
    updateSelectInput(
      session,
      "group_selected3",
      choices=select_group())
  })
  
  vals13 <- reactiveValues(heatmap_clonal2=NULL)
  vals14 <- reactiveValues(heatmap_clonal3=NULL)
  
  heatmap_matrix2 <- function() {
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    dat <- subset(dat, get(input$group_column)==input$group_selected3)
    df <- as.data.frame(ddply(dat,(c(input$heatmap_2,input$group.heatmap)),numcolwise(sum)))
    
    names(df) <- c("group","chain","cloneCount")
    head(df)
    min.FC <- 0
    med.FC <- median(df$cloneCount)
    max.FC <- max(df$cloneCount)
    
    df.1 <- acast(df, group~chain, value.var="cloneCount")
    df.1
    head(df.1)
    df.1[is.na(df.1)] <- 0
    dim(df.1)
    
    # ha = HeatmapAnnotation(text = anno_text(df.1), which = "row", gp = gpar(fontfamily = input$font_type, fontface = "bold"))
    ht <- Heatmap(df.1,
                  heatmap_legend_param = list(title = "count",
                                              title_gp = gpar(fontsize = 10, 
                                                              fontface = "bold",fontfamily=input$font_type),
                                              labels_gp = gpar(fontsize = 10,fontfamily=input$font_type)),
                  col = colorRamp2(c(min.FC,max.FC), c("white",input$col.heatmap)),
                  row_names_gp = grid::gpar(fontsize = input$heat.font.size.row,fontfamily=input$font_type),
                  column_names_gp = grid::gpar(fontsize = input$heat.font.size.col,fontfamily=input$font_type),
                  
    )
    
    draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
    
    
    
    
    
  }
  heatmap_matrix3 <- function() {
    # group by chain dat <- read.csv("test-data/Group/Merged_A-B_2021.08.02.csv")
    dat <- input.data2();
    validate(
      need(nrow(dat)>0,
           error_message_val1)
    )
    dat <- as.data.frame(dat)
    
    df <- as.data.frame(ddply(dat,(c(input$heatmap_2,input$group.heatmap)),numcolwise(sum)))
    names(df) <- c("group","chain","cloneCount")
    head(df)
    
    min.FC <- 0
    med.FC <- median(df$cloneCount)
    max.FC <- max(df$cloneCount)
    
    df.1 <- acast(df, group~chain, value.var="cloneCount")
    df.1
    head(df.1)
    df.1[is.na(df.1)] <- 0
    dim(df.1)
    par(family=input$font_type)
    ht <- Heatmap(df.1,
                  heatmap_legend_param = list(title = "count"),
                  col = colorRamp2(c(min.FC,max.FC), c("white",input$col.heatmap)),
                  row_names_gp = grid::gpar(fontsize = input$heat.font.size.row,fontfamily=input$font_type),
                  column_names_gp = grid::gpar(fontsize = input$heat.font.size.col,fontfamily=input$font_type)
                  
    )
    
    draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
    
  }
  
  output$heatmap_out2 <- renderPlot({
    
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    if (input$group_hm == "yes") {
      
      heatmap_matrix2()
    }
    
    else {
      
      heatmap_matrix3()
      
    }
  })
  
  
  output$downloadPlot_heatmap <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_heatmap_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_heatmap,height=input$height_heatmap, onefile = FALSE) # open the pdf device
      if (input$group_hm == "yes") {
        print(heatmap_matrix2())
      }
      
      else {
        print(heatmap_matrix3())
        
      }
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_heatmap <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_heatmap_", gsub("/", "-", x), ".png", sep = "")
    },    
    content = function(file) {
      png(file, width = input$width_png_heatmap, 
          height = input$height_png_heatmap, 
          res = input$resolution_PNG_heatmap)
      if (input$group_hm == "yes") {
        print(heatmap_matrix2())
      }
      
      else {
        print(heatmap_matrix3())
        
      }
      dev.off()},   contentType = "application/png" # MIME type of the image
  )
  
  # simpson calc -----
  vals11 <- reactiveValues(Simp1=NULL)
  vals12 <- reactiveValues(Simp2=NULL)
  
  observe({
    updateSelectInput(
      session,
      "group_column_simp",
      choices=names(input.data_old2()),
      selected = "Indiv")
    
  })
  observe({
    updateSelectInput(
      session,
      "group_column_simp3",
      choices=names(input.data_old2()),
      selected = "group")
    
  })
  observe({
    updateSelectInput(
      session,
      "group_column_simp2",
      choices=names(input.data_old2()),
      selected = names(input.data_old2())[dim(input.data_old2())[2]] 
      )
  })
  
  inv.simpson.index <- function() {
    dataframe = input.data_old2()
    
    test <- as.data.frame(dataframe)
    validate(
      need(nrow(test)>0,
           "Upload file")
    )
    
    df1 <- ddply(dataframe,c(input$group_column_simp,input$group_column_simp2,input$group_column_simp3),numcolwise(sum))
    df1 <- df1[order(df1$cloneCount, decreasing = T),]
    V1 <- df1[names(df1) %in% input$group_column_simp]
    V2 <- df1[names(df1) %in% input$group_column_simp3]
    V1V2 <- cbind(V1,V2)
    names(V1V2) <- c("V1","V2")
    df1$selected.groups <- paste(V1V2$V1,V1V2$V2,sep="_")
    df.group <- unique(df1[names(df1) %in% "selected.groups"])
    names(df.group) <- "V1"
    
    column.length <- length(df.group$V1)
    df.group2 <- unique(df1[names(df1) %in% input$group_column_simp2])
    names(df.group2) <- "V1"
    row.length <- length(df.group2$V1)
    
    m = matrix(NA,ncol=column.length, nrow=row.length)
    samps <- df.group$V1
    
    for (j in 1:column.length){
      df2 <- subset(df1,df1$selected.groups==samps[j])
      m[,j] <- c(df2$cloneCount, rep(NA, row.length - length(df2$cloneCount)))
    }
    m <- as.data.frame(m)
    names(m) <- samps
    
    richness_samp <- matrix(nrow=1,ncol=dim(m)[2])
    inverse_simp <- matrix(nrow=1,ncol=dim(m)[2])
    inverse_simp_corrected <- matrix(nrow=1,ncol=dim(m)[2])
    shannon_div <- matrix(nrow=1,ncol=dim(m)[2])
    chao1_div <- matrix(nrow=1,ncol=dim(m)[2])
    Pielou_even <- matrix(nrow=1,ncol=dim(m)[2])
    b <- matrix(nrow=1,ncol=dim(m)[2])
    d <- matrix(nrow=1,ncol=dim(m)[2])
    
    for( i in 1:dim(m)[2]) {
      
      samp <- m[,i]
      samp <- na.omit(samp)
      richness_samp[,i] <-specnumber(samp)
      inverse_simp[,i] <- diversity(samp,"invsimpson")
      inverse_simp_corrected[,i] <- diversity(samp,"invsimpson")/(specnumber(samp))
      shannon_div[,i] <- diversity(samp,"shannon")
      Pielou_even[,i] <- diversity(samp,"shannon")/log(specnumber(samp))
      chao1_div[,i] <-chao1(samp, taxa.row = TRUE)
      b[,i] <- sum(samp)
      d[,i] <- nrow(as.data.frame(samp))
    }
    
    a1 <- rbind(b, richness_samp, shannon_div, Pielou_even, inverse_simp,inverse_simp_corrected,chao1_div)  
    a1 <- as.data.frame(a1)
    names(a1) <- names(m)

    
    df_name <- as.data.frame(do.call(rbind, strsplit(as.character(names(m)), "_")))
    df_name$both <- paste(df_name$V1,df_name$V2,sep = "_")
    head(df_name) 
    names(df_name) <- c(input$group_column_simp,input$group_column_simp3,"both")
    
    head(df_name) 
    
    a2 <- as.data.frame(t(a1))
    names(a2) <- c("total # clones","unique # clones", "shannon_div", "Pielou_even", "inverse_simp","inverse_simp_corrected","chao1_div")
    a2
    both <- cbind(a2,df_name)
    # both$get(input$index.type) <- both$inv.simpson.index/both$`total # clones`
    as.data.frame(both)
    
  }
  
  output$table_display.diversity <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    dat <- inv.simpson.index()
    dat <- as.data.frame(dat)
    dat
  })
  
  observe({
    df <- as.list(c(input$group_column_simp,input$group_column_simp3,"both"))
    updateSelectInput(
      session,
      "group.index",
      choices=df,
      selected= "group"
    )
  })
  
  observe({
    df <- as.list(c(input$group_column_simp,input$group_column_simp3,"both"))
    
    updateSelectInput(
      session,
      "group2.index",
      choices=df,
      selected = "Indiv"
    )
    
  })
  
  
  # simp cal other plots ------
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
        colourInput(paste("col.inv.simpson", i, sep="_"), paste(num[i,]), input$one.colour.default)
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
    
    vals11$Simp1 <- ggplot(both,aes(x=get(input$group.index),y=get(input$index.type)))+
      stat_summary(fun.data = quantiles_95, geom = "errorbar", position = position_dodge(1)) +       
      stat_summary(fun.data = middle, geom = "boxplot", position = position_dodge(1))  +
            theme_classic() +
            geom_dotplot(aes(fill=get(input$group2.index)),binaxis = 'y',
                         dotsize = 1,
                         show.legend = T,
                         stackdir = "center", binpositions="all", stackgroups=TRUE
            ) +
      
      scale_fill_manual(values = c(unique.col$simp.inv_palette)) +
            theme(text=element_text(size=20,family=input$font_type),
                  axis.title = element_text(colour="black", size=20,family=input$font_type),
                  axis.text.x = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
                  axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
                  axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
                  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
                  legend.title  =element_blank(),
                  legend.position = input$legend.placement.simp,
                  legend.text = element_text(colour="black", size=input$legend.text.simp,family=input$font_type)) +
            guides(fill=guide_legend(ncol=input$col.num.simp)) +
            xlab("") +
            ylab(input$index.type)
    vals11$Simp1
    
    
      # if (input$IQR_95=="95% CI" ) {
      # 
      # 
      #   
      #   if (input$show_all_dots=="Colour all dots") {
      #     
      #     vals11$Simp1 <- ggplot(both,aes(x=get(input$group.index),y=get(input$index.type)))+
      #       stat_summary(fun.data = quantiles_95, geom = "errorbar", position = position_dodge(1)) +       
      #       stat_summary(fun.data = middle, geom = "boxplot", position = position_dodge(1))  
      #     
      #     # vals11$Simp1 +
      #     #   theme_classic() +
      #     #   geom_dotplot(aes(fill=get(input$group2.index)),binaxis = 'y',
      #     #                dotsize = 1,
      #     #                show.legend = T,
      #     #                stackdir = "center", binpositions="all", stackgroups=TRUE
      #     #   ) +
      #     #   scale_fill_manual(values = c(unique.col$simp.inv_palette)) +
      #     #   theme(text=element_text(size=20,family=input$font_type),
      #     #         axis.title = element_text(colour="black", size=20,family=input$font_type),
      #     #         axis.text.x = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
      #     #         axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
      #     #         axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
      #     #         axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
      #     #         legend.title  =element_blank(),
      #     #         legend.position = input$legend.placement.simp,
      #     #         legend.text = element_text(colour="black", size=input$legend.text.simp,family=input$font_type)) +
      #     #   guides(fill=guide_legend(ncol=input$col.num.simp)) +
      #     #   xlab("")+
      #     #   ylab("Inverse SDI (corrected)")
      #     
      #   }
      #   
      # 
      # 
      # else {
      #   vals11$Simp1 <- ggplot(both,aes(x=get(input$group.index),y=get(input$index.type)))+
      #     stat_summary(fun.data = quantiles_95, geom = "errorbar", position = position_dodge(1)) +       
      #     stat_summary(fun.data = middle, geom = "boxplot", position = position_dodge(1))  
      #   # vals11$Simp1 +
      #   #   
      #   #   stat_summary(fun.y = o, geom="point", position = position_dodge(1), shape = 1) +
      #   #   theme_classic() +
      #   #   scale_fill_manual(values = c(unique.col$simp.inv_palette)) +
      #   #   theme(text=element_text(size=20,family=input$font_type),
      #   #         axis.title = element_text(colour="black", size=20,family=input$font_type),
      #   #         axis.text.x = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
      #   #         axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
      #   #         axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
      #   #         axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
      #   #         legend.title  =element_blank(),
      #   #         legend.position = input$legend.placement.simp,
      #   #         legend.text = element_text(colour="black", size=input$legend.text.simp,family=input$font_type)) +
      #   #   guides(fill=guide_legend(ncol=input$col.num.simp)) +
      #   #   xlab("")+
      #   #   ylab("Inverse SDI (corrected)")
      # }
      # 
      # 
      # }
      # 
      # else {
      #   
      #   vals11$Simp1 <- ggplot(both,aes(x=get(input$group.index),y=get(input$index.type)))+
      #     geom_boxplot(show.legend = F)+
      #     geom_dotplot(aes(fill=get(input$group2.index)),binaxis = 'y',
      #                  dotsize = 1,
      #                  show.legend = T,
      #                  stackdir = "center", binpositions="all", stackgroups=TRUE
      #     ) +
      #     theme_classic() +
      #     scale_fill_manual(values = c(unique.col$simp.inv_palette)) +
      #     theme(text=element_text(size=20,family=input$font_type),
      #           axis.title = element_text(colour="black", size=20,family=input$font_type),
      #           axis.text.x = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
      #           axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
      #           axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
      #           axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
      #           legend.title  =element_blank(),
      #           legend.position = input$legend.placement.simp,
      #           legend.text = element_text(colour="black", size=input$legend.text.simp,family=input$font_type)) +
      #     guides(fill=guide_legend(ncol=input$col.num.simp)) +
      #     xlab("")+
      #     ylab("Inverse SDI (corrected)")
      #   
      #   vals11$Simp1
      #   
      # }
    
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
    
    vals12$Simp2 <- ggplot(both,aes(x=get(input$x.axis.index), y=get(input$index.type),color=get(input$group2.index)))+
      geom_point(size =3, alpha =1, show.legend =T)+
      theme_bw() +
      scale_color_manual(values=unique.col$simp.inv_palette) +
      theme(text=element_text(size=20,family=input$font_type),
            axis.title = element_text(colour="black", size=20,family=input$font_type),
            axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
            axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
            legend.title  =element_blank(),
            legend.position = "none",
            legend.text = element_text(colour="black", size=8,family=input$font_type)) +
      # scale_alpha(guide = 'none') +
      scale_x_continuous(labels = function(x) format(x, scientific = input$scale_x_continuous_x))+
      labs(x=input$x.axis.index,
           y=input$index.type)
    
    vals12$Simp2
    
  }
  
  # group.diversity2 <- function() {
  #   both <- inv.simpson.index()
  #   cols <- unlist(colors_inv.simp())
  #   validate(
  #     need(nrow(both)>0,
  #          "select correct chain")
  #   )
  #   
  #   both <- as.data.frame(both)
  #   
  #   selected.col <- both[names(both) %in% input$group2.index]
  #   names(selected.col) <- "V1"
  #   both[names(both) %in% input$group2.index] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))
  #   
  #   unique.col <- as.data.frame(unique(both[names(both) %in% input$group2.index]))
  #   names(unique.col) <- "V1"
  #   unique.col$simp.inv_palette <- cols
  #   if (input$index.type == "Sample size corrected Inverse SDI") {
  #     if (input$scale_x_continuous_x=="scientific") {
  #       vals12$Simp2 <- ggplot(both,aes(x=get(input$x.axis.index), y=get(input$index.type),color=get(input$group2.index)))+
  #         geom_point(size =3, alpha =1, show.legend =T)+
  #         theme_bw() +
  #         scale_color_manual(values=unique.col$simp.inv_palette) +
  #         theme(text=element_text(size=20,family=input$font_type),
  #               axis.title = element_text(colour="black", size=20,family=input$font_type),
  #               axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
  #               axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               legend.title  =element_blank(),
  #               legend.position = "none",
  #               legend.text = element_text(colour="black", size=8,family=input$font_type)) +
  #         # scale_alpha(guide = 'none') +
  #         scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  #         labs(x="number of clones",
  #              y="Inverse SDI (corrected)")
  #       
  #       vals12$Simp2
  #     }
  #     else {
  #       vals12$Simp2 <- ggplot(both,aes(x=get(input$x.axis.index), y=get(input$index.type),color=get(input$group2.index)))+
  #         geom_point(size =3, alpha =1, show.legend =T)+
  #         theme_bw() +
  #         scale_color_manual(values=unique.col$simp.inv_palette) +
  #         theme(text=element_text(size=20,family=input$font_type),
  #               axis.title = element_text(colour="black", size=20,family=input$font_type),
  #               axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
  #               axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               legend.title  =element_blank(),
  #               legend.position = "none",
  #               legend.text = element_text(colour="black", size=8,family=input$font_type)) +
  #         labs(x="number of clones",
  #              y="Inverse SDI (corrected)")
  #       
  #       vals12$Simp2
  #     }}
  #   else {
  #     if (input$scale_x_continuous_x=="scientific") {
  #       vals12$Simp2 <- ggplot(both,aes(x=get(input$x.axis.index), y=inv.simpson.index,color=get(input$group2.index)))+
  #         geom_point(size =3, alpha =1, show.legend =T)+
  #         theme_bw() +
  #         scale_color_manual(values=unique.col$simp.inv_palette) +
  #         theme(text=element_text(size=20,family=input$font_type),
  #               axis.title = element_text(colour="black", size=20,family=input$font_type),
  #               axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
  #               axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               legend.title  =element_blank(),
  #               legend.position = "none",
  #               legend.text = element_text(colour="black", size=8,family=input$font_type)) +
  #         # scale_alpha(guide = 'none') +
  #         scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  #         labs(x="number of clones",
  #              y="Inverse SDI")
  #       
  #       vals12$Simp2
  #     }
  #     else {
  #       vals12$Simp2 <- ggplot(both,aes(x=get(input$x.axis.index), y=inv.simpson.index,color=get(input$group2.index)))+
  #         geom_point(size =3, alpha =1, show.legend =T)+
  #         theme_bw() +
  #         scale_color_manual(values=unique.col$simp.inv_palette) +
  #         theme(text=element_text(size=20,family=input$font_type),
  #               axis.title = element_text(colour="black", size=20,family=input$font_type),
  #               axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type),
  #               axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type),
  #               legend.title  =element_blank(),
  #               legend.position = "none",
  #               legend.text = element_text(colour="black", size=8,family=input$font_type)) +
  #         labs(x="number of clones",
  #              y="Inverse SDI")
  #       
  #       vals12$Simp2
  #     }}
  #   
  #   
  # }
  
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
  
  select_group2 <- function () {
    df <- input.data2();
    
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    
    df2 <- as.data.frame(unique(df[names(df) %in% input$group1_column]))
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
  observe({
    df <- as.list(c(input$group_column_simp,input$group_column_simp3,"both"))
    
    updateSelectInput(
      session,
      "group1_column",
      choices=df,
      selected = "group"
    )
    
  })
  ttestout <- reactive({
    dat <- table.inv.simpson()
    conf <- input$conf
    dat <- dat[order(dat$both),]
    ve <- ifelse(input$varequal == 'y', TRUE, FALSE)
    pair_samp <- ifelse(input$paired == 'y', TRUE, FALSE)
    group1 <- subset(dat, get(input$group1_column)==input$group1_selected) # group 1
    group2 <- subset(dat, get(input$group1_column)==input$group2_selected) # group 2
    
if (input$test_ttest == "parametric") {
  t.test(group1[,names(group1) %in% input$index.type],group2[,names(group2) %in% input$index.type], 
         paired = pair_samp, var.equal = ve, alternative = input$tail,conf.level = conf)
}
      
      else {
        wilcox.test(group1[,names(group1) %in% input$index.type],group2[,names(group2) %in% input$index.type], 
               paired = pair_samp)
      }
    
  })
  
  AOVtestout <- reactive({
    dat <-  table.inv.simpson()
    dat2 <- dat[names(dat) %in% c(input$index.type,input$group1_column)]
    
    aov(get(input$index.type) ~ get(input$group1_column),data = dat2)
  })
  
  postAOVtestout <- reactive({
    dat <-  table.inv.simpson()
    dat2 <- dat[names(dat) %in% c(input$index.type,input$group1_column)]
    
    TukeyHSD(aov(get(input$index.type) ~ get(input$group1_column),data = dat2))
  })
  
  output$tvalue <- renderPrint({
    vals <- ttestout()
    if (is.null(vals)){return(NULL)}
    vals
  })
  
  output$aov_test <- renderPrint({
    vals <- AOVtestout()
    if (is.null(vals)){return(NULL)}
    summary(vals)
  })
  output$postaov_test <- renderPrint({
    vals <- postAOVtestout()
    if (is.null(vals)){return(NULL)}
    vals
  })

  output$downloadTABLE_simpson.inv <- downloadHandler(
    filename = function(){
      paste("inv.simpson.index",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      write.csv(table.inv.simpson(),file, row.names = FALSE)
    })
  
  output$downloadPlot_simpson.inv <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("inv.simpson.index.",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_simpson.inv,height=input$height_simpson.inv, onefile = FALSE) # open the pdf device
      grid.arrange(print(group.diversity1()),print(group.diversity2()),ncol=2)
      # group.diversity1()
      dev.off()},
    contentType = "application/pdf" )
  
  output$downloadPlotPNG_simpson.inv <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("inv.simpson.index.", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_simpson.inv, height = input$height_png_simpson.inv, res = input$resolution_PNG_simpson.inv)
      grid.arrange(print(group.diversity1()),print(group.diversity2()),ncol=2)
      # group.diversity1()
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  # FACS index data -----
  input.data_FACS <- reactive({switch(input$dataset3,"test-FACS" = test.data_FACS(), "own_FACS" = own.data_FACS())})
  test.data_FACS <- reactive({
    read.FCS("test-data/Index/Murine Lymph Node_INX_780 Fib index 2_001_018.fcs")
  })
  own.data_FACS <- reactive({
    input$file_FACS
    if (is.null(input$file_FACS)) return(NULL)
    
    else {
      read.FCS(input$file_FACS$datapath)}
    
  })
  output$FACS_to_CSV <- DT::renderDataTable(escape = FALSE,{
    
    samp <-  input.data_FACS()
    validate(
      need(nrow(samp)>0,
           error_message_val1)
    )
    req(samp)
    samp_index <- getIndexSort(samp)
    samp_index
  })
  output$FACS.CSV <- DT::renderDataTable( {
    samp <-  input.data_FACS()
    validate(
      need(nrow(samp)>0,
           error_message_val1)
    )
    
    req(samp)
    samp_index <- getIndexSort(samp)
    samp_index <- as.data.frame(samp_index[1:dim(samp_index)[1],])
    samp_index$group <- input$group_FACS
    samp_index$Indiv <- input$indiv_FACS
    samp_index$plate <- input$Plate_FACS
    datatable(samp_index[1:dim(samp_index)[1],], extensions = "Buttons", options = list(searching = TRUE,
                                                                                        ordering = TRUE,
                                                                                        buttons = c('copy','csv', 'excel'),
                                                                                        dom = 'Bfrtip',
                                                                                        pageLength=5, 
                                                                                        lengthMenu=c(2,5,10,20,50,100), 
                                                                                        scrollX = TRUE
    ))
  }, server = FALSE)  
  
  
  select_group_FACS <- function () {
    df <- input.data.clone.file()
    df <- as.data.frame(df)
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    
    df2 <- as.data.frame(unique(df[names(df) %in% "group"]))
    df2 <- as.data.frame(df2)
    #names(df2) <- "V1"
    df2
  }
  
  
  observe({
    updateSelectInput(
      session,
      "group_FACS",
      choices=select_group_FACS(),
      selected = c("other"))
  })
  
  select_group_FACS2 <- function () {
    df <- input.data.clone.file()
    df <- as.data.frame(df)
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    
    df2 <- as.data.frame(unique(df[names(df) %in% "Indiv"]))
    df2 <- as.data.frame(df2)
    #names(df2) <- "V1"
    df2
  }
  
  
  observe({
    updateSelectInput(
      session,
      "indiv_FACS",
      choices=select_group_FACS2(),
      selected = c("780"))
  })
  
  
  
  
  merged.index <- function(){
    samp <-  input.data_FACS()
    req(samp)
    samp_index <- getIndexSort(samp)
    samp_index <- as.data.frame(samp_index[1:dim(samp_index)[1],])
    samp_index$group <- input$group_FACS
    samp_index$Indiv <- input$indiv_FACS
    samp_index$plate <- input$Plate_FACS
    
    replace_ID <- read.csv("test-data/Index/Loc_to_ID.csv")
    head(replace_ID)
    index_updated_ID <- merge(samp_index,replace_ID,by=c("XLoc","YLoc"))
    index_updated_ID
    index_updated_ID$plate.well <- paste(index_updated_ID$plate,index_updated_ID$well,sep="")
    index_updated_ID
    
  }
  
  output$merged.clone <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    merged.index()
  })
  input.data.clone.file <- reactive({switch(input$data_clone.index, "ab.test.clone3" = test.data.gd.index.csv2(),"own.clone.file" = own.data.clone.file.csv())})
  test.data.gd.index.csv2 <- reactive({
    dataframe = read.csv("test-data/Index/DR4-780 TCR sequence data.csv",header=T) 
  })
  
  own.data.clone.file.csv <- reactive({
    inFile21 <- input$file_diversity.index.2
    if (is.null(inFile21)) return(NULL)
    
    else {
      dataframe <- read.csv(
        inFile21$datapath)}
    
  })
  
  with.clone.data <- function () {
    index_updated_ID <- as.data.frame(merged.index())
    clonal.file <-  as.data.frame(input.data.clone.file())
    clonal.file <- as.data.frame(clonal.file)
    
    if (input$multiple_plates == T) {
      index_updated_ID$well <- index_updated_ID$plate.well
      index.clonal.file <- merge(clonal.file,index_updated_ID,by=c("well","Indiv","group"))
      index.clonal.file <- index.clonal.file[!names(index.clonal.file) %in% c("row","column","XLoc","YLoc",'name',"plate","plate.well")]
      index.clonal.file <- index.clonal.file[ , -which(names(index.clonal.file) %in% c(names(index.clonal.file)[grep("Sequence",names(index.clonal.file))],
                                                                                       names(index.clonal.file)[grep("*allele*",names(index.clonal.file))],
                                                                                       names(index.clonal.file)[grep("JUNCTION_",names(index.clonal.file))],
                                                                                       names(index.clonal.file)[grep("IMGT",names(index.clonal.file))],
                                                                                       names(index.clonal.file)[grep("GENE",names(index.clonal.file))] ))]
      index.clonal.file
    }
    else {
      index.clonal.file <- merge(clonal.file,index_updated_ID,by=c("well","Indiv","group"))
      index.clonal.file <- index.clonal.file[!names(index.clonal.file) %in% c("row","column","XLoc","YLoc",'name',"plate","plate.well")]
      index.clonal.file <- index.clonal.file[ , -which(names(index.clonal.file) %in% c(names(index.clonal.file)[grep("Sequence",names(index.clonal.file))],names(index.clonal.file)[grep("*allele*",names(index.clonal.file))],names(index.clonal.file)[grep("JUNCTION_",names(index.clonal.file))],names(index.clonal.file)[grep("IMGT",names(index.clonal.file))] ,names(index.clonal.file)[grep("GENE",names(index.clonal.file))] ))]
      
      index.clonal.file
      
    }
    
  }
  
  output$merged.index.clone <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    samp <-  input.data_FACS()
    validate(
      need(nrow(samp)>0,
           "upload FACS file")
    )
    clonal.file <-  input.data.clone.file();
    
    validate(
      need(nrow(clonal.file)>0,
           "Upload clone file")
    )
    with.clone.data()
    
    
  })
  
  output$downloadTABLE_FACS <- downloadHandler(
    filename = function(){
      paste(input$name.colour2,"TCR_Explore_index.clonal.",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      write.csv(with.clone.data(),file, row.names = FALSE)
    }
  )
  
  
  # colouring columns -----
  input.data_CSV1 <-  reactive({switch(input$dataset7,"test-csv"=test.data_csv1(),"own_csv" = own.data_CSV1())})
  test.data_csv1 <- reactive({
    dataframe = read.csv("test-data/Index/TCR_Explore_index.clonal.2021.11.19.csv",header = T)
  })
  own.data_CSV1 <- reactive({
    inFile_CSV1 <- input$file_FACS.csv1
    if (is.null(inFile_CSV1)) return(NULL)
    
    else {
      dataframe <- read.csv(inFile_CSV1$datapath, header=T)}
    
  })
  
  vals15 <- reactiveValues(complex_dot=NULL)
  vals16 <- reactiveValues(complex_dot2=NULL)
  
  output$names.in.file <- renderPrint( {
    df <- input.data_CSV1();
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    df <- as.data.frame(df)
    names(df)
    
  })
  
  observe({
    updateSelectInput(
      session,
      "group.col.dot",
      choices=names(input.data_CSV1()),
      selected = "Indiv")
  })
  observe({
    updateSelectInput(
      session,
      "V.gene.1",
      choices=names(input.data_CSV1()),
      selected = "AJ")
  })
  observe({
    updateSelectInput(
      session,
      "CDR3.1",
      choices=names(input.data_CSV1()),
      selected = "CDR3a.Sequence")
  })
  observe({
    updateSelectInput(
      session,
      "V.gene.2",
      choices=names(input.data_CSV1()),
      selected = "BJ")
  })
  observe({
    updateSelectInput(
      session,
      "CDR3.2",
      choices=names(input.data_CSV1()),
      selected = "CDR3b.Sequence")
  })
  observe({
    updateSelectInput(
      session,
      "string.data",
      choices=names(input.data_CSV1()),
      selected = c("Indiv","group","TRBV","CDR3b.Sequence","TRBJ","TRAV","CDR3a.Sequence", "TRAJ","AJ", "BJ","AJBJ")) 
  }) 
  
  index.cleaning1 <- reactive({
    df <- input.data_CSV1();
    validate(
      need(nrow(df)>0,
           error_message_val2)
    )
    df <- as.data.frame(df)
    head(df)
    
    
    df2 <- df[,c("cloneCount",input$string.data)] 
    df2
    df3 <- as.data.frame(ddply(df2,input$string.data,numcolwise(sum)))
    df3
    df1 <- subset(df3,df3$cloneCount>input$numeric.cloneCount)
    df1$clonal <- "yes"
    colnames(df1)[which(names(df1) == "cloneCount")] <- "# of clones"
    a <- df1[names(df1) %in% input$V.gene.1]
    names(a) <- "V1"
    b <- df1[names(df1) %in% input$CDR3.1]
    names(b) <- "V1"
    group.CDR <- df1[names(df1) %in% input$group.col.dot]
    names(group.CDR) <- "V1"
    d <- df1[names(df1) %in% input$V.gene.2]
    names(d) <- "V1"
    e  <- df1[names(df1) %in% input$CDR3.2]
    names(e) <- "V1"
    df1$gene.CDR3.1 <- paste(a$V1,b$V1,group.CDR$V1,sep="_")
    df1$gene.CDR3.2 <- paste(d$V1,e$V1,group.CDR$V1,sep = "_")
    df1$gene.CDR3.both <- paste(a$V1,b$V1,d$V1,e$V1,group.CDR$V1,sep = "_")
    df1
    a2 <- merge(df,df1,by=input$string.data,all=T)
    #a[is.na(a)] <- 'unique'
    a2 <- as.data.frame(a2)
    a2[a2< -10000] <- 0.0001
    a2[a2< -9000] <- 0.0002
    a2[a2< -8000] <- 0.0003
    a2[a2< -7000] <- 0.0004
    a2[a2< -6000] <- 0.0005
    a2[a2< -5000] <- 0.0006
    a2[a2< -4000] <- 0.0007
    a2[a2< -3000] <- 0.0008
    a2[a2< -2000] <- 0.0009
    a2[a2< -1000] <- 0.0010
    a2[a2< -900] <- 0.0020
    a2[a2< -800] <- 0.0030
    a2[a2< -700] <- 0.0040
    a2[a2< -600] <- 0.0050
    a2[a2< -500] <- 0.0060
    a2[a2< -400] <- 0.0070
    a2[a2< -300] <- 0.0080
    a2[a2< -200] <- 0.0090
    a2[a2< -190] <- 0.0091
    a2[a2< -180] <- 0.0092
    a2[a2< -170] <- 0.0093
    a2[a2< -160] <- 0.0094
    a2[a2< -150] <- 0.0095
    a2[a2< -140] <- 0.0096
    a2[a2< -130] <- 0.0097
    a2[a2< -120] <- 0.0098
    a2[a2< -110] <- 0.0099
    a2[a2< -100] <- 0.0100
    a2[a2< -90] <- 0.020
    a2[a2< -80] <- 0.030
    a2[a2< -70] <- 0.040
    a2[a2< -60] <- 0.050
    a2[a2< -50] <- 0.060
    a2[a2< -40] <- 0.070
    a2[a2< -30] <- 0.080
    a2[a2< -20] <- 0.090
    a2[a2< -10] <- 0.1
    a2[a2< -9] <- 0.2
    a2[a2< -8] <- 0.3
    a2[a2< -7] <- 0.4
    a2[a2< -6] <- 0.5
    a2[a2< -5] <- 0.6
    a2[a2< -4] <- 0.7
    a2[a2< -3] <- 0.8
    a2[a2< -2] <- 0.9
    a2[a2< 0] <- 1
    a2
  })
  
  output$table.index.1 <- DT::renderDataTable(escape = FALSE,options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    df <- index.cleaning1()
    df
    
  })
  
  observe({
    updateSelectInput(
      session,
      "string.data_UMAP",
      choices=names(input.data_CSV1()),
      selected = c("dump.fitc","CD8.pb","LD.aqua","CD4.buv395","tcr.beta.APC.Cy7",	"Tetramer.1.APC","Tetramer.2.PE","CD62L.bv605")) 
  }) 
  
  index.cleaning.UMAP <- reactive({
    UMAP <- index.cleaning1()
    df2 <- UMAP[,c(input$string.data_UMAP)] 
    mydata <- na.omit(df2) # listwise deletion of missing
    mydata <- log10(mydata) # standardize variables
    
    umap.data <- umap(mydata)
    fit <- pamk(mydata, krange= input$lower_range:input$upper_range, ns=input$upper_range)
    # append cluster assignment
    mydata <- data.frame(UMAP, fit$pamobject$clustering,umap.data$layout)
    aggregate(df2,by=list(fit$pamobject$clustering),FUN=mean)
    
    y = dim(mydata)[2]
    x = y-2
    names(mydata)[x:y] <- c("cluster","umap.1","umap.2")
    names(mydata) <- ifelse(grepl("umap",names(mydata)),toupper(names(mydata)),names(mydata))
    mydata
  }) 
  
  index.cleaning.UMAP2 <- reactive({
    UMAP <- index.cleaning1()
    df2 <- UMAP[,c(input$string.data_UMAP)] 
    mydata <- na.omit(df2) # listwise deletion of missing
    mydata <- log10(mydata) # standardize variables
    
    umap.data <- umap(mydata)
    fit <- pamk(mydata, krange= input$lower_range:input$upper_range, ns=input$upper_range)
    # append cluster assignment
    mydata <- data.frame(UMAP, fit$pamobject$clustering,umap.data$layout)
    names(mydata) <- ifelse(grepl("umap",names(mydata)),toupper(names(mydata)),names(mydata))
    
    aggregate(df2,by=list(fit$pamobject$clustering),FUN=mean)
  }) 
  
  output$table.index.UMAP <- DT::renderDataTable(escape = FALSE,options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
    df <- index.cleaning.UMAP()
    df
    
  })
  
  output$table.index.UMAP2 <- DT::renderDataTable( {
    df <- log10(index.cleaning.UMAP2())
    df$Group.1 <- index.cleaning.UMAP2()$Group.1
    
    datatable(as.data.frame(df), extensions = "Buttons", options = list(searching = TRUE,
                                                                            ordering = TRUE,
                                                                            buttons = c('copy','csv', 'excel'),
                                                                            dom = 'Bfrtip',
                                                                            pageLength=20, 
                                                                            lengthMenu=c(2,5,10,20,50,100), 
                                                                            scrollX = TRUE
    ))
  }, server = FALSE)
    
  
  output$downloadTABLE_cleaning <- downloadHandler(
    filename = function(){
      paste(input$name.colour,"colouring column",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      write.csv(index.cleaning.UMAP(),file, row.names = FALSE)
    })
  
  
  # creating the dot plot ----
  input.data_CSV2 <-  reactive({switch(input$dataset_index.2,"test-csv"=test.data_csv2(),"own_csv_file" = own.data_CSV2())})
  test.data_csv2<- reactive({
    dataframe = read.csv("test-data/Index/ID.780_colouring column2022.12.29.csv",header = T)
  })
  own.data_CSV2 <- reactive({
    inFile_CSV2 <- input$file_FACS.csv2
    if (is.null(inFile_CSV2)) return(NULL)
    
    else {
      dataframe <- read.csv(inFile_CSV2$datapath, header=T)}
    
  })
  
  observe({
    dat <- input.data_CSV2()
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
            )
    names(dat) <- gsub("\\.", " ", names(dat))
    
    if (input$plot_type_umap == "overlaid dot plot") {
      updateSelectInput(
        session,
        "x.axis2",
        choices=names(dat),
        selected = "Tetramer 2 PE")
    }
    
   
    else if (input$plot_type_umap == "Ridge plot") {
      updateSelectInput(
        session,
        "x.axis2",
        choices=names(dat),
        selected = "Tetramer 2 PE")
    }
    
    else {
      updateSelectInput(
        session,
        "x.axis2",
        choices=names(dat),
        selected = "UMAP 1")
    
    }

  })
  
  observe({
    dat <- input.data_CSV2()
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
    )
    names(dat) <- gsub("\\.", " ", names(dat))
    
    if (input$plot_type_umap == "overlaid dot plot") {
      updateSelectInput(
        session,
        "y.axis2",
        choices=names(dat),
        selected = "Tetramer 1 APC")
    }
    
    else if (input$plot_type_umap == "Ridge plot") {
      updateSelectInput(
        session,
        "y.axis2",
        choices=names(dat),
        selected = "Tetramer 1 APC")
    }
    
    else {
      updateSelectInput(
        session,
        "y.axis2",
        choices=names(dat),
        selected = "UMAP 2")
    }
  })
  
  observe({
    dat <- input.data_CSV2()
    
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
    )
    
    names(dat) <- gsub("\\.", " ", names(dat))
    updateSelectInput(
      session,
      "group_complex_dot",
      choices=names(dat),
      selected = "AJBJ")
  })
  
  cols.FACS.index <- reactive({
    set.seed(123)
    dat <- input.data_CSV2()
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
    )
    
    names(dat) <- gsub("\\.", " ", names(dat))
    
    
    dat <- as.data.frame(dat)
    dat[is.na(dat)] <- "not_clonal"
    num <- unique(dat[names(dat) %in% input$group_complex_dot])
    col.gg <- gg_fill_hue(dim(num)[1])
    unique.col <- as.data.frame(unique(dat[grep(input$group_complex_dot,names(dat))]))
    
    
    if (input$FACS.index_colour.choise == "default") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.FACS.index", i, sep="_"), paste(num[i,]), col.gg[i])        
      })
    }
    else if (input$FACS.index_colour.choise == "random") {
      palette1 <- distinctColorPalette(dim(unique.col)[1])
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.FACS.index", i, sep="_"), paste(num[i,]), palette1[i])        
      })
      
    }
    
    else {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.FACS.index", i, sep="_"), paste(num[i,]), "grey")        
      })
      
      
    }
    
  })
  shape.FACS.index <- reactive({
    dat <- input.data_CSV2()
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
    )
    names(dat) <- gsub("\\.", " ", names(dat))
    
    validate(
      need(nrow(dat)>0,
           "Upload file")
    )
    dat <- as.data.frame(dat)
    dat[is.na(dat)] <- "not_clonal"
    num <- unique(dat[names(dat) %in% input$group_complex_dot])
    lapply(1:dim(num)[1], function(i) {
      numericInput(paste("shape.FACS.index", i, sep="_"), paste(num[i,]), 19)        
    })
  })
  size.FACS.index <- reactive({
    dat <- input.data_CSV2()
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
    )
    names(dat) <- gsub("\\.", " ", names(dat))
    
    validate(
      need(nrow(dat)>0,
           "Upload file")
    )
    dat <- as.data.frame(dat)
    dat[is.na(dat)] <- "not_clonal"
    num <- unique(dat[names(dat) %in% input$group_complex_dot])
    lapply(1:dim(num)[1], function(i) {
      numericInput(paste("size.FACS.index", i, sep="_"), paste(num[i,]), 3)        
    })
    
  })
  alpha.FACS.index <- reactive({
    dat <- input.data_CSV2();
    
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
    )
    dat <- as.data.frame(dat)
    dat[is.na(dat)] <- "not_clonal"
    num <- unique(dat[names(dat) %in% input$group_complex_dot])
    lapply(1:dim(num)[1], function(i) {
      numericInput(paste("alpha.FACS.index", i, sep="_"), paste(num[i,]), 1)        
    })
  })
  
  output$myPanel.FACS.index.shape <- renderUI({shape.FACS.index()})
  output$myPanel.FACS.index <- renderUI({cols.FACS.index()})
  output$myPanel.FACS.index.size <- renderUI({size.FACS.index()})
  
  colors.FACS.index <- reactive({
    set.seed(123)
    dat <- input.data_CSV2()
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
    )
    names(dat) <- gsub("\\.", " ", names(dat))
    validate(
      need(nrow(dat)>0,
           "Upload file")
    )
    dat <- as.data.frame(dat)
    dat[is.na(dat)] <- "not_clonal"
    
    num <- unique(dat[names(dat) %in% input$group_complex_dot])
    lapply(1:dim(num)[1], function(i) {
      input[[paste("col.FACS.index", i, sep="_")]]
    })
  })
  shape.FACS.index2 <- reactive({
    dat <- input.data_CSV2()
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
    )
    names(dat) <- gsub("\\.", " ", names(dat))
    validate(
      need(nrow(dat)>0,
           "Upload file")
    )
    dat <- as.data.frame(dat)
    dat[is.na(dat)] <- "not_clonal"
    
    num <- unique(dat[names(dat) %in% input$group_complex_dot])
    lapply(1:dim(num)[1], function(i) {
      input[[paste("shape.FACS.index", i, sep="_")]]
    })
  })
  size.FACS.index2 <- reactive({
    dat <- input.data_CSV2()
    validate(
      need(nrow(dat)>0,
           "Upload file for dotplot")
    )
    names(dat) <- gsub("\\.", " ", names(dat))
    validate(
      need(nrow(dat)>0,
           "Upload file")
    )
    dat <- as.data.frame(dat)
    dat[is.na(dat)] <- "not_clonal"
    
    num <- unique(dat[names(dat) %in% input$group_complex_dot])
    lapply(1:dim(num)[1], function(i) {
      input[[paste("size.FACS.index", i, sep="_")]]
    })
  })
  
  dot_plot.complex <- reactive({
    set.seed(123)
    index <- input.data_CSV2();
    validate(
      need(nrow(index)>0,
           "Upload file")
    )
    names(index) <- gsub("\\.", " ", names(index))
    
    
    index <- as.data.frame(index)
    y_lable1 <- bquote(.(input$y.axis2))
    x_lable1 <-  bquote(.(input$x.axis2))
    
    index[is.na(index)] <- "not_clonal"
    selected.col <- index[names(index) %in% input$group_complex_dot]
    names(selected.col) <- "V1"
    index[names(index) %in% input$group_complex_dot] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))
    palette.complex <- unlist(colors.FACS.index())
    shape.ggplot <- unlist(shape.FACS.index2())
    size.ggplot <- unlist(size.FACS.index2())
    
    names_unique <- as.data.frame(unique(selected.col))
    names_unique_size <- cbind(names_unique,as.data.frame(size.ggplot))
    
    names(names_unique_size)[1] <- input$group_complex_dot
    
    index2 <- merge(index,names_unique_size,by=input$group_complex_dot)
    
    vals15$complex_dot <- ggplot(index2, aes(x=get(input$x.axis2), y=get(input$y.axis2),
                                            colour = get(input$group_complex_dot),
                                            fill = get(input$group_complex_dot),
                                            shape = get(input$group_complex_dot), 
                                            size = get(input$group_complex_dot), 
                                            ))+
      geom_point() + 
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    limits = c(input$min.x,10^input$max.x),
                    labels = trans_format("log10", math_format(10^.x))) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    limits = c(input$min.y,10^input$max.y),
                    labels = trans_format("log10", math_format(10^.x))) +
      theme_bw() +
      scale_color_manual(values=palette.complex) + 
      scale_shape_manual(values=shape.ggplot)+
      scale_size_manual(values=size.ggplot)+
      scale_fill_manual(values=palette.complex) +
      geom_hline(yintercept = input$yintercept,colour=input$intercept.col,linetype=input$int.type)+
      geom_vline(xintercept = input$xintercept,colour=input$intercept.col, linetype=input$int.type)+
      annotation_logticks()  +
      theme(text=element_text(size=20,family=input$font_type2),
            plot.margin = margin(20, 20, 20, 20),
            axis.text.x = element_text(colour="black",size=input$axis.numeric.size,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type2),
            axis.text.y = element_text(colour="black",size=input$axis.numeric.size,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type2),
            axis.title.x=element_text(colour="black",size=input$axis.title.size,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type2),
            axis.title.y = element_text(colour="black",size=input$axis.title.size,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type2),
            legend.title  =element_blank(),
            legend.position = input$legend.dot,
            legend.text = element_text(colour="black",size=input$legend.size.cd,hjust=.5,vjust=.5,face="plain",family=input$font_type2)) +
      scale_alpha(guide = 'none') +
      guides(size=FALSE, col = guide_legend(ncol=input$legend.column,override.aes = list(size=input$leg.dot.size)))+
      labs(x=x_lable1,
           y=y_lable1)
      
    
    vals15$complex_dot
  })
  
  dot_plot.complex_UMAP <- reactive({
    set.seed(123)
    index <- input.data_CSV2();
    validate(
      need(nrow(index)>0,
           "Upload file")
    )
    names(index) <- gsub("\\.", " ", names(index))
    
    
    index <- as.data.frame(index)
    y_lable1 <- bquote(.(input$y.axis2))
    x_lable1 <-  bquote(.(input$x.axis2))
    
    index[is.na(index)] <- "not_clonal"
    selected.col <- index[names(index) %in% input$group_complex_dot]
    names(selected.col) <- "V1"
    index[names(index) %in% input$group_complex_dot] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))
    palette.complex <- unlist(colors.FACS.index())
    shape.ggplot <- unlist(shape.FACS.index2())
    size.ggplot <- unlist(size.FACS.index2())
    
    names_unique <- as.data.frame(unique(selected.col))
    names_unique_size <- cbind(names_unique,as.data.frame(size.ggplot))
    
    names(names_unique_size)[1] <- input$group_complex_dot
    
    index2 <- merge(index,names_unique_size,by=input$group_complex_dot)
    
    vals15$complex_dot <- ggplot(index2, aes(x=get(input$x.axis2), y=get(input$y.axis2),
                                             colour = get(input$group_complex_dot),
                                             fill = get(input$group_complex_dot),
                                             shape = get(input$group_complex_dot), 
                                             size = get(input$group_complex_dot), 
    ))+
      geom_point() + 
      theme_bw() +
      scale_color_manual(values=palette.complex) + 
      scale_shape_manual(values=shape.ggplot)+
      scale_size_manual(values=size.ggplot)+
      scale_fill_manual(values=palette.complex) +
      theme(text=element_text(size=20,family=input$font_type2),
            plot.margin = margin(20, 20, 20, 20),
            axis.text.x = element_text(colour="black",size=input$axis.numeric.size,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type2),
            axis.text.y = element_text(colour="black",size=input$axis.numeric.size,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type2),
            axis.title.x=element_text(colour="black",size=input$axis.title.size,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type2),
            axis.title.y = element_text(colour="black",size=input$axis.title.size,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font_type2),
            legend.title  =element_blank(),
            legend.position = input$legend.dot,
            legend.text = element_text(colour="black",size=input$legend.size.cd,hjust=.5,vjust=.5,face="plain",family=input$font_type2)) +
      scale_alpha(guide = 'none') +
      guides(size=FALSE, col = guide_legend(ncol=input$legend.column,override.aes = list(size=input$leg.dot.size)))+
      labs(x=x_lable1,
           y=y_lable1)
    
    vals15$complex_dot
    
    
  })
  
  output$Add_size <- DT::renderDataTable( {
      
    set.seed(123)
    index <- input.data_CSV2();
    validate(
      need(nrow(index)>0,
           "Upload file")
    )
    names(index) <- gsub("\\.", " ", names(index))
    
    
    index <- as.data.frame(index)
    y_lable1 <- bquote(.(input$y.axis2))
    x_lable1 <-  bquote(.(input$x.axis2))
    
    index[is.na(index)] <- "not_clonal"
    selected.col <- index[names(index) %in% input$group_complex_dot]
    names(selected.col) <- "V1"
    index[names(index) %in% input$group_complex_dot] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))
    palette.complex <- unlist(colors.FACS.index())
    shape.ggplot <- unlist(shape.FACS.index2())
    size.ggplot <- unlist(size.FACS.index2())
    
    names_unique <- as.data.frame(unique(selected.col))
    names_unique_size <- cbind(names_unique,as.data.frame(size.ggplot))
    
    names(names_unique_size)[1] <- input$group_complex_dot
    
    index2 <- merge(index,names_unique_size,by=input$group_complex_dot)
    
    # rownames(size.ggplot) <- unique(selected.col$V1)

    datatable(as.data.frame(index2), extensions = "Buttons", options = list(searching = TRUE,
                                                                           ordering = TRUE,
                                                                           buttons = c('copy','csv', 'excel'),
                                                                           dom = 'Bfrtip',
                                                                           pageLength=2, 
                                                                           lengthMenu=c(2,5,10,20,50,100), 
                                                                           scrollX = TRUE
      ))
    }, server = FALSE)
  
  dot_plot.complex1 <- function () {
    
    if(input$plot_type_umap == "overlaid dot plot") {
      if (input$density_dotplot =="no" & input$grid.lines.dot =='no') {
        
        if (input$Add_ellipse =="no") {
          dot_plot.complex() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }
       else {
         dot_plot.complex() +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
           stat_ellipse(geom="polygon",level=0.8,alpha=0.2,lwd = 0.8)
       }
      }
      else if (input$density_dotplot =="no" & input$grid.lines.dot =='yes') {
        if (input$Add_ellipse =="no") {
          dot_plot.complex() 
        }
        else {
          dot_plot.complex() +
            stat_ellipse(geom="polygon",level=0.8,alpha=0.2,lwd = 0.8)
          
        }
        
        
      }
      else if (input$density_dotplot =="yes" & input$grid.lines.dot =='no') {
        
        if (input$Add_ellipse =="no") {
          dot_plot <- dot_plot.complex() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
          ggExtra::ggMarginal(dot_plot,groupColour = TRUE, groupFill = TRUE)
        }
        else {
          dot_plot <- dot_plot.complex() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + stat_ellipse(geom="polygon",level=0.8,alpha=0.2,lwd = 0.8)
          ggExtra::ggMarginal(dot_plot,groupColour = TRUE, groupFill = TRUE)

        }
      }
      else {
        if (input$Add_ellipse =="no") {
          ggExtra::ggMarginal(dot_plot.complex(),groupColour = TRUE, groupFill = TRUE)
        }
        else {
          dot_plot <- dot_plot.complex() +
            stat_ellipse(geom="polygon",level=0.8,alpha=0.2,lwd = 0.8)
          
          ggExtra::ggMarginal(dot_plot,groupColour = TRUE, groupFill = TRUE)
          
        }
        
        
        
      }
    }
    else if (input$plot_type_umap == "Ridge plot") {
      dot_plot_ridge_df()
      
    }
    else{
      
      if (input$Add_ellipse =="no" & input$grid.lines.dot =='no') {
        
        dot_plot.complex_UMAP() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
      else if (input$Add_ellipse =="no" & input$grid.lines.dot =='yes') {
        
        dot_plot.complex_UMAP() 
      }
      else if (input$Add_ellipse =="yes" & input$grid.lines.dot =='no') {
        dot_plot.complex_UMAP() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
          stat_ellipse(geom="polygon",level=0.8,alpha=0.2,lwd = 0.8)
        
      }
      else {
        dot_plot.complex_UMAP() +
          stat_ellipse(geom="polygon",level=0.8,alpha=0.2,lwd = 0.8)
      }
      
      
    }
    
    
    
  }
  output$dot_plot.complex2 <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    
    dot_plot.complex1()
  })
  
  
  dot_plot_ridge_df  <- reactive({
    set.seed(123)
    index <- input.data_CSV2();
    validate(
      need(nrow(index)>0,
           "Upload file")
    )
    names(index) <- gsub("\\.", " ", names(index))
    
    title_axis <- bquote(.(input$x.axis2))
    index <- as.data.frame(index)
    index[is.na(index)] <- "not_clonal"
    selected.col <- index[names(index) %in% input$group_complex_dot]
    names(selected.col) <- "V1"
    index[names(index) %in% input$group_complex_dot] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))
    palette.complex <- unlist(colors.FACS.index())
    index$selected <- index[,names(index) %in% input$x.axis2]
    index_clone <- index[,names(index) %in% c("cloneCount",input$group_complex_dot)]
    index_clone <- index_clone %>%
      select(cloneCount, everything())
    df3 <- as.data.frame(ddply(index_clone,names(index_clone[-c(1)]),numcolwise(sum)))
    names(df3)[2] <- "Sum_count"
    df4 <- merge(index,df3,by=input$group_complex_dot)
    index2 <- subset(df4,df4$Sum_count>2)
    num <- unique(index[names(index) %in% input$group_complex_dot])
    num$palette.complex <- palette.complex
    
    df <- as.data.frame(unique(index2[,names(index2) %in% c("cloneCount",input$group_complex_dot)]))
    # names(df) <- input$group_complex_dot
    df.col1 <- merge(df,num,by=input$group_complex_dot)
    x_lable1 <- bquote(Log[10]~(.(input$x.axis2)))
    x_lable1
    
   ggplot(index, aes(x = log10(selected), y = as.character(get(input$group_complex_dot)), fill = as.character(get(input$group_complex_dot)))) +
      geom_density_ridges() +
      theme_ridges() + 
      scale_fill_manual(values=df.col1$palette.complex)+
      theme(legend.position = "none") +
      labs(x=x_lable1) +
     theme(text=element_text(size=20,family=input$font_type2),
           plot.margin = margin(20, 20, 20, 20),
           axis.text.x = element_text(colour="black",size=input$axis.numeric.size,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type2),
           axis.text.y = element_text(colour="black",size=input$axis.numeric.size,angle=0,hjust=1,vjust=0,face="plain",family=input$font_type2),
           axis.title.x = element_text(colour="black",size=input$axis.title.size,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font_type2),
           axis.title.y = element_blank(),
           # panel.grid.major.x = element_blank() ,
           
           )
           # legend.title  =element_blank(),
           
    
  })

  output$dot_plot_ridge_tab <- DT::renderDataTable( {
    index <- input.data_CSV2();
    validate(
      need(nrow(index)>0,
           "Upload file")
    )
    names(index) <- gsub("\\.", " ", names(index))
    
    index <- as.data.frame(index)
    index[is.na(index)] <- "not_clonal"
    
    index$selected <- log10(index[,names(index) %in% input$x.axis2])
    
    tab <-as.data.frame(TukeyHSD( aov(selected ~ as.character(get(input$group_complex_dot)),data = index))[1])
    names(tab) <- c("diff" ,"lwr","upr","p.adj")

    tab$stat <- ifelse(tab$p.adj<0.0001,"****",
                       ifelse(tab$p.adj<0.001,"***",
                              ifelse(tab$p.adj<0.01,"**",
                                     ifelse(tab$p.adj<0.05,"*","NS"
                                     ))))
    

    
    datatable(tab, extensions = "Buttons", options = list(searching = TRUE,
                                                                     ordering = TRUE,
                                                                     buttons = c('copy','csv', 'excel'),
                                                                     dom = 'Bfrtip',
                                                                     pageLength=10, 
                                                                     lengthMenu=c(2,5,10,20,50,100), 
                                                                     scrollX = TRUE
    ))
  }, server = FALSE)
  
  output$dot_plot_ridge <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    
    dot_plot_ridge_df()
  })
  
  output$downloadPlot_complex.dotplot <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste(input$name.colour3,"complex.dotplot_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_complex.dotplot,height=input$height_complex.dotplot, onefile = FALSE) # open the pdf device
      print(dot_plot.complex1())
      dev.off()}, 
    contentType = "application/pdf" )
  
  output$downloadPlotPNG_complex.dotplot <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste(input$name.colour3,"complex.dotplot_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_complex.dotplot, height = input$height_png_complex.dotplot, res = input$resolution_PNG_complex.dotplot)
      
      
      print(dot_plot.complex1())
      
      
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  # upset diagram -----
  vals23 <- reactiveValues(upset.plot=NULL)
  observe({
    updateSelectInput(
      session,
      "upset.select",
      choices=names(input.data_old2()),
      selected = names(input.data_old2())[dim(input.data_old2())[2]] 
      )
    
  })
  observe({
    updateSelectInput(
      session,
      "upset.group.select",
      choices=names(select_group()),
      selected = names(select_group())
      )
    
  })
  observe({
    updateSelectInput(
      session,
      "order.of.group",
      choices=select_group(),
      selected = c("CD8","IFN"))
  })
  # df <-df[df$group %in% input$string.data2,]
  upset.parameters <- function () {
    df <- input.data_old2();
    
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    
    
    df <- as.data.frame(df)
    head(df)
    
    
    
    unique.df <- unique(df[c(input$upset.select,input$upset.group.select)])
    names(unique.df) <- c("chain","group")
    unique.df <-unique.df[unique.df$group %in% input$order.of.group,]
    
    
    unique.df$cloneCount <- 1
    mat <- acast(unique.df, chain~group, value.var="cloneCount")
    mat[is.na(mat)] <- 0
    mat <- as.data.frame(mat)
    # a <- as.data.frame(unique(names(mat)))
    #   a$V1 <- distinctColorPalette(dim(a)[1])
    
    df.x <- make_comb_mat(mat)
    
    
    if (input$upset_anno == "Colour by degree") {
      
      your_list <- c(input$upset_colours_list1)
      your_list_df <- as.data.frame((unlist(strsplit(your_list, ','))))
      names(your_list_df) <- "ID"
      
      ht = draw(UpSet(df.x,
                      pt_size = unit(input$upset.point.size, "mm"),
                      lwd = input$upset.lwd,
                      row_names_gp =  gpar(fontfamily = input$font_type,fontsize = input$upset.text.size),#changes font size of "set size" labels
                      column_names_gp = gpar(fontfamily = input$font_type),
                      comb_col = c(your_list_df$ID)[comb_degree(df.x)],
                      top_annotation = upset_top_annotation(df.x,
                                                            numbers_gp = gpar(fontfamily = input$font_type,fontsize = input$font.size.anno.upset),
                                                            annotation_name_gp = gpar(fontfamily = input$font_type,fontsize=input$font.size.anno.upset),
                                                            gp = gpar(fill = input$top_annotation_colour),
                                                            
                      ),
                      right_annotation = upset_right_annotation(df.x,
                                                                add_numbers = T,
                                                                numbers_gp = gpar(fontfamily = input$font_type,fontsize = input$font.size.anno.upset),
                                                                annotation_name_gp = gpar(fontfamily = input$font_type,fontsize=input$font.size.anno.upset),
                                                                gp = gpar(fill = input$right_annotation_colour),
                      ),
                      
                      
                      set_order  = c(input$order.of.group)
      ), padding = unit(c(20, 20, 20, 20), "mm"))
      
      
    }
    
    else {
      
      ht = draw(UpSet(df.x,
                      pt_size = unit(input$upset.point.size, "mm"),
                      lwd = input$upset.lwd,
                      row_names_gp =  gpar(fontfamily = input$font_type, fontsize = input$upset.text.size),
                      column_names_gp = gpar(fontfamily = input$font_type),
                      top_annotation = upset_top_annotation(df.x,
                                                            annotation_name_gp = gpar(fontfamily = input$font_type),
                                                            gp = gpar(fill = input$top_annotation_colour),
                      ),
                      
                      
                      right_annotation = upset_right_annotation(df.x,
                                                                add_numbers = T,
                                                                numbers_gp = gpar(fontfamily = input$font_type,fontsize = input$font.size.anno.upset),
                                                                annotation_name_gp = gpar(fontfamily = input$font_type,fontsize=input$font.size.anno.upset),
                                                                gp = gpar(fill = input$right_annotation_colour),
                      ),
                      
                      
                      set_order  = c(input$order.of.group)
                      
      ), padding = unit(c(20, 20, 20, 20), "mm"))
    }
    

    od = column_order(ht)
    cs = comb_size(df.x)
    decorate_annotation("intersection_size", {
      grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
                default.units = "native", just = "bottom", gp = gpar(fontsize = input$upset.font.size, fontfamily = input$font_type)
      ) })
    
  }
  output$UpSet.plot <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    upset.parameters()
  })
  upset.datatable1 <- reactive({
    df <- input.data2();
    
    validate(
      need(nrow(df)>0,
           error_message_val1)
    )
    
    df <- as.data.frame(df)
    head(df)
    unique.df <- unique(df[c(input$upset.select,input$upset.group.select)])
    names(unique.df) <- c("chain","group")
    unique.df$cloneCount <- 1
    mat <- acast(unique.df, chain~group, value.var="cloneCount")
    mat[is.na(mat)] <- 0
    mat <- as.data.frame(mat)
    mat
  })
  
  output$upset.datatable <- DT::renderDataTable( {
    datatable(upset.datatable1(), extensions = "Buttons", options = list(searching = TRUE,
                                                                         ordering = TRUE,
                                                                         buttons = c('copy','csv', 'excel'),
                                                                         dom = 'Bfrtip',
                                                                         pageLength=2, 
                                                                         lengthMenu=c(2,5,10,20,50,100), 
                                                                         scrollX = TRUE
    ))
  }, server = FALSE)
  
  output$downloadPlot_upset <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("upset_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_upset,height=input$height_upset, onefile = FALSE) # open the pdf device
      print(upset.parameters())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_upset <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("upset_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_upset, height = input$height_png_upset, res = input$resolution_PNG_upset)
      print(upset.parameters())
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
}

shinyApp(ui, server)
