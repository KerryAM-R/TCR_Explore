
index <- c("shannon", "simpson","invsimpson")

draw_colnames_rotate <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, hjust = .5, rot = 0, gp = gpar(...)) # rot = rotation for # degrees
  return(res)}
gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Nucleotide <- function (Nucleotide, seqlength) 
{
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
require("treemap")
require("treemapify")
require("circlize")
require("motifStack")
require("VLF")
require("gplots")
require("MASS") # to access Animals data sets
require("scales") # to access break formatting functions
require("flowCore")
library("readxl")
require("RColorBrewer")
require("vegan")
require("ggheatmap")
library("randomcoloR")
library("colourpicker")
require("ComplexHeatmap")

test_fun <- function()
{
  for (i in 1:15) {
    incProgress(1/15)
    sum(runif(1000000,0,1))
  }
}


graph_type <- c("histogram","density")
axis_density_group <- c("x-axis","y-axis")

angle <- c(0,90,180,270)

error_message_val1 <- "No data found"
error_message_val2 <- "Suggest uploading file"
error_message_val3 <- "No data found\n \nSuggest uploading file\nheaders=ID, logFC, Pvalue"
error_message_val4 <- "no own list found\n \nSuggest uploading file\nheaders=ID"
# user interface  ----
ui <- navbarPage("TCR_explore", position = "fixed-top",collapsible = TRUE,
                 
                 
                 navbarMenu("Workflow",
                 tabPanel("Overview",
                          fluidRow(includeMarkdown("README.md"))
                            ),     
                 tabPanel("QC",
                          fluidRow(includeMarkdown("READMEQC.md"))
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
                           
                 # UI IMGT only ----
                            tabPanel("IMGT",
                                     sidebarLayout(
                                       sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                                    selectInput("dataset_IMGT3", "Choose a dataset:", choices = c("ab-test-data1","gd-test-data1", "own_data")),
                                                    fileInput('file_IMGT3', 'Select file for IMGT datafile',
                                                              accept=c('xls/xlsx', '.xls')),
                                                    selectInput("IMGT_chain2","alpha-beta or gamma-delta",choices = c("ab","gd")),
                                                    selectInput("dataset_IMGT_afterQC", "Choose a dataset:", choices = c("ab-test-data1","gd-test-data1", "own1")),
                                                    
                                                    fileInput('file_IMGT_afterQC', 'Select chromatogram checked file (.csv)',
                                                             accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
                                       ),
                                       mainPanel(
                                         tabsetPanel(
                                             tabPanel("IMGT check uploaded",
                                                      tags$head(tags$style("#IMGT2_out  {white-space: nowrap;  }")),
                                             div(DT::dataTableOutput("IMGT2_out")),
                                             downloadButton('downloadTABLE_IMGTonly','Download table')
                                             ),
                                             tabPanel("QC file unsummarised",
                                                      tags$head(tags$style("#chain_table_IMGT.QC1  {white-space: nowrap;  }")),
                                                      div(DT::dataTableOutput("chain_table_IMGT.QC1")),
                                                      downloadButton('downloadTABLE.QC1','Download table')
                                             ), 
                                             tabPanel("Summary table",
                                                      verbatimTextOutput("names.in.file3"),
                                                      fluidRow(column(12, textInput("string.data3","column names for summary","group, Indiv, V.GENE.and.allele_A, J.GENE.and.allele_A, JUNCTION..AA._A, V.GENE.and.allele_B, J.GENE.and.allele_B, D.GENE.and.allele_B, JUNCTION..AA._B", width = "1200px") )),
                                                      tags$head(tags$style("#chain_table_IMGT.QC3  {white-space: nowrap;  }")),
                                                      div(DT::dataTableOutput("chain_table_IMGT.QC3")),
                                                      downloadButton('downloadTABLE.QC3','Download table')
                                                      
                                             ) 
                                         
                                       )
                                       )
                                     
                                     )
                                     ),
                            
                            
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
                 
                 tabPanel("scTCR analysis",
                          sidebarLayout(
                            sidebarPanel(id = "tPanel",style = "overflow-y:scroll; max-height: 700px; position:relative;", width=3,
                                         tags$style(type="text/css", "body {padding-top: 70px; padding-left: 10px;}"),
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
                                         
                                         selectInput("group_column",label = h5("column of group"), ""),
                                         
                                         tags$hr(),
                                         selectInput("shannon.index", "Choose a dataset:", choices = c("gd.test.index", "own.index")),
                                         fileInput('file_diversity.index', 'diversity statistics file',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                         ),
                            
                            
                            mainPanel(tabsetPanel(
                 # UI Treemap -----
                              
                              tabPanel("Treemap",
                                                                              
                                         fluidRow( 
                                           column(3, selectInput("tree_colour.choise",label = h5("colour default choises"), choices =  c("default","random","grey"))),
                                           column(3, selectInput( "fill2",label = h5("Colour treemap by"),"" )),
                                           column(3, selectInput( "sub_group2",label = h5("separate panels by"),"" )),
                                           column(3,selectInput( "wrap",label = h5("Group"),"sub" )),
                                           column(3,selectInput( "count2",label = h5("Count column"),"")),
                                           column(3, selectInput("tree.lab",label = h5 ("add label"),choices = c("yes","no")))
                                           
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
                              tabPanel("Circular plot", 
                                       h5("Circular plot requires an unsummarised dataset"),
                                       
                                       fluidRow(

                                         column(2,selectInput( "group_selected2",label = h5("group"),"" )),
                                         column(2,selectInput( "chain1",label = h5("chain one"),"" )),
                                         column(2,selectInput( "chain2",label = h5("chain two"),"" )),
                                        
                                         column(2,tableOutput("table_display")),
                                         
                                       ),
                                       fluidRow(
                                          column(3, selectInput("circ.lab",label = h5("add label"),choices = c("yes","no"))),
                                          column(3,selectInput( "colour_cir",label = h5("colour"),choices = c("rainbow","random","grey"))),    
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
                                       
                              ),
                 # UI bar graphs ----- 
                              tabPanel("Chain bar graph",
                                       h5("Chain plot requires an unsummarised dataset"),
                                       fluidRow(
                                         column(2,selectInput( "variable_chain",label = h5("y-axis"),"" )),
                                         column(2,selectInput( "graph_bar_type",label = h5("x-axis"),choices = c("count","percentage"))),
                                         column(2,selectInput( "selected_group_chain",label = h5("group"),"" )),
                                         column(2,style = "margin-top: 10px;", textInput("colour_bar.usage","colour of graph", value = "black"))

                                       ),
                                       plotOutput("Chain1_usage",height="800px"),
                                       
                                       
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
                 # UI CDR3 length distribution graphs ----- 
                              tabPanel("CDR3 length distribution",
                                       h5("length distribution plot requires an unsummarised dataset"),
                                       fluidRow(column(4, h6("please select V1 with the default data")),
                                                column(4,selectInput( "aa.or.nt",label = h5("amino acid or nucleotide sequence column"),"" )),
                                                column(4,selectInput( "selected_group_len",label = h5("group"),"" ))
                                       ),
                                       fluidRow(
                                         column(4, numericInput("bin","binds of histogram",value=30)),
                                         column(4, textInput("hist_col","colour of histogram","red")),
                                         column(4, selectInput('graph_type', 'Type of graph', graph_type))
                                       ),
                                       plotOutput("Chain1_length"),
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
                                         column(3,selectInput( "aa.or.nt2",label = h5("amino acid CDR3 column"),"" )),
                                         column(3, numericInput("len","CD3R amino acid length", value = 15)),                               
                                         column(3,selectInput( "group_selected_motif",label = h5("group"),"" ))
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
                                         column(3,selectInput( "aa.or.nt3",label = h5("nucleotide CDR3 column"),"")),
                                         column(3,selectInput( "group_selected",label = h5("group"),"" )),
                                         column(3, numericInput("len_nt","CD3R nucleotide length", value = 30))
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
                 # UI Pie ----
                              tabPanel("pie chart of clonal distribution",
                                       fluidRow(column(3,selectInput("pie_chain",label = h5("Colour by this group"),"")),
                                                column(3,selectInput("pie_colour.choise",label = h5("colour default choises"), choices =  c("default","random","grey")))),
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
                              
                 # UI heatmap -----
                              tabPanel("ggplot heatmap groups",
                                       selectInput("group_hm", "Select specific groups", choices = c("yes", "no")),
                                       fluidRow(
                                         column(3,selectInput("group_selected3",label = h5("select group"),"" )),
                                         column(3,selectInput( "heatmap_2",label = h5("x-axis"),"" )),
                                         column(3,selectInput("group.heatmap",label = h5("y-axis"),"" ))),
                                                                
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
          
                 # UI inverse simpson index -----
                              tabPanel("inverse Simpson Index",
                                       p("Inverse simpson index; ∞=infinite diversity and 1=limited diversity"),
                                       
                                       fluidRow(column(9, div(DT::dataTableOutput("table_display.diversity")))),
                                       downloadButton('downloadTABLE_simpson.inv','Download table'),
                                       
                                       fluidRow(
                                                column(3,selectInput("inv.simp_colour.choise",label = h5("colour default choises"), choices =  c("default","random","grey"))),
                                                column(3,selectInput("group.index",label = h5("Group by on x-axis"),"")),
                                                column(3,selectInput("group2.index",label = h5("Fill by this group"),"")),
                                                column(3,selectInput("x.axis.index",label = h5("x-axis (total or unique clones"),""))),
                                       
                                       fluidRow(
                                                column(3,
                                                wellPanel(id = "tPanel22",style = "overflow-y:scroll; max-height: 400px",
                                                          uiOutput('myPanel.inv.simp'))),
                                         
                                         
                                                column(4,plotOutput("simpson.index1", height="400px")),
                                                column(4,plotOutput("simpson.index2", height="400px"))),
                                       
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
                                       ),
                                       
                                       # gini index is created from Lorentz Surface Calculation pone.0125373.s004.xlsx
                                       
                              ),
                 # upset plot -----
                 tabPanel("upset plot",
                          fluidRow(
                            column(3,selectInput("upset.select",label = h5("select chain"), choices = "", selected = "")),
                            column(3,selectInput("upset.group.select",label = h5("select group column (max 31 groups)"), choices = "",selected= ""))),
                          
                          tags$head(tags$style("#upset.datatable  {white-space: nowrap;  }")),
                          div(DT::dataTableOutput("upset.datatable")),
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
                          )
                          
                          
                          )

                                  
                            )
                            )
                          )
                          
                 ),
                 
                 
                 # UI Index data graphs -----


                 
                 tabPanel("FACS Index data",
                          sidebarLayout(
                            sidebarPanel(id = "tPanel3",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                         
                                         selectInput("dataset3", "FACS file:", choices = c("test-FACS", "own_FACS")),
                                         fileInput('file_FACS', 'Raw index FACS file',
                                                   accept=c('FACS files', '.fcs')),
                                         
                                         selectInput("data_clone.index", "Unsummarised clone file source:", choices = c("gd.test.clone" ,"own.clone.file")),
                                         fileInput('file_diversity.index.2', 'upload unsummarised clone file',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),

                                         selectInput("dataset7", "Merged FACS and clone file (make colours", choices = c("test-csv" ,"own_csv")),
                                         fileInput('file_FACS.csv1', 'several groups .csv file',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                         selectInput("dataset_index.2", "Choose a dataset for complex plot:", choices = c("test-csv" ,"own_csv")),
                                         fileInput('file_FACS.csv2', 'several groups .csv file',
                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),

                            ),
                            mainPanel(tabsetPanel(
                              tabPanel("FACS to .csv",
                                       div(DT::dataTableOutput("FACS.CSV")),
                              ),
                              
                              tabPanel("merging scTCR with Index data",
                              div(DT::dataTableOutput("merged.clone")),
                              div(DT::dataTableOutput("merged.index.clone")),
                              
                              
                              
                              
                              downloadButton('downloadTABLE_FACS','Download table')
                              
                              
                              
                              ),
  
                 # UI complex dotplot add columns if needed -----
                              tabPanel("Adding columns to complex datasheet for colouring",
                              
                                verbatimTextOutput("names.in.file"),
                                textInput("string.data","column names for selecting clonal","TRDV, dCDR3, TRGV, gCDR3, group"),
                                fluidRow(
                                  column(2,selectInput("V.gene.1",label = h5("V Gene 1"),"")),
                                  column(2,selectInput("CDR3.1",label = h5("CDR3 1"),"")),
                                  column(2,selectInput('V.gene.2', label = h5("V Gene 2"), "")),
                                  column(2,selectInput("CDR3.2",label = h5("CDR3 1"),"")),
                                  column(2,numericInput("numeric.cloneCount","select clones greater than",value=1))
                                ),
                                div(DT::dataTableOutput("table.index.1")),
                                verbatimTextOutput("NAMES.df"),
                                downloadButton('downloadTABLE_cleaning','Download table')),
                 # UI complex dotplot -----
                                tabPanel("complex dotplot",
                                         fluidRow(
                                           column(3,selectInput("x.axis2",label = h5("x-axis label"),"")),
                                           column(3,selectInput("y.axis2",label = h5("y-axis label"),"")),
                                           column(3,selectInput("density_dotplot",label = h5("add histogram"), choices = c("no","yes"))),
                                           column(3,selectInput("group_complex_dot",label = h5("create col column"),"")),
                                         ),
                                         fluidRow(
                                           column(3,numericInput("legend.size.cd","legend text size",value=12)),
                                           column(3, selectInput("legend.dot",label=h5("legend position"),choices = c("top","bottom","left","right","none"),selected = "right")),
                                           column(4,numericInput("legend.column", "# of legend columns", value=1)),
                                           column(2,selectInput( "FACS.index_colour.choise",label = h5("colour"),choices = c("default","random","grey"))),
                                           
                                           
                                         ),
                                         
                                         fluidRow(column(3,
                                                         wellPanel(id = "tPanel222",style = "overflow-y:scroll; max-height: 600px",
                                                                   uiOutput('myPanel.FACS.index'))),
                                                  column(9, plotOutput("dot_plot.complex2",height = "600px"))),

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
1200/8*6  
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
  
  # QC MiXCR and IMGT files -----
  input.data_IMGT.xls <- reactive({switch(input$dataset_IMGT,"ab-test-data" = test.data_ab.xls(), "gd-test-data" = test.data_gd.xls(),"own" = own.data.IMGT())})
  test.data_ab.xls <- reactive({
    dataframe = read_excel("test-data/IMGT/T00016.xls") 
  })
  test.data_gd.xls <- reactive({
    dataframe = read_excel("test-data/IMGT/T00016.xls") 
  })
  own.data.IMGT <- reactive({
    inFile4 <- input$file_IMGT1
    if (is.null(inFile4)) return(NULL)
    
    else {
      dataframe <- read_excel(
        inFile4$datapath
        
      )}
    
  })
  
  input.data_MiXCR.txt <- reactive({switch(input$dataset_IMGT,"ab-test-data" = test.data_ab.txt(),"gd-test-data" = test.data_gd.txt(),"own" = own.data.MiXCR())})
  test.data_ab.txt <- reactive({
    dataframe = read.table("test-data/MiXCR/SJS_TEN/T00016.txt",sep="\t",header=T) 
  })
  test.data_gd.txt <- reactive({
    dataframe = read.table("test-data/MiXCR/SJS_TEN/T00016.txt",sep="\t",header=T) 
  })
  own.data.MiXCR <- reactive({
    inFile5 <- input$file_IMGT2
    if (is.null(inFile5)) return(NULL)
    
    else {
      dataframe <- read.table(
        inFile5$datapath,
        header=inFile5$header_IMGT2,
        sep=input$sep_IMGT2,
        quote=input$quote_IMGT2)}
    
  })
  
  # combining MiXCR and IMGT -----
  input.data.IMGT.MiXCR <- reactive({switch(input$dataset_IMGT,"ab-test-data" = test.data.ab.csv(), "gd-test-data" = test.data.gd.csv(),"own" = own.data.csv())})
  test.data.ab.csv <- reactive({
    dataframe = read.csv("test-data/QC/merged_IMGT_table2021.08.01.csv",header=T) 
  })
  test.data.gd.csv <- reactive({
    dataframe = read.csv("test-data/QC/merged_IMGT_table2021.08.01.csv",header=T) 
  })
  own.data.csv <- reactive({
    inFile6 <- input$file_IMGT.MiXCR
    if (is.null(inFile6)) return(NULL)
    
    else {
      dataframe <- read.csv(
        inFile6$datapath)}
    
  })
  
  # unmerged IMGT df
  
  output$unmerged_chain1 <- DT::renderDataTable(escape = FALSE,options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
    df <- input.data_IMGT.xls();
    
    
    
    df_chain1 <- as.data.frame(df)
    df_chain1$`J-GENE and allele` <- gsub('Homsap ','',df_chain1$`J-GENE and allele`)
    df_chain1$`V-GENE and allele` <- gsub('Homsap ','',df_chain1$`V-GENE and allele`)
    df_chain1$`D-GENE and allele` <- gsub('Homsap ','',df_chain1$`D-GENE and allele`)
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
    df_chain1$`J-GENE and allele` <- gsub('TR','',df_chain1$`J-GENE and allele`)
    df_chain1$`V-GENE and allele` <- gsub('TR','',df_chain1$`V-GENE and allele`)
    df_chain1$`D-GENE and allele` <- gsub('TR','',df_chain1$`D-GENE and allele`)
    df_chain1$JUNCTION <- toupper(df_chain1$JUNCTION) 
    df_chain1
    
  })
  output$unmerged_chain2 <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
    df <- input.data_MiXCR.txt();
    df <- as.data.frame(df)
    df
    
  })
  IMGT <- reactive({
    df_chain1 <- input.data_IMGT.xls();
    df_chain1 <- as.data.frame(df_chain1)
    df_chain1$`J-GENE and allele` <- gsub('Homsap ','',df_chain1$`J-GENE and allele`)
    df_chain1$`V-GENE and allele` <- gsub('Homsap ','',df_chain1$`V-GENE and allele`)
    df_chain1$`D-GENE and allele` <- gsub('Homsap ','',df_chain1$`D-GENE and allele`)
    df_chain1$`J-GENE and allele` <- gsub('[(]','',df_chain1$`J-GENE and allele`)
    df_chain1$`J-GENE and allele` <- gsub('[)]','',df_chain1$`J-GENE and allele`)
    df_chain1$`J-GENE and allele` <- gsub(' F','',df_chain1$`J-GENE and allele`)
    df_chain1$`J-GENE and allele` <- gsub(' or ',', ',df_chain1$`J-GENE and allele`)
    
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
    

    df_chain1$JUNCTION <- toupper(df_chain1$JUNCTION) 
    
    # MiXCR -----
    df_chain2 <- input.data_MiXCR.txt();
    df_chain2 <- as.data.frame(df_chain2)
    df <- merge(df_chain1,df_chain2, by.x="Sequence ID",by.y="ID",all=T)
    df$V.Gene.IMGT <- df$`V-GENE and allele`
    df$V.Gene.IMGT <- gsub("[*]0.","",df$V.Gene.IMGT)
    df$V.Gene.IMGT <- gsub("[*]0..","",df$V.Gene.IMGT)
    df$V.Gene.IMGT <- gsub("/","",df$V.Gene.IMGT)
    df$V.Gene.IMGT <- gsub(",.*", "", df$V.Gene.IMGT)
    df$V.Gene.IMGT <- gsub(" ", "", df$V.Gene.IMGT)
    df$V.Gene.IMGT <- gsub("TR", "", df$V.Gene.IMGT)
    df$`J-GENE and allele` <- gsub('TR','',df$`J-GENE and allele`)
    df$`V-GENE and allele` <- gsub('TR','',df$`V-GENE and allele`)
    df$`D-GENE and allele` <- gsub('TR','',df$`D-GENE and allele`)
    
    df <- mutate(df, CDR3_sequence_match_nt = ifelse(df$JUNCTION==df$nSeqCDR3,"matched","different"),
                 CDR3_sequence_match_aa = ifelse(df$`JUNCTION (AA)`==df$aaSeqCDR3,"matched","different"),
                 Matched_TRV = ifelse(df$V.Gene.IMGT==df$bestVHit,"matched","different"))
    
    df$clone_quality <- NA 
    df$comments <- NA
    df
    
  })
  output$merged_IMGT <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    IMGT()
  })
  output$downloadTABLE2 <- downloadHandler(
    filename = function(){
      paste("merged_IMGT_table",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      write.csv(IMGT(),file, row.names = FALSE)
    })
  # IMGT only  -----
  input.data_IMGT.xls3 <- reactive({switch(input$dataset_IMGT3,"ab-test-data1" = test.data_ab.xls3(), "gd-test-data1" = test.data_gd.xls3(),"own_data" = own.data.IMGT3())})
  test.data_ab.xls3 <- reactive({
    dataframe = read_excel("Raw_data/vquest-2.xls") 
  })
  test.data_gd.xls3 <- reactive({
    dataframe = read_excel("Raw_data/vquest-2.xls") 
  })
  own.data.IMGT3 <- reactive({
    inFile_IMGT3 <- input$file_IMGT3
    if (is.null(inFile_IMGT3)) return(NULL)
    
    else {
      dataframe <- read_excel(
        inFile_IMGT3$datapath
        
      )}
    
  })
  
  
  input.data_IMGT.xls4 <- reactive({switch(input$dataset_IMGT3,"ab-test-data1" = test.data_ab.xls4(), "gd-test-data1" = test.data_gd.xls4(),"own_data" = own.data.IMGT4())})
  test.data_ab.xls4 <- reactive({
    dataframe = read_xls("Raw_data/vquest-2.xls",sheet = 2) 
  })
  test.data_gd.xls4 <- reactive({
    dataframe = read_xls("Raw_data/vquest-2.xls",sheet = 2) 
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
    df2 <- input.data_IMGT.xls4();
    
    df3 <- df1[names(df1) %in% c("Sequence number","Sequence ID","V-DOMAIN Functionality", "V-GENE and allele","V-REGION identity %","J-GENE and allele","J-REGION identity %","D-GENE and allele","JUNCTION frame","JUNCTION (with frameshift)","CDR3-IMGT (with frameshift)")]
    
    df4 <- df2[names(df2) %in% c("Sequence number","Sequence ID","JUNCTION","JUNCTION (AA)","JUNCTION (with frameshift)","JUNCTION (AA) (with frameshift)","CDR3-IMGT","CDR3-IMGT (AA)")]
    
    df_chain1 <- merge(df3,df4,by=c("Sequence number","Sequence ID"))
    df_chain1 <- as.data.frame(df_chain1)
    df_chain1$`J-GENE and allele` <- gsub('Homsap ','',df_chain1$`J-GENE and allele`)
    df_chain1$`V-GENE and allele` <- gsub('Homsap ','',df_chain1$`V-GENE and allele`)
    df_chain1$`D-GENE and allele` <- gsub('Homsap ','',df_chain1$`D-GENE and allele`)
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
    df_chain1$`V-GENE and allele` <- gsub('TR','',df_chain1$`V-GENE and allele`)
    df_chain1$`D-GENE and allele` <- gsub('TR','',df_chain1$`D-GENE and allele`)
    
    df_chain1$JUNCTION <- toupper(df_chain1$JUNCTION) 
    df_chain1$V.sequence.quality.check <- ifelse(df_chain1$`V-REGION identity %`>=90,"quality sequence alignment","check chromatogram")
    df_chain1$clone_quality <- NA 
    df_chain1$comments <- NA
    df_chain1
    
  })
  
  output$IMGT2_out <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    IMGT2()
  })
  
  output$downloadTABLE_IMGTonly <- downloadHandler(
    filename = function(){
      paste("IMGT_only.QC",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- IMGT2()
      write.csv(df,file, row.names = FALSE)
    } )
  
  # QC after merged MiXCR.IMGT ---- 
  chain_merge <- reactive({
    df1 <- input.data.IMGT.MiXCR();
    df1 <- as.data.frame(df1)
    df <- subset(df1,df1$clone_quality=="pass")
    df <- as.data.frame(df)
    
    df2 <- df[c(1:7,11,22,13)]
    y = dim(df2)[2]-1
    y
    df2$V.Gene.IMGT <- gsub("TR","",df2$V.Gene.IMGT)
    df2$cloneCount <- 1
    
    if (input$chain=="ab") {

      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      head(df_name2)
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "-")))
      df_name4 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), "-")))
      df_name3$V1 <- gsub("A","",df_name3$V1)
      df_name3$V1 <- gsub("B","",df_name3$V1)
      
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      head(df_name5)
      df2$ID <- df_name2$V1
      
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$clone <- df_name4$V2
      names(df2)
      
      dim(df)
      
      chain1 <- df2[grep("A",df2$V.Gene.IMGT),]
      chain2 <- df2[grep("B",df2$V.Gene.IMGT),]
      # paste chain into name
      chain1$ID <- gsub("A-","-",chain1$ID)
      names(chain1)[1:y] <- paste(names(chain1)[1:y],"A",sep="_")
      head(chain1)
      chain2$ID <- gsub("B-","-",chain2$ID)
      names(chain2)[1:y] <- paste(names(chain2)[1:y],"B",sep="_")
      
      x <- names(chain2)[10:dim(chain2)[2]]
      
      merged_chain <- merge(chain1,chain2,by =x)
      head(merged_chain)
      merged_chain2 <- merged_chain[ , -which(names(merged_chain) %in% c("ID","Sequence.ID_A","Sequence.ID_B","V.DOMAIN.Functionality_A","V.DOMAIN.Functionality_B","D.GENE.and.allele_A","JUNCTION.frame_A","JUNCTION.frame_B"))]
      names(merged_chain2)
      dat <- merged_chain2
      dat$AJ <- paste(dat$V.Gene.IMGT_A,".",dat$J.GENE.and.allele_A,sep="")
      dat$BJ <- paste(dat$V.Gene.IMGT_B,".",dat$J.GENE.and.allele_B,sep="")
      dat$AJ <- gsub("[*]0.","",dat$AJ)
      dat$BJ <- gsub("[*]0.","",dat$BJ)
      dat$AJ <- gsub("TR","",dat$AJ)
      dat$AJBJ <- paste(dat$AJ,"_",dat$BJ,sep="")
      dat$AJ_aCDR3 <- paste(dat$AJ,dat$JUNCTION..AA._A,sep="_")
      dat$BJ_bCDR3 <- paste(dat$BJ,dat$JUNCTION..AA._B,sep="_")
      dat$AJ_aCDR3_BJ_bCDR3 <- paste(dat$AJ_aCDR3,dat$BJ_bCDR3,sep=" & ")
      head(dat)
      
      dat
      
    }
    else {
      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "-")))
      df_name4 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), "-")))
      df_name3$V1 <- gsub("G","",df_name4$V1)
      df_name3$V1 <- gsub("D","",df_name4$V1)
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      df2$ID <- df_name2$V1
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$clone <- df_name4$V2
      chain1 <- df2[grep("G",df2$V.Gene.IMGT),]
      chain2 <- df2[grep("D",df2$V.Gene.IMGT),]
      # paste chain into name
      chain1$ID <- gsub("G-",".",chain1$ID)
      names(chain1)[1:y] <- paste(names(chain1)[1:y],"G",sep="_")
      chain2$ID <- gsub("D-",".",chain2$ID)
      names(chain2)[1:y] <- paste(names(chain2)[1:y],"D",sep="_")
      x <- names(chain2)[10:dim(chain2)[2]]
      merged_chain <- merge(chain1,chain2,by=x)
      merged_chain2 <- merged_chain[ , -which(names(merged_chain) %in% c("ID","Sequence.ID_G","Sequence.ID_D","V.DOMAIN.Functionality_G","V.DOMAIN.Functionality_D","D.GENE.and.allele_G","JUNCTION.frame_G","JUNCTION.frame_D"))]
      dat <- merged_chain2
      head(dat)
      dat$GJ <- paste(dat$V.Gene.IMGT_G,"-",dat$J.GENE.and.allele_G,sep="")
      
      dat$DJ <- paste(dat$V.Gene.IMGT_D,"-",dat$J.GENE.and.allele_D,sep="")
      dat$GJ <- gsub("[*]0.","",dat$GJ)
      dat$DJ <- gsub("[*]0.","",dat$DJ)
      dat$GJDJ <- paste(dat$GJ,"_",dat$DJ,sep="")
      dat$GJ_gCDR3 <- paste(dat$GJ,dat$JUNCTION..AA._G,sep="_")
      dat$DJ_dCDR3 <- paste(dat$DJ,dat$JUNCTION..AA._D,sep="_")
      dat$GJ_gCDR3_DJ_dCDR3 <- paste(dat$AJ_aCDR3,dat$BJ_bCDR3,sep=" & ")
      
      
      dat
      
    }
  })
  output$chain_table <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    df1 <- input.data.IMGT.MiXCR();
    df1 <- as.data.frame(df1)
    
    a <- subset(df1 ,is.na(df1$clone_quality)==TRUE)
    if (dim(a)[1]>0) {
      df <- as.data.frame("please complete QC analysis")
      names(df) <- " "
      df
    }
    else {
      df <- chain_merge()
      df <- as.data.frame(df)
      df
    }
  })
  
  output$downloadTABLE3 <- downloadHandler(
    filename = function(){
      paste("Merged_A-B.or.G-D_",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- chain_merge()
      df <- as.data.frame(df)
      write.csv(df,file, row.names = FALSE)
    } )
  
  output$names.in.file2 <- renderPrint( {
    df <- chain_merge()
    df <- as.data.frame(df)
    names(df)
    
  })
  output$chain_table2 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    df1 <- input.data.IMGT.MiXCR();
    df1 <- as.data.frame(df1)
    
    a <- subset(df1 ,is.na(df1$clone_quality)==TRUE)
    if (dim(a)[1]>0) {
      df <- as.data.frame("please complete QC analysis")
      names(df) <- " "
      df
      
    }
    else {
      df <- chain_merge()
      df <- as.data.frame(df)
      
      your_list <- c("cloneCount",input$string.data2)
      your_list_df <- as.data.frame((unlist(strsplit(your_list, ', '))))
      names(your_list_df) <- "V1"
      your_list_df
      df.your_list <- df[names(df) %in% your_list_df$V1]
      names.df <- names(df.your_list[ , -which(names(df.your_list) %in% "cloneCount")])
      
      df2 <- as.data.frame(ddply(df.your_list,names.df,numcolwise(sum)))
      df2
      
    }
  })
  output$downloadTABLE4 <- downloadHandler(
    filename = function(){
      paste("summarised_CDR3_chain",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- chain_merge()
      df <- as.data.frame(df)
      
      your_list <- c("cloneCount",input$string.data2)
      your_list_df <- as.data.frame((unlist(strsplit(your_list, ', '))))
      names(your_list_df) <- "V1"
      your_list_df
      df.your_list <- df[names(df) %in% your_list_df$V1]
      names.df <- names(df.your_list[ , -which(names(df.your_list) %in% "cloneCount")])
      
      df2 <- as.data.frame(ddply(df.your_list,names.df,numcolwise(sum)))
      df2
      write.csv(df2,file, row.names = FALSE)
    } )
  
  
  # summary table of IMGT only -----
  input.data.IMGT_afterQC <- reactive({switch(input$dataset_IMGT_afterQC,"ab-test-data1" = test.data.ab.csv3(), "gd-test-data1" = test.data.gd.csv3(),"own1" = own.data.csv3())})
  test.data.ab.csv3 <- reactive({
    dataframe = read.csv("test-data/QC/IMGT_only.QC2021.08.29.csv",header=T) 
  })
  test.data.gd.csv3 <- reactive({
    dataframe = read.csv("test-data/QC/IMGT_only.QC2021.08.29.csv",header=T) 
  })
  own.data.csv3 <- reactive({
    inFile12 <- input$file_IMGT_afterQC
    if (is.null(inFile12)) return(NULL)
    
    else {
      dataframe <- read.csv(
        inFile12$datapath)}
    
  })
  
  
  chain_merge_IMGTonly <- reactive({
    df1 <- input.data.IMGT_afterQC();
    df1 <- as.data.frame(df1)
    df <- subset(df1,df1$clone_quality=="pass")
    df <- as.data.frame(df)
    df2 <- df[!names(df) %in% c("V.sequence.quality.check","clone_quality","comments","JUNCTION..with.frameshift.","CDR3.IMGT..with.frameshift.","JUNCTION..AA...with.frameshift.")]
    
    df.Vgene <- as.data.frame(do.call(rbind, strsplit(as.character(df2$V.GENE.and.allele), ",")))
    df2$V.GENE <- df.Vgene$V1
    y = dim(df2)[2]
    y
    df2$V.GENE <- gsub(" ","",df2$V.GENE)
    df2$cloneCount <- 1
    if (input$IMGT_chain2 =="ab") {
      
      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      head(df_name2)
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "-")))
      df_name4 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), "-")))
      df_name3$V1 <- gsub("A","",df_name3$V1)
      df_name3$V1 <- gsub("B","",df_name3$V1)
      
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      head(df_name5)
      df2$ID <- df_name2$V1
      
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$clone <- df_name4$V2
      names(df2)
      
      dim(df)
      
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
      merged_chain2 <- merged_chain[ , -which(names(merged_chain) %in% c("ID","Sequence.ID_A","Sequence.ID_B","V.DOMAIN.Functionality_A","V.DOMAIN.Functionality_B","D.GENE.and.allele_A","JUNCTION.frame_A","JUNCTION.frame_B"))]
      names(merged_chain2)
      dat <- merged_chain2
      dat$AJ <- paste(dat$V.GENE_A,".",dat$J.GENE.and.allele_A,sep="")
      dat$BJ <- paste(dat$V.GENE_B,".",dat$J.GENE.and.allele_B,sep="")
      dat$AJ <- gsub("[*]0.","",dat$AJ)
      dat$BJ <- gsub("[*]0.","",dat$BJ)
      dat$AJ <- gsub("TR","",dat$AJ)
      dat$AJBJ <- paste(dat$AJ,"_",dat$BJ,sep="")
      dat$AJ_aCDR3 <- paste(dat$AJ,dat$JUNCTION..AA._A,sep="_")
      dat$BJ_bCDR3 <- paste(dat$BJ,dat$JUNCTION..AA._B,sep="_")
      dat$AJ_aCDR3_BJ_bCDR3 <- paste(dat$AJ_aCDR3,dat$BJ_bCDR3,sep=" & ")
      head(dat)
      
      dat
      
      }
    else {
      df_name2 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "_")))
      head(df_name2)
      df_name3 <- as.data.frame(do.call(rbind, strsplit(as.character(df2$Sequence.ID), "-")))
      df_name4 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name2$V1), "-")))
      df_name3$V1 <- gsub("G","",df_name3$V1)
      df_name3$V1 <- gsub("D","",df_name3$V1)
      
      df_name5 <- as.data.frame(do.call(rbind, strsplit(as.character(df_name3$V1), "[.]")))
      head(df_name5)
      df2$ID <- df_name2$V1
      
      df2$Indiv.group <- df_name3$V1
      df2$Indiv <-df_name5$V1
      df2$group <- df_name5$V2
      df2$clone <- df_name4$V2

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
      dat$GJ <- paste(dat$V.GENE_G,".",dat$J.GENE.and.allele_G,sep="")
      dat$DJ <- paste(dat$V.GENE_D,".",dat$J.GENE.and.allele_D,sep="")
      dat$GJ <- gsub("[*]0.","",dat$GJ)
      dat$DJ <- gsub("[*]0.","",dat$DJ)
      dat$GJ <- gsub(", or GJ","/",dat$GJ)
      dat$GJDJ <- paste(dat$GJ,".",dat$DJ,sep="")
      dat$GJ_gCDR3 <- paste(dat$GJ,dat$JUNCTION..AA._G,sep="_")
      dat$DJ_dCDR3 <- paste(dat$DJ,dat$JUNCTION..AA._D,sep="_")
      dat$GJ_gCDR3_DJ_dCDR3 <- paste(dat$DJ_dCDR3,dat$GJ_gCDR3,sep=" & ")
      dat
    }
      
  })
  
  output$chain_table_IMGT.QC1 <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
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
     # df2 <- as.data.frame(ddply(df,names(df)[2:14],numcolwise(sum)))
     # df2
      
    }
  })
  output$downloadTABLE.QC1 <- downloadHandler(
    filename = function(){
      paste("paired_unsummarised",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- chain_merge_IMGTonly()
      df <- as.data.frame(df)
      write.csv(df,file, row.names = FALSE)
    } )
  
  
  # summarised IMGT only -----
  output$names.in.file3 <- renderPrint( {
    df <- chain_merge_IMGTonly()
    df <- as.data.frame(df)
    names(df)
    
  })
  
  
  output$chain_table_IMGT.QC3 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
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
      
      your_list <- c("cloneCount",input$string.data3)
      your_list_df <- as.data.frame((unlist(strsplit(your_list, ', '))))
      names(your_list_df) <- "V1"
      your_list_df
      df.your_list <- df[names(df) %in% your_list_df$V1]
      names.df <- names(df.your_list[ , -which(names(df.your_list) %in% "cloneCount")])
      
      df2 <- as.data.frame(ddply(df.your_list,names.df,numcolwise(sum)))
      df2
      
    }
  })

  
  output$downloadTABLE.QC3 <- downloadHandler(
    filename = function(){
      paste("paired_chain_CDR3",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- chain_merge_IMGTonly()
      df <- as.data.frame(df)
      
      your_list <- c("cloneCount",input$string.data3)
      your_list_df <- as.data.frame((unlist(strsplit(your_list, ', '))))
      names(your_list_df) <- "V1"
      your_list_df
      df.your_list <- df[names(df) %in% your_list_df$V1]
      names.df <- names(df.your_list[ , -which(names(df.your_list) %in% "cloneCount")])
      
      df2 <- as.data.frame(ddply(df.your_list,names.df,numcolwise(sum)))
      df2
        write.csv(df2,file, row.names = FALSE)
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
  
  observe({
    updateSelectInput(
      session,
      "count2",
      choices=names(input.data2()),
      selected = "cloneCount")
    
  })
  observe({
    updateSelectInput(
      session,
      "fill2",
      choices=names(input.data2()),
      selected = "AJ" )
    
  })
  observe({
    updateSelectInput(
      session,
      "sub_group2",
      choices=names(input.data2()),
      selected = "AJBJ")
    
  })
  observe({
    updateSelectInput(
      session,
      "wrap",
      choices=names(input.data2()),
      selected = "group")
    
  })
  
  
  cols <- reactive({
    dat <- input.data2();
    dat <- as.data.frame(dat)
    
    num <- unique(dat[names(dat) %in% input$fill2])
    col.gg <- gg_fill_hue(dim(num)[1])
    unique.col <- as.data.frame(unique(dat[grep(input$fill2,names(dat))]))
    unique.col <- unique.col; 
    palette1 <- distinctColorPalette(dim(unique.col)[1])
    
    if (input$tree_colour.choise == "default") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col", i, sep="_"), paste(num[i,]), col.gg[i])        
      })
    }
    else if (input$tree_colour.choise == "random") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col", i, sep="_"), paste(num[i,]), palette1[i])        
      })
      
    }
    
    else {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col", i, sep="_"), paste(num[i,]), "grey")        
      })
      
      
    }
    
  })
  
  
  output$myPanel <- renderUI({cols()})
  
  colors <- reactive({
    dat <- input.data2();
    dat <- as.data.frame(dat)
    num <- unique(dat[names(dat) %in% input$fill2])
    lapply(1:dim(num)[1], function(i) {
      input[[paste("col", i, sep="_")]]
    })
  })
  
  
  
  tree_plot_dynamic <- reactive({
    dat <- input.data2();
    dat <- as.data.frame(dat)
    
    if (is.null(input$col_1)) {
      cols <- rep("#000000", ncol(dat))
    } else {
      cols <- unlist(colors())
    }
    
    df1 <- dat[names(dat) %in% c(input$count2,input$fill2,input$sub_group2,input$wrap)]
    df2 <- as.data.frame(ddply(dat,names(df1)[-c(1)],numcolwise(sum)))
    unique.col <- as.data.frame(unique(dat[names(dat) %in% input$fill2]))
    names(unique.col) <- "V1"
    unique.col$tree_palette <- cols
    df3 <- as.data.frame(merge(df2,unique.col,by.x=input$fill2,by.y = "V1"))
    
    if (input$tree.lab == "yes") {
      vals22$Treemap22 <- ggplot(df3, aes(area = get(input$count2),
                                          fill = get(input$fill2),
                                          subgroup = get(input$sub_group2))) +
        geom_treemap(aes(alpha = 1),colour="white",show.legend = F, fill = df3$tree_palette) +
        geom_treemap_subgroup_border(colour = "white", show.legend = F,size=12) +
        geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 1, family = "serif",
                                   colour = "black", fontface = "italic", min.size = 0,show.legend = F) +
        facet_wrap(~get(input$wrap)) +
        theme(strip.text = element_text(size = 20, family = "serif"))
      vals22$Treemap22

    }
    
    
    else {
      vals22$Treemap22 <- ggplot(df3, aes(area = get(input$count2),
                                          fill = get(input$fill2),
                                          subgroup = get(input$sub_group2))) +
        geom_treemap(aes(alpha = 1),colour="white",show.legend = F, fill = df3$tree_palette) +
        geom_treemap_subgroup_border(colour = "white", show.legend = F,size=12) +
        facet_wrap(~get(input$wrap))+
        theme(strip.text = element_text(size = 20, family = "serif"))
      vals22$Treemap22

    }

    
  })
  ?facet_wrap
  
  output$Treemap2 <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    print(tree_plot_dynamic())
  })
  
  output$downloadPlot_scTREE <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_explore_treemap_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_tree,height=input$height_tree, onefile = FALSE) # open the pdf device
      print(tree_plot_dynamic())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_scTREE <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_explore_treemap_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_tree, height = input$height_png_tree, res = input$resolution_PNG_tree)
      print(tree_plot_dynamic())
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  # circular plot =====
  observe({
    updateSelectInput(
      session,
      "group_column",
      choices=names(input.data2()),
      selected = "group")
    
  }) # group 
  select_group <- function () {
    df <- input.data2();
    df2 <- as.data.frame(unique(df[names(df) %in% input$group_column]))
    df2 <- as.data.frame(df2)
    #names(df2) <- "V1"
    df2
  }
  observe({
    updateSelectInput(
      session,
      "group_selected2",
      choices=select_group())
    
  }) # group 
  observe({
    updateSelectInput(
      session,
      "chain1",
      choices=names(input.data2()),
      selected = "AJ")
    
  }) # chain 1
  observe({
    updateSelectInput(
      session,
      "chain2",
      choices=names(input.data2()),
      selected = "BJ")
    
  }) # chain 2
  
  output$table_display <- renderTable({
    dat <- input.data2();
    dat <- as.data.frame(dat)
    
    
    
    dat <- subset(dat, get(input$group_column)==input$group_selected2)
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    head(hierarchy, n=2)
  })
  cols_circ <- reactive({
    dat <- input.data2();
    dat <- as.data.frame(dat)
    dat <- subset(dat, get(input$group_column)==input$group_selected2)
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    df.col1 <- as.data.frame(unique(hierarchy[,1]))
    names(df.col1) <- "V1"
    
    df.col.j <- as.data.frame(unique(hierarchy[,2]))
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    df.col.2
    palette2 <- distinctColorPalette(dim(df.col.2)[1])
    palette_rainbow <- rev(rainbow(length(t(df.col.2))))
    
    if (input$colour_cir == "rainbow") {
    lapply(1:dim(df.col.2)[1], function(i) {
      colourInput(paste("col.cir", i, sep="_"), paste(df.col.2[i,]), palette_rainbow[i])        
    }) }
    
    else if (input$colour_cir == "random") {
      lapply(1:dim(df.col.2)[1], function(i) {
        colourInput(paste("col.cir", i, sep="_"), paste(df.col.2[i,]), palette2[i])        
      }) }
    else  {
      lapply(1:dim(df.col.2)[1], function(i) {
        colourInput(paste("col.cir", i, sep="_"), paste(df.col.2[i,]), "grey")        
      }) }
    
  })
  output$myPanel_circ <- renderUI({cols_circ()})
  colors_cir <- reactive({
    dat <- input.data2();
    dat <- as.data.frame(dat)
    dat <- subset(dat, get(input$group_column)==input$group_selected2)
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
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
  Circular_plot2 <- function () {
    dat <- input.data2();
    dat <- as.data.frame(dat)
    dat <- subset(dat, get(input$group_column)==input$group_selected2)
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    df.col1 <- as.data.frame(unique(hierarchy[,1]))
    names(df.col1) <- "V1"
    
    df.col.j <- as.data.frame(unique(hierarchy[,2]))
    names(df.col.j) <- "V1"
      df.col.2 <- rbind(df.col1,df.col.j)
      col2 = unlist(colors_cir())
      col2
      df.col.2$colour <- col2
      df.col.2
      

      
      par(mar = rep(0, 4), cex=0.8)
      
      grid.col <- as.data.frame(as.matrix(t(as.data.frame(df.col.2$colour))))
      names(grid.col) <- df.col.2$V1
      grid.col <- as.matrix(grid.col)
      
      if (input$circ.lab=="yes") {
        
        
        circos.clear()
        #par(new = TRUE) # <- magic
        circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-1, 1))
        chordDiagram(hierarchy, annotationTrack = "grid", grid.col = grid.col,
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
          
          
          #circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
          #            facing = "clockwise", niceFacing = TRUE, adj = c(0, 0))
          
          
        }, bg.border = NA)
        
      }

      else {
        
        circos.clear()
        #par(new = TRUE) # <- magic
        circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-1, 1))
        chordDiagram(hierarchy, annotationTrack = "grid", grid.col = grid.col,
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
      paste("TCR_explore_circular_plot_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_circ,height=input$height_circ, onefile = FALSE) # open the pdf device
      print(Circular_plot2())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_circ <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_explore_circular_plot_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_circ, height = input$height_png_circ, res = input$resolution_PNG_circ)
      print(Circular_plot2())
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  
  # bar graph =====
  observe({
    updateSelectInput(
      session,
      "aa.or.nt",
      choices=names(input.data2()),
      selected = "JUNCTION_B") 
  }) # amino acid or nucleotides column
  observe({
    updateSelectInput( 
      session,
      "selected_group_len",
      choices=select_group()) }) # group
  
  Chain1_length <- function () {
    df <- input.data2(); 
    # not summarised dataset
    df <- as.data.frame(df)
    df$len1 <- nchar(df[,grep(input$aa.or.nt,names(df))])
    df <- subset(df,get(input$group_column)==input$selected_group_len)
    
    if (input$graph_type == "histogram") {
      vals4$bar.len <- ggplot(df,aes(x=len1)) +
        geom_histogram(fill = input$hist_col, bins = input$bin) + 
        theme_bw()  +
        theme(legend.title = element_blank(),
              legend.position = "none") +
        
        labs(y="Count",
             x="",
             title="") +
        theme(
          axis.title.y = element_text(colour="black",family="serif"),
          axis.text.y = element_text(colour="black",family="serif"),
          axis.text.x = element_text(colour="black",family="serif",angle=90),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
          legend.text = element_blank()
        ) 
      vals4$bar.len
    }
    else {
      vals4$bar.len <- ggplot(df,aes(x=len1)) +
        #geom_histogram(fill = input$hist_col, bins = input$bin) + 
        geom_density(color = input$hist_col) +
        #facet_wrap(~sub) +
        theme_bw()  +
        theme(legend.title = element_blank(),
              legend.position = "none") +
        #scale_fill_manual(values=c() +
        
        labs(y="CDF",
             x="",
             title="") +
        # Add a line showing the alpha = 0.01 level
        theme(
          axis.title.y = element_text(colour="black",family="serif"),
          axis.text.y = element_text(colour="black",family="serif"),
          axis.text.x = element_text(colour="black",family="serif",angle=90),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
          legend.text = element_blank()
        )
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
      paste("TCR_explore_length_plot_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_length,height=input$height_length, onefile = FALSE) # open the pdf device
      print(Chain1_length())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_length <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_explore_length_plot_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_length, height = input$height_png_length, res = input$resolution_PNG_length)
      print(Chain1_length())
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  # bar plot -----
  observe({
    updateSelectInput(
      session,
      "variable_chain",
      choices=names(input.data2()),
      selected = "AJ")
    
  })
  
  observe({
    updateSelectInput( 
      session,
      "selected_group_chain",
      choices=select_group()) }) # group
  Chain1_usage <- function () {
    df <- input.data2(); 
    df <- as.data.frame(df)
    df <- subset(df, get(input$group_column)==input$selected_group_chain)
    df2 <- as.data.frame(ddply(df,c(input$variable_chain),numcolwise(sum)))[1:2]
    names(df2) <- c("chain","cloneCount")
    
    df2 <- df2[order(df2$cloneCount),]
    df2$chain <- factor(df2$chain, levels = unique(df2$chain),labels = df2$chain)
    
    if (input$graph_bar_type == "count") {
      
      vals5$bar.usage <- ggplot(df2,aes(x=chain,y=cloneCount)) +
        geom_bar(stat="identity", position = "dodge",fill=input$colour_bar.usage) +
        theme_bw()  +
        theme(legend.title = element_blank(),
              legend.position = "bottom") +
        #scale_fill_manual(values=c("red","blue","darkgreen")) +
        labs(y="count",
             x="",
             title="") +
        theme(
          axis.title.y = element_text(colour="black",family="serif",size = 12),
          axis.text.y = element_text(colour="black",family="serif",size = 12),
          axis.text.x = element_text(colour="black",family="serif",angle=0,size = 12),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family="serif",size = 12),
          legend.text = element_text(colour="black",family="serif")
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
          axis.title.y = element_text(colour="black",family="serif",size = 12),
          axis.text.y = element_text(colour="black",family="serif",size = 12),
          axis.text.x = element_text(colour="black",family="serif",angle=0,size = 12),
          axis.title.x = element_text(colour="black",angle=0,hjust=.5,vjust=.5,face="plain",family="serif",size = 12),
          legend.text = element_text(colour="black",family="serif")
        ) +
        coord_flip()
      vals5$bar.usage
      
    }
    

    
  }
  output$Chain1_usage <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    Chain1_usage()
  }) 
  output$downloadPlot_chain.usage <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("bar.plot_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_chain.usage,height=input$height_chain.usage, onefile = FALSE) # open the pdf device
      print(Chain1_usage())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_chain.usage <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("bar.plot_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png_chain.usage, height = input$height_png_chain.usage, res = input$resolution_PNG_chain.usage)
      print(Chain1_usage())
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  
  # motif -----
  observe({
    updateSelectInput(
      session,
      "aa.or.nt2",
      choices=names(input.data2()),
      selected = "JUNCTION..AA._B")})
  observe({
    updateSelectInput(
      session,
      "group_selected_motif",
      choices=select_group()) })
  
  
  output$Motif <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE), {
    df <- input.data2();
    df <- as.data.frame(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt2)),numcolwise(sum)))
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt2])
    df_subset <- subset(df_unique,df_unique$len1==input$len)
    df_subset <- subset(df_subset,get(input$group_column)==input$group_selected_motif)
    
    motif <- as.data.frame(t(as.data.frame(strsplit(df_subset[,grep(input$aa.or.nt2,names(df_subset))], ""))))
    rownames(motif) <- 1:dim(motif)[1]
    motif
    
  })
  
  output$length <- renderPrint( {
    df <- input.data2();
    df$len1 <- nchar(df[,grep(input$aa.or.nt2,names(df))])
    df <- df[order(df$len1),]
    df <- subset(df,get(input$group_column)==input$group_selected_motif)
    df_len <- unique(df$len1)
    cat("The dataset contains CDR3 lengths of:",  df_len)
  })
  
  output$length.table <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE), {
    df <- input.data2();
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt2)),numcolwise(sum)))
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt2])
    df_unique
    

  })
  
  Motif_plot2 <- reactive({
    df <- input.data2();
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
    motif<-new("pfm", mat=motif_count, name="",
               color=colorset(alphabet="AA",
                              colorScheme="chemistry"))
    motif
    })
  
  output$Motif_plot <- renderPlot( {
    motif <- Motif_plot2()
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    plot(motif)
  })
  
  output$downloadPlot_motif <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_explore_motif_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_motif,height=input$height_motif, onefile = FALSE) # open the pdf device
      plot(Motif_plot2())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_motif <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_explore_motif_", gsub("/", "-", x), ".png", sep = "")
    },    
    content = function(file) {
      png(file, width = input$width_png_motif, 
          height = input$height_png_motif, 
          res = input$resolution_PNG_motif)
      plot(Motif_plot2())
      dev.off()},   contentType = "application/png" # MIME type of the image
  )

  # motif NT ------
  observe({
    updateSelectInput(
      session,
      "aa.or.nt3",
      choices=names(input.data2()),
      selected = "JUNCTION_A")
    
  })
  
  observe({
    updateSelectInput(
      session,
      "group_selected",
      choices=select_group())
    
  })
  
  output$Motif_nt <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5), {
    df <- input.data2();
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
    df <- input.data2();
    df <- as.data.frame(df)
    df$len1 <- nchar(df[,grep(input$aa.or.nt3,names(df))])
    df <- subset(df,get(input$group_column)==input$group_selected)
    df <- df[order(df$len1),]
    df_len <- unique(df$len1)
    cat("The dataset contains nucleotide CDR3 lengths of:",  df_len)
    })
  output$length.table_nt <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE), {
    df <- input.data2();
    df <- as.data.frame(df)
    df_unique <- as.data.frame(ddply(df,(c(input$group_column,input$aa.or.nt3)),numcolwise(sum)))
    df_unique$len1 <- nchar(df_unique[,names(df_unique) %in% input$aa.or.nt3])
    df_unique
    
    
  })
  Motif_plot2_nt <- reactive( {
    df <- input.data2();
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
      paste("TCR_explore_motif_",gsub("/", "-", x), ".pdf", sep = "")
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
      paste("TCR_explore_motif_", gsub("/", "-", x), ".png", sep = "")
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
  # pie graph -----
  observe({
    updateSelectInput(
      session,
      "pie_chain",
      choices=names(input.data2()),
      selected = "AJBJ")
    
  })
  
  
  cols_pie <- reactive({
    dat <- input.data2();
    dat <- as.data.frame(dat)
    df <- as.data.frame(ddply(dat,(c(input$group_column,input$pie_chain)),numcolwise(sum)))
    names(df) <- c("group","chain","cloneCount")
    num <- unique(df$chain)
    col.gg <- gg_fill_hue(length(num))
    length(num)
    palette1 <- distinctColorPalette(length(num))
    
    if (input$pie_colour.choise == "default") {
      lapply(1:length(num), function(i) {
        colourInput(paste("col.pie", i, sep="_"), paste(num[i]), col.gg[i])        
      })
    }
    else if (input$pie_colour.choise == "random") {
      lapply(1:length(num), function(i) {
        colourInput(paste("col.pie", i, sep="_"), paste(num[i]), palette1[i])        
      })
      
    }
    
    else {
      lapply(1:length(num), function(i) {
        colourInput(paste("col.pie", i, sep="_"), paste(num[i]), "grey")        
      })
      
      
    }
    
  })
  output$myPanel_pie <- renderUI({cols_pie()})
  
  colors_pie <- reactive({
    dat <- input.data2();
    dat <- as.data.frame(dat)
    df <- as.data.frame(ddply(dat,(c(input$group_column,input$pie_chain)),numcolwise(sum)))
    names(df) <- c("group","chain","cloneCount")
    
    
    num <- unique(df$chain)
    lapply(1:length(num), function(i) {
      input[[paste("col.pie", i, sep="_")]]
    })
  })
  
  pie_chart <- function() {
    set.seed(123)
    dat <- input.data2();
    dat <- as.data.frame(dat)
    cols <- unlist(colors_pie())
    
    df <- as.data.frame(ddply(dat,(c(input$group_column,input$pie_chain)),numcolwise(sum)))
    names(df) <- c("group","chain","cloneCount")
    palette <- cols
    a <- unique(df$chain)
    df$chain <- factor(df$chain,levels = a,labels=a)
    
   vals9$pie <- ggplot(df, aes(x="", y=cloneCount, fill=chain)) +
      geom_bar(width = 1, stat = "identity",aes(colour = "black")) +
      scale_fill_manual(values=palette) +
      scale_color_manual(values = "black") +
      coord_polar("y", start=1) + facet_wrap(~group) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid  = element_blank(),
            axis.title.y= element_blank(),
            legend.position = "bottom") +
      guides(color = "none", size = "none")
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
      paste("TCR_explore_pie_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_pie,height=input$height_pie, onefile = FALSE) # open the pdf device
      plot(pie_chart())
      dev.off()}, contentType = "application/pdf" )
  
  output$downloadPlotPNG_pie <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_explore_pie_", gsub("/", "-", x), ".png", sep = "")
    },    
    content = function(file) {
      png(file, width = input$width_png_pie, 
          height = input$height_png_pie, 
          res = input$resolution_PNG_pie)
      plot(pie_chart())
      dev.off()},   contentType = "application/png" # MIME type of the image
  )
  
  
  
  # heatmap ----
  observe({
    updateSelectInput(
      session,
      "heatmap_1",
      choices=names(input.data2()),
      selected = "AJ")
    
  })
  

  
  
  heatmap_matrix <- function() {
    # group by chain dat <- read.csv("test-data/Group/Merged_A-B_2021.08.02.csv")
    dat <- input.data2();
    dat <- as.data.frame(dat)
    
    df <- as.data.frame(ddply(dat,(c("group",input$heatmap_2)),numcolwise(sum)))
    names(df) <- c("group","heatmap_2","cloneCount")
    upper <- max(df$cloneCount)
    
   # df <- scale(t(df))
    ord <- hclust( dist(df, method = "euclidean"), method = "ward.D" )$order
    df2 <- df[order(ord),]
    a <- unique(df2$heatmap_2)
    df$heatmap_2 <- factor(df$heatmap_2,levels = a,labels=a)
    vals10$heatmap_clonal <-  ggplot(df, aes(x=factor(group), y=factor(heatmap_2)))+ 
      geom_tile(aes(fill=cloneCount), color=NA)+
      scale_fill_gradient(low="#FFFFFF", high="#FF0000", limits=c(0,upper))+
      coord_fixed() +
      theme_bw(base_size = 16)+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),axis.text.x=element_text(angle=90,hjust=1),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_blank(), axis.title.y=element_blank()) 
      
      
      vals10$heatmap_clonal
    
  }
  # y.axis heatmap
  observe({
    updateSelectInput(
      session,
      "group.heatmap",
      choices=names(input.data2()),
      selected = "AJ")
    
  })
  # x-axis heatmap
  observe({
    updateSelectInput(
      session,
      "heatmap_2",
      choices=names(input.data2()),
      selected = "BJ")
    
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
    # group by chain dat <- read.csv("test-data/Group/Merged_A-B_2021.08.02.csv")
    dat <- input.data2();
    dat <- as.data.frame(dat)
    dat <- subset(dat, get(input$group_column)==input$group_selected3)
    df <- as.data.frame(ddply(dat,(c(input$heatmap_2,input$group.heatmap)),numcolwise(sum)))
    names(df) <- c("group","chain","cloneCount")
    head(df)
    
    mat <- acast(df, chain~group, value.var="cloneCount")
    mat[is.na(mat)] <- 0

    
    vals13$heatmap_clonal2 <- ggheatmap(mat,color=colorRampPalette(c( "white","#efc000","red"))(100), scale = "none",cluster_rows  = TRUE,cluster_cols = TRUE,
                                        text_show_rows = row.names(mat))%>%
      ggheatmap_theme(1,theme =list(
        theme(axis.text.x = element_text(angle = 90,face = "bold"),
              axis.text.y = element_text(colour = "black",face = "bold"))
      ))
    
    
    vals13$heatmap_clonal2
    
  }
  heatmap_matrix3 <- function() {
    # group by chain dat <- read.csv("test-data/Group/Merged_A-B_2021.08.02.csv")
    dat <- input.data2();
    dat <- as.data.frame(dat)
  
    df <- as.data.frame(ddply(dat,(c(input$heatmap_2,input$group.heatmap)),numcolwise(sum)))
    names(df) <- c("group","chain","cloneCount")
    head(df)
    
    mat <- acast(df, chain~group, value.var="cloneCount")
    mat[is.na(mat)] <- 0
    
    
    vals14$heatmap_clonal3 <- ggheatmap(mat,color=colorRampPalette(c( "white","#efc000","red"))(100),cluster_rows  = TRUE,cluster_cols = TRUE,
                                        text_show_rows = row.names(mat))%>%
      ggheatmap_theme(1,theme =list(
        theme(axis.text.x = element_text(angle = 90,face = "bold"),
              axis.text.y = element_text(colour = "black",face = "bold"))
      ))
    
    
    vals14$heatmap_clonal3
    
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
      paste("TCR_explore_heatmap_",gsub("/", "-", x), ".pdf", sep = "")
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
      paste("TCR_explore_heatmap_", gsub("/", "-", x), ".png", sep = "")
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
  
  input.data.diversity.index <- reactive({switch(input$shannon.index,"gd.test.index" = test.data.gd.index.csv(),"own.index" = own.data.index.csv())})
  
  test.data.gd.index.csv <- reactive({
  dataframe = read.csv("../Sanger/diversity.csv",header=T)   
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
  inv.simpson.index <- function() {
    plots <- input.data.diversity.index();
    
    validate(
      need(nrow(plots)>0,
           "Upload clone file")
    )
    
    a <- matrix(nrow=1,ncol=dim(plots)[2])
    b <- matrix(nrow=1,ncol=dim(plots)[2])
    d <- matrix(nrow=1,ncol=dim(plots)[2])
    
    for( i in 1:dim(plots)[2]) {
      
      samp <- plots[,i]
      samp <- na.omit(samp)
      a[,i] <- diversity(samp,"invsimpson")
      b[,i] <- sum(samp)
      d[,i] <- nrow(as.data.frame(samp))
    }
    
    a1 <- rbind(a,b,d)  
    a1 <- as.data.frame(a1)
    names(a1) <- names(plots)
    a1 <- rbind(a,b,d)  
    a1 <- as.data.frame(a1)
    names(a1) <- names(plots)
    
    df_name <- as.data.frame(do.call(rbind, strsplit(as.character(names(a1)), "_")))
    head(df_name) 
    
    a2 <- as.data.frame(t(a1))
    names(a2) <- c("inv.simpson.index","total # clones","unique # clones")
    both <- cbind(a2,df_name)
    
    both
    
    
  }
  
  
  observe({
    updateSelectInput(
      session,
      "group.index",
      choices=names(inv.simpson.index()),
      selected = "V2")
    
  })
  observe({
    updateSelectInput(
      session,
      "group2.index",
      choices=names(inv.simpson.index()),
      selected = "V1")
    
  })
  observe({
    updateSelectInput(
      session,
      "x.axis.index",
      choices=names(inv.simpson.index()),
      selected = "total # clones")
    
  })
  
  cols_simp.index <- reactive({
    dat <- inv.simpson.index();
    dat <- as.data.frame(dat)
    
    selected.col <- dat[names(dat) %in% input$group2.index]
    names(selected.col) <- "V1"
    dat[names(dat) %in% input$group2.index] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))
    
    
    num <- unique(dat[names(dat) %in% input$group2.index])
    col.gg <- gg_fill_hue(dim(num)[1])
    palette1 <- distinctColorPalette(dim(num)[1])
    
    if (input$inv.simp_colour.choise == "default") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.inv.simpson", i, sep="_"), paste(num[i,]), col.gg[i])        
      })
    }
    else if (input$inv.simp_colour.choise == "random") {
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
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    limits = c(1,10^6),
                    labels = trans_format("log10", math_format(10^.x))) +
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

  
  # 

  
  DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    
    clonal.file <-  input.data.clone.file();
    clonal.file <- as.data.frame(clonal.file)
    with.clone.data()
  })
  
  output$table_display.diversity <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      
    dat <- inv.simpson.index()
    dat <- as.data.frame(dat)
    dat
    })

  table.inv.simpson <- function () {
    dat <- inv.simpson.index()
    dat <- as.data.frame(dat)
    dat
    
  }
  
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
      dev.off()}, contentType = "application/png" # MIME type of the image
  )
  
  
  
  # FACS index data -----
  input.data_FACS <- reactive({switch(input$dataset3,"test-FACS" = test.data_FACS(), "own_FACS" = own.data_FACS())})
  test.data_FACS <- reactive({
    read.FCS("../../Results_2020/Flow Cytometry/FACS Aria/AFKM200923_PBMC/Specimen_001_INX_PBMCs 4905 CBZ sort IFN 1_014.fcs")
  })
  own.data_FACS <- reactive({
    input$file_FACS
    if (is.null(input$file_FACS)) return(NULL)
    
    else {
      read.FCS(input$file_FACS$datapath)}
    
  })
  output$FACS_to_CSV <- DT::renderDataTable(escape = FALSE,{
    
    samp <-  input.data_FACS()
    req(samp)
    samp_index <- getIndexSort(samp)
    samp_index
  })
  output$FACS.CSV <- DT::renderDataTable( {
    samp <-  input.data_FACS()
    req(samp)
    samp_index <- getIndexSort(samp)
    datatable(samp_index[1:dim(samp_index)[1],], extensions = "Buttons", options = list(searching = TRUE,
                                                                                        ordering = TRUE,
                                                                                        buttons = c('copy','csv', 'excel'),
                                                                                        dom = 'Bfrtip',
                                                                                        pageLength=5, 
                                                                                        lengthMenu=c(2,5,10,20,50,100), 
                                                                                        scrollX = TRUE
    ))
  }, server = FALSE)  
  merged.index <- function(){
    samp <-  input.data_FACS()
    req(samp)
    samp_index <- getIndexSort(samp)
    samp_index <- getIndexSort(samp)
    head(samp_index)
    replace_ID <- read.csv("test-data/Index/Loc_to_ID.csv")
    head(replace_ID)
    
    index_updated_ID <- merge(samp_index,replace_ID,by=c("XLoc","YLoc"))
    index_updated_ID
    
  }

  output$merged.clone <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    merged.index()
  })
  input.data.clone.file <- reactive({switch(input$data_clone.index, "gd.test.clone" = test.data.gd.index.csv2(),"own.clone.file" = own.data.clone.file.csv())})
  test.data.gd.index.csv2 <- reactive({
    dataframe = read.csv("test-data/Index/4905_VD1_paired_ex_vivo1.csv",header=T,row.names = 1) 
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
    clonal.file <-  input.data.clone.file();
    clonal.file <- as.data.frame(clonal.file)
    
    validate(
      need(nrow(clonal.file)>0,
           as.data.frame("Upload clone file"))
    )
    
    index.clonal.file <- merge(clonal.file,index_updated_ID,by="clone")
    index.clonal.file
    
    
    
    
    
    
  }
  
  output$merged.index.clone <- DT::renderDataTable(escape = FALSE, options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    
    clonal.file <-  input.data.clone.file();
    clonal.file <- as.data.frame(clonal.file)
    with.clone.data()
  })
  
  output$downloadTABLE_FACS <- downloadHandler(
    filename = function(){
      paste("TCR_explore_index.clonal.",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      write.csv(with.clone.data(),file, row.names = FALSE)
    }
  )
  
  
  
  input.data_CSV1 <-  reactive({switch(input$dataset7,"test-csv"=test.data_csv1(),"own_csv" = own.data_CSV1())})
  test.data_csv1 <- reactive({
    dataframe = read.csv("test-data/Index/4905.csv",header = T,row.names = 1)
  })
  own.data_CSV1 <- reactive({
    inFile_CSV1 <- input$file_FACS.csv1
    if (is.null(inFile_CSV1)) return(NULL)
    
    else {
      dataframe <- read.csv(inFile_CSV1$datapath, header=T)}
    
  })


  

  vals15 <- reactiveValues(complex_dot=NULL)
  
  output$names.in.file <- renderPrint( {
    df <- input.data_CSV1();
    df <- as.data.frame(df)
    names(df)
    
  })
  
  observe({
    updateSelectInput(
      session,
      "V.gene.1",
      choices=names(input.data_CSV1()),
      selected = "TRGV")
  })
  observe({
    updateSelectInput(
      session,
      "CDR3.1",
      choices=names(input.data_CSV1()),
      selected = "gCDR3")
  })
  observe({
    updateSelectInput(
      session,
      "V.gene.2",
      choices=names(input.data_CSV1()),
      selected = "TRDV")
  })
  observe({
    updateSelectInput(
      session,
      "CDR3.2",
      choices=names(input.data_CSV1()),
      selected = "dCDR3")
  })
  
  index.cleaning1 <- reactive({
    df <- input.data_CSV1();
    validate(
      need(nrow(df)>0,
           "Upload file")
    )
    df <- as.data.frame(df)
    head(df)
    df$cloneCount <- 1
    
    your_list <- c("cloneCount",input$string.data)
    your_list_df <- as.data.frame((unlist(strsplit(your_list, ', '))))
    names(your_list_df) <- "V1"
    your_list_df
    df.your_list <- df[names(df) %in% your_list_df$V1]
    head(df.your_list)
    
    names.df <- names(df.your_list[ , -which(names(df.your_list) %in% "cloneCount")])
    
    df2 <- as.data.frame(ddply(df.your_list,names.df,numcolwise(sum)))
    
    df1 <- subset(df2,df2$cloneCount>input$numeric.cloneCount)
    
    df1$clonal <- "yes"
    colnames(df1)[which(names(df1) == "cloneCount")] <- "# of clones"
    a <- df1[names(df1) %in% input$V.gene.1]
    names(a) <- "V1"
    b <- df1[names(df1) %in% input$CDR3.1]
    names(b) <- "V1"
    group.CDR <- df1[names(df1) %in% "group"]
    names(group.CDR) <- "V1"
    d <- df1[names(df1) %in% input$V.gene.2]
    names(d) <- "V1"
    e  <- df1[names(df1) %in% input$CDR3.2]
    names(e) <- "V1"
    df1$gene.CDR3.1 <- paste(a$V1,b$V1,group.CDR$V1,sep="_")
    df1$gene.CDR3.2 <- paste(d$V1,e$V1,group.CDR$V1,sep = "_")
    df1$gene.CDR3.both <- paste(a$V1,b$V1,d$V1,e$V1,group.CDR$V1,sep = "_")
    df1
    a2 <- merge(df,df1,by=names.df,all=T)
    #a[is.na(a)] <- 'unique'
    a2 <- as.data.frame(a2)
    a2
    
  })
  
  
  output$downloadTABLE_cleaning <- downloadHandler(
    filename = function(){
      paste("colouring column",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      write.csv(index.cleaning1(),file, row.names = FALSE)
    })
  output$table.index.1 <- DT::renderDataTable(escape = FALSE,options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    df <- index.cleaning1()
    df
    
  })
  
  
  # creating the dot plot ----
  
  input.data_CSV2 <-  reactive({switch(input$dataset_index.2,"test-csv"=test.data_csv2(),"own_csv" = own.data_CSV2())})
  test.data_csv2<- reactive({
    dataframe = read.csv("test-data/Index/colouring column2021.09.08 (1).csv",header = T)
  })
  own.data_CSV2 <- reactive({
    inFile_CSV2 <- input$file_FACS.csv2
    if (is.null(inFile_CSV2)) return(NULL)
    
    else {
      dataframe <- read.csv(inFile_CSV2$datapath, header=T)}
    
  })
  
  
  
  observe({
    updateSelectInput(
      session,
      "x.axis2",
      choices=names(input.data_CSV2()),
      selected = "CD69")
  })
  
  observe({
    updateSelectInput(
      session,
      "y.axis2",
      choices=names(input.data_CSV2()),
      selected = "CD38")
  })
  
  observe({
    updateSelectInput(
      session,
      "group_complex_dot",
      choices=names(input.data_CSV2()),
      selected = "group")
  })
  
  cols.FACS.index <- reactive({
    dat <- input.data_CSV2();
    
    validate(
      need(nrow(dat)>0,
           "Upload file")
    )
    dat <- as.data.frame(dat)
    dat[is.na(dat)] <- "not_clonal"
    num <- unique(dat[names(dat) %in% input$group_complex_dot])
    col.gg <- gg_fill_hue(dim(num)[1])
    unique.col <- as.data.frame(unique(dat[grep(input$group_complex_dot,names(dat))]))
    palette1 <- distinctColorPalette(dim(unique.col)[1])
    
    if (input$FACS.index_colour.choise == "default") {
      lapply(1:dim(num)[1], function(i) {
        colourInput(paste("col.FACS.index", i, sep="_"), paste(num[i,]), col.gg[i])        
      })
    }
    else if (input$FACS.index_colour.choise == "random") {
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
  
  
  output$myPanel.FACS.index <- renderUI({cols.FACS.index()})
  
  colors.FACS.index <- reactive({
    dat <- input.data_CSV2();
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
  
  
  
  
  
  
  dot_plot.complex <- reactive({
    index <- input.data_CSV2();
    validate(
      need(nrow(index)>0,
           "Upload file")
    )
    
    
    index <- as.data.frame(index)
    y_lable1 <- bquote(.(input$y.axis2))
    x_lable1 <-  bquote(.(input$x.axis2))
    
    index[is.na(index)] <- "not_clonal"
    selected.col <- index[names(index) %in% input$group_complex_dot]
    names(selected.col) <- "V1"
    index[names(index) %in% input$group_complex_dot] <- factor(selected.col$V1, levels = unique(selected.col$V1),labels = unique(selected.col$V1))
    palette.complex <- unlist(colors.FACS.index())
    
    vals15$complex_dot <- ggplot(index, aes(x=get(input$x.axis2), y=get(input$y.axis2),colour = get(input$group_complex_dot)))+
      geom_point(size =3, alpha =1) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    limits = c(1,10^6),
                    labels = trans_format("log10", math_format(10^.x))) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    limits = c(1,10^6),
                    labels = trans_format("log10", math_format(10^.x))) +
      theme_bw() +
      scale_color_manual(values=palette.complex) + 
      geom_hline(yintercept = 1000,colour="grey")+
      geom_vline(xintercept = 1000,colour="grey")+
      annotation_logticks()  +
      theme(text=element_text(size=20,family="serif"),
            axis.title = element_text(colour="black", size=40,family="serif"),
            axis.text.x = element_text(colour="black",size=28,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
            axis.text.y = element_text(colour="black",size=28,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
            axis.title.x=element_text(colour="black",size=40,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
            axis.title.y = element_text(colour="black",size=40,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
            legend.title  =element_blank(),
            legend.position = input$legend.dot,
            legend.text = element_text(colour="black",size=input$legend.size.cd,hjust=.5,vjust=.5,face="plain",family="serif")) +
      scale_alpha(guide = 'none') +
      guides(size=FALSE, col = guide_legend(ncol=input$legend.column))+
      labs(x=x_lable1,
           y=y_lable1)
    
    vals15$complex_dot
  })
  dot_plot.complex1 <- function () {
      
      
      if (input$density_dotplot =="no") {
        
        dot_plot.complex()
        
      }
      
      else {
        ggExtra::ggMarginal(dot_plot.complex(),groupColour = TRUE, groupFill = TRUE)
      }
      
    }
  output$dot_plot.complex2 <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })

    dot_plot.complex1()

    
  })
  
  output$downloadPlot_complex.dotplot <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("complex.dotplot_",gsub("/", "-", x), ".pdf", sep = "")
    }, content = function(file) {
      pdf(file, width=input$width_complex.dotplot,height=input$height_complex.dotplot, onefile = FALSE) # open the pdf device
      print(dot_plot.complex1())
      dev.off()}, 
    contentType = "application/pdf" )
  
  output$downloadPlotPNG_complex.dotplot <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("complex.dotplot_", gsub("/", "-", x), ".png", sep = "")
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
      choices=names(input.data2()),
      selected = "AJBJ")
    
  })
  
  observe({
    updateSelectInput(
      session,
      "upset.group.select",
      choices=names(input.data2()),
      selected = "group")
    
  })
  
  
  upset.parameters <- function () {
    df <- input.data2();
    df <- as.data.frame(df)
    head(df)
    unique.df <- unique(df[c(input$upset.select,input$upset.group.select)])
    names(unique.df) <- c("chain","group")
    unique.df$cloneCount <- 1
    mat <- acast(unique.df, chain~group, value.var="cloneCount")
    mat[is.na(mat)] <- 0
    mat <- as.data.frame(mat)
   # a <- as.data.frame(unique(names(mat)))
 #   a$V1 <- distinctColorPalette(dim(a)[1])
    
    df.x <- make_comb_mat(mat)
    
    vals23$upset.plot <- UpSet(df.x, comb_order = order(comb_size(df.x)))
    vals23$upset.plot
    
    
  }

  output$UpSet.plot <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    upset.parameters()
  })
  
  
  
  output$upset.datatable <- DT::renderDataTable(escape = FALSE,filter = 'top', options = list(lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
    df <- input.data2();
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
