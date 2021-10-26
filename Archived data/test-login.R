library(shiny)
library(shinydashboard)
library(DT)
library(shinyjs)
library(sodium)
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
require("muscle")
require("DiffLogo")
library("shinyauthr")




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


# Main login screen ----
loginpage <- div(id = "loginpage", style = "width: 500px; max-width: 100%; margin: 0 auto; padding: 20px;",
                 wellPanel(
                   tags$h2("LOG IN", class = "text-center", style = "padding-top: 0;color:#333; font-weight:600;"),
                   textInput("userName", placeholder="Username", label = tagList(icon("user"), "Username")),
                   passwordInput("passwd", placeholder="Password", label = tagList(icon("unlock-alt"), "Password")),
                   br(),
                   div(
                     style = "text-align: center;",
                     actionButton("login", "SIGN IN", style = "color: white; background-color:#3c8dbc;
                                 padding: 10px 15px; width: 150px; cursor: pointer;
                                 font-size: 18px; font-weight: 600;"),
                     shinyjs::hidden(
                       div(id = "nomatch",
                           tags$p("Oops! Incorrect username or password!",
                                  style = "color: red; font-weight: 600; 
                                            padding-top: 5px;font-size:16px;", 
                                  class = "text-center"))),
                     br(),
                     br(),
                     tags$code("Username: myuser  Password: mypass"),
                     br(),
                     tags$code("Username: myuser1  Password: mypass1")
                   ))
)

credentials = data.frame(
  username_id = c("admin"),
  passod   = sapply(c("123"),password_store),
  stringsAsFactors = F
)

header <- dashboardHeader(title = "TCR_Explore", uiOutput("logoutbtn"))

sidebar <- dashboardSidebar(uiOutput("sidebarpanel")) 
body <- dashboardBody(shinyjs::useShinyjs(), uiOutput("body"))
ui<-dashboardPage(header, sidebar, body, skin = "blue")

server <- function(input, output, session) {
  
  login = FALSE
  USER <- reactiveValues(login = login)
  observe({ 
    if (USER$login == FALSE) {
      if (!is.null(input$login)) {
        if (input$login > 0) {
          Username <- isolate(input$userName)
          Password <- isolate(input$passwd)
          if(length(which(credentials$username_id==Username))==1) { 
            pasmatch  <- credentials["passod"][which(credentials$username_id==Username),]
            pasverify <- password_verify(pasmatch, Password)
            if(pasverify) {
              USER$login <- TRUE
            } else {
              shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
            }
          } else {
            shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
            shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
          }
        } 
      }
    }
  })
  

  
 
  
  
  
  output$logoutbtn <- renderUI({
    req(USER$login)
    tags$li(a(icon("th"), "Logout", 
              href="javascript:window.location.reload(true)"),
            class = "dropdown", 
            style = "background-color: #eee !important; border: 0;
                    font-weight: bold; margin:5px; padding: 10px;")
  })
  
  output$sidebarpanel <- renderUI({
    if (USER$login == TRUE ){ 
      sidebarMenu(id = 'tab',
        menuItem("Workflow", tabName = "1", icon = icon("list-alt")),
        menuItem("QC", tabName = "2", icon = icon("th")),
        menuItem("scTCR", tabName = "3", icon = icon("th"))
      )
    }
  })
  
  output$body <- renderUI({
    if (USER$login == TRUE ) {
        
        if (input$tab == "1") {
          
          dyn_ui <- tabsetPanel(id = "tabset_id", selected = "t1", 
                                tabPanel("Overview", value = "t1",
                                         fluidRow(includeMarkdown("README.md"))),
                                tabPanel("Quality control", value = "t2",
                                         fluidRow(includeMarkdown("READMEQC.md"))),
                                tabPanel("scTCR analysis", value = "t3",
                                         fluidRow(includeMarkdown("README.scTCR.md"))),
                                tabPanel("FACS index analysis", value = "t4",
                                         fluidRow(includeMarkdown("README.FACS.md"))),
                                        )
        } 
        if (input$tab == "2") {
          
          dyn_ui <- tabPanel("IMGT",
                             sidebarLayout(
                               sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                            selectInput("dataset_IMGT3", "Choose a dataset:", choices = c("ab-test-data1","gd-test-data1", "own_data")),
                                            fileInput('file_IMGT3', 'Select file for IMGT datafile',
                                                      accept=c('xls/xlsx', '.xls')),
                                            selectInput("IMGT_chain2","Alpha-beta or gamma-delta",choices = c("ab","gd")),
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
                                            actionButton("update.string","update string"),
                                            fluidRow(column(12, selectInput("string.data3","column names for summary","",multiple = T, width = "1200px") )),
                                            tags$head(tags$style("#chain_table_IMGT.QC3  {white-space: nowrap;  }")),
                                            div(DT::dataTableOutput("chain_table_IMGT.QC3")),
                                            downloadButton('downloadTABLE.QC3','Download table')
                                            
                                   ) 
                                   
                                 )
                               )
                               
                             )
          )
        }
      
      
       if (input$tab == "3") {
         dyn_ui <- tabPanel("scTCR")
         
         
       }
      
        return(dyn_ui)
      
      
      
    }
    else {
      loginpage
    }
  })
  
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

  chain_merge_IMGTonly <- function (){
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

  }
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
    df <- input.data.IMGT_afterQC()
    df <- as.data.frame(df)
    names(df)

  })
  #
  observeEvent( input$update.string, {
    df <- isolate(chain_merge_IMGTonly())

        updateSelectInput(
          session,
          "string.data3",
          choices=names(df),
          selected = c("group","JUNCTION_A"))
  })
chain_table_summary <- reactive({
  df <- chain_merge_IMGTonly()
  df <- as.data.frame(df)
  df2 <- df[,c("cloneCount",input$string.data3)]
  df2
  df3 <- as.data.frame(ddply(df2,input$string.data3,numcolwise(sum)))
  df3
})




  output$chain_table_IMGT.QC3 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    df1 <- input.data.IMGT_afterQC()
    df1 <- as.data.frame(df1)

    a <- subset(df1 ,is.na(df1$clone_quality)==TRUE)
    if (dim(a)[1]>0) {
      df <- as.data.frame("please complete QC analysis")
      names(df) <- " "
      df

    }
    else {

       chain_table_summary()
    }
  })


  output$downloadTABLE.QC3 <- downloadHandler(
    filename = function(){
      paste("paired_chain_CDR3",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- chain_merge_IMGTonly()
      df <- as.data.frame(df)
      df2 <- df[,c("cloneCount",input$string.data3)]
      df2
      df3 <- as.data.frame(ddply(df2,input$string.data3,numcolwise(sum)))
      df3
      write.csv(df3,file, row.names = FALSE)
    })
  
  
  
  
  
  
}

runApp(list(ui = ui, server = server))
