library(shiny)
library(dplyr)

ui <- fluidPage(
  fluidPage(
    titlePanel("MY SEQ FILES MERGER"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file1_seq.file",
                  "Choose .seq files from directory",
                  multiple = TRUE,
                  accept=c('.seq')),
        h5("Add Indiv and group/chain name"),
        h6("IndividualID.groupChain-initialwell"),
        fluidRow(
          column(4,selectInput("indiv_miss","Add Indiv label", choices = c("No","Yes"))),
          column(4,selectInput("group_miss","Add Group label", choices = c("No","Yes"))),
          
                        ),
        
       
        fluidRow(column(4,
                        conditionalPanel(
                          condition = "input.indiv_miss == 'Yes'",
                          textInput("indiv.miss.name","Individual ID","Other"),
                          
                        )),
                column(4,
                       conditionalPanel(
                         condition = "input.group_miss == 'Yes'",
                         textInput("group.miss.name","Group ID","Other"),
                         
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
  )
)

server <-  function(input, output) {
  getData <- reactive({
    inFile.seq <- input$file1_seq.file
    if (is.null(inFile.seq)){
      return(NULL)
    }
    
    else {

      numfiles = nrow(inFile.seq) 

      
      df_total = data.frame()
      
      
      for (i in 1:numfiles) { 
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
}

shinyApp(ui = ui, server = server)
