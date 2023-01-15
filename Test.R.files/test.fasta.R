runApp_TCR_EXPLORE_V1()

# QC ----

ui <- navbarPage("page",
          navbarMenu("QC",
           tags$head(
           
             )
           ),
           # Create and merge 50 fasta files -----
           tabPanel("Making fasta files",
                    directoryInput('directory', label = 'select a directory'),
                    verbatimTextOutput("dir", placeholder = TRUE),
                    actionButton("do", "Click Me to make fasta file (50 per file)"),
                    h4('Merging statistics'),
                    uiOutput('textWithHTML') # ui output as a list of HTML p() tags
           ))

server  <- function(input, output, session) {


# making 50 fasta files -----
path = reactive({readDirectoryInput(session, 'directory')})
observeEvent(
  ignoreNULL = TRUE,
  eventExpr = {
    input$directory
  },
  handlerExpr = {
    if (input$directory > 0) {
      # condition prevents handler execution on initial app launch
      
      # launch the directory selection dialog with initial path read from the widget
      path = choose.dir(default = readDirectoryInput(session, 'directory'))
      # update the widget value
      updateDirectoryInput(session, 'directory', value = path)
      setwd(path)
      observeEvent(input$do, {
        
      })
      
    }
  }
)

observeEvent(input$do, { output$textWithHTML <- renderUI({
  system(command = system.file("extdata","test-data/scripts/alignment_211230.sh", package = "TCR.Explore"))
  rawText <- readLines('time.txt') # get raw text
  # split the text into a list of character vectors
  #   Each element in the list contains one line
  splitText <- stringi::stri_split(str = rawText, regex = '\\n')
  
  # wrap a paragraph tag around each element in the list
  replacedText <- lapply(splitText, p)
  
  return(replacedText)
  
})
})

} # final bracket

shinyApp(ui, server)


