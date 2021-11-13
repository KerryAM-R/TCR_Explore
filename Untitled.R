library(shiny)
library(shinydashboard)

ui <- fluidPage(
  titlePanel( div(column(width = 6, h2("My Header")), 
                  column(width = 6, tags$img(src = "Logo.png"))),
              windowTitle="MyPage"
  )
)


server <- function(input, output) {
  
}

shinyApp(ui = ui, server = server)
