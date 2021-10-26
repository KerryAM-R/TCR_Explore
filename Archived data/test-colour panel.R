library(shiny)

ui <- pageWithSidebar(
  
  headerPanel("Dynamic number of plots"),
  
  sidebarPanel(
    sliderInput("n", "Number of plots", value=1, min=1, max=5),
    sliderInput("max", "Max number of plots", value=5, min=1, max=25,step=1)
  ),
  
  mainPanel(
    # This is the dynamic UI for the plots
    uiOutput("plots")
  )
)



server = function(input, output,session) {
  
  
  
  # Insert the right number of plot output objects into the web page
  output$plots <- renderUI({
    plot_output_list <- lapply(1:input$n, function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname, height = 280, width = 250)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  
  #making max_plot variable reactiveVal instead of global, input$max would also do, but that would be too easy...
  max_plots = reactiveVal(5)
  observe({
    max_plots(input$max)
    updateSliderInput(session,"n",min=1, max=input$max)
  })
  
  ##wrapping the, for me strange(but nice) for-local-loop, seems to restore normal shiny reactive behavior
  observeEvent({c(input$max,input$n)},{
    mp= max_plots()
    
    for (i in 1:mp) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      local({
        my_i <- i
        plotname <- paste("plot", my_i, sep="")
        
        output[[plotname]] <- renderPlot({
          plot(1:my_i, 1:my_i,
               xlim = c(1, mp),
               ylim = c(1, mp),
               main = paste("1:", my_i, ".  n is ", input$n, sep = "")
          )
        })
      })
    }
    
  })
  
}

shinyApp(ui, server)
