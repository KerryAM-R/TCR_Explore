library(shiny)
library(ggplot2)

ui <- shinyUI(fluidPage(
  navbarPage(
    windowTitle = "Automated EDA",
    title = "Automated EDA",
    collapsible = TRUE,
    
    tabPanel(title = "EDA",
             icon = icon("edit"),
             
             fluidRow(column(
               width = 12,
               uiOutput('univariate_plots')
             )))
  )
))

server <- function(input, output,session) {
  
  rv <- reactiveValues()
  
  observe({
    isolate({
      
      dat <- read.csv("../ShinyTCR/test-data/Group/Merged_A-B_2021.08.02.csv",header=T) 
      dat <- subset(dat, dat$group=="IFN")
      hierarchy <- dat[names(dat) %in% c("AJ","BJ")]
      hierarchy$cloneCount <- 1
      hierarchy
      
      chain1 <- as.data.frame(ddply(hierarchy,names(hierarchy)[-c(2,3)],numcolwise(sum)))
      chain1 <- chain1[order(chain1$cloneCount, decreasing = T),]
      chain1
      
      chain2 <- as.data.frame(ddply(hierarchy,names(hierarchy)[-c(1,3)],numcolwise(sum)))
      chain2 <- chain2[order(chain2$cloneCount, decreasing = T),]
      chain2
      
      
      rv$data_filtered <- unique(dat$group)
      
    })
  })
  
  observeEvent(rv$data_filtered, {
    
    
    
    
    
    
    uniPlot <- lapply(names(rv$data_filtered), function(i){
      
      dat <- read.csv("../ShinyTCR/test-data/Group/Merged_A-B_2021.08.02.csv",header=T) 
      hierarchy <- dat[names(dat) %in% c("AJ","BJ")]
      df.col1 <- as.data.frame(unique(hierarchy[,1]))
      names(df.col1) <- "V1"
      df.col.j <- as.data.frame(unique(hierarchy[,2]))
      names(df.col.j) <- "V1"
      df.col.2 <- rbind(df.col1,df.col.j)
      length(t(df.col.2))
      col2 <- distinctColorPalette(dim(df.col.2)[1])
      df.col.2$colour <- col2
      grid.col <- as.data.frame(as.matrix(t(as.data.frame(df.col.2$colour))))
      names(grid.col) <- df.col.2$V1
      grid.col <- as.matrix(grid.col)
      
      a <- rv$data_filtered
      
      dat2 <- subset(dat, dat$group==a[i])
      hierarchy <- dat2[names(dat2) %in% c("AJ","BJ")]
      df.col1 <- as.data.frame(unique(hierarchy[,1]))
      names(df.col1) <- "V1"
      df.col.j <- as.data.frame(unique(hierarchy[,2]))
      names(df.col.j) <- "V1"
      df.col.2 <- rbind(df.col1,df.col.j)
      df.col.2
      grid.col1 <- as.data.frame(grid.col)
      grid.col2 <- grid.col1[names(grid.col1) %in% df.col.2$V1]
      grid.col2
      grid.col <- as.matrix(grid.col2)
      
      chordDiagram(hierarchy, annotationTrack = "grid", grid.col = grid.col,
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(hierarchy))))))
      title(a[i], cex = 0.6)

    })
    
    output$univariate_plots <- renderUI({
      
      plot_output_list <- lapply(seq_along(1:length(uniPlot)),function(plot_i){
        column(width = 6,
               tags$div(style = "margin-top: 10px; margin-bottom: 10px;",
                        plotOutput(outputId = paste0("uni_", plot_i))
               ))
      })
      
      # Convert the list to a tagList - this is necessary for the list of items to display properly.
      do.call(tagList, plot_output_list)
      
      ## either works
      # plot_output_list
      
    })
    
    rv$uniPlot <- uniPlot
  })
  
  observeEvent(rv$uniPlot,{
    
    lapply(seq_along(1:length(rv$uniPlot)), function(plot_i) {
      output[[paste0("uni_", plot_i)]] <- renderPlot({
        rv$uniPlot[[plot_i]]
      })
    })
    
  }, ignoreInit = FALSE)
  
}

shinyApp(ui, server)
