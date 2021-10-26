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
  
  chard.data <- function () {
    
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
  a <- unique(dat$group)
    
  }
  

    uniPlot <- lapply(a, function(i){
      a <- unique(dat$group)
      dat <- chard.data()
      # 
      # dat2 <- subset(dat, dat$group==a[i])
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
      grid.col3 <- as.matrix(grid.col2)
      canvas.x.lim <- 1
      canvas.y.lim <-1
      chordDiagram(hierarchy, annotationTrack = "grid", grid.col = grid.col3,
                   preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(hierarchy))))))
      
      
      
      # title(a[i])
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
      circos.clear()
      # we go back to the first track and customize sector labels
    
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
