max_plots <- 5

ui <- fluidPage(
  
  headerPanel("Dynamic number of plots"),
  
  sidebarPanel(
    sliderInput("n", "Number of plots", value=1, min=1, max=5)
  ),
  
  mainPanel(
    # This is the dynamic UI for the plots
    uiOutput("plots.circular.multiple")
  )
)

server <- function(input, output) {
  options(error=NULL)
  # Insert the right number of plot output objects into the web page
  output$plots.circular.multiple <- renderUI({
    set.seed(123)
    plot_output_list <- lapply(1:input$n, function(i) {
      dat <- read.csv("../ShinyTCR/test-data/Group/Merged_A-B_2021.08.02.csv",header=T) 
      a <- unique(dat$group)
      
      plotname <- paste("plot", a[i], sep="")
      plotOutput(plotname, height = 500, width = 500)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  for (i in 1:max_plots) {
    set.seed(123)
    
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      
      dat <- read.csv("../ShinyTCR/test-data/Group/Merged_A-B_2021.08.02.csv",header=T) 
      a <- unique(dat$group)
      my_i <- i
      plotname <- paste("plot", a[i], sep="")
      
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
      b <- a[i]
      output[[plotname]] <- renderPlot({
        par(cex=0.5)
        circos.par("canvas.xlim" = c(-1.25, 1.25), "canvas.ylim" = c(-1, 1))
        chordDiagram(hierarchy, annotationTrack = "grid", grid.col = grid.col,
                     preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(hierarchy))))))
        title(main = b, cex = 0.6)
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
      })
    })
  }
}

shinyApp(ui, server)
