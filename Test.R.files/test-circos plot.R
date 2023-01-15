gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

ui <- navbarPage("ShinyTCR", position = "fixed-top",collapsible = TRUE,
                 
                 tabPanel("Workflow",
                          sidebarLayout(
                            sidebarPanel(id = "tPanel",style = "overflow-y:scroll; max-height: 700px; position:relative;", width=3,
                                         tags$style(type="text/css", "body {padding-top: 70px; padding-left: 10px;}"),
                                         #textInput(inputId = "lab1", label = "Group label of file 1",value = "Ex.vivo"),
                                         tags$head(tags$style(HTML(".shiny-notification {position:fixed;top: 50%;left: 30%;right: 30%;}"))),
                                         tags$head(tags$style(HTML('.progress-bar {background-color: purple;}'))),
                                         selectInput("app.colour","app to selectively colour",choices = c("Treemap","circular"))
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Circular plot", 
                                         h5("Circular plot requires an unsummarised dataset"),
                                         
                                         fluidRow(
                                           
                                           column(2,selectInput( "chain1",label = h5("chain one"),"" )),
                                           column(2,selectInput( "chain2",label = h5("chain two"),"" )),
                                           column(2, selectInput("cir_colour.choise",label = h5("colour default choises"), choices =  c("default","random","grey"))),
                                           column(3,tableOutput("table_display"))
                                         ),
                                         
                                         fluidRow(column(3,wellPanel(id = "tPanel21",style = "overflow-y:scroll; max-height: 600px",uiOutput("myPanel_circ")))),
                                         renderPlot("Circular"),
                                         
                                         
                                         
                                         
                                         
                                         
                                ),
                                tabPanel("bulk TCR: analysis"
                                         
                                )
                              )
                              
                            )
                          )
                          
                 ))

server <- function(input, output, session) {
  
  observe({
    updateSelectInput(
      session,
      "chain1",
      choices=names(input.data2()))
    
  }) # chain 1
  observe({
    updateSelectInput(
      session,
      "chain2",
      choices=names(input.data2()))
    
  }) # chain 2
  
  output$table_display <- renderTable({
    dat <- input.data2();
    dat <- as.data.frame(dat)
    hierarchy <- dat[names(dat) %in% c(input$chain1,input$chain2)]
    
    
    head(hierarchy)
  })
  
  input.data2 <- reactive({
    dat <- read.csv("../ShinyTCR/test-data/Group/Merged_A-B_2021.08.02.csv",header=T) 
    dat <- as.data.frame(dat)
    dat
  })
  
  
  
  cols_circ <- reactive({
    dat <- read.csv("../ShinyTCR/test-data/Group/Merged_A-B_2021.08.02.csv",header=T) 
    dat <- as.data.frame(dat)
    unique.col <- as.data.frame(unique(unlist(dat[names(dat) %in% c("AJ","BJ")])))
    dim(unique.col)
    palette2 <- distinctColorPalette(dim(unique.col)[1])
    
    lapply(1:dim(unique.col)[1], function(i) {
      colourInput(paste("col", i, sep="_"), paste(unique.col[i,]), palette2[i])        
    })
  })
  
  output$myPanel_circ <- renderUI({cols_circ()})
  
  colors_cir <- reactive({
    
    dat <- read.csv("../ShinyTCR/test-data/Group/Merged_A-B_2021.08.02.csv",header=T) 
    hierarchy <- dat[names(dat) %in% c("AJ","BJ")]
    df.col1 <- as.data.frame(unique(hierarchy[,1]))
    names(df.col1) <- "V1"
    df.col.j <- as.data.frame(unique(hierarchy[,2]))
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    dim(df.col.2)
    
    lapply(1:dim(df.col.2)[1], function(i) {
      input[[paste("col", i, sep="_")]]
    })
  })
  
  
  Circular_plot2 <-  reactive( {
    par(mfrow = c(2,2))
    
    dat <- read.csv("../ShinyTCR/test-data/Group/Merged_A-B_2021.08.02.csv",header=T) 
    hierarchy <- dat[names(dat) %in% c("AJ","BJ")]
    df.col1 <- as.data.frame(unique(hierarchy[,1]))
    names(df.col1) <- "V1"
    
    df.col.j <- as.data.frame(unique(hierarchy[,2]))
    names(df.col.j) <- "V1"
    df.col.2 <- rbind(df.col1,df.col.j)
    length(t(df.col.2))
    
    
    col2 <- cols
    col2 <- distinctColorPalette(dim(df.col.2)[1])
    df.col.2$colour <- col2
    grid.col <- as.data.frame(as.matrix(t(as.data.frame(df.col.2$colour))))
    names(grid.col) <- df.col.2$V1
    grid.col <- as.matrix(grid.col)
    
    for (i in 1:2) {
      dat <- read.csv("../ShinyTCR/test-data/Group/Merged_A-B_2021.08.02.csv",header=T) 
      a <- unique(dat$group)
      dat <- subset(dat, dat$group==a[i])
      hierarchy <- dat[names(dat) %in% c("AJ","BJ")]
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
      par(cex=0.5)
      # circos.par("canvas.xlim" = c(-canvas.x.lim, canvas.x.lim), "canvas.ylim" = c(-canvas.y.lim, canvas.y.lim))
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
    }
    
    
  })
  output$Circular <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    Circular_plot2()
  })
  
  
  
  
  
}

shinyApp(ui, server)

##### 
