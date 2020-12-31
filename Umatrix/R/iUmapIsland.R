iUmapIsland <- function(Umatrix=NULL, BestMatches=NULL, Cls=NULL){
  # Imx = iUmapIsland(Umatrix, BestMatches, Cls)
  #
  # INPUT
  # Umatrix
  # BestMatches
  # Cls
  # OUTPUT
  # island          the generated Imx
  # Author: FL

  #requireRpackage("shiny")
  #requireRpackage("shinyjs")
  #requireRpackage("png")
  #requireRpackage("tcltk")
  FilePath<-getwd() #MT: sonst funktionierts nicht
  DataName = NULL
  UmatrixSize = 1

  QuitButtonText = "Quit"

  ########
  ########

  ############
  # Fehler abfangen
  ############
  if(is.null(Umatrix))
    stop("Missing Umatrix")
  if(!is.null(Cls)){
    if(is.list(Cls)) stop("Cls is not a vector")
    if(nrow(BestMatches)!=length(Cls))
      stop('The length of the given classification does not equal the row length of the bestmatches')
  }
  # umatrix normieren
  Umatrix = Umatrix/max(Umatrix)

  if(!is.null(BestMatches))
    if(is.null(Cls)) Cls=rep(1,nrow(BestMatches))

  Imx = matrix(0,nrow=nrow(Umatrix)*2,ncol=ncol(Umatrix)*2)
  islandLeft = 1
  islandRight = ncol(Umatrix)*2
  islandTop = 1
  islandBottom = nrow(Umatrix)*2

  Width = islandRight - islandLeft + 1
  Height = islandBottom - islandTop + 1

  # Vorschlag f?r eine Insel berechnen
  idealIsland = NULL
  if(!is.null(BestMatches))
    idealIsland = bestUmatrixTranslation(Umatrix, BestMatches)
  ##########
  # Shiny Fenster
  ##########
  UmatrixUi = fluidPage(
    useShinyjs(),
    sidebarLayout(position="right",
      mainPanel(
        #div(plotOutput("UmatrixPlot", width=as.character(Width), height=as.character(Height), click = "clickOnUmatrix"), style="margin-left:-70px")
        div(plotOutput("UmatrixPlot", click = "clickOnUmatrix"))
      ),
      div(style="max-width:1150px", # die sidebar ist dann auf 33% davon beschraenkt
      sidebarPanel(
        # busy balken
        tags$head(tags$style(type="text/css", "
                             #loadmessage {
                             position: fixed;
                             top: 0px;
                             left: 0px;
                             width: 100%;
                             padding: 5px 0px 5px 0px;
                             text-align: center;
                             font-weight: bold;
                             font-size: 100%;
                             color: #000000;
                             background-color: #CCFF66;
                             z-index: 105;}")),
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),

        fluidRow(
        actionButton("createIsland", "Create Island"),
        actionButton("resetIsland", "Reset Island")),


        fluidRow(
          column(12, selectInput("BmSize", "Bestmatchsize:",
                                 c(
                                   "1x" = 1,
                                   "2x" = 2,
                                   "3x" = 3,
                                   "4x" = 4,
                                   "6x" = 6,
                                   "8x" = 8),
                                 selected = 2))
        ),
        fluidRow(
          column(12, selectInput("UmxSize", "Umatrixsize:",
                                 c(
                                   "0.5x" = 0.5,
                                   "1x" = 1,
                                   "1.5x" = 1.5,
                                   "2x" = 2,
                                   "3x" = 3,
                                   "4x" = 4),
                                 selected = UmatrixSize))
        ),

        checkboxInput("ShowWarnings", "Show Warnings", value = TRUE),

        htmlOutput("warnings"),

        br(),

        fluidRow(
          actionButton("quit", QuitButtonText)
          )
      )
    )))


  BestMatches = CheckBestMatches(BestMatches, Cls, shiny=T)
  CheckUmatrix(Umatrix, shiny=T)

  outputApp=runApp(list(ui = UmatrixUi, server = function(input, output, session){

    updateSelectInput(session, "UmxSize", selected = UmatrixSize)
    observe({UmatrixSize <<- as.numeric(input$UmxSize)})

    val = reactiveValues()
    val$UmatrixImg = NULL
    val$currentLines = c()
    val$BestMatches = BestMatches
    val$Cls = Cls
    val$Imx = matrix(0,nrow=nrow(Umatrix)*2,ncol=ncol(Umatrix)*2)

    val$Width <- Width
    val$Height <- Height

    val$Stretchfactor = 800/Width

    observe({
      if(input$ShowWarnings)
        output$warnings <- renderText(val$warnings)
      else
        output$warnings <- renderText(" ")
    })


    val$warnings = "<font color=\"#00FF00\"><b>No Warnings</b></font>"


    ##########
    # rerender Umatrix if island has changed
    ##########
    observe({
      CurrentDirectory <- getwd()
      setwd(tempdir())
      png("tmpUmatrix.png", width= val$Stretchfactor*val$Width,
          height = val$Stretchfactor*val$Height  )


      islandExists = !(sum(val$Imx)==0)


      umx <- plotMatrix(Umatrix, BmSize = as.numeric(input$BmSize)*2, BestMatches = val$BestMatches, Cls = val$Cls, Clean = T,
                        Imx=val$Imx, RemoveOcean = T, TransparentContours = F,
                        Toroid = T, MarkDuplicatedBestMatches=(islandExists&input$ShowWarnings))

      if(!is.null(val$BestMatches))
        idealIsland = bestUmatrixTranslation(Umatrix, val$BestMatches)

      if((!is.null(idealIsland))&&(!islandExists)){
        test <- data.frame(x=c(idealIsland$cols,
                               idealIsland$cols,
                               idealIsland$cols+ncol(Umatrix),
                               idealIsland$cols+ncol(Umatrix),
                               idealIsland$cols),
                           y=c(idealIsland$lines,
                               idealIsland$lines+nrow(Umatrix),
                               idealIsland$lines+nrow(Umatrix),
                               idealIsland$lines,
                               idealIsland$lines))

        # zeichne standard insel drueber
        umx = umx + geom_polygon(data=test,
                              aes(x=x, y=y),
                              size=0.5, alpha=0.01,
                              color="yellow", fill=NA)
      }

      print(umx)

      dev.off()
      val$UmatrixImg <- readPNG("tmpUmatrix.png")
      setwd(CurrentDirectory)
    })

    ##########
    # plot Umatrix
    ##########
    observe({
      output$UmatrixPlot <- renderPlot({
        if(!is.null(val$UmatrixImg)){
          par(mar=c(0,0,0,0))
          par(mai=c(0,0,0,0))
          par(xaxs = 'i',yaxs='i')
          rasterImage(val$UmatrixImg,0,0,1,1)

          lines(val$currentLines[,1], val$currentLines[,2],col="red")
        }
      },
      width = function(){val$Width * as.numeric(input$UmxSize) * val$Stretchfactor},
      height = function(){val$Height * as.numeric(input$UmxSize) * val$Stretchfactor})
    })

    ##########
    # react on click
    ##########
    observeEvent(input$clickOnUmatrix,{
      x <- isolate(input$clickOnUmatrix$x)
      y <- isolate(input$clickOnUmatrix$y)

      val$currentLines = rbind(val$currentLines,c(x,y))
    })

    #########
    # react on button: create island
    #########
    observeEvent(input$createIsland,{
      # check if there are already enough points to draw a complete polygon
      if(is.null(val$currentLines)) return()
      if(nrow(val$currentLines) < 3) return()

      Polygon = cbind( val$currentLines[,1]*ncol(Umatrix)*2+0.5, (1-val$currentLines[,2])*nrow(Umatrix)*2+0.5 )
      val$Imx <- createMaskFromPolygon(nrow(Umatrix)*2, ncol(Umatrix)*2, CutoutPol = Polygon)
      Imx <<- val$Imx
      val$currentLines = c()
      islandBM = UmxBestMatchesFromIsland(Umatrix, val$BestMatches, Imx = val$Imx, Toroid = T)$BestMatches

      if(!is.null(val$BestMatches)){
        if(length(unique(islandBM[,1])) < nrow(val$BestMatches)){
          val$warnings = paste0("<font color=\"#FF0000\"><b>",nrow(val$BestMatches)-length(unique(islandBM[,1]))," bestmatches are missing!</b></font>")
        }
        else if(length(islandBM[,1]) > nrow(val$BestMatches)){
          val$warnings = "<font color=\"#FF0000\"><b>Some bestmatches are duplicated! (marked with red triangles)</b></font>"
        }
        else{
          val$warnings = "<font color=\"#00FF00\"><b>No Warnings</b></font>"
        }
      }
    })

    #######
    # react on button: reset island
    #######
    observeEvent(input$resetIsland,{
      val$Imx = matrix(0,nrow=nrow(Umatrix)*2, ncol=ncol(Umatrix)*2)
      val$currentLines = c()
      val$warnings = "<font color=\"#00FF00\"><b>No Warnings</b></font>"
    })


    ##########
    ##########

    session$onSessionEnded(function() {
      print("program closed")
      stopApp(list(Imx = Imx))
    })

    observeEvent(input$quit, {
      stopApp(list(Imx = Imx))
    })

  }))

  return(outputApp)
}
