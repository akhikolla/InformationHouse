iUstarmatrix <- function(Weights = NULL, Lines = NULL, Columns = NULL, Data=NULL ,Imx = NULL, Cls=NULL, Radius = NULL,
                         Toroid = T){
  # V = iUstarmatrix(Umatrix, BestMatches, Cls)
  # UstarMatrix = V$UstarMatrix
  # PMatrix = V$PMatrix
  #
  # INPUT
  # Weights     The Weights in List Format
  # Data            Data that was used to learn the EsomNeurons
  # BestMatches
  # Cls
  # Radius          Radius for selecting neighbours and generating the pmatrix
  # OUTPUT
  # list of:
  # UstarMatrix
  # PMatrix/
  # Author: FL & MT

  DataName = NULL
  UmatrixSize = 1

  QuitButtonText = "Quit"

  ########
  ########


  EsomNeurons = Weights

  ############
  # Fehler abfangen
  ############
  if(is.null(Weights))
    stop("Missing Weights")
  if(is.null(Data))
    stop("Missing Data (lrn)")
  if(is.null(Lines))
    stop("Missing Lines")
  if(is.null(Columns))
    stop("Missing Columns")

  ##############
  # Bestmatches aus Weights und Data generieren
  ##############
  V = esomProjection(Weights, Data, Columns, Lines)
  BestMatches = cbind(nrow(V), V)

  if(!is.null(Cls)){
    if(!is.vector(Cls)) stop("Cls is not a vector")
    if(nrow(BestMatches)!=length(Cls))
      stop('The length of the given classification does not equal the number of data rows')
  }

  # if(!is.null(BestMatches))
  #   if(is.null(Cls)) Cls=rep(1,nrow(BestMatches))
  # BestMatches = CheckBestMatches(BestMatches, Cls, shiny=T)

  origData = nrow(Data)

  Width = Columns
  Height = Lines

  ##########
  # Bestimme alle Distanzen und schlage radius durch EM vor
  ###########
  distancesAsMatrix = as.matrix(stats::dist(Data))
  distances=distancesAsMatrix[lower.tri(distancesAsMatrix, diag = FALSE)]
  PradiusMinimum = signif(min(distances), digits=2)
  PradiusMaximum = signif(max(distances), digits=2)
  RadiusByEM = NULL

  V = Delaunay4BestMatches(BestMatches, MatrixOrSize = c(Lines,Columns), IsToroid = Toroid)
  Delaunay = V$Delaunay
  ToroidDelaunay = V$ToroidDelaunay

  DelaunayDistances = Delaunay*distancesAsMatrix # nur distanzen zwischen nachbarn
  neighbourDistances = DelaunayDistances = DelaunayDistances[DelaunayDistances!=0]

  RadiusByEM = 1
  tryCatch({
    V = EMGauss(neighbourDistances, 2)
    V <- BayesDecisionBoundaries(V$Means, V$SDs, V$Weights)
    if(length(V)>1){
      print("Multiple Decision boundaries found. max(Distances)/2 chosen a Radius")
      RadiusByEM <- signif(max(distances)/2,2)
    }
    else RadiusByEM <- signif(V * 0.8, 2)
  }, error = function(e){
      print("No Decision boundaries found. max(Distances)/2 chosen a Radius")
      RadiusByEM <- signif(max(distances)/2,2)
  })


  #########
  # initialisiere alles
  ###########
  if((!is.null(Radius))){
    if(PradiusMinimum > Radius) stop(paste0("Minimum for Radius is ", PradiusMinimum,". Given Radius is ", Radius))
    if(PradiusMaximum < Radius) stop(paste0("Maximum for Radius is ", PradiusMaximum,". Given Radius is ", Radius))
    Pradius <- Radius
  }
  else{
    Pradius <- RadiusByEM
  }


  Umatrix <- umatrixForEsom(Weights, Lines, Columns, Toroid=Toroid)

  UstarMatrix = NULL
  PMatrix = NULL

  mapGridText = 'Planar'
  if(Toroid) mapGridText = 'Toroid'


  ##########
  # Shiny Fenster
  ##########
  UmatrixUi = fluidPage(
    useShinyjs(),
    sidebarLayout(position="right",
      mainPanel(
        div("U*Matrix", style="font-size:16pt"),
        plotOutput("UStarmatrixPlot", height = "100%"),
        div("UMatrix", style="font-size:16pt;"),
        plotOutput("UmatrixPlot", height = "100%"),

        style = "text-align: left; width = 100%"
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


        tabsetPanel(
          tabPanel("Radius",
            fluidRow(
              column(6,numericInput("newPradius", "Pradius", value = Pradius)),
              column(6,actionButton("submitPradius", "Submit"), style="margin-top:24px")),

            fluidRow(
              column(6,div(style="font-weight: bold;",
                "Suggestion By EM: ", textOutput("RadiusByEMText", inline=T))),
                column(6, actionButton("resetRadiusByEM", "Reset"))),

            plotOutput("DistNeighHistPlot", click="distneighhistclick"),
            "PMatrix",
            plotOutput("PmatrixPlot")
            #plotOutput("DistNeighHistPlot", click="distneighhistclick")
            #plotOutput("DistHistPlot", click="disthistclick"),
          ),

          tabPanel("View",
            selectInput('mapGrid','Grid Topology',choices=c('Toroid','Planar'),selected=mapGridText),
            fluidRow(
              column(12, selectInput("BmSize", "Bestmatchsize:",
                                     c(
                                       "1x" = 1,
                                       "2x" = 2,
                                       "3x" = 3,
                                       "4x" = 4,
                                       "6x" = 6,
                                       "8x" = 8),
                                     selected = 2))),
            fluidRow(
              column(12, selectInput("UmxSize", "Umatrixsize:",
                                     c(
                                       "0.5x" = 0.5,
                                       "1x" = 1,
                                       "1.5x" = 1.5,
                                       "2x" = 2,
                                       "3x" = 3,
                                       "4x" = 4),
                                     selected = UmatrixSize))),

            fluidRow(column(12,checkboxInput("TransparentContours", "transparent contours", value=F))),
            fluidRow(column(12,checkboxInput("showBMU", "show bestmatches", value=T)))
          ) # tabPanel View
        ), # tabSetPanel

        br(),br(),br(),

        fluidRow(
          actionButton("quit", QuitButtonText))
        )# sidebarpanel
    )))

  outputApp=runApp(list(ui = UmatrixUi, server = function(input, output, session){

    updateSelectInput(session, "UmxSize", selected = UmatrixSize)
    observe({UmatrixSize <<- as.numeric(input$UmxSize)
    })

    val = reactiveValues()
    val$UmatrixImg = NULL
    val$PmatrixImg = NULL
    val$PMatrix <- NULL
    val$BestMatches <- BestMatches
    val$Cls <- Cls
    val$PradiusLimitsCalculated <- F
    val$Imx = NULL
    val$Distances = NULL
    val$Pradius = Pradius

    val$Toroid = Toroid
    observe({val$Toroid = (input$mapGrid == "Toroid")})

    val$Width = Width
    val$Height = Height
    val$Stretchfactor = 1

    # initialize empty global variables
    PradiusOut <- 0

    ###
    # update size
    ###
    observe({
      if(!is.null(val$PMatrix)){
        if(val$Toroid){
          val$Width = ncol(val$PMatrix) * 2
          val$Height = nrow(val$PMatrix) * 2
        }
        else{
          val$Width = ncol(val$PMatrix)
          val$Height = nrow(val$PMatrix)
        }
      }

      # zooming
      val$Stretchfactor = 800/val$Width
    })

    ##########
    # recalculate and rerender u*matrix
    ##########
    observe({
      output$UStarmatrixPlot <- renderPlot({
        if(is.null(val$PMatrix)) return()

        UstarMatrix <<- ustarmatrixCalc(Umatrix,val$PMatrix)

        if(input$showBMU) BMU = val$BestMatches
        else BMU = NULL

        plotMatrix(UstarMatrix, TransparentContours = input$TransparentContours, BmSize = as.numeric(input$BmSize), BestMatches = BMU,
                      Toroid = val$Toroid, Cls = val$Cls, Clean = T,DrawLegend = T, Imx = Imx, RemoveOcean=T)

      }, width=val$Width*val$Stretchfactor*as.numeric(input$UmxSize), height=val$Height*val$Stretchfactor*as.numeric(input$UmxSize))
    })

    #########
    # render umatrix
    #########
    observe({
      output$UmatrixPlot <- renderPlot({
        if(input$showBMU) BMU = val$BestMatches
        else BMU = NULL
        plotMatrix(Umatrix, TransparentContours = input$TransparentContours, BmSize = as.numeric(input$BmSize), BestMatches = BMU,
                      Toroid = val$Toroid, Cls = val$Cls, Clean = T,DrawLegend = T, Imx = Imx, RemoveOcean=T)
      }, width=val$Width*val$Stretchfactor*as.numeric(input$UmxSize), height=val$Height*val$Stretchfactor*as.numeric(input$UmxSize))
    })

    ######
    # recalculate pmatrix
    ######
    observe({
      if(val$Pradius < 0.1)
        val$PMatrix <- pmatrixForEsom(Data, Weights=Weights,Lines=Lines, Columns=Columns, Radius=0.1)
      else
        val$PMatrix <- pmatrixForEsom(Data, Weights=Weights,Lines=Lines, Columns=Columns, Radius=val$Pradius)
      PMatrix <<- isolate(val$PMatrix)
      PradiusOut<<-isolate(val$Pradius)
    })

    ##########
    # rerender pmatrix
    ##########
    output$PmatrixPlot <- renderPlot({
       if(input$showBMU) BMU = val$BestMatches
       else BMU = NULL
       V <- plotMatrix(val$PMatrix, ColorStyle = "Pmatrix", BmSize = as.numeric(input$BmSize), Toroid=val$Toroid, Nrlevels = 10,
                 DrawLegend=F, Imx = Imx, TransparentContours=F, FixedRatio=T, Cls=val$Cls, Clean=T, BestMatches = BMU, RemoveOcean = T)
       return(V)
    })


    # output$DistHistPlot <- renderPlot({
    #   #distances <- pdist(Hepta$Data)
    #   GGDistancesPlot <- PDEplot(as.vector(distances), title="PDE Plot of distances")$ggPlot
    #   GGDistancesPlot = GGDistancesPlot + geom_vline(xintercept = val$Pradius) +
    #     geom_vline(xintercept = RadiusByEM, colour="violet")
    #   GGDistancesPlot
    # })

    output$DistNeighHistPlot <- renderPlot({
      if(is.null(BestMatches)) return()
      
      Data = neighbourDistances
      # simple pde plot
      # radius
      distvec = (stats::dist(as.matrix(Data)))
      paretoRadius = quantile(distvec, 18/100, type = 8, na.rm = TRUE)
      paretoRadius = paretoRadius * 4 / length(Data)^0.2
      # no of bins
      iqr = quantile(Data, 0.75) - quantile(Data, 0.25)
      obw = 3.49 * (min(sd(Data), iqr/1.349) / nrow(Data)^(1/3))
      noOfBins = max((max(Data) - min(Data)) / obw, 10)
      # density
      kernels = as.vector(seq(min(Data), max(Data), length.out = noOfBins))
      inRad = sapply(kernels, function(k)   (Data >= (k - paretoRadius)) & (Data <= (k + paretoRadius)))
      nrInRad = as.vector(colSums(inRad))
      normv = pracma::trapz(kernels, nrInRad)
      dens = nrInRad / normv
      # plot itself
      GGDistancesPlot = ggplot2::ggplot(data  = data.frame(x = kernels, y = dens), ggplot2::aes_string("x", "y")) + ggplot2::geom_line() + ggplot2::ggtitle("Distances of Delaunay Graph")
      
      GGDistancesPlot = GGDistancesPlot + geom_vline(xintercept = val$Pradius) +
        geom_vline(xintercept = RadiusByEM, colour="violet")
      GGDistancesPlot
    })

    output$RadiusByEMText <- renderText(RadiusByEM)
    observeEvent(input$resetRadiusByEM, {
      val$Pradius <- RadiusByEM
      #updateSliderInput(session, "Pradius", value=RadiusByEM)
    })
    observeEvent(input$submitPradius, {
      #updateSliderInput(session, "Pradius", value=input$newPradius)
      val$Pradius <- input$newPradius
    })
    observe({
      updateNumericInput(session, "newPradius", value=val$Pradius)
    })

    #########
    # clicks on histograms (new pradius)
    #########
    # observeEvent(input$disthistclick,{
    #   if(!is.null(input$disthistclick$x))
    #     val$Pradius <- signif(input$disthistclick$x,2)
    # })
    observeEvent(input$distneighhistclick,{
      if(!is.null(input$distneighhistclick$x))
        val$Pradius <- signif(input$distneighhistclick$x,2)
    })

    ############
    ############

    session$onSessionEnded(function() {
      print("program closed")
      stopApp(list(UstarMatrix = UstarMatrix, PMatrix = PMatrix))
    })

    observeEvent(input$quit, {
      stopApp(list(UstarMatrix = UstarMatrix, PMatrix = PMatrix))
    })

  }))


  return(outputApp)
}
