iEsomTrain <- function(Data=NULL, BestMatches=NULL, Cls=NULL, Key = NULL, Toroid = T){
  # V = iEsomTrain()
  # Umatrix <- V$Umatrix
  # BestMatches <- V$BestMatches
  # Lines <- V$Lines
  # Columns <- V$Columns
  # Epochs <- V$Epochs
  # Weights <- V$Weights
  #
  # INPUT
  # Data            Data that will be used to learn
  # BestMatches
  # Cls
  # Key            names of the datapoints
  # OUTPUT
  # list of:
  # Umatrix
  # BestMatches
  # Lines           Lines that where used for training
  # Columns
  # Epochs
  # Weights         The EsomNeurons in Listformat
  # Author: FL


  QuitButtonText = "Quit"

  fixedBestMatches = FALSE
  Umatrix = NULL
  BestMatches = NULL
  FilePath = NULL
  DataName = NULL
  EsomDetails<-list()

  if(!is.null(BestMatches)) fixedBestMatches = TRUE

  if(!is.null(Data)){
    if(!is.matrix(Data)) stop("Data is not given as a matrix!")
    if(is.null(Key)) warning("Data but no keys given. \n1:nrow(Data) will be used as key.\nClassification may be broken.\nProceed with care!")
  }

  ########
  ########

  ########
  # Fehler abfangen
  ########
  if(is.null(FilePath))
    FilePath = getwd()
  if(!is.null(Cls)){
    if(nrow(Data)!=length(Cls))
      stop('The length of the given classification does not equal the row length of the data set')
  }

  if(!is.null(BestMatches)){
    if(!is.matrix(BestMatches)){
      Data=as.matrix(BestMatches)
      warning('BestMatches should be a matrix, Umatrix() called as.matrix(). If an follwing error occurs, please change BestMatches to matrix format manually.')
  }
   if(nrow(Data)!=nrow(BestMatches))
      stop('The row length of the BestMatches does not equal the row length of the data set')
  }

  mapGridText = 'Planar'
  if(Toroid) mapGridText = 'Toroid'

  ##########
  # Shiny Fenster
  ##########
  UmatrixUi = fluidPage(
    sidebarLayout(position="right",
	    mainPanel(
	      plotOutput("plotMatrix", width = "800", height = "500")),
	    div(style="max-width:1150px", # die sidebar ist dann auf 33% davon beschraenkt
	    sidebarPanel(      # U-Matrix training Panel # busy balken

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


        id="trainingPanel",
	        numericInput('trainingEpochs','Number of training epochs',value=20),
	        fluidRow(
              selectInput('mapGrid','Grid Topology',choices=c('Toroid','Planar'),selected=mapGridText),
              fluidRow(
                column(6,numericInput('numRows','Lines',value=50)),
                column(6,numericInput('numColumns','Columns',value=80))
              )
          ),
		     checkboxInput("excalc", "Expert Modus"),

		     conditionalPanel(condition = "input.excalc==true",
  	       fluidRow(
  	         selectInput('neighborFunction','Neighborhood Function',choices=c('gauss','cone','mexicanhat'), selected='gauss')
  	       ),
           fluidRow(
  	         selectInput('initFunction','Initialization',choices=c('uniform','norm'), selected='uniform')
  	       ),
  		     h4("Radius", align="center"),
  		     fluidRow(
  		      column(6,numericInput('startRadius','Start',value=46)),
  		      column(6,numericInput('finalRadius','End',value=1))
  		      ),
  		     h4("Learningrate", align="center"),
  		     fluidRow(
  		       column(6,numericInput('startLearning','Start',value=0.5)),
  		       column(6,numericInput('finalLearning','End',value=0.1))
  		     ),
  		     h4("Cooling Strategy", align="center"),
  		     fluidRow(
  		       column(6,selectInput('coolingRadius','Radius',choices=c('linear','Lead In Lead Out'))),
  		       column(6,selectInput('coolingLearning','Learningrate',choices=c('linear','Lead In Lead Out')))
  		     )
		     ), # conditionalPanel




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
                                 selected = 1))
        ),

          fluidRow(column(12,checkboxInput("TransparentContours", "transparent contours", value=F))),
          fluidRow(column(12,checkboxInput("showBMU", "show bestmatches", value=T))),

actionButton("trainEsom", "Train"),


actionButton("quit", QuitButtonText)


      )
	  ))
	 )



  outputApp=runApp(list(ui = UmatrixUi,
      server = function(input, output, session){
      # initialize container for reactive values
      val = reactiveValues()
      val$Umatrix = NULL
      val$Cls = Cls
      val$Width = 80  # anzahl sichtbare spalten
      val$Height = 50 # anzahl sichtbare zeilen
      val$Stretchfactor = 10 # faktor, so dass val$Width * val$Stretchfactor = 800

      # observe({
      #   UmatrixSize <<- as.numeric(input$UmxSize)
      #   #UmatrixSize = as.numeric(isolate(input$UmxSize))
      # })

      # initialize empty variables
      Epochs <- 0
      Lines <- 0
      Columns <- 0
      Weights <- 0

      progress <- Progress$new()




      #####################
      # Training der Esom / Berechnung der Umatrix (Button)
      #####################
    observeEvent(input$trainEsom, {
      # lade konfiguration f?r die Esom
      epochs <- isolate(input$trainingEpochs)
      k <- isolate(input$numRows)
      m <- isolate(input$numColumns)
      lradius0 <- isolate(input$startRadius)
      lradius1 <- isolate(input$finalRadius)
      Toroid <<- (isolate(input$mapGrid) == "Toroid")

      # FL: maxRadius existiert nichtmehr, da jetzt Ueberlappungen auf dem Torus erlaubt sind
      #maxRadius <- floor(sqrt(k^2 + m^2)/2)
      #if(lradius0 >= maxRadius){
      #  print(paste0("Radius too big. Please use max. Radius of ",(maxRadius-1)))
      #  return(NULL)
      #}

      # verstecke Trainingpanel
      toggle("trainingPanel")

      #if (isolate(input$coolingLearning)=='linear'){
      #  lcoolfun <- 'linear'
      #} else {
      #  lcoolfun <- 'exp'
      #}

      #if (isolate(input$coolingRadius)=='linear'){
      #  rcoolfun <- 'linear'
      #} else {
      #  rcoolfun <- 'exp'
      #}
      lrate0 <- isolate(input$startLearning)
      lrate1 <- isolate(input$finalLearning)

      neighbourhoodFunction <- isolate(input$neighborFunction)
      res <- isolate(input$initFunction)
      if(res == "uniform") initFunction = "uni_min_max"
	  else if(res == "norm") initFunction = "norm_mean_2std"
  #    else if(res == "ica") initFunction = "uni_min_max+ica"

      # eigentliches training der Esom
      print("Initializing ESOM algorithm")
      #MT: Fuers Ausschreiben
      Epochs<<-epochs
      Lines <<- k
      Columns <<- m

      trainSOM <- esomTrain(Data,k,m,Epochs=epochs, StartLearningRate=lrate0, EndLearningRate=lrate1,StartRadius=lradius0,
                       EndRadius=lradius1, NeighbourhoodCooling=input$coolingRadius, LearningRateCooling=input$coolingLearning,
                       NeighbourhoodFunction=neighbourhoodFunction, Toroid=Toroid, ShinyProgress = progress, InitMethod=initFunction, Key = Key, UmatrixForEsom=F)
      EsomDetails<<-list(StartLearningRate=lrate0, EndLearningRate=lrate1,StartRadius=lradius0, EndRadius=lradius1,NeighbourhoodCooling="linear", LearningRateCooling="linear", NeighbourhoodFunction=neighbourhoodFunction)

      if(!fixedBestMatches){
        BestMatches <<- trainSOM$BestMatches
      }
      Weights <<- trainSOM$Weights

      # berechne die Umatrix
      print("Begin calculation of Umatrix")
      val$Umatrix <<- umatrixForEsom(Weights, Lines, Columns,Toroid = Toroid)

      Umatrix <<- val$Umatrix
      BestMatches <<- BestMatches

      # show trainingpanel
      toggle("trainingPanel")
      print("Calculation of Umatrix finished")

   })


    ###############################
    # Plotte Umatrix
    ###############################
    observe({
      output$plotMatrix <- renderPlot({
        if(!is.null(val$Umatrix)){
          if(input$showBMU)
            u <- plotMatrix(val$Umatrix, BestMatches, val$Cls, Toroid = Toroid, TransparentContours
                             = input$TransparentContours, BmSize = as.numeric(input$BmSize), Clean=T)
          else
            u <- plotMatrix(val$Umatrix, Cls = val$Cls, Toroid = Toroid, TransparentContours
                             = input$TransparentContours, BmSize = as.numeric(input$BmSize), Clean=T)

          if(!is.null(DataName)) return(u + ggtitle(DataName))
          else return(u + ggtitle('U-Matrix'))
        }
      },
      height = as.numeric(input$UmxSize) * val$Height * val$Stretchfactor,
      width = as.numeric(input$UmxSize) * val$Width * val$Stretchfactor)
    })

    # update size
    observe({
      if(!is.null(val$Umatrix)){
        if(input$mapGrid == "Toroid"){
          val$Width = ncol(val$Umatrix) * 2
          val$Height = nrow(val$Umatrix) * 2
        }
        else{
          val$Width = ncol(val$Umatrix)
          val$Height = nrow(val$Umatrix)
        }
      }

      # zooming
      val$Stretchfactor = 800/val$Width
    })

    #########
    ############

    session$onSessionEnded(function() {
      stopApp(list(Umatrix=Umatrix, BestMatches=BestMatches, Lines = Lines,
                   Columns = Columns, Epochs = Epochs, Weights=Weights, Toroid = Toroid,
                   EsomDetails = EsomDetails))
    })

    observeEvent(input$quit, {
      stopApp(list(Umatrix=Umatrix, BestMatches=BestMatches, Lines = Lines,
                   Columns = Columns, Epochs = Epochs, Weights=Weights, Toroid = Toroid,
                   EsomDetails = EsomDetails))
    })
  }))

  invisible(outputApp)
}

