iClassification <- function(Umatrix=NULL, BestMatches=NULL, Cls=NULL, Imx = NULL, Toroid = T){
# V = iClassification(Umatrix, BestMatches, Cls)
# Cls = V$Cls
# Key = V$ClsKey
#
  # INPUT
  # Umatrix     The EsomNeurons in Array Format
  # BestMatches
  # Cls
  # Imx             Island Mask
  # OUTPUT
  # Cls
  # Author: FL

  DataName = NULL
  UmatrixSize = 1

  QuitButtonText = "Quit"

  ########
  ###########

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

  if(is.null(BestMatches))
    stop('No BestMatches available')

  if(!is.matrix(BestMatches)){
    Data=as.matrix(BestMatches)
    warning('BestMatches should be a matrix, Umatrix() called as.matrix(). If an follwing error occurs, please change BestMatches to matrix format manually.')
  }
  if(is.null(Cls)) Cls=rep(1,nrow(BestMatches))
  if(ncol(BestMatches) == 2) BestMatches = cbind(1:nrow(BestMatches),BestMatches)

  islandLeft = 1
  islandRight = ncol(Umatrix)*2
  islandTop = 1
  islandBottom = nrow(Umatrix)*2

  if(!is.null(Imx)){
    # start und endpunkte der IMX bestimmen
    lineHeights = rowSums(Imx,na.rm=T)
    islandTop = min(which(lineHeights<ncol(Imx)),na.rm=T)
    islandBottom = length(lineHeights) - min(which(rev(lineHeights)<ncol(Imx)),na.rm=T) + 1

    colHeights = colSums(Imx,na.rm=T)
    islandLeft = min(which(colHeights<nrow(Imx)),na.rm=T)
    islandRight = length(colHeights) - min(which(rev(colHeights)<nrow(Imx)),na.rm=T) + 1
  }

  updateIslandProportions <- function(Imx){
    lineHeights = rowSums(Imx,na.rm=T)
    islandTop <<- min(which(lineHeights<ncol(Imx)),na.rm=T)
    islandBottom <<- length(lineHeights) - min(which(rev(lineHeights)<ncol(Imx)),na.rm=T) + 1

    colHeights = colSums(Imx,na.rm=T)
    islandLeft <<- min(which(colHeights<nrow(Imx)),na.rm=T)
    islandRight <<- length(colHeights) - min(which(rev(colHeights)<nrow(Imx)),na.rm=T) + 1
  }

  if(is.null(Imx)){
    Imx = matrix(0,nrow=nrow(Umatrix)*2, ncol=ncol(Umatrix)*2)
  }

  Width = islandRight - islandLeft + 1
  Height = islandBottom - islandTop + 1

  # zooming
  # if(Width>Height) zoomFactor = 800/Width
  # else zoomFactor = 800/Height
  # Width=Width*zoomFactor
  # Height=Height*zoomFactor
  #
  # Width = round(Width)
  # Height = round(Height)

  mapGridText = 'Planar'
  if(Toroid) mapGridText = 'Toroid'


  ##########
  # Shiny Fenster
  ##########
  UmatrixUi = fluidPage(
    useShinyjs(),
    sidebarLayout( position="right",
      mainPanel(
        #div(plotOutput("UmatrixPlot", Width=as.character(Width), Height=as.character(Height), click = "clickOnUmatrix"),
        #    style = "margin-left:-70px")
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



        uiOutput("buttons"),
        actionButton("addCluster", "Add Cluster"),

        selectInput('mapGrid','Grid Topology',choices=c('Toroid','Planar'),selected=mapGridText),
        uiOutput("clusterOverview"),
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



        fluidRow(column(12,checkboxInput("TransparentContours", "transparent contours", value=F))),
        fluidRow(column(12,checkboxInput("showBMU", "show bestmatches", value=T))),

        fluidRow(
          
          actionButton("quit", QuitButtonText))

      )
    )))

  BestMatches = CheckBestMatches(BestMatches, Cls, shiny=T)
  CheckUmatrix(Umatrix, shiny=T)
  CheckImx(Imx, shiny=T)

  outputApp=runApp(list(ui = UmatrixUi, server = function(input, output, session){

    updateSelectInput(session, "UmxSize", selected = UmatrixSize)
    observe({UmatrixSize <<- as.numeric(input$UmxSize)})

    val = reactiveValues()
    val$UmatrixImg = NULL
    val$Cls = Cls
    val$BestMatches = BestMatches
    val$Imx = Imx
    val$currentLines = c()
    val$existingClusters = 1
    # val$imgWidth = Width
    # val$imgHeight = Height
    val$Width = Width
    val$Height = Height
    val$Stretchfactor = 1

    val$Toroid = Toroid
    observe({val$Toroid = (input$mapGrid == "Toroid")})

 

    ##########
    # button observer
    ##########
    classes <- c(unique(Cls), setdiff(1:100,unique(Cls)))
    lapply(classes, function(i){
      observeEvent(input[[paste0('deleteCluster', i)]],{
        val$Cls[val$Cls == i] <<- 0
        Cls <<- val$Cls
      })
    })

    #########
    # create all buttons
    #########
    output$buttons = renderUI({
      cls <- val$Cls[val$Cls != 0]
      classes <- sort(na.last=T,c(unique(cls), setdiff(1:100,unique(cls))))
      if(any(val$Cls == 0))
        classes = c(classes,0)
      #print(classes)
	  
      colors = DefaultColorSequence()
      colors[colors == "green"] = "lightgreen"


      lapply(classes, function(x){
        if(x == 0)
          fluidRow(id=paste0("row",x), column(8, span(paste0("No Cluster (",sum(val$Cls == x),")"),style=paste0("color:","darkgreen",";font-size:1.3em"))))
        else if(x %in% val$Cls)
          fluidRow(id=paste0("row",x), column(8, span(paste0('Cluster ', x, "(",sum(val$Cls==x),")"),style=paste0("color:",colors[which(classes==x)],";font-size:1.3em"))),
                 column(4, actionButton(paste0("deleteCluster",x),"delete")),
                 style="position:relative;top:50%")
        else
          hidden(fluidRow(id=paste0("row",x), column(5, span(paste0('Cluster ', x),style=paste0("color:",colors[which(classes==x)],";font-size:1.3em"))),
                          column(5, actionButton(paste0("deleteCluster",x),"delete")),
                          style="position:relative;top:50%"))
      })
    })



    ##########
    # rerender Umatrix if classification has changed
    ##########
    observe({
      # update proportions
      updateIslandProportions(val$Imx)

      val$Width = (islandRight - islandLeft + 1)
      val$Height = (islandBottom - islandTop + 1)

      val$Stretchfactor = 800/val$Width

      CurrentDirectory <- getwd()
      setwd(tempdir())
      png("tmpUmatrix.png", width=val$Stretchfactor*val$Width,
          height = val$Stretchfactor*val$Height )
sample(1:212, 50)
      if(input$showBMU)
        print(plotMatrix(Umatrix, BmSize = as.numeric(input$BmSize)*2, BestMatches = val$BestMatches, Cls = val$Cls, Clean = T,
                          Imx=val$Imx, RemoveOcean = T, TransparentContours = input$TransparentContours,
                          Toroid=val$Toroid))
      else
        print(plotMatrix(Umatrix, BmSize = as.numeric(input$BmSize)*2, Cls = val$Cls, Clean = T,
                          Imx=val$Imx, RemoveOcean = T, TransparentContours = input$TransparentContours,
                          Toroid=val$Toroid))

      dev.off()
      val$UmatrixImg <- readPNG("tmpUmatrix.png")
      setwd(CurrentDirectory)
    })

    ##########
    # plot Umatrix
    ##########
    output$UmatrixPlot <- renderPlot({

      # plot.new()
      #plot.window(c(0,1),c(0,1))
      par(mar=c(0.01,0.01,0.01,0.01))
      par(mai=c(0,0,0,0))
      par(xaxs = 'i',yaxs='i')

      rasterImage(val$UmatrixImg,0,0,1,1)

      lines(val$currentLines[,1], val$currentLines[,2],col="red")
    }, width = function(){val$Width*val$Stretchfactor*as.numeric(input$UmxSize)},
    height = function(){val$Height*val$Stretchfactor*as.numeric(input$UmxSize)})#, Height=val$imgHeight)

    ##########
    # react on click
    ##########
    observeEvent(input$clickOnUmatrix,{
      x <- isolate(input$clickOnUmatrix$x)
      y <- isolate(input$clickOnUmatrix$y)

      #print(paste(x ,y))

      val$currentLines = rbind(val$currentLines,c(x,y))
    })


    ################
    # reaction for button "Add Cluster"
    ################
    observeEvent(input$addCluster,{
        # check if there are already enough points to draw a complete polygon
        if(is.null(val$currentLines)) return()
        if(nrow(val$currentLines) < 3) return()

        print("Add new Cluster")

        # define unused new class
        newClassId = sort(na.last=T,setdiff(1:100,val$Cls))[1] # first unused number

        if(val$Toroid){
          Height = (islandBottom-islandTop+1)
          Width = (islandRight - islandLeft +1)

          Polygon = cbind( val$currentLines[,1]*Width + (islandLeft)-0.5, (1-val$currentLines[,2])*Height + (islandTop) - 0.5)

          # BestMatches vervierfachen
          res = ToroidUmatrix(Umatrix, BestMatches)
          BestMatches4X = res$BestMatches

          # filter out everything thats not on the Imx
          t <- c()
          for(i in 1:nrow(BestMatches4X))
            t[i] = (Imx[BestMatches4X[i,2], BestMatches4X[i,3]] == 1)
          #BestMatches4X = BestMatches4X[t,]

          # find all BestMatches within polygon
          FilteredBm <- in.poly(BestMatches4X[,c(3,2)]*1000, Polygon*1000)

          # filtered Key of BestMatches
          FilteredKey = BestMatches4X[FilteredBm,1,drop=F]
        }
        else{ # kein toroid
          Polygon = cbind( val$currentLines[,1]*ncol(Umatrix)+0.5, (1-val$currentLines[,2])*nrow(Umatrix) + 0.5)
          FilteredBm <- in.poly(BestMatches[,c(3,2)]*1000, Polygon*1000)
          FilteredKey = BestMatches[FilteredBm,1,drop=F]
        }
        if(sum(FilteredBm)>=1)
          val$Cls[BestMatches[,1] %in% FilteredKey] <<- newClassId

        val$currentLines <<- c()
        Cls <<- val$Cls

        print("New class added")
    })

    ###########
    ###########
    
    session$onSessionEnded(function() {
      print("program closed")
      Cls = list(Cls = Cls, ClsKey=BestMatches[,1], UniqueClasses = unique(Cls))
      stopApp(Cls)
    })

    observeEvent(input$quit, {
      Cls = list(Cls = Cls, ClsKey=BestMatches[,1], UniqueClasses = unique(Cls))
      stopApp(Cls)
    })

  }))

  return(outputApp)
}
