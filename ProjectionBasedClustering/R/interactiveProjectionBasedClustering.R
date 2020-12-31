IPBC=interactiveProjectionBasedClustering = function(Data, Cls=NULL) {

    ## Schritt 0: Check and sanetize input ####
    if (missing(Data))
      stop("No data was given as input")
    
    if (!inherits(Data,"matrix")){
      warning("Input data must be a matrix, assuming data frame and calling as.matrix()")
      Data=as.matrix(Data)
    }
    
    
    Key=rownames(Data)
    if(is.null(Key)){
      warning('"Data" has no rownames. Setting Key to 1:NoOfCases and Naming Output "Cls" accordingly.')
      Key=1:nrow(Data)
      
    }
    # Parameter initialisieren
    LastProjectionMethodUsed=NULL
    BestMatchesNotExtended=NULL
    UmatrixNotExtended=NULL
   
    Umatrix = NULL 
    bestmatches = NULL
    Umatrix4Toroid =NULL
    bestmatches4Toroid =NULL
    bmus = NULL 
    mergei = 1:nrow(Data) #give back all cases in one group, if function exists before clustering and no cls given
    ubmus = NULL 
    udim = NULL 
    uniq = NULL 
    outplot = NULL 
    oubmus = NULL 
    extendBorders=10
    LC=NULL
    imx<-NA
    
    if(is.null(Cls))#give back all cases in one group, if function exists before clustering and no cls given
      Cls=rep(1,nrow(Data))
     
    # Original CLS speichern.
    ClsO=Cls 
    ## Helper functions ####
    
    # Normalizes the U-Matrix. As used in plotmatrix
    matrixnormalization <- function(Matrix){
      ## MT: Normalization der Umatrix werte
      # Milligan, Copper 1988
      # A Study of Standadization of Variables in Cluster Analysis,
      # robust Normalization Z_5 :"Z_5 is bounded by 0.0 and 1.0
      # with at least one observed value at each of these end points"
      quants <- quantile(as.vector(Matrix), c(0.01, 0.5, 0.99), na.rm = T)
      
      minu <- quants[1]
      maxu <- quants[3]
      # Verhaeltnis zwischen minhoehe/maxHoehe = 1/HeightScale
      
      Matrix <- (Matrix - minu) / (maxu - minu)
      
      #MT: Die Level muessen vor der Begrenzung der Werte
      # auf 0 und 1 gesetz werden, aber nachdem der Wertebereich umgeschoben
      # wurde, damit die Dichte in den Grenzbereichen abgeschetzt werden kann
      indmax <- which(Matrix > 1, arr.ind = T)
      indmin <- which(Matrix < 0, arr.ind = T)
      
      if (length(indmax) > 0)
        Matrix[indmax] <- 1
      if (length(indmin) > 0)
        Matrix[indmin] <- 0
      
      return(Matrix)
    }

    colormap=GeneralizedUmatrix::UmatrixColormap
    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      #showticklabels = FALSE,
      showgrid = FALSE
    )
    ay <- list(
      autorange = "reversed",
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      #showticklabels = FALSE,
      showgrid = FALSE
    )
    
    createParams = function (Umatrix, bestmatches, Cls){
      
      udim <<- dim(Umatrix)
      bmus  <<- cbind(bestmatches[, 1], bestmatches[, 2])
      uniq  <<- getUniquePoints(bmus)
      ubmus  <<- uniq$unique
      mergei <<- uniq$mergeind
      
      if (missing(Cls) || length(Cls) != nrow(bestmatches)){
        
        Cls <<- rep(1, nrow(ubmus))
      } else {
        Cls <<- Cls[uniq$sortind]
      }
      oubmus <<- cbind(uniq$sortind, uniq$unique[, 1], uniq$unique[, 2])
      
    }
    
    ## Schritt 2: shiny shiny interactive tool ####
    ui <- shiny::shinyUI(
      
      tagList(
        ## Script for resizing plotly plots
        tags$head(tags$script('
                              var dimension = [0, 0];
                              $(document).on("shiny:connected", function(e) {
                              var plot = document.getElementsByClassName("col-sm-9")[0];
                              dimension[0] = plot.clientWidth*2;
                              dimension[1] = plot.innerHeight*2;
                              Shiny.onInputChange("dimension", dimension);
                              });
                              $(window).resize(function(e) {
                              var plot = document.getElementsByClassName("col-sm-9")[0];
                              dimension[0] = plot.clientWidth;
                              dimension[1] = plot.innerHeight;
                              Shiny.onInputChange("dimension", dimension);
                              });
                              ')),

        navbarPage(
          title = "Interactive Projection-Based Clustering",
          theme = shinythemes::shinytheme("flatly"),
            tabPanel(title = "2D                ",
            sidebarPanel(
              width = 3,
              height="100%",
              tabsetPanel(
                type = "tabs",
              tabPanel(
                title = "Projection",
                h3(" "),
                selectInput('projections','Choose Projection Method',choices=c('PCA','CCA','ICA','MDS','NeRV','ProjectionPursuit',"Pswarm",'SammonsMapping','UniformManifoldApproximationProjection','tSNE'),selected='NeRV'),
                
                # Further Parameter-Querys for Projections
                #PCA
                
                conditionalPanel(condition ="input.projections=='PCA'", 
                                 checkboxInput("PCAscale", "Scale", value = FALSE, width = NULL)),
                conditionalPanel(condition ="input.projections=='PCA'", 
                                 checkboxInput("PCACenter", "Center", value = FALSE, width = NULL)),
                
                #CCA
                
                conditionalPanel(condition ="input.projections=='CCA'", selectInput(
                  "CCAMethod", "Method",
                  c('euclidean','cosine','chebychev','jaccard','manhattan','binary'))),
                conditionalPanel(condition ="input.projections=='CCA'",
                                 numericInput("CCAEpochs", "Epochs", value=20, min = 1, max = NA, step = 1)),
                conditionalPanel(condition ="input.projections=='CCA'",
                                 numericInput("CCASteps", "initial Stepsize", value=0.5, min = 0.001, max = NA, step = 0.0001)),
                
                
                
                #ICA
                
                conditionalPanel(condition ="input.projections=='ICA'",
                                 numericInput("ICAIterations", "Iterations", value=100, min = 1, max = NA, step = 1)),
                conditionalPanel(condition ="input.projections=='ICA'",
                                 numericInput("ICASteps", "Alpha", value=1, min = 1, max = 2, step = 0.0001)),
                
                #MDS
                
                conditionalPanel(condition ="input.projections=='MDS'", selectInput(
                  "MDSMethod", "Distance Method", selected="euclidan",
                  c('euclidean','maximum','canberra','manhattan'))),
                
                
                
                #NERV
                
                conditionalPanel(condition ="input.projections=='NeRV'",
                                 numericInput("NERVIterations", "Iterations", value=20, min = 1, max = NA, step = 1)),
                conditionalPanel(condition ="input.projections=='NeRV'",
                                 numericInput("NERVLambda", "Lambda", value=0.1, min = 0.0001, max = NA, step = 0.0001)),
                conditionalPanel(condition ="input.projections=='NeRV'",
                                 numericInput("NERVNeighbors", "Neighbors", value=20, min = 1, max = NA, step = 1)),
                
                #ProjectionPursuit
                
                conditionalPanel(condition ="input.projections=='ProjectionPursuit'", selectInput(
                  "PPMethod", "Indexfunction", c('logcosh','exp'))),
                
                conditionalPanel(condition ="input.projections=='ProjectionPursuit'",
                                 numericInput("PPIterations", "Iterations", value=200, min = 1, max = NA, step = 1)),
                conditionalPanel(condition ="input.projections=='ProjectionPursuit'",
                                 numericInput("PPAlpha", "Alpha", value=1, min = 0.0001, max = NA, step = 0.0001)),
                
                #Polar Swarm
                
                conditionalPanel(condition ="input.projections=='Pswarm'", selectInput(
                  "PswarmDistMethod", "Distance Method", selected="euclidan",
                  c('euclidean','maximum','canberra','manhattan','chord','hellinger','geodesic','kullback','mahalanobis'))),
                
                #SammonsMapping
                
                conditionalPanel(condition ="input.projections=='SammonsMapping'", selectInput(
                  "SMMethod", "Method",
                  c('euclidean','maximum','canberra','manhattan'))),
                
                #TSNE
                
                conditionalPanel(condition ="input.projections=='tSNE'", selectInput(
                  "tSNEMethod", "Method"  , selected = 'logcosh', c('euclidean','maximum','canberra','manhattan'))),
                
                conditionalPanel(condition ="input.projections=='tSNE'",
                                 numericInput("tSNEIterations", "Iterations", value=1000, min = 1, max = NA, step = 1)),
                conditionalPanel(condition ="input.projections=='tSNE'", 
                                 checkboxInput("tSNEWhite", "Whitening", value = FALSE, width = NULL)),
                
                # UniformManifoldApproximationProjection
                
                 conditionalPanel(condition ="input.projections=='UniformManifoldApproximationProjection'",
                                  numericInput("knn", "k nearest neighbors", value=15, min = 2, max = NA, step = 1)),
                 conditionalPanel(condition ="input.projections=='UniformManifoldApproximationProjection'",
                                  numericInput("Epochs", "training length, number of epochs", value=200, min = 2, max = NA, step = 1)),
                
                #Other things
                
                actionButton("generate", HTML("Visual Analytics <br/> with g. U-Matrix"), icon = icon("calculator")),
                tags$hr(),
                numericInput("markerSize", "Bestmatch Size", value=7, min = 1, max = 30, step = 0.1),
                numericInput("ExtendBorders", 'Set Extend Borders "x"', value=10, min = 0, max = 100, step = 1),
                checkboxInput(inputId = 'ExtendBorders10',label = 'Extend Borders by "x" Lines & Columns',value = FALSE),
                checkboxInput(inputId = 'Toroid',label = 'Toroidal',value = FALSE)
              ),
  
              tabPanel(
                title = "Clustering",
                h3(" "),
                
                                h5("Execute Interactive",align="center"),
                actionButton(inputId = "Clust", label = "Create Cluster", icon = icon("location-arrow")),
               
                h4("--"),
                h5("Modify Existing",align="center"),
                selectInput(
                  inputId = "ClsSelect",
                  label = "Select Class",
                  choices = unique(Cls)
                ),
                actionButton(inputId = "AddToCls", label = "Add Selected", icon = icon("plus-circle")),
                actionButton(inputId = "ClearAll", label = "Clear All", icon = icon("undo-alt")),
             
                h4("--"),
                h5("Execute Automatic",align="center"),
                numericInput(inputId = "NumberClusters", 
                             label ="No. Clusters", value = 1),
                checkboxInput(inputId = 'StructureType',label = 'Change Structure Type',value = TRUE),
                actionButton(inputId = "Clustering", label = "Clustering", icon = icon("calculator")),

                tags$hr()
                
              ),
              tabPanel(
                title = "Info",
                h3(" "),
                htmlOutput("umxinfo"),
                tags$h5("Corresponding Data"),
                actionButton("pointinfo", "View selected"),
                tags$hr()
              )
            ),
            actionButton(
              inputId = "Exit",
              label = "Exit",
              class = "btn-primary", icon = icon("window-close")
            )
            )#end sidebar pannel
            ),
              mainPanel(
              
              plotly::plotlyOutput("Plot",
                                   # height =  
                                   #   udim[1] * 5 * (quadtile + 1),
                                   # width = 
                                   #   udim[2] * 5 * (quadtile + 1)),
                                   width = "auto"),
          
              width = 9,
              height = "100%"
            )
            

        )#end navbarPage
      )# end tagList
      )# end shiny::shinyUI
    
    server <- shiny::shinyServer(function(input, output, session) {
      
      ## Helper functions:
      selectedids <- function() {
        
        # Beware of off by one errors here
        # Keep in mind: shiny and plotly use javascript ennumeration
        selected <- plotly::event_data("plotly_selected")
        if (is.null(selected) || length(selected) == 0) {
          return(NA)
        }
        classes  <- unique(Cls)[selected$curveNumber]
        indicies <- selected$pointNumber
        d        <- cbind(classes, indicies)
        # give toroid copys the id of the original
        # and convert to R indicies
        d[, 2]   <- (d[, 2] %% as.vector(table(Cls))[d[, 1]]) + 1
        # now transfer from inner-cluster indicies to indicies over Cls.
        Clsids   <- apply(
          d,
          MARGIN = 1,
          FUN = function(row) {
            return(which(Cls == row[1])[row[2]])
          }
        )
        return(Clsids)
      }#end selectedids
      
      TopographicMapTopView_hlp=function(Cls,Tiled){
        #@TIM: die Parameteruebergabe Umatrix,ubmus muss du besser loesen. aus dem globalen workspace das zuu nehmen ist sehr schlecht
        # das interakCtive anpassen an fenster groesse erfordert nun ein button click, weis nicht wieso
        # pruef auch das mal bitte 
     
        plt=GeneralizedUmatrix::TopviewTopographicMap(GeneralizedUmatrix = Umatrix,BestMatchingUnits = ubmus,
                   Cls=Cls,Tiled=Tiled,BmSize = input$markerSize,ShinyBinding=TRUE,ShinyDimension=input$dimension[1],Session=session)
        
        shiny::updateSelectInput(session,
                                 "ClsSelect",
                                 label = "Select Class",
                                 choices = unique(Cls))
        
        requireNamespace("plotly")
        Rendered <- plotly::renderPlotly({
          plt
        })
        output$Plot=Rendered
        outplot<<-plt
       
      }#end TopographicMapTopView_hlp
     
      selectedpoints <- function() {
        
        d <- plotly::event_data("plotly_selected")
        if (length(d) == 0)
          "Nothing selected!"
        else{
          classes  <- unique(Cls)[d$curveNumber]
          indicies <- d$pointNumber
          d        <- cbind(classes, indicies)
          # give toroid copys the id of the original
          # and convert to R indicies
          d[, 2]    <- (d[, 2] %% as.vector(table(Cls))[d[, 1]]) + 1
          # now transfer from inner-cluster indicies to indicies over Cls.
          Clsids   <- apply(
            d,
            MARGIN = 1,
            FUN = function(row) {
              return(which(Cls == row[1])[row[2]])
            }
          )
          # to get original ids:
          unipointids <- uniq$sortind[Clsids]
          # and now the Points which sit on the same gridpostions:
          # if(length(unipointids) > 1) # TODO
          opoints <- bestmatches[colSums((diag(nrow(
            uniq$IsDuplicate
          )) + uniq$IsDuplicate)[unipointids,]) > 0,]
          if (all(is.na(Data)))
            return(opoints) # Beaware: opoints often more than unipointsid
          Data[opoints[, 1],]
        }
      } #end helper fun selectedpoints
      
      # HTML Output: Umatrix Properties
     #  observeEvent(LastProjectionMethodUsed,{ 
     # if(!is.null(LastProjectionMethodUsed))
      #ToDo: Currently just works on the first click of info but does not actualize its data.
      output$umxinfo <- renderUI({
        HTML(paste0("
                    <table>
                    <tr>
                    <th>Lines: </th>
                    <td>",  udim[1], "</td>
                    </tr>
                    <tr>
                    <th>Columns: </th>
                    <td>",  udim[2], "</td>
                    </tr>
                    <tr>
                    <th>Bestmatches: </th>
                    <td>", nrow(bestmatches), "</td>
                    </tr>
                    <tr>
                    <th>Last DR method used: </th>
                    <td>", LastProjectionMethodUsed, "</td>
                    </tr>
                    <tr>
                    <th>No. Clusters: </th>
                    <td>", length(unique(Cls)), "</td>
                    </tr>                    
                    <tr>
                    <th>Last Structure type used: </th>
                    <td>",  input$StructureType, "</td>
                    </tr>  
                    </table>
                    "))
      })
     # })#end observe info
      observeEvent(input$dimension,{ 
      })
      
      # Button: View Selected
      observeEvent(input$pointinfo,{
        seldat <- selectedpoints()
        output$dataout <- renderDataTable({seldat})
        shiny::showModal(shiny::modalDialog(title ="Data corresponding to selected points", shiny::dataTableOutput("dataout"), easyClose = TRUE )
      )})

      
      # Button: Merge Clusters
      observeEvent(input$Merge, {
        selected <- plotly::event_data("plotly_selected")
        if (is.null(selected) || length(selected) == 0 ){
          return()
        }
        points <- (selected[, 2] %% nrow(ubmus)) + 1
        Clss <- unique(Cls[points])
        ids <- which(Cls %in% Clss)
        Cls[ids] <- min(Clss)
        Cls <<- getUniquePoints(Cls)$mergeind
        TopographicMapTopView_hlp(Cls=Cls,Tiled=input$Toroid)
      })
      
      ##Simple Extend ----
      observeEvent(input$ExtendBorders10, {
        if(!is.null(Umatrix)){ #wird direkt bei start augeprueft
          if(input$ExtendBorders10 & isFALSE(input$Toroid)){
            V=GeneralizedUmatrix::ExtendToroidalUmatrix(Umatrix4Toroid,bestmatches4Toroid,extendBorders)
            Umatrix<<-V$Umatrix
            bestmatches<<-V$Bestmatches
            ClsTemp=Cls#create params overwrites cls as it is not the length of bestmachtes
            createParams(Umatrix, bestmatches,Cls = Cls)
            Cls<<-ClsTemp #but in this cases cls is only the length of unique bestmatches, ToDo: kill createParams function anywhere
            TopographicMapTopView_hlp(Cls,input$Toroid) #shortcut: Umatrix muss existieren
          }else{
            bestmatches<<-bestmatches4Toroid
            Umatrix<<-Umatrix4Toroid
            ClsTemp=Cls
            createParams(Umatrix4Toroid, bestmatches4Toroid,Cls = Cls)
            Cls<<-ClsTemp
            TopographicMapTopView_hlp(Cls,input$Toroid) #shortcut: Umatrix muss existieren
          }
        }
      })
      
      ##Toroidal ----
      observeEvent(input$Toroid, {
        if(!is.null(Umatrix)){ #wir direkt bei start augeprueft
          bestmatches<<-bestmatches4Toroid
          Umatrix<<-Umatrix4Toroid
          # print(length(Cls))
          # createParams(Umatrix4Toroid, bestmatches4Toroid,Cls = Cls)
          # print(length(Cls))
          # Cls <<- getUniquePoints(Cls)$mergeind
          TopographicMapTopView_hlp(Cls,input$Toroid) #shortcut: Umatrix muss existieren
        }
      })
      
      # Button: Create new cluster
      observeEvent(input$Clust, {
        Clsids <- selectedids()
        if (all(!is.na(Clsids))) { 
          id <- max(Cls) + 1
          points <- Clsids
          Cls[points] <- id
          Cls <<- getUniquePoints(Cls)$mergeind
          TopographicMapTopView_hlp(Cls=Cls,Tiled=input$Toroid)
        }
      })
      
      # Button: Add to cluster
      observeEvent(input$AddToCls, {
        
        id <- input$ClsSelect
        Clsids <- selectedids()
        
        if (all(!is.na(Clsids))){ 
          Cls[Clsids] <- id
          Cls <- getUniquePoints(as.numeric(Cls))$mergeind
          Cls <<- Cls
          TopographicMapTopView_hlp(Cls=Cls,Tiled=input$Toroid)
        }
      })
      
      observeEvent(input$ClearAll, {
        ind=getUniquePoints(data = bestmatches)$sortind
        Cls <<-rep(1,length(ind))
        TopographicMapTopView_hlp(Cls=Cls,Tiled=input$Toroid)

      })
      
      ##AutomaticClustering ----
      observeEvent(input$Clustering, {
        if(!is.null(Umatrix)){
          #wir direkt bei start augeprueft
          #shortcut: Umatrix muss existieren
          Clsfull <- ProjectionBasedClustering::ProjectionBasedClustering(input$NumberClusters,
                                            Data = Data,
                                            BestMatches = bestmatches,
                                            LC = LC,
                                            StructureType = input$StructureType,
                                            PlotIt = FALSE
                                            )
          ind=getUniquePoints(data = bestmatches)$sortind
          Cls <<-Clsfull[ind] #intern unique cls mitfueren
          TopographicMapTopView_hlp(Cls=Cls,Tiled=input$Toroid)
        } 
 
      })
      observeEvent(input$ExtendBorders, {
        if(!is.null(input$ExtendBorders)){
                if(is.finite(input$ExtendBorders)){
                  extendBorders<<-input$ExtendBorders
                }
        }
      })
      ## ExtendBorders----
      #ToDo: extend more than once
      # observeEvent(input$ExtendBorders, {
      #   if(!is.null(Umatrix)){
      #     if(!is.null(input$ExtendBorders)){
      #       if(is.finite(input$ExtendBorders)){
      #     if(input$ExtendBorders>0){
      #     extendBorders<<-input$ExtendBorders
      #     Umatrix<<-cbind(Umatrix[,c((ncol(Umatrix)-extendBorders-1):ncol(Umatrix))],Umatrix,Umatrix[,1:extendBorders])
      #     Umatrix<<-rbind(Umatrix[c((nrow(Umatrix)-extendBorders-1):nrow(Umatrix)),],Umatrix,Umatrix[1:extendBorders,])
      #     Umatrix<<-matrixnormalization(Umatrix)
      #     bmu=bestmatches
      #     bmu[,1]=bmu[,1]+length(c((nrow(Umatrix)-extendBorders-1):nrow(Umatrix)))
      #     bmu[,2]=bmu[,2]+length(c((ncol(Umatrix)-extendBorders-1):ncol(Umatrix)))
      #     bestmatches<<-bmu
      #     createParams(Umatrix, bestmatches,Cls = Cls)
      #     extendBorders<<-0
      #     }
      #     if(input$ExtendBorders<0){
      #       extendBorders<<- -(input$ExtendBorders)
      #       Umatrix<<-Umatrix[-c((nrow(Umatrix)-extendBorders+1):nrow(Umatrix)),]
      #       Umatrix<<-Umatrix[-c(1:(extendBorders+1)),]
      #       
      #       Umatrix<<-Umatrix[,-c((ncol(Umatrix)-extendBorders+1):ncol(Umatrix))]
      #       Umatrix<<-Umatrix[,-c(1:(extendBorders+1))]#
      #       Umatrix<<-matrixnormalization(Umatrix)
      #       bmu=bestmatches
      #       bmu[,1]=bmu[,1]-length(c((nrow(Umatrix)-extendBorders+1):nrow(Umatrix)))
      #       bmu[,2]=bmu[,2]-length(c(ncol(Umatrix)-extendBorders+1):ncol(Umatrix))
      #       bestmatches<<-bmu
      #       createParams(Umatrix, bestmatches,Cls = Cls)
      #       extendBorders<<-0
      #     }
      #       
      #   
      #     TopographicMapTopView_hlp(Cls=Cls,Tiled=input$Toroid)
      #     } 
      #     }
      #   }
      # })
      
      ##Select Projection ----
      # Calculate Projection and GeneralizedUmatrix and plot result.
      observeEvent(input$generate,{
        
        shiny::showModal(shiny::modalDialog("Please wait while the Generalized U-matrix is being calculated",style = "font-size:20px", easyClose = TRUE))
        
        type=input$projections 
   
        project=function(type){
          switch(type,
                 
                 PCA = ProjectionBasedClustering::PCA(Data,Center=input$PCACenter,Scale=input$PCAscale,OutputDimension = 2,Cls=Cls)$ProjectedPoints,
                 CCA = ProjectionBasedClustering::CCA(Data,OutputDimension = 2,Cls=Cls,Epochs =input$CCAEpochs, method = input$CCAMethod, alpha0=input$CCASteps)$ProjectedPoints,
                 ICA = ProjectionBasedClustering::ICA(Data,OutputDimension = 2,Cls=Cls, Iterations =input$ICAIterations,Alpha = input$ICASteps )$ProjectedPoints,
                 MDS = ProjectionBasedClustering::MDS(Data,OutputDimension = 2,Cls=Cls,  method=input$MDSMethod)$ProjectedPoints,
                 NeRV = ProjectionBasedClustering::NeRV(Data= Data,OutputDimension = 2,Cls=Cls, iterations = input$NERVIterations, lambda = input$NERVLambda, neighbors =input$NERVNeighbors),
                 ProjectionPursuit= ProjectionBasedClustering::ProjectionPursuit(Data,OutputDimension = 2,Cls=Cls, Iterations = input$PPIterations, Indexfunction =input$PPMethod , Alpha =input$PPAlpha )$ProjectedPoints,
                 SammonsMapping = ProjectionBasedClustering::SammonsMapping(Data,OutputDimension = 2,Cls=Cls,  method=input$MDSMethod)$ProjectedPoints,
                 Pswarm= ProjectionBasedClustering::PolarSwarm(Data,Cls=Cls, method = input$PswarmDistMethod)$ProjectedPoints,
                 tSNE = ProjectionBasedClustering::tSNE(Data,OutputDimension = 2,Cls=Cls, method = input$tSNEMethod, Iterations = input$tSNEIterations, Whitening = input$tSNEWhite)$ProjectedPoints,
                 UniformManifoldApproximationProjection = ProjectionBasedClustering::UniformManifoldApproximationProjection(Data,Cls=Cls, k = input$knn, Epochs = input$Epochs)$ProjectedPoints
          )
        }
        LastProjectionMethodUsed<<-type
        pData=project(type)
        
        # construct full CLS from unique bestmatches
        if(length(Cls)!=length(ClsO)){
          reCls =1:(length(uniq$mergeind))
          Cls <<- Cls[uniq$mergeind]
        }
        
        ## Generalized Umatrix ----
       
        newU=GeneralizedUmatrix(Data=Data,ProjectedPoints = pData,Cls=Cls,Toroid = TRUE)
        #only for output/external usage
        UmatrixNotExtended<<-newU$Umatrix
        BestMatchesNotExtended<<- newU$Bestmatches#
        #for internal usage and extentend bordes
        Umatrix <<- newU$Umatrix
        bestmatches <<- newU$Bestmatches
        
        Umatrix4Toroid <<- newU$Umatrix
        bestmatches4Toroid <<- newU$Bestmatches
        
        #To be deleted later, if extend borders works interactively ----
        if(input$ExtendBorders10){
          
          V=GeneralizedUmatrix::ExtendToroidalUmatrix(Umatrix,bestmatches,extendBorders)
          Umatrix<<-V$Umatrix
          bestmatches<<-V$Bestmatches
            #Umatrix<<-cbind(Umatrix[,c((ncol(Umatrix)-extendBorders-1):ncol(Umatrix))],Umatrix,Umatrix[,1:extendBorders])
            #Umatrix<<-rbind(Umatrix[c((nrow(Umatrix)-extendBorders-1):nrow(Umatrix)),],Umatrix,Umatrix[1:extendBorders,])

            #bmu=bestmatches
            #bmu[,1]=bmu[,1]+length(c((nrow(Umatrix)-extendBorders-1):nrow(Umatrix)))
            #bmu[,2]=bmu[,2]+length(c((ncol(Umatrix)-extendBorders-1):ncol(Umatrix)))

        }
        ## until here
        createParams(Umatrix, bestmatches,Cls = Cls)
        LC <<-c(newU$Lines,newU$Columns)
        TopographicMapTopView_hlp(Cls=Cls,Tiled=input$Toroid)
        
        shiny::removeModal()
        
      })
      
      
      # Button: Exit ----
      #MT hier fehlt noch names(Cls)=Key
      observeEvent(input$Exit, {
        Cls=normCls(Cls[mergei])$normalizedCls
        names(Cls)=Key
        stopApp(list(Cls = Cls,Umatrix=UmatrixNotExtended, Bestmatches=BestMatchesNotExtended,LastProjectionMethodUsed=LastProjectionMethodUsed,TopView_TopographicMap = outplot))
      })
    })
    
    output <- runApp(list(ui = ui, server = server))
    return(output)
  }

getUniquePoints <- function(data,mindist = 1e-10){
  # U <- getUniquePoints(data)
  # return only the unique points in data
  #
  # INPUT
  # data[1:n,1:d]   The vector/matrix with n data points of dimension d
  #				          the points are in the  rows
  # mindist					everything with less distance is considered equal
  #
  # OUTPUT
  # a list U containg:
  # UniqueData  <- U$unique[1:u,1:d]       the data points  without duplicate points
  # UniqueInd   <- U$sortind[1:u]		       an index vector such that unique ==  data[sortind,]
  # Uniq2DataInd<- U$mergeind[1:n] 	       an index vector such that   data ==  unique[mergeind,]
  # IsDuplicate <- U$IsDuplicate[1:n,1:n]  for i!=j IsDuplicate[i,j]== 1  iff data[i,] == data[j,] ;  IsDuplicate[i,i]==0
  #
  # Complete Code rewrite by FP 03/17 
  # 1.Editor: MT 03/17: mindist als Parameter eingefuehrt
  
  
  # when data is a vector, convert to matrix
  if (methods::is(data,"numeric") || methods::is(data,"complex")) {
    data <- matrix(data, ncol = 1)
  } else if (is.data.frame(data)) {
    data <- as.matrix(data)
  } else if (!is.matrix(data)) {
    stop("getUniquePoints input is neither a (numeric or complex) vector, matrix or data.frame.")
  }
  
  NumData = nrow(data)
  
  distsmat = as.matrix(dist(data))
  # * 1 is to "hack" TRUE to 1 and FALSE to 0
  IsDuplicate = ((distsmat) < mindist) * 1 - diag(NumData)
  
  # Hack: set distance of every node to itself to NA (we're not interested in these distances)
  # otherwise every element would be a duplicate of itself later.
  distsmat = distsmat + diag(rep(NA, nrow(distsmat)))
  
  # No duplicates found
  if (length(which((distsmat) < mindist)) == 0) {
    return(list(
      "unique" = unique(data),
      "sortind" = c(1:NumData),
      "mergeind" = c(1:NumData),
      IsDuplicate = IsDuplicate
    ))
  }
  
  # save rownames
  origrows = rownames(data)
  # use rownames to save original index
  rownames(data) = c(1:NumData)
  
  #  find out which distances are smaller than mindist (and as such duplicates)
  # (except the diagonal, so we artificially increase it)
  dups = which((distsmat) < mindist, arr.ind = T)
  
  # indicies of rows which will be deleted
  delinds = c()
  while (dim(dups)[1] > 0) {
    ind = dups[1, 1]
    delinds = c(delinds, ind)
    # remove every row in which the duplicate also occurs.
    # this is to prevent to also delete the "original"
    # if not save, encapuslate right hand side in matrix( ... , nrow = 2)
    dups = dups[-(which(dups == ind, arr.ind = T)[, 1]), ]
  }
  # delete duplicates, stay as matrix (important for vector input) and keep rownames
  uniquedata = data[-delinds, ]
  uniquedata = matrix(uniquedata, ncol = ncol(data))
  rownames(uniquedata) = rownames(data)[-delinds]
  
  # the rownames contain the indicies in the original data
  sortind = as.numeric(rownames(uniquedata))
  
  # calculate the mergind. At the end there should not be a NA in mergeind left
  mergeind = rep(NA, NumData)
  mergeind[sortind] = 1:nrow(uniquedata)
  # this works due to the fact that:
  #   union(sortind, delinds) == 1:NumData (if everything is sorted)
  for (i in delinds) {
    candidates = which(IsDuplicate[i, ] == 1)
    original = intersect(candidates, sortind)
    mergeind[i] = which(sortind == original)
  }
  
  # restore rownames of the original data
  rownames(uniquedata) = origrows[-delinds]
  
  return(
    list(
      "unique" = as.matrix(uniquedata),
      "sortind" = sortind,
      "mergeind" = mergeind,
      IsDuplicate = IsDuplicate
    )
  )
  
} # end fun bestmatches

normCls <- function(Cls) {
  #E<-normCls(Cls);
  #NormalizedCls    <- E$normalizedCls      #    Cls consistently recoded to positive consecutive integers
  #NormalizedClasses<- E$normalizedClasses  #    the different class numbers in NormalizedCls
  #UniqueCls        <- E$uniqueClasses      #    the different class numbers in Cls such that 
  #AnzClasses       <- E$numberOfClasses    #    the number of different classes
  # 
  # Values in Cls are consistently recoded to positive consecutive integers
  # INPUT
  # Cls                  vector of class identifiers can be integers or
  #                      NaN's, need not be consecutive nor positive
  # OUTPUT list of 
  # normalizedCls           Cls consistently recoded to positive consecutive integers
  # normalizedClasses        the different class numbers in NormalizedCls
  # uniqueClasses            the different class numbers in Cls such that 
  #                           NormalizedCls(i) <-> UniqueCls(i)
  # numberOfClasses           the number of different classes
  
  # ALU 2014
  # angepasst an Mdbt und Doku standards
  
  
  uniqueClasses <- sort(na.last = T, unique(Cls))
  numberOfClasses <- length(uniqueClasses)
  unique2Cls <- NULL #  initializing the vector
  
  for (i in 1:length(Cls)) {
    # calculating the indexes of elements of Cls in uniqueClasses
    unique2Cls <- c(unique2Cls, which(uniqueClasses == Cls[i]))
  }
  
  if (numberOfClasses > 0) {
    normalizedClasses <- c(1:numberOfClasses)
    normalizedCls <- normalizedClasses[unique2Cls]
  }
  else {
    normalizedClasses <- Cls
  }
  
  return(
    list(
      normalizedCls = normalizedCls,
      normalizedClasses = normalizedClasses,
      uniqueClasses = uniqueClasses,
      numberOfClasses = numberOfClasses
    )
  )
}
