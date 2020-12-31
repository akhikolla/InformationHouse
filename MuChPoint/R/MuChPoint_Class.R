##' Class "MuChPoint"
##'
##' Class of object returned by the \code{MuChPoint} function.
##'
##' @slot S a vector object of type \code{numeric},
##' giving the values of the statistics S_n(n_1,...,n_L) following the number L.
##'
##' @slot N a numeric vector with the position of the different break-points.
##'
##' @slot bt an inferior triangular matrix containing the positions of break-points following the
##' number of break-points (in rows).
##'
##' @param object an object with class \code{MuChPoint}
##'
##' @seealso See also \code{\link{plot,MuChPoint-method}} and \code{\link{MuChPoint}}.
##'
##' @references Article: BRAULT V., OUADAH S., SANSONNET L. and LEVY-LEDUC C. Nonparametric
##' homogeneity tests and multiple change-point estimation for analyzing large Hi-C data matrices.
##' Journal of Multivariate Analysis, 2017
##'
##' @rdname MuChPoint-class
##'
##' @exportClass MuChPoint
##'
setClass(
  Class="MuChPoint",
  representation=representation(
    S="numeric",
    N="numeric",
    bt="matrix")
)

### Fonction affichage
##' Print for the class of object returned by the \code{MuChPoint} function.
##'
##' @param x an object with class \code{MuChPoint}
##' @param N a numeric between 1 and length(x@@N) for the number of break-points desired.
##' @exportMethod print
##' @rdname print-MuChPoint
setMethod("print", "MuChPoint", definition =
            function(x, N=NULL) {
              cat("Multiple Change-Point (MuChPoint)\n")
              cat("- max number of change-points :", length(methods::slot(x,"N")),"\n")
              if (!is.null(N)){
                if (is.numeric(N)){
                  if ((N<1)||(N>nrow(methods::slot(x,"bt")))){
                    warning("N must be a positive integer less than ",nrow(methods::slot(x,"bt")))
                  }else{
                    cat("- position of change-points for N equal to",as.character(N),":", methods::slot(x,"bt")[N,1:N], "\n")
                  }
                }else{
                  warning("N must be a positive integer less than ",nrow(methods::slot(x,"bt")))
                }
              }
              cat("Use 'plot' to represent the results with 'shiny=TRUE' to select the number of break points \n")
              invisible(x)
            }
)

##' Summary of a \code{MuChPoint} object.
##'
##' @param object an object of class \code{MuChPoint}.
##'
##' @seealso \code{\linkS4class{MuChPoint}}.
##' @rdname summary-MuChPoint
##'
##' @examples
##' require(MuChPoint)
##' mu=c(rep(c(rep(1,25),rep(0,25)),3))%*%t(rep(c(rep(0,25),rep(1,25)),3))
##' Y=matrix(rnorm(150^2,0,2),150)+mu+t(mu)
##' Y=as.matrix(Matrix::forceSymmetric(Y))
##' res=MuChPoint(Y)
##' summary(res)
##'
##' @exportMethod summary
setMethod("summary", "MuChPoint", definition =
            function(object) {
              cat("Multiple Change-Point (MuChPoint)\n")
              cat("- max number of change-points :", length(methods::slot(object,"N")),"\n")
              cat("- slot '@bt' contains the break-points for each L (use 'print' with 'L=...' to display the positions)\n")
              cat("Use 'plot' to represent the results with 'shiny=TRUE' to select the number of break points \n")
              invisible(object)
            }
)

##' @rdname MuChPoint-class
setMethod("show", "MuChPoint", definition =
            function(object) {print(object)}
)

##' Produce a plot of two-dimensional segmentation of a \code{MuChPoint} fit.
##'
##' @param x an object of class \code{MuChPoint}.
##' @param y used  for S4 compatibility represented the matrix (typically,
##' the matrix used in the program \code{\link{MuChPoint}}).
##' @param shiny for a representation with a shiny application.
##' @param col for the colors of the representations.
##' @param L the summarized matrix with L break-points (L can be a vector).
##' @param ask If \code{TRUE}, to hit will be necessary to see next plot.
##'
##' @seealso \code{\linkS4class{MuChPoint}}, \code{\link{capushe}}.
##' @rdname plot-MuChPoint
##'
##' @references Article: BRAULT V., OUADAH S., SANSONNET L. and LEVY-LEDUC C. Nonparametric
##' homogeneity tests and multiple change-point estimation for analyzing large Hi-C data matrices.
##' Journal of Multivariate Analysis, 2017
##'
##' @examples
##' require(MuChPoint)
##' mu=c(rep(c(rep(1,25),rep(0,25)),3))%*%t(rep(c(rep(0,25),rep(1,25)),3))
##' Y=matrix(rnorm(150^2,0,2),150)+mu+t(mu)
##' Y=as.matrix(Matrix::forceSymmetric(Y))
##' res=MuChPoint(Y)
##' plot(res,Y,L=5,shiny=FALSE)
##' plot(res,Y,L=1:5,shiny=FALSE,ask=FALSE)
##'
##' @exportMethod plot

setMethod(
  f="plot",
  signature="MuChPoint",
  definition=function(x,y,shiny=TRUE,col="Color",L=NULL,ask=TRUE){

    if (!is.matrix(y)){
      y=as.matrix(y)
    }

    if (nrow(y)!=ncol(y)){
      stop("y must have the same number of lines and columns")
    }

    if (!is.character(col)){
      stop("col must be a character")
    }else{
      if (length(col)==1){
        if (col=="GrayLevel"){
          mypalette=grDevices::gray(seq(0,1, length=256))
        }else{
          mypalette=c(grDevices::rgb((0:201)/201*200,(0:201)/201*200,255,maxColorValue = 255),"white",
                      grDevices::rgb(255,(0:201)/201*200,(0:201)/201*200,maxColorValue = 255)[201:0])
        }
      }else{
        mypalette=col
      }
    }
    brek=seq(min(y),max(y),length=length(mypalette)+1)
    if (!is.null(L)){
      if (!is.numeric(L)){
        stop("L must be a positive ingeter vector between 1 and Lmax")
      }
      if ((length(methods::slot(x,"bt"))==0)&&(any(L!=length(x@N)))){
        stop("bt is required for the evolution. In the function MuChPoint, use
             the option Lmax instead of the option N")
      }
      bt=x@bt
      if ((any(floor(L)!=L))||(any((L<1)||(L>nrow(x@bt))))){
        stop("L must be a positive ingeter vector between 1 and Lmax")
      }
      }else{
        L=max(1,floor(length(x@S)/10))
        bt=x@bt
      }
    if (!is.logical(ask)){
      warning("ask must be a logical. Enforcing ask=TRUE")
      ask=TRUE
    }
    if (!is.logical(shiny)){
      warning("shiny must be a logical. Enforcing shiny=TRUE")
      shiny=TRUE
    }
    if (shiny){
      ############### Version Shiny
      ui_MuChPoint = shiny::fluidPage(
        shiny::fluidRow(
          shiny::column(12,shiny::sliderInput("L","Number of breaks",min=1,max=length(x@S),step=1,value=L))),
        shiny::textOutput("result2"),
        shiny::tabsetPanel(type = "tabs",
                           shiny::tabPanel('Out of MuChPoint',
                                           shiny::checkboxInput("Gray", "Gray version", col=="GrayLevel"),
                                           shiny::plotOutput(outputId="MuChPoint", height="1000px")),
                           shiny::tabPanel('Capushe',shiny::fluidRow(

                             shiny::sidebarPanel("For DDSE and Djump",
                                                 shiny::sliderInput("L2","Number of penalties",min=1,
                                                                    max=length(x@S),
                                                                    value=c(1,min(length(x@S),
                                                                                  max(1,floor(nrow(y)*0.45)))),
                                                                    step=1)
                             ),
                             shiny::sidebarPanel("For DDSE",
                                                 shiny::sliderInput("pct","pct",min=0,max=1,value=0.1,
                                                                    step=0.01)
                             )),

                             shiny::mainPanel("Out of Capushe",
                                              shiny::column(12,shiny::plotOutput(outputId="Capushe")),
                                              shiny::column(12,shiny::plotOutput(outputId="Capushe2")))))

      )
      ########### Partie Server
      server_MuChPoint = function(input, output) {
        output$result2 <- shiny::renderText({
          n=nrow(y)*(nrow(y)+1)
          K=(2:(length(x@S)+1))*(3:(length(x@S)+2))
          AIC=which.max(x@S[1:length(x@S)]-K)
          BIC=which.max(x@S[1:length(x@S)]-K/2*n)
          paste0(AIC," break-points chosen by AIC and ",BIC," break-points chosen by BIC")

        })
        output$MuChPoint<- shiny::renderPlot({
          if (input$Gray){
            plot(x,y,shiny=FALSE,L=input$L,col="GrayLevel",ask=FALSE)
          }else{
            plot(x,y,shiny=FALSE,L=input$L,col=col,ask=FALSE)
          }
        })
        output$Capushe<- shiny::renderPlot({
          data=data.frame(matrix(c(input$L2[1]:input$L2[2],input$L2[1]:input$L2[2],
                                   input$L2[1]:input$L2[2],-x@S[input$L2[1]:input$L2[2]]),
                                 ncol=4))
          cap=capushe::DDSE(data,pct =input$pct)
          plot(cap,newwindow=FALSE)
        })
        output$Capushe2<- shiny::renderPlot({
          data=data.frame(matrix(c(input$L2[1]:input$L2[2],input$L2[1]:input$L2[2],
                                   input$L2[1]:input$L2[2],-x@S[input$L2[1]:input$L2[2]]),
                                 ncol=4))
          n=nrow(y)*(nrow(y)+1)
          cap=capushe::Djump(data,Careajump=sqrt(log(n)/n))
          plot(cap,newwindow=FALSE)
        })
      }
      ########### Affichage
      shiny::shinyApp(ui = ui_MuChPoint,server=server_MuChPoint)

    }else{
      par_temp=graphics::par()
      par_old=list(mfrow=par_temp$mfrow,oma=par_temp$oma)
      graphics::par(mfrow=c(1,2),oma=c(0,0,3,0))
      for (liter in L){
        N=bt[liter,1:liter]
        ty=sort(nrow(y)+1-N)
        tx=N
        ## Diplay of the original matrix
        graphics::image(1:nrow(y),1:nrow(y),t(y)[,nrow(y):1],xlab="",ylab="",xaxt="n",
                        yaxt="n",main="Original data",col=mypalette,breaks=brek)

        graphics::abline(v=tx+0.5,col="purple")
        graphics::abline(h=ty-0.5,col="purple")
        # Diplay of the estimated matrix
        emplz=diff(unique(c(0,N,nrow(y))))
        z=as.matrix(Matrix::bdiag(lapply(emplz,function(i){rep(1,i)})))
        nz=length(emplz)
        resum=as.matrix(t(z)%*%y%*%z/(emplz%*%t(emplz)))
        graphics::image(y = unique(c(1,ty,nrow(y)+1)),
                        x = unique(c(1,tx,nrow(y)+1)),
                        z = as.matrix(resum[,nz:1]),xlab="",ylab="",xaxt="n",yaxt="n",
                        main="Estimated matrix",col=mypalette)
        graphics::abline(v=tx,col="purple")
        graphics::abline(h=ty,col="purple")
        graphics::title(main=paste(as.character(length(N))," break-points", sep=""),outer=TRUE)
        graphics::par(ask=ask)
      }
      graphics::par(par_old)
    }
  })

