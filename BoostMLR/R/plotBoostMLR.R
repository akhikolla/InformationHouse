#--------------------------------------------------
# Colors for plotting
ColorPlot <- c("black","red3","green4","blue","pink4","magenta2","orange")

#--------------------------------------------------
# Plot rho, phi, or error rate
plotBoostMLR <- function(Result,xlab = "",ylab = "",legend_fraction_x = 0.10,legend_fraction_y = 0,...){
  if(is.vector(Result)){
    Result <- matrix(Result,ncol = 1)
  }
  nr <- nrow(Result)
  nc <- ncol(Result)
  Range_x <- c(0,nr)
  x_lim <- c(0, (nr + round(legend_fraction_x*nr,0)) )
  Range_y <- range(c(Result),na.rm = TRUE)
  y_lim <- c(Range_y[1], ifelse(Range_y[2] > 0, Range_y[2] + legend_fraction_y*Range_y[2], Range_y[2] + legend_fraction_y*abs(Range_y[2])))
  plot(1:nr,Result[,1],type = "n",xlab = xlab,ylab = ylab,xlim = x_lim,ylim = y_lim,cex.lab = 1.5,frame.plot=FALSE,tck = 0,cex.axis = 1,xaxt = "n",yaxt = "n")
  axis(1,lwd = 2,at = seq(Range_x[1],Range_x[2],length.out = 5),tick = TRUE,cex = 0.5,font.axis = 2)
  axis(2,lwd = 2,at = seq(Range_y[1],Range_y[2],length.out = 5),labels = round(seq(Range_y[1],Range_y[2],length.out = 5),2),tick = TRUE,cex = 0.5,font.axis = 2)
  for(i in 1:nc){
    lines(1:nr,Result[,i],type = "l",lwd = 2,col = ColorPlot[i])
  }
# q3 <-  quantile(Result,probs = 0.75)
# legend(nr,q3,col = ColorPlot[1:nc],lwd = 3,lty = 1,box.col = "white",colnames(Result))
  legend(Range_x[2],Range_y[2],col = ColorPlot[1:nc],lwd = 3,lty = 1,box.col = "white",colnames(Result))
}
