## 1v PCA for quantiles -------
# Principal components analysis  of histogram variable based on Wasserstein distance----
#' Principal components analysis  of histogram variable based on Wasserstein distance
#' @description The function implements a Principal components analysis of histogram variable 
#' based on Wasserstein distance. It performs a centered (not standardized) PCA on a set of quantiles of a variable.
#' Being a distribution a multivalued description, the analysis performs a dimensional reduction and a visualization of distributions. 
#' It is a 1d (one dimension) becuse it is considered just one histogram variable.
#' @param data A MatH object (a matrix of distributionH).
#' @param var An integer, the variable number.
#' @param quantiles An integer, it is the number of quantiles used in the analysis.
#' @param plots a logical value. Default=TRUE plots are drawn.
#' @param listaxes A vector of integers listing the axis for the 2d factorial reperesntations.
#' @param axisequal A logical value. Default TRUE, the plot have the same scale for the x and the y axes.
#' @param qcut  a number between 0.5 and 1, it is used for the plot of densities, and avoids very peaked densities. 
#' Default=1, all the densities are considered.
#' @param outl a number between 0 (default)  and 0.5. For each distribution, is the amount of mass removed from the tails 
#' of the distribution. For example, if 0.1, from each distribution is cut away a left tail and a right one each containing 
#' the 0.1 of mass.   

#' @return a list with the results of the PCA in the MFA format of package \pkg{FactoMineR} for function MFA
#' @references Verde, R.; Irpino, A.; Balzanella, A., "Dimension Reduction Techniques for Distributional Symbolic Data," Cybernetics, IEEE Transactions on , vol.PP, no.99, pp.1,1
#' doi: 10.1109/TCYB.2015.2389653
#' keywords: {Correlation;Covariance matrices;Distribution functions;Histograms;Measurement;Principal component analysis;Shape;Distributional data;Wasserstein distance;principal components analysis;quantiles},
#' \url{http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7024099&isnumber=6352949}
#' @details In the framework of symbolic data analysis (SDA), distribution-valued data 
#' are defined as multivalued data, where each unit is described by a distribution 
#' (e.g., a histogram, a density, or a quantile function) of a quantitative variable.
#'  SDA provides different methods for analyzing multivalued data. Among them, the most relevant 
#'  techniques proposed for a dimensional reduction of multivalued quantitative variables is principal 
#'  component analysis (PCA). This paper gives a contribution in this context of analysis. 
#'  Starting from new association measures for distributional variables based on a peculiar metric 
#'  for distributions, the squared Wasserstein distance, a PCA approach is proposed for distribution-valued data,
#'   represented by quantile-variables. 
#' @examples
#' results=WH.1d.PCA(data = BLOOD,var = 1, listaxes=c(1:2))
#' @importFrom FactoMineR PCA
#' @importFrom graphics plot abline axis hist plot.new plot.window polygon segments text title
#' @importFrom grDevices dev.new rgb
#' @importFrom stats density quantile
#' @export
WH.1d.PCA=function(data,var, quantiles=10, plots=TRUE, listaxes=c(1:4),axisequal=FALSE,qcut=1,outl=0){
  if (is(data)[1]!="MatH"){stop("Input data must be a MatH object (A matrix of histograms)")}
  #Build the matrix of Quantiles
  VARS=ncol(data@M)
  INDIV=nrow(data@M)
  if ((qcut<0.5)||(qcut>1)){qcut=1}
  if ((outl<0)||(outl>0.5)){outl=0}
  if ((outl==0.5)){outl=0.49}
  if (missing(var)){
    var=1
    varname=colnames(data@M)[1]
    cat(paste("Var is missing, We do a PCA on variable ", varname, "\n"))
  } else{
    if (length(var)>1){
      varname=colnames(data@M)[var[1]]
      cat(paste("Var has several values, We do a PCA on variable -->", var[1], "\n"))
      var=var[1]
    }
    else{
      if (var>VARS){stop(paste("The variables are less than ",var))}
      else{
        varname=colnames(data@M)[var]
        cat(paste("We do a PCA on variable ---> ", varname, "\n"))
      }
    }
  }
  data=data[,var]
  #check indiv and remove empty rows
  toremove=numeric(0)
  for (i in 1:INDIV){
    if (length(data@M[i,1][[1]]@x)<2){
      toremove=c(toremove,i)
    }
  }
  if (length(toremove)>0){
    data@M=as.matrix(data@M[-toremove,1])
  }
  colnames(data@M)=varname
  INDIV=nrow(data@M)
  if (INDIV<2){stop("At least two individuals are necessary, please check data")}
  
  #Create the quantile matrix
  MATQ=matrix(0,nrow = INDIV,ncol = (quantiles+1))
  if (outl==0){
    p=c(0, 1:quantiles)/quantiles
    namec=list("Min")
    for (j in 2:(quantiles)){
      namec=c(namec,paste0("Q(",format(p[j], digits=2, nsmall=2),")"))
    }
    namec=c(namec, "Max")
  }
  else{
    
    p=c(0, 1:quantiles)/quantiles*(1-2*outl)+outl 
    namec=list(paste0("Q(",format(p[1], digits=2, nsmall=2),")"))
    for (j in 2:length(p)){
      namec=c(namec,paste0("Q(",format(p[j], digits=2, nsmall=2),")"))
    }
  }
  
  for (i in 1:INDIV)
  {
    for (j in 1:(quantiles+1))
    {MATQ[i,j]=compQ(data@M[i,1][[1]],p[j])}
    
  }
  rownames(MATQ)=rownames(data@M)
  colnames(MATQ)=namec
  
  # do a PCA
  vmeans=numeric(0)
  vstds=numeric(0)
  
  for (ind in 1:INDIV){
    vmeans=c(vmeans,data@M[ind,1][[1]]@m)
    vstds=c(vstds,data@M[ind,1][[1]]@s)
  }
  minM=min(vmeans)
  maxM=max(vmeans)
  minS=min(vstds)
  maxS=max(vstds)
  
  MATF=cbind(MATQ/sqrt((quantiles+1)),vmeans,vstds)
  res.pca = PCA(MATF,quanti.sup = c((ncol(MATF)-1),ncol(MATF)), scale.unit=FALSE,  graph=F)
  
  VARW=0#WH.var.covar(data)
  TOTINE=sum(res.pca$eig[,1])
  ## plotting PCA results ----------------
  if (plots){
    #par(mfrow=c(1,1))
    # lt's define the couple of axes
    dev.new(noRStudioGD = FALSE )
    planes=matrix(0,1,2)
    for (cc in 1:floor(length(listaxes)/2)){
      planes=rbind(planes,c(listaxes[(cc*2-1)],listaxes[(cc*2)]))
    }
    planes=as.matrix(planes[-1,])
    dim(planes)=c(floor(length(listaxes)/2),2)
    if ((length(listaxes)%%2)==1){
      planes=rbind(planes, c(listaxes[(length(listaxes)-1)],
                             listaxes[length(listaxes)]))
    }
    
    #plot Spanish-fan plot
    for (pl in 1:nrow(planes)){
      axe1=planes[pl,1];axe2=planes[pl,2]
      
      labX=paste("Axis ",axe1," (", format(res.pca$eig[axe1,2], digits=2, nsmall=2), "%)")
      labY=paste("Axis ",axe2," (", format(res.pca$eig[axe2,2], digits=2, nsmall=2), "%)")
      CVAR=res.pca$var$coord[,c(axe1,axe2)]
      plot.new()
      
      if (axisequal){
        xrange=c(min(0,min(cbind(CVAR[,1],CVAR[,2]))), max(0,max(cbind(CVAR[,1],CVAR[,2]))))
        yrange=xrange
        plot.window(xrange, yrange)
        title(xlab=labX)
        title(ylab=labY)
        
        #         plot(type="n",
        #              xlab=labX,ylab=labY)
        segments(xrange[1],0,xrange[2],0)
        segments(yrange[1],0,yrange[2],0)
      }
      else{
        plot.window(c(min(0,CVAR[,1]),max(0,CVAR[,1])),
                    c(min(0,CVAR[,2]),max(0,CVAR[,2])))
        title(xlab=labX)
        title(ylab=labY)
        #       plot(c(min(0,CVAR[,1]),max(0,CVAR[,1])), c(min(0,CVAR[,2]),max(0,CVAR[,2])),type="n",
        #            xlab=labX,ylab=labY)Y)
        segments(min(0,CVAR[,1]),0,max(0,CVAR[,1]),0)
        segments(0,min(0,CVAR[,2]),0,max(0,CVAR[,2]))}
      title("PCA Variable plot (Spanish-fan plot)")
      if ((nrow(CVAR)%%2)==0){centr=nrow(CVAR)/2} else{centr=(nrow(CVAR)-1)/2}
      centr=nrow(CVAR)/2
      cxl=0.8
      for (tr in 2:nrow(CVAR)){
        centrality=(abs(tr-centr-1)/(centr-1))
        
        #cat(centrality,"\n")
        x=c(0, CVAR[(tr-1),1], CVAR[tr,1], 0)
        y=c(0, CVAR[(tr-1),2], CVAR[tr,2], 0)
        red=1
        green=centrality
        #cat(red,green,"\n")
        polygon(x, y, col=rgb(red, green, 0,0.7), lty = 1, lwd = 1,  border = "black")
      }
      text(x = CVAR[1,1],y= CVAR[1,2],labels = colnames(MATQ)[1],cex=cxl)
      if (nrow(CVAR)<8){
        for (tr in 2:(nrow(CVAR)-1)){
          text(x = CVAR[tr,1],y= CVAR[tr,2],labels = colnames(MATQ)[tr],cex=cxl)
        }}
      else{
        
        text(x = CVAR[ceiling(nrow(CVAR)/8),1],y= CVAR[ceiling(nrow(CVAR)/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)/8)],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*2/8),1],y= CVAR[ceiling(nrow(CVAR)*2/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*2/8)],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*3/8),1],y= CVAR[ceiling(nrow(CVAR)*3/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*3/8)],cex=cxl)
        text(x = CVAR[centr+1,1],y= CVAR[centr+1,2],labels = colnames(MATQ)[centr+1],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*5/8),1],y= CVAR[ceiling(nrow(CVAR)*5/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*5/8)],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*6/8),1],y= CVAR[ceiling(nrow(CVAR)*6/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*6/8)],cex=cxl)
        text(x = CVAR[ceiling(nrow(CVAR)*7/8),1],y= CVAR[ceiling(nrow(CVAR)*7/8),2],
             labels = colnames(MATQ)[ceiling(nrow(CVAR)*7/8)],cex=cxl)
        
      }
      text(x = CVAR[nrow(CVAR),1],y= CVAR[nrow(CVAR),2],
           labels = colnames(MATQ)[ncol(MATQ)],cex=cxl)
      axis(1, pos = 0, cex.axis=0.7)#, at=xM,labels=labX)
      axis(2, pos = 0, cex.axis=0.7)#, at=yM,labels=labY)
      dev.new(noRStudioGD = FALSE )
      
      #plot individuals
      
      #layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
      
      xm=min(res.pca$ind$coord[,axe1])
      xM=max(res.pca$ind$coord[,axe1])
      ym=min(res.pca$ind$coord[,axe2])
      yM=max(res.pca$ind$coord[,axe2])
      Xlen=xM-xm
      Ylen=yM-ym
      
      #Compute scaling factor for domain
      matX=matrix(0,128,INDIV)
      matD=matrix(0,128,INDIV)
      for (ind in 1:INDIV){
        distr=data@M[ind,1][[1]]
        #generate 200 random points according to the QF
        rn=200
        xn=c(rep(0,rn))
        random_no=c(0:rn)/rn
        
        for (i in 1:rn){
          xn[i]=compQ(distr,random_no[i])
        }
        d <- density(xn,n=128)
        matX[,ind]=d$x
        matD[,ind]=d$y
      }
      MinX=min(matX)
      MaxX=max(matX)
      RX=MaxX-MinX
      q95D=quantile(matD,probs = qcut)
      matD[matD>q95D]=as.numeric(q95D)
      MaxD=max(matD)
      
      meanXRange=diff(range(apply(matX,2,FUN = range)))
      xfact=0.3
      yfact=0.1
      if (axisequal){xm1=min(xm,ym)
      xM1=max(xM,yM)
      ym1=xm1
      yM1=xM1
      Xlen=xM1-xm1
      Ylen=yM1-ym1
      xm=(xM1+xm1)/2-Xlen/2*(1+1.5*xfact)
      xM=(xM1+xm1)/2+Xlen/2*(1+1.5*xfact)
      ym=(yM1+ym1)/2-Ylen/2*(1+0.5*yfact) 
      yM=(yM1+ym1)/2+Ylen/2*(1+1.5*yfact)
      
      }else{
        xm1=xm
        xM1=xM
        ym1=ym
        yM1=yM
        xm=(xM1+xm1)/2-Xlen/2*(1+1.5*xfact)
        xM=(xM1+xm1)/2+Xlen/2*(1+1.5*xfact)
        ym=(yM1+ym1)/2-Ylen/2*(1+0.5*yfact) 
        yM=(yM1+ym1)/2+Ylen/2*(1+1.5*yfact)     
      }
      dev.new(noRStudioGD = FALSE )
      plot.new()
      plot.window(c(xm,xM),c(ym,yM))
      title(main="PCA plot of distributions", sub="Smoothed distributions")
      for (ind in 1:INDIV){
        x=matX[,ind]
        y=matD[,ind]
        x=(x-data@M[ind,1][[1]]@m)/meanXRange*xfact*Xlen+res.pca$ind$coord[ind,axe1]
        x=c(x,x[length(x)],x[1])
        y=c(y,0,0)/q95D*yfact*Ylen+res.pca$ind$coord[ind,axe2]
        polygon(x,y, col=rgb(1, 1,0,0.6))
        
      }
      for (ind in 1:INDIV){
        text(res.pca$ind$coord[ind,axe1],res.pca$ind$coord[ind,axe2],
             label=rownames(MATQ)[ind], pos=1, cex=0.7)
      }
      axis(1, pos = 0, cex.axis=0.7)#, at=xM,labels=labX)
      axis(2, pos = 0, cex.axis=0.7)#, at=yM,labels=labY)
      segments(xm,0,xM,0)
      segments(0,ym,0,yM)
      
      title(xlab=labX)
      title(ylab=labY)
      dev.new(noRStudioGD = FALSE )
      # plot.new()
      plot.new()
      plot.window(c(xm,xM),c(ym,yM))
      title(main="Coloured using the mean values")
      trasp=0.6
      
      for (ind in 1:INDIV){
        x=matX[,ind]
        y=matD[,ind]
        x=(x-data@M[ind,1][[1]]@m)/meanXRange*xfact*Xlen+res.pca$ind$coord[ind,axe1]
        x=c(x,x[length(x)],x[1])
        y=c(y,0,0)/q95D*yfact*Ylen+res.pca$ind$coord[ind,axe2]
        red=(data@M[ind,1][[1]]@m-minM)/(maxM-minM)
        if (red>0.5){
          green=1-red
          blue=0
        }
        else{green=1-red
        red=0
        blue=1-green
        }
        polygon(x,y, col=rgb(red,green,blue,trasp))
        
      }
      for (ind in 1:INDIV){
        text(res.pca$ind$coord[ind,axe1],res.pca$ind$coord[ind,axe2],
             label=rownames(MATQ)[ind], pos=1, cex=0.7)
      }
      axis(1, pos = 0, cex.axis=0.7)#, at=xM,labels=labX)
      axis(2, pos = 0, cex.axis=0.7)#, at=yM,labels=labY)
      segments(xm,0,xM,0)
      segments(0,ym,0,yM)
      
      abline(0,res.pca$quanti.sup$coord[1,2]/res.pca$quanti.sup$coord[1,1],lty=3,lwd=2, col="red")
      title(xlab=labX)
      title(ylab=labY)
      dev.new(noRStudioGD = FALSE )
      plot.new()
      plot.window(c(xm,xM),c(ym,yM))
      title(main="Coloured using std values")
      for (ind in 1:INDIV){
        
        x=matX[,ind]
        y=matD[,ind]
        x=(x-data@M[ind,1][[1]]@s)/meanXRange*xfact*Xlen+res.pca$ind$coord[ind,axe1]
        x=c(x,x[length(x)],x[1])
        y=c(y,0,0)/q95D*yfact*Ylen+res.pca$ind$coord[ind,axe2]
        red=(data@M[ind,1][[1]]@s-minS)/(maxS-minS)
        if (red>0.5){
          green=1-red
          blue=0
        }
        else{green=1-red
        red=0
        blue=1-green
        }
        polygon(x,y, col=rgb(red,green,blue,trasp))
        
      }
      for (ind in 1:INDIV){
        text(res.pca$ind$coord[ind,axe1],res.pca$ind$coord[ind,axe2],
             label=rownames(MATQ)[ind], pos=1, cex=0.7)
      }
      axis(1, pos = 0, cex.axis=0.7)#, at=xM,labels=labX)
      axis(2, pos = 0, cex.axis=0.7)#, at=yM,labels=labY)
      segments(xm,0,xM,0)
      segments(0,ym,0,yM)
      
      abline(0,res.pca$quanti.sup$coord[2,2]/res.pca$quanti.sup$coord[2,1],lty=3,lwd=2, col="red")
      
      title(xlab=labX)
      title(ylab=labY)
    }
  }
  return(list(PCAout=res.pca, WASSVARIANCE=VARW, INERTIA=TOTINE))
}
# Principal components analysis  of a set of histogram variables based on Wasserstein distance----
#' Principal components analysis  of a set of histogram variable based on Wasserstein distance
#' @description (Beta version) The function implements a Principal components analysis of a set of histogram variables 
#' based on Wasserstein distance. It performs a centered (not standardized) PCA on a set of quantiles of a variable.
#' Being a distribution a multivalued description, the analysis performs a dimensional reduction and a visualization of distributions. 
#' It is a 1d (one dimension) becuse it is considered just one histogram variable.
#' @param data A MatH object (a matrix of distributionH).
#' @param list.of.vars A list of  integers, the active variables.
#' @param quantiles An integer, it is the number of quantiles used in the analysis.
#' Default=10.
#' @param outl a number between 0 (default)  and 0.5. For each distribution, is the amount of mass removed from the tails 
#' of the distribution. For example, if 0.1, from each distribution is cut away a left tail and a right one each containing 
#' the 0.1 of mass.
#' @return a list with the results of the PCA in the MFA format of package \pkg{FactoMineR} for function MFA
#' @details It is an extension of WH.1d.PCA to the multiple case. 
#' @importFrom FactoMineR MFA
#' @importFrom stats rnorm
#' @export
## MFA PCA for quantiles -------
WH.MultiplePCA=function(data, list.of.vars, quantiles=10, outl=0){
  data=data[,list.of.vars]
  VARS=ncol(data@M)
  INDIV=nrow(data@M)
  MATQ=matrix(0,nrow = INDIV,ncol = VARS*(quantiles+1))
  if ((outl<0)||(outl>0.5)){outl=0}
  if ((outl==0.5)){outl=0.49}
 # browser()
  #create the names of the variables
  #################################################
  COLnames=list()
  if (outl==0){
    p=c(0, 1:quantiles)/quantiles
    for (v in 1:VARS){
      namec=list("Min")
      for (j in 2:(quantiles)){
        namec=c(namec,paste0("Q(",format(p[j], digits=2, nsmall=2),")"))
      }
      namec=c(namec, "Max")
      COLnames=c(COLnames, paste(substr(colnames(data@M)[v],1,4),namec,sep="."))
      for (i in 1:INDIV)
      {
       # tmp=compQ_vect(data@M[i,v][[1]],p)
        for (j in 1:(quantiles+1)){
          MATQ[i,((v-1)*(quantiles+1)+j)]=compQ(data@M[i,v][[1]],p[j])
        }
      }
    }
    
  }
  else{
    
    p=c(0, 1:quantiles)/quantiles*(1-2*outl)+outl 
    namec=list(paste0("Q(",format(p[1], digits=2, nsmall=2),")"))
    for (j in 2:length(p)){
      namec=c(namec,paste0("Q(",format(p[j], digits=2, nsmall=2),")"))
    }
    for (v in 1:VARS){
      namec=list(paste0("Q(",format(p[1], digits=2, nsmall=2),")"))
      for (j in 2:length(p)){
        namec=c(namec,paste0("Q(",format(p[j], digits=2, nsmall=2),")"))
      }
      COLnames=c(COLnames, paste(substr(colnames(data@M)[v],1,4),namec,sep="."))
      for (i in 1:INDIV)
      {
        for (j in 1:(quantiles+1)){
          MATQ[i,((v-1)*(quantiles+1)+j)]=compQ(data@M[i,v][[1]],p[j])
        }
      }
    }
  }
  #################################################
  # COLnames=list();
  # for (v in 1:VARS){
  #   namec=list("Min")
  #   for (j in 2:(quantiles)){
  #     namec=c(namec,paste0("Q(",format(p[j], digits=2, nsmall=2),")"))
  #   }
  #   namec=c(namec, "Max")
  #   COLnames=c(COLnames, paste(substr(colnames(data@M)[v],1,4),namec,sep="."))
  #   for (i in 1:INDIV)
  #   {
  #     for (j in 1:(quantiles+1)){
  #       MATQ[i,((v-1)*(quantiles+1)+j)]=compQ(data@M[i,v][[1]],p[j])
  #     }
  #   }
  # }
  
  rownames(MATQ)=rownames(data@M)
  colnames(MATQ)=COLnames
  
  # do a MFA
  MATF=cbind(MATQ/sqrt((quantiles+1)))
  #check all columns
  for(cc in 1:ncol(MATF)){
    if (sd(MATF[,cc])<1e-100){
      MATF[,cc]=MATF[,cc]+rnorm(nrow(MATF),0,1e-20)
    }
  }
  MFA.res= MFA(MATF,group=c(rep(quantiles+1,length(list.of.vars))),type=c(rep("c",length(list.of.vars))),
               ncp=3,name.group=colnames(data@M),graph=TRUE)
  # do Plots!!!!
  #plot Spanish-fan plot
  #plot individuals
  return(MFA.res)
}

# Plotting Spanish fun plots for Multiple factor analysis of Histogram Variables ----
#' Plotting Spanish fun plots for Multiple factor analysis of Histogram Variables
#' @description The function plots the circle of correlation of the quantiles
#' of the histogrma variables after a Multiple factor analysis.
#' @param res Results from WH.MultiplePCA, or WH.1D.PCA.
#' @param axes A list of  integers, the new factorial axes c(1,2) are the default.
#' @param var A list of integers are the variables to plot.
#' @param LABS Logical, if TRUE graph is labeled, otherwise it does not.
#' @param multi Logical, if TRUE (default) results come from a WH.MultiplePCA,
#'  if FALSE results come from WH.1D.PCA.
#' @param corplot Logical, if TRUE (default) the plot reports correlations, if FALSE
#' the coordinates of quantiles on the factorial plane
#' @return a plot of class ggplot
#' @examples
#' #Do a MultiplePCA on the BLOOD dataset
#' \dontrun{
#' res=WH.MultiplePCA(BLOOD,list.of.vars = c(1:3)) 
#' }
#' #Plot results
#' \dontrun{
#' WH.plot_multiple_Spanish.funs(res,axes=c(1,2),var=c(1:3))
#' }
#' @export
WH.plot_multiple_Spanish.funs=function(res,
                                       axes=c(1,2),var=1,LABS=TRUE,
                                       multi=TRUE,corplot=TRUE){
  #require(ggplot2)
  if (multi){
    labX=paste("Axis ",axes[1]," (", format(res$eig[axes[1],2], 
                                            digits=2, nsmall=2), "%)")
    labY=paste("Axis ",axes[2]," (", format(res$eig[axes[2],2], 
                                            digits=2, nsmall=2), "%)")
  }
  else
  {
    labX=paste("Axis ",axes[1]," (", format(res$PCAout$eig[axes[1],2], 
                                            digits=2, nsmall=2), "%)")
    labY=paste("Axis ",axes[2]," (", format(res$PCAout$eig[axes[2],2], 
                                            digits=2, nsmall=2), "%)")
  }
  v1=numeric()
  v2=numeric()
  tt=character()
  tt2=numeric()
  VA=numeric()
  x1=numeric()
  x2=numeric()
  hj=numeric()
  vj=numeric()
  if (!multi){var=1}
  for (j in var){
    
    #plot Spanish-fan plot
    if (multi){
      els=which(res$summary.quanti$group==j)
    } 
    else{
      els=c(1:nrow(res$PCAout$var$coord))
    }
    tmp_v1=numeric()
    tmp_v2=numeric()
    tmp_tt=character()
    tmp_tt2=numeric()
    tmp_VA=numeric()
    if (multi){
      if (corplot){
        CVAR=res$global.pca$var$cor[els,c(axes[1],axes[2])]
      }
      else{
        CVAR=res$global.pca$var$coord[els,c(axes[1],axes[2])]
      }
    }
    else{
      
      if (corplot){
        CVAR=res$PCAout$var$cor[els,c(axes[1],axes[2])]}
      else{
        CVAR=res$PCAout$var$coord[els,c(axes[1],axes[2])]
      }
    }
    
    for (pp in 1:(nrow(CVAR)-1)){
      tmp_v1=c(tmp_v1,CVAR[pp,1],0,CVAR[pp+1,1],CVAR[pp,1])
      tmp_v2=c(tmp_v2,CVAR[pp,2],0,CVAR[pp+1,2],CVAR[pp,2])
      tmp_tt=c(tmp_tt,paste(as.character(j),rep(as.character(pp),4),sep="_"))
      tmp_tt2=c(tmp_tt2,rep(pp,4))
      tmp_VA=c(tmp_VA,rep(j,4))
    }
    tmp_tt2=abs(tmp_tt2-mean(tmp_tt2))
    tmp_tt2=1-tmp_tt2/max(tmp_tt2)
    v1=c(v1,tmp_v1,0,0,0,0)
    v2=c(v2,tmp_v2,0,0,0,0)
    tt=c(tt,tmp_tt,paste(as.character(j),rep(as.character(pp),4),sep="_"))
    tt2=c(tt2,tmp_tt2,1,1,1,1)
    tmp_VA=c(tmp_VA,rep(j,4))
    VA=c(VA,tmp_VA)
    x1=c(x1,CVAR[,1])
    x2=c(x2,CVAR[,2])
  }
  
  d=data.frame(v1=v1,v2=v2,tt=tt,tt2=tt2,VA=VA)
  #CVAR=res$global.pca$var$cor[,c(axes[1],axes[2])]
  txt_lab=data.frame(v1=x1,v2=x2,hj=(-sign(x1)+1)/2,vj=(-sign(x2)+1)/2)
  ### the circle of correlation
  angle <- seq(-pi, pi, length = 250)
  x_c = sin(angle)
  y_c = cos(angle)
  df <- data.frame(x_c = x_c, y_c = y_c)
  ##
  l1=data.frame(x=c(-1.2,1.2),y=c(0,0))
  l2=data.frame(x=c(0,0),y=c(-1.2,1.2))
  
  
  p <- ggplot(d, aes(x=v1, y=v2)) + 
    geom_polygon(aes(group=tt, fill=as.factor(VA),alpha=tt2),colour='grey30')
  if (LABS==TRUE){
    p=p+geom_text(data=txt_lab,aes(x=v1,y=v2,
                                   label=rownames(txt_lab),
                                   hjust=hj,vjust=vj))
  }
  if (corplot){
    p=p+xlim(-1.2,1.2)+ylim(-1.2,1.2)+
      ggtitle("Correlation plot of variables (Spanish-fan plot)")
  }
  else{
    XL=range(d$v1)
    XL=(XL-mean(XL))*1.2+mean(XL)
    YL=range(d$v2)
    YL=(YL-mean(YL))*1.2+mean(YL)
    p=p+xlim(XL)+ylim(YL)+
      ggtitle("Plot of variables (Spanish-fan plot)")
  }
  
  p=p+theme_bw()+
    xlab(labX)+ylab(labY)+
    theme(legend.position="none")
  if (corplot){
    p=p+coord_fixed(ratio = 1)}
  p=p+theme(plot.title = element_text(hjust = 0.5))
  if (corplot){
    p=p+geom_path( data = df, aes(x=x_c, y=y_c), inherit.aes = F,linetype = 2)
  }
  p=p+geom_vline(xintercept = 0)+geom_hline(yintercept = 0)
  
  
  print(p)    
  
  return(p)
}

# Plot histograms of individuals after a  Multiple factor analysis of Histogram Variables ----
#'  Plot histograms of individuals after a  Multiple factor analysis of Histogram Variables
#' @description (Beta version) The function plots histogram data of the individuals for a particular
#' variable on a factorial palne after a Multiple factor analysis.
#' @param data a MatH object
#' @param res Results from WH.MultiplePCA.
#' @param axes A list of  integers, the new factorial axes c(1,2) are the default.
#' @param indiv A list of objects (rows) of data to plot. Default=0 all the objects of data.
#' @param var An integer indicating an original histogrma variable to plot.
#' @param strx a resizing factor for the domain of histograms (default=0.1 means 
#' that each distribution has a support that is one tenth of the spread of the x axis)
#' @param stry a resizing factor for the density of histograms (default=0.1 means 
#' that each distribution has a density that is one tenth of the spread of the y axis)
#' @param HISTO a logical value. Default=TRUE plots histograms, FALSE plot smooth densities.
#' @param coor (optional) if 0 (Default) takes the coordinates in res, if a a matrix is passed the coordinates are those passed
#' @param stat (optional) if 'mean'(Default) a plot of individuals labeled by the means is produced.
#' Otherwise if 'std', 'skewness' or 'kurtosis', data are labeled with this statistic. 
#' @return a plot of class ggplot
#' @examples
#' #Do a MultiplePCA on the BLOOD dataset
#' \dontrun{
#' #' results=WH.MultiplePCA(BLOOD,list.of.vars = c(1:3)) 
#' #Plot histograms of variable 1 of BLOOD dataset on the first 
#' #factorial plane showing histograms
#' WH.plot_multiple_indivs(BLOOD,results,axes=c(1,2),var=1,strx=0.1,
#'  stry=0.1, HISTO=TRUE)
#' #Plot histograms of variable 1 of BLOOD dataset on the first 
#' #factorial plane showing densities

#' WH.plot_multiple_indivs(BLOOD,results,axes=c(1,2),var=1,strx=0.1,
#'  stry=0.1, HISTO=FALSE)
#'  }

#' @export 
WH.plot_multiple_indivs=function(data,res,axes=c(1,2),indiv=0,var=1, 
                                 strx=0.1, stry=0.1, HISTO=TRUE, coor=0,
                                 stat="mean"){
  #require(ggplot2)
  
  if (indiv[1]==0){
    el_ind=c(1:get.MatH.nrows(data))
  } else{
    el_ind=indiv
  }
  if (length(coor)==1){
    coord=res$ind$coord[el_ind,axes]}else{
      coord=coor
    }
  # browser()
  
  if (HISTO){
    x=numeric()
    y=numeric()
    ID_O=numeric()  
    minx=numeric()
    mea=numeric()
    rmax=-Inf
    for (i in el_ind){
      
      tmp_r=max(data@M[i,var][[1]]@x)-data@M[i,var][[1]]@x[1]+1e-100
      mm=data@M[i,var][[1]]@m
      if (tmp_r>rmax){rmax=tmp_r}
      vv=COMP_POLY_H(data@M[i,var][[1]])
      x=c(x,min(vv$x),vv$x,max(vv$x))
      y=c(y,0,vv$y,0)
      ID_O=c(ID_O,rep(i,length(vv$x)+2))
      mea=c(mea,rep(mm,length(vv$x)+2))
      minx=c(minx,rep(min(vv$x),length(vv$x)+2))
    }
    
  }else
  { 
    # compute densities
    pt=300
    x=numeric()
    y=numeric()
    ID_O=numeric()  
    minx=numeric()
    mea=numeric()
    rmax=-Inf
    
    for (i in el_ind){
      vals=COMP_MQ(data@M[i,var][[1]],Mp=c(0:pt)/pt)
      mm=data@M[i,var][[1]]@m
      tmp_r=max(vals)-min(vals)
      if (tmp_r>rmax){rmax=tmp_r}
      vv=density(vals,n = 128,from=min(vals), to=max(vals))
      x=c(x,min(vv$x),vv$x,max(vv$x))
      y=c(y,0,vv$y,0)
      ID_O=c(ID_O,rep(i,length(vv$x)+2))
      mea=c(mea,rep(mm,length(vv$x)+2))
      minx=c(minx,rep(min(vals),length(vv$x)+2))
    }
    
    ### end of compute density
  }
  ID_o=ID_O
  x_o=x
  y_o=y
  ALL_OBJ=data.frame(x_o=x_o,y_o=y_o,ID_o=ID_o,minx=minx,mea=mea)
  limX=range(ALL_OBJ$x_o)
  
  EXTRY=quantile(ALL_OBJ$y_o,probs = 0.995)
  #bisogna traslare ogni oggetto e scalarlo
  #scaliamo gli oggetti tra 0 e 1
  ALL_OBJ$x=(ALL_OBJ$x_o-minx)/(rmax)
  #tagliamo picchi eccessivi
  ALL_OBJ$y[ALL_OBJ$y_o>EXTRY]=EXTRY
  limY=quantile(ALL_OBJ$y_o,probs = 0.99)
  ALL_OBJ$y=ALL_OBJ$y_o/limY
  #ogni oggetto va traslato e rimpicciolito di un fattore strech x strech y
  dim1MAX=max(coord[,1])
  dim1Min=min(coord[,1])
  rX=(dim1MAX-dim1Min)
  dim2MAX=max(coord[,2])
  dim2Min=min(coord[,2])
  rY=(dim2MAX-dim2Min)
  cc=0
  for (i in el_ind){
    cc=cc+1
    xc=coord[cc,1]
    yc=coord[cc,2]
    tmp_x=ALL_OBJ$x_o[which(ALL_OBJ$ID_o==i)]
    tmp_y=ALL_OBJ$y_o[which(ALL_OBJ$ID_o==i)]
    tmp_x=tmp_x*strx*rX
    tmp_y=tmp_y*stry*rY
    tmp_x=(tmp_x-mean(tmp_x))+xc
    tmp_y=tmp_y+yc
    ALL_OBJ$x_o[which(ALL_OBJ$ID_o==i)]=tmp_x
    ALL_OBJ$y_o[which(ALL_OBJ$ID_o==i)]=tmp_y
  }
  
  
  title=paste("Plot of individuals: distributions for variable",
              colnames(data@M)[var])
  labX=paste("Axis ",axes[1]," (", format(res$eig[axes[1],2], 
                                          digits=2, nsmall=2), "%)")
  labY=paste("Axis ",axes[2]," (", format(res$eig[axes[2],2], 
                                          digits=2, nsmall=2), "%)")
  
  colnames(coord)=c('x','y')
  coord=as.data.frame(coord)
  p=ggplot(data=ALL_OBJ,aes(x=x_o,y=y_o))+
    geom_polygon(aes(group=ID_o, fill=as.factor(ID_o),alpha=0.4),colour='black')+
    geom_point(data=coord,aes(x=x,y=y))+
    geom_text(data=coord,
              aes(x=x,y=y,label=rownames(coord),
                  hjust=0.5, vjust=1))+theme_bw()+
    xlab(labX)+ylab(labY)+
    ggtitle(title)+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_vline(xintercept = 0)+geom_hline(yintercept = 0)
  ### add a line for mean snd sd
  Ms=get.MatH.stats(data[el_ind,var],stat=stat)$mat
  # SD=get.MatH.stats(data[el_ind,var],stat="std")$mat
  # Skew=get.MatH.stats(data[el_ind,var],stat="std")$mat
  # kurtH()=get.MatH.stats(data[el_ind,var],stat="std")$mat
  ### print a plot with means
  
  coord=cbind(coord,Ms)
  colnames(coord)[[3]]=stat
  
  p2=ggplot(ALL_OBJ,aes(x=x_o,y=y_o))+
    scale_fill_gradient(low = "yellow", high = "red")+
    geom_polygon(aes(group=ID_o, fill=mea,alpha=0.4),colour='black')+
    geom_point(data=coord,aes(x=x,y=y))+theme_bw()+
    geom_text(data=coord,
              aes(x=x,y=y,
                  label=format(mean,digits=1,nsmall=2),
                  hjust=0.5, vjust=1))+
    xlab(labX)+ylab(labY)+
    ggtitle(paste(title,"each label is the ",stat,"." ))+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_vline(xintercept = 0)+geom_hline(yintercept = 0)
  
  #+ geom_abline(intercept = 0, slope = cor(coord)[1,3]*sd(coord$x),linetype=2)
  # browser()
  print(p)
  print(p2)
  
  return(pp=list(p,p2))
}

COMP_MQ=function(object,Mp){
  res=numeric()
  n=length(Mp)
  
  for (i in 1:n){
    p=Mp[i]
    #Check for errors
    if (p<0 || p>1) stop("p must be a value between 0 and 1")
    
    if (p<=0) {
      q=object@x[1]
    }else{
      if (p>=1) {q=object@x[length(object@x)]}
      else{
        ini=max(object@x[object@p<=p])
        pos1=which.max(object@x[object@p<=p])
        pos2=pos1+1
        fin=object@x[pos2];
        if (ini==fin){
          q=ini
        }
        else{
          q=ini+(object@x[pos2]-object@x[pos1])*(p-object@p[pos1])/(object@p[pos2]-object@p[pos1])
        }
      }
      res=c(res,q)
    }
  }
  return(res)
}

COMP_POLY_H=function(object){
  #browser()
  x=object@x[1]
  y=0
  if (length(object@x)<2){
    x=c(x,x,x)
    y=c(y,1,0)}
  else{
    for (i in 2:length(object@x)){
      tmp_dens=(object@p[i]-object@p[i-1])/(object@x[i]-object@x[i-1]+1e-100)
      x=c(x,object@x[i-1],object@x[i])
      y=c(y,tmp_dens,tmp_dens)
    }
    x=c(x,object@x[length(object@x)])
    y=c(y,0)
  }
  
  return(res=data.frame(x=x,y=y))
}