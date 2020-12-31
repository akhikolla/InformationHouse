# $Id: sdea.R 229 2020-07-04 13:39:18Z lao $


# Calculates super efficiency
sdea <- function(X,Y, RTS="vrs", ORIENTATION="in", DIRECT=NULL, param=NULL,
                 TRANSPOSE=FALSE, LP=FALSE, CONTROL=NULL)
{
   # Input is as for the method eff

    # Hvis data er en data.frame saa tjek om det er numerisk data og lav
    # dem i saa fald om til en matrix
    X <- tjek_data(X)
    Y <- tjek_data(Y)

   # Antal firmaer i data er K
   if ( TRANSPOSE )  {
      X <- t(X)
      Y <- t(Y)
      if ( !is.null(DIRECT) & is(DIRECT, "matrix") )
         DIRECT <- t(DIRECT)
   }
   K = dim(X)[1]
   if (LP)  {
      print(paste("K =",K))
      print(dim(X))
      print(dim(Y))
   }
   lambda <- matrix(0,K,K)
   # Superefficiens skal gemmes i en K-vektor som vi til en start 
   # saetter til NA.
   supereff = rep(NA,K)
   if ( !is.null(dimnames(X)[[1]]) )  {
      names(supereff) <- dimnames(X)[[1]]
   }
   rts <- c("fdh","vrs","drs","crs","irs","irs","add","fdh+","fdh++","fdh0")
   if ( missing(RTS) ) RTS <- "vrs" 
   if ( is.numeric(RTS) )  {
      RTS_ <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTS_
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) ) {
      print(paste("Unknown value for RTS:",RTS),quote=F)
      RTS <- "vrs"
      print(paste("Continues with RTS =",RTS),quote=F)
   }
   orientation <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   if ( !(ORIENTATION %in% orientation) ) {
      print(paste("Unknown value for ORIENTATION:",ORIENTATION),quote=F)
      ORIENTATION <- "in"
      print(paste("Continues with ORIENTATION =",ORIENTATION),quote=F)
   }

   direct <- NULL
   directMatrix <- FALSE
   if ( !is.null(DIRECT) )  {
      if ( is(DIRECT, "matrix") )
         directMatrix <- TRUE
      else
         direct <- DIRECT
   }

   for ( i in 1:K ) {
      if (LP) print(paste("=====>> Unit",i),quote=F)
      # For hver enhed i laver vi beregningen og saetter resultat paa
      # den i'te plads i supereff
      # Den forste brug af X og Y er data for den enhed der skal 
      # beregnes efficiens for, den i'te.
      # Den anden brug er ved XREF og YREF for at angive hvilken teknologi
      # der skal bruges, de definerer teknologien.
      if ( directMatrix ) direct <- DIRECT[i,]
      e <- dea(X[i,,drop=FALSE], Y[i,,drop=FALSE],
         RTS,ORIENTATION, XREF=X[-i,,drop=FALSE], YREF=Y[-i,,drop=FALSE],
         TRANSPOSE=FALSE, DIRECT=direct, param=param, LP=LP, CONTROL=CONTROL)
      supereff[i] <- e$eff
      # print(dim(lambda))
      if (LP) print(dim(e$lambda))
      lambda[i,-i] <- e$lambda[1,]
   }
   if (LP) print("sdea: faerdig med gennemloeb")

# print(colnames(X))
   if ( is.null(rownames(X)) )  {
      colnames(lambda) <- paste("L",1:K,sep="")
   } else {
       colnames(lambda) <- paste("L",rownames(X),sep="_")
   }
   rownames(lambda) <- rownames(X)
   if (TRANSPOSE)  {
      lambda <- t(lambda)
      if ( !is.null(DIRECT) & is(DIRECT, "matrix") )
         DIRECT <- t(DIRECT)
   }

   if (LP) {
      print("sdea: Om lamda, dim og lambda")
      print(dim(lambda))
      print(rownames(lambda))
      print(colnames(lambda))
      print(lambda)
   }

   # return(supereff)
   objval <- NULL
   oe <- list(eff=supereff, lambda=lambda, objval=objval, RTS=RTS,
              ORIENTATION=ORIENTATION, TRANSPOSE=TRANSPOSE,
              slack=NULL, sx=NULL, sy=NULL)
   class(oe) <- "Farrell"
   return(oe)
}  # sdea

