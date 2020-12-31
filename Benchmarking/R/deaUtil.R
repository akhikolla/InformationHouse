# $Id: deaUtil.R 218 2020-05-21 21:28:28Z lao $

# Naesten alle funktioner transponerer lambda i forhold til normal.
# Normal er lambda K x Kr, men i naesten alle funktioner laves den
# til en Kr x K matrix; uvist af hvilken grund

efficiencies <- function( object, ... )  {
    UseMethod( "efficiencies" )
}
eff <- function( object, ... )  {
    UseMethod( "efficiencies" )
}
eff.add <- function( object, ... )  {
    UseMethod( "efficiencies" )
}


# default method
efficiencies.default <- function( object, ... )  {
   return( object$eff )
}
eff.default <- function( object, ... )  {
   return( object$eff )
}



efficiencies.Farrell <- function(object, type="Farrell", ...)  {
# Returnerer efficencer som et array
   if ( type == "Farrell" )
     return(object$eff)
   else if ( type == "Shephard" )
     return(1/object$eff)
   else
     warning("Unknown type:", type)
   #   e <- as.matrix(object$eff)
   #   if ( object$ORIENTATION == "in" )  {
   #      colnames(e) <- "E"
   #   } else if ( object$ORIENTATION == "out" )  {
   #      colnames(e) <- "F"
   #   } else if ( object$ORIENTATION == "graph" )  {
   #      colnames(e) <- "G"
   #   }
   #   if ( !is.null(names(object$eff)) )  {
   #      rownames(e) <- names(object$eff)
   #   }
   #   if ( object$TRANSPOSE == TRUE ) {
   #      e <- t(e)
   #   }
   #   return(e)
} ## efficiency



eff.Farrell <- function(object, type="Farrell", ...)  {
      return( efficiencies.Farrell(object, type, ...) )
}





print.Farrell  <- function(x, digits=4, ...)  {
#   a <- cbind("Efficiens"=x$eff)
    a <- x$eff
#   a <- x@eff
   print(a, digits=digits, ...)
   invisible(a)
} ## print.Farrell



summary.Farrell <- function(object, digits=4, ...)  {
   eps <- 1e-7
   eff <- object$eff
   cat("Summary of ", ifelse(is.null(object$direct),"","directional "), 
       "efficiencies\n", sep="")
   cat(toupper(object$RTS)," technology and ",object$ORIENTATION,
       "put orientated efficiency\n", sep="")
   if ( sum(is.infinite(eff)) | sum(is.nan(eff)) )  {
    cat("Number of firms with infinite efficiency are ", 
        sum(is.infinite(eff)),"; removed below\n", sep="")
    eff <- eff[is.finite(eff)]
   }
   if ( is.null(object$direct) ) 
      cat("Number of firms with efficiency==1 are",
         sum(abs(eff-1) < eps, na.rm=TRUE), "out of", length(eff), 
         "\nMean efficiency:", format(mean(eff, na.rm=TRUE),digit=3), "\n---" )
   if ( object$ORIENTATION!="out" && is.null(object$direct) )  {
      # Orientering er in eller graph
      minE <- min(eff, na.rm=TRUE)
      minE <- floor( 10 * minE ) / 10
      dec <- seq(from=minE, to=1, by=.1)
      Estr <- "<= E <"
      Eeff <- "      E ==1   "
      n <- length(dec)
      estr <- rep(NA,n)
      for ( i in 1:(n-1) )
         estr[i] <- paste(dec[i],Estr,dec[i+1],"  ",sep="")
      estr[n-1] <- paste(estr[n-1]," ")
      estr[n] <- Eeff
      antal <- rep(NA,n)
      for ( i in 1:(n-1) )
         antal[i] <- sum(dec[i]-eps <= eff & eff < dec[i+1]-eps, na.rm=TRUE)
      antal[n] <- sum(abs(eff-1) < eps, na.rm=TRUE)
   } else if ( is.null(object$direct) )  {
      # Orientering er out
      maxF <- max(eff, na.rm=TRUE)
      maxF <- ceiling( 10 * maxF ) / 10
      dec <- seq(from=1, to=maxF, by=.1)
      Estr <- "< F =<"
      Eeff <- "F ==1   "
      n <- length(dec)
      if ( n > 10 )  {
         dec_  <- c(1,1.1,1.2,1.3,1.5,2.0,5.0,10.0,100.0,Inf)
         n <- length(dec_) 
         while ( n>1 && dec_[n-1] > maxF )
            n <- n - 1 
         dec <- dec_[1:n]
      }
      estr <- rep(NA,n)
      estr[1] <- paste("    ",Eeff)
      for ( i in 2:n )
         estr[i] <- paste(format(dec[i-1],digits=2,width=3),Estr,
                          format(dec[i],digits=3,width=3),"  ",sep="")
      his <- hist(object$eff, breaks=dec, plot=FALSE)
      antal <- his$counts
      # Foerste gruppe skal vaere eff==1; fra dist er foerte gruppe eff mellem 1 og 1.1
      antal[1] <- antal[1] - sum(abs(eff-1) < eps, na.rm=TRUE)
      antal <- c(sum(abs(eff-1) < eps, na.rm=TRUE), antal)
   } else {
      # directional er det saa
      cat("Number of firms with directional efficiency==0 are",
         sum(abs(eff) < eps, na.rm=TRUE), 
         "\nMean efficiency:", format(mean(eff, na.rm=TRUE),digit=3), "\n---" )
      his <- hist(object$eff, breaks=7, plot=FALSE)
      antal <- his$counts
      antal[1] <- antal[1] - sum( abs(eff) < eps , na.rm=TRUE)
      antal <- c(sum( abs(eff) < eps , na.rm=TRUE), antal)
      dec <- his$breaks
      Estr <- "< D =<"
      Eeff <- "D ==0   "
      estr <- rep(NA,length(his$counts)+1)
      estr[1] <- paste("    ",Eeff)
      for ( i in 1:length(his$counts) )  {
         estr[1+i] <- paste(format(dec[i],digits=2,width=3),Estr,
                          format(dec[i+1],digits=3,width=3),"  ",sep="")
      }
   }
   andel <- antal/sum(!is.na(eff))

   a <- cbind(antal , 100*andel)
   dimnames(a) <- list("  Eff range"=estr,c( "#", "%"))
   print(a,digits=c(2,3),quote=F,...)
   print(summary(eff))
    #   if ( SLACK & !is.null(object$slack) )  {
    #       sl <- object
    #       class(sl) <- "slack"
    #       summary(sl)
    #   }
   invisible(object)
}  ## summary.Farrell




# returns peers, i.e. numbers for units with positive lambda,
# efficient units to be compared to
peers <- function(object, NAMES=FALSE, N=1:dim(object$lambda)[1], LAMBDA=0)  {
   #  if ( object$TRANSPOSE ) {
   #    print("Colnames i lambda")
   #    print(colnames(object$lambda))
   #  } else {
   #    print("Rownames i lambda")
   #    print(rownames(object$lambda))
   #  }
   if ( class(object) != "Farrell" && class(object) != "slack" )
       stop("'Object' is not of class 'Farrell' (or 'slack');",
            " you might have used FAST=TRUE in 'dea'.")
   if ( object$TRANSPOSE ) {
        # Gem kun de raekker/firms i lambda, der skal laves peers for
        lam <- t(object$lambda[,N,drop=FALSE])
   } else {
        lam <- object$lambda[N,,drop=FALSE]   #'lam' er NxKr matrix
   }

   # Fjern foranstillet L_ eller L i soejlenavne for lambda
   if ( all("L_" == substr(colnames(lam),1,2)) )  {
      colnames(lam) <- substring(colnames(lam),3)
   }
   if ( all("L" == substr(colnames(lam),1,1)) )  {
      colnames(lam) <- substring(colnames(lam),2)
   }
   # 'lam' er N x Kr matrix

    # Liste, for hver firm et array af peers
    pt_ <- apply(lam, 1, function(x){which(x>LAMBDA)})
    if ( dim(lam)[1] == 1 )  {
        # Kun een firm; problem at 'pt_' bliver vektor, men
        # skal vaere liste derfor denne omskrivning
        pt_ <- list(c(pt_))
   }
   # Lav liste om til matrix og transponer den saa firms er raekker
   # Bemaerk pt_ er indeks er bench ogsaa indeks i reference matrix,
   # og ikke navne.
   # Herunder er 'sapply(pt_, length)' antal positive i hver raekke
   bench <- t(mapply(function(x) x[1:max(sapply(pt_, length))], pt_))
   # Hvis der kun er en peer for hver, bliver 'bench' et array,
   # en raekke-vektor, og ikke en soejle-vektor; derfor saettes
   # dim eksplicit.
    # Foerst finder vi det stoerste antal peers
    maxp <- max(sapply(pt_, length))
    # Hvis er ingen loesning er paa LP problemet er lambda's elementer
    # alle NA, og saa skal peers vaere NA, dvs. mindst een peer.
    if ( maxp==0 | is.na(maxp) ) maxp <- 1
    dim(bench) <- c(dim(lam)[1], maxp)
    rownames(bench) <- rownames(lam)
   
    # Skal der navne i matricen bench med peers i steder for blot indeks
    if ( is.logical(NAMES) && NAMES & 
        (!is.null(colnames(lam)) || !is.null(names(object$eff))) )  {
            # der skal navne, og enten er der names paa lambda eller paa eff.
        bench_ <- matrix(colnames(lam)[bench], nrow=dim(bench)[1])
        rownames(bench_) <- rownames(bench)
        bench <- bench_
    } else if (is(NAMES, "character") | is(NAMES, "integer") & 
                length(NAMES)==dim(bench)[1])  {
        # NAMES er et array med navne der bruges, dvs. character ell. 
        # integer array
        bench_ <- matrix(NAMES[bench], nrow=dim(bench)[1])
        rownames(bench_) <- rownames(bench)
        bench <- bench_
    }

    if ( object$TRANSPOSE ) {
        bench <- t(bench)
    }

    colnames(bench) <- paste("peer",1:dim(bench)[2],sep="")
    return(bench)
} ## peers



# For each unit return lambda-values for peers 
get.peers.lambda <- function(object, N=1:dim(object$lambda)[1], LAMBDA=0)  {
   lambda <- object$lambda
   if (object$TRANSPOSE) {
       lambda <- t(lambda)
   }
   if ( is.null(rownames(lambda)) ) {
    colnames(lambda) <- rownames(lambda) <- 1:dim(lambda)[1]
   }
   bench <- apply(lambda[N,,drop=FALSE], 1, function(x) {x[x>LAMBDA]})
   if (0)  {
        bench <- array(NA,dim=c(2,dim(lambda)[2],dim(lambda)[1]))
        maxj = 0  # Stoerste antal peers for en unit
        for ( i in 1:dim(lambda)[1] ) {  # For hver unit
            j = 0  # hvor mange peers for unit 'i'
            for ( h in 1:dim(lambda)[2] ) {  # For hver reference unit
                if ( lambda[i,h] > 0 ) {
                        j = j+1    # Har fundet en peer
                        bench[1,j,i] = h
                        bench[2,j,i] = lambda[i,h]
                }
                maxj = max(maxj,j)
            }
        }
        bench <- bench[,1:maxj,]
   }  ## if (0)
   return(bench)
} ## get.peers.lambda



print.peers  <- function(x, ...)  {
   a <- peers(x)
   print(a,...)
   invisible(a)
} ## print.peers



get.number.peers  <-  function(object, NAMES=FALSE, N=1:dim(object$lambda)[2], LAMBDA=0)  {
    # For hver peer get antal units den er peer for
   if ( object$TRANSPOSE ) {
       lam <- t(object$lambda)
   } else {
       lam <- object$lambda
   }
   # Kun de raekker af lambda som svarer til de oenskede peers
   lam <- lam[,N,drop=FALSE]

   # Fjern foranstillet L i soejlenavne for lambda
   # if ( "L" %in% substr(rownames(lam),1,1) )  {
   #    rownames(lam) <- substring(rownames(lam),2)
   # }
   # rownames er overfloedige da de fremgaar af foerste soejle, index er
   # ofte lettere for at kunne udpege raekker
   # rownames(lam) <- NULL
   if ( all("L_" == substr(colnames(lam),1,2)) )  {
       colnames(lam) <- substring(colnames(lam),3)
    }

    # Hvem er overhovedet peer for firms
   peer <- which(colSums(lam, na.rm=TRUE)>LAMBDA)
   # names(peer) <- NULL

    # Find hvor mange 'peer' er forbillede for
   count <- colSums(lam[,peer,drop=FALSE]>LAMBDA, na.rm=TRUE)
   np <- data.frame(peer,count)
   # if (NAMES) rownames(np) <- colnames(lam)[peer]
   if (NAMES) np$peer <- colnames(lam)[np$peer]
   return(np)
}  # get.number.peers



get.which.peers <- function(object, N=1:dim(object$lambda)[2], LAMBDA=0 )  {
   # Hvilke units en peer er peer for
   # Problem: Hvis N ikke er blandt peers gives en R fejl.
   # Jeg har ikke fundet en maade at teste for dette for at undgaa
   # R fejlen.
   if ( object$TRANSPOSE ) {
       lam <- t(object$lambda)
   } else {
      lam <- object$lambda
   }
   # Fjern foranstillet L_ eller L i soejlenavne for lambda
   if ( all("L_" == substr(colnames(lam),1,2)) )  {
        colnames(lam) <- substring(colnames(lam),3)
   }
   # Nu er lambda en K x Kr matrix
   # Er peer for en unit, naar units lambda er positivt,
   # positivt element i soejlen for N
   p <- apply(lam[,N,drop=FALSE] > LAMBDA, 2, which)
   p0 <- p[lapply(p,length) > 0]
   if ( length(p0) > 0 )
    return(p0)
   else
    return(NULL)
}  # get.which.peers




lambda.print  <- function(x, KEEPREF=FALSE, ...)  {
   if ( x$TRANSPOSE ) {
      lam <- t(x$lambda)
   } else {
      lam <- x$lambda
   }
   # print(class(lam))
   if (!KEEPREF && dim(lam)[1]>1 ) {
      lam <- lam[,rowSums(as.matrix(lam))>0]
   }
   xx <- round(unclass(lam)*1000)/1000
   if (any(ina <- is.na(lam))) 
      xx[ina] <- ""
   if ( any(i0 <- !ina & abs(lam) < 1e-5) ) 
      xx[i0] <- sub("0.0000", ".", xx[i0])
   if ( x$TRANSPOSE )
      xx <- t(as.matrix(xx))
   print(xx, quote=FALSE, rigth=TRUE, ...)
   invisible(xx)
   # printSpMatrix(Matrix(lam),digits=4, col.names=T,...)
   # invisible(lam)
} ## print.lambda



lambda <- function(object, KEEPREF=FALSE)  {
   if ( object$TRANSPOSE ) {
      lam <- object$lambda
   } else {
      lam <- t(object$lambda)
   }
   if (!KEEPREF && dim(lam)[2]>1 ) {
      lam <- lam[rowSums(lam, na.rm=TRUE)>0,,drop=FALSE]
   } else if (!KEEPREF && dim(lam)[2]==1 ) {
      lam <- lam[lam>0,,drop=FALSE]
   }

   if ( !object$TRANSPOSE )  lam <- t(lam)
   return(lam)
}



# Calculate excess input or output
excess <- function(object, X=NULL, Y=NULL)  {
   if ( class(object) != "Farrell" )
      stop("Only works for object of class/type 'Farrell'",
           " as output from dea and like functions")
   if ( is.null(object$direct) && is.null(X) && is.null(Y) )
      stop("Either X or Y is needed in the arguments for",
             " objects with no direction")

   e <- object$objval

   if ( is.null(object$direct) )  {
      # no direction, must be Farrell so direction is set by X or Y
      if ( object$ORIENTATION == "in" && !is.null(X) )
         ex <-  X * (1-e)
      else if ( object$ORIENTATION == "out" && !is.null(Y) )
         ex <-  Y * (e-1)
      else if ( object$ORIENTATION == "graph"&& !is.null(X)&& !is.null(Y) )
         ex <- cbind((1-e)*X, (1/e-1)*Y )
      else # ( is.null(dir) )
         stop("X/Y missing for ORIENTATION = ", object$ORIENTATION )
   } else {
       if ( is(object$direct, "matrix") )  {
          ex <- apply(object$direct,2,"*",e)
       } else {
           dir <- matrix(object$direct, nrow=length(e), 
                           ncol=length(object$direct), byrow=TRUE )
           ex <- e * dir
      }
   }
   # Afrund til 0 hvis ex er naer 0
   eps <- sqrt(.Machine$double.eps)
   ex[abs(ex)<eps] <- 0
   return(ex)
} # excess




eladder <- function(n, X, Y, RTS="vrs", ORIENTATION="in", 
                    XREF=NULL, YREF=NULL, DIRECT=NULL, param=NULL, MAXELAD=NULL)  {

   if ( is.null(XREF) )  {
      XREF <- X
   }
   if ( is.null(YREF) )  {
      YREF <- Y
   }

   idx <- NULL
   if ( missing(MAXELAD) || is.null(MAXELAD) ) {
    MAXELAD <- dim(XREF)[1]
   } else {
    if( !is.numeric(MAXELAD) ) stop("MAXELAD must be an integer")
        MAXELAD <- min(abs(MAXELAD), dim(XREF)[1])
   }
   elad <- rep(NA, MAXELAD)
   for ( i in 1:MAXELAD )  {
    # if (LP) print(paste("===> ",i),quote=FALSE)
    # if (LP) print("FRONT.IDX")
    # if (LP) print(idx)
    # print( length(idx) == dim(X)[1] )  
    if ( length(idx) == MAXELAD )  {
       break  
    }
    # Brug FRONT.IDX for at kunne bruge de oprindelige indeks i X og Y
    e <- dea(X[n,,drop=F],Y[n,,drop=F], RTS=RTS, ORIENTATION=ORIENTATION,
             XREF=XREF, YREF=YREF, FRONT.IDX=idx, DIRECT=DIRECT, 
             param=param, LP=FALSE)
    # if (LP) print(paste("Eff =",eff(e)), quote=FALSE)
    # if (LP) print(paste("Peers =",peers(e)), quote=FALSE)
    # if (LP) print(lambda(e))
    # Er der nogen peers overhovedet ellers kan vi bare slutte nu
    if ( abs(eff(e)) == Inf ) break
    if ( is.na(peers(e)[1]) ) break
    elad[i] <- e$eff
    # Array nr. for den stoerste vaerdi af lambda
    # Bruger kun den foerste hvis der er flere
    p <- which.max(e$lambda)
    # if (LP) print(paste("p =",p))
    # firm number for array number p, firm numbers follow L in the colnames
    str <- substring(colnames(e$lambda)[p],2)
    # if (LP) print(str)
    suppressWarnings(ip <- as.integer(str))
    # if (LP) print(paste("    ", ip),quote=FALSE)
    # nok en firm der ikke laengere skal indgaa i referenceteknologien
    if ( is.na(ip) )  {
       # det er et navn/en streng saa den skal laves om til et indeks,
       # et heltal, for elleres kan den ikke bruges som indeks til FRONT.IDX
       str0 <- substring(str,2)
       # if (LP) print(str0)
       navne <- rownames(XREF)
       # if (LP) print(navne)
       if ( is.null(navne) )
           navne <- rownames(YREF) 
       # Find indeks for placering af str0 i rownames(XREF)/rownames(YREF)
       ip <- which( navne %in% str0 )
    }
    # saa er ip et tal
    idx <- c(idx,-ip)
  }
  elad <- na.omit(elad)
  elad <- elad[1:length(elad)]
  if ( is.null(idx) )  {
     idx <- NA
  } else {
     idx <- -idx
  }
  return(list(eff=elad, peer=idx, lastp=peers(e)))
}  ## eladder



eladder.plot <- function(elad, peer, TRIM=NULL, ...)  {
   if ( all(is.na(elad)) )
      stop("All values of first argument are NA")
   if ( !is.null(TRIM) & !is.numeric(TRIM) )
       stop("TRIM must be an integer")
   if ( is.null(TRIM) )  {
      TRIM <- 0
      for ( i in 1:length(peer) )  {
         TRIM <- max(TRIM, nchar(toString(peer[i])))
      }
   }
   linje <- ifelse(TRIM==1,2,TRIM^(1/1.3))
   opar <- 
       par(mar=c(linje+2,4.1,4.1,2.1))
   plot(elad, xaxt="n", xlab="", ylab="Efficiency", ...)
   mtext("Most influential peers", side=1, line=linje+.5)
   if ( class(peer) == "character" || class(peer) == "factor" )  {
      axis(1, at=1:length(peer),
            labels=strtrim(peer,TRIM), las=ifelse(TRIM>1,2,0) )
   } else {
      axis(1, at=1:length(peer), labels=peer, las=ifelse(TRIM>1,2,0) )
   }
   abline(v=which(elad==1), lty=3)
   abline(h=1, lty=3)
   par(opar)
}  ## eladder.plot


# Funktion til at droppe en eller flere units fra et
# Farrrell objekt
#dropUnit <- function(E, dmu)  {
#   # Det foerste element er typisk eff og kan give
#   # antal units.
#   K <- length(E[[1]])
#   for (n in 1:length(names(E)))  {
#       if ( class(E[[n]]) == "matrix" ) 
#           E[[n]] <- E[[n]][-dmu,,drop=FALSE]
#       else if ( class(E[[n]]) == "numeric" )
#           E[[n]] <- E[[n]][-dmu]
#   }
#   return(E)
#}
