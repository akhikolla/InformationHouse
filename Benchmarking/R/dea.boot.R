# $Id: dea.boot.R 229 2020-07-04 13:39:18Z lao $

# Boot No FEAR: boot.nf
# Bootstrap af dea model a la Simar Wilson 1998.



dea.boot <- function(X,Y, NREP=200, EFF=NULL, RTS="vrs", 
       ORIENTATION="in", alpha=0.05, XREF=NULL, YREF=NULL, FRONT.IDX=NULL,
       EREF=NULL, DIRECT=NULL, TRANSPOSE=FALSE, SHEPHARD.INPUT=TRUE, LP=FALSE, 
       CONTROL=NULL)
{

   if ( !is.null(DIRECT) )
      stop("Use of DIRECT does not yet work in dea.boot")

   rts_ <- c("fdh","vrs","drs","crs","irs","irs","add")
   if ( is.numeric(RTS) )  {
      RTStemp <- rts_[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
   }
   RTS <- tolower(RTS)
   rts <- which(RTS == rts_) -1 # Foerste i RTS 'fdh' er 0
 
   if ( rts < 1 || rts > 4 )
      stop("Invalid value of RTS in call to boot:",RTS)


   orientation <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   orientation <- which(ORIENTATION == orientation) -1
   if ( orientation < 1 || 3 < orientation )
      stop("Invalid value of ORIENTATION in call to boot")

   if (LP) print(paste("rts =",RTS,"   orientation=",ORIENTATION), quote=FALSE)
   if (LP) print(paste("rts =",rts,"     orientation=",orientation), quote=FALSE)

   if ( !is(X, "matrix") )
      stop("X is not a matrix")
   if ( !is(Y, "matrix") )
      stop("Y is not a matrix")
   if ( !is.null(XREF) && !is(XREF, "matrix") )
      stop("XREF is not a matrix")
   if ( !is.null(YREF) && !is(YREF, "matrix") )
      stop("YREF is not a matrix")

   if ( !is.null(EREF) && (is.null(XREF) || is.null(YREF)) )
      stop("When EREF is present then XREF and YREF must also be present") 

   .xyref.missing <- FALSE
   if ( missing(XREF) || is.null(XREF) )  {
      .xyref.missing <- TRUE
      XREF <- X
   }
   if ( missing(YREF) || is.null(YREF) )  {
      .xyref.missing <- TRUE && .xyref.missing
      YREF <- Y
   }
   
   if ( TRANSPOSE )  {
      X <- t(X)
      Y <- t(Y)
      XREF <- t(XREF)
      YREF <- t(YREF)
      if ( !is.null(DIRECT) & is(DIRECT, "matrix") )
         DIRECT <- t(DIRECT)
   }

   if ( length(FRONT.IDX) > 0 )  {
      if ( !is.vector(FRONT.IDX)  & !(ifelse(is.matrix(FRONT.IDX), dim(FRONT.IDX)[2]==1, FALSE) ))
         stop("FRONT.IDX is not a vector or collumn matrix in 'dea'")
      XREF <- XREF[FRONT.IDX,, drop=FALSE]
      YREF <- YREF[FRONT.IDX,, drop=FALSE]
   }

   rNames <- rownames(XREF)
   if ( is.null(rNames) & !is.null(rownames(YREF)) )
      rNames <- rownames(YREF)

   K <- dim(X)[1]  # number of units, firms, DMUs
   Kr <- dim(XREF)[1] # number of units,units in the reference technology
 
   if ( is.null(EFF) )  
      eff <- dea(X,Y, RTS, ORIENTATION, XREF, YREF, FAST=TRUE, CONTROL=CONTROL)
   else if ( is(EFF, "Farrell") )
      eff <- EFF$eff
   else 
      eff <- EFF

   if ( length(eff) != K ) 
      stop("The length of EFF must be the number of firms in X and Y")

   if ( is.null(EREF) ) {
    if ( identical(X, XREF) & identical(Y, YREF) )  {
        EREF <- eff
    } else {
        EREF <- dea(XREF, YREF, RTS, ORIENTATION, FAST=TRUE, CONTROL=CONTROL)
      }
    }


    # Hvis eff eller EREF har NA elementer vil X eller XREF faa NA 
    # elemnter og dermed kommer der NA i elementer i LP problemet og
    # det giver fejl i 'set.row' ved kald af dea. Derfor tjek af om 
    # der er NA.
    if (!all(!is.na(eff))) {
        nr <- which(is.na(eff))
        msg <- paste0("The dea calculation has NA values for firms ", 
            paste(nr, collapse=", "), ".\n",
            "For bootstrap either remove these firms from the data, use a (another) scaling\n",
            "of the data, or use a scaling option using CONTROL, cf. the manual for 'dea'.\n")
        stop(msg)
    }
    if (!all(!is.na(EREF))) {
        nr <- which(is.na(EREF))
        msg <- paste0("The dea calculation has NA values for firms ", 
            paste(nr, collapse=", "), " in the reference technology.\n",
            "For bootstrap either remove these firms from the reference data, use a (another) scaling\n",
            "of the data, or use a scaling option using CONTROL, cf. the manual for 'dea'.\n")
        stop(msg)
    }

    # Bootstrap noget der er over 1, dvs. Farrel output eller Shephard input
   if (ORIENTATION=="out")  farrell<-TRUE else farrell<-FALSE


   # Beregn vinduesbredden til udglatning af efficiensers fordeling.
   # Lav det alene for efficencer paa 1 eller over, dvs. for Farrell
   # output eller Shephard input orienterede efficiencer fordi det er
   # dem der samples paa.
   # Fjern 1-erne saa de ikke dominerer og spejl efficiencer i 1,
   # Daraio and Simar (2007, 61) ligning (3.23) og (3.26)
   dist <- eff
   if ( !farrell )  dist <- 1/dist
   # Behold efficiencer over 1, drop 1-erne naar bredden beregnes
   zeff <- eff[ dist > 1 + 1e-6 ]
   if ( length(zeff) == 0 )  {
      cat("No unit with efficiency different from 1.0000.\n", quote=FALSE)
      stop("The range of efficiencies is degenrate, 'dea.boot' stops.")
   }
   # Spejling om 1
   neff <- c(zeff,2-zeff)
   # Ratio |adjust| som bredden skal ganges med fordi der er for mange
   # elementer i den spejlede vektor som baandbredden bregnes ud fra
   adjust <- sd(dist)/sd(neff) * (length(neff)/length(dist))^(1/5)
   # Ligning (3.28), (3.30) og (3.31) i Silverman (1986, 45--47)
   std <- sd(neff) 
   iqr <- IQR(neff)/1.349
   # Hvis bredden er for lille kan IQR vaere 0 og saa skal std bruges
   # Hvorfor '0.9' og ikke '1.06' som staar i Daraio and  Simar (2007, 61)
   # ligning (3.23)?  Fordi Silverman side 48 anbefaler 0.9.
   if ( iqr > 1e-6 && iqr < std  ) std <- iqr
   h0 <- .9 * std * length(neff)^(-1/5)
   # h0 <- 1.06 * std * length(neff)^(-1/5)
   h <- adjust * h0
   if (LP) cat("Bandwidth =",h,"\n")

    # matrix til at gemme de bootstrappede efficiencer
    boot <- matrix(NA, nrow=K, ncol=NREP)


if ( ORIENTATION == "in" )  {
   for ( b in 1:NREP )  {
      # if (LP) print(paste(b," -----"))
      estar <- dea.sample(eff, h, Kr)
      # estar/eff ganges paa hver soejle i X; benytter at R har data
      # efter soejler, byrow=FALSE er default i R.
      # eff*XREF er paa randen, og eff/estar*XREF bliver indre punkt
      xstar <- EREF/estar * XREF
      #### xstar <- eff/estar * XREF
      boot[,b] <- dea(X,Y, RTS, ORIENTATION, XREF=xstar, YREF, FAST=TRUE, CONTROL=CONTROL)
   }  # b in 1:NREP
} else if ( ORIENTATION == "out" )  {   # farrel er TRUE
   for ( b in 1:NREP )  {
      # Farrell output efficiencer
      estar <- dea.sample(eff, h, Kr)
      ystar <- EREF/estar * YREF
      boot[,b] <- dea(X,Y, RTS, ORIENTATION, XREF, YREF=ystar, FAST=TRUE, CONTROL=CONTROL)
   }  # b in 1:NREP
} else if ( ORIENTATION == "graph" )  {
   for ( b in 1:NREP )  {
      estar <- dea.sample(eff, h, Kr)
      xstar <- EREF/estar * XREF
      ystar <- estar/EREF * YREF
      boot[,b] <- dea(X,Y, RTS, ORIENTATION, XREF=xstar, YREF=ystar,
                      FAST=TRUE, CONTROL=CONTROL)
   }  # b in 1:NREP
} else {
   stop("Unknown ORIENTATION for dea.boot")
}


    if ( SHEPHARD.INPUT & ORIENTATION != "out" )  {
        # Lav input efficiencer storre end 1 hvis den er under 1
        eff <- 1/eff
        boot <- 1/boot
    }
    # Efficiencer stoerre end 1
    bias <- rowMeans(boot, na.rm=TRUE) - eff
    eff.bc <- eff - bias

    #ci <- t(apply(boot,1, quantile, 
    #           probs=c(0.5*alpha, 1-0.5*alpha), type=9, na.rm=TRUE))
    # ci_ <- t(apply(eff - boot, 1, quantile, 
    #            probs=c(0.5*alpha, 1-0.5*alpha), na.rm=TRUE, type=9))
    ci_ <- t(apply(eff-boot, 1, function(x) {
        quantile(x, probs=c(0.5*alpha, 1-0.5*alpha), type=9, na.rm=TRUE)}))
    ci <- eff + ci_

    if (SHEPHARD.INPUT & ORIENTATION!="out")  {
        # Lav efficiencer mm. tilbage til under 1
        eff <- 1/eff
        boot <- 1/boot
        eff.bc <- 1/eff.bc
        ci <- t(apply(1/ci, 1, rev))  # Byt om paa de lave og hoeje graenser
        bias <- eff - eff.bc
    }
    var <- apply(boot,1, function(x) {var(x, na.rm=TRUE)})

   bb <- list(eff=eff, eff.bc=eff.bc,
           bias=bias,
           var=var,
           conf.int=ci,
           eref=EREF, boot=boot)

   return( bb )

}  # boot




# Bootstrap sample fra efficiencer over 1, 
# Shephard input eller Farrell output
dea.sample <- function(e, h, K=NULL)  {
   # sample for Shaphard input afstandsfunktion, dvs. vaerdier over 1
   if ( min(e) < 1-1e-6 )  {
      e <- 1/e   
      farrell <- TRUE
   } else {
      farrell <- FALSE  # Vaerdi over 1 er ikke Farrel, men Shephard
   }
   if ( is.null(K) )
      K <- length(e)
   # spejling omkring 1
   espejl <- c(e, 2-e)
   beta <- sample(espejl, K, replace=TRUE)
   etilde <- beta + h * rnorm(K)
   # Juster saa varians bliver rigtig, Daraio and Simar (2007, 62) [3]
   estar <- mean(beta) + (etilde-mean(beta)) / sqrt(1+h^2/var(espejl))
   # Genspejl for at faa alle vaerdier storre end 1
   estar <- ifelse(estar < 1, 2-estar, estar)
   if ( farrell )  {
      # Omsaet til Farrell
      estar <- 1/estar
   }
   return(estar)
}


