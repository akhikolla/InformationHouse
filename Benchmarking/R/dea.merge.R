# $Id: dea.merge.R 229 2020-07-04 13:39:18Z lao $


dea.merge <- function(X, Y, M, RTS = "vrs", ORIENTATION = "in", 
    XREF = NULL, YREF = NULL, FRONT.IDX = NULL, TRANSPOSE = FALSE, 
    CONTROL = NULL)
{

   rts <- c("fdh","vrs","drs","crs","irs","irs","add")
   if ( is.numeric(RTS) )  {
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) )  stop("Unknown scale of returns: ", RTS)

   orientation <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   if ( !(ORIENTATION %in% orientation) ) {
      stop("Unknown value for ORIENTATION: ", ORIENTATION)
   }


# Mergers should be measured against the originale technology set,
# therefore we must use XREF og YREF otherwise it will be Xmerger and
# Ymerger that determines the technology.

.xyref.missing <- FALSE
if ( missing(XREF) || is.null(XREF) )  {
   .xyref.missing <- TRUE
   XREF <- X
}
if ( missing(YREF) || is.null(YREF) )  {
   .xyref.missing <- TRUE && .xyref.missing
   YREF <- Y
}

Xmerger <- M %*% X
Ymerger <- M %*% Y

# Potential gains
E <- dea(Xmerger, Ymerger, RTS=RTS, ORIENTATION=ORIENTATION, 
    XREF = XREF, YREF = YREF, CONTROL = CONTROL,
    FRONT.IDX = FRONT.IDX, TRANSPOSE = TRANSPOSE, FAST = TRUE)


# Individual efficiencies
e <- dea(X, Y, RTS=RTS, ORIENTATION=ORIENTATION, 
      XREF = XREF, YREF = YREF, FRONT.IDX = FRONT.IDX, 
      TRANSPOSE = TRANSPOSE, FAST = TRUE, CONTROL = CONTROL)


# Inputs of individual firms projected on efficient frontier and
# inputs of merged firms after elimination of individual inefficiency
if ( ORIENTATION == "in" )  {
   Xeff <- diag(e) %*% X
   Xmerger_proj <- M %*% Xeff
   Ymerger_proj <- Ymerger
} else if ( ORIENTATION == "out" )  {
   Yeff <- diag(e) %*% Y
   Xmerger_proj <- Xmerger
   Ymerger_proj <- M %*% Yeff
} else if ( ORIENTATION == "graph" )  {
   Xeff <- diag(e) %*% X
   Yeff <- diag(e) %*% Y
   Xmerger_proj <- M %*% Xeff
   Ymerger_proj <- M %*% Yeff
} else  {
   stop("Unknown ORIENTATION: ", ORIENTATION, ". The function dea.merge stops.")
}

# Pure gains from mergers
Estar <- dea(Xmerger_proj, Ymerger_proj, RTS=RTS, ORIENTATION=ORIENTATION, 
      XREF = XREF, YREF = YREF, FRONT.IDX = FRONT.IDX, 
      TRANSPOSE = TRANSPOSE, FAST = TRUE, CONTROL = CONTROL)


# Learning effect
LE <- E/Estar 

# Inputs and outputs for merged firms in harmony calculation
Xharm <- diag(1/rowSums(M), nrow=dim(M)[1]) %*% Xmerger_proj
Yharm <- diag(1/rowSums(M), nrow=dim(M)[1]) %*% Ymerger_proj

# Harmony effect
HA <- dea(Xharm,Yharm, RTS=RTS, ORIENTATION=ORIENTATION, 
      XREF = XREF, YREF = YREF, FRONT.IDX = FRONT.IDX, 
      TRANSPOSE = TRANSPOSE, FAST = TRUE, CONTROL = CONTROL) 


# Size effect
SI <- Estar/HA

em <- list(Eff=E, Estar=Estar, learning=LE, harmony=HA, size=SI)

return(em)

} 
