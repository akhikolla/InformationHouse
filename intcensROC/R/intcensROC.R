# function to compute the knots
.knotT <- function(U, V, delta, dim) {
    size <- dim - 2
    knot_pre_t = c(U[delta != 1], U[delta != 2], V[delta != 2], V[delta != 3])
    qt = rep(0, size + 1)
    qt[1] = 0
    for (i in 2:(size + 1)) qt[i] = quantile(knot_pre_t, (i - 1)/size, name = F, 
        na.rm = TRUE)
    knots = c(qt[1], qt[1], qt, qt[size + 1], qt[size + 1])
}

.knotM <- function(marker, dim) {
    size <- dim - 2
    knot_pre_m = marker
    qt = rep(0, size + 1)
    qt[1] = 0
    for (i in 2:(size + 1)) qt[i] = quantile(knot_pre_m, (i - 1)/size, name = F, 
        na.rm = TRUE)
    qt[size + 1] = max(marker + 0.1)
    knots = c(qt[1], qt[1], qt, qt[size + 1], qt[size + 1])
}

intcensROC <- function(U, V, Marker, Delta, PredictTime, gridNumber = 500) {
    
    if (any(U < 0)) 
        stop(paste0("Negative U found!"))
    if (any(V - U < 0)) 
        stop(paste0("V - U < 0 found!"))
    if (any(Marker < 0)) 
        stop(paste0("Negative marker value found!"))
    if (PredictTime < min(U) || PredictTime > max(V)) 
        stop(paste0("Predict time out of range min(U) and max(v)!"))
    
    # detemine the dimension of spline function
    size <- length(U)
    cadSize <- size^(1/3)
    if (cadSize - floor(cadSize) < ceiling(cadSize) - cadSize) {
        Dim <- floor(cadSize) + 2
    } else {
        Dim <- ceiling(cadSize) + 2
    }
    
    # compute the knots for time and marker
    knotT <- .knotT(U, V, Delta, Dim)
    knotM <- .knotM(Marker, Dim)
    
    # compute the thetas and produce the AUC curve
    theta = .Call("_intcensROC_sieve", PACKAGE = "intcensROC", U, V, Marker, Delta, 
        knotT, knotM, Dim)
    max_m = max(Marker)
    TPt = NULL
    FPt = NULL
    markerGrid = max_m * seq(1:gridNumber)/gridNumber
    for (i in 1:gridNumber) {
        TPt[i] = .Call("_intcensROC_truePos", PACKAGE = "intcensROC", theta, markerGrid[i], 
            max_m, PredictTime, knotT, knotM)
        FPt[i] = .Call("_intcensROC_falsePos", PACKAGE = "intcensROC", theta, markerGrid[i], 
            max_m, PredictTime, knotT, knotM)
    }
    
    # check if marker value and survivual time are positve related.
    ROCdata <- data.frame(fp = FPt, tp = TPt)
    FP = ROCdata$fp
    TP = ROCdata$tp
    
    tmpAUC <- pracma::trapz(sort(FP, decreasing = FALSE), sort(TP, decreasing = FALSE))
    if (tmpAUC < 0.5) {
        TPt = NULL
        FPt = NULL
        Marker = max(Marker) + min(Marker) - Marker
        knotM <- .knotM(Marker, Dim)
        
        theta = .Call("_intcensROC_sieve", PACKAGE = "intcensROC", U, V, Marker, 
            Delta, knotT, knotM, Dim)
        max_m = max(Marker)
        TPt = NULL
        FPt = NULL
        markerGrid = max_m * seq(1:gridNumber)/gridNumber
        for (i in 1:gridNumber) {
            TPt[i] = .Call("_intcensROC_truePos", PACKAGE = "intcensROC", theta, 
                markerGrid[i], max_m, PredictTime, knotT, knotM)
            FPt[i] = .Call("_intcensROC_falsePos", PACKAGE = "intcensROC", theta, 
                markerGrid[i], max_m, PredictTime, knotT, knotM)
        }
        ROCdata <- data.frame(fp = FPt, tp = TPt)
    }
    ROCdata
}

