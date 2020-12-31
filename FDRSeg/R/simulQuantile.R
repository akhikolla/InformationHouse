# simulate quantiles
simulQuantile <- function(alpha, n, r=round(50/min(alpha, 1-alpha)), type=c("smuce","fdrseg","dfdrseg"), convKern, pos = .GlobalEnv)
{
  type <- match.arg(type)
  if (type == "smuce") {                # SMUCE
    ret <- .scalar_quantile(alpha, n, r)
  } else if (type == "fdrseg") {        # FDRSeg
    ret <- .vector_quantile(alpha, n, r)
  } else { # if (type == "dfdrseg")     # D-FDRSeg
    if (missing(convKern))
      stop("Simulate D-FDRSeg quatiles: 'convKern' has to be specified!")
    simData <- NULL    
    ret <- .dependent_vector_quantile(alpha, n, r, convKern, pos)
  }
  ret
}

# local function hidden from users
.dependent_vector_quantile <- function(alpha, n, r, convKern, pos) 
{
  isdone = 0
  lenK   = length(convKern)
  if (exists(".simulData", where = pos, inherits = TRUE)) {
    simData = get(".simulData", pos = pos, inherits = TRUE)
    n0      = simData$n
    r0      = simData$r
    lenK0   = length(simData$convKern)
    if (n0 >= n && r0 >= r && lenK0 == lenK) {
      if (all.equal(convKern, simData$convKern)) {
        isdone = 1
      }
    }
  } 
  if (isdone == 1) {
    data = simData$data[1:n, ]
  } else {
    data  <- matrix(nrow=n, ncol=r)
    fKern <- fft(c(convKern, rep(0,n-1))) / (n+lenK-1)
    sd    <- 1 / sqrt(sum(convKern^2))
    for (i in 1:r) {
      Y <- rnorm(n+lenK-1, 0, sd) 
      Y <- Re(fft(fft(Y) * fKern, inverse=TRUE))[lenK:(lenK+n-1)] 
      # compute multiresolution statistics via c code
      data[,i] <- .mrstatvec_cpp(Y)
    }
    simData <- list(n=n, r=r, convKern=convKern, data=data)
    assign(".simulData", simData, pos = pos, inherits = TRUE)
  }
  qnt <- numeric(n)
  for (i in 1:n) {
    qnt[i] <- quantile(data[i,], alpha)
  }
  qnt 
} 

# # local function hidden from users
# .dependent_vector_quantile <- function(alpha, n, r, convKern) 
# {
#   isdone = 0
#   lenK   = length(convKern)
#   if (exists("simData")) {
#     n0    = simData$n
#     r0    = simData$r
#     lenK0 = length(simData$convKern)
#     if (n0 >= n && r0 >= r && lenK0 == lenK) {
#       if (all.equal(convKern, simData$convKern)) {
#         isdone = 1
#       }
#     }
#   } 
#   if (isdone == 1) {
#     data = simData$data[1:n, ]
#   } else {
#     data  <- matrix(nrow=n, ncol=r)
#     fKern <- fft(c(convKern, rep(0,n-1))) / (n+lenK-1)
#     sd    <- 1 / sqrt(sum(convKern^2))
#     for (i in 1:r) {
#       Y <- rnorm(n+lenK-1, 0, sd) 
#       Y <- Re(fft(fft(Y) * fKern, inverse=TRUE))[lenK:(lenK+n-1)] 
#       # compute multiresolution statistics via c code
#       data[,i] <- .mrstatvec_cpp(Y)
#     }
#     simData <<- list(n=n, r=r, convKern=convKern, data=data)
#   }
#   qnt <- numeric(n)
#   for (i in 1:n) {
#     qnt[i] <- quantile(data[i,], alpha)
#   }
#   qnt 
# } 
