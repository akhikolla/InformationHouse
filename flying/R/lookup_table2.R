# load table 2 generate. Description in Pennycuick earlier version
# @name .gen.table2
# @author Brian Masinde
# @return table2
# @description Pennycuick's table II aids in calculation of D factor for finding
#              lift:drag ratio. C factor for finding power required at maximum
#              range speed and B a ratio of maximum range speed and minimum
#              power speed.
#


.gen.table2 <- function() {
  x1Plusx2 <- c(0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00,
                2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25,
                4.50, 4.75, 5.00)
  B <- c(1.360, 1.386, 1.452, 1.515, 1.574, 1.631, 1.684, 1.735, 1.784, 1.830,
         1.875, 1.918, 1.959, 1.999, 2.038, 2.075, 2.111, 2.146, 2.180, 2.213,
         2.246)
  C <- c(1.140, 1.458, 1.783, 2.115, 2.453, 2.795, 3.141, 3.490, 3.841, 4.195,
         4.550, 4.907, 5.266, 5.625, 5.986, 6.348, 6.711, 7.074, 7.438,
         7.803, 8.168)
  D <- c(1.000, 0.824, 0.706, 0.621, 0.556, 0.506, 0.465, 0.431, 0.402, 0.378,
         0.357, 0.339, 0.322, 0.308, 0.295, 0.283, 0.273, 0.263, 0.254, 0.246,
         0.238)
  table2 <- as.data.frame(cbind(x1Plusx2, B, C, D))

  return(table2)
}

#### Inerpolate D values in table 2 ####
#' @author Brian Masinde
#' @name .interpolate
#' @param x1plusx2 sum of metabolic power profile power ratio
#' @param table2 interpolation table II for B, C D factors based on sum of
#'               metabolic power ratio and profile power ratio
#' @return D factor


.interpolate <- function(x1plusx2, table2) {
  if (length(which(table2$x1Plusx2 == x1plusx2)) == 0) {
    upId <- which(table2$x1Plusx2 > x1plusx2)[1]
    lowId <- tail(which(table2$x1Plusx2 < x1plusx2), 1)

    if ((abs(x1plusx2 - table2$x1Plusx2[upId])) >=
        (abs(x1plusx2 - table2$x1Plusx2[lowId]))) {
      D <- table2$D[lowId]
    } else {
      D <- table2$D[upId]
    }
  } else {
    D <- table2$D[which(table2$x1Plusx2 == x1plusx2)]
  }

  return(D)
}

