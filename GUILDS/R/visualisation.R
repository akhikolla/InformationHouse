octave_index <- function(ab_in) {
  min_ab <- 1
  max_ab <- 2
  result <- 1
  while (!((ab_in < max_ab) && (ab_in >= min_ab))) {
    min_ab <- min_ab * 2
    max_ab <- max_ab * 2
    result <- result + 1
  }
  return(result)
}


preston_sort <- function(abund) {
  output <- rep(0, pracma::ceil(log(max(abund)) / log(2)))
  for (i in seq_along(abund)) {
    index <- octave_index(abund[i])
    output[index] <- output[index] + 1
  }
  return(output)
}

preston_plot <- function(abund, expected, ...) {
  v <- preston_sort(abund)
  maxy <- max(v)
  plot_expected <- FALSE
  if(!missing(expected)) {
    plot_expected <- TRUE
  }
  
  if(plot_expected) {
    maxy <- max(maxy,max(expected))
  }
  
  df.bar <- graphics::barplot(v, names.arg = 0:(length(v)-1), 
                                      cex.names = 0.7, 
                                      space = 0, 
                                      ylim = c(0,maxy), ...)
  if(plot_expected) {
    xx <- c()
    yy <- seq_along(df.bar)
    lmmm <- stats::lm(df.bar~yy)
    slope <- lmmm$coefficients[[2]]
    intercept <- lmmm$coefficients[[1]]
    xvals <- seq_along(expected)
    xvals <- xvals * slope + intercept
    graphics::lines(expected~xvals, lwd = 2)
  }
}






