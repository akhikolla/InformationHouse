## File Name: plot.immer_hrm.R
## File Version: 0.07

plot.immer_hrm <- function( x, ... ){
    class(x) <- "mcmc.sirt"
    graphics::plot( x, ... )
}

