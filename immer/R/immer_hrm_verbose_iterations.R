## File Name: immer_hrm_verbose_iterations.R
## File Version: 0.05

immer_hrm_verbose_iterations <- function(iter, it, print_iter, mcmc_start_time)
{
    if ( it %% print_iter==0 ){
        s_temp <- Sys.time()
        t0 <- difftime( s_temp, mcmc_start_time )
        time_diff <- ( iter - it ) * t0 / it
        units <- "mins"
        time_diff <- as.numeric( time_diff, units=units )
        time_diff <- round( time_diff, digits=2 )
        t0 <- round( as.numeric( t0, units=units ), digits=2 )
        p1 <- paste0("*** Iteration ", it, " | ",
                     "Remaining: ", time_diff, " ", units,
                    " | Elapsed: ", t0, " ", units,
                    " \n" )
        cat(p1)
        utils::flush.console()
    }
}
