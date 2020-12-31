## File Name: immer_trim_increment.R
## File Version: 0.04

immer_trim_increment <- function(incr, max_incr)
{
    do_cut <- TRUE
    while( do_cut ){
        incr <- ifelse( abs(incr) > max_incr, incr / 2, incr )
        do_cut <- ( sum( abs(incr) > max_incr ) > 0 )
    }
    return(incr)
}
