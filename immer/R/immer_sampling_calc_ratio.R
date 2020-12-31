## File Name: immer_sampling_calc_ratio.R
## File Version: 0.06

immer_sampling_calc_ratio <- function(ll_old, ll_new, p_old, p_new, eps=1E-20)
{
    max_ll_diff <- 200
    ll_diff <- ll_new - ll_old
    ll_diff <- ifelse( ll_diff > max_ll_diff, max_ll_diff, ll_diff)
    ratio <- p_new * exp( ll_diff ) / ( p_old  + eps )
    ratio[ is.na(ratio) ] <- 0
    return(ratio)
}
