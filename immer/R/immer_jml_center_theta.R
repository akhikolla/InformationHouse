## File Name: immer_jml_center_theta.R
## File Version: 0.04

immer_jml_center_theta <- function( theta, center_theta)
{
    if (center_theta){
        theta <- theta - mean(theta)
    }
    return(theta)
}
