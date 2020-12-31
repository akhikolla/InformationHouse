## File Name: harm_mean.R
## File Version: 0.01


#***********************************************
# compute harmonic mean
harm_mean <- function(x)
{
    exp( mean( log(x) ) )
}
#***********************************************
