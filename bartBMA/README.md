# bartBMA

<!-- badges: start -->
<!-- badges: end -->

The goal of bartBMA is to provide an implementation of Bayesian Additive Regression Trees Using Bayesian Model Averaging (BART-BMA) (Hernandez et al. 2018)

Hernández, B., Raftery, A. E., Pennington, S. R., & Parnell, A. C. (2018). Bayesian additive regression trees using Bayesian model averaging. Statistics and computing, 28(4), 869-890.


## Installation

``` r
library(devtools)
install_github("bartBMA")
```

## Example


``` r
library(bartBMA)
## basic example code

N <- 100
p<- 100
set.seed(100)

epsilon <- rnorm(N)
xcov <- matrix(runif(N*p), nrow=N)
y <- sin(pi*xcov[,1]*xcov[,2]) + 20*(xcov[,3]-0.5)^2+10*xcov[,4]+5*xcov[,5]+epsilon

epsilontest <- rnorm(N)
xcovtest <- matrix(runif(N*p), nrow=N)
ytest <- sin(pi*xcovtest[,1]*xcovtest[,2]) + 20*(xcovtest[,3]-0.5)^2+10*xcovtest[,4]+5*xcovtest[,5]+epsilontest


bart_bma_example <- bartBMA(x.train = xcov,y.train=y,x.test=xcovtest, 
                    zero_split = 1, only_max_num_trees = 1,split_rule_node = 0)


```

