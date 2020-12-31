---
output:
  pdf_document: default
  html_document: default
---
# ALassoSurvIC

The `ALassoSurvIC` package provides penalized variable selection tools for the Cox proportional hazards model with interval censored and possibly left truncated data. The main function `alacoxIC` performs the variable selection via penalized nonparametric maximum likelihood estimation with an adaptive lasso penalty. The function also finds the optimal thresholding parameter automatically by minimizing the Bayesian information criterion (BIC). The unpenalized Non-Parametric Maximum Likelihood Estimate (NPMLE) for interval censored and possibly left truncated data is also available with another main function `unpencoxIC`. The asymptotic validity of the methodology is established in Li et al. (2019). 

## Installation
``` r
install.packages("ALassoSurvIC")
intsall.packages("parallel") # required for parallel computing
```

## Overview
The package contains two main functions (`alacoxIC` and `unpencoxIC`) and two methods (`baseline` and `plot`) for the objects returned by the main functions. The cluster object, created by `makeCluster` in the `parallel` package, can be supplied with the `cl` argument in the main functions to reduce computation time via parallel computing. The parallel computing will be used when searching the optimal thresholding parameter and calculating the hessian matrix of the log profile likelihood. How to use the parallel computing is illustrated in one of the examples given below.

* `alacoxIC` : The function performs variable selection for interval censored data or for interval censored and left truncated data. The users can supply the value of a theresholding parameter with the argument `theta` in the function. If `theta` is not supplied by users, the function will automatically find the optimal thresholding parameter using a grid search algorithm, based on the Bayesian information criterion (BIC).

* `unpencoxIC`: The function allows users to get unpenalized NPMLEs along with standard errors and 95\% confidence intervals. 

* `basline`: The method to extract the NPMLEs for the baseline cumulative hazard function from an object returned by the `alacoxIC` function or the `unpencoxIC` function.

* `plot`: The method to plot the estimated baseline cumulative hazard function or the estimated baseline survival function from an object returned by the `alacoxIC` function or the `unpencoxIC` function.

## Example
The examples below show how to use the main functions and the methods with two virtual data sets; `ex_IC` is interval censored data and `ex_ICLT` is interval censored and left truncated data. Any inference cannot be drawn from these data sets.

### Example 1: Interval censored data
```
library(ALassoSurvIC)

data(ex_IC) # 'ex_IC' is interval censored data
lowerIC <- ex_IC$lowerIC
upperIC <- ex_IC$upperIC
X <- ex_IC[, -c(1:2)]

## Performing the variable selection algorithm using a single core
## Use the `cl` argument to reduce computation time.
res <- alacoxIC(lowerIC, upperIC, X)
res           # main result
baseline(res) # obtaining the baseline cumulative hazard estimate
plot(res)     # plotting the estimated baseline cumulative hazard function by default
plot(res, what = "survival")  # plotting the estimated baseline survival hazard function

## Getting the unpenalized NPMLEs for interval censored data
res2 <- unpencoxIC(lowerIC, upperIC, X)
res2
```

### Example 2: Interval censored and left truncated data
```
data(ex_ICLT) # 'ex_ICLT' is interval censored and left truncated data
lowerIC <- ex_ICLT$lowerIC
upperIC <- ex_ICLT$upperIC
trunc <- ex_ICLT$trunc
X <- ex_ICLT[, -c(1:3)]

## Performing the variable selection algorithm using a single core
## Use the `cl` argument to reduce computation time.
res3 <- alacoxIC(lowerIC, upperIC, X, trunc)
res3
baseline(res3)
plot(res3)
plot(res3, what = "survival")

## Getting the unpenalized NPMLEs for interval censored data
res4 <- unpencoxIC(lowerIC, upperIC, X, trunc)
res4
```

### Example 3: Reducing computation time using parallel computing
```
data(ex_IC) # 'ex_IC' is interval censored data
lowerIC <- ex_IC$lowerIC
upperIC <- ex_IC$upperIC
X <- ex_IC[, -c(1:2)]

library(parallel)
cl <- makeCluster(2L)  # making the cluster object 'cl' with two CPU cores
# cl <- makeCluster(detectCores()) # run this code instead to use all available CPU cores

## Compare two computation times
## Note that the `unpencoxIC` function also allows users to use the `cl` argument.
system.time(res_parallel <- alacoxIC(lowerIC, upperIC, X, cl = cl)) # Use two cores
system.time(res <- alacoxIC(lowerIC, upperIC, X)) # Use a single core
```

## Reference
Li, C., Pak, D., & Todem, D. (2019). Adaptive lasso for the Cox regression with interval censored and possibly left truncated data. Statistical methods in medical research, in press.

