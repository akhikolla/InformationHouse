README.Rmd
================

Tools for Individual-Based Models in Infectious Disease
=======================================================

I have been developing an individual-based model to derive the cost-effective strategies to target malaria hotspots and eliminate malaria in Myanmar. In order to explore as many model structures as possible, I'm developing this tools which are generic enough to be used in any individual-based model for any infectious disease. At this moment, the package has 2 generic functions.

1. Create a synthetic population having several states.
-------------------------------------------------------

This function populates a matrix in which columns represent the states of the individuals and rows represent the individuals. Making it a generic function will let you explore as many disease state as you want. This is expecially useful when you're comparing your IBM with your ODE model.

``` r
library(ibmcraftr)
syn_pop(c(3,2,1)) # will populate 3 individuals in state 1, 2 in state 2 and 1 in state 3.
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    1    0    0
#> [3,]    1    0    0
#> [4,]    0    1    0
#> [5,]    0    1    0
#> [6,]    0    0    1
```

2. Make state transitions.
--------------------------

Using the state matrix of a population created previously, calculate the transitions from one state to other state(s) using the transition rate(s). This version has **stRCPP** function which is based on the codes in C++ to make it run faster.

``` r
pop <- syn_pop(c(19,1,0,0))
state_trans(1,2,.1,pop) #state transition from 1 to 2, at rate .1
#>       [,1] [,2] [,3] [,4]
#>  [1,]    0    0    0    0
#>  [2,]    0    0    0    0
#>  [3,]    0    0    0    0
#>  [4,]    0    0    0    0
#>  [5,]    0    0    0    0
#>  [6,]   -1    1    0    0
#>  [7,]    0    0    0    0
#>  [8,]    0    0    0    0
#>  [9,]    0    0    0    0
#> [10,]    0    0    0    0
#> [11,]    0    0    0    0
#> [12,]    0    0    0    0
#> [13,]    0    0    0    0
#> [14,]    0    0    0    0
#> [15,]    0    0    0    0
#> [16,]    0    0    0    0
#> [17,]    0    0    0    0
#> [18,]   -1    1    0    0
#> [19,]    0    0    0    0
#> [20,]    0    0    0    0
stRCPP(1,4,100,pop) #state transition from 1 to 4, at rate 100
#>       [,1] [,2] [,3] [,4]
#>  [1,]    0    0    0    0
#>  [2,]   -1    0    0    1
#>  [3,]   -1    0    0    1
#>  [4,]   -1    0    0    1
#>  [5,]   -1    0    0    1
#>  [6,]   -1    0    0    1
#>  [7,]    0    0    0    0
#>  [8,]    0    0    0    0
#>  [9,]   -1    0    0    1
#> [10,]    0    0    0    0
#> [11,]   -1    0    0    1
#> [12,]    0    0    0    0
#> [13,]   -1    0    0    1
#> [14,]   -1    0    0    1
#> [15,]   -1    0    0    1
#> [16,]    0    0    0    0
#> [17,]   -1    0    0    1
#> [18,]    0    0    0    0
#> [19,]   -1    0    0    1
#> [20,]    0    0    0    0
```

3. Make state transitions for 'n' number of timesteps.
------------------------------------------------------

**run\_state\_trans** function organizes how the transitions are calculated for the specified number of timesteps.

``` r
 pop <- syn_pop(c(19,1,0,0,0)) #synthesizing population
 b <- 2 #effective contact rate
 param <- list(
 list(1,c(2,5),c(NA,.1)), #transition from state 1 to 2 using FOI lambda
 list(2,3,100), #transition from state 2 to 3,
 list(3,4,100)  #the 3rd term ensures the transition to the next stage
 )

 timesteps <- 10
 transient <- c("param[[1]][[3]][1] <- rate2prob(b*sum(pop[,2],pop[,3])/sum(pop))")
 eval(parse(text=transient))

 run_state_trans(timesteps, param, pop, transient)
#>       [,1] [,2] [,3] [,4] [,5]
#>  [1,]   18    1    1    0    0
#>  [2,]   13    2    1    0    4
#>  [3,]   10    4    1    0    5
#>  [4,]    6    5    4    0    5
#>  [5,]    4    4    4    2    6
#>  [6,]    3    2    4    5    6
#>  [7,]    1    2    5    6    6
#>  [8,]    0    3    1   10    6
#>  [9,]    0    2    1   11    6
#> [10,]    0    2    0   12    6
 run_state_trans(timesteps, param, pop, transient, useC = FALSE)
#>       [,1] [,2] [,3] [,4] [,5]
#>  [1,]   16    1    1    0    2
#>  [2,]   10    3    1    1    5
#>  [3,]    5    6    2    2    5
#>  [4,]    2    5    5    3    5
#>  [5,]    0    2    8    5    5
#>  [6,]    0    0    6    9    5
#>  [7,]    0    0    2   13    5
#>  [8,]    0    0    2   13    5
#>  [9,]    0    0    1   14    5
#> [10,]    0    0    0   15    5
```
