
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CGGP

[![Travis-CI Build
Status](https://travis-ci.org/CollinErickson/CGGP.svg?branch=master)](https://travis-ci.org/CollinErickson/CGGP)
[![Coverage
Status](https://img.shields.io/codecov/c/github/CollinErickson/CGGP/master.svg)](https://codecov.io/github/CollinErickson/CGGP?branch=master)

The goal of CGGP is to provide a sequential design of experiment
algorithm that can efficiently use many points and interpolate exactly.

## Installation

You can install CGGP from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("CollinErickson/CGGP")
```

## Example

To create a CGGP object:

``` r
## basic example code
library(CGGP)
d <- 4
CG <- CGGPcreate(d=d,200)
print(CG)
#> CGGP object
#>    d = 4
#>    output dimensions = 1
#>    CorrFunc = CauchySQ
#>    number of design points             = 193
#>    number of unevaluated design points = 193
#>    Available functions:
#>      - CGGPfit(CGGP, Y) to update parameters with new data
#>      - CGGPpred(CGGP, xp) to predict at new points
#>      - CGGPappend(CGGP, batchsize) to add new design points
#>      - CGGPplot<name>(CGGP) to visualize CGGP model
```

A new `CGGP` object has design points that should be evaluated next,
either from `CG$design` or `CG$design_unevaluated`.

``` r
f <- function(x) {x[1]^2*cos(x[3]) + 4*(0.5-x[2])^3*(1-x[1]/3) + x[1]*sin(2*2*pi*x[3]^2)}
Y <- apply(CG$design, 1, f)
```

Once you have evaluated the design points, you can fit the object with
`CGGPfit`.

``` r
CG <- CGGPfit(CG, Y)
CG
#> CGGP object
#>    d = 4
#>    output dimensions = 1
#>    CorrFunc = CauchySQ
#>    number of design points             = 193
#>    number of unevaluated design points = 0
#>    Available functions:
#>      - CGGPfit(CGGP, Y) to update parameters with new data
#>      - CGGPpred(CGGP, xp) to predict at new points
#>      - CGGPappend(CGGP, batchsize) to add new design points
#>      - CGGPplot<name>(CGGP) to visualize CGGP model
```

If you want to use the model to make predictions at new input points,
you can use `CGGPpred`.

``` r
xp <- matrix(runif(10*CG$d), ncol=CG$d)
CGGPpred(CG, xp)
#> $mean
#>             [,1]
#>  [1,] -0.5797854
#>  [2,] -0.3749652
#>  [3,]  0.1128483
#>  [4,]  0.9209868
#>  [5,]  0.9907518
#>  [6,] -0.2265420
#>  [7,]  0.1801634
#>  [8,] -0.3477451
#>  [9,]  0.7743188
#> [10,] -0.0355732
#> 
#> $var
#>               [,1]
#>  [1,] 4.354093e-03
#>  [2,] 1.928985e-02
#>  [3,] 3.206114e-03
#>  [4,] 4.775720e-03
#>  [5,] 1.564388e-04
#>  [6,] 4.594134e-06
#>  [7,] 1.772568e-04
#>  [8,] 4.485993e-03
#>  [9,] 1.386959e-02
#> [10,] 1.293479e-02
```

To add new design points to the already existing design, use
`CGGPappend`. It will use the data already collected to find the most
useful set of points to evaluate next.

``` r
# To add 100 points
CG <- CGGPappend(CG, 100)
```

Now you will need to evaluate the points added to `CG$design`, and refit
the model.

``` r
ynew <- apply(CG$design_unevaluated, 1, f)
CG <- CGGPfit(CG, Ynew=ynew)
```

### Plot functions

There are a few functions that will help visualize the CGGP design.

#### `CGGPplotblocks`

`CGGPplotblocks` shows the block structure when projected down to all
pairs of two dimensions. The plot is symmetric. The facet labels be a
little bit confusing. The first column has the label 1, and it looks
like it is saying that the x-axis for each plot in that column is for
`X1`, but it is actually the y-axis that is `X1` for each plot in that
column.

``` r
CGGPplotblocks(CG)
#> Registered S3 methods overwritten by 'ggplot2':
#>   method         from 
#>   [.quosures     rlang
#>   c.quosures     rlang
#>   print.quosures rlang
```

![](tools/README-plotblocks-1.png)<!-- -->

#### `CGGPplotheat`

`CGGPplotheat` is similar to `CGGPplotblocks` and can be easier to read
since it is only a single plot. The \((i,j)\) entry shows the maximum
value for which a block was selected with \(X_i\) and \(X_j\) at least
that large. The diagonal entries, \((i, i)\), show the maximum depth for
\(X_i\). A diagonal entry must be at least as large as any entry in its
column or row. This plot is also symmetric.

``` r
CGGPplotheat(CG)
```

![](tools/README-heat-1.png)<!-- -->

#### `CGGPhist`

`CGGPhist` shows histograms of the block depth in each direction. The
dimensions that have more large values are dimensions that have been
explored more. These should be the more active dimensions.

``` r
CGGPplothist(CG)
#> Warning: Transformation introduced infinite values in continuous y-axis
#> Warning: Removed 8 rows containing missing values (geom_bar).
```

![](tools/README-hist-1.png)<!-- -->

#### `CGGPplotcorr`

`CGGPplotcorr` gives an idea of what the correlation structure in each
dimension is. The values plotted do not represent the actual data given
to CGGP. Each wiggly line represents a random Gaussian process drawn
using the correlation parameters for that dimension from the given CGGP
model. Dimensions that are more wiggly and have higher variance are the
more active dimensions. Dimensions with nearly flat lines mean that the
corresponding input dimension has a relatively small effect on the
output.

``` r
CGGPplotcorr(CG)
```

![](tools/README-corrplot-1.png)<!-- -->

#### `CGGPplotvariogram`

`CGGPplotvariogram` shows something similar to the semi-variogram for
the correlation parameters found for each dimension. Really it is just
showing how the correlation function decays for points that are further
away. It should always start at `y=1` for `x=0` and decrease in `y` as
`x` gets larger

``` r
CGGPplotvariogram(CG)
```

![](tools/README-vario-1.png)<!-- -->

#### `CGGPplotslice`

`CGGPplotslice` shows what the predicted model along each individual
dimension when the other input dimensions are held constant, i.e., a
slice along a single dimension. By default the slice is done holding all
other inputs at 0.5, but this can be changed by changing the argument
`proj`. The black dots are the data points that are in that slice If you
change `proj` to have values not equal to 0.5, you probably won’t see
any black dots. The pink regions are the 95% prediction intervals.

This plot is the best for giving an idea of what the higher dimension
function look like. You can see how the output changes as each input is
varied.

``` r
CGGPplotslice(CG)
```

![](tools/README-plotslice-1.png)<!-- -->

The next plot changes so that all the other dimensions are held constant
at 0.15 for each slice plot. When moving from the center line, the error
bounds generally should be larger since it is further from the data, but
we should see similar patterns unless the function is highly nonlinear.

``` r
CGGPplotslice(CG, proj = rep(.15, CG$d))
```

![](tools/README-plotslice2-1.png)<!-- -->

#### `CGGPplottheta`

`CGGPplottheta` is useful for getting an idea of how the samples for the
correlation parameters (theta) vary compared to the maximum a posteriori
(MAP). This may be helpful when using `UCB` or `TS` in `CGGPappend` to
get an idea of how much uncertainty there is in the parameters. Note
that there are likely multiple parameters for each input dimension.

``` r
CGGPplottheta(CG)
```

![](tools/README-plottheta-1.png)<!-- -->

#### `CGGPplotsamplesneglogpost`

`CGGPplotsamplesneglogpost` shows the negative log posterior for each of
the different samples for theta. The value for the MAP is shown as a
blue line. It should be at the far left edge if it is the true MAP.

``` r
CGGPplotsamplesneglogpost(CG)
```

![](tools/README-samplesneglogpost-1.png)<!-- -->
