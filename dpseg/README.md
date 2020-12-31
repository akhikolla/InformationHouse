[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/dpseg)](https://cran.r-project.org/package=dpseg)
[![Downloads](https://cranlogs.r-pkg.org/badges/dpseg)](https://cran.r-project.org/package=dpseg)


# R package `dpseg`: piecewise linear segmentation by a simple dynamic programing algorithm

authors: "Rainer Machne, Peter F. Stadler"

This package performs piecewise linear segmentation of ordered data by
a dynamic programing algorithm, provided via the function `dpseg`. It
was developed for time series data, eg. growth curves, and for
genome-wide read-count data from next generation sequencing.

The package and its documentation are also intended to serve as a tutorial
on dynamic programing and the segmentation problem. A `movie` function
visualizes the progress of the algorithm through the data.

Moreover, the package features generic implementations of dynamic
programing routines, where new segmentation criteria ("scoring
functions") can be tested in base `R` and efficient versions
implemented in `Rcpp`.

## Documentation

See the package vignette (`vignette("dpseg")`) for details.


## Installation

Install package from within `R` via `cran`: 


```R
install.packages("dpseg")
```



### Developtment Version

```{r}
library(devtools)
install_gitlab("raim/dpseg")
```

## Basic Usage

```{r}
library(dpseg)

# get example data `oddata` - bacterial growth measured as optical density OD
x <- oddata$Time
y <- log(oddata[,"A2"]) # note: exponential growth -> log(y) is linear

segs <- dpseg(x=x, y=y, jumps=FALSE, P=0.0004)

## inspect resulting segments
print(segs)

## plot results
plot(segs)

## use predict method
lines(predict(segs),lty=2, lwd=3, col="yellow")

## view the algorithm in action
movie(segs)
```

## Theory

### Piecewise Linear Segmentation

The problem is to find break-points in 2-dimensional data, eg. timeseries, that 
split the data into linear segments. This can be formulated as an optimization 
problem that can be solved by `dynamic programing`:

```math
\newcommand{\Ell}{\mathcal{L}}
\newcommand{\jump}{\mathcal{J}}
\newcommand{\Var}{\mathrm{Var}}
\newcommand{\Cov}{\mathrm{Cov}}
\newcommand{\lmax}{\ell_\text{max}}
\newcommand{\lmin}{\ell_\text{min}}
S_j = \max_{i\le j} (S_{i-\jump} + \text{score}(i,j)) - P\;,
```

where the $`\text{score}`$ quantifies how well a segment between points
$`i`$ and $`j`$ is defined, eg. some goodness-of-fit measure such as the
negative variance of the residuals

```math
\text{score}(i,j) = -s_r^2
```

of a straight line fitted through data points from points $`i`$ to
$`j`$. $`P`$ is a penalty term, and $`P>0`$ allows to fine-tune segment
lengths. At constant scores it will accumulate in $`S_j`$ (subtracted
for each $`i`$), forcing the algorithm to "wait" until a higher score
is reached. Similarly, initialization of $`S_0>0`$ to a relatively
large value avoids spurious segments of length 1 at $`j=1`$ by
enforcing a break-point before $`i=1`$. $`S_1`$ has to be initialized
to $`S_1=-P`$.

Discontinuous jumps between adjacent segments can be allowed with
$`\newcommand{\jump}{\mathcal{J}} \jump =1`$, while segment borders
(break-points) are part of both left and right segments with
$`\newcommand{\jump}{\mathcal{J}} \jump =0`$, (in which case $`S_0`$
has no effect).


See the package vignette (`vignette("dpseg")`) for details.
