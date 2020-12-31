
<!-- README.md is generated from README.Rmd. Please edit that file -->

# backbone

<!-- badges: start -->

<!-- badges: end -->

The `backbone` package provides methods for extracting from a weighted
graph a binary or signed backbone that retains only the significant
edges. The user may input a weighted graph, or a bipartite graph from
which a weighted graph is first constructed via projection. Backbone
extraction methods include the stochastic degree sequence model [(Neal,
Z. P. (2014))](https://doi.org/10.1016/j.socnet.2014.06.001),
hypergeometric model [(Neal, Z.
(2013))](https://doi.org/10.1007/s13278-013-0107-y), the fixed degree
sequence model [(Zweig, K. A., and Kaufmann, M.
(2011))](https://doi.org/10.1007/s13278-011-0021-0), as well as a
universal threshold method.

In a graph `G`, edges are either present (i.e. `G_{ij}=1`) or absent
(i.e. `G_{ij}=0`). However in a weighted or valued graph, edges can take
a range of values that may capture such properties as the strength or
capacity of the edge. Although weighted graphs contain a large amount of
information, there are some cases (e.g. visualization, application of
statistical models not developed for weighted graphs) where it is useful
to reduce this information by focusing on an unweighted subgraph that
contains only the most important edges. We call this subgraph the
backbone of `G`, which we denote as `G’`. Extracting `G’` from `G`
requires deciding which edges to preserve. This usually involves
selecting a threshold `T_{ij}` such that edges are preserved if they are
above the threshold (i.e. `G_{ij}’=1` if `G_{ij} > T_{ij}`), and omitted
if they are below the threshold (i.e. `G_{ij}’=0` if `G_{ij} < T_{ij}`).
It is also possible to extract a signed backbone by selecting upper
`T_{ij}` and lower `T’_{ij}` thresholds such that `G_{ij}’=1` if
`G_{ij}>T_{ij}`, `G_{ij}’=-1` if `G_{ij} < T’_{ij}`, and `G_{ij}’=0` if
`G_{ij} > T’_{ij}` and `G_{ij} < T_{ij}`. The key to all backbone
extraction methods lies in the selection of `T`. The `backbone` package
provides several different methods for selecting `T` and thus extracting
`G’` from `G`.

## Installation

You can install the released version of backbone from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("backbone")
```

You can install from GitHub with:

``` r
library(devtools)
install_github("domagal9/backbone", build_vignettes = TRUE)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(backbone)
data(davis)
sdsm_props <- sdsm(davis)
#> Finding the distribution using SDSM with polytope model.
sdsm_bb <- backbone.extract(sdsm_props, signed = TRUE, alpha = 0.05)
```

For more detailed examples and background on the topic, see
`vignette("backbone_introduction", package = "backbone")` or our
manuscript on the backbone package: <https://arxiv.org/abs/1912.12779>
