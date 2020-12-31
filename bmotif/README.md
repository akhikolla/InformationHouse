
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/SimmonsBI/bmotif.svg?branch=master)](https://travis-ci.org/SimmonsBI/bmotif)
[![codecov](https://codecov.io/gh/SimmonsBI/bmotif/branch/master/graph/badge.svg)](https://codecov.io/gh/SimmonsBI/bmotif)

## Overview

`bmotif` is software for motif analyses of bipartite networks. It can
count occurrences of motifs in bipartite networks, as well as the number
of times each node or link appears in each unique node or link position
within motifs (a node or link’s structural role). `bmotif` supports
weighted as well as unweighted networks: the mean weight of motifs can
be calculated, as well as the standard deviation of motifs mean weights;
weighted versions of node and link position counts are also supported
(see below). As well as R, core functionality is also available in
[MATLAB](https://github.com/SimmonsBI/bmotif-matlab) and
[Python](https://github.com/SimmonsBI/bmotif-python). `bmotif` was
originally developed to analyse bipartite species interaction networks
in ecology but its methods are general and can be applied to any
bipartite graph.

Note that the released version on CRAN currently only has support for
counting motif occurrences and node positions. It does not have support
for link positions or weighted networks, which are currently only
available in the development version (these will be added to the CRAN
version when development is finalised).

## Installation

To install the released version from CRAN:

``` r
install.packages("bmotif")
```

To install the development version from GitHub:

``` r
install.packages("devtools") # install the devtools package
devtools::install_github("SimmonsBI/bmotif", build_vignettes = TRUE) # install bmotif
```

### Common issue when installing

Some users have reported issues installing `bmotif` without the `Rcpp`
package being installed. Therefore we recommend installing `Rcpp` before
installing `bmotif`:

``` r
install.packages("Rcpp")
```

### Common issue when installing on Windows

If installing on Windows, you might get an error saying something like
‘Rtools is required to build R packages, but no compatible version was
found’. Even after you install the latest version of Rtools you may get
‘ERROR: compilation failed for package ’bmotif’. If you get these
errors, try these steps:

1.  Install the latest version of the ‘devtools’ package:

<!-- end list -->

``` r
install.packages("devtools")
```

2.  Download the latest version of Rtools and follow the guidelines
    [here](https://thecoatlessprofessor.com/programming/installing-rtools-for-compiled-code-via-rcpp/)
    (archived version
    [here](https://web.archive.org/web/20180814151143/https://thecoatlessprofessor.com/programming/installing-rtools-for-compiled-code-via-rcpp/))
    to install. It describes how to set some PATH variables.
    Essentially, when given the option to modify the Window PATH
    variables during the Rtools installation, choose to do this, then
    add the following to the top of the PATH variables text field:
      - `c:\Rtools\bin;`
      - `c:\Rtools\mingw_32\bin;`
3.  Now you can use `Sys.getenv('PATH')` in R to make sure the correct
    PATH variables appear.
4.  If the PATH variables are there, run `devtools::find_rtools()` which
    should now return `TRUE`
5.  Try installing bmotif again and it should work\!

## Dictionary

`bmotif` considers all 44 unique bipartite motifs up to six nodes.
Within these motifs there are 148 unique node positions and 106 unique
link positions. All motifs, node positions and link positions considered
by `bmotif` are shown in the Figure below. This is the ‘dictionary’ used
by bmotif: the canonical reference for all motif, node position and link
position IDs used by the package and returned by the functions.

![Motif dictionary](./man/figures/dictionary.png?raw=true
"Motif dictionary")

Large numbers above and to the left of each motif represent the ID of a
motif. Small numbers at the end of links represent node position IDs
within motifs. Small numbers to the left of each motif represent link
position IDs within motifs: the colour of the link position number
corresponds to the colour of the links in the motif. Colours are
colourblind safe following the palette proposed by Wong et al (2011).
Motif IDs and node position IDS can also be found in **Simmons, B. I.,
Sweering, M. J. M., Schillinger, M., Dicks, L. V., Sutherland W. J., Di
Clemente, R. bmotif: a package for motif analyses of bipartite networks.
Methods in Ecology and Evolution, 10(5), 695-701.**. Link position IDs
are only in the Figure above. Node positions were defined following
Baker et al (2015) Appendix 1 Figure A27.

Consider the example of motif 5:

![Motif 5](./man/figures/motif5.png?raw=true "Motif 5")

We know this is motif 5 because the large number above and to the left
of the motif gives the motif ID as 5. This motif contains four unique
node positions, given by the numbers at the ends of each link: 9, 10, 11
and 12. This motif also contains three unique link positions given by
the three coloured numbers to the left of the motif. The leftmost link
between positions 9 and 12 is in link position 5, the rightmost link
between positions 10 and 11 is in link position 6, and the middle
diagonal link between positions 10 and 12 is in link position 7.

## Use

`bmotif` has three functions, with can all be used with binary or
weighted networks:

1.  `mcount`
      - **Unweighted**: counts how many times each motif occurs in a
        bipartite network
      - **Weighted**: calculates the mean weight of motifs and the
        standard deviation of their weights
2.  `node_positions`:
      - **Unweighted**: counts the number of times each node occurs in
        each of the unique node positions within the motifs.  
      - **Weighted**: calculates a range of weighted metrics, such as
        the mean link strength of a each node in each position or the
        contribution of each node to each motif’s total weight
3.  `link_positions`:
      - **Unweighted**: counts the number of times each link occurs in
        each of the unique link positions within the motifs.
      - **Weighted**: calculates the number of times each link occurs in
        each unique link position, multiplied by each link’s strength

Weighted methods for `mcount`, `link_positions` and the ‘mora’ method in
`node_positions` were originally defined by Mora et al. (2018).

## License

The code is released under the MIT license (see LICENSE file).

## Citation

If you use the package in your work, please cite: Simmons, B. I.,
Sweering, M. J. M., Schillinger, M., Dicks, L. V., Sutherland W. J., Di
Clemente, R. bmotif: a package for motif analyses of bipartite networks.
Methods in Ecology and Evolution 10(5), 695-701.

If you use any of the weighted analyses originally defined by Mora et
al. (2018) (weighted `mcount`, weighted `link_positions` and the ‘mora’
method in `node_positions`), please additionally cite: Mora, B.B.,
Cirtwill, A.R. and Stouffer, D.B., 2018. pymfinder: a tool for the motif
analysis of binary and quantitative complex networks. bioRxiv, 364703.

## References

Baker, N.J., Kaartinen, R., Roslin, T. and Stouffer, D.B., 2015.
Species’ roles in food webs show fidelity across a highly variable oak
forest. Ecography, 38(2), pp.130-139.

Mora, B.B., Cirtwill, A.R. and Stouffer, D.B., 2018. pymfinder: a tool
for the motif analysis of binary and quantitative complex networks.
bioRxiv, 364703.

Wong, B. Points of view: Color blindness. Nat. Methods 8, 441 (2011).
