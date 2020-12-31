# gretel: Generalized Path Analysis for Social Networks

Features methods for quantifying path values and identifying optimal paths under a 
variety of modeling assumptions. Intended to be used in service of other structural 
analyses.

## Getting Started

Installation Instructions
```
library(devtools)
install_github("davidbuch/gretel")
```
## Example

```
# Identify the path of optimal conductivity between nodes 1 and 5 of a sociomatrix
# Using example data from *Yang, Knoke* (2001) <DOI: 10.1016/S0378-8733(01)00043-0>

best_path <- opt_gpv(YangKnoke01, source = 1, target = 5, alpha = 1)

# Compare the conductivity of this path to that of an inferior path

gpv(YangKnoke01, path = best_path, alpha = 1)

gpv(YangKnoke01, path = c(1,2,3,4,5), alpha = 1)

```
Please see the package vignette for more information and examples.

## Overview

This package contains two categories of functions. The first category is concerned
with assigning values to user specified paths, while the second identifies
paths of optimal value.

Key functions in the path value calculation category are

* ```gpv```, which calculates Generalized Path Value
* ```ppv```, which calculates Probabilistic Path Value
* ```binary_distance```, ```peay_path_value```, ```flament_path_length```,
  ```peay_average_path_value```, and ```flament_average_path_length```, which
  calculate path value measures described in *Yang, Knoke* (2001).
* ```generate_proximities```, which generates a matrix of values representing the
  measures of optimal paths from each source node (row index) to each target node
  (column index).

Key functions in the optimal path identification category are

* ```opt_gpv```, which identifies the path of optimal Generalized Path Value from
  a particular source node to a particular target node
* ```opt_ppv```, which identifies the path of optimal Probabilistic Path Value from
  a particular source node to a particular target node
* ```all_opt_gpv```, which identifies the 'gpv'-optimal paths from every source node
  to every target node
* ```all_opt_ppv```, which identifies the 'ppv'-optimal paths from every source node
  to every target node
* ```unpack```, which unpacks the Dijkstra-format encoded shortest paths returned by
  ```all_opt_gpv``` and ```all_opt_ppv```. See their help pages for details.


### Author

David A. Buch

### Citation

TBA

### Acknowledgements

To my Dad, on his birthday.

### License

GPL-3