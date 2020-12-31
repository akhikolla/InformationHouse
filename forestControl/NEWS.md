### NEWS

#### V0.2.1

- Support of `randomForest` and `ranger` models which have been created using `parsnip`

- Fix class checking for when formula notation is used in randomForest

#### V0.2.0

- Significant reduction in compute time for calculating false positive rates by sampling only unique selection frequencies

- Addition of `tidy` tools (dplyr, tibble, magrittr) 

#### V0.1.1

- `internals` now implemented in C++ _via_ `Rcpp` thanks to Dr Jasen Finch (@jasenfinch)

#### v0.1.0

- Implemented Strategy-1 from _Konukoglu,E. and Ganz,M.,2014_.  __Approximate false positive rate control in selection frequency for random forest__

- Support for `randomForest` and `ranger` forest objects

- Calculate selection frequency threshold for a given false positive rate (alpha)

- False positive rate feature selection

- Wrapper for selection frequencies extract from objects
