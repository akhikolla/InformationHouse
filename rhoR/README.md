# rhoR <img src="man/figures/logo.png" align="right" alt="" width="120" />

[![pipeline status](https://gitlab.com/epistemic-analytics/qe-packages/rhoR/badges/master/pipeline.svg)](https://gitlab.com/epistemic-analytics/qe-packages/rhoR/pipelines?page=1&scope=branches&ref=master)
[![coverage report](https://gitlab.com/epistemic-analytics/qe-packages/rhoR/badges/master/coverage.svg)](https://gitlab.com/epistemic-analytics/qe-packages/rhoR/-/commits/master)
[![cran status](https://www.r-pkg.org/badges/version-ago/rhoR)](https://CRAN.R-project.org/package=rhoR) 
[![cran downloads](https://cranlogs.r-pkg.org/badges/grand-total/rhoR)](https://cranlogs.r-pkg.org/badges/grand-total/rhoR) 

## What is Rho
Rho is used to test the generalization of inter rater reliability
(IRR) statistics. Calculating rho starts by generating a large number of
simulated, fully-coded data sets: a sizable collection of hypothetical
populations, all of which have a kappa value below a given threshold -- which
indicates unacceptable agreement. Then kappa is calculated on a sample from
each of those sets in the collection to see if it is equal to or higher than
the kappa in then real sample. If less than five percent of the distribution
of samples from the simulated data sets is greater than actual observed kappa,
the null hypothesis is rejected and one can conclude that if the two raters had
coded the rest of the data, we would have acceptable agreement (kappa above the
threshold).

---

## Installation

### Install release version from CRAN

```
install.packages("rhoR")
```

### Install development version

```
devtools::install_git(url = "https://gitlab.com/epistemic-analytics/qe-packages/rhoR/", ref="master")
```

---
## Resources

To learn more about Rho, visit the [resources page](http://www.epistemicnetwork.org/resources/).
