# frailtyEM 1.0.1 (CRAN release)
- fixed some interals to work with `survival` package version 3.0
- added a `conf_level` argument that can set the width of the confidence intervals for the predicted cumulative hazards

# frailtyEM 1.0.0 (CRAN release)
- release version
- updated CITATION info and documentation

## frailtyEM 0.8.10
- fixed a bug that would occur in predict() when strata() was also used. 

## frailtyEM 0.8.9
- minor fixes, and by default predict() now doesn't calculate CI for the hazard function (it's slow and most people don't need that anyway)

## frailtyEM 0.8.8
- fixed an error that would happen in recognizing the `autoplot()` method on Linux systems

## frailtyEM 0.8.7
- fixed the vignette and made the plots nicer

## frailtyEM 0.8.6
- Nicer print in summary() and print()
- Added a sort of safety net when the likelihood is very flat. The program will switch then from `nlm()` to `optimize()`, which is generally more stable in these cases. The reason why I do not switch all the time top `optimize()` is because that one cannot be programmed with. Then, in combination with `numDeriv::hessian`, estimation would take forever and be basically impossible to check. 

# frailtyEM 0.8.4 (CRAN release)
- Fixed a bug where se = FALSE would break the predict() method

# frailtyEM 0.8.3 (CRAN release)
- Added the `zph` option in `emfrail_control()` so that the result of the `cox.zph` for the frailty model is also returned. This can be used to for goodness of fit. A guide on that soon to come!
- Bugfix (when empty strata was part of the input)

### frailtyEM 0.8.1
Major update. Now stratified models are supported! 
Several improvements in the documentation and in the performance section. 

Smaller fixes, as compared to the previoius CRAN release:

- removed `rev(cumsum(rev(rowsum)))` statement and replaced with an Rcpp function `rowsum_vec`
- using the cholesky decomposition instead of `solve`, seems that this is way better for symmetric matrices (0.7.16)
- simplified a bit the `emfrail_control()` function (0.7.15)
- documented options in `summary()` that control what is printed (if you want the output of a package to fit on one slide, for example) (0.7.15)
- fixed a bug in `ca_test` that was not reading correctly the input because of the partial matching of arguments in R (0.7.14)
- fixed some inconsistencies in summary and predict methods, when certain options are passed (such as no standard errors, for example) (0.7.13)
- improved print method for the summary object with some options to make things shorter (0.7.12)
- improved behaviour when frailty variance is actually 0 (0.7.11)
- fixed some warnings that were actually expected behaviour (0.7.10)
- some improvements in the limiting case where there is no frailty (0.7.10)
- more consistent notation and less vague passing of limits for the confidence intervals (now it's clear whether it's log scale or not) (0.7.10)

# frailtyEM 0.7.9 (CRAN release)

As compared to the previous CRAN release, 0.7.2:
- fixed a bug where the estimation would go wrong when the data set was not ordered according to the cluster
- fixed a bug where `emfrail` would crash when the cluster colum would be a character vector
- fixed a bug where the test for dependent censoring would not work
- part of the output is now nicer (e.g. the `frail` vector is named, the `autoplot.emfrail()` gives a nicer plot)
- removed a bunch of redundant calculations and old pieces of code
- minor corrections in the vignette

# frailtyEM 0.7.2
As compared to the previous CRAN release, 0.7.0:

- `ca_test()` now provides an interface to use the Commenges-Andersen test for heterogeneity outside the `emfrail()` function. It takes as input a `coxph` object. Therefore, it can work with other baseline hazard estimators and with strata. 
- Various fixes in the documentation and vignette, mostly typos. 

As usual, feedback is welcome. 

### frailtyEM 0.7.1-6
- Various fixes for `ca_test()`: no more model frame needed, works well with strata.
- fixed some small things in vignette

### frailtyEM 0.7.1-4
- fixed some comments and some documentation
- fixed the `ca_test()`, a small bug that was leading to wrong answers sometimes. Now it should give the sam result as the one in `emfrail`.

### frailtyEM 0.7.1-3
- `ca_test()` now works for `coxph` models properly as long as they have covariates
- fixed a bug where the CA test would not give the correct results in `emfrail`. 

### frailtyEM 0.7.1-0
- added a ca_test() function for `coxph` objects. Basically this is also done in `emfrail()`, but now you can also use `strata` or other things that are not supported by `emfrail().`

### frailtyEM 0.7.0-2
- added a warning for when the limits for searching of the likelihood based confidence interval are reached.
- removed that message with calculating information matrix

# frailtyEM 0.7.0
- big update comprising all previous changes: many new methods, organized plot methdos, speed improvements.
- updated documentation
- minor bug fixes

# frailtyEM 0.6.8
- now it's `emfrail_dist()` rather than `emfrail_distribution()`
- a bunch of small fixes and improvements

# frailtyEM 0.6.8
- added a larger number of methods for `emfrail` objects.

# frailtyEM 0.6.7
- the `predict.emfrail` method suffered some alterations: first of all, it now gives predictions for each
`lp` or each row of `newdata`, and it also gained the argumnet `individual`. If true, then the `newdata` argument
is taken as coming from the same individual. This can be used with time-dependent covariates and adjusting
the time at risk.

# frailtyEM 0.6.7
- the `emfrail` object type has been re-vamped into a more conventional object
- cleaned up the code of the methods 

# frailtyEM 0.6.6
- now no more arguments starting with dots and a more conventional `emfrail(formula, data, stuff)` phrasing of the main fitting function.
- added more checks of the input and warnings that try to tell the user whether the old `.formula` or `.data` arguments are still used.

# frailtyEM 0.6.5
- removed all the plot functions and replaced them by methods with `plot.emfrail()` and `autoplot.emfrail()` (for `ggplot2`).

# frailtyEM 0.6.4
- massive performance improvement, when there are a lot of distinct event time points. moved part of the calculation of the information matrix to c++

# frailtyEM 0.6.3
- fixed intervals for calculating confidence based intervals. seems there is a problem when the frailty variance is very large (e.g. 30, 40) 
- fixed an issue where nlm was not taking the parameters from the .control argument
- set the step size smaller for the nlm maximizer so that it doesn't overshoot (see 1st issue)


# frailtyEM 0.6.2
- big overhaul of the `control` argument and the `emfrail_control()` function

# frailtyEM 0.6.1 
- removed some old dependencies in the documentation and DESCRIPTION

# frailtyEM 0.6.0 (release)
- overall, numerous improvements compared to the previous release. Key new features include likelihood based confidence interval for the frailty parameter, more measures of dependence calculated with `summary()`, plots using `ggplot2`, and numerous bug fixes. 

# frailtyEM 0.5.13
- now the call is printed also when the summary is printed

# frailtyEM 0.5.12
- performance improvements. Now the likelihood-based confidence intervals should take less time as they know better where to look. 

# frailtyEM 0.5.11
- moved from `optimize` + `numDeriv` to `nlm`

# frailtyEM 0.5.11
- added a number of dependence measures that can be compared across distributions such as Kendall's tau, median concordance. 
- changed quite a lot in the structure of the summary object and the print method to make it more consistent and easier to develop in the future

# frailtyEM 0.5.10
- added score test for dependent censoring

# frailtyEM 0.5.9
- `ggplot_emfrail()`  added! Now the same plots (and more) can be done with the good looking `ggplot2` engine. 

# frailtyEM 0.5.8
- `summary.emfrail()` now has a new argument `print_opts` that is used in `print.emfrail_summary()`; if the output becomes too big, then some parts of the output may be ommitted

# frailtyEM 0.5.7

- The optimization now is regulated by search intervals described in the `emfrail_control()` and the `.control` argument. 
- The parametrization of the stable distribution has been changed, just removed the $1-$ in the beginning (why did I have that there again?)
- There are different intervals for the gamma/pvf and stable distributions. That's roughly because the stable chokes with small values of `theta`. This should be tuned somehow in the future. The problem lies in the M step where `agreg.fit` can't deal with large offset values.
- Likelihood confidence based intervals now do the correct thing when the estimate is close to the parameter space but not quite there
- Eliminated the fast fit for the inverse gaussian, this also seems to choke (the fast E step, can't figure out why), while the slow fit in C++ works fine...
- A slight update in documentation. 

TODO: 
- recover lost features in this update: measures of dependence in `summary.emfrail`, first of all
- bring back the fast fit for inverse gaussian or... who knows, maybe now
- document `emfrail_control` properly
- update vignette

# frailtyEM 0.5.6
Likelihood based confidence intervals are here! 

# frailtyEM 0.5.5
Removed the maximization by `optimx` and doing it with `optimize()`, since it's one dimensional. 
A hessian estimate is obtained from `numDeriv()`.

# frailtyEM 0.5.4
Minor bug fixes 

# frailtyEM 0.5.3
Some big changes in how the confidence intervals are constructed in predict.emfrail. Now - they are first constructed with the delta method for the log(cumulative hazard) and then exponentiated, so they do not have to be truncated at 0 or 1 any more. 

# frailtyEM 0.5.2
Further improved compatibility with CRAN policies and added a bunch of stuff in the examples in `\dontrun` statements (now they should be less than 5 seconds runtime)

# frailtyEM 0.5.1
Improved compatibility with R-devel 3.4.0. Registered C++ files to get rid of an R CMD check NOTE. Small modifications in the C++ file - for some reason a segfault started happening out of nowhere, think it's fixed now.

# frailtyEM 0.5.0
Added vignette, fixed small things for R CMD check
R CMD check: PASS, 0 warnings, 1 note / about new developer, that's ok.

# frailtyEM 0.4.9
Added the Commenges-Andersen test for heterogeneity. 
The test is implemented in a pretty non-efficient way, and it can be skipped with a proper `emfrail_control()` call, see `?emfrail_control`. Also there I added an option to *just* calculate the test, instead of doing anything else, and then just that is returned. A nice idea would be to implement this as a post-hoc calculation for `coxph` objects but that seems like another project atm.

R CMD check: PASS, 0 warnings, 0 notes.

# frailtyEM 0.4.8
Changed name to the more professional `frailtyEM`.
Added CI and SE for Kendall's tau with gamma

bugfixes: CI for tau with stable is now ok

# frailtoys 0.4.7
Added a `newdata` option for the `predict` method and for the `plot` methods. 
This can be used instead of `lp`, and basically calculates the corresponding linear predictor for the 
given covariate values.

bugfixes

# frailtoys 0.4.6
There is an option now to calculate the unadjusted SE or no SE at all

# frailtoys 0.4.5

* There are now plot methods available! Check out `?plot_emfrail`\
* Documentation updated accordingly





# frailtoys 0.4.3

* Added a `NEWS.md` file to track changes to the package.



