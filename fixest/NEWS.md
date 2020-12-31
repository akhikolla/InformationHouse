
# News for the R Package `fixest`

## Changes in version 0.8.0

#### Bugs

 - Major bug when predict was used in the presence of fixed-effects (thanks to @jurojas5, [#54](https://github.com/lrberge/fixest/issues/54)). Introduced in version 0.7.

  - When using variable names to cluster the standard-errors inside functions, summary may not fetch the data in the right frame (thanks to @chenwang, [#52](https://github.com/lrberge/fixest/issues/52)). Now a completely new internal mechanic is in place.
  
  - When using variables with varying slopes and the number of iterations is greater than 300, a bug occurred in the function checking the convergence was right (thanks to @kendonB, [#53](https://github.com/lrberge/fixest/issues/53)). 
  
  - Fix bug in the demeaning algorithm when two variables with varying slopes were identical.
  
  - Fix bug in femlm/feNmlm when factor variables are removed due to the removal of some observations.
  
  - In `summary`, fix bug when the argument `cluster` was equal to a formula with expressions and not a variable name (thanks to @edrubin [#55](https://github.com/lrberge/fixest/issues/55)).
  
  - Fix bug when integers are present in the RHS (thanks to @zozotintin [#56](https://github.com/lrberge/fixest/issues/56)).
  
  - Fix bug when nb_FE >= 2 and the data was large (thanks to @zozotintin [#56](https://github.com/lrberge/fixest/issues/56)). 
  
  - Fix bug display of how the standard-errors were clustered in `etable`.
  
  - Fix bug occurring when lags were used in combination with combined fixed-effects (i.e. fe1 ^ fe2) (thanks to @SuperMayo [#59](https://github.com/lrberge/fixest/issues/59)).
  
  - Fix bug `coefplot` when representing multiple estimations and coefficient names are numbers.
 
#### IV

 - IV estimations are now supported. It is summoned by adding a formula defining the endogenous regressors and the instruments after a pipe.
```R
base = iris
names(base) = c("y", "x1", "x_endo", "x_inst", "species")
base$endo_bis = 0.5 * base$y + 0.3 * base$x_inst + rnorm(150)
base$inst_bis = 0.2 * base$x_endo + 0.3 * base$endo_bis + rnorm(150)

# The endo/instrument is defined in a formula past a pipe
res_iv1 = feols(y ~ x1 | x_endo ~ x_inst, base)

# Same with the species fixed-effect
res_iv2 = feols(y ~ x1 | species | x_endo ~ x_inst, base)

# To add multiple endogenous regressors: embed them in c()
res_iv3 = feols(y ~ x1 | c(x_endo, x_endo_bis) ~ x_inst + x_inst_bis, base)

```

#### fit statistics

  - The `fitstat` function has been significantly enhanced. 
  
  - Now the following types are supported:
  
    * Likelihood ratios
    
    * F-tests
    
    * Wald tests
    
    * IV related tests (F/Wald/Sargan)
    
    * common stats like the R2s, the RMSE, Log-likelihood, etc
    
  - You can register your own fit statistics. These can then be seamlessly summoned in `etable` via the argument `fitstat`.
  
  - The `print.fixest` function now supports the `fitstat` argument. This means that you can display your own desired fit statistics when printing `fixest` objects. This is especially useful in combination with the `setFixest_print` function that allows to define the default fit statistics to display once and for all. See the example in the "Instrumental variables" section of the Walkthrough vignette.
  
  - The new function `wald` computes basic Wald tests.
 
#### Multiple estimations

  - New arguments `split` and `fsplit`: you can now perform split sample estimations (`fsplit` adds the full sample).
    
  - Estimations for multiple left-hand-sides can be done at once by wrapping the variables in `c()`.
    
  - In the right-hand-side and the fixed-effects parts of the formula, stepwise estimations can be performed with the new stepwise functions (`sw`, `sw0`, `csw` and `csw0`).
    
  - The object returned is of class `fixest_multi`. You can easily navigate through the results with its subset methods.

```R
aq = airquality[airquality$Month %in% 5:6, ]
est_split = feols(c(Ozone, Solar.R) ~ sw(poly(Wind, 2), poly(Temp, 2)),
                 aq, split = ~ Month)
                 
# By default: sample is the root
etable(est_split)

# Let's reorder, by considering lhs the root
etable(est_split[lhs = TRUE])

# Selecting only one LHS and RHS
etable(est_split[lhs = "Ozone", rhs = 1])

# Taking the first root (here sample = 5)
etable(est_split[I = 1])

# The first and last estimations
etable(est_split[i = c(1, .N)])
```

#### Formula macros

  - The algorithm now accepts regular expressions with the syntax `..("regex")`:
  ```R
  data(longley)
  # All variables containing "GNP" or "ployed" in their names are fetched
  feols(Armed.Forces ~ Population + ..("GNP|ployed"), longley)
  ```
  
    
#### New features in etable

  - New `style.tex` and `style.df` arguments that define the look of either Latex tables or the output data.frames. 
  
    * it can be set with the new functions `style.tex` and `style.df` that contain their own documentation. 
    
    * some `etable` arguments have been ported to the `style` functions (`yesNo`, `tablefoot`).
    
  - New `postprocess.tex` and `postprocess.df` arguments which allow the automatic postprocessing of the outputs. See the dedicated vignette on exporting tables for an illustration.
    
  - new `tabular` arguments which allows to create `tabular*` tables (suggestion by @fostermeijer [#51](https://github.com/lrberge/fixest/issues/51)).

  - polynomials and powers are automatically renamed to facilitate comparison across models. You can set their style with the argument `poly_dict`.
  
  - the labeling of models is enhanced when `rep.fixest` is used with different standard-errors (the model names are now "model INDEX.SUB-INDEX").
  
  - the argument `subtitles` has been improved, and now automatically displays the samples when split sample estimations are performed.
 
#### Other new features

 - In all estimations:
    
    * `subset`: regular subset (long overdue).
    
    * `split`, `fsplit`: to perform split sample estimations.
    
    * `se`, `cluster`: to cluster the standard-errors during the call.
    
    * `lean`: if `TRUE`, then summary is applied and any large object is removed from the result. To save memory => but many methods won't work afterwards.
    
    * `fixef.rm`: argument that accepts `none`, `perfect`, `singleton`, `both`. Controls the removal of fixed-effects from the observation.
    
    * auto parsing of powers. Now you don't need to use `I()` to have powers of variables in the RHS, it is automatically done for you (i.e. `x^3` becomes `I(x^3)`):
    ```R
    base = iris
    names(base) = c("y", "x1", "x2", "x3", "species")

    # The multiple estimation below works just fine
    feols(y ~ csw(x, x^2, x^3), base)
    ```
  - Estimation options can be set globally with `setFixest_estimation()`.
  
  - The `demean` function has been enhanced (with the contribution of Sebastian Krantz).
  
#### Improvements of the internal algorithm

 - Internal demeaning algorithm: some copies of the data are avoided when using `feglm`.
 
 - Internal algorithm of `to_integer` (used in all estimations): one copy of the input data is now avoided.
 
 - All estimations: smarter handling of the intercept, thus avoiding the reconstruction of the design matrix.

## Changes in version 0.7.1 (27-10-2020)

#### Hotfixes

 - Fix bug int overflow in estimations with only one variable.
 
 - Fix bug in tests occurring in R old release.
 
 - Fix bug in examples occurring in R old release.

#### Improvements

 - Function `i()` now behaves as `factor()`, setting automatically a reference when appropriate.
 
 - Internal algorithm of `i()` is much faster.
 
#### New features

 - In `etable`, the user can now provide a type of clustering for each model.
 
 - New method `rep.fixest` to replicate fixest objects, mostly useful in `etable` when several SEs for the same models are to be reported.
 
 - Automatic fix when the variance is not positive definite.


## Changes in version 0.7.0 (24-10-2020)

#### Bugs

 - Major bug when fixed-effects were combined with `^` and they contained NAs (thanks to @poliquin [#35](https://github.com/lrberge/fixest/issues/35)).
 
 - Bug when using lead/lags in estimations. The bug was due to a bug in a dependency ([dreamerr](https://cran.r-project.org/package=dreamerr)) and was fixed. Now fixest requires **dreamerr** version >= 1.2.1. Bug spotted by @seunghoon001 ([#44](https://github.com/lrberge/fixest/issues/44)).
 
 - Major bug when n_obs x n_vars > 2B or n_obs x n_fixed-effects > 2B. In such cases estimations could just not be done, even leading R to crash when using nthreads > 1. The algorithm was fixed to allow datasets with up to 2B observations to be estimated in all circumstances. Bug reported, and many help for checking provided, by Howard Zihao Zhang.
 
 - `coefplot`: Problem regarding interactions when observations, and hence coefficients, were removed from the estimation. Now the coefficients are removed from the plot. Bug reported by @phisherblack [#45](https://github.com/lrberge/fixest/issues/45).
 
 - `coefplot`: Corrected various bugs when asked for the plotting of several estimations. 
 
 - Fix the stack imbalance warning (report by @shoonlee, [#46](https://github.com/lrberge/fixest/issues/46)).
 
#### Internal improvements

 - Brand new internal algorithm which now uses closed form solutions when dealing with variables with varying slopes. This means that when variables with varying slopes are present, the algorithm is incomparably faster and more accurate.

 - Two deep copies of some data are now avoided in the demeaning function. This improves the performance in terms of memory footprint, and also makes the algorithm faster. 
 
#### Standard-errors, important changes
  
 - New default values for standard-errors (only concerns multiway clustering). They become similar to `reghdfe` to increase cross-software comparability. Computing the standard-errors the old way is still possible using the argument `dof`. See the dedicated vignette: [On standard errors](https://cran.r-project.org/package=fixest/vignettes/standard_errors.html).
 
 - Name change in `summary`/`vcov`/`etable`: To get heteroskedasticity-robust standard-errors, `se = "hetero"` now replaces `se = "white"` to enhance clarity. Note that `se = "white"` still works.
 
#### New function: `fitstat` 

  - New function `fitsat` that computes various fit statistics. It is integrated with `etable` and can be invoked with the argument `fitstat`. So far only two fit statistics are included, but more will come.
 
#### New features in `interact()`

  - You can now use `i(var)` to treat the variable `var` as a factor. You can select which values to drop/keep with the respective arguments. 
  
  - Using `i(var)` leads to a special treatment of these variables in the functions `coefplot` and `etable`.
 
#### New features in `etable`

  - New argument `placement` to define the position of the float in Latex (suggestion by Caleb Kwon).
  
  - New argument `drop.section`, with which you can drop a) the fixed-effects, b) the variables with varying slopes, or c) the statistics, sections (suggestion by Caleb Kwon).
  
  - Fix glitch in help pages regarding the use of the '%' (percentage) character in regular expressions.
  
  - Two new arguments `.vcov` and `.vcov_args` to compute the standard-errors with custom functions.
  
  - The number of observations (`n`) is now treated as a regular statistic and can be placed where one wants.
  
  - The statistics can now have custom aliases using the argument `dict`.
  
  - The overdispersion becomes a regular fit statistic that can be included (or not) using `fitstat`.
  
  - The dictionnary now applies to the factors of interactions, and the values of factors.

#### User visible changes

 - Argument `nthreads`:
 
    * The new default of argument `nthreads` is 50% of all available threads. 
  
    * Accepts new values: a) 0 means all available threads, b) a number strictly between 0 and 1 will represent the fraction of all threads to use.
  
  - When setting formula macros:
  
    * The functions `xpd` and `setFixest_fml` now accept character vectors and numeric scalars on top of formulas.
  
  - `demean`:
  
    * speed improvement.
  
  - `coefplot`:
  
    * The argument `group` now accepts a special character `"^^"`, when used, it cleans the beginning of the coefficient name. Very useful for, e.g., factors although factors created with `i()` need not that.
    
     * When `horiz = TRUE`, the order of the coefficients is not reversed any more.

 - Improved display of numbers in `print` method.
 
 - Added variables names to `X_demeaned` from `feols`.
 
 - Lagging functions:
 
    * Now `time.step = NULL` by default, which means that the choice of how to lag is automatically set. This means that the default behavior for time variables equal to Dates or character values should be appropriate.
  
    * New operator `d` which is the difference operator.
 
 - In all estimations: 
 
    * new argument `mem.clean`: internally, intermediary objects are removed as much as possible and `gc()` is called before each memory intensive C++ section. Only useful when you're at the edge of reaching the memory limit.
    * new output: `collin.min_norm`, this value informs on the possible presence of collinearity in the system of variables.
    * new arguments `only.env` and `env`:
      * The first, `only.env`, allows to recover only the environment used to perform the estimation (i.e. all the preprocessing done before the estimation).
      * The second, `env`, accepts a fixest environment created by `only.env`, and performs the estimation using this environment--all other arguments are ignored. 
      * These changes are a prerequisite to the efficient implementation of bootstraping (since, by applying modifications directly in `env`, we cut all preprocessing).
    

 - In non-linear estimations: 
 
    * non-numeric variables can now be used. 
    * argument `NL.start` now accepts numeric scalars, initializing all coefficients to the same value (avoids the use of the other argument `NL.start.init`).
    
  - `summary.fixest`:
  
    * argument `.vcov` now accepts functions that compute the vcov. This ensures convenient compatibility with the `sandwich` package (compatibility is still not full though: bootstraped SEs don't work yet).
    
  - `update.fixest`:
  
    * new argument `evaluate` to ensure consistency with the `update` method from stats.
    
  - `feols` & `feglm`:
  
    * the Cholesky decomposition now checks for user interrupts (matters for models with MANY variables to estimate).
  
#### Deprecation

 - Argument `na_inf.rm` has been removed. It was present for historical reasons, and is removed to increase code clarity.

## Changes in version 0.6.0 (13-07-2020)

#### Bugs

 - In `vcov`, the degree-of-freedom in the small sample correction correction was fixed to "nested" and couldn't be modified, now corrected. Further, "nested" was not properly accounted for, now corrected.
 
 - In `etable`, `fitsat = FALSE` or `fitsat = NA` led to a bug.
 
 - `r2`: bug when the estimation contained only fixed effects (thanks to Luis Fonseca [#27](https://github.com/lrberge/fixest/issues/27)).
 
 - Now the `BIC` of `feglm` is similar to the one of `glm`.
 
 - Bug in the log-likelihood in the presence of weights, now corrected.
 
 - Bug in `coefplot` when some interacted variables were removed because of collinearity. Now corrected.
 
#### New vignettes

 - On standard-errors: how are the SEs computed in fixest and how to replicate the SEs from other software.
 
 - Exporting estimation tables: how to use efficiently `etable`, in particular how to customize the tables.
 
#### Major changes: etable

  - New arguments: `group`, `extraline`, `notes`, `tablefoot`. 
  
    - `group` allows to eliminate variables (like `drop`) and adds an extra line with TRUE/FALSE if the model contained those variables.
    
    - `extraline` allows to add extra lines with any content.
    
    - `notes` allows to add notes after the table (suggestion by @bgchamps [#25](https://github.com/lrberge/fixest/issues/25)).
    
    - `tablefoot` controls whether the table footer, containing the type of standard-errors and the significance codes, should be displayed.
    
  - Renaming: `yesNoFixef` => `yesNo`.
  
  - Most default values can be set globally with the new function `setFixest_etable`.
 
#### Major changes: dof

 - Function `dof`, used to adjust the small sample corrections, is now much more complete and allows to replicate a large set of estimation results from alternative software.

#### User visible changes

 - You can now provide custom VCOVs to summary by using the argument `.vcov`.
 
 - A warning is now prompted when the maximum number of iterations of the algorithm is reached (suggestion by @clukewatson  [#24](https://github.com/lrberge/fixest/issues/24)]). 
 
 - The types of standard-errors can now be set globally with the function `setFixest_se` (suggestion by @dlindzee [#28](https://github.com/lrberge/fixest/issues/28))
 
 - New `feols` argument `demeaned`. If `TRUE`, then the centered variables are returned (`y_demeaned` and `X_demeaned`). (Suggestion by Linus Holtermann.)
 
 - `interact` gains two new arguments: `drop` and `keep` (suggestion by @SuperMayo [#23](https://github.com/lrberge/fixest/issues/23)).
 
#### New methods

 - `hatvalues` has been implemented for feols and feglm estimations.
 
 - the `estfun` from `sandwich` has been implemented.

## Changes in version 0.5.1 (18-06-2020)

#### Hotfix
 
 - Fixed bug introduced in the previous update (memory access error). Does not affect any of the results but could lead R to crash unexpectedly (odds were low though since access was adjacent).
 
#### Bugs

 - Fix image link to the equation in the README.md.
 - Fix bug R2 and logLik when observations were removed because of NA values. Due to the update in `residuals.fixest`.
 
#### User visible change

 - Rewriting of the internal algorithm computing the VCOV. 1) About 30% performance gain for estimations with many variables. 2) The code is much less memory hungry. 

#### Major update of `etable`

 - New argument `style` which allows to set many elements of the output table. 
 - (minor) `signifCode` can be equal to `"letters"` to display letters instead of stars.
 
#### Other

 - `setFixest_nthreads` now respects the `OMP_THREAD_LIMIT` environment variable.
 - Rd links are now made to the proper htlm files.


## Changes in version 0.5.0 (10-06-2020)

#### Bug fixes
        
 - Bug with estimations with varying slopes if the fixed-effect relative to the slope is not in its decreasing order (thanks to Davide Proserpio).
 - Bug when interacting two variables with the `var::fe` syntax with `confirm = TRUE` and no reference.
 - Bug in `etable` when the standard-errors where `NA`.
 - Fixed very minor bug when computing the SEs (1e-6 difference).
 - Standard-errors in `feglm` for non-poisson, non-binomial families, are now correct (minor differences).
 - `fixef` did not work when the slope was an integer, now corrected (thanks to @clerousset [#20](https://github.com/lrberge/fixest/issues/20)).

#### New functionality: formula macros
        
 - You can use macros in formulas.
 - To set a macro variable, use e.g., `setFixest_fml(..ctrl = ~ var1 + var2)`. Here the macro variable `..ctrl` has been set to the value `"var1 + var2"`.
 - Now you can use this macro variable in any `fixest` estimation: e.g. `data(airquality) ; setFixest_fml(..ctrl = ~ Temp + Day) ; feols(Ozone ~ Wind + ..ctrl, airquality)`.
 - You can use macros in non-fixest estimations with `xpd`, which expands formulas. E.g. `lm(xpd(Ozone ~ Wind + ..ctrl), airquality)`.
 
#### New functions

 - `to_integer`: user-level version of the internal algorithm transforming any kind of vector (or combination of vectors) into an integer ranging from 1 to the number of unique elements of the vector. Very fast.
 - `demean`: user-level version of the demeaning algorithm used in `feols`.

#### Major user-visible changes
        
 - New internal algorithm to estimate OLS (applies to both `feols` and `feglm`):
 
    1. It is numerically more stable.
    
    2. Incomparably faster when factors are to be estimated (and not explicitly used as fixed-effects).
    
    3. Collinear variables are removed on the fly.

#### User-visible changes
        
 - Interactions in `var::fe(ref)` now accept multiple references (i.e. `ref` can be a vector).
 - In `etable`, the variable names of non-Latex output can now be changed.
 - You can use the argument `n` when applying summary to choose the number of coefficients to display.
 - Argument `confirm` has been removed from the function `interact`.
 - `r2` allows more flexibility in the keywords it accepts.
 - Function `dof` gains a new argument `adj` which allows to make different types of common small sample corrections. Its other arguments have been renamed for clarity (`fixef` => `fixef.K`, `exact` => `fixef.exact`, `cluster` => `cluster.adj`).
 - Now t-statistics are used for `feols` and non-poisson, non-binomial models in `feglm`. For all other models, z-statistics are used. This complies with the default's R-stats behavior.

#### New Methods
        
 - The `residuals` method has been substantially improved, now allowing different types.
 - New stats methods: sigma, deviance, weights.

#### Vignette and Readme
        
 - Typos corrected.
 - Images in the Readme set to 1200px.
 
#### Issue found: convergence problems with multiples variables with varying slopes
        
 - Convergence problems may arise in the presence of **multiple** variables with varying slopes. Theoretical work helped find a solution to this problem, but the implementation in R is proving not instantaneous.
 - In the meantime, now a warning is prompted when the algorithm suspects a convergence problem leading to poor precision of the estimated coefficients.

#### Error-handling
        
 - Improved error-handling with [dreamerr](https://cran.r-project.org/package=dreamerr)'s functions.
 
#### Other

 -  Dependency to MASS has been removed.

## Changes in version 0.4.1 (2020-04-13)

#### Bug fixes
  
 - Major bug leading R to crash when using non-linear-in-parameters right-hand-sides in feNmlm. Only occured when some observations were removed from the data set (due to NAness or to perfect fit). [Thanks to @marissachilds, GH issue [#17]](https://github.com/lrberge/fixest/issues/17).]
 - In the `collinearity` help pages: an example could lead to an error (due to random data generation). It has been removed.
 - In `collinearity`, corrected the problem of display of the intercept in some situations.
 - Defaults for the arguments `cex` and `lwd` in `coefplot` have been changed to 1 and 1 (instead of par("cex") and par("lwd")). Otherwise this led to the creation of `Rplots.pdf` in the working directory (thanks to Kurt Hornik).
 - Corrected a typo in the article's title in the vignette.

#### Help
  
 - Rewriting of sections, correction of small mistakes (wrong argument names), dropping completely the 'cluster' terminology (meant for fixed-effects), addition of what is contained in fixest objects.

#### Other
  
 - Adding a README.md.
 - Small corrections in the vignette.

## Changes in version 0.4.0 (2020-03-27)

#### User visible changes: Latex export
        
 - Better Latex special character escaping (errors reported by @dlindzee [#15](https://github.com/lrberge/fixest/issues/15)).
 - New argument `fixef_sizes.simplify`, which provides the sizes of the fixed-effects in parentheses when there is no ambiguity.
 - You can suppress the line with the significance codes with `signifCode = NA`.
 - New argument `float` which decides whether to embed the table into a table environment. By default it is set to `TRUE` if a `title` or `label` is present.
 - New argument `keep` to select the variables to keep in the table.
 - New way to keep/drop/order variables with the special argument "\%". If you use "\%var", then it makes reference to the original variable name, not the aliased one (which is the default).
 - New argument `coefstat` defining what should be shown below the coefficients (standard-errors, t-stats or confidence intervals). Suggestion by @d712 [#16](https://github.com/lrberge/fixest/issues/16).
 - Better rendering of significant digits.

#### User visible changes: coefplot
        
 - Argument `horiz`. The coefficients can now be displayed horizontally instead of vertically.
 - The coefficient labels, when in the x-axis, can now be displayed in three different ways thanks to the new argument `lab.fit`: "simple", the classic axis, "multi", the labels appear across multiple lines to avoid collision, and "tilted" for tilted labels.
 - The margins now automatically fit.
 - Argument `style` allows you to set styles with the function `setFixest_coefplot`, you can then summon the style in `coefplot` with this argument.
 - Use the ampersand to set dictionnary variables specific to `coefplot` (e.g. setFixest_dict(c(var1 = "a $\\times$ b", "&var1" = "&a\%*\%b")), the variable `var1` will be rendered differently in `etable` and in `coefplot`).
 - Better display of groups (with the arguments `group` and `group.par`).

#### New methods
        
 - `terms.fixest` giving the terms of the estimation.

#### Other
        
 - All `donttest` sections were removed from help pages.


## Changes in version 0.3.1 (2020-02-09)

#### Major bug fix
        
 - [panel] Fixed faulty memory access when taking the lead of a variable.

#### Other bug fixes
        
 - [esttable/esttex] These two functions were replaced by the function `etable`. In the process, some of their arguments were "lost", this is now corrected.
 - [etable] Better escaping of special characters.
 - [estimations] Bug when particular non-numeric vectors were used in explanatory variables.

#### New features
        
 - [coefplot] The function `coefplot` now accepts lists of estimations.


## Changes in version 0.3.0 (2020-02-01)

#### New feature: Lagging
         
 - You can now add lags and leads in any `fixest` estimations. You only need to provide the panel identifiers with the new argument `panel.id`, then you're free to use the new functions `l()` for lags and `f()` for leads.
 
 - You can also set up a panel data set using the function `panel` which allows you to use the lagging functions without having to provide the argument `panel.id`, and which dispose of more options for setting the panel.

#### New feature: Interactions
         
 - You can now add interactions in formulas with a new syntax: `var::fe(ref)`
 
 - The command `var::fe(ref)` interacts the variable `var` with each value of `fe` and sets `ref` as a reference. Note that if you don't use the argument `ref`, the command `var::fe` is identical to `var:factor(fe)`.
 
 - Using `var::fe(ref)` to write interactions opens up a special treatment of such variables in the exporting function `etable` and in the coefficient plotting function `coefplot`.

#### New feature: `coefplot`
         
 - You can plot coefficients and their associated confidence intervals with the function `coefplot`.
 
 - `coefplot` dispose of many options, whose default values can be set with the function `setFixest_coefplot`.
 
 - As for the function `etable`, you can easily rename/drop/order the coefficients.
 
 - `coefplot` detects when interactions have been used and offers a special display for it.

#### New functions
         
 - [etable] Estimations table: new function to export the results of multiple estimations. Replaces the two functions `esttex` and `esttable` (the two functions still exist but they will be deprecated in the future).
 
 - [Lagging] New functions related to lagging: `l`, `f`, `panel`, `unpanel` and `[.fixest_panel`.
 
 - [Utilities] A set of small utility functions has been added. They allow to extract part a coefficient table or parts of it (like the t-statistics of the standard-error) from an estimation. These functions are `coeftable`, `ctable` (an alias to `coeftable`), `se`, `tstat` and `pvalue`.
 
 - [coefplot] The functions `coefplot` and `setFixest_coefplot`.
 
 - [dof] New function to set the type of degree of freedom adjustment when computing the variance-covariance matrix. You can permanently set the type of DoF adjustment with the new function setFixest_dof().

#### User visible changes
         
 - [all estimations] A key pre-processing step has been paralellized => algorithm faster in general and much faster for multi-FEs.
 - [predict & fitted] Predict and fitted now returns vectors of the length equal to the one of original data.
 - [standard-errors] New ways to compute the standard-errors have been implemented. In particular, now it account for the "nestedness" of the fixed-effects in the clusters by default. You can freely change how to compute the degrees of freedom correction with the function dof().
 - [r2] Computation of the within-R2 for feglm models is now self-contained.
 - [all estimations] New, more accurate, stopping criterion for 2+ fixed-effects.
 - [feols] Estimations are slightly faster.
 - [etable/esttex] When there are interactions, R may change the order of the interactions, making two interactions in two different estimations look different while they are in fact the same (e.g. x3:x2 and x2:x3). Now esstable automatically reorders the interactions when needed for comparison across estimations.
 - [etable/esttable] The type of standard errors is now always shown.
 - [etable/esttex] The aliases provided by 'dict' are also applied within interactions. For example: `dict=c(x1="Wind", x2="Rain")`, with an estimation with the following variables 'x1', 'x2', 'x1:x2' will lead to the following aliases in Latex 'Wind', 'Rain' and 'Wind $\\times$ Rain'.
 - [etable/esttex] Interactions of similar values but of different order (e.g. x1:x2 and x2:x1) are reorderd to appear in the same lines.
 - [etable/esttex] The i) type of standard errors and ii) the significance codes, are now displayed in two separate lines (otherwise the line would be too wide).
 - [etable/esttex] Argument `yesNoFixef` can be of length one, defaulting the second element to the empty string.
 - [etable/esttex] Escaping of Latex special characters is now much more robust.

#### Bug correction
         
 - [all estimations] Fixed: bug when functions in the formula returned matrices.
 - [update] Fixed: error message when the data is missing.
 - [feglm] Fixed: bug double estimation when family not equal to poisson or logit
 - [feglm] Fixed: severe bug occurring for families not equal to poisson or logit
 - [predict] Fixed: bug when the estimation contained combined FEs.
 - [summary] Regarding small sample only: now Student t distribution is used instead of the Normal to compute the pvalue.
 - [esttex] Different variables with the same aliases (given by the argument 'dict') now appear in the same row.
 - [esttex] Arguments 'drop' and 'order' are now applied post aliasing (alias given by the argument 'dict').
 - [esttex] But when exporting multi-way standard errors.
 - [r2] Small bug regarding objects obtained with `did_estimate_yearly_effects`.
 - [estimations] bug when using weights in `feglm`.

## Changes in version 0.2.1 (2019-11-22)

#### Major bug correction
        
 - lag.formula: Bug introduced from previous update which could lead to wrong results. Now fixed.

#### Major user visible changes
         
 - [All estimation methods] Significant speed improvement when the fixed-effects variables (i.e. the identifiers) are string vectors.


## Changes in version 0.2.0 (2019-11-19)

#### New function
         
 -[did_means] New function `did_means` to conveniently compare means of groups of observations (both treat/control and pre/post). Contains tools to easily export in Latex.



#### Major user visible changes
        
 - [All estimation methods] Significant speed improvement when the fixed-effects variables (i.e. the identifiers) are of type integer or double.
 - [esttex, esttable] New argument 'fitstat' to select which fit statistic to display. The default adapts to the models. Old arguments (loglik, bic, aic, sq.cor) are dropped.
 - [esttable] Significantly better rendering of SE types.
 - [r2] Now NA is returned for R2s that have no theoretical justification (e.g. within R2 when no FEs, or 'regular' R2 for ML models).



#### Minor user visible changes
        
 - [did_plot_yearly_effects] Now the name of the dependent variable appears on the y-axis.
 - [esttex] Usage of the `sym` macro in Latex is dropped.



#### Bug correction
        
 - [fixef.fixest] bug could appear when using varying slopes coefficients in some specific circumstances (when the slope FEs were different from the regular FEs).
 - [fixef.fixest] bug when many regular FEs jointly with varying slopes.
 - [fixef.fixest] regarding slope coefficients: now the algorithm also evaluates functions of variables.
 - [esttable] Width of the "separating lines" now appropriately set for long dependent variable names.
 - [esttex] Spelling mistake corrected.
 - [estimations] Bug could occur for extremely small data sets (< 10 observations).



#### Error handling
        
 - [esttex, esttable] More informative error messages in functions esttex and esttable.


## Changes in version 0.1.2 (2019-10-04)
#### Major bug correction
        
 - lag.formula: When the data was not in a particular format, the results could be wrong. Now corrected.


## Changes in version 0.1.1 (2019-09-20)
#### Major bug correction
	    
 - feglm: bug when a) the deviance at initialization was higher than the deviance of the first iteration of the IRWLS and b) the step-halving was unable to find a lower deviance. This led the estimation to fail with an error although it should have been performed properly.
 - did_estimate_yearly_effects: bug when the estimation involved periods with negative values.

#### Minor bug correction
        
 - esttex: bug regarding the number of digits of negative coefficients to be displayed
 - esttex: now properly escaping the percentage and the underscore for exports in Latex
 - esttex: bug when changing the names of the dependent variables using a dictionnary
 - vcov: some warning messages were misleading
 - update: bug update when using the argument nframes
 - update: bug when updating the function fepois

#### Error handling
        
 - Better error messages for: did_estimate_yearly_effects, main estimation functions, setFixest_dict, fepois and fenegbin.

## Version 0.1.0 (2019-09-03)

#### First version
	    
 - This package is an effort to create a family of fast and user-friendly functions to perform estimations with multiple fixed-effects (F.E.).

 - Estimations with fixed-effects (or call it factor variables) is a staple in social science. Hence having a package gathering many methods with fast execution time is of prime importance. At the time of this version, this is the fastest existing method to perform F.E. estimations (often by orders of magnitude, compared to the most efficient alternative methods [both in R and Stata]). The underlying method to obtain the F.E. is based on Berge 2018, and the workhorse of the code is in c++ parallelized via OpenMP (btw thanks Rcpp for simplifying coders' life!).

 - This package is the follow up of the (now deprecated) package `FENmlm` which performed fixed-effects estimations but for only four likelihood families. Package `fixest` completely supersedes `FENmlm` by extending the method to regular OLS and all GLM families, and adding new utility functions. Further, the design of the functions has been completely overhauled and extended towards much more user-friendliness. Massive effort has been put into providing a set of informative error messages to the user for quick debugging of her workflow (e.g. one of the functions contains over 100 different errors messages).


