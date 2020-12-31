# TauStar Package

## Purpose

This package allows you to efficiently compute, and perform tests of
independence with, the U/V-statistic corresponding to the tau* coefficient
described in the paper:

Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based on a
sign covariance related to Kendall's tau. Bernoulli 20 (2014), no. 2, 1006–1028.

The tau* statistic has the special property that it is 0 if and only if the
bivariate distribution it is computed upon is independent (under some weak
conditions on the bivariate distribution) and is positive otherwise. Since t*, 
the U-statistic corresponding to tau*, is an unbiased estimator of tau* this 
gives a consistent test of independence. Computing t* naively results an 
algorithm that takes O(n^4) time where n is the sample size. Luckily it is 
possible to compute t* much faster (in O(n^2) time) using the algorithm 
described in:

Heller, Yair and Heller, Ruth. "Computing the Bergsma Dassios sign-covariance."
arXiv preprint arXiv:1605.08732 (2016).

building off of the O(n^2*log(n)) algorithm of:

Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the
Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).

This fast algorithm is implemented in this package. Moreover, the package also
uses the results of Nandy, Weihs, and Drton (2016) to allow the use of t* in
performing tests of independence. In particular, we provide the function
tauStarTest which automates tests of independence using the asymptotic null
distribution of t*.

## Example

A simple example of computing t* on a independent bivariate normal distribution
follows:

```
> set.seed(2342)
> x = rnorm(1000)
> y = rnorm(1000)
> tStar(x, y)
[1] 0.0003637266
```

Similarly, we may obtain the asymptotic p-value corresponding to a test of
independence as follows:

```
> set.seed(2341)
> x = rnorm(1000)
> y = rnorm(1000)
> tauStarTest(x, y)$pVal
[1] 0.5692797
```

## Where to go

The main functionality of this package is currently included in the functions
`tStar` (which computes the t* statistic on two input vectors) and `tauStarTest`
(which performs tests of independence using t*). One may also be interested in
the functions

1. `pHoeffInd`, `dHoeffInd`, `rHoeffInd`, `qHoeffInd`
2. `pDisHoeffInd`, `dDisHoeffInd`, `rDisHoeffInd`, `qDisHoeffInd`
3. `pMixHoeffInd`, `dMixHoeffInd`, `rMixHoeffInd`, `qMixHoeffInd`

which compute distribution functions, densities, random samples, and quantiles
for the asymptotic distribution of t* in different cases.