---
title: "NEWS"
author: "Erik Peter"
date: "Saturday, 24th March, 2019"
output: html_document
---

## fbroc 0.4.1

### Other Changes

Update for compability with newer R Versions like switch to roxygen2.

## fbroc 0.4.0

### New features

* Partial AUCs over both TPR and FPR ranges can be calculated
* You can now adjust text size for plots
* In the ROC plot the (partial) AUC can now optionally be shown instead of confidence regions

### Other Changes

* The location of the text showing the performance in the ROC plot has been shifted downwards and
to the left

## fbroc 0.3.1

### Bugfixed

* Fixed broken plot function for single ROC curves when showing a metric

## fbroc 0.3.0

### New features

* Allows the comparison of two paired ROC curves

### Bugfixes

* Bad defaults caused plotting to fail with a large number of negative samples

### Other changes

* perf.roc is now deprecated. Use the new S3 generic perf instead
* conf.roc is now deprecated. Use the new S3 generic conf instead

## fbroc 0.2.1

### Bugfixes

* fixed a off-by-one pointer error

## fbroc 0.2.0

### New features

* Allow uncached bootstrap of the ROC curve to avoid memory issues, this now the new default
* New performance metrices: TPR at fixed FPR and FPR at fixed TPR

### Other changes

* Stand-alone function to find thresholds `calculate.thresholds` was removed. To calculate thresholds
please call boot.roc and look at list item `roc` of the outpot
* Smarter default for the number of steps in `conf.roc`
* Smarter default for the number of bins in `plot.fbroc.perf`

### Internal Changes

* Completely refactored C++ code for improved maintability

### Bugfixes

* Function `boot.tpr.at.fpr` now works properly 
* For duplicated predictions not all relevant thresholds were found reliably, this was fixed

## fbroc 0.1.0

* Initial release