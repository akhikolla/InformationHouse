# PAutilities 1.0.1

## What's new:

* Some minor modifications have been made to clean up internal code.

* The transition pairing method unit tests have been temporarily suspended.

# PAutilities 1.0.0

## What's new:

* Added a function (weight_status) to classify body mass index for adults
    using CDC cutoffs. (This can also be used as a wrapper for the existing
    function that classifies youth BMI based on CDC percentiles.)
    
* Added Rcpp functions related to rolling windows (`get_indices` to extract
    indices for rolling windows; `rolling_groups` to extract rolling
    subsequences of raw data).
    
* Removed dependency on orphaned `clues` package (resulting in changes to
    unit test cache)

# PAutilities 0.3.1

## What's new:

* Added a vignette for the Transition Pairing Method

# PAutilities 0.3.0

## What's new:

* Added paired equivalence testing functionality, including plot method
* Added a function to classify activity intensity from metabolic equivalents and
    (if available) posture
* Added a wrapper for `base::rle` that gives a data frame with original `rle`
    information, plus start/stop indices for each run
* Added various performance indicators to `summary.transition`
* Converted `summary.transition` output to an S4 framework, and added
    addition/subraction methods
* Added spurious curve generation functions
* Adjusted Transition Pairing Method to make it more robust (e.g. to deal with
    missing values)

# PAutilities 0.2.0

## What's new:

* Removed unused dependency on package `AGread`
* Added support for calculation of RMR via sliding window analysis
* Added process management utility for easy messaging and timekeeping
* Added tools to compare objects more specifically than via `all.equal`
* Added `na.rm = TRUE` to the mean bias line on Bland-Altman plots
* Added shape control as an option for Bland-Altman plots
* Cleaned up TPM summary method
* Added rejection of non-consecutive pairings to the TPM

# PAutilities 0.1.2

## Summary

This is a second resubmission of version 0.1.0. The following changes
    have been made:
    
* Changed title to eliminate redundancy
* Expanded the `description` field in DESCRIPTION to provide more
    detail about what the package does.
* Swapped donttest{} in for dontrun{} in documentation where appropriate
* Added detail to package documentation

# PAutilities 0.1.1

## Summary

This is a resubmission of version 0.1.0. The extraneous LICENSE file
    has been removed.

# PAutilities 0.1.0

## Summary

This is the intial release of PAutilities. Current features include:

* Bland-Altman plots
* Bouted moderate-to-vigorous physical activity analysis
* Formatted descriptive statistics
* Demographic calculations (age, BMI-for-age percentile)
* Metabolic calculations (energy expenditure conversions,
    basal metabolic rate predictions)
* Analysis of bout detection algorithm performance
