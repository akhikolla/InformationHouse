# Change log

**Version 0.0.1**

* First submission to CRAN.
* Remove plotting functions. Users can create their own plots using examples in the vignette.
* Make the API of PCR and LDS reconstructions consistent.

**Version 0.0.0.9003**

* Added ensemble model function

**Version 0.0.0.9002**

* Added utility function `water_to_calendar_year()`.
* Added example dataset `P1monthly`, which is monthly streamflow at P1 following the Thai water year.
* Modified cross-validation function `cvLDS()` to use one core less than the total number of cores.

**Version 0.0.0.9001**

Added experimental features: learing with genetic algorithm and with L-BFGS-B. These features are not fully tested and should not be used. Users are strongly advised to only use the EM method.
