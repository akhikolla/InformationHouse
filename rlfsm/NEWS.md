# rlfsm 1.0.1

* References and DESCRIPTION were updated

* a_tilda was rewritten in C++. The old function was left for compatibility and tests under the name a_tilda_R.

# rlfsm 1.0.0

* path function deals with index shift differently now. On top of that, the default value for seed is NULL.

* More autotests were added

* CLT was deprecated and replaced by MCestimLFSM

* Some bugs associated with plot functions were fixed

# rlfsm 0.3.1

* Automated testing was implemented via testthat package

* Documentation slightly improved

# rlfsm 0.2.1

* Various optimizations were implemented in functions related to increment().

* Parameter estimator from Ljungdahl and Podolskij, 2019 was added.

* CLT() was largely amended.

* New plotting function for rlfsm sample paths was added.

* Slight changes to documentation were introduced.

# rlfsm 0.2.0

* CLT function was equipped with a faster version of the sample path generator which is unavailable for users. User version (function path) remains intact.

* CLT documentation was amended.

* The developer repository was moved to Gitlab.com.
