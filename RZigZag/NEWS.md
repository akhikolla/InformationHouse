
# RZigZag 0.2.0
* general refactoring which is unfortunately not downwards compatible
* R helper functions to obtain ESS and samples from skeleton
* transposed dataX matrix in logistic regression: now the different rows of dataX represent different observations 
* allow user to provide starting direction and custom reference points for control variates (ZigZagLogistic)
* BPS and Zig-Zag routinesStudent distributions

# RZigZag 0.1.6
* Allows for fixed time horizon (in additional to fixed number of iterations)
* Preprocessing of input to logistic regression disabled. Adding intercept to design matrix and possible recentering is up to user.

# RZigZag 0.1.5
* Implements Bouncy Particle Sampler for Gaussian target
* Allows to specify starting value for ZigZagGaussian

# RZigZag 0.1.4
* Changed from n_epochs argument to n_iter, representing number of iterations, as this is more flexible. This update is not downwards compatible, but if you want to run ZigZagLogistic when subsampling (with or without control variates) with n_epochs gradient evaluations, simply pick n_iter = n_epochs * (number of observations).
* Resolved ambiguities in overloaded functions which previously resulted in a compile error on Solaris.

# RZigZag 0.1.3
* Corrected arXiv reference

# RZigZag 0.1.1
* This is the first version, which introduces the R functions `ZigZagLogistic` and `ZigZagGaussian`.
