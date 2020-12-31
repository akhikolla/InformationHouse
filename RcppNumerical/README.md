## Rcpp Integration for Numerical Computing Libraries <img src="https://statr.me/images/sticker-rcppnumerical.png" alt="RcppNumerical" height="150px" align="right" />

- [Introduction](#introduction)
- [Numerical Integration](#numerical-integration)
  - [One-dimensional](#one-dimensional)
  - [Multi-dimensional](#multi-dimensional)
- [Numerical Optimization](#numerical-optimization)
- [A More Interesting Example](#a-more-interesting-example)

### Introduction

[Rcpp](https://CRAN.R-project.org/package=Rcpp) is a
powerful tool to write fast C++ code to speed up R programs. However,
it is not easy, or at least not straightforward, to compute numerical
integration or do optimization using pure C++ code inside Rcpp.

**RcppNumerical** integrates a number of open source numerical computing
libraries into Rcpp, so that users can call convenient functions to
accomplish such tasks.

- To use **RcppNumerical** with `Rcpp::sourceCpp()`, add
```cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
```
in the C++ source file.
- To use **RcppNumerical** in your package, add `Imports: RcppNumerical`
and `LinkingTo: Rcpp, RcppEigen, RcppNumerical` to the `DESCRIPTION` file,
and `import(RcppNumerical)` to the `NAMESPACE` file.

### Numerical Integration

#### One-dimensional

The one-dimensional numerical integration code contained in **RcppNumerical**
is based on the [NumericalIntegration](https://github.com/tbs1980/NumericalIntegration)
library developed by [Sreekumar Thaithara Balan](https://github.com/tbs1980),
[Mark Sauder](https://github.com/mcsauder), and Matt Beall.

To compute integration of a function, first define a functor derived from
the `Func` class (under the namespace `Numer`):

```cpp
class Func
{
public:
    virtual double operator()(const double& x) const = 0;
    virtual void eval(double* x, const int n) const
    {
        for(int i = 0; i < n; i++)
            x[i] = this->operator()(x[i]);
    }
    
    virtual ~Func() {}
};
```

The first function evaluates one point at a time, and the second version
overwrites each point in the array by the corresponding function values.
Only the second function will be used by the integration code, but usually it
is easier to implement the first one.

**RcppNumerical** provides a wrapper function for the **NumericalIntegration**
library with the following interface:

```cpp
inline double integrate(
    const Func& f, const double& lower, const double& upper,
    double& err_est, int& err_code,
    const int subdiv = 100, const double& eps_abs = 1e-8, const double& eps_rel = 1e-6,
    const Integrator<double>::QuadratureRule rule = Integrator<double>::GaussKronrod41
)
```

- `f`: The functor of integrand.
- `lower`, `upper`: Limits of integral.
- `err_est`: Estimate of the error (output).
- `err_code`: Error code (output). See `inst/include/integration/Integrator.h`
[Line 676-704](https://github.com/yixuan/RcppNumerical/blob/master/inst/include/integration/Integrator.h#L676).
- `subdiv`: Maximum number of subintervals.
- `eps_abs`, `eps_rel`: Absolute and relative tolerance.
- `rule`: Integration rule. Possible values are
`GaussKronrod{15, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 201}`. Rules with
larger values have better accuracy, but may involve more function calls.
- Return value: The final estimate of the integral.

See a full example below, which can be compiled using the `Rcpp::sourceCpp`
function in Rcpp.

```cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// P(0.3 < X < 0.8), X ~ Beta(a, b)
class BetaPDF: public Func
{
private:
    double a;
    double b;
public:
    BetaPDF(double a_, double b_) : a(a_), b(b_) {}

    double operator()(const double& x) const
    {
        return R::dbeta(x, a, b, 0);
    }
};

// [[Rcpp::export]]
Rcpp::List integrate_test()
{
    const double a = 3, b = 10;
    const double lower = 0.3, upper = 0.8;
    const double true_val = R::pbeta(upper, a, b, 1, 0) -
                            R::pbeta(lower, a, b, 1, 0);

    BetaPDF f(a, b);
    double err_est;
    int err_code;
    const double res = integrate(f, lower, upper, err_est, err_code);
    return Rcpp::List::create(
        Rcpp::Named("true") = true_val,
        Rcpp::Named("approximate") = res,
        Rcpp::Named("error_estimate") = err_est,
        Rcpp::Named("error_code") = err_code
    );
}
```

Runing the `integrate_test()` function in R gives

```r
integrate_test()
## $true
## [1] 0.2528108
##
## $approximate
## [1] 0.2528108
##
## $error_estimate
## [1] 2.806764e-15
##
## $error_code
## [1] 0
```

#### Multi-dimensional

Multi-dimensional integration in **RcppNumerical** is done by the
[Cuba](http://www.feynarts.de/cuba/) library developed by
[Thomas Hahn](http://wwwth.mpp.mpg.de/members/hahn/).

To calculate the integration of a multivariate function, one needs to define
a functor that inherits from the `MFunc` class:

```cpp
class MFunc
{
public:
    virtual double operator()(Constvec& x) = 0;
    
    virtual ~MFunc() {}
};
```

Here `Constvec` represents a read-only vector with the definition

```cpp
// Constant reference to a vector
typedef const Eigen::Ref<const Eigen::VectorXd> Constvec;
```

(Basically you can treat `Constvec` as a `const Eigen::VectorXd`. Using
`Eigen::Ref` is mainly to avoid memory copy. See the explanation
[here](http://eigen.tuxfamily.org/dox/classEigen_1_1Ref.html).)

The function provided by **RcppNumerical** for multi-dimensional
integration is

```cpp
inline double integrate(
    MFunc& f, Constvec& lower, Constvec& upper,
    double& err_est, int& err_code,
    const int maxeval = 1000,
    const double& eps_abs = 1e-6, const double& eps_rel = 1e-6
)
```

- `f`: The functor of integrand.
- `lower`, `upper`: Limits of integral. Both are vectors of the same
dimension of `f`.
- `err_est`: Estimate of the error (output).
- `err_code`: Error code (output). Non-zero values indicate failure of
convergence.
- `maxeval`: Maximum number of function evaluations.
- `eps_abs`, `eps_rel`: Absolute and relative tolerance.
- Return value: The final estimate of the integral.

See the example below:

```cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// P(a1 < X1 < b1, a2 < X2 < b2), (X1, X2) ~ N([0], [1   rho])
//                                            ([0], [rho   1])
class BiNormal: public MFunc
{
private:
    const double rho;
    double const1;  // 2 * (1 - rho^2)
    double const2;  // 1 / (2 * PI) / sqrt(1 - rho^2)
public:
    BiNormal(const double& rho_) : rho(rho_)
    {
        const1 = 2.0 * (1.0 - rho * rho);
        const2 = 1.0 / (2 * M_PI) / std::sqrt(1.0 - rho * rho);
    }

    // PDF of bivariate normal
    double operator()(Constvec& x)
    {
        double z = x[0] * x[0] - 2 * rho * x[0] * x[1] + x[1] * x[1];
        return const2 * std::exp(-z / const1);
    }
};

// [[Rcpp::export]]
Rcpp::List integrate_test2()
{
    BiNormal f(0.5);  // rho = 0.5
    Eigen::VectorXd lower(2);
    lower << -1, -1;
    Eigen::VectorXd upper(2);
    upper << 1, 1;
    double err_est;
    int err_code;
    const double res = integrate(f, lower, upper, err_est, err_code);
    return Rcpp::List::create(
        Rcpp::Named("approximate") = res,
        Rcpp::Named("error_estimate") = err_est,
        Rcpp::Named("error_code") = err_code
    );
}
```

We can test the result in R:

```r
library(mvtnorm)
trueval = pmvnorm(c(-1, -1), c(1, 1), sigma = matrix(c(1, 0.5, 0.5, 1), 2))
integrate_test2()
## $approximate
## [1] 0.4979718
##
## $error_estimate
## [1] 4.612333e-09
##
## $error_code
## [1] 0
trueval - integrate_test2()$approximate
## [1] 2.893336e-11
```

### Numerical Optimization

Currently **RcppNumerical** contains the L-BFGS algorithm for unconstrained
minimization problems based on the
[LBFGS++](https://github.com/yixuan/LBFGSpp) library.

Again, one needs to first define a functor to represent the multivariate
function to be minimized.

```cpp
class MFuncGrad
{
public:
    virtual double f_grad(Constvec& x, Refvec grad) = 0;
    
    virtual ~MFuncGrad() {}
};
```

Same as the case in multi-dimensional integration, `Constvec` represents a
read-only vector and `Refvec` a writable vector. Their definitions are

```cpp
// Reference to a vector
typedef Eigen::Ref<Eigen::VectorXd>             Refvec;
typedef const Eigen::Ref<const Eigen::VectorXd> Constvec;
```

The `f_grad()` member function returns the function value on vector `x`,
and overwrites `grad` by the gradient.

The wrapper function for **LBFGS++** is

```cpp
inline int optim_lbfgs(
    MFuncGrad& f, Refvec x, double& fx_opt,
    const int maxit = 300, const double& eps_f = 1e-6, const double& eps_g = 1e-5
)
```

- `f`: The function to be minimized.
- `x`: In: the initial guess. Out: best value of variables found.
- `fx_opt`: Out: Function value on the output `x`.
- `maxit`: Maximum number of iterations.
- `eps_f`: Algorithm stops if `|f_{k+1} - f_k| < eps_f * |f_k|`.
- `eps_g`: Algorithm stops if `||g|| < eps_g * max(1, ||x||)`.
- Return value: Error code. Negative values indicate errors.

Below is an example that illustrates the optimization of the Rosenbrock function
`f(x1, x2) = 100 * (x2 - x1^2)^2 + (1 - x1)^2`:

```cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>

using namespace Numer;

// f = 100 * (x2 - x1^2)^2 + (1 - x1)^2
// True minimum: x1 = x2 = 1
class Rosenbrock: public MFuncGrad
{
public:
    double f_grad(Constvec& x, Refvec grad)
    {
        double t1 = x[1] - x[0] * x[0];
        double t2 = 1 - x[0];
        grad[0] = -400 * x[0] * t1 - 2 * t2;
        grad[1] = 200 * t1;
        return 100 * t1 * t1 + t2 * t2;
    }
};

// [[Rcpp::export]]
Rcpp::List optim_test()
{
    Eigen::VectorXd x(2);
    x[0] = -1.2;
    x[1] = 1;
    double fopt;
    Rosenbrock f;
    int res = optim_lbfgs(f, x, fopt);
    return Rcpp::List::create(
        Rcpp::Named("xopt") = x,
        Rcpp::Named("fopt") = fopt,
        Rcpp::Named("status") = res
    );
}
```

Calling the generated R function `optim_test()` gives

```r
optim_test()
## $xopt
## [1] 1 1
##
## $fopt
## [1] 3.12499e-15
##
## $status
## [1] 0
```

### A More Practical Example

It may be more meaningful to look at a real application of the **RcppNumerical**
package. Below is an example to fit logistic regression using the L-BFGS
algorithm. It also demonstrates the performance of the library.

```cpp
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>

using namespace Numer;

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

class LogisticReg: public MFuncGrad
{
private:
    const MapMat X;
    const MapVec Y;
public:
    LogisticReg(const MapMat x_, const MapVec y_) : X(x_), Y(y_) {}

    double f_grad(Constvec& beta, Refvec grad)
    {
        // Negative log likelihood
        //   sum(log(1 + exp(X * beta))) - y' * X * beta

        Eigen::VectorXd xbeta = X * beta;
        const double yxbeta = Y.dot(xbeta);
        // X * beta => exp(X * beta)
        xbeta = xbeta.array().exp();
        const double f = (xbeta.array() + 1.0).log().sum() - yxbeta;

        // Gradient
        //   X' * (p - y), p = exp(X * beta) / (1 + exp(X * beta))

        // exp(X * beta) => p
        xbeta.array() /= (xbeta.array() + 1.0);
        grad.noalias() = X.transpose() * (xbeta - Y);

        return f;
    }
};

// [[Rcpp::export]]
Rcpp::NumericVector logistic_reg(Rcpp::NumericMatrix x, Rcpp::NumericVector y)
{
    const MapMat xx = Rcpp::as<MapMat>(x);
    const MapVec yy = Rcpp::as<MapVec>(y);
    // Negative log likelihood
    LogisticReg nll(xx, yy);
    // Initial guess
    Eigen::VectorXd beta(xx.cols());
    beta.setZero();

    double fopt;
    int status = optim_lbfgs(nll, beta, fopt);
    if(status < 0)
        Rcpp::stop("fail to converge");

    return Rcpp::wrap(beta);
}
```

Here is the R code to test the function:

```r
set.seed(123)
n = 5000
p = 100
x = matrix(rnorm(n * p), n)
beta = runif(p)
xb = c(x %*% beta)
p = exp(xb) / (1 + exp(xb))
y = rbinom(n, 1, p)

system.time(res1 <- glm.fit(x, y, family = binomial())$coefficients)
##   user  system elapsed
##  0.229   0.006   0.234
system.time(res2 <- logistic_reg(x, y))
##   user  system elapsed
##  0.005   0.000   0.006
max(abs(res1 - res2))
## [1] 0.0001873564
```

It is much faster than the standard `glm.fit()` function in R! (Although
`glm.fit()` calculates some other quantities besides beta.)

**RcppNumerical** also provides the `fastLR()` function to run fast logistic
regression, which is a modified and more stable version of the code above.

```r
system.time(res3 <- fastLR(x, y)$coefficients)
##   user  system elapsed
##  0.007   0.001   0.008
max(abs(res1 - res3))
## [1] 7.066969e-06
```
