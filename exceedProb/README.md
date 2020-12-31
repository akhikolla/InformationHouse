# exceedProb
This R package computes exceedance probabilities and associated confidence intervals. It currently only supports general linear models.

## Installation

```{r}
install.packages("exceedProb")
```

## Examples

```{r}
library(exceedProb)

# Sample mean -----------------------------------------------------------------
n <- 100
x <- rnorm(n = n)

theta_hat <- mean(x)
sd_hat <- sd(x)

cutoff <- seq(from = theta_hat - 0.5, to = theta_hat + 0.5, by = 0.1)

exceedProb(cutoff = cutoff, 
           theta_hat = theta_hat, 
           sd_hat = sd_hat, 
           alpha = 0.05, 
           d = 1,
           n = n,
           m = n)

# Linear regression -----------------------------------------------------------
n <- 100
beta <- c(1, 2)
x <-runif(n = n, min = 0, max = 10)
y <- rnorm(n = n, mean = cbind(1, x) %*% beta, sd = 1)

j <- 2
fit <- lm(y ~ x)
theta_hat <- coef(fit)[j]
sd_hat <- sqrt(n * vcov(fit)[j, j])

cutoff <- seq(from = theta_hat - 0.5, to = theta_hat + 0.5, by = 0.1)

exceedProb(cutoff = cutoff, 
           theta_hat = theta_hat, 
           sd_hat = sd_hat, 
           alpha = 0.05, 
           d = length(beta),
           n = n,
           m = n)
```

## References

Segal, B. D. (submitted). Towards replicability with confidence intervals for the exceedance probability.
