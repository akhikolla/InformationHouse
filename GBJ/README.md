
<!-- README.md is generated from README.Rmd. Please edit that file -->
What is GBJ?
------------

The Generalized Berk-Jones statistic was developed to perform set-based inference in genetic association studies. It is an alternative to tests such as the Sequence Kernel Association Test (SKAT), Generalized Higher Criticism (GHC), and Minimum p-value (minP).

Why use GBJ?
------------

GBJ is a generalization of the Berk-Jones (BJ) statistic, which offers - in a certain sense - asymptotic power guarantees for detection of rare and weak signals. GBJ modifies BJ to account for correlation between factors in a set. GBJ has been demonstrated to outperform other tests when signals are moderately sparse (more precisely, when the number of signals is between *d*<sup>1/4</sup> and *d*<sup>1/2</sup>, where *d* is the number of factors in the set).

Other advantages include:
1. Analytic p-value calculation (no need for permutation inference).
2. Can be applied to individual-level genotype data or GWAS summary statistics.
3. No tuning parameters. Accepts standard inputs (similar to glm() function).

Example
-------

We show a simple example for testing the association between a set of 50 SNPs (which could be, for example, from the same gene or pathway) and a binary outcome.

``` r
library(GBJ)
set.seed(1000)

# Case-control study, 1000 subjects
cancer_status <- c(rep(1,500), rep(0,500))

# We have 50 SNPs each with minor allele frequency of 0.3 in this example
genotype_data <- matrix(data=rbinom(n=1000*50, size=2, prob=0.3), nrow=1000)
age <- round( runif(n=1000, min=30, max=80) )
gender <- rbinom(n=1000, size=1, prob=0.5)     

# Fit the null model, calculate marginal score statistics for each SNP
# (asymptotically equivalent to those calculated by, for example, PLINK)
null_mod <- glm(cancer_status~age+gender, family=binomial(link="logit"))
log_reg_stats <- calc_score_stats(null_model=null_mod, factor_matrix=genotype_data, link_function="logit")

# Run the test
GBJ(test_stats=log_reg_stats$test_stats, cor_mat=log_reg_stats$cor_mat)
#> $GBJ
#> [1] 1.43984
#> 
#> $GBJ_pvalue
#> [1] 0.330911
#> 
#> $err_code
#> [1] 0
```

What else is in here?
---------------------

We may not have convinced you that GBJ is the best option for your application. If that is the case, then you may still be interested in trying the Berk-Jones (BJ), Generalized Higher Criticism (GHC), Higher Criticism (HC), or Minimum p-value (minP) tests, which can be run with the same inputs, i.e. GHC(test\_stats=score\_stats, cor\_mat=cor\_Z) to run the GHC. We also have developed an omnibus test which information from multiple different methods. Please see the vignette for more details.
