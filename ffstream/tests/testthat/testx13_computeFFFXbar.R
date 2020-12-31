context("Test 13: direct vector computation of FFF xbar")

test_that("make computeFFFMean produce errors", {
        #check for non-numerics
        x1 <- c(1, "a", 3)
        #interestingly, it just matches a regular expression, 
        #so does not need to be exact
        expect_error(computeFFFMean(x1), "x contains non-numeric")

        #check for AN
        x2 <- c(1, NA, 3)
        expect_error(computeFFFMean(x2), "contains NA")

        #check for NaNs
        x3 <- c(1, NaN, 3)
        expect_error(computeFFFMean(x3), "contains NA")

        #check for Inf
        x4 <- c(1, Inf, 4)
        expect_error(computeFFFMean(x4), "contains non-finite")

        #now check lambda
        x5 <- c(1, 2, 4)
        expect_error(computeFFFMean(x5, lambda=NA), "not a finite")
        expect_error(computeFFFMean(x5, lambda=NULL), "is NULL")
        expect_error(computeFFFMean(x5, lambda=0), NA) #success
        expect_error(computeFFFMean(x5, lambda=1), NA) #success
        expect_error(computeFFFMean(x5, lambda=0.9), NA) #success
        expect_error(computeFFFMean(x5, lambda=-0.1), "not in range") 
        expect_error(computeFFFMean(x5, lambda=1.1), "not in range") 
        
        })


test_that("check computeFFFMean works", {

        #by hand
        x1 <- c(1, 2, 3)
        lambda1 <- 0.9
        m1 <- (0.9)^2 * 1 + (0.9) * 2 + 3
        w1 <- (0.9)^2 + 0.9 + 1
        expect_equal(computeFFFMean(x1, lambda1), m1/w1)

        #by hand
        x2 <- c(5)
        lambda2 <- 0.9
        expect_equal(computeFFFMean(x2, lambda2), 5)
        
        })


test_that("long check computeFFFMean works", {

        #by hand
        set.seed(1)
        x1 <- rnorm(100)
        lambda1 <- 0.9

        #compute with R
        m <- 0
        w <- 0
        for (i in seq_len(length(x1))){
            m <- update_m(m, lambda1, x1[i])
            w <- update_w(w, lambda1)
        }
        fffxbar <- get_xbar(m, w)

        expect_equal(computeFFFMean(x1, lambda1), fffxbar)
        
        })

