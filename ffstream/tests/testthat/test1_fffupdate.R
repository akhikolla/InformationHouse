context("Test  1: computation of FFF w, xbar, u and v and s2")

#will be used in the tests

library(ffstream)
library(Rcpp)

test_that("w is computed properly", {
        fff1 <- initFFFMean(lambda=0.9)

        x <- c(1, 2, 3)
        expect_equal(fff1$w, 0)

        fff1$update(x[1])
        expect_equal(fff1$w, 1)

        fff1$update(x[2])
        expect_equal(fff1$w, 1.9)

        fff1$update(x[3])
        expect_equal(fff1$w, 2.71)
        
        })


test_that("xbar is computed properly", {
#        fff1 <- new (FFF, 0.9)
        fff1 <- initFFFMean(lambda=0.9)
        x <- c(1, 2, 3)
        
        expect_equal(fff1$xbar, 0)

        fff1$update(x[1])
        expect_equal(fff1$xbar, 1)
        
        fff1$update(x[2])
        expect_equal(fff1$xbar, 2.9/1.9)

        fff1$update(x[3])
        expect_equal(fff1$xbar, 5.61/2.71)
        
        })



test_that("u is computed properly", {
        fff1 <- new (FFF, 0.9)
#        fff1 <- initFFFMean(lambda=0.9)
        x <- c(1, 2, 3)
        
        expect_equal(fff1$u, 0)

        fff1$update(x[1])
        expect_equal(fff1$u, 1)

        fff1$update(x[2])
        expect_equal(fff1$u, 0.5013850416)
        
        fff1$update(x[3])
        expect_equal(fff1$u,  0.3357933579)
        
        })


test_that("v is computed properly", {
#        fff1 <- new (FFF, 0.9)
        fff1 <- initFFFMean(lambda=0.9)
        x <- c(1, 2, 3)

        expect_equal(fff1$v, 0)

        fff1$update(x[1])
        expect_equal(fff1$v, 0)

        fff1$update(x[2])
        expect_equal(fff1$v, 0.947368421)
        #expect_equal(fff1$xbar, 3.9/1.9)

        fff1$update(x[3])
        #crazy! how did this work out...?
        expect_equal(fff1$v, 1.8)

        })


test_that("s2 is computed properly", {
#        fff1 <- new (FFF, 0.9)
        fff1 <- initFFFMean(lambda=0.9)
        x <- c(1, 2, 3)

        expect_equal(fff1$s2, 0)

        fff1$update(x[1])
        expect_equal(fff1$s2, 0)

        fff1$update(x[2])
        expect_equal(fff1$s2, 0.50)

        fff1$update(x[3])
        expect_equal(fff1$s2, 0.9981549815)

        })



test_that("non-FFF s2 is computed properly", {
        fff1 <- new (FFF, 1)
#        fff1 <- initFFFMean(lambda=1)
        set.seed(1)
        x <- rnorm(10)

        for (i in 1:length(x)){
            fff1$update(x[i])
        }

       s2 <- var(x)

        expect_equal(fff1$s2, s2)

        })

test_that("processVector works", {
        set.seed(2)
        x <- rnorm(10)
        #two forgetting factors
        fff1 <- initFFFMean(lambda=0.9)
        fff2 <- initFFFMean(lambda=0.9)

        #sequential
        for (i in seq_len(length(x))){
            obs <- x[i]
            fff1$update(obs)
        }

        #update by vector
        fff2$processVector(x)

        expect_equal(fff1$xbar, fff2$xbar)
        expect_equal(fff1$s2, fff2$s2)

        })
