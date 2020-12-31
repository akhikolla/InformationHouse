context("Test  7: initialisation of EWMA")

library(Rcpp)
test_that("init EWMA (r fails)", {

        expect_error(ewmacd1 <- initEWMAMeanCD(r="A", L=3.00, BL=50), " r ")
        expect_error(ewmacd1 <- initEWMAMeanCD(r="A", L=3.00, BL=50), " not a finite ")

        expect_error(ewmacd1 <- initEWMAMeanCD(r=NaN, L=3.00, BL=50), " r ")
        expect_error(ewmacd1 <- initEWMAMeanCD(r=NaN, L=3.00, BL=50), " not a finite ")
        
        expect_error(ewmacd1 <- initEWMAMeanCD(r=Na, L=3.00, BL=50), "not found")

        expect_error(ewmacd1 <- initEWMAMeanCD(r=NA, L=3.00, BL=50), " r ")
        expect_error(ewmacd1 <- initEWMAMeanCD(r=NA, L=3.00, BL=50), " not a finite ")

        expect_error(ewmacd1 <- initEWMAMeanCD(r=-0.01, L=3.00, BL=50), " r ")
        expect_error(ewmacd1 <- initEWMAMeanCD(r=-0.01, L=3.00, BL=50), " not in range")

        expect_error(ewmacd1 <- initEWMAMeanCD(r=1.01, L=3.00, BL=50), " r ")
        expect_error(ewmacd1 <- initEWMAMeanCD(r=1.01, L=3.00, BL=50), " not in range")

        })


test_that("init EWMA (L fails)", {

        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L="A", BL=50), " L ")
        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L="A", BL=50), "not a finite")

        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L=NaN, BL=50), " L ")
        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L=NaN, BL=50), "not a finite")

        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L=Na, BL=50), "not found")

        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L=NA, BL=50), " L ")
        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L=NA, BL=50), "not a finite")

        })


test_that("EWMA init (BL)", {
        #should FAIL - BL not in range 
        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L=3.00, BL=-1), "BL")
        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L=3.00, BL=-1), "not greater than")
        })



test_that("EWMA init (passes)", {
        #should pass
        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L=3.00, BL=50), NA)

        #should pass
        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, L=3.00), NA)

        #should pass
        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25, BL=50), NA)

        #should pass
        expect_error(ewmacd1 <- initEWMAMeanCD(L=3.00, BL=50), NA)

        #should pass
        expect_error(ewmacd1 <- initEWMAMeanCD(r=0.25), NA)

        #should pass
        expect_error(ewmacd1 <- initEWMAMeanCD(L=3.00), NA)

        #should pass
        expect_error(ewmacd1 <- initEWMAMeanCD(BL=50), NA)

        #should pass
        expect_error(ewmacd1 <- initEWMAMeanCD(), NA)

        })
