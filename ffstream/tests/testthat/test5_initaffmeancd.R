context("Test  5: initialisation of AFF mean cd")

#will be used in the tests

test_that("AFF change detector is initialised correctly (alpha)", {

        #----------------------------------#
        ###TESTS THAT FAIL 
        #----------------------------------#
        #initAFFMeanCD <- function(alpha, eta=0.01, BL=50){

        #should FAIL - alpha not numeric
        expect_error(fffcd1 <- initAFFMeanCD("A", 0.9), "alpha")
        expect_error(fffcd1 <- initAFFMeanCD("A", 0.9), "not a finite")

        #should FAIL - alpha not finite 
        expect_error(fffcd1 <- initAFFMeanCD(NaN, 0.9), "alpha")
        expect_error(fffcd1 <- initAFFMeanCD(NaN, 0.9), "not a finite")
        
        #should FAIL - alpha not finite 
        expect_error(fffcd1 <- initAFFMeanCD(Na, 0.9), "not found")
        
        #should FAIL - alpha not finite 
        expect_error(fffcd1 <- initAFFMeanCD(NA, 0.9), "alpha")
        expect_error(fffcd1 <- initAFFMeanCD(NA, 0.9), "not a finite")

        #should FAIL - alpha not in range 
        expect_error(fffcd1 <- initAFFMeanCD(-1, 0.9), "alpha")
        expect_error(fffcd1 <- initAFFMeanCD(-1, 0.9), "not in range")
        
        #should FAIL - alpha not in range 
        expect_error(fffcd1 <- initAFFMeanCD(1.1, 0.9), "alpha")
        expect_error(fffcd1 <- initAFFMeanCD(1.1, 0.9), "not in range")

        #should FAIL - alpha not in range 
        expect_error(fffcd1 <- initAFFMeanCD(1, 0), "alpha")
        expect_error(fffcd1 <- initAFFMeanCD(1, 0), "not in range")

        #should FAIL - alpha not in range 
        expect_error(fffcd1 <- initAFFMeanCD(1, 1), "alpha")
        expect_error(fffcd1 <- initAFFMeanCD(1, 1), "not in range")

        })


test_that("AFF change detector is initialised correctly (eta)", {
        #should FAIL - lambda not in range 
        expect_error(fffcd1 <- initAFFMeanCD(eta=1.1, alpha=0.01), "eta")
        expect_error(fffcd1 <- initAFFMeanCD(eta=1.1, alpha=0.01), "not in range")

        })


test_that("AFF change detector is initialised correctly (BL)", {
        #should FAIL - BL not in range 
        expect_error(affcd1 <- initAFFMeanCD(0.01, 0.01, BL=-1), "BL")
        expect_error(affcd1 <- initAFFMeanCD(0.01, 0.01, BL=-1), "not greater than")


        #SHOULD FAIL - BL not in range
        expect_error(affcd1 <- initAFFMeanCD(0.01, 0.01, BL=1), "BL")
        expect_error(affcd1 <- initAFFMeanCD(0.01, 0.01, BL=1), "not greater than")

        #no more lambdaBL
        #should FAIL - lambdaBL not numeric 
        #expect_error(fffcd1 <- initAFFMeanCD(1, 0.9, BL=50, lambdaBL="A"))

        #should FAIL - lambdaBL not finite 
        #expect_error(fffcd1 <- initAFFMeanCD(1, 0.9, BL=50, lambdaBL=NA))

        #should FAIL - lambdaBL not numeric 
        #expect_error(fffcd1 <- initAFFMeanCD(1, 0.9, BL=50, lambdaBL=NaN))

        #should FAIL - lambdaBL not numeric 
        #expect_error(fffcd1 <- initAFFMeanCD(1, 0.9, BL=50, lambdaBL=Na))

        #should FAIL - lambdaBL not in range 
        #expect_error(fffcd1 <- initAFFMeanCD(1, 0.9, BL=50, lambdaBL=-1))
        
        
        #should FAIL - lambdaBL not in range 
        #expect_error(fffcd1 <- initAFFMeanCD(1, 0.9, BL=50, lambdaBL=1.1))
        })



test_that("AFF change detector is initialised correctly (pass)", {

        #----------------------------------#
        ###TESTS THAT PASS
        #----------------------------------#

        #SHOULD PASS
        expect_error(fffcd1 <- initAFFMeanCD(eta=0.01, alpha=0.01), NA) 

        #SHOULD PASS - prints nothing
        expect_error(fffcd1 <- initAFFMeanCD(eta=0.01, alpha=0.01, BL=50), NA)
        
        #SHOULD PASS - prints nothing
        expect_error(fffcd1 <- initAFFMeanCD(eta=0.01, alpha=0.01, BL=50.1), NA)

        #SHOULD PASS - prints nothing
        expect_error(fffcd1 <- initAFFMeanCD(eta=0.01, alpha=0.01, BL=50.1), NA)

        #SHOULD PASS - prints nothing
        expect_error(fffcd1 <- initAFFMeanCD(eta=0.01, alpha=0.01, BL=2), NA)

        #NOTE: case BL=1 must fail, but BL=0 is okay, the user must just set
        #prechange values for mean and variance

        #SHOULD PASS - prints nothing
        expect_error(fffcd1 <- initAFFMeanCD(eta=0.01, alpha=0.01, BL=0), NA)

        #SHOULD PASS
        expect_error(fffcd1 <- initAFFMeanCD(), NA) 

        })

