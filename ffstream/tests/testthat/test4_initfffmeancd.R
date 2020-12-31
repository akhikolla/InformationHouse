context("Test  4: initialisation of FFF mean cd")

#will be used in the tests

#test 1: alpha
test_that("fff change detector is initialised correctly (alpha)", {

        #----------------------------------#
        ###TESTS THAT FAIL 
        #----------------------------------#

        #should FAIL - alpha not numeric
        expect_error(fffcd1 <- initFFFMeanCD("A", 0.9), "alpha")
        expect_error(fffcd1 <- initFFFMeanCD("A", 0.9), "not a finite")

        #should FAIL - alpha not finite 
        expect_error(fffcd1 <- initFFFMeanCD(NaN, 0.9), "alpha")
        expect_error(fffcd1 <- initFFFMeanCD(NaN, 0.9), "not a finite")
        
        #should FAIL - alpha not finite 
        expect_error(fffcd1 <- initFFFMeanCD(Na, 0.9), "not found")
        
        #should FAIL - alpha not finite 
        expect_error(fffcd1 <- initFFFMeanCD(NA, 0.9), "alpha")
        expect_error(fffcd1 <- initFFFMeanCD(NA, 0.9), "not a finite")

        #should FAIL - alpha not in range 
        expect_error(fffcd1 <- initFFFMeanCD(-1, 0.9), "alpha")
        expect_error(fffcd1 <- initFFFMeanCD(-1, 0.9), "not in range")
        
        #should FAIL - alpha lambda not in range 
        expect_error(fffcd1 <- initFFFMeanCD(1.1, 0.9), "alpha")
        expect_error(fffcd1 <- initFFFMeanCD(1.1, 0.9), "not in range")

        #should FAIL - alpha not in range 
        expect_error(fffcd1 <- initFFFMeanCD(1, 0), "alpha")
        expect_error(fffcd1 <- initFFFMeanCD(1, 0), "not in range")

        #should FAIL - alpha not in range 
        expect_error(fffcd1 <- initFFFMeanCD(1, 1), "alpha")
        expect_error(fffcd1 <- initFFFMeanCD(1, 1), "not in range")
        })


#test 2: lambda
test_that("fff change detector is initialised correctly (lambda)", {

        #should FAIL - lambda not in range 
        expect_error(fffcd1 <- initFFFMeanCD(lambda=1.1, alpha=0.01), "lambda")
        expect_error(fffcd1 <- initFFFMeanCD(lambda=1.1, alpha=0.01), "not in range")
        })




#test 3: BL 
test_that("fff change detector is initialised correctly (BL)", {
        #should FAIL - BL not in range 
        expect_error(fffcd1 <- initFFFMeanCD(0.01, 0.9, BL=-1), "BL")
        expect_error(fffcd1 <- initFFFMeanCD(0.01, 0.9, BL=-1), "not greater than")

        #no more lambdaBL
        #should FAIL - lambdaBL not numeric 
        #expect_error(fffcd1 <- initFFFMeanCD(1, 0.9, BL=50, lambdaBL="A"))

        #should FAIL - lambdaBL not finite 
        #expect_error(fffcd1 <- initFFFMeanCD(1, 0.9, BL=50, lambdaBL=NA))

        #should FAIL - lambdaBL not numeric 
        #expect_error(fffcd1 <- initFFFMeanCD(1, 0.9, BL=50, lambdaBL=NaN))

        #should FAIL - lambdaBL not numeric 
        #expect_error(fffcd1 <- initFFFMeanCD(1, 0.9, BL=50, lambdaBL=Na))

        #should FAIL - lambdaBL not in range 
        #expect_error(fffcd1 <- initFFFMeanCD(1, 0.9, BL=50, lambdaBL=-1))
        
        
        #should FAIL - lambdaBL not in range 
        #expect_error(fffcd1 <- initFFFMeanCD(1, 0.9, BL=50, lambdaBL=1.1))

        })


#test 4: pass
test_that("fff change detector is initialised correctly (pass)", {
        #----------------------------------#
        ###TESTS THAT PASS
        #----------------------------------#

        #SHOULD PASS
        expect_error(fffcd1 <- initFFFMeanCD(lambda=0.99, alpha=0.01), NA) 

        #SHOULD PASS - prints nothing
        expect_error(fffcd1 <- initFFFMeanCD(lambda=0.99, alpha=0.01, BL=50), NA)
        
        #SHOULD PASS - prints nothing
        expect_error(fffcd1 <- initFFFMeanCD(lambda=0.99, alpha=0.01, BL=50.1), NA)

        #SHOULD PASS - prints nothing
        expect_error(fffcd1 <- initFFFMeanCD(lambda=0.99, alpha=0.01, BL=50.1), NA)

        #SHOULD PASS - prints nothing
        expect_error(fffcd1 <- initFFFMeanCD(lambda=0.99, alpha=0.01, BL=0), NA)

        #SHOULD PASS - no commands, all defaults
        expect_error(fffcd1 <- initFFFMeanCD(), NA) 
        })

