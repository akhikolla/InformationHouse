context("Test  3: initialisation of FFF, AFF")

#will be used in the tests


test_that("FFF mean is initialised correctly", {

        #----------------------------------#
        ###TESTS THAT FAIL 
        #----------------------------------#

        #need to do two, because error message cannot handle
        #parentheses, e.g. (value=1), specifically the open 
        #parenthesis, for some reason.

        #should not work - not numeric
        expect_error(fff1 <- initFFFMean("A"), "lambda")
        expect_error(fff1 <- initFFFMean("A"), "not a finite")
       #lambda (value=A)  
        #should not work - Na does not exist
        expect_error(fff1 <- initFFFMean(Na), "not found")
        
        #should not work - NaN not finite
        expect_error(fff1 <- initFFFMean(NaN), "lambda")
        expect_error(fff1 <- initFFFMean(NaN), "not a finite")

        #should not work - NA not available
        expect_error(fff1 <- initFFFMean(NA), "lambda")
        expect_error(fff1 <- initFFFMean(NA), "not a finite")

        #should not work - 1.1 not in range
        expect_error(fff1 <- initFFFMean(1.1), "lambda")
        expect_error(fff1 <- initFFFMean(1.1), "not in range")

        #should not work - -1 not in range
        expect_error(fff1 <- initFFFMean(-1), "lambda")
        expect_error(fff1 <- initFFFMean(-1), "not in range")


        #----------------------------------#
        ###TESTS THAT PASS
        #----------------------------------#

        #old test assumed nothign was printed
        #expect_that(fff1 <- initFFFMean(), prints_text(""))
        #same as expect_error(..., NA)
        expect_error(fff1 <- initFFFMean(), NA)

        #SAME AS ABOVE
        #SHOULD PASS - prints nothing
        #expect_output(fff1 <- initFFFMean(), "")

        #SHOULD PASS - prints nothing
        expect_error(fff1 <- initFFFMean(0), NA)
        
        #SHOULD PASS - prints nothing
        expect_error(fff1 <- initFFFMean(1), NA)
        
        #SHOULD PASS - prints nothing
        expect_error(fff1 <- initFFFMean(lambda=0.5), NA)
        })






test_that("AFF mean is initialised correctly", {

        #----------------------------------#
        ###TESTS THAT FAIL 
        #----------------------------------#

        #should not work - not numeric
        expect_error(aff1 <- initAFFMean("A"), "eta")
        expect_error(aff1 <- initAFFMean("A"), "not a finite")
        
        #should not work - Na does not exist
        expect_error(aff1 <- initAFFMean(Na), "not found")
        
        #should not work - NaN not finite
        expect_error(aff1 <- initAFFMean(NaN), "eta")
        expect_error(aff1 <- initAFFMean(NaN), "not a finite")

        #should not work - NA not available
        expect_error(aff1 <- initAFFMean(NA), "eta")
        expect_error(aff1 <- initAFFMean(NA), "not a finite")

        #should not work - 1.1 not in range
        expect_error(aff1 <- initAFFMean(1.1), "eta")
        expect_error(aff1 <- initAFFMean(1.1), "not in range")

        #should not work - -1 not in range
        expect_error(aff1 <- initAFFMean(-1), "eta")
        expect_error(aff1 <- initAFFMean(-1), "not in range")


        #----------------------------------#
        ###TESTS THAT PASS
        #----------------------------------#

        #very simple - prints nothing
        #expect_that(aff1 <- initAFFMean(), prints_text(""))
        expect_error(aff1 <- initAFFMean(), NA)

        #SAME AS ABOVE
        #SHOULD PASS - prints nothing
        #expect_output(aff1 <- initAFFMean(), "")

        #SHOULD PASS - prints nothing
        #expect_output(aff1 <- initAFFMean(0), "")
        expect_error(aff1 <- initAFFMean(0), NA)
        
        #SHOULD PASS - prints nothing
        #expect_output(aff1 <- initAFFMean(1), "")
        expect_error(aff1 <- initAFFMean(1), NA)

        #SHOULD PASS - prints nothing
        #expect_output(aff1 <- initAFFMean(eta=0.1), "")
        expect_error(aff1 <- initAFFMean(eta=0.1), NA)
        })

