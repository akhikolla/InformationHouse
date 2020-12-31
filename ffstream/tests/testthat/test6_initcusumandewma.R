context("Test  6: initialisation of CUSUM")

test_that("init CUSUM (k fails)", {

        expect_error(cusumcd1 <- initCUSUMMeanCD(k="A", h=8.00, BL=50), " k ")
        expect_error(cusumcd1 <- initCUSUMMeanCD(k="A", h=8.00, BL=50), " not a finite ")

        expect_error(cusumcd1 <- initCUSUMMeanCD(k=NaN, h=8.00, BL=50), " k ")
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=NaN, h=8.00, BL=50), " not a finite ")
        
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=Na, h=8.00, BL=50), "not found")

        expect_error(cusumcd1 <- initCUSUMMeanCD(k=NA, h=8.00, BL=50), " k ")
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=NA, h=8.00, BL=50), " not a finite ")


        })


test_that("init CUSUM (h fails)", {

        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h="A", BL=50), " h ")
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h="A", BL=50), "not a finite")

        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h=NaN, BL=50), " h ")
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h=NaN, BL=50), "not a finite")

        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h=Na, BL=50), "not found")

        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h=NA, BL=50), " h ")
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h=NA, BL=50), "not a finite")

        })


test_that("CUSUM init (BL)", {
        #should FAIL - BL not in range 
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h=8.00, BL=-1), "BL")
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h=8.00, BL=-1), "not greater than")
        })


test_that("CUSUM init (pass)", {

        #----------------------------------#
        ###TESTS THAT PASS
        #----------------------------------#
        #should pass
        expect_error(cusumcd1 <- initCUSUMMeanCD(BL=50), NA)

        #should pass
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25), NA)
        
        #should pass
        expect_error(cusumcd1 <- initCUSUMMeanCD(h=8.00), NA)

        #should pass
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, BL=50), NA)
        
        #should pass
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h=8.00), NA)
        
        #should pass
        expect_error(cusumcd1 <- initCUSUMMeanCD(h=8.00, BL=50), NA)

        #should pass
        expect_error(cusumcd1 <- initCUSUMMeanCD(k=0.25, h=8.00, BL=50), NA)

        })

