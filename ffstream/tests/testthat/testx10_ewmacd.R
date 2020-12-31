context("Test 10: EWMA change detection")


library(Rcpp)
changepointStr <<- "tauhat"

test_that("EWMA change detection on stream 1", {
        stream1 <- makeStreamMeanChangeR(seednum=1, numChanges=3)

        r_1 <- 0.20
        L_1 <- 3.00
        BL1 <- 50
        ewmaparams1 <- c(r_1, L_1)
        changesDetectedOriginal1 <- EWMA_stream_jumpdetect(stream1, BL1, ewmaparams1)

        ewma1 <- initEWMAMeanCD(r = r_1, L = L_1, BL=BL1)
        returnList <- ewma1$processVectorSave(stream1)
        changesDetectedRcpp1 <- returnList[[changepointStr]]

        #need to have -1 for R version (change signalled at start of new regime in R
        #rather than at the end of the old regime...
        expect_equal(changesDetectedOriginal1-1, changesDetectedRcpp1)
        })


test_that("EWMA change detection on stream 2", {
        stream2 <- makeStreamMeanChangeR(seednum=2, numChanges=5)
        r_2 <- 0.05
        L_2 <- 2.615
        BL2 <- 50
        ewmaparams2 <- c(r_2, L_2)
        changesDetectedOriginal2 <- EWMA_stream_jumpdetect(stream2, BL2, ewmaparams2)

        ewma2 <- initEWMAMeanCD(r = r_2, L = L_2, BL=BL2)
        returnList <- ewma2$processVectorSave(stream2)
        changesDetectedRcpp2 <- returnList[[changepointStr]]
        expect_equal(changesDetectedOriginal2-1, changesDetectedRcpp2)
        })



#now checking getting/setting BL (derived class!)
test_that("checking BL 1", {
        r_3 <- 0.05
        L_3 <- 2.615
        BL3 <- 87
        ewma3 <- initEWMAMeanCD(r = r_3, L = L_3, BL=BL3)

        #default values
        expect_equal(ewma3$BL, 87)
        })

test_that("checking BL 2", {
        r_3 <- 0.05
        L_3 <- 2.615
        BL3 <- 87
        ewma3 <- initEWMAMeanCD(r = r_3, L = L_3, BL=BL3)
        #now change the BL
        ewma3$BL <- 100
        expect_equal(ewma3$BL, 100)

        })




#now checking pval (getter, only)
test_that("checking pval", {
        r_4 <- 0.05
        L_4 <- 2.615
        BL4 <- 50
        ewma4 <- initEWMAMeanCD(r = r_4, L = L_4, BL=BL4)

        #default value is 0.5
        expect_equal(ewma4$pval, 0.5)

        })




test_that("EWMA change detection on stream 1", {
        stream1 <- makeStreamMeanChangeR(seednum=3, numChanges=3)

        r_1 <- 0.20
        L_1 <- 3.00
        BL1 <- 50
        ewmaparams1 <- c(r_1, L_1)
        changesDetectedOriginal1 <- EWMA_stream_jumpdetect(stream1, BL1, ewmaparams1)

        ewma1 <- initEWMAMeanCD(r = r_1, L = L_1, BL=BL1)
        returnList <- ewma1$processVectorSave(stream1)
        changesDetectedRcpp1 <- returnList[[changepointStr]]

        #need to have -1 for R version (change signalled at start of new regime in R
        #rather than at the end of the old regime...
        expect_equal(changesDetectedOriginal1-1, changesDetectedRcpp1)
        })



test_that("checking detectVector multiple", {
        stream1 <- makeStreamMeanChangeR(seednum=1, numChanges=3)

        r_1 <- 0.20
        L_1 <- 3.00
        BL1 <- 50
        ewmaparams1 <- c(r_1, L_1)
        changesDetectedOriginal1 <- EWMA_stream_jumpdetect(stream1, BL1, ewmaparams1)

        returnList <- detectEWMAMean(x=stream1, r=r_1, L=L_1, BL=BL1, multiple=TRUE)
        changesDetectedRcpp1 <- returnList[[changepointStr]]

        #need to have -1 for R version (change signalled at start of new regime in R
        #rather than at the end of the old regime...
        expect_equal(changesDetectedOriginal1-1, changesDetectedRcpp1)
        })



#testing single change
test_that("checking detectorVector single", {
        stream1 <- makeStreamMeanChangeR(seednum=2, numChanges=3)

        r_1 <- 0.20
        L_1 <- 3.00
        BL1 <- 50
        ewmaparams1 <- c(r_1, L_1)
        changesDetectedOriginal1 <- EWMA_stream_jumpdetect(stream1, BL1, ewmaparams1)
    changesDetectedOriginal1 <- changesDetectedOriginal1[1]

        returnList<- detectEWMAMean(x=stream1, r=r_1, L=L_1, BL=BL1, single=TRUE)
        changesDetectedRcpp1 <- returnList[[changepointStr]]
        

        #only check for first element, because single change
        expect_equal(changesDetectedOriginal1-1, changesDetectedRcpp1)
        
        })


#testing single change
test_that("checking detectorVector single (v2)", {
        stream1 <- makeStreamMeanChangeR(seednum=2, numChanges=2, mu0=2, sigma=3)

        r_1 <- 0.3
        L_1 <- 2.54
        BL1 <- 50
        ewmaparams1 <- c(r_1, L_1)
        changesDetectedOriginal1 <- EWMA_stream_jumpdetect(stream1, BL1, ewmaparams1)
    changesDetectedOriginal1 <- changesDetectedOriginal1[1]

        returnList<- detectEWMAMean(x=stream1, r=r_1, L=L_1, BL=BL1, single=TRUE)
        changesDetectedRcpp1 <- returnList[[changepointStr]]
        

        #only check for first element, because single change
        expect_equal(changesDetectedOriginal1-1, changesDetectedRcpp1)
        
        })


#testing single change WITH PRECHANGE KNOWN
test_that("checking detectorVector single with prechange", {
        mu0 <- 0
        sigma0 <- 1
        stream1 <- makeStreamMeanChangeR(seednum=2, numChanges=3, mu0=mu0, sigma0=sigma0)

        r_1 <- 0.20
        L_1 <- 3.00
        BL1 <- 0
        ewmaparams1 <- c(r_1, L_1)

        changesDetectedOriginal1prechange <- EWMA_stream_jumpdetect_prechange(stream1, 0, ewmaparams1, mu0, sigma0)

        #must subtract 1
    changesDetectedOriginal1prechange <- changesDetectedOriginal1prechange - 1



        returnList<- detectEWMAMean(x=stream1, r=r_1, L=L_1, BL=0, single=TRUE, usePrechange=TRUE, prechangeMean=mu0, prechangeSigma=sigma0)
        changesDetectedRcpp1 <- returnList[[changepointStr]]


        #only check for first element, because single change
        expect_equal(changesDetectedOriginal1prechange, changesDetectedRcpp1)

        })



#testing single change WITH PRECHANGE KNOWN
test_that("checking detectorVector single with prechange (v2)", {
        mu0 <- 0
        sigma0 <- 1
        stream1 <- makeStreamMeanChangeR(seednum=5, numChanges=2, mu0=mu0, sigma0=sigma0)

        r_1 <- 0.3
        L_1 <- 2.54
        BL1 <- 70
        ewmaparams1 <- c(r_1, L_1)

        changesDetectedOriginal1prechange <- EWMA_stream_jumpdetect_prechange(stream1, BL1, ewmaparams1, mu0, sigma0)

        #must subtract 1
    changesDetectedOriginal1prechange <- changesDetectedOriginal1prechange - 1



        returnList<- detectEWMAMean(x=stream1, r=r_1, L=L_1, BL=BL1, single=TRUE, usePrechange=TRUE, prechangeMean=mu0, prechangeSigma=sigma0)
        changesDetectedRcpp1 <- returnList[[changepointStr]]


        #only check for first element, because single change
        expect_equal(changesDetectedOriginal1prechange, changesDetectedRcpp1)

        })
