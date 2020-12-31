context("Test  9: CUSUM change detection")

changepointStr <<- "tauhat"

test_that("CUSUM change detection on stream 1", {
        stream1 <- makeStreamMeanChangeR(seednum=1, numChanges=3)

        d_1 <- 0.25
        B_1 <- 8.01
        BL1 <- 50
        cusumparams1 <- c(d_1, B_1)
        changesDetectedOriginal1 <- CUSUM_stream_jumpdetect(stream1, BL1,cusumparams1)

        cusum1 <- initCUSUMMeanCD(k=d_1, h=B_1, BL=BL1)
        returnList <- cusum1$processVectorSave(stream1)
        changesDetectedRcpp1 <- returnList[[changepointStr]]

        #need to have -1 for R version (change signalled at start of new regime in R
        #rather than at the end of the old regime...
        expect_equal(changesDetectedOriginal1-1, changesDetectedRcpp1)
        })


test_that("CUSUM change detection on stream 2", {
        stream2 <- makeStreamMeanChangeR(seednum=2, numChanges=5)
        d_2 <- 0.50
        B_2 <- 4.77
        BL2 <- 50
        cusumparams2 <- c(d_2, B_2)
        changesDetectedOriginal2 <- CUSUM_stream_jumpdetect(stream2, BL2,cusumparams2)

        cusum2 <- initCUSUMMeanCD(k=d_2, h=B_2, BL=BL2)
        returnList <- cusum2$processVectorSave(stream2)
        changesDetectedRcpp2 <- returnList[[changepointStr]]

        expect_equal(changesDetectedOriginal2-1, changesDetectedRcpp2)
        })



#now checking getting/setting BL (derived class!)
test_that("checking BL", {
        d_3 <- 0.50
        B_3 <- 4.77
        BL3 <- 87
        cusum3 <- initCUSUMMeanCD(k=d_3, h=B_3, BL=BL3)

        #default values
        expect_equal(cusum3$BL, 87)

        #now change the BL
        cusum3$BL <- 100
        expect_equal(cusum3$BL, 100)

        })




#now checking pval (getter, only)
test_that("checking pval", {
        d_4 <- 0.50
        B_4 <- 4.77
        BL4 <- 87
        cusum4 <- initCUSUMMeanCD(k=d_4, h=B_4, BL=BL4)

        #default value is 0.5
        expect_equal(cusum4$pval, 0.5)

        })



#testing multiple change
test_that("checking detectorVector multiple", {
        stream1 <- makeStreamMeanChangeR(seednum=1, numChanges=3)

        d_1 <- 0.25
        B_1 <- 8.01
        BL1 <- 50
        cusumparams1 <- c(d_1, B_1)
        changesDetectedOriginal1 <- CUSUM_stream_jumpdetect(stream1, BL1,cusumparams1)

        returnList<- detectCUSUMMean(x=stream1, k=d_1, h=B_1, BL=BL1, multiple=TRUE)
        changesDetectedRcpp1 <- returnList[[changepointStr]]

        expect_equal(changesDetectedOriginal1-1, changesDetectedRcpp1)
        
        })



#testing single change
test_that("checking detectorVector single", {
        stream1 <- makeStreamMeanChangeR(seednum=2, numChanges=3)

        d_1 <- 0.25
        B_1 <- 8.01
        BL1 <- 50
        cusumparams1 <- c(d_1, B_1)
        changesDetectedOriginal1 <- CUSUM_stream_jumpdetect(stream1, BL1,cusumparams1)
    changesDetectedOriginal1 <- changesDetectedOriginal1[1]

        returnList<- detectCUSUMMean(x=stream1, k=d_1, h=B_1, BL=BL1, single=TRUE)
        changesDetectedRcpp1 <- returnList[[changepointStr]]
        

#        cat("\nSingle change prechange unknown: \n")
#        cat("Rcpp: single: ")
#        print(changesDetectedRcpp1)
#        cat("R original: ")
#        print(changesDetectedOriginal1)
#        cat("\n here\n\n ")

        #only check for first element, because single change
        expect_equal(changesDetectedOriginal1-1, changesDetectedRcpp1)
        
        })

#testing single change
test_that("checking detectorVector single (v2)", {
        stream1 <- makeStreamMeanChangeR(seednum=2, numChanges=2, mu0=2, sigma0=3)

        d_2 <- 0.5
        B_2 <- 4.77
        BL2 <- 30
        cusumparams1 <- c(d_2, B_2)
        changesDetectedOriginal2 <- CUSUM_stream_jumpdetect(stream1, BL2,cusumparams1)
    changesDetectedOriginal2 <- changesDetectedOriginal2[1]

        returnList<- detectCUSUMMean(x=stream1, k=d_2, h=B_2, BL=BL2, single=TRUE)
        changesDetectedRcpp2 <- returnList[[changepointStr]]
        
        #only check for first element, because single change
        expect_equal(changesDetectedOriginal2-1, changesDetectedRcpp2)
        
        })



#testing single change WITH PRECHANGE KNOWN
test_that("checking detectorVector single with prechange", {
        mu0 <- 0
        sigma0 <- 1
        stream1 <- makeStreamMeanChangeR(seednum=2, numChanges=3, mu0=mu0, sigma0=sigma0)

        d_1 <- 0.25
        B_1 <- 8.01
        BL1 <- 50

        cusumparams1 <- c(d_1, B_1)
        changesDetectedOriginal1prechange <- CUSUM_stream_jumpdetect_prechange(stream1, BL1,cusumparams1, mu0, sigma0)

        #must subtract 1
    changesDetectedOriginal1prechange <- changesDetectedOriginal1prechange - 1


        returnList<- detectCUSUMMean(x=stream1, k=d_1, h=B_1, BL=BL1, single=TRUE, usePrechange=TRUE, prechangeMean=mu0, prechangeSigma=sigma0)
        changesDetectedRcpp1 <- returnList[[changepointStr]]


        #only check for first element, because single change
        expect_equal(changesDetectedOriginal1prechange, changesDetectedRcpp1)

        
        })



#testing single change WITH PRECHANGE KNOWN
test_that("checking detectorVector single with prechange (v2)", {
        mu0 <- 2
        sigma0 <- 2
        stream2 <- makeStreamMeanChangeR(seednum=3, numChanges=2, mu0=mu0, sigma0=sigma0)

        d_2 <- 0.5
        B_2 <- 4.77

        cusumparams2 <- c(d_2, B_2)
        changesDetectedOriginal2prechange <- CUSUM_stream_jumpdetect_prechange(stream2, 0,cusumparams2, mu0, sigma0)

        #must subtract 1
    changesDetectedOriginal2prechange <- changesDetectedOriginal2prechange - 1


        returnList<- detectCUSUMMean(x=stream2, k=d_2, h=B_2, BL=0, single=TRUE, usePrechange=TRUE, prechangeMean=mu0, prechangeSigma=sigma0)
        changesDetectedRcpp2 <- returnList[[changepointStr]]

        #only check for first element, because single change
        expect_equal(changesDetectedOriginal2prechange, changesDetectedRcpp2)

        
        })
