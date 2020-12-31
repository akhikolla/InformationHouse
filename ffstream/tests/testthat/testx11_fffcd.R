context("Test 11: FFF change detection")




test_that("test checkIfChange", {
#will be used in the tests
#        fffcd1 <- new(FFFChangeDetector, 0.9, 0.05)
        #BL is default
        fffcd1 <- initFFFMeanCD(lambda=0.9, alpha=0.05)
        fffcd1$streamEstMean <- 0
        fffcd1$streamEstSigma <- 1

        #nothing has happened, just checking values
        expect_equal(fffcd1$fffxbar, 0)
        expect_equal(fffcd1$changeDetected, FALSE)

        #there should not be a change
        fffcd1$checkIfChange()
        expect_equal(fffcd1$changeDetected, FALSE)
        
        #artificially forcing xbar - just BELOW change threshold
        fffcd1$fffxbar <- 1.95
        fffcd1$checkIfChange()
        expect_equal(fffcd1$changeDetected, FALSE)
        
        #artificially forcing xbar - just ABOVE change threshold
        fffcd1$fffxbar <- 1.96
        fffcd1$checkIfChange()
        expect_equal(fffcd1$changeDetected, TRUE)
        })

test_that("burn in estimator works", {
        set.seed(1)
        BL <- 60
        N <- 100
        x <- rnorm(N)
#        fffcd2 <- new(FFFChangeDetector, 0.95, 0.01, BL)
        fffcd2 <- initFFFMeanCD(lambda=0.95, BL=BL, alpha=0.01)
#        fffcd2$print()
        
        expect_equal(fffcd2$streamEstMean, 0)
        expect_equal(fffcd2$streamEstSigma, 0)
        
        #to do it 'manually' with a for loop
#        for (i in 1:length(x)){
#            fffcd2$update(x[i])
#        }
        fffcd2$processVector(x)

        burnMean <- mean(x[1:BL]) 
        burnSigma <- sqrt(var(x[1:BL]))
        expect_equal(fffcd2$streamEstMean, burnMean)
        expect_equal(fffcd2$streamEstSigma, burnSigma)
        
        })


test_that("test it is possible to set streamEstMean", {
        fffcd3 <- initFFFMeanCD(lambda=0.95, alpha=0.01, BL=50)
        fffcd3$streamEstSigma <- 3
        expect_equal(fffcd3$streamEstSigma, 3)
        })


#now checking getting/setting BL (derived class!)
test_that("checking BL", {
        alpha4 <- 0.05
        BL <- 93
        lambda4 <- 0.99
        fffcd4 <- initFFFMeanCD(BL=BL, alpha=alpha4, lambda=lambda4)

        #default values
        expect_equal(fffcd4$BL, 93)

        #now change the BL
        fffcd4$BL <- 100
        expect_equal(fffcd4$BL, 100)

        })




#now checking pval (getter, only)
test_that("checking pval", {
        alpha5 <- 0.05
        BL <- 93
        lambda5 <- 0.99
        fffcd5 <- initFFFMeanCD(BL=BL, alpha=alpha5, lambda=lambda5)

        #default value is 0.5
        expect_equal(fffcd5$pval, 0.5)

        })





test_that("test 1 of FFF CD on stream - compare to R version", {
        #will be used in the tests
        stream1 <- makeStreamMeanChangeR(seednum=1)
        BL1 <- 50

        #set params
        lambda1 <- 0.95
        p1 <- 0.99
        #this is a number, if zero, that will reset u, v, w, ffmean after a change.
        #if not zero, then will continue with old ffmean
        resettozero <- 1
        u_init <- 0
        v_init <- 0
        w_init <- 0
        ffmean_init <- 0
        ffvar_init <- 0
        initlist <- list(resettozero, u_init, v_init, w_init, ffmean_init, ffvar_init)
        fffparams1 <- c( list(lambda1, p1), initlist )

        #R version:
        changepointsR <- FFF_stream_jumpdetect(stream1, BL1, fffparams1)

        #now cpp:
        fffcd1 <- new(FFFChangeDetector, 0.95, 0.01)
        returnList <- fffcd1$processVectorSave(stream1)
        changepointsCpp <- returnList$tauhat

        expect_equal(changepointsR, changepointsCpp)
        })



test_that("test 2 of FFF CD on stream - compare to R version", {
        #will be used in the tests
        stream2 <- makeStreamMeanChangeR(seednum=2, numChanges=5)
        BL2 <- 50

        #set params
        lambda2 <- 0.99
        p2 <- 0.95
        #this is a number, if zero, that will reset u, v, w, ffmean after a change.
        #if not zero, then will continue with old ffmean
        resettozero <- 1
        u_init <- 0
        v_init <- 0
        w_init <- 0
        ffmean_init <- 0
        ffvar_init <- 0
        initlist <- list(resettozero, u_init, v_init, w_init, ffmean_init, ffvar_init)
        fffparams2 <- c( list(lambda2, p2), initlist )

        #R version:
        changepointsR <- FFF_stream_jumpdetect(stream2, BL2, fffparams2)

        #now cpp:
        fffcd2 <- initFFFMeanCD(lambda=0.99, alpha=0.05)
        returnList <- fffcd2$processVectorSave(stream2)
        changepointsCpp <- returnList$tauhat

        expect_equal(changepointsR, changepointsCpp)
        })





test_that("detectMultiple works for FFFMeanCD", {
        stream1 <- makeStreamMeanChangeR(seednum=1, numChanges=3)
        lambda1 <- 0.95
        alpha1 <- 0.01
        BL1 <- 50


        #set params
        p1 <- 1-alpha1
        #this is a number, if zero, that will reset u, v, w, ffmean after a change.
        #if not zero, then will continue with old ffmean
        resettozero <- 1
        u_init <- 0
        v_init <- 0
        w_init <- 0
        ffmean_init <- 0
        ffvar_init <- 0
        initlist <- list(resettozero, u_init, v_init, w_init, ffmean_init, ffvar_init)
        fffparams1 <- c( list(lambda1, p1), initlist )

        #R version:
        changepointsR <- FFF_stream_jumpdetect(stream1, BL1, fffparams1)


        #now cpp:
        #fffcd1 <- new(FFFChangeDetector, lambda1, alpha1, BL1)
        #returnList <- fffcd1$processVectorSave(stream1)
        #changepointsCpp <- returnList$changepoints
        #expect_equal(tauhatList$tauhat, changepointsCpp)

        tauhatList <- cpp_detectFFFMeanMultiple(stream1, lambda1, alpha1, BL1)

        expect_equal(tauhatList$tauhat, changepointsR)
        #print(tauhatList)

        tauhatList2 <- detectFFFMean(x=stream1, lambda=lambda1, alpha=alpha1, single=FALSE)

        #now using wrapper in R
        expect_equal(tauhatList2$tauhat, changepointsR)
        
        })



test_that("detectSingle works for FFFMeanCD", {
        stream1 <- makeStreamMeanChangeR(seednum=1, numChanges=2)
        lambda1 <- 0.95
        alpha1 <- 0.01
        BL1 <- 50


        #set params
        p1 <- 1-alpha1
        #this is a number, if zero, that will reset u, v, w, ffmean after a change.
        #if not zero, then will continue with old ffmean
        resettozero <- 1
        u_init <- 0
        v_init <- 0
        w_init <- 0
        ffmean_init <- 0
        ffvar_init <- 0
        initlist <- list(resettozero, u_init, v_init, w_init, ffmean_init, ffvar_init)
        fffparams1 <- c( list(lambda1, p1), initlist )

        #R version:
        changepointsR <- FFF_stream_jumpdetect(stream1, BL1, fffparams1)
        #only keep first changepoint
        changepointsR <- changepointsR[1]

        tauhatList <- cpp_detectFFFMeanSingle(stream1, lambda1, alpha1, BL1)

        expect_equal(tauhatList$tauhat, changepointsR)
        #print(tauhatList)
        
        tauhatList2 <- detectFFFMean(x=stream1, lambda=lambda1, alpha=alpha1, single=TRUE, usePrechange=FALSE)

        #now using wrapper in R
        expect_equal(tauhatList2$tauhat, changepointsR)
        })


test_that("detectSingle works for FFFMeanCD with prechange known", {
        stream1 <- makeStreamMeanChangeR(seednum=5, numChanges=5)
        lambda1 <- 0.95
        alpha1 <- 0.01
        BL1 <- 50

        #set params
        p1 <- 1-alpha1
        #this is a number, if zero, that will reset u, v, w, ffmean after a change.
        #if not zero, then will continue with old ffmean
        resettozero <- 1
        u_init <- 0
        v_init <- 0
        w_init <- 0
        ffmean_init <- 0
        ffvar_init <- 0
        initlist <- list(resettozero, u_init, v_init, w_init, ffmean_init, ffvar_init)
        fffparams1 <- c( list(lambda1, p1), initlist )

        prechangeMean <- 0
        prechangeSigma <- 1

        #no BL
        fffcd1 <- initFFFMeanCD(lambda=lambda1, alpha=alpha1, BL=0)
        fffcd1$streamEstMean <- prechangeMean
        fffcd1$streamEstSigma <- prechangeSigma
        changeNotFound <- TRUE
        index <- 0
        tauhatSingle <- 0
        while ( (changeNotFound) & (index < length(stream1)) ) {
            index <- index + 1
            x <- stream1[index]
            fffcd1$update(x)
            if (fffcd1$changeDetected){
                tauhatSingle <- index
                changeNotFound <- FALSE
            }
        }

        tauhatList <- cpp_detectFFFMeanSinglePrechange(stream1, lambda1, alpha1, prechangeMean, prechangeSigma)

        #compare 
        expect_equal(tauhatList$tauhat, tauhatSingle)
        #print(tauhatList)

        tauhatList2 <- detectFFFMean(x=stream1, lambda=lambda1, alpha=alpha1, multiple=FALSE, usePrechange=TRUE, prechangeMean=prechangeMean, prechangeSigma=prechangeSigma)


        #now using wrapper in R
        expect_equal(tauhatList2$tauhat, tauhatSingle)

        tauhatR <- FFF_stream_jumpdetect_prechange(stream1, 0, fffparams1, 
                                                prechangeMean, prechangeSigma)

        #now using wrapper in R
        expect_equal(tauhatR, tauhatSingle)
        expect_equal(tauhatR, tauhatList2$tauhat)
#        print(length(tauhatR))
        
        })
