context("Test 12: AFF change detection")




#very similar to fff tests
#just checking if change occurs
test_that("test checkIfChange", {
#will be used in the tests
        #new initialisation
        #affcd1 <- new(AFFChangeDetector, 0.05)
        affcd1 <- initAFFMeanCD(alpha=0.05)

        affcd1$streamEstMean <- 0
        affcd1$streamEstSigma <- 1
#        affcd1$setStreamEstMean(0)
#        affcd1$setStreamEstSigma(1)

        #nothing has happened, just checking values
        expect_equal(affcd1$affxbar, 0)
        expect_equal(affcd1$changeDetected, FALSE)

        #there should NOT be a change
        affcd1$checkIfChange()
        expect_equal(affcd1$changeDetected, FALSE)
        
        #artificially forcing xbar - just BELOW change threshold
        #there should NOT be a change
        affcd1$affxbar <- 1.95
        affcd1$checkIfChange()
        expect_equal(affcd1$changeDetected, FALSE)
        
        #artificially forcing xbar - just ABOVE change threshold
        #there SHOULD BE a change
        affcd1$affxbar <- 1.96
        affcd1$checkIfChange()
        expect_equal(affcd1$changeDetected, TRUE)
        })


#now checking the burn-in estimator
test_that("burn in estimator works", {
        set.seed(1)
        eta1 <- 0.01
        alpha1 <- 0.05
        lambdaDuringBL <- 1
        BL <- 50
        N <- 100
        x <- rnorm(N)
        affcd2 <- initAFFMeanCD(alpha1, eta1, BL)
        
        expect_equal(affcd2$streamEstMean, 0)
        expect_equal(affcd2$streamEstSigma, 0)



        xstart <- x[1:BL]
        xend <- x[(BL+1):N]

        #do it 'manually' with a for loop
#        for (i in 1:length(xstart)){
#            affcd2$update(x[i])
#        }

        affcd2$processVector(xstart)
        burnMean <- mean(x[1:BL]) 
        burnSigma <- sqrt(var(x[1:BL]))
        
        })


#now checking the set/get of streamEst mean and sigma
test_that("checking set/get streamEst mean and sigma", {
        eta3 <- 0.01
        alpha3 <- 0.05
        BL <- 50
        affcd3 <- initAFFMeanCD(BL=BL, alpha=alpha3, eta=eta3)

        #default values
        expect_equal(affcd3$streamEstMean, 0)
        expect_equal(affcd3$streamEstSigma, 0)

        #now set them
        affcd3$streamEstMean <- 3
        expect_equal(affcd3$streamEstMean, 3)

        affcd3$streamEstSigma <- 5
        expect_equal(affcd3$streamEstSigma, 5)
        })



#now checking getting/setting BL (derived class!)
test_that("checking BL", {
        eta4 <- 0.01
        alpha4 <- 0.05
        BL <- 93
        affcd4 <- initAFFMeanCD(BL=BL, alpha=alpha4, eta=eta4)

        #default values
        expect_equal(affcd4$BL, 93)

        #now change the BL
        affcd4$BL <- 100
        expect_equal(affcd4$BL, 100)

        })




#now checking pval (getter, only)
test_that("checking pval", {
        eta5 <- 0.01
        alpha5 <- 0.05
        BL <- 50
        affcd5 <- initAFFMeanCD(BL=BL, alpha=alpha5, eta=eta5)

        #default value is 0.5
        expect_equal(affcd5$pval, 0.5)

        })







test_that("lambda is computed properly (ignore changepoints)", {
        #now using burn-in and change detection method, but just monitoring the AFF
        stream3 <- makeStreamMeanChangeR(seednum=5, numChanges=1)
        BL3 <- 50

        lambdatol <- 1e-10

        #aff R:
        #---------AFF Scaled params--------------#
        lambda_init <- 1
        p_saff1 <- 0.99
        resettozero <- 1
        u_init <- 0 
        v_init <- 0
        w_init <- 0
        affmean_init <- 0
        affvar_init <- 0
        low_bound <- 0.6
        up_bound  <- 1
        #should be minus 1, I think
        signchosen <- -1
        #make alpha=0.1 here (eta)
        alpha <- 0.01	
        
        initsafflist <- list(resettozero, u_init, v_init, w_init, affmean_init, affvar_init)
        safflist <- list(low_bound, up_bound, signchosen, alpha)
        affparams <- c(list(lambda_init, p_saff1), initsafflist, safflist)

        #R version:
        returnListR <- AFF_scaled_stream_jumpdetect(stream3, BL3, affparams)
        lambdaR <- returnListR$lambda
        tauR <- returnListR$tau

        #Cpp version:
        alpha3 <- 0.01
        eta3 <- 0.01
        affcd3 <- initAFFMeanCD(alpha=alpha3, eta=eta3, BL=BL3)
        returnListCpp <- affcd3$processVectorSave(stream3)
        lambdaCpp <- returnListCpp$lambda
        tauCpp <- returnListCpp$tauhat
#        vecR <- AFF_scaled_stream_no_change_detection(stream3, BL3, affparams, FALSE)

        #check lambdas and taus are equal
        expect_equal(lambdaCpp, lambdaR, tolerance=lambdatol)
        expect_equal(tauCpp, tauR)

        })




test_that("testing AFF multiple change detection", {
        #now using burn-in and change detection method, but just monitoring the AFF
        stream3 <- makeStreamMeanChangeR(seednum=1, numChanges=5)
        BL3 <- 50


        #aff R:
        #---------AFF Scaled params--------------#
        lambda_init <- 1
        p_saff1 <- 0.99
        resettozero <- 1
        u_init <- 0 
        v_init <- 0
        w_init <- 0
        affmean_init <- 0
        affvar_init <- 0
        low_bound <- 0.6
        up_bound  <- 1
        #should be minus 1, I think
        signchosen <- -1
        #make alpha=0.1 here (eta)
        alpha <- 0.01	
        
        initsafflist <- list(resettozero, u_init, v_init, w_init, affmean_init, affvar_init)
        safflist <- list(low_bound, up_bound, signchosen, alpha)
        affparams <- c(list(lambda_init, p_saff1), initsafflist, safflist)

        #R version:
        returnListR <- AFF_scaled_stream_jumpdetect(stream3, BL3, affparams)
        tauhatR <- returnListR$tau

        #Cpp version:
        alpha3 <- 0.01
        eta3 <- 0.01
        affcd3 <- initAFFMeanCD(alpha=alpha3, eta=eta3, BL=BL3)
        #R version:

        #now cpp:
        returnList <- affcd3$processVectorSave(stream3)
        tauhatCpp <- returnList$tauhat

        expect_equal(tauhatR, tauhatCpp)
        })






test_that("detectMultiple works for AFFMeanCD", {
        stream4 <- makeStreamMeanChangeR(seednum=2, numChanges=5)
        BL4 <- 50


        #aff R:
        #---------AFF Scaled params--------------#
        lambda_init <- 1
        p_saff1 <- 0.99
        resettozero <- 1
        u_init <- 0 
        v_init <- 0
        w_init <- 0
        affmean_init <- 0
        affvar_init <- 0
        low_bound <- 0.6
        up_bound  <- 1
        #should be minus 1, I think
        signchosen <- -1
        #make alpha=0.1 here (eta)
        alpha <- 0.01	
        
        initsafflist <- list(resettozero, u_init, v_init, w_init, affmean_init, affvar_init)
        safflist <- list(low_bound, up_bound, signchosen, alpha)
        affparams <- c(list(lambda_init, p_saff1), initsafflist, safflist)

        #R version:
        returnListR <- AFF_scaled_stream_jumpdetect(stream4, BL4, affparams)
        tauhatR <- returnListR$tau

        #Cpp version:
        alpha4 <- 0.01
        eta4 <- 0.01
        tauhatList <- cpp_detectAFFMeanMultiple(stream4, alpha=alpha4, 
                                                eta=eta4, BL=BL4)

        expect_equal(tauhatList$tauhat, tauhatR)
        #print(tauhatList)

        tauhatList2 <- detectAFFMean(x=stream4, alpha=alpha4, eta=eta4, 
                                     BL=BL4, single=FALSE)

        #now using wrapper in R
        expect_equal(tauhatList2$tauhat, tauhatR)
        
        })




test_that("detectSingle works for AFFMeanCD", {
        stream5 <- makeStreamMeanChangeR(seednum=7, numChanges=2)
        BL5 <- 50


        #aff R:
        #---------AFF Scaled params--------------#
        lambda_init <- 1
        p_saff1 <- 0.99
        resettozero <- 1
        u_init <- 0 
        v_init <- 0
        w_init <- 0
        affmean_init <- 0
        affvar_init <- 0
        low_bound <- 0.6
        up_bound  <- 1
        #should be minus 1, I think
        signchosen <- -1
        #make alpha=0.1 here (eta)
        alpha <- 0.01	
        
        initsafflist <- list(resettozero, u_init, v_init, w_init, affmean_init, affvar_init)
        safflist <- list(low_bound, up_bound, signchosen, alpha)
        affparams <- c(list(lambda_init, p_saff1), initsafflist, safflist)

        #R version:
        returnListR <- AFF_scaled_stream_jumpdetect(stream5, BL5, affparams)
        tauhatR <- returnListR$tau[1]


        #Cpp version:
        alpha5 <- 0.01
        eta5 <- 0.01
        tauhatList <- cpp_detectAFFMeanSingle(stream5, alpha=alpha5, 
                                              eta=eta5, BL5)

        expect_equal(tauhatList$tauhat, tauhatR)
        #print(tauhatList)
        
        tauhatList2 <- detectAFFMean(x=stream5, alpha=alpha5, eta=eta5, single=TRUE, usePrechange=FALSE)

        #now using wrapper in R
        expect_equal(tauhatList2$tauhat, tauhatR)
        })





test_that("detectSingle works for FFFMeanCD with prechange known", {
        #using seednum=7 results in changepoint being found at tauhat=15
        #an example of a lack of burn-in hurting the estimator
        stream6 <- makeStreamMeanChangeR(seednum=8, numChanges=2)
        stream7 <- makeStreamMeanChangeR(seednum=7, numChanges=2)
        BL6 <- 0


        #aff R:
        #---------AFF Scaled params--------------#
        lambda_init <- 1
        p_saff1 <- 0.99
        resettozero <- 1
        u_init <- 0 
        v_init <- 0
        w_init <- 0
        affmean_init <- 0
        affvar_init <- 0
        low_bound <- 0.6
        up_bound  <- 1
        #should be minus 1, I think
        signchosen <- -1
        #make alpha=0.1 here (eta)
        alpha <- 0.01	
        
        initsafflist <- list(resettozero, u_init, v_init, w_init, affmean_init, affvar_init)
        safflist <- list(low_bound, up_bound, signchosen, alpha)
        affparams <- c(list(lambda_init, p_saff1), initsafflist, safflist)

        #Cpp version:
        alpha6 <- 0.01
        eta6 <- 0.01

        prechangeMean <- 0
        prechangeSigma <- 1
        returnList <- AFF_scaled_stream_jumpdetect_prechange(stream6, 0, 
                                     affparams, prechangeMean, prechangeSigma)
        tauhatR <- returnList$tau

        returnList7 <- AFF_scaled_stream_jumpdetect_prechange(stream7, 0, 
                                     affparams, prechangeMean, prechangeSigma)
        tauhatR7 <- returnList7$tau

        #no BL
        affcd6 <- initAFFMeanCD(alpha=alpha6, eta=eta6, BL=0)
        affcd6$streamEstMean  <- prechangeMean
        affcd6$streamEstSigma <- prechangeSigma
        changeNotFound <- TRUE
        index <- 0
        tauhatSingle <- 0
        lambdavec <- rep(0, length(stream6))
        while ( (changeNotFound) & (index < length(stream6)) ) {
            index <- index + 1
            x <- stream6[index]
            affcd6$update(x)
            if (affcd6$changeDetected){
                tauhatSingle <- index
                changeNotFound <- FALSE
            }
            lambdavec[index] <- affcd6$lambda
        }

        tauhatList2 <- detectAFFMean(x=stream6, alpha=alpha6, eta=eta6, multiple=FALSE, usePrechange=TRUE, prechangeMean=prechangeMean, prechangeSigma=prechangeSigma)

        tauhatList7 <- detectAFFMean(x=stream7, alpha=alpha6, eta=eta6, multiple=FALSE, usePrechange=TRUE, prechangeMean=prechangeMean, prechangeSigma=prechangeSigma)


        #now using wrapper in R
        expect_equal(tauhatList2$tauhat, tauhatSingle)
        expect_equal(tauhatList2$tauhat, tauhatR)
        
        expect_equal(tauhatList7$tauhat, tauhatR7)
        })




#TODO: Need to add tests for getLambda and getLderiv to reconcile ideas
#test_that("checking getLderiv for BL=0", {
#        #create aff and affcd objects
#        etaval <- 0.01
#        alphaval <- 0.05
#        BL <- 0
#
#        aff1 <- initAFFMean(eta=etaval)
#        affcd1 <- initAFFMeanCD(eta=etaval, alpha=alphaval)
#        
#
#        #create stream
#        set.seed(1)
#        N <- 50
#        x <- rnorm(2*N) + rep(c(0, 3), each=N)
#
#        lastK <- 60
#        cat("\n\nhello\n")
#        for ( i in seq_along(x) ){
#            obs <- x[i]
#            aff1$update(obs)
#            affcd1$update(obs)
#
#            if (i > N - lastK){
#                cat(i, ": " , aff1$lambda, " ", affcd1$lambda, "\n ")
#            }
#
#        }
#        expect_equal(aff1$lambda, affcd1$lambda)
#        cat("\n")
#        
#
#        })
