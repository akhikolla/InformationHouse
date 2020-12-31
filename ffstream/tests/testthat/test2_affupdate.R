context("Test  2: AFF; lambda, xbar, Omega, Delta, etc")

#will be used in the tests


test_that("lambda is computed properly", {
        #eta specified
        eta1 <- 0.01
        #aff1 <- new(AFF, eta1)
        aff1 <- initAFFMean(eta=eta1)
        x <- c(1, 2, 3, 4)
        # initial values
        expect_equal(aff1$lambda, 1)
        expect_equal(aff1$Omega, 0)
        expect_equal(aff1$Delta, 0)
        expect_equal(aff1$xbarDeriv, 0)
        expect_equal(aff1$Lderiv, 0)
        #not including this one for display purposes, but it is true
        expect_equal(aff1$xbar, 0)


        aff1$update(x[1])
        #after 1 tick, w_0 = 0, 
        expect_equal(aff1$Omega, 0)
        #after 1 tick, Delta = m_0 = 0
        expect_equal(aff1$Delta, 0)
        Delta1 <- 0
        #xbarDeriv is 0, because Omega and Delta are both 0
        expect_equal(aff1$xbarDeriv, 0)
        #xbarDeriv is 0, so Lderiv is 0
        expect_equal(aff1$Lderiv, 0)
        lderiv1 <- 0
        #no change in lambda, because Lderiv from before was zero
        #actually lambda_1 is computed after Lderiv1
        lambda1 <- 1
        expect_equal(aff1$lambda, lambda1)
        expect_equal(aff1$xbar, x[1])
        


        aff1$update(x[2])
        #no change in lambda, because Lderiv from before was zero
#        expect_equal(aff1$lambda, 1)
        #after 2 ticks, w_1 = 0, 
        #after 2 ticks, Delta = m = x[1]
        # after 2 ticks:
        # m = 1 * x[2] + x[1]
        m2 <- x[1] + x[2]
        # w = 2
        w2 <- 2
        #Delta2 <- 0 + x[1]
        Delta2 <- x[1]
        #Omega2 <- 0 + w1 = 1
        Omega2 <- 1
        xbarderiv2 <- (Delta2*w2 - m2*Omega2)/(w2*w2)
        xbar1 <- x[1]
        xbarderiv1 <- 0
        lderiv2 <- 2*(xbar1 - x[2]) * xbarderiv1
        lambda2 <- 1

        #now can test
        expect_equal(aff1$Lderiv, lderiv2)
        expect_equal(aff1$Delta, Delta2)
        expect_equal(aff1$Omega, Omega2)
        expect_equal(aff1$xbarDeriv, xbarderiv2)
        #still no change for lambda, because Lderiv is zero
        expect_equal(aff1$lambda, lambda2)
        expect_equal(aff1$xbar, m2/w2)


        #first case that uses eta
        aff1$update(x[3])
        #after 3 ticks, w_2 = 2, Omega_2 = 1 so = 3
        w3 <- 3
        Omega3 <- 3
        #after 3 ticks, m_2 = 3 (x[1] + x[2]), Delta_2 = [1] so =4
        m3 <- 6
        Delta3 <- 4
        xbarderiv3 <- (Delta3*w3 - m3*Omega3)/(w3*w3)
        xbar2 <- (1+2)/2
        lderiv3 <- 2*(xbar2 - x[3]) * xbarderiv2 
        #lambda3 <- 1 - eta1 * lderiv3
        lambda3 <- lambda2 - eta1 * lderiv3

        #now can test
        expect_equal(aff1$Lderiv, lderiv3)
        expect_equal(aff1$Omega, Omega3)
        expect_equal(aff1$Delta, Delta3)
        expect_equal(aff1$xbarDeriv, xbarderiv3)
        expect_equal(aff1$lambda, lambda3)
        expect_equal(aff1$xbar, m3/w3)

        #first case that uses eta
        aff1$update(x[4])
        #after 3 ticks, w_2 = 2, Omega_2 = 1 so = 3
        w4 <- w3 * lambda3 + 1 
        Omega4 <- Omega3 * lambda3 + w3 
        #after 3 ticks, m_2 = 3 (x[1] + x[2]), Delta_2 = [1] so =4
        m4 <- m3 * lambda3 + x[4]
        Delta4 <- Delta3 * lambda3 + m3

        xbarderiv4 <- (Delta4*w4 - m4*Omega4)/(w4*w4)
        xbar3 <- m3/w3
        lderiv4 <- 2*(xbar3 - x[4]) * xbarderiv3 
        lambda4 <- lambda3 - eta1 * lderiv4
        xbar4 <- m4/w4

        expect_equal(aff1$xbar, m4/w4)
        expect_equal(aff1$Omega, Omega4)
        expect_equal(aff1$Delta, Delta4)
        expect_equal(aff1$xbarDeriv, xbarderiv4)
        expect_equal(aff1$Lderiv, lderiv4)
        expect_equal(aff1$eta, eta1)
        expect_equal(aff1$lambda, lambda4)

#        cat("lambdas:\n")
#        cat(lambda1, ", ", lambda2, ", ", lambda3, ", ", lambda4, "\n")
#        cat("Deltas: ", Delta1, ", ", Delta2, ", ", Delta3, ", ", Delta4, "\n", sep="")
#        cat("lderivs: ", lderiv1, ", ", lderiv2, ", ", lderiv3, ", ", lderiv4, "\n", sep="")


        })


test_that("lambda is computed properly for long vector, even after change", {

        N <- 10
        stream1 <- c(1:N)


        #aff R:
        #---------AFF Scaled params--------------#
        lambda_init <- 1
        p_saff1 <- 0.995
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

        BL1 <- 0
        vecR <- AFF_scaled_stream_no_change_detection(stream1, BL1, affparams)
        #do not need to set showValues anymore - this is removed
#        vecR <- AFF_scaled_stream_no_change_detection(stream1, BL1, affparams, FALSE)

        aff1 <- initAFFMean(eta=0.01)
        vecCpp <- aff1$processVectorSave(stream1)$lambda
        expect_equal(vecR, vecCpp)



        #now with a longer stream
        stream2 <- makeStreamMeanChangeR(seednum=5, numChanges=1)
        BL2 <- 0
        #do not show values
        vecR2 <- AFF_scaled_stream_no_change_detection(stream2, BL2, affparams)
        
        aff2 <- initAFFMean(eta=0.01)
        vecCpp <- aff2$processVectorSave(stream2)$AFF
        
        })



test_that("processVector works", {
        set.seed(2)
        x <- rnorm(10)
        #two forgetting factors
        aff1 <- initAFFMean(eta=0.01)
        aff2 <- initAFFMean(eta=0.01)

        #sequential
        for (i in seq_len(length(x))){
            obs <- x[i]
            aff1$update(obs)
        }

        #update by vector
        aff2$processVector(x)

        expect_equal(aff1$xbar, aff2$xbar)
        expect_equal(aff1$s2, aff2$s2)

        })
