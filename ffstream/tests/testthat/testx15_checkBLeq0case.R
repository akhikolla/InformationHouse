context("Test 15: check BL=0 case for affcd is same as for aff")

test_that("check 1, setting streamEstSigma", {
    #num checks
    numChecks <- 5
    showChecks <- F

    #estimators/change detector
    etaval <- 0.01
    alphaval <- 0.01
    BLval <- 0
    aff1 <- initAFFMean(eta=etaval)
    affcd1 <- initAFFMeanCD(eta=etaval, alpha=alphaval, BL=BLval)
    affcd1$streamEstMean <- 0

    #NOTE: setting streamEstSigma here
    affcd1$streamEstSigma <- 1

    #stream variables
    seednum <- 3
    regimeLength <- 100
    numChanges <- 3
    delta <- 3
    sigma <- 1
    numChanges <- 3
    N <- (regimeLength) * (numChanges + 1)
    stepChanges <- rep(c(0:numChanges), each=regimeLength) * delta * sigma

    set.seed(seednum)
    stream <- rnorm(N, mean=0, sd=sigma)  +  stepChanges 


    #save results
    lambda <- rep(0, length(stream))
    lambdaCD <- rep(0, length(stream))
    lambdaDeriv <- rep(0, length(stream))
    lambdaCDDeriv <- rep(0, length(stream))
    affmean <- rep(0, length(stream))
    affmeanCD <- rep(0, length(stream))
    #derivSign <- rep(0, length(stream))

    #if same as last
    #derivSignNoChange <- rep(0, length(stream))


    t <- seq_along(stream)

    ##NOW ADD GRAPH FOR AFFCD
    for (i in t ){
        obs <- stream[i]
        aff1$update(obs)
        affcd1$update(obs)
        lambda[i] <- aff1$lambda
        lambdaCD[i] <- affcd1$lambda
        affmean[i] <- aff1$xbar
        affmeanCD[i] <- affcd1$affxbar
        lambdaDeriv[i] <- aff1$Lderiv
        lambdaCDDeriv[i] <- affcd1$Lderiv
    }


    #check last values
    endCheck <- length(stream)
    startCheck <- endCheck - numChecks + 1
    for ( j in seq(from=startCheck, to=endCheck, by=1) ){
        expect_equal(lambda[j], lambdaCD[j])
        expect_equal(lambdaDeriv[j], lambdaCDDeriv[j])
        expect_equal(affmean[j], affmeanCD[j])

        if (showChecks){
            cat("\n\n Test1 \n")
            cat(j, ": ", lambda[j], ", ", lambdaCD[j], ", ", 
                lambda[j]==lambdaCD[j],  "\n")
            cat(j, ": ", lambdaDeriv[j], ", ", lambdaCDDeriv[j], ", ", 
                lambdaDeriv[j]==lambdaCDDeriv[j],  "\n")
            cat(j, ": ", affmean[j], ", ", affmeanCD[j], ", ", 
                affmean[j]==affmeanCD[j], "\n")
        }
    }


})


test_that("check 2, NOT setting streamEstSigma", {
    #num checks
    numChecks <- 5
    showChecks <- F

    #estimators/change detector
    etaval <- 0.01
    alphaval <- 0.01
    BLval <- 0
    aff1 <- initAFFMean(eta=etaval)
    affcd1 <- initAFFMeanCD(eta=etaval, alpha=alphaval, BL=BLval)
    affcd1$streamEstMean <- 0

    #NOTE:  NOT setting streamEstSigma 
    ##affcd1$streamEstSigma <- 1


    #stream variables
    seednum <- 3
    regimeLength <- 100
    numChanges <- 3
    delta <- 3
    sigma <- 1
    numChanges <- 3
    N <- (regimeLength) * (numChanges + 1)
    stepChanges <- rep(c(0:numChanges), each=regimeLength) * delta * sigma

    set.seed(seednum)
    stream <- rnorm(N, mean=0, sd=sigma)  +  stepChanges 


    #save results
    lambda <- rep(0, length(stream))
    lambdaCD <- rep(0, length(stream))
    lambdaDeriv <- rep(0, length(stream))
    lambdaCDDeriv <- rep(0, length(stream))
    affmean <- rep(0, length(stream))
    affmeanCD <- rep(0, length(stream))

    #if same as last
    #derivSignNoChange <- rep(0, length(stream))


    t <- seq_along(stream)

    ##NOW ADD GRAPH FOR AFFCD
    for (i in t ){
        obs <- stream[i]
        aff1$update(obs)
        affcd1$update(obs)
        lambda[i] <- aff1$lambda
        lambdaCD[i] <- affcd1$lambda
        affmean[i] <- aff1$xbar
        affmeanCD[i] <- affcd1$affxbar
        lambdaDeriv[i] <- aff1$Lderiv
        lambdaCDDeriv[i] <- affcd1$Lderiv
    }


    #check last values
    endCheck <- length(stream)
    startCheck <- endCheck - numChecks + 1
    for ( j in seq(from=startCheck, to=endCheck, by=1) ){
        expect_equal(lambda[j], lambdaCD[j])
        expect_equal(lambdaDeriv[j], lambdaCDDeriv[j])
        expect_equal(affmean[j], affmeanCD[j])

        if (showChecks){
            cat("\n\n Test2 \n")
            cat(j, ": ", lambda[j], ", ", lambdaCD[j], ", ", 
                lambda[j]==lambdaCD[j],  "\n")
            cat(j, ": ", lambdaDeriv[j], ", ", lambdaCDDeriv[j], ", ", 
                lambdaDeriv[j]==lambdaCDDeriv[j],  "\n")
            cat(j, ": ", affmean[j], ", ", affmeanCD[j], ", ", 
                affmean[j]==affmeanCD[j], "\n")
        }
    }


})

