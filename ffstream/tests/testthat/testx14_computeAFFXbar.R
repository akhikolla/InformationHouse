context("Test 14: direct vector computation of AFF xbar")

test_that("make computeAFFMean produce errors", {
        #check for non-numerics
        x1 <- c(1, "a", 3)
        #interestingly, it just matches a regular expression, 
        #so does not need to be exact
        expect_error(computeAFFMean(x1), "x contains non-numeric")

        #check for AN
        x2 <- c(1, NA, 3)
        expect_error(computeAFFMean(x2), "contains NA")

        #check for NaNs
        x3 <- c(1, NaN, 3)
        expect_error(computeAFFMean(x3), "contains NA")

        #check for Inf
        x4 <- c(1, Inf, 4)
        expect_error(computeAFFMean(x4), "contains non-finite")

        #now check lambda
        x5 <- c(1, 2, 4)
        expect_error(computeAFFMean(x5, eta=NA), "not a finite")
        expect_error(computeAFFMean(x5, eta=NULL), "is NULL")
        expect_error(computeAFFMean(x5, eta=0), NA) #success
        expect_error(computeAFFMean(x5, eta=1), NA) #success
        expect_error(computeAFFMean(x5, eta=0.9), NA) #success
        expect_error(computeAFFMean(x5, eta=-0.1), "not in range") 
        expect_error(computeAFFMean(x5, eta=1.1), "not in range") 
        
        })


#unnecessary, really, since test 2 already covers this
test_that("check computeAFFMean works", {
        #eta specified
        eta1 <- 0.01
        #aff1 <- new(AFF, eta1)
        x <- c(1, 2, 3, 4)


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
        xbar2 <- m2/w2

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

        #xbarvec <- c(xbar1, xbar2, xbar3, xbar4)

        expect_equal(computeAFFMean(x, eta=eta1), xbar4)

        })

