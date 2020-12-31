context("wpp functionality")


test_that("wasserstein for balanced wpp objects", {

  balanced_random <- function(nrep, method, p) {
    m <- 30
    n <- 40
    res <- rep(0,nrep)
    for (i in 1:nrep) {
      set.seed(190611+p*100+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      massx <- rexp(m)
      massx <- massx/sum(massx)
      set.seed(190611+p*400+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      massy <- rexp(n)
      massy <- massy/sum(massy)
      set.seed(190611+p*700+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      x <- wpp(matrix(runif(2*m),m,2),massx)
      set.seed(190611+p*1000+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      y <- wpp(matrix(runif(2*n),n,2),massy)
      res[i] <- wasserstein(x,y,method=method,p=p)
    }
    return(res)
  }

  # 20 each, p1 checked for revsimplex, shortsimplex, primaldual (exactly equal!)
  precomp = data.frame(
    p1 = c( 0.2204111080, 0.1844101315, 0.2971009077, 0.1886649427, 0.2790101370, 0.2546042195, 0.1611890814, 0.2074879170,
            0.2688192942, 0.2676081249, 0.2734300136, 0.2042014557, 0.1737081418, 0.4119017624, 0.1586022259, 0.2059735692,
            0.2162859769, 0.1514188963, 0.2699822686, 0.1859990710),
    p2 = c(0.4047509243, 0.2178269822, 0.1914773931, 0.2721772849, 0.2302878681, 0.3246175470, 0.3040419122, 0.2468305838,
           0.1838419748, 0.2363982393, 0.2071055851, 0.2193084666, 0.1995616592, 0.2743842796, 0.1683908060, 0.2665572555,
           0.4027684820, 0.2988371087, 0.2199253236, 0.2574573042)
  )

  expect_equal(balanced_random(nrep=20, method="revsimplex", p=1), precomp$p1)
  expect_equal(balanced_random(nrep=20, method="shortsimplex", p=1), precomp$p1)
  expect_equal(balanced_random(nrep=20, method="primaldual", p=1), precomp$p1)
  expect_equal(balanced_random(nrep=20, method="networkflow", p=1), precomp$p1)
  expect_equal(balanced_random(nrep=20, method="revsimplex", p=2), precomp$p2)
  expect_equal(balanced_random(nrep=20, method="shortsimplex", p=2), precomp$p2)
  expect_equal(balanced_random(nrep=20, method="primaldual", p=2), precomp$p2)
  expect_equal(balanced_random(nrep=20, method="networkflow", p=2), precomp$p2)
})




test_that("wasserstein for unbalanced wpp objects", {
  
  skip_if(getRversion() < "3.6.0", "discrete uniform generation was different for R versions prior to 3.6.0")
  
  unbalanced_random <- function(nrep, method, p) {
    m <- 30
    n <- 40
    res <- rep(0,nrep)
    for (i in 1:nrep) {
      set.seed(190611+p*100+i,kind="Mersenne-Twister",normal.kind = "Inversion",sample.kind="Rejection")
      temp <- sample(30,m,replace=TRUE)
      massx <- runif(m)*10^temp
      massx <- massx/sum(massx)
      set.seed(190611+p*400+i,kind="Mersenne-Twister",normal.kind = "Inversion",sample.kind="Rejection")
      temp <- sample(30,n,replace=TRUE)
      massy <- runif(n)*10^temp
      massy <- massy/sum(massy)
      set.seed(190611+p*700+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      x <- wpp(matrix(runif(2*m),m,2),massx)
      set.seed(190611+p*1000+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      y <- wpp(matrix(runif(2*n),n,2),massy)
      res[i] <- wasserstein(x,y,method=method,p=p)
    }
    return(res)
  }

  # 20 each
  precomp = data.frame(
      p1 = c(0.408712501991695, 0.265729468180307, 0.516829831105873, 0.882572832422254, 
             0.362460109564705, 0.44352849614055, 0.524396418818738, 0.567140446957284, 
             0.199605145672852, 0.500836289468377, 0.608465149250878, 0.526292289496797, 
             0.295525516429156, 0.75368039100222, 0.539846547719169, 0.578900881963893, 
             0.512086614860834, 0.460683637826441, 0.364161062855492, 0.641518533510968),
      p2 = c(0.469550449073627, 0.490510114076552, 0.846518466977827, 0.504892622550457, 
             0.347462537086477, 0.335342245751588, 0.504182118170615, 0.449929012598699, 
             0.311890213777329, 0.225573702901285, 0.260565764246565, 0.372631582632826, 
             0.340233336357178, 0.785096673418506, 0.483406022771683, 0.855022095898563, 
             0.522595036560708, 0.600555445021545, 0.666915983059887, 0.276794690989451)
  )

  # shortlist creates warnings that are no problem, primaldual deprecated
  # --> we don't test these currently (but tests would be ok otherwise)
  expect_equal(unbalanced_random(nrep=20, method="revsimplex", p=1), precomp$p1)
  #expect_equal(unbalanced_random(nrep=20, method="shortsimplex", p=1), precomp$p1)
  #expect_equal(unbalanced_random(nrep=20, method="primaldual", p=1), precomp$p1)
  expect_equal(unbalanced_random(nrep=20, method="revsimplex", p=2), precomp$p2)
  #expect_equal(unbalanced_random(nrep=20, method="shortsimplex", p=2), precomp$p2)
  #expect_equal(unbalanced_random(nrep=20, method="primaldual", p=2), precomp$p2)
})



