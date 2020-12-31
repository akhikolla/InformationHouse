# Begin Exclude Linting

# function created by Lance Waller as part of his book
# Applied Spatial Statistics for Public Health Data (2004)
cepp <- function(varx = 0,
                 vary = 0,
                 cases = 0,
                 pop = 0,
                 pop.radius = 0,
                 numsim = 9,
                 type = "poisson") {
  # Monte Carlo test of Turnbull et al's Cluster Evaluation Permutation
  # procedure (CEPP) for cases under
  #   constant risk hypothesis (Poisson)
  #   constant risk hypothesis (multinomial)

  dist <- matrix(0, length(varx), length(varx))
  wt.cepp <- matrix(0, length(varx), length(varx))
  turn.stat.obs <- rep(0, length(varx))
  turn.stat.sim <- matrix(999, length(varx), numsim)
  p.val <- rep(0, length(varx))
  max.sim <- rep(0, numsim)
  ind.max.count <- rep(0, length(varx))
  test <- rep(0, length(varx))
  ind <- 1:length(varx)

  for (i in 1:length(varx)) {
    dist[i, ] <- sqrt((varx[i] - varx) ^ 2 + (vary[i] - vary) ^ 2)
    pop.ord <- pop[order(dist[i, ])]
    cas.ord <- cases[order(dist[i, ])]
    ind.ord <- order(dist[i, ])

    cpop <- cumsum(pop.ord)
    ccas <- cumsum(cas.ord)

    frac.ind <- min(ind[cpop > pop.radius])
    if (frac.ind == 1) {
      frac <- pop.radius / pop.ord[frac.ind]
      wt.cepp[i, ind.ord[frac.ind]] <- frac
      turn.stat.obs[i] <- frac * cas.ord[frac.ind]
      test.pop <- frac * pop.ord[frac.ind]
    }
    if (frac.ind > 1) {
      frac <- (pop.radius - cpop[frac.ind - 1]) / pop.ord[frac.ind]
      wt.cepp[i, ind.ord[cpop <= pop.radius]] <- 1
      wt.cepp[i, ind.ord[frac.ind]] <- frac
      test.pop <- max(cpop[cpop <= pop.radius]) +
        frac * pop.ord[frac.ind]
    }

    if (abs(test.pop - pop.radius) > 1e-10) {
      print(paste("test.pop = ", test.pop, " i = ", i))
    }
  }

  turn.stat.obs <- wt.cepp %*% cases

  # Find Monte Carlo p-values

  for (sim in 1:numsim) {
    if (type == "multinomial") {
      test <- runif(sum(cases))
      rand.var <- hist(test,
                       breaks = c(0, (cumsum(cases) / sum(cases))), plot =
                         F)$counts
    }

    if (type == "poisson") {
      rand.var <- rpois(length(cases),
                        lambda = (sum(cases) / sum(pop)) * pop)
    }

    test <- wt.cepp %*% rand.var
    turn.stat.sim[, sim] <- test
    max.sim[sim] <- max(test)
    ind.max.count[ind[test[ind] == max(test)]] <-
      ind.max.count[ind[test[ind] == max(test)]] + 1

  }

  for (i in 1:length(cases)) {
    p.val[i] <-
      length(turn.stat.sim[i, ][turn.stat.sim[i, ] > turn.stat.obs[i]]) /
      (numsim + 1)
  }

  # make histogram
  overall.pval <-
    length(max.sim[max.sim > max(turn.stat.obs)]) / (numsim + 1)

  return(
    list(
      turn.stat.obs = turn.stat.obs,
      turn.stat.sim =
        turn.stat.sim,
      p.val = p.val,
      wt.cepp = wt.cepp,
      overall.pval =
        overall.pval,
      ind.max.count = ind.max.count
    )
  )
}

data(nydf)
data(nyw)
coords = with(nydf, cbind(x, y))
cases = nydf$cases
pop = nydf$pop

cepp1 = cepp(varx = coords[,1], vary = coords[,2],
             cases = cases, pop = pop, pop.radius = 15000,
             numsim = 9, type = "poisson")

cepp2 = suppressWarnings(cepp.test(coords = coords, cases = cases, pop = pop,
                  nstar = 15000, alpha = 1, nsim = 0))

d = sp::spDists(coords)
# find smallest windows with cumulative population of
# at least n* = 15000
nn = casewin(d, pop, 15000)
# compute weights
w = cepp.weights(nn, pop, 15000)

context("check accuracy of cepp.test")
test_that("check accuracy for cepp.test for NY data", {
  for (i in seq_along(w)) {
    wts1i = cepp1$wt.cepp[i, ]
    expect_equal(sort(wts1i[wts1i > 0]), sort(w[[i]]))
  }
  expect_equal(max(cepp1$turn.stat.obs),
               cepp2$clusters[[1]]$test_statistic)
})
# End Exclude Linting
