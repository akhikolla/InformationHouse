context("check scan.sim accuracy for different distributions")
set.seed(2)
nsim = 499
data(nydf)
coords = nydf[, c("x", "y")]
nn = nnpop(as.matrix(dist(coords)), pop = nydf$pop, ubpop = 0.1)
cases = floor(nydf$cases)
ty = sum(cases)
e = ty / sum(nydf$population) * nydf$population
ein = nn.cumsum(nn, e)
tpop = sum(nydf$population)
popin = nn.cumsum(nn, nydf$population)
sa = scan.sim(nsim, nn, ty = ty, ex = e, type = "poisson",
              ein = ein, eout = ty - ein, simdist = "multinomial",
              pop = nydf$pop)

sb = scan.sim(nsim, nn, ty = ty, ex = e, type = "poisson",
              ein = ein, eout = ty - ein, simdist = "poisson",
              pop = nydf$pop)

sc = scan.sim(nsim, nn, ty = ty, ex = e, type = "binomial",
              ein = ein, eout = ty - ein, simdist = "binomial",
              tpop = tpop, popin = popin,
              popout = tpop - popin, pop = nydf$pop)
summa = summary(sa)
summb = summary(sb)
summc = summary(sc)

test_that("check accuracy for scan.sim", {
  expect_true(round(summa[2], 1) - round(summb[2], 1) <= 0.1)
  expect_true(round(summb[2], 1) - round(summc[2], 1) <= 0.1)
  expect_true(round(summa[3], 1) - round(summb[3], 1) <= 0.1)
  expect_true(round(summb[3], 1) - round(summc[3], 1) <= 0.1)
  expect_true(round(summa[4], 1) - round(summb[4], 1) <= 0.1)
  expect_true(round(summb[4], 1) - round(summc[4], 1) <= 0.1)
})
