context("test scan_test related functions")
data(nydf)
cases = floor(nydf$cases)
pop = nydf$pop
ty = sum(cases)
tpop = sum(pop)
ex = ty / tpop * pop

# number of cases and population
# poisson test statistic
zone = c(52, 50, 53, 38, 49, 48, 15, 39, 37, 1, 16, 44,
         47, 40, 14, 2, 51, 13, 43, 45, 17, 55, 11, 3,
         12, 46, 36, 35, 54, 10, 5)
yin = sum(cases[zone])
ein = sum(ex[zone])
yout = ty - yin
eout = ty - ein
t1 = scan.stat(yin = yin, ty = ty, ein = ein)
t1_ = scan_stat(yin = yin, ty = ty, ein = ein)
t2 = stat.poisson(yin = yin, yout = yout, ein = ein, eout = eout)
t2_ = stat_poisson(yin = yin, yout = yout, ein = ein, eout = eout)

zone = c(1, 2, 3, 12, 13, 14, 15, 35, 47, 49)
yin = sum(cases[zone])
popin = sum(pop[zone])
yout = ty - yin
popout = tpop - popin
t3 = scan.stat(yin = yin, ty = ty, popin = popin, tpop = tpop,
               type = "binomial")
t3_ = scan_stat(yin = yin, ty = ty, popin = popin, tpop = tpop,
               type = "binomial")

t4 = stat.binom(yin = yin, yout = yout, ty = ty,
                popin = popin, popout = popout, tpop = tpop)
t4_ = stat_binom(yin = yin, yout = yout, ty = ty,
                popin = popin, popout = popout, tpop = tpop)

test_that("check accuracy for scan_test for NY data", {
  # taken from satscannyoutpoisson
  expect_equal(round(t1, 6), 14.780276)
  expect_equal(round(t2, 6), 14.780276)
  # taken from scan_binomial
  expect_equal(round(t3, 5), 8.47836)
  expect_equal(round(t4, 5), 8.47836)

  # compare stat.poisson, stat_poisson, etc
  expect_equal(t1, t1_)
  expect_equal(t2, t2_)
  expect_equal(t3, t3_)
  expect_equal(t4, t4_)
})
