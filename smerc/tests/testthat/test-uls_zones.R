context("test uls.zones matches PULSE example")

# Recreate examples from PULSE paper
# PULSE, progressive upper level set scan statistic for
# geospatial hotspot detection
# GP Patil, SW Joshi, RE Koli

# Table 1
cases = c(90, 80, 80, 80, 80, 80, 70, 70, 60, 50, 40, 30)
pop = rep(100, times = 12)
w = matrix(0, nrow = 12, ncol = 12)
ubpop = 1
w[1, c(2, 3, 8)] = 1
w[2, c(1, 3, 7, 4)] = 1
w[3, c(1, 2, 7, 8)] = 1
w[4, c(2, 7)] = 1
w[5, c(7, 6, 8)] = 1
w[6, c(5, 12, 11)] = 1
w[7, c(4, 2, 3, 5)] = 1
w[8, c(1, 3, 5, 11)] = 1
w[9, c(11, 10, 12)] = 1
w[10, c(11, 9)] = 1
w[11, c(8, 6, 9, 10)] = 1
w[12, c(6, 9)] = 1

z = vector("list", 8)
z[[1]] = 1
z[[2]] = 1:4
z[[3]] = 5:6
z[[4]] = 1:8
z[[5]] = 9
z[[6]] = 9:10
z[[7]] = 1:11
z[[8]] = 1:12

uz = uls.zones(cases, pop, w, ubpop = 1, check.unique = TRUE)

# example 9.1.1
cases2 = c(50, 60, 40, 50, 30, 50)
pop2 = rep(100, times = 6)
w2 = matrix(0, nrow = 6, ncol = 6)
ubpop = 1
w2[1, c(4, 2)] = 1
w2[2, c(1, 3, 5)] = 1
w2[3, c(2, 6)] = 1
w2[4, c(1, 5)] = 1
w2[5, c(4, 2, 6)] = 1
w2[6, c(3, 5)] = 1

z2 = vector("list", 5)
z2[[1]] = 2
z2[[2]] = c(1, 2, 4)
z2[[3]] = 6
z2[[4]] = c(1:4, 6)
z2[[5]] = 1:6

uz2 = uls.zones(cases2, pop2, w2, ubpop = 1, check.unique = TRUE)
uz2 = sapply(uz2, sort)

test_that("check accuracy of uls.zones", {
  expect_equal(z[1:4], uz[1:4])
  expect_equal(z2, uz2)
})
