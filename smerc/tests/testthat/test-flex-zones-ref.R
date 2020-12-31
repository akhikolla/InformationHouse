fpath = system.file("testdata",  package = "smerc")
fname = paste(fpath, "/flex5_zones_ref.rda", sep = "")
load(fname)

data(nydf)
coords = as.matrix(nydf[, c("x", "y")])
data(nyw)
flex5_zones = flex.zones(coords, nyw, k = 5)

context("check flex5_zones_ref w/ flex5_zones")
test_that("flex.zones result matches reference", {
  expect_equal(flex5_zones_ref, flex5_zones)
})
