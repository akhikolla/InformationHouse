data(meuse, package = "sp")
sp::coordinates(meuse) = ~ x + y
maxd = max(dist(sp::coordinates(meuse)))/2
# gmeuse = gstat::gstat(id = "cadmium", formula = cadmium ~ 1, data = meuse)
test_that("check args of evgram function", {
  expect_error(evgram(1:10)) # not formula
  f = ~ q + p # bad formula, no response
  expect_error(evgram(f))
  f = j ~ q + p
  expect_error(evgram(f, 1:10)) # not data frame or similar
  expect_error(evgram(f, meuse)) # bad formula
  f = cadmium ~ x + y # good formula
  meuse = as.data.frame(meuse)
  expect_error(evgram(f, meuse)) # need to specify coords
  expect_error(evgram(f, meuse, coords = 1:2)) # bad coords
  expect_error(evgram(f, meuse, coords = ~ x + q)) # bad coords
  cf = ~ x + y # good coords
  # bad nbins
  expect_error(evgram(f, meuse, coords = cf, nbins = -1))
  expect_error(evgram(f, meuse, coords = cf, nbins = 1:2))
  expect_error(evgram(f, meuse, coords = cf, nbins = "b"))
  # bad maxd
  expect_error(evgram(f, meuse, coords = cf, maxd = -1))
  expect_error(evgram(f, meuse, coords = cf, maxd = 1:2))
  expect_error(evgram(f, meuse, coords = cf, maxd = "b"))
  
  # bad angle
  expect_error(evgram(f, meuse, coords = cf, angle = -1))
  expect_error(evgram(f, meuse, coords = cf, angle = 1:2))
  expect_error(evgram(f, meuse, coords = cf, angle = "b"))
  
  # bad ndir
  expect_error(evgram(f, meuse, coords = cf, ndir = 0.99))
  expect_error(evgram(f, meuse, coords = cf, ndir = 1:2))
  expect_error(evgram(f, meuse, coords = cf, ndir = "b"))
  
  # bad type
  expect_error(evgram(f, meuse, coords = cf, type = 1:2))
  expect_error(evgram(f, meuse, coords = cf, type = 1))
  
  # bad npmin
  expect_error(evgram(f, meuse, coords = cf, npmin = 0.99))
  expect_error(evgram(f, meuse, coords = cf, npmin = 1:2))
  expect_error(evgram(f, meuse, coords = cf, npmin = "b"))
})
