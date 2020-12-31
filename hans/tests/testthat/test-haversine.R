test_that("haversine works on dataframe", {

  set.seed(42)
  lon1 <- runif(-160, -60, n = 3)
  lat1 <- runif(40, 60, n = 3)
  lon2 <- runif(-160, -60, n = 3)
  lat2 <- runif(40, 60, n = 3)

  df <- data.frame(lat1, lon1, lat2, lon2)

  df$havers <- haversine(df$lat1, df$lon1, df$lat2, df$lon2)

  rows <- nrow(df)

  expect_equal(rows, 3)
})
