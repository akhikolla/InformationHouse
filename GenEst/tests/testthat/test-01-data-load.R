context("Check load of data sets")

test_that("Check mock data", {
  data(mock)
  expect_equal(names(mock), c("SE", "CP", "SS", "DWP", "CO"))
  expect_is(mock$SE, "data.frame")
  expect_is(mock$CP, "data.frame")
  expect_is(mock$SS, "data.frame")
  expect_is(mock$DWP, "data.frame")
  expect_is(mock$CO, "data.frame")
})

test_that("Check wind_cleared data", {
  data(wind_cleared)
  expect_equal(names(wind_cleared), c("SE", "CP", "SS", "DWP", "CO"))
  expect_is(wind_cleared$SE, "data.frame")
  expect_is(wind_cleared$CP, "data.frame")
  expect_is(wind_cleared$SS, "data.frame")
  expect_is(wind_cleared$DWP, "data.frame")
  expect_is(wind_cleared$CO, "data.frame")
})

test_that("Check wind_RP data", {
  data(wind_RP)
  expect_equal(names(wind_RP), c("SE", "CP", "SS", "DWP", "CO"))
  expect_is(wind_RP$SE, "data.frame")
  expect_is(wind_RP$CP, "data.frame")
  expect_is(wind_RP$SS, "data.frame")
  expect_is(wind_RP$DWP, "data.frame")
  expect_is(wind_RP$CO, "data.frame")
})

test_that("Check wind_RPbat data", {
  data(wind_RPbat)
  expect_equal(names(wind_RPbat), c("SE", "CP", "SS", "DWP", "CO"))
  expect_is(wind_RPbat$SE, "data.frame")
  expect_is(wind_RPbat$CP, "data.frame")
  expect_is(wind_RPbat$SS, "data.frame")
  expect_is(wind_RPbat$DWP, "data.frame")
  expect_is(wind_RPbat$CO, "data.frame")
})

test_that("Check solar_powerTower data", {
  data(solar_powerTower)
  expect_equal(names(solar_powerTower), c("SE", "CP", "SS", "DWP", "CO"))
  expect_is(solar_powerTower$SE, "data.frame")
  expect_is(solar_powerTower$CP, "data.frame")
  expect_is(solar_powerTower$SS, "data.frame")
  expect_is(solar_powerTower$DWP, "data.frame")
  expect_is(solar_powerTower$CO, "data.frame")
})

test_that("Check solar_PV data", {
  data(solar_PV)
  expect_equal(names(solar_PV), c("SE", "CP", "SS", "DWP", "CO"))
  expect_is(solar_PV$SE, "data.frame")
  expect_is(solar_PV$CP, "data.frame")
  expect_is(solar_PV$SS, "data.frame")
  expect_is(solar_PV$DWP, "data.frame")
  expect_is(solar_PV$CO, "data.frame")
})

test_that("Check solar_trough data", {
  data(solar_trough)
  expect_equal(names(solar_trough), c("SE", "CP", "SS", "DWP", "CO"))
  expect_is(solar_trough$SE, "data.frame")
  expect_is(solar_trough$CP, "data.frame")
  expect_is(solar_trough$SS, "data.frame")
  expect_is(solar_trough$DWP, "data.frame")
  expect_is(solar_trough$CO, "data.frame")
})

