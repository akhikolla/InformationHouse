context("Path values")
library(gretel)


path1 <- c(1,2,3)
path2 <- c(1,5,2,3)
path3 <- c(1,5,4,3)
not_a_path <- c(1,5,3)
junk <- matrix(c(1,"a",3,3,2,1), nrow = 2)

test_that("Historical measures recover Yang, Knoke (2001)", {
  expect_equal(binary_distance(YangKnoke01, path1), 2)
  expect_equal(binary_distance(YangKnoke01, path2), 3)
  expect_equal(binary_distance(YangKnoke01, path3), 3)
  
  expect_equal(peay_path_value(YangKnoke01, path1), 1)
  expect_equal(peay_path_value(YangKnoke01, path2), 2)
  expect_equal(peay_path_value(YangKnoke01, path3), 3)
  
  expect_equal(flament_path_length(YangKnoke01, path1), 3)
  expect_equal(flament_path_length(YangKnoke01, path2), 8)
  expect_equal(flament_path_length(YangKnoke01, path3), 13)
  
  expect_equal(peay_average_path_value(YangKnoke01, path2), 2/3)
  expect_equal(flament_average_path_length(YangKnoke01, path2), 8/3)
})

test_that("Historical measures detect bad inputs", {
  # Not a real path
  # Not the right format
  expect_equal(binary_distance(YangKnoke01, path = not_a_path), Inf)
  expect_error(binary_distance(YangKnoke01, path = junk))
  
  expect_equal(peay_path_value(YangKnoke01, path = not_a_path), 0)
  expect_error(peay_path_value(YangKnoke01, path = junk))
  
  expect_equal(flament_path_length(YangKnoke01, path = not_a_path), NA)
  expect_error(flament_path_length(YangKnoke01, path = junk))
})

test_that("'gpv' and 'ppv' evaluate correctly", {
  #two values of p
  expect_equal(gpv(YangKnoke01, path = path1, p = Inf), 1)
  expect_equal(gpv(YangKnoke01, path = path1, p = 1), 2/3)
  #one value of p with node costs
  expect_equal(gpv(YangKnoke01, path = path1, p = 1, node_costs = c(1,2,1,2,1)), 0.2857143)
  #two values of odds_scale
  expect_equal(ppv(YangKnoke01, path = path1, odds_scale = 1), 1/2)
  expect_equal(ppv(YangKnoke01, path = path1, odds_scale = 2), 4/10)
  #odds scale by node
  expect_equal(ppv(YangKnoke01, path = path1, odds_scale_by_node = c(2,1,1,2,3)), 0.2857143)
})

test_that("'gpv' and 'ppv' detect bad inputs", {
  # Not the right format path
  expect_error(ppv(YangKnoke01, path = junk))
  # not the right format p
  expect_error(gpv(YangKnoke01, path = path1, p = -1))
  expect_error(gpv(YangKnoke01, path = path1, p = "a"))
  # not the right format node costs
  expect_error(gpv(YangKnoke01, path = path1, node_costs = c(-5,2,3,1,2)))
  # not the right format odds_scale
  expect_error(ppv(YangKnoke01, path = path1, odds_scale = "a"))
  # not the right format odds_scale_by_node
  expect_error(ppv(YangKnoke01, path = path1 , odds_scale_by_node = c(1,2)))
  # Not a real path
  expect_equal(gpv(YangKnoke01, path = not_a_path, p = 2), 0)
  expect_equal(ppv(YangKnoke01, path = not_a_path, odds_scale = 3), 0)
})