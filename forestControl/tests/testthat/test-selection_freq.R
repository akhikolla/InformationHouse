context("selection-freqs")

test_that("selection-freqs",{
library(parsnip)

  rf_no_forest <- randomForest::randomForest(iris[,-5], factor(iris$Species), ntree = 100, keep.forest = FALSE)
  set.seed(1234)
  rf_w_forest <- randomForest::randomForest(iris[,-5], factor(iris$Species), ntree = 100, keep.forest = TRUE)

  ranger_no_forest <- ranger::ranger(factor(iris$Species) ~., iris[,-5], write.forest = FALSE,num.trees = 100)
  set.seed(1234)
  ranger_w_forest <- ranger::ranger(factor(iris$Species) ~., iris[,-5], write.forest = TRUE, num.trees = 100)

  model_parsnip_ranger <-
    rand_forest(mode = 'classification') %>% set_engine("ranger", importance = 'impurity') %>%
    fit(Species ~ ., iris)

  model_parsnip_rf <-
    rand_forest(mode = 'classification') %>% set_engine("randomForest") %>%
    fit(Species ~ ., iris)

  expect_true(is.data.frame(selection_freqs(rf_w_forest)))
  expect_true(is.data.frame(selection_freqs(ranger_w_forest)))

  expect_true(is.data.frame(selection_freqs(model_parsnip_ranger)))
  expect_true(is.data.frame(selection_freqs(model_parsnip_rf)))

  expect_error(selection_freqs(rf_no_forest))
  expect_error(selection_freqs(ranger_no_forest))


})
