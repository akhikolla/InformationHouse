## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("superml")

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github("saraswatmks/superml")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("superml", dependencies=TRUE)

## -----------------------------------------------------------------------------
library(superml)

# should be a vector of texts
sents <-  c('i am going home and home',
          'where are you going.? //// ',
          'how does it work',
          'transform your work and go work again',
          'home is where you go from to work')

# generate more sentences
n <- 10
sents <- rep(sents, n) 
length(sents)

## -----------------------------------------------------------------------------
# initialise the class
tfv <- TfIdfVectorizer$new(max_features = 10, remove_stopwords = FALSE)

# generate the matrix
tf_mat <- tfv$fit_transform(sents)

head(tf_mat, 3)


## -----------------------------------------------------------------------------
# initialise the class
tfv <- TfIdfVectorizer$new(min_df = 0.4, remove_stopwords = FALSE, ngram_range = c(1, 3))

# generate the matrix
tf_mat <- tfv$fit_transform(sents)

head(tf_mat, 3)


## ---- warning=FALSE-----------------------------------------------------------

library(data.table)
library(superml)

# use sents from above
sents <-  c('i am going home and home',
          'where are you going.? //// ',
          'how does it work',
          'transform your work and go work again',
          'home is where you go from to work',
          'how does it work')

# create dummy data
train <- data.table(text = sents, target = rep(c(0,1), 3))
test <- data.table(text = sample(sents), target = rep(c(0,1), 3))

## -----------------------------------------------------------------------------
head(train, 3)


## -----------------------------------------------------------------------------
head(test, 3)

## -----------------------------------------------------------------------------
# initialise the class
tfv <- TfIdfVectorizer$new(min_df = 0.3, remove_stopwords = FALSE, ngram_range = c(1,3))

# we fit on train data
tfv$fit(train$text)

train_tf_features <- tfv$transform(train$text)
test_tf_features <- tfv$transform(test$text)

dim(train_tf_features)
dim(test_tf_features)


## -----------------------------------------------------------------------------
head(train_tf_features, 3)

## -----------------------------------------------------------------------------
head(test_tf_features, 3)

## -----------------------------------------------------------------------------

# ensure the input to classifier is a data.table or data.frame object
x_train <- data.table(cbind(train_tf_features, target = train$target))
x_test <- data.table(test_tf_features)


xgb <- XGBTrainer$new(n_estimators = 10, objective = "binary:logistic")
xgb$fit(x_train, "target")

predictions <- xgb$predict(x_test)
predictions


