#' A list consisting of Training and Test data along with corresponding class labels.
#' @format A list consisting of:
#' \describe{
#'   \item{TrainData}{ (179 x 4) Matrix of training data features. the first two features satisfy sqrt(x_{i1}^2 + x_{i2}^2) > 2/3 if the ith sample is in class 1. 
#' Otherwise, they satisfy sqrt(x_{i1}^2 + x_{i2}^2) < 2/3 - 1/10 if the ith sample is in class 2. 
#' The third and fourth features are generated as independent N(0, 1/2) noise.}
#'   \item{TestData}{ (94 x 4) Matrix of test data features. the first two features satisfy sqrt(x_{i1}^2 + x_{i2}^2) > 2/3 if the ith sample is in class 1. 
#' Otherwise, they satisfy sqrt(x_{i1}^2 + x_{i2}^2) < 2/3 - 1/10 if the ith sample is in class 2. 
#' The third and fourth features are generated as independent N(0, 1/2) noise.}
#'   \item{CatTrain}{ (179 x 1) Vector of class labels for the training data.}
#'   \item{CatTest}{ (94 x 1) Vector of class labels for the test data.}
#'   ...
#' }
#' @source Simulation model 1 from [Lapanowski and Gaynanova, preprint].
#' @references Lapanowski, Alexander F., and Gaynanova, Irina. ``Sparse Feature Selection in Kernel Discriminant Analysis via Optimal Scoring'', preprint.
"KOS_Data"

#' A list consisting of Training and Test data along with corresponding class labels.
#' @format A list consisting of:
#' \describe{
#'   \item{TrainData}{ (10000 x 10) Matrix of independent normally-distributed training samples conditioned on class membership.
#'   There are 7000 samples belonging to class 1, and 3000 samples belonging to class 2.
#'   The class 1 mean vector is the vector of length 10 consisting only of -2. Likewise,
#'   the class 2 mean vector is the vector of length 10 consisting only of 2.
#'   The shared covariance matrix has (i,j) entry (0.5)^|i-j|.}
#'   \item{TestData}{ (1000 x 10) Matrix of independenttest data features with the same distributions and class proportions as \code{TrainData}}.
#'   \item{Train}{ (10000 x 1) Vector of class labels for the samples in \code{TrainData}.}
#'   \item{TestCat}{ (1000 x 1) Vector of class labels for the samples in \code{TestData}.}
#'   ...
#' }
#' @references Lapanowski, Alexander F., and Gaynanova, Irina. ``Compressing large-sample data for discriminant analysis'', preprint.
"LDA_Data"

#' A list consisting of Training and Test data along with corresponding class labels.
#' @format A list consisting of:
#' \describe{
#'   \item{TrainData}{ (10000 x 10) Matrix of independent normally-distributed training samples conditioned on class membership.
#'   There are 7000 samples belonging to class 1, and 3000 samples belonging to class 2.
#'   The class 1 mean vector is the vector of length 10 consisting only of -2. Likewise,
#'   the class 2 mean vector is the vector of length 10 consisting only of 2.
#'   The class 1 covariance matrix has (i,j) entry (0.5)^|i-j|.
#'   The class 2 covariance matrix has (i,j) entry (-0.5)^|i-j|.}
#'   \item{TestData}{ (1000 x 10) Matrix of independenttest data features with the same distributions and class proportions as \code{TrainData}}.
#'   \item{Train}{ (10000 x 1) Vector of class labels for the samples in \code{TrainData}.}
#'   \item{TestCat}{ (1000 x 1) Vector of class labels for the samples in \code{TestData}.}
#'   ...
#' }
#' @references Lapanowski, Alexander F., and Gaynanova, Irina. ``Compressing large-sample data for discriminant analysis'', preprint.
"QDA_Data"