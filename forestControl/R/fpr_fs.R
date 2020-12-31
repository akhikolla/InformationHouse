#' False Postivie Rate Feature Selection
#'
#' Calculate the False Positive Rate (FPR) for each feature using it's selection frequency
#'
#' @param x a `randomForest` or `ranger` object
#' @return a `tibble` of selection frequencies and their false positive rate
#'
#' @author Jasen Finch \email{jsf9@@aber.ac.uk}
#' @importFrom purrr map_dbl
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#' @export
#' @examples
#' library(randomForest)
#' data(iris)
#' iris.rf <- randomForest(iris[,-5], iris[,5], forest = TRUE)
#'
#' iris.features <- fpr_fs(iris.rf)
#' print(iris.features)

fpr_fs <- function(x)
{
  params <- extract_params(x)
  freq <- selection_freqs(x)
  fpr <-
    map_dbl(
      unique(freq$freq),
      fpr_fs_calc,
      Ft = params$Ft,
      Fn = params$Fn,
      Tr = params$Tr,
      K = params$K
    ) %>%
    {
      tibble(freq = unique(freq$freq), fpr = .)
    }

  feat_sel <- left_join(freq, fpr, by = "freq")

  return(feat_sel)
}
