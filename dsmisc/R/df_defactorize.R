
#' df_defactorize
#'
#' @param df a data.frame like object 
#'
#' @return returns the same data.frame except that factor columns have been transformed into 
#'         character columns
#' @export
#'
#' @examples
#' 
#' df <- 
#'   data.frame(
#'     a = 1:2, 
#'     b = factor(c("a", "b")), 
#'     c = as.character(letters[3:4]), 
#'     stringsAsFactors = FALSE
#'  )
#' vapply(df, class, "")
#' 
#' df_df <- df_defactorize(df)
#' vapply(df_df, class, "")
#' 
df_defactorize <- function(df){
  iffer <- vapply(df, is.factor, TRUE)
  df[iffer] <- 
    lapply(
      X = df[iffer], 
      FUN = as.character
    ) 
  df
}

