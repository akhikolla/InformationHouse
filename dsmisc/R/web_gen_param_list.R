#' URL Parameter Combinations
#' 
#' Generate URL parameter combinations from sets of parameter values.
#'
#' @param ... multiple vectors passed on as named arguments or a single list or a data.frame 
#' @param sep_1 first separator to use between key and value 
#' @param sep_2 second separator to use between key-value pairs
#'
#' @return string vector with assembled query string parameter combinations
#'
#' @export
#'
#' @examples
#' 
#' web_gen_param_list_expand(q = "beluga", lang = c("de", "en"))
#' 
web_gen_param_list_expand <- 
  function(..., sep_1 = "=", sep_2 = "&"){
    
    # handling edge cases
    lst <- list(...)
    lst <- 
      lapply(
        lst, 
        function(x){
          if ( length(x)==0 ){
            x <- ""
          }
          x
        }
      )
    
    # make combinations
    df <- expand.grid(lst)
    
    # concatenation 1
    for(i in seq_along(df)){
      concat  <- paste0(names(df)[i], sep_1, as.character(df[[i]]))
      df[[i]] <- concat
    }
    
    # concatenation 2
    df <- apply(df, 1, paste0, collapse = sep_2)
    
    # return
    unlist(df)
  }