#' Identify Gaps in Therapy
#'
#' @description Compute gaps in a patient's therapy from the end of their prior fill to the beginning of the next. This function assumes that one has
#' arranged the dates and grouped appropriately outside of the function. The length of any gap will be appended 
#' to the row after the gap has occurred.
#' 
#' @param .data data frame
#' @note This function relies an \code{adjusted_date} column to identify gaps in therapy. So, if you don't want to use \code{propagate_date()} beforehand,
#' you'll need to rename the date variable you wish to use to \code{adjusted_date}.
#'
#' @return A new claims tibble with an appended column, \code{gap}
#' @export
#'
#' @examples 
#' library(adheRenceRX)
#' library(dplyr)
#' 
#' toy_claims %>% 
#'   filter(ID == "D") %>% 
#'   propagate_date(.date_var = date, .days_supply_var = days_supply) %>% 
#'   identify_gaps()
#'
#'
identify_gaps <- function(.data){
  
  .data %>% 
    mutate(
      gap = as.numeric(.data$adjusted_date - lag(.data$adjusted_date)),
      gap = .data$gap - lag(.data$days_supply),
      gap = if_else(is.na(.data$gap), 0, .data$gap)
    )
}

#' Summarise Gaps in Therapy
#'
#' This function serves as a convenience wrapper of \code{dplyr::summarise()}, which takes the grouped variables and
#' summarises their gaps in therapy. This function is to be used after \code{propagate_date()}.
#'
#' @param .data Data to be piped into the function
#' @note This function relies an \code{adjusted_date} column to identify gaps in therapy. So, if you don't want to use \code{propagate_date()} beforehand,
#' you'll need to rename the date variable you wish to use to \code{adjusted_date}.
#' @return A summary of gaps in therapy
#' @export
#'
#' @examples 
#' library(adheRenceRX)
#' library(dplyr)
#' 
#' toy_claims %>% 
#'   filter(ID == "D") %>% 
#'   propagate_date(.date_var = date, .days_supply_var = days_supply) %>% 
#'   summarise_gaps()
#'
#'
#'
summarise_gaps <- function(.data){

  .data %>%
    mutate(gap = as.numeric(.data$adjusted_date - lag(.data$adjusted_date)),
           gap = .data$gap - lag(.data$days_supply)) %>%
    summarise(Sum_Of_Gaps = sum(.data$gap, na.rm = TRUE))
}
