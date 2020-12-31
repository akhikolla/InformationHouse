#' Restructuring Dates to Remove Overlap
#'
#' This is a function meant to be utilized within propagate_date() in order to adjust pharmaceutical claims dates to prevent
#' overlapping in adherence calculations, per Canfield SL, Zuckerman A, Anguiano RH, Jolly JA, DeClercq J.
#' Navigating the wild west of medication adherence reporting in specialty pharmacy. J Manag Care Spec Pharm. 2019;25(10):1073-77.
#'
#' @param df a claims data frame with a date of a medication claim and the corresponding days supply
#'
#' @import Rcpp
#' @return A new claims data frame with an appended column, \code{adjusted_date}
#'
date_check <- function(df){
  adjusted_date <- df$date
  days_supply <- df$days_supply
  # length to loop over
  l_date <- length(adjusted_date)
  # Ensure multiple rows (else, no need dates to compare to)
  if(l_date < 2){
    return(df)
  } else {
    adjusted_date <- .date_checkCpp(adjusted_date, days_supply)
    df$adjusted_date <- as_date(adjusted_date)
    return(df)
  }
}

#' Adjust Overlapping Fill Dates
#'
#' When assessing pharmaceutical adherence, one should adjust overlapping dates forward for a specified group (e.g. patient ids or medication classes) 
#' so that there is no overlap in days supply. For example, if a patient receives a 30 days supply on January 1st, and another 15 days later, the next fill date
#' should be moved up 15 days. This function is modeled after recommendations from Canfield SL, Zuckerman A, Anguiano RH, Jolly JA, DeClercq J.
#' Navigating the wild west of medication adherence reporting in specialty pharmacy. J Manag Care Spec Pharm. 2019;25(10):1073-77.
#'
#' @param .data Data to be piped into the function
#' @param .date_var Date, column indicating the date of a given fill
#' @param .days_supply_var Integer, column indicating the days supply of a given fill
#' @note This function relies on \link[anytime]{anydate} to parse the users date variable into a date class. So, for most columns passed to .date_var, this function will run without warning or error.
#' For example, \code{anydate(30)} will return "1970-01-31" even though 30 is most likely a days supply. If strange results are produced, double check that the 
#' date variable being specified is indeed a fill date.
#' @rawNamespace import(dplyr, except = data_frame)
#' @importFrom anytime anydate
#' @import tidyr
#' @importFrom purrr map
#' @import lubridate
#' @importFrom rlang enquo !! quo_is_null
#'
#' @return The initial claims data frame with an appended column,  \code{adjusted_date}
#' @export
#'
#' @examples 
#' library(adheRenceRX)
#' library(dplyr)
#' 
#' toy_claims %>% 
#'   filter(ID == "D") %>% 
#'   propagate_date(.date_var = date, .days_supply_var = days_supply)

propagate_date <- function(.data, .date_var = NULL, .days_supply_var = NULL){

  
  
  # naming 'date_var' and 'days_supply_var' for consistent application
  .date <- enquo(.date_var)
  .days_supply <- enquo(.days_supply_var)
  
  
  # Determine whether data is grouped to regroup at the end
  if(is_grouped_df(.data)){
    
    # determine the grouping
    groupingVariable <- syms(group_vars(.data))
    
    .data %>%
      # rename data appropriately to standardize with the rest of the functions
      # will consider more flexible applications at some point
      rename(date = !!.date,
             days_supply = !!.days_supply) %>%
      # guessing date if given in weird format
      mutate(date = anydate(.data$date)) %>% 
      arrange(date, .by_group = TRUE) %>%
      group_nest() %>% 
      # apply date_check() to all groups
      mutate(propagated_date = map(.data$data, date_check)) %>%
      unnest(.data$propagated_date) %>%
      select(-.data$data) %>% 
      mutate(adjusted_date = if_else(is.na(.data$adjusted_date), .data$date, .data$adjusted_date)) %>% 
      # regroup by the grouping variables (to prevent constant grouping)
      # most of the time you're going to want to keep these groups
      group_by(!!!groupingVariable)
    
  } else {
    
    .data %>%
      rename(date = !!.date,
             days_supply = !!.days_supply) %>%
      mutate(date = anydate(.data$date)) %>% 
      arrange(date, .by_group = TRUE) %>%
      group_nest() %>% 
      mutate(propagated_date = map(.data$data, date_check)) %>%
      unnest(.data$propagated_date) %>%
      select(-.data$data) %>% 
      mutate(adjusted_date = if_else(is.na(.data$adjusted_date), .data$date, .data$adjusted_date))
  
  }
  

}
