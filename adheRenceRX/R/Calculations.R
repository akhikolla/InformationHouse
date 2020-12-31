#' Calculate Proportion Days Covered
#'
#' Calculate the proportion of days covered (PDC) from a pharmaceutical claims database. This function is suggested only  
#' after one has properly adjusted their dates (\code{propagate_date()}) and identified gaps in therapy
#' (\code{identify_gaps()}). This function calculates a length of total therapy as the first fill date to the last for a given grouping. 
#' Finally, if you'd like to view adherence by episodes after you have used \code{rank_episodes()}, the function will re-adjust gaps for you so that the gap that defined the episode isn't included.
#' 
#'
#' @param .data data frame 
#' @param .summarise Logical value (defaulting to TRUE) indicating whether the output should be summarised or not
#' @return a summarised tibble, by default, with proportion of days covered calculated
#' @export
#'
#' @examples
#' library(adheRenceRX)
#' library(dplyr)
#' 
#' toy_claims %>% 
#'   group_by(ID) %>% 
#'   propagate_date(.date = date, .days_supply = days_supply) %>% 
#'   identify_gaps() %>% 
#'   calculate_pdc()
#'   
#' #OR, one could group by the ID and episode of care like...
#' toy_claims %>% 
#'   group_by(ID) %>% 
#'   propagate_date(.date = date, .days_supply = days_supply) %>% 
#'   identify_gaps() %>% 
#'   rank_episodes(.permissible_gap = 30) %>% 
#'   ungroup() %>% 
#'   group_by(ID, episode) %>% 
#'   calculate_pdc()

calculate_pdc <- function(.data, .summarise = TRUE){
  
  
  if(.summarise){
    .data %>% 
      mutate(
        # If one wants to group by episode ranking, the minimum date will have
        # a gap (since the episode was defined by that gap)
        # so, we need to re-adjust gap here to be 0 in this case
        gap = if_else(.data$adjusted_date == min(.data$adjusted_date), 0, .data$gap),
        max_date = max(.data$adjusted_date),
        min_date = min(.data$adjusted_date),
        overlap = as.numeric(.data$max_date - .data$adjusted_date),
        adjusted_gap = if_else(.data$overlap < .data$gap, .data$overlap, .data$gap)
      ) %>% 
      summarise(total_gaps = sum(.data$adjusted_gap),
                total_days = max(as.numeric(.data$max_date-.data$min_date)),
                adherence = 1 - .data$total_gaps/.data$total_days)
  } else {
    .data %>% 
      mutate(
        gap = if_else(.data$adjusted_date == min(.data$adjusted_date), 0, .data$gap),
        max_date = max(.data$adjusted_date),
        min_date = min(.data$adjusted_date),
        overlap = as.numeric(.data$max_date - .data$adjusted_date),
        adjusted_gap = if_else(.data$overlap < .data$gap, .data$overlap, .data$gap)
      ) %>% 
      mutate(total_gaps = sum(.data$adjusted_gap),
             total_days = max(as.numeric(.data$max_date-.data$min_date)),
             adherence = 1 - .data$total_gaps/.data$total_days)
  }
}




