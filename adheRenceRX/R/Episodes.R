#' Ranking Episodes of Care
#' 
#' This is a helper function to assist \code{rank_episodes}
#'
#' @param df a data frame with "gap", "initial_rank", and "permi_gap" columns appended from \code{identify_gaps()}
#' @return a data frame with an "episode" column appended, which ranks episodes of care in time

episode_check <- function(df){
  gap <- df$gap
  epi_rank <- unique(df$initial_rank)
  perm_gap <- unique(df$permi_gap)
  l_df <- nrow(df)
  epi_vec <- c()
  i <- 1
  if(l_df < 2){
    return(df)
  } else {
    epi_vec <- .episode_checkCpp(gap, perm_gap, epi_rank)
    df$episode <- epi_vec
    return(df)
  }
}

#' Rank Episodes of Care
#'
#' This function identifies and labels all episodes of care for a given patient in chronological order. A new episode begins after a specified gap in therapy has occurred. 
#' It is meant to be used after one has appropriately adjusted dates (\code{propagate_date()}) and identified gaps in therapy (\code{identify_gaps()}).
#' 
#' @param .data Data frame with a "gap" column appended from \code{identify_gaps()}
#' @param .permissible_gap Integer value suggesting the maximum gap allowed before labeling a new episode of care
#' @param .initial_rank Integer value to identify what the indexing rank should be (defaults to 1). 
#' @note This function assumes an \code{adjusted_date} column, which is produced by the \code{propagate_date()} function and a
#' \code{gap} column, which is produced by \code{identify_gaps()}. If you would like to rank episodes of care using other dates and a separate
#' column for gaps, you'll need to rename those columns before passing the frame to \code{rank_episodes()}. Notably, this is on purpose as this step should
#' almost always come after the former two.
#' @return The initial claims data frame with an \code{episode} column appended, which ranks episodes of care in time
#' @export
#'
#' @examples 
#' library(adheRenceRX)
#' library(dplyr)
#' 
#' toy_claims %>% 
#'   filter(ID == "D") %>% 
#'   propagate_date() %>% 
#'   identify_gaps() %>% 
#'   rank_episodes(.permissible_gap = 20, .initial_rank = 1)
#'   
rank_episodes <- function(.data, .permissible_gap = NULL, .initial_rank = 1){
  
  # Check to see if user has entered a permissible gap
  if(is.null(.permissible_gap)){
    stop("'.permissible_gap' cannot be NULL; check ?rank_episodes for more information")
  }
  
  perm_gap <- .permissible_gap
  epi_rank <- .initial_rank
  
  
  if(is_grouped_df(.data)){
    # determine the grouping
    groupingVariable <- syms(group_vars(.data))
    
    .data %>% 
      arrange(.data$adjusted_date, .by_group = TRUE) %>%
      # Add columns to pass to episode_check 
      # -- Not as familiar with map below but could probably do that another way
      mutate(permi_gap = perm_gap,
             initial_rank = epi_rank) %>% 
      # Group nest the grouped df
      group_nest() %>% 
      # apply episode_check() to all groups
      mutate(episode = map(.data$data, episode_check)) %>% 
      unnest(.data$episode) %>%
      select(-.data$data) %>% 
      mutate(episode = if_else(is.na(.data$episode), .data$initial_rank, .data$episode)) %>% 
      # regroup by the grouping variables (to prevent constant grouping)
      # most of the time you're going to want to keep these groups
      group_by(!!!groupingVariable) %>% 
      select(-.data$permi_gap, -.data$initial_rank)
      
    
  } else {
    
    .data %>% 
      arrange(.data$adjusted_date, .by_group = TRUE) %>%
      # Add columns to pass to episode_check 
      # -- Not as familiar with map below but could probably do that another way
      mutate(permi_gap = perm_gap,
             initial_rank = epi_rank) %>% 
      nest(data = everything()) %>% 
      # apply episode_check() to all groups
      mutate(episode = map(.data$data, episode_check)) %>% 
      unnest(.data$episode) %>%
      select(-.data$data) %>% 
      mutate(episode = if_else(is.na(.data$episode), .data$initial_rank, .data$episode)) %>% 
      select(-.data$permi_gap, -.data$initial_rank)
    
    
  }
}

  






