qqplot_masterList <- function() {
  
  variables <- c("P", "PET", "Temp", "Radn", "wetDay", "1dayMaxP", "2dayMaxP")
  
  units <- list("mm",  "mm",   "K",  "W/sq.m", "frac",  "mm/day",   "mm/day")
  
  # Is the variable calculated from P/PET/Temp/Radn?
  calcFrom <- list(NULL, NULL, NULL,   NULL,     "P",       "P",       "P")
  
  # Calculation functions: How should the variable be calculated? (eg: identify wet days, take a moving average)
  calc_FUN_list <- list(NULL, NULL, NULL, NULL, wetdays,  NULL, ma2)
  
  # Aggregation functions: How should the data for the year be aggregated?
  agg_FUN_list <- list(list(FUN = sum, na.rm = FALSE), 
                       list(FUN = sum, na.rm = FALSE), 
                       list(FUN = mean, na.rm = FALSE), 
                       list(FUN = mean, na.rm = FALSE), 
                       list(FUN = mean, na.rm = FALSE), 
                       list(FUN = max, na.rm = FALSE),
                       list(FUN = max, na.rm = TRUE)
  )
  
  # CDF Functions: normal or extreme?
  cdf_FUN_list <- list(qnorm, qnorm, qnorm, qnorm, qnorm, evd::qgumbel, evd::qgumbel)
  
  qqplot_labels <- list("Annual Total P",
                        "Annual Total PET",
                        "Annual Mean Temp",
                        "Annual Mean Radn",
                        "Annual Wet Day Proportion",
                        "1-Day Annual Maximum P",
                        "2-Day Annual Maximum P"
  )
  
  names(units) <- names(calcFrom) <- names(calc_FUN_list) <- names(agg_FUN_list) <- names(cdf_FUN_list) <- names(qqplot_labels) <- variables
  
  return(list(variables = variables,
              units = units,
              calcFrom = calcFrom,
              calc_FUN_list = calc_FUN_list,
              agg_FUN_list = agg_FUN_list,
              cdf_FUN_list = cdf_FUN_list,
              qqplot_labels = qqplot_labels))
  
}
