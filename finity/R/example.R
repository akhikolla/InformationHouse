.example_finite_moment_test <- function() {
  
  # Generate sample
  rvs <- stabledist::rstable(100000, 1.9, 0.5, 1, 0, pm = 0)
  
  # Perform test
  result <- finite_moment_test(rvs, 2)
  
  # Print results
  message(paste("Test statistic:", result[1], "p-value:", result[2], "\n\n"))
}
