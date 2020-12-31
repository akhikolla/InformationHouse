
test_that(
  "simple graphs work",
  {
    expect_true(
      graphs_find_subgraphs(1,1, verbose = 0) == 1
    )
    
    expect_true(
      all(
        graphs_find_subgraphs(1:3,1:3, verbose = 0) == 1:3
      )
    )
    
    expect_true(
      all(
        graphs_find_subgraphs( c(33:30, 1:5), c(34:31, 2:6), verbose = 0) %in% 1:2 
      ) & 
      all(
        1:2 %in% graphs_find_subgraphs( c(33:30, 1:5), c(34:31, 2:6), verbose = 0)
      )
    )
  }
)