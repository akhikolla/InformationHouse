test_that("Proper frame for eval.parent(sexpr[[2]])", {
  # We need to make sure that when and3s enters the dots, the symbols it uses
  # on the second run are the symbols provided by the user, not the symbols
  # in the frame of the function (e.g. ox)

  ox <- oy <- oz <- 2L  # name of symbol conflicts with op2M(">") internally
  y <- 1:5
  B <- C <- D <- (y > 3L)
  A <- !B



  and4b <- function(a, b, c, d) a & b & c & d
   or4b <- function(a, b, c, d) a | b | c | d

  expect_equal(and3s(y > 0, y > 0, y > 0, y > oy),
               and4b(y > 0, y > 0, y > 0, y > oy))

  expect_equal(and3s(y > 0, y > 0, A, y > oy),
               and4b(y > 0, y > 0, A, y > oy))

  expect_equal( or3s(y < 0, y < 0, A, y < oy),
                or4b(y < 0, y < 0, A, y < oy))

  expect_equal(and3s(y > 0, y > 0, A, B),
               and4b(y > 0, y > 0, A, B))

  expect_equal( or3s(y < 0, y < 0, A, B),
                or4b(y < 0, y < 0, A, B))

  expect_equal(and3s(y > 0, y > 0, A, !B),
               and4b(y > 0, y > 0, A, !B))

  expect_equal( or3s(y < 0, y < 0, A, !B),
                or4b(y < 0, y < 0, A, !B))

  test_multi <-
    vapply(0:9, function(i) {
      ox <- i %% 5L
      oy <- i %/% 5L
      oz <- 3L
      A <- 1:5
      !identical(and3s(y >= 0, y >= 0, y >= 0, y >= oy),
                 and4b(y >= 0, y >= 0, y >= 0, y >= oy))
    },
    logical(1))

  expect_equal(which_first(test_multi), 0L)
})
