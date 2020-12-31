# if (requireNamespace("lintr", quietly = TRUE)) {
#   context("lints")
#   mylints = lintr::with_defaults(assignment_linter = NULL,
#                                  closed_curly_linter = NULL,
#                                  commas_linter = NULL,
#                                  object_name_linter = NULL)
#   test_that("Package style", {
#     lintr::expect_lint_free(linters = mylints)
#   })
# }
