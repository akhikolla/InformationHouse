library(testthat)
library(ldsr)
#
# force tests to be executed if in dev release which we define as
## having a sub-release, eg 0.9.15.5 is one whereas 0.9.16 is not
# if (length(strsplit(packageDescription("ldsr")$Version, "\\.")[[1]]) > 3) {
#   Sys.setenv("RunAllTests" = "yes")
# }
test_check("ldsr")
