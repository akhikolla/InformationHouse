
context("Pstatus_WGEN check")


test_that("Pstatus_WGEN output correctly simulates Markov chain", {
  
  parPwd <- rep(0.4,4)
  parPdd <- rep(0.7,4)
  ndays <- 4
  randomVector <- c(0.08770825, 0.32354876, 0.73957656, 0.69327496)
  
  expect_equal(Pstatus_WGEN(parPwd, parPdd, ndays, randomVector), c(0,0,1,1)) # Compared to hand-calculated transition probabilities
})




# test_that("Pstatus_WGEN correctly signals erroneous or dubious inputs", {
#   
#   expect_error(Pstatus_WGEN(parPwd = c(0.4, 0.4, 0.4, -0.1), 
#                             parPdd = c(0.7, 0.7, 0.7, 0.7), 
#                             ndays = 4, 
#                             randomVector = c(0.08770825, 0.32354876, 0.73957656, 0.69327496)), 
#                "Either parPwd or parPdd \\(or both\\) contain values outside range \\[0,1\\]") 
#   expect_error(Pstatus_WGEN(parPwd = c(0.4, 0.4, 0.4, 0.4), 
#                             parPdd = c(-0.7, 0.7, 0.7, 0.7), 
#                             ndays = 4, 
#                             randomVector = c(0.08770825, 0.32354876, 0.73957656, 0.69327496)), 
#                "Either parPwd or parPdd \\(or both\\) contain values outside range \\[0,1\\]") 
#   expect_error(Pstatus_WGEN(parPwd = c(0.4, 1.4, 0.4, 0.4), 
#                             parPdd = c(0.7, 0.7, 0.7, 0.7), 
#                             ndays = 4, 
#                             randomVector = c(0.08770825, 0.32354876, 0.73957656, 0.69327496)), 
#                "Either parPwd or parPdd \\(or both\\) contain values outside range \\[0,1\\]") 
#   expect_error(Pstatus_WGEN(parPwd = c(0.4, 0.4, 0.4, 0.4), 
#                             parPdd = c(0.7, 1.7, 1.7, 0.7), 
#                             ndays = 4, 
#                             randomVector = c(0.08770825, 0.32354876, 0.73957656, 0.69327496)), 
#                "Either parPwd or parPdd \\(or both\\) contain values outside range \\[0,1\\]") 
#   expect_error(Pstatus_WGEN(parPwd = c(-0.4, 0.4, 0.4, 0.4), 
#                             parPdd = c(0.7, -0.7, 0.7, 0.7), 
#                             ndays = 4, 
#                             randomVector = c(0.08770825, 0.32354876, 0.73957656, 0.69327496)), 
#                "Either parPwd or parPdd \\(or both\\) contain values outside range \\[0,1\\]") 
#   expect_warning(Pstatus_WGEN(parPwd = c(0, 0.4, 0.4, 0.4), 
#                               parPdd = c(0.7, 0.7, 0.7, 0.7), 
#                               ndays = 4, 
#                               randomVector = c(0.08770825, 0.32354876, 0.73957656, 0.69327496)), 
#                "Either parPwd or parPdd \\(or both\\) contain values of either 0 or 1 \\(implausible to occur in real rainfall series\\)")   
#   expect_warning(Pstatus_WGEN(parPwd = c(0.4, 0.4, 0.4, 0.4), 
#                               parPdd = c(0.7, 1, 0.7, 0.7), 
#                               ndays = 4, 
#                               randomVector = c(0.08770825, 0.32354876, 0.73957656, 0.69327496)), 
#                  "Either parPwd or parPdd \\(or both\\) contain values of either 0 or 1 \\(implausible to occur in real rainfall series\\)") 
#   expect_error(Pstatus_WGEN(parPwd = c(0.4, 0.4, 0.4), 
#                             parPdd = c(0.7, 0.7, 0.7, 0.7), 
#                             ndays = 4, 
#                             randomVector = c(0.08770825, 0.32354876, 0.73957656, 0.69327496)), 
#                "parPwd and\\/or parPdd and\\/or randomVector are not of length") 
#   expect_error(Pstatus_WGEN(parPwd = c(0.4, 0.4, 0.4, 0.4), 
#                             parPdd = c(0.7, 0.7, 0.7, 0.7, 0.7), 
#                             ndays = 4, 
#                             randomVector = c(0.08770825, 0.32354876, 0.73957656, 0.69327496)), 
#                "parPwd and\\/or parPdd and\\/or randomVector are not of length") 
#   expect_error(Pstatus_WGEN(parPwd = c(0.4, 0.4, 0.4, 0.4), 
#                             parPdd = c(0.7, 0.7, 0.7, 0.7), 
#                             ndays = 4, 
#                             randomVector = c(0.08770825, 0.32354876, 0.73957656)), 
#                "parPwd and\\/or parPdd and\\/or randomVector are not of length")   
#   })

