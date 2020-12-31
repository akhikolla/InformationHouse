context("testcreate")

test_that("CGGPcreate works", {
  expect_error(SG <- CGGPcreate())
  expect_error(SG <- CGGPcreate(d=3))
  expect_error(SG <- CGGPcreate(batchsize=20))
  expect_error(SG <- CGGPcreate(d=1, batchsize=20))
  
  
  SG <- CGGPcreate(d=3, batchsize=20, corr = "MATern32")
  expect_is(SG, "list")
  expect_is(SG, "CGGP")
  # Can give in function directly
  expect_error(CGGPcreate(3, 30, corr=CGGP_internal_CorrMatMatern32), NA)
  # Get error if bad string for corr
  expect_error(CGGPcreate(3, 30, corr="notrealcorrfunc"))
  
  # Create a big one that forces it to allocate more memory
  expect_error(CGGPcreate(d=4, 7000, corr="powerEXP"), NA)
})
