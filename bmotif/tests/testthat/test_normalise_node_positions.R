context("normalise_node_positions")

pc <- matrix(1, 10, 148, dimnames = list(c(paste0("r",1:5),paste0("c",1:5)),paste0("np",1:148)))

test_that("'none' doesn't change anything",{
  nnp <- normalise_node_positions(pc = pc, type = "none", six_node = TRUE)
  expect_identical(pc, nnp)
})

test_that("'sum' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "sum", six_node = TRUE)
  rsnpp <- sapply(rowSums(nnp), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
})

test_that("'position' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "position", six_node = TRUE)
  csnpp <- sapply(colSums(nnp), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(csnpp), TRUE)
})

test_that("'sizeclass' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "sizeclass", six_node = TRUE)
  rsnpp <- sapply(rowSums(nnp[,1:2]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,3:6]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,7:16]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,17:46]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,47:148]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
})

test_that("'levelsize' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "levelsize", six_node = TRUE)
  rsnpp <- sapply(rowSums(nnp[,1:2]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,3:4]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,5:6]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,7:8]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,9:14]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,15:16]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,17:18]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,19:31]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,32:44]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,45:46]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,47:48]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,49:70]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,71:124]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,125:146]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
  rsnpp <- sapply(rowSums(nnp[,147:148]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
  expect_identical(all(rsnpp), TRUE)
})

test_that("'motif' works correctly when pc is all 1s",{
  nnp <- normalise_node_positions(pc = pc, type = "motif", six_node = TRUE)
  mps <- lapply(1:44, function(x) motif_info(x, link = FALSE)) # motif positions
  for(i in mps){
    rsnpp <- sapply(rowSums(nnp[,i]), function(x) all.equal(x, 1, tolerance = sqrt(.Machine$double.eps)))
    expect_identical(all(rsnpp), TRUE)
  }
})

### ------ Testing _NAzero
test_that("'sizeclass_NAzero' replaces all NAs/NaNs with zero",{
  for(i in 1:10){
    m <- motifs[[sample(17,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "sizeclass", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "sizeclass_NAzero", six_node = TRUE)
    expect_identical(all(np_NAzero[is.na(np_normal)] == 0), TRUE)
  }
})

test_that("'levelsize_NAzero' replaces all NAs/NaNs with zero",{
  for(i in 1:10){
    m <- motifs[[sample(44,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "levelsize", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "levelsize_NAzero", six_node = TRUE)
    expect_identical(all(np_NAzero[is.na(np_normal)] == 0), TRUE)
  }
})

test_that("'motif_NAzero' replaces all NAs/NaNs with zero",{
  for(i in 1:10){
    m <- motifs[[sample(44,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "motif", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "motif_NAzero", six_node = TRUE)
    expect_identical(all(np_NAzero[is.na(np_normal)] == 0), TRUE)
  }
})

### ------ Testing _plus1
test_that("'sizeclass_plus1' removes all NAs/NaNs",{
  for(i in 1:10){
    m <- motifs[[sample(17,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "sizeclass", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "sizeclass_plus1", six_node = TRUE)
    expect_identical(all(!is.na(np_NAzero[is.na(np_normal)])), TRUE)
  }
})

test_that("'levelsize_plus1' removes all NAs/NaNs",{
  for(i in 1:10){
    m <- motifs[[sample(44,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "levelsize", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "levelsize_plus1", six_node = TRUE)
    expect_identical(all(!is.na(np_NAzero[is.na(np_normal)])), TRUE)
  }
})

test_that("'motif_plus1' removes all NAs/NaNs",{
  for(i in 1:10){
    m <- motifs[[sample(44,size = 1)]]
    np <- node_positions(M = m, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    np_normal <- normalise_node_positions(np, type = "motif", six_node = TRUE)
    np_NAzero <- normalise_node_positions(np, type = "motif_plus1", six_node = TRUE)
    expect_identical(all(!is.na(np_NAzero[is.na(np_normal)])), TRUE)
  }
})

### ------ Make sure we never have NaNs
test_that("Check that we never have NaNs, sum-normalisation", {
  M <- matrix(1,1,1)
  npc <- node_positions(M, normalisation = "none", six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none")
  np <- normalise_node_positions(npc, type = "sum", six_node = TRUE)
  l <- lapply(np, is.nan)
  for (item in l) {
    expect_true(!any(item))
  }
})

test_that("Check that we never have NaNs, sizeclass-normalisation", {
  M <- matrix(1,1,1)
  npc <- node_positions(M, normalisation = "none", six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none")
  np <- normalise_node_positions(npc, type = "sizeclass", six_node = TRUE)
  l <- lapply(np, is.nan)
  for (item in l) {
    expect_true(!any(item))
  }
})

test_that("Check that we never have NaNs, levelsize-normalisation", {
  M <- matrix(1,1,1)
  npc <- node_positions(M, normalisation = "none", six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none")
  np <- normalise_node_positions(npc, type = "levelsize", six_node = TRUE)
  l <- lapply(np, is.nan)
  for (item in l) {
    expect_true(!any(item))
  }
})

test_that("Check that we never have NaNs, motif-normalisation", {
  M <- matrix(1,1,1)
  npc <- node_positions(M, normalisation = "none", six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none")
  np <- normalise_node_positions(npc, type = "motif", six_node = TRUE)
  l <- lapply(np, is.nan)
  for (item in l) {
    expect_true(!any(item))
  }
})

