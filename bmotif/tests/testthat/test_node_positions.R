context("node_positions")

test_that("In a motif, each node is only in one position",{
  mp <- lapply(1:44, function(x) motif_info(x, link = FALSE)) # motif positions
  for(i in 1:length(mp)){
    ps <- mp[[i]]
    mot <- motifs[[i]]
    nps <- node_positions(M = mot, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    nps <- nps[,ps]
    expect_identical(all(rowSums(nps) == 1), TRUE)
  }

  # tests with block matrices
  mlist2 <- motifs[1]
  bm2 <- block_matrix(mlist2)
  np <- node_positions(M = bm2, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,1:2]
  expect_identical(all(rowSums(np) == 1), TRUE)

  mlist3 <- motifs[2:3]
  bm3 <- block_matrix(mlist3)
  np <- node_positions(M = bm3, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,3:6]
  expect_identical(all(rowSums(np) == 1), TRUE)

  mlist4 <- motifs[4:7]
  bm4 <- block_matrix(mlist4)
  np <- node_positions(M = bm4, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,7:16]
  expect_identical(all(rowSums(np) == 1), TRUE)

  mlist5 <- motifs[8:17]
  bm5 <- block_matrix(mlist5)
  np <- node_positions(M = bm5, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,17:46]
  expect_identical(all(rowSums(np) == 1), TRUE)

  mlist6 <- motifs[18:44]
  bm6 <- block_matrix(mlist6)
  np <- node_positions(M = bm6, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")[,47:148]
  expect_identical(all(rowSums(np) == 1), TRUE)
})

test_that("Test counts for motif 1", {
  np <- node_positions(M = motifs[[1]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,1], c(0,1))
  expect_equal(np[,2], c(1,0))
  expect_identical(all(np[,3:148] == 0), TRUE)
})

test_that("Test counts for motif 2", {
  np <- node_positions(M = motifs[[2]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,1], c(0,1,1))
  expect_equal(np[,2], c(2,0,0))
  expect_equal(np[,3], c(0,1,1))
  expect_equal(np[,4], c(1,0,0))
  expect_identical(all(np[,5:148] == 0), TRUE)
})

test_that("Test counts for motif 3", {
  np <- node_positions(M = motifs[[3]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,1], c(0,0,2))
  expect_equal(np[,2], c(1,1,0))
  expect_equal(np[,3], c(0,0,0))
  expect_equal(np[,4], c(0,0,0))
  expect_equal(np[,5], c(0,0,1))
  expect_equal(np[,6], c(1,1,0))
  expect_identical(all(np[,7:148] == 0), TRUE)
})

test_that("Test counts for motif 4", {
  np <- node_positions(M = motifs[[4]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,1], c(0,0,0,3))
  expect_equal(np[,2], c(1,1,1,0))
  expect_equal(np[,3], c(0,0,0,0))
  expect_equal(np[,4], c(0,0,0,0))
  expect_equal(np[,5], c(0,0,0,3))
  expect_equal(np[,6], c(2,2,2,0))
  expect_equal(np[,7], c(0,0,0,1))
  expect_equal(np[,8], c(1,1,1,0))
  expect_identical(all(np[,9:148] == 0), TRUE)
})

test_that("Test counts for motif 5", {
  np <- node_positions(M = motifs[[5]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,1], c(0,0,1,2))
  expect_equal(np[,2], c(2,1,0,0))
  expect_equal(np[,3], c(0,0,1,1))
  expect_equal(np[,4], c(1,0,0,0))
  expect_equal(np[,5], c(0,0,0,1))
  expect_equal(np[,6], c(1,1,0,0))
  expect_equal(np[,7], c(0,0,0,0))
  expect_equal(np[,8], c(0,0,0,0))
  expect_equal(np[,9], c(0,0,1,0))
  expect_equal(np[,10], c(0,0,0,1))
  expect_equal(np[,11], c(0,1,0,0))
  expect_equal(np[,12], c(1,0,0,0))
  expect_identical(all(np[,13:148] == 0), TRUE)
})

test_that("Test counts for motif 6", {
  np <- node_positions(M = motifs[[6]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,1], c(0,0,2,2))
  expect_equal(np[,2], c(2,2,0,0))
  expect_equal(np[,3], c(0,0,2,2))
  expect_equal(np[,4], c(1,1,0,0))
  expect_equal(np[,5], c(0,0,1,1))
  expect_equal(np[,6], c(2,2,0,0))
  expect_equal(np[,7], c(0,0,0,0))
  expect_equal(np[,8], c(0,0,0,0))
  expect_equal(np[,9], c(0,0,0,0))
  expect_equal(np[,10], c(0,0,0,0))
  expect_equal(np[,11], c(0,0,0,0))
  expect_equal(np[,12], c(0,0,0,0))
  expect_equal(np[,13], c(0,0,1,1))
  expect_equal(np[,14], c(1,1,0,0))
  expect_identical(all(np[,15:148] == 0), TRUE)
})

test_that("Test counts for motif 7", {
  np <- node_positions(M = motifs[[7]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,1], c(0,1,1,1))
  expect_equal(np[,2], c(3,0,0,0))
  expect_equal(np[,3], c(0,2,2,2))
  expect_equal(np[,4], c(3,0,0,0))
  expect_equal(np[,5], c(0,0,0,0))
  expect_equal(np[,6], c(0,0,0,0))
  expect_equal(np[,7], c(0,0,0,0))
  expect_equal(np[,8], c(0,0,0,0))
  expect_equal(np[,9], c(0,0,0,0))
  expect_equal(np[,10], c(0,0,0,0))
  expect_equal(np[,11], c(0,0,0,0))
  expect_equal(np[,12], c(0,0,0,0))
  expect_equal(np[,13], c(0,0,0,0))
  expect_equal(np[,14], c(0,0,0,0))
  expect_equal(np[,15], c(0,1,1,1))
  expect_equal(np[,16], c(1,0,0,0))
  expect_identical(all(np[,17:148] == 0), TRUE)
})

test_that("Test counts for motif 8", {
  np <- node_positions(M = motifs[[8]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,17], c(0,0,0,0,1))
  expect_equal(np[,18], c(1,1,1,1,0))
  expect_identical(all(np[,19:148] == 0), TRUE)
})

test_that("Test counts for motif 9", {
  np <- node_positions(M = motifs[[9]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,19], c(0,0,0,0,1))
  expect_equal(np[,20], c(0,0,0,1,0))
  expect_equal(np[,21], c(1,1,0,0,0))
  expect_equal(np[,22], c(0,0,1,0,0))
  expect_identical(all(np[,23:148] == 0), TRUE)
})

test_that("Test counts for motif 10", {
  np <- node_positions(M = motifs[[10]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,23], c(0,0,0,1,1))
  expect_equal(np[,24], c(1,0,1,0,0))
  expect_equal(np[,25], c(0,1,0,0,0))
  expect_identical(all(np[,26:148] == 0), TRUE)
})

test_that("Test counts for motif 11", {
  np <- node_positions(M = motifs[[11]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,26], c(0,0,0,0,1))
  expect_equal(np[,27], c(0,0,0,1,0))
  expect_equal(np[,28], c(1,0,0,0,0))
  expect_equal(np[,29], c(0,1,1,0,0))
  expect_identical(all(np[,30:148] == 0), TRUE)
})

test_that("Test counts for motif 12", {
  np <- node_positions(M = motifs[[12]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,30], c(0,0,0,1,1))
  expect_equal(np[,31], c(1,1,1,0,0))
  expect_identical(all(np[,32:148] == 0), TRUE)
})

test_that("Test counts for motif 13", {
  np <- node_positions(M = motifs[[13]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,32], c(0,0,0,1,1))
  expect_equal(np[,33], c(0,0,1,0,0))
  expect_equal(np[,34], c(1,0,0,0,0))
  expect_equal(np[,35], c(0,1,0,0,0))
  expect_identical(all(np[,36:148] == 0), TRUE)
})

test_that("Test counts for motif 14", {
  np <- node_positions(M = motifs[[14]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,36], c(0,0,1,0,1))
  expect_equal(np[,37], c(0,0,0,1,0))
  expect_equal(np[,38], c(1,1,0,0,0))
  expect_identical(all(np[,39:148] == 0), TRUE)
})

test_that("Test counts for motif 15", {
  np <- node_positions(M = motifs[[15]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,39], c(0,0,0,0,1))
  expect_equal(np[,40], c(0,0,1,1,0))
  expect_equal(np[,41], c(1,0,0,0,0))
  expect_equal(np[,42], c(0,1,0,0,0))
  expect_identical(all(np[,43:148] == 0), TRUE)
})

test_that("Test counts for motif 16", {
  np <- node_positions(M = motifs[[16]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,43], c(0,0,1,1,1))
  expect_equal(np[,44], c(1,1,0,0,0))
  expect_identical(all(np[,45:148] == 0), TRUE)
})

test_that("Test counts for motif 17", {
  np <- node_positions(M = motifs[[17]], six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
  expect_equal(np[,45], c(0,1,1,1,1))
  expect_equal(np[,46], c(1,0,0,0,0))
  expect_identical(all(np[,47:148] == 0), TRUE)
})

test_that("Error messages work",{
  expect_error(object = node_positions("a"), "'M' must be an object of class 'matrix'")
  expect_error(object = node_positions(matrix("a",3,3)), "Elements of 'M' must be numeric")
  expect_error(object = node_positions(matrix(-1,3,3)), "Elements of 'M' must be greater than or equal to zero")
  expect_error(object = node_positions(matrix(1,3,3), six_node = TRUE, level = 7, weights_method = "none", weights_combine = "none", normalisation = "none"), "'level' must be of class 'character'")
  expect_error(object = node_positions(matrix(1,3,3), six_node = TRUE, level = "wrong", weights_method = "none", weights_combine = "none", normalisation = "none"), "'level' must equal 'rows', 'columns' or 'all'")
  expect_error(object = node_positions(matrix(1,3,3), six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = 7), "'normalisation' must be of class 'character'")
  expect_error(object = node_positions(matrix(1,3,3), six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "wrong"), "'normalisation' must equal 'none','sum','sizeclass', 'sizeclass_plus1', 'sizeclass_NAzero', 'position','levelsize','levelsize_plus1', 'levelsize_NAzero', 'motif','motif_plus1','motif_NAzero'")
})

test_that("Matches validated results",{
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none"),positions_T_all_none)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_T_all_sum)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "sizeclass"), positions_T_all_sc)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_T_rows_none)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_T_rows_sum)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "sizeclass"), positions_T_rows_sc)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_T_columns_none)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_T_columns_sum)
  expect_equal(node_positions(M = mat, six_node = TRUE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "sizeclass"), positions_T_columns_sc)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_F_all_none)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_F_all_sum)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "sizeclass"), positions_F_all_sc)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_F_rows_none)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_F_rows_sum)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "rows", weights_method = "none", weights_combine = "none", normalisation = "sizeclass"), positions_F_rows_sc)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "none"), positions_F_columns_none)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "sum"), positions_F_columns_sum)
  expect_equal(node_positions(M = mat, six_node = FALSE, level = "columns", weights_method = "none", weights_combine = "none", normalisation = "sizeclass"), positions_F_columns_sc)
})

test_that("Correct output class for all argument combinations", {
  param_comb <- as.matrix(expand.grid(c('none', 'mean_motifweights', 'total_motifweights', 'mean_nodeweights', 'total_nodeweights', 'contribution', 'mora', 'all'),
                                      c('none', 'mean', 'sum'),
                                      c('none','sum','sizeclass', 'sizeclass_plus1', 'sizeclass_NAzero', 'position','levelsize','levelsize_plus1','levelsize_NAzero','motif','motif_plus1','motif_NAzero'),
                                      c(TRUE, FALSE)))

  M <- rbm(5,5)
  n <- nrow(param_comb)
  for (i in 1:n) {
    wm <- param_comb[i,1]
    wc <- param_comb[i,2]
    norm <- param_comb[i,3]
    sn <- ifelse(param_comb[i,4] == " TRUE" | param_comb[i,4] == "TRUE", yes = TRUE, no = FALSE)
    found_err <- FALSE
    tryCatch({
      np <- node_positions(M, level = "all", six_node = sn, weights_method = wm, weights_combine = wc, normalisation = norm)
      if (wm != 'all') {
        expect_equal(class(np), "data.frame")
      } else {
        expect_equal(class(np), "list")
        k <- length(np)
        for (j in 1:k) {
          expect_equal(class(np[[j]]), "data.frame")
        }
      }
    },
    error = function(e) {}
    )
  }
})

test_that("Default dimension names are not used for weighted node positions when the original matrix has dimension names",{
  param_comb <- as.matrix(expand.grid(c("rows","columns","all"),
                                      c('mean', 'sum'),
                                      c('mean_motifweights', 'total_motifweights', 'mean_nodeweights', 'total_nodeweights', 'contribution', 'mora', 'all'),
                                      c('none','sum','sizeclass', 'sizeclass_plus1', 'sizeclass_NAzero', 'position','levelsize','levelsize_plus1','levelsize_NAzero','motif','motif_plus1','motif_NAzero')))
  param_comb <- param_comb[!(param_comb[,2] == "mean" & param_comb[,4] != "none"),]
  M <- rwm(5,5)
  dimnames(M) <- list(LETTERS[1:5], LETTERS[6:10])
  n <- nrow(param_comb)
  for (i in 1:n) {
    l <- param_comb[i,1]
    wc <- param_comb[i,2]
    wm <- param_comb[i,3]
    norm <- param_comb[i,4]

    npw <- node_positions(M, level = l, six_node = FALSE, weights_method = wm, weights_combine = wc, normalisation = norm)

    if(wm != "all"){
      if(l == "rows"){
        expect_false(any(paste0("r", 1:nrow(M)) %in% rownames(npw)))
      } else if(l == "columns"){
        expect_false(any(paste0("c", 1:ncol(M)) %in% rownames(npw)))
      } else if(l == "all"){
        expect_false(any(paste0("r", 1:nrow(M)) %in% rownames(npw)))
        expect_false(any(paste0("c", 1:ncol(M)) %in% rownames(npw)))
      }
    } else if(wm == "all"){
      for(i in 1:length(npw)){
        if(l == "rows"){
          expect_false(any(paste0("r", 1:nrow(M)) %in% rownames(npw[[i]])))
        } else if(l == "columns"){
          expect_false(any(paste0("c", 1:ncol(M)) %in% rownames(npw[[i]])))
        } else if(l == "all"){
          expect_false(any(paste0("r", 1:nrow(M)) %in% rownames(npw[[i]])))
          expect_false(any(paste0("c", 1:ncol(M)) %in% rownames(npw[[i]])))
        }
      }
    }
  }
})

test_that("Dimension names are passed through correctly for weighted node positions when the original matrix has dimension names",{
  param_comb <- as.matrix(expand.grid(c("rows","columns","all"),
                                      c('mean', 'sum'),
                                      c('mean_motifweights', 'total_motifweights', 'mean_nodeweights', 'total_nodeweights', 'contribution', 'mora', 'all'),
                                      c('none','sum','sizeclass', 'sizeclass_plus1', 'sizeclass_NAzero', 'position','levelsize','levelsize_plus1','levelsize_NAzero','motif','motif_plus1','motif_NAzero')))
  param_comb <- param_comb[!(param_comb[,2] == "mean" & param_comb[,4] != "none"),]
  M <- rwm(5,5)
  n <- nrow(param_comb)
  for(j in 1:4){
    if(j == 1){
      dimnames(M) <- list(LETTERS[1:5], LETTERS[6:10])
    } else if(j == 2){
      dimnames(M) <- list(LETTERS[1:5], LETTERS[6:10])
      rownames(M) <- NULL
    } else if(j == 3){
      dimnames(M) <- list(LETTERS[1:5], LETTERS[6:10])
      colnames(M) <- NULL
    } else if(j == 4){
      dimnames(M) <- list(LETTERS[1:5], LETTERS[6:10])
      dimnames(M) <- NULL
    }
    for (i in 1:n) {
      l <- param_comb[i,1]
      wc <- param_comb[i,2]
      wm <- param_comb[i,3]
      norm <- param_comb[i,4]

      npw <- node_positions(M, level = l, six_node = FALSE, weights_method = wm, weights_combine = wc, normalisation = norm)
      np <- node_positions(M, level = l, six_node = FALSE, weights_method = "none", weights_combine = "none", normalisation = "none")

      if(wm != "all"){
        expect_identical(rownames(npw), rownames(np))
        if(l == "rows"){
          if(j == 1){
            expect_identical(rownames(npw), rownames(M))
          } else if(j == 2){
            expect_identical(rownames(npw), paste0("r",1:nrow(M)))
          } else if(j == 3){
            expect_identical(rownames(npw), rownames(M))
          } else if(j == 4){
            expect_identical(rownames(npw), paste0("r",1:nrow(M)))
          }
        } else if(l == "columns"){
          if(j == 1){
            expect_identical(rownames(npw), colnames(M))
          } else if(j == 2){
            expect_identical(rownames(npw), colnames(M))
          } else if(j == 3){
            expect_identical(rownames(npw), paste0("c",1:nrow(M)))
          } else if(j == 4){
            expect_identical(rownames(npw), paste0("c",1:nrow(M)))
          }
        } else if(l == "all"){
          if(j == 1){
            expect_identical(rownames(npw), c(rownames(M), colnames(M)))
          } else if(j == 2){
            expect_identical(rownames(npw), c(paste0("r",1:nrow(M)), colnames(M)))
          } else if(j == 3){
            expect_identical(rownames(npw), c(rownames(M), paste0("c",1:nrow(M))))
          } else if(j == 4){
            expect_identical(rownames(npw), c(paste0("r",1:nrow(M)), paste0("c",1:nrow(M))))
          }
        }
      } else if(wm == "all"){
        for(k in 1:length(npw)){
          expect_identical(rownames(npw[[k]]), rownames(np))
          if(l == "rows"){
            if(j == 1){
              expect_identical(rownames(npw[[k]]), rownames(M))
            } else if(j == 2){
              expect_identical(rownames(npw[[k]]), paste0("r",1:nrow(M)))
            } else if(j == 3){
              expect_identical(rownames(npw[[k]]), rownames(M))
            } else if(j == 4){
              expect_identical(rownames(npw[[k]]), paste0("r",1:nrow(M)))
            }
          } else if(l == "columns"){
            if(j == 1){
              expect_identical(rownames(npw[[k]]), colnames(M))
            } else if(j == 2){
              expect_identical(rownames(npw[[k]]), colnames(M))
            } else if(j == 3){
              expect_identical(rownames(npw[[k]]), paste0("c",1:nrow(M)))
            } else if(j == 4){
              expect_identical(rownames(npw[[k]]), paste0("c",1:nrow(M)))
            }
          } else if(l == "all"){
            if(j == 1){
              expect_identical(rownames(npw[[k]]), c(rownames(M), colnames(M)))
            } else if(j == 2){
              expect_identical(rownames(npw[[k]]), c(paste0("r",1:nrow(M)), colnames(M)))
            } else if(j == 3){
              expect_identical(rownames(npw[[k]]), c(rownames(M), paste0("c",1:nrow(M))))
            } else if(j == 4){
              expect_identical(rownames(npw[[k]]), c(paste0("r",1:nrow(M)), paste0("c",1:nrow(M))))
            }
          }
        }
      }
    }
  }
})
