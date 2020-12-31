context("weighted node positions")

# test for binary matrices
# motif weights are always 1

# save number of edges and vertices for each motif
motif_n_edges <- c(1,2,2,3,3,4,3,4,4,4,5,6,4,4,5,6,4)
motif_n_vertices <- c(2,3,3,4,4,4,4,rep(5, 10))
# save degree of each position
pos_degrees <- c(c(1,1), c(1,2,2,1), c(3,1,1,2,1,2,2,2,1,3), c(4,1,1,3,1,2,2,1,2,2,3,1,2), c(3,2,1,2,1,3,1,2,2,1,2,2,3,2,3,1,4))


test_that("Test motif weights for binary, weights_combine = mean", {
  for (i in 1:10) {
    W <- rbm(5,5)
    np <- node_positions(W, six_node = FALSE, level = "all", weights_method = "all", weights_combine = "mean", normalisation = "none")
   # print(np)

    npc <- as.matrix(node_positions(W, six_node = FALSE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none"))
    mmw <- as.matrix(np[[1]])
    npc_zero <- which(npc == 0)
    expect_true(all(is.na(mmw[npc_zero]))) # when position count is zero, mean motif weights are NA
    expect_true(all(mmw[-npc_zero] == 1)) # all other motif weights are one
    mmw <- np[[1]]

    mnw <- as.matrix(np[[3]])
    expect_true(all(is.na(mnw[npc_zero]))) # when position count is zero,mean node weights are NA
    expect_true(all(mnw[-npc_zero] == 1)) # all other node weights are one
    mnw <- np[[3]]

    # --- more stuff needed for testing tmw, tnw, con and py ---

    # create a matrix with number of edges for each motif
    nedges <- matrix(NA, nrow(W) + ncol(W), 46)
    for (i in 1:17) {
      positions <- motif_info(i, node = TRUE, link = FALSE)
      nedges[,positions] <- motif_n_edges[i]
    }

    degrees <- matrix(NA, nrow(W) + ncol(W), 46)
    for (i in 1:46) {
      degrees[,i] <- pos_degrees[i]
    }

    nvertices <- matrix(NA, nrow(W) + ncol(W), 46)
    for (i in 1:17) {
      positions <- motif_info(i, node = TRUE, link = FALSE)
      nvertices[,positions] <- motif_n_vertices[i]
    }

    # -------- actual testing ------------------



    tmw <- np[[2]]
    expect_equal(tmw / nedges, mmw)

    tnw <- np[[4]]
    expect_equal(tnw / degrees, mnw)

    # con is degree / number of edges in the motif
    con <- as.matrix(np[[5]])
    con_res <- degrees / nedges
    con_not_na <- which(!is.na(con))
    expect_equal(con[con_not_na], con_res[con_not_na])

    con_na <- which(is.na(con))
    expect_true(all(npc[con_na] == 0)) # check if we really only get NAs if the count is zero
    expect_true(all(is.na(con[npc_zero]))) # reverse: if the count is zero, we expect NAs
    con <- np[[5]]

    # py is 1 / number of edges in the motif
    py <- as.matrix(np[[6]])
    py_res <- 1 / nvertices

    py_not_na <- which(!is.na(py))
    expect_equal(py[py_not_na], py_res[py_not_na])

    py_na <- which(is.na(py))
    expect_true(all(npc[py_na] == 0))
    expect_true(all(is.na(py[npc_zero])))
    py <- np[[6]]

    # NOW TEST SINGLE METHODS
    np_mmw <- node_positions(W, six_node = FALSE, level = "all", weights_method = "mean_motifweights", weights_combine = "mean", normalisation = "none")
    expect_equal(np_mmw, mmw)

    np_tmw <- node_positions(W, six_node = FALSE, level = "all", weights_method = "total_motifweights", weights_combine = "mean", normalisation = "none")
    expect_equal(np_tmw, tmw)

    np_mnw <- node_positions(W, six_node = FALSE, level = "all", weights_method = "mean_nodeweights", weights_combine = "mean", normalisation = "none")
    expect_equal(np_mnw, mnw)

    np_tnw <- node_positions(W, six_node = FALSE, level = "all", weights_method = "total_nodeweights", weights_combine = "mean", normalisation = "none")
    expect_equal(np_tnw, tnw)

    np_con <- node_positions(W, six_node = FALSE, level = "all", weights_method = "contribution", weights_combine = "mean", normalisation = "none")
    expect_equal(np_con, con)

    np_py <- node_positions(W, six_node = FALSE, level = "all", weights_method = "mora", weights_combine = "mean", normalisation = "none")
    expect_equal(np_py, py)
}
})

test_that("Test motif weights for binary, weights_combine = sum", {
  for (i in 1:10) {
    W <- rbm(5,5)
    np <- node_positions(W, six_node = FALSE, level = "all", weights_method = "all", weights_combine = "sum", normalisation = "none")
    # print(np)

    npc <- node_positions(W, six_node = FALSE, level = "all", weights_method = "none", weights_combine = "none", normalisation = "none")
    mmw <- np[[1]]
    expect_equal(mmw, npc) # mean motif weight is always one, so sum will equal the count

    mnw <- np[[3]]
    expect_equal(mnw, npc) # mean node weight is always one, so sum will equal the count

    # --- more stuff needed for testing tmw, tnw, con and py ---

    # create a matrix with number of edges for each motif
    nedges <- matrix(NA, nrow(W) + ncol(W), 46)
    for (i in 1:17) {
      positions <- motif_info(i, node = TRUE, link = FALSE)
      nedges[,positions] <- motif_n_edges[i]
    }

    degrees <- matrix(NA, nrow(W) + ncol(W), 46)
    for (i in 1:46) {
      degrees[,i] <- pos_degrees[i]
    }

    nvertices <- matrix(NA, nrow(W) + ncol(W), 46)
    for (i in 1:17) {
      positions <- motif_info(i, node = TRUE, link = FALSE)
      nvertices[,positions] <- motif_n_vertices[i]
    }

    # -------- actual testing ------------------

    tmw <- np[[2]]
    expect_equal(tmw / nedges, mmw)

    tnw <- np[[4]]
    expect_equal(tnw / degrees, mnw)

    # con is count * degree / number of edges in the motif
    con <- np[[5]]
    con_res <- npc * degrees / nedges
    expect_equal(con, con_res)

    # py is count / number of edges in the motif
    py <- np[[6]]
    py_res <- npc / nvertices

    expect_equal(py, py_res)

    # NOW TEST SINGLE METHODS
    np_mmw <- node_positions(W, six_node = FALSE, level = "all", weights = "mean_motifweights", weights_combine = "sum", normalisation = "none")
    expect_equal(np_mmw, mmw)

    np_tmw <- node_positions(W, six_node = FALSE, level = "all", weights = "total_motifweights", weights_combine = "sum", normalisation = "none")
    expect_equal(np_tmw, tmw)

    np_mnw <- node_positions(W, six_node = FALSE, level = "all", weights = "mean_nodeweights", weights_combine = "sum", normalisation = "none")
    expect_equal(np_mnw, mnw)

    np_tnw <- node_positions(W, six_node = FALSE, level = "all", weights = "total_nodeweights", weights_combine = "sum", normalisation = "none")
    expect_equal(np_tnw, tnw)

    np_con <- node_positions(W, six_node = FALSE, level = "all", weights = "contribution", weights_combine = "sum", normalisation = "none")
    expect_equal(np_con, con)

    np_py <- node_positions(W, six_node = FALSE, level = "all", weights = "mora", weights_combine = "sum", normalisation = "none")
    expect_equal(np_py, py)
  }
})



# test_that("We test one weighted version of motif 9",{
#   W <- matrix(c(1,2,3,0,0,4), byrow = FALSE, nrow = 3)
#   np <- node_positions(W, six_node = FALSE, level = "all", weights_method = "all", weights_combine = "mean", normalisation = "none")
#
#   # If I have time this week, will add this example also
# })
