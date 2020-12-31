context("link_positions")

# -------------- TEST SPECIAL CASES -------------------------------------

test_that("Test for matrix with a single link", {
  lp <- link_positions(matrix(1,1,1), six_node = FALSE, weights = FALSE, normalisation = "none")
  expect_true(class(lp) == "data.frame")
  expect_true(nrow(lp) == 1)
  expect_equivalent(lp[1,1], 1)
  expect_equivalent(as.vector(lp[,2:29]), rep(0,28))
})

test_that("Test for matrix with a single link, six_node is true", {
  lp <- link_positions(matrix(1,1,1), six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_true(class(lp) == "data.frame")
  expect_true(nrow(lp) == 1)
  expect_equivalent(lp[1,1], 1)
  expect_equivalent(as.vector(lp[,2:106]), rep(0,105))
})

test_that("Test for matrix without a link", {
  expect_error(link_positions(matrix(0,1,1), six_node = FALSE, weights = FALSE, normalisation = "none"))
})

# --------------- TEST ON MOTIFS ----------------------------------------

test_that("Test counts for the motif 1", {
    lr01 <- link_positions(motifs[[1]], six_node = FALSE, weights = FALSE, normalisation = "none")
    expect_equal( as.vector(lr01[,1]), c(1))
    expect_equal( as.vector(lr01[,2]), c(0))
})

test_that("Test counts for the motif 2", {
  lr02 <- link_positions( motifs[[2]], six_node = FALSE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr02[,1]), c(1,1))
  expect_equal(as.vector(lr02[,2]), c(1,1))
  expect_equal(as.vector(lr02[,3]), c(0,0))
})

test_that("Test counts for the motif 3", {
  lr <- link_positions( motifs[[3]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,1]), c(1,1))
  expect_equal(as.vector(lr[,2]), c(0,0))
  expect_equal(as.vector(lr[,3]), c(1,1))
})

test_that("Test counts for the motif 4", {
  lr <- link_positions( motifs[[4]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,1]), c(1,1,1))
  expect_equal(as.vector(lr[,2]), c(0,0,0))
  expect_equal(as.vector(lr[,3]), c(2,2,2))
  expect_equal(as.vector(lr[,4]), c(1,1,1))
})

test_that("Test counts for the motif 5", {
  lr <- link_positions( motifs[[5]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,1]), c(1,1,1))
  expect_equal(as.vector(lr[,2]), c(1,1,0))
  expect_equal(as.vector(lr[,3]), c(0,1,1))
  expect_equal(as.vector(lr[,4]), c(0,0,0))
  expect_equal(as.vector(lr[,5]), c(1,0,0))
  expect_equal(as.vector(lr[,6]), c(0,0,1))
  expect_equal(as.vector(lr[,7]), c(0,1,0))
})

test_that("Test counts for the motif 6", {
  lr <- link_positions( motifs[[6]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,1]), c(1,1,1,1))
  expect_equal(as.vector(lr[,2]), rep(1,4))
  expect_equal(as.vector(lr[,3]), rep(1,4))
  expect_equal(as.vector(lr[,4]), rep(0,4))
  expect_equal(as.vector(lr[,5]), rep(0,4))
  expect_equal(as.vector(lr[,6]), rep(0,4))
  expect_equal(as.vector(lr[,7]), rep(0,4))
  expect_equal(as.vector(lr[,8]), rep(1,4))
})

test_that("Test counts for the motif 7", {
  lr <- link_positions( motifs[[7]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,1]), c(1,1,1))
  expect_equal(as.vector(lr[,2]), rep(2,3))
  expect_equal(as.vector(lr[,3]), rep(0,3))
  expect_equal(as.vector(lr[,4]), rep(0,3))
  expect_equal(as.vector(lr[,5]), rep(0,3))
  expect_equal(as.vector(lr[,6]), rep(0,3))
  expect_equal(as.vector(lr[,7]), rep(0,3))
  expect_equal(as.vector(lr[,8]), rep(0,3))
  expect_equal(as.vector(lr[,9]), rep(1,3))
})

# now only test selected columns

test_that("Test counts for the motif 8", {
  lr <- link_positions( motifs[[8]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,3]), rep(3,4))
  expect_equal(as.vector(lr[,4]), rep(3,4))
})

test_that("Test counts for the motif 9", {
  lr <- link_positions( motifs[[9]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,5]), c(0,0,0,2))
  expect_equal(as.vector(lr[,6]), c(1,1,0,0))
  expect_equal(as.vector(lr[,7]), c(0,0,2,0))
})

# skip to complicated motifs
test_that("Test counts for the motif 19", {
  lr <- link_positions( motifs[[19]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,5]), c(3,0,0,0,0))
  expect_equal(as.vector(lr[,6]), c(0,0,1,1,1))
  expect_equal(as.vector(lr[,7]), c(0,3,0,0,0))
  expect_equal(as.vector(lr[,11]), c(3,0,0,0,0))
  expect_equal(as.vector(lr[,12]), c(0,0,2,2,2))
  expect_equal(as.vector(lr[,13]), c(0,3,0,0,0))
})

test_that("Test counts for the motif 20", {
  lr <- link_positions( motifs[[20]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,14]), c(2,0,0,1,1))
  expect_equal(as.vector(lr[,15]), c(0,2,2,0,0))
})

test_that("Test counts for the motif 21", {
  lr <- link_positions( motifs[[21]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,10]), c(1,1,1,0,1,0))
  expect_equal(as.vector(lr[,11]), c(0,0,0,1,0,1))
  expect_equal(as.vector(lr[,12]), c(2,2,0,0,0,0))
  expect_equal(as.vector(lr[,13]), c(0,0,1,0,1,0))
  expect_equal(as.vector(lr[,14]), rep(0,6))
  expect_equal(as.vector(lr[,16]), c(0,0,0,2,0,2))
  expect_equal(as.vector(lr[,17]), c(1,1,0,0,0,0))
  expect_equal(as.vector(lr[,18]), c(0,0,2,0,2,0))
})

test_that("Test counts for the motif 22", {
  lr <- link_positions( motifs[[22]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,10]), rep(0,6))
  expect_equal(as.vector(lr[,11]), rep(0,6))
  expect_equal(as.vector(lr[,12]), rep(0,6))
  expect_equal(as.vector(lr[,13]), rep(0,6))
  expect_equal(as.vector(lr[,14]), c(2,0,0,0,0,2))
  expect_equal(as.vector(lr[,15]), c(0,1,1,1,1,0))
  expect_equal(as.vector(lr[,16]), c(0,1,1,1,1,0))
  expect_equal(as.vector(lr[,17]), c(1,0,0,0,0,1))
  expect_equal(as.vector(lr[,18]), c(0,1,1,1,1,0))
})

test_that("Test counts for the motif 23", {
  lr <- link_positions( motifs[[23]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,10]), c(0,1,0,1,0,1,1))
  expect_equal(as.vector(lr[,12]), rep(0,7))
  expect_equal(as.vector(lr[,14]), rep(0,7))
  expect_equal(as.vector(lr[,16]), c(2,0,2,0,2,0,0))
  expect_equal(as.vector(lr[,17]), c(rep(0,6), 3))
  expect_equal(as.vector(lr[,18]), c(0,2,0,2,0,2,0))
  expect_equal(as.vector(lr[,19]), c(rep(1,6),0))
})


test_that("Test counts for the motif 24", {
  lr <- link_positions( motifs[[24]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,19]), rep(3,8))
})

# tests edges 20-22
test_that("Test counts for the motif 25", {
  lr <- link_positions( motifs[[25]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,20]), c(0,0,0,2,2))
  expect_equal(as.vector(lr[,21]), c(1,1,0,0,0))
  expect_equal(as.vector(lr[,22]), c(0,0,2,0,0))
})
test_that("Test counts for the motif 26", {
  lr <- link_positions( motifs[[26]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,20]), c(0,1,2,1,0))
  expect_equal(as.vector(lr[,21]), c(1,0,0,0,1))
  expect_equal(as.vector(lr[,22]), c(0,1,0,1,0))
})

test_that("Test counts for the motif 29", {
  lr <- link_positions( motifs[[29]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,20]), c(rep(0,4), 1, 1))
  expect_equal(as.vector(lr[,21]), c(1,rep(0,5)))
  expect_equal(as.vector(lr[,22]), c(0,0,0,1,0,0))
})

# tests edges 23-24
test_that("Test counts for the motif 30", {
  lr <- link_positions( motifs[[30]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,23]), c(2,0,0,1,0,1))
  expect_equal(as.vector(lr[,24]), c(0,2,1,0,1,0))
})

test_that("Test counts for the motif 28", {
  lr <- link_positions( motifs[[28]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,23]), c(1,0,0,1,0))
  expect_equal(as.vector(lr[,24]), c(0,1,1,0,0))
})

# tests edges 25-27
test_that("Test counts for the motif 40", {
  lr <- link_positions( motifs[[40]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,25]), c(rep(0,4),1,1))
  expect_equal(as.vector(lr[,26]), c(2,2,rep(0,4)))
  expect_equal(as.vector(lr[,27]), c(0,0,2,2,0,0))
})

test_that("Test counts for the motif 41", {
  lr <- link_positions( motifs[[41]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,25]), c(1, rep(0,4), 1))
  expect_equal(as.vector(lr[,26]), c(0,1,1,1,1,0))
  expect_equal(as.vector(lr[,27]), c(0,1,1,1,1,0))
})


# tests edge 28
test_that("Test counts for the motif 42", {
  lr <- link_positions( motifs[[42]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,28]), c(rep(1,6),0))
})

test_that("Test counts for the motif 43", {
  lr <- link_positions( motifs[[43]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,28]), c(rep(3,8)))
})

# test edge 29
test_that("Test counts for the motif 44", {
  lr <- link_positions( motifs[[44]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,29]), c(rep(4,5)))
})

test_that("Test counts for the motif 43", {
  lr <- link_positions( motifs[[43]], six_node = FALSE, weights = FALSE, normalisation = "none" )
  expect_equal(as.vector(lr[,29]), c(rep(1,8)))
})

#####################################
# 6 vertices

# motif 18
test_that("Test counts for the motif 18", {
  lr <- link_positions( motifs[[18]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,30]), rep(1,5))
})

# motif 19
test_that("Test counts for the motif 19", {
  lr <- link_positions( motifs[[19]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,31]), c(1, rep(0,4)))
  expect_equal(as.vector(lr[,32]), c(0,0,1,1,1))
  expect_equal(as.vector(lr[,33]), c(0,1,0,0,0))
})

# motif 20
test_that("Test counts for the motif 20", {
  lr <- link_positions( motifs[[20]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,34]), c(1, rep(0,4)))
  expect_equal(as.vector(lr[,35]), c(0,1,0,0,0))
  expect_equal(as.vector(lr[,36]), c(0,0,0,1,1))
  expect_equal(as.vector(lr[,37]), c(0,0,1,0,0))
})

# motif 21
test_that("Test counts for the motif 21", {
  lr <- link_positions( motifs[[21]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,38]), c(0,0,0,1,0,1))
  expect_equal(as.vector(lr[,39]), c(1,1,rep(0,4)))
  expect_equal(as.vector(lr[,40]), c(0,0,1,0,1,0))
})

# motif 22
test_that("Test counts for the motif 22", {
  lr <- link_positions( motifs[[22]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,41]), c(1,0,0,0,0,1))
  expect_equal(as.vector(lr[,42]), c(0,1,1,1,1,0))
})

# motif 23
test_that("Test counts for the motif 23", {
  lr <- link_positions( motifs[[23]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,43]), c(1,0,1,0,1,0,0))
  expect_equal(as.vector(lr[,44]), c(rep(0,6),1))
  expect_equal(as.vector(lr[,45]), c(0,1,0,1,0,1,0))
})

# motif 24
test_that("Test counts for the motif 24", {
  lr <- link_positions( motifs[[24]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,46]), c(rep(1,8)))
})


# motif 25
test_that("Test counts for the motif 25", {
  lr <- link_positions( motifs[[25]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,47]), c(0,0,0,1,1))
  expect_equal(as.vector(lr[,48]), c(1,1,0,0,0))
  expect_equal(as.vector(lr[,49]), c(0,0,1,0,0))
})

test_that("Test counts for the motif 26", {
  lr <- link_positions( motifs[[26]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,50]), c(0,1,0,1,0))
  expect_equal(as.vector(lr[,51]), c(1,0,0,0,1))
  expect_equal(as.vector(lr[,52]), c(0,0,1,0,0))
})

test_that("Test counts for the motif 27", {
  lr <- link_positions( motifs[[27]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,53]), c(0,1,0,1,0))
  expect_equal(as.vector(lr[,54]), c(0,0,1,0,0))
  expect_equal(as.vector(lr[,55]), c(1,0,0,0,1))
})

test_that("Test counts for the motif 28", {
  lr <- link_positions( motifs[[28]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,56]), c(1,rep(0,4)))
  expect_equal(as.vector(lr[,57]), c(rep(0,4),1))
  expect_equal(as.vector(lr[,58]), c(0,0,0,1,0))
  expect_equal(as.vector(lr[,59]), c(0,1,0,0,0))
  expect_equal(as.vector(lr[,60]), c(0,0,1,0,0))
})

test_that("Test counts for the motif 29", {
  lr <- link_positions( motifs[[29]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,61]), c(rep(0,5),1))
  expect_equal(as.vector(lr[,62]), c(0,0,1,0,0,0))
  expect_equal(as.vector(lr[,63]), c(0,0,0,0,1,0))
  expect_equal(as.vector(lr[,64]), c(0,1,rep(0,4)))
  expect_equal(as.vector(lr[,65]), c(0,0,0,1,0,0))
  expect_equal(as.vector(lr[,66]), c(1,rep(0,5)))
})

test_that("Test counts for the motif 30", {
  lr <- link_positions( motifs[[30]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,67]), c(1,rep(0,5)))
  expect_equal(as.vector(lr[,68]), c(0,0,0,1,0,1))
  expect_equal(as.vector(lr[,69]), c(0,1,0,0,0,0))
  expect_equal(as.vector(lr[,70]), c(0,0,1,0,1,0))
})

test_that("Test counts for the motif 31", {
  lr <- link_positions( motifs[[31]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,71]), c(1,rep(0,6)))
  expect_equal(as.vector(lr[,72]), c(0,0,0,1,1,1,1))
  expect_equal(as.vector(lr[,73]), c(0,1,1,0,0,0,0))
})

test_that("Test counts for the motif 32", {
  lr <- link_positions( motifs[[32]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,74]), c(rep(0,5),1))
  expect_equal(as.vector(lr[,75]), c(0,0,0,0,1,0))
  expect_equal(as.vector(lr[,76]), c(1,1,0,0,0,0))
  expect_equal(as.vector(lr[,77]), c(0,0,1,1,0,0))
})

test_that("Test counts for the motif 33", {
  lr <- link_positions( motifs[[33]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,78]), c(0,0,1,1,0,1,1))
  expect_equal(as.vector(lr[,79]), c(1,rep(0,6)))
  expect_equal(as.vector(lr[,80]), c(0,1,0,0,1,0,0))
})

test_that("Test counts for the motif 34", {
  lr <- link_positions( motifs[[34]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,81]), rep(1,6))
})

test_that("Test counts for the motif 35", {
  lr <- link_positions( motifs[[35]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,82]), c(1,0,0,0,0,0,1))
  expect_equal(as.vector(lr[,83]), c(0,0,1,0,1,0,0))
  expect_equal(as.vector(lr[,84]), c(0,1,0,0,0,1,0))
  expect_equal(as.vector(lr[,85]), c(0,0,0,1,0,0,0))
})

test_that("Test counts for the motif 36", {
  lr <- link_positions( motifs[[36]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,86]), c(1,0,0,1,0,0,0,0))
  expect_equal(as.vector(lr[,87]), c(rep(0,6), 1,1))
  expect_equal(as.vector(lr[,88]), c(0,1,1,0,1,1,0,0))
})

test_that("Test counts for the motif 37", {
  lr <- link_positions( motifs[[37]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,89]), rep(1,9))
})

# motif 38
test_that("Test counts for the motif 38", {
  lr <- link_positions( motifs[[38]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,90]), c(0,0,1,1,1))
  expect_equal(as.vector(lr[,91]), c(1,rep(0,4)))
  expect_equal(as.vector(lr[,92]), c(0,1,0,0,0))
})

# motif 39
test_that("Test counts for the motif 39", {
  lr <- link_positions( motifs[[39]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,93]), c(1,rep(0,4)))
  expect_equal(as.vector(lr[,94]), c(0,0,0,1,1))
  expect_equal(as.vector(lr[,95]), c(0,1,0,0,0))
  expect_equal(as.vector(lr[,96]), c(0,0,1,0,0))
})

# motif 40
test_that("Test counts for the motif 40", {
  lr <- link_positions( motifs[[40]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,97]), c(rep(0,4),1,1))
  expect_equal(as.vector(lr[,98]), c(1,1,rep(0,4)))
  expect_equal(as.vector(lr[,99]), c(0,0,1,1,0,0))
})

# motif 41
test_that("Test counts for the motif 41", {
  lr <- link_positions( motifs[[41]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,100]), c(1,rep(0,4),1))
  expect_equal(as.vector(lr[,101]), c(0,rep(1,4),0))
})

# motif 42
test_that("Test counts for the motif 42", {
  lr <- link_positions( motifs[[42]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,102]), c(rep(0,6),1))
  expect_equal(as.vector(lr[,103]), c(rep(1,3), rep(0,4)))
  expect_equal(as.vector(lr[,104]), c(rep(0,3), rep(1,3), 0))
})

# motif 43
test_that("Test counts for the motif 43", {
  lr <- link_positions( motifs[[43]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,105]), c(rep(1,8)))
})

# motif 44
test_that("Test counts for the motif 44", {
  lr <- link_positions( motifs[[44]], six_node = TRUE, weights = FALSE, normalisation = "none")
  expect_equal(as.vector(lr[,106]), c(rep(1,5)))
})


# ---------- Do test counts for big block matrices --------------------

mlist4 <- motifs[4:7]
bm4 <- block_matrix(mlist4)

test_that("Test counts for the motifs with 4 nodes", {
  lr <- link_positions(bm4, six_node = TRUE, weights = FALSE, normalisation = "none")
  rownames(lr) <- c()
  colnames(lr) <- c()
  expect_equal(rowSums(lr[,4:9]), rep(1,nrow(lr)))
  expect_equal(colSums(lr[,4:9]), c(3,1,1,1,4,3))
  expect_equal(rowSums(lr[,10:106]), rep(0,nrow(lr)))
})

mlist5 <- motifs[8:17]
bm5 <- block_matrix(mlist5)

test_that("Test counts for the motifs with 5 nodes", {
  lr <- link_positions(bm5, six_node = TRUE, weights = FALSE, normalisation = "none")
  rownames(lr) <- c()
  colnames(lr) <- c()
  expect_equal(rowSums(lr[,10:29]), rep(1,nrow(lr)))
  v <- c(4, 1, 2, 1, 2, 2, 2, 1, 2, 6, 2, 1, 1, 2, 2, 1, 2, 2, 6, 4)
  expect_equal(colSums(lr[,10:29]), v)
  expect_equal(rowSums(lr[,30:106]), rep(0,nrow(lr)))
})


mlist6 <- motifs[18:44]
bm6 <- block_matrix(mlist6)

test_that("Test counts for the motifs with 6 nodes", {
  lr <- link_positions(bm6, six_node = TRUE, weights = FALSE, normalisation = "none")
  rownames(lr) <- c()
  colnames(lr) <- c()
  expect_equal(rowSums(lr[,30:106]), rep(1,nrow(lr)))
  expect_equal(sum(lr[,30]), sum(motifs[[18]]))
  expect_equal(colSums(lr[,31:33]), c(1,3,1))
  expect_equal(colSums(lr[,34:37]), c(1,1,2,1))
  expect_equal(colSums(lr[,38:40]), c(2,2,2))
  expect_equal(colSums(lr[,41:50]), c(2, 4, 3, 1, 3, 8, 2, 2, 1, 2))
  expect_equal(colSums(lr[,51:60]), c(2, 1, 2, 1, 2, 1, 1, 1, 1, 1))
  expect_equal(colSums(lr[,61:70]), c(1, 1, 1, 1, 1, 1, 1, 2, 1, 2))
  expect_equal(colSums(lr[,71:80]), c(1, 4, 2, 1, 1, 2, 2, 4, 1, 2))
  expect_equal(colSums(lr[,81:90]), c(6, 2, 2, 2, 1, 2, 2, 4, 9, 3))
  expect_equal(colSums(lr[,91:100]), c(1, 1, 1, 2, 1, 1, 2, 2, 2, 2))
  expect_equal(colSums(lr[,101:106]), c(4, 1, 3, 3, 8, 5))
})


# ------------------ TEST NORMALISATION ---------------------------

test_that("Check correct output formats", {
  lp <- link_positions(matrix(1,5,5), six_node = FALSE, weights = FALSE, normalisation = "none")
  lp_sum <- link_positions(matrix(1,5,5), six_node = FALSE, weights = FALSE, normalisation = "sum")
  lp_sizeclass <- link_positions(matrix(1,5,5), six_node = FALSE, weights = FALSE, normalisation = "sizeclass")
  lp_pos <- link_positions(matrix(1,5,5), six_node = FALSE, weights = FALSE, normalisation = "position")
  lp_level <- link_positions(matrix(1,5,5), six_node = FALSE, weights = FALSE, normalisation = "levelsize")
  lp_motif <- link_positions(matrix(1,5,5), six_node = FALSE, weights = FALSE, normalisation = "motif")
  expect_equal(rownames(lp_sum), rownames(lp))
  expect_equal(rownames(lp_sizeclass), rownames(lp))
  expect_equal(rownames(lp_pos), rownames(lp))
  expect_equal(rownames(lp_level), rownames(lp))
  expect_equal(rownames(lp_motif), rownames(lp))
  expect_equal(ncol(lp_sum), 29)
  expect_equal(ncol(lp_sizeclass), 29)
  expect_equal(ncol(lp_pos), 29)
  expect_equal(ncol(lp_level), 29)
  expect_equal(ncol(lp_motif), 29)
})

test_that("Check correct output formats", {
  lp <- link_positions(matrix(1,5,5), six_node = TRUE, weights = FALSE, normalisation = "none")
  lp_sum <- link_positions(matrix(1,5,5), six_node = TRUE, weights = FALSE, normalisation = "sum")
  lp_sizeclass <- link_positions(matrix(1,5,5), six_node = TRUE, weights = FALSE, normalisation = "sizeclass")
  lp_pos <- link_positions(matrix(1,5,5), six_node = TRUE, weights = FALSE, normalisation = "position")
  lp_level <- link_positions(matrix(1,5,5), six_node = TRUE, weights = FALSE, normalisation = "levelsize")
  lp_motif <- link_positions(matrix(1,5,5), six_node = TRUE, weights = FALSE, normalisation = "motif")
  expect_equal(rownames(lp_sum), rownames(lp))
  expect_equal(rownames(lp_sizeclass), rownames(lp))
  expect_equal(rownames(lp_pos), rownames(lp))
  expect_equal(rownames(lp_level), rownames(lp))
  expect_equal(rownames(lp_motif), rownames(lp))
  expect_equal(ncol(lp_sum), 106)
  expect_equal(ncol(lp_sizeclass), 106)
  expect_equal(ncol(lp_pos), 106)
  expect_equal(ncol(lp_level), 106)
  expect_equal(ncol(lp_motif), 106)
})


test_that("Test sum-normalisation for single link matrix", {
  lp <- link_positions(matrix(1,1,1), six_node = FALSE, weights = FALSE, normalisation = "sum")
  expect_equivalent(lp[1,1], 1)
  expect_equal(nrow(lp), 1)
  expect_equivalent(as.vector(lp[,2:29]), rep(0,28))
  expect_equal(rownames(lp), "r1 -- c1")
})

test_that("Test sum-normalisation for single link matrix, six_node = TRUE", {
  lp <- link_positions(matrix(1,1,1), six_node = TRUE, weights = FALSE, normalisation = "sum")
  expect_equivalent(lp[1,1], 1)
  expect_equal(nrow(lp), 1)
  expect_equal(ncol(lp), 106)
  expect_equivalent(as.vector(lp[,2:106]), rep(0,105))
  expect_equal(rownames(lp), "r1 -- c1")
})

test_that("Test sizeclass-normalisation for single link matrix", {
  lp <- link_positions(matrix(1,1,1), six_node = FALSE, weights = FALSE, normalisation = "sizeclass")
  expect_equivalent(lp[1,1], 1)
  expect_equal(nrow(lp), 1)
  expect_equal(rownames(lp), "r1 -- c1")
  expect_true(all(is.na(lp[,2:29])))
})

test_that("Test sizeclass-normalisation for single link matrix, six_node = TRUE", {
  lp <- link_positions(matrix(1,1,1), six_node = TRUE, weights = FALSE, normalisation = "sizeclass")
  expect_equivalent(lp[1,1], 1)
  expect_equal(nrow(lp), 1)
  expect_equal(ncol(lp), 106)
  expect_equal(rownames(lp), "r1 -- c1")
  expect_true(all(is.na(lp[,2:106])))
})

test_that("Test levelsize-normalisation for single link matrix", {
  lp <- link_positions(matrix(1,1,1), six_node = FALSE, weights = FALSE, normalisation = "levelsize")
  expect_equivalent(lp[1,1], 1)
  expect_equal(nrow(lp), 1)
  expect_equal(rownames(lp), "r1 -- c1")
  expect_true(all(is.na(lp[,2:29])))
})

test_that("Test levelsize-normalisation for single link matrix, six_node = TRUE", {
  lp <- link_positions(matrix(1,1,1), six_node = TRUE, weights = FALSE, normalisation = "levelsize")
  expect_equivalent(lp[1,1], 1)
  expect_equal(nrow(lp), 1)
  expect_equal(ncol(lp), 106)
  expect_equal(rownames(lp), "r1 -- c1")
  expect_true(all(is.na(lp[,2:106])))
})

test_that("Test motif-normalisation for single link matrix", {
  lp <- link_positions(matrix(1,1,1), six_node = FALSE, weights = FALSE, normalisation = "motif")
  expect_equivalent(lp[1,1], 1)
  expect_equal(nrow(lp), 1)
  expect_equal(rownames(lp), "r1 -- c1")
  expect_true(all(is.na(lp[,2:29])))
})

test_that("Test motif-normalisation for single link matrix, six_node = TRUE", {
  lp <- link_positions(matrix(1,1,1), six_node = TRUE, weights = FALSE, normalisation = "motif")
  expect_equivalent(lp[1,1], 1)
  expect_equal(nrow(lp), 1)
  expect_equal(ncol(lp), 106)
  expect_equal(rownames(lp), "r1 -- c1")
  expect_true(all(is.na(lp[,2:106])))
})


test_that("Test position-normalisation for single link matrix", {
  lp <- link_positions(matrix(1,1,1), six_node = FALSE, weights = FALSE, normalisation = "sizeclass")
  expect_equivalent(lp[1,1], 1)
  expect_equal(nrow(lp), 1)
  expect_equal(rownames(lp), "r1 -- c1")
  expect_true(all(is.na(lp[,2:29])))
})


test_that("Test sum-normalised counts for the motif 4", {
  lr <- link_positions( motifs[[4]], six_node = FALSE, weights = FALSE, normalisation = "sum" )
  expect_equal( as.vector(lr[,1]), rep(0.25,3))
  expect_equal( as.vector(lr[,2]), rep(0,3))
  expect_equal( as.vector(lr[,3]), rep(0.5,3))
  expect_equal( as.vector(lr[,4]), rep(0.25,3))
  expect_equal( as.vector(lr[,5]), rep(0,3))
})

test_that("Test sizeclass-normalised counts for the motif 4", {
  lr <- link_positions( motifs[[4]], six_node = FALSE, weights = FALSE, normalisation = "sizeclass" )
  expect_equal( as.vector(lr[,1]), rep(1,3))
  expect_equal( as.vector(lr[,2]), rep(0,3))
  expect_equal( as.vector(lr[,3]), rep(1,3))
  expect_equal( as.vector(lr[,4]), rep(1,3))
  expect_true(all(lr[,5:9] == 0))
  expect_true(all(is.na(lr[,10:29])))
})

test_that("Test position-normalised counts for the motif 4", {
  lr <- link_positions( motifs[[4]], six_node = FALSE, weights = FALSE, normalisation = "position" )
  expect_equal( as.vector(lr[,1]), rep(1/3,3))
  expect_true(all(is.na(lr[,2])))
  expect_equal( as.vector(lr[,3]), rep(1/3,3))
  expect_equal( as.vector(lr[,4]), rep(1/3,3))
  expect_true(all(is.na(lr[,5:29])))
})

test_that("Test levelsize-normalised counts for the motif 4", {
  lr <- link_positions( motifs[[4]], six_node = FALSE, weights = FALSE, normalisation = "levelsize" )
  expect_equal( as.vector(lr[,1]), rep(1,3))
  expect_true(all(is.na(lr[,2])))
  expect_equal( as.vector(lr[,3]), rep(1,3))
  expect_equal( as.vector(lr[,4]), rep(1,3))
  expect_true(all(is.na(lr[,5:29])))
})

test_that("Test motif-normalised counts for the motif 4", {
  lr <- link_positions( motifs[[4]], six_node = FALSE, weights = FALSE, normalisation = "motif" )
  expect_equal( as.vector(lr[,1]), rep(1,3))
  expect_true(all(is.na(lr[,2])))
  expect_equal( as.vector(lr[,3]), rep(1,3))
  expect_equal( as.vector(lr[,4]), rep(1,3))
  expect_true(all(is.na(lr[,5:29])))
})

test_that("Test sum-normalised counts for the motif 6", {
  lr <- link_positions( motifs[[6]], six_node = FALSE, weights = FALSE, normalisation = "sum" )
  expect_equal( as.vector(lr[,1]), rep(0.25,4))
  expect_equal( as.vector(lr[,2]), rep(0.25,4))
  expect_equal( as.vector(lr[,3]), rep(0.25,4))
  expect_equal( as.vector(lr[,4]), rep(0,4))
  expect_equal( as.vector(lr[,5]), rep(0,4))
  expect_equal( as.vector(lr[,8]), rep(0.25,4))
  expect_true(all(lr[,9:29] == 0))
})

test_that("Test sizeclass-normalised counts for motif 6", {
  lr <- link_positions( motifs[[6]], six_node = FALSE, weights = FALSE, normalisation = "sizeclass" )
  expect_equal( as.vector(lr[,1]), rep(1,4))
  expect_equal( as.vector(lr[,2]), rep(0.5,4))
  expect_equal( as.vector(lr[,3]), rep(0.5,4))
  expect_equal( as.vector(lr[,4]), rep(0,4))
  expect_equal( as.vector(lr[,5]), rep(0,4))
  expect_equal( as.vector(lr[,8]), rep(1,4))
  expect_equal( as.vector(lr[,9]), rep(0,4))
  expect_true(all(is.na(lr[,10:29])))
})

test_that("Test levelsize-normalised counts for motif 6", {
  lr <- link_positions( motifs[[6]], six_node = FALSE, weights = FALSE, normalisation = "levelsize" )
  expect_equal( as.vector(lr[,1]), rep(1,4))
  expect_equal( as.vector(lr[,2]), rep(1,4))
  expect_equal( as.vector(lr[,3]), rep(1,4))
  expect_true(all(is.na(lr[,4])))
  expect_true(all(lr[,5:7] == 0))
  expect_equal( as.vector(lr[,8]), rep(1,4))
  expect_true(all(is.na(lr[,9:29])))
})

test_that("Test motif-normalised counts for motif 6", {
  lr <- link_positions( motifs[[6]], six_node = FALSE, weights = FALSE, normalisation = "motif" )
  expect_equal( as.vector(lr[,1]), rep(1,4))
  expect_equal( as.vector(lr[,2]), rep(1,4))
  expect_equal( as.vector(lr[,3]), rep(1,4))
  expect_true(all(is.na(lr[,4])))
  expect_true(all(is.na(lr[,5:7])))
  expect_equal( as.vector(lr[,8]), rep(1,4))
  expect_true(all(is.na(lr[,9:29])))
})

test_that("Test sum-normalised counts for the motif 9", {
  lr <- link_positions( motifs[[9]], six_node = TRUE, weights = FALSE, normalisation = "sum")
  expect_equal( as.vector(lr[,1]), c(1/6, 1/6, 1/8, 1/5))
  expect_equal( as.vector(lr[,2]), c(0, 0, 1/8, 1/5))
  expect_equal( as.vector(lr[,3]), c(2/6, 2/6, 2/8, 0))
  expect_equal( as.vector(lr[,30]), rep(0,4))
})

test_that("Test sizeclass-normalised counts for the motif 9", {
  lr <- link_positions( motifs[[9]], six_node = TRUE, weights = FALSE, normalisation = "sizeclass" )
  expect_equal( as.vector(lr[,1]), rep(1,4))
  expect_equal( as.vector(lr[,2]), c(0,0, 1/3, 1))
  expect_equal( as.vector(lr[,3]), c(1,1, 2/3, 0))
  expect_equal( as.vector(lr[,4]), c(1/2, 1/2, 1/3, 0))
  expect_equal( as.vector(lr[,5]), c(0,0, 0, 1))
  expect_true(all(lr[,14:29] == 0))
  expect_true(all(is.na(lr[,30:106])))
})

test_that("Test positions-normalised counts for the motif 9", {
  lr <- link_positions( motifs[[9]], six_node = FALSE, weights = FALSE, normalisation = "position")
  expect_equal( as.vector(lr[,1]), rep(1/4,4))
  expect_equal( as.vector(lr[,2]), c(0,0, 1/2, 1/2))
  expect_equal( as.vector(lr[,3]), c(1/3,1/3, 1/3, 0))
  expect_equal( as.vector(lr[,4]), c(1/3, 1/3, 1/3, 0))
  expect_equal( as.vector(lr[,5]), c(0,0, 0, 1))
  expect_true(all(is.na(lr[,9:10])))
})

test_that("Test levelsize-normalised counts for the motif 9", {
  lr <- link_positions( motifs[[9]], six_node = FALSE, weights = FALSE, normalisation = "levelsize")
  expect_equal( as.vector(lr[,1]), rep(1,4))
  expect_equal( as.vector(lr[,5]), c(0,0,0,1))
  expect_equal( as.vector(lr[,6]), c(1,1,0,0))
  expect_equal( as.vector(lr[,7]), c(0,0,1,0))
  expect_equal( as.vector(lr[,8]), c(0,0,0,0))
  expect_true(all(is.na(lr[,9:10])))
})

test_that("Test motif-normalised counts for the motif 9", {
  lr <- link_positions( motifs[[9]], six_node = FALSE, weights = FALSE, normalisation = "motif")
  expect_equal( as.vector(lr[,1]), rep(1,4))
  expect_equal( as.vector(lr[,5]), c(0,0,0,1))
  expect_equal( as.vector(lr[,6]), c(1,1,0,0))
  expect_equal( as.vector(lr[,7]), c(0,0,1,0))
  expect_true(all(is.na(lr[,8])))
  expect_true(all(is.na(lr[,9:10])))
})

test_that("Test levelsize-normalised counts for the motif 10", {
  lr <- link_positions( motifs[[10]], six_node = FALSE, weights = FALSE, normalisation = "levelsize")
  expect_equal( as.vector(lr[,5]), c(0,0.5,0.5,0))
  expect_equal( as.vector(lr[,6]), c(1,0,0,1))
  expect_equal( as.vector(lr[,7]), c(0,0.5,0.5,0))
  expect_equal( as.vector(lr[,8]), c(0,0,0,0))
  expect_true(all(is.na(lr[,9:10])))
})

test_that("Test motif-normalised counts for the motif 10", {
  lr <- link_positions( motifs[[10]], six_node = FALSE, weights = FALSE, normalisation = "motif")
  expect_equal( as.vector(lr[,1]), rep(1,4))
  expect_equal( as.vector(lr[,5]), c(0,0.5,0.5,0))
  expect_equal( as.vector(lr[,6]), c(1,0,0,1))
  expect_equal( as.vector(lr[,7]), c(0,0.5,0.5,0))
  expect_true(all(is.na(lr[,8])))
  expect_true(all(is.na(lr[,9:10])))
})

# ------------------- TEST EXTRA NORMALISATION METHODS ----------------------

test_that("Test sizeclass_plus1", {
  lr <- link_positions(rbm(5,5), six_node = FALSE, weights = FALSE, normalisation = "sizeclass_plus1")
  expect(!any(is.na(lr)), failure_message = "failed")
})

test_that("Test sizeclass_NAzero", {
  lr <- link_positions(rbm(5,5), six_node = FALSE, weights = FALSE, normalisation = "sizeclass_NAzero")
  expect(!any(is.na(lr)), failure_message = "failed")
})

test_that("Test levelsize_plus1", {
  lr <- link_positions(rbm(5,5), six_node = FALSE, weights = FALSE, normalisation = "levelsize_plus1")
  expect(!any(is.na(lr)), failure_message = "failed")
})

test_that("Test levelsize_NAzero", {
  lr <- link_positions(rbm(5,5), six_node = FALSE, weights = FALSE, normalisation = "levelsize_NAzero")
  expect(!any(is.na(lr)), failure_message = "failed")
})

test_that("Test motif_plus1", {
  lr <- link_positions(rbm(5,5), six_node = FALSE, weights = FALSE, normalisation = "motif_plus1")
  expect(!any(is.na(lr)), failure_message = "failed")
})

test_that("Test motif_NAzero", {
  lr <- link_positions(rbm(5,5), six_node = FALSE, weights = FALSE, normalisation = "motif_NAzero")
  expect(!any(is.na(lr)), failure_message = "failed")
})
