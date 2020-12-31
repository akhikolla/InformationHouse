node_sets <- function(M, six_node){
  # clean matrix
  M[M > 0] <- 1 # ensure M is binary
  dimnames(M) <- NULL # strip row and column names

  # create results container
  sets <- stats::setNames(object = rep(NA, 17), nm = paste0("m",1:17))

  # record rows and columns
  nr <- nrow(M)
  nc <- ncol(M)

  # calculate possible node sets, given in the format m_rows_columns in the motif
  m1_1 <- choose(nr,1) * choose(nc,1)

  m1_2 <- choose(nr,1) * choose(nc,2)
  m2_1 <- choose(nr,2) * choose(nc,1)

  m3_1 <- choose(nr,3) * choose(nc,1)
  m2_2 <- choose(nr,2) * choose(nc,2)
  m1_3 <- choose(nr,1) * choose(nc,3)

  m4_1 <- choose(nr,4) * choose(nc,1)
  m3_2 <- choose(nr,3) * choose(nc,2)
  m2_3 <- choose(nr,2) * choose(nc,3)
  m1_4 <- choose(nr,1) * choose(nc,4)

  if(six_node == TRUE){
    m5_1 <- choose(nr,5) * choose(nc,1)
    m4_2 <- choose(nr,4) * choose(nc,2)
    m3_3 <- choose(nr,3) * choose(nc,3)
    m2_4 <- choose(nr,2) * choose(nc,4)
    m1_5 <- choose(nr,1) * choose(nc,5)
  }

  if(six_node == TRUE){
    sets <- stats::setNames(c(m1_1,
                       m1_2,
                       m2_1,
                       m3_1,
                       rep(m2_2,2),
                       m1_3,
                       m4_1,
                       rep(m3_2,4),
                       rep(m2_3,4),
                       m1_4,
                       m5_1,
                       rep(m4_2,6),
                       rep(m3_3,13),
                       rep(m2_4,6),
                       m1_5),
                     nm = paste0("m", 1:44))
  } else {
    sets <- stats::setNames(c(m1_1,
                       m1_2,
                       m2_1,
                       m3_1,
                       rep(m2_2,2),
                       m1_3,
                       m4_1,
                       rep(m3_2,4),
                       rep(m2_3,4),
                       m1_4),
                     nm = paste0("m", 1:17))
  }

  # output
  return(sets)
}
