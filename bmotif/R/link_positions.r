#' Calculate link position vectors
#'
#' Counts the number of times each link in a network occurs in each unique link position within the motifs
#' @param M A numeric matrix representing interactions between two groups of nodes. Each row corresponds to a node in one level
#' and each column corresponds to a node in the other level. Elements of M are positive numbers if nodes interact, and 0
#' otherwise. Formally, M is a biadjacency matrix. When nodes i and j interact, m_ij > 0; if they do not interact, m_ij = 0.
#'
#' @param six_node Logical; should positions in six node motifs be counted? Defaults to FALSE.
#' @param weights Logical; Should weights of the links be taken into account?
#' @param normalisation Which normalisation should be used: 'none','sum', 'position', 'sizeclass', 'sizeclass_plus1', 'sizeclass_NAzero', 'levelsize', 'levelsize_plus1', 'levelsize_NAzero','motif', 'motif_plus1' or 'motif_NAzero'?  Defaults to "none". (see details)
#' @details Counts the number of times each link in a network occurs in each of the 29 (if \code{six_node} = FALSE) or 106 (if \code{six_node} = TRUE) unique link positions within motifs (to quantify a link's structural role).
#' If \code{six_node} = FALSE, link positions in all motifs containing between 2 and 5 nodes are counted. If \code{six_node} = TRUE, link positions in all motifs containing between 2 and 6 nodes are counted. Analyses where \code{six_node} = FALSE are substantially faster
#' than when \code{six_node} = TRUE, especially for large networks. For large networks, counting six node motifs is also memory intensive. In some cases, R can crash if there is not enough memory.
#'
#' If interactions are weighted (non-zero matrix elements take values other than 1), these can be incorporated by setting \code{weights = TRUE}.
#' If \code{weights = TRUE}, the function will return the number of times each link occurs in each position,
#' multiplied by the weight of the link, following Mora et al. (2018).
#'
#' Links between nodes with more interactions will tend to appear in more positions. Normalisation helps control for this effect. bmotif include four types of normalisation:
#' \itemize{
#'  \item{\strong{"none"}: performs no normalisation and will return the raw position counts.}
#'  \item{\strong{"sum"}: divides the position measure for each link by the total number of times that link appears in any position (divides each element in a row by the row sum).}
#'  \item{\strong{"position"}: divides the position measure for each link by the total number of times any link occurs in that link position (divides each element in a column by the column sum). This gives a measure of how often a link occurs in a position relative to the other links in the network.}
#'  \item{\strong{Size class normalisation}
#'  \itemize{
#'  \item{\strong{"sizeclass"}: divides the position measure for each link by the total number of times that link appears in any position within the same motif size class (the number of nodes a motif contains).}
#'  \item{\strong{"sizeclass_plus1"}: same as 'sizeclass' but adds one to all position measure values. If a link does not occur in any motifs in a given size class, 'sizeclass' normalisation
#'  will return NAs. 'sizeclass_plus1' avoids this by adding one to all counts.}
#'  \item{\strong{"sizeclass_NAzero"}: same as 'sizeclass' but replaces all NA values with 0. If a link does not occur in any motifs in a given size class, 'sizeclass' normalisation
#'  will return NAs. 'sizeclass_NAzero' avoids this by replacing NAs with zero.}
#'  }
#'  }
#'  \item{\strong{Levelsize normalisation}
#'  \itemize{
#'  \item{\strong{"levelsize"}: divides the position measure for each link by the total number of times that link appears in any position within a motif with a given number of nodes in the top level and the bottom level.
#'  For example, the relative frequencies of all position measures in motifs with three nodes in the top level and two nodes in the bottom level will sum to one, as will the relative frequency of all position measures in motifs with 2 nodes in the top level and
#'  two nodes in the bottom level, and so on.}
#'  \item{\strong{"levelsize_plus1"}: same as 'levelsize' but adds one to all position measure values. If a link does not occur in any motifs with a given number of nodes in the top level and the bottom level, 'levelsize' normalisation
#'  will return NAs. 'levelsize_plus1' avoids this by adding one to all counts.}
#'  \item{\strong{"levelsize_NAzero"}: same as 'levelsize' but replaces all NA values with 0. If a link does not occur in any motifs with a given number of nodes in the top level and the bottom level, 'levelsize' normalisation
#'  will return NAs. 'levelsize_NAzero' avoids this by replacing NAs with zero.}
#'  }
#'  }
#'  \item{\strong{Motif normalisation}
#'  \itemize{
#'  \item{\strong{"motif"}: divides the position measure for each link by the total number of times that link appears in any position within the same motif.
#'  For example, the relative frequencies of all position measures in motif 5 will sum to one, as will the relative frequency of all position measures in motif 10, and so on.}
#'  \item{\strong{"motif_plus1"}: same as 'motif' but adds one to all position measure values. If a link does not occur in a particular motif, 'motif' normalisation
#'  will return NAs. 'motif_plus1' avoids this by adding one to all counts.}
#'  \item{\strong{"motif_NAzero"}: same as 'motif' but replaces all NA values with 0. If a link does not occur in a particular motif, 'levelsize' normalisation
#'  will return NAs. 'motif_NAzero' avoids this by replacing NAs with zero.}
#'  }
#'  }
#'  }
#'
#' If a matrix is provided without row or column names, default names will be assigned: the first row will be called called 'r1', the second row will be called 'r2' and so on. Similarly, the first column will be called 'c1', the second column will be called 'c2' and so on.
#'
#' @return  Returns a data frame with one column for each link position: 29 columns if \code{six_node} is FALSE, and 106 columns if \code{six_node} is TRUE.
#' Columns names are given as "lpx" where x is the ID of the position as described in the motif dictionary. \strong{To view the 'motif dictionary' showing
#' which link position a given ID corresponds to, enter \code{vignette("bmotif-dictionary")}.}
#'
#' Each row corresponds to one link in the network. Row names are gives as "x -- y", where x is the species in the first level (rows) and y is the species in the second level (columns).
#'
#' By default, the elements of the data frame will be the raw link position counts.
#' If \code{weight = TRUE}, link position counts will be multiplied by the link weight.
#' If \code{normalisation} is set to "sum", "sizeclass" or "position", the elements will be
#' normalised position counts as described above.
#'
#' @export
#' @references
#' Mora, B.B., Cirtwill, A.R. and Stouffer, D.B., 2018. pymfinder: a tool for the motif analysis of binary and quantitative complex networks. bioRxiv, 364703.
#'
#' Simmons, B. I., Sweering, M. J. M., Dicks, L. V., Sutherland, W. J. and Di Clemente, R. bmotif: a package for counting motifs in bipartite networks. bioRxiv. doi: 10.1101/302356
#' @examples
#' set.seed(123)
#' row <- 10
#' col <- 10
#'
#' # link positions in a binary network
#' m <- matrix(sample(0:1, row*col, replace=TRUE), row, col)
#' link_positions(M = m, six_node = TRUE, weights = FALSE, normalisation = "none")
#'
#' # link positions in a weighted network
#' m[m>0] <- stats::runif(sum(m), 0, 100)
#' link_positions(M = m, six_node = TRUE, weights = TRUE, normalisation = "none")

link_positions <- function(M, six_node = FALSE, weights, normalisation = "none") {
  # M is the adjacence matrix

  # check inputs
  if(inherits(M, "matrix") != TRUE){stop("'M' must be an object of class 'matrix'")} # make sure M is a matrix
  if(!all(apply(M, 1:2, is.numeric))){stop("Elements of 'M' must be numeric")} # make sure all elements of M are numbers
  if(!all(apply(M, 1:2, function(x) length(x) > 0))){stop("Elements of 'M' cannot have length zero")} # make sure no elements of M have 0 length e.g. numeric(0)
  if(!all(apply(M, 1:2, function(x) x >= 0))){stop("Elements of 'M' must be greater than or equal to zero")} # make sure all elements of M are >= zero
  if(all(apply(M, 1:2, function(x) x == 0))){stop("The matrix has no links (all elements are zero)")} # make sure M has some links
  if(inherits(six_node,"logical") != TRUE){stop("'six_node' must be of class 'logical' i.e. TRUE or FALSE")} # make sure six_node is logical i.e. TRUE or FALSE
  if(inherits(normalisation,"character") != TRUE){stop("'normalisation' must be of class 'character'")} # make sure 'normalisation' is a character
  if(!normalisation %in% c("none","sum","sizeclass", "sizeclass_plus1", "sizeclass_NAzero","levelsize", "levelsize_plus1", "levelsize_NAzero", "motif", "motif_plus1", "motif_NAzero","position")){
    stop("'normalisation' must equal 'none','sum', 'position', 'sizeclass', 'sizeclass_plus1', 'sizeclass_NAzero', 'levelsize', 'levelsize_plus1', 'levelsize_NAzero','motif', 'motif_plus1' or 'motif_NAzero'")
  }
  if(any(duplicated(rownames(M))) | any(duplicated(colnames(M)))){stop("Input matrix must not have duplicate row or column names")}
  if((normalisation == "sum" | normalisation == "sizeclass" | normalisation == "levelsize" | normalisation == "motif") & weights == TRUE){warning("Please note that with 'sum', 'sizeclass', 'levelsize' or 'motif' normalisation, the results won't change if weights are taken into account. This is because when weights = TRUE, each row of the output is multiplied by a fixed factor (the link's weight) and therefore the relative proportions are the same as if weights were not considered. Consider setting normalisation to 'position' or 'none'.")}

  W <- M # store copy of weighted matrix
  M[M > 0] <- 1 # ensure M is binary

  # compute auxiliary stuff
  NZ = nrow(M) # this is just Z in the latex file, but Z is also a matrix later
  NP = ncol(M) # this is just P in the latex file, but P is also a matrix later
  jZ = as.matrix(rep(1, NZ))
  jP = as.matrix(rep(1, NP))

  J <- matrix(rep(1, NZ * NP), nrow = NZ, ncol = NP) # Z x P matrix of ones
  JP <- matrix(rep(1, NP * NP), nrow = NP, ncol = NP) # P x P matrix of ones
  JZ <- matrix(rep(1, NZ * NZ), nrow = NZ, ncol = NZ) # Z x Z matrix of ones

  JT <- t(J) # transpose of J, P x Z matrix of ones
  MT <- t(M) # transpose of M
  N <- J - M # complement of M
  NT <- t(N) # transpose of the complement of M

  dZ = M %*% jP
  dP = MT %*% jZ

  V = matrix(rep(dZ, each = NP), nrow = NP, ncol = NZ, byrow = FALSE)
  U = matrix(rep(dP, times = NZ), nrow = NP, ncol = NZ, byrow = FALSE)

  # matrices like in the formulae
  Z <- M %*% MT
  Y <- M %*% NT
  X <- N %*% MT

  P <- MT %*% M
  Q <- MT %*% N
  R <- NT %*% M

  ##########################

  # now compute positions
  # possibly put this in a separate function later

  # create a data frame
  # each row corresponds to one edge in the graph
  # order: (1,1), (1,2), (1,3), ...(1,P), (2,1), ..., (2,P), ..., (Z,1), ..., (Z,P)
  # each column corresponds to an edge position, i.e. column p3 gives count for link position 3

  # note: the as.vector() command is not strictly needed
  # since R even converts the matrix automatically into a vector if we assign a matrix to a column like below

  k = 29 + (six_node) * (106 - 29) # k depends on whether 6 species are wanted or not

  lp <- matrix(rep (-1, k* NZ * NP), ncol = k)

  ######

  names <- rep("", NZ * NP)

  if (is.null(rownames(M))) {
    class1_names <- paste0("r", 1:NZ)
  } else {
    class1_names <- rownames(M)
  }
  if (is.null(colnames(M))) {
    class2_names <- paste0("c", 1:NP)
  } else {
    class2_names <- colnames(M)
  }

  for (i in 1:NZ) {
    for (j in 1:NP) {
      names[(i - 1)* NP + j] <- name_edge(class1_names[i], class2_names[j])
    }
  }
  rownames(lp) <- names

  if (six_node) {
    colnames(lp) <- paste('lp', 1:106, sep = '')
  }
  if (!six_node) {
    colnames(lp) <- paste('lp', 1:29, sep = '')
  }

  # want to compute only for present edges
  # so make a list of all the edges present in the network
  # column one is just a numbering
  # column 2 is the index of the row vertex, column 3 the index of the column vertex
  # for later use, if weights = TRUE, also store the weight of the link in column 4
  ct <- 0
  E <- sum(M) # number of edges

  if (six_node & !weights) {
    edges <- matrix(rep(0,3*E), nrow = E)
    for (i in 1:NZ) {
      for (j in 1:NP) {
        if (M[i,j] == 1) {
          ct <- ct + 1
          edges[ct,1] <- ct
          edges[ct,2:3] <- c(i, j)
        }
      }
    }
  }

  if (weights) {
    edges <- matrix(rep(0,4*E), nrow = E)
    for (i in 1:NZ) {
      for (j in 1:NP) {
        if (M[i,j] == 1) {
          ct <- ct + 1
          edges[ct,1] <- ct
          edges[ct,2:3] <- c(i, j)
          edges[ct, 4] <- W[i,j]
        }
      }
    }
  }

  # ------- PUT IN VALUES
  # up to 4 species

  lp[,1] <- as.vector(t(M))
  lp[,2] <- as.vector(V - JT)

  lp[,3] <- as.vector(U - JT)
  lp[,4] <- 0.5*as.vector((U - JT)*(U-2*JT))
  lp[,5] <- as.vector(NT %*% Z)
  lp[,6] <- as.vector(MT %*% Y)
  lp[,7] <- as.vector(MT %*% X)
  lp[,8] <- as.vector(MT %*% (Z - JZ) - V + JT)
  lp[,9] <- 0.5*as.vector((V - JT)*(V - 2*JT))

  # 5 species

  lp[,10] <- 1/6 * as.vector((U - JT)*(U - 2*JT) * (U - 3*JT))
  lp[,11] <- 0.5 * as.vector(t((M %*% (Q * (Q - JP)))))
  lp[,12] <- as.vector(t(N %*% (P * (R - JP))))
  lp[,13] <- 0.5 * as.vector(t(M %*% (R * (R - JP))))
  lp[,14] <- as.vector(t(N %*% (P * Q)))
  lp[,15] <- as.vector(t(M %*% (Q * R)))
  lp[,16] <- as.vector(t(M %*% (Q * (P - JP))))
  lp[,17] <- 0.5 * as.vector(t(N %*% (P * (P - JP))))
  lp[,18] <- as.vector(t(M %*% (R * (P - JP))))
  lp[,19] <- 0.5 * as.vector(t(M %*% ( (P - JP) * (P - 2*JP))) - (U - JT) * (U - 2*JT) )
  lp[,20] <- as.vector(NT %*% (Z * (X - JZ)))
  lp[,21] <- 0.5 * as.vector(MT %*% (Y * (Y - JZ)))
  lp[,22] <- 0.5 * as.vector(MT %*% (X * (X - JZ)))
  lp[,23] <- as.vector(NT %*% (Z * Y))
  lp[,24] <- as.vector(MT %*% (X * Y))
  lp[,25] <- 0.5 * as.vector(NT %*% (Z * (Z - JZ)))
  lp[,26] <- as.vector(MT %*% (Y * (Z - JZ)))
  lp[,27] <- as.vector(MT %*% (X * (Z - JZ)))
  lp[,28] <- 0.5 * as.vector(MT %*% ( (Z - JZ) * (Z - 2*JZ) ) - (V - JT) * (V - 2*JT))
  lp[,29] <- 1/6 * as.vector((V - JT) * (V - 2*JT) * (V - 3*JT))

  # 6 species
  if (six_node) {
    lp[,30] <- 1/24 * as.vector((U - JT) * (U - 2*JT) * (U - 3*JT) * (U - 4*JT))
    lp[,31] <- 1/6 * as.vector(t(M %*% (Q * (Q - JP) * (Q - 2*JP))))
    lp[,32] <- 1/2 * as.vector(t(N %*% ((R - JP) * (R - 2*JP) * P)))
    lp[,33] <- 1/6 * as.vector(t(M %*% (R * (R - JP) * (R - 2*JP))))
    lp[,34] <- 1/2 * as.vector(t(N %*% (Q * (Q - JP) * P)))
    lp[,35] <- 1/2 * as.vector(t(M %*% (Q * (Q - JP) * R)))
    lp[,36] <- as.vector(t(N %*% (Q * P * (R - JP))))
    lp[,37] <- 1/2 * as.vector(t(M %*% (R * (R - JP) * Q)))
    lp[,38] <- 1/2 * as.vector(t(M %*% (Q * (Q - JP) * (P - JP))))
    lp[,39] <- 1/2 * as.vector(t(N %*% (P * (P - JP) * (R - JP))))
    lp[,40] <- 1/2 * as.vector(t(M %*% (R * (R - JP) * (P - JP))))
    lp[,41] <- 1/2 * as.vector(t(N %*% (P * (P - JP) * Q)))
    lp[,42] <- as.vector(t(M %*% ((P - JP) * Q * R)))
    lp[,43] <- 1/2 * as.vector(t(M %*% ((P- JP) * (P - 2*JP) * Q)))
    lp[,44] <- 1/6 * as.vector(t(N %*% (P* (P - JP) * (P - 2*JP))))
    lp[,45] <- 1/2 * as.vector(t(M %*% ((P- JP) * (P - 2*JP) * R)))
    lp[,46] <- 1/6 * as.vector(t(M %*% ((P- JP) * (P - 2*JP) * (P - 3* JP))) - ((U - JT) * (U - 2*JT) * (U - 3*JT)))

    lp[,90] <- 1/2 * as.vector(NT %*% ((X - JZ) * (X - 2*JZ) * Z))
    lp[,91] <- 1/6 * as.vector(MT %*% (Y * (Y - JZ) * (Y - 2*JZ)))
    lp[,92] <- 1/6 * as.vector(MT %*% (X * (X - JZ) * (X - 2*JZ)))
    lp[,93] <- 1/2 * as.vector(NT %*% (Z * Y * (Y - JZ)))
    lp[,94] <- as.vector(NT %*% (Z * Y * (X - JZ)))
    lp[,95] <- 0.5 * as.vector(MT %*% (X * Y * (Y - JZ)))
    lp[,96] <- 0.5 * as.vector(MT %*% (Y * X * (X - JZ)))
    lp[,97] <- 0.5 * as.vector(NT %*% ((X - JZ) * Z * (Z - JZ)))
    lp[,98] <- 0.5 * as.vector(MT %*% ((Z - JZ) * Y * (Y - JZ)))
    lp[,99] <- 0.5 * as.vector(MT %*% ((Z - JZ) * X * (X - JZ)))
    lp[,100] <- 0.5 * as.vector(NT %*% (Z * (Z - JZ) * Y))
    lp[,101] <- as.vector(MT %*% (X * Y * (Z - JZ)))
    lp[,102] <- 1/6 * as.vector(NT %*% (Z * (Z - JZ) * (Z - 2*JZ)))
    lp[,103] <- 0.5 * as.vector(MT %*% (Y * (Z - JZ) * (Z - 2*JZ)))
    lp[,104] <- 0.5 * as.vector(MT %*% (X * (Z - JZ) * (Z - 2*JZ)))
    lp[,105] <- 1/6 * as.vector(MT %*% ((Z - JZ) * (Z - 2*JZ) * (Z - 3*JZ)) - ((V - JT) * (V - 2*JT) * (V - 3*JT)))
    lp[,106] <- 1/24 * as.vector((V - JT) * (V - 2*JT) * (V - 3*JT) * (V - 4*JT))

    # ---------- TENSOR CASES

    # now clean rows where the edge is not present
    lp <- subset(lp, lp[,1] == 1)

    if (NP >= NZ) { # NP > NZ
      JZ3 <- array(rep(1, NZ * NZ * NZ), c(NZ, NZ, NZ)) # Z x Z x Z array of ones
      KZ3 <- JZ3 # Z x Z x Z matrix, 0 if two indices are equal, 1 otherwise
      for (i in 1 : NZ){
        for (j in 1 : NZ){
          KZ3[i,j,j] <- 0
          KZ3[j,i,j] <- 0
          KZ3[j,j,i] <- 0
        }
      }

      # auxiliary matrices to create A, ..., G
      AP <- maketensor(MT, MT)
      BP <- maketensor(MT, NT)
      CP <- maketensor(NT, MT)
      DP <- maketensor(NT, NT)

      # the following are like the matrices A to G in the Latex File
      # multiply with KZ3 accounts for the terms (1 - delta_{ij}) etc.
      A_ <- tensor::tensor(M, AP, 2, 1) * KZ3
      B_ <- tensor::tensor(M, BP, 2, 1) * KZ3
      C_ <- tensor::tensor(M, CP, 2, 1) * KZ3
      D_ <- tensor::tensor(M, DP, 2, 1) * KZ3
      E_ <- tensor::tensor(N, AP, 2, 1) * KZ3
      F_ <- tensor::tensor(N, BP, 2, 1) * KZ3
      G_ <- tensor::tensor(N, CP, 2, 1) * KZ3

      ####################

      lp[,47] <- 0.5 * large_tensor(DP, A_ * D_ - A_, edges)
      lp[,48] <- 0.5 * large_tensor(AP, G_ * (G_ - JZ3), edges)
      lp[,49] <- 0.25 * large_tensor(AP, D_ * (D_ - JZ3), edges)
      lp[,50] <- large_tensor(BP, D_ * C_, edges)
      lp[,51] <- large_tensor(BP, E_ * F_, edges)
      lp[,52] <- 0.5 * large_tensor(DP, B_ * C_, edges)
      lp[,53] <- large_tensor(AP, D_ * G_, edges)
      lp[,54] <- 0.5 * large_tensor(AP, F_ * G_, edges)
      lp[,55] <- large_tensor(DP, A_ * G_, edges)
      lp[,56] <- large_tensor(DP, B_ * E_, edges)
      lp[,57] <- large_tensor(BP, E_ * G_, edges)
      lp[,58] <- large_tensor(BP, C_ * G_, edges)
      lp[,59] <- large_tensor(BP, D_ * E_, edges)
      lp[,60] <- large_tensor(BP, C_ * F_, edges)
      lp[,61] <- large_tensor(DP, A_ * B_, edges)
      lp[,62] <- large_tensor(CP, A_ * G_, edges)
      lp[,63] <- large_tensor(BP, A_ * D_, edges)
      lp[,64] <- large_tensor(AP, C_ * G_, edges)
      lp[,65] <- large_tensor(AP, B_ * D_, edges)
      lp[,66] <- large_tensor(AP, E_ * G_, edges)
      lp[,67] <- 0.5 * large_tensor(DP, A_ * E_, edges)
      lp[,68] <- large_tensor(BP, A_ * G_, edges)
      lp[,69] <- 0.5 * large_tensor(AP, D_ * E_, edges)
      lp[,70] <- large_tensor(AP, C_ * F_, edges)
      lp[,71] <- 0.25 * large_tensor(DP, A_ * (A_ - JZ3), edges)
      lp[,72] <- large_tensor(AP, (A_ - JZ3) * F_, edges)
      lp[,73] <- 0.5 * large_tensor(AP, (A_ - JZ3) * D_, edges)
      lp[,74] <- 0.5 * large_tensor(BP, E_ * (E_ - JZ3), edges)
      lp[,75] <- 0.5 * large_tensor(BP, C_ * (C_ - JZ3), edges)
      lp[,76] <- large_tensor(BP, (B_ - JZ3) * E_, edges)
      lp[,77] <- large_tensor(BP, (B_ - JZ3) * C_, edges)
      lp[,78] <- large_tensor(BP, A_ * (B_ - JZ3), edges)
      lp[,79] <- 0.25 * large_tensor(AP, E_ * (E_ - JZ3), edges)
      lp[,80] <- 0.5 * large_tensor(AP, B_ * (B_ - JZ3), edges)
      lp[,81] <- large_tensor(BP, C_ * E_, edges)
      lp[,82] <- large_tensor(BP, A_ * E_, edges)
      lp[,83] <- large_tensor(BP, A_ * C_, edges)
      lp[,84] <- large_tensor(AP, B_ * E_, edges)
      lp[,85] <- 0.5 * large_tensor(AP, B_ * C_, edges)
      lp[,86] <- 0.5 * large_tensor(BP, A_ * (A_ - JZ3), edges)
      lp[,87] <- 0.5 * large_tensor(AP, (A_ - JZ3) * E_, edges)
      lp[,88] <- large_tensor(AP, (A_ - JZ3) * B_, edges)
      lp[,89] <- 0.25 * large_tensor(AP, (A_ - JZ3) * (A_ - 2*JZ3) * KZ3, edges)

    }
    if (NP < NZ) { # NP <= NZ
      JP3 <- array(rep(1, NP * NP * NP), c(NP, NP, NP)) # P x P x P array of ones
      KP3 <- JP3 # P x P x P matrix, 0 if two indices are equal, 1 otherwise
      for (i in 1 : NP){
        for (j in 1 : NP){
          KP3[i,j,j] <- 0
          KP3[j,i,j] <- 0
          KP3[j,j,i] <- 0
        }
      }

      ####
      # those are needed if z > p
      AZ <- maketensor(M, M)
      BZ <- maketensor(M, N)
      CZ <- maketensor(N, M)
      DZ <- maketensor(N, N)

      # the following are like the matrices A' to G' in the Latex File
      # multiply with KP3 accounts for the terms (1 - delta_{ij}) etc.
      A_2 <- tensor::tensor(M, AZ, 1, 1) * KP3
      B_2 <- tensor::tensor(M, BZ, 1, 1) * KP3
      C_2 <- tensor::tensor(M, CZ, 1, 1) * KP3
      D_2 <- tensor::tensor(M, DZ, 1, 1) * KP3
      E_2 <- tensor::tensor(N, AZ, 1, 1) * KP3
      F_2 <- tensor::tensor(N, BZ, 1, 1) * KP3
      G_2 <- tensor::tensor(N, CZ, 1, 1) * KP3

      ####
      lp[,47] <- 0.5 * large_tensor(F_2 * (F_2 - JP3), AZ, edges)
      lp[,48] <- 0.5 * large_tensor(A_2 * (D_2 - JP3), DZ, edges)
      lp[,49] <- 0.25 * large_tensor(D_2 * (D_2 - JP3), AZ, edges)
      lp[,50] <- large_tensor(D_2 * G_2, AZ, edges)
      lp[,51] <- large_tensor(A_2 * G_2, DZ, edges)
      lp[,52] <- 0.5 * large_tensor(F_2 * G_2, AZ, edges)
      lp[,53] <- large_tensor(C_2 * D_2, BZ, edges)
      lp[,54] <- 0.5 * large_tensor(B_2 * C_2, DZ, edges)
      lp[,55] <- large_tensor(E_2 * F_2, BZ, edges)
      lp[,56] <- large_tensor(E_2 * G_2, BZ, edges)
      lp[,57] <- large_tensor(B_2 * E_2, DZ, edges)
      lp[,58] <- large_tensor(D_2 * E_2, BZ, edges)
      lp[,59] <- large_tensor(C_2 * G_2, BZ, edges)
      lp[,60] <- large_tensor(C_2 * F_2, BZ, edges)
      lp[,61] <- large_tensor(E_2 * G_2, AZ, edges)
      lp[,62] <- large_tensor(A_2 * F_2, BZ, edges)
      lp[,63] <- large_tensor(C_2 * G_2, AZ, edges)
      lp[,64] <- large_tensor(A_2 * D_2, BZ, edges)
      lp[,65] <- large_tensor(B_2 * D_2, AZ, edges)
      lp[,66] <- large_tensor(A_2 * B_2, DZ, edges)
      lp[,67] <- 0.5 * large_tensor(E_2 * (E_2 - JP3), BZ, edges)
      lp[,68] <- large_tensor((B_2 - JP3) * E_2, BZ, edges)
      lp[,69] <- 0.5 * large_tensor(C_2 * (C_2 - JP3), BZ, edges)
      lp[,70] <- large_tensor(B_2 * (C_2 - JP3), CZ, edges)
      lp[,71] <- 0.25* large_tensor(E_2 * (E_2 - JP3), AZ, edges)
      lp[,72] <- large_tensor(A_2 * (C_2 - JP3), CZ, edges)
      lp[,73] <- 0.5 * large_tensor(C_2 * (C_2 - JP3), AZ, edges)
      lp[,74] <- 0.5 * large_tensor(A_2 * E_2, DZ, edges)
      lp[,75] <- 0.5 * large_tensor(D_2 * E_2, AZ, edges)
      lp[,76] <- large_tensor(A_2 * G_2, BZ, edges)
      lp[,77] <- large_tensor(B_2 * G_2, AZ, edges)
      lp[,78] <- large_tensor((A_2 - JP3) * G_2, AZ, edges)
      lp[,79] <- 0.25 * large_tensor(A_2* (A_2 - JP3), DZ, edges)
      lp[,80] <- 0.5 * large_tensor((A_2 - JP3) * D_2, AZ, edges)
      lp[,81] <- large_tensor(C_2 * E_2, BZ, edges)
      lp[,82] <- large_tensor(A_2 * E_2, BZ, edges)
      lp[,83] <- large_tensor(B_2 * E_2, AZ, edges)
      lp[,84] <- large_tensor(A_2 * C_2, BZ, edges)
      lp[,85] <- 0.5 * large_tensor(B_2 * C_2, AZ, edges)
      lp[,86] <- 0.5 * large_tensor((A_2 - JP3) * E_2, AZ, edges)
      lp[,87] <- 0.5 * large_tensor(A_2 * (A_2 - JP3), BZ, edges)
      lp[,88] <- large_tensor((A_2 - JP3) * C_2, AZ, edges)
      lp[,89] <- 0.25 * large_tensor((A_2 - JP3) * (A_2 - 2*JP3) * KP3, AZ, edges)
    }
  }

  #######################
  # now clean rows where the edge is not present
  lp <- subset(lp, lp[,1] == 1)

  # --------------- CONSIDER WEIGHTS -------------------------------------
  if (weights) {
    for (k in 1:nrow(lp)) {
      lp[k,] <- lp[k,] * edges[k,4]
    }
  }

  # --------------- NORMALISE --------------------------------------------

  if (normalisation == "none") {
    return(as.data.frame(lp))
  } else {
    norm_lp <- normalise_link_positions(lp, type = normalisation, six_node = six_node)
    return(as.data.frame(norm_lp))
  }

}
