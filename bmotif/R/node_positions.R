#' Calculate node position vectors
#'
#' For binary networks, counts the number of times each node appears in each unique node position within motifs; for weighted networks calculates a range of weighted node position measures.
#' @param M A numeric matrix representing interactions between two groups of nodes. Each row corresponds to a node in one level
#' and each column corresponds to a node in the other level. Elements of M are positive numbers if nodes do interact, and 0
#' otherwise. Formally, M is a biadjacency matrix. When nodes i and j interact, m_ij > 0; if they do not interact, m_ij = 0.
#' If interactions are weighted (non-zero matrix elements take values other than 1), the function will automatically convert the matrix to a binary
#' matrix.
#' @param six_node Logical; should six node motifs be counted? Defaults to FALSE.
#' @param level Which node level should positions be calculated for: "rows", "columns" or "all"?  Defaults to "all".
#' @param weights_method The method for calculating weighted positions; must be one of 'none', 'mean_motifweights', 'total_motifweights', 'mean_nodeweights', 'total_nodeweights', 'contribution', 'mora' or 'all' (see details).
#' @param weights_combine Method for combining weighted position measures; must be one of 'none', 'mean' or 'sum' (see details). Defaults to 'none'.
#' @param normalisation Which normalisation should be used: "none","sum","sizeclass", "sizeclass_plus1", "sizeclass_NAzero", "position","levelsize","levelsize_plus1","levelsize_NAzero","motif","motif_plus1" or "motif_NAzero" (see details)?  Defaults to "none".
#' @details For binary networks, counts the number of times each node occurs in each unique position within motifs (to quantify a node's structural role).
#' If networks are weighted, \code{node_positions} can also calculate various weighted node position measures.
#'
#' If a matrix is provided without row or column names, default names will be assigned: the first row will be called called 'r1', the second row will be called 'r2' and so on. Similarly, the first column will be called 'c1', the second column will be called 'c2' and so on.
#'
#' \strong{Six node}
#'
#' If \code{six_node} = FALSE, node positions in all motifs containing between 2 and 5 nodes are counted. If \code{six_node} = TRUE, node positions in all motifs containing between 2 and 6 nodes are counted.
#' Analyses where \code{six_node} = FALSE are substantially faster than when \code{six_node} = TRUE, especially for large networks. For large networks, counting six node motifs is also memory intensive.
#' In some cases, R can crash if there is not enough memory.
#'
#' \strong{Level}
#'
#' The \code{level} argument controls which level of nodes positions are calculated for: "rows" returns position counts for all nodes in rows, "columns"
#' returns position counts for all nodes in columns, and "all" return counts for all nodes in the network.
#'
#' \strong{Weighted networks}
#'
#' \code{node_positions} also supports weighted networks for motifs up to five nodes. Weighted analyses are controlled using two arguments: \code{weights_method} and \code{weights_combine}. These are described in detail below:
#'
#' \itemize{\item{\strong{'weights_method'}: determines how the weighted position of a node is calculated in each motif occurrence.
#'
#' \itemize{
#' \item{\strong{'none'}: weights are ignored and \code{node_positions} returns the frequency with which each node occurs in each unique position within motifs. \code{'weights_combine'} must also be 'none'.}
#' \item{\strong{'mean_motifweights'}: for a given node in a given position in a motif occurrence (formally a subgraph isomorphic to a particular motif), returns the mean weight of that motif occurrence
#'  i.e. the mean of all link strengths in that motif occurrence.}
#' \item{\strong{'total_motifweights'}: for a given node in a given position in a motif occurrence (formally a subgraph isomorphic to a particular motif), returns the total weight of that motif occurrence
#'  i.e. the sum of all link strengths in that motif occurrence.}
#' \item{\strong{'mean_nodeweights'}: for a given node in a given position in a motif occurrence (formally a subgraph isomorphic to a particular motif), returns the mean weight of that focal node's links.}
#' \item{\strong{'total_nodeweights'}: for a given node in a given position in a motif occurrence (formally a subgraph isomorphic to a particular motif), returns the total weight of that focal node's links.}
#' \item{\strong{'contribution'}: for a given node in a given position in a motif occurrence (formally a subgraph isomorphic to a particular motif), returns the total weight of that focal node's links as a proportion of the total
#' weight of that motif occurrence. i.e. the sum of the focal node's links divided by the sum of all link strengths in that motif occurrence.}
#' \item{\strong{'mora'}: calculates a contribution measure following Mora et al. (2018).}
#' \item{\strong{'all'}: calculates all the above measures (except 'none') and returns them as a list of length five.}
#' }
#' }
#' \item{\strong{'weights_combine'}: determines how weighted position measures are combined across motif occurrences to give an overall measure for a each node in a each position.
#' \itemize{
#' \item{\strong{'none'}: weights are ignored and \code{node_positions} returns the frequency with which each node occurs in each unique positions within motifs. \code{'weights_method'} must also be 'none'.}
#' \item{\strong{'sum'}: weighted measures are summed across occurrences.}
#' \item{\strong{'mean'}: the mean of the weighted measure across occurrences is calculated.}
#' }
#' }
#' }
#'
#' \strong{Normalisation}
#'
#' Nodes with more interactions will tend to appear in more positions. Normalisation helps control for this effect. bmotif include six main types of normalisation:
#' \itemize{
#'  \item{\strong{"none"}: performs no normalisation and will return the raw position measure}
#'  \item{\strong{"sum"}: divides the position measure for each node by the total number of times that node appears in any position (divides each element in a row by the row sum).}
#'  \item{\strong{"position"}: divides the position measure for each node by the total number of times any node occurs in that node position (divides each element in a column by the column sum). This gives a measure of how often a node occurs in a position relative to the other nodes in the network.}
#'  \item{\strong{Size class normalisation}
#'  \itemize{
#'  \item{\strong{"sizeclass"}: divides the position measure for each node by the total number of times that node appears in any position within the same motif size class (the number of nodes a motif contains).}
#'  \item{\strong{"sizeclass_plus1"}: same as 'sizeclass' but adds one to all position measures If a species does not occur in any motifs in a given size class, 'sizeclass' normalisation
#'  will return NAs. 'sizeclass_plus1' avoids this by adding one to all counts.}
#'  \item{\strong{"sizeclass_NAzero"}: same as 'sizeclass' but replaces all NA values with 0. If a species does not occur in any motifs in a given size class, 'sizeclass' normalisation
#'  will return NAs. 'sizeclass_NAzero' avoids this by replacing NAs with zero.}
#'  }
#'  }
#'  \item{\strong{Levelsize normalisation}
#'  \itemize{
#'  \item{\strong{"levelsize"}: divides the position measure for each node by the total number of times that node appears in any position within a motif with a given number of nodes in the top level and the bottom level.
#'  For example, the relative frequencies of all position measures in motifs with three nodes in the top level and two nodes in the bottom level will sum to one, as will the relative frequency of all position counts in motifs with 2 nodes in the top level and
#'  two nodes in the bottom level, and so on.}
#'  \item{\strong{"levelsize_plus1"}: same as 'levelsize' but adds one to all position measures If a species does not occur in any motifs with a given number of nodes in the top level and the bottom level, 'levelsize' normalisation
#'  will return NAs. 'levelsize_plus1' avoids this by adding one to all counts.}
#'  \item{\strong{"levelsize_NAzero"}: same as 'levelsize' but replaces all NA values with 0. If a species does not occur in any motifs with a given number of nodes in the top level and the bottom level, 'levelsize' normalisation
#'  will return NAs. 'levelsize_NAzero' avoids this by replacing NAs with zero.}
#'  }
#'  }
#'  \item{\strong{Motif normalisation}
#'  \itemize{
#'  \item{\strong{"motif"}: divides the position measure for each node by the total number of times that node appears in any position within the same motif.
#'  For example, the relative frequencies of all position measures in motif 5 will sum to one, as will the relative frequency of all position counts in motif 10, and so on.}
#'  \item{\strong{"motif_plus1"}: same as 'motif' but adds one to all position measures. If a species does not occur in a particular motif, 'motif' normalisation
#'  will return NAs. 'motif_plus1' avoids this by adding one to all counts.}
#'  \item{\strong{"motif_NAzero"}: same as 'motif' but replaces all NA values with 0. If a species does not occur in a particular motif, 'levelsize' normalisation
#'  will return NAs. 'motif_NAzero' avoids this by replacing NAs with zero.}
#'  }
#'  }
#'  }
#'
#' @return
#' Returns a data frame with one column for each node position: 46 columns if \code{six_node} is FALSE, and 148 columns if \code{six_node} is TRUE.
#' Columns names are given as "npx" where x is the ID of the position as described in Simmons et al. (2017) (and originally in Appendix 1 of Baker et al. (2015)). \strong{To view the 'motif dictionary' showing
#' which node position a given ID corresponds to, enter \code{vignette("bmotif-dictionary")}.}
#'
#' For a network with A rows and P columns, by default (where \code{level} = "all") the data frame has A + P rows, one for each node. If \code{level} = "rows", the data frame will have A rows, one for each row node;
#' if \code{level} = "columns", it will have P rows, one for each column node.
#'
#' By default, the elements of this data frame will be the raw binary or weighted position measures (depending on which was requested). If \code{normalisation} is set to something other than "none", the elements will be
#' normalised position counts as described above.
#'
#' If \code{weights_method} is set to 'all', \code{node_positions} instead returns a list of length five, each containing a data.frame corresponding to
#' one of the five weighting methods described above.
#'
#' @export
#' @references
#' Baker, N., Kaartinen, R., Roslin, T., and Stouffer, D. B. (2015). Species’ roles in food webs show fidelity across a highly variable oak forest. Ecography, 38(2):130–139.
#'
#' Mora, B.B., Cirtwill, A.R. and Stouffer, D.B., 2018. pymfinder: a tool for the motif analysis of binary and quantitative complex networks. bioRxiv, 364703.
#'
#' Simmons, B. I., Sweering, M. J. M., Dicks, L. V., Sutherland, W. J. and Di Clemente, R. bmotif: a package for counting motifs in bipartite networks. bioRxiv. doi: 10.1101/302356
#' @examples
#' set.seed(123)
#' row <- 10
#' col <- 10
#'
#' # node positions in a binary network
#' m <- matrix(sample(0:1, row*col, replace=TRUE), row, col)
#' node_positions(M = m, six_node = TRUE, weights_method = "none", weights_combine = "none")
#'
#' # node positions in a weighted network
#' m[m>0] <- stats::runif(sum(m), 0, 100)
#' node_positions(M = m, six_node = FALSE, weights_method = "all", weights_combine = "sum")


node_positions <- function(M, six_node = FALSE, level = "all", weights_method, weights_combine = "none", normalisation = "none"){
  # check inputs
  if(inherits(M, "matrix") != TRUE){stop("'M' must be an object of class 'matrix'")} # make sure M is a matrix
  if(!all(apply(M, 1:2, is.numeric))){stop("Elements of 'M' must be numeric")} # make sure all elements of M are numbers
  if(!all(apply(M, 1:2, function(x) length(x) > 0))){stop("Elements of 'M' cannot have length 0")} # make sure no elements of M have 0 length e.g. numeric(0)
  if(!all(apply(M, 1:2, function(x) x >= 0))){stop("Elements of 'M' must be greater than or equal to zero")} # make sure all elements of M are >= zero
  if(inherits(level,"character") != TRUE){stop("'level' must be of class 'character'")} # make sure level is a character
  if(!level %in% c("rows","columns","all")){stop("'level' must equal 'rows', 'columns' or 'all'")} # make sure level equals 'rows', 'columns' or 'all'
  if(inherits(normalisation,"character") != TRUE){stop("'normalisation' must be of class 'character'")} # make sure 'normalisation' is a character
  if(!normalisation %in% c("none","sum","sizeclass", "sizeclass_plus1", "sizeclass_NAzero", "position","levelsize","levelsize_plus1","levelsize_NAzero","motif","motif_plus1","motif_NAzero")){stop("'normalisation' must equal 'none','sum','sizeclass', 'sizeclass_plus1', 'sizeclass_NAzero', 'position','levelsize','levelsize_plus1', 'levelsize_NAzero', 'motif','motif_plus1','motif_NAzero'")} # make sure normalisation equals 'none', 'sum' or 'sizeclass'
  if(any(duplicated(rownames(M))) | any(duplicated(colnames(M)))){stop("Input matrix must not have duplicate row or column names")}
  if(!weights_method %in% c('none', 'mean_motifweights', 'total_motifweights', 'mean_nodeweights', 'total_nodeweights', 'contribution', 'mora', 'all')) {
    stop("weights_method must be one of 'none', 'mean_motifweights', 'total_motifweights', 'mean_nodeweights', 'total_nodeweights', 'contribution', 'mora', 'all'.")
  }
  if(!weights_combine %in% c('none', 'mean', 'sum')) {
    stop("weights_combine must be one of 'none', 'sum', 'mean'.")
  }

  # Give warnings or errors if two arguments can't be used together
  if(weights_method != "none" & six_node == TRUE){
    stop("Sorry, weighted methods are not available for six node motifs. Set six_node = FALSE if you want to use weighted position measures.")
  }

  if(weights_method == "none" & weights_combine != "none"){
    stop("Cannot set a weights_combine method when weights_method = 'none'.\nIf you wish to return binary position counts, set weights_combine to 'none'.\nIf you want weighted position measures, set weights_method to one of 'mean_motifweights', 'total_motifweights', 'mean_nodeweights', 'total_nodeweights', 'contribution', 'mora', 'all'.")
  }

  if(weights_method != "none" & weights_combine == "none"){
    stop("Need to set a weights_combine method. Cannot have weights_method != 'none' and weights_combine = 'none'. Set weights_combine to 'sum' or 'mean' if you want to use weighted position measures.")
  }

  if(weights_combine == "mean" & normalisation != "none"){
    stop("Taking the mean is already a form of normalisation. Set normalisation = 'none' if you want to use weights_combine = 'mean'. Alternatively, set weights_combine to something other than 'mean' if you want to use one of the normalisation methods.")
  }

  # clean matrix
  W <- M # save copy of weighted matrix
  M[M > 0] <- 1 # now M is the binarised version
  if(is.null(rownames(M))){ # if M has no row names, give it some
    rownames(M) <- paste0("r", 1:nrow(M))
  }
  if(is.null(colnames(M))){ # if M has no column names, give it some
    colnames(M) <- paste0("c", 1:ncol(M))
  }
  rn <- rownames(M) # record row names
  cn <- colnames(M) # record column names
  dimnames(M) <- NULL # strip row and column names


  # ------------ BINARY COUNT ---------------------------------------

  # ACTUAL COUNT, NEEDED FOR weights_combine = MEAN OR weights_method= NONE

  if (weights_method == "none" | weights_method == "all" | weights_method == "contribution" | weights_combine == "mean") {
    # calculate inputs
    Z <- nrow(M) # number of row species
    P <- ncol(M) # number of column species

    jZ <- matrix(rep(1, Z)) # Z x 1 vector of ones
    jP <- matrix(rep(1, P)) # P x 1 vector of ones
    J <- matrix(rep(1, Z * P), nrow = Z, ncol = P) # Z x P matrix of ones
    JP <- matrix(rep(1, P * P), nrow = P, ncol = P) # P x P matrix of ones
    JZ <- matrix(rep(1, Z * Z), nrow = Z, ncol = Z) # Z x Z matrix of ones

    if(six_node == TRUE){
      JP3 <- array(rep(1, P * P * P), c(P, P, P)) # P x P x P array of ones
      JZ3 <- array(rep(1, Z * Z * Z), c(Z, Z, Z)) # Z x Z x Z array of ones
      KP3 <- JP3 # P x P x P matrix, 0 if two indices are equal, 1 otherwise
      for (i in 1 : P){
        for (j in 1 : P){
          KP3[i,j,j] <- 0
          KP3[j,i,j] <- 0
          KP3[j,j,i] <- 0
        }
      }
      KZ3 <- JZ3 # Z x Z x Z matrix, 0 if two indices are equal, 1 otherwise
      for (i in 1 : Z){
        for (j in 1 : Z){
          KZ3[i,j,j] <- 0
          KZ3[j,i,j] <- 0
          KZ3[j,j,i] <- 0
        }
      }
    }

    MT <- t(M) # transpose of M
    N <- J - M # complement of M
    NT <- t(N) # transpose of the complement of M

    dZ <- M %*% jP # degrees of the row species
    dP <- MT %*% jZ # degrees of the column species

    Z <- M %*% MT
    Y <- M %*% NT
    X <- N %*% MT

    P <- MT %*% M
    Q <- MT %*% N
    R <- NT %*% M

    if(six_node == TRUE){
      AZ <- maketensor(M, M)
      BZ <- maketensor(M, N)
      CZ <- maketensor(N, M)
      DZ <- maketensor(N, N)

      AP <- maketensor(MT, MT)
      BP <- maketensor(MT, NT)
      CP <- maketensor(NT, MT)
      DP <- maketensor(NT, NT)

      MTA <- multtensor(MT, AZ)
      MTB <- multtensor(MT, BZ)
      MTC <- multtensor(MT, CZ)
      MTD <- multtensor(MT, DZ)

      MA <- multtensor(M, AP)
      MB <- multtensor(M, BP)
      MC <- multtensor(M, CP)
      MD <- multtensor(M, DP)

      NTA <- multtensor(NT, AZ)
      NTB <- multtensor(NT, BZ)
      NTC <- multtensor(NT, CZ)

      Na <- multtensor(N, AP) # because NA already means something
      NB <- multtensor(N, BP)
      NC <- multtensor(N, CP)
    }

    # create results containers
    if(six_node == FALSE){
      pos_row <- matrix(0, ncol = 46, nrow = nrow(M), dimnames = list(rn,paste0("np",1:46)))
      pos_col <- matrix(0, ncol = 46, nrow = ncol(M), dimnames = list(cn,paste0("np",1:46)))
    } else {
      pos_row <- matrix(0, ncol = 148, nrow = nrow(M), dimnames = list(rn,paste0("np",1:148)))
      pos_col <- matrix(0, ncol = 148, nrow = ncol(M), dimnames = list(cn,paste0("np",1:148)))
    }

    # count positions
    if(six_node == FALSE){
      for(i in 1:46){
        rc <- rowcolumn(i)
        f <- countposition(M = M, p = i, jZ = jZ, jP = jP, JP = JP, JZ = JZ, MT = MT, N = N, NT = NT, dZ = dZ, dP = dP, Z = Z, Y = Y, X = X, P = P, Q = Q, R = R)
        if(rc == "row"){
          pos_row[,i] <- f
        } else {
          pos_col[,i] <- f
        }
      }
    } else {
      for(i in 1:148){
        rc <- rowcolumn(i)
        f <- countposition(M = M, p = i, jZ = jZ, jP = jP, JP = JP, JZ = JZ, JP3 = JP3, JZ3 = JZ3, KP3 = KP3, KZ3 = KZ3, MT = MT, N = N, NT = NT, dZ = dZ, dP = dP, Z = Z, Y = Y, X = X, P = P, Q = Q, R = R, MTA = MTA, MTB = MTB, MTC = MTC, MTD = MTD, MA = MA, MB = MB, MC = MC, MD = MD, NTA = NTA, NTB = NTB, NTC = NTC, Na = Na, NB = NB, NC = NC)
        if(rc == "row"){
          pos_row[,i] <- f
        } else {
          pos_col[,i] <- f
        }
      }
    }

  } # end of calculation of binary count



  # -------------------------- weights_method------------------------------

  # ------------------------- SOME AUXILIARY STUFF ------------------------

  NZ <- nrow(W)
  NP <- ncol(W)

  # save number of edges and vertices for each motif
  motif_n_edges <- c(1,2,2,3,3,4,3,4,4,4,5,6,4,4,5,6,4)
  motif_n_vertices <- c(2,3,3,4,4,4,4,rep(5, 10))
  # save degree of each position
  pos_degrees <- c(c(1,1), c(1,2,2,1), c(3,1,1,2,1,2,2,2,1,3), c(4,1,1,3,1,2,2,1,2,2,3,1,2), c(3,2,1,2,1,3,1,2,2,1,2,2,3,2,3,1,4))

  # create a matrix with number of edges for each motif
  nedges <- matrix(NA, NZ + NP, 46)
  for (i in 1:17) {
    positions <- motif_info(i, node = TRUE, link = FALSE)
    nedges[,positions] <- motif_n_edges[i]
  }

  degrees <- matrix(NA, nrow(W) + ncol(W), 46)
  for (i in 1:46) {
    degrees[,i] <- pos_degrees[i]
  }

  # ------------------- WEIGHTS = "NONE" ------------------------------------

  if (weights_method== "none") {
    # output
    if(level == "all"){
      out <- rbind(pos_row, pos_col)
      out <- normalise_node_positions(pc = out, type = normalisation, six_node = six_node)
      out <- as.data.frame(out)
      return(out)
    } else if(level == "rows"){
      pos_row <- normalise_node_positions(pc = pos_row, type = normalisation, six_node = six_node)
      pos_row <- as.data.frame(pos_row)
      return(pos_row)
    } else if(level == "columns"){
      pos_col <- normalise_node_positions(pc = pos_col, type = normalisation, six_node = six_node)
      pos_col <- as.data.frame(pos_col)
      return(pos_col)
    }
  }

  # ------------------ WEIGHTS = ALL --------------------------

  if (weights_method== "all") {

    out <- list()
    # mmw gives back the sum of the mean motif weights whenever species x is in position p
    mmw <- matrix(NA, nrow = NZ + NP, ncol = 46) # mean motif weights
    colnames(mmw) <- paste('np', 1:46, sep = '')
    rownames(mmw) <- weighted_node_positions_output_row_names(x = W, NZ = NZ, NP = NP) #c(paste('r', 1:NZ, sep = ''), paste('c', 1:NP, sep = ''))
    # mnw gives back the sum of the mean weights that species x is in when in pos p (mean node weights)
    mnw <- matrix(NA, nrow = NZ + NP, ncol = 46) # mean node weights
    colnames(mnw) <- paste('np', 1:46, sep = '')
    rownames(mnw) <- weighted_node_positions_output_row_names(x = W, NZ = NZ, NP = NP) #c(paste('r', 1:NZ, sep = ''), paste('c', 1:NP, sep = ''))
    # con gives back the sum of contributions of species x to the motif weights
    con <- matrix(NA, nrow = NZ + NP, ncol = 46) # mean contribution
    colnames(con) <- paste('np', 1:46, sep = '')
    rownames(con) <- weighted_node_positions_output_row_names(x = W, NZ = NZ, NP = NP) #c(paste('r', 1:NZ, sep = ''), paste('c', 1:NP, sep = ''))
    # py gives back the summed pymfinder measure
    py <- matrix(NA, nrow = NZ + NP, ncol = 46) # pymfinder contribution
    colnames(py) <- paste('np', 1:46, sep = '')
    rownames(py) <- weighted_node_positions_output_row_names(x = W, NZ = NZ, NP = NP) #c(paste('r', 1:NZ, sep = ''), paste('c', 1:NP, sep = ''))

    # I also need the node_position count
    np_count <- rbind(pos_row, pos_col)

    # position 1 and 2
    # motif weights
    mmw[1:NZ,1] <- 0
    mmw[(NZ + 1):(NZ + NP),1] <- apply(W, 2, sum)
    mmw[1:NZ,2] <- apply(W, 1, sum)
    mmw[(NZ + 1):(NZ + NP), 2] <- 0

    # node weights (the same)
    mnw[1:NZ,1] <- 0
    mnw[(NZ + 1):(NZ + NP),1] <- apply(W, 2, sum)
    mnw[1:NZ,2] <- apply(W, 1, sum)
    mnw[(NZ + 1):(NZ + NP), 2] <- 0

    # contribution (always 1, so it is just the count)
    con[,1:2] <- np_count[,1:2]

    # pymfinder measure
    py[1:NZ, 1] <- 0
    py[(NZ + 1):(NZ + NP),1] <- 0.5 * apply(W, 2, sum)
    py[1:NZ, 2] <- 0.5 * apply(W, 1, sum)
    py[(NZ + 1):(NZ + NP),2] <- 0

    # for pos 3-46, calculation is done in C++
    pos <- 3
    for (i in 2:17) {
      funname <- paste('np_m', i, '_c', sep ='')
      lst <- get(funname)(NZ, NP, W)
      len <- length(lst)
      for (j in 1:len) {
        rc <- rowcolumn(pos)
        # print(paste('top of loop j', j, 'pos ', pos))
        # The 1st, 5th, 9th etc entry of the list is for mean weights
        # The 2nd, 6th, 8th etc. entry of the list is for node weights
        # the 3rd, 7th etc entry of the list is for contribution
        # the 4th, 8th entry of the list is for pymfinder measure
        # after this, we move on to a new position

        if ((j %% 4 == 1) & (rc == "row")) {
          # print(paste('j ', j, 'pos', pos, 'rc', rc, 'list '))
          # print(lst[[j]])
          mmw[1:NZ,pos] <- lst[[j]]
          mmw[(NZ + 1):(NZ + NP),pos] <- 0
        } else if ((j %% 4 == 1) & (rc == "column")) {
          # print(paste('j ', j, 'pos', pos, 'rc', rc, 'list '))
          # print(lst[[j]])
          mmw[1:NZ,pos] <- 0
          mmw[(NZ + 1):(NZ + NP),pos] <- lst[[j]]
        } else if ((j %% 4 == 2) & (rc == "row")) {
          mnw[1:NZ,pos] <- lst[[j]]
          mnw[(NZ + 1):(NZ + NP),pos] <- 0
        } else if ((j %% 4 == 2) & (rc == "column")) {
          mnw[1:NZ,pos] <- 0
          mnw[(NZ + 1):(NZ + NP),pos] <- lst[[j]]
        } else if ((j %% 4 == 3) & (rc == "row")) {
          con[1:NZ,pos] <- lst[[j]]
          con[(NZ + 1):(NZ + NP),pos] <- 0
        } else if ((j %% 4 == 3) & (rc == "column")) {
          con[1:NZ,pos] <- 0
          con[(NZ + 1):(NZ + NP),pos] <- lst[[j]]
        }  else if ((j %% 4 == 0) & (rc == "row")) {
          py[1:NZ,pos] <- lst[[j]]
          py[(NZ + 1):(NZ + NP),pos] <- 0
          pos <- pos + 1
        } else if ((j %% 4 == 0) & (rc == "column")) {
          py[1:NZ,pos] <- 0
          py[(NZ + 1):(NZ + NP),pos] <- lst[[j]]
          pos <- pos + 1
        }

      } # close for j in 1:len
    } # close for i in 2:17

    # now compute total motif weight and total node weight
    # total motif weight: need to multiply each position-column by the number of edges in the motif
    tmw <- mmw * nedges
    # total node weight: need to multiply each position-column by the degree
    tnw <- mnw * degrees

    # so far, we have the sums, now we act according to the weights_combine argument
    if (weights_combine == "mean") {
      # divide by position count
      mmw <- mmw / np_count
      mnw <- mnw / np_count
      tmw <- tmw / np_count
      tnw <- tnw / np_count
      con <- con / np_count
      py <- py / np_count

      # replace NaNs by NAs
      mmw <- replace(mmw, which(is.nan(mmw)), NA)
      mnw <- replace(mnw, which(is.nan(mnw)), NA)
      tmw <- replace(tmw, which(is.nan(tmw)), NA)
      tnw <- replace(tnw, which(is.nan(tnw)), NA)
      con <- replace(con, which(is.nan(con)), NA)
      py <- replace(py, which(is.na(py)), NA)

      # no need to normalise here
      # now set output list
      out[[1]] <- mmw
      out[[2]] <- tmw
      out[[3]] <- mnw
      out[[4]] <- tnw
      out[[5]] <- con
      out[[6]] <- py
      names(out) <- c("mean_motifweights", "total_motifweights", "mean_nodeweights", "total_nodeweights", "contribution", "mora")

      # depending on level, delete unused rows
      if(level == "rows"){
        out <- lapply(out, function(x) {x[1:NZ, ]})
      } else if(level == "columns"){
        out <- lapply(out, function(x) {x[(NZ + 1):(NZ + NP), ]})
      }
      out <- lapply(out, function(x) {as.data.frame(x)})
      return(out)
    } # end weights_combine = "mean"

    if (weights_combine == "sum") {
      # mw, nw and con already contain the correct results
      # normalisation is possible now
      if(normalisation != "none"){
        mmw <- normalise_node_positions(pc = mmw, type = normalisation, six_node = six_node)
        tmw <- normalise_node_positions(pc = tmw, type = normalisation, six_node = six_node)
        mnw <- normalise_node_positions(pc = mnw, type = normalisation, six_node = six_node)
        tnw <- normalise_node_positions(pc = tnw, type = normalisation, six_node = six_node)
        con <- normalise_node_positions(pc = con, type = normalisation, six_node = six_node)
        py <- normalise_node_positions(pc = py, type = normalisation, six_node = six_node)
      }

      # now set output list
      out[[1]] <- mmw
      out[[2]] <- tmw
      out[[3]] <- mnw
      out[[4]] <- tnw
      out[[5]] <- con
      out[[6]] <- py
      names(out) <- c("mean_motifweights", "total_motifweights", "mean_nodeweights", "total_nodeweights", "contribution", "mora")

      # depending on level, delete unused rows
      if(level == "rows"){
        out <- lapply(out, function(x) {x[1:NZ, ]})
      } else if(level == "columns"){
        out <- lapply(out, function(x) {x[(NZ + 1):(NZ + NP), ]})
      }
      out <- lapply(out, function(x) {as.data.frame(x)})
      return(out)
    } # end weights_combine = "sum"
  } # end weights_method= "all"

  # --------------- weights_method = MOTIFWEIGHTS -----------------------------

  if (weights_method== "mean_motifweights" | weights_method== "total_motifweights") {
    # motifweights gives back the sum of the mean motif weights whenever species x is in position p

    mmw <- matrix(NA, nrow = NZ + NP, ncol = 46)
    colnames(mmw) <- paste('np', 1:46, sep = '')
    rownames(mmw) <- weighted_node_positions_output_row_names(x = W, NZ = NZ, NP = NP) #c(paste('r', 1:NZ, sep = ''), paste('c', 1:NP, sep = ''))

    # position 1

    mmw[1:NZ,1] <- 0
    mmw[(NZ + 1):(NZ + NP),1] <- apply(W, 2, sum)
    mmw[1:NZ,2] <- apply(W, 1, sum)
    mmw[(NZ + 1):(NZ + NP), 2] <- 0

    # for pos 3-46, calculation is done in C++
    pos <- 3
    for (i in 2:17) {
      funname <- paste('np_m', i, '_mw', sep ='')
      lst <- get(funname)(NZ, NP, W)
      len <- length(lst)
      for (j in 1:len) {
        if (rowcolumn(pos) == "row") {
          mmw[1:NZ,pos] <- lst[[j]]
          mmw[(NZ + 1):(NZ + NP),pos] <- 0
        } else {
          mmw[1:NZ,pos] <- 0
          mmw[(NZ + 1):(NZ + NP),pos] <- lst[[j]]
        }
        pos <- pos + 1
      }
    }

    # at this point we have all information from the C++ file we need, now adjust

    if (weights_method== "mean_motifweights") {
      if (weights_combine == "mean") {
        # I also need the node_position count
        np_count <- rbind(pos_row, pos_col)

        mmw <- mmw / np_count
        mmw <- replace(mmw, which(is.nan(mmw)), NA)

        # no normalisation necessary

        # convert to data frame
        mmw <- as.data.frame(mmw)

        # depending on level, delete unused rows
        if (level == "all") {
          return(mmw)
        } else if(level == "rows"){
          return(mmw[1:NZ,])
        } else if(level == "columns"){
          return(mmw[(NZ + 1):(NZ + NP),])
        }
      } # end weights_combine == "mean"

      if (weights_combine == "sum") {

        # can have normalisation
        if(normalisation != "none"){
          mmw <- normalise_node_positions(pc = mmw, type = normalisation, six_node = six_node)
        }

        # convert to data frame
        mmw <- as.data.frame(mmw)

        # depending on level, delete unused rows
        if (level == "all") {
          return(mmw)
        } else if(level == "rows"){
          return(mmw[1:NZ,])
        } else if(level == "columns"){
          return(mmw[(NZ + 1):(NZ + NP),])
        }
      } # end weights_combine == "sum"
    } # end weights_method= "mean_motifweights"


    if (weights_method== "total_motifweights") {

      tmw <- mmw * nedges # multiply mean motif weights by number of edges

      if (weights_combine == "mean") {
        # I also need the node_position count
        np_count <- rbind(pos_row, pos_col)

        tmw <- tmw / np_count
        tmw <- replace(tmw, which(is.nan(tmw)), NA)

        # no normalisation necessary

        # convert to data frame
        tmw <- as.data.frame(tmw)

        # depending on level, delete unused rows
        if (level == "all") {
          return(tmw)
        } else if(level == "rows"){
          return(tmw[1:NZ,])
        } else if(level == "columns"){
          return(tmw[(NZ + 1):(NZ + NP),])
        }
      } # end weights_combine == "mean"

      if (weights_combine == "sum") {

        # can have normalisation
        if(normalisation != "none"){
          tmw <- normalise_node_positions(pc = tmw, type = normalisation, six_node = six_node)
        }

        # convert to data frame
        tmw <- as.data.frame(tmw)

        # depending on level, delete unused rows
        if (level == "all") {
          return(tmw)
        } else if(level == "rows"){
          return(tmw[1:NZ,])
        } else if(level == "columns"){
          return(tmw[(NZ + 1):(NZ + NP),])
        }
      } # end weights_combine == "sum"
    } # end weights_method= "total_motifweights"


  } # end weights_method= "mean_motifweights" or "total_motifweights"

  # -------- weights_method= NODEWEIGHTS ------------------------

  if (weights_method== "mean_nodeweights" | weights_method== "total_nodeweights") {
    # nodeweights gives back the sum of the mean of the nodeweights (i.e. the weights of the links that species x is in whenever in pos p)

    mnw <- matrix(NA, nrow = NZ + NP, ncol = 46)
    colnames(mnw) <- paste('np', 1:46, sep = '')
    rownames(mnw) <- weighted_node_positions_output_row_names(x = W, NZ = NZ, NP = NP) #c(paste('r', 1:NZ, sep = ''), paste('c', 1:NP, sep = ''))

    # position 1
    mnw[1:NZ,1] <- 0
    mnw[(NZ + 1):(NZ + NP),1] <- apply(W, 2, sum)
    mnw[1:NZ,2] <- apply(W, 1, sum)
    mnw[(NZ + 1):(NZ + NP), 2] <- 0

    # for pos 3-46, calculation is done in C++
    pos <- 3
    for (i in 2:17) {
      funname <- paste('np_m', i, '_nw', sep ='')
      lst <- get(funname)(NZ, NP, W)
      len <- length(lst)
      for (j in 1:len) {
        if (rowcolumn(pos) == "row") {
          mnw[1:NZ,pos] <- lst[[j]]
          mnw[(NZ + 1):(NZ + NP),pos] <- 0
        } else {
          mnw[1:NZ,pos] <- 0
          mnw[(NZ + 1):(NZ + NP),pos] <- lst[[j]]
        }
        pos <- pos + 1
      }
    }

    if (weights_method== "mean_nodeweights") {
      if (weights_combine == "mean") {
        # I also need the node_position count
        np_count <- rbind(pos_row, pos_col)

        mnw <- mnw / np_count
        mnw <- replace(mnw, which(is.nan(mnw)), NA)

        # convert to data frame
        mnw <- as.data.frame(mnw)

        # again, normalisation does not make sense
        # depending on level, delete unused rows
        if (level == "all") {
          return(mnw)
        } else if(level == "rows"){
          return(mnw[1:NZ,])
        } else if(level == "columns"){
          return(mnw[(NZ + 1):(NZ + NP),])
        }
      } # end weights_combine == "mean"

      if (weights_combine == "sum") {
        # can have normalisation
        if(normalisation != "none"){
          mnw <- normalise_node_positions(pc = mnw, type = normalisation, six_node = six_node)
        }

        # convert to data frame
        mnw <- as.data.frame(mnw)

        # depending on level, delete unused rows
        if (level == "all") {
          return(mnw)
        } else if(level == "rows"){
          return(mnw[1:NZ,])
        } else if(level == "columns"){
          return(mnw[(NZ + 1):(NZ + NP),])
        }
      } # end weights_combine == "sum"
    } # end weights_method== "mean_nodeweights"

    if (weights_method== "total_nodeweights") {

      tnw <- mnw * degrees

      if (weights_combine == "mean") {
        # I also need the node_position count
        np_count <- rbind(pos_row, pos_col)

        tnw <- tnw / np_count
        tnw <- replace(tnw, which(is.nan(tnw)), NA)

        # convert to data frame
        tnw <- as.data.frame(tnw)

        # again, normalisation does not make sense
        # depending on level, delete unused rows
        if (level == "all") {
          return(tnw)
        } else if(level == "rows"){
          return(tnw[1:NZ,])
        } else if(level == "columns"){
          return(tnw[(NZ + 1):(NZ + NP),])
        }
      } # end weights_combine == "mean"

      if (weights_combine == "sum") {
        # can have normalisation
        if(normalisation != "none"){
          tnw <- normalise_node_positions(pc = tnw, type = normalisation, six_node = six_node)
        }

        # convert to data frame
        tnw <- as.data.frame(tnw)

        # depending on level, delete unused rows
        if (level == "all") {
          return(tnw)
        } else if(level == "rows"){
          return(tnw[1:NZ,])
        } else if(level == "columns"){
          return(tnw[(NZ + 1):(NZ + NP),])
        }
      } # end weights_combine == "sum"
    } # end weights_method== "mean_nodeweights"


  }

  # ------------- weights_method= CONTRIBUTION -------------------

  if (weights_method== "contribution") {

    con <- matrix(NA, nrow = NZ + NP, ncol = 46)
    colnames(con) <- paste('np', 1:46, sep = '')
    rownames(con) <- weighted_node_positions_output_row_names(x = W, NZ = NZ, NP = NP) #c(paste('r', 1:NZ, sep = ''), paste('c', 1:NP, sep = ''))

    # position 1
    np_count <- rbind(pos_row, pos_col)
    con[,1:2] <- np_count[,1:2]

    # for pos 3-46, calculation is done in C++
    pos <- 3
    for (i in 2:17) {
      funname <- paste('np_m', i, '_con', sep ='')
      lst <- get(funname)(NZ, NP, W)
      len <- length(lst)
      for (j in 1:len) {
        if (rowcolumn(pos) == "row") {
          con[1:NZ,pos] <- lst[[j]]
          con[(NZ + 1):(NZ + NP),pos] <- 0
        } else {
          con[1:NZ,pos] <- 0
          con[(NZ + 1):(NZ + NP),pos] <- lst[[j]]
        }
        pos <- pos + 1
      }
    }

    if (weights_combine == "mean") {
      con <- con / np_count
      con <- replace(con, which(is.nan(con)), NA)

      # convert to data frame
      con <- as.data.frame(con)

      # again, normalisation does not make sense
      # depending on level, delete unused rows
      if (level == "all") {
        return(con)
      } else if(level == "rows"){
        return(con[1:NZ,])
      } else if(level == "columns"){
        return(con[(NZ + 1):(NZ + NP),])
      }

    } # end weights_combine == "mean"

    if (weights_combine == "sum") {

      # can have normalisation
      if(normalisation != "none"){
        con <- normalise_node_positions(pc = con, type = normalisation, six_node = six_node)
      }

      con <- as.data.frame(con)

      if (level == "all") {
        return(con)
      } else if(level == "rows"){
        return(con[1:NZ,])
      } else if(level == "columns"){
        return(con[(NZ + 1):(NZ + NP),])
      }

    } # end weights_combine == "sum"
  } # end weights_method == "contribution"

  # -------------- weights_method = STOUFFER --------------------------------------

  if (weights_method == "mora") {

    py <- matrix(NA, nrow = NZ + NP, ncol = 46)
    colnames(py) <- paste('np', 1:46, sep = '')
    rownames(py) <- weighted_node_positions_output_row_names(x = W, NZ = NZ, NP = NP) #c(paste('r', 1:NZ, sep = ''), paste('c', 1:NP, sep = ''))

    # position 1

    py[1:NZ, 1] <- 0
    py[(NZ + 1):(NZ + NP),1] <- 0.5 * apply(W, 2, sum)
    py[1:NZ, 2] <- 0.5 * apply(W, 1, sum)
    py[(NZ + 1):(NZ + NP),2] <- 0

    # for pos 3-46, calculation is done in C++
    pos <- 3
    for (i in 2:17) {
      funname <- paste('np_m', i, '_py', sep ='')
      lst <- get(funname)(NZ, NP, W)
      len <- length(lst)
      for (j in 1:len) {
        if (rowcolumn(pos) == "row") {
          py[1:NZ,pos] <- lst[[j]]
          py[(NZ + 1):(NZ + NP),pos] <- 0
        } else {
          py[1:NZ,pos] <- 0
          py[(NZ + 1):(NZ + NP),pos] <- lst[[j]]
        }
        pos <- pos + 1
      }
    }

    if (weights_combine == "mean") {

      np_count <- rbind(pos_row, pos_col)
      py <- py / np_count
      py <- replace(py, which(is.nan(py)), NA)

      # convert to data frame
      py <- as.data.frame(py)

      # again, normalisation does not make sense
      # depending on level, delete unused rows
      if (level == "all") {
        return(py)
      } else if(level == "rows"){
        return(py[1:NZ,])
      } else if(level == "columns"){
        return(py[(NZ + 1):(NZ + NP),])
      }

    } # end weights_combine == "mean"

    if (weights_combine == "sum") {

      # can have normalisation
      if(normalisation != "none"){
        py <- normalise_node_positions(pc = py, type = normalisation, six_node = six_node)
      }

      # convert to data frame
      py <- as.data.frame(py)

      if (level == "all") {
        return(py)
      } else if(level == "rows"){
        return(py[1:NZ,])
      } else if(level == "columns"){
        return(py[(NZ + 1):(NZ + NP),])
      }

    } # end weights_combine == "sum"
  }


  # -------------------------- END WEIGHTS ----------------------------


}
