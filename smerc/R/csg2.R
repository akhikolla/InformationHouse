#' Construct connected subgraphs
#'
#' \code{csg2}, \code{lcsg2}, and \code{scsg2} construct
#' connected subgraphs. These functions are not intended
#' for users.
#' \code{nn} contains a list of nearest
#' neighbors for each region. \code{idx} is a
#' vector of possible vertices being considered as a
#' subgraph. \code{w} is a connectivity matrix relating the
#' N vertices. \code{w[i,j] = 1} if vertices i and j are
#' connected, i.e., if they share an edge. The dimensions of
#' \code{w} are \eqn{N times k}, where \code{k =
#' length(idx)}. While the rows of \code{w} contain
#' adjacency information for all N vertices, only the
#' \code{idx} columns of the complete adjacency matrix are
#' used in \code{w}.  See Details for discussion of
#' \code{scsg}.
#'
#' \code{scsg2} performs a sequence of \code{lcsg2} calls.
#' Starting with \code{lcz == list(idx[1])}, \code{scsg}
#' keeps iteratively building more connected subsgraphs by
#' perfoming something like:  lcz1 = list(idx[1]).  lcz2 =
#' lcsg2(lcz1, ...). lcz3 = lcsg2(lcz2, ...).  This is
#' done until there are no more connected subgraphs among
#' the elements of \code{idx}.
#'
#' @param nn A list of the nearest neighbors for each vertex (region).
#' @param w A binary adjacency matrix indicating connected neighbors.
#' @param idx A vector of vertices for which to construct the set of connected subgraphs.
#' @param nlevel The maximum size of each subgraph.
#' @param verbose A logical value indicating whether descriptive messages should be provided.  Default is
#'   \code{FALSE}.  If \code{TRUE}, this can be useful for
#'   diagnosing where the sequences of connected subgraphs
#'   is slowing down/having problems.
#' @param logical A logical value indicating whether a list of logical vectors should be returned. The default is \code{FALSE},
#' indicating that the \code{scsg} function should return a list of vectors with
#' each vector containing the vertex indices included in each subgraph.
#' @param cz A logical vector representing the current subgraph.
#' @param cw A binary adjacency matrix for the neighbors of the current vertex.
#' @param cnn The indices of the neighbors of the current vertex.
#' @param lcz A list of current zones (in the form of logical vectors).
#' @return A list with all possible connected subgraphs based on the
#' user-provided parameters.
#' @export
#'
#' @examples
#' data(nydf)
#' data(nyw)
#' # determine 50 nn of region 1 for NY data
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' nn3 = knn(coords, longlat = TRUE, k = 3)
#' z1 = scsg2(nn3, nyw)
#' z2 = flex.zones(coords, nyw, k = 3, longlat = TRUE)
#' all.equal(z1, z2)
#' @export
#' @rdname csg
csg2 = function(cz, cnn, cw) {

  # determine neighbors of current set
  nb = which(((cz + colSums(cw[cz, , drop = FALSE])) >= 1) - cz == 1)

  # if there are some neighbors
  if (length(nb) > 0) {
    # return new set of candidate zones
    lapply(nb, function(i, tcz) {
      tcz[i] = TRUE
      tcz
    }, tcz = cz)
  } else {# return NULL if there are no neighbors
    return(logical(0))
  }
}

#' @export
#' @rdname csg
lcsg2 = function(lcz, cnn, cw) {
  x = unlist(lapply(lcz, function(cz) {
    csg2(cz, cnn, cw)
  }), recursive = FALSE)

  if (is.null(x)) {
    return(NULL)
  } else {
    x[smerc::distinct(x, N = length(cnn))]
  }
}

#' @export
#' @rdname csg
scsg2 = function(nn, w, idx = seq_along(nn),
                 nlevel = NULL, verbose = FALSE,
                 logical = FALSE) {
  if (is.null(nlevel)) {
    nlevel = max(sapply(nn, length))
  }
  # list to store all results
  lcz_all = vector("list", length(idx))

  # set stopping conditions
  # stop if nidx == 1, otherwise, consider expanding subgraph
  for (i in idx) {
    # create temporary list to store results for region i
    temp = vector("list", nlevel)
    temp[[1]] = list(logical(length(nn[[i]])))
    temp[[1]][[1]][1] = TRUE
    clevel =  1
    # while we haven't reached the max level
    while (clevel < nlevel) {
      # increment current level
      clevel = clevel + 1
      # print message, if necessary
      if (verbose) {
        txt = paste("region ", i, ": ", clevel, "/", nlevel,
                    " at ", Sys.time(), ".", sep = "")
        message(txt)
      }
      # determine next set of candidate zones from
      # previous set of candidate zones
      temp[[clevel]] = lcsg2(lcz = temp[[clevel - 1]],
                             cnn = nn[[i]],
                             cw = w[nn[[i]], nn[[i]], drop = FALSE])
      # if there are no more connected neighbors for
      # any of the previous level of zones, end algorithm
      if (length(temp[[clevel]]) < 1) {
        clevel = nlevel
      }
    }

    if (!logical) {
    lcz_all[[i]] = lapply(unlist(temp, recursive = FALSE),
                          function(x) {
                            z = logical(nrow(w))
                            z[nn[[i]][x]] = TRUE
                            return(z)
                          })
    } else {
      lcz_all[[i]] = unlist(temp, recursive = FALSE)
    }
  }

  if (!logical) {
    lcz_all = unlist(lcz_all, recursive = FALSE)
    lcz_all = lapply(lcz_all, which)
    return(lcz_all[smerc::distinct(lcz_all, N = nrow(w))])
  } else {
    return(lcz_all)
  }
}
