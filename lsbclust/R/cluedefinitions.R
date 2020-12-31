#' S3 export
#' 
#' These export
#' into the framework set out in package \pkg{clue}.
#' 
#' @rdname lsbclusttoclue
#' @aliases cl_class_ids.int.lsbclust is.cl_partition.int.lsbclust is.cl_hard_partition.int.lsbclust
#' @param x An object of class \code{int.lsclust}
#' @export
cl_class_ids.int.lsbclust <- function(x) as.cl_class_ids(x$cluster)

#' @rdname lsbclusttoclue
#' @export
is.cl_partition.int.lsbclust <- function(x) TRUE

#' @rdname lsbclusttoclue
#' @export
is.cl_hard_partition.int.lsbclust <- function(x) TRUE

#' @rdname lsbclusttoclue
#' @export
cl_class_ids.lsbclust_sim_part <- function(x) as.cl_class_ids(x$cluster)

#' @rdname lsbclusttoclue
#' @export
is.cl_partition.lsbclust_sim_part <- function(x) TRUE

#' @rdname lsbclusttoclue
#' @export
is.cl_hard_partition.lsbclust_sim_part <- function(x) TRUE

#' @rdname lsbclusttoclue
#' @export
cl_class_ids.T3Clusf <- function(x) as.cl_class_ids(x$cluster)

#' @rdname lsbclusttoclue
#' @export
is.cl_partition.T3Clusf <- function(x) TRUE

#' @rdname lsbclusttoclue
#' @export
is.cl_hard_partition.T3Clusf <- function(x) TRUE

#' @rdname lsbclusttoclue
#' @export
cl_class_ids.akmeans <- function(x) as.cl_class_ids(x$cluster)

#' @rdname lsbclusttoclue
#' @export
is.cl_partition.akmeans <- function(x) TRUE

#' @rdname lsbclusttoclue
#' @export
is.cl_hard_partition.akmeans <- function(x) TRUE