
redundant <- function(cond, x = NULL, simplify = TRUE){
  cond <- noblanks(cond)
  .redund(x, cond, simplify)
}

# ==== .redund() ====
# Switch order of first 2 args to provide dispatching on x
# Generic function
.redund <- function(x, cond, ...) UseMethod(".redund")

# ==== Method for class 'cti' ====
# identifies the asf that are redundant within some csf
#   x         cti
#   cond      character vector with the csf
#   simplify  output matrix instead of list if all csf have the same number of asf
# value: A list of logical vectors (lengths corresponding to the number of asf), 
#        or a matrix if simplify=TRUE and all csf have the same number of asf
.redund.cti <- function(x, cond, simplify = TRUE, full = FALSE, names = TRUE,
                        qc_full = qcond_csf(cond, sc, flat = TRUE)){
  if (!full) x <- full.ct(x)
  sc <- x$scores
  ok <- .qcondType(cond, colnames(sc), x$type, 
                   stdComplex.multiple.only = FALSE) %in% c("stdAtomic", "stdComplex")
  if (any(!ok)) 
    stop("Invalid input to redundant:\n", paste0("  ", cond[!ok], collapse = "\n"),
         call. = FALSE)
  if (!length(cond)) return(logical(0))
  asfs <- extract_asf(cond)
  uasfs <- unique(unlist(asfs))
  hmatches <- happly(asfs, match, table = uasfs)
  qc <- qcond_asf(uasfs, sc, force.bool = TRUE)
  mode(qc) <- "logical"
  out <- lapply(hmatches, function(x) C_redund(qc[, x, drop = FALSE]))
  names(out) <- cond
  if (simplify && length(ul <- unique(lengths(out, use.names = FALSE))) == 1L){
    nms <- names(out)
    out <- matrix(unlist(out, use.names = FALSE), ncol = ul, byrow = TRUE)
    if (names) rownames(out) <- nms
  }
  out
}

# ==== Method for class 'configTable' ====
# Function suited for interactive use
.redund.configTable <- function(x, cond, simplify = TRUE, full = FALSE){
  .redund.cti(ctInfo(x), cond, simplify = simplify, full = full) 
}

# ==== Default Method (for matrix or data.frame) ====
# builds full.ct if x is NULL
#   x       configTable or NULL
# value:    configTable, mv if original is mv, cs else
.redund.default  <- function(x, cond, simplify = TRUE, ...){
  if (is.null(x)){
    x <- full.ct.default(cond)
  } else {
    x <- full.ct(x)
  }
  .redund.cti(ctInfo(x), cond, simplify = simplify, full = TRUE)
}
