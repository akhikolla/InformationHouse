
# function cyclic():
# Determine if a csf is cyclic
#   x must be a character vector
# Returns a logical vector
cyclic <- function(x, cycle.type = c("factor", "value"), use.names = TRUE, verbose = FALSE){
  stopifnot(is.character(x))
  if (length(x) == 0) return(logical(0))
  x <- noblanks(x)
  cycle.type <- match.arg(cycle.type)

  # Check input conditions [identical code in is.submodel()]
  px <- lapply(x, tryparse)
  ok <- !vapply(px, is.null, logical(1)) 
  ct_type <- if (any(grepl("=", x, fixed = T))) "mv" else "cs"
  if (ct_type == "cs"){
    vals <- unique.default(unlist(lapply(px, all.vars)))
  } else {
    vals <- rapply(px, .call2list, how = "unlist",
                   stopOps = c("==", "<", ">", "<=", ">="), 
                   validOps = c("<-", "<<-", "=", "&", "|", "(", "-"))
    vals <- unique.default(vals[!vapply(vals, is.symbol, FUN.VALUE = TRUE)])
    vals <- sub(" == ", "=", vapply(vals, deparse, character(1)))
  }
  cond_type <- .qcondType(x, values = vals, ct_type = ct_type, stdComplex.multiple.only = FALSE)
  ok <- ok & cond_type %in% c("stdAtomic", "stdComplex")
  if (any(!ok)) 
    stop("Invalid input to cyclic:\n", paste0("  ", x[!ok], collapse = "\n"),
         call. = FALSE)

  # check cyclicity
  fStr <- factorStruct(x, cycle.type = cycle.type, ct_type = ct_type)
  fStrGr <- factorStructGroups(fStr)
  # .cyclic1() executed only once per fStr
  ii <- which(!duplicated(fStrGr))
  outCycl <- logical(length(ii))
  cntr <- 0L
  for (i in ii){
    if (verbose) cat("---", x[i], "---\n")
    outCycl[[cntr <- cntr+1L]] <- .cyclic1(fStr[[i]], verbose = verbose)
  }
  # Expand to length(x)
  out <- outCycl[fStrGr]
  if (use.names) names(out) <- x
  out
}


# .cyclic1
# Processes one csf, coded as factorStruct (see below)
.cyclic1 <- function(fstr, verbose = FALSE){
  r <- names(fstr)
  l <- fstr
  ul <- unlist(l, use.names = FALSE)
  
  # "causal successors" of all factors
  successors <- split.default(rep(r, lengths(l)), ul)
  successors <- append(successors, rep_names(list(""), setdiff(r, ul)))
  successors[] <- lapply(successors, unique.default)
  # Auxiliary function step()
  .step <- function(s){
    current <- vapply(s, "[", 1, FUN.VALUE = character(1))
    succ <- successors[current]
    out <- mapply(c, unlist(succ, use.names = FALSE), rep(s, lengths(succ)), 
                  SIMPLIFY = FALSE, USE.NAMES = FALSE)
    out <- unique(out)
    out
  }
    
  # starting points
  to_check <- intersect(r, ul)  
  found_cycle <- FALSE
  
  # Expand "causal paths"
  while (length(to_check) && !found_cycle){

    states <- to_check[[1]]
    visited <- states  # Filter(nzchar, unique(unlist(s2)))
    
    repeat {
      states_new <- .step(states)
      to_follow <- nzchar(current <- vapply(states_new, "[", 1, FUN.VALUE = character(1)))
      
      if (!any(to_follow)){  # -> ok so far, continue with to_check
        break
      }
      if (any(unlist(lapply(states_new, duplicated), use.names = FALSE))){ # -> cyclic!!
        found_cycle <- TRUE
      }
      visited_new <- Filter(nzchar, union(visited, unlist(states_new, use.names = FALSE)))

      states <- states_new
      visited <- visited_new
      
      if (found_cycle)
        break
    }
    
    if (!found_cycle) to_check <- setdiff(to_check, visited)

    if (verbose){
      writeLines(vapply(states, function(x) C_concat(rev(Filter(nzchar, x)), " > "),
                        FUN.VALUE = character(1)))
    }
  }
  if (verbose) cat("\n")
  found_cycle
}


# extract factor and outcome configurations from csfs
factorStruct <- function(x, cycle.type = "factor", ct_type = "cs"){
  if (cycle.type == "factor" && ct_type == "cs")
    x[] <- toupper(x)
  if (cycle.type == "factor" && ct_type == "mv")
    x[] <- gsub("=[0-9]+", "", x)
  # extract asf
  asfs0 <- extract_asf(x)
  # Remove duplicate factor (values)
  asfs1 <- happly(asfs0, function(x) strsplit(lhs(x), "\\+|\\*"))
  asfs1[] <- rapply(asfs1, unique.default, how = "replace", classes = "character")
  # Sort factor (values) within asf
  ll <- hlengths(asfs1)
  ul_factors <- unlist(asfs1, use.names = FALSE, recursive = TRUE)
  asfs1 <- hrelist(ul_factors[order(rep(seq_along(ll[[2]]), ll[[2]]), ul_factors)], ll)
  # Set outcomes as names of asf
  rr <- happly(asfs0, rhs)
  relst <- relist1(setNames(unlist(asfs1, recursive = FALSE), unlist(rr)),
                   ll[[1]])
  # Sort asf within csf
  ul_asfs <- unlist(relst, use.names = TRUE, recursive = FALSE)
  ord <- order(rep(seq_along(ll[[1]]), ll[[1]]), names(ul_asfs))
  C_relist_List(ul_asfs[ord], ll[[1]])
}


# categorize factorStructs into equivalent groups
# csfs within these groups will be identical wrt cyclicity
factorStructGroups <- function(x){
  ux <- unlist(x, recursive = FALSE)
  uxchar <- paste(C_mconcat(ux, ","), names(ux), sep = ">")
  xchar <- C_mconcat(relist1(uxchar, lengths(x)), ";")
  as.integer(factor(xchar, levels = unique(xchar)))
}


# Aux fun
rep_names <- function(x, nms){
  out <- rep_len(x, length(nms))
  names(out) <- nms
  out
}
