# RPP Algorithm


#' Cluster and Subcluster effects
#' 
#' Partitions model residuals into cluster and subcluster effects using RPP
#' algorithm.
#' 
#' 
#' @param f A model formula which specifices the random effects (see example)
#' @param resd The residuals from the fitted model
#' @param data The data the model was fitted on
#' @param rprpair Character string indicating the location and scale parameters
#' to use. Default to "hl-disp", but may also be "med-mad". See Bilgic (2012).
#' @return
#'   \item{siga2}{ Variance from cluster }
#'   \item{sigw2}{ Variance from subcluster }
#'   \item{sigmae2}{ Remaining variance not accounted for by variance of cluster and subcluster }
#'   
#' @author J. W. McKean and Y. K. Bilgic
#' 
#' @seealso rprmeddis, dispvar
#' 
#' @references Y. K. Bilgic. Rank-based estimation and prediction for mixed
#' effects models in nested designs. 2012. URL
#' http://scholarworks.wmich.edu/dissertations/40. Dissertation.
#' 
#' @examples
#' 
#' # Load school data
#' data(schools)
#' 
#' # Fit fixed effects model with lmr
#' lmr.fit = lmr(y ~ age + sex, data=schools)
#' 
#' # Three level design
#' # Partition residuals into school and region effects with rpp algorithm
#' rpr(y ~ age + sex + (1 | school) + (1 | school:region), lmr.fit$ehat, schools)
#' 
#' # Two level design
#' # Estimate variance in residuals from school
#' rpr(y ~ age + sex + (1 | school), lmr.fit$ehat, schools)
#' 
#' @importFrom stringr str_extract_all
#' @export
rpr <- function(f, resd, data, rprpair='hl-disp') {
  random_terms = str_extract_all(as.character(f)[[3]], "\\(([0-9 a-zA-Z.:|]+)\\)")
  random_terms = lapply(random_terms[[1]], function(str) substr(str, 
                        2, nchar(str) - 1))
  
  rprpair = tolower(rprpair)
  location = scale = 2
  if (rprpair == "med-mad") {
      location = scale = 1
  }
  
  if (length(random_terms) == 2) {
    levels = 3
    school_name = tail(strsplit(random_terms[[1]], "\\W")[[1]], 
                       1)
    section_name = tail(strsplit(random_terms[[2]], "\\W")[[1]], 
                        1)
    I = length(unique(factor(data[[school_name]])))
    sec = as.vector(sec_vec(data[[school_name]], data[[section_name]]))
    mat = mat_vec(data[[school_name]], data[[section_name]])
    
    return(rprmeddis(I, sec, mat, resd, location, scale))
  }
  if (length(random_terms) == 1) {
    levels = 2
    school_name = tail(strsplit(random_terms[[1]], "\\W")[[1]], 
                       1)
    I = length(unique(factor(data[[school_name]])))
    mat = mat_vec(data[[school_name]], rep(1, length(data[[school_name]])))
    
    return(rprmeddis2(I, c(), mat, resd, location, scale))
  }
}
