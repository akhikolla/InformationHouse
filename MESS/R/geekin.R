#' Fit a generalized estimating equation (GEE) model with fixed additive
#' correlation structure
#'
#' The geekin function fits generalized estimating equations but where the
#' correlation structure is given as linear function of (scaled) fixed
#' correlation structures.
#'
#' The geekin function is essentially a wrapper function to \code{geeglm}.
#' Through the varlist argument, it allows for correlation structures of the
#' form
#'
#' R = sum_i=1^k alpha_i R_i
#'
#' where alpha_i are(nuisance) scale parameters that are used to scale the
#' off-diagonal elements of the individual correlation matrices, R_i.
#'
#' @aliases geekin print.geekin
#' @param formula See corresponding documentation to \code{glm}.
#' @param family See corresponding documentation to \code{glm}.
#' @param data See corresponding documentation to \code{glm}.
#' @param weights See corresponding documentation to \code{glm}.
#' @param subset See corresponding documentation to \code{glm}.
#' @param id a vector which identifies the clusters. The length of \code{id}
#' should be the same as the number of observations. Data must be sorted so
#' that observations on a cluster are contiguous rows for all entities in the
#' formula. If not the function will give an error
#' @param na.action See corresponding documentation to \code{glm}.
#' @param control See corresponding documentation to \code{glm}.
#' @param varlist a list containing one or more matrix or bdsmatrix objects
#' that represent the correlation structures
#' @param \dots further arguments passed to or from other methods.
#' @return Returns an object of type \code{geeglm}.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{lmekin}, \code{geeglm}
#' @keywords models
#' @examples
#'
#'
#'  # Get dataset
#'  library(kinship2)
#'  library(mvtnorm)
#'  data(minnbreast)
#'
#'  breastpeda <- with(minnbreast[order(minnbreast$famid), ], pedigree(id,
#'                    fatherid, motherid, sex,
#'                    status=(cancer& !is.na(cancer)), affected=proband,
#'                    famid=famid))
#'
#' set.seed(10)
#'
#' nfam <- 6
#' breastped <- breastpeda[1:nfam]
#'
#'  # Simulate a response
#'
#' # Make dataset for lme4
#' df <- lapply(1:nfam, function(xx) {
#'             as.data.frame(breastped[xx])
#'             })
#'
#' mydata <- do.call(rbind, df)
#' mydata$famid <- rep(1:nfam, times=unlist(lapply(df, nrow)))
#'
#' y <- lapply(1:nfam, function(xx) {
#'             x <- breastped[xx]
#'             rmvtnorm.pedigree(1, x, h2=0.3, c2=0)
#'             })
#' yy <- unlist(y)
#'
#' library(geepack)
#'
#' geekin(yy ~ 1, id=mydata$famid, varlist=list(2*kinship(breastped)))
#'
#' # lmekin(yy ~ 1 + (1|id), data=mydata, varlist=list(2*kinship(breastped)),method="REML")
#'
#'
#'
#'
#' @export geekin
geekin <- function(formula,
                   family=gaussian,
                   data,
                   weights,
                   subset,
                   id,
                   na.action,
                   control=geepack::geese.control(...),
                   varlist,
#                   type=c("pearson", "OR"),
                   ...
                   ) {

  # TODO
  # Pass everything correctly to geeglm

  # Check that the input is correct
  if (is.matrix(varlist))
    varlist <- list(varlist)
  if (!is.list(varlist))
    stop("varlist must be a list")

#  type <- match.arg(type)

  # Setup the initial glm call for a model without the correlation structure taken into account
  # This is used to check for missing data and is not necessary otherwise
  Call <- match.call()
  glmcall <- Call
  glmcall$id <- NULL
  glmcall$id <- glmcall$control <- glmcall$corstr <- glmcall$zcor <- glmcall$type <- NULL
  glmcall$varlist <- NULL
  glmcall[[1]] <- as.name("glm")
  glmFit <- eval(glmcall, parent.frame())
  mf <- model.frame(glmFit)


  # Extract clusters from formula
  # Can only accept ONE simple random effect
#  modterm <- attr(terms(formula), "term.labels")
#  clusters <- modterm[grep("\\|", attr(terms(formula), "term.labels"))]
#  if (length(clusters) != 1) {
#    stop("must provide exactly one random effect that determines the clusters")
#  }
#  if (length(grep("^1 \\| ", clusters))!=1) {
#    stop("only accepts a simple random effect to identify the clusters")
#  }

#  id <- clusters
  id <- id

  # Eval det rigtige sted

  # Check subset på variansmatricerne

#  print(get(id))
#  print(id)

  # Check that data are in the correct order
  if (!ordered.clusters(id)) {
    stop("the clusters must appear as contiguous blocks")
  }

  # The updated cluster index
  id.fixed <- id

  # The vector of observations to be removed
  missing.index <- as.numeric(glmFit$na.action)
  if (length(missing.index>0))
    id.fixed <- id.fixed[-missing.index]

  # for each element k in the varlist
  # extract the lower triangular matrix for use with the zcorr specification in geeglm
  varelements <- lapply(varlist, function(x) {
    if (length(missing.index>0)) {
        x <- x[-missing.index,-missing.index]
    }
    lower.tri.vector(x, cluster=id.fixed, diag=FALSE)
  })

  #  print(id.fixed)

  # Setup the proper call to geeglm
  geecall <- Call
  geecall$varlist <- NULL
  geecall$zcor <- do.call("cbind", varelements)
  geecall$corstr <- "userdefined"
  geecall[[1]] <- as.name("geeglm")
  geecall$control <- control
  geecall$id <- id.fixed
  geecall$data <- mf

  geeFit <- eval(geecall, parent.frame())

  # Cheat the output so the varlist and data don't look too horrible
  geeFit$call <- Call
  geeFit$na.action <- glmFit$na.action
  class(geeFit) <- c("geekin", class(geeFit))
  geeFit
}



print.geekin <- function (x, digits = NULL, quote = FALSE, prefix = "", ...)
{
    xg <- x$geese
    if (is.null(digits)) {
        digits <- options()$digits
    } else { options(digits = digits) }
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(unclass(x$coefficients), digits = digits)
    cat("\nDegrees of Freedom:", length(x$y), "Total (i.e. Null); ",
        x$df.residual, "Residual\n")
    if (!xg$model$scale.fix) {
        cat("\nScale Link:                  ", xg$model$sca.link)
        cat("\nEstimated Scale Parameters:  ")
        print(as.numeric(unclass(xg$gamma)), digits = digits)
    } else cat("\nScale is fixed.\n")
    cat("\nCorrelation:  Structure =", xg$model$corstr, " ")
    if (pmatch(xg$model$corstr, "independence", 0) == 0) {
        cat("  Link =", xg$model$cor.link, "\n")
        cat("Estimated Correlation Parameters:\n")
        print(unclass(xg$alpha), digits = digits)
    }
    cat("\nNumber of clusters:  ", length(xg$clusz), "  Maximum cluster size:",
        max(xg$clusz), "\n")
    if (nzchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")

    cat("\n")
    invisible(x)
}
