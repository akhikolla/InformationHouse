#' @title Stopover mass calculator
#' @description During stop-overs birds replenish fat mass. Using simplifications
#'               from Lindström 1991. The implementation here is simplistic in that
#'               muscle mass is not restored as theory and field experiments have
#'               shown.
#' @name stopover.mass.calculator
#' @param bodyMass left-over after running function migrate
#' @param fatMass left-over after running function migrate
#' @param taxon (or order) classified into two categories (passerines and non-passerines)
#' @param duration number of hours spent at stop-over site. This must be an integer see example
#' @return fat_mass, body_mass

#' @export stopover.mass.calculator
#' @examples
#' stopover.mass.calculator(bodyMass = c(2.2, 3.4), fatMass = c(0.34, 0.42),
#' taxon = c(1,2), duration = 36L)


# SOURCE #######################################################################
# Lindström, Å. (1991). Maximum Fat Deposition Rates in Migrating Birds.
# Ornis Scandinavica (Scandinavian Journal of Ornithology), 22(1),
# 12-19. doi:10.2307/3676616
################################################################################

stopover.mass.calculator <-
  function(bodyMass, fatMass, taxon, duration) {

     # check factors in taxon
    if (any(taxon != 1 & taxon != 2)) {
      stop("taxon column values 1 or 2", call. = TRUE)
    }

    if (length(bodyMass) != length(fatMass) || length(bodyMass) != length(taxon)) {
      stop("number of observations doesn't match", call. = TRUE)
    }

    # check duration is an non-zero integer
    if (is.integer(duration) == FALSE | duration < 0) {
      stop("duration is a non-integer, see documentation \n", call. = TRUE)
    }

    results <- list(
      bodyMass = vector(length = length(bodyMass)),
      fatMass = vector(length = length(bodyMass))
    )

    # Get lean mass
    leanMass <- bodyMass - fatMass

    # maximum fat deposition rate based on order
    maxFatDepositionRate <-
      ifelse(taxon == 1, 2.22 * leanMass ^ (-0.27), 2.80 * leanMass ^ (-0.27))

    # check maxFatDepositionRate is not less than zero and not greater than hundred
    if (any(maxFatDepositionRate < 0) |
        any(maxFatDepositionRate > 100)) {
      warning("Maximum fat deposition rate contains values out of bounds [0,1]")
    }

    fatMassGained <-
      leanMass * (maxFatDepositionRate / 100) * duration

    results$bodyMass <- bodyMass + fatMassGained
    results$fatMass <- fatMass + fatMassGained
    # return object of class migrate
    class(results) <- append(class(results), "stopOverMass")
    return(results)
  }
