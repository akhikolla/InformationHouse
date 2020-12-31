# Method 1 on practical range calculation based on Breguets equations
#
# @author Brian Masinde
# @param bodyMass all up mass
# @param wingSpan wing span of bird in metres
# @param fatMass fat mass of bird
# @param ordo Passerine (1) or non-passerine (2)
# @param wingArea area of wing
# @param constants A list of re-definition of constants (i.e *airDensity*,
#             *consume*, *enegry e*, *mechanical mce n*).
# @importFrom utils tail
# @return List with range (in km), constants used and fat fraction
# @include misc_functions.R lookup_table2.R
#

#' @importFrom utils tail

.breguet <- function(bodyMass, wingSpan, fatMass, ordo, wingArea, constants) {

  ##############################################################################
  # fat fraction
  fatFrac <- fatMass/bodyMass

  # metabolic power ratio metPowRatio
  metPowRatio <- .met.pow.ratio(constants, bodyMass, wingSpan, ordo)

  # x1:ppcons/Aspect ratio + metPowRatio:mpratio check for Drag
  # Aspect ratio = wingSpan^2 / wingArea
  # drag is the effective drag force found by interpolation (table 2)
  # add ppratio to metPowRatio and interpolate
  # round off to 2 digits
  table2 <- .gen.table2()

  dFactor <-
    sapply(round((
      .prof.pow.ratio(ws = wingSpan, wa = wingArea, constants) + metPowRatio
    ),
    2), .interpolate, table2)


  ##############################################################################
  # Effective lift:drag ratio
  # Disk area diskArea
  diskArea <- 0.25 * pi * (wingSpan ^ 2)

  # flat-plate area
  flatPlateArea <- 0.00813 * (bodyMass ^ 0.666) * constants$bdc

  # lift drag ratio at beginning of flight
  liftDragRatio <- (dFactor / ((constants$ipf ^ 0.5) * constants$vcp)) *
    ((diskArea / flatPlateArea) ^ 0.5)

  # increase by 10F%
  liftDragRatio <- liftDragRatio + (liftDragRatio * (10 * fatFrac) / 100)

  # range in kilometres
  kmRange <-
    ((constants$fed * constants$mce) / constants$g) * liftDragRatio *
    log(1 / (1 - fatFrac))/1000

  return(round(kmRange, 1))
}
