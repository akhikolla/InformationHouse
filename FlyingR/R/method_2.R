# Method 2 practical range calculation based on Breguets equations with mean
#  of effective lift: drag ratio
# @author Brian Masinde
# @param bodyMass All up mass. Including fuel, crop contents and equipment in Kg
# @param wingSpan Wing span in metres
# @param fatMass Fat mass in kg (fuel)
# @param ordo Passerine (1) or non-passerine (2)
# @param wingArea Wing area
# @param ctrl A list of re-definition of constants (i.e *airDensity*,
#             *consume*, *enegry e*, *mechanical mce n*).
# @importFrom utils tail
# @include misc_functions.R lookup_table2.R
# @description Practical range estimation using Breguet equation for fixed wing
#              with crude adjustments. Mean lift:drag ratio between start and
#              end of flight is used as proxy for engine burn.

#' @importFrom utils tail
.breguet_adj <- function(bodyMass, wingSpan, fatMass, ordo, wingArea, constants) {
  ##############################################################################
  # fat fraction
  fatFrac <- fatMass/bodyMass

  ## lift:drag end of flight ###################################################
  # m2 mass end of flight
  #bodyMassEnd <- bodyMass - (fatMass * constants$consume)
  bodyMassEnd <- bodyMass - fatMass

  # x2
  metPowRatioEnd <- .met.pow.ratio(constants, bodyMassEnd, wingSpan, ordo)

  # x1:ppcons/Aspect ratio + x2:mpratio check for D ----------------------------
  # Aspect ratio = wingSpan^2 / wingArea
  # D is the effective drag force found by interpolation (table 2)
  # add ppratio to x2 and interpolate
  # round off to 2 digits
  table2 <- .gen.table2()

  # dFactorEnd <- sapply(round((constants$ppcons / (wingSpan^2/wingArea)) +
  #                              metPowRatioEnd, 2), .interpolate, table2)
  dFactorEnd <-
    sapply(round((
      .prof.pow.ratio(ws = wingSpan, wa = wingArea, constants) + metPowRatioEnd
    ),
    2), .interpolate, table2)

  # dFactorEnd <- sapply(round(1.2 +
  #                              metPowRatioEnd, 2), .interpolate, table2)

  ### Ask if we should round off when interpolating

  # disk area
  diskArea <- 0.25 * pi * (wingSpan ^ 2)

  # flat-plate area
  flatPlateAreaEnd <- 0.00813 * (bodyMassEnd ^ 0.666) * constants$bdc

  # lift drag ratio at begining of flight
  liftDragEnd <-
    dFactorEnd / (constants$ipf ^ 0.5 * constants$vcp) * ((diskArea / flatPlateAreaEnd) ^ 0.5)


  ## lift:drag ratio start of flight ###########################################
  # why not calculate metPowRatioStart using the funciton but with bodymass at
  # start
  metPowRatioStart <- metPowRatioEnd / ((1 / (1 - fatFrac)) ^ 1.75)

  # dFactorStart <- sapply(round((constants$ppcons / (wingSpan^2/wingArea)) +
  #                                metPowRatioStart, 2), .interpolate, table2)

  dFactorStart <-
    sapply(round((
      .prof.pow.ratio(ws = wingSpan, wa = wingArea, constants) + metPowRatioStart
    ),
    2), .interpolate, table2)

  # dFactorStart <- sapply(round(1.2 +
  #                                metPowRatioStart, 2), .interpolate, table2)

  liftDragStart <-
    (dFactorStart / ((constants$ipf ^ 0.5) * constants$vcp)) *
    (((diskArea / flatPlateAreaEnd) ^ 0.5) / ((bodyMass / bodyMassEnd) ^ 0.5))


  ## Range in km ###############################################################
  kmRange <-
    ((constants$fed * constants$mce) / constants$g) *
    apply(cbind(liftDragStart, liftDragEnd), 1, mean) *
    log(1 / (1 - fatFrac)) / 1000

  return(round(kmRange, 1))
}

