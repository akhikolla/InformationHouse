
################################################################################
#' @name .met.pow.ratio
#' @author Brian Masinde
#' @param cons Constants defined
#' @param m body mass at start/end of flight
#' @param ws wing span
#' @param ord Order (passerine or non-passerine)
#' @return x2 (metabolic power ratio)
#' @description Metabolic power ratio is a ratio between mechanical power and
#'              the absolute minimum power. The mechanical power is a factor of
#'              basal metabolic rate, which differs between passerines and
#'              non-passerines, and the mechanical conversion efficiency $\etta$.
#'              Basal metabolism is assumed to be needed at all times irrespective
#'              of what the bird is doing.


.met.pow.ratio <- function(cons, m, ws, ordo) {

  a <- ifelse(ordo == 1, cons$alpha[[1]], cons$alpha[[2]])

  d <- ifelse(ordo == 1, cons$delta[[1]], cons$delta[[2]])

  num <-
    6.023 * a * cons$mce * cons$airDensity ^ 0.5 * ws ^ 1.5 * m ^(d - (5 / 3))
  den <-  cons$ipf ^ (3 / 4) * cons$g ^ (5 / 3)

  x2 <- num / den

  return(x2)
}

################################################################################
#' @name .basal.met
#' @author Brian Masinde
#' @param cons Constants defined
#' @param mFrame mass of frame
#' @param mMusc mass of muscle
#' @param ordo (passerine or non-passerine)
#' @return Pbmr (basal metabolism)
#' @description Rate at which fuel energy (not mechanical work) is needed
#' irrespective of what the bird is doing.

.basal.met <- function(cons, mFrame, mMusc, ordo) {
    pbmr <-
      ifelse(ordo == 1, cons$alpha[[1]], cons$alpha[[2]]) * (mFrame + mMusc) ^
      ifelse(ordo == 1, cons$delta[[1]], cons$delta[[2]])

  return(pbmr)
}

################################################################################
#' @name .basal.met2
#' @author Brian Masinde
#' @param cons Constants defined
#' @param bodyMass mass of frame
#' @param fatMass mass of muscle
#' @param ordo Order (passerine or non-passerine)
#' @return Pbmr (basal metabolism)
#' @description Rate at which fuel energy (not mechanical work) is needed
#' irrespective of what the bird is doing.
.basal.met2 <- function(cons, bodyMass, fatMass, ordo) {
  pbmr <-
    ifelse(ordo == 1, cons$alpha[[1]], cons$alpha[[2]]) * (bodyMass - fatMass) ^
    ifelse(ordo == 1, cons$delta[[1]], cons$delta[[2]])

  return(pbmr)
}

################################################################################
#' @name .min.pow.speed
#' @author Brian Masinde
#' @param m all body mass
#' @param ws wing span
#' @param cons constants
#' @return Vmp (minimum power speed)
#' @description Minimum power speed is the minimum speed required for a bird to
#'              maintain horizontal flight. This formula estimates a starting
#'              value but the actual value is estimated by simulation.

.min.pow.speed <- function(m, ws, cons) {
  Vmp <- (0.807 * cons$ipf ^ 0.25 * m ^ 0.5 * cons$g ^ 0.5) /
            (
              cons$airDensity ^ 0.5 * ws ^ 0.5 * .body.front.area(m) ^ 0.25 *
                cons$bdc ^ 0.25
            )

  return(Vmp)
}

################################################################################
#' @name .body.front.area
#' @author Brian Masinde
#' @param  m body mass
#' @return body frontal area
#' @description Body frontal area is cross-sectional area of the widest point.
#'              Used in calculating parasite power of a bird (power required)
#'              to overcome drag by multiplying it by body drag coefficient.
#'              Pennycuick provides a formula that estimates this other than
#'              repetitive measurements. Formula solely depends on body mass.

.body.front.area <- function(m) {
  Sb <- 0.00813 * m ^ 0.666

  return(Sb)
}

################################################################################
#' @name .disc.area
#' @author Brian Masinde
#' @param ws wing span
#' @return disk area of a bird
#' @description Area of circle whose diameter is the wing span.
#'
.disc.area <- function(ws) {
  Sd <- (pi*(ws^2))/4
}

################################################################################
#' @name .wing.frequency
#' @author Brian Masinde
#' @param bm all-up mass
#' @param g gravity
#' @param ws wing span
#' @param wa wing area
#' @param rho air density
#' @description wing beat frequency at specified air density
#'
.wingbeat.freq <- function(bm, ws, wa, cons){
  f <- bm^0.375*cons$g^0.5*ws^(-23/24)*wa^(-1/3)*cons$airDensity^(-3/8)
}

################################################################################
#' @name .max.range.speed
#' @author Brian Masinde
#' @param m all up mass
#' @param ws wing span
#' @param cons constants
#' @description Speed at which the bird covers the most range. However, it is not
#'              necessarily the optimal speed because it requires more power to
#'              maintain this speed. This function gives a lower estimate of the
#'              speed. In flight this speed is identified when effective lift drag
#'              ratio ceases to increases.

.max.range.speed <- function(bm, ws, cons) {
  num <- cons$ipf ^ 0.25 * bm ^ 0.5 * cons$g ^ 0.5

  den <-
    cons$airDensity ^ 0.5 * (cons$bdc * .body.front.area(bm)) ^ 0.25 *
    .disc.area(ws) ^ 0.25

  Vmr <- num / den

  return(Vmr)
}

################################################################################
#' @name .induced.pow
#' @author Brian Masinde
#' @param m all up mass
#' @param ws wing span
#' @param Vt true airspeed
#' @param cons constants
#' @return induced power in horizontal flight
#' @description Induced power is rate at which flight muscles have to provide work.
#'              In the this case, induced power is not considered during hovering
#'              (true airspeed = 0). Instead the minimum power speed is used as
#'              a starting point.


.induced.pow <- function(m, ws, Vt, cons) {
  # sapply(Vt, function(x)
  #   if (x <= 0) {
  #     x  <-  0.1
  #     cat("minimum true airspeed zero for some obesrvations")
  #   })

  pind <- (2 * cons$ipf * (m * cons$g) ^ 2) / (Vt * pi * (ws ^ 2) * cons$airDensity)

  return(pind)
}

################################################################################
#' @name .parasite.pow
#' @author Brian Masinde
#' @param Vt true airspeed
#' @param m all up mass
#' @param cons constants
#' @return ppar parasite power
#' @description Parasite power is the rate at which work must be done in order
#'              to overcome drag of the body. Found by multiplying the speed and
#'              the body drag coefficient.

.parasite.pow <- function(m, Vt, cons) {
  ppar <- (cons$airDensity * Vt^3 * .body.front.area(m) * cons$bdc) / 2

  return(ppar)
}

################################################################################
#' @name .abs.min.pow
#' @author Brian Masinde
#' @param m all up mass
#' @param ws wing span
#' @param cons constants
#' @return  pam (absolute minimum power for an ideal bird to fly at Vmp)
#' @description Power required for a bird to fly at minimum power speed.

.abs.min.pow <- function(m, ws, cons) {
  pam <-
    (1.05 * (cons$ipf ^ 0.75) * (m ^ 1.5) * (cons$g ^ 1.5) * (.body.front.area(m) ^ 0.25) *
       cons$bdc^0.25) / ((cons$airDensity) ^ 0.5 * (ws ^ 1.5))

  return(pam)
}

################################################################################
#' @name .prof.pow.ratio
#' @author Brian Masinde
#' @param ws wing span
#' @param wa wing area
#' @return profile power ratio
#' @description A ratio of profile power to absolute minimum power. Initially,
#'              assigned a value of 1.2 in earlier values but then later noted
#'              that profile power would be proportional to wing area, other
#'              things being constant. Therefore, it calculated as profile power
#'              constant (8.4) divided by the aspect ratio.

.prof.pow.ratio <- function(ws, wa, cons) {
  # ws = wing span
  # wa = wing area
  X1 <- cons$ppc / (ws ^ 2 / wa)

  return(X1)
}

################################################################################
#' @name .prof.pow
#' @author Brian Masinde
#' @param x1 profile power ratio
#' @param amp absolute minimum power
#' @return profile power
#' @description Power needed to overcome the drag of the wings.

.prof.pow <- function(x1, amp) {
  # x1 = profile power ratio
  # amp = results of absolute minimum power
  profilePower <- x1 * amp

  return(profilePower)
}

################################################################################
#' @name .total.mech.power
#' @author Brian Masinde
#' @inheritParams .induced.pow
#' @inheritParams .parasite.pow
#' @inheritParams .prof.pow
#' @return total power
#' @description This is a procedure for calculating the total power of a bird.
#'              Induced, parasite, and profile power are summed to get the total
#'              power at a certain speed. Usually starts from minimum power
#'               speed.


.total.mech.power <- function(m, ws, wa, Vt, cons) {

  # induced power at starting speed
  pind <- .induced.pow(m, ws, Vt, cons)

  # parasite power
  ppar <- .parasite.pow(m, Vt, cons)

  # profile power
  x1 <- .prof.pow.ratio(ws, wa, cons)
  amp <- .abs.min.pow(m, ws, cons)
  ppr <- .prof.pow(x1, amp)

  # total power
  tpow <- pind + ppar + ppr

  return(tpow)
}

################################################################################
#' @name .pow.curve
#' @author Brian Masinde
#' @param bm all-up body mass
#' @param ws wing span
#' @param wa wing area
#' @param tas true airspeed
#' @param cons settings
#' @return Total mechanical power using true airspeed
#' @description In calculating the power curve, the true air-speed is
#'              estimated where total power is stable. Since we are only interested
#'              in the mechanical power where it is stable the estimation ends slightly
#'              above the total mechanical power estimated at start of flight.
#'

.pow.curve <- function(bm, ws, wa, tas, cons) {
  #cons <- cons
  # starting mechanical power
  start_mech_pow <- .total_Mech_Pow_cpp(
    bm = bm,
    ws = ws,
    wa = wa,
    vt = tas, # holding speed constant
    g = cons$g,
    airDensity = cons$airDensity,
    ipf = cons$ipf,
    bdc = cons$bdc,
    ppc = cons$ppc
  )

  # reduce tas slightly
  true_speed <- tas * 0.8
  current_mech_pow <- 0  # power curve estimate of mechanical power
  repeat {
    if ((current_mech_pow > start_mech_pow) == TRUE) {
      break
    } else{
      true_speed <- true_speed + 0.1

      current_mech_pow <- .total_Mech_Pow_cpp(
        bm = bm,
        ws = ws,
        wa = wa,
        vt = true_speed,
        # holding speed constant
        g = cons$g,
        airDensity = cons$airDensity,
        ipf = cons$ipf,
        bdc = cons$bdc,
        ppc = cons$ppc
      )
    }
  }
  return(current_mech_pow)
}
