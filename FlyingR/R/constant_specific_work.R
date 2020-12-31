# constant specific work time marching computation
# @author Brain Masinde.
# @name .constant.specific.work
# @param data Data as output from .colnames.match
# @param constants
# @param speed_control speed control as either

.constant.specific.work <- function(data, constants, speed_control, min_energy_protein) {

# defensives ###################################################################
  if (missing(data) == TRUE) {
    stop("Missing data argument", call. = FALSE)
  }

  if (missing(constants) == TRUE) {
    stop("Missing constants", call. = FALSE)
  }

  if (missing(speed_control) == TRUE) {
    stop("Missing speed control method", call. = FALSE)
  }

  if (speed_control != 1 && speed_control != 0) {
    stop("speed control should either be 1 or 0")
  }

  # muscle mass is a must for this function
  if (is.null(data$muscleMass)) {
    stop("Muscle mass column missing", call.FALSE = TRUE)
  }
################################################################################
  # number of observations
  n <- nrow(data)

  allMass <- data$allMass

  fatMass <- data$fatMass

  wingSpan <- data$wingSpan

  wingArea <- data$wingArea

  muscleMass <- data$muscleMass

  taxon <- data$taxon

# time marching ################################################################

  # basal metabolic constants
  alphaPasserines <- constants$alpha[1]
  alphaNonPasserines <- constants$alpha[2]
  deltaPasserines <- constants$delta[1]
  deltaNonPasserines <- constants$delta[2]

  results <- list(
    distance = vector(length = n),
    allUpMass = vector(length = n),
    fatMass = vector(length = n),
    muscleMass = vector(length = n),
    startMinSpeed = vector(length = n),
    endMinSpeed = vector(length = n)
  )

  if (speed_control == 1 && min_energy_protein == 0) {
    for (i in seq_len(nrow(data))) {
      # things to keep track of ################################################
      bm <- allMass[i]
      fm <- fatMass[i]
      mm <- muscleMass[i]
      airframeMass <- allMass[i] - (fatMass[i] + muscleMass[i]) # airframe mass
      dist <- 0

      # true start speed
      startMinSpeed <- .minpowspeed_cpp(
        bm = bm,
        ws = wingSpan[i],
        ipf = constants$ipf,
        g = constants$g,
        airDensity = constants$airDensity,
        bdc = constants$bdc
      )

      # we want to hold true speed constant
      trueSpeed <- startMinSpeed * constants$speedRatio

      mechPower <-
        .mechanical_power(
          bm = bm,
          ws = wingSpan[i],
          wa = wingArea[i],
          tas =  trueSpeed,
          g = constants$g,
          airDensity = constants$airDensity,
          ipf = constants$ipf,
          bdc = constants$bdc,
          ppc = constants$ppc
        )

      wingFreq <-
        .wingbeat.freq(bm = bm,
                       ws = wingSpan[i],
                       wa = wingArea[i],
                       constants)
      # subdivide muscle mass into respective components #######################
      myofibrils <-
        muscleMass[i] * (1 - constants$mipd * (mechPower /
                                                 muscleMass[i]) * constants$muscDensity)
      mitochondria <- muscleMass[i] - myofibrils

      # specific work at start of flight #######################################
      specWorkStart <- mechPower / (myofibrils * wingFreq)

      # return starting minimum power speed ####################################
      results$startMinSpeed[i] <- startMinSpeed

      while (fm > 0.000001) {
        # mechanical power from power curve holding true air-speed constant ####
        mechPower <-
          .mechanical_power(
            bm = bm,
            # changes from previous J iteration
            ws = wingSpan[i],
            wa = wingArea[i],
            tas =  trueSpeed,
            g = constants$g,
            airDensity = constants$airDensity,
            ipf = constants$ipf,
            bdc = constants$bdc,
            ppc = constants$ppc
          )

        # wing frequency
        wingFreq <-
          .wingbeat.freq(bm = bm,
                         # changes from previous J iteration
                         ws = wingSpan[i],
                         wa = wingArea[i],
                         constants)

        # chemical power #######################################################
        chemPower <- constants$vcp * (mechPower + constants$mce * .basal_metabolic_pow(
          airframeMass,
          # doesn't change from previous J iteration
          mm,
          # changes from previous J iteration
          taxon[i],
          alphaPasserines,
          alphaNonPasserines,
          deltaPasserines,
          deltaNonPasserines
        )) / constants$mce

        # fat consumed in the interval?
        usedFat <- chemPower / constants$fed * 360

        #Reduce body composition by this consumed fat and determine what amount
        # of protein to consume to achieve specific power at start of flight
        bmDummy <- bm - usedFat

        # because of this reduction minimum speed, mechanical power, and
        # wing frequency reduce
        trueSpeedDummy <-
          .minpowspeed_cpp(
            bm = bmDummy,
            ws = wingSpan[i],
            ipf = constants$ipf,
            g = constants$g,
            airDensity = constants$airDensity,
            bdc = constants$bdc
          ) * constants$speedRatio

        mechPowerDummy <-
          .mechanical_power(
            bm = bmDummy,
            ws = wingSpan[i],
            wa = wingArea[i],
            tas = trueSpeedDummy,
            g = constants$g,
            airDensity = constants$airDensity,
            ipf = constants$ipf,
            bdc = constants$bdc,
            ppc = constants$ppc
          )

        # wing frequency
        wingFreqDummy <-
          .wingbeat.freq(bm = bmDummy,
                         ws = wingSpan[i],
                         wa = wingArea[i],
                         constants)

        # amount of Myofibrils consumed that restores specific work to initial value
        usedMyofibrils <-
          -(mechPowerDummy / (wingFreqDummy * specWorkStart)) + myofibrils

        # amount of fuel energy released is found by multiplying the mass of dry protein
        # removed by the energy density of dry protein. Efficiency is involved in
        # this conversion.
        usedFatEquiv <-
          (usedMyofibrils * constants$ped * constants$mce) / constants$fed

        # update body measurements #############################################
        fm <- fm - (usedFat - usedFatEquiv)
        myofibrils <- myofibrils - usedMyofibrils
        mm <- mitochondria + myofibrils - (usedMyofibrils * constants$phr)
        bm <- airframeMass + mm + fm
        # distance increment ###################################################
        dist <- dist + trueSpeed * 360
      }
      results$distance[i] <- dist
      results$allUpMass[i] <- bm
      results$fatMass[i] <- fm
      results$muscleMass[i] <- mm
      results$endMinSpeed[i] <- .minpowspeed_cpp(
        bm = bm,
        ws = wingSpan[i],
        ipf = constants$ipf,
        g = constants$g,
        airDensity = constants$airDensity,
        bdc = constants$bdc
      )
    }
  }  else if (speed_control == 1 && min_energy_protein > 0) {
    for (i in seq_len(nrow(data))) {
      # things to keep track of ################################################
      bm <- allMass[i]
      fm <- fatMass[i]
      mm <- muscleMass[i]
      airframeMass <- allMass[i] - (fatMass[i] + muscleMass[i])
      dist <- 0

      # true start speed
      startMinSpeed <- .minpowspeed_cpp(
        bm = bm,
        ws = wingSpan[i],
        ipf = constants$ipf,
        g = constants$g,
        airDensity = constants$airDensity,
        bdc = constants$bdc
      )

      # we want to hold true speed constant
      trueSpeed <- startMinSpeed * constants$speedRatio

      mechPower <-
        .mechanical_power(
          bm = bm,
          ws = wingSpan[i],
          wa = wingArea[i],
          tas =  trueSpeed,
          g = constants$g,
          airDensity = constants$airDensity,
          ipf = constants$ipf,
          bdc = constants$bdc,
          ppc = constants$ppc
        )

      wingFreq <-
        .wingbeat.freq(bm = bm,
                       ws = wingSpan[i],
                       wa = wingArea[i],
                       constants)
      # subdivide muscle mass into respective components #######################
      myofibrils <-
        muscleMass[i] * (1 - constants$mipd * (mechPower /
                                                 muscleMass[i]) * constants$muscDensity)
      mitochondria <- muscleMass[i] - myofibrils

      # specific work at start of flight #######################################
      specWorkStart <- mechPower / (myofibrils * wingFreq)

      # return starting minimum power speed ####################################
      results$startMinSpeed[i] <- startMinSpeed

      while (fm > 0.000001) {
        # mechanical power from power curve holding true air-speed constant
        mechPower <-
          .mechanical_power(
            bm = bm,
            ws = wingSpan[i],
            wa = wingArea[i],
            tas =  trueSpeed,
            g = constants$g,
            airDensity = constants$airDensity,
            ipf = constants$ipf,
            bdc = constants$bdc,
            ppc = constants$ppc
          )
        #cat("mech power within while", mechPower, sep = " ", "\n")
        # wing frequency
        wingFreq <-
          .wingbeat.freq(bm = bm,
                         ws = wingSpan[i],
                         wa = wingArea[i],
                         constants)

        # chemical power #######################################################
        chemPower <- constants$vcp * (mechPower + constants$mce * .basal_metabolic_pow(
          airframeMass,
          # doesn't change from previous J iteration
          mm,
          # changes from previous J iteration
          taxon[i],
          alphaPasserines,
          alphaNonPasserines,
          deltaPasserines,
          deltaNonPasserines
        )) / constants$mce


        # fat consumed in the interval?
        usedFat <- (chemPower  / constants$fed) * 360
        #Reduce body composition by this consumed fat and determine what amount
        # of protein to consume to achieve specific power at start of flight
        bmDummy <- bm - usedFat

        # because of this reduction minimum speed, mechanical power, and
        # wing frequency reduce

        trueSpeedDummy <-
          .minpowspeed_cpp(
            bm = bmDummy,
            ws = wingSpan[i],
            ipf = constants$ipf,
            g = constants$g,
            airDensity = constants$airDensity,
            bdc = constants$bdc
          ) * constants$speedRatio

        mechPowerDummy <-
          .mechanical_power(
            bm = bmDummy,
            ws = wingSpan[i],
            wa = wingArea[i],
            tas = trueSpeedDummy,
            g = constants$g,
            airDensity = constants$airDensity,
            ipf = constants$ipf,
            bdc = constants$bdc,
            ppc = constants$ppc
          )

        wingFreqDummy <-
          .wingbeat.freq(bm = bmDummy,
                         ws = wingSpan[i],
                         wa = wingArea[i],
                         constants)

        # amount of protein consumed that restores specific work to initial value
        usedMyofibrils <-
          -(mechPowerDummy / (wingFreqDummy * specWorkStart)) + myofibrils

        # amount of fuel energy released is found by multiplying the mass of dry protein
        # removed by the energy density of dry protein
        usedFatEquiv <-
          (usedMyofibrils * constants$ped * constants$mce) / constants$fed

        # total energy produced in the interval
        totalEnergy <- chemPower * 360 + usedMyofibrils * constants$ped
        # % of protein energy remaining for minimum energy from protein to be met
        # this will come from airframe mass.
        meetProtein <-
          min_energy_protein - (usedMyofibrils * constants$ped * constants$mce) / totalEnergy

        # adjust body components ###############################################
        fm <- fm - (usedFat - usedFatEquiv)
        myofibrils <- myofibrils - usedMyofibrils
        mm <- mitochondria + myofibrils - (usedMyofibrils * constants$phr)
        airframeMass <- airframeMass - (((totalEnergy * meetProtein) / constants$ped) * constants$phr) - ((totalEnergy * meetProtein)/ constants$ped)
        bm <- airframeMass + mm + fm
        # distance increment ###################################################
        dist <- dist + trueSpeed * 360
      }
      results$distance[i] <- dist
      results$allUpMass[i] <- bm
      results$fatMass[i] <- fm
      results$muscleMass[i] <- mm
      results$endMinSpeed[i] <- .minpowspeed_cpp(
        bm = bm,
        ws = wingSpan[i],
        ipf = constants$ipf,
        g = constants$g,
        airDensity = constants$airDensity,
        bdc = constants$bdc
      )
    } # closes for loop
  } else if (speed_control == 0 && min_energy_protein == 0) {
    for (i in seq_len(nrow(data))) {
      # things to keep track of ################################################
      bm <- allMass[i]
      fm <- fatMass[i]
      mm <- muscleMass[i]
      airframeMass <- allMass[i] - (fatMass[i] + muscleMass[i]) # airframe mass
      dist <- 0

      # true start speed
      startMinSpeed <- .minpowspeed_cpp(
        bm = bm,
        ws = wingSpan[i],
        ipf = constants$ipf,
        g = constants$g,
        airDensity = constants$airDensity,
        bdc = constants$bdc
      )

      # we want to hold ratio of minimum power speed and true speed constant
      trueSpeed <- startMinSpeed * constants$speedRatio

      mechPower <-
        .mechanical_power(
          bm = bm,
          ws = wingSpan[i],
          wa = wingArea[i],
          tas =  trueSpeed,
          g = constants$g,
          airDensity = constants$airDensity,
          ipf = constants$ipf,
          bdc = constants$bdc,
          ppc = constants$ppc
        )

      wingFreq <-
        .wingbeat.freq(bm = bm,
                       ws = wingSpan[i],
                       wa = wingArea[i],
                       constants)
      # subdivide muscle mass into respective components #######################
      myofibrils <-
        muscleMass[i] * (1 - constants$mipd * (mechPower /
                                                 muscleMass[i]) * constants$muscDensity)
      mitochondria <- muscleMass[i] - myofibrils

      # specific work at start of flight #######################################
      specWorkStart <- mechPower / (myofibrils * wingFreq)

      # return starting minimum power speed ####################################
      results$startMinSpeed[i] <- startMinSpeed

      while (fm > 0.000001) {
        # hold ratio of minimum power speed and true airspeed constant
        # true start speed
        minSpeed <- .minpowspeed_cpp(
          bm = bm, # bm changes after each J iteration
          ws = wingSpan[i],
          ipf = constants$ipf,
          g = constants$g,
          airDensity = constants$airDensity,
          bdc = constants$bdc
        )

        # we want to hold ratio of minimum power speed and true speed constant
        trueSpeed <- minSpeed * constants$speedRatio

        # mechanical power from power curve holding true air-speed constant ####
        mechPower <-
          .mechanical_power(
            bm = bm,
            # changes from previous J iteration
            ws = wingSpan[i],
            wa = wingArea[i],
            tas =  trueSpeed,
            g = constants$g,
            airDensity = constants$airDensity,
            ipf = constants$ipf,
            bdc = constants$bdc,
            ppc = constants$ppc
          )

        # wing frequency
        wingFreq <-
          .wingbeat.freq(bm = bm,
                         # changes from previous J iteration
                         ws = wingSpan[i],
                         wa = wingArea[i],
                         constants)

        # chemical power #######################################################
        chemPower <- constants$vcp * (mechPower + constants$mce * .basal_metabolic_pow(
          airframeMass,
          # doesn't change from previous J iteration
          mm,
          # changes from previous J iteration
          taxon[i],
          alphaPasserines,
          alphaNonPasserines,
          deltaPasserines,
          deltaNonPasserines
        )) / constants$mce

        # fat consumed in the interval?
        usedFat <- chemPower / constants$fed * 360

        #Reduce body composition by this consumed fat and determine what amount
        # of protein to consume to achieve specific power at start of flight
        bmDummy <- bm - usedFat

        # because of this reduction minimum speed, mechanical power, and
        # wing frequency reduce
        trueSpeedDummy <-
          .minpowspeed_cpp(
            bm = bmDummy,
            ws = wingSpan[i],
            ipf = constants$ipf,
            g = constants$g,
            airDensity = constants$airDensity,
            bdc = constants$bdc
          ) * constants$speedRatio

        mechPowerDummy <-
          .mechanical_power(
            bm = bmDummy,
            ws = wingSpan[i],
            wa = wingArea[i],
            tas = trueSpeedDummy,
            g = constants$g,
            airDensity = constants$airDensity,
            ipf = constants$ipf,
            bdc = constants$bdc,
            ppc = constants$ppc
          )

        # wing frequency
        wingFreqDummy <-
          .wingbeat.freq(bm = bmDummy,
                         ws = wingSpan[i],
                         wa = wingArea[i],
                         constants)

        # amount of protein consumed that restores specific work to initial value
        usedMyofibrils <-
          -(mechPowerDummy / (wingFreqDummy * specWorkStart)) + myofibrils

        # amount of fuel energy released is found by multiplying the mass of dry protein
        # removed by the energy density of dry protein
        usedFatEquiv <-
          (usedMyofibrils * constants$ped * constants$mce) / constants$fed

        # update body measurements #############################################
        fm <- fm - (usedFat - usedFatEquiv)
        myofibrils <-  myofibrils - usedMyofibrils
        mm <- mitochondria + myofibrils - (usedMyofibrils * constants$phr)
        bm <- airframeMass + mm + fm

        # distance increment ###################################################
        dist <- dist + trueSpeed * 360
      }
      results$distance[i] <- dist
      results$allUpMass[i] <- bm
      results$fatMass[i] <- fm
      results$muscleMass[i] <- mm
      results$endMinSpeed[i] <- .minpowspeed_cpp(
        bm = bm,
        ws = wingSpan[i],
        ipf = constants$ipf,
        g = constants$g,
        airDensity = constants$airDensity,
        bdc = constants$bdc
      )
    } # closes for loop (iterates over rows of species)
  } else if (speed_control == 0 && min_energy_protein > 0) {
    for (i in seq_len(nrow(data))) {
      # things to keep track of ################################################
      bm <- allMass[i]
      fm <- fatMass[i]
      mm <- muscleMass[i]
      airframeMass <- allMass[i] - (fatMass[i] + muscleMass[i]) # airframe mass
      dist <- 0

      # true start speed
      startMinSpeed <- .minpowspeed_cpp(
        bm = bm,
        ws = wingSpan[i],
        ipf = constants$ipf,
        g = constants$g,
        airDensity = constants$airDensity,
        bdc = constants$bdc
      )

      # we want to hold true speed constant
      trueSpeed <- startMinSpeed * constants$speedRatio

      mechPower <-
        .mechanical_power(
          bm = bm,
          ws = wingSpan[i],
          wa = wingArea[i],
          tas =  trueSpeed,
          g = constants$g,
          airDensity = constants$airDensity,
          ipf = constants$ipf,
          bdc = constants$bdc,
          ppc = constants$ppc
        )

      wingFreq <-
        .wingbeat.freq(bm = bm,
                       ws = wingSpan[i],
                       wa = wingArea[i],
                       constants)
      # subdivide muscle mass into respective components #######################
      myofibrils <-
        muscleMass[i] * (1 - constants$mipd * (mechPower /
                                                 muscleMass[i]) * constants$muscDensity)
      mitochondria <- muscleMass[i] - myofibrils

      # specific work at start of flight #######################################
      specWorkStart <- mechPower / (myofibrils * wingFreq)

      # return starting minimum power speed ####################################
      results$startMinSpeed[i] <- startMinSpeed

      while (fm > 0.000001) {
        # mechanical power from power curve holding true air-speed constant
        # hold ratio of minimum power speed and true airspeed constant
        # true start speed
        minSpeed <- .minpowspeed_cpp(
          bm = bm, # bm changes after each J iteration
          ws = wingSpan[i],
          ipf = constants$ipf,
          g = constants$g,
          airDensity = constants$airDensity,
          bdc = constants$bdc
        )

        # we want to hold ratio of minimum power speed and true speed constant
        trueSpeed <- minSpeed * constants$speedRatio
        mechPower <-
          .mechanical_power(
            bm = bm,
            ws = wingSpan[i],
            wa = wingArea[i],
            tas =  trueSpeed,
            g = constants$g,
            airDensity = constants$airDensity,
            ipf = constants$ipf,
            bdc = constants$bdc,
            ppc = constants$ppc
          )

        # wing frequency
        wingFreq <-
          .wingbeat.freq(bm = bm,
                         ws = wingSpan[i],
                         wa = wingArea[i],
                         constants)

        # chemical power #######################################################
        chemPower <- constants$vcp * (mechPower + constants$mce * .basal_metabolic_pow(
          airframeMass,
          # doesn't change from previous J iteration
          mm,
          # changes from previous J iteration
          taxon[i],
          alphaPasserines,
          alphaNonPasserines,
          deltaPasserines,
          deltaNonPasserines
        )) / constants$mce

        # fat consumed in the interval?
        usedFat <- chemPower / constants$fed * 360

        #Reduce body composition by this consumed fat and determine what amount
        # of protein to consume to achieve specific work at start of flight
        bmDummy <- bm - usedFat

        # because of this reduction minimum speed, mechanical power, and
        # wing frequency reduce

        trueSpeedDummy <-
          .minpowspeed_cpp(
            bm = bmDummy,
            ws = wingSpan[i],
            ipf = constants$ipf,
            g = constants$g,
            airDensity = constants$airDensity,
            bdc = constants$bdc
          ) * constants$speedRatio

        mechPowerDummy <-
          .mechanical_power(
            bm = bmDummy,
            ws = wingSpan[i],
            wa = wingArea[i],
            tas = trueSpeedDummy,
            g = constants$g,
            airDensity = constants$airDensity,
            ipf = constants$ipf,
            bdc = constants$bdc,
            ppc = constants$ppc
          )

        wingFreqDummy <-
          .wingbeat.freq(bm = bmDummy,
                         ws = wingSpan[i],
                         wa = wingArea[i],
                         constants)

        # amount of protein consumed that restores specific work to initial value
        usedMyofibrils <-
          -(mechPowerDummy / (wingFreqDummy * specWorkStart)) + myofibrils

        # amount of fuel energy released is found by multiplying the mass of dry protein
        # removed by the energy density of dry protein
        usedFatEquiv <- (usedMyofibrils * constants$ped * constants$mce) / constants$fed

        # total energy produced in the interval
        totalEnergy <- chemPower * 360 + usedMyofibrils * constants$ped

        # % protein energy remaining for minimum % energy from protein to be achieved
        meetProtein <-
          min_energy_protein - (usedMyofibrils * constants$ped * constants$mce) / totalEnergy

        #cat("energy met", meetProtein, sep = " ", "\n")

        # adjust body components ###############################################
        fm <- fm - (usedFat - usedFatEquiv)
        myofibrils <- myofibrils - usedMyofibrils
        mm <- mitochondria + myofibrils - (usedMyofibrils * constants$phr)
        airframeMass <- airframeMass - (((totalEnergy * meetProtein) / constants$ped) * constants$phr) - ((totalEnergy * meetProtein)/ constants$ped)
        bm <- airframeMass + mm + fm

        # distance increment ###################################################
        dist <- dist + trueSpeed * 360
        # increase counter
      }
      results$distance[i] <- dist
      results$allUpMass[i] <- bm
      results$fatMass[i] <- fm
      results$muscleMass[i] <- mm
      results$endMinSpeed[i] <- .minpowspeed_cpp(
        bm = bm,
        ws = wingSpan[i],
        ipf = constants$ipf,
        g = constants$g,
        airDensity = constants$airDensity,
        bdc = constants$bdc
      )
    }
  } # closes conditioning speed_control and min_energy_protein
  return(results)
}
