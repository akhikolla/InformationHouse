# constant specific power time marching
# @author Brain Masinde.
# @name .constant.specific.power
# @param data Data as output from .colnames.match
# @param cons
# @param speed_control speed control as either

.constant.specific.power <- function(data, cons, speed_control, protein_met) {
  if (missing(data) == TRUE) {
    stop("Missing data argument", call. = FALSE)
  }

  if(missing(cons) == TRUE) {
    stop("Missing constants", call. = FALSE)
  }

  if(missing(speed_control) == TRUE) {
    stop("Missing speed control method", call. = FALSE)
  }

  # muscle mass is a must for this function
  if(is.null(data$muscleMass)) {
    stop("Muscle mass column missing", call.FALSE = TRUE)
  }

  all_mass <- data$allMass

  fat_mass <- data$fatMass

  wing_span <- data$wingSpan

  wing_area <- data$wingArea

  muscle_mass <- data$muscleMass

  id_col <- data$id

  name_col <- data$name

  taxon <- data$taxon

  # get fractions to get frame mass
  fat_frac <- fat_mass/all_mass

  muscle_frac <- muscle_mass/all_mass

  frame_mass <- all_mass*(1 - fat_frac - muscle_frac)

  # time marching

  sim_results <- list(
    bm = rep(list(vector()), nrow(data)),
    fm = rep(list(vector()), nrow(data)),
    mm = rep(list(vector()), nrow(data)),
    afm = rep(list(vector()), nrow(data)),
    mitochondria = rep(list(vector()), nrow(data)),
    myofibrils = rep(list(vector()), nrow(data)),
    dist = rep(list(vector()), nrow(data)),
    #delta_m = rep(list(vector()), nrow(data)),
    mechPow = rep(list(vector()), nrow(data)),
    chemPow = rep(list(vector()), nrow(data)),
    min_speed = rep(list(vector()), nrow(data)),
    true_speed = rep(list(vector()), nrow(data)),
    speed_ratio = rep(list(vector()), nrow(data)),
    wing_freq = rep(list(vector()), nrow(data)),
    spec_work = rep(list(vector()), nrow(data))
  )

  if(speed_control == "constant_speed") {
    for (i in 1:nrow(data)) {
      sim_results$bm[[i]][1] <- all_mass[i]
      sim_results$fm[[i]][1] <- fat_mass[i]
      sim_results$mm[[i]][1] <-  muscle_mass[i]
      sim_results$afm[[i]][1] <- all_mass[i] - (fat_mass[i] + muscle_mass[i])
      current_fm <- sim_results$fm[[i]][1]

      j <- 1
      while (sim_results$fm[[i]][j] > 0.000001) {
        # calculate speed
        sim_results$min_speed[[i]][j] <-
          .minpowspeed_cpp(
            bm = sim_results$bm[[i]][j],
            ws = wing_span[i],
            ipf = cons$ipf,
            g = cons$g,
            airDensity = cons$airDensity,
            bdc = cons$bdc
          )
        # true speed (constant speed always)
        sim_results$true_speed[[i]][j] <- sim_results$min_speed[[i]][1] * cons$speedRatio

        # mechanical power from power curve holding true air-speed constant
        sim_results$mechPow[[i]][j] <-
          .pow.curve(bm = sim_results$bm[[i]][j],
                     ws = wing_span[i],
                     wa = wing_area[i],
                     tas =  sim_results$true_speed[[i]][j], cons = cons)

        # wing frequency
        sim_results$wing_freq[[i]][j] <-
          .wingbeat.freq(bm = sim_results$bm[[i]][j],
                         ws = wing_span[i],
                         wa = wing_area[i],
                         cons)

        # subdivide muscle mass into mitochondria and myofibrils
        sim_results$myofibrils[[i]][1] <-
          sim_results$mm[[i]][1] * (1 - cons$invPower * (sim_results$mechPow[[i]][1] /
                                                           sim_results$mm[[i]][1]) * cons$muscDensity)
        sim_results$mitochondria[[i]][j] <-
          sim_results$mm[[i]][1] - sim_results$myofibrils[[i]][1]

        # mass specific power is at beginning of flight; mechanical power is divided
        # by the mass of the myofibrils (page 225 Pennycuick 2008)
        spec_pow_start <- sim_results$mechPow[[i]][1] / sim_results$myofibrils[[i]][1]

        if (protein_met == 0) {
          # chemical power
          deltaE <- (sim_results$mechPow[[i]][j] / cons$mce) +
            .basal.met(
              cons = cons,
              mFrame = frame_mass[i],
              mMusc = sim_results$mm[[i]][j],
              ordo = taxon[i]
            )

          # increase deltaE by 10% respiratory
          sim_results$chemPow[[i]][j] <- deltaE + (deltaE * 0.1)

          # fat consumed in the interval?
          used_fat <- sim_results$chemPow[[i]][j] / cons$eFat * 360

          #Reduce body composition by this consumed fat and determine what amount
          # of protein to consume to achieve specific power at start of flight
          dummy_bm <- sim_results$bm[[i]][j] - used_fat
          dummy_fm <- sim_results$fm[[i]][j] - used_fat

          # because of this reduction minimum speed, mechanical power, and
          # wing frequency reduce

          dummy_true_speed <-
            .minpowspeed_cpp(
              bm = dummy_bm,
              ws = wing_span[i],
              ipf = cons$ipf,
              g = cons$g,
              airDensity = cons$airDensity,
              bdc = cons$bdc
            ) * cons$speedRatio

          dummy_mechPow <- .pow.curve(bm = dummy_bm, ws = wing_span[i],
                                      wa = wing_area[i], tas = dummy_true_speed, cons = cons)

          # wing frequency
          dummy_wing_freq <-
            .wingbeat.freq(bm = dummy_bm, ws = wing_span[i], wa = wing_area[i], cons)

          #dummy_spec_pow as mechpow/mass myofibrils
          #dummy_spec_pow <- dummy_mechPow/sim_results$myofibrils[[i]][j]

          # amount of protein consumed that restores specific power to initial value
          used_protein <- - (dummy_mechPow/spec_pow_start) + sim_results$myofibrils[[i]][j]

          # checking if specific work has been restored
          sim_results$spec_work[[i]][j] <-
            dummy_mechPow / ((sim_results$myofibrils[[i]][j] - used_protein) * dummy_wing_freq)

          # amount of fule energy released is found by multiplying the mass of dry protein
          # removed by the energy density of dry protein
          used_fat_equiv <- (used_protein * cons$eProtein) / cons$eFat

          # adjust body components
          # new fat mass
          sim_results$fm[[i]][j+1] <-  sim_results$fm[[i]][j] - (used_fat - used_fat_equiv)
          sim_results$myofibrils[[i]][j+1] <- sim_results$myofibrils[[i]][j] - (used_protein * cons$phr)
          sim_results$mm[[i]][j+1] <- sim_results$mitochondria[[i]][1] + (sim_results$myofibrils[[i]][j] - (used_protein *cons$phr))
          sim_results$bm[[i]][j+1] <- sim_results$bm[[i]][j] - (used_fat - used_fat_equiv) - (used_protein *cons$phr)
          sim_results$dist[[i]][j] <- sim_results$true_speed[[i]][j] * 360
          j <-  j + 1
        }else {
          # chemical power
          deltaE <- (sim_results$mechPow[[i]][j] / cons$mce) +
            .basal.met(
              cons = cons,
              mFrame = frame_mass[i],
              mMusc = sim_results$mm[[i]][j],
              ordo = taxon[i]
            )

          # increase deltaE by 10% respiratory
          sim_results$chemPow[[i]][j] <- deltaE + (deltaE * 0.1)

          # protein requried for metabolic pathways
          EFromProtein <- sim_results$chemPow[[i]][j] * protein_met

          # fat consumed in the interval?
          used_fat <- (sim_results$chemPow[[i]][j] - EFromProtein)  / cons$eFat * 360

          #Reduce body composition by this consumed fat and determine what amount
          # of protein to consume to achieve specific power at start of flight
          dummy_fm <- sim_results$fm[[i]][j] - used_fat
          dummy_bm <-
            (sim_results$afm[[i]][j] - (EFromProtein / cons$eProtein)) + dummy_fm + sim_results$mm[[i]][j]


          # because of this reduction minimum speed, mechanical power, and
          # wing frequency reduce

          dummy_true_speed <-
            .minpowspeed_cpp(
              bm = dummy_bm,
              ws = wing_span[i],
              ipf = cons$ipf,
              g = cons$g,
              airDensity = cons$airDensity,
              bdc = cons$bdc
            ) * cons$speedRatio

          dummy_mechPow <- .pow.curve(bm = dummy_bm, ws = wing_span[i],
                                      wa = wing_area[i], tas = dummy_true_speed, cons = cons)

          # wing frequency
          dummy_wing_freq <-
            .wingbeat.freq(bm = dummy_bm, ws = wing_span[i], wa = wing_area[i], cons)

          #dummy_spec_pow as mechpow/mass myofibrils
          #dummy_spec_pow <- dummy_mechPow/sim_results$myofibrils[[i]][j]

          # amount of protein consumed that restores specific power to initial value
          used_protein <- - (dummy_mechPow/spec_pow_start) + sim_results$myofibrils[[i]][j]

          # checking if specific work has been restored
          sim_results$spec_work[[i]][j] <-
            dummy_mechPow / ((sim_results$myofibrils[[i]][j] - used_protein) * dummy_wing_freq)

          # amount of fule energy released is found by multiplying the mass of dry protein
          # removed by the energy density of dry protein
          used_fat_equiv <- (used_protein * cons$eProtein) / cons$eFat

          # adjust body components
          # new fat mass
          sim_results$fm[[i]][j+1] <-  sim_results$fm[[i]][j] - (used_fat - used_fat_equiv)
          sim_results$myofibrils[[i]][j+1] <- sim_results$myofibrils[[i]][j] - (used_protein * cons$phr)
          sim_results$mm[[i]][j+1] <- sim_results$mitochondria[[i]][1] + (sim_results$myofibrils[[i]][j] - (used_protein *cons$phr))
          sim_results$afm[[i]][j+1] <- sim_results$afm[[i]][j] - (EFromProtein / cons$eProtein)
          sim_results$bm[[i]][j + 1] <-
            sim_results$bm[[i]][j] - (used_fat - used_fat_equiv) - (used_protein *
                                                                      cons$phr) - (EFromProtein / cons$eProtein)
          sim_results$dist[[i]][j] <- sim_results$true_speed[[i]][j] * 360
          j <-  j + 1
        }
      }
    }
  } else {
    for (i in 1:nrow(data)) {
      sim_results$bm[[i]][1] <- all_mass[i]
      sim_results$fm[[i]][1] <- fat_mass[i]
      sim_results$mm[[i]][1] <-  muscle_mass[i]
      current_fm <- sim_results$fm[[i]][1]

      j <- 1
      while (sim_results$fm[[i]][j] > 0.000001) {
        # calculate speed
        sim_results$min_speed[[i]][j] <-
          .minpowspeed_cpp(
            bm = sim_results$bm[[i]][j],
            ws = wing_span[i],
            ipf = cons$ipf,
            g = cons$g,
            airDensity = cons$airDensity,
            bdc = cons$bdc
          )
        # hold speed ratio constant
        sim_results$true_speed[[i]][j] <- sim_results$min_speed[[i]][j] * cons$speedRatio

        # mechanical power from power curve holding true air-speed constant
        sim_results$mechPow[[i]][j] <-
          .pow.curve(bm = sim_results$bm[[i]][j],
                     ws = wing_span[i],
                     wa = wing_area[i],
                     tas =  sim_results$true_speed[[i]][j], cons = cons)

        # wing frequency
        sim_results$wing_freq[[i]][j] <-
          .wingbeat.freq(bm = sim_results$bm[[i]][j],
                         ws = wing_span[i],
                         wa = wing_area[i],
                         cons)

        # subdivide muscle mass into mitochondria and myofibrils
        sim_results$myofibrils[[i]][1] <-
          sim_results$mm[[i]][1] * (1 - cons$invPower * (sim_results$mechPow[[i]][1] /
                                                           sim_results$mm[[i]][1]) * cons$muscDensity)
        sim_results$mitochondria[[i]][j] <-
          sim_results$mm[[i]][1] - sim_results$myofibrils[[i]][1]

        # mass specific power is at beginning of flight; mechanical power is divided
        # by the mass of the myofibrils (page 225 Pennycuick 2008)
        spec_pow_start <- sim_results$mechPow[[i]][1] / sim_results$myofibrils[[i]][1]

        # chemical power
        deltaE <- (sim_results$mechPow[[i]][j] / cons$mce) +
          .basal.met(
            cons = cons,
            mFrame = frame_mass[i],
            mMusc = sim_results$mm[[i]][j],
            ordo = taxon[i]
          )

        # increase deltaE by 10% respiratory
        sim_results$chemPow[[i]][j] <- deltaE + (deltaE * 0.1)

        # fat consumed in the interval?
        used_fat <- sim_results$chemPow[[i]][j] / cons$eFat * 360

        #Reduce body composition by this consumed fat and determine what amount
        # of protein to consume to achieve specific power at start of flight
        dummy_bm <- sim_results$bm[[i]][j] - used_fat
        dummy_fm <- sim_results$fm[[i]][j] - used_fat

        # because of this reduction minimum speed, mechanical power, and
        # wing frequency reduce

        dummy_true_speed <-
          .minpowspeed_cpp(
            bm = dummy_bm,
            ws = wing_span[i],
            ipf = cons$ipf,
            g = cons$g,
            airDensity = cons$airDensity,
            bdc = cons$bdc
          ) * cons$speedRatio

        dummy_mechPow <- .pow.curve(bm = dummy_bm, ws = wing_span[i],
                                    wa = wing_area[i], tas = dummy_true_speed, cons = cons)

        # wing frequency
        dummy_wing_freq <-
          .wingbeat.freq(bm = dummy_bm, ws = wing_span[i], wa = wing_area[i], cons)

        #dummy_spec_pow as mechpow/mass myofibrils
        #dummy_spec_pow <- dummy_mechPow/sim_results$myofibrils[[i]][j]

        # amount of protein consumed that restores specific power to initial value
        used_protein <- - (dummy_mechPow/spec_pow_start) + sim_results$myofibrils[[i]][j]

        # checking if specific work has been restored
        sim_results$spec_work[[i]][j] <-
          dummy_mechPow / ((sim_results$myofibrils[[i]][j] - used_protein) * dummy_wing_freq)

        # amount of fule energy released is found by multiplying the mass of dry protein
        # removed by the energy density of dry protein
        used_fat_equiv <- (used_protein * cons$eProtein) / cons$eFat

        # adjust body components
        # new fat mass
        sim_results$fm[[i]][j+1] <-  sim_results$fm[[i]][j] - (used_fat - used_fat_equiv)
        sim_results$myofibrils[[i]][j+1] <- sim_results$myofibrils[[i]][j] - (used_protein * cons$phr)
        sim_results$mm[[i]][j+1] <- sim_results$mitochondria[[i]][1] + (sim_results$myofibrils[[i]][j] - (used_protein *cons$phr))
        sim_results$bm[[i]][j+1] <- sim_results$bm[[i]][j] - (used_fat - used_fat_equiv) - (used_protein * cons$phr)
        sim_results$dist[[i]][j] <- sim_results$true_speed[[i]][j] * 360
        j <-  j + 1

      }
    }
  }
  return(sim_results)
}
