# function takes user settings and generates constants to be used for simulation
.control <- function(settings) {

  if (missing(settings) == FALSE &&
      is.list(settings) == FALSE) {
    stop("contorl must be a list")
  }

  cons <- list(
    # profile power constant
    ppc = 8.4,

    # energy content of fuel per kg
    eFat = 3.9 * 10^7,

    # energy content of
    eProtein = 1.8 * 10^7,

    # accelaration due to gravity
    g = 9.81,

    # mechanical conversion efficiency  [0,1]
    mce = 0.23,

    # induced power factor
    ipf = 1.20,

    # ventilation and circulation power (Tucker's data)
    vcp =  1.10,

    # air density at flight height
    airDensity = 1.00,

    # body drag coefficient
    bdc = 0.10,

    # constant varies btw passerines and non-passerines
    alpha = c(6.25, 3.79),

    delta = c(0.724, 0.723),

    # inverse power density of mitochondria
    invPower= 1.2 * 10 ^-6,

    # ratio V:Vmp
    speedRatio = 1.2,

    # density of muscle
    muscDensity = 1060,

    # protein hydration ratio
    phr = 2.2
  )

  if (missing(settings) == TRUE) {
    # use team of default parameters
    message("## settings not defined. Using default constants.
            \nDefault airDensity = 1.00 kg m^3 \n")
  }else{
    extArgs <- c(
      "ppc",
      "eFat",
      "eProtein",
      "g",
      "mce",
      "ipf",
      "vcp",
      "airDensity",
      "bdc",
      "alpha",
      "delta",
      "invPower",
      "speedRatio",
      "muscDensity",
      "phr"
    )

    # match extArgs to user provided
    given <- which(extArgs %in% names(settings) == TRUE)


    # extract names
    consGive <- extArgs[given]
    for (i in 1:length(consGive)) {
      cons[consGive[i]] <- settings[consGive[i]]
    }

    # throw error wrong argument in settings
    if(length(cons) > 15){
      stop("Wrong argument in settings", call. = FALSE)
    }

    if(length(cons$delta) != 2 || length(cons$alpha) != 2) {
      stop("In settings, alpha and delta as vectors of length == 2", call. = FALSE)
    }

  }

  return(cons)
}
