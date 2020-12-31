# function takes user settings and generates constants to be used for simulation
.control <- function(settings) {

  if (missing(settings) == FALSE &&
      is.list(settings) == FALSE) {
    stop("control must be a list")
  }

  # default constants
  constants <- list(
    # profile power constant
    ppc = 8.4,

    # energy content of fuel per kg
    fed = 3.9 * 10^7,

    # energy content of protein
    ped = 1.8 * 10^7,

    # acceleration due to gravity
    g = 9.81,

    # mechanical conversion efficiency  [0,1]
    mce = 0.23,

    # induced power factor
    ipf = 1.2,

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
    mipd = 1.2 * 10^-6,

    # ratio V:Vmp
    speedRatio = 1.2,

    # density of muscle
    muscDensity = 1060,

    # protein hydration ratio
    phr = 2.2
  )

  if (missing(settings) == TRUE) {
    # use team of default constants
    message("## settings not defined. Using default constants.
            \nDefault airDensity = 1.00 kg m^3 \n")
  }else{
    args <- c(
      "ppc",
      "fed",
      "ped",
      "g",
      "mce",
      "ipf",
      "vcp",
      "airDensity",
      "bdc",
      "alpha",
      "delta",
      "mipd",
      "speedRatio",
      "muscDensity",
      "phr"
    )

    # match args to user provided
    usrGiven <- which(args %in% names(settings) == TRUE)


    # extract names
    given <- args[usrGiven]
    for (i in seq_len(length(given))) {
      constants[given[i]] <- settings[given[i]]
    }

    # throw error wrong argument in settings
    if (length(constants) > 15) {
      stop("Wrong argument in settings", call. = FALSE)
    }

    if (length(constants$delta) != 2 || length(constants$alpha) != 2) {
      stop("In settings, alpha and delta as vectors of length == 2", call. = FALSE)
    }

  }
  # return object constants and class control
  class(constants) <- append(class(constants), "control")
  return(constants)
}

# print function for control
# print.control <- function(constants) {
#   cat("Constant    ",  "    Value", "\n")
#   for (i in seq_along(constants)) {
#     cat(names(constants)[i], "\t ", constants[[i]], "\n")
#   }
# }

