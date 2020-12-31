#' @title Range Estimation
#' @description  Practical range estimation of birds using methods from Pennycuick (2008).
#'
#' @author Brian Masinde
#'
#' @name migrate
#' @param file Path to file where data resides.
#' @param header Logical. If TRUE use first row as column headers
#' @param sep separator
#' @param quote The set of quoting characters. see read.csv
#' @param dec The character used in the file for decimal points
#' @param fill See read.csv
#' @param comment.char For more details see read.csv
#' @param ... further arguments see read.csv
#' @param data A data frame with required columns: body mass, fat mass, etc.
#' @param settings A list for re-defining constants. See details for these with
#' default values from Pennycuick(2008) and Pennycuick(1998).
#' @param method Methods for protein energy consumption from muscle mass
#' @param speed_control One of two speed control methods. By default
#'        \emph{1} is used. \emph{0} is the alternative.
#'        The former holds the true airspeed constant while the latter holds the
#'        ratio of true airspeed to the minimum power speed constant (V:Vmp constant).
#' @param min_energy_protein Percentage of energy attributed to protein due to
#'        metabolism. Default value is 5 percent (0.05). If method  "csw" or "csp"
#'        is chosen, 2% -3% of the total energy is accounted for. The remainder
#'        would be attained from consuming protein in the airframe mass.
#'
#' @details The option *control takes the following arguments
#' \itemize{
#'    \item ppc: Profile power constant (8.4).
#'    \item fed: Energy content of fuel from fat (3.9E+07).
#'    \item ped: Energy content of protein (1.8E+07).
#'    \item g: Acceleration due to gravity (9.81).
#'    \item mce: Mechanical conversion efficiency [0,1]. Efficiency at which
#'    mechanical power is converted to chemical power (0.23).
#'    \item ipf: Induced power factor (1.2).
#'    \item vcp: Ventilation and circulation power (1.1).
#'    \item airDensity: Air density at cruising altitude (1.00).
#'    \item bdc: Body drag coefficient (0.1).
#'    \item alpha: Basal metabolism factors in passerines and non passerines (6.25, 3.79).
#'    \item delta: Basal metabolism factors in passerines and non passerines (0.724, 0.723)
#'    alpha*bodyMass^delta.
#'    \item mipd: Inverse power density of the mitochondria (1.2E-06).
#'    \item speedRatio: True air speed to minimum power speed ratio (1.2).
#'    \item muscDensity: Density of the flight muscles (1060).
#'    \item phr: Protein hydration ratio (2.2). Whenever protein is consumed from
#'    the muscle mass some amount of water is lost in the process. This water
#'    is estimated as the mass of dry protein time the phr.
#'}
#' @include control.R input_match.R misc_functions.R constant_muscle_mass.R
#' @return S3 class object with range estimates based on methods defined and
#'        settings
#' \itemize{
#'    \item range estimates (Km)
#'    \item remaining body mass (Kg)
#'    \item remaining fat mass (Kg)
#'    \item remaining muscle mass (Kg)
#'    \item minimum power speed at start of flight (m/s)
#'    \item minimum power speed at end of flight (m/s)
#'    \item taxon (order)
#' }
#'
#' @import Rcpp
#' @export migrate
#'
#' @examples
#' migrate(data = birds, settings = list(fed = 3.89*10^7))
#' migrate(data = birds,  method = "cmm", settings = list(airDensity = 0.905))
#'
#'
#' @usage migrate(file, header = TRUE, sep = ",", quote = "\"", dec = ".",
#'                fill = TRUE, comment.char = "", ...,
#'                data = NULL, settings = list(), method = "cmm",
#'                speed_control = 1, min_energy_protein = 0.05)
#'


migrate <- function(file, header = TRUE, sep = ",", quote = "\"", dec = ".",
                    fill = TRUE, comment.char = "", ...,
                    data = NULL, settings = list(), method = "cmm",
                    speed_control = 1, min_energy_protein = 0.05) {

  # both file and data not given throw an error
  if (missing(file) == TRUE & is.null(data) == TRUE) {
    stop("Data not found \n", call. = TRUE)
  }

  if (missing(file) == FALSE & is.null(data) == FALSE) {
    stop("Both path and data given. Function needs only one of the two \n", call. = TRUE)
  }

  if (speed_control != 1 && speed_control != 0) {
    stop("speed control should either be 1 or 0")
  }

  # min_protein < 0 or > 1 throw error
  if (min_energy_protein < 0 || min_energy_protein > 1) {
    stop("min_protein within [0,1] \n", call. = TRUE)
  }

  # process data ###############################################################
  if (missing(file) == FALSE) {
    dataRead <- read.csv(file = file, header = header, sep = sep, quote = quote,
                         dec = dec, fill = fill, comment.char,
                         stringsAsFactors = FALSE, ...)
    data <- .colnames.match(dataRead)
  } else {
    data <- .colnames.match(data)
  }

  # if any fat mass 0 zero stop (needs fuel)
  if (any(data$fatMass == 0.00)) {
    stop("Observations with zero fat mass. Birds need fat mass to migrate \n")
  }

  if (any(data$muscleMass == 0.00)) {
    stop("Observations with zero muscle mass \n")
  }

  if (any(data$bodyMass == 0.00)) {
    stop("Observations with zero body mass \n")
  }

  if (any(data$wingSpan == 0.00)) {
    stop("Observations with zero wingspan \n")
  }

  if (any(data$wingArea == 0.00)) {
    stop("Observations with zero wing area \n")
  }

  n <- nrow(data)
  # object with results ########################################################
  results <- list(
    range = vector(length = n),
    bodyMass = vector(length = n),
    fatMass = vector(length = n),
    muscleMass = vector(length = n),
    startMinSpeed = vector(length = n),
    endMinSpeed = vector(length = n),
    taxon = vector(length = n)
  )

  # control constants using default vs supplied
  if (missing(settings) == TRUE) {
    constants <- .control()
  } else {
    constants <- .control(settings)
  }

  if (method == "cmm") {
    simulation <-
      .constant.muscle.mass(data, constants, speed_control, min_energy_protein)
  } else if (method == "csw") {
    simulation <-
      .constant.specific.work(data, constants, speed_control, min_energy_protein)
  } else if (method == "csp") {
    simulation <-
      .constant.specific.power(data, constants, speed_control, min_energy_protein)
  }

  # aggregate dist from simulation to get range in Km
  results$range <- sapply(simulation$distance, function(x) sum(x) / 1000)
  results$bodyMass <- simulation$allUpMass
  results$fatMass <- round(simulation$fatMass,2)
  results$muscleMass <- simulation$muscleMass
  results$startMinSpeed <- simulation$startMinSpeed
  results$endMinSpeed <- simulation$endMinSpeed
  results$taxon <- data$taxon
  # range as named vectors
  if (!is.null(data$name)) {
    #names(results$range) <- as.vector(data$name)
    for (i in seq_len(length(simulation))) {
        names(results[[i]]) <- as.vector(data$name)
    }
  } else {
    #names(results$range) <- as.vector(data$ID)
    for (i in seq_len(length(simulation))) {
      names(results[[i]]) <- as.vector(data$ID)
    }
  }

  # return object of class migrate
  class(results) <- append(class(results), "migrate")
  return(results)
}


