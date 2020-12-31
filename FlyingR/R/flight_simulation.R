#' @title Range Estimation
#' @description  Practical range estimation of birds using methods in Pennycuick (1975)
#' Mechanics of Flight. These methods are based on Breguet equations.
#'
#' @author Brian Masinde
#'
#' @name flysim
#' @param file Arguments for path to data.
#' @param header Logical. If TRUE use first row as column headers
#' @param sep separator
#' @param quote The set of quoting characters. see read.csv
#' @param dec The character used in the file for decimal points.
#' @param fill See read.csv
#' @param comment.char For more details see read.csv
#' @param ... further arguments see read.csv
#' @param data A data frame.
#' @param settings A list for re-defining constants. See details.
#'
#' @details The option *settings takes the arguments (those particularly
#' required by this function)
#' \itemize{
#'    \item ppc: Profile power constant
#'    \item fed: Energy content of fuel from fat
#'    \item g: Acceleration due to gravity
#'    \item mce: Mechanical conversion efficiency [0,1]
#'    \item ipf: Induced power factor
#'    \item vcp: Ventilation and circulation power
#'    \item airDensity: Air density at cruising altitude
#'    \item bdc: Body drag coefficient
#'    \item alpha: Basal metabolism factors in passerines and non passerines
#'    \item delta: Basal metabolism factors in passerines and non passerines
#'    alpha*bodyMass^delta
#'}
#' @include control.R misc_functions.R lookup_table2.R input_match.R method_1.R method_2.R
#' @return S3 class object with range estimates based on methods defined and
#'        settings used
#' \itemize{
#'    \item range estimates (Km)
#'    \item settings used
#'    \item data
#' }
#'
#' @importFrom  utils read.csv
#' @export flysim
#'
#' @examples
#' flysim(data = birds, settings = list(fatEnergy = 3.89*10^7))
#' flysim(data = birds,  settings = list(airDensity = 0.905))
#'
#' @usage flysim(file, header = TRUE, sep = ",", quote = "\"", dec = ".",
#'              fill = TRUE, comment.char = "", ..., data = NULL,
#'              settings = list())


flysim <- function(file, header = TRUE, sep = ",", quote = "\"", dec = ".",
                   fill = TRUE, comment.char = "", ..., data = NULL, settings = list()) {

  ##  Error check data #########################################################
  # missing file and data should throw an error
  if (missing(file) == TRUE & is.null(data) == TRUE) {
    stop("Data not found \n", call. = TRUE)
  }

  if (missing(file) == FALSE & is.null(data) == FALSE) {
    stop("Both path and data given. Function needs only one of the two \n", call. = TRUE)
  }

  # if (is.data.frame(data) == FALSE & is.list(data) == FALSE) {
  #   stop(">> data as list or data frame <<")
  # }
  #
  # # check number of columns
  # if (is.data.frame(data) == TRUE && ncol(data) < 6) {
  #   stop(">> at least 5 columns for data <<")
  # }
  # # check number of fields
  # if (is.list(data) == TRUE && length(data) < 6) {
  #   stop("data list should have at least 4 fields")
  # }

  # data processing
  if (missing(file) == FALSE) {
    # read data from file path
    dataRead <- read.csv(file = file, header = header, sep = sep,
                         quote = quote, dec = dec, fill = fill, comment.char,
                         stringsAsFactors = FALSE, ...)
    # column identification
    data <- .colnames.match(dataRead)
  } else if (missing(file) == TRUE) {
    data <- .colnames.match(data)
  }

  # control check
  if (missing(settings) == TRUE) {
    constants <- .control()
  } else {
    constants <- .control(settings)
  }

  results <- list("range" = vector(),
                  #"fuel" = vector(),
                  #"Vmp" = vector(),
                  #"Vmr" = vector(),
                  "constants" = unlist(constants, use.names = TRUE),
                  "data" = data
                  )

  results$range <- ifelse( data$taxon == 1, .breguet(bodyMass = data$allMass,
                                                    wingSpan = data$wingSpan,
                                                    fatMass = data$fatMass,
                                                    ordo = data$taxon,
                                                    wingArea = data$wingArea,
                                                    constants = constants),
      .breguet_adj(bodyMass = data$allMass,
                   wingSpan = data$wingSpan,
                   fatMass = data$fatMass,
                   ordo = data$taxon,
                   wingArea = data$wingArea,
                   constants = constants)
    )

  # results should be named vectors
  if (!is.null(data$name)) {
    names(results$range) <- as.vector(data$name)
  } else {
    names(results$range) <- as.vector(data$ID)
  }


  # return object of class flysim
  class(results) <- append(class(results), "flysim")
  # return class object
  return(results)
}



