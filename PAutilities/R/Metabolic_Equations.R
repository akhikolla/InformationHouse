#' Retrieve estimated basal metabolic rate for an individual
#'
#' @param Sex The individual's sex
#' @param Ht The individual's height, in meters
#' @param Wt The individual's weight, in kilograms
#' @param Age The individual's age, in years
#' @param verbose Logical. Should processing updates be printed?
#' @param RER numeric. The respiratory exchange ratio
#' @param equation The equation to apply
#' @param kcal_table The table to reference for converting kilocalories to
#'   oxygen consumption. See \code{\link{get_kcal_vo2_conversion}}
#' @param method The calculation method to use
#' @param MJ_conversion The value to use for converting megajoules to
#'   kilocalories. Defaults to thermochemical.
#' @param kcal_conversion numeric. If RER is NULL (default), the factor to use
#'   for converting kilocalories to oxygen consumption
#'
#' @references
#' Schofield, W. N. (1985). Predicting basal metabolic rate, new standards and
#' review of previous work. \emph{Human nutrition. Clinical nutrition}, 39,
#' 5-41.
#'
#' @export
#' @examples
#'
#' # Get BMR for an imaginary 900-year-old person (Age is only
#' # used to determine which equations to use, in this case the
#' # equations for someone older than 60)
#'
#' get_bmr(
#'   Sex = "M", Ht = 1.5, Wt = 80, Age = 900, equation = "both",
#'   method = "both", RER = 0.865, kcal_table = "both",
#'   MJ_conversion = c("all")
#' )
#'
#' get_bmr(
#'   Sex = "M", Ht = 1.5, Wt = 80, Age = 900, MJ_conversion = "all",
#'   kcal_conversion = 4.86
#' )
#'
#' get_bmr(
#'   Sex = "M", Ht = 1.5, Wt = 80, Age = 900, method = "FAO",
#'   kcal_conversion = 4.86
#' )
#'
get_bmr <- function(
  Sex = c("M", "F"), Ht = NULL, Wt, Age,
  verbose = FALSE, RER = NULL,
  equation = c("ht_wt", "wt", "both"),
  kcal_table = c("Lusk", "Peronnet", "both"),
  method = c("Schofield", "FAO", "both"),
  MJ_conversion = c("thermochemical", "dry", "convenience", "all"),
  kcal_conversion = 5
){

  if (verbose) cat(
    "\nCalculating Schofield predicted",
    "basal metabolic rate"
  )

  ## Set up arguments

    Sex <- match.arg(Sex, c("M", "F", "ERROR"))
    Sex <- switch(Sex, "M" = "male", "F" = "female")

    if (!is.null(Ht)) {
      if (Ht > 3){
        message("Detected height in cm. Converting to M.")
        Ht <- Ht / 100
      }
    }

    equation <- match.arg(
      equation, c("both", "ht_wt", "wt", "ERROR"), TRUE
    )
    if ("both" %in% equation) equation <- c("wt", "ht_wt")

    method <- match.arg(
      method, c("Schofield", "FAO", "both", "ERROR"), TRUE
    )
    if ("both" %in% method) method <- c("Schofield", "FAO")

    if (!is.null(RER)) {
      kcal_table <- match.arg(
        kcal_table, c("Lusk", "Peronnet", "both", "ERROR"), TRUE
      )
      if ("both" %in% kcal_table) kcal_table <- c("Lusk", "Peronnet")
      kcal_conversion <- get_kcal_vo2_conversion(RER, kcal_table)
    } else {
      kcal_table <- "NA"
    }

    MJ_conversion <- match.arg(
      MJ_conversion,
      c("thermochemical", "dry", "convenience", "all", "ERROR"),
      TRUE
    )

    if (identical(
      c("thermochemical", "dry", "convenience", "all"),
      MJ_conversion
    )) MJ_conversion <- "thermochemical"

    if ("all" %in% MJ_conversion) MJ_conversion <- c(
      "thermochemical", "dry", "convenience"
    )

    MJ_conversion <- sapply(
      MJ_conversion,
      function(x) {
        switch(
          x,
          "thermochemical" = 239.006,
          "dry" = 238.846,
          "convenience" = 239
        )
      }
    )

  ## Match participant to proper equation

    labels <- c(
      "less3", "3to10", "10to18",
      "18to30", "30to60", "over60"
    )

    agegroup <- as.character(cut(
      Age,
      c(-Inf, 3, 10, 18, 30, 60, Inf),
      labels,
      right = FALSE
    ))

    name_test <- paste("", Sex, agegroup, sep = "_")

    weights_schofield <- lapply(
      schofield_weights,
      function(x) {
        x[grepl(name_test, rownames(x)), ]
      }
    )

    weights_fao <- fao[rownames(weights_schofield[[1]]), ]

    all_weights <- c(
      weights_schofield,
      list(FAO = weights_fao)
    )

    stopifnot(all(sapply(all_weights, nrow) == 1))

  ## Determine and dispatch calculations

    calculations <- stats::setNames(
      expand.grid(
        equation, kcal_table, method,
        MJ_conversion, stringsAsFactors = FALSE
      ),
      c("equation", "kcal_table", "method", "MJ_conversion")
    )
    calculations$Wt <- Wt
    calculations$Ht <- Ht
    calculations$MJ_conversion_char <- names(MJ_conversion)[match(
      calculations$MJ_conversion,
      MJ_conversion
    )]

    if (!is.null(RER)) {
      calculations$kcal_conversion <- kcal_conversion[
        match(calculations$kcal_table, names(kcal_conversion))
      ]
    } else {
      calculations$kcal_conversion <- kcal_conversion
    }

    values <- do.call(
      rbind,
      lapply(
        split(calculations, seq(nrow(calculations))),
        metabolic_row_wise,
        weights = all_weights
      )
    )
    values <- values[!duplicated(values), ]

    if (verbose) cat(
      "... Done.\n"
    )

    return(values)

}

#' Calculate energy expenditure using the Weir equation
#'
#' @param VO2 Oxygen consumption
#' @param VCO2 Carbon dioxide production
#' @param epochSecs The averaging window of the metabolic data, in seconds
#'
#' @export
#'
#' @examples
#' weir_equation(3.5, 3.1, 60)
#'
weir_equation <- function(VO2, VCO2, epochSecs){
  kcal <- 1.44*(3.94*VO2 + 1.11*VCO2)
  kcal <- (kcal / (24*60*60))*epochSecs
  return(kcal)
}
