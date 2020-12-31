
#' Generate a matrix of two-phase relative permeability data for the water-oil system using the modified Brooks-Corey model
#'
#' The 'kr2p_ow()' creates a table of two-phase water and oil relative permeability data for water and oil saturation values between zero and one.
#' @param SWCON connate water saturation, fraction
#' @param SWCRIT critical water saturation, fraction
#' @param SOIRW irreducible oil saturation, fraction
#' @param SORW residual oil saturation, fraction
#' @param KRWIRO water relative permeability at irreducible oil
#' @param KROCW oil relative permeability at connate water
#' @param NW exponent term for calculating krw
#' @param NOW exponent term for calculating krow
#' @param NP number of saturation points in the table, the maximum acceptable value is 501
#' @return A matrix with water saturation, oil saturation, water relative permeability, and oil relative permeability values, respectively
#' @examples
#' rel_perm_wo <- kr2p_ow(0.15, 0.2, 0.15, 0.15, 0.4, 1, 3, 2, 101)
#'
#' @references
#' \insertRef{Brooks1964}{Rrelperm}
#'
#' @export
kr2p_ow <- function(SWCON, SWCRIT, SOIRW, SORW, KRWIRO, KROCW, NW, NOW, NP) {

   if (is.null(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (SWCON < 0 | SWCON > 1) stop("'SWCON' must be a numeric between zero and one.")

   if (is.null(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (SWCRIT < 0 | SWCRIT > 1) stop("'SWCRIT' must be a numeric between zero and one.")
   if (SWCRIT < SWCON) stop("'SWCRIT' must be equal to or greater than 'SWCON'.")

   if (is.null(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (SOIRW < 0 | SOIRW > 1) stop("'SOIRW' must be a numeric between zero and one.")

   if (is.null(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (SORW < 0 | SORW > 1) stop("'SORW' must be a numeric between zero and one.")
   if (SORW < SOIRW) stop("'SORW' must be equal to or greater than 'SOIRW'.")

   if (is.null(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (!is.numeric(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (KRWIRO < 0 | KRWIRO > 1) stop("'KRWIRO' must be a numeric between zero and one.")

   if (is.null(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (!is.numeric(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (KROCW < 0 | KROCW > 1) stop("'KROCW' must be a numeric between zero and one.")

   if (is.null(NW)) stop("'NW' must be a numeric between zero and one.")
   if (!is.numeric(NW)) stop("'NW' must be a numeric between zero and one.")
   if (NW < 0) stop("'NW' must be a numeric greater than zero.")

   if (is.null(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (!is.numeric(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (NOW < 0) stop("'NOW' must be a numeric greater than zero.")

   if (is.null(NP)) stop("'NP' must be a numeric between zero and one.")
   if (!is.numeric(NP)) stop("'NP' must be a numeric between zero and one.")
   if (NP < 1 | NP > 501) stop("'NP' must be a numeric between one and 500.")

   results <- kr2p_ow_cpp(SWCON, SWCRIT, SOIRW, SORW, KRWIRO, KROCW, NW, NOW, NP)
   return(results)
}


#' Generate a matrix of two-phase relative permeability data for the gas-liquid system using the modified Brooks-Corey model
#'
#' The 'kr2p_gl()' creates a table of two-phase gas and liquid relative permeability data for gas and liquid saturation values between zero and one.
#' @param SWCON connate water saturation, fraction
#' @param SOIRG irreducible oil saturation, fraction
#' @param SORG residual oil saturation, fraction
#' @param SGCON connate gas saturation, fraction
#' @param SGCRIT critical gas saturation, fraction
#' @param KRGCL gas relative permeability at connate liquid
#' @param KROGCG oil relative permeability at connate gas
#' @param NG exponent term for calculating krg
#' @param NOG exponent term for calculating krog
#' @param NP number of saturation points in the table, the maximum acceptable value is 501
#' @return A matrix with gas saturation, liquid saturation, gas relative permeability, and oil relative permeability values, respectively
#' @examples
#' rel_perm_gl <- kr2p_gl(0.15, 0.1, 0.1, 0.05, 0.05, 0.3, 1, 4, 2.25, 101)
#'
#' @references
#' \insertRef{Brooks1964}{Rrelperm}
#'
#' @export
kr2p_gl <- function(SWCON, SOIRG, SORG, SGCON, SGCRIT, KRGCL, KROGCG, NG, NOG, NP) {

   if (is.null(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (SWCON < 0 | SWCON > 1) stop("'SWCON' must be a numeric between zero and one.")

   if (is.null(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (SOIRG < 0 | SOIRG > 1) stop("'SOIRG' must be a numeric between zero and one.")

   if (is.null(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (SORG < 0 | SORG > 1) stop("'SORG' must be a numeric between zero and one.")
   if (SORG < SOIRG) stop("'SORG' must be equal to or greater than 'SOIRG'.")

   if (is.null(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (SGCON < 0 | SGCON > 1) stop("'SGCON' must be a numeric between zero and one.")

   if (is.null(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (SGCRIT < 0 | SGCRIT > 1) stop("'SGCRIT' must be a numeric between zero and one.")
   if (SGCRIT < SGCON) stop("'SGCRIT' must be equal to or greater than 'SGCON'.")

   if (is.null(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (!is.numeric(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (KRGCL < 0 | KRGCL > 1) stop("'KRGCL' must be a numeric between zero and one.")

   if (is.null(KROGCG)) stop("'KROGCG' must be a numeric between zero and one.")
   if (!is.numeric(KROGCG)) stop("'KROGCG' must be a numeric between zero and one.")
   if (KROGCG < 0 | KROGCG > 1) stop("'KROGCG' must be a numeric between zero and one.")

   if (is.null(NG)) stop("'NG' must be a numeric between zero and one.")
   if (!is.numeric(NG)) stop("'NG' must be a numeric between zero and one.")
   if (NG < 0) stop("'NG' must be a numeric greater than zero.")

   if (is.null(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (!is.numeric(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (NOG < 0) stop("'NOG' must be a numeric greater than zero.")

   if (is.null(NP)) stop("'NP' must be a numeric between zero and one.")
   if (!is.numeric(NP)) stop("'NP' must be a numeric between zero and one.")
   if (NP < 1 | NP > 501) stop("'NP' must be a numeric between one and 500.")

   results <- kr2p_gl_cpp(SWCON, SOIRG, SORG, SGCON, SGCRIT, KRGCL, KROGCG, NG,
                          NOG, NP)
   return(results)
}


#' Generate a matrix of three-phase relative permeability data for the water-gas-oil system using the modified Stone I model
#'
#' The 'kr3p_StoneI_So()' creates a table of three-phase oil relative permeability data for water, gas, and oil saturation values between zero and one. This model reads the oil saturation in the three-phase region as oil saturation input into two-phase relative permeability models
#' @param SWCON connate water saturation, fraction
#' @param SWCRIT critical water saturation, fraction
#' @param SOIRW irreducible oil saturation, fraction
#' @param SORW residual oil saturation, fraction
#' @param SOIRG irreducible oil saturation, fraction
#' @param SORG residual oil saturation, fraction
#' @param SGCON connate gas saturation, fraction
#' @param SGCRIT critical gas saturation, fraction
#' @param KRWIRO water relative permeability at irreducible oil
#' @param KROCW oil relative permeability at connate water
#' @param KRGCL gas relative permeability at connate liquid
#' @param NW exponent term for calculating krw
#' @param NOW exponent term for calculating krow
#' @param NG exponent term for calculating krg
#' @param NOG exponent term for calculating krog
#' @param NP number of saturation points in the two-phase relative permeability tables, the maximum acceptable value is 501. The number of data points in the three-phase relative permeability table is (0.5 * NP * (NP + 1))
#' @return A matrix with water saturation, gas saturation, oil saturation, and oil relative permeability values, respectively.
#' @examples
#' rel_perm_wgo <- kr3p_StoneI_So(0.15, 0.2, 0.15, 0.15, 0.2, 0.2, 0.05, 0.05,
#' 0.4, 1, 0.3, 3, 2, 4, 2.5, 101)
#'
#' @references
#' \insertRef{Stone1970}{Rrelperm}
#'
#' \insertRef{Fayers1984}{Rrelperm}
#'
#' @export

kr3p_StoneI_So <- function(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON, SGCRIT,
                           KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP) {

   if (is.null(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (SWCON < 0 | SWCON > 1) stop("'SWCON' must be a numeric between zero and one.")

   if (is.null(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (SWCRIT < 0 | SWCRIT > 1) stop("'SWCRIT' must be a numeric between zero and one.")
   if (SWCRIT < SWCON) stop("'SWCRIT' must be equal to or greater than 'SWCON'.")

   if (is.null(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (SOIRW < 0 | SOIRW > 1) stop("'SOIRW' must be a numeric between zero and one.")

   if (is.null(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (SORW < 0 | SORW > 1) stop("'SORW' must be a numeric between zero and one.")
   if (SORW < SOIRW) stop("'SORW' must be equal to or greater than 'SOIRW'.")

   if (is.null(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (SOIRG < 0 | SOIRG > 1) stop("'SOIRG' must be a numeric between zero and one.")

   if (is.null(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (SORG < 0 | SORG > 1) stop("'SORG' must be a numeric between zero and one.")
   if (SORG < SOIRG) stop("'SORG' must be equal to or greater than 'SOIRG'.")

   if (is.null(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (SGCON < 0 | SGCON > 1) stop("'SGCON' must be a numeric between zero and one.")

   if (is.null(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (SGCRIT < 0 | SGCRIT > 1) stop("'SGCRIT' must be a numeric between zero and one.")
   if (SGCRIT < SGCON) stop("'SGCRIT' must be equal to or greater than 'SGCON'.")

   if (is.null(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (!is.numeric(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (KRWIRO < 0 | KRWIRO > 1) stop("'KRWIRO' must be a numeric between zero and one.")

   if (is.null(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (!is.numeric(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (KROCW < 0 | KROCW > 1) stop("'KROCW' must be a numeric between zero and one.")

   if (is.null(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (!is.numeric(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (KRGCL < 0 | KRGCL > 1) stop("'KRGCL' must be a numeric between zero and one.")

   if (is.null(NW)) stop("'NW' must be a numeric between zero and one.")
   if (!is.numeric(NW)) stop("'NW' must be a numeric between zero and one.")
   if (NW < 0) stop("'NW' must be a numeric greater than zero.")

   if (is.null(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (!is.numeric(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (NOW < 0) stop("'NOW' must be a numeric greater than zero.")

   if (is.null(NG)) stop("'NG' must be a numeric between zero and one.")
   if (!is.numeric(NG)) stop("'NG' must be a numeric between zero and one.")
   if (NG < 0) stop("'NG' must be a numeric greater than zero.")

   if (is.null(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (!is.numeric(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (NOG < 0) stop("'NOG' must be a numeric greater than zero.")

   if (is.null(NP)) stop("'NP' must be a numeric between zero and one.")
   if (!is.numeric(NP)) stop("'NP' must be a numeric between zero and one.")
   if (NP < 1 | NP > 501) stop("'NP' must be a numeric between one and 500.")

   results <- kr3p_StoneI_So_cpp(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON,
                                 SGCRIT, KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP)
   return(results)
}



#' Generate a matrix of three-phase relative permeability data for the water-gas-oil system using the modified Stone I model
#'
#' The 'kr3p_StoneI_SwSg()' creates a table of three-phase oil relative permeability data for water, gas, and oil saturation values between zero and one. This model reads the water, and gas saturation values in the three-phase region as water, and gas saturation inputs into two-phase relative permeability models
#' @param SWCON connate water saturation, fraction
#' @param SWCRIT critical water saturation, fraction
#' @param SOIRW irreducible oil saturation, fraction
#' @param SORW residual oil saturation, fraction
#' @param SOIRG irreducible oil saturation, fraction
#' @param SORG residual oil saturation, fraction
#' @param SGCON connate gas saturation, fraction
#' @param SGCRIT critical gas saturation, fraction
#' @param KRWIRO water relative permeability at irreducible oil
#' @param KROCW oil relative permeability at connate water
#' @param KRGCL gas relative permeability at connate liquid
#' @param NW exponent term for calculating krw
#' @param NOW exponent term for calculating krow
#' @param NG exponent term for calculating krg
#' @param NOG exponent term for calculating krog
#' @param NP number of saturation points in the two-phase relative permeability tables, the maximum acceptable value is 501. The number of data points in the three-phase relative permeability table is (0.5 * NP * (NP + 1))
#' @return A matrix with water saturation, gas saturation, oil saturation, and oil relative permeability values, respectively.
#' @examples
#' rel_perm_wgo <- kr3p_StoneI_SwSg(0.15, 0.2, 0.15, 0.15, 0.2, 0.2, 0.05, 0.05,
#' 0.4, 1, 0.3, 3, 2, 4, 2.5, 101)
#'
#' @references
#' \insertRef{Stone1970}{Rrelperm}
#'
#' \insertRef{Fayers1984}{Rrelperm}
#'
#' @export
kr3p_StoneI_SwSg <- function(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON, SGCRIT,
                             KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP) {

   if (is.null(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (SWCON < 0 | SWCON > 1) stop("'SWCON' must be a numeric between zero and one.")

   if (is.null(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (SWCRIT < 0 | SWCRIT > 1) stop("'SWCRIT' must be a numeric between zero and one.")
   if (SWCRIT < SWCON) stop("'SWCRIT' must be equal to or greater than 'SWCON'.")

   if (is.null(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (SOIRW < 0 | SOIRW > 1) stop("'SOIRW' must be a numeric between zero and one.")

   if (is.null(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (SORW < 0 | SORW > 1) stop("'SORW' must be a numeric between zero and one.")
   if (SORW < SOIRW) stop("'SORW' must be equal to or greater than 'SOIRW'.")

   if (is.null(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (SOIRG < 0 | SOIRG > 1) stop("'SOIRG' must be a numeric between zero and one.")

   if (is.null(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (SORG < 0 | SORG > 1) stop("'SORG' must be a numeric between zero and one.")
   if (SORG < SOIRG) stop("'SORG' must be equal to or greater than 'SOIRG'.")

   if (is.null(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (SGCON < 0 | SGCON > 1) stop("'SGCON' must be a numeric between zero and one.")

   if (is.null(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (SGCRIT < 0 | SGCRIT > 1) stop("'SGCRIT' must be a numeric between zero and one.")
   if (SGCRIT < SGCON) stop("'SGCRIT' must be equal to or greater than 'SGCON'.")

   if (is.null(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (!is.numeric(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (KRWIRO < 0 | KRWIRO > 1) stop("'KRWIRO' must be a numeric between zero and one.")

   if (is.null(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (!is.numeric(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (KROCW < 0 | KROCW > 1) stop("'KROCW' must be a numeric between zero and one.")

   if (is.null(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (!is.numeric(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (KRGCL < 0 | KRGCL > 1) stop("'KRGCL' must be a numeric between zero and one.")

   if (is.null(NW)) stop("'NW' must be a numeric between zero and one.")
   if (!is.numeric(NW)) stop("'NW' must be a numeric between zero and one.")
   if (NW < 0) stop("'NW' must be a numeric greater than zero.")

   if (is.null(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (!is.numeric(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (NOW < 0) stop("'NOW' must be a numeric greater than zero.")

   if (is.null(NG)) stop("'NG' must be a numeric between zero and one.")
   if (!is.numeric(NG)) stop("'NG' must be a numeric between zero and one.")
   if (NG < 0) stop("'NG' must be a numeric greater than zero.")

   if (is.null(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (!is.numeric(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (NOG < 0) stop("'NOG' must be a numeric greater than zero.")

   if (is.null(NP)) stop("'NP' must be a numeric between zero and one.")
   if (!is.numeric(NP)) stop("'NP' must be a numeric between zero and one.")
   if (NP < 1 | NP > 501) stop("'NP' must be a numeric between one and 500.")

   results <- kr3p_StoneI_SwSg_cpp(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON,
                                   SGCRIT, KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP)
   return(results)
}



#' Generate a matrix of three-phase relative permeability data for the water-gas-oil system using the modified Stone II model
#'
#' The 'kr3p_StoneII_So()' creates a table of three-phase oil relative permeability data for water, gas, and oil saturation values between zero and one. This model reads the oil saturation in the three-phase region as oil saturation input into two-phase relative permeability models
#' @param SWCON connate water saturation, fraction
#' @param SWCRIT critical water saturation, fraction
#' @param SOIRW irreducible oil saturation, fraction
#' @param SORW residual oil saturation, fraction
#' @param SOIRG irreducible oil saturation, fraction
#' @param SORG residual oil saturation, fraction
#' @param SGCON connate gas saturation, fraction
#' @param SGCRIT critical gas saturation, fraction
#' @param KRWIRO water relative permeability at irreducible oil
#' @param KROCW oil relative permeability at connate water
#' @param KRGCL gas relative permeability at connate liquid
#' @param NW exponent term for calculating krw
#' @param NOW exponent term for calculating krow
#' @param NG exponent term for calculating krg
#' @param NOG exponent term for calculating krog
#' @param NP number of saturation points in the two-phase relative permeability tables, the maximum acceptable value is 501. The number of data points in the three-phase relative permeability table is (0.5 * NP * (NP + 1))
#' @return A matrix with water saturation, gas saturation, oil saturation, and oil relative permeability values, respectively.
#' @examples
#' rel_perm_wgo <- kr3p_StoneII_So(0.15, 0.2, 0.15, 0.15, 0.2, 0.2, 0.05, 0.05,
#' 0.4, 1, 0.3, 3, 2, 4, 2.5, 101)
#'
#' @references
#' \insertRef{Stone1970}{Rrelperm}
#'
#' \insertRef{Fayers1984}{Rrelperm}
#'
#' @export
kr3p_StoneII_So <- function(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON, SGCRIT,
                            KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP) {

   if (is.null(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (SWCON < 0 | SWCON > 1) stop("'SWCON' must be a numeric between zero and one.")

   if (is.null(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (SWCRIT < 0 | SWCRIT > 1) stop("'SWCRIT' must be a numeric between zero and one.")
   if (SWCRIT < SWCON) stop("'SWCRIT' must be equal to or greater than 'SWCON'.")

   if (is.null(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (SOIRW < 0 | SOIRW > 1) stop("'SOIRW' must be a numeric between zero and one.")

   if (is.null(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (SORW < 0 | SORW > 1) stop("'SORW' must be a numeric between zero and one.")
   if (SORW < SOIRW) stop("'SORW' must be equal to or greater than 'SOIRW'.")

   if (is.null(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (SOIRG < 0 | SOIRG > 1) stop("'SOIRG' must be a numeric between zero and one.")

   if (is.null(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (SORG < 0 | SORG > 1) stop("'SORG' must be a numeric between zero and one.")
   if (SORG < SOIRG) stop("'SORG' must be equal to or greater than 'SOIRG'.")

   if (is.null(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (SGCON < 0 | SGCON > 1) stop("'SGCON' must be a numeric between zero and one.")

   if (is.null(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (SGCRIT < 0 | SGCRIT > 1) stop("'SGCRIT' must be a numeric between zero and one.")
   if (SGCRIT < SGCON) stop("'SGCRIT' must be equal to or greater than 'SGCON'.")

   if (is.null(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (!is.numeric(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (KRWIRO < 0 | KRWIRO > 1) stop("'KRWIRO' must be a numeric between zero and one.")

   if (is.null(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (!is.numeric(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (KROCW < 0 | KROCW > 1) stop("'KROCW' must be a numeric between zero and one.")

   if (is.null(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (!is.numeric(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (KRGCL < 0 | KRGCL > 1) stop("'KRGCL' must be a numeric between zero and one.")

   if (is.null(NW)) stop("'NW' must be a numeric between zero and one.")
   if (!is.numeric(NW)) stop("'NW' must be a numeric between zero and one.")
   if (NW < 0) stop("'NW' must be a numeric greater than zero.")

   if (is.null(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (!is.numeric(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (NOW < 0) stop("'NOW' must be a numeric greater than zero.")

   if (is.null(NG)) stop("'NG' must be a numeric between zero and one.")
   if (!is.numeric(NG)) stop("'NG' must be a numeric between zero and one.")
   if (NG < 0) stop("'NG' must be a numeric greater than zero.")

   if (is.null(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (!is.numeric(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (NOG < 0) stop("'NOG' must be a numeric greater than zero.")

   if (is.null(NP)) stop("'NP' must be a numeric between zero and one.")
   if (!is.numeric(NP)) stop("'NP' must be a numeric between zero and one.")
   if (NP < 1 | NP > 501) stop("'NP' must be a numeric between one and 500.")

   results <- kr3p_StoneII_So_cpp(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON,
                                  SGCRIT, KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP)
   return(results)
}




#' Generate a matrix of three-phase relative permeability data for the water-gas-oil system using the modified Stone II model
#'
#' The 'kr3p_StoneII_SwSg()' creates a table of three-phase oil relative permeability data for water, gas, and oil saturation values between zero and one. This model reads the water, and gas saturation values in the three-phase region as water, and gas saturation inputs into two-phase relative permeability models
#' @param SWCON connate water saturation, fraction
#' @param SWCRIT critical water saturation, fraction
#' @param SOIRW irreducible oil saturation, fraction
#' @param SORW residual oil saturation, fraction
#' @param SOIRG irreducible oil saturation, fraction
#' @param SORG residual oil saturation, fraction
#' @param SGCON connate gas saturation, fraction
#' @param SGCRIT critical gas saturation, fraction
#' @param KRWIRO water relative permeability at irreducible oil
#' @param KROCW oil relative permeability at connate water
#' @param KRGCL gas relative permeability at connate liquid
#' @param NW exponent term for calculating krw
#' @param NOW exponent term for calculating krow
#' @param NG exponent term for calculating krg
#' @param NOG exponent term for calculating krog
#' @param NP number of saturation points in the two-phase relative permeability tables, the maximum acceptable value is 501. The number of data points in the three-phase relative permeability table is (0.5 * NP * (NP + 1))
#' @return A matrix with water saturation, gas saturation, oil saturation, and oil relative permeability values, respectively.
#' @examples
#' rel_perm_wgo <- kr3p_StoneII_SwSg(0.15, 0.2, 0.15, 0.15, 0.2, 0.2, 0.05, 0.05,
#' 0.4, 1, 0.3, 3, 2, 4, 2.5, 101)
#'
#' @references
#' \insertRef{Stone1970}{Rrelperm}
#'
#' \insertRef{Fayers1984}{Rrelperm}
#'
#' @export
kr3p_StoneII_SwSg <- function(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON, SGCRIT,
                              KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP) {

   if (is.null(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (SWCON < 0 | SWCON > 1) stop("'SWCON' must be a numeric between zero and one.")

   if (is.null(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (SWCRIT < 0 | SWCRIT > 1) stop("'SWCRIT' must be a numeric between zero and one.")
   if (SWCRIT < SWCON) stop("'SWCRIT' must be equal to or greater than 'SWCON'.")

   if (is.null(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (SOIRW < 0 | SOIRW > 1) stop("'SOIRW' must be a numeric between zero and one.")

   if (is.null(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (SORW < 0 | SORW > 1) stop("'SORW' must be a numeric between zero and one.")
   if (SORW < SOIRW) stop("'SORW' must be equal to or greater than 'SOIRW'.")

   if (is.null(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (SOIRG < 0 | SOIRG > 1) stop("'SOIRG' must be a numeric between zero and one.")

   if (is.null(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (SORG < 0 | SORG > 1) stop("'SORG' must be a numeric between zero and one.")
   if (SORG < SOIRG) stop("'SORG' must be equal to or greater than 'SOIRG'.")

   if (is.null(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (SGCON < 0 | SGCON > 1) stop("'SGCON' must be a numeric between zero and one.")

   if (is.null(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (SGCRIT < 0 | SGCRIT > 1) stop("'SGCRIT' must be a numeric between zero and one.")
   if (SGCRIT < SGCON) stop("'SGCRIT' must be equal to or greater than 'SGCON'.")

   if (is.null(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (!is.numeric(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (KRWIRO < 0 | KRWIRO > 1) stop("'KRWIRO' must be a numeric between zero and one.")

   if (is.null(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (!is.numeric(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (KROCW < 0 | KROCW > 1) stop("'KROCW' must be a numeric between zero and one.")

   if (is.null(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (!is.numeric(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (KRGCL < 0 | KRGCL > 1) stop("'KRGCL' must be a numeric between zero and one.")

   if (is.null(NW)) stop("'NW' must be a numeric between zero and one.")
   if (!is.numeric(NW)) stop("'NW' must be a numeric between zero and one.")
   if (NW < 0) stop("'NW' must be a numeric greater than zero.")

   if (is.null(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (!is.numeric(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (NOW < 0) stop("'NOW' must be a numeric greater than zero.")

   if (is.null(NG)) stop("'NG' must be a numeric between zero and one.")
   if (!is.numeric(NG)) stop("'NG' must be a numeric between zero and one.")
   if (NG < 0) stop("'NG' must be a numeric greater than zero.")

   if (is.null(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (!is.numeric(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (NOG < 0) stop("'NOG' must be a numeric greater than zero.")

   if (is.null(NP)) stop("'NP' must be a numeric between zero and one.")
   if (!is.numeric(NP)) stop("'NP' must be a numeric between zero and one.")
   if (NP < 1 | NP > 501) stop("'NP' must be a numeric between one and 500.")

   results <- kr3p_StoneII_SwSg_cpp(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON,
                                    SGCRIT, KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP)
   return(results)
}


#' Generate a matrix of three-phase relative permeability data for the water-gas-oil system using Baker's linear model
#'
#' The 'kr3p_Baker()' creates a table of three-phase oil relative permeability data for water, gas, and oil saturation values between zero and one. This model reads the water, and gas saturation values in the three-phase region as water, and gas saturation inputs into two-phase relative permeability models
#' @param SWCON connate water saturation, fraction
#' @param SWCRIT critical water saturation, fraction
#' @param SOIRW irreducible oil saturation, fraction
#' @param SORW residual oil saturation, fraction
#' @param SOIRG irreducible oil saturation, fraction
#' @param SORG residual oil saturation, fraction
#' @param SGCON connate gas saturation, fraction
#' @param SGCRIT critical gas saturation, fraction
#' @param KRWIRO water relative permeability at irreducible oil
#' @param KROCW oil relative permeability at connate water
#' @param KRGCL gas relative permeability at connate liquid
#' @param NW exponent term for calculating krw
#' @param NOW exponent term for calculating krow
#' @param NG exponent term for calculating krg
#' @param NOG exponent term for calculating krog
#' @param NP number of saturation points in the two-phase relative permeability tables, the maximum acceptable value is 501. The number of data points in the three-phase relative permeability table is (0.5 * NP * (NP + 1))
#' @return A matrix with water saturation, gas saturation, oil saturation, and oil relative permeability values, respectively.
#' @examples
#' rel_perm_wgo <- kr3p_Baker(0.15, 0.2, 0.15, 0.15, 0.2, 0.2, 0.05, 0.05, 0.4,
#' 1, 0.3, 3, 2, 4, 2.5, 101)
#'
#' @references
#' \insertRef{Baker1988}{Rrelperm}
#'
#' @export

kr3p_Baker <- function(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON, SGCRIT,
                       KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP) {

   if (is.null(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCON)) stop("'SWCON' must be a numeric equal to or greater than zero.")
   if (SWCON < 0 | SWCON > 1) stop("'SWCON' must be a numeric between zero and one.")

   if (is.null(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SWCRIT)) stop("'SWCRIT' must be a numeric equal to or greater than zero.")
   if (SWCRIT < 0 | SWCRIT > 1) stop("'SWCRIT' must be a numeric between zero and one.")
   if (SWCRIT < SWCON) stop("'SWCRIT' must be equal to or greater than 'SWCON'.")

   if (is.null(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRW)) stop("'SOIRW' must be a numeric equal to or greater than zero.")
   if (SOIRW < 0 | SOIRW > 1) stop("'SOIRW' must be a numeric between zero and one.")

   if (is.null(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORW)) stop("'SORW' must be a numeric equal to or greater than zero.")
   if (SORW < 0 | SORW > 1) stop("'SORW' must be a numeric between zero and one.")
   if (SORW < SOIRW) stop("'SORW' must be equal to or greater than 'SOIRW'.")

   if (is.null(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SOIRG)) stop("'SOIRG' must be a numeric equal to or greater than zero.")
   if (SOIRG < 0 | SOIRG > 1) stop("'SOIRG' must be a numeric between zero and one.")

   if (is.null(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SORG)) stop("'SORG' must be a numeric equal to or greater than zero.")
   if (SORG < 0 | SORG > 1) stop("'SORG' must be a numeric between zero and one.")
   if (SORG < SOIRG) stop("'SORG' must be equal to or greater than 'SOIRG'.")

   if (is.null(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCON)) stop("'SGCON' must be a numeric equal to or greater than zero.")
   if (SGCON < 0 | SGCON > 1) stop("'SGCON' must be a numeric between zero and one.")

   if (is.null(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (!is.numeric(SGCRIT)) stop("'SGCRIT' must be a numeric equal to or greater than zero.")
   if (SGCRIT < 0 | SGCRIT > 1) stop("'SGCRIT' must be a numeric between zero and one.")
   if (SGCRIT < SGCON) stop("'SGCRIT' must be equal to or greater than 'SGCON'.")

   if (is.null(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (!is.numeric(KRWIRO)) stop("'KRWIRO' must be a numeric between zero and one.")
   if (KRWIRO < 0 | KRWIRO > 1) stop("'KRWIRO' must be a numeric between zero and one.")

   if (is.null(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (!is.numeric(KROCW)) stop("'KROCW' must be a numeric between zero and one.")
   if (KROCW < 0 | KROCW > 1) stop("'KROCW' must be a numeric between zero and one.")

   if (is.null(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (!is.numeric(KRGCL)) stop("'KRGCL' must be a numeric between zero and one.")
   if (KRGCL < 0 | KRGCL > 1) stop("'KRGCL' must be a numeric between zero and one.")

   if (is.null(NW)) stop("'NW' must be a numeric between zero and one.")
   if (!is.numeric(NW)) stop("'NW' must be a numeric between zero and one.")
   if (NW < 0) stop("'NW' must be a numeric greater than zero.")

   if (is.null(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (!is.numeric(NOW)) stop("'NOW' must be a numeric between zero and one.")
   if (NOW < 0) stop("'NOW' must be a numeric greater than zero.")

   if (is.null(NG)) stop("'NG' must be a numeric between zero and one.")
   if (!is.numeric(NG)) stop("'NG' must be a numeric between zero and one.")
   if (NG < 0) stop("'NG' must be a numeric greater than zero.")

   if (is.null(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (!is.numeric(NOG)) stop("'NOG' must be a numeric between zero and one.")
   if (NOG < 0) stop("'NOG' must be a numeric greater than zero.")

   if (is.null(NP)) stop("'NP' must be a numeric between zero and one.")
   if (!is.numeric(NP)) stop("'NP' must be a numeric between zero and one.")
   if (NP < 1 | NP > 501) stop("'NP' must be a numeric between one and 500.")

   results <- kr3p_Baker_R(SWCON, SWCRIT, SOIRW, SORW, SOIRG, SORG, SGCON, SGCRIT,
                           KRWIRO, KROCW, KRGCL, NW, NOW, NG, NOG, NP)
   return(results)
}

