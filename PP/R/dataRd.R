#' Simulated data set
#' 
#' This dataset contains data simulated with the \code{sim_4pl()} function.
#'
#' @format A data.frame with 60 rows and 14 columns.
#' @source simulation
#' @seealso \link{sim_4pl}
#' @name fourpl_df
NULL

#' Adaptive Matrices Test data
#' 
#' This dataset contains real data from the 'Adaptive Matrices Test' (AMT), which is a computer-administered test.
#' This power test assesses logical reasoning as an indicator of general intelligence. The ability to identify regularities and draw logical conclusions is a very good predictor of long-term success at work. The dataset is sparse, because the test tailores a specific set of items for each examinee's ability level. (More information about adaptive testing: \url{https://en.wikipedia.org/wiki/Computerized_adaptive_testing})
#'
#' 
#' The data are provided from the Unitersity of Vienna, Faculty of Psychology, Department of Psychological Assessment.
#' Thanks to Schuhfried \url{https://www.schuhfried.at/test/AMT}.
#' 
#' @format A list with two data.frames. The first data.frame 'daten_amt' contains 298 columns and 710 rows. Each row contains responses from on examinee. The second data.frame 'betas' contains the difficulty parameter (1PL) (These parameters came with the raw-score extraction.).
#' @source Division of Psychological Assessment and Applied Psychometrics, Fakulty of Psychology, University of Vienna
#' \itemize{
#' \item ID: id of person
#' \item AGE: age in years (with ages below 18 and above 34 are collapsed)
#' \item TE_GA:
#' \itemize{
#' \item TE: self-assessment. To pass a psychological assessment course, the students have to complete several hours self assessment on a bunch of tests, to get familiar with them.
#' \item GA: testing for an assessment report. The students also have to test other people (not psychologists nor psychology students) in order to write an assessment report.
#' }
#' 
#' 
#' \item FORM: There are several different versions of this test, which differ in test length, duration etc \ldots
#' \item TIME1: start time
#' \item TIME2: end time
#' \item REL: reliability for each person
#' \item i: items
#' }
#' @references 
#' \itemize{
#' \item Hornke, L. F., Etzel, S., & Rettig, K. (2003). Manual Adaptive Matrices Test (AMT). \emph{Moedling: SCHUHFRIED GmbH}.
#' }
#' 
#' @seealso \link{PPass}
#' @name pp_amt
NULL
