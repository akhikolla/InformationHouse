#' Gastric Cancer Data.
#'
#' Survival data from a trial comparing chemotherapy versus combined chemotherapy plus radiotherapy in the treatment of gastric cancer.
#'
#' @usage data(gastric)
#'
#' @format A data frame with 90 observations (45 in each treatment group) with the following 3 columns:
#' \describe{
#'   \item{time}{the observed follow-up times in days.}
#'   \item{status}{the event indicators, 0=right censored, 1= event.}
#'   \item{group}{ the group labels, 1 = chemotherapy, 2 = chemotherapy plus radiotherapy.}
#' }
#' 
#' @source Stablein, D. M. and Koutrouvelis, I. A. (1985) A two-sample test sensitive to crossing hazards in uncensored and singly
#' censored data. Biometrics 41, 643–652. (Page 649).
#' 
#' @references Gastrointestinal Tumor Study Group: Schein, P. D., Stablein, D. M., Bruckner, H. W., Douglass, H.
#' O., Mayer, R., et al. (1982). A comparison of combination chemotherapy and combined modality
#' therapy for locally advanced gastricarcinoma. Cancer 49, 1771-1777.
"gastric"



#' Urothelial carcinoma.
#'
#' Survival data from a trial comparing chemotherapy versus atezolizumab in the treatment of Urothelial carcinoma.
#'
#' @usage data(carcinoma)
#'
#' @format A data frame with 625 observations (316 in the atezolizumab group and 309 chemotherapy group) with the following 3 columns:
#' \describe{
#'   \item{time}{the observed follow-up times in days.}
#'   \item{status}{the event indicators, 0=right censored, 1= event.}
#'   \item{group}{ the group labels, 1 = atezolizumab, 2 = chemotherapy.}
#' }
#' 
#' @references Powles T, Dura?n I, van der Heijden MS, et al. Atezolizumab versus chemotherapy in patients with platinum-treated locally
#' advanced or metastatic urothelial carcinoma (IMvigor211): a multicentre, open-label, phase 3 randomised controlled trial.
#' Lancet 2018; 391: 748-757.
"carcinoma"