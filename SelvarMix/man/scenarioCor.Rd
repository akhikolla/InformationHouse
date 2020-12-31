\docType{data}
\name{scenarioCor}
\alias{scenarioCor}
\title{Simulated quantitative data according SRUW modeling}
\format{A data matrix with 2000 observations on 14 variables and the last column contains the labels.

\describe{

  \item{\code{scenarioCor[,1:14]}}{a numeric matrix containing the observations}

  \item{\code{scenarioCor[,15]}}{an integer vector containing the labels}

}}
%\source{
%  \url{}
%}
\description{
 The dataset consists of 2000 data points in \eqn{R^{14}}. 
 On the subset of relevant clustering variables \eqn{S = \{1, 2\}}, 
 data are distributed from a mixture of four equiprobable 
 spherical Gaussian distributions with means \eqn{(0,0), (4,0)
 (0,2)} and \eqn{(4,2)}. The subset of redundant variables 
 is \eqn{U =\{3-11\}} that are explained by the subset of 
 predictor variables \eqn{R = \{1,2\}}. 
 The last three variables are independent 
 \eqn{W = \{11, 12, 13\}}.
 
}
\details{
The subset \eqn{U} of redundant variables is simulated as follows :

 \eqn{x^{U} =  (0,0, 0.4, 0.8, ..., 2) + x^{S} b + \varepsilon}, with \eqn{\varepsilon \sim N(0_9, \Omega)}
 
The subset \eqn{W} of independent variables is simulated as follows :

\eqn{x^{W} \sim N((3.2, 3.6, 4), I_3)}

For more details on the regression coefficients \eqn{b} and the covariance
matrix \eqn{\Omega} see Maugis et al.(2009).
}
\examples{
data(scenarioCor)
}
\references{
 Maugis, C., Celeux, G., and Martin-Magniette, M. L., 2009. 
  "Variable selection in model-based clustering: 
  A general variable role modeling". Computational 
  Statistics and Data Analysis, vol. 53/11, pp. 3872-3882.

%   Celeux, G., Martin-Magniette, M. L., Maugis, C., and Raftery, A. E., 2014. 
%  "Comparing Model Selection and Regularization Approaches to Variable
%   Selection in Model-Based Clustering". Journal de la Soci\'et\'e Fran\,caise
%   de Statistique, vol. 155/2, pp. 57-71.

}
\keyword{datasets}
