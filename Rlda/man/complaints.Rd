\name{complaints}
\alias{complaints}
\docType{data}
\title{
Complaints received for the **Bureau of Consumer Financial Protection** in US about financial products and services.
}
\description{
Specifically in this dataset we work with only credit card complaint's for the 2015 year.
}
\usage{data("complaints")}
\format{
  A data frame with 17301 observations on the following 3 variables.
  \describe{
    \item{\code{Product}}{a factor with levels \code{Credit card}}
    \item{\code{Issue}}{a factor with levels \code{Advertising and marketing} \code{Application processing delay} \code{APR or interest rate} \code{Arbitration} \code{Balance transfer} \code{Balance transfer fee} \code{Bankruptcy} \code{Billing disputes} \code{Billing statement} \code{Cash advance} \code{Cash advance fee} \code{Closing/Cancelling account} \code{Convenience checks} \code{Credit card protection / Debt protection} \code{Credit determination} \code{Credit line increase/decrease} \code{Customer service / Customer relations} \code{Delinquent account} \code{Forbearance / Workout plans} \code{Identity theft / Fraud / Embezzlement} \code{Late fee} \code{Other} \code{Other fee} \code{Overlimit fee} \code{Payoff process} \code{Privacy} \code{Rewards} \code{Sale of account} \code{Transaction issue} \code{Unsolicited issuance of credit card}}
    \item{\code{Company}}{a factor variable describing the companies available}
  }
}
\examples{
data(complaints)
}
\keyword{datasets}
