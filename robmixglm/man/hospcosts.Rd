\name{hospcosts}
\alias{hospcosts}
\docType{data}
\title{Hospital Costs data}
\description{
Data for the analysis in Beath (2018),  previously analysed in Marazzi and Yohai (2004),  Cantoni and Ronchetti (2006) and Heritier et al (2009). The data is for 100 patients hospitalised at the Centre Hospitalier Universitaire Vaudois in Lausanne,  Switzerland for "medical back problems" (APDRG 243).
}
\usage{hospcosts}
\format{
  A data frame with 100 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{patient id}
    \item{\code{costs}}{cost of stay in Swiss francs}
    \item{\code{los}}{length of stay in days}
    \item{\code{adm}}{admission type,  0 = planned,  1 = emergency}
    \item{\code{ins}}{insurance type,  0 = regular,  1 = private}
    \item{\code{age}}{age in years}
    \item{\code{sex}}{sex,  0 = female,  1 = male}
    \item{\code{dest}}{discharge destination,  0 = another health institution,  1 = home}
    \item{\code{loglos}}{log of length of stay}
  }
}
\source{
  Heritier et al (2009)
  }
\references{
Cantoni,  E.,  & Ronchetti,  E. (2006). A robust approach for skewed and heavy-tailed outcomes in the analysis of health care expenditures. Journal of Health Economics,  25(2),  198213. http://doi.org/10.1016/j.jhealeco.2005.04.010

Heritier,  S.,  Cantoni,  E.,  Copt,  S. and Victoria-Feser,  M-P (2009). Robust Methods in Biostatistics. Wiley.

Marazzi,  A.,  & Yohai,  V. J. (2004). Adaptively truncated maximum likelihood regression with asymmetric errors. Journal of Statistical Planning and Inference,  122(12),  271291. http://doi.org/10.1016/j.jspi.2003.06.011
}
\examples{
\donttest{
hospcosts.robustmix <- robmixglm(costs~adm+age+dest+ins+loglos+sex, family = "gamma", 
    data = hospcosts, cores = 1)
summary(hospcosts.robustmix)
}}
\keyword{datasets}
