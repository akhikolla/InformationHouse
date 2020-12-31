\name{diabdata}
\alias{diabdata}
\docType{data}
\title{Diabetes data}
\description{
Data from Heritier et al (2009),  originally from Harrell (2001, p379). This data was from a study of the prevalence of cardiovascular risk factors such as obesity and diabetes for African Americans. (Willems et al,  19997) Data was available for 403 subjects screened for diabetes,  reduced to 372 after removal of cases with missing data.   
}
\usage{diabdata}
\format{
  A data frame with 372 observations on the following 8 variables.
  \describe{
    \item{\code{glyhb}}{Glycosated haemoglobin (values above 7.0 are usually taken as a positive diagnosis of diabetes)}
    \item{\code{age}}{age in years}
    \item{\code{gender}}{male or female}
    \item{\code{bmi}}{body mass index in kg/m^2}
    \item{\code{waisthip}}{ratio of waist to hip measurement}
    \item{\code{frame}}{body frame,  small,  medium or large}
    \item{\code{stab.glu}}{glucose}
    \item{\code{location}}{location,  Buckingham or Louisa}
  }
}
\source{
  Heritier et al (2009)
  }
\references{
Harrell,  F.E. (2001). Regression Modeling Strategies: With Applications to Linear Models,  Logistic Regression and Survival Analysis. Springer.

Heritier,  S.,  Cantoni,  E.,  Copt,  S. and Victoria-Feser,  M-P (2009). Robust Methods in Biostatistics. Wiley.

Willems,  J.P.,  Saunders,  J.T.,  Hunt,  D.E. and Schorling,  J.B. (1997) Prevalence of coronary heart disease risk factors among rural blacks: A community-based study. Southern Medical Journal,  90:814-820.
}
\examples{
\donttest{
diabdata.robustmix <- robmixglm(glyhb~age+gender+bmi+waisthip+frame+location, 
         data = diabdata, cores = 1)
summary(diabdata.robustmix)

diabdata.step <- step(diabdata.robustmix, glyhb~age+gender+bmi+waisthip+frame+location, cores = 1)
summary(diabdata.step)
}}
\keyword{datasets}
