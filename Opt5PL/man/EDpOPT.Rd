\name{EDpOPT}
\alias{EDpOPT}

\title{Search c-optimal designs for estimating the EDp under the 5-parameter logistic model

}
\description{
Obtain c-optimal design (i.e., optimal dose levels and corresponding optimal weights) for estimating the EDp under the 5-parameter logistic model. The EDp is the dose level that achieve the 100p\% of the difference between the maximum and the minimum responses.
}
\usage{
EDpOPT(LB,UB,P,EDp,grid,r,epsilon,N_dose,log_scale)
}

\arguments{
  \item{LB}{Numeric. Predetermined lower bound of the dose range for the original dose level.

}
  \item{UB}{Numeric. Predetermined lower bound of the dose range for the original dose level.

}
  \item{P}{A numeric vector. Solicited information on nominal values for the vector for the 5PL model. P=(p1, p2, p3, p4, p5), where p1 is the upper limit of the response, p2 is the slope that control the stiffness of the response curve, p3 is the position of the transition region in dose, p4 is the lower limit of the response, and p5 is the asymmetric factor.

}
  \item{EDp}{Numeric. An user select value to define the EDp. For example, EDp=.05 to define the ED50.

}
 \item{grid}{Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.

}
  \item{r}{Numeric. The number of iterations to set an initial design to search the c-optimal design. Default is 30 and needed to be increased (for example, r=50) if the searched optimal design is not a true optimal. It can be verified by the verification plot.

}
  \item{epsilon}{Numeric. Stopping criterion for the algorithm to search the c-optimal design. Default is 0.001.

}
  \item{N_dose}{Set to TRUE if the original dose levels take negative values. Set to FALSE if the original dose levels take positive values. Default is FALSE

}
  \item{log_scale}{Set to TRUE if a user wants to set the X-axis(dose level) on log scale in the verification plot of the optimal design. Otherwise, set to FALSE. Default is TRUE. Note: this need to set FALSE when the original dose levels take negative values.

}
}

\references{
Hyun, S. W., Wong, W.K., Yang, Y. Optimal designs for asymmetric sigmoidal response curves in bioassays and immunoassays. (submitted to Statistical Methods in Medical Research)

Hyun, S.W., Wong, W.K, Yang, Y. VNM: An R Package for Finding Multiple-Objective Optimal
Designs for the 4-Parameter Logistic Model. (Journal of Statistical Software, 83, 1-19, 2018, doi:10.18637/jss.v083.i05.)

}
\author{
Seung Won Hyun, Weng Kee Wong, and Yarong Yang
}





\examples{
##Search the c-optimal design for estimating the ED50 under the 5PL model

##Parameter values for the 5PL model
T5=c(30000,0.5,800,0.5,2)

##Dose range=[1.95,32000]

EDpOPT(LB=1.95,UB=32000, P=T5,EDp=.5,grid=.01, r =30,epsilon=.001)

}
