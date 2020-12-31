\name{DsOPT}
\alias{DsOPT}

\title{Search Ds-optimal design for estimating the asymmetric factor under the 5-parameter logistic model.

}
\description{Obtain Ds-optimal design (i.e., optimal dose levels and corresponding optimal weights) for estimating the asymmetric factor under the 5-parameter logistic model.

}
\usage{
DsOPT(LB,UB,P,grid,r,epsilon,N_dose,log_scale)

}

\arguments{
    \item{LB}{Numeric. Predetermined lower bound of the dose range for the original dose level.

}
  \item{UB}{Numeric. Predetermined lower bound of the dose range for the original dose level.

}
  \item{P}{A numeric vector. Solicited information on nominal values for the vector for the 5PL model. P=(p1, p2, p3, p4, p5), where p1 is the upper limit of the response, p2 is the slope that control the stiffness of the response curve, p3 is the position of the transition region in dose, p4 is the lower limit of the response, and p5 is the asymmetric factor.

}
 \item{grid}{Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.

}
  \item{r}{Numeric. The number of iterations to set an initial design to search the Ds-optimal design. Default is 30 and needed to be increased (for example, r=50) if the searched optimal design is not a true optimal. It can be verified by the verification plot.

}
  \item{epsilon}{Numeric. Stopping criterion for the algorithm to search the Ds-optimal design. Default is 0.001.

}
  \item{N_dose}{Set to TRUE if the original dose levels take negative values. Set to FALSE if the original dose levels take positive values. Default is FALSE

}
  \item{log_scale}{Set to TRUE if a user wants to set the X-axis(dose level) on log scale in the verification plot of the optimal design. Otherwise, set to FALSE. Default is TRUE. Note: this need to set FALSE when the original dose levels take negative values.

}
}

\references{Hyun, S. W., Wong, W.K., Yang, Y. Optimal designs for asymmetric sigmoidal response curves in bioassays and immunoassays. (submitted to Statistical Methods in Medical Research)

Hyun, S.W., Wong, W.K, Yang, Y. VNM: An R Package for Finding Multiple-Objective Optimal
Designs for the 4-Parameter Logistic Model. (Journal of Statistical Software, 83, 1-19, 2018, doi:10.18637/jss.v083.i05.)

}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang

}

\examples{

##Find the Ds-optimal design for estimating theta_5 under the 5PL model
##with the parameter values T5=(30000,0.5,800,0.5,2) and the dose range [1.95,32000].

##Parameter values for the 5PL model
T5=c(30000,0.5,800,0.5,2)

##Find the Ds-optimal design
DsOPT(LB=1.95, UB=32000, P=T5,grid=.01,epsilon=.001)

}


