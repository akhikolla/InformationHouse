\name{RDOPT}
\alias{RDOPT}

\title{Search the robust D-optimal designs for estimating model parameters

}
\description{Obtain robust D-optimal design (i.e., optimal dose levels and corresponding optimal weights) for estimating the model parameters that works well under the 3, 4, 5-parameter logistic models. In addition, the function can be used to obtain D-optimal design for each of the three models.

}
\usage{
RDOPT(LB,UB, P3,P4,P5,q,grid,r,epsilon,N_dose,log_scale)
}

\arguments{
  \item{LB}{Numeric. Predetermined lower bound of the dose range for the original dose level.

}
  \item{UB}{Numeric. Predetermined lower bound of the dose range for the original dose level.

}
  \item{P3}{A numeric vector. Solicited information on nominal values for the vectors for the 3PL model. P3=(p1, p2, p3), where p1 is the upper limit of the response, p2  is the slope that control the stiffness of the response curve, and p3 is the ED50.

}
  \item{P4}{A numeric vector. Solicited information on nominal values for the vectors for the 4PL model. P4=(p1, p2, p3, p4),  where p1 is the upper limit of the response, p2 is the slope that control the stiffness of the response curve, p3 is the ED50, and p4 is the lower limit of the response.

}
  \item{P5}{A numeric vector. Solicited information on nominal values for the vectors for the 5PL model. P5=(p1, p2, p3, p4, p5), where p1 is the upper limit of the response, p2 is the slope that control the stiffness of the response curve, p3 is the position of the transition region in dose, p4 is the lower limit of the response, and p5 is the asymmetric factor. For searching the robust D-optimal design, the nominal values for P3, P4, P5 need to be specified. However, when a user searches the D-optimal design for the 3PL, 4PL, or 5PL model, only P5 needs to be specified. Details are given in the examples.

}
  \item{q}{A numeric vector. User select weights. q=c(q1, q2, q3), where q1, q2, q3 represent the relative importance of the 3, 4, 5PL models respectively. They are non-negative and q1+q2<=1. For example, q=c(0,0,1) is used to find D-optimal design for the 5PL model only and q=c(1/3,1/3,1/3) is used to find robust D-optimal design for the three models under the assumption that the three models are equally important.

}
  \item{grid}{Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.

}
  \item{r}{Numeric. The number of iterations to set an initial design to search the D-optimal design. Default is 30 and needed to be increased (for example, r=50) if the searched optimal design is not a true optimal. It can be verified by the verification plot.

}
  \item{epsilon}{Numeric. Stopping criterion for the algorithm to search the D-optimal design. Default is 0.001.

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

##1. Search the robust D-optimal design for estimating the model parameters
##under the 3,4,5PL models

##Parameter values for the 5PL model
T5=c(30000,0.5,800,0.5,2)
##Parameter values for the 4PL model
T4=c(27264.92,0.67,3340.95,-225.55)
##Parameter values for the 3PL model
T3=c(26715.52,0.70,3204.92)

##Dose range=[1.95,32000]
##Assuming the three models are equally important
RDOPT(LB=1.95,UB=32000,P3=T3,P4=T4,P5=T5,q=c(1/3,1/3,1/3))

##2. Search the D-optimal design for the 3PL model using
##the parameter values T3=(26715.52,0.70,3204.92)
##and the dose range [1.95,32000].

#Set the parameter values for the 3PL model in the form of P5
P5=c(26715.52,0.70,3204.92,0,1)
RDOPT(LB=1.95, UB=32000,P5=P5,q=c(1,0,0))

##3. Search the D-optimal design for the 4PL model using
##the parameter values T4=(27264.92,0.67,3340.95,-225.55)
##and the dose range [1.95,32000].

#Set the parameter values for the 4PL model in the form of P5.
P5=c(27264.92,0.67,3340.95,-225.55,1)
RDOPT(LB=1.95, UB=32000,P5=P5,q=c(0,1,0))

}

