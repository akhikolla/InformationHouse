\name{Deff}
\alias{Deff}

\title{Obtaining D-efficiency for estimating model parameters

}
\description{Obtaining D-efficiency for estimating model parameters for a given design under the 3, 4, 5-parameter logistic models. For the output, the function shows three things: 1. the obtained D-optimal design for the given parameter values; 2. the D-efficiency of the given design; and 3. the number of subjects required for the given design in order to provide the same performance as the D-optimal design does.

}
\usage{
Deff(weight,dose,model,P,LB,UB,grid,N_dose)

}

\arguments{
  \item{weight}{A numeric vector. The first K-1 weights for a given design with K design points(K dose levels). The weights represent the proportional allocations of subjects to the dose levels in a given design.

}
  \item{dose}{A numeric vector. K dose levels for a given design.

}
  \item{model}{Numeric. The number of parameters in the model. This is used to specify the model to check the D-efficiency.

}
  \item{P}{A numeric vector. Solicited information on nominal values for the vector for the 5PL model. P=( p1, p2, p3, p4, p5), where p1 is the upper limit of the response, p2 is the slope that control the stiffness of the response curve, p3 is the position of the transition region in dose(this is the ED50 under the 3PL and the 4PL model), p4 is the lower limit of the response, and p5 is the asymmetric factor. When a user computes the D-efficiency under the 3PL or 4PL model, P is the vector of model parameter values for each model with setting p4 and p5 differently. Details are given in the examples.

}
\item{LB}{Numeric. Predetermined lower bound of the dose range for the original dose level.

}
  \item{UB}{Numeric. Predetermined lower bound of the dose range for the original dose level.

}
\item{grid}{Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.

}
 \item{N_dose}{Set to TRUE if the original dose levels take negative values. Set to FALSE if the original dose levels take positive values. Default is FALSE

}
}

\references{Hyun, S. W., Wong, W.K., Yang, Y. Optimal designs for asymmetric sigmoidal response curves in bioassays and immunoassays. (submitted to Statistical Methods in Medical Research)

Hyun, S.W., Wong, W.K, Yang, Y. VNM: An R Package for Finding Multiple-Objective Optimal
Designs for the 4-Parameter Logistic Model. (Journal of Statistical Software, 83, 1-19, 2018, doi:10.18637/jss.v083.i05.)

}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang

}

\examples{

##Under the 5PL model with the parameter values
##T5=(30000,0.5,800,0.5,2) and the dose range [1.95,32000],
##find the D-efficiency of the broad range design.

##The broad range design
dose=c(1.95,7.8,31.25,125,500,2000,8000,32000)
dlength=length(dose)
weight=rep(1/dlength,dlength-1)

##Parameter values
T5=c(30000,0.5,800,0.5,2)

##Compute D-efficiency of the broad range design
Deff(weight,dose,model=5,P=T5,LB=1.95,UB=32000,grid=.01)


}

