#' Dataset funcdep_data
#'
#' This dependence structure is constructed by applying the function :
#' \deqn{h(u_1,u_2,u_3) = (u_{1},\sin(2\pi u_{1})-\frac{u_{2}}{\pi},(1+\frac{u_{3}}{\pi^{2}})(\frac{u_{3}}{2} I_{\frac{1}{4}\ge u_1}-\sin(\pi^{x_{1}}) I_{\frac{1}{4} < u_{1}}))}
#' to uniformly drawn 3-dimensional random vectors. The dataset is the ranks of \eqn{h(u)}.
#'
#' This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
#'
#' @format A matrix with 500 rows and 3 columns
#'
#' The example section below gives the code to re-generate this data if needed.
#'
#' @examples
#' set.seed(seed = 12,kind = "Mersenne-Twister",normal.kind = "Inversion")
#' x = matrix(runif(1500),500,3)
#' x[,2] = sin(2*pi*x[,1])-x[,2]/pi
#' x[,3] = (x[,3]*(x[,1]<1/4)/2 - sin(pi**(x[,1]))*(x[,1]>1/4))*(1+x[,3]/(pi^2))
#' funcdep_data = apply(x,2,function(x){return(rank(x,ties.method = "max"))})/(501)
#'
#' @references
#' \insertRef{laverny2020}{cort}
"funcdep_data"

#' Dataset impossible_data
#'
#' We simulate from a density inside the piecewise linear copula class, by applying the function:
#' \deqn{h(u) = (u_1,          \frac{u_2}{2} + \frac{1}{2}I_{u_1 \notin (\frac{1}{3}, \frac{2}{3})})}
#' to a 200x2 uniform sample, and taking ranks.
#'
#' This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
#'
#' @format A matrix with 200 rows and 2 columns
#'
#' The example section below gives the code to re-generate this data if needed.
#'
#' @examples
#' set.seed(seed = 12, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' x = matrix(runif(400),200,2)
#' x = t(apply(x, 1,function(u){
#'   if(u[1]< 1/3){
#'     u[2] = 1/2 + u[2]/2
#'   } else{ if(u[1]<2/3){
#'     u[2] = u[2]/2
#'   } else {
#'     u[2] = 1/2 + u[2]/2
#'   }}
#'   return(u)
#' }))
#' impossible_data = apply(x,2,function(x){return(rank(x,ties.method = "max"))})/(201)
#'
#' @references
#' \insertRef{laverny2020}{cort}
"impossible_data"

#' Dataset recoveryourself_data
#'
#' This dataset is a simple test: we simulate random samples from a density inside the piecewise copula class,
#' and test whether or not the estimator can recover it. For that, we will use a 2-dimensional sample with 500
#' observations, uniform on the unit hypercube, and apply the following function:
#' \deqn{h(u) = (u_1, \frac{u_2 + I_{u_1 \le \frac{1}{4}} + 2I_{u_1 \le \frac{1}{2}} + I_{\frac{3}{4} \le u_1}}{4})}
#'
#' This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
#'
#'
#' @format A matrix with 500 rows and 2 columns
#'
#' The example section below gives the code to re-generate this data if needed.
#'
#' @examples
#' set.seed(seed = 12, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' x = matrix(runif(1000),500,2)
#' recoveryourself_data = t(apply(x, 1,function(u){
#'   if(u[1]< 1/4){
#'     u[2] = 3/4 + u[2]/4
#'   } else{ if(u[1]<1/2){
#'     u[2] = 1/2 + u[2]/4
#'   } else { if(u[1]<3/4){
#'     u[2] = u[2]/4
#'   } else {
#'     u[2] = 1/4 + u[2]/4
#'   }}}
#'   return(u)
#' }))
#'
#' @references
#' \insertRef{laverny2020}{cort}
"recoveryourself_data"

#' Dataset clayton_data
#'
#' This dataset is a simulation of 200 points from a 3-dimensional clayton copula with \eqn{\theta = 7},
#' hence highly dependent, for the first, third and fourth marginals. The second marginal is added
#' as independent uniform draws. Lastly, the third marginal is flipped,
#' inducing a negative dependence structure.
#'
#' This dataset is studied in O. Laverny, V. Maume-Deschamps, E. Masiello and D. Rullière (2020).
#'
#' @format A matrix with 200 rows and 4 columns
#'
#' The example section below gives the code to re-generate this data if needed.
#'
#' @examples
#' psi <- function(t,alpha) (1 + sign(alpha)*t) ^ (-1/alpha) # generator
#' rClayton <- function(n,dim,alpha){
#'   val <- matrix(runif(n * dim), nrow = n)
#'   gam <- rgamma(n, shape = 1/alpha, rate = 1)
#'   gam <- matrix(gam, nrow = n, ncol = dim)
#'   psi(- log(val) / gam,alpha)
#' }
#' set.seed(12,kind = "Mersenne-Twister",normal.kind = "Inversion")
#' clayton_data <- matrix(nrow=200,ncol=4)
#' clayton_data[,c(1,4,3)] = rClayton(n=200,dim=3,alpha=7)
#' clayton_data[,2] = runif(200)
#' clayton_data[,3] <- 1 - clayton_data[,3]
#'
#' @references
#' \insertRef{laverny2020}{cort}
"clayton_data"
