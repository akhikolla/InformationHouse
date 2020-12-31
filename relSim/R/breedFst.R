#' Breed a population with an approximate level of
#' \eqn{\theta}{theta}(\eqn{F_{ST}}{Fst})
#' 
#' This function simulates a population with an approximate level of population
#' substructure. This is achieved by subdividing a population into equal sized
#' subpopulations and allowing them to breed within themselves for \deqn{t =
#' }{t = ceiling(log(1-theta)/log(1-1/(2*Ns)))}\deqn{
#' \lceil{\frac{\log_e(1-\theta)}{\log\left(1-\frac{1}{2N_s}\right)}}\rceil}{t
#' = ceiling(log(1-theta)/log(1-1/(2*Ns)))} generations, where \eqn{N_s}{Ns} is
#' the number of individuals in each subpopulation. This will produce a
#' population with an estimated coancestry coefficient approximately equal to
#' \eqn{\theta}{theta}
#' 
#' 
#' @param Freqs A list with an element, \code{freqs} which contains a list of
#' vectors, where each vector is a set of allele frequencies for a locus
#' @param theta A desired level of inbreeding, where \eqn{0 < \theta < 0.5}{0 < theta < 0.5}
#' @param N Total population size
#' @param ns The number of subpopulations. \eqn{N/n_s}{N/ns} needs to be
#' greater than 100
#' @param DNAtools If \code{TRUE} then the profiles in the return population
#' will be formatted as a data frame with an id column and two columns per
#' locus.
#' @return An object of class 'population' which is a list with the following
#' elements \itemize{ \item \code{profiles} - a vector of profiles where the
#' level of inbreeding is approximately equal to \eqn{\theta}{theta} \item
#' \code{nProfiles} - the total number of individuals in the population \item
#' \code{nSubpops} - the number of sub-populations in the population \item
#' \code{nLoci} - the number of loci each individual is typed at \item
#' \code{theta} - the desired level of substructure in the population. The
#' actual value will be near to this.  \item \code{Freqs} - a Freq object
#' representing the ancestral frequencies of the population }
#' @author James M. Curran
#' @examples
#' 
#' data(USCaucs)
#' pop = breedFst(USCaucs)
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export breedFst
breedFst = function(Freqs, theta = 0.01, N = 10000, ns = 10,
                    DNAtools = FALSE){
    if(N<1000){
        stop("N must be >= 1000")
    }

    if(N %% ns != 0){
        stop("ns must divide N into a whole number\tThat is the subpopulation sizes must be integers")
    }

    if(theta <=0 || theta >= 0.5){
        stop("0 < theta < 0.5")
    }

    Ns = N / ns

    nGen = ceiling(log(1 - theta) / log(1 - 1/(2*Ns)))
    cat(paste("Breeding for", nGen, "generations\n"))

    nLoci = length(Freqs$freqs)

    ## generate the parental population

    parents = .randomProfiles(Freqs$freqs, N)
    pb = txtProgressBar(min = 0, max = nGen, style = 3)

    for(t in 1:nGen){
        parents = .breed(parents, ns, Ns, nLoci)
        setTxtProgressBar(pb, t)
    }
    setTxtProgressBar(pb, nGen)

    if(DNAtools){
        parents = data.frame(cbind(1:N, matrix(parents, nrow = N, byrow = T)))
        colNames = paste(rep(Freqs$loci, rep(2, nLoci)), c('a','b'), sep = '_')
        names(parents) = c('id', colNames)
    }

    pop = list(profiles = parents, nProfiles = N,  nSubpops = ns,
               nLoci = nLoci, theta = theta, Freqs = Freqs)
    class(pop) = "population"

    return(pop)
}
