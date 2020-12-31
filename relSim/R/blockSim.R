#' Perform relatives simulations using large memory blocks in C
#' 
#' Generate N pairs with a given relationship, calculate the LR for sibs,
#' parent-child and the number of matching alleles and count the number of
#' pairs that meet the threshold criteria.
#' 
#' This function is used for fast accurate estimation of false positive and
#' false negative rates. It achieves part of its speed by block exectution in
#' C, and part by not saving the LR or IBS results. It can do 1 billion
#' iterations in about an hour.
#' 
#' @param N The number of iterations to carry out
#' @param Freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex.
#' @param rel generate unrelated (\code{rel = 'UN'}), full-sibs (\code{rel =
#' 'FS'}), or parent child (\code{rel = 'PC'}) pairs
#' @param ibsthresh A vector of one or more IBS thresholds
#' @param kithresh A vector of one or more KI/LR thresholds
#' @param code A code from 1 to 6 which dictates the events that will be
#' counted.  \enumerate{ \item the LR for siblings will be compared to the
#' values in \code{kithresh} and incremented if the LR is greater than the
#' threshold \item the LR for parent/child will be compared to the values in
#' \code{kithresh} and incremented if the LR is greater than the threshold
#' \item the number of matching alleles (IBS) will be compared to the values in
#' \code{ibsthresh} and incremented if the IBS is greater than the threshold
#' \item the LR for siblings and the number of matching alleles will be
#' compared to the values in \code{kithresh} and \code{ibsthresh} and
#' incremented if both the LR and IBS is greater than the thresholds.
#' \code{ibsthresh} and \code{kithresh} must be of equal length for this option
#' to work \item the LR for parent/child and the number of matching alleles
#' will be compared to the values in \code{kithresh} and \code{ibsthresh} and
#' incremented if both the LR and IBS is greater than the thresholds.
#' \code{ibsthresh} and \code{kithresh} must be of equal length for this option
#' to work \item this option is equivalent to performing code 4 and 5
#' simulataneously. It is not currently implemented }
#' @param falseNeg if TRUE then the number of results that DO NOT satisfy the
#' conditions are counted, otherwise the number of results DO satisfy the
#' conditions are counted
#' @param BlockSize Sets the number of random profiles to be generated in each
#' iteration. By default the block size is set to 10 percent of the total
#' sample size. It is unclear whether the procedure is more efficient if a
#' bigger percentage of the total is used. Users must take care to make sure
#' that the block size evenly divides \code{N} otherwise the procedure will
#' exit. Users must also make sure that they have enough memory.
#' @param showProgress If \code{TRUE} then a progress bar will be displayed in
#' the console showing the progress of the simulation.
#' @return A vector containing the number of profile pairs that satisfied the
#' threshold conditions
#' @author James M. Curran
#' @seealso sim
#' @examples
#' 
#' ## not run
#' ## this counts the number of unrelated pairs that are falsely identified
#' ## as siblings using the policy that there are 16 or more matching
#' ## alleles, and the LR/KI is greater than 100,000
#' ## this is a very rare event for the FBI Caucasians with a frequency of
#' ## about 4-5 times in 10 million pairs
#' \dontrun{
#' data(fbiCaucs)
#' N = 1e8
#' ki = 1e5
#' ibs = 16
#' code = 5
#' BlockSize = 1e6
#' blockSim(N, fbiCaucs, rel = "UN", ibsthresh = ibs, kithresh = ki,
#'          code = code, falseNeg = FALSE, BlockSize = BlockSize)
#' }
#' @importFrom utils txtProgressBar setTxtProgressBar
blockSim = function(N, Freqs, rel = "UN", ibsthresh = NULL, kithresh = NULL,
                    code = 1, falseNeg = TRUE, BlockSize = N/10, showProgress = FALSE){
  
  rel = toupper(rel)
  if(!grepl("(UN|FS|PC)", rel)){
    stop("rel must be one of 'UN', 'FS' or 'PC'")
  }
  
  nBlocks = N / BlockSize
  
  if(is.null(ibsthresh) & is.null(kithresh))
    stop("You must specify one or both of ibsthresh or kithresh")
  
  nResults = 0
  if(is.null(ibsthresh) & !is.null(kithresh)){
    nResults = length(kithresh)
    ibsthresh = rep(0, nResults) ## dummy vals
  }else if(is.null(kithresh) & !is.null(ibsthresh)){
    nResults = length(ibsthresh)
    kithresh = rep(0, nResults)
  }else{
    if(length(ibsthresh) != length(kithresh)){
      stop("ibsthresh and kithresh must be the same length")
    }else{
      nResults = length(ibsthresh)
    }
  }
  
  if(nResults == 0)
    stop("Nothing to count")
  
  nTotal = rep(0, nResults)
  
  if(showProgress){
    pb = txtProgressBar(min = 0, max = nBlocks, style = 3)
  }
  
  if(rel == "UN"){
    for(block in 1:nBlocks){
      Prof1 = randomProfiles(Freqs$freqs, BlockSize)
      Prof2 = randomProfiles(Freqs$freqs, BlockSize)
      
      count = blockStatCounts(Prof1, Prof2, BlockSize, Freqs$freqs, 
                              code, falseNeg, ibsthresh, kithresh, nResults)
      nTotal = nTotal + count
      if(showProgress){
        setTxtProgressBar(pb, block)
      }
    }
  }else if(rel == "FS"){
    for(block in 1:nBlocks){
      Prof1 = randomProfiles(Freqs$freqs, BlockSize)
      Prof2 = randomSibs(Prof1, Freqs$freqs, BlockSize)
      
      count = blockStatCounts(Prof1, Prof2, BlockSize, Freqs$freqs, 
                              code, falseNeg, ibsthresh, kithresh, nResults)
      
      nTotal = nTotal + count
      if(showProgress){
        setTxtProgressBar(pb, block)
      }
    }
  }else if(rel == "PC"){
    for(block in 1:nBlocks){
      Prof1 = randomProfiles(Freqs$freqs, BlockSize)
      Prof2 = randomChildren(Prof1, Freqs$freqs, BlockSize)
      
      count = blockStatCounts(Prof1, Prof2, BlockSize, Freqs$freqs, 
                              code, falseNeg, ibsthresh, kithresh, nResults)
      
      nTotal = nTotal + count
      if(showProgress){
        setTxtProgressBar(pb, block)
      }
    }
  }
  
  if(showProgress){
    close(pb)
  }
  
  return(list(nTotal = nTotal, N = N, p = nTotal / N))
}
