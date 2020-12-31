#' Perform the relatives simulation
#' 
#' Generate N pairs with a given relationship and calculate the LR for sibs,
#' parent-child and the number of matching alleles
#' 
#' This is the function that generates all the data for the results in the
#' paper. WARNING: this function is not especially fast. To achieve the 100
#' million iterations used in the paper, 30 instances of R were launched on a
#' multicore server. Each instance represented one relationship with 10 million
#' iterations. The compute time for this arrangement was approximately 1 hours,
#' meaning a full serial run would have taken over 30 hours to achieve the same
#' result.
#' 
#' @param N The number of iterations to carry out
#' @param Freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex.
#' @param rel generate unrelated (\code{rel = 'UN'}), full-sibs (\code{rel =
#' 'FS'}), or parent child (\code{rel = 'PC'}) pairs
#' @param save Write the results to disk if \code{TRUE}
#' @param strPath Optional prefix to add to the results file path so that the
#' output location can be specified
#' @param strVer Optional suffix for the results file. This is useful when
#' running multiple instances of R
#' @param BlockSize Sets the number of random profiles to be generated in each
#' iteration. By default the block size is set to 1 percent of the total sample
#' size. It is unclear whether the procedure is more efficient if a bigger
#' percentage of the total is used. Users must take care to make sure that the
#' block size evenly divides \code{N} otherwise the procedure will exit
#' @param fileName This argument lets the user override the default result file
#' naming scheme
#' @return a data frame with three columns: sib, pc, ibs containing the LRs for
#' full-siblings, parent-child, and the number of matching alleles for each
#' generated pair of profiles.
#' @author James M. Curran
#' @seealso readResults, errorRate
#' @examples
#' 
#' ## not run
#' ## this replicates Ge et al.'s experiment and takes about 45 minutes
#' ## to run (I think)
#' \dontrun{
#' data(fbiCaucs)
#' N = 1000000
#' sim(N, fbiCaucs, save = T)
#' sim(N, fbiCaucs, 'FS', save = T)
#' sim(N, fbiCaucs, 'PC', save = T)
#' }
#'
#' @importFrom utils setTxtProgressBar 
#' @export sim
sim = function(N, Freqs, rel = "UN", save = FALSE, strPath = "",
               strVer = "", BlockSize = N/100, fileName = NULL){

    rel = toupper(rel)
    if(!grepl("(UN|FS|PC)", rel, ignore.case = T, perl = TRUE)){
        stop("Unrecognized relationship. Must be one of 'U', 'FS', 'PC'")
    }


    f1 = NULL

    if(save){
        if(is.null(fileName)){
            if(nchar(strVer) > 0){
                fileName = paste(strPath, sprintf("results-sim-%s-%d-%s.csv.gz",
                                 rel, N, strVer),sep = "")
            }else{
                fileName = paste(strPath, sprintf("results-sim-%s-%d.csv.gz",
                                 rel, N),sep = "")
            }
        }

        f1 = gzfile(fileName, 'w')
        if(!isOpen(f1)){
            stop(paste("Couldn't open", fileName, "for writing"))
        }
    }

    ## make space for results
    lrsibs = rep(0, N)
    lrpc = lrsibs
    ibs = lrsibs

    ## progressBar

    pb = txtProgressBar(min = 0, max = 100, style = 3)
    step = N/100

    ## number of blocks (usually 100), but...
    numBlocks = N/BlockSize

    if(numBlocks - floor(numBlocks) > 0){
        stop("BlockSize must be a whole integer factor of sample size")
    }


    j = 0
    k = 0
    setTxtProgressBar(pb, k)

    ## predetermine some inputs
    nLoci = length(Freqs$loci)
    f = unlist(Freqs$freqs)
    n = sapply(Freqs$freqs, length)

    if(rel == "UN"){ ## put this out here so that it is faster
        for(block in 1:numBlocks){
            P = randomProfilePairs(Freqs, BlockSize)

            for(i in 1:BlockSize){
                i1 = (block-1)*BlockSize + i

                p1 = P[[i]]$prof1
                p2 = P[[i]]$prof2
                ibs[i1] = IBS(p1, p2, nLoci)
                lrsibs[i1] = lrSib(p1,p2, NULL, nLoci, f = f, n = n)
                lrpc[i1] = lrPC(p1, p2, NULL, nLoci, f = f, n = n)

                if(j == step){
                    k = k + 1
                    setTxtProgressBar(pb, k)
                    j = 0
                }

                j = j + 1

                if(save){ ## bit expensive to check every iteration
                    line = sprintf("%10.8E,%10.8E,%d", lrsibs[i1],
                                                       lrpc[i1], ibs[i1])
                    writeLines(line, f1)
                }
            }
         }
    }else if(rel == "FS"){
        for(block in 1:numBlocks){
            P = randomSibPairs(Freqs, BlockSize)

            for(i in 1:BlockSize){
                i1 = (block-1)*BlockSize + i

                p1 = P[[i]]$sib1
                p2 = P[[i]]$sib2
                ibs[i1] = IBS(p1, p2, nLoci)
                lrsibs[i1] = lrSib(p1,p2, NULL, nLoci, f = f, n = n)
                lrpc[i1] = lrPC(p1, p2, NULL, nLoci, f = f, n = n)

                if(j == step){
                    k = k + 1
                    setTxtProgressBar(pb, k)
                    j = 0
                }

                j = j + 1

                if(save){ ## bit expensive to check every iteration
                    line = sprintf("%10.8E,%10.8E,%d", lrsibs[i1],
                                                       lrpc[i1], ibs[i1])
                    writeLines(line, f1)
                }
            }
         }
    }else if(rel == "PC"){
        for(block in 1:numBlocks){
            P = randomPCPairs(Freqs, BlockSize)

            for(i in 1:BlockSize){
                i1 = (block-1)*BlockSize + i

                p1 = P[[i]]$parent
                p2 = P[[i]]$child
                ibs[i1] = IBS(p1, p2, nLoci)
                lrsibs[i1] = lrSib(p1,p2, NULL, nLoci, f = f, n = n)
                lrpc[i1] = lrPC(p1, p2, NULL, nLoci, f = f, n = n)

                if(j == step){
                    k = k + 1
                    setTxtProgressBar(pb, k)
                    j = 0
                }

                j = j + 1

                if(save){ ## bit expensive to check every iteration
                    line = sprintf("%10.8E,%10.8E,%d", lrsibs[i1],
                                                       lrpc[i1], ibs[i1])
                    writeLines(line, f1)
                }
            }
        }
    }

    if(save){
        close(f1)
        cat(paste("Results written to", fileName, "\n"))
    }

    invisible(data.frame(sib = lrsibs, pc = lrpc, ibs = ibs))
}



