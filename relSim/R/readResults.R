#' Read a simulation result set from file
#' 
#' This function will read the output from \code{sim} that has been saved to
#' disk
#' 
#' The arguments to this file are used to generate the input file name. The
#' format is very rigid, being 'results-sim-rel-N(-strVer).csv(.gz)' That is,
#' if strVer is something than an empty string then it is included after the
#' number of interations. Similarly if \code{gzip == TRUE} then the filename is
#' assumed to end with '.gz'
#' 
#' @param N The number of iterations in the simulation
#' @param rel 'UN' = unrelated, 'FS' = full-sib, 'PC' = parent-child
#' @param gzip If \code{TRUE} then it is assumed that the file is compressed
#' @param strPath Optional location of files. Must terminate with / otherwise
#' it will not work
#' @param strVer A version string, useful if more than simulation has been run
#' @param fileName This argument allows the user to override the default file
#' naming conventions of the result file
#' @return a data frame with three columns labelled sib, pc, and ibs. These
#' represent the LRs for sibs and parent-child calculated on each simulated
#' profile pair, and the number of matching alleles (IBS).
#' @author James M. Curran
#' @seealso sim
#' @examples
#' 
#' data(fbiCaucs)
#' ## not run
#' ## write the results of 100 unrelated profile pairs to
#' ## results-sim-UN-100.csv.gz
#' ## and read it back in
#' \dontrun{
#' sim(100, save = T)
#' unrel = readResults(100)
#' sim(100, rel = "FS", strVer = "01", save = T)
#' sibs = readResults(100, rel = "FS", strVer = "01")
#' }
#' 
#' @export readResults
readResults = function(N = 0, rel = "UN", gzip = TRUE,
                       strPath = "", strVer = "", fileName = NULL){

    if(is.null(fileName)){
        if(nchar(strVer)==0){
            fileName = sprintf("results-sim-%s-%d.csv", rel, N)
        }else{
            fileName = sprintf("results-sim-%s-%d-%s.csv", rel, N, strVer)
        }
    }

    if(nchar(strPath)>0){
        fileName = paste(strPath, fileName, sep = "")
    }

    if(gzip & !grepl("gz$", fileName))
        fileName = paste(fileName, ".gz", sep = "")

    res = matrix(scan(fileName, sep = ','), ncol = 3, byrow = T)

    return(data.frame(sib = res[,1], pc = res[,2], ibs = res[,3]))
}
