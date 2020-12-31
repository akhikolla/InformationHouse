#' Convert Frequency Table into Raw Data
#' 
#' Data with unique rows and a frequency column is converted into
#' data with duplicate rows.
#' 
#' The ouput data frame can be used as input to \code{\link{bbl}}.
#' 
#' @param data Data frame with factors in columns 
#' @param freq Vector of frequency of each row in \code{data}
#' @return Data frame with one row per instances
#' @examples
#' Titanic
#' x <- as.data.frame(Titanic)
#' head(x)
#' titanic <- freq2raw(data=x[,1:3], freq=x$Freq)
#' head(titanic)
#' @export
freq2raw <- function(data,freq){
  
  if(length(freq)!=NROW(data)) 
     stop('Frequency length does not match data')
  n <- nrow(data)
  dat <- NULL
  for(i in 1:n){
    w <- data[rep(i,freq[i]),]
    dat <- rbind(dat, w)
  }
  rownames(dat) <- seq_len(NROW(dat))
  return(dat)
}


#' Read FASTA File
#' 
#' Read nucleotide sequence files in FASTA format
#' 
#' Sequence data in FASTA files are converted into data frame
#' suitable as input to \code{\link{bbl}}. If sequence lengths are different,
#' instances longer than those already read will be truncated. Empty sequences
#' are skipped.
#' 
#' @param file File name of FASTA input.
#' @param rownames Use the sequence annotation line in file (starts with
#'        \code{'>'}) as the row names. Will fail if there are duplicate items.
#' @return Data frame of each sequence in rows.
#' @examples
#' file <- tempfile('data')
#' write('>seq1', file)
#' write('atgcc', file, append=TRUE)
#' write('>seq2', file, append=TRUE)
#' write('gccaa', file, append=TRUE)
#' system(paste0('cat ',file))
#' x <- readFasta(file)
#' x
#' @export
readFasta <- function(file, rownames=FALSE){
  
  if(!file.exists(file)) stop(paste0(file,' does not exist'))
  fl <- readLines(file)
  if(length(fl)==0) stop(paste0(file, ' is empty'))
  label <- dat <- NULL
  
  for(i in seq_along(fl)){
    if(substr(fl[[i]],start=1,stop=1)=='>'){
      x <- strsplit(fl[[i]],split='\t')[[1]]
      x <- paste(x,collapse=' ')
      label <- c(label, x)
    }
    else{
      seqs <- strsplit(fl[[i]],split='')[[1]]
      if(NROW(dat)>0){ 
        if(length(seqs)<NCOL(dat)) next()
        if(length(seqs)>NCOL(dat)) seqs <- seqs[seq_len(NCOL(dat))]
      }
      dat <- rbind(dat, as.data.frame(matrix(seqs,nrow=1)))
    }
  }
  if(rownames){
    if(length(label)!=NROW(dat)) 
      stop('Annotation lines do not match sequences.')
    rownames(dat) <- label
  }
  colnames(dat) <- seq_len(NCOL(dat))
  
  return(dat)
}

#' Remove Non-varying Predictors
#' 
#' Constant predictor is identified and removed
#' 
#' Variables with only one factor level is removed from data. Intended
#' for use before calling \code{\link{bbl}}.
#' @param x Data frame containing discrete factor variables in each column
#' @return Data frame omitting non-varying variables from \code{x}.
#' @examples
#' set.seed(351)
#' nt <- c('a','c','g','t')
#' x <- data.frame(v1=sample(nt,size=50,replace=TRUE),
#'                 v2=rep('a',50),v3=sample(nt,size=50,replace=TRUE))
#' y <- sample(c('case','ctrl'),size=50,replace=TRUE)
#' dat <- cbind(data.frame(y=y), x)
#' summary(dat)
#' dat <- removeConst(dat)
#' summary(dat)
#' @export
removeConst <- function(x){
  
  bad <- lapply(x, function(x){length(levels(factor(x)))==1})
  bad <- unlist(bad)
  return(x[,!bad])
}
