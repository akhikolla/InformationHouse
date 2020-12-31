freq_to_unit <- function(freq_distr){
	
#' Individual rankings/orderings from the frequency distribution
#' 
#' Construct the dataset of individual rankings/orderings from the frequency distribution of the distinct observed sequences.
#' 
#' @param freq_distr Numeric matrix of the distinct observed sequences with the corresponding frequencies indicated in the last \eqn{(K+1)}-th column. 
#' 
#' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of observed individual sequences.
#' 
#' @author Cristina Mollica and Luca Tardella
#' @examples
#' 
#' library(gtools)
#' K <- 4
#' perm_matrix <- permutations(n=K, r=K)
#' freq_data <- cbind(perm_matrix, sample(1:factorial(K)))
#' freq_data
#' freq_to_unit(freq_distr=freq_data)
#' 
#' @export 
  
  if(is.vector(freq_distr)){
    freq_distr <- t(freq_distr)
  }
  
	K <- ncol(freq_distr)-1
	
	r_seq <- fill_single_entries(data=freq_distr[,-(K+1)])
	out <- r_seq[rep(1:nrow(r_seq),freq_distr[,(K+1)]),]
	rownames(out) <- NULL
	return(out)
	
	######### TUTTE LE DIRETTIVE PER CREARE IL FILE NAMESPACE 
	######### LE INSERIAMO QUI SOTTO 
	
	#'@useDynLib PLMIX, .registration = TRUE
	#'@importFrom stats median 
	#'@importFrom stats var
	#'@importFrom stats rgamma
	#'@importFrom stats dgamma
	#'@importFrom stats na.omit
	#'@importFrom utils getFromNamespace
	#'@importFrom abind adrop
	#'@importFrom coda as.mcmc
	#'@importFrom coda HPDinterval
	#'@importFrom foreach foreach
	#'@importFrom foreach %dopar%
	#'@importFrom graphics plot
	#'@importFrom gtools ddirichlet
	#'@importFrom gtools permutations 
	#'@importFrom gridExtra grid.arrange
	#'@importFrom ggmcmc ggmcmc
	#'@importFrom ggmcmc ggs
	#'@importFrom ggplot2 ggplot
	#'@importFrom ggplot2 aes
	#'@importFrom ggplot2 aes_string
	#'@importFrom ggplot2 position_stack
	#'@importFrom ggplot2 geom_bar
	#'@importFrom ggplot2 coord_polar
	#'@importFrom ggplot2 labs
	#'@importFrom ggplot2 geom_tile
	#'@importFrom ggplot2 scale_fill_brewer
	#'@importFrom ggplot2 theme_void
	#'@importFrom ggplot2 xlim
	#'@importFrom ggplot2 scale_fill_gradient
	#'@importFrom ggplot2 element_blank
	#'@importFrom ggplot2 theme
	#'@importFrom label.switching pra
	#'@importFrom label.switching permute.mcmc
	#'@importFrom MCMCpack rdirichlet
	#'@importFrom reshape2 melt
	#'@importFrom rcdd makeH
	#'@importFrom rcdd scdd 
	#'@importFrom radarchart chartJSRadar
	#'@importFrom Rcpp evalCpp
	#'
	
	
}

unit_to_freq <- function(data){

#' Frequency distribution from the individual rankings/orderings
#' 
#' Construct the frequency distribution of the distinct observed sequences from the dataset of individual rankings/orderings.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of observed individual sequences.
#' @return Numeric matrix of the distinct observed sequences with the corresponding frequencies indicated in the last \eqn{(K+1)}-th column. 
#' 
#' @author Cristina Mollica and Luca Tardella
#' @examples
#'
#' ## Frequency distribution for the APA top-ordering dataset
#' data(d_apa)
#' unit_to_freq(data=d_apa)
#' @export 

data <- fill_single_entries(data=data)
  
K <- ncol(data)
freq <- table(apply(data,1,paste,collapse="-"))

obs_seq <- matrix(as.numeric(unlist(strsplit(names(freq),split="-"))),nrow=length(freq),ncol=K,byrow=TRUE)
rownames(obs_seq) <- NULL
out <- cbind(obs_seq,freq=freq,deparse.level=0)
rownames(out) <- NULL
return(out)
}

fill_single_entries <- function(data){
  
  #/' Utility to fill in single missing entries of top-(K-1) sequences in partial ordering/ranking datasets
  #/'
  #/' @param data Numeric data matrix of partial sequences.
  #/'
  #/' @return Numeric data matrix of partial sequences in the same format of the input \code{data} with possible single missing entries filled.
  #/' @author Cristina Mollica and Luca Tardella

  if(is.vector(data)){
    data <- t(data)
  }
  
  K=ncol(data)
  r_single_miss <- (rowSums(data==0)==1)
  
  if(any(r_single_miss)){  
    w_row <- which(r_single_miss)
    w_col <- apply(data[w_row,,drop=FALSE],1,function(x)which(x==0))
    w_item <- apply(data[w_row,,drop=FALSE],1,setdiff,x=1:K)
    data[cbind(w_row,w_col)] <- w_item
    warning(paste(paste0("Top-",K-1,""),"sequencies correspond to full orderings. Single missing entries filled."), call. = FALSE)
  }
  
  return(data)
  
}

is.top_ordering <- function(data,...){
  
  #' Top-ordering datasets
  #' 
  #' Check the consistency of partial ordering data with a top-ordering dataset.
  #' 
  #' The argument \code{data} requires the partial sequences expressed in ordering format. When the value of \code{is.top-ordering} is \code{FALSE}, the membership function returns also a message with the conditions that are not met for the \code{data} to be a top-ordering dataset. \code{NA}'s in the input \code{data} are tacitly converted into zero entries. 
  #' 
  #' @param data An object containing the partial orderings whose consistency with a top-ordering dataset has to be tested. The following classes are admissible for \code{data}: numeric \code{matrix}, \code{data.frame}, \code{RandData} from the \code{rankdist} package and \code{rankings} from the \code{PlackettLuce} package.
  #' @param ... Further arguments passed to or from other methods (not used).
  #' 
  #' @return Logical: \code{TRUE} if the \code{data} argument is consistent with a top-ordering dataset (with a possible warning message if the supplied data need a further treatment with the coercion function \code{\link{as.top_ordering}} before being processed with the core functions of \pkg{PLMIX}) and \code{FALSE} otherwise.
  #' 
  #' @references 
  #' Turner, H., Kormidis, I. and Firth, D. (2018). PlackettLuce: Plackett-Luce Models for Rankings. R package version 0.2-3. \url{https://CRAN.R-project.org/package=PlackettLuce}
  #' 
  #' Qian, Z. (2018). rankdist: Distance Based Ranking Models. R package version 1.1.3. \url{https://CRAN.R-project.org/package=rankdist}
  #' 
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @seealso \code{\link[PlackettLuce]{rankings}} and \code{\link[PlackettLuce]{rankings}}
  #' 
  #' @examples
  #' 
  #' ## A toy example of data matrix not satisfying the conditions to be a top-ordering dataset
  #' toy_data=rbind(1:5,
  #' c(0,4,3,2,1),
  #' c(4,3.4,2,1,5),
  #' c(2,3,0,0,NA),
  #' c(4,4,3,2,5),
  #' c(3,5,4,2,6),
  #' c(2,-3,1,4,5),
  #' c(2,0,1,4,5),
  #' c(2,3,1,1,1),
  #' c(2,3,0,4,0))
  #' 
  #' is.top_ordering(data=toy_data) 
  #' 
  #' ## A dataset from the StatRank package satisfying the conditions to be a top-ordering dataset
  #' library(StatRank) 
  #' data(Data.Election9)
  #' is.top_ordering(data=Data.Election9)
  #' 
  #' @export is.top_ordering
  #' @export 
  
  if(!(class(data)[1]%in%c("top_ordering","RankData","rankings","matrix","data.frame"))){
    stop("Invalid 'type' of data argument.")
  }
  
  if(any(class(data)=="top_ordering")){
    out=TRUE
  }
  
  if(class(data)[1]=="RankData"){
    warning("Objects of class 'RankData' are compatible with top-ordering datasets, but need to be coerced with as.top_ordering() before using the other functions of the PLMIX package.")
    out=TRUE
  }
  
  if(class(data)[1]=="rankings"){
    ttt=try(as.top_ordering(data=data, format_input="ranking"),silent=TRUE)
    if(class(ttt)=="try-error"){
      warning("The supplied data of class 'rankings' is not compatible with a top-ordering dataset because all rankings contain ties.")
      out=FALSE
    }else{
      warning("The supplied data of class 'rankings' is compatible with a top-ordering dataset, but needs to be coerced with as.top_ordering() before using the other functions of the PLMIX package.")
      out=TRUE
    }
  }
  
  if(class(data)[1]=="matrix" | class(data)[1]=="data.frame"){
    if(class(data)[1]=="data.frame"){
      data <- as.matrix(data)
    }
    if(is.vector(data)){
      data <- t(data)
    }
    
    data[which(is.na(data))]=0
    
    K=ncol(data)
    
    if(any(!(data%in%(0:K)))){
      check1=FALSE
      message(paste0("->> Only integers {", paste(0:K,collapse=", "), "} are allowed as entries of the top-ordering dataset:"))
      if(any(data<0)){
        message("* Some entries are negative.")
      }
      if(any(data>K)){
        message(paste0("* Some entries are > ", K,"."))
      }
      if(any((data%%1)>0)){
        message(paste("* Some entries are not integer."))
      }
    }else{
      check1=TRUE
    }
    
    data_dupl <- t(apply(data,1,duplicated))
    data_dupl[data==0] <- NA  
    
    if(any(data_dupl,na.rm=TRUE)){
      check2=FALSE
      message("->> Ties are not allowed.")
    }else{
      check2=TRUE
    }
    
    non_cons_zeros=apply(data,1,function(x) if(0%in%x) length(setdiff(min(which(x==0)):K,which(x==0)))>0 else FALSE )
    if(any(non_cons_zeros)){
      check3=FALSE
      message("->> Non-consecutive zero are not allowed in the rows of a top-ordering dataset.")
    }else{
      check3=TRUE
    }
    
    if(any(data[,1]==0)){
      check4=FALSE
      message("->> Rows starting with zero entries are not allowed in a top-ordering dataset.")
    }else{
      check4=TRUE
    }
    
    out=all(c(check1,check2,check3,check4))
    
  }   
  
  
  return(out)
  
}

as.top_ordering <- function(data,format_input=NULL,aggr=NULL,freq_col=NULL,ties_method="random",...){
  
  #' Coercion into top-ordering datasets
  #' 
  #' Attempt to coerce the input data into a top-ordering dataset.
  #' 
  #' The coercion function \code{as.top_ordering} tries to coerce the input data into an object of class \code{top_ordering} after checking for possible partial sequences that do not satisfy the top-ordering requirements. If none of the supplied sequences satisfies the top-ordering conditions, an error message is returned. \code{NA}'s in the input \code{data} are tacitly converted into zero entries.
  #' 
  #' @param data An object containing the partial sequences to be coerced into an object of class \code{top_ordering}. The following classes are admissible for \code{data}: numeric \code{matrix}, \code{data.frame}, \code{RandData} from the \code{rankdist} package and \code{rankings} from the \code{PlackettLuce} package.
  #' @param format_input Character string indicating the format of the \code{data} input, namely \code{"ordering"} or \code{"ranking"}. Used only when the class of the \code{data} argument is matrix or data frame. Default is \code{NULL}.
  #' @param aggr Logical: whether the \code{data} argument collects the distinct observed sequences with the corresponding frequencies (aggregated format). Used only when the class of the \code{data} aargument is matrix or data frame. Default is \code{NULL}.
  #' @param freq_col Integer indicating the column of the \code{data} argument containing the frequencies of the distinct observed sequences. Used only when the class of the \code{data} argument is matrix or data frame and \code{aggr} argument is \code{TRUE}. Default is \code{NULL}.
  #' @param ties_method Character string indicating the treatment of sequences with ties (not used for data of class \code{RankData}). If \code{"remove"}, the sequences with ties are removed before acting the coercion; if \code{"random"} (default), tied positions are re-assigned at random before acting the coercion.
  #' @param ... Further arguments passed to or from other methods (not used).
  #' 
  #' @return An object of S3 class \code{c("top_ordering","matrix")}.
  #' 
  #' @references 
  #' Turner, H., Kormidis, I. and Firth, D. (2018). PlackettLuce: Plackett-Luce Models for Rankings. R package version 0.2-3. \url{https://CRAN.R-project.org/package=PlackettLuce}
  #' 
  #' Qian, Z. (2018). rankdist: Distance Based Ranking Models. R package version 1.1.3. \url{https://CRAN.R-project.org/package=rankdist}
  #' 
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @seealso \code{\link{is.top_ordering}}, \code{\link[PlackettLuce]{as.rankings}} and \code{\link[PlackettLuce]{rankings}}
  #' 
  #' @examples
  #' 
  #' ## Coerce an object of class 'rankings' into an object of class 'top_ordering'
  #' library(PlackettLuce)
  #' RR <- matrix(c(1, 2, 0, 0,
  #' 4, 1, 2, 3,
  #' 2, 1, 1, 1,
  #' 1, 2, 3, 0,
  #' 2, 1, 1, 0,
  #' 1, 0, 3, 2), nrow = 6, byrow = TRUE) 
  #' RR_rank=as.rankings(RR)
  #' RR_rank
  #' as.top_ordering(RR_rank, ties_method="random")
  #' 
  #' ## Coerce an object of class 'RankData' into an object of class 'top_ordering'
  #' library(rankdist) 
  #' data(apa_partial_obj)
  #' d_apa_top_ord=as.top_ordering(data=apa_partial_obj)
  #' identical(d_apa,d_apa_top_ord)
  #' 
  #' ## Coerce a data frame from the package prefmod into an object of class 'top_ordering'
  #' library(prefmod)
  #' data(carconf)
  #' carconf_rank=carconf[,1:6]
  #' carconf_top_ord=as.top_ordering(data=carconf_rank,format_input="ranking",aggr=FALSE)
  #' identical(d_carconf,carconf_top_ord)
  #' 
  #' ## Coerce a data frame from the package pmr into an object of class 'top_ordering'
  #' library(pmr)
  #' data(big4)
  #' head(big4)
  #' big4_top_ord=as.top_ordering(data=big4,format_input="ranking",aggr=TRUE,freq_col=5)
  #' head(big4_top_ord)
  #' 
  #' @export as.top_ordering
  #' @export 
  
  
  
  if(!(class(data)[1]%in%c("top_ordering","RankData","rankings","matrix","data.frame"))){
    stop("Invalid 'type' of data argument (see 'Details').")
  }
  
  if(any(class(data)=="top_ordering")){
    out=data
  }
  
  if(class(data)[1]=="RankData"){
    K=data@nobj                
    dist_rankings=data@ranking
    tied_rows=apply(dist_rankings,1,function(x)any(duplicated(x)))
    dist_rankings[tied_rows,]=t(apply(dist_rankings[tied_rows,],1,rank,ties.method="max"))
    temp_tied=dist_rankings[tied_rows,]
    temp_tied[which(temp_tied==K)]=0
    dist_rankings[tied_rows,]=temp_tied
    dist_orderings=rank_ord_switch(data=dist_rankings,format_input="ranking")
    n_dist=data@ndistinct
    out=dist_orderings[rep(1:n_dist,times=data@count),]
  }
  
  if(class(data)[1]=="rankings"){
    temp_rankings=unclass(data)          
    N=nrow(temp_rankings)
    temp_rankings[temp_rankings==0] <- NA  
    tied_rows=which(apply(temp_rankings,1,function(x)any(duplicated(na.omit(x)))))
    if(length(tied_rows)>0){
      if(ties_method=="remove"){
        if(length(tied_rows)<N){
          warning("Rankings with ties are removed from the supplied dataset.")
          temp_rankings=temp_rankings[-tied_rows,,drop=FALSE]
        }else{
          if(length(tied_rows)==N){
            stop("Supplied data cannot be coerced into a top-ordering dataset because all rankings contain ties.")
          } 
        }
      }else{
        warning("Tied positions are re-assigned at random to satisfy the top-ordering requirements.")
        temp_rankings[tied_rows,]=t(apply(temp_rankings[tied_rows,],1,rank,na.last="keep",ties.method="random"))
      }
    }
    temp_rankings[is.na(temp_rankings)] <- 0  
    out=rank_ord_switch(temp_rankings,format_input="ranking")
  }
  
  if(class(data)[1]=="matrix" | class(data)[1]=="data.frame"){
    if(class(data)[1]=="data.frame"){
      data <- as.matrix(data)
    }
    if(is.vector(data)){
      data <- t(data)
    }
    if(aggr){
      NN=nrow(data)
      data_aggr=data[,-freq_col,drop=FALSE]
      freq=data[,freq_col]
      data=data_aggr[rep(1:NN,times=freq),]
    }
    K=ncol(data)
    N=nrow(data)
    data[which(is.na(data))]=0
    if(format_input=="ordering"){
      
      check1=which(apply(data,1,function(x)any(!(x%in%(0:K)))))
      
      data_dupl <- t(apply(data,1,duplicated))
      data_dupl[data==0] <- NA
      check2=which(apply(data_dupl,1,any,na.rm=TRUE))
      
      non_cons_zeros=apply(data,1,function(x) if(0%in%x) length(setdiff(min(which(x==0)):K,which(x==0)))>0 else FALSE )
      check3=which(non_cons_zeros)
      
      check4=which(data[,1]==0)
      
      checks=unique(c(check1,check2,check3,check4))
      
      if(length(checks)>0 & length(checks)<N){
        warning("Rows not satisfying the requirements of a top-ordering dataset have been removed. Please, apply the function is.top_ordering() to the supplied data for more information.")
        data=data[-checks,,drop=FALSE]
      }else{
        if(length(checks)==N){
          stop("Supplied data cannot be coerced because the provided orderings do not satisfy the requirements of a top-ordering dataset.")
        }
      }  
    }else{
      
      check1=which(apply(data,1,function(x)any(!(x%in%(0:K)))))
      
      data_dupl <- t(apply(data,1,duplicated))
      data_dupl[data==0] <- NA
      check2=which(apply(data_dupl,1,any,na.rm=TRUE))
      
      check3=which(apply(data,1,function(x)any(diff(sort(x))>1)))
      
      if(ties_method=="remove"){
        
        checks=unique(c(check1,check2,check3))
        
      }else{
        
        checks=unique(c(check1,check3))
        if(length(check2)>0){
          warning("Tied positions are re-assigned at random to satisfy the top-ordering requirements.")
          temp=data[check2,]
          temp[temp==0]=NA
          data[check2,]=t(apply(temp,1,rank,na.last="keep",ties.method="random"))
        }
        
      }
      
      if(length(checks)>0 & length(checks)<N){
        warning("Rows not satisfying the requirements of a top-ordering dataset have been removed. Please, apply the function is.top_ordering() to the supplied data for more information.")
        data=data[-checks,,drop=FALSE]
      }else{
        if(length(checks)==N){
          stop("Supplied data cannot be coerced because the provided rankings do not satisfy the requirements of a top-ordering dataset.")
        }
      }
      
      data=rank_ord_switch(data,format_input="ranking")
    }  
    out=data
  } 
  
  class(out) <- c("top_ordering","matrix")
  return(out)
  
} 



myorder <- function(x){
  
#/' Utility to switch from a partial ranking to a partial ordering (missing positions denoted with zero)
#/' @param x Numeric integer vector
#/' 
#/' @author Cristina Mollica and Luca Tardella

  k <- sum(is.na(x))
  out <- c(order(x,na.last=NA),rep(0,k))
  return(out)
}


rank_ord_switch <- function(data,format_input,nranked=NULL){

#' Switch from orderings to rankings and vice versa
#' 
#' Convert the format of the input dataset from orderings to rankings and vice versa.
#' 
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences whose format has to be converted.
#' @param format_input Character string indicating the format of the \code{data} input, namely \code{"ordering"} or \code{"ranking"}.
#' @param nranked Optional numeric vector of length \eqn{N} with the number of items ranked by each sample unit. 
#' 
#' @return Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences with inverse format.
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' ## From orderings to rankings for the Dublin West dataset
#' data(d_dublinwest)
#' head(d_dublinwest)
#' rank_ord_switch(data=head(d_dublinwest), format_input="ordering")
#' @export 

    
  data <- fill_single_entries(data=data)
     
  K <- ncol(data)
  
  if(any(data==0)){
    	
    	data[data==0] <- NA
    	
    	if(format_input=="ranking"){
    		out <- t(apply(data,1,myorder))
    		colnames(out) <- paste0("Rank_",1:K)
    	}else{
    		N <- nrow(data)
    		if(is.null(nranked)) nranked=rowSums(!is.na(data))
    		out <- matrix(0,nrow=N,ncol=K)
    		out[cbind(rep(1:N,nranked),na.omit(c(t(data))))] <- unlist(sapply(nranked,seq,from=1))
    	}
    	
    }else{ 

    out <- t(apply(data,1,order))
    
    }

    if(format_input=="ranking"){
      colnames(out) <- paste0("Rank_",1:K)
    }else{
      colnames(out) <- paste0("Item_",1:K)
    }
    
    return(out)
}


rank_summaries <- function(data,format_input,mean_rank=TRUE,marginals=TRUE,pc=TRUE){

#' Descriptive summaries for a partial ordering/ranking dataset
#' 
#' Compute rank summaries and censoring patterns for a partial ordering/ranking dataset.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences.
#' @param format_input Character string indicating the format of the \code{data} input, namely \code{"ordering"} or \code{"ranking"}.
#' @param mean_rank Logical: whether the mean rank vector has to be computed. Default is \code{TRUE}.
#' @param marginals Logical: whether the marginal rank distributions have to be computed. Default is \code{TRUE}.
#' @param pc Logical: whether the paired comparison matrix has to be computed. Default is \code{TRUE}.
#' 
#' @return A list of named objects:
#' 
#'  \item{\code{nranked}}{ Numeric vector of length \eqn{N} with the number of items ranked by each sample unit.}
#'  \item{\code{nranked_distr}}{ Frequency distribution of the \code{nranked} vector.}
#'  \item{\code{na_or_not}}{ Numeric \eqn{3}\eqn{\times}{x}\eqn{K} matrix with the counts of sample units that ranked or not each item. The last row contains the total by column, corresponding to the sample size \eqn{N}.}
#'  \item{\code{mean_rank}}{ Numeric vector of length \eqn{K} with the mean rank of each item.}
#'  \item{\code{marginals}}{ Numeric \eqn{K}\eqn{\times}{x}\eqn{K} matrix of the marginal rank distributions: the \eqn{(i,j)}-th entry indicates the number of units that ranked item \eqn{i} in the \eqn{j}-th position.}
#'  \item{\code{pc}}{ Numeric \eqn{K}\eqn{\times}{x}\eqn{K} paired comparison matrix: the \eqn{(i,i')}-th entry indicates the number of sample units that preferred item \eqn{i} to item \eqn{i'}.}
#' 
#' 
#' @references 
#' Marden, J. I. (1995). Analyzing and modeling rank data. \emph{Monographs on Statistics and Applied Probability} (64). Chapman & Hall, ISSN: 0-412-99521-2. London.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#' 
#' data(d_carconf)
#' rank_summaries(data=d_carconf, format_input="ordering")
#' @export 

data <- fill_single_entries(data=data)
N <- nrow(data)
K <- ncol(data)	
if(format_input=="ordering"){
	 data <- rank_ord_switch(data=data,format_input=format_input,nranked=NULL)
	 format_input <- "ranking"
}
data[data==0] <- NA
isna_data <- is.na(data)
nranked <- rowSums(!isna_data) 
#nranked_distr <- table(nranked,dnn=NULL,deparse.level=0) 
nranked_distr <- table(factor(nranked,levels=1:K)) 
#names(nranked_distr) <- paste0("Top-",1:(K-1))
names(nranked_distr) <- paste0("Top-",names(nranked_distr))
na <- colSums(isna_data) 
na_or_not <- rbind(na, N-na, rep(N, K))
dimnames(na_or_not) <- list(c("n.a.","not n.a.","total"),paste0("Item_",1:K))
if(mean_rank){
    mean_rank <- colMeans(data,na.rm=TRUE)  
    names(mean_rank) <- paste0("Item_",1:K)
}else{
	mean_rank <- NULL
}	
if(marginals){
    marginals <- apply(data,2,tabulate,nbins=K)
    dimnames(marginals) <- list(paste0("Rank_",1:K),paste0("Item_",1:K))
}else{
	marginals <- NULL
}	
if(pc){
    data[is.na(data)] <- 0
    pc <- paired_comparisons(data=data,format_input=format_input,nranked=nranked)
    rownames(pc) <- colnames(pc) <- paste0("Item_",1:K)
}else{
	pc <- NULL
}	
out <- list(nranked=nranked,nranked_distr=nranked_distr,
         na_or_not=na_or_not,mean_rank=mean_rank,
         marginals=marginals,pc=pc)
return(out)
}

paired_comparisons <- function(data,format_input,nranked=NULL){

#' Paired comparison matrix for a partial ordering/ranking dataset
#' 
#' Construct the paired comparison matrix for a partial ordering/ranking dataset.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences.
#' @param format_input Character string indicating the format of the \code{data} input, namely \code{"ordering"} or \code{"ranking"}.
#' @param nranked Optional numeric vector of length \eqn{N} with the number of items ranked by each sample unit. 
#' 
#' @return Numeric \eqn{K}\eqn{\times}{x}\eqn{K} paired comparison matrix: the \eqn{(i,i')}-th entry indicates the number of sample units that preferred item \eqn{i} to item \eqn{i'}.
#' 
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#' 
#' @author Cristina Mollica and Luca Tardella
#'
#' @seealso \code{\link{rank_summaries}} 
#' 
#' @examples
#' 
#' data(d_dublinwest)
#' paired_comparisons(data=d_dublinwest, format_input="ordering")
#' @export 
	
  data <- fill_single_entries(data=data)
	N <- nrow(data)
	K <- ncol(data)
    
    if(format_input=="ranking"){
    	if(is.null(nranked)) nranked <- rowSums(data!=0)
    	data <- rank_ord_switch(data,format_input=format_input,nranked=nranked)
    } 
    pc <- tau(pi_inv=data)
    rownames(pc) <- colnames(pc) <- paste0("Item_",1:K)
    return(pc)
}   # K*K matrix


make_partial <- function(data,format_input,nranked=NULL,probcens=rep(1,ncol(data)-1)){

#' Censoring of complete rankings/orderings
#' 
#' Return partial top rankings/orderings from complete sequences obtained either with user-specified censoring patterns or with a random truncation.
#' 
#' The censoring of the complete sequences can be performed in: (i) a deterministic way, by specifying the number of top positions to be retained for each sample unit in the \code{nranked} argument; (ii) a random way, by sequentially specifying the probabilities of the top-1, top-2, \eqn{...}, top-\eqn{(K-1)} censoring patterns in the \code{probcens} argument. Recall that a top-\eqn{(K-1)} sequence corresponds to a complete ordering/ranking.
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of complete sequences to be censored.
#' @param format_input Character string indicating the format of the \code{data} input, namely \code{"ordering"} or \code{"ranking"}.
#' @param nranked Numeric vector of length \eqn{N} with the desired number of items ranked by each sample unit after censoring. If not supplied (\code{NULL}), the censoring patterns are randomly generated according to the probabilities in the \code{probcens} argument. 
#' @param probcens Numeric vector of length \eqn{(K-1)} with the probability of each censoring pattern to be employed for the random truncation of the complete sequences (normalization is not necessary). It works only if \code{nranked} argument is \code{NULL} (see 'Details'). Default is equal probabilities.
#' 
#' @return A list of two named objects:
#' 
#'  \item{\code{partialdata}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial (censored) sequences with the same format of the input \code{data} and missing positions/items denoted with zero entries.}	
#'  \item{\code{nranked}}{ Numeric vector of length \eqn{N} with the number of items ranked by each sample unit after censoring.}
#' 
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#' 
#' data(d_german)
#' head(d_german)
#' d_german_cens <- make_partial(data=d_german, format_input="ordering", 
#'                               probcens=c(0.3, 0.3, 0.4))  
#' head(d_german_cens$partialdata)
#' 
#' ## Check consistency with the nominal censoring probabilities
#' round(prop.table(table(d_german_cens$nranked)), 2)
#' 
#' @export 
  
  data <- fill_single_entries(data=data)
  K <- ncol(data)

  if(format_input=="ranking"){
    	data <- rank_ord_switch(data,format_input=format_input)
    } 

  if(is.null(nranked)){
    N <- nrow(data)	
    nranked <- sample(c(1:(K-2),K),size=N,replace=TRUE,prob=probcens)	
  }

  out <- data*t(sapply(nranked,function(x)rep(c(1,0),c(x,K-x))))

  if(format_input=="ranking"){
    	out <- rank_ord_switch(out,format_input="ordering",nranked=nranked)
    } 

  return(list(partialdata=out,nranked=nranked))	
} # N*K censored data matrix

make_complete <- function(data,format_input,nranked=NULL,probitems=rep(1,ncol(data))){

#' Completion of partial rankings/orderings
#' 
#' Return complete rankings/orderings from partial sequences relying on a random generation of the missing positions/items.
#' 
#' The completion of the partial top rankings/orderings is performed according to the Plackett-Luce scheme, that is, with a sampling without replacement of the not-ranked items by using the positive values in the \code{probitems} argument as support parameters (normalization is not necessary).
#' 
#' @param data Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial sequences to be completed.
#' @param format_input Character string indicating the format of the \code{data} input, namely \code{"ordering"} or \code{"ranking"}.
#' @param nranked Optional numeric vector of length \eqn{N} with the number of items ranked by each sample unit. 
#' @param probitems Numeric vector with the \eqn{K} item-specific probabilities to be employed for the random generation of the missing positions/items (see 'Details'). Default is equal probabilities.
#' 
#' @return A list of two named objects:
#' 
#'  \item{\code{completedata}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of complete sequences with the same format of the input \code{data}.}	
#'  \item{\code{nranked}}{ Numeric vector of length \eqn{N} with the number of items ranked by each sample unit of the input \code{data}.}
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @examples
#'
#' ## Completion based on the top item frequencies
#' data(d_dublinwest)
#' head(d_dublinwest)
#' top_item_freq <- rank_summaries(data=d_dublinwest, format_input="ordering", mean_rank=FALSE, 
#'                                 pc=FALSE)$marginals["Rank_1",]
#' 
#' d_dublinwest_compl <- make_complete(data=d_dublinwest, format_input="ordering", 
#'                                     probitems=top_item_freq)
#' head(d_dublinwest_compl$completedata)
#' 
#' @export 
  
  data <- fill_single_entries(data=data)
  K <- ncol(data)

  if(is.null(nranked)){
		nranked <- rowSums(data!=0) 
	}

  if(format_input=="ranking"){
    	data <- rank_ord_switch(data,format_input=format_input,nranked=nranked)
    } 

  data[data==0] <- NA	
  out <- data
  partialdata <- out[which(nranked!=K),]
	
  out[which(nranked!=K),] <- t(apply(partialdata,1,function(x){ notrankeditems=setdiff(1:K,x); c(na.omit(x),sample(notrankeditems,prob=probitems[notrankeditems]))}))

  if(format_input=="ranking"){
    	out <- rank_ord_switch(out,format_input="ordering")
    } 

  return(list(completedata=out,nranked=nranked))	
	
}

### Utility to simulate from a EPL

mysample <- function(support,pr){ 
	 sample(x=support,prob=pr)
}


rPLMIX <- function(n=1,K,G,p=t(matrix(1/K,nrow=K,ncol=G)),ref_order=t(matrix(1:K,nrow=K,ncol=G)),weights=rep(1/G,G),format_output="ordering"){

#' Random sample from a mixture of Plackett-Luce models
#' 
#' Draw a random sample of complete orderings/rankings from a \eqn{G}-component mixture of Plackett-Luce models.
#' 
#' Positive values are required for \code{p} and \code{weights} arguments (normalization is not necessary).
#' 
#' The \code{ref_order} argument accommodates for the more general mixture of Extended Plackett-Luce models (EPL), involving the additional reference order parameters (Mollica and Tardella 2014). A permutation of the first \eqn{K} integers can be specified in each row of the \code{ref_order} argument to generate a sample from a \eqn{G}-component mixture of EPL. Since the Plackett-Luce model is a special instance of the EPL with the reference order equal to the identity permutation \eqn{(1,\dots,K)}, the default value of the \code{ref_order} argument is forward orders. 
#' 
#' @param n Number of observations to be sampled. Default is 1.
#' @param K Number of possible items.
#' @param G Number of mixture components. 
#' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters. Default is equal support parameters (uniform mixture components).
#' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders. Default is forward orders (identity permutations) in each row, corresponding to Plackett-Luce mixture components (see 'Details').
#' @param weights Numeric vector of \eqn{G} mixture weights. Default is equal weights.
#' @param format_output Character string indicating the format of the returned simulated dataset (\code{"ordering"} or \code{"ranking"}). Default is \code{"ordering"}.
#' 
#' @return If \eqn{G=1}, a numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of simulated complete sequences. If \eqn{G>1}, a list of two named objects:
#' 
#'  \item{\code{comp}}{ Numeric vector of \eqn{N} mixture component memberships.}
#'  \item{\code{sim_data}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{K} matrix of simulated complete sequences.}
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' K <- 6
#' G <- 3
#' support_par <- matrix(1:(G*K), nrow=G, ncol=K)
#' weights_par <- c(0.50, 0.25, 0.25)
#' 
#' set.seed(47201)
#' simulated_data <- rPLMIX(n=5, K=K, G=G, p=support_par, weights=weights_par)
#' simulated_data$comp
#' simulated_data$sim_data
#' 
#' @export 
	if(G==1){
		if(is.matrix(p)) p <- c(p)
		if(is.matrix(ref_order)) ref_order <- c(ref_order)
		p_par <- p/sum(p)
		perm_par <- matrix(p_par,nrow=K,ncol=n)
		out <- t(apply(perm_par,2,mysample,support=1:K))
		out <- out[,order(ref_order)]		
		if(format_output=="ranking") out <- rank_ord_switch(out,format_input="ordering",nranked=rep(K,n))
		return(out)
	}else{
		p_par <- p/rowSums(p)
		comp <- sample(x=1:G,size=n,replace=T,prob=weights)
		perm_par <- p[comp,]
		out <- t(apply(perm_par,1,mysample,support=1:K))
		for(g in 1:G){
			out[comp==g,] <- out[comp==g,order(ref_order[g,])]
		}
		if(format_output=="ranking"){
			out <- rank_ord_switch(out,format_input="ordering",nranked=rep(K,n))
    	}
		return(list(comp=comp,sim_data=out))	
    }  
}


likPLMIX <- function(p,ref_order,weights,pi_inv){
#' @rdname loglikelihood
#' @name Loglikelihood
#' @aliases likPLMIX loglikPLMIX loglikelihood Likelihood likelihood
#' @title Likelihood and log-likelihood evaluation for a mixture of Plackett-Luce models
#' 
#' @description Compute either the likelihood or the log-likelihood of the Plackett-Luce mixture model parameters for a partial ordering dataset.
#' @details The \code{ref_order} argument accommodates for the more general mixture of Extended Plackett-Luce models (EPL), involving the additional reference order parameters (Mollica and Tardella 2014). A permutation of the first \eqn{K} integers can be specified in each row of the \code{ref_order} argument. Since the Plackett-Luce model is a special instance of the EPL with the reference order equal to the identity permutation, the \code{ref_order} argument must be a matrix with \eqn{G} rows equal to \eqn{(1,\dots,K)} when dealing with Plackett-Luce mixtures.
#' @param p Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters.
#' @param ref_order Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific reference orders.
#' @param weights Numeric vector of \eqn{G} mixture weights.
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @return Either the likelihood or the log-likelihood value of the Plackett-Luce mixture model parameters for a partial ordering dataset.
#'
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' data(d_apa)
#' K <- ncol(d_apa)
#' G <- 3
#' support_par <- matrix(1:(G*K), nrow=G, ncol=K)
#' weights_par <- c(0.50, 0.25, 0.25)
#' loglikPLMIX(p=support_par, ref_order=matrix(1:K, nrow=G, ncol=K, byrow=TRUE), 
#'             weights=weights_par, pi_inv=d_apa)
#' 
#' @export
  
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  
  pi_inv <- fill_single_entries(data=pi_inv)
  lik <- exp(loglikPLMIX(p,ref_order,weights,pi_inv))
  return(lik)
}



bicPLMIX <- function(max_log_lik,pi_inv,G,ref_known=TRUE,ref_vary=FALSE){
	
#' BIC for the MLE of a mixture of Plackett-Luce models
#' 
#' Compute BIC value for the MLE of a mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' The \code{max_log_lik} and the BIC values can be straightforwardly obtained from the output of the \code{\link{mapPLMIX}} and \code{\link{mapPLMIX_multistart}} functions when the default noninformative priors are adopted in the MAP procedure. So, the \code{bicPLMIX} function is especially useful to compute the BIC value from the output of alternative MLE methods for mixtures of Plackett-Luce models implemented, for example, with other softwares. 
#' 
#' The \code{ref_known} and \code{ref_vary} arguments accommodate for the more general mixture of Extended Plackett-Luce models (EPL), involving the additional reference order parameters (Mollica and Tardella 2014). Since the Plackett-Luce model is a special instance of the EPL with the reference order equal to the identity permutation \eqn{(1,\dots,K)}, the default values of \code{ref_known} and \code{ref_vary} are set equal, respectively, to \code{TRUE} and \code{FALSE}.
#' 
#' @param max_log_lik Maximized log-likelihood value.
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @param G Number of mixture components.
#' @param ref_known Logical: whether the component-specific reference orders are known (not to be estimated). Default is \code{TRUE}.
#' @param ref_vary Logical: whether the reference orders vary across mixture components. Default is \code{FALSE}.
#' 
#' @return A list of two named objects:
#' 
#'  \item{\code{max_log_lik}}{ The \code{max_log_lik} argument.}
#'  \item{\code{bic}}{ BIC value.}	
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{Ann. Statist.}, \bold{6}(2), pages 461--464, ISSN: 0090-5364, DOI: 10.1002/sim.6224.
#'
#' @author Cristina Mollica and Luca Tardella
#'
#' @seealso \code{\link{mapPLMIX}} and \code{\link{mapPLMIX_multistart}} 
#'
#' @examples
#' 
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' MAP_mult <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=3, n_start=2, n_iter=400*3)
#' bicPLMIX(max_log_lik=MAP_mult$mod$max_objective, pi_inv=d_carconf, G=3)$bic
#' 
#' ## Equivalently
#' MAP_mult$mod$bic
#' 
#' @export 
	
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  
  pi_inv <- fill_single_entries(data=pi_inv)
  N <- nrow(pi_inv)					 
	K <- ncol(pi_inv)
	if(!ref_known){
		if(ref_vary){
   	       bic <- -2*max_log_lik+(G*(K-1)+G+(G-1))*log(N)
   	    }else{
   	       bic <- -2*max_log_lik+(G*(K-1)+1+(G-1))*log(N)
        }
	}else{
	    bic <- -2*max_log_lik+(G*(K-1)+(G-1))*log(N)
	}
    return(list(max_log_lik=max_log_lik,bic=bic))
    }


gammamat <- function(u_bin,z_hat){
	 gam <- t(z_hat)%*%u_bin
	 return(gam)
}

binary_group_ind <- function(class,G){

#' Binary group membership matrix
#'
#' Construct the binary group membership matrix from the multinomial classification vector.
#'
#' @param class Numeric vector of class memberships.
#' @param G Number of possible different classes.
#'
#' @return Numeric \code{length(class)}\eqn{\times}{x}\eqn{G} matrix of binary group memberships.
#' @author Cristina Mollica and Luca Tardella
#' @examples
#' 
#' binary_group_ind(class=c(3,1,5), G=6)
#' 
#' @export 

	N <- length(class)
	temp <- (rep(1:G,length(class))==rep(class,each=G))*1
	out <- matrix(temp,nrow=N,ncol=G,byrow=TRUE)
	return(out)
	}  # N*G matrix


##########################################################       
############# EM for MAP estimation #############################


mapPLMIX <- function(pi_inv,K,G,
                    init=list(p=NULL,omega=NULL),n_iter=1000,
                    hyper=list(shape0=matrix(1,nrow=G,ncol=K),rate0=rep(0,G),alpha0=rep(1,G)),  
                    eps=10^(-6),
					          centered_start=FALSE,
                    plot_objective=FALSE){
#' MAP estimation for a Bayesian mixture of Plackett-Luce models
#' 
#' Perform MAP estimation via EM algorithm for a Bayesian mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' Under noninformative (flat) prior setting, the EM algorithm for MAP estimation corresponds to the EMM algorithm described by Gormley and Murphy (2006) to perform frequentist inference. In this case, the MAP solution coincides with the MLE and the output vectors \code{log_lik} and \code{objective} coincide as well. 
#' 
#' The \code{\link{mapPLMIX}} function performs the MAP procedure with a single starting value. To address the issue of local maxima in the posterior distribution, see the \code{\link{mapPLMIX_multistart}} function.
#' 
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param init List of named objects with initialization values: \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters; \code{omega} is a numeric vector of \eqn{G} mixture weights. If starting values are not supplied (\code{NULL}), they are randomly generated with a uniform distribution. Default is \code{NULL}. 
#' @param n_iter Maximum number of EM iterations.
#' @param hyper List of named objects with hyperparameter values for the conjugate prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is noninformative (flat) prior setting.
#' @param eps Tolerance value for the convergence criterion.
#' @param centered_start Logical: whether a random start whose support parameters and weights should be centered around the observed relative frequency that each item has been ranked top. Default is \code{FALSE}. Ignored when \code{init} is not \code{NULL}.
#' @param plot_objective Logical: whether the objective function (that is the kernel of the log-posterior distribution) should be plotted. Default is \code{FALSE}.
#' 
#' @return A list of S3 class \code{mpPLMIX} with named elements:
#' 
#'  \item{\code{W_map}}{ Numeric vector with the MAP estimates of the \eqn{G} mixture weights.}
#'  \item{\code{P_map}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters.}
#'  \item{\code{z_hat}}{ Numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of estimated posterior component membership probabilities.}
#'  \item{\code{class_map}}{ Numeric vector of \eqn{N} mixture component memberships based on MAP allocation from the \code{z_hat} matrix.}
#'  \item{\code{log_lik}}{ Numeric vector of the log-likelihood values at each iteration.}
#'  \item{\code{objective}}{ Numeric vector of the objective function values (that is the kernel of the log-posterior distribution) at each iteration.}
#'  \item{\code{max_objective}}{ Maximized objective function value.}
#'  \item{\code{bic}}{ BIC value (only for the default flat priors, otherwise \code{NULL}).}
#'  \item{\code{conv}}{ Binary convergence indicator: 1 = convergence has been achieved, 0 = otherwise.}
#'  \item{\code{call}}{ The matched call.}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Gormley, I. C. and Murphy, T. B. (2006). Analysis of Irish third-level college applications data. \emph{Journal of the Royal Statistical Society: Series A}, \bold{169}(2), pages 361--379, ISSN: 0964-1998, DOI: 10.1111/j.1467-985X.2006.00412.x.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link{mapPLMIX_multistart}} 
#'
#' @examples
#' 
#' data(d_carconf)
#' MAP <- mapPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_iter=400*3)
#' str(MAP)
#' MAP$P_map
#' MAP$W_map
#'
#' @export 
 
    cl <- match.call()
    
    if(class(pi_inv)[1]!="top_ordering"){
      if(class(pi_inv)[1]=="RankData"){
        pi_inv=as.top_ordering(data=pi_inv)
      }
      if(class(pi_inv)[1]=="rankings"){
        pi_inv=as.top_ordering(data=pi_inv)
      }
      if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
        pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
      }
    }
    pi_inv <- fill_single_entries(data=pi_inv)
    N <- nrow(pi_inv)
    n_rank <-  howmanyranked(pi_inv)

    rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
    
    ref_known <- TRUE
    ref_vary <- FALSE
    	
	if(is.null(init$omega)){
	#omega <- runif(G)
	#omega <- omega/sum(omega)	
        omega <- rdirichlet(1,rep(1,G))
        }else{
	omega <- init$omega
      if(sum(omega)!=1){
        warning("Initial mixture weights must add to one ==> input arguments has been normalized!") 
        omega <- omega/sum(omega)    
      }
    }


 		if(is.null(init$p)){

                    if(centered_start){

# print("CENTERED START !!")

            mle1comp <- matrix(prop.table(table(factor(pi_inv[,1],levels=1:K))),nrow=1)
            p <- random_start(mlesupp=mle1comp, givenweights=omega)
            p <- p/rowSums(p)

                    }else{

# print("COMPLETELY RANDOM (uniform support, rescaled) START")                    
# p <- matrix(runif(G*K),nrow=G,ncol=K)
# p <- p/rowSums(p)

            p <- rdirichlet(G,rep(1,K))
        }
	    }else{
	    p <- init$p
	      if(is.vector(p)){
          p <- t(p)
	      }
          if(!all(rowSums(p)==1)){
          warning("Initial support parameters for each mixture component must 
                       add to one ==> input arguments has been normalized!") 
          p <- p/rowSums(p)    
          }
        }

		
	init <- list(p=p,omega=omega)
	
    
	shape0 <- hyper$shape0
	rate0 <- hyper$rate0
	alpha0 <- hyper$alpha0
	
	u_bin <- umat(pi_inv=pi_inv)
			
	log_lik <- rep(NA,n_iter)
	
	if(!(all(shape0==1) & all(rate0==0) & all(shape0==1))){
#		print("Non-flat prior input")
		log_prior <- log_lik
	}
	
	objective <- log_lik
	conv <- 0
	l <- 1
		
	
	while(l<=n_iter){
		
	z_hat <- Estep(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv)
	
	omega <- UpWhet(z_hat=z_hat,alpha0=alpha0)

if(any(is.na(omega))){
print("==>  PROBLEM WITH *omega* update")
print(omega)
}

	p <- UpPhetpartial(p=p,ref_order=rho,pi_inv=pi_inv,z_hat=z_hat,shape0=shape0,
	                  rate0=rate0,n_rank=n_rank,u_bin=u_bin)

if(any(is.na(p))){
print("==>  PROBLEM WITH *p* update")
print(p)
}

    
    log_lik[l] <- loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv)
    if(is.na(log_lik[l])){
        print(p)
        print(omega)
        threshold <- -17
    while(is.na(log_lik[l]) & threshold<(-3)){
        p[p<=(10^threshold)] <- 10^threshold
        threshold  <-  threshold+1
        log_lik[l] <- loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv)
        print(paste0("Likelihood/parameter approximation for support parameter <=10^(-",threshold,")"))
}
    }
    
    if(!(all(shape0==1) & all(rate0==0) & all(shape0==1))){
          log_prior[l] <- log(ddirichlet(omega,alpha0))+sum(dgamma(p,shape=shape0,rate=rate0,log=TRUE))
          objective[l] <- log_lik[l]+log_prior[l]

        }else{
    	  objective[l] <- log_lik[l]
    }
    
    
    if(l>=2){

	    if((objective[l]-objective[l-1])/abs(objective[l-1])<eps |
               ((objective[l]-objective[l-1])==0 & objective[l-1]==0)){
		conv <- 1
		l <- n_iter+1
	       }
      }
    l <- l+1
    }
    
    P_map=p/rowSums(p)
    dimnames(P_map)=list(paste0("g_",1:G),paste0("p_",1:K))
    
    names(omega)=paste0("w_",1:G)
    
    log_lik <- log_lik[!(is.na(log_lik))]
    max_log_lik <- max(log_lik)

    objective <- objective[!(is.na(objective))]
    max_objective <- max(objective)
    if(all(shape0==1) & all(rate0==0) & all(shape0==1)){
      bic <- bicPLMIX(max_log_lik=max_log_lik,pi_inv=pi_inv,
                   G=G,ref_known=ref_known,
                   ref_vary=ref_vary)$bic
    }else{
      bic <- NULL
    }
    if(plot_objective){
     	plot(objective,ylab="Log-joint distribution",xlab="Iteration",
     	        main=paste("MAP estimation for PL mixture with",G,"components"),type="l")
     }
     
    out=list(W_map=omega,P_map=P_map,z_hat=z_hat,class_map=apply(z_hat,1,which.max),
               log_lik=log_lik,objective=objective,max_objective=max_objective,bic=bic,conv=conv,call=cl)

    class(out)="mpPLMIX"
    
    return(out)
    

}


mapPLMIX_multistart <- function(pi_inv,K,G,n_start=1,
							     init=rep(list(list(p=NULL,omega=NULL)),times=n_start),
                                 n_iter=200,
                                 hyper=list(shape0=matrix(1,nrow=G,ncol=K),rate0=rep(0,G),alpha0=rep(1,G)),
                                 eps=10^(-6),
                                 plot_objective=FALSE,
                                 init_index=1:n_start,
                                 parallel=FALSE,
                                 centered_start=FALSE){
#' MAP estimation for a Bayesian mixture of Plackett-Luce models with multiple starting values
#' 
#' Perform MAP estimation via EM algorithm with multiple starting values for a Bayesian mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' Under noninformative (flat) prior setting, the EM algorithm for MAP estimation corresponds to the EMM algorithm described by Gormley and Murphy (2006) to perform frequentist inference. In this case the MAP solution coincides with the MLE. The best model in terms of maximized posterior distribution is returned.
#' 
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param n_start Number of starting values.
#' @param init List of \code{n_start} lists of named objects with initialization values: \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters; \code{omega} is a numeric vector of \eqn{G} mixture weights. If starting values are not supplied (\code{NULL}), they are randomly generated with a uniform distribution. Default is \code{NULL}.
#' @param n_iter Maximum number of EM iterations.
#' @param hyper List of named objects with hyperparameter values for the conjugate prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is noninformative (flat) prior setting.
#' @param eps Tolerance value for the convergence criterion.
#' @param plot_objective Logical: whether the objective function (that is the kernel of the log-posterior distribution) should be plotted. Default is \code{FALSE}.
#' @param init_index Numeric vector indicating the positions of the starting values in the \code{init} list to be actually launched. Useful to launch the most promising starting values identified after a preliminary run. Default is run all the starting points in the \code{init} list.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#' @param centered_start Logical: whether a random start whose support parameters and weights should be centered around the observed relative frequency that each item has been ranked top. Default is \code{FALSE}. Ignored when \code{init} is not \code{NULL}.
#' 
#' @return A list of S3 class \code{mpPLMIX} with named elements:
#' 
#'  \item{\code{mod}}{ List of named objects describing the best model in terms of maximized posterior distribution. See output values of the single-run \code{\link{mapPLMIX}} function for a detailed explanation of the list elements.}
#'  \item{\code{max_objective}}{ Numeric vector of the maximized objective function values for each initialization.}
#'  \item{\code{convergence}}{ Binary vector with \code{length(init_index)} convergence indicators for each initialization: 1 = convergence has been achieved, 0 = otherwise.}
#'  \item{\code{call}}{ The matched call.}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Gormley, I. C. and Murphy, T. B. (2006). Analysis of Irish third-level college applications data. \emph{Journal of the Royal Statistical Society: Series A}, \bold{169}(2), pages 361--379, ISSN: 0964-1998, DOI: 10.1111/j.1467-985X.2006.00412.x.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link{mapPLMIX}} 
#'
#' @examples
#' 
#' data(d_carconf)
#' MAP_mult <- mapPLMIX_multistart(pi_inv=d_carconf, K=ncol(d_carconf), G=3, 
#'                                             n_start=2, n_iter=400*3)
#' str(MAP_mult)
#' MAP_mult$mod$P_map
#' MAP_mult$mod$W_map
#'
#' @export 
	
  cl <- match.call()
  
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
	for(i in 1:n_start){
#            print(paste0("Multiple starting point #",i))

            	if(is.null(init[[i]]$omega)){
	#omega <- runif(G)
	#omega <- omega/sum(omega)	
        omega <- rdirichlet(1,rep(1,G))
        }else{
	omega <- init[[i]]$omega
      if(sum(omega)!=1){
        warning("Initial mixture weights must add to one ==> input arguments has been normalized!") 
        omega <- omega/sum(omega)    
      }
	}

 		if(is.null(init[[i]]$p)){

                    if(centered_start){

#            print("CENTERED START !!")

            mle1comp <- matrix(prop.table(table(factor(pi_inv[,1],levels=1:K))),nrow=1)
            p <- random_start(mlesupp=mle1comp, givenweights=omega)
            p <- p/rowSums(p)

                    }else{

#            print("COMPLETELY RANDOM (uniform support, rescaled) START")                    
	    #p <- matrix(runif(G*K),nrow=G,ncol=K)
	    #p <- p/rowSums(p)

            p <- rdirichlet(G,rep(1,K))
        }
	    }else{
	    p <- init[[i]]$p
	      if(is.vector(p)){
          p <- t(p)
	      }
          if(!all(rowSums(p)==1)){
          warning("Initial support parameters for each mixture component must 
                       add to one ==> input arguments has been normalized!") 
          p <- p/rowSums(p)    
          }
        }

	
	
	init[[i]] <- list(p=p,omega=omega)
	}
	
	if(!parallel){
		
    mod <- vector(mode="list",length=length(init_index))
	max_objective <- rep(NA,length(init_index))     
	convergence <- rep(NA,length(init_index))
	record <- rep(NA,length(init_index))

	l <- 0

	for(i in init_index){
		l <- l+1
#		print(paste("INITIALIZATION",l))
	    mod[[l]] <- mapPLMIX(pi_inv=pi_inv,K=K,G=G,init=init[[i]],n_iter=n_iter,hyper=hyper,
                         eps=eps,centered_start=centered_start,plot_objective=plot_objective)
      max_objective[l] <- mod[[l]]$max_objective
      convergence[l] <- mod[[l]]$conv
      record[l] <- max(max_objective[1:l])
		print(paste("Starting value #",l," => best objective function value so far =",record[l]))
	}
    mod <- mod[[which.max(max_objective)]]
    class(mod) <- "list"
    mod <- mod[-length(mod)]
    out=list(mod=mod,max_objective=max_objective,convergence=convergence,call=cl)

		
	}else{
		
		
	mod <- foreach(i=init_index) %dopar%{   
    tempmod <- mapPLMIX(pi_inv=pi_inv,K=K,G=G,init=init[[i]],n_iter=n_iter,hyper=hyper,
                         eps=eps,centered_start=centered_start,plot_objective=plot_objective)
           }
    max_objective <- sapply(mod,"[[","max_objective")
    convergence <- sapply(mod,"[[","conv")
  
    outmod <- mod[[which.max(max_objective)]]
    class(outmod) <- "list"
    outmod <- outmod[-length(outmod)]
    out=list(mod=outmod,max_objective=max_objective,convergence=convergence,call=cl)
    
	}
  
  class(out)="mpPLMIX"
  
  return(out)
  
	
}

print.mpPLMIX <- function(x,...){
  #' Print of the MAP estimation algorithm for a Bayesian mixture of Plackett-Luce models
  #' 
  #' \code{print} method for class \code{mpPLMIX}. It shows some general information on the MAP estimation procedure for a Bayesian mixture of Plackett-Luce models.
  #'
  #'  
  #' @param x Object of class \code{mpPLMIX} returned by the \code{mapPLMIX} or \code{mapPLMIX_multistart} function.
  #' @param ... Further arguments passed to or from other methods (not used).
  #'
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @seealso \code{\link{mapPLMIX}} and \code{\link{mapPLMIX_multistart}} 
  #' 
  #' @examples
  #' 
  #' ## Print of the MAP procedure with a single starting point
  #' data(d_carconf)
  #' MAP <- mapPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3)
  #' print(MAP)
  #' 
  #' ## Print of the MAP procedure with 5 starting points
  #' MAP_multi <- mapPLMIX_multistart(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_start=5)
  #' print(MAP_multi)
  #' @export print.mpPLMIX  
  #' @export 
  
  mpPLMIX_out=x
  
  if(class(mpPLMIX_out)!="mpPLMIX"){
    stop("The function requires an object of S3 class 'mpPLMIX' as its first argument.")
  }
  
  cat("\nCall:\n", paste(deparse(mpPLMIX_out$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if(is.null(mpPLMIX_out$convergence)){
    G=length(mpPLMIX_out$W_map)
    K=ncol(mpPLMIX_out$P_map)
    L=length(mpPLMIX_out$log_lik)
    cat("MAP estimation procedure for a Bayesian mixture of Plackett-Luce models:\n")
    if(!is.null(mpPLMIX_out$bic)){
      cat("Prior distribution used: flat (default) ====> MAP = MLE\n")
    }else{
      cat("Prior distribution used: subjective (see 'Call')\n")
    }
    cat("\n")
    cat("No. of items:",K,"\n")
    cat("No. of mixture components:",G,"\n")
    cat("No. of iterations:",L,"\n")
    cat("\n")
    cat("Max. log-likelihood:",max(mpPLMIX_out$log_lik,na.rm=TRUE),"\n")
    cat("Max. objective function:",mpPLMIX_out$max_objective,"\n")
    if(!is.null(mpPLMIX_out$bic)){
      cat("BIC:",mpPLMIX_out$bic,"\n")
    }
    cat("\n")
    cat("Algorithm convergence check:",ifelse(mpPLMIX_out$conv,"Ok.","Failed."),"\n")
  }else{
    G=length(mpPLMIX_out$mod$W_map)
    K=ncol(mpPLMIX_out$mod$P_map)
    L=length(mpPLMIX_out$mod$log_lik)
    cat("MAP estimation procedure for a Bayesian mixture of Plackett-Luce models with",length(mpPLMIX_out$convergence),"starting values ====> best solution reported:\n")
    if(!is.null(mpPLMIX_out$mod$bic)){
      cat("Prior distribution used: flat (default) ====> MAP = MLE\n")
    }else{
      cat("Prior distribution used: subjective (see 'Call')\n")
    }
    cat("\n")
    cat("No. of items:",K,"\n")
    cat("No. of mixture components:",G,"\n")
    cat("No. of iterations:",L,"\n")
    cat("\n")
    cat("Max. log-likelihood:",max(mpPLMIX_out$mod$log_lik,na.rm=TRUE),"\n")
    cat("Max. objective function:",mpPLMIX_out$mod$max_objective,"\n")
    if(!is.null(mpPLMIX_out$mod$bic)){
      cat("BIC:",mpPLMIX_out$mod$bic,"\n")
    }
    cat("\n")
    cat("Algorithm convergence check:",ifelse(mpPLMIX_out$mod$conv,"Ok.","Failed."),"\n")
    
  }
  cat("\n")
  cat("Use functions summary() and plot() to summarize and visualize the object of class 'mpPLMIX'.")
  
}

print.summary.mpPLMIX <- function(x,...){
  #/' Print of the summary of MAP estimation for a Bayesian mixture of Plackett-Luce models
  #/' 
  #/' \code{print} method for class \code{summary.mpPLMIX}. It provides summaries for the MAP estimation of a Bayesian mixture of Plackett-Luce models.
  #/'
  #/'  
  #/' @param x Object of class \code{summary.mpPLMIX} returned by the \code{summary.mpPLMIX} function.
  #/' @param ... Further arguments passed to or from other methods (not used).
  #/'
  #/'
  #/' @references 
  #/' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
  #/'
  #/' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
  #/' 
  #/' @author Cristina Mollica and Luca Tardella
  
  summary.mpPLMIX_out=x
  
  if(class(summary.mpPLMIX_out)!="summary.mpPLMIX"){
    stop("The function requires an object of S3 class 'summary.mpPLMIX' as its first argument.")
  }
  
  G=length(summary.mpPLMIX_out$MAP_w)
  cat("\nCall:\n", paste(deparse(summary.mpPLMIX_out$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if(G>1){
    cat("MAP estimates of the mixture weghts:\n")
    print(summary.mpPLMIX_out$MAP_w)
    cat("\n")
  }  
  cat("MAP estimates of the support parameters:\n")
  print(summary.mpPLMIX_out$MAP_p)
  cat("\n")  
  cat("Estimated component-specific modal orderings:\n")
  print(summary.mpPLMIX_out$Modal_orderings)
  cat("\n")
  if(G>1){
    cat("Relative frequency distribution of group memberships:\n")
    print(summary.mpPLMIX_out$group_distr)
    cat("\n")
  }
  if(!is.null(summary.mpPLMIX_out$perc_conv_rate)){
    cat("Convergence percentage over multiple initialization:\n")
    print(c(Conv=summary.mpPLMIX_out$perc_conv_rate))
  }
  
} 


summary.mpPLMIX <- function(object,digits=2,...){
  #' Summary of the MAP estimation for a Bayesian mixture of Plackett-Luce models
  #' 
  #' \code{summary} method for class \code{mpPLMIX}. It provides summaries for the MAP estimation of a Bayesian mixture of Plackett-Luce models.
  #'
  #' @param object Object of class \code{mpPLMIX} returned by the \code{mapPLMIX} or \code{mapPLMIX_multistart} function.
  #' @param digits Number of decimal places for rounding the summaries.
  #' @param ... Further arguments passed to or from other methods (not used). 
  #'
  #' @return A list of summaries for the \code{mpPLMIX} class object:
  #' 
  #'  \item{\code{MAP_w}}{ Numeric vector with the MAP estimates of the \eqn{G} mixture weights. Returned only when when \eqn{G>1}.}
  #'  \item{\code{MAP_p}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters.}
  #'  \item{\code{MAP_modal_orderings}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the estimated modal orderings of each mixture component.}
  #'  \item{\code{group_distr}}{ Numeric vector with the relative frequency distribution of the mixture component memberships based on MAP allocation. Returned only when when \eqn{G>1}.}
  #'  \item{\code{perc_conv_rate}}{ Numeric scalar with the percentage of MAP algorithm convergence over the multiple starting points. Returned only when \code{summary.mpPLMIX} is applied to the output of the \code{mapPLMIX_multistart} function.}
  #'
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @examples
  #' 
  #' ## Summary of the MAP procedure with a single starting point
  #' data(d_carconf)
  #' MAP <- mapPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3)
  #' summary(MAP)
  #' 
  #' ## Summary of the MAP procedure with 5 starting points
  #' MAP_multi <- mapPLMIX_multistart(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_start=5)
  #' summary(MAP_multi)
  #' @export summary.mpPLMIX  
  #' @export 
  
  mpPLMIX_out=object
  
  if(class(mpPLMIX_out)!="mpPLMIX"){
    stop("The function requires an object of S3 class 'mpPLMIX' as its first argument.")
  }
  
  cl=mpPLMIX_out$call
  if(is.null(mpPLMIX_out$convergence)){
    G=length(mpPLMIX_out$W_map)
    K=ncol(mpPLMIX_out$P_map)
    out=list(MAP_w=mpPLMIX_out$W_map,MAP_p=mpPLMIX_out$P_map,
             Modal_orderings=t(apply(matrix(mpPLMIX_out$P_map,nrow=G,ncol=K),1,order,decreasing=TRUE)),
             group_distr=prop.table(table(factor(mpPLMIX_out$class_map,levels=1:G))),
             call=cl)
    out[c(1:2,4)]=lapply(out[c(1:2,4)],round,digits=digits)
    
    dimnames(out$Modal_orderings) <- list(paste0("g_",1:G),paste0("Rank_",1:K))
    
  }else{
    G=length(mpPLMIX_out$mod$W_map)
    K=ncol(mpPLMIX_out$mod$P_map)
    out=list(MAP_w=mpPLMIX_out$mod$W_map,MAP_p=mpPLMIX_out$mod$P_map,
             Modal_orderings=t(apply(matrix(mpPLMIX_out$mod$P_map,nrow=G,ncol=K),1,order,decreasing=TRUE)),
             group_distr=prop.table(table(factor(mpPLMIX_out$mod$class_map,levels=1:G))),
             perc_conv_rate=100*mean(mpPLMIX_out$convergence),
             call=cl)
    out[c(1:2,4:5)]=lapply(out[c(1:2,4:5)],round,digits=digits)
    
  }
  
  dimnames(out$Modal_orderings) <- list(paste0("g_",1:G),paste0("Rank_",1:K))
  
  class(out)="summary.mpPLMIX"
  
  out
  
}

plot.mpPLMIX <- function(x,max_scale_radar=NULL,...){
  #' Plot the MAP estimates for a Bayesian mixture of Plackett-Luce models
  #' 
  #' \code{plot} method for class \code{mpPLMIX}.
  #'
  #' By recalling the \code{chartJSRadar} function from the \code{radarchart} package and the routines of the \code{ggplot2} package, \code{plot.mpPLMIX} produces a radar plot of the support parameters and, when \eqn{G>1}, a donut plot of the mixture weights and a heatmap of the component membership probabilities based on the MAP estimates.  The radar chart is returned in the Viewer Pane.
  #'  
  #' @param x Object of class \code{mpPLMIX} returned by the \code{mpPLMIX} function.
  #' @param max_scale_radar Numeric scalar indicating the maximum value on each axis of the radar plot for the support parameter point estimates. Default is \code{NULL} meaning that the maximum of the estimated support parameters is used.
  #' @param ... Further arguments passed to or from other methods (not used). 
  #'
  #'
  #' @references 
  #' Ashton, D. and Porter, S. (2016). radarchart: Radar Chart from 'Chart.js'. R package version 0.3.1. \url{https://CRAN.R-project.org/package=radarchart}
  #' 
  #' Wickham, H. (2009). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.
  #' 
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @seealso \code{\link[radarchart]{chartJSRadar}} and \code{\link[ggplot2]{ggplot}}
  #' 
  #' @examples
  #' 
  #' # Not run:
  #' data(d_carconf)
  #' MAP <- mapPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3)
  #' plot(MAP)
  #' 
  #' # Not run:
  #' MAP_multi <- mapPLMIX_multistart(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_start=5)
  #' plot(MAP_multi)
  #' @export plot.mpPLMIX  
  #' @export 
  
  mpPLMIX_out=x
  
  if(class(mpPLMIX_out)!="mpPLMIX"){
    stop("The function requires an object of S3 class 'mpPLMIX' as its first argument.")
  }
  
  if(is.null(mpPLMIX_out$convergence)){
    G=length(mpPLMIX_out$W_map)
    K=ncol(mpPLMIX_out$P_map)
    N=nrow(mpPLMIX_out$z_hat)
  }else{
    G=length(mpPLMIX_out$mod$W_map)
    K=ncol(mpPLMIX_out$mod$P_map)
    N=nrow(mpPLMIX_out$mod$z_hat)
  }
  
  labs <- paste("Item",1:K)
  if(is.null(mpPLMIX_out$convergence)){
    scores=as.list(as.data.frame(t(mpPLMIX_out$P_map)))
  }else{
    scores=as.list(as.data.frame(t(mpPLMIX_out$mod$P_map)))
  }
  names(scores)=paste("Group",1:G)
  main_radar="MAP estimates of the support parameters"
  oo=chartJSRadar(scores = scores, labs = labs, main=main_radar,maxScale = ifelse(is.null(max_scale_radar),max(unlist(scores)),max_scale_radar))
  print(oo)

  if(G>1){
    if(is.null(mpPLMIX_out$convergence)){
      df_w <- data.frame(Composition = paste0(paste("Group",1:G),":"),value=mpPLMIX_out$W_map,label=paste(paste0(paste("Group",1:G),":"),paste0(round(mpPLMIX_out$W_map*100), "%")))
      df_z=as.data.frame(mpPLMIX_out$z_hat)
    }else{
      df_w <- data.frame(Composition = paste0(paste("Group",1:G),":"),value=mpPLMIX_out$mod$W_map,label=paste(paste0(paste("Group",1:G),":"),paste0(round(mpPLMIX_out$mod$W_map*100), "%")))
      df_z=as.data.frame(mpPLMIX_out$mod$z_hat)
    }
    pp=ggplot(df_w, aes_string(x=2, y = "value", fill = "label")) +
       geom_bar(stat = "identity", color = "white") +
       coord_polar(theta = "y", start = 0)+
       labs(x = NULL, y = NULL, fill = NULL, title = "Sample composition by group membership")+
       scale_fill_brewer(palette="Blues")+
       theme_void()+
       xlim(0.5, 2.5)
      
    names(df_z)=paste("Group",1:G)
    df_z=data.frame(Unit=paste("Unit",1:N),df_z)
    df_z_m=melt(df_z,id.vars="Unit")
    # zz=ggplot(df_z_m, aes_string("variable", "Unit")) +
    #    geom_tile(aes_string(fill = "value"), colour = "white") +
    #    labs(x = "", y = "Sample units", fill = NULL, title = "Component membership probabilities")+
    #    theme(axis.text.y = element_blank(),
    #    axis.ticks = element_blank())+
    #    scale_fill_gradient(low = "white", high = "steelblue")
      
    zz=ggplot(df_z_m, aes_string("Unit", "variable")) +
      geom_tile(aes_string(fill = "value"), colour = "white") +
      labs(x = "Sample units", y = "", fill = NULL, title = "Component membership probabilities")+
      theme(axis.text.x = element_blank(),
            axis.ticks = element_blank())+
      scale_fill_gradient(low = "white", high = "steelblue")

    grid.arrange(pp,zz,nrow=2)
    
  }
  
}



##########################################################       
############# GIBBS SAMPLING #############################

gibbsPLMIX <- function(pi_inv,K,G,
						   init=list(z=NULL,p=NULL),
						   n_iter=1000,
						   n_burn=500,
               hyper=list(shape0=matrix(1,nrow=G,ncol=K),rate0=rep(0.001,G),alpha0=rep(1,G)),
               centered_start=FALSE){
#' Gibbs sampling for a Bayesian mixture of Plackett-Luce models
#' 
#' Perform Gibbs sampling simulation for a Bayesian mixture of Plackett-Luce models fitted to partial orderings.
#' 
#' The size \eqn{L} of the final MCMC sample is equal to \code{n_iter}-\code{n_burn}.
#' 
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @param K Number of possible items. 
#' @param G Number of mixture components.
#' @param init List of named objects with initialization values: \code{z} is a numeric \eqn{N}\eqn{\times}{x}\eqn{G} matrix of binary mixture component memberships; \code{p} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of component-specific support parameters. If starting values are not supplied (\code{NULL}), they are randomly generated with a uniform distribution. Default is \code{NULL}.
#' @param n_iter Total number of MCMC iterations.
#' @param n_burn Number of initial burn-in drawings removed from the returned MCMC sample.
#' @param hyper List of named objects with hyperparameter values for the conjugate prior specification: \code{shape0} is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of shape hyperparameters; \code{rate0} is a numeric vector of \eqn{G} rate hyperparameters; \code{alpha0} is a numeric vector of \eqn{G} Dirichlet hyperparameters. Default is vague prior setting.
#' @param centered_start Logical: whether a random start whose support parameters and weights should be centered around the observed relative frequency that each item has been ranked top. Default is \code{FALSE}. Ignored when \code{init} is not \code{NULL}.
#'
#' @return A list of S3 class \code{gsPLMIX} with named elements:
#' 
#'  \item{\code{W}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with MCMC samples of the mixture weights.}
#'  \item{\code{P}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with MCMC samples of the component-specific support parameters.}
#'  \item{\code{log_lik}}{ Numeric vector of \eqn{L} posterior log-likelihood values.}
#'  \item{\code{deviance}}{ Numeric vector of \eqn{L} posterior deviance values (\eqn{-2 * }\code{log_lik}).}
#'  \item{\code{objective}}{ Numeric vector of \eqn{L} objective function values (that is the kernel of the log-posterior distribution).}
#'  \item{\code{call}}{ The matched call.}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#' 
#' @author Cristina Mollica and Luca Tardella
#' 
#' @examples
#' 
#' data(d_carconf)
#' GIBBS <- gibbsPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_iter=30, n_burn=10)
#' str(GIBBS)
#' GIBBS$P
#' GIBBS$W
#' 
#' @export 

    cl=match.call()
    
    if(class(pi_inv)[1]!="top_ordering"){
      if(class(pi_inv)[1]=="RankData"){
        pi_inv=as.top_ordering(data=pi_inv)
      }
      if(class(pi_inv)[1]=="rankings"){
        pi_inv=as.top_ordering(data=pi_inv)
      }
      if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
        pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
      }
    }
    pi_inv <- fill_single_entries(data=pi_inv)
    N <- nrow(pi_inv)
    n_rank <-  howmanyranked(pi_inv)
    rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
	
	if(is.null(init$z)){
	z <- binary_group_ind(class=sample(1:G,size=N,replace=TRUE),G=G)
	}else{
	z <- init$z
	}


	omega <- colMeans(z)
    
	if(is.null(init$p)){
         if(centered_start){

            print("CENTERED START !!")

            # omega <- rdirichlet(1,rep(1,G))
            mle1comp <- matrix(prop.table(table(factor(pi_inv[,1],levels=1:K))),nrow=1)
            p <- random_start(mlesupp=mle1comp, givenweights=omega)
            # p <- p/rowSums(p)

          }else{

            print("COMPLETELY RANDOM (uniform support, rescaled) START")                    
		
			p <- matrix(rgamma(n=G*K,shape=1,rate=1),nrow=G,ncol=K)
    
          }
	
	}else{
	
	    p <- init$p
	}

	shape0 <- hyper$shape0
	rate0 <- hyper$rate0
	alpha0 <- hyper$alpha0
	
	u_bin <- umat(pi_inv=pi_inv)
    
	
    log_lik <- c(loglikPLMIX(p=p,ref_order=rho,weights=omega,pi_inv=pi_inv),
              rep(NA,n_iter))
    
    log_prior <- c(log(ddirichlet(omega,alpha0))+sum(dgamma(p,shape=shape0,rate=rate0,log=TRUE)),
                 rep(NA,n_iter))
    
    objective <- log_lik+log_prior

    Pi <- array(NA,dim=c(G,K,n_iter+1))
    Pi[,,1] <- p
    Zeta <- z
    Omega <- matrix(NA,nrow=n_iter+1,ncol=G)
    Omega[1,] <- omega
    
    
for(l in 1:n_iter){

    if(l%%500==0){
    print(paste("GIBBS iteration",l))
}
    
Omega[l+1,] <- rdirichlet(n=1,alpha=alpha0+colSums(Zeta))

temprate <- CompRateYpartial(p=adrop(Pi[,,l,drop=FALSE],3),pi_inv=pi_inv,ref_order=rho,z=Zeta,n_rank=n_rank)
Ypsilon <- SimYpsilon(rate=temprate,n_rank=n_rank)
    
Pi[,,l+1] <- matrix(rgamma(n=G*K,shape=shape0+gammamat(u_bin=u_bin,z_hat=Zeta),
          rate <- CompRateP(pi_inv=pi_inv, Y=Ypsilon, z=Zeta, u_bin=u_bin, n_rank=n_rank, rate0=rate0)),nrow=G,ncol=K)


Zeta <- binary_group_ind(apply(CompProbZpartial(p=adrop(Pi[,,l+1,drop=FALSE],3),pi_inv=pi_inv,Y=Ypsilon, u_bin=u_bin,n_rank,omega=Omega[l+1,]),1,FUN=sample,x=1:G,replace=TRUE,size=1),G=G)

log_lik[l+1] <- loglikPLMIX(p=adrop(Pi[,,l+1,drop=FALSE],3),ref_order=rho,weights=Omega[l+1,],
                                     pi_inv=pi_inv)

log_prior[l+1] <- log(ddirichlet(Omega[l+1,],alpha0))+sum(dgamma(adrop(Pi[,,l+1,drop=FALSE],3),shape=shape0,rate=rate0,log=TRUE))

objective[l+1] <- log_lik[l+1]+log_prior[l+1]

    }

log_lik <- log_lik[-c(1:(n_burn+1))]

objective <- objective[-c(1:(n_burn+1))]

Omega <- Omega[-c(1:(n_burn+1)),,drop=FALSE]
colnames(Omega) <- paste0("w_",1:G)


Pi <- array(apply(Pi,3,FUN=function(x)x/rowSums(x)),c(G,K,n_iter+1))	

Pi=t(apply(Pi,3,c))[-c(1:(n_burn+1)),]
colnames(Pi) <- paste0("p_",rep(1:G,K),rep(1:K,each=G))

out=list(W=Omega,P=Pi,log_lik=log_lik,deviance=-2*log_lik,objective=objective,call=cl)

class(out)="gsPLMIX"

return(out)

}

gsPLMIX_to_mcmc <- function(gsPLMIX_out){
  #' MCMC class objects from the Gibbs sampling simulations of a Bayesian mixture of Plackett-Luce models
  #' 
  #' Coerce the Gibbs sampling simulations for a Bayesian mixture of Plackett-Luce models into an \code{mcmc} class object.
  #' 
  #' \code{gsPLMIX_to_mcmc} attemps to coerce its argument by recalling the \code{as.mcmc} function of the \code{coda} package.
  #' 
  #' @param gsPLMIX_out Object of class \code{gsPLMIX} returned by the \code{gibbsPLMIX} function.
  #'
  #' @return An \code{mcmc} class object.
  #' 
  #' @references 
  #' Plummer, M., Best, N., Cowles, K. and Vines, K. (2006). CODA: Convergence Diagnosis and Output Analysis for MCMC, \emph{R News}, \bold{6}, pages 7--11, ISSN: 1609-3631.
  #' 
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @seealso \code{\link[coda]{as.mcmc}} 
  #' 
  #' @examples
  #' 
  #' data(d_carconf)
  #' GIBBS <- gibbsPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_iter=30, n_burn=10)
  #' 
  #' ## Coerce the posterior samples into an mcmc class object
  #' gsPLMIX_to_mcmc(GIBBS)
  #' 
  #' @export 
  
  if(class(gsPLMIX_out)!="gsPLMIX"){
    stop("The function requires an object of S3 class 'gsPLMIX' as its first argument.")
  }

  G=ncol(gsPLMIX_out$W)
  K=ncol(gsPLMIX_out$P)/G
  class(gsPLMIX_out)="list"
  gsPLMIX_out=gsPLMIX_out[-length(gsPLMIX_out)] # to remove call element
  out=as.data.frame(gsPLMIX_out)
  colnames(out)=c(paste0("w_",1:G),paste0("p_",rep(1:G,K),rep(1:K,each=G)),"log_lik","deviance","objective")
  out=as.mcmc(out)
  
  return(out)
  
}


summary.gsPLMIX <- function(object,quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975),hpd_prob=0.95,digits=2,...){
  #' Summary of the Gibbs sampling procedure for a Bayesian mixture of Plackett-Luce models
  #' 
  #' \code{summary} method for class \code{gsPLMIX}. It provides summary statistics and credible intervals for the Gibbs sampling simulation of a Bayesian mixture of Plackett-Luce models.
  #'
  #' Posterior summaries include means, standard deviations, naive standard errors of the means (ignoring autocorrelation of the chain) and time-series standard errors based on an estimate of the spectral density at 0. They correspond to the \code{statistics} element of the output returned by the \code{summary.mcmc} function of the \code{coda} package. Highest posterior density (HPD) intervals are obtained by recalling the \code{HPDinterval} function of the \code{coda} package.
  #'  
  #' @param object Object of class \code{gsPLMIX} returned by the \code{gibbsPLMIX} function.
  #' @param quantiles Numeric vector of quantile probabilities.
  #' @param hpd_prob Numeric scalar in the grid of values spanning the interval (0,1) by 0.05, giving the posterior probability content of the HPD intervals. Supplied values outside the grid are rounded.
  #' @param digits Number of decimal places for rounding the posterior summaries.
  #' @param ... Further arguments passed to or from other methods (not used). 
  #'
  #' @return A list of summary statistics for the \code{gsPLMIX} class object:
  #' 
  #'  \item{\code{statistics}}{ Numeric matrix with posterior summaries in each row (see 'Details').}
  #'  \item{\code{quantiles}}{ Numeric matrix with posterior quantiles at the given \code{quantiles} probabilities in each row.}
  #'  \item{\code{HPDintervals}}{ Numeric matrix with 100\eqn{*}\code{hpd_prob}\% HPD intervals in each row.}
  #'  \item{\code{Modal_orderings}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the estimated posterior modal orderings of each mixture component.}
  #'  \item{\code{call}}{ The matched call.}
  #'
  #' @references 
  #' Plummer, M., Best, N., Cowles, K. and Vines, K. (2006). CODA: Convergence Diagnosis and Output Analysis for MCMC, \emph{R News}, \bold{6}, pages 7--11, ISSN: 1609-3631.
  #'
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @seealso \code{\link[coda]{summary.mcmc}} and \code{\link[coda]{HPDinterval}} 
  #' 
  #' @examples
  #' 
  #' data(d_carconf)
  #' GIBBS <- gibbsPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_iter=30, n_burn=10)
  #' 
  #' ## Summary of the Gibbs sampling procedure
  #' summary(GIBBS)
  #' @export summary.gsPLMIX  
  #' @export 
  
  gsPLMIX_out=object
  G=ncol(gsPLMIX_out$W)
  K=ncol(gsPLMIX_out$P)/G
  cl=gsPLMIX_out$call
  mcmc_obj=gsPLMIX_to_mcmc(gsPLMIX_out)
  p_idx=grep(pattern="p_",x=colnames(mcmc_obj))
  
  temp=getFromNamespace("summary.mcmc",ns="coda")(object=mcmc_obj,quantiles=quantiles)[1:2]
  hpd_int=HPDinterval(mcmc_obj,prob=hpd_prob)
#  attr(hpd_int,"Probability")=NULL
  out=list(statistics=temp[[1]],quantiles=as.matrix(temp[[2]]),HPD_intervals=hpd_int,
           Modal_orderings=t(apply(matrix(temp$statistics[p_idx,"Mean"],nrow=G,ncol=K),1,order,decreasing=TRUE)),
           call=cl)
  out[1:3]=lapply(out[1:3],round,digits=digits)
  
  names(out)[3]=paste0(100*hpd_prob,"%_HPD_intervals")
  if(length(quantiles)==1){
    colnames(out$quantiles)=paste0(quantiles*100,"%")
  }
  dimnames(out$Modal_orderings) <- list(paste0("g_",1:G),paste0("Rank_",1:K))
  
  class(out)="summary.gsPLMIX"

  out
  
}

plot.gsPLMIX <- function(x,file="ggmcmc-output.pdf",family=NA,plot=NULL,param_page=5,width=7,height=10,dev_type_html="png",post_est="mean",max_scale_radar=NULL,...){
  #' Plot the Gibbs sampling simulations for a Bayesian mixture of Plackett-Luce models
  #' 
  #' \code{plot} method for class \code{gsPLMIX}. It builds a suite of plots, visual convergence diagnostics and credible intervals for the MCMC samples of a Bayesian mixture of Plackett-Luce models. Graphics can be plotted directly into the current working device or stored into an external file placed into the current working directory.
  #'
  #' Plots of the MCMC samples include histograms, densities, traceplots, running means plots, overlapped densities comparing the complete and partial samples, autocorrelation functions, crosscorrelation plots and caterpillar plots of the 90 and 95\% equal-tails credible intervals. Note that the latter are created for the support parameters (when either \code{family=NA} or \code{family="p"}), for the mixture weights in the case \eqn{G>1} (when either \code{family=NA} or \code{family="w"}), for the log-likelihood values (when \code{family="log_lik"}), for the deviance values (when \code{family="deviance"}). Convergence tools include the potential scale reduction factor and the Geweke z-score. These functionalities are implemented with a call to the \code{ggs} and \code{ggmcmc} functions of the \code{ggmcmc} package (see 'Examples' for the specification of the \code{plot} argument) and for the objective function values (when \code{family="objective"}). 
  #' 
  #' By recalling the \code{chartJSRadar} function from the \code{radarchart} package and the routines of the \code{ggplot2} package, \code{plot.gsPLMIX} additionally produces a radar plot of the support parameters and, when \eqn{G>1}, a donut plot of the mixture weights based on the posterior point estimates. The radar chart is returned in the Viewer Pane.
  #'  
  #' @param x Object of class \code{gsPLMIX} returned by the \code{gibbsPLMIX} function.
  #' @param file Character vector with the name of the file to be created in the current working directory. Defaults is "ggmcmc-output.pdf". When NULL, plots are directly returned into the current working device (not recommended). This option allows also the user to work with an opened pdf (or other) device. When the file has an html file extension, the output is an Rmarkdown report with the figures embedded in the html file.
  #' @param family Character string indicating the name of the family of parameters to be plotted. A family of parameters is considered to be any group of parameters with the same name but different numerical values (for example \code{w[1]}, \code{w[2]}, etc). Default is \code{NA} meaning that all the parameters in the chain are plotted. Alternatively, one can choose \code{"w"}, \code{"p"}, \code{"log_lik"}, \code{"deviance"} or \code{"objective"}.
  #' @param plot Character vector containing the names of the desired plots. Default is \code{NULL} meaning that all the plots and convergence diagnostics are built (see 'Details'). 
  #' @param param_page Number of parameters to be plotted in each page. Defaults is 5.
  #' @param width Numeric scalar indicating the width of the pdf display in inches. Defaults is 7.
  #' @param height Numeric scalar indicating the height of the pdf display in inches. Defaults is 10.
  #' @param dev_type_html Character vector indicating the type of graphical device for the html output. Default is \code{"png"}. Alternatively, one can choose \code{"svg"}.
  #' @param post_est Character string indicating the point estimates of the Plackett-Luce mixture parameters to be computed from the \code{gsPLMIX} class object and then plotted in the current working device. Default is \code{"mean"}. Alternatively, one can choose \code{"median"}.
  #' @param max_scale_radar Numeric scalar indicating the maximum value on each axis of the radar plot for the support parameter point estimates. Default is \code{NULL} meaning that the maximum of the estimated support parameters is used.
  #' @param ... Further arguments passed to or from other methods (not used). 
  #'
  #'
  #' @references 
  #' Ashton, D. and Porter, S. (2016). radarchart: Radar Chart from 'Chart.js'. R package version 0.3.1. \url{https://CRAN.R-project.org/package=radarchart}
  #' 
  #' Wickham, H. (2009). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.
  #' 
  #' Fernandez-i-Marin, X. (2006). ggmcmc: Analysis of MCMC Samples and Bayesian Inference, \emph{Journal of Statistical Software}, \bold{70}(9), pages 1--20, DOI: 10.18637/jss.v070.i09.
  #' 
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @seealso \code{\link[ggmcmc]{ggs}}, \code{\link[ggmcmc]{ggmcmc}}, \code{\link[radarchart]{chartJSRadar}} and \code{\link[ggplot2]{ggplot}}
  #' 
  #' @examples
  #' 
  #' # Not run:
  #' data(d_carconf)
  #' GIBBS <- gibbsPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=5, n_iter=30, n_burn=10)
  #' 
  #' # Not run:
  #' # Plot posterior samples supplied as an gsPLMIX class object
  #' # plot(GIBBS)
  #' 
  #' # Selected plots of the posterior samples of the support parameters
  #' # plot(GIBBS, family="p", plot=c("compare_partial","Rhat","caterpillar"), param_page=6)
  #' 
  #' # Selected plots of the posterior samples of the mixture weights
  #' # plot(GIBBS, family="w", plot=c("histogram","running","crosscorrelation","caterpillar"))
  #' 
  #' # Selected plots of the posterior log-likelihood values
  #' # plot(GIBBS, family="log_lik", plot=c("autocorrelation","geweke"), param_page=1)
  #' 
  #' # Selected plots of the posterior deviance values
  #' # plot(GIBBS, family="deviance", plot=c("traceplot","density"), param_page=1)

  #' @export plot.gsPLMIX  
  #' @export 
  
  gsPLMIX_out=x
  mcmc_obj=gsPLMIX_to_mcmc(gsPLMIX_out=gsPLMIX_out)
  G=ncol(gsPLMIX_out$W)
  K=ncol(gsPLMIX_out$P)/G
  n_par=G+G*K
  colnames(mcmc_obj)[1:n_par]=gsub("_","[",colnames(mcmc_obj)[1:n_par])
  colnames(mcmc_obj)[1:n_par]=paste0(colnames(mcmc_obj)[1:n_par],"]")
  if(G==1){
    mcmc_obj=mcmc_obj[,-1]
  }
  tbl_obj = ggs(S=mcmc_obj)
  
  if(!is.na(family) & family=="w" & G==1){
    message(paste("No. of mixture components:",G, "====> w[1] = 1"))
  }else{
    simplify_traceplot=NULL
    mcmc_plot=ggmcmc(tbl_obj,file=file,family=family,plot=plot,param_page=param_page,width=width,height=height,
                     simplify_traceplot=simplify_traceplot,dev_type_html=dev_type_html)
    print(mcmc_plot)
  }
  
#  if(!is.na(family) & family=="p"){
    labs <- paste("Item",1:K)
    temp_radar=summary(object=gsPLMIX_out,quantiles=0.5)
    if(post_est=="mean"){
      scores=as.list(as.data.frame(t(matrix(temp_radar$statistics[grep("p",rownames(temp_radar$quantiles)),"Mean"],nrow=G,ncol=K))))
      main_radar="Posterior means of the support parameters"
    }else{
      scores=as.list(as.data.frame(t(matrix(temp_radar$quantiles[grep("p",rownames(temp_radar$quantiles)),],nrow=G,ncol=K))))
      main_radar="Posterior medians of the support parameters"
    }
    names(scores)=paste("Group",1:G)
    oo=chartJSRadar(scores = scores, labs = labs, main=main_radar,maxScale = ifelse(is.null(max_scale_radar),max(unlist(scores)),max_scale_radar))
    print(oo)
#  }  
  
#  if(!is.na(family) & family=="w" & G>1){
  if(G>1){  
    temp_radar=summary(object=gsPLMIX_out,quantiles=0.5)
    if(post_est=="mean"){
      temp_value=temp_radar$statistics[grep("w",rownames(temp_radar$quantiles)),"Mean"]
    }else{
      temp_value=temp_radar$quantiles[grep("w",rownames(temp_radar$quantiles)),]
    }
    df_w <- data.frame(Composition = paste0(paste("Group",1:G),":"),value=temp_value,label=paste(paste0(paste("Group",1:G),":"),paste0(round(temp_value*100), "%")))
    pp=ggplot(df_w, aes_string(x = 2, y = "value", fill = "label")) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar(theta = "y", start = 0)+
      labs(x = NULL, y = NULL, fill = NULL, title = "Sample composition by group membership")+
      #      geom_text(aes(label = paste0(round(value*100), "%")), position = position_stack(vjust = 0.5))+
      scale_fill_brewer(palette="Blues")+
      theme_void()+
      xlim(0.5, 2.5)
    print(pp)
  }
  
}

print.gsPLMIX <- function(x,...){
  #' Print of the Gibbs sampling simulation of a Bayesian mixture of Plackett-Luce models
  #' 
  #' \code{print} method for class \code{gsPLMIX}. It shows some general information on the Gibbs sampling simulation for a Bayesian mixture of Plackett-Luce models.
  #'
  #'  
  #' @param x Object of class \code{gsPLMIX} returned by the \code{gibbsPLMIX} function.
  #' @param ... Further arguments passed to or from other methods (not used).
  #'
  #' @author Cristina Mollica and Luca Tardella
  #' 
  #' @seealso \code{\link{gibbsPLMIX}} 
  #' 
  #' @examples
  #' 
  #' ## Print of the Gibbs sampling procedure
  #' data(d_carconf)
  #' GIBBS <- gibbsPLMIX(pi_inv=d_carconf, K=ncol(d_carconf), G=3, n_iter=30, n_burn=10)
  #' print(GIBBS)
  #' @export print.gsPLMIX  
  #' @export 
  
  gsPLMIX_out=x
  
  if(class(gsPLMIX_out)!="gsPLMIX"){
    stop("The function requires an object of S3 class 'gsPLMIX' as its first argument.")
  }

  G=ncol(gsPLMIX_out$W)
  K=ncol(gsPLMIX_out$P)/G
  L=nrow(gsPLMIX_out$W)
  cat("\nCall:\n", paste(deparse(gsPLMIX_out$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Gibbs sampling procedure for a Bayesian mixture of Plackett-Luce models:\n")
  cat("\n")
  cat("No. of items:",K,"\n")
  cat("No. of mixture components:",G,"\n")
  cat("No. of saved MCMC samples:",L,"\n")
  cat("\n")
  cat("Max. posterior log-likelihood:",max(gsPLMIX_out$log_lik,na.rm=TRUE),"\n")
  cat("Min. posterior deviance:",min(gsPLMIX_out$deviance,na.rm=TRUE),"\n")
  cat("Max. objective function:",max(gsPLMIX_out$objective,na.rm=TRUE),"\n")
  cat("\n")
  cat("Use functions summary() and plot() to summarize and visualize the object of class 'gsPLMIX'.")

}

print.summary.gsPLMIX <- function(x,...){
  #/' Print of the summary of Gibbs sampling simulation of a Bayesian mixture of Plackett-Luce models.
  #/' 
  #/' \code{print} method for class \code{summary.gsPLMIX}. It shows some general information on the Gibbs sampling simulation of a Bayesian mixture of Plackett-Luce models.
  #/'
  #/'  
  #/' @param x Object of class \code{summary.gsPLMIX} returned by the \code{summary.gibbsPLMIX} function.
  #/' @param ... Further arguments passed to or from other methods (not used).
  #/'
  #/'
  #/' @references 
  #/' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
  #/'
  #/' Mollica, C. and Tardella, L. (2014). Epitope profiling via mixture modeling for ranked data. \emph{Statistics in Medicine}, \bold{33}(21), pages 3738--3758, ISSN: 0277-6715, DOI: 10.1002/sim.6224.
  #/' 
  #/' @author Cristina Mollica and Luca Tardella
  
  summary.gsPLMIX_out=x
  
  if(class(summary.gsPLMIX_out)!="summary.gsPLMIX"){
    stop("The function requires an object of S3 class 'summary.gsPLMIX' as its first argument.")
  }
  
  cat("\nCall:\n", paste(deparse(summary.gsPLMIX_out$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Posterior statistics:\n")
  print(summary.gsPLMIX_out$statistics)
  cat("\n")
  cat("Quantiles:\n")
  print(summary.gsPLMIX_out$quantiles)
  cat("\n")
  pr=paste0(100*attr(summary.gsPLMIX_out[[grep("HPD",names(summary.gsPLMIX_out))]],"Probability"),"%")
  cat(pr,"HPD intervals:\n")
  attr(summary.gsPLMIX_out[[grep("HPD",names(summary.gsPLMIX_out))]],"Probability")=NULL
  print(summary.gsPLMIX_out[[grep("HPD",names(summary.gsPLMIX_out))]])
  cat("\n")
  cat("Estimated component-specific modal orderings:\n")
  print(summary.gsPLMIX_out$Modal_orderings)
  
} 



random_start <- function(mlesupp, givenweights, alpha=rep(1,G)){

#/' Appropriate simulation of starting values for tandom initialization of Gibbs Sampling. It start from the mle corresponding to no-group structure and then it randomly selects rescaled random support points (with sum 1) of G mixture components such that the marginal support coincides with the mle support for G=1

#/' Random generation of starting values of the component-specific support parameters for Gibbs sampling
#/' 
#/' @param mlesupp MLE of support parameters
#/' @param givenweights A numeric vector of \eqn{G} mixture weights
#/' @param alpha A numeric vector of \eqn{G} positive reals to be used as Dirichlet parameters for the random start which corresponds to a convex combination of \eqn{G} support parameter vertices
#/'
#/' @return \code{out} A numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with starting values of the component-specific support parameters  
#/' 
#/' @author Cristina Mollica and Luca Tardella
           K <- length(mlesupp)
           G <- length(givenweights)
           out <- matrix(NA,nrow=G,ncol=K)

if(G==1){
    out[1,] <- mlesupp
}else{
# for each component
# compute the H-representation
# transform it into the V-representation
# draw a random sample from the symplex

           for( j in 1:K ) {

               Aineq <- rbind(-diag(G))
               bineq <- c(rep(0, G))

               Aeq <- matrix(givenweights,nrow=1)
               beq <- mlesupp[j]

               hr <- makeH(Aineq,bineq,Aeq,beq)
               vr <- scdd(hr)

               Vertexes <- t(vr$output[,-c(1,2)]) # as column vectors
               myrandomcomponentwithconstrainedmean <- Vertexes%*%t(rdirichlet(1,alpha))
               out[,j] <- myrandomcomponentwithconstrainedmean
               
           }
       }
           
           return(out)

       }

#### Selection criteria

selectPLMIX_single <- function(pi_inv,G,
			      MCMCsampleP=NULL,
			      MCMCsampleW=NULL,
			      MAPestP,
			      MAPestW,
            deviance,
            post_est="mean"){
#/' Bayesian selection criteria for mixtures of Plackett-Luce models
#/' 
#/' Compute Bayesian comparison criteria for mixtures of Plackett-Luce models with a different number of components.
#/' 
#/' Two versions of DIC and BPIC are returned corresponding to two alternative ways of computing the penalty term: the former was proposed by Spiegelhalter et al. (2002) and is denoted with \code{pD}, whereas the latter was proposed by Gelman et al. (2004) and is denoted with \code{pV}. DIC2 coincides with AICM, that is, the Bayesian counterpart of AIC introduced by Raftery et al. (2007). 
#/' 
#/' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the MCMC samples of the component-specific support parameters.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#/' @param MAPestP Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of MAP component-specific support parameter estimates.
#/' @param MAPestW Numeric vector of the \eqn{G} MAP estimates of the mixture weights.
#/' @param deviance Numeric vector of posterior deviance values.
#/' @param post_est Character string indicating the  point estimates of the Plackett-Luce mixture parameters to be computed from the MCMC sample. This argument is ignored when MAP estimates are supplied in the \code{MAPestP} and \code{MAPestW} arguments. Default is \code{"mean"}. Alternatively, one can choose \code{"median"}.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{point_estP}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{(K+1)} matrix with the point estimates of the Plackett-Luce mixture parameters. The \eqn{(K+1)}-th column contains estimates of the mixture weights.}
#/'  \item{\code{point_estW}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{(K+1)} matrix with the point estimates of the Plackett-Luce mixture parameters. The \eqn{(K+1)}-th column contains estimates of the mixture weights.}
#/'  \item{\code{D_bar}}{ Posterior expected deviance.}
#/'  \item{\code{D_hat}}{ Deviance function evaluated at \code{point_est}.}
#/'  \item{\code{pD}}{ Effective number of parameters computed as \code{D_bar}-\code{D_hat}.}
#/'  \item{\code{pV}}{ Effective number of parameters computed as half the posterior variance of the deviance.}
#/'  \item{\code{DIC1}}{ Deviance Information Criterion with penalty term equal to \code{pD}.}
#/'  \item{\code{DIC2}}{ Deviance Information Criterion with penalty term equal to \code{pV}.}
#/'  \item{\code{BPIC1}}{ Bayesian Predictive Information Criterion obtained from \code{DIC1} by doubling its penalty term.}
#/'  \item{\code{BPIC2}}{ Bayesian Predictive Information Criterion obtained from \code{DIC2} by doubling its penalty term.}
#/'  \item{\code{BICM1}}{ Bayesian Information Criterion-Monte Carlo.}
#/'  \item{\code{BICM2}}{ Bayesian Information Criterion-Monte Carlo based on the actual MAP estimate given in the \code{MAPestP} and \code{MAPestW} arguments (unlike \code{BICM1}, no approximation of the MAP estimate from the MCMC sample).}
#/' 
#/' 
#/' @references 
#/' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#/'
#/' Ando, T. (2007). Bayesian predictive information criterion for the evaluation of hierarchical Bayesian and empirical Bayes models. \emph{Biometrika}, \bold{94}(2), pages 443--458.
#/'
#/' Raftery, A. E, Satagopan, J. M., Newton M. A. and Krivitsky, P. N. (2007). BAYESIAN STATISTICS 8. \emph{Proceedings of the eighth Valencia International Meeting 2006}, pages 371--416. Oxford University Press.
#/' 
#/' Gelman, A., Carlin, J. B., Stern, H. S. and Rubin, D. B. (2004). Bayesian data analysis. Chapman & Hall/CRC, Second Edition, ISBN: 1-58488-388-X. New York.
#/' 
#/' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{64}(4), pages 583--639.
#/' 
#/' @author Cristina Mollica and Luca Tardella

  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
  N <- nrow(pi_inv)
  K <- ncol(pi_inv)

  D_bar <- mean(deviance)

  if(!is.null(MAPestP) & !is.null(MAPestW)){
  	   point_estP <- MAPestP
  	   point_estW <- MAPestW
  }else{
     if(post_est=="mean"){
	         point_estP <- matrix(colMeans(MCMCsampleP),G,K)  	
   	       point_estW <- colMeans(MCMCsampleW)
       }else{
           point_estP <- matrix(apply(MCMCsampleP,2,FUN=median),G,K)  	
           point_estW <- apply(MCMCsampleW,2,FUN=median)
  	   }
  }
  
  	rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)
  	D_hat <- -2*loglikPLMIX(p=point_estP,weights=point_estW,ref_order=rho,pi_inv=pi_inv)
  
  pD <- D_bar-D_hat
  pV <- var(deviance)/2

  return(list(point_estP=point_estP,point_estW=point_estW,D_bar=D_bar,D_hat=D_hat,pD=pD,pV=pV,DIC1=D_bar+pD,DIC2=D_bar+pV,
              BPIC1=D_bar+2*pD,BPIC2=D_bar+2*pV,BICM1=D_bar+pV*(log(x=N)-1),BICM2=D_hat+pV*log(x=N)))
}

selectPLMIX <- function(pi_inv,seq_G,
			      MCMCsampleP=vector(mode="list",length=length(seq_G)),
			      MCMCsampleW=vector(mode="list",length=length(seq_G)),
			      MAPestP,
			      MAPestW,
            deviance,
            post_est="mean",
            parallel=FALSE){
#' Bayesian selection criteria for mixtures of Plackett-Luce models
#' 
#' Compute Bayesian comparison criteria for mixtures of Plackett-Luce models with a different number of components.
#'
#' The \code{selectPLMIX} function privileges the use of the MAP point estimates to compute the Bayesian model comparison criteria, since they are not affected by the label switching issue. By setting both the \code{MAPestP} and \code{MAPestW} arguments equal to NULL, the user can alternatively compute the selection measures by relying on a different posterior summary (\code{"mean"} or \code{"median"}) specified in the \code{post_est} argument. In the latter case, the MCMC samples for each Plackett-Luce mixture must be supplied in the lists \code{MCMCsampleP} and \code{MCMCsampleW}. The drawback when working with point estimates other than the MAP is that the possible presence of label switching has to be previously removed from the traces to obtain meaningful results. See the \code{\link{label_switchPLMIX}} function to perfom label switching adjustment of the MCMC samples. 
#' 
#' Several model selection criteria are returned. The two versions of DIC correspond to alternative ways of computing the effective number of parameters: DIC1 was proposed by Spiegelhalter et al. (2002) with penalty named \code{pD}, whereas DIC2 was proposed by Gelman et al. (2004) with penalty named \code{pV}. The latter coincides with the AICM introduced by Raftery et al. (2007), that is, the Bayesian counterpart of AIC. BPIC1 and BPIC2 are obtained from the two DIC by simply doubling the penalty term, as suggested by Ando (2007) to contrast DIC's tendency to overfitting. BICM1 is the Bayesian variant of the BIC, originally presented by Raftery et al. (2007) and entirely based on the MCMC sample. The BICM2, instead, involved the MAP estimate without the need of its approximation from the MCMC sample as for the BICM1.
#'   
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @param seq_G Numeric vector with the number of components of the Plackett-Luce mixtures to be compared.
#' @param MCMCsampleP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the MCMC samples of the component-specific support parameters. Default is list of \code{NULL} elements.
#' @param MCMCsampleW List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights. Default is list of \code{NULL} elements.
#' @param MAPestP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters.
#' @param MAPestW List of size \code{length(seq_G)}, whose generic element is a numeric vector with the MAP estimates of the \eqn{G} mixture weights.
#' @param deviance List of size \code{length(seq_G)}, whose generic element is a numeric vector of posterior deviance values.
#' @param post_est Character string indicating the  point estimates of the Plackett-Luce mixture parameters to be computed from the MCMC sample. This argument is ignored when MAP estimates are supplied in the \code{MAPestP} and \code{MAPestW} arguments. Default is \code{"mean"}. Alternatively, one can choose \code{"median"} (see 'Details').
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{point_estP}}{ List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the point estimates of the component-specific support parameters employed for the computation of the criteria.}
#'  \item{\code{point_estW}}{ List of size \code{length(seq_G)}, whose generic element is a numeric vector with the \eqn{G} point estimates of the mixture weights employed for the computation of the criteria.}
#'  \item{\code{fitting}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix with the fitting terms of the comparison measures, given by the posterior expected deviance \code{D_bar} and the deviance \code{D_hat} evaluated at the point estimate.}
#'  \item{\code{penalties}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix with the penalty terms \code{pD} and \code{pV} (effective number of parameters).}
#'  \item{\code{criteria}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{6} matrix of Bayesian model selection criteria: \code{DIC1}, \code{DIC2}, \code{BPIC1}, \code{BPIC2}, \code{BICM1} and \code{BICM2} (see 'Details').}
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' Ando, T. (2007). Bayesian predictive information criterion for the evaluation of hierarchical Bayesian and empirical Bayes models. \emph{Biometrika}, \bold{94}(2), pages 443--458.
#'
#' Raftery, A. E, Satagopan, J. M., Newton M. A. and Krivitsky, P. N. (2007). BAYESIAN STATISTICS 8. \emph{Proceedings of the eighth Valencia International Meeting 2006}, pages 371--416. Oxford University Press.
#' 
#' Gelman, A., Carlin, J. B., Stern, H. S. and Rubin, D. B. (2004). Bayesian data analysis. Chapman & Hall/CRC, Second Edition, ISBN: 1-58488-388-X. New York.
#' 
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{64}(4), pages 583--639.
#'
#' @author Cristina Mollica and Luca Tardella
#' @examples
#' 
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' 
#' ## Fit 1- and 2-component PL mixtures via MAP estimation 
#' MAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=1, 
#'                                    n_start=2, n_iter=400*1)
#' 
#' MAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=2, 
#'                                    n_start=2, n_iter=400*2)
#' 
#' mcmc_iter <- 30
#' burnin <- 10
#' 
#' ## Fit 1- and 2-component PL mixtures via Gibbs sampling procedure
#' GIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_1$mod$P_map,
#'                       z=binary_group_ind(MAP_1$mod$class_map,G=1)))
#' GIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_2$mod$P_map,
#'                       z=binary_group_ind(MAP_2$mod$class_map,G=2)))

#' ## Select the optimal number of components 
#' SELECT <- selectPLMIX(pi_inv=d_carconf, seq_G=1:2, 
#'                       MAPestP=list(MAP_1$mod$P_map, MAP_2$mod$P_map), 
#'                       MAPestW=list(MAP_1$mod$W_map, MAP_2$mod$W_map), 
#'                       deviance=list(GIBBS_1$deviance, GIBBS_2$deviance))
#' SELECT$criteria
#' 
#' @export 
  
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
  ncomp <- length(seq_G)

	if(!parallel){
		
	  selection <- vector(mode="list",length=ncomp)
		
	  for(l in 1:ncomp){
		
		print(paste("SELECTION CRITERIA FOR G=",seq_G[l]))
	    selection[[l]] <- selectPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],deviance=deviance[[l]],post_est=post_est)
	  }

		
	}else{
		
		
	  selection <- foreach(l=1:ncomp) %dopar%{   
      tempselection <- selectPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
               MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
               MAPestW=MAPestW[[l]],deviance=deviance[[l]],post_est=post_est)
          }
		
  }
  
  
  point_estP <- sapply(selection,"[[","point_estP")
  point_estW <- sapply(selection,"[[","point_estW")
  fitting <- t(sapply(lapply(selection,"[",c("D_bar","D_hat")),unlist))
  effective_numer_of_parameters <- t(sapply(lapply(selection,"[",c("pD","pV")),unlist))
  criteria <- t(sapply(lapply(selection,"[",c("DIC1","DIC2","BPIC1","BPIC2","BICM1","BICM2")),unlist))

  names(point_estP) <- names(point_estW) <- rownames(fitting) <- rownames(effective_numer_of_parameters) <- rownames(criteria) <- paste0("G_",seq_G)                           
    
  out <- list(point_estP=point_estP,point_estW=point_estW,fitting=fitting,
           effective_numer_of_parameters=effective_numer_of_parameters,criteria=criteria)
           
  return(out)
              
}

#### Label switching adjustment

label_switchPLMIX_single <- function(pi_inv,G,
                              MCMCsampleP,
                              MCMCsampleW,
                              MAPestP,
                              MAPestW){ 
#/' Label switching adjustment for mixtures of Plackett-Luce models
#/' 
#/' Remove the label switching phenomenon from the MCMC samples of Bayesian mixtures of Plackett-Luce models with a different number of components.
#/' 
#/' The \code{label_switchPLMIX} function performs the label switching adjustment of the MCMC samples via the Pivotal Reordering Algorithm (PRA) described in Marin et al (2005), by recalling the \code{\link[label.switching]{pra}} function from the \code{\link[label.switching]{label.switching}} package.
#/' 
#/' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the MCMC samples of the component-specific support parameters to be processed.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights to be processed.
#/' @param MAPestP Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of MAP component-specific support parameter estimates to be used as pivot in the PRA method.
#/' @param MAPestW Numeric vector of the \eqn{G} MAP estimates of the mixture weights as pivot in the PRA method.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{final_sampleP}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K}\eqn{\times}{x}\eqn{L} array MCMC samples of the component-specific support parameters adjusted for label switching.}
#/'  \item{\code{final_sampleW}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix of MCMC samples of the mixture weights adjusted for label switching.}
#/' 
#/' @author Cristina Mollica and Luca Tardella
 
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
  N <- nrow(pi_inv)
  K <- ncol(pi_inv)
  L <- nrow(MCMCsampleW)
  
  mcmc.sample <- array(cbind(MCMCsampleP,MCMCsampleW),c(L,G,(K+1)))
  
  if(G==1){
    reordered.pra <- list(output=NULL)
    reordered.pra$output <- mcmc.sample
  }else{
    print("LABEL SWITCHING ADJUSTMENT WITH PIVOTAL REORDERING ALGORITHM")
    pivot.input <- cbind(MAPestP,MAPestW)
    lab.pra <- pra(mcmc.pars=mcmc.sample,pivot=pivot.input)
    reordered.pra <- permute.mcmc(mcmc=mcmc.sample,permutations=lab.pra$permutations)
  }
  
  final.sample <- matrix(reordered.pra$output,nrow=L,ncol=G*(K+1))
  final_sampleP <- array(t(final.sample[,1:(G*K)]),c(G,K,L))
  final_sampleW <- final.sample[,-c(1:(G*K)),drop=FALSE]
  
  out <- list(final_sampleP=final_sampleP,final_sampleW=final_sampleW)
  
  return(out)
  
}


label_switchPLMIX <- function(pi_inv,seq_G,
                       MCMCsampleP,
                       MCMCsampleW,
                       MAPestP,
                       MAPestW,
                       parallel=FALSE){
#' Label switching adjustment  of the Gibbs sampling simulations for Bayesian mixtures of Plackett-Luce models
#' 
#' Remove the label switching phenomenon from the MCMC samples of Bayesian mixtures of Plackett-Luce models with \eqn{G>1} components.
#' 
#' The \code{label_switchPLMIX} function performs the label switching adjustment of the MCMC samples via the Pivotal Reordering Algorithm (PRA) described in Marin et al (2005), by recalling the \code{\link[label.switching]{pra}} function from the \code{\link[label.switching]{label.switching}} package.
#' 
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @param seq_G Numeric vector with the number of components of the Plackett-Luce mixtures to be assessed.
#' @param MCMCsampleP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the MCMC samples of the component-specific support parameters to be processed.
#' @param MCMCsampleW List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights to be processed.
#' @param MAPestP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix with the MAP estimates of the component-specific support parameters to be used as a pivot in the PRA method (see 'Details').
#' @param MAPestW List of size \code{length(seq_G)}, whose generic element is a numeric vector with the MAP estimates of the \eqn{G} mixture weights to be used as a pivot in the PRA method (see 'Details').
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list of named objects:
#' 
#'  \item{\code{final_sampleP}}{ List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{G}\eqn{\times}{x}\eqn{K}\eqn{\times}{x}\eqn{L} array with the MCMC samples of the component-specific support parameters adjusted for label switching.}
#'  \item{\code{final_sampleW}}{ List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights adjusted for label switching.}
#' 
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#' 
#' Papastamoulis, P. (2016). label.switching: An R Package for Dealing with the Label Switching Problem in MCMC Outputs. \emph{Journal of Statistical Software}, \bold{69}(1), pages 1--24, DOI: 10.18637/jss.v069.c01.
#'   
#' Marin, J. M., Mengersen, K. and Robert, C.P. (2005). Bayesian modelling and inference on mixtures of distributions. \emph{Handbook of Statistics} (25), D. Dey and C.R. Rao (eds). Elsevier-Sciences.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link[label.switching]{pra}} 
#' 
#' @examples
#' 
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' 
#' ## Fit 1- and 2-component PL mixtures via MAP estimation
#' MAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=1, 
#'                              n_start=2, n_iter=400*1)
#' 
#' MAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=2, 
#'                              n_start=2, n_iter=400*2)
#'                                    
#' MAP_3 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=3, 
#'                              n_start=2, n_iter=400*3)
#'                                    
#' mcmc_iter <- 30
#' burnin <- 10
#' 
#' ## Fit 1- and 2-component PL mixtures via Gibbs sampling procedure
#' GIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_1$mod$P_map,
#'                       z=binary_group_ind(MAP_1$mod$class_map,G=1)))
#' GIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_2$mod$P_map,
#'                       z=binary_group_ind(MAP_2$mod$class_map,G=2)))
#' GIBBS_3 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=3, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_3$mod$P_map,
#'                       z=binary_group_ind(MAP_3$mod$class_map,G=3)))
#'                             
#' ## Adjusting the MCMC samples for label switching
#' LS <- label_switchPLMIX(pi_inv=d_carconf, seq_G=1:3, 
#'                    MCMCsampleP=list(GIBBS_1$P, GIBBS_2$P, GIBBS_3$P), 
#'                    MCMCsampleW=list(GIBBS_1$W, GIBBS_2$W, GIBBS_3$W), 
#'                    MAPestP=list(MAP_1$mod$P_map, MAP_2$mod$P_map, MAP_3$mod$P_map), 
#'                    MAPestW=list(MAP_1$mod$W_map, MAP_2$mod$W_map, MAP_3$mod$W_map))
#' str(LS)       
#' @export 
  
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
  ncomp <- length(seq_G)
  
  if(!parallel){
    
    adjust <- vector(mode="list",length=ncomp)
    
    
    for(l in 1:ncomp){
      adjust[[l]] <- label_switchPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                       MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
                                       MAPestW=MAPestW[[l]])
    }
    
    
  }else{
    
    
    adjust <- foreach(l=1:ncomp) %dopar%{   
      tempadjust <- label_switchPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                      MCMCsampleW=MCMCsampleW[[l]],MAPestP=MAPestP[[l]],
                                      MAPestW=MAPestW[[l]])
    }
    
  }
  
  # OLD final_sampleP <- sapply(adjust,"[[","final_sampleP")   
  # OLD final_sampleW <- sapply(adjust,"[[","final_sampleW")     
  
  final_sampleP <-  drop(simplify2array(simplify2array(lapply(adjust,function(x){lapply("final_sampleP",function(y)do.call("[[",list(x,y)))}))))
  
  final_sampleW <-  drop(simplify2array(simplify2array(lapply(adjust,function(x){lapply("final_sampleW",function(y)do.call("[[",list(x,y)))}))))
  
  if(length(seq_G)>1){
    names(final_sampleP) <- names(final_sampleW) <- paste0("G_",seq_G)
  }else{
    final_sampleP <- list(final_sampleP)
    final_sampleW <- list(final_sampleW)
    names(final_sampleP) <- names(final_sampleW) <- paste0("G_",seq_G)
  }
  
  out <- list(final_sampleP=final_sampleP,final_sampleW=final_sampleW)
  return(out)
  
}


#### Posterior predictive check

ppcheckPLMIX_single <- function(pi_inv,G,
			               MCMCsampleP,
			               MCMCsampleW,
						         top1=TRUE,			                 
						         paired=TRUE){ 
#/' Posterior predictive check for a mixture of Plackett-Luce models
#/' 
#/' Compute predictive posterior \eqn{p}-values based on top item and paired comparison frequencies to assess the goodness-of-fit of a Bayesian mixtures of Plackett-Luce models for partial orderings. 
#/' 
#/' In the case of partial orderings, the same missingness patterns of the observed dataset, i.e., the number of items ranked by each sample unit, are reproduced on the replicated datasets.
#/' 
#/' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the MCMC samples of the component-specific support parameters.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#/' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on top frequencies has to be computed. Default is \code{TRUE}.
#/' @param paired Logical: whether the posterior predictive \eqn{p}-value based on paired comparison frequencies has to be computed. Default is \code{TRUE}.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{post_pred_pvalue_top1}}{ If \code{top1} is \code{TRUE}, posterior predictive \eqn{p}-value based on top frequencies, otherwise \code{NULL}.}
#/'  \item{\code{post_pred_pvalue_paired}}{ If \code{paired} is \code{TRUE}, posterior predictive \eqn{p}-value based on paired comparison frequencies, otherwise \code{NULL}.}
#/' 
#/' @author Cristina Mollica and Luca Tardella
 
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
  N <- nrow(pi_inv)
  K <- ncol(pi_inv)
  L <- nrow(MCMCsampleW)
  

	final.sample <- cbind(MCMCsampleP,MCMCsampleW)
	final_sampleP <- array(c(t(MCMCsampleP)),c(G,K,L))
	final_sampleW <- MCMCsampleW

  pi_inv_int <- pi_inv
  mode(pi_inv_int) <- "integer"
  rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)

  if(top1){
  	
    print(paste("POSTERIOR PREDICTIVE CHECK FOR G=",G))
  	print("Top1 frequencies-based posterior predictive p-value")
  	chi.obs.top1 <- rep(NA,L)
	  chi.rep.top1 <- rep(NA,L)

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.top1[l] <- chisqmeasureobs1dim(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.top1[l] <- chisqmeasuretheo1dim(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	}

	post_pred_pvalue_top1 <- mean(chi.rep.top1 >= chi.obs.top1)	

  }else{
  	
	post_pred_pvalue_top1 <- NA

  }

  if(paired){

    print(paste("POSTERIOR PREDICTIVE CHECK FOR G=",G))
  	print("Paired comparison frequencies-based posterior predictive p-value")  	
  	chi.obs.paired <- rep(NA,L)
	  chi.rep.paired <- rep(NA,L)

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.paired[l] <- chisqmeasureobs(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.paired[l] <- chisqmeasuretheo(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	}

	post_pred_pvalue_paired <- mean(chi.rep.paired >= chi.obs.paired)	

  }else{
  	
	post_pred_pvalue_paired <- NA

  }
  	

  out <- list(post_pred_pvalue_top1=post_pred_pvalue_top1,post_pred_pvalue_paired=post_pred_pvalue_paired)
  
  return(out)
    
}


ppcheckPLMIX <- function(pi_inv,seq_G,
                         MCMCsampleP,
                         MCMCsampleW,
                         top1=TRUE,			                 
                         paired=TRUE,
                         parallel=FALSE){ 
#' Posterior predictive check for Bayesian mixtures of Plackett-Luce models
#' 
#' Perform posterior predictive check to assess the goodness-of-fit of Bayesian mixtures of Plackett-Luce models with a different number of components.
#' 
#' The \code{ppcheckPLMIX} function returns two posterior predictive \eqn{p}-values based on two chi squared discrepancy variables involving: (i) the top item frequencies and (ii) the paired comparison frequencies. In the presence of partial sequences in the \code{pi_inv} matrix, the same missingness patterns observed in the dataset (i.e., the number of items ranked by each sample unit) are reproduced on the replicated datasets from the posterior predictive distribution. 
#' 
#' 
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @param seq_G Numeric vector with the number of components of the Plackett-Luce mixtures to be assessed.
#' @param MCMCsampleP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the MCMC samples of the component-specific support parameters.
#' @param MCMCsampleW List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on the top item frequencies has to be computed. Default is \code{TRUE}.
#' @param paired Logical: whether the posterior predictive \eqn{p}-value based on the paired comparison frequencies has to be computed. Default is \code{TRUE}.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list with a named element:
#' 
#'  \item{\code{post_pred_pvalue}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix of posterior predictive \eqn{p}-values based on the top item and paired comparison frequencies. If either \code{top1} or \code{paired} argument is \code{FALSE}, the corresponding matrix entries are \code{NA}.}
#' 
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link{ppcheckPLMIX_cond}} 
#' 
#' @examples
#' 
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' 
#' ## Fit 1- and 2-component PL mixtures via MAP estimation
#' MAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=1, 
#'                              n_start=2, n_iter=400*1)
#' 
#' MAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=2, 
#'                              n_start=2, n_iter=400*2)
#'                                    
#' MAP_3 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=3, 
#'                              n_start=2, n_iter=400*3)
#'                                    
#' mcmc_iter <- 30
#' burnin <- 10
#' 
#' ## Fit 1- and 2-component PL mixtures via Gibbs sampling procedure
#' GIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_1$mod$P_map,
#'                       z=binary_group_ind(MAP_1$mod$class_map,G=1)))
#' GIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_2$mod$P_map,
#'                       z=binary_group_ind(MAP_2$mod$class_map,G=2)))
#' GIBBS_3 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=3, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_3$mod$P_map,
#'                       z=binary_group_ind(MAP_3$mod$class_map,G=3)))
#'                             
#' ## Checking goodness-of-fit of the estimated mixtures
#' CHECK <- ppcheckPLMIX(pi_inv=d_carconf, seq_G=1:3, 
#'                       MCMCsampleP=list(GIBBS_1$P, GIBBS_2$P, GIBBS_3$P), 
#'                       MCMCsampleW=list(GIBBS_1$W, GIBBS_2$W, GIBBS_3$W))
#' CHECK$post_pred_pvalue
#' 
#' @export 
  
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
  ncomp <- length(seq_G)
  
  if(!parallel){
    
    fitting <- vector(mode="list",length=ncomp)
    
    
    for(l in 1:ncomp){
      fitting[[l]] <- ppcheckPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                       MCMCsampleW=MCMCsampleW[[l]],top1=top1,paired=paired)
    }
    
    
  }else{
    
    
    fitting <- foreach(l=1:ncomp) %dopar%{   
      tempfitting <- ppcheckPLMIX_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                      MCMCsampleW=MCMCsampleW[[l]],top1=top1,paired=paired)
    }
    
  }
  
  post_pred_pvalue <- t(sapply(lapply(fitting,"[",c("post_pred_pvalue_top1","post_pred_pvalue_paired")),unlist))
  

  if(!is.numeric(post_pred_pvalue)){
    post_pred_pvalue <- matrix(NA,nrow=length(seq_G),ncol=2)
  }
  attributes(post_pred_pvalue) <- attributes(post_pred_pvalue)[c("dim","dimnames")]
  post_pred_pvalue <- as.matrix(post_pred_pvalue)
  
  rownames(post_pred_pvalue) <- paste0("G_",seq_G)
  
 
  out <- list(post_pred_pvalue=post_pred_pvalue)
  return(out)
  
}


ppcheckPLMIX_cond_single <- function(pi_inv,G,
			               MCMCsampleP,
			               MCMCsampleW,
						         top1=TRUE,			                 
						         paired=TRUE){ 
#/' Conditional predictive posterior \eqn{p}-values
#/' 
#/' Compute conditional predictive posterior \eqn{p}-values based on top paired comparison frequencies to assess the goodness-of-fit of a Bayesian mixtures of Plackett-Luce models for partial orderings. 
#/' 
#/' In the case of partial orderings, the same missingness patterns of the observed dataset, i.e., the number of items ranked by each sample unit, are reproduced on the replicated datasets.
#/' 
#/' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#/' @param G Number of mixture components.
#/' @param MCMCsampleP Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the MCMC samples of the component-specific support parameters.
#/' @param MCMCsampleW Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#/' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on top frequencies has to be computed. Default is \code{TRUE}.
#/' @param paired Logical: whether the posterior predictive \eqn{p}-value based on paired comparison frequencies has to be computed. Default is \code{TRUE}.
#/'
#/' @return A list of named objects:
#/' 
#/'  \item{\code{post_pred_pvalue_top1}}{ If \code{top1} is \code{TRUE}, posterior predictive \eqn{p}-value based on top frequencies, otherwise \code{NULL}.}
#/'  \item{\code{post_pred_pvalue_paired}}{ If \code{paired} is \code{TRUE}, posterior predictive \eqn{p}-value based on paired comparison frequencies, otherwise \code{NULL}.}
#/' 
#/' @author Cristina Mollica and Luca Tardella
 
  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
  N <- nrow(pi_inv)
  K <- ncol(pi_inv)
  L <- nrow(MCMCsampleW)
  
	final.sample <- cbind(MCMCsampleP,MCMCsampleW)
	final_sampleP <- array(c(t(MCMCsampleP)),c(G,K,L))
	final_sampleW <- MCMCsampleW

  pi_inv_int <- pi_inv
  mode(pi_inv_int) <- "integer"
  rho <- matrix(1:K,nrow=G,ncol=K,byrow=TRUE)

  if(top1){
  	
    print(paste("CONDITIONAL POSTERIOR PREDICTIVE CHECK FOR G=",G))
  	print("Conditional top1 frequencies-based posterior predictive p-value")
  	chi.obs.top1.cond <- rep(NA,L)
	  chi.rep.top1.cond <- rep(NA,L)

  	chi.obs.top1.mat <- array(NA,dim=c(K,K,L))
	  chi.rep.top1.mat <- array(NA,dim=c(K,K,L))

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.top1.mat[,,l] <- chisqmeasureobsmatrix1dim(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.top1.mat[,,l] <- chisqmeasuretheomatrix1dim(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	 chi.obs.top1.cond[l] <- sum(chi.obs.top1.mat[,,l])
  	 chi.rep.top1.cond[l] <- sum(chi.rep.top1.mat[,,l])
  	}

	post_pred_pvalue_top1_cond <- mean(chi.rep.top1.cond >= chi.obs.top1.cond)	

  }else{
  	
	post_pred_pvalue_top1_cond <- NA

  }

  if(paired){

    print(paste("CONDITIONAL POSTERIOR PREDICTIVE CHECK FOR G=",G))
  	print("Conditional paired comparison frequencies-based posterior predictive p-value")  	
  	chi.obs.paired.cond=rep(NA,L)
	  chi.rep.paired.cond=rep(NA,L)

  	for(l in 1:L){
     (if((l%%200)==0) print(l))
     chi.obs.paired.cond[l] <- chisqmeasureobscond(pi_inv_int, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,])
     chi.rep.paired.cond[l] <- chisqmeasuretheocond(N,ref_order=rho, p=matrix(final_sampleP[,,l],nrow=G), weights=final_sampleW[l,],pi_inv_int)
  	}

	post_pred_pvalue_paired_cond <- mean(chi.rep.paired.cond >= chi.obs.paired.cond)	

  }else{
  	
	post_pred_pvalue_paired_cond <- NA

  }

  out <- list(post_pred_pvalue_top1_cond=post_pred_pvalue_top1_cond,post_pred_pvalue_paired_cond=post_pred_pvalue_paired_cond)
  
  return(out)
    
}


ppcheckPLMIX_cond <- function(pi_inv,seq_G,
			               MCMCsampleP,
			               MCMCsampleW,
						         top1=TRUE,			                 
						         paired=TRUE,
						         parallel=FALSE){ 
#' Conditional posterior predictive check for Bayesian mixtures of Plackett-Luce models
#' 
#' Perform conditional posterior predictive check to assess the goodness-of-fit of Bayesian mixtures of Plackett-Luce models with a different number of components. 
#' 
#' The \code{ppcheckPLMIX_cond} function returns two posterior predictive \eqn{p}-values based on two chi squared discrepancy variables involving: (i) the top item frequencies and (ii) the paired comparison frequencies. In the presence of partial sequences in the \code{pi_inv} matrix, the same missingness patterns observed in the dataset (i.e., the number of items ranked by each sample unit) are reproduced on the replicated datasets from the posterior predictive distribution. Differently from the \code{ppcheckPLMIX} function, the condional discrepancy measures are obtained by summing up the chi squared discrepancies computed on subsamples of observations with the same number of ranked items.
#' 
#' 
#' @param pi_inv An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.
#' @param seq_G Numeric vector with the number of components of the Plackett-Luce mixtures to be assessed.
#' @param MCMCsampleP List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{(G*K)} matrix with the MCMC samples of the component-specific support parameters.
#' @param MCMCsampleW List of size \code{length(seq_G)}, whose generic element is a numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights.
#' @param top1 Logical: whether the posterior predictive \eqn{p}-value based on the top item frequencies has to be computed. Default is \code{TRUE}.
#' @param paired Logical: whether the posterior predictive \eqn{p}-value based on the paired comparison frequencies has to be computed. Default is \code{TRUE}.
#' @param parallel Logical: whether parallelization should be used. Default is \code{FALSE}.
#'
#' @return A list with a named element:
#' 
#'  \item{\code{post_pred_pvalue_cond}}{ Numeric \code{length(seq_G)}\eqn{\times}{x}\eqn{2} matrix of posterior predictive \eqn{p}-values based on the top item and paired comparison frequencies. If either \code{top1} or \code{paired} argument is \code{FALSE}, the corresponding matrix entries are \code{NA}.}
#' 
#' 
#' @references 
#' Mollica, C. and Tardella, L. (2017). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika}, \bold{82}(2), pages 442--458, ISSN: 0033-3123, DOI: 10.1007/s11336-016-9530-0.
#'
#' @author Cristina Mollica and Luca Tardella
#' 
#' @seealso \code{\link{ppcheckPLMIX}} 
#' 
#' @examples
#' 
#' data(d_carconf)
#' K <- ncol(d_carconf)
#' 
#' ## Fit 1- and 2-component PL mixtures via MAP estimation
#' MAP_1 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=1, 
#'                              n_start=2, n_iter=400*1)
#' 
#' MAP_2 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=2, 
#'                              n_start=2, n_iter=400*2)
#'                                    
#' MAP_3 <- mapPLMIX_multistart(pi_inv=d_carconf, K=K, G=3, 
#'                              n_start=2, n_iter=400*3)
#'                                    
#' mcmc_iter <- 30
#' burnin <- 10
#' 
#' ## Fit 1- and 2-component PL mixtures via Gibbs sampling procedure
#' GIBBS_1 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=1, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_1$mod$P_map,
#'                       z=binary_group_ind(MAP_1$mod$class_map,G=1)))
#' GIBBS_2 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=2, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_2$mod$P_map,
#'                       z=binary_group_ind(MAP_2$mod$class_map,G=2)))
#' GIBBS_3 <- gibbsPLMIX(pi_inv=d_carconf, K=K, G=3, n_iter=mcmc_iter, 
#'                       n_burn=burnin, init=list(p=MAP_3$mod$P_map,
#'                       z=binary_group_ind(MAP_3$mod$class_map,G=3)))
#'                             
#' ## Checking goodness-of-fit of the estimated mixtures
#' CHECKCOND <- ppcheckPLMIX_cond(pi_inv=d_carconf, seq_G=1:3, 
#'                                MCMCsampleP=list(GIBBS_1$P, GIBBS_2$P, GIBBS_3$P), 
#'                                MCMCsampleW=list(GIBBS_1$W, GIBBS_2$W, GIBBS_3$W))
#' CHECKCOND$post_pred_pvalue
#' 
#' @export 

  if(class(pi_inv)[1]!="top_ordering"){
    if(class(pi_inv)[1]=="RankData"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="rankings"){
      pi_inv=as.top_ordering(data=pi_inv)
    }
    if(class(pi_inv)[1]=="matrix" | class(pi_inv)[1]=="data.frame"){
      pi_inv=as.top_ordering(data=pi_inv,format_input="ordering",aggr=FALSE)
    }
  }
  pi_inv <- fill_single_entries(data=pi_inv)
  ncomp <- length(seq_G)

	if(!parallel){
		
	  fitting <- vector(mode="list",length=ncomp)
		

	  for(l in 1:ncomp){
	    fitting[[l]] <- ppcheckPLMIX_cond_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                            MCMCsampleW=MCMCsampleW[[l]],top1=top1,paired=paired)
	  }

		
	}else{
		
		
	  fitting <- foreach(l=1:ncomp) %dopar%{   
      tempfitting <- ppcheckPLMIX_cond_single(pi_inv=pi_inv,G=seq_G[l],MCMCsampleP=MCMCsampleP[[l]],
                                           MCMCsampleW=MCMCsampleW[[l]],top1=top1,paired=paired)
          }
		
  }

	
  post_pred_pvalue_cond <- t(sapply(lapply(fitting,"[",c("post_pred_pvalue_top1_cond","post_pred_pvalue_paired_cond")),unlist))


  if(!is.numeric(post_pred_pvalue_cond)){
    post_pred_pvalue_cond <- matrix(NA,nrow=length(seq_G),ncol=2)
  }
  attributes(post_pred_pvalue_cond) <- attributes(post_pred_pvalue_cond)[c("dim","dimnames")]
  post_pred_pvalue_cond <- as.matrix(post_pred_pvalue_cond)
  

  rownames(post_pred_pvalue_cond) <- paste0("G_",seq_G)                           
    
  out <- list(post_pred_pvalue_cond=post_pred_pvalue_cond)
  
  return(out)

}


	

