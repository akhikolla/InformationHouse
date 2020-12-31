%% File Name: immer-package.Rd
%% File Version: 0.55

\name{immer-package}
\alias{immer-package}
\alias{immer}
\docType{package}

\title{
\packageTitle{immer}
}

\description{
\packageDescription{immer}
}

\author{
\packageAuthor{immer}

Maintainer: \packageMaintainer{immer}
}


\details{

% URL: \tab \url{https://sites.google.com/site/alexanderrobitzsch/software} \cr

The \pkg{immer} package has following features: 

\itemize{   

%% Hierarchical rater model
\item Estimation of the hierarchical rater model (Patz et al., 2002) with 
\code{\link{immer_hrm}} and simulation of it with \code{\link{immer_hrm_simulate}}.

%% conditional maximum likelihood estimation
\item The linear logistic partial credit model as an extension to the 
linear logistic test model (LLTM) for dichotomous data can be estimated with
conditional maximum likelihood (Andersen, 1995) using \code{\link{immer_cml}}.

%% Composite conditional maximum likelihood estimation 
\item The linear logistic partial credit model can be estimated with
composite conditional maximum likelihood (Varin, Reid & Firth, 2011) using the
\code{\link{immer_ccml}} function. 

%% joint maximum likelihood estimation
\item The linear logistic partial credit model can be estimated with a bias-corrected 
joint maximum likelihood method (Bertoli-Barsotti, Lando & Punzo, 2014)
using the \code{\link{immer_jml}} function.  

%% wrapper function
\item Wrapper function \code{\link{immer_FACETS}} to the commercial
program FACETS (Linacre, 1999) for analyzing multi-faceted Rasch models.

\item \dots

}

}


%Federal Institute for Education Research, Innovation and Development of the 
%Austrian School System (BIFIE Salzburg), Austria


\references{
Andersen, E. B. (1995). Polytomous Rasch models and their estimation.
In G. H. Fischer & I. W. Molenaar (Eds.). \emph{Rasch Models} (pp. 39-52).
New York: Springer.    

Bertoli-Barsotti, L., Lando, T., & Punzo, A. (2014). Estimating a Rasch Model via 
fuzzy empirical probability functions. In D. Vicari, A. Okada, G. Ragozini & 
C. Weihs (Eds.). 
\emph{Analysis and Modeling of Complex Data in Behavioral and Social Sciences}, 
Springer.

Linacre, J. M. (1999). \emph{FACETS} (Version 3.17)[Computer software]. Chicago: MESA.

Patz, R. J., Junker, B. W., Johnson, M. S., & Mariano, L. T. (2002). 
The hierarchical rater model for rated test items and its 
application to large-scale educational assessment data. 
\emph{Journal of Educational and Behavioral Statistics, 27}(4), 341-384.

Robitzsch, A., & Steinfeld, J. (2018). Item response models for human ratings: Overview, 
estimation methods, and implementation in R. 
\emph{Psychological Test and Assessment Modeling, 60}(1), 101-139.

Varin, C., Reid, N., & Firth, D. (2011). An overview of composite likelihood methods. 
\emph{Statistica Sinica, 21}, 5-42.

Wang, W. C., Su, C. M., & Qiu, X. L. (2014). Item response models for local dependence 
among multiple ratings. \emph{Journal of Educational Measurement, 51}(3), 260-280.
}



\seealso{
For estimating the Rasch multi-facets model with marginal
maximum likelihood see also the
\code{\link[TAM:tam.mml.mfr]{TAM::tam.mml.mfr}} and
\code{\link[sirt:rm.facets]{sirt::rm.facets}} functions.

For estimating the hierarchical rater model based on signal
detection theory see \code{\link[sirt:rm.sdt]{sirt::rm.sdt}}.

For conditional maximum likelihood estimation of linear logistic
partial credit models see the \pkg{eRm} (e.g. \code{\link[eRm:LPCM]{eRm::LPCM}})
and the \pkg{psychotools} (e.g. \code{\link[psychotools:pcmodel]{psychotools::pcmodel}})
packages.
}


\examples{   
  ##  
  ##  
  ##  immer 0.0-16 (2015-06-27)
  ##  
  ##   
  ##      #@#################################################################@   
  ##     :#@##################################################################`  
  ##     #####################################################################+  
  ##     #####################################################################@  
  ##    .##########@###########################################################` 
  ##    +######@+;'@##########################################################@: 
  ##    ######'     +@########################################################## 
  ##    #####+ ##'  `##########################################################@ 
  ##    ####@ ####  +@#@##################################################@##### 
  ##   .#####.#####@######@#########@#@######@@####@##@############@######@#####`
  ##   :####+:###@:,@##@,;##@+@##@+'###@;.'##@##@##@''@#####@@'+@###,;##@;######.
  ##   ;####@ @##   ##@  ,;    :,    .#    ':    :.    ,@@#'     ';   ,+   @####.
  ##   '#####`:##   #@   : ''    ++   @@;  ` +;    +'   #@' `+@   #@  .    #####,
  ##   +######, #   : ,   @#@   +#@  `##,   @#@   @##  `#@  `##   @#`  +,  @####:
  ##   +#######@@  :+@`  ;##'  `##@  .@#   +##:  `###  .@.  @#+   #:  .##@@#####:
  ##   '########:  ##@   @#@,  ;##;  ;#@   @##`  +##,  +#       `@#`  ##########,
  ##   ;########   @##   ###`  @##`  @#+  `###   @##   ##   ,;'@##@   @#########,
  ##   :#######@  `##+  ,##@   @#@   @#,  '##@   ##@   @@  `#######  `##########.
  ##   .#######+  ,##.  @##+  `##@  .##`  @##;  ,##@  `##   ###; #;  :#@,+###### 
  ##    #######'    #   @##:  ,##;  ;##   @##`  +###   `@:     `;#.  @@:  ###### 
  ##    #######@,`,##.`,###,..###,  #@#..,##@,..####:` +#@'  ``#+#` `@#+ :#####@ 
  ##    #########################:   ,###################@#####@###@#@#########@ 
  ##    '########################@+;+##########################################: 
  ##    ,######################################################################` 
  ##     ##@##################################################################@  
}

\keyword{package}
