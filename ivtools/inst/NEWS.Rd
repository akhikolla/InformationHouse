\name{NEWS}
\encoding{UTF-8}
\title{News for Package \pkg{ivtools}}
\section{Changes in version 2.3.0 (2019-06-11)}{
  \subsection{Visible changes}{
    \itemize{
      \item Correction of the help page for \code{ivglm}: 
      \eqn{h_i(\psi)=\textrm{expit}\{\hat{E}(Y|L_i,Z_i,X_i)-m^T(L_i)\psi X_i\}}
      replaced by
      \eqn{h_i(\psi)=\textrm{expit}[\textrm{logit}\{\hat{E}(Y|L_i,Z_i,X_i)\}-m^T(L_i)\psi X_i]}.
    }
  }
  \subsection{Non-visible changes}{
    \itemize{
      \item Less internal calls to \code{coxph.detail} in \code{stdCoxph},
      which makes the function slightly faster.
      \item Replaced \code{class(fit)=="name"} with 
      \code{inherits(x=fit, what="name")} in places, which should be safer. 
      \item Minor changes in the internal handling of formula arguments,
      to avoid errors if formulas are given implicitly to models,
      e.g. \code{f <- y~x}, \code{fit <- glm(formula=f)}
    }
  } 
}
