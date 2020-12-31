\name{NEWS}
\alias{NEWS}

\title{News for \R Package \pkg{FENmlm}}

\section{Changes in version 2.4.2 (2019-05-16)}{
	\subsection{Bug fixes}{
		\itemize{
			\item When logit/negbin estimation with 2+ clusters and 2+ cores was made, if the user asked to stop the process, R could crash: now corrected.
			\item Now dates can be used as cluster variables.
			\item Small bug that prevented the obtention of the right cluster-coefficients for the Poisson and Gaussian models with 2 clusters (important only for small samples).
			\item For 2+ cluster convergence, one rare problem could appear because of numerical precision. Corrected.
		}
	}
	\subsection{Other}{
		\itemize{
			\item Better OpenMP handling.
			\item New convergence criterion that better handles large number cases.
			\item More precise warning messages.
			\item Modification of the trade data: now the variable Euros is scaled in million euros.
		}
	}
}


\section{Changes in version 2.4.1 (2019-02-25)}{
	\subsection{Vignette}{
		\itemize{
			\item Correction of a bug introduced after the modification of the \code{getFE} function.
		}
	}
}

\section{Changes in version 2.4 (2019-02-07)}{

	\subsection{Major changes}{
		\itemize{
			\item Major improvement of the internal code. All functions relating to the fixed-effects have been ported to c++, significantly improving performance (both in computing time and memory usage).
		}
	}

	\subsection{Other improvements}{
		\itemize{
			\item In the \code{print} method, now the BIC and likelihood numbers are properly aligned.
			\item[res2table] The titles of each column are more consistent.
			\item[femlm] New controls to check that the formula given by the user is valid and can be evaluated.
			\item[femlm] Better controls when using \code{na.rm = TRUE}.
			\item[res2tex] Argument \dQuote{append} has been superseded by \dQuote{replace} for improved clarity (the former can still be used though).
		}
	}

	\subsection{Function getFE}{
		\itemize{
			\item Now \code{getFE} returns a list of simple vectors (the fixed-effects on each dimension) instead of \code{femlm.cluster} objects.
			\item The \code{print} method has been superseded by a \code{summary} method.
		}
	}

	\subsection{Bug fixes}{
		\itemize{
			\item[femlm] For estimations with fixed-effects only (no variable), the print method did not work. Now corrected.
			\item[femlm] For estimations in Poisson with 2 fixed-effects, when convergence was "difficult", a bug could emerge. Now corrected.
			\item[femlm] Bug when there were NAs in the data and a factor variable in the right hand side.
			\item[print.femlm] For estimations without variables (only the fixed-effects), the print method now works properly.
			\item[diagnostic] Now function \code{diagnostic} works properly when there is a variable that is factored in the formula (i.e. as.factor(var) is applied) and not in the fixed-effect part.
		}
	}
}

\section{Changes in version 2.3.1 (2018-12-06)}{
	\subsection{New argument}{
		\itemize{
			\item[femlm] Function \code{femlm()} has the new argument \code{showWarning} whose default is \code{TRUE}. If \code{FALSE}, then warnings relating to a) the number of NA values, b) the clusters removed from estimation and c) the convergence state, are not prompted anymore.
			\item[res2tex] The new argument \code{powerBelow} controls the display of coefficients/standard-errors of low value in Latex. For values lower than \code{10**(powerBelow+1)} the display in Latex is now \code{x.xx $\\times 10^{pow}$}, for example \code{0.0000456} would be displayed as \code{4.56 $\\times 10^{-5}$}.
		}
	}

	\subsection{Improvements}{
		\itemize{
			\item Function \code{diagnostic()} now finds more involved collinearity cases.
			\item Function res2tex now provides a \dQuote{nicer} display of coefficients/standard-errors.
		}
	}
	\subsection{Bug fixes}{
		\itemize{
			\item Function \code{diagnostic()} did not work when some observations were removed from the estimation. This is now corrected.
			\item[femlm] For 3+ clusters, when finding the cluster coefficients was \sQuote{difficult}, there was a bug. Now corrected.
			\item[femlm] Now whenever the data of an estimation contains NA, \code{na.rm=TRUE} and no fixed-effects are used, the observations removed are correctly returned in the object \code{obsRemoved}. This was causing a bug in \code{summary} with custom clustered standard errors.
		}
	}
}

\section{Changes in version 2.3 (2018-11-15)}{
	\subsection{New arguments}{
		\itemize{
			\item Argument \dQuote{na.rm} is introduced. If the variables necessary for the estimation contain NAs and \code{na.rm = TRUE}, then all observations containing NA are removed prior to estimation and a warning message is raised detailing the number of observations removed. By default, if the necessary variables contain NA, an error is raised.
		}
	}

	\subsection{Deprecated arguments}{
		\itemize{
			\item Argument \dQuote{useAcc} is removed since now the acceleration algorithm is controlled dynamically.
		}
	}

	\subsection{Improvements}{
		\itemize{
			\item Improved algorithm to obtain the cluster coefficients for 3+ clusters. For difficult cases, convergence is faster.
			\item Now the acceleration algorithm (to get the cluster coefficients) is set up dynamically. Increases convergence speed overall.
			\item When there is no acceleration, the algorithm converges faster (for 2+ clusters).
			\item Important improvement of the algorithm getting the derivatives of the cluster coefficients. Convergence is fatser overall, especially for difficult cases.
			\item Algorithm has been improved when factor variables are used as clusters.
			\item Update method is now faster when using the same clusters and family.
		}
	}

	\subsection{Control}{
		\itemize{
			\item Better errors messages regarding standard-errors clustering.
			\item Better error messages when a cluster variable is not in the data.
			\item Amelioration of function \code{diagnostic}.
		}
	}

	\subsection{Vignette}{
		\itemize{
			\item New examples regarding the clustering of standard-errors.
		}
	}
}

\section{Changes in version 2.2.2 (2018-09-24)}{
	\subsection{Bug Fixes}{
		\itemize{
			\item [femlm/cores, reported by Valerio Leone Sciabolazza] Multi-core estimations did not work anymore since version 2.2.0. This has been corrected.
		}
	}
	\subsection{New arguments}{
		\itemize{
			\item Arguments \sQuote{itermax.cluster} and \sQuote{itermax.deriv} are added to the function femlm. They control the number of iterations in the step used to obtain the cluster coefficients and the step used to obtain the derivatives of the cluster coefficients.
		}
	}
	\subsection{Improvements}{
		\itemize{
			\item For 2+ clusters, the algorithm to obtain the cluster coefficients has been improved. Convergence is usually faster.
			\item For the Gaussian and Poisson models, and in the case of two clusters, the algorithm getting the cluster coefficients has been improved. Now it is much more efficient when dealing with difficult cases.
			\item For 2 clusters, the algorithm to obtain the derivative of the cluster coefficients has been improved. Convergence is usually faster for difficult cases.
			\item Information displayed when verbose > 0 is more detailed.
			\item Better handling of the dynamic precision leading to (usually) faster convergence.
		}
	}
	\subsection{Other}{
		\itemize{
			\item Dependence to package Matrix has been removed.
			\item Slight modifications in the vignette.
		}
	}
}

\section{Changes in version 2.2.1 (2018-09-15)}{
	\subsection{New features}{
		\itemize{
			\item Function diagnostic() is created to quickly spot where possible estimation problem come from. Ususally, perfect collinearity (between the explanatory variables or with the cluster variables) leads to false convergence of the optimization algorithm. This function helps figure out which variables that are collinear.
			\item Function obs2remove() is created to obtain the observations that need to be withdrawn from a Poisson, Negative Binomial or Logit estimation (because of only 0 or 1 outcomes for some clusters). The function femlm automatically removes such observations before performing the estimation; this function simply allows the user to clean the data before the estimation is done, if needed.
		}
	}
	\subsection{Bug fix}{
		\itemize{
			\item [femlm] Correction of a bug regarding the constuction of the cluster matrix when the clusters were complex variables in the formula and some observations had to be removed.
		}
	}
	\subsection{Vignette}{
		\itemize{
			\item Bug in displaying one table of the vignette has been corrected.
			\item Some parts of the vignette have been rewritten.
		}
	}
}

\section{Changes in version 2.2.0 (2018-08-03)}{
	\subsection{New features}{
		\itemize{
			\item [stats.methods] Several statistics methods have been implemented: nobs, AIB, BIC, logLik, coef, fitted, residuals, formula, model.matrix, predict, vcov, confint and update.
			\item [res2table] New arguments \sQuote{depvar} and \sQuote{family}: to control whether the dependent variables (resp. model families) should be displayed.
			\item [res2table] Now the \code{call} is used to produce more informative model names by default.
			\item [res2tex] New argument \sQuote{family} to control whether the families of the models should be displayed.
			\item [summary/vcov] Added the possibility to cluster the standard-errors with direct reference to a cluster variable used in the estimation. Example: \code{summary(est, se = "cluster", cluster = "var_2")}.
		}
	}
	\subsection{Improvements}{
		\itemize{
			\item [femlm] Formulas are more robustly handled with the package Formula.
			\item [femlm--negative binomial] Added a warning message when there is no overdispersion.
			\item [print.femlm] Information is now displayed more compactly
			\item Slightly faster internal functions.
		}
	}
	\subsection{Changing arguments}{
		\itemize{
			\item Some argument names in femlm have been modified for clarification:
			\itemize{
				\item [start => NL.start] Argument \code{start} now becomes \code{NL.start} to emphasize that it concerns the non-linear in parameters part only. This is the same for argument \code{start.init} which becomes \code{NL.start.init}.
				\item [debug => verbose] Argument \code{debug} of type logical now becomes \code{verbose} of type integer which can take values 0, 1 or 2 depending of the level of information to be reported.
			}
		}
	}
	\subsection{Documentation}{
		\itemize{
			\item Improved examples.
			\item A new reference is added which details the theory behind the method.
		}
	}
	\subsection{Bug}{
		\itemize{
			\item [femlm--data] Now the function properly handles matrices in argument \sQuote{data}.
		}
	}
	\subsection{Miscellaneous}{
		\itemize{
			\item OpenMP is handled differently (now used conditionnaly).
			\item CITATION is added.
			\item Some items in the femlm output have changed names to fit the \code{lm/glm} vocabuary more closely. This concerns: (old name => new name) coef => coefficients, expected.predictor => fitted.values, k => nparams, resids => residuals.
			\item A demo has been added.
			\item Warnings have been suppressed when the clustering of standard-errors was equivocal. They have been replaced by explicit mention of the cluster variables in the type of the standard-errors.
		}
	}
}

\section{Changes in version 2.1.2}{
	\subsection{Improvements}{
		\itemize{
			\item [femlm] Improved warning messages for overfitting problems.
			\item [femlm] Improved error-handling.
		}
	}
	\subsection{Documentation}{
		\itemize{
			\item Amelioration of the Vignette.
		}
	}
}

\section{Changes in version 2.1.1}{
	\subsection{New features}{
		\itemize{
			\item [femlm] Added argument \dQuote{precision.cluster} so that the user can set the level of precision manually (concerns the algorithm to obtain the cluster coefficients).
			\item [summary.femlm] For one-way clustering, the argument \dQuote{cluster} can now be a vector.
			\item A print method for the fixed-effects obtained with function \code{\link[FENmlm]{getFE}} has been implemented.
		}
	}
	\subsection{Documentation}{
		\itemize{
			\item A vignette has been created.
			\item Documentation has been added for fixed-effects plotting and for a new print method.
			\item Some examples have been updated.
		}
	}
	\subsection{Bug}{
		\itemize{
			\item [getFE] Corrected a small glitch that impeded the algorithm to find the right solution in presence of  complicated cases.
			\item [femlm] The function could not be evaluated at all when parallel::detectCores() gives NA. This has been corrected. Now when parallel::detectCores() gives NA, the number of cores cannot be set greater than 1.
		}
	}
	\subsection{Miscellaneous}{
		\itemize{
			\item [res2tex] In the exported Latex table the word \sQuote{Cluster} is changed to \sQuote{Fixed-Effects}.
			\item [summary.femlm] Using cluster=NULL does not cause problems any more.
			\item [print.femlm] Corrected a spelling mistake.
		}
	}

}

\section{Changes in version 2.1.0}{
	\subsection{New Features}{
		\itemize{
			\item New argument 'cores' => added multi cores with openMP. Particularly efficient for the Negative Binomial and the Logit models.
			\item Added a fixed-point acceleration algorithm for 3+ clusters and also to obtain the derivatives.
			\item Cluster coefficient precision is dynnamically adjusted for complex problems.
			\item  Faster matrix creation.
		}
	}
  \subsection{Other}{
		\itemize{
			\item [Design] Argument noWarning is dropped.
			\item [Design] clearer message when clusters are removed.
			\item [Bug] corrected undue information reporting.
			\item [Bug] corrected important bug that made some estimations fail.
			\item [Help] Added more accurate information on the returned objects.
		}
	}

}


\section{Changes in version 2.0.1}{
	\itemize{
	  \item [c++ code] Fixed overloading problem.
	  \item Improved help pages.
	  \item Code more robust to overfitting situations (in particular for the Logit case).
	  \item Added an acceleration algorithm to get the fixed-point for 2 clusters.
	  \item New collinearity warnings.
    }
  }

\section{Changes in version 2.0}{
  \subsection{New Features}{
    \itemize{
	  \item Algorithm is greatly improved.
      \item Two new likelihood models: Negative Binomial and Gaussian.
      \item Any number of clusters is supported.
      \item The functions res2table/res2tex to export the results are greatly enhanced.
      \item cluster-robust covariance: the method to compute it has been improved (more robust to rounding problems).
	  \item Computing standard errors: up to four-way clustering is supported.
    }
  }
  \subsection{New design}{
  	\itemize{
  		\item The function design has changed.
		\item There is only one function called femlm that performs \dQuote{normal} estimations and \dQuote{non.linear} in parameters estimations when adding the argument NL.fml.
		\item clusters can now directly be inserted in the formula without using the arument \dQuote{cluster}.
  }
  }
  \subsection{Argument name change}{
  	\itemize{
  		\item Name change: the argument \dQuote{dummy} becomes \dQuote{cluster}.
  }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item Bug related to the handling of the upper and lower bounds is now corrected.
    }
  }
}

