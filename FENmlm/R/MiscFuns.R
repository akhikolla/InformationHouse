



#' A print facility for \code{femlm} objects. It can compute different types of standard errors.
#'
#' This function is very similar to usual \code{summary} functions as it provides the table of coefficients along with other information on the fit of the estimation.
#'
#' @method print femlm
#'
#' @param x A femlm object. Obtained using \code{\link[FENmlm]{femlm}}.
#' @param n Integer, number of coefficients to display. By default, all estimated coefficients are displayed.
#' @param ... Other arguments to be passed to \code{\link[FENmlm]{vcov.femlm}}.
#'
#' @seealso
#' See also the main estimation functions \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade => we account for 3 cluster effects
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # displaying the results
#' print(est_pois)
#'
#' # with other type of standard error:
#' print(est_pois, se = "c")
#'
#'
print.femlm <- function(x, n, ...){

	x = summary(x, fromPrint = TRUE, ...)

	if(missing(n)){
		n = Inf
	} else if(!length(n) == 1 || !is.numeric(n) || n<=0){
		stop("Argument 'n' must be a single positive integer.")
	}

	if(!x$convStatus){
		message("The optimization algorithm did not converge, the results are not reliable.")
	}

	coeftable = x$coeftable

	# The type of SE
	se.type = attr(coeftable, "type")
	family_format = c(poisson="Poisson", negbin="Negative Binomial", logit="Logit", gaussian="Gaussian")

	msg = ifelse(is.null(x$call$NL.fml), "", "Non-linear ")
	cat(msg, "ML estimation, family = ", family_format[x$family], ", Dep. Var.: ", as.character(x$fml)[[2]], "\n", sep="")
	cat("Observations:", addCommas(x$n), "\n")
	if(!is.null(x$clusterSize)) cat("Cluster sizes: ", paste0(x$clusterNames, ": ", x$clusterSize, collapse=",  "), "\n", sep="")

	if(is.null(x$onlyCluster)){

		cat("Standard-errors type:", se.type, "\n")

		# The matrix of coefficients
		if(x$family=="negbin"){
			if(nrow(coeftable) == 2){
				new_table = coeftable[1, , FALSE]
			} else {
				new_table = coeftable[-nrow(coeftable), ]
			}

			myPrintCoefTable(head(new_table, n))

			theta = coeftable[".theta", 1]
			noDispInfo = ifelse(theta > 1000, "(theta >> 0, no sign of overdispersion, you may consider a Poisson model)", "")
			cat("Over-dispersion parameter: theta =", theta, noDispInfo, "\n")
		} else {
			myPrintCoefTable(head(coeftable, n))
		}
	}

	bic_ll = formatBicLL(BIC(x), x$loglik)
	cat("           BIC:", bic_ll$ll, "   Pseudo-R2:", round(x$pseudo_r2, 5), "\n")
	cat("Log-likelihood:", bic_ll$bic, "Squared Cor.:", round(x$sq.cor, 5), "\n")

	if(!x$convStatus && is.null(x$onlyCluster)){
		cat("# Evaluations:", x$iterations, "--", x$message, "\n")
	}

}

##

#' Summary of a \code{femlm} object. Computes different types of standard errors.
#'
#' This function is similar to \code{print.femlm}. It provides the table of coefficients along with other information on the fit of the estimation. It can compute different types of standard errors. The new variance covariance matrix is an object returned.
#'
#' @method summary femlm
#'
#' @param se Character scalar. Which kind of standard error should be computed: \dQuote{standard} (default), \dQuote{White}, \dQuote{cluster}, \dQuote{twoway}, \dQuote{threeway} or \dQuote{fourway}?
#' @param cluster A list of vectors. Used only if \code{se="cluster"}, \dQuote{se=twoway}, \dQuote{se=threeway} or \dQuote{se=fourway}. The vectors should give the cluster of each observation. Note that if the estimation was run using \code{cluster}, the standard error is automatically clustered along the cluster given in \code{\link[FENmlm]{femlm}}. For one-way clustering, this argument can directly be a vector (instead of a list). If the estimation has been done with cluster variables, you can give a character vector of the dimensions over which to cluster the SE.
#' @param object A femlm object. Obtained using \code{\link[FENmlm]{femlm}}.
#' @param dof_correction Logical, default is \code{FALSE}. Should there be a degree of freedom correction to the standard errors of the coefficients?
#' @param forceCovariance (Advanced users.) Logical, default is \code{FALSE}. In the peculiar case where the obtained Hessian is not invertible (usually because of collinearity of some variables), use this option force the covariance matrix, by using a generalized inverse of the Hessian. This can be useful to spot where possible problems come from.
#' @param keepBounded (Advanced users.) Logical, default is \code{FALSE}. If \code{TRUE}, then the bounded coefficients (if any) are treated as unrestricted coefficients and their S.E. is computed (otherwise it is not).
#' @param ... Not currently used.
#'
#' @return
#' It returns a \code{femlm} object with:
#' \item{cov.scaled}{The new variance-covariance matrix (computed according to the argument \code{se}).}
#' \item{se}{The new standard-errors (computed according to the argument \code{se}).}
#' \item{coeftable}{The table of coefficients with the new standard errors.}
#'
#' @seealso
#' See also the main estimation function \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade (with 3 cluster effects)
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # Comparing different types of standard errors
#' sum_white = summary(est_pois, se = "white")
#' sum_oneway = summary(est_pois, se = "cluster")
#' sum_twoway = summary(est_pois, se = "twoway")
#' sum_threeway = summary(est_pois, se = "threeway")
#'
#' res2table(sum_white, sum_oneway, sum_twoway, sum_threeway)
#'
#' # Alternative ways to cluster the SE:
#' \dontrun{
#' # two-way clustering: Destination and Product
#' summary(est_pois, se = "twoway", cluster = c("Destination", "Product"))
#' summary(est_pois, se = "twoway", cluster = list(trade$Destination, trade$Product))
#' }
#'
#'
summary.femlm <- function(object, se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), cluster, dof_correction=FALSE, forceCovariance = FALSE, keepBounded = FALSE, ...){
	# computes the clustered SD and returns the modified vcov and coeftable

	sd.val = match.arg(se)

	if(!is.null(object$onlyCluster)){
		# means that the estimation is done without variables
		return(object)
	}

	dots = list(...)

	# If cov.scaled exists => means that it has already been computed
	if(!is.null(object$cov.scaled) && "fromPrint" %in% names(dots)) return(object)

	# The new VCOV
	vcov = vcov(object, se=se, cluster=cluster, dof_correction=dof_correction, forceCovariance = forceCovariance, keepBounded = keepBounded, ...)

	sd2 = diag(vcov)
	sd2[sd2 < 0] = NA
	se = sqrt(sd2)

	# used to handle the case of bounded parameters
	params = names(object$coefficients)
	if(length(se) != length(params)){
		se = se[params]
	}
	names(se) = params

	# The coeftable is modified accordingly
	coeftable = object$coeftable

	# th z & p values
	zvalue <- object$coefficients/se
	pvalue <- 2*pnorm(-abs(zvalue))

	# update of se if bounded
	se_format = se
	isBounded = object$isBounded
	if(!is.null(isBounded) && any(isBounded)){
		if(!keepBounded){
			se_format[!isBounded] = decimalFormat(se_format[!isBounded])
			se_format[isBounded] = attr(isBounded, "type")
		}
	}

	# modifs of the table
	coeftable[, 2] = se_format
	coeftable[, 3] = zvalue
	coeftable[, 4] = pvalue

	attr(coeftable, "type") = attr(se, "type") = attr(vcov, "type")

	object$cov.scaled = vcov
	object$coeftable = coeftable
	object$se = se

	return(object)
}

#' Facility to export the results of multiple \code{femlm} estimations in a Latex table.
#'
#' This function aggregates the results of multiple estimations and display them in the form of  one Latex table whose rownames are the variables and the columns contain the coefficients and standard-errors.
#'
#' @inheritParams summary.femlm
#'
#' @param ... Used to capture different \code{\link[FENmlm]{femlm}} objects. Note that any other type of element is discarded. Note that you can give a list of \code{\link[FENmlm]{femlm}} objects.
#' @param digits Integer, default is 4. The number of digits to be displayed.
#' @param pseudo Logical, default is \code{TRUE}. Should the pseudo R2 be displayed?
#' @param title Character scalar. The title of the Latex table.
#' @param sdBelow Logical, default is \code{TRUE}. Should the standard-errors be displayed below the coefficients?
#' @param drop Character vector. This element is used if some variables are not to be displayed. This should be a regular expression (see \code{\link[base]{regex}} help for more info). There can be more than one regular expression. Each variable satisfying the regular expression will be discarded.
#' @param order Character vector. This element is used if the user wants the variables to be ordered in a certain way. This should be a regular expression (see \code{\link[base]{regex}} help for more info). There can be more than one regular expression. The variables satisfying the first regular expression will be placed first, then the order follows the sequence of regular expressions.
#' @param dict A named character vector. If provided, it changes the original variable names to the ones contained in the \code{dict}. Example: I want to change my variable named "a" to "$log(a)$" and "b3" to "$bonus^3$", then I used \code{dict=c(a="$log(a)$",b3="$bonus^3$")}.
#' @param file A character scalar. If provided, the Latex table will be saved in a file whose path is \code{file}.
#' @param replace Logical, default is \code{FALSE}. Only used if option \code{file} is used. Should the Latex table be written in a new file that replaces any existing file?
#' @param convergence Logical, default is missing. Should the convergence state of the algorithm be displayed? By default, convergence information is displayed if at least one model did not converge.
#' @param signifCode Named numeric vector, used to provide the significance codes with respect to the p-value of the coefficients. Default is \code{c("***"=0.01, "**"=0.05, "*"=0.10)}.
#' @param label Character scalar. The label of the Latex table.
#' @param aic Logical, default is \code{FALSE}. Should the AIC be displayed?
#' @param sqCor Logical, default is \code{FALSE}. Should the squared correlation be displayed?
#' @param subtitles Character vector of the same lenght as the number of models to be displayed. If provided, subtitles are added underneath the dependent variable name.
#' @param showClusterSize Logical, default is \code{FALSE}. If \code{TRUE} and clusters were used in the models, then the number "individuals" of per cluster is also displayed.
#' @param keepFactors Logical, default is \code{FALSE}. By default, when factor variables are contained in the estimation, they are printed as if they were a cluster variable. Put to \code{TRUE} to display all the coefficients of the factor variables.
#' @param bic Logical, default is \code{TRUE}.Should the BIC be reported?
#' @param loglik Logical, default is \code{TRUE}. Should the log-likelihood be reported?
#' @param yesNoCluster A character vector of lenght 2. Default is \code{c("Yes", "No")}. This is the message displayed when a given cluster is (or is not) included in a regression.
#' @param family A logical, default is missing. Whether to display the families of the models. By default this line is displayed when at least two models are from different families.
#' @param powerBelow Integer, default is -5. A coefficient whose value is below \code{10**(powerBelow+1)} is written with a power in Latex. For example \code{0.0000456} would be written \code{4.56$\\times 10^{-5}$} by default. Setting \code{powerBelow = -6} would lead to \code{0.00004} in Latex.
#'
#' @return
#' There is nothing returned, the result is only displayed on the console or saved in a file.
#'
#' @seealso
#' See also the main estimation function \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#'# two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#' res2 = femlm(Sepal.Length ~ Petal.Width | Species, iris)
#'
#' # We export the three results in one Latex table,
#' # with clustered standard-errors:
#' res2tex(res1, res2, se = "cluster")
#'
#' # Changing the names & significance codes
#' res2tex(res1, res2, dict = c(Sepal.Length = "The sepal length", Sepal.Width = "SW"),
#'         signifCode = c("**" = 0.1, "*" = 0.2, "n.s."=1))
#'
res2tex <- function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), cluster, digits=4, pseudo=TRUE, title, sdBelow=TRUE, drop, order, dict, file, replace=FALSE, convergence, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), label, aic=FALSE, sqCor = FALSE, subtitles, showClusterSize = FALSE, bic = TRUE, loglik = TRUE, yesNoCluster = c("Yes", "No"), keepFactors = FALSE, family, powerBelow = -5){
	# drop: a vector of regular expressions
	# order: a vector of regular expressions
	# dict: a 'named' vector
	# file: a character string

	if(missing(se)){
		useSummary = FALSE
	} else {
		useSummary = TRUE
		sdType = match.arg(se)
	}

	# to get the model names
	dots_call = match.call(expand.dots = FALSE)[["..."]]

	if("append" %in% names(dots_call)){
		value = dots_call$append
		replace = !value
		warning("Argument 'append' is deprecated, use argument 'replace' instead.")
		dots_call$append = NULL
	}

	info = results2formattedList(..., se=se, cluster=cluster, digits=digits, sdBelow=sdBelow, signifCode=signifCode, subtitles=subtitles, yesNoCluster=yesNoCluster, keepFactors=keepFactors, isTex = TRUE, useSummary=useSummary, sdType=sdType, dots_call=dots_call, powerBelow=powerBelow)

	n_models = length(info$depvar_list)
	# Getting the information
	se_type_list = info$se_type_list
	var_list = info$var_list
	coef_list = info$coef_list
	coef_below = info$coef_below
	sd_below = info$sd_below
	depvar_list = info$depvar_list
	obs_list = info$obs_list
	r2_list = info$r2_list
	aic_list = info$aic_list
	bic_list = info$bic_list
	loglik_list = info$loglik_list
	convergence_list = info$convergence_list
	sqCor_list = info$sqCor_list
	factorNames = info$factorNames
	isFactor = info$isFactor
	nbFactor = info$nbFactor
	family_list = info$family_list
	theta_list = info$theta_list

	if(!missing(subtitles)){
		isSubtitles = TRUE
	} else {
		isSubtitles = FALSE
	}

	#
	# prompting the infos gathered
	#

	# Starting the table
	myTitle = ifelse(!missing(title), title, "no title")
	if(!missing(label)) myTitle = paste0("\\label{", label, "} ", myTitle)
	start_table = paste0("\\begin{table}[htbp]\\centering\n\\caption{",  .cleanPCT(myTitle), "}\n")
	end_table = "\\end{table}"

	# intro and outro Latex tabular
	myAmpLine = paste0(paste0(rep(" ", length(depvar_list)+1), collapse="&"), "\\tabularnewline\n")
	intro_latex <- paste0("\\begin{tabular}{l", paste0(rep("c", n_models), collapse=""), "}\n",
								 myAmpLine,
								 "\\hline\n",
								 "\\hline\n")

	outro_latex <- "\\end{tabular}\n"

	# 1st lines => dep vars
	# first_line <- paste0("Variables&", paste0(depvar_list, collapse="&"), "\\\\\n\\hline\n\\hline\n")
	depvar_list = c(depvar_list, recursive = TRUE)
	if(!missing(dict)){
		if(!is.character(dict)|| is.null(names(dict))) stop("the arg. 'dict' must be a named character vector.")
		depvar_list = c(depvar_list, recursive = TRUE)
		qui = which(depvar_list%in%names(dict))
		who = depvar_list[qui]
		depvar_list[qui] = dict[who]
	}

	# We write the dependent variables properly, with multicolumn when necessary
	# to do that, we count the number of occurences of each variable (& we respect the order provided by the user)
	nb_multi = 1
	names_multi = depvar_list[1]

	if(n_models > 1){
		k = 1
		old_dep = depvar_list[1]
		for(i in 2:length(depvar_list)){
			if(depvar_list[i] == old_dep){
				nb_multi[k] = nb_multi[k] + 1
			} else {
				k = k + 1
				nb_multi[k] = 1
				names_multi[k] = old_dep = depvar_list[i]
			}
		}
	}

	# now the proper format
	first_line <- "Dependent Variables:"
	if(length(nb_multi) == 1) first_line = "Dependent Variable:"
	for(i in 1:length(nb_multi)){
		if(nb_multi[i] == 1){
			# no multi column
			first_line = paste0(first_line, "&", names_multi[i])
		} else {
			first_line = paste0(first_line, "&\\multicolumn{", nb_multi[i], "}{c}{", names_multi[i], "}")
		}
	}
	first_line = paste0(first_line, "\\\\\n")

	# Model line
	model_line = paste0("Model:&", paste0("(", 1:n_models, ")", collapse = "&"), "\\\\\n")

	# a simple line with only "variables" written in the first cell
	variable_line = "\\hline\n\\emph{Variables}\\tabularnewline\n"


	# Coefficients,  the tricky part
	coef_lines <- list()
	all_vars <- unique(c(var_list, recursive=TRUE))

	# dropping some coefs
	if(!missing(drop)){
		if(!is.character(drop)) stop("the arg. 'drop' must be a character vector of regular expression (see help regex).")
		for(var2drop in drop) all_vars = all_vars[!grepl(var2drop, all_vars)]
	}

	# ordering the coefs
	if(!missing(order)){
		if(!is.character(order)) stop("the arg. 'order' must be a character vector of regular expression (see help regex).")
		for(var2order in rev(order)){
			who = grepl(var2order, all_vars)
			all_vars = c(all_vars[who], all_vars[!who])
		}
	}

	# changing the names of the coefs
	aliasVars = all_vars
	names(aliasVars) = all_vars
	if(!missing(dict)){

		if(!is.character(dict)|| is.null(names(dict))){
			stop("the arg. 'dict' must be a named character vector.")
		}

		qui = all_vars %in% names(dict)
		who = aliasVars[qui]
		aliasVars[qui] = .cleanPCT(dict[who])
	}


	coef_mat <- all_vars
	for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
	coef_mat[is.na(coef_mat)] <- "  "
	if(sdBelow){
		coef_lines = c()
		for(v in all_vars){
			myCoef = mySd= myLine = c()
			for(m in 1:n_models){
				myCoef = c(myCoef, coef_below[[m]][v])
				mySd = c(mySd, sd_below[[m]][v])
			}

			myCoef[is.na(myCoef)] = "  "
			mySd[is.na(mySd)] = "  "
			myCoef = paste0(aliasVars[v], "&", paste0(myCoef, collapse="&"))
			mySd = paste0("  &", paste0(mySd, collapse="&"))
			myLines = paste0(myCoef, "\\\\\n", mySd, "\\\\\n")
			coef_lines = c(coef_lines, myLines)
		}
		coef_lines = paste0(coef_lines, collapse="")
	} else {
		coef_lines = paste0(paste0(apply(coef_mat, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")
	}

	# Factors (if needed)
	if(length(factorNames)>0){
		dumIntro = paste0("\\hline\n\\emph{Fixed-Effects}& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")

		for(m in 1:n_models) {
			quoi = isFactor[[m]][factorNames]
			quoi[is.na(quoi)] = yesNoCluster[2]
			isFactor[[m]] = quoi

			# We do the same for the number of items
			quoi = nbFactor[[m]][factorNames]
			quoi[is.na(quoi)] = "--"
			nbFactor[[m]] = quoi
		}

		allFactors = matrix(c(isFactor, recursive=TRUE), nrow = length(factorNames))
		# We change the names of the factors
		if(!missing(dict)){
			qui = which(factorNames %in% names(dict))
			factorNames[qui] = dict[factorNames[qui]]
		}

		allFactors = cbind(factorNames, allFactors)
		factor_lines <- paste0(paste0(apply(allFactors, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

		# For the number of items
		all_nb_Factors = matrix(c(nbFactor, recursive=TRUE), nrow = length(factorNames))
		factorNames_nbItems = paste0("# ", factorNames)
		all_nb_Factors = cbind(factorNames_nbItems, all_nb_Factors)
		nb_factor_lines <- paste0(paste0(apply(all_nb_Factors, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")


	} else {
		factor_lines = NULL
		dumIntro = NULL
	}

	# Subtitles
	if(isSubtitles){
		info_subtitles = paste0("  & ", paste(subtitles, collapse="&"), "\\\\\n")
	} else {
		info_subtitles = ""
	}

	# Fit statistics
	fit_part <- paste0("\\hline\n\\emph{Fit statistics}& ", paste(rep(" ", n_models), collapse="&"), "\\\\\n")
	# Misc
	info_aic <- paste0("AIC & ", paste(numberFormatLatex(aic_list), collapse="&"), "\\\\\n")
	info_loglik <- paste0("Log-Likelihood & ", paste(numberFormatLatex(loglik_list), collapse="&"), "\\\\\n")
	info_bic <- paste0("BIC & ", paste(numberFormatLatex(bic_list), collapse="&"), "\\\\\n")

	info_obs <- paste0("Observations& ", paste(addCommas(obs_list), collapse="&"), "\\\\\n")
	info_r2 <- paste0("Adj-pseudo $R^2$ &", paste(r2_list, collapse="&"), "\\\\\n")
	info_sqCor <- paste0("$R^2$ &", paste(sqCor_list, collapse="&"), "\\\\\n")

	# Convergence information
	info_convergence = ""
	if((missing(convergence) && any(convergence_list == FALSE)) || (!missing(convergence) && convergence)){
		info_convergence = paste0("Convergence &", paste(convergence_list, collapse="&"), "\\\\\n")
	}

	info_theta <- paste0("Overdispersion& ", paste(theta_list, collapse="&"), "\\\\\n")

	# information on family
	if((!missing(family) && family) || (missing(family) && length(unique(family_list)) > 1)){
		info_family <- paste0("Family& ", paste(family_list, collapse="&"), "\\\\\n")
	} else {
		info_family = ""
	}


	# The standard errors
	isUniqueSD = length(unique(unlist(se_type_list))) == 1
	if(isUniqueSD){
		my_se = unique(unlist(se_type_list)) # it comes from summary
		# every model has the same type of SE
		if(my_se == "Standard") my_se = "Normal"
		if(my_se == "White") my_se = "White-corrected"

		# Now we modify the names of the clusters if needed
		if(!missing(dict) && grepl("\\(", my_se)){
			# we extract the clusters
			se_cluster = strsplit(gsub("(^.+\\()|(\\))", "", my_se), " & ")[[1]]
			qui = se_cluster %in% names(dict)
			se_cluster[qui] = dict[se_cluster[qui]]
			new_se = gsub("\\(.+", "", my_se)
			my_se = paste0(new_se, "(", paste0(se_cluster, collapse = " & "), ")")
		}

		nb_col = length(obs_list)+1
		info_SD = paste0("\\hline\n\\hline\n\\multicolumn{", nb_col, "}{l}{\\emph{", my_se, " standard-errors in parenthesis. Signif Codes: ", paste(names(signifCode), signifCode, sep=": ", collapse = ", "), "}}\\\\\n")
		info_muli_se = ""
	} else {
		info_muli_se = paste0("Standard-Error type& ", paste(se_type_list, collapse="&"), "\\\\\n")
		info_SD = "\\hline\n\\hline\n\\\\\n"
	}

	# Information on number of items

	supplemental_info = "\\global\\long\\def\\sym#1{\\ifmmode^{#1}\\else\\ensuremath{^{#1}}\\fi}\n"

	if(!pseudo) info_r2 <- ""
	if(!sqCor) info_sqCor <- ""
	if(!aic) info_aic = ""
	if(!bic) info_bic = ""
	if(!loglik) info_loglik = ""
	if(!showClusterSize) nb_factor_lines = ""
	if(all(theta_list == "")) info_theta = ""

	if(!missing(file)) sink(file = file, append = !replace)

	cat(paste0(supplemental_info,
				  start_table,
				  intro_latex,
				  first_line,
				  info_subtitles,
				  model_line,
				  info_family,
				  variable_line,
				  coef_lines,
				  info_theta,
				  dumIntro,
				  factor_lines,
				  fit_part,
				  info_obs,
				  nb_factor_lines,
				  info_convergence,
				  info_muli_se,
				  info_r2,
				  info_sqCor,
				  info_aic,
				  info_loglik,
				  info_bic,
				  info_SD,
				  outro_latex,
				  end_table))

	if(!missing(file)) sink()

}

#' Facility to display the results of multiple \code{femlm} estimations.
#'
#' This function aggregates the results of multiple estimations and display them in the form of only one table whose rownames are the variables and the columns contain the coefficients and standard-errors.
#'
#' @inheritParams res2tex
#' @inheritParams summary.femlm
#'
#' @param depvar Logical, default is missing. Whether a first line containing the dependent variables should be shown. By default, the dependent variables are shown only if they differ across models.
#'
#' @return
#' Returns a data.frame containing the formatted results.
#'
#' @seealso
#' See also the main estimation function \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#' # estimation without clusters
#' res2 = update(res1, . ~ Sepal.Width | 0)
#'
#' # We export the two results in one Latex table:
#' res2table(res1, res2)
#'
#' # With clustered standard-errors + showing the dependent variable
#' res2table(res1, res2, se = "cluster", cluster = iris$Species, depvar = TRUE)
#'
#' # Changing the model names + the order of the variables
#' # + dropping the intercept.
#' res2table(model_1 = res1, res2,
#'           order = c("Width", "Petal"), drop = "Int",
#'           signifCode = c("**" = 0, "*" = 0.2, "n.s."=1))
#'
#'
#'
res2table <- function(..., se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), cluster, depvar, drop, order, digits=4, pseudo=TRUE, convergence, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), subtitles, keepFactors = FALSE, family){

	# Need to check for the presence of the se
	if(missing(se)){
		useSummary = FALSE
	} else {
		useSummary = TRUE
		sdType = match.arg(se)
	}

	# to get the model names
	dots_call = match.call(expand.dots = FALSE)[["..."]]

	info = results2formattedList(..., se=se, cluster=cluster, digits=digits, signifCode=signifCode, subtitles=subtitles, keepFactors=keepFactors, useSummary=useSummary, sdType=sdType, dots_call=dots_call)

	n_models = length(info$depvar_list)
	# Getting the information
	se_type_list = info$se_type_list
	var_list = info$var_list
	coef_list = info$coef_list
	coef_below = info$coef_below
	sd_below = info$sd_below
	depvar_list = info$depvar_list
	obs_list = info$obs_list
	r2_list = info$r2_list
	aic_list = info$aic_list
	bic_list = info$bic_list
	loglik_list = info$loglik_list
	convergence_list = info$convergence_list
	sqCor_list = info$sqCor_list
	factorNames = info$factorNames
	isFactor = info$isFactor
	nbFactor = info$nbFactor
	useSummary = info$useSummary
	depvar_list = info$depvar_list
	model_names = info$model_names
	family_list = info$family_list
	theta_list = info$theta_list


	if(!missing(subtitles)){
		isSubtitles = TRUE
	} else {
		isSubtitles = FALSE
	}

	# The coefficients

	all_vars <- unique(c(var_list, recursive=TRUE))

	# dropping some coefs
	if(!missing(drop)){
		if(!is.character(drop)) stop("the arg. 'drop' must be a character vector of regular expression (see help regex).")
		for(var2drop in drop) all_vars = all_vars[!grepl(var2drop, all_vars)]
	}

	# ordering the coefs
	if(!missing(order)){
		if(!is.character(order)) stop("the arg. 'order' must be a character vector of regular expression (see help regex).")
		for(var2order in rev(order)){
			who = grepl(var2order, all_vars)
			all_vars = c(all_vars[who], all_vars[!who])
		}
	}

	coef_mat <- all_vars
	for(m in 1:n_models) coef_mat <- cbind(coef_mat, coef_list[[m]][all_vars])
	coef_mat[is.na(coef_mat)] <- "  "
	res = coef_mat

	if("Neg. Bin." %in% family_list){
		theta_line = c("Overdispersion:", unlist(theta_list))
		res = rbind(res, theta_line)
	}

	# Used to draw a line
	myLine = "______________________________________"
	longueur = apply(res, 2, function(x) max(nchar(as.character(x))))
	theLine = sapply(longueur, function(x) sprintf("%.*s", x, myLine))
	theLine[1] = sprintf("%.*s", max(nchar(theLine[1]), 19), myLine)

	if(length(factorNames)>0){

		for(m in 1:n_models) {
			quoi = isFactor[[m]][factorNames]
			quoi[is.na(quoi)] = "No"
			isFactor[[m]] = quoi
		}
		allFactors = matrix(c(isFactor, recursive=TRUE), nrow = length(factorNames))
		allFactors = cbind(factorNames, allFactors)
		factor_lines <- paste0(paste0(apply(allFactors, 1, paste0, collapse="&"), collapse="\\\\\n"), "\\\\\n")

		myLine = "-------------------------------"

		res = rbind(res, c("Fixed-Effects:", sprintf("%.*s", longueur[-1], myLine)))
		factmat = matrix(c(strsplit(strsplit(factor_lines, "\n")[[1]], "&"), recursive = TRUE), ncol=n_models+1, byrow=TRUE)
		factmat[, ncol(factmat)]=gsub("\\", "", factmat[, ncol(factmat)], fixed = TRUE)
		res = rbind(res, factmat)
	}

	# The line with the dependent variable
	preamble = c()
	if((missing(depvar) && length(unique(unlist(depvar_list))) > 1) || (!missing(depvar) && depvar)){
		preamble = rbind(c("Dependent Var.:", unlist(depvar_list)), preamble)
	}



	if(length(preamble) > 0){
		# preamble = rbind(preamble, c("  ", theLine[-1]))
		preamble = rbind(preamble, rep("   ", length(theLine)))
		res <- rbind(preamble, res)
	}

	res <- rbind(res, theLine)

	# the line with the families
	if((missing(family) && length(unique(unlist(family_list))) > 1) || (!missing(family) && family)){
		# preamble = rbind(c("Family:", unlist(family_list)), preamble)
		res = rbind(res, c("Family:", unlist(family_list)))
	}

	res <- rbind(res, c("Observations", addCommas(obs_list)))
	if(!useSummary){
		se_type_format = c()
		for(m in 1:n_models) se_type_format[m] = charShorten(se_type_list[[m]], longueur[[1+m]])
		res <- rbind(res, c("S.E. type", c(se_type_format, recursive = TRUE)))
	}

	# convergence status
	if((missing(convergence) && any(convergence_list == FALSE)) || (!missing(convergence) && convergence)){
		res <- rbind(res, c("Convergence", c(convergence_list, recursive = TRUE)))
	}

	res <- rbind(res, c("Squared-Correlation", c(sqCor_list, recursive = TRUE)))
	res <- rbind(res, c("Adj-pseudo R2", c(r2_list, recursive = TRUE)))
	# res <- rbind(res, c("AIC", c(aic_list, recursive = TRUE)))
	res <- rbind(res, c("Log-Likelihood", numberFormatNormal(loglik_list)))
	res <- rbind(res, c("BIC", numberFormatNormal(bic_list)))

	# if subtitles
	if(isSubtitles){
		modelNames = subtitles
	} else {
		# modelNames = paste0("model ", 1:n_models)
		modelNames = model_names
	}

	# we shorten the model names to fit the width
	for(m in 1:n_models) modelNames[m] = charShorten(modelNames[[m]], longueur[[1+m]])

	res <- as.data.frame(res)
	names(res) <- c("variables", modelNames)
	row.names(res) = res$variables
	res$variables = NULL

	# We rename theta when NB is used
	quiTheta = which(row.names(res) == ".theta")
	row.names(res)[quiTheta] = "Dispersion Parameter"

	return(res)
}

results2formattedList = function(..., se=c("standard", "white", "cluster", "twoway"), cluster, digits=4, pseudo=TRUE, sdBelow=TRUE, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10), label, subtitles, yesNoCluster = c("Yes", "No"), keepFactors = FALSE, isTex = FALSE, useSummary, sdType, dots_call, powerBelow){
	# This function is the core of the functions res2table and res2tex

	signifCode = sort(signifCode)
	if(any(signifCode<0) | any(signifCode>1)) stop("The argument 'signifCode' must lie between 0 and 1.")

	if(length(yesNoCluster) != 2) stop("The argument 'yesNoCluster' must be of length 2.")

	# To take care of old verions:
	allowed_types = c("femlm", "feNmlm")

	# We get all the models
	dots <- list(...)

	# for retro-compatibility:
	if("sd" %in% names(dots)){
		warning("The use of the argument 'sd' is deprecated, it is now replaced by the argument 'se'.")
		se = dots$sd
	}

	# formatting the names of the models
	dots_names = names(dots_call)
	if(!is.null(dots_names)){

		for(i in 1:length(dots_call)){
			if(dots_names[i] != ""){
				dots_call[[i]] = dots_names[i]
			} else {
				dots_call[[i]] = deparse(dots_call[[i]])
			}
		}
	}

	n = length(dots)
	all_models = list()
	model_names = list()
	k = 1
	for(i in 1:n){
		di = dots[[i]]

		if(any(allowed_types %in% class(di))){
			all_models[[k]] = di
			if(any(class(dots_call[[i]]) %in% c("call", "name"))){
				model_names[[k]] = deparse(dots_call[[i]])
			} else {
				model_names[[k]] = as.character(dots_call[[i]])
			}

			k = k+1
		} else if(length(class(di))==1 && class(di)=="list"){
			# we get into this list to get the FENmlm objects
			types = sapply(di, class)
			qui = which(types %in% allowed_types)
			for(m in qui){
				all_models[[k]] = di[[m]]

				# handling names
				if(n > 1){
					if(is.null(names(di)[m]) || names(di)[m]==""){
						model_names[[k]] = paste0(dots_call[[i]], "[[", m, "]]")
					} else {
						model_names[[k]] = paste0(dots_call[[i]], "$", names(di)[m])
					}
				} else {
					model_names[[k]] = as.character(names(di)[m])
				}

				k = k+1
			}
		}

	}

	if(length(all_models)==0) stop("Not any proper model (femlm) as argument!")

	n_models <- length(all_models)

	# formatting the names (if needed)
	alternative_names = paste0("model ", 1:n_models)
	who2replace = sapply(model_names, function(x) length(x) == 0 || x == "")
	model_names[who2replace] = alternative_names[who2replace]

	# we keep track of the SEs
	se_type_list = list()

	var_list <- coef_list <- coef_below <- sd_below <- list()
	depvar_list <- obs_list <- list()
	r2_list <- aic_list <- bic_list <- loglik_list <- convergence_list <- list()
	sqCor_list = family_list = theta_list = list()

	# To take care of factors
	factorNames = c()
	isFactor = vector(mode = "list", n_models)
	nbFactor = vector(mode = "list", n_models) # the number of items per factor

	# if there are subtitles
	if(!missing(subtitles)){
		if(length(subtitles) != n_models){
			stop("If argument 'subtitles' is provided, it must be of the same length as the number of models.")
		} else {
			isSubtitles = TRUE
		}
	} else {
		isSubtitles = FALSE
	}

	for(m in 1:n_models){

		# If se is provided, we use summary
		if(useSummary){
			x = summary(all_models[[m]], se=sdType, cluster)
		} else {
			x = all_models[[m]]
		}
		se_type_list[[m]] = attr(x$se, "type")

		# family
		family = all_models[[m]]$family
		family_list[[m]] = switch(family, poisson = "Poisson", negbin = "Neg. Bin.", gaussian = "Gaussian", logit = "Logit")

		# Negbin parameter
		theta = all_models[[m]]$theta
		theta_list[[m]] = ifelse(is.null(theta), "", numberFormatNormal(theta))

		# variable dependante:
		depvar <- gsub(" ", "", as.character(x$fml)[[2]])

		a <- x$coeftable
		if(!is.data.frame(a)){
			class(a) <- NULL
			a = as.data.frame(a)
		}

		# We drop the .theta coefficient
		if(family == "negbin"){
			quiTheta = rownames(a) == ".theta"
			a = a[!quiTheta, ]
		}

		#
		# Formatting of the factors
		#

		# on enleve les facteurs des variables a garder
		if(!keepFactors){
			fact = rownames(a)
			qui_drop = grepl("factor(", fact, fixed = TRUE)
			a = a[!qui_drop, , FALSE]
			b = fact[qui_drop]
			c = sapply(b, function(x) strsplit(x, "factor(", fixed=TRUE)[[1]][2])
			d = sapply(c, function(x) strsplit(x, ")", fixed=TRUE)[[1]][1])
			factor_var = unique(d)

			# Now the number of items per factor
			if(length(factor_var) == 0){
				nbItems = character(0)
			} else {
				nbItems = addCommas(sapply(factor_var, function(x) 1+sum(grepl(x, b))))
			}
		} else {
			factor_var = c()
			nbItems = character(0)
		}

		# now the normal clusters
		if(!is.null(x$clusterNames)){
			factor_var = c(factor_var, x$clusterNames, recursive=TRUE)

			if(!is.null(x$clusterSize)){
				new_items = addCommas(as.vector(x$clusterSize))
				names(new_items) = names(x$clusterSize)
			} else {
				# for retro compatibility
				new_items = addCommas(sapply(x$id_dummies, function(y) length(unique(y))))
			}

			nbItems = c(nbItems, new_items)
		}

		nbFactor[[m]] = nbItems

		# Formatting

		lFactor = rep(yesNoCluster[1], length(factor_var))
		names(lFactor) = factor_var
		isFactor[[m]] = lFactor

		factorNames = unique(c(factorNames, factor_var, recursive=TRUE))

		#
		# END: cluster formatting
		#

		# on enleve les espaces dans les noms de variables
		var <- c(gsub(" ", "", row.names(a)))

		if(isTex){
			coef = coefFormatLatex(a[, 1], digits = digits, power = abs(powerBelow))
			se = coefFormatLatex(a[, 2], digits = digits, power = abs(powerBelow))
		} else {
			coef = as.character(round(a[, 1], digits))
			se = as.character(myRound(a[, 2], digits))
		}

		if(isTex){
			pval = cut(a[, 4], breaks = c(-1, signifCode, 100), labels = c(paste0("\\sym{",names(signifCode),"}"), ""))
		} else {
			pval = cut(a[, 4], breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
		}

		# If the coefficient is bounded, we supress the 'stars'
		isBounded = grepl("bounded", se)
		if(any(isBounded)){
			pval[isBounded] = ""
		}

		structured_coef = c(paste0(coef, pval, " (", se, ")"))

		# saving the infos
		var_list[[m]] <- var
		names(structured_coef) <- var
		coef_list[[m]] <- structured_coef
		if(sdBelow){
			cb = c(paste0(coef, pval))
			sb = c(paste0("(", se, ")"))
			names(cb) = names(sb) = var
			coef_below[[m]] = cb
			sd_below[[m]] = sb
		}

		# La depvar
		depvar_list[[m]] <- depvar

		# statistics
		# Pseudo-R2 // AIC // BIC // N
		n <- x$n
		obs_list[[m]] <- n
		convergence_list[[m]] = x$convStatus

		K <- x$nparams
		ll <- x$loglik
		bic_list[[m]] <- round(-2*ll+K*log(n), 3)
		aic_list[[m]] <- round(-2*ll+2*K, 3)
		loglik_list[[m]] <- round(ll, 3)
		r2_list[[m]] <- round(x$pseudo_r2, 5)
		sqCor_list[[m]] <- round(x$sq.cor, 3)

	}

	res = list(se_type_list=se_type_list, var_list=var_list, coef_list=coef_list, coef_below=coef_below, sd_below=sd_below, depvar_list=depvar_list, obs_list=obs_list, r2_list=r2_list, aic_list=aic_list, bic_list=bic_list, loglik_list=loglik_list, convergence_list=convergence_list, sqCor_list=sqCor_list, factorNames=factorNames, isFactor=isFactor, nbFactor=nbFactor, useSummary=useSummary, model_names=model_names, family_list=family_list, theta_list=theta_list)

	return(res)
}

.cleanPCT = function(x){
	# changes % into \% => to escape that character in Latex
	res = gsub("%", "\\%", x, fixed = TRUE)
	res = gsub("_", "\\_", res, fixed = TRUE)
	res
}

myPrintCoefTable = function(coeftable){
	# Simple function that does as the function coeftable but handles special cases
	# => to take care of the case when the coefficient is bounded

	if(!is.data.frame(coeftable)){
		class(coeftable) = NULL
		ct = as.data.frame(coeftable)
	} else {
		ct = coeftable
	}

	signifCode = c("***"=0.001, "** "=0.01, "*  "=0.05, ".  "=0.1)

	pvalues = ct[, 4]

	stars = cut(pvalues, breaks = c(-1, signifCode, 100), labels = c(names(signifCode), ""))
	stars[is.na(stars)] = ""

	whoIsLow = !is.na(pvalues) & pvalues < 2.2e-16

	for(i in 1:4){
		ct[, i] = decimalFormat(ct[, i])
	}

	# browser()

	ct[whoIsLow, 4] = "< 2.2e-16"
	ct[is.na(ct[, 4]), 4] = "NA"

	ct[, 5] = stars
	names(ct)[5] = ""

	print(ct)

	cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

}


#' Summary method for cluster coefficients
#'
#' This function summarizes the main characteristics of the cluster coefficients. It shows the number of fixed-effects that have been set as references and the first elements of the fixed-effects.
#'
#' @method summary femlm.allClusters
#'
#' @param object An object returned by the function \code{\link[FENmlm]{getFE}}.
#' @param n Positive integer, defaults to 5. The \code{n} first fixed-effects for each cluster are reported.
#' @param ... Not currently used.
#'
#' @return
#' It prints the number of fixed-effect coefficients per cluster, as well as the number of fixed-effects used as references for each cluster, and the mean and variance of the cluster coefficients. Finally it reports the first 5 elements of each cluster.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{getFE}}, \code{\link[FENmlm]{plot.femlm.allClusters}}.
#'
#' @examples
#'
#' data(trade)
#'
#' # We estimate the effect of distance on trade
#' # => we account for 3 cluster effects
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # obtaining the cluster coefficients
#' fe_trade = getFE(est_pois)
#'
#' # printing some summary information on the cluster coefficients:
#' fe_trade
#'
#'
summary.femlm.allClusters = function(object, n=5, ...){
	# This function shows some generic information on the clusters

	Q = length(object)
	clustNames = names(object)

	isRegular = TRUE
	if(Q > 1){
		nb_ref = attr(object, "References")
		nb_per_cluster = sapply(object, length)
		mean_per_cluster = var_per_cluster = c()
		for(i in 1:Q){
			mean_per_cluster[i] = as.character(signif(mean(object[[i]]), 3))
			var_per_cluster[i] = as.character(signif(var(object[[i]]), 3))
		}
		res = as.data.frame(rbind(nb_per_cluster, nb_ref, mean_per_cluster, var_per_cluster))
		rownames(res) = c("Number of fixed-effects", "Number of references", "Mean", "Variance")

		if(sum(nb_ref) > Q-1){
			isRegular = FALSE
		}
	}

	# The message
	cat("Fixed-effects coefficients.\n")
	if(Q == 1){
		x1 = object[[1]]
		cat("Number of fixed-effects for variable ", clustNames, " is ", length(x1), ".\n", sep = "")
		cat("\tMean = ", signif(mean(x1), 3), "\tVariance = ", signif(var(x1), 3), "\n", sep = "")
	} else {
		print(res)
		if(!isRegular){
			cat("NOTE: The fixed-effects are NOT regular, so cannot be straighforwardly interpreted.\n")
		}
	}

	# We print the first 5 elements of each cluster
	cat("\nCOEFFICIENTS:\n")
	for(i in 1:Q){
		m = head(object[[i]], n)

		m_char = as.data.frame(t(as.data.frame(c("", as.character(signif(m, 4))))))
		names(m_char) = c(paste0(clustNames[i], ":"), names(m))
		rownames(m_char) = " "

		n_cluster = length(object[[i]])
		if(n_cluster > n){
			m_char[["   "]] = paste0("... ", addCommas(n_cluster - n), " remaining")
		}

		print(m_char)
		if(i != Q) cat("-----\n")
	}

}


#' Extract the Fixed-Effects from a \code{femlm} estimation.
#'
#' This function retrives the fixed effects from a femlm estimation. It is useful only when there are more than one cluster.
#'
#' @param x A \code{\link[FENmlm]{femlm}} object.
#'
#' If the cluster coefficients not regular, then several reference points need to be set, leading to the coefficients to be NOT interpretable. If this is the case, then a warning is raised.
#'
#' @return
#' A list containig the vectors of the fixed effects.
#'
#' If there is more than 1 cluster, then the attribute \dQuote{References} is created. This is a vector of length the number of clusters, each element contains the number of fixed-effects set as references. By construction, the elements of the first clusters are never set as references. In the presence of regular clusters, there should be Q-1 references (with Q the number of clusters).
#'
#' @seealso
#' \code{\link[FENmlm]{plot.femlm.allClusters}}. See also the main estimation function \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#'
#' @examples
#'
#' data(trade)
#'
#' # We estimate the effect of distance on trade => we account for 3 cluster effects
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # obtaining the cluster coefficients
#' fe_trade = getFE(est_pois)
#'
#' # plotting them
#' plot(fe_trade)
#'
#' # plotting only the Products fixed-effects & showing more of them
#' plot(fe_trade$Product, n=8)
#'
getFE = function(x){
	# x is a femlm object
	# This function retrieves the dummies

	if(class(x) != "femlm"){
		stop("Argument 'x' mus be a femlm object.")
	}

	# Preliminary stuff
	if("dummies" %in% names(x)){
		# for retro-compatibility
		S = x$dummies
	} else {
		S = x$sumFE
	}

	if(is.null(S)){
		stop("There is no cluster to be retrieved.")
	}

	family = x$family
	clustNames = x$clusterNames

	id_dummies = x$id_dummies

	Q = length(id_dummies)
	N = length(S)

	# either (we need to clean its attributes for unlist to be efficient)
	id_dummies_vect = list()
	for(i in 1:Q) id_dummies_vect[[i]] = as.vector(id_dummies[[i]])

	if(Q == 1){
		# This is the simplest case
		id = id_dummies_vect[[1]]

		myOrder = order(id)
		myDiff = c(1, diff(id[myOrder]))

		select = myOrder[myDiff == 1]

		cluster_values = list(S[select])

		# There are no references => no need to set nb_ref
	} else if(Q == 2) {
		# specific method that is faster but specific to the case of 2 FE
		dum_i = id_dummies_vect[[1]] - 1L
		dum_j = id_dummies_vect[[2]] - 1L

		order_ij = order(dum_i, dum_j)
		i_sorted_index_j = dum_j[order_ij]

		order_ji = order(dum_j, dum_i)
		j_sorted_index_i = dum_i[order_ji]

		i_sorted_sumFE = S[order_ij]
		j_sorted_sumFE = S[order_ji]

		fe <- cpp_get_fe_2(clusterSize = x$clusterSize, i_sorted_index_j = i_sorted_index_j, i_sorted_sumFE = i_sorted_sumFE, j_sorted_index_i = j_sorted_index_i, j_sorted_sumFE = j_sorted_sumFE, r_cumtable_i = cumsum(table(dum_i)), r_cumtable_j = cumsum(table(dum_j)))

		# cpp_get_fe_2 returns a matrix (we will update cpp_get_fe_gnl to return a matrix later)
		# so we put it into a list (for now)

		cluster_values = list()
		cluster_values[[1]] = fe[fe[, 1] == 1, 3]
		cluster_values[[2]] = fe[fe[, 1] == 2, 3]

		# the number of references needed
		nb_ref = c(0, sum(fe[, 4]))

	} else {
		# We apply a Rcpp script to handle complicated cases (and we don't know beforehand if the input is one)

		dumMat <- matrix(unlist(id_dummies_vect), N, Q) - 1
		orderCluster <- matrix(unlist(lapply(id_dummies_vect, order)), N, Q) - 1

		nbCluster = sapply(id_dummies, max)

		cluster_values <- cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)

		# the information on the references
		nb_ref = cluster_values[[Q+1]]
		cluster_values[[Q+1]] = NULL
	}

	# now saving & adding information
	all_clust = list()
	for(i in 1:Q){
		cv = cluster_values[[i]]
		names(cv) = attr(id_dummies[[i]], "clust_names")
		all_clust[[clustNames[i]]] = cv
	}

	class(all_clust) = c("femlm.allClusters", "list")

	# Dealing with the references
	if(Q > 1){
		names(nb_ref) = clustNames
		attr(all_clust, "References") = nb_ref

		# warning if unbalanced
		if(sum(nb_ref) > Q-1){
			warning("The fixed-effects are not regular, they cannot be straightforwardly interpreted.", call. = FALSE)
		}
	}

	return(all_clust)
}


plot_single_cluster = function(x, n=5, ...){
	# It plots the n first and last most notable FEs

	# we compare with the average of the coefficients
	x = sort(x) - mean(x)
	x_name = names(x)

	k = length(x)
	nb_show = min(k, 2*n+1)
	mid_value = floor(nb_show/2)
	xlim = c(1, nb_show)
	xlim = xlim + c(-1, 1) * diff(xlim)/30

	if(k <= 2*n+1){
		# no need of space to split the data
		x_values = 1:k
		y_values = x
		isSplit = FALSE
	} else {
		nb_show = nb_show - 1 # because we don't want to show the middle point
		x_values = c(1:n, (n+2):(2*n+1))
		y_values = c(head(x, n), tail(x, n))
		isSplit = TRUE
	}

	# very specific case where the axis confonds with the boxing
	ylim = range(y_values)
	if(pretty(range(y_values))[1] < min(y_values) || tail(pretty(range(y_values)), 1) > max(y_values)){
		ylim = range(y_values) + c(-1, 1) * diff(range(y_values))/30
	}

	plot(x_values, y_values, ann = FALSE, axes = FALSE, xlim=xlim, ylim = ylim, col = 0)

	# display
	box()
	y = axis(2)
	abline(h = y, lty = 4, col = "lightgrey")

	# name information & points
	points(x_values, y_values)
	text(head(x_values, mid_value), head(y_values, mid_value), head(x_name, mid_value), pos = 3)
	text(tail(x_values, nb_show-mid_value), tail(y_values, nb_show-mid_value), tail(x_name, nb_show-mid_value), pos = 1)

	axis(4, at=y, labels = signif(exp(y), 2))

	title(xlab = "Centered Fixed-Effects", line = 1)

	if(isSplit){
		abline(v = c(n+0.75, n+1.25), lty = 2)

		axis(1, at = c(n+0.75, n+1.25), col = "grey", labels = NA, lwd.ticks = 0)
		axis(3, at = c(n+0.75, n+1.25), col = "grey", labels = NA, lwd.ticks = 0)
	}

	coord = par("usr")
	axis(1, at = coord[1:2], labels = c("coef", "exp(coef)"), lwd = 0, lwd.ticks = 1)

}


#' Displaying the most notable fixed-effects
#'
#' This function plots the 5 fixed-effects with the highest and lowest values, for each of the clusters. It takes as an argument the fixed-effects obtained from the function \code{\link{getFE}} after and estimation using \code{\link{femlm}}.
#'
#' @method plot femlm.allClusters
#'
#' @param x An object obtained from the function \code{\link{getFE}}.
#' @param n The number of fixed-effects to be drawn. Defaults to 5.
#' @param ... Not currently used.
#'
#' Note that the fixed-effect coefficients might NOT be interpretable. This function is useful only for fully regular panels.
#'
#' If the data are not regular in the cluster coefficients, this means that several \sQuote{reference points} are set to obtain the fixed-effects, thereby impeding their interpretation. In this case a warning is raised.
#'
#' @seealso
#' \code{\link[FENmlm]{getFE}} to extract clouster coefficients. See also the main estimation function \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' data(trade)
#'
#' # We estimate the effect of distance on trade
#' # => we account for 3 cluster effects
#' est_pois = femlm(Euros ~ log(dist_km)|Origin+Destination+Product, trade)
#'
#' # obtaining the cluster coefficients
#' fe_trade = getFE(est_pois)
#'
#' # plotting them
#' plot(fe_trade)
#'
#'
plot.femlm.allClusters = function(x, n=5, ...){

	Q = length(x)

	mfrow = as.character(c(11, 12, 22, 22, 32, 32, 33, 33))

	if(Q > 1 && sum(attr(x, "References")) > Q-1){
		warning("The fixed-effects are not regular, they cannot be straightforwardly interpreted.", call. = FALSE)
	}

	op = par(mfrow = as.numeric(strsplit(mfrow[Q], "")[[1]]), mar = c(3, 3, 2.5, 3))

	for(i in 1:Q){
		plot_single_cluster(x[[i]], n=n)
		title(main = names(x)[i])
	}

	par(op)
}

#
# To compute clustered standard errors
#

vcovClust <- function (cluster, myBread, scores, dof_correction=FALSE, do.unclass=TRUE){
	# Internal function: no need for controls, they come beforehand
	#	- cluster: the vector of dummies
	#	- myBread: the naive variance covariance matrix
	# - scores
	#Note: if length(unique(cluster)) == n (i.e. White correction), then the dof are such that vcovClust is equivalent to vcovHC(res, type="HC1")

	n <- NROW(scores)
	k <- NCOL(scores)

	# Control for cluster type
	if(do.unclass){
		cluster <- unclass(as.factor(cluster))
	}

	Q <- max(cluster)
	RightScores = cpp_tapply_sum(Q, scores, cluster)

	# Finite sample correction:
	if(dof_correction) dof  <- Q / (Q-1) * (n-1) / (n-k)
	else dof = 1

	return(crossprod(RightScores%*%myBread) * dof)
}


formatBicLL = function(bic, ll){
	# entry: bic and ll

	bic = numberFormatNormal(bic)
	ll = numberFormatNormal(ll)

	bic_split = strsplit(bic, "\\.")[[1]]
	ll_split = strsplit(ll, "\\.")[[1]]

	n = max(nchar(bic_split[1]), nchar(ll_split[1]))

	bic_new = sprintf("% *s.%s", n, bic_split[1], ifelse(length(bic_split) == 2, bic_split[2], 0))
	ll_new = sprintf("% *s.%s", n, ll_split[1], ifelse(length(ll_split) == 2, ll_split[2], 0))

	sep = "  "
	myWidth = max(nchar(c(bic_new, ll_new))) + length(sep) + 1

	bic_format = paste0(bic_new, sprintf("% *s", myWidth - nchar(bic_new), sep))
	ll_format = paste0(ll_new, sprintf("% *s", myWidth - nchar(ll_new), sep))

	list(bic = bic_format, ll = ll_format)
}

addCommas_single = function(x){

	if (!is.finite(x)) return(as.character(x))

	s = sign(x)
	x = abs(x)
	decimal = x - floor(x)
	if (decimal > 0){
		dec_string = substr(decimal, 2, 4)
	} else {
		dec_string = ""
	}

	entier = sprintf("%.0f", floor(x))
	quoi = rev(strsplit(entier, "")[[1]])
	n = length(quoi)
	sol = c()
	for (i in 1:n) {
		sol = c(sol, quoi[i])
		if (i%%3 == 0 && i != n) sol = c(sol, ",")
	}
	res = paste0(ifelse(s == -1, "-", ""), paste0(rev(sol), collapse = ""),
					 dec_string)
	res
}

addCommas = function(x){
	sapply(x, addCommas_single)
}

myRound_single = function(x, digits=5){
	# There can be non numeric values...
	# we give away the non numeric ones and round the others

	if(is.na(x)){
		return(NA)
	}

	if(is.numeric(x)){
		res = round(x, digits)
	} else {

		if(!grepl("[[:digit:]]", x)){
			# means it is a character
			res = x
		} else {
			res = round(as.numeric(x), digits)
		}
	}

	res
}

myRound = function(x, digits=5){
	sapply(x, myRound_single, digits = digits)
}

decimalFormat_single = function(x){
	# for very small numbers: format 5.2e-08

	if(is.na(x) || !is.numeric(x)) return(x)

	xPower = log10(abs(x))

	if(xPower < -5){
		res = signif(x, 3)
	} else if(xPower < 0){
		res = round(x, 6)
	} else {
		res = round(x, max(1, 5 - ceiling(xPower)))
	}

	res
}

decimalFormat = function(x){
	sapply(x, decimalFormat_single)
}

coefFormatLatex_single = function(x, digits, power){
	# format decimals: 5.3 10**-7 instead of 0.00000053
	# format large numbers 6356516.12464 => 6356516.1

	nbSignif = 3

	if(is.na(x)) return(x)

	if(!is.numeric(x)){
		if(grepl("[^[:digit:]e\\.-]", x)){
			return(x)
		} else {
			x = as.numeric(x)
		}
	}

	exponent = floor(log10(abs(x)))

	if(exponent > 0){
		return(sprintf("%.*f", max(1, digits - abs(exponent)), x))
	}

	if(abs(exponent) >= power){
		left_value = round(x*10**-exponent, 3)
		res = paste0("$", left_value, "\\times 10^{", exponent, "}$")
	} else if(x > 10**(-digits)){
		res = sprintf("%.*f", digits, x)
	} else {
		res = sprintf("%.*f", abs(exponent), x)
	}

	res
}

coefFormatLatex = function(x, digits = 4, power = 5){
	sapply(x, coefFormatLatex_single, digits = digits, power = power)
}

numberFormat_single = function(x, type = "normal"){
	# For numbers higher than 1e9 => we apply a specific formatting
	# idem for numbers lower than 1e-4

	if(x == 0) return("0")

	exponent = floor(log10(abs(x)))

	if(-4 < exponent && exponent < 9){
		if(exponent > 0){
			return(addCommas(x))
		} else {
			return(decimalFormat(x))
		}

	}

	left_value = round(x*10**-exponent, 3)

	if(type == "latex"){
		res = paste0("$", left_value, "\\times 10^{", exponent, "}$")
	} else {
		res = paste0(left_value, "e", ifelse(exponent > 0, "+", ""), exponent)
	}

	res
}

numberFormatLatex = function(x){
	sapply(x, numberFormat_single, type = "latex")
}

numberFormatNormal = function(x){
	sapply(x, numberFormat_single)
}

enumerate_items = function (x, endVerb = c("is", "no", "contain"), addS = FALSE, past = FALSE){
	# function that enumerates items and add verbs
	endVerb = match.arg(endVerb)
	n = length(x)

	if(past){
		endWord = switch(endVerb, is = ifelse(n == 1, " was", " were"), no = "", contain = "contained")
	} else {
		endWord = switch(endVerb, is = ifelse(n == 1, " is", " are"), no = "", contain = ifelse(n == 1, " contains", " contain"))
	}

	if (addS) {
		startWord = ifelse(n == 1, " ", "s ")
	} else {
		startWord = ""
	}

	if (n == 1) {
		res = paste0(startWord, x, endWord)
	} else {
		res = paste0(startWord, paste0(x[-n], collapse = ", "), " and ", x[n], endWord)
	}

	res
}

charShorten = function(x, width){
	# transforms characters => they can't go beyond a certain width
	# two dots are added to suggest longer character
	# charShorten("bonjour", 5) => "bon.."
	n = nchar(x)

	if(n > width){
		res = substr(x, 1, width - 2)
		res = paste0(res, "..")
	} else {
		res = x
	}

	res
}


show_vars_limited_width = function(charVect, nbChars = 60){
	# There are different cases

	n = length(charVect)

	if(n==1){
		text = paste0(charVect, ".")
		return(text)
	}

	nb_char_vect = nchar(charVect)
	sumChar = cumsum(nb_char_vect) + (0:(n-1))*2 + 3 + 1

	if(max(sumChar) < nbChars){
		text = paste0(paste0(charVect[-n], collapse = ", "), " and ", charVect[n])
		return(text)
	}

	qui = max(which.max(sumChar > nbChars - 8) - 1, 1)

	nb_left = n - qui

	if(nb_left == 1){
		text = paste0(paste0(charVect[1:qui], collapse = ", "), " and ", nb_left, " other.")
	} else {
		text = paste0(paste0(charVect[1:qui], collapse = ", "), " and ", nb_left, " others.")
	}

	return(text)
}


char2num = function(x){
	# we transform the data to numeric => faster analysis

	# special case
	qui = which(x == "")
	if(length(qui) > 0){
		x[qui] = "xxEMPTYxx"
	}

	x_unik = unique(x)
	dict = 1:length(x_unik)
	names(dict) = x_unik
	x_num = dict[x]

	names(x_num) = NULL

	x_num
}

quickUnclassFactor = function(x){
	# does as unclass(as.factor(x))
	# but waaaaay quicker

	if(!is.numeric(x)){
		x = as.character(x)
		res = char2num(x)
		return(res)
	}

	myOrder = order(x)
	x_sorted = x[myOrder]
	g = cpp_unclassFactor(x_sorted)
	res = g[order(myOrder)]

	return(res)
}


getItems = function(x){
	# to get the unique elements of x before quickunclassfactor
	# needs to be done because differs depending on the type of x

	if(!is.numeric(x)){
		x = as.character(x)
		res = unique(x)
	} else {
		res = sort(unique(x))
	}

	return(res)
}

prepare_matrix = function(fml, base){
	# This function is faster than model.matrix when there is no factor, otherwise model.matrix is faster
	# The argument fml **MUST** not have factors!

	rhs = fml[c(1,3)]
	t = terms(rhs, data = base)

	all_var_names = attr(t, "term.labels")
	all_vars = gsub(":", "*", all_var_names)

	all_vars_call = parse(text = paste0("list(", paste0(all_vars, collapse = ", "), ")"))
	data_list <- eval(all_vars_call, base)
	names(data_list) = all_var_names

	if(attr(t, "intercept") == 1){
		data_list <- c(list("(Intercept)" = rep(1, nrow(base))), data_list)
	}

	res = do.call("cbind", data_list)

	res
}

#' Finds observations to be removed from ML estimation with factors/clusters
#'
#' For Poisson, Negative Binomial or Logit estimations with fixed-effects, when the dependent variable is only equal to 0 (or 1 for Logit) for one cluster value this leads to a perfect fit for that cluster value by setting its associated cluster coefficient to \code{-Inf}. Thus these observations need to be removed before estimation. This function gives the observations to be removed. Not that by default the function \code{\link[FENmlm]{femlm}} drops them before performing the estimation.
#'
#' @param fml A formula contaning the dependent variable and the clusters. It can be of the type: \code{y ~ cluster_1 + cluster_2} or \code{y ~ x1 | cluster_1 + cluster_1} (in which case variables before the pipe are ignored).
#' @param data A data.frame containing the variables in the formula.
#' @param family Character scalar: either \dQuote{poisson} (default), \dQuote{negbin} or \dQuote{logit}.
#'
#' @return
#' It returns an integer vector of observations to be removed. If no observations are to be removed, an empty integer vector is returned. In both cases, it is of class \code{femlm.obs2remove}.
#' The vector has an attribut \code{cluster} which is a list giving the IDs of the clusters that have been removed, for each cluster dimension.
#'
#' @examples
#'
#' base = iris
#' # v6: Petal.Length with only 0 values for 'setosa'
#' base$v6 = base$Petal.Length
#' base$v6[base$Species == "setosa"] = 0
#'
#' (x = obs2remove(v6 ~ Species, base))
#' attr(x, "cluster")
#'
#' # The two results are identical:
#' res_1 = femlm(v6 ~ Petal.Width | Species, base)
#' # => warning + obsRemoved is created
#'
#' res_2 = femlm(v6 ~ Petal.Width | Species, base[-x, ])
#' # => no warning because observations are removed before
#'
#' res2table(res_1, res_2)
#'
#' all(res_1$obsRemoved == x)
#'
obs2remove = function(fml, data, family = c("poisson", "negbin", "logit")){
	# in the formula, the clusters must be there:
	# either y ~ cluster_1 + cluster_2
	# either y ~ x1 + x2 | cluster_1 + cluster_2

	#
	# CONTROLS
	#

	# FAMILY

	family = match.arg(family)

	# FML

	if(!"formula" %in% class(fml) || length(fml) != 3){
		stop("Argument 'fml' must be a formula of the type: 'y ~ x1 | cluster_1 + cluster_1' or of the type 'y ~ cluster_1 + cluster_2'.")
	}

	FML = Formula::Formula(fml)
	n_rhs = length(FML)[2]

	if(n_rhs > 2){
		stop("Argument 'fml' must be a formula of the type: 'y ~ x1 | cluster_1 + cluster_1' or of the type 'y ~ cluster_1 + cluster_2'.")
	}

	# DATA

	if(is.matrix(data)){
		if(is.null(colnames(data))){
			stop("If argument data is to be a matrix, its columns must be named.")
		}
		data = as.data.frame(data)
	}
	# The conversion of the data (due to data.table)
	if(!"data.frame" %in% class(data)){
		stop("The argument 'data' must be a data.frame or a matrix.")
	}
	if("data.table" %in% class(data)){
		# this is a local change only
		class(data) = "data.frame"
	}

	dataNames = names(data)

	# Extracting the variables
	vars_left = all.vars(formula(FML, lhs=1, rhs=0))
	cluster_fml = formula(FML, lhs=0, rhs=n_rhs)
	vars_clusters = all.vars(cluster_fml)

	if(length(left_missing <- setdiff(vars_left, dataNames)) > 0){
		stop("Left hand side could not be evaluated, following variables are missing from the data: ", paste0(left_missing, collapse = ", "), ".")
	}

	if(length(right_missing <- setdiff(vars_clusters, dataNames)) > 0){
		stop("The clsuters could not be evaluated, following variables are missing from the data: ", paste0(right_missing, collapse = ", "), ".")
	}

	# Evaluation variables
	lhs = as.vector(eval(fml[[2]], data))
	cluster_mat = model.frame(cluster_fml, data)
	cluster_name = names(cluster_mat)

	#
	# -- CORE --
	#

	Q = length(cluster_name)
	dummyOmises = list()
	obs2remove = c()
	for(q in 1:Q){

		dum_raw = cluster_mat[, q]

		thisNames = getItems(dum_raw)
		dum = quickUnclassFactor(dum_raw)
		k = length(thisNames)

		# We delete "all zero" outcome
		sum_y_clust = cpp_tapply_vsum(k, lhs, dum)
		n_perClust = cpp_table(k, dum)

		if(family %in% c("poisson", "negbin")){
			qui = which(sum_y_clust == 0)
		} else if(family == "logit"){
			qui = which(sum_y_clust == 0 | sum_y_clust == n_perClust)
		}

		if(length(qui > 0)){
			# We first delete the data:
			dummyOmises[[q]] = thisNames[qui]
			obs2remove = unique(c(obs2remove, which(dum %in% qui)))
		} else {
			dummyOmises[[q]] = character(0)
		}
	}

	names(dummyOmises) = cluster_name

	if(length(obs2remove) == 0){
		print("No observation to be removed.")
		obs2remove = integer(0)
		class(obs2remove) = "femlm.obs2remove"
		return(invisible(obs2remove))
	}

	class(obs2remove) = "femlm.obs2remove"
	attr(obs2remove, "family") = family
	attr(obs2remove, "cluster") = dummyOmises

	return(obs2remove)
}



#' Print method for femlm.obs2remove objects
#'
#' This function show synthetizes the information of function \code{\link[FENmlm]{obs2remove}}. It reports the number of observations to be removed as well as the number of clusters removed per cluster dimension.
#'
#' @method print femlm.obs2remove
#'
#' @param x A \code{femlm.obs2remove} object obtained from function \code{\link[FENmlm]{obs2remove}}.
#' @param ... Not currently used.
#'
#'
#' @examples
#' base = iris
#' # v6: Petal.Length with only 0 values for 'setosa'
#' base$v6 = base$Petal.Length
#' base$v6[base$Species == "setosa"] = 0
#'
#' (x = obs2remove(v6 ~ Species, base))
#' attr(x, "cluster")
#'
print.femlm.obs2remove = function(x, ...){

	if(length(x) == 0){
		print("No observation to be removed.")
	} else {
		cat(length(x), " observations removed because of only zero", ifelse(attr(x, "family") == "logit", ", or only one,", ""), " outcomes.\n", sep = "")
		cluster = attr(x, "cluster")
		cat("# clusters removed: ", paste0(names(cluster), ": ", lengths(cluster), collapse = ", "), ".", sep = "")
	}

}

#' Collinearity diagnostics for femlm objects
#'
#' In some occasions, the optimization algorithm of \code{\link[FENmlm]{femlm}} may fail to converge, or the variance-covariance matrix may not be available. The most common reason of why this happens is colllinearity among variables. This function helps to find out which variable is problematic.
#'
#' @param x A \code{femlm} object obtained from function \code{\link[FENmlm]{femlm}}.
#'
#' @details
#' This function tests: 1) collinearity with the cluster variables, 2) perfect multi-collinearity between the variables, and 3) identification issues when there are non-linear in parameters parts.
#'
#' @return
#' It returns a text message with the identified diagnostics.
#'
#' @examples
#'
#' # Creating an example data base:
#' cluster_1 = sample(3, 100, TRUE)
#' cluster_2 = sample(20, 100, TRUE)
#' x = rnorm(100, cluster_1)**2
#' y = rnorm(100, cluster_2)**2
#' z = rnorm(100, 3)**2
#' dep = rpois(100, x*y*z)
#' base = data.frame(cluster_1, cluster_2, x, y, z, dep)
#'
#' # creating collinearity problems:
#' base$v1 = base$v2 = base$v3 = base$v4 = 0
#' base$v1[base$cluster_1 == 1] = 1
#' base$v2[base$cluster_1 == 2] = 1
#' base$v3[base$cluster_1 == 3] = 1
#' base$v4[base$cluster_2 == 1] = 1
#'
#' # Estimations:
#'
#' # Collinearity with the cluster variables:
#' res_1 = femlm(dep ~ log(x) + v1 + v2 + v4 | cluster_1 + cluster_2, base)
#' diagnostic(res_1)
#' # => collinearity with cluster identified, we drop v1 and v2
#' res_1bis = femlm(dep ~ log(x) + v4 | cluster_1 + cluster_2, base)
#' diagnostic(res_1bis)
#'
#' # Multi-Collinearity:
#' res_2 =  femlm(dep ~ log(x) + v1 + v2 + v3 + v4, base)
#' diagnostic(res_2)
#'
#' # In non-linear part:
#' res_3 = femlm(dep ~ log(z), base, NL.fml = ~log(a*x + b*y),
#'               NL.start = list(a=1, b=1), lower = list(a=0, b=0))
#' diagnostic(res_3)
#'
#'
diagnostic = function(x){
	# x: femlm estimation

	if(class(x) != "femlm"){
		stop("Argument 'x' must be a femlm object.")
	}

	# I) (linear) collinearity with clusters
	# II) (linear) multi collinearity
	# III) (non-linear) overidentification

	# flags
	isCluster = !is.null(x$id_dummies)
	linear_fml = formula(Formula(formula(x, "linear")), lhs=0, rhs=1)
	isLinear = length(all.vars(linear_fml)) > 0
	NL_fml = x$NL.fml
	isNL = !is.null(NL_fml)
	coef = x$coefficients

	# Getting the data
	data = NULL
	try(data <- eval(x$call$data, parent.frame()))

	if(is.null(data)){
		dataName = x$call$data
		stop("To apply 'diagnostic', we fetch the original database in the parent.frame -- but it doesn't seem to be there anymore (btw it was ", deparse(dataName), ").")
	}

	if(!is.null(x$obsRemoved)){
		data = data[-x$obsRemoved, ]
	}

	if(isCluster){
		linear_fml = update(linear_fml, ~ . + 1)
	}

	if(isLinear || isCluster || "(Intercept)" %in% names(coef)){
		linear.matrix = model.matrix(linear_fml, data)
	}

	#
	# I) collinearity with clusters
	#

	if(isCluster && isLinear){
		# We project each variable onto the cluster subspace
		linear_mat_noIntercept = linear.matrix[, -1, drop = FALSE]

		cluster = x$id_dummies
		Q = length(cluster)
		for(q in 1:Q){
			dum = cluster[[q]]
			k = max(dum)
			value = cpp_tapply_sum(Q = k, x = linear_mat_noIntercept, dum = dum)
			nb_per_cluster = cpp_table(Q = k, dum = dum)

			# residuals of the linear projection on the cluster space
			residuals = linear_mat_noIntercept - (value/nb_per_cluster)[dum, ]

			sum_residuals = colSums(abs(residuals))

			if(any(sum_residuals < 1e-4)){
				varnames = colnames(linear_mat_noIntercept)
				collin_var = varnames[sum_residuals < 1e-4]
				if(length(collin_var) == 1){
					message = paste0("Variable '", collin_var, "' is collinear with cluster ", names(cluster)[q], ".")
				} else {
					message = paste0("Variables ", show_vars_limited_width(collin_var), " are collinear with cluster '", names(cluster)[q], "'.")
				}

				print(message)
				return(invisible(message))

			}

		}
	}

	#
	# Variables are equal to 0
	#

	if(isLinear){
		sum_all = colSums(abs(linear.matrix))
		if(any(sum_all == 0)){
			var_problem = colnames(linear.matrix)[sum_all == 0]
			message = paste0("Variable", enumerate_items(var_problem, addS = TRUE), " constant and equal to 0.")

			print(message)
			return(invisible(message))
		}
	}

	#
	# II) perfect multicollinearity
	#

	name2change = grepl("\\)[[:alnum:]]", colnames(linear.matrix))
	if(any(name2change)){
		linbase = as.data.frame(linear.matrix)
		names(linbase)[name2change] = gsub("(\\(|\\))", "_", names(linbase)[name2change])
		linearVars = setdiff(names(linbase), "(Intercept)")
	} else {
		linbase = as.data.frame(linear.matrix)
		linearVars = setdiff(colnames(linear.matrix), "(Intercept)")
	}

	# we add possibly missing variables
	varmiss = setdiff(all.vars(linear_fml), names(linbase))
	for(v in varmiss) linbase[[v]] = data[[v]]

	# browser()

	# linearVars = setdiff(colnames(linear.matrix), "(Intercept)")
	if(isLinear && length(linearVars) >= 2){

		for(v in linearVars){
			fml2estimate = as.formula(paste0(v, "~", paste0(setdiff(linearVars, v), collapse = "+")))

			# res = lm(fml2estimate, data)
			# res = lm(fml2estimate, as.data.frame(linear.matrix))
			res = lm(fml2estimate, linbase)

			sum_resid = sum(abs(resid(res)))
			if(sum_resid < 1e-4){
				coef_lm = coef(res)
				collin_var = names(coef_lm)[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
				message = paste0("Variable '", v, "' is collinear with variable", ifelse(length(collin_var)>=2, "s", ""), ": ", paste0(collin_var, collapse = ", "), ".")

				print(message)
				return(invisible(message))
			}
		}
	}

	#
	# II.b) perfect multicollinearity + cluster
	#

	# linearVars = setdiff(colnames(linear.matrix), "(Intercept)")
	if(isLinear && length(linearVars) >= 2 && isCluster){

		dum_names = names(x$id_dummies)
		n_clust = length(dum_names)
		new_dum_names = paste0("__DUM_", 1:n_clust)
		for(i in 1:n_clust){
			# data[[paste0("__DUM_", i)]] = x$id_dummies[[i]]
			linbase[[paste0("__DUM_", i)]] = x$id_dummies[[i]]
		}

		for(v in linearVars){
			fml2estimate = as.formula(paste0(v, "~", paste0(setdiff(linearVars, v), collapse = "+")))

			for(id_cluster in 1:n_clust){
				# res = femlm(fml2estimate, data, cluster = new_dum_names[id_cluster], family = "gaussian", showWarning = FALSE)
				res = femlm(fml2estimate, linbase, cluster = new_dum_names[id_cluster], family = "gaussian", showWarning = FALSE)

				sum_resid = sum(abs(resid(res)))
				if(sum_resid < 1e-4){
					coef_lm = coef(res)
					collin_var = names(coef_lm)[!is.na(coef_lm) & abs(coef_lm) > 1e-6]
					message = paste0("Variable '", v, "' is collinear with variable", ifelse(length(collin_var)>=2, "s", ""), ": ", paste0(collin_var, collapse = ", "), " and cluster ", dum_names[id_cluster], ".")

					print(message)
					return(invisible(message))
				}
			}
		}
	}

	#
	# III) NL problem
	#

	if(isNL){
		NL_vars = all.vars(NL_fml)
		varNotHere = setdiff(NL_vars, c(names(coef), names(data)))
		if(length(varNotHere) > 0){
			stop("Some variables used to estimate the model (in the non-linear formula) are missing from the original data: ", paste0(varNotHere, collapse = ", "), ".")
		}

		var2send = intersect(NL_vars, names(data))
		env = new.env()
		for(var in var2send){
			assign(var, data[[var]], env)
		}
	}

	if(isNL && (length(coef) >= 2 || isCluster)){

		if(isCluster){
			# we add the constant
			coef["CONSTANT"] = 1
			data$CONSTANT = 1
			linear.matrix = model.matrix(update(linear_fml, ~.-1+CONSTANT), data)
		}

		coef_name = names(coef)

		NL_coef = setdiff(all.vars(NL_fml), names(data))
		# L_coef = setdiff(coef_name, NL_coef)
		L_coef = colnames(linear.matrix)

		#
		# We compute mu:
		#

		mu = 0
		if(length(L_coef) > 0){
			mu = linear.matrix %*% coef[L_coef]
		}

		for(iter_coef in NL_coef){
			assign(iter_coef, coef[iter_coef], env)
		}

		# Evaluation of the NL part
		value_NL = eval(NL_fml[[2]], env)
		mu = mu + value_NL
		data$mu = mu

		if(var(mu) == 0){
			message = "Variance of the NL part is 0."
			print(message)
			return(invisible(message))
		}

		#
		# The loop
		#

		for(var in coef_name){
			# we modify the coef of var, then we fix it:
			# if we can find a set of other parameters that give the same fit,
			# then there is an identification issue

			if(var %in% L_coef){
				# if linear coef: we modify the fml and add an offset
				if(length(L_coef) == 1){
					fml = mu ~ 0
				} else {
					fml = as.formula(paste0("mu ~ 0+", paste0(setdiff(L_coef, var), collapse = "+")))
				}

				# offset with the new value
				offset = as.formula(paste0("~", coef[var] + 1, "*", var))

				res = femlm(fml, data = data, family = "gaussian", NL.fml = NL_fml, NL.start = as.list(coef), offset = offset)

			} else {
				# NL case:
				# we modify both fml and NL.fml

				if(isCluster){
					fml = update(x$fml, mu ~ . - 1 + CONSTANT)
				} else {
					fml = update(x$fml, mu ~ .)
				}

				# the new formula with the new variable
				NL_fml_char = as.character(NL_fml)[2]
				NL_fml_new = as.formula(paste0("~", gsub(paste0("\\b", var, "\\b"), 1 + coef[var], NL_fml_char)))

				if(length(NL_coef) == 1){
					# there is no more parameter to estimate => offset
					res = femlm(fml, data = data, family = "gaussian", offset = NL_fml_new)
				} else {
					res = femlm(fml, data = data, family = "gaussian", NL.fml = NL_fml_new, NL.start = as.list(coef))
				}

			}

			sum_resids = sum(abs(resid(res)))
			if(sum_resids < 1e-4){
				coef_esti = coef(res)
				coef_diff = abs(coef_esti - coef[names(coef_esti)])
				collin_var = names(coef_diff)[coef_diff > 1e-3]
				message = paste0("Coefficients ", show_vars_limited_width(c(var, collin_var)), " are not uniquely identifed.")

				print(message)
				return(invisible(message))
			}

		}

	}

	message = "No visible collinearity problem. (Doesn't mean there's none!)"
	print(message)
	return(invisible(message))

}


#### ................. ####
#### Aditional Methods ####
####

# Here we add common statistical functions

#' Extract the number of observations form a femlm object
#'
#' This function simply extracts the number of obsrvations used to estimate a \code{\link[FENmlm]{femlm}} model.
#'
#'
#' @param object An object of class \code{femlm}. Typically the result of a \code{\link[FENmlm]{femlm}} estimation.
#' @param ... Not currently used.
#'
#' @seealso
#' See also the main estimation functions \code{\link[FENmlm]{femlm}}. Use \code{\link[FENmlm]{summary.femlm}} to see the results with the appropriate standard-errors, \code{\link[FENmlm]{getFE}} to extract the cluster coefficients, and the functions \code{\link[FENmlm]{res2table}} and \code{\link[FENmlm]{res2tex}} to visualize the results of multiple estimations.
#'
#' @author
#' Laurent Berge
#'
#' @return
#' It returns an interger.
#'
#' @examples
#'
#' # simple estimation on iris data, clustering by "Species"
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' nobs(res)
#' logLik(res)
#'
#'
nobs.femlm = function(object, ...){
	object$n
}

#' Aikake's an information criterion
#'
#' This function computes the AIC (Aikake's, an information criterion) from a \code{\link[FENmlm]{femlm}} estimation.
#'
#' @inheritParams nobs.femlm
#'
#' @param ... Optionally, more fitted objects.
#' @param k A numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC (i.e. \code{AIC=-2*LL+k*nparams}).
#'
#' @details
#' The AIC is computed as:
#' \deqn{AIC = -2\times LogLikelihood + k\times nbParams}
#' with k the penalty parameter.
#'
#' You can have more information on this crtierion on \code{\link[stats]{AIC}}.
#'
#' @return
#' It return a numeric vector, with length the same as the number of objects taken as arguments.
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{AIC.femlm}}, \code{\link[FENmlm]{logLik.femlm}}, \code{\link[FENmlm]{nobs.femlm}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#' res2 = femlm(Sepal.Length ~ Petal.Width | Species, iris)
#'
#' AIC(res1, res2)
#' BIC(res1, res2)
#'
#'
AIC.femlm = function(object, ..., k = 2){

	dots = list(...)
	if(length(dots) > 0){
		# we check consistency with observations
		nobs_all = c(nobs(object), sapply(dots, nobs))

		if(any(diff(nobs_all) != 0)){
			warning("Models are not all fitted to the same number of observations.")
		}

		otherAIC = sapply(dots, AIC)
	} else {
		otherAIC = c()
	}

	all_AIC = c(-2*object$loglik + k*object$nparams, otherAIC)

	all_AIC
}

#' Bayesian information criterion
#'
#' This function computes the BIC (Bayesian information criterion) from a \code{\link[FENmlm]{femlm}} estimation.
#'
#'
#' @inheritParams nobs.femlm
#'
#' @param ... Optionally, more fitted objects.
#'
#' @details
#' The BIC is computed as follows:
#' \deqn{BIC = -2\times LogLikelihood + \log(nobs)\times nbParams}
#' with k the penalty parameter.
#'
#' You can have more information on this crtierion on \code{\link[stats]{AIC}}.
#'
#' @return
#' It return a numeric vector, with length the same as the number of objects taken as arguments.
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{AIC.femlm}}, \code{\link[FENmlm]{logLik.femlm}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # two fitted models with different expl. variables:
#' res1 = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#' res2 = femlm(Sepal.Length ~ Petal.Width | Species, iris)
#'
#' AIC(res1, res2)
#' BIC(res1, res2)
#'
BIC.femlm = function(object, ...){

	dots = list(...)
	if(length(dots) > 0){
		# we check consistency with observations
		nobs_all = c(nobs(object), sapply(dots, nobs))

		if(any(diff(nobs_all) != 0)){
			warning("Models are not all fitted to the same number of observations.")
		}

		otherBIC = sapply(dots, BIC)
	} else {
		otherBIC = c()
	}

	all_BIC = c(-2*object$loglik + 2*object$nparams*log(object$n), otherBIC)

	all_BIC
}

#' Extracts the log-likelihood
#'
#' This function extracts the log-likelihood from a \code{\link[FENmlm]{femlm}} estimation.
#'
#' @inheritParams nobs.femlm
#'
#' @param ... Not currently used.
#'
#' @details
#' This function extracts the log-likelihood based on the model fit. You can have more information on the likelihoods in the details of the function \code{\link[FENmlm]{femlm}}.
#'
#' @return
#' It returns a numeric scalar.
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{AIC.femlm}}, \code{\link[FENmlm]{BIC.femlm}}, \code{\link[FENmlm]{nobs.femlm}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # simple estimation on iris data, clustering by "Species"
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' nobs(res)
#' logLik(res)
#'
#'
logLik.femlm = function(object, ...){
	ll = object$loglik
	ll
}

#' Extracts the coefficients from a femlm fit
#'
#' This function extracts the coefficients obtained from a model estimated with \code{\link[FENmlm]{femlm}}.
#'
#' @inheritParams nobs.femlm
#'
#' @param ... Not currently used.
#'
#' @details
#' The coefficients are the ones that have been found to maximize the log-likelihood of the specified model. More information can be found on \code{\link[FENmlm]{femlm}} help page.
#'
#' Note that if the model has been estimated with clusters, to obtain the cluster coefficients, you need to use the function \code{\link[FENmlm]{getFE}}.
#'
#' @return
#' This function returns a named numeric vector.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{confint.femlm}}, \code{\link[FENmlm]{vcov.femlm}}, \code{\link[FENmlm]{res2table}}, \code{\link[FENmlm]{res2tex}}, \code{\link[FENmlm]{getFE}}.
#'
#' @examples
#'
#' # simple estimation on iris data, clustering by "Species"
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # the coefficients of the variables:
#' coef(res)
#'
#' # the cluster coefficients:
#' getFE(res)
#'
#'
coef.femlm = coefficients.femlm = function(object, ...){
	object$coefficients
}

#' @rdname coef.femlm
"coefficients.femlm"


#' Extracts fitted values from a femlm fit
#'
#' This function extracts the fitted values from a model estimated with \code{\link[FENmlm]{femlm}}. The fitted values that are returned are the \emph{expected predictor}.
#'
#' @inheritParams nobs.femlm
#'
#' @param type Character either equal to \code{"response"} (default) or \code{"link"}. If \code{type="response"}, then the output is at the level of the response variable, i.e. it is the expected predictor \eqn{E(Y|X)}. If \code{"link"}, then the output is at the level of the explanatory variables, i.e. the linear predictor \eqn{X\cdot \beta}.
#' @param ... Not currently used.
#'
#' @details
#' This function returns the \emph{expected predictor} of a \code{\link[FENmlm]{femlm}} fit. The likelihood functions are detailed in \code{\link[FENmlm]{femlm}} help page.
#'
#' @return
#' It returns a numeric vector of length the number of observations used to estimate the model.
#'
#' If \code{type = "response"}, the value returned is the expected predictor, i.e. the expected value of the dependent variable for the fitted model: \eqn{E(Y|X)}.
#' If \code{type = "link"}, the value returned is the linear predictor of the fitted model, that is \eqn{X\cdot \beta} (remind that \eqn{E(Y|X) = f(X\cdot \beta)}).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{resid.femlm}}, \code{\link[FENmlm]{predict.femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{vcov.femlm}}, \code{\link[FENmlm]{getFE}}.
#'
#' @examples
#'
#' # simple estimation on iris data, clustering by "Species"
#' res_poisson = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris)
#'
#' # we extract the fitted values
#' y_fitted_poisson = fitted(res_poisson)
#'
#' # Same estimation but in OLS (Gaussian family)
#' res_gaussian = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris, family = "gaussian")
#'
#' y_fitted_gaussian = fitted(res_gaussian)
#'
#' # comparison of the fit for the two families
#' plot(iris$Sepal.Length, y_fitted_poisson)
#' points(iris$Sepal.Length, y_fitted_gaussian, col = 2, pch = 2)
#'
#'
fitted.femlm = fitted.values.femlm = function(object, type = c("response", "link"), ...){

	type = match.arg(type)

	if(type == "response"){
		res = object$fitted.values
	} else if(!is.null(object$mu)){
		res = object$mu
	} else {
		family = object$family
		famFuns = switch(family,
							  poisson = ml_poisson(),
							  negbin = ml_negbin(),
							  logit = ml_logit(),
							  gaussian = ml_gaussian())

		res = famFuns$linearFromExpected(object$fitted.values)
	}

	res
}

#' @rdname fitted.femlm
"fitted.values.femlm"

#' Extracts residuals from a femlm object
#'
#' This function extracts residuals from a fitted model estimated with \code{\link[FENmlm]{femlm}}.
#'
#' @inheritParams nobs.femlm
#'
#' @param ... Not currently used.
#'
#' @details
#' The residuals returned are the difference between the dependent variable and the expected predictor.
#'
#' @return
#' It returns a numeric vector of the length the number of observations used for the estimation.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{fitted.femlm}}, \code{\link[FENmlm]{predict.femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{vcov.femlm}}, \code{\link[FENmlm]{getFE}}.
#'
#' @examples
#'
#' # simple estimation on iris data, clustering by "Species"
#' res_poisson = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'                     Petal.Width | Species, iris)
#'
#' # we plot the residuals
#' plot(resid(res_poisson))
#'
resid.femlm = residuals.femlm = function(object, ...){
	object$residuals
}

#' @rdname resid.femlm
"residuals.femlm"

#' Predict method for femlm fits
#'
#' This function obtains prediction from a fitted model estimated with \code{\link[FENmlm]{femlm}}.
#'
#' @inheritParams nobs.femlm
#' @inheritParams fitted.femlm
#'
#' @param newdata A data.frame containing the variables used to make the prediction. If not provided, the fitted expected (or linear if \code{type = "link"}) predictors are returned.
#' @param ... Not currently used.
#'
#' @return
#' It returns a numeric vector of length equal to the number of observations in argument \code{newdata}.
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{update.femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{vcov.femlm}}, \code{\link[FENmlm]{getFE}}.
#'
#' @examples
#'
#' # Estimation on iris data
#' res = femlm(Sepal.Length ~ Petal.Length | Species, iris)
#'
#' # what would be the prediction if the data was all setosa?
#' newdata = data.frame(Petal.Length = iris$Petal.Length, Species = "setosa")
#' pred_setosa = predict(res, newdata = newdata)
#'
#' # Let's look at it graphically
#' plot(c(1, 7), c(3, 11), type = "n", xlab = "Petal.Length",
#'      ylab = "Sepal.Length")
#'
#' newdata = iris[order(iris$Petal.Length), ]
#' newdata$Species = "setosa"
#' lines(newdata$Petal.Length, predict(res, newdata))
#'
#' # versicolor
#' newdata$Species = "versicolor"
#' lines(newdata$Petal.Length, predict(res, newdata), col=2)
#'
#' # virginica
#' newdata$Species = "virginica"
#' lines(newdata$Petal.Length, predict(res, newdata), col=3)
#'
#' # The original data
#' points(iris$Petal.Length, iris$Sepal.Length, col = iris$Species, pch = 18)
#' legend("topleft", lty = 1, col = 1:3, legend = levels(iris$Species))
#'
predict.femlm = function(object, newdata, type = c("response", "link"), ...){

	# Controls
	type = match.arg(type)

	# if newdata is missing
	if(missing(newdata)){
		if(type == "response"){
			return(object$fitted.values)
		} else {
			family = object$family
			famFuns = switch(family,
								  poisson = ml_poisson(),
								  negbin = ml_negbin(),
								  logit = ml_logit(),
								  gaussian = ml_gaussian())

			return(famFuns$linearFromExpected(object$fitted.values))
		}
	}

	if(!is.matrix(newdata) && !"data.frame" %in% class(newdata)){
		stop("Argument 'newdata' must be a data.frame.")
	}

	# we ensure it really a clean data.frame
	newdata = as.data.frame(newdata)

	# We deconstruct it in four steps:
	# 1) cluster
	# 2) linear
	# 3) non-linear
	# 4) offset

	n = nrow(newdata)

	# 1) Cluster

	# init cluster values
	value_cluster = 0

	clusterNames = object$clusterNames
	if(!is.null(clusterNames)){

		n_cluster = length(clusterNames)

		# Extraction of the clusters
		id_cluster = list()
		for(i in 1:n_cluster){
			# checking if the variable is in the newdata
			variable = all.vars(parse(text = clusterNames[i]))
			isNotHere = !variable %in% names(newdata)
			if(any(isNotHere)){
				stop("The variable ", variable[isNotHere][1], " is absent from the 'newdata' but is needed for prediction (it is a cluster variable).")
			}

			# Obtaining the unclassed vector of clusters
			cluster_current = eval(parse(text = clusterNames[i]), newdata)
			cluster_current_unik = unique(cluster_current)

			cluster_values_possible = attr(object$id_dummies[[i]],"clust_names")
			valueNotHere = setdiff(cluster_current_unik, cluster_values_possible)
			if(length(valueNotHere) > 0){
				stop("The cluster value ", valueNotHere[1], " (cluster ", clusterNames[i], ") was not used in the initial estimation, prediction cannot be done for observations with that value. Prediction can be done only for cluster values present in the main estimation.")
			}

			cluster_current_num = unclass(factor(cluster_current, levels = cluster_values_possible))
			id_cluster[[i]] = cluster_current_num
		}

		# Value of the cluster coefficients
		cluster_coef = getFE(object)
		for(i in 1:n_cluster){
			cluster_current_num = id_cluster[[i]]
			cluster_coef_current = cluster_coef[[i]]

			value_cluster = value_cluster + cluster_coef_current[cluster_current_num]
		}
	}

	# 2) Linear values

	coef = object$coefficients

	# intercept
	if("(Intercept)" %in% names(coef)){
		value_linear = rep(coef["(Intercept)"], n)
	} else {
		value_linear = 0
	}

	rhs_fml = formula(Formula(object$fml), lhs = 0, rhs = 1)
	if(length(all.vars(rhs_fml)) > 0){
		# Checking all variables are there
		varNotHere = setdiff(all.vars(rhs_fml), names(newdata))
		if(length(varNotHere) > 0){
			stop("Some variables used to estimate the model (in fml) are missing from argument 'newdata': ", paste0(varNotHere, collapse = ", "), ".")
		}

		# we create the matrix
		matrix_linear = model.matrix(rhs_fml, newdata)

		keep = intersect(names(coef), colnames(matrix_linear))
		value_linear = value_linear + as.vector(matrix_linear[, keep, drop = FALSE] %*% coef[keep])
	}

	# 3) Non linear terms

	value_NL = 0
	NL_fml = object$NL.fml
	if(!is.null(NL_fml)){
		# controlling that we can evaluate that
		NL_vars = all.vars(NL_fml)
		varNotHere = setdiff(NL_vars, c(names(coef), names(newdata)))
		if(length(varNotHere) > 0){
			stop("Some variables used to estimate the model (in the non-linear formula) are missing from argument 'newdata': ", paste0(varNotHere, collapse = ", "), ".")
		}

		var2send = intersect(NL_vars, names(newdata))
		env = new.env()
		for(var in var2send){
			assign(var, newdata[[var]], env)
		}

		coef2send = setdiff(NL_vars, names(newdata))
		for(iter_coef in coef2send){
			assign(iter_coef, coef[iter_coef], env)
		}

		# Evaluation of the NL part
		value_NL = eval(NL_fml[[2]], env)
	}

	# 4) offset value

	value_offset = 0
	offset = object$offset
	if(!is.null(offset)){
		# evaluation of the offset
		varNotHere = setdiff(all.vars(offset), names(newdata))
		if(length(varNotHere) > 0){
			stop("Some variables used to estimate the model (in the offset) are missing from argument 'newdata': ", paste0(varNotHere, collapse = ", "), ".")
		}

		value_offset = eval(offset[[length(offset)]], newdata)
	}

	value_predicted = value_cluster + value_linear + value_NL + value_offset

	if(type == "link"){
		return(value_predicted)
	}

	# Now the expected predictor
	family = object$family
	famFuns = switch(family,
						  poisson = ml_poisson(),
						  negbin = ml_negbin(),
						  logit = ml_logit(),
						  gaussian = ml_gaussian())

	if(family == "gaussian"){
		exp_value = 0
	} else {
		exp_value = exp(value_predicted)
	}

	expected.predictor = famFuns$expected.predictor(value_predicted, exp_value)

	expected.predictor
}

#' Extract the variance/covariance of a femlm fit
#'
#' This function extracts the variance-covariance of estimated parameters from a model estimated with \code{\link[FENmlm]{femlm}}.
#'
#' @inheritParams summary.femlm
#' @inheritParams nobs.femlm
#'
#' @param ... Other arguments to be passed to \code{\link[FENmlm]{summary.femlm}}.
#'
#' The computation of the VCOV matrix is first done in \code{\link[FENmlm]{summary.femlm}}.
#'
#' @return
#' It returns a \eqn{N\times N} square matrix where \eqn{N} is the number of variables of the fitted model.
#' This matrix has an attribute \dQuote{type} specifying how this variance/covariance matrix has been commputed (i.e. was it created using White correction, or was it clustered along a specific factor, etc).
#'
#' @author
#' Laurent Berge
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{confint.femlm}}, \code{\link[FENmlm]{resid.femlm}}, \code{\link[FENmlm]{predict.femlm}}, \code{\link[FENmlm]{getFE}}.
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade (with 3 fixed-effects)
#' est_pois = femlm(Euros ~ log(dist_km) + log(Year) | Origin + Destination +
#'                  Product, trade)
#'
#' # "normal" VCOV
#' vcov(est_pois)
#'
#' # "white" VCOV
#' vcov(est_pois, se = "white")
#'
#' # "clustered" VCOV (with respect to the Origin factor)
#' vcov(est_pois, se = "cluster")
#'
#' # "clustered" VCOV (with respect to the Product factor)
#' vcov(est_pois, se = "cluster", cluster = trade$Product)
#' # another way to make the same request:
#' vcov(est_pois, se = "cluster", cluster = "Product")
#'
#' # Another estimation without cluster:
#' est_pois_simple = femlm(Euros ~ log(dist_km) + log(Year), trade)
#'
#' # We can still get the clustered VCOV,
#' # but we need to give the cluster-vector:
#' vcov(est_pois_simple, se = "cluster", cluster = trade$Product)
#'
#'
vcov.femlm = function(object, se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), cluster, dof_correction=FALSE, forceCovariance = FALSE, keepBounded = FALSE, ...){
	# computes the clustered vcov

	if(!is.null(object$onlyCluster)){
		# means that the estimation is done without variables
		stop("No explanatory variable was used: vcov cannot be applied.")
	}

	sd.val = match.arg(se)

	#
	# non-linear: handling bounded parameters
	#

	# We handle the bounded parameters:
	isBounded = object$isBounded
	if(is.null(isBounded)){
		isBounded = rep(FALSE, length(object$coefficients))
	}

	if(any(isBounded)){
		if(keepBounded){
			# we treat the bounded parameters as regular variables
			myScore = object$score
			object$cov.unscaled = solve(object$hessian)
		} else {
			myScore = object$score[, -which(isBounded), drop = FALSE]
		}
	} else {
		myScore = object$score
	}


	#
	# Core function
	#

	n = object$n
	k = object$nparams

	correction.dof = n / (n-k*dof_correction)

	VCOV_raw = object$cov.unscaled

	# information on the variable used for the clustering
	type_info = ""

	if(anyNA(VCOV_raw)){

		if(!forceCovariance){
			# warning("Standard errors are NA because of likely presence of collinearity. You can use option 'forceCovariance' to try to force the computation of the vcov matrix (to see what's wrong).", call. = FALSE)
			warning("Standard errors are NA because of likely presence of collinearity. Use function diagnostic() to detect collinearity problems.", call. = FALSE)
			return(VCOV_raw)
		} else {
			VCOV_raw_forced = MASS::ginv(object$hessian)
			if(anyNA(VCOV_raw_forced)) {
				stop("The covariance matrix could not be 'forced'.")
			}

			object$cov.unscaled = VCOV_raw_forced
			return(vcov(object, se=sd.val, cluster=cluster, dof_correction=dof_correction))
		}

	} else if(sd.val == "standard"){

		vcov = VCOV_raw * correction.dof

	} else if(sd.val == "white"){

		vcov = crossprod(myScore %*% VCOV_raw) * correction.dof

	} else {
		# Clustered SD!
		nway = switch(sd.val, cluster=1, twoway=2, threeway=3, fourway=4)

		#
		# Controls
		#

		# Controlling the clusters
		do.unclass = TRUE
		if(missing(cluster) || is.null(cluster)){

			if(is.null(object$id_dummies)){
				stop("To display clustered standard errors, you must provide the argument 'cluster'.")

			} else if(length(object$id_dummies) < nway) {
				stop("Since the result was not clustered with ", nway, " clusters, you need to provide the argument 'cluster' with ", nway, "clusters.")

			} else {
				cluster = object$id_dummies[1:nway]

				type_info = paste0(" (", paste0(object$clusterNames[1:nway], collapse = " & "), ")")

				# in that specific case, there is no need of doing unclass.factor because already done
				do.unclass = FALSE
			}

		} else {

			isS = ifelse(nway>1, "s, each", "")
			if(length(cluster) == nway && is.character(cluster)){
				if(any(!cluster %in% object$clusterNames)){
					var_problem = setdiff(cluster, object$clusterNames)
					stop("Cannot apply ", nway, "-way clustering with current 'cluster' argument. Variable", enumerate_items(var_problem, past = TRUE, addS = TRUE), " not used as clusters in the estimation.\nAlternatively, use a matrix for argument cluster (ex: cluster = base[, c('", paste0(cluster, collapse = "', '"), "')]")
				}

				type_info = paste0(" (", paste0(cluster, collapse = " & "), ")")
				cluster = object$id_dummies[cluster]

			} else if(nway == 1){
				if(!is.list(cluster) && (is.vector(cluster) || is.factor(cluster))){
					cluster = list(cluster)

				} else if(! (is.list(cluster) && length(cluster) == 1)){
					stop("For one way clustering, the argument 'cluster' must be either the vector of cluster ids, or a list containing the vector of cluster ids.")

				}
			} else if(! (is.list(cluster) && length(cluster)==nway) ){
				stop("The 'cluster' must be a list containing ", nway, " element", isS, " being the vector of IDs of each observation.")

			}

			cluster = as.list(cluster)
		}

		# now we check the lengths:
		n_per_cluster = sapply(cluster, length)
		if(!all(diff(n_per_cluster) == 0)){
			stop("The vectors of the argument 'cluster' must be of the same length.")
		}

		# Either they are of the same length of the data
		if(n_per_cluster[1] != object$n){
			# Then two cases: either the user introduces the original data and it is OK
			if(n_per_cluster[1] == (object$n + length(object$obsRemoved))){
				# We modify the clusters
				for(i in 1:nway) cluster[[i]] = cluster[[i]][-object$obsRemoved]
			} else {
				# If this is not the case: there is a problem
				stop("The length of the cluster does not match the original data.")
			}
		}

		#
		# Calculus
		#

		# initialisation
		vcov = VCOV_raw * 0

		if(do.unclass){
			for(i in 1:nway){
				cluster[[i]] = quickUnclassFactor(cluster[[i]])
			}
		}

		for(i in 1:nway){

			myComb = combn(nway, i)

			power = floor(1 + log10(sapply(cluster, max)))

			for(j in 1:ncol(myComb)){

				if(i == 1){
					index = cluster[[myComb[, j]]]
				} else if(i > 1){

					vars = myComb[, j]

					if(sum(power[vars]) > 14){
						myDots = cluster[vars]
						myDots$sep = "_"
						index = do.call("paste", myDots)
					} else {
						# quicker, but limited by the precision of integers
						index = cluster[[vars[1]]]
						for(k in 2:length(vars)){
							index = index + cluster[[vars[k]]]*10**sum(power[vars[1:(k-1)]])
						}
					}

					index = quickUnclassFactor(index)

				}

				vcov = vcov + (-1)**(i+1) * vcovClust(index, VCOV_raw, myScore, dof_correction, do.unclass=FALSE)

			}
		}
	}

	if(any(diag(vcov)<0)){
		warning("Some variances are negative (likely problem in the model).")
	}

	sd.dict = c("standard" = "Standard", "white"="White", "cluster"="Clustered", "twoway"="Two-way", "threeway"="Three-way", "fourway"="Four-way")
	attr(vcov, "type") = paste0(as.vector(sd.dict[sd.val]), type_info)

	vcov
}


#' Confidence interval for parameters estimated with femlm
#'
#' This function computes the confidence interval of parameter estimates obtained from a model estimated with \code{\link[FENmlm]{femlm}}.
#'
#' @inheritParams nobs.femlm
#' @inheritParams vcov.femlm
#'
#' @param parm The parameters for which to compute the confidence interval (either an integer vector OR a character vector with the parameter name). If missing, all parameters are used.
#' @param level The confidence level. Default is 0.95.
#'
#' @return
#' Returns a data.frame with two columns giving respectively the lower and upper bound of the confidence interval. There is as many rows as parameters.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Load trade data
#' data(trade)
#'
#' # We estimate the effect of distance on trade (with 3 cluster effects)
#' est_pois = femlm(Euros ~ log(dist_km) + log(Year) | Origin + Destination +
#'                  Product, trade)
#'
#' # confidence interval with "normal" VCOV
#' confint(est_pois)
#'
#' # confidence interval with "clustered" VCOV (w.r.t. the Origin factor)
#' confint(est_pois, se = "cluster")
#'
#'
confint.femlm = function(object, parm, level = 0.95, se=c("standard", "white", "cluster", "twoway", "threeway", "fourway"), cluster, dof_correction=FALSE, ...){

	# Control
	if(!is.numeric(level) || !length(level) == 1 || level >= 1 || level <= .50){
		stop("The argument 'level' must be a numeric scalar greater than 0.50 and strictly lower than 1.")
	}

	# the parameters for which we should compute the confint
	all_params = names(object$coefficients)

	if(missing(parm)){
		parm_use = all_params
	} else if(is.numeric(parm)){
		if(any(parm %% 1 != 0)){
			stop("If the argument 'parm' is numeric, it must be integers.")
		}

		parm_use = unique(na.omit(all_params[parm]))
		if(length(parm_use) == 0){
			stop("There are ", length(all_params), " coefficients, the argument 'parm' does not correspond to any of them.")
		}
	} else if(is.character(parm)){
		parm_pblm = setdiff(parm, all_params)
		if(length(parm_pblm) > 0){
			stop("some parameters of 'parm' have no estimated coefficient: ", paste0(parm_pblm, collapse=", "), ".")
		}

		parm_use = intersect(parm, all_params)
	}

	# The proper SE
	sum_object = summary(object, se = se, cluster = cluster, dof_correction = dof_correction, ...)

	se_all = sum_object$se
	coef_all = object$coefficients

	# multiplicative factor
	val = (1 - level) / 2
	fact <- abs(qnorm(val))

	# The confints
	lower_bound = coef_all[parm_use] - fact * se_all[parm_use]
	upper_bound = coef_all[parm_use] + fact * se_all[parm_use]

	res = data.frame(lower_bound, upper_bound, row.names = parm_use)
	names(res) = paste0(round(100*c(val, 1-val), 1), " %")

	res
}

#' Updates a femlm estimation
#'
#' Updates and re-estimates a \code{\link[FENmlm]{femlm}} model. This function updates the formulas and use previous starting values to estimate a new \code{\link[FENmlm]{femlm}} model. The data is obtained from the original \code{call}.
#'
#' @method update femlm
#'
#' @inheritParams nobs.femlm
#'
#' @param fml.update Changes to be made to the original argument \code{fml}. See more information on \code{\link[stats]{update.formula}}. You can add/withdraw both variables and clusters. E.g. \code{. ~ . + x2 | . + z2} would add the variable \code{x2} and the cluster \code{z2} to the former estimation.
#' @param ... Other arguments to be passed to the function \code{\link[FENmlm]{femlm}}.
#'
#' @return
#' It returns a \code{\link[FENmlm]{femlm}} object (see details in \code{\link[FENmlm]{femlm}}.
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{predict.femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{vcov.femlm}}, \code{\link[FENmlm]{getFE}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # Example using trade data
#' data(trade)
#'
#' # main estimation
#' est_pois <- femlm(Euros ~ log(dist_km) | Origin + Destination, trade)
#'
#' # we add the variable log(Year)
#' est_2 <- update(est_pois, . ~ . + log(Year))
#'
#' # we add another cluster: "Product"
#' est_3 <- update(est_2, . ~ . | . + Product)
#'
#' # we remove the cluster "Origin" and the variable log(dist_km)
#' est_4 <- update(est_3, . ~ . - log(dist_km) | . - Origin)
#'
#' # Quick look at the 4 estimations
#' res2table(est_pois, est_2, est_3, est_4)
#'
update.femlm = function(object, fml.update, ...){
	# Update method
	# fml.update: update the formula
	# If 1) SAME DATA and 2) SAME dep.var, then we make initialisation


	if(missing(fml.update)){
		fml.update = . ~ .
	} else {
		if(!"formula" %in% class(fml.update)){
			stop("The argument 'fml.update' is required.")
		}
	}

	FML = Formula(fml.update)
	call_new = match.call()
	dots = list(...)

	dot_names = names(dots)
	if("cluster" %in% dot_names){
		stop("Argument 'cluster' is not accepted in the 'update.femlm' method. Please make modifications to clusters directly in the argument 'fml.update'. (E.g. .~.|.+v5 to add variable v5 as a cluster.)")
	}

	if(any(dot_names == "")){
		call_new_names = names(call_new)
		problems = call_new[call_new_names == ""][-1]
		stop("In 'update.femlm' the arguments of '...' are passed to the function femlm, and must be named. Currently there are un-named arguments (e.g. '", deparse(problems[[1]]), "').")
	}



	#
	# I) Linear formula update
	#

	fml_old = object$fml
	fml = update(fml_old, formula(FML, lhs = 1, rhs = 1))
	fml_char = as.character(fml)

	useInit = TRUE
	if(fml[[2]] != fml_old[[2]]){
		# means different dependent variables
		# 	=> initialisation with past parameters is useless
		useInit = FALSE
	}

	if(!is.null(dots$family)){
		family_new = match.arg(dots$family, c("poisson", "negbin", "gaussian", "logit"))
		if(family_new != object$family){
			# if different families: initialisation is useless
			useInit = FALSE
		}
	}

	if(!is.null(dots$na.rm) && dots$na.rm){
		# Too complicated to initialize with na.rm
		# I would have to make tons of controls, and it would work
		# only in some cases...
		useInit = FALSE
	}

	#
	# II) evaluation data
	#

	# We find out if it is the same data
	#	=> only to find out if init is useful
	if(useInit){
		# we test if we can use the initialisation of the parameters

		if(is.null(dots$data)){
			# evaluation
			data = NULL
			try(data <- eval(object$call$data, parent.frame()))

			if(is.null(data)){
				dataName = object$call$data
				stop("To apply 'update.femlm', we fetch the original database in the parent.frame -- but it doesn't seem to be there anymore (btw it was ", deparse(dataName), ").")
			}
		} else {
			if(!is.matrix(data) && !"data.frame" %in% class(dots$data)){
				stop("The argument 'data' must be a data.frame.")
			}
			data = dots$data
		}

		# if same data => we apply init
		n_old = object$n + length(object$obsRemoved)
		n_new = nrow(data)
		if(n_old != n_new){
			useInit = FALSE
		}
	}

	#
	# III) cluster updates
	#

	clusterNames = object$clusterNames
	clusterStart = NULL
	clusterFromUpdate = FALSE
	if(length(FML)[2] > 1){
		# modification of the clusters
		cluster_old = as.formula(paste0("~", paste0(c(1, clusterNames), collapse = "+")))
		cluster_new = update(cluster_old, formula(FML, lhs = 0, rhs = 2))

		if(useInit){
			# Only if we use the init => the starting cluster values
			isThere = sapply(clusterNames, function(x) grepl(x, as.character(cluster_new)[2]))
			if(is.null(object$obsRemoved)){
				# ONLY when there is no cluster removed (otherwise computationaly too complex to be worth)
				if(all(isThere)){
					# we just have to put the old
					clusterStart = object$sumFE
				} else if(any(isThere)){
					# we use the dummies only for the ones that are there
					my_fe = getFE(object)
					clusterStart = 0
					for(i in which(isThere)){
						clusterStart = clusterStart + my_fe[[i]][object$id_dummies[[i]]]
					}
				}
			}
		}

		if(length(all.vars(cluster_new)) > 0){
			# means there is a cluster
			fml_new = as.formula(paste0(fml_char[2], "~", fml_char[3], "|", as.character(cluster_new)[2]))
		} else {
			# there is no cluster
			fml_new = fml
		}

	} else if(!is.null(clusterNames)){
		# Means we keep the same clusters
		clusterFromUpdate = TRUE

		# the starting value:
		clusterStart = object$sumFE

		# the formula updated:
		fml_new = as.formula(paste0(fml_char[2], "~", fml_char[3], "|", paste0(clusterNames, collapse = "+")))

	} else {
		# there is no cluster in the initial model
		fml_new = fml
	}


	#
	# The call
	#

	call_old = object$call

	# we drop the argument cluster from old call (now it's in the fml_new)
	call_old$cluster = NULL

	# new call: call_clear
	call_clear = call_old
	for(arg in names(call_new)[-1]){
		call_clear[[arg]] = call_new[[arg]]
	}

	call_clear$fml = as.call(fml_new)

	if(useInit){
		# we can use the initialisation of parameters
		if(is.null(dots$linear.start) &&  is.null(object$onlyCluster)){
			# I do that to inform the "future" call about the init
			linear.start = c(list(name="c"), lapply(object$coefficients, signif))
			call_clear$linear.start = do.call("call", linear.start)
		}

		if(object$family == "negbin"){
			if(is.null(dots$theta.init)){
				theta.init = object$theta.init
				call_clear$theta.init = theta.init
			}
		}

		call_clear$clusterFromUpdate = clusterFromUpdate
		call_clear$clusterStart = as.name("clusterStart")
	}

	# The variable "clusterStart" must be evaluated here!
	res = eval(call_clear, list(clusterStart=clusterStart), parent.frame())

	res
}


#' Extract the formula of a femlm fit
#'
#' This function extracts the formula from a \code{\link[FENmlm]{femlm}} estimation. If the estimation was done with fixed-effects, they are added in the formula after a pipe (\dQuote{|}). If the estimation was done with a non linear in parameters part, then this will be added in the formula in between \code{I()}.
#'
#' @inheritParams nobs.femlm
#'
#' @param x An object of class \code{femlm}. Typically the result of a \code{\link[FENmlm]{femlm}} estimation.
#' @param type A character scalar. Default is \code{type = "full"} which gives back a formula containing the linear part of the model along with the clusters (if any) and the non-linear in parameters part (if any). If \code{type = "linear"} then only the linear formula is returned. If \code{type = "NL"} then only the non linear in parameters part is returned.
#' @param ... Not currently used.
#'
#' @return
#' It returns a formula.
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{model.matrix.femlm}}, \code{\link[FENmlm]{update.femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{vcov.femlm}}.
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # simple estimation on iris data, clustering by "Species"
#' res = femlm(Sepal.Length ~ Sepal.Width + Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' # formula with the cluster variable
#' formula(res)
#' # linear part without the cluster variable
#' formula(res, "linear")
#'
#'
formula.femlm = function(x, type = c("full", "linear", "NL"), ...){
	# Extract the formula from the object
	# we add the clusters in the formula if needed

	type = match.arg(type)

	if(type == "linear"){
		return(x$fml)
	} else if(type == "NL"){
		NL.fml = x$NL.fml
		if(is.null(NL.fml)){
			stop("There was no nonlinear part estimated, option type = 'NL' cannot be used.")
		}

		return(NL.fml)
	}

	res = x$fml

	NL.fml = x$NL.fml
	if(!is.null(NL.fml)){
		fml_char = as.character(res)
		nl_char = as.character(NL.fml)
		res = as.formula(paste(fml_char[2], "~", fml_char[3], "+ I(", nl_char[2], ")"))
	}

	clusterNames = x$clusterNames
	if(length(clusterNames) > 0){
		fml_char = as.character(res)
		res = as.formula(paste(fml_char[2], "~", fml_char[3], "|", paste0(clusterNames, collapse = "+")))
	}

	res
}


#' Design matrix of a femlm model
#'
#' This function creates a design matrix of the linear part of a \code{\link[FENmlm]{femlm}} estimation. Note that it is only the linear part and the cluster variables (which can be considered as factors) are excluded from the matrix.
#'
#' @method model.matrix femlm
#'
#' @inheritParams nobs.femlm
#'
#' @param data If missing (default) then the original data is obtained by evaluating the \code{call}. Otherwise, it should be a \code{data.frame}.
#' @param ... Not currently used.
#'
#' @return
#' It returns a design matrix.
#'
#' @seealso
#' \code{\link[FENmlm]{femlm}}, \code{\link[FENmlm]{formula.femlm}}, \code{\link[FENmlm]{update.femlm}}, \code{\link[FENmlm]{summary.femlm}}, \code{\link[FENmlm]{vcov.femlm}}.
#'
#'
#' @author
#' Laurent Berge
#'
#' @examples
#'
#' # simple estimation on iris data, clustering by "Species"
#' res = femlm(Sepal.Length ~ Sepal.Width*Petal.Length +
#'             Petal.Width | Species, iris)
#'
#' head(model.matrix(res))
#'
#'
#'
model.matrix.femlm = function(object, data, ...){
	# We evaluate the formula with the past call

	# I) we obtain the right formula
	fml = object$fml

	# we kick out the intercept if there is presence of clusters
	if(attr(terms(fml), "intercept") == 1 && !is.null(object$clusterNames)){
		fml = update(fml, . ~ . - 1)
	}

	# II) evaluation with the data
	if(missing(data)){
		call_old = object$call

		data = NULL
		try(data <- eval(object$call$data, parent.frame()))

		if(is.null(data)){
			dataName = deparse(object$call$data)
			stop("To apply 'model.matrix.femlm', we fetch the original database in the parent.frame -- but it doesn't seem to be there anymore (btw it was ", dataName, ").")
		}

	}

	# control of the data
	if(is.matrix(data)){
		if(is.null(colnames(data))){
			stop("If argument data is to be a matrix, its columns must be named.")
		}
		data = as.data.frame(data)
	}
	# The conversion of the data (due to data.table)
	if(!"data.frame" %in% class(data)){
		stop("The argument 'data' must be a data.frame or a matrix.")
	}

	data = as.data.frame(data)

	res = model.matrix(fml, data)

	res
}



#### .................. ####
#### DOCUMENTATION DATA ####
####



#' Trade data sample
#'
#' This data reports trade information between countries of the European Union (EU15).
#'
#' @usage
#' data(trade)
#'
#' @format
#' \code{trade} is a data frame with 38,325 observations and 6 variables named \code{Destination}, \code{Origin}, \code{Product}, \code{Year}, \code{dist_km} and \code{Euros}.
#'
#' \itemize{
#' \item{Origin: 2-digits codes of the countries of origin of the trade flow.}
#' \item{Destination: 2-digits codes of the countries of destination of the trade flow.}
#' \item{Products: Number representing the product categories (from 1 to 20).}
#' \item{Year: Years from 2007 to 2016}
#' \item{dist_km: Geographic distance in km between the centers of the countries of origin and destination.}
#' \item{Euros: The total amount of trade flow in million euros for the specific year/product category/origin-destination country pair.}
#'
#' }
#'
#' @source
#' This data has been extrated from Eurostat on October 2017.
#'
#'
#'
"trade"


































































