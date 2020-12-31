#' Plot rlme Fit
#' 
#' Generates Normal Q-Q plot of residuals from rlme fit
#' 
#' 
#' @param x A list of class rlme. Store as fit.rlme.
#' @param \dots not used
#' @examples
#' 
#' data(schools)
#' rlme.fit = rlme(y ~ 1 + sex + age + (1 | region) + (1 | region:school), schools, method="gr")
#' plot(rlme.fit)
#'
#' @export
plot.rlme <- function(x, ...) {
    if (x$method == "GR") {
        getgrstplot(x)
    }
    else if (x$method == "REML") {
        getlmestplot(x)
    }
    else {
        cat("rlme only supports plotting fits from GR and REML methods\n")
    }
    
    invisible()
}



#' Rank-based Estimates for Mixed-Effects Nested Models
#' 
#' This function estimates fixed effects and predicts random effects in two-
#' and three-level random effects nested models using three rank-based fittings
#' (GR, GEER, JR) via the prediction method algorithm RPP.
#' 
#' The iterative methods GR and GEER can be quite slow for large datasets; try
#' JR for faster analysis. If you want to use the GR method, try using
#' rprpair='med-mad'. This method avoids building a NxN covariance matrix which
#' can quickly become unwieldly with large data.
#' 
#' 
#' @param f An object of class formula describing the mixed effects model. The
#' syntax is same as in the lme4 package. Example: y ~ 1 + sex + age + (1 |
#' region) + (1 | region:school) - sex and age are the fixed effects, region
#' and school are the nested random effects, school is nested within region.
#' @param data The dataframe to analyze. Data should be cleaned prior to
#' analysis: cluster and subcluster columns are expected to be integers and in
#' order (e.g. all clusters and subclusters )
#' @param method string indicating the method to use (one of "gr", "jr",
#' "reml", and "geer"). defaults to "gr".
#' @param print Whether or not to print a summary of results. Defaults to
#' false.
#' @param na.omit Whether or not to omit rows containing NA values. Defaults to
#' true.
#' @param weight When weight="hbr", it uses hbr weights in GEE weights. By
#' default, ="wil", it uses Wilcoxon weights. See the theory in the references.
#' @param rprpair By default, it uses "hl-disp" in the random prediction
#' procedure (RPP). Also, "med-mad" would be an alternative.
#' @param verbose Boolean indicating whether to print out diagnostic messages.
#' 
#' @return The function returns a list of class "rlme". Use summary.rlme to see
#' a summary of the fit.
#' 
#'   \item{formula}{ The model formula. }
#'   \item{method}{ The method used. }
#'   \item{fixed.effects}{ Estimate of fixed effects. }
#'   \item{random.effects}{Estimate of random effects. }
#'   \item{standard.residual}{ Residuals. }
#'   \item{intra.class.correlations}{ Intra/inter-class correlationa estimates obtained from RPP. }
#'   \item{t.value}{ t-values. }
#'   \item{p.value}{ p-values. }
#'   \item{location}{ Location.  }
#'   \item{scale}{ Scale.  }
#'   \item{y}{ The response variable y.  }
#'   \item{num.obs}{ Number of observations in provided dataset.}
#'   \item{num.clusters}{ The number of clusters. }
#'   \item{num.subclusters}{ The number of subclusters. }
#'   \item{effect.err}{ Effect from error. }
#'   \item{effect.cluster}{ Effect from cluster. }
#'   \item{effect.subcluster}{Effect from subcluster. }
#'   \item{var.b}{ Variances of fixed effects estimate (Beta estimates). }
#'   \item{xstar}{ Weighted design matrix with error covariance matrix. }
#'   \item{ystar}{ Weighted response vector with its covariance matrix.}
#'   \item{ehat}{ The raw residual.  }
#'   \item{ehats}{ The raw residual after weighted step. Scaled residual. }
#'   
#' @author Yusuf Bilgic \email{yekabe@@hotmail.com} and Herb Susmann \email{hps1@@geneseo.edu}
#' 
#' @seealso summary.rlme, plot.rlme, compare.fits
#' 
#' @references Y. K. Bilgic. Rank-based estimation and prediction for mixed
#' effects models in nested designs. 2012. URL
#' http://scholarworks.wmich.edu/dissertations/40. Dissertation.
#' 
#' T. P. Hettmansperger and J. W. McKean. Robust Nonparametric Statistical
#' Methods. Chapman Hall, 2012.
#' @keywords models
#' 
#' @examples
#' 
#' data(schools)
#' 
#' rlme.fit = rlme(y ~ 1 + sex + age + (1 | region) + (1 | region:school), schools, method="gr")
#' summary(rlme.fit)
#' 
#' # Try method="geer", "reml", "ml" and "jr" along with 
#' # rprpair="hl-disp" (not robust), and "med-mad" (robust),
#' # weight="hbr" is for the gee method.
#' 
#' @importFrom utils tail
#' @importFrom stringr str_extract_all
#' @export
rlme <- function(f, data, method = "gr", print = FALSE, na.omit = TRUE,  weight = "wil", rprpair = "hl-disp", verbose=FALSE)  {
    location = 2
    method = tolower(method)
    fit = list(formula = f, location = location, scale = scale)
    if (na.omit == TRUE) {
        data = na.omit(data)
    }
    response_name = as.character(as.formula(f)[[2]])
    random_terms = str_extract_all(as.character(f)[[3]], "\\(([0-9 a-zA-Z.:|]+)\\)")
    random_terms = lapply(random_terms[[1]], function(str) substr(str, 
        2, nchar(str) - 1))
    covariate_names = setdiff(attributes(terms(f))$term.labels, 
        random_terms)
    school_name = section_name = ""
    levels = 0
    if (length(random_terms) == 2) {
        levels = 3
        school_name = tail(strsplit(random_terms[[1]], "\\W")[[1]], 
            1)
        section_name = tail(strsplit(random_terms[[2]], "\\W")[[1]], 
            1)
        I = length(unique(factor(data[[school_name]])))
        sec = as.vector(sec_vec(data[[school_name]], data[[section_name]]))
        
        #ss = rhosch(data[[school_name]], data[[school_name]], 
        #    data[[section_name]])$nvec
        
        
        fit$num.clusters = I
        fit$num.subclusters = sum(sec)
        mat = mat_vec(data[[school_name]], data[[section_name]])
        
        ss = rowSums(mat)
        
        J = sum(sec)
        n = sum(mat)
        one = rep(1, length(data[[response_name]]))
        schsize = aggregate(one, by = list(data[[response_name]]), 
            FUN = sum)[, 2]
        school1 = factor(rep(1:length(schsize), schsize))
        sss = aggregate(one, by = list(data[[school_name]], data[[section_name]]), 
            FUN = sum)[, 3]
        section = factor(rep(1:length(ss), ss))
    }
    if (length(random_terms) == 1) {
        levels = 2
        school_name = tail(strsplit(random_terms[[1]], "\\W")[[1]], 
            1)
        I = length(unique(factor(data[[school_name]])))
        mat = mat_vec(data[[school_name]], rep(1, length(data[[school_name]])))
        n = sum(mat)
    }
    x = as.matrix(data[covariate_names])
    x = apply(x, 2, function(x) {
        x - mean(x)
    })
    y = data[[response_name]]
    fit$num.obs = length(y)
    fit$y = y
    if (length(random_terms) == 0) {
        cat("You have entered an independent linear model.\n")
        cat("The function 'lmr' can be used to fit these models.\n")
        cat("Continuing using lmr.\n")
        fitw = lmr(f, data=data)
        return(fitw)
    }
    if (method == "reml" || method == "ml") {
        if (levels == 3) {
            REML = LM_est(x, y, dat = data.frame(y = y, x, school = data[[school_name]], 
                section = data[[section_name]]), method = toupper(method))
        }
        else if (levels == 2) {
            REML = LM_est2(x, y, dat = data.frame(y = y, x, school = data[[school_name]]), 
                method = toupper(method))
        }
        REMLb = REML$theta
        REMLs = REML$sigma
        REMLe = REML$ses
        REML$ehat
        tvalue = REMLb/REMLe
        pvalue = 2 * pnorm(-abs(tvalue))
        intracoeffs = c(REMLs[1]/sum(REMLs), (REMLs[1] + REMLs[2])/sum(REMLs))
        fit$method = "REML"
        fit$ehat = REML$ehat
        fit$effect.err = REML$effect_err
        fit$effect.cluster = REML$effect_sch
        if (levels == 3) {
            fit$effect.subcluster = REML$effect_sec
        }
        fit$fixed.effects = data.frame(RowNames = c("(Intercept)", 
            covariate_names), Estimate = REMLb, StdError = REMLe, 
            tvalue = tvalue, pvalue = pvalue)
        fit$var.b = REML$varb
        fit$intra.class.correlations = intracoeffs
        fit$t.value = tvalue
        fit$p.value = pvalue
        fit$standr.lme = REML$standr.lme
    }
    else if (method == "jr") {
        if (levels == 3) {
            JR = JR_est(x, y, I, sec, mat, data[[school_name]], 
                data[[section_name]], rprpair = rprpair, verbose=verbose)
        }
        else if (levels == 2) {
            JR = JR_est2(x, y, I, 1, mat, data[[school_name]], 
                rep(1, length(data[[school_name]])), rprpair = rprpair, verbose=verbose)
        }
        JRb = JR$theta
        JRe = JR$ses
        JRs = JR$sigma
        tvalue = JRb/JRe
        pvalue = 2 * pnorm(-abs(tvalue))
        intracoeffs = c(JRs[1]/sum(JRs), (JRs[1] + JRs[2])/sum(JRs))
        fit$method = "JR"
        fit$ehat = JR$ehat
        fit$fixed.effects = data.frame(RowNames = c("(Intercept)", 
            covariate_names), Estimate = JRb, StdError = JRe, 
            tvalue = tvalue, pvalue = pvalue)
        fit$effect.err = JR$effect_err
        fit$effect.cluster = JR$effect_sch
        if (levels == 3) {
            fit$random.effects = data.frame(Groups = c(paste(school_name, 
                ":", section_name, sep = ""), school_name, "Residual"), 
                Name = c("(Intercept)", "(Intercept)", ""), Variance = JRs)
            fit$effect.subcluster = JR$effect_sec
        }
        else if (levels == 2) {
            fit$random.effects = data.frame(Groups = c(school_name, 
                "Residual"), Name = c("(Intercept)", ""), Variance = JRs)
        }
        fit$intra.class.correlations = intracoeffs
        fit$var.b = JR$varb
        fit$t.value = tvalue
        fit$p.value = pvalue
    }
    else if (method == "gr") {
        if (levels == 3) {
            GR = GR_est(x, y, I, sec, mat, data[[school_name]], 
                data[[section_name]], rprpair = rprpair, verbose=verbose)
        }
        else if (levels == 2) {
            GR = GR_est2(x, y, I, rep(1, length(data[[school_name]])), 
                mat, data[[school_name]], rep(1, length(data[[school_name]])), 
                rprpair = rprpair, verbose=verbose)
        }
        GRb = GR$theta
        GRe = GR$ses
        tvalue = GRb/GRe
        pvalue = 2 * pnorm(-abs(tvalue))
        varb = GR$varb
        GRs = GR$sigma
        intracoeffs = c(GRs[1]/sum(GRs), (GRs[1] + GRs[2])/sum(GRs))
        coll.stres <- stanresidgr(GR$xstar, GR$ystar, resid = GR$ehats, 
            delta = 0.8, param = 2, conf = 0.95)
        stresgr <- coll.stres$stanr
        fit$method = "GR"
        fit$ehat = GR$ehat
        fit$ehats = GR$ehats
        fit$xstar = GR$xstar
        fit$ystar = GR$ystar
        fit$fixed.effects = data.frame(RowNames = c("(Intercept)", 
            covariate_names), Estimate = GRb, StdError = GRe, 
            tvalue = tvalue, pvalue = pvalue)
        fit$effect.err = GR$effect_err
        fit$effect.cluster = GR$effect_sch
        if (levels == 3) {
            fit$random.effects = data.frame(Groups = c(school_name, paste(school_name, 
                ":", section_name, sep = ""), "Residual"), 
                Name = c("(Intercept)", "(Intercept)", ""), Variance = GRs)
            fit$effect.subcluster = GR$effect_sec
        }
        else if (levels == 2) {
            fit$random.effects = data.frame(Groups = c(school_name, 
                "Residual"), Name = c("(Intercept)", ""), Variance = GRs)
        }
        fit$standard.residual = stresgr
        fit$intra.class.correlations = intracoeffs
        fit$var.b = varb
        fit$t.value = tvalue
        fit$p.value = pvalue
    }
    else if (method == "geer") {
        tol = 1.1e-25
        if (levels == 3) {
            GEER = GEER_est(x, y, I, sec, mat, data[[school_name]], 
                data[[section_name]], weight = weight, rprpair = rprpair, verbose=verbose)
        }
        else if (levels == 2) {
            GEER = GEER_est2(x, y, I, rep(1, length(data[[school_name]])), 
                mat, data[[school_name]], rep(1, length(data[[school_name]])), 
                weight = weight, rprpair = rprpair, verbose=verbose)
        }
        GEERb = GEER$theta
        GEERs = GEER$sigma
        GEER_e2 = GEER$ses_AP
        intracoeffs = c(GEERs[1]/sum(GEERs), (GEERs[1] + GEERs[2])/sum(GEERs))
        GEERe = GEER$ses_AP
        tvalue = GEERb/GEERe
        pvalue = 2 * pnorm(-abs(tvalue))
        fit$method = "GEER"
        fit$ehat = GEER$ehat
        fit$fixed.effects = data.frame(RowNames = c("(Intercept)", 
            covariate_names), Estimate = GEERb, StdError = GEERe, 
            tvalue = tvalue, pvalue = pvalue)
        fit$effect.err = GEER$effect_err
        fit$effect.cluster = GEER$effect_sch
        if (levels == 3) {
            fit$random.effects = data.frame(Groups = c(school_name, paste(school_name, 
                ":", section_name, sep = ""), "Residual"), 
                Name = c("(Intercept)", "(Intercept)", ""), Variance = GEERs)
            fit$effect.subcluster = GEER$effect_sec
        }
        else if (levels == 2) {
            fit$random.effects = data.frame(Groups = c(school_name, 
                "Residual"), Name = c("(Intercept)", ""), Variance = GEERs)
        }
        fit$intra.class.correlations = intracoeffs
        fit$var.b = GEER$varb
        fit$t.value = tvalue
        fit$p.value = pvalue
    }
    class(fit) = "rlme"
    if (print == TRUE) {
        summary(fit)
    }
    return(fit)
}

#' rlme Summary
#' 
#' Summarizes a model fit from the rmle function
#' 
#' 
#' @param object A list of class rlme
#' @param \dots not used
#' @author Herb Susmann \email{hps1@@geneseo.edu}
#' @seealso \code{\link{rlme}} \code{\link{plot.rlme}}
#' @export
summary.rlme <- function(object, ...) {
    fit = object
    cat("Linear mixed model fit by ", fit$method, "\n")
    cat("Formula: ", deparse(fit$formula), "\n")
    
    if ("random.effects" %in% attributes(fit)$names) {
        cat("Random effects:\n")
        random.effects = fit$random.effects
        names(random.effects) = c("Groups", "Name", "Variance")
        print(random.effects, row.names = FALSE, right = FALSE)
    }
    
    cat("\nNumber of obs:\n")
    cat(fit$num.obs)
    if ("num.clusters" %in% attributes(fit)$names && "num.subclusters" %in% 
        attributes(fit)$names) {
        cat(" observations,", fit$num.clusters, "clusters,", 
            fit$num.subclusters, "subclusters")
    }
    cat("\n")
    cat("\nFixed effects:\n")
    fixed.effects = fit$fixed.effects
    
    if(length(fixed.effects) == 5) {
      names(fixed.effects) = c("", "Estimate", "Std. Error", "t value", 
          "p value")
    } else if(length(fixed.effects) == 3) {
      names(fixed.effects) = c("", "Estimate", "Std. Error")
    } else if(length(fixed.effects) == 2) {
      names(fixed.effects) = c("", "Estimate")
    }
    
    print(fixed.effects, row.names = FALSE, right = FALSE)
    if ("intra.class.correlations" %in% attributes(fit)$names) {
        cat("\nIntra-class correlation coefficients\n")
        intra.coeffs = data.frame(names = c("intra-cluster", 
            "intra-subcluster"), Estimates = fit$intra.class.correlations)
        names(intra.coeffs) = c("", "Esimates")
        print(intra.coeffs, row.names = FALSE, right = FALSE)
    }
    
    if(!is.null(fit$var.b)) {
      cat("\ncov-var (fixed effects)\n")
      print(fit$var.b)
    }
}