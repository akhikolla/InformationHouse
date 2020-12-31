#' Generic autoplot function
#'
#' The following is imported and then re-exported to avoid conflicts with \code{ggplot2}
#' @importFrom ggplot2 autoplot labs
#' @name autoplot
#' @export
NULL

#' Plots for emfrail objects using \code{ggplot2}
#' @importFrom ggplot2 autoplot theme_minimal ggplot geom_step geom_path geom_point geom_histogram geom_abline geom_errorbar aes_string ylim xlab ylab scale_colour_manual aes_string geom_hline scale_x_continuous
#' @param object \code{emfrail} object, typically result of \code{emfrail()}
#' @param type One (or more) of \code{hist} for a histogram of the estimated frailty values,
#'  \code{hr} for a plot of the conditional and marginal hazard ratio between two cases,
#'  \code{pred} for the predicted conditional and marginal cumulative hazard or survival for one case,
#'  \code{frail} for a caterpillar plot of the ordered frailty estimates with confidence intervals, where available.
#' @param newdata A \code{data.frame} with values of the covariates. For \code{type == "hr"} the hazard ratio
#' between the first two rows of \code{newdata} is calculated. For \code{type == "pred"} the prediction
#' for the first row of \code{newdata} is calculated.
#' @param lp A numeric vector of values of the linear predictor, each corresponding to a case. For \code{type == "hr"} the hazard ratio
#' between the first two values of \code{lp} is calculated. For \code{type == "pred"} the prediction
#' for the first value of \code{lp} is calculated.
#' @param strata The name of the strata (if applicable) for which the prediction should be made.
#' @param quantity One of \code{c("cumhaz", "survival")} for \code{type == "pred"}; see \code{quantity} in \code{\link{predict.emfrail}}
#' @param type_pred One of \code{c("conditional", "marginal")} for \code{type == "pred"}; see \code{type} in \code{\link{predict.emfrail}}
#' @param conf_int One of \code{c("regular", "adjusted")} for \code{type == "pred"}; see \code{conf_int} in \code{\link{predict.emfrail}}
#' @param conf_level The width of the confidence interval for \code{type == "pred"}; see \code{conf_level} in \code{\link{predict.emfrail}}
#' @param individual Logical, for \code{type == "pred"} to be used for drawing a curve when the rows of \code{newdata} refer to the same individual; see
#' \code{individual} in \code{\link{predict.emfrail}}
#' @param ... Further arguments to be passed on to `ggplot` (ignored)
#'
#' @return A list of \code{ggplot2} objects corresponding to the required plots, or one \code{ggplot2} if only one plot is selected
#'
#' For the catterpillar plot, in the case of the gamma frailty model, the vertical lines represent the 0.025 and 0.975 quantiles of the
#' posterior gamma distribution. For other distributions, this quantity is not easy to calculate (at least not in closed form)
#' and only the frailty estimates are shown.
#'
#' @note It's normal for \code{autoplot} to give a warning of the type \code{Warning: Ignoring unknown aesthetics: id
#' }. This is because, in \code{ggplot2} terms, the \code{id} aesthetic is not recognized. This is correct, and for any
#' practical purpose this will not make a difference (you can safely ignore the warnings). However, this makes it
#' easier to create an interactive plot out of the resulting object.
#'
#'
#' @export
#' @seealso \code{\link{predict.emfrail}}, \code{\link{summary.emfrail}}, \code{\link{plot.emfrail}}.
#'
#' @examples
#' mod_rec <- emfrail(Surv(start, stop, status) ~ treatment + number + cluster(id), bladder1,
#' control = emfrail_control(ca_test = FALSE, lik_ci = FALSE))
#'
#' # Histogram of the estimated frailties
#' autoplot(mod_rec, type = "hist")
#'
#' # Ordered estimated frailties (with confidence intervals, for gamma distribution)
#' autoplot(mod_rec, type = "frail")
#'
#' # hazard ratio between placebo and pyridoxine
#' newdata1 <- data.frame(treatment = c("placebo", "pyridoxine"),
#'                        number = c(1, 3))
#'
#' autoplot(mod_rec, type = "hr", newdata = newdata1)
#'
#' # predicted cumulative hazard for placebo, and number = 1
#' autoplot(mod_rec, type = "pred", newdata = newdata1[1,])
#'
#' # predicted survival for placebo, and number = 1
#' autoplot(mod_rec, type = "pred", quantity = "survival", newdata = newdata1[1,])
#'
#' # predicted survival for an individual that switches from
#' # placebo to pyridoxine at time = 15
#' \dontrun{
#' newdata2 <- data.frame(treatment = c("placebo", "pyridoxine"),
#'                        number = c(1, 3),
#'                        tstart = c(0, 15),
#'                        tstop = c(15, Inf))
#'
#' autoplot(mod_rec, type = "pred", quantity = "survival", newdata = newdata2, individual = TRUE)
#' }
autoplot.emfrail <- function(object,
                             type = c("hist", "hr", "pred", "frail"),
                             newdata = NULL, lp = NULL,
                             strata = NULL,
                             quantity = "cumhaz",
                             type_pred = c("conditional", "marginal"),
                             conf_int = "adjusted",
                             conf_level = 0.95,
                             individual = FALSE,
                             ...) {

  # browser()
  if(any(!(type %in% c("hist", "hr", "pred", "frail"))))
    stop("type misspecified, check ?autoplot.emfrail")

  res <- vector("list", length(type))
  i <- 1

  if("hist" %in% type) {
    sobj <- summary.emfrail(object)
    plot1 <- sobj$frail %>%
      ggplot(aes_string(x = "z")) +
      geom_histogram() +
      labs(x = "Frailties", y = "Frequency")
    res[[i]] <- plot1
    i <- i+1
  }

  if("hr" %in% type) {


    if(is.null(newdata))
      if(missing(lp)) stop("lp or newdata should be specified") else
        if(length(lp) < 2) stop("lp must be at least of length 2 or newdata must have 2 rows with type = hr") else {
          if(length(lp) > 2) warning("just the first 2 lp values are used for type = hr")
          if(lp[1] == lp[2]) stop("values of lp should be different")
        }


    if(!is.null(newdata))
      if(!inherits(newdata, "data.frame")) stop("newdata must be a data.frame") else
        if(nrow(newdata) < 2) stop("newdata must have 2 rows with type = hr") else {
          if(nrow(newdata) > 2) warning("just the first 2 newdata rows are used for type = hr")
          if(isTRUE(all.equal(newdata[1,], newdata[2,]))) stop("the rows of newdata must be different")
        }


    p <- predict.emfrail(object,
                         lp = lp,
                         newdata = newdata,
                         strata = strata,
                         quantity = "cumhaz",
                         conf_int = NULL,
                         conf_level = conf_level)


    hr_cond <- p[[1]]$cumhaz / p[[2]]$cumhaz

    hr_mar <- p[[1]]$cumhaz_m / p[[2]]$cumhaz_m

    hr_cond[1] <- hr_cond[2]  # that's cause it's 0/0
    # the hr_mar in the beginning is the same for gamma / pvf
    # except with the PS. then it is the same all the time
    # this part here is simply a cosmetic change to get nice plots
    if(object$distribution$dist == "stable")
      hr_mar[1] <- hr_mar[2] else
        hr_mar[1] <- hr_cond[1]


    plot2 <- p[[1]] %>%
      ggplot(aes_string(x = "time")) +
      geom_hline(aes_string("yintercept" = 1), col = "gray") +
      geom_step(aes_string(y = "hr_cond", col = shQuote("1"))) +
      geom_step(aes_string(y = "hr_mar", col = shQuote("2"))) +
      theme_minimal() +
      scale_colour_manual(name = "type", values = c("#fc8d59", "#000000" ), labels = c("conditional", "marginal") ) +
      ylab("hazard ratio")

    res[[i]] <- plot2
    i <- i+1
  }

  if("pred" %in% type) {


    if(is.null(newdata) & is.null(lp)) stop("newdata or lp have to be specified for type = pred")

    if(!is.null(newdata))
      if(!inherits(newdata, "data.frame")) stop("newdata must be a data.frame") else
        if(nrow(newdata) < 1) stop("newdata must have eat least 1 row for type = pred") else
          if(nrow(newdata) > 1 & !isTRUE(individual)) {
            warning("just the first row is used for type = pred")
            newdata <- newdata[1,,drop = FALSE]
          }

    if(is.null(newdata) & !is.null(lp))
      if(length(lp) < 1) stop("lp should have at least length 1") else
        if(length(lp) > 1) {
          warning("just the first value of lp is used for type = pred")
          lp <- lp[1]
        }

    if(!(quantity %in% c("cumhaz", "survival")) | length(quantity) != 1)
      stop("quantity must be cumhaz or survival for type = pred")

    if("both" %in% type_pred)
      type_pred <- c("conditional", "marginal")
    # browser()
    p1 <- predict.emfrail(object,
                          lp = lp,
                          strata = strata,
                          newdata = newdata,
                          quantity = quantity,
                          type = type_pred,
                          conf_int = conf_int,
                          individual = individual)


    if("conditional" %in% type_pred & "marginal" %in% type_pred) type_pred <- "both"

    if(type_pred == "conditional") {
      if(quantity == "cumhaz") {
        # if(conf_int == "adjusted") ylim <- c(0, max(p1$cumhaz_r_a)) else
        #   if(conf_int == "regular") ylim <- c(0, max(p1$cumhaz_r)) else
        #     ylim <- c(0, max(p1$cumhaz_m))

        plot1 <- p1 %>%
          ggplot(aes_string(x = "time", y = "cumhaz")) +
          geom_step()

        if(conf_int == "adjusted") {
          plot1 <-
            plot1 +
            geom_path(aes_string(y = "cumhaz_l_a"), lty = 2) +
            geom_path(aes_string(y = "cumhaz_r_a"), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <-
            plot1 +
            geom_path(aes_string(y = "cumhaz_l"), lty = 2) +
            geom_path(aes_string(y = "cumhaz_r"), lty = 2)
        }
      }

      if(quantity == "survival") {
        plot1 <- p1 %>%
          ggplot(aes_string(x = "time", y = "survival")) +
          geom_step() +
          ylim(c(0,1))

        if(conf_int == "adjusted") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_l_a"), lty = 2) +
            geom_path(aes_string(y = "survival_r_a"), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_l"), lty = 2) +
            geom_path(aes_string(y = "survival_r"), lty = 2)
        }
      }
    }

    if(type_pred == "marginal") {
      if(quantity == "cumhaz") {

        plot1 <- p1 %>%
          ggplot(aes_string(x = "time", y = "cumhaz_m")) +
          geom_step()

        if(conf_int == "adjusted") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "cumhaz_m_l_a"), lty = 2) +
            geom_path(aes_string(y = "cumhaz_m_r_a"), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "cumhaz_m_l"), lty = 2) +
            geom_path(aes_string(y = "cumhaz_m_r"), lty = 2)
        }
      }

      if(quantity == "survival") {
        plot1 <- p1 %>%
          ggplot(aes_string(x = "time", y = "survival_m")) +
          geom_step() +
          ylim(c(0,1)) +
          ylab("S(t)")

        if(conf_int == "adjusted") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_m_l_a"), lty = 2) +
            geom_path(aes_string(y = "survival_m_r_a"), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_m_l"), lty = 2) +
            geom_path(aes_string(y = "survival_m_r"), lty = 2)
        }
      }


    }

    if(type_pred == "both") {
      if(quantity == "cumhaz") {

        plot1 <- p1 %>%
          ggplot(aes_string(x = "time")) +
          geom_step(aes_string(y = "cumhaz", col = shQuote("1"))) +
          geom_step(aes_string(y = "cumhaz_m", col = shQuote("2"))) +
          theme_minimal() +
          ylab("cumulative hazard") +
          scale_colour_manual(name = "type", values = c("#fc8d59", "#000000" ), labels = c("conditional", "marginal") )

        if(conf_int == "adjusted") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "cumhaz_m_l_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_m_r_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_l_a", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_r_a", col = shQuote("1")), lty = 2)

          # with(p1, lines(time, cumhaz_m_l_a, lty = 2, col = 2))
          # with(p1, lines(time, cumhaz_m_r_a, lty = 2, col = 2))
          # with(p1, lines(time, cumhaz_l_a, lty = 2))
          # with(p1, lines(time, cumhaz_r_a, lty = 2))
        }

        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "cumhaz_m_l", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_m_r", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_l", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "cumhaz_r", col = shQuote("1")), lty = 2)

        }

      }

      if(quantity == "survival") {
        plot1 <- p1 %>%
          ggplot(aes_string(x = "time")) +
          geom_step(aes_string(y = "survival", col = shQuote("1"))) +
          geom_step(aes_string(y = "survival_m", col = shQuote("2"))) +
          ylim(c(0,1)) +
          ylab("S(t)") +
          theme_minimal() +
          scale_colour_manual(name = "type", values = c("#fc8d59", "#000000" ), labels = c("conditional", "marginal") )

        # with(p1, plot(time, survival_m,
        #               type = "s",
        #               #main = "Survival",
        #               ylab = "S(t)",
        #               xlab = "time",
        #               col = 2,
        #               ylim = c(0,1),
        #               ...))
        # with(p1, lines(time, survival,
        #                type = "s",
        #                col = 1,
        #                ...))
        #
        if(conf_int == "adjusted") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_m_l_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_m_r_a", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_l_a", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "survival_r_a", col = shQuote("1")), lty = 2)
        }
        if(conf_int == "regular") {
          plot1 <- plot1 +
            geom_path(aes_string(y = "survival_m_l", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_m_r", col = shQuote("2")), lty = 2) +
            geom_path(aes_string(y = "survival_l", col = shQuote("1")), lty = 2) +
            geom_path(aes_string(y = "survival_r", col = shQuote("1")), lty = 2)
        }

        # legend(x = "topright", legend = c("conditional", "marginal"), col = 1:2, lty = 1)
      }


    }

    res[[i]] <- plot1
    i <- i+1

  }

  if("frail" %in% type) {

    sobj <- summary.emfrail(object)

    frdat <- sobj$frail[order(sobj$frail$z),]

    # it's normal that id gives a warning, that's just there
    # in case you want to make an interactive plot with ggplotly

    plot1 <- frdat %>%
      ggplot(aes_string(x = seq_along(frdat$z), y = "z")) +
      suppressWarnings(geom_point(aes_string(id = "id"))) +
      scale_x_continuous(labels = as.character(frdat$id), breaks = seq_along(frdat$z)) +
      xlab("cluster")

    if(sobj$est_dist$dist == "gamma")
      plot1 <- plot1 +
      suppressWarnings(geom_errorbar(aes_string(ymin = "lower_q", ymax = "upper_q", id = "id")))

    res[[i]] <- plot1
  }

  if(length(res) == 1) res <- res[[1]]

  res
}
