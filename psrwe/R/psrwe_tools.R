## --------------------------------------------------------------------
## --------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## --------------------------------------------------------------------
## --------------------------------------------------------------------

##
## RWE related class names
##
get_rwe_class <- function(c.str) {
    switch(c.str,
           DWITHPS   = "RWE_DWITHPS",
           PSDIST    = "RWE_PSDIST",
           D_GPS     = "RWE_D_GPS",
           GPSDIST   = "RWE_GPSDIST",
           CLRST     = "RWE_CLRST",
           PPRST     = "RWE_POWERPRST"
           )
}


## compute propensity scores
get_ps <- function(dta, ps_method = c("logistic", "randomforest"),
                   ps_fml, ntree = 5000,
                   ..., grp = NULL,
                   ps_cov = NULL) {

    type <- match.arg(ps_method)

    ## generate formula
    if (is.null(ps_fml))
        ps_fml <- as.formula(paste(grp, "~",
                                   paste(ps_cov, collapse = "+"),
                                   sep = ""))

    ## identify grp if passed from formula
    grp <- all.vars(ps_fml)[1]

    ## fit model
    switch(type,
           logistic = {
               glm_fit <- glm(ps_fml, family = binomial, data = dta, ...)
               est_ps  <- glm_fit$fitted
           },
           randomforest = {
               dta[[grp]] <- as.factor(dta[[grp]])
               rf_fit     <- randomForest(ps_fml, data = dta,
                                          ntree = ntree, ...)
               est_ps     <- predict(rf_fit, type = "prob")[, 2]
           })
    est_ps
}


## Get number of subjects borrowed
##
## @param A target number of subjects to be borrowed
## @param m.lambda method to split A. rs: by overlapping coefficient; even: by
##     minimizing trt and control imbalance in numbers
##
## @return power parameter before standardization
##
##
##
get_aborrow <- function(total_borrow, ns0, rs = NULL,
                        ns1_trt = NULL, ns1_ctl = NULL,
                        m_lambda = c("dist", "inverse"), ...) {

    m_lambda <- match.arg(m_lambda);
    rst      <- switch(m_lambda,
                       dist    = {
                           apply(cbind(ns0, total_borrow * rs / sum(rs)),
                                 1, min)
                       },
                       inverse = {
                           mrs <- 1 / (1 - rs)
                           apply(cbind(ns0, total_borrow * mrs / sum(mrs)),
                                 1, min)})
    rst
}


##
## Get current stratum
##
##
get_cur_d <- function(data, i, v_outcome) {
    cur_d1 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 1,
                   v_outcome]
    cur_d0 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 0,
                   v_outcome]

    if (is.data.frame(cur_d1)) {
        flag <- 0 == nrow(cur_d1) | 0 == nrow(cur_d0)
    } else {
        flag <- 0 == length(cur_d1) | 0 == length(cur_d0)
    }

    if (flag) {
        warning(paste("Stratum ", i,
                   " contains no subjects from group 1 or 0",
                   sep = ""))
    }

    list(cur_d1 = cur_d1,
         cur_d0 = cur_d0)
}


## Generate frequency table for factor columns
##
## @return a vector with the number of samples in group 0, the number of samples
##     in group 1, and the KL divergence from group 0 to group 1
##
get_freq_tbl <- function(data, var.groupby, vars = NULL) {

    if (is.null(vars))
        vars <- colnames(data);

    rst <- NULL;
    for (v in vars) {
        if (!is.factor(data[[v]]))
            next

        cur.freq <- data %>%
          count_(c(var.groupby, v)) %>%
          group_by_(.dots = var.groupby) %>%
          mutate(Sum  = sum(.data$n),
                 Freq = .data$n / sum(.data$n)) %>%
          mutate_if(is.factor, as.character) %>%
          mutate(Cov = v) %>%
          rename_(Value = v)

        rst <- rbind(rst,
                     data.frame(cur.freq))
    }

    rst
}



##
## plot density of propensity score
##
plot_ps <- function(data.withps, overall.inc = TRUE, add.text = TRUE,
                    facet.scales = "free_y", ...) {

    N0 <- N1 <- Dist <- Ps <- Group <- NULL

    stopifnot(inherits(data.withps,
                       what = get_rwe_class("DWITHPS")))

    pskl        <- rwe_ps_dist(data.withps, ...)
    nstrata     <- data.withps$nstrata;
    dtaps       <- data.withps$data;

    pskl$Strata <- as.factor(c(paste("Stratum ",
                                     1:nstrata, sep = ""),
                               "Overall"))

    xlim        <- range(dtaps[which(!is.na(dtaps[["_strata_"]])),"_ps_"],
                         na.rm = TRUE)

    all.data <- NULL;
    for (i in 1:nstrata) {
        cur.sub  <- dtaps[which(i == dtaps[["_strata_"]]),];
        cur.data <- data.frame(Strata = paste("Stratum ", i, sep = ""),
                               Ps     = cur.sub[["_ps_"]],
                               Group  = cur.sub[["_grp_"]]);
        all.data <- rbind(all.data, cur.data);
    }

    if (overall.inc) {
        cur.data <-  data.frame(Strata = "Overall",
                                Ps     = dtaps[["_ps_"]],
                                Group  = dtaps[["_grp_"]]);

        all.data <- rbind(all.data, cur.data);
    } else {
        pskl <- pskl %>% dplyr::filter(.data$Strata != "Overall");
    }

    all.data$Group <- as.factor(all.data$Group);
    rst <- ggplot(data = all.data, aes(x = Ps)) +
        geom_density(alpha = 0.2,
                     aes(group = Group,
                         fill  = Group,
                         linetype = Group),
                     trim  = TRUE,
                     na.rm = TRUE) +
        labs(x = "Propensity Score", y = "Density") +
        scale_y_continuous(breaks = NULL) +
        scale_x_continuous(limits = xlim) +
        scale_fill_manual(values=c("gray20", "gray80")) +
        theme_bw() +
        theme(strip.background = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black"),
              panel.spacing = unit(0, "lines")) +
        facet_grid(Strata ~ ., scales = facet.scales);

    if (add.text) {
        rst <- rst +
            geom_text(x = Inf, y = Inf, hjust = 1, vjust = 1,
                      aes(label = paste('n0=', N0,
                                        ", n1=", N1,
                                        ", Distantce=",
                                        format(Dist, digits = 3),
                                        sep = "")),
                      data = pskl, size = 4);
    }
    rst
}

plot_balance_fac <- function(dtaps, v, overall.inc = TRUE) {
    cur.d <- get_freq_tbl(dtaps,
                          var.groupby = c("Strata", "Group"),
                          vars = v)

    cur.d <- cur.d %>%
      dplyr::filter(!is.na(.data$Strata))

    cur.d$Strata <- paste("Stratum ", cur.d$Strata, sep = "")
    if (overall.inc) {
        cur.overall <- get_freq_tbl(dtaps,
                                    var.groupby = "Group",
                                    vars = v);
        cur.overall$Strata <- "Overall";
        cur.d <- rbind(cur.d, cur.overall);
    }

    cur.d$Group <- as.factor(cur.d$Group);
    cur.d$Value <- as.factor(cur.d$Value);

    rst <- ggplot(data = cur.d, aes(x = .data$Value, y = .data$Freq)) +
        geom_bar(alpha = 0.4,
                 stat = "identity",
                 position = "dodge",
                 color = "black",
                 aes(group = .data$Group,
                     fill  = .data$Group)) +
        scale_fill_manual(values = c("gray20", "gray80")) +
        scale_y_continuous(breaks = NULL, limits = c(0,1)) +
        labs(x = "", y = "") +
        facet_grid(Strata ~ .);
    rst
}

plot_balance_cont <- function(dtaps, v, nstrata,
                              overall.inc = TRUE,
                              facet.scales = "free_y") {

    Value <- Group <- NULL
    cur.d <- NULL;
    for (i in 1:nstrata) {
        cur.sub      <- dtaps[which(i == dtaps[["_strata_"]]),];
        cur.v        <- data.frame(Cov    = v,
                                   Value  = cur.sub[[v]],
                                   Group  = cur.sub[["_grp_"]]);
        cur.v$Strata <- paste("Stratum ", i, sep = "");
        cur.d        <- rbind(cur.d, cur.v);
    }

    if (overall.inc) {
        cur.sub      <- dtaps;
        cur.v        <- data.frame(Cov    = v,
                                   Value  = cur.sub[[v]],
                                   Group  = cur.sub[["_grp_"]]);
        cur.v$Strata <- paste("Overall");
        cur.d        <- rbind(cur.d, cur.v);
    }
    cur.d$Group <- as.factor(cur.d$Group);

    rst <- ggplot(data = cur.d, aes(x = Value)) +
        geom_density(alpha = 0.2,
                     aes(group = Group,
                         fill  = Group,
                         linetype = Group),
                     na.rm = TRUE) +
        scale_y_continuous(breaks = NULL) +
        scale_fill_manual(values=c("gray20", "white")) +
        labs(x = "", y = "") +
        facet_grid(Strata ~ ., scales = facet.scales);
    rst
}


plot_balance <- function(data.withps, overall.inc = TRUE, v.cov = NULL,
                         facet.scales = "free_y", label.cov = v.cov,
                         legend.width = 0.08,
                         ...) {

    if (is.null(v.cov))
        v.cov <- all.vars(data.withps$ps_fml)[-1];

    if (is.null(label.cov))
        label.cov <- v.cov;

    nstrata      <- data.withps$nstrata;
    dtaps        <- data.withps$data;
    dtaps$Strata <- dtaps[["_strata_"]];
    dtaps$Group  <- dtaps[["_grp_"]];

    rst <- list();
    for (v in v.cov) {
        if (is.factor(dtaps[[v]])) {
            cur.p <- plot_balance_fac(dtaps, v, overall.inc = overall.inc);
        } else {
            cur.p <- plot_balance_cont(dtaps, v, nstrata = nstrata,
                                       overall.inc = overall.inc,
                                       facet.scales = facet.scales);
        }

        cur.p <- cur.p +
            labs(title = label.cov[v == v.cov]) +
            theme_bw() +
            theme(strip.background = element_blank(),
                  strip.placement  = "right",
                  strip.text       = element_blank(),
                  panel.grid       = element_blank(),
                  panel.border     = element_blank(),
                  panel.spacing    = unit(0, "lines"),
                  plot.title = element_text(hjust = 0.5),
                  legend.position  = "none",
                  plot.margin      = unit(c(1,0,1,-0.5), "lines"));

        rst[[v]] <- cur.p;
    }

    rst[[length(rst)]] <- rst[[length(rst)]] +
      theme(strip.text = element_text(size = 8),
            legend.position = "right")

    rst$nrow       <- 1;
    rst$rel_widths <- c(rep(1, length(v.cov)-1),
                        1 + legend.width * length(v.cov))
    do.call(plot_grid, rst)
}
