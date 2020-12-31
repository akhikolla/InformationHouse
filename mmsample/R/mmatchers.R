mmatcher <- function(ds, group_var, x_vars = "_all_", id_var = NA,
                     distance = "mahal", caliper = 0.10, seed = 12345,
                     max_candidates = 1000, n_per_match = 1, loud = TRUE) {
  # check ds is a dataframe
  if (!is.data.frame(ds)) { stop("ds argument must be a data.frame.") }

  # check distance name is valid, default (mahalanobia) otherwise
  if(!(distance %in% c("mahal", "euclid", "norm_euclid", "sad"))) {
    distance <- "mahal"
    warning("Invalid distance, defaulting to 'mahal' (Mahalanobis).")
  }

  # check id and group variables appear in ds
  if ((!is.na(id_var) & !(id_var %in% names(ds))) | !(group_var %in% names(ds))) {
    stop("id_var and group_var must appear in ds.")
  }

  # get x_var names if "_all_", otherwise stop if any not present
  if (length(x_vars) == 1 & x_vars[1] == "_all_") {
    if (!is.na(id_var)) {
      x_vars <- names(ds)[!(names(ds) %in% c(id_var, group_var))]
    } else {
      x_vars <- names(ds)[names(ds) != group_var]
    }
  } else if (!(all(x_vars %in% names(ds)))) {
    stop("All variables in x_vars must appear in ds.")
  }

  # remove any rows containing NAs on required variables
  req_vars <- c(group_var, x_vars)
  if (!is.na(id_var)) req_vars <- c(id_var, req_vars)
  for (v in req_vars) { ds <- ds[!is.na(ds[[v]]), ] }

  # check sufficient data remains
  if (min(table(ds[[group_var]])) < 2) {
    stop("Insufficient data. Try specifying fewer variables using x_vars argument.")
  }

  # propensity scores and caliper width
  prp_var <- "propensity_score"
  if (prp_var %in% names(ds)) {
    i = 1 ; while (prp_var %in% names(ds)) { prp_var <- paste0("propensity_score", i) ; i = i+1 }
  }
  ds[[prp_var]] <- round(stats::predict(stats::glm(stats::as.formula(paste(group_var, "~", paste(x_vars, collapse = "+"))),
                                     ds, family = "binomial"), type = "response", newdata = ds), 6)
  ds <- ds[order(ds[[prp_var]]), ]
  x_vars <- c(x_vars, prp_var)
  d_cal <- diff(range(ds[[prp_var]])) * caliper

  # inverted covariance matrix
  if(distance == "mahal") {
    cov_inv <- MASS::ginv(stats::cov(ds[, x_vars]))
  } else if (distance == "euclid") {
    cov_inv <- diag(length(x_vars))
  } else if (distance == "norm_euclid") {
    cov_inv <- diag(MASS::ginv(stats::cov(ds[, x_vars]))) * diag(length(x_vars))
  }

  # initialise
  set.seed(seed)
  treatment_n <- sum(ds[[group_var]])
  five_pc <- floor(0.05 * treatment_n)
  if (five_pc == 0) five_pc <- 1
  control_flag <- ds[[group_var]] == 0

  dist_var <- paste0(distance, "_dist")
  if (dist_var %in% names(ds)) {
    i = 1 ; while (dist_var %in% names(ds)) { dist_var <- paste0(distance, "_dist", i) ; i = i+1 }
  }

  pair_var <- "mmatch_pair"
  if (pair_var %in% names(ds)) {
    i = 1 ; while (pair_var %in% names(ds)) { pair_var <- paste0("mmatch_pair", i) ; i = i+1 }
  }

  ds[[pair_var]] <- ds[[dist_var]] <- 0
  ds[[pair_var]][ds[[group_var]] == 1] <- sample(treatment_n)

  # start outer assignment loop (once per required number of matches)
  for (k in 1:n_per_match) {
    # start greedy matching loop (once per target)
    for (target in 1:treatment_n) {
      # apply calipers
      tgt_prp = ds[[prp_var]][ds[[pair_var]] == target & ds[[group_var]] == 1]
      candidates <- ds[(control_flag &
                          ds[[pair_var]] == 0 &
                          ds[[prp_var]] <= (tgt_prp + d_cal) &
                          ds[[prp_var]] >= (tgt_prp - d_cal)), ]

      # proceed if any candidates
      if (nrow(candidates) >= 1) {
        u <- as.numeric(ds[ds[[pair_var]] == target & ds[[group_var]] == 1, x_vars])

        # cut down list to closest on propensity if too many
        if (nrow(candidates) > max_candidates) {
          candidates[["abs_gap"]] <- abs(tgt_prp - candidates[[prp_var]])
          candidates <- candidates[order(candidates[["abs_gap"]]), ]
          candidates <- candidates[1:max_candidates, ]
        }

        # calculate distances
        if (distance == "sad") {
          candidates[[dist_var]] <- round(colSums(abs(t(candidates[, x_vars]) - u)), 8)
        } else {
          candidates[[dist_var]] <- round(sqrt(as.vector(ruler(as.matrix(candidates[, x_vars]), u, cov_inv))), 8)
        }

        # get selected match and update pair number
        selected <- candidates[candidates[[dist_var]] == min(candidates[[dist_var]]), ]
        selected <- selected[1,]
        ds[row.names(ds) == row.names(selected), c(pair_var, dist_var)] <- c(target, selected[[dist_var]])
      } else {
        ds[[pair_var]][ds[[pair_var]] == target] <- 0
      }

      if (loud) {
        if (target %% five_pc == 0) {
          pc_done <- 5 * target / five_pc
          if (pc_done <= 100) cat("|", pc_done, sep = "")
        }
      }

    }
  }

  paired_ds <- ds[order(ds[[pair_var]]), ]
  paired_ds <- paired_ds[paired_ds[[pair_var]] != 0, ]

  if (loud) {
    fmt_rnd <- function(x, d = 3) format(round(x, 3), nsmall = 3)
    distance_values <- paired_ds[[dist_var]][paired_ds[[group_var]] == 0]
    distance_stats_string <- paste0(fmt_rnd(mean(distance_values), 3),
           " (", paste(fmt_rnd(as.vector(stats::quantile(distance_values, c(0.025, 0.975))), 3), collapse = ", "), ")")
    cat("\nMatched: ", sum(paired_ds[[group_var]]), "/", treatment_n ,
        "\tMean (95% CI) distance: ", distance_stats_string, "\n", sep = "")
  }

  return(paired_ds)
}


