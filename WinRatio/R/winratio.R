winratio <- function(id = NULL, trt, active = NULL, outcomes, fu, data, keep.matrix = FALSE){
  
  old <- options("warning.length")
  on.exit(options(old))
  options(warning.length = 8170)
  
  ################################################################################;
  # CHECK TYPE OF ARGUMENTS ####
  ################################################################################;
  
  err <- NULL
  
  # Check if id is a string
  if (!is.vector(id) | !is.character(id) | length(id) != 1){
    err <- c(err, "The id argument must be a string indicating the patient ID variable.")
  }
  
  # Check if trt is a string
  if (!is.vector(trt) | !is.character(trt) | length(trt) != 1){
    err <- c(err, "The trt argument must be a string indicating the treatment variable.")
  }
  
  # Check if active is a vector of length 1 (if not NULL)
  if (!is.null(active)){
    if (!is.vector(active) | length(active) != 1){
      err <- c(err, "The active argument must be a numeric or string value used to define the active treatment group for the calculation of the win ratio.")
    }  
  }
  
  # Check if outcomes is a list of character which each element of the list outcomes must be of length 3
  if (!is.list(outcomes)){
    err <- c(err, "The outcomes argument must be a list used to define all outcomes in order of priority.")
  } else if (!is.character(unlist(outcomes)) | !(min(sapply(outcomes, length)) == 3 & max(sapply(outcomes, length)) == 3)){
    err <- c(err, "Each element of the outcomes argument must be a character vector of length 3 or a character list of length 3.")
  } else {
    types <- sapply(outcomes, "[[", 2)
    if (any(!types %in% c("s", "r", "c"))){
      err <- c(err, "For each element in the outcomes argument, the 2nd sub-element must be 's' for survival 'failure-time' events, 'r' for repeated survival 'failure-time' events or 'c' for continuous or ordinal 'non-failure time' events.")
    } else {
      if (any(types == "s")){
        if (sum(types == "s") == 1){
          if (!(length(sapply(outcomes[which(types=="s")], "[[", 1)) == 1 & length(sapply(outcomes[which(types=="s")], "[[", 3)) == 1)){
            err <- c(err, "Each element in the outcomes argument referring to a survival 'failure-time' event must be a character vector of length 3, where the 1st sub-element is a string indicating the name of the event variable and the 3rd sub-element is a string indicating the name of the time-to-event variable.")
          }
        } else {
          if (any(!sapply(sapply(outcomes[which(types=="s")], "[[", 1), length) == 1) | any(!sapply(sapply(outcomes[which(types=="s")], "[[", 3), length) == 1)){
            err <- c(err, "Each element in the outcomes argument referring to a survival 'failure-time' event must be a character vector of length 3, where the 1st sub-element is a string indicating the name of the event variable and the 3rd sub-element is a string indicating the name of the time-to-event variable.")
          }
        }
      }
      if (any(types == "r")){
        if (sum(types == "r") == 1){
          if ((length(sapply(outcomes[which(types=="r")], "[[", 1)) != length(sapply(outcomes[which(types=="r")], "[[", 3))) | 
              (!is.list(outcomes[which(types=="r")]))){
            err <- c(err, "Each element in the outcomes argument referring to a repeated survival 'failure-time' event must be a character list of length 3, where the 1st sub-element is a a character vector indicating the name of event variables and the 3rd sub-element is a character vector indicating the name of time-to-event variables corresponding to the event variables.")
          }
        } else {
          if (any(!sapply(sapply(outcomes[which(types=="r")], "[[", 1), length) == sapply(sapply(outcomes[which(types=="r")], "[[", 3), length)) |
              any(!(sapply(outcomes[which(types=="r")], is.list)))){
            err <- c(err, "Each element in the outcomes argument referring to a repeated survival 'failure-time' event must be a character list of length 3, where the 1st sub-element is a a character vector indicating the name of event variables and the 3rd sub-element is a character vector indicating the name of time-to-event variables corresponding to the event variables.")
          }
        }
      }
      if (any(types == "c")){
        if (sum(types == "c") == 1){
          if (!(length(sapply(outcomes[which(types=="c")], "[[", 1)) == 1)){
            err <- c(err, "Each element in the outcomes argument referring to a continuous or ordinal 'non-failure time' event must be a character vector of length 3, where the 1st sub-element is a string indicating the name of the continuous or ordinal variable and the 3rd sub-element is a string with two options '<' or '>' indicating the direction of deterioration or worsening.")
          } else {
            dir <- sapply(outcomes[which(types == "c")], "[[", 3)
            if (any(!dir %in% c(">","<"))){
              err <- c(err, "Each element in the outcomes argument referring to a continuous or ordinal 'non-failure time' event must be a character vector of length 3, where the 1st sub-element is a string indicating the name of the continuous or ordinal variable and the 3rd sub-element is a string with two options '<' or '>' indicating the direction of deterioration or worsening.")
            }
          }
        } else {
          if (any(!sapply(sapply(outcomes[which(types=="c")], "[[", 1), length) == 1)){
            err <- c(err, "Each element in the outcomes argument referring to a continuous or ordinal 'non-failure time' event must be a character vector of length 3, where the 1st sub-element is a string indicating the name of the continuous or ordinal variable and the 3rd sub-element is a string with two options '<' or '>' indicating the direction of deterioration or worsening.")
          } else {
            dir <- sapply(outcomes[which(types == "c")], "[[", 3)
            if (any(!unlist(dir) %in% c(">","<"))){
              err <- c(err, "Each element in the outcomes argument referring to a continuous or ordinal 'non-failure time' event must be a character vector of length 3, where the 1st sub-element is a string indicating the name of the continuous or ordinal variable and the 3rd sub-element is a string with two options '<' or '>' indicating the direction of deterioration or worsening.")
            }
          }
        }
      }
    }
  }
  
  # Check if fu is a character vector of length 1
  if (!is.vector(fu) | !is.character(fu) | length(fu) != 1){
    err <- c(err, "The fu argument must be a string indicating the name of the follow-up time variable.")
  }
  
  # Check if data is a data frame
  if (!is.data.frame(data)) {
    err <- c(err, "The data argument must be a data frame containing all the variables listed in the id, trt, outcomes and fu arguments.")
  }
  
  # Check if keep.matrix is a logical vector of length 1
  if (!is.vector(keep.matrix) | !is.logical(keep.matrix) | length(keep.matrix) != 1) {
    err <- c(err, "The keep.matrix argument must be a logical value indicating if the 'win-loss' matrix is kept.")
  }
  
  # Stop if at least one error
  if (!is.null(err)){
    stop(paste0("\n", paste0(err, collapse = "\n")), call. = F)
  }
  
  ################################################################################;
  # CHECKS IF VARIABLES EXIST IN THE DATAFRAME ####
  ################################################################################;
  
  vars <- unique(as.character(unlist(c(id, trt, 
                                       sapply(outcomes[which(types != "c")], "[[", 1),
                                       sapply(outcomes[which(types != "c")], "[[", 3),
                                       sapply(outcomes[which(types == "c")], "[[", 1), fu))))
  vars_miss <- setdiff(vars, names(data))
  if (length(vars_miss) > 0){
    stop(paste0("\nThe following variables are missing in the data frame: ", paste0("'", vars_miss, collapse = "', "), "'."), call. = F)
  }
  
  ################################################################################;
  # CHECKS IF VARIABLES ARE MISSING VALUES ####
  ################################################################################;
  
  # Check id variable (missing ID, duplicate ID)
  if (sum(is.na(data[[id]])) > 0){
    stop(paste0("\nMissing values for ID variable '", id,"'."), call. = F)
  } else if (nrow(data) > length(unique(data[[id]]))) {
    stop(paste0("\nDuplicate values for ID variable '", id,"'."), call. = F)
  }
  
  all_miss <- sapply(data[vars], FUN = function(x){all(is.na(x))})
  if (any(all_miss)){
    stop(paste0("\nThe following variables have only missing values: ", paste0("'", vars[all_miss == T], collapse = "', "), "'."), call. = F)
  }
  
  nb_miss <- sapply(data[vars], FUN = function(x){sum(is.na(x))})
  if (any(nb_miss > 0)){
    N <- nrow(data)
    data <- data %>% dplyr::select(all_of(vars)) %>% na.omit() %>% data.frame()
    warning("\nSome variables have missing values: ",
            paste0(paste0("'", vars[nb_miss > 0], "' (", nb_miss[nb_miss > 0], " NAs)"), collapse = ", "),".\n",
            paste0(N - nrow(data), " rows containing missing values were dropped"), call. = F, immediate. = T)
  }
  
  ################################################################################;
  # CHECKS TYPE OF EACH VARIABLE ####
  ################################################################################;
  
  # Check treatment variable + active argument
  if (length(setdiff(unique(data[[trt]]), c(NA, NaN))) != 2){
    stop(paste0("\nTreatment variable '", trt, "' must have 2 unique values/levels."), call. = F)
  }
  
  err <- NULL
  
  if (is.binary(data[[trt]])){
    if (is.null(active)){
      active <- 1
      warning(paste0("\nActive group = ", active), call. = F, immediate. = T)
    } else if (!active %in% unique(data[[trt]])){
      err <- c(err, "The active argument must be a binary value (0 or 1).")
    }
  } else if (is.numeric(data[[trt]])){
    if (is.null(active)){
      active <- max(data[[trt]])
      warning(paste0("\nActive group = ", active), call. = F, immediate. = T)
    } else if (!active %in% unique(data[[trt]])){
      err <- c(err, paste0("The active argument must be a numeric value (",
                           paste0(c(min(data[[trt]]), max(data[[trt]])), collapse=" or "), ")."))
    }
  } else if (is.character(data[[trt]])){
    if (is.null(active)){
      active <- sort(unique(data[[trt]]))[2]
      warning(paste0("\nActive group = '", active, "'"), call. = F, immediate. = T)
    } else if (!active %in% unique(data[[trt]])){
      err <- c(err, paste0("The active argument must be a string (", 
                           paste0(paste0("'", sort(unique(data[[trt]])), "'"), collapse=" or "), ")."))
    }
  } else if (is.factor(data[[trt]])){
    if (is.null(active)){
      tab <- table(data[[trt]])
      active <- names(tab)[tab>0][2]
      warning(paste0("\nActive group = '", active, "'"), call. = F, immediate. = T)
    } else if (!active %in% unique(data[[trt]])){
      err <- c(err, paste0("The active argument must be a string (", 
                           paste0(paste0("'", sort(unique(data[[trt]])), "'"), collapse=" or "), ")."))
    }
  } 
  
  # Check fu variable
  if (!is.numeric(data[[fu]])){
    err <- c(err, paste0("Follow-up time variable '", fu, "' must be a numeric variable with all values >= 0."))
  } else if (min(data[[fu]]) < 0){
    err <- c(err, paste0("Follow-up time variable '", fu, "' must be a numeric variable with all values >= 0."))
  }
  
  if (!is.null(err)){
    stop(paste0("\n",paste0(err, collapse = "\n")), call. = F)
  }
  
  # Check outcomes variables
  err <- vector("list", length(outcomes))
  for (k in seq_along(outcomes)){
    
    if (types[k] == "s"){
      if (!is.binary(data[[outcomes[[k]][1]]])){
        err[[k]] <- c(err[[k]], paste0("'", outcomes[[k]][1], "' must be a binary variable."))
      }
      if (!is.numeric(data[[outcomes[[k]][3]]])){
        err[[k]] <- c(err[[k]], paste0("'", outcomes[[k]][3], "' must be a numeric variable with all values >= 0."))
      } else if (min(data[[outcomes[[k]][3]]]) < 0){
        err[[k]] <- c(err[[k]], paste0("'", outcomes[[k]][3], "' must be a numeric variable with all values >= 0."))
      } else if (any(data[[outcomes[[k]][3]]] > data[[fu]])){
        err[[k]] <- c(err[[k]], paste0("'", outcomes[[k]][3], "' must be <= '", fu, "'"))
      }
    } else if (types[k] == "r"){
      events <- unlist(outcomes[[k]][1])
      test_bin <- sapply(data[, events], function(x){is.binary(x)})
      times <- unlist(outcomes[[k]][3])
      test_num <- sapply(data[, times], function(x){is.numeric(x)})
      if (any(!test_bin)){
        err[[k]] <- c(err[[k]], paste0("The following variables must be binary: ", 
                                       paste0(paste0("'",events[which(!test_bin)], "'"), collapse = ", "), "."))
      }
      if (any(!test_num)){
        err[[k]] <- c(err[[k]], paste0("The following variables must be numeric with all values >= 0: ", 
                                       paste0(paste0("'",times[which(!test_num)], "'"), collapse = ", "), "."))
      } else {
        test_pos <- sapply(data[, times], function(x){min(x) >=0 })
        if (any(!test_pos)){
          err[[k]] <- c(err[[k]], paste0("The following variables must be numeric with all values >= 0: ", 
                                         paste0(paste0("'",times[which(!test_pos)], "'"), collapse = ", "), "."))
        } 
      }
      if (all(test_num) & all(test_bin)){
        if (all(test_pos)){
          if (length(events) == 1){
            if (any(data[[times]] > data[[fu]])){
              err[[k]] <- c(err[[k]], paste0("'", times, "' must be <= '", fu, "'."))
            }
          } else {
            for (j in 1:(length(events))){
              if (j < length(events)){
                if (any(!data[events[j]] >= data[events[j+1]])){
                  err[[k]] <- c(err[[k]], paste0("'", events[j], " must be >= '", events[j+1], "'."))
                }
                if (any(!data[times[j]] <= data[times[j+1]])){
                  err[[k]] <- c(err[[k]], paste0("'", times[j], " must be <= '", times[j+1], "'."))
                }
              }
              if (any(!data[[times[j]]] <= data[[fu]])){
                err[[k]] <- c(err[[k]], paste0("'", times[j], " must be <= '", fu, "'."))
              }
            }
          }
        }
      }
      
      if (is.null(err[[k]])){
        for (j in 1:(length(events))){
          index <- which(data[[events[j]]] == 0 & data[[times[j]]] != data[[fu]])
          if (length(index) > 0){
            err[[k]] <- c(err[[k]], paste0("'", times[j], "' must be == '", fu, "' for rows with '", events[j], "' == 0."))
          }
          if (j<length(events)){
            index11 <- which(data[[events[j]]] == 1 & data[[events[j+1]]] == 1 & data[[times[j]]] >= data[[times[j+1]]])
            if (length(index11) > 0){
              err[[k]] <- c(err[[k]], paste0("'", times[j], "' must be < '", times[j+1], "' for rows with '", events[j], "' == 1 and '" , events[j+1], "' == 1."))
            }
            index10 <- which(data[[events[j]]] == 1 & data[[events[j+1]]] == 0 & data[[times[j]]] > data[[times[j+1]]])
            if (length(index10) > 0){
              err[[k]] <- c(err[[k]], paste0("'", times[j], "' must be <= '", times[j+1], "' for rows with '", events[j], "' == 1 and '" , events[j+1], "' == 0."))
            }
            index00 <- which(data[[events[j]]] == 0 & data[[events[j+1]]] == 0 & data[[times[j]]] != data[[times[j+1]]])
            if (length(index00) > 0){
              err[[k]] <- c(err[[k]], paste0("'", times[j], "' must be == '", times[j+1], "' for rows with '", events[j], "' == 0 and '" , events[j+1], "' == 0."))
            }
          }
        }
      }
    } else if (types[k] == "c"){
      event <- outcomes[[k]][1] 
      if (!is.numeric(data[[event]])){
        err[[k]] <- c(err[[k]], paste0("'", event, "' must be a numeric variable."))
      }
    }
  }
  
  
  index <- which(!sapply(err, function(x){is.null(x)})) 
  if (length(index) > 0){
    stop(paste0("\n", paste0(paste0("Outcome ", index, ":\n", sapply(err, function(x){paste0(x, collapse = "\n")})[index]), collapse = "\n")), call. = F)
  }
  
  ################################################################################;
  # Generate matrix for calculation of win ratio
  ################################################################################;
  index1 <- which(data[[trt]]== active)
  n1 <- length(index1)
  fu1 <- data[index1, fu]
  index0 <- which(data[[trt]] != active)
  n0 <- length(index0)
  fu0 <- data[index0, fu]
  mat <- matrix(data = 0,  nrow = n1, ncol = n0)
  mat_list <- NULL
  for (k in 1:length(outcomes)){
    if (outcomes[[k]][2] == "s"){
      mat_list[[k]] <- .Call(`_WinRatio_mat_comp_surv_cpp`, data[[outcomes[[k]][3]]][index1],
                             data[[outcomes[[k]][1]]][index1],
                             data[[outcomes[[k]][3]]][index0],
                             data[[outcomes[[k]][1]]][index0]) * k
    } else if (outcomes[[k]][2] == "r"){
      mat_list[[k]] <- .Call(`_WinRatio_mat_comp_repeated_cpp`, as.matrix(data[unlist(outcomes[[k]][3])][index1,]), 
                                             as.matrix(data[unlist(outcomes[[k]][1])][index1,]),
                                             data[[fu]][index1],
                                             as.matrix(data[unlist(outcomes[[k]][3])][index0,]), 
                                             as.matrix(data[unlist(outcomes[[k]][1])][index0,]),
                                             data[[fu]][index0]) * k
    } else if (outcomes[[k]][2] == "c"){
      mat_list[[k]] <- .Call(`_WinRatio_mat_comp_cont_cpp`, data[[outcomes[[k]][1]]][index1],
                                         data[[outcomes[[k]][1]]][index0],
                                         direction = outcomes[[k]][3]) * k
    } 
    mat[which(mat == 0)] <- mat_list[[k]][which(mat == 0)]
  }
  
  if (keep.matrix){
    wr.matrix <- mat
  } else {
    wr.matrix <- NULL
  }
  
  # Calculate win ratio statistic and asymptotic confidence interval (p-value)
  n <- n1+n0
  
  # Number of wins and losses for active group
  wins <- loss <- rep(0, length(outcomes))
  for (k in 1:max(abs(mat))){
    wins[k] <- sum(mat == k)
    loss[k] <- sum(mat == -k)
  }
  ttw <- sum(wins)/(n^2)
  ttl <- sum(loss)/(n^2)
  
  # Win ratio statistics
  if (ttw>0 & ttl>0){
    wr <- ttw/ttl
    v <- (sum(((apply(mat>0, 1, sum) - apply(mat<0, 1, sum) * wr)/n0)^2)*(n0/n)^2/n +
      sum(((apply(mat>0, 2, sum) - apply(mat<0, 2, sum) * wr)/n1)^2)*(n1/n)^2/n)/(ttl^2)
    z <- sqrt(n)*log(wr)*wr/sqrt(v)
    p.value <- pnorm(abs(z), lower.tail = F) * 2
    wr.lower <- wr * exp(- qnorm(0.975) * sqrt(v/wr^2/n))
    wr.upper <- wr * exp(qnorm(0.975) * sqrt(v/wr^2/n))
  } else {
    wr <- v <- z <- p.value <- wr.lower <- wr.upper <- NA
    warning("The win ratio is not defined.", call. = F, immediate. = T)
  }
  
  res <- list(call = list(id = id, trt = trt, active = active, outcomes = outcomes, fu = fu, data = data, keep.matrix = keep.matrix), 
              group1 = active, group0 = setdiff(unique(data[[trt]]), active), n1 = n1, n0 = n0, n = n, 
              wins = wins, loss = loss, total.wins = sum(wins), total.loss = sum(loss), total.ties = n1 * n0 - (sum(wins) + sum(loss)), 
              wr = wr, v = v, z = z, p.value = p.value, wr.lower = wr.lower, wr.upper = wr.upper, wr.matrix = wr.matrix)
  class(res) <- "WinRatio"
  return(res)
}
