partial.auc.index <- function(partial.auc, area, give.warning, tpr.area) {
  # case area = FPR, assume TPR = FPR
  auc.min <- 0.5 * (area[2]^2 - area[1]^2)
  auc.max <- area[2] - area[1]
  
  if (tpr.area) auc.min <- auc.max - auc.min
  
  # care area = TPR, assume FPR = TPR
  
  
  
  partial.auc <- 0.5*(1 + ((partial.auc - auc.min) / (auc.max - auc.min))) #McClish
  if (any(partial.auc < 0.5) & give.warning) warning("Partial AUCs under 0.5 generated.")
  
  return(partial.auc)
}



# transforms list to roc data.frame
c.list.to.roc <- function(input) {
  roc.frame <- as.data.frame(input[1:3])
  names(roc.frame) <- c("TPR", "FPR", "threshold")
  return(roc.frame)
}

# Check if number is a single numeric between 0 and 1
validate.single.numeric <- function(number, var.name) {
  if (is.null(number)) stop(paste("Please pass ", var.name, " to perf.roc!", sep = ""))
  if (class(number) != "numeric") stop(paste(var.name, " must be numeric!", sep = ""))
  if (length(number) != 1) stop(paste(var.name, " must have length 1!", sep = ""))
  if (is.na(number)) stop(paste(var.name, " is NA!", sep = ""))
  if ((number < 0) | (number > 1)) stop(paste(var.name, " must be in [0, 1]!", sep = ""))
  return(number)
}


fbroc.plot.base <- function(plot.frame) {
  roc.plot <- ggplot(data = plot.frame, aes(x = FPR, y = TPR)) +               
              ggtitle("ROC Curve") + xlab("False Positive Rate") +
              ylab("True Positive Rate") + theme_bw() +
              theme(title = element_text(size = 22),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_text(size = 16))
  return(roc.plot)
}

# adds confidence region to roc curve plot
fbroc.plot.add.conf <- function(roc1, conf.level = conf.level, steps = steps, fill = fill) {
  conf.frame <- conf(roc1, conf.level = conf.level, steps = steps)
  conf.frame$Segment <- 1
 
  geom_ribbon(data = conf.frame, fill = fill, alpha = 0.5,
              aes(y = NULL, ymin = Lower.TPR, ymax = Upper.TPR))
}




# add performance metric visualization to roc plot (paired roc curve)
fbroc.plot.add.metric.paired <- function(x, roc.plot,show.metric, show.area, perf, 
                                         col1, fill1, col2, fill2) {
  if (show.metric == "tpr") {
    extra.frame <- data.frame(FPR = perf$params, 
                              TPR = c(perf$Observed.Performance.Predictor1, 
                                      perf$Observed.Performance.Predictor2),
                              Segment = 1,
                              lower = c(perf$CI.Performance.Predictor1[1], 
                                        perf$CI.Performance.Predictor2[1]),
                              upper = c(perf$CI.Performance.Predictor1[2],
                                        perf$CI.Performance.Predictor2[2]))
  
    roc.plot <- roc.plot + geom_errorbar(data = extra.frame, width = 0.02, size = 1.25,
                              aes(ymin = lower, ymax = upper), col = c(col1, col2), alpha = 0.7) + 
      geom_point(data = extra.frame, size = 4, col = c(col1, col2))
  }
  if (show.metric == "fpr") {
    extra.frame <- data.frame(TPR = perf$params, 
                              FPR = c(perf$Observed.Performance.Predictor1, 
                                      perf$Observed.Performance.Predictor2),
                              Segment = 1,
                              lower = c(perf$CI.Performance.Predictor1[1], 
                                        perf$CI.Performance.Predictor2[1]),
                              upper = c(perf$CI.Performance.Predictor1[2],
                                        perf$CI.Performance.Predictor2[2]))
    roc.plot <- roc.plot + geom_errorbarh(data = extra.frame, height = 0.02, size = 1.25,
                              aes(xmin = lower, xmax = upper), col = c(col1, col2), alpha = 0.7) +
      geom_point(data = extra.frame, size = 4)
  }
  if (show.metric == "auc" & show.area) {
    roc.plot <- fbroc.plot.add.auc(extract.roc(x, 1), roc.plot, fill1)
    roc.plot <- fbroc.plot.add.auc(extract.roc(x, 2), roc.plot, fill2)
  }
  if (show.metric == "partial.auc") { # This function does the main work
    roc.plot <- fbroc.plot.add.partial.auc(extract.roc(x, 1), roc.plot, show.metric, show.area, 
                                           perf, fill1)
    roc.plot <- fbroc.plot.add.partial.auc(extract.roc(x, 2), roc.plot, show.metric, show.area, 
                                           perf, fill2)
  }
  
  return(roc.plot)
}

remove.tpr.shift <- function(roc.data) {
  shift <- c(1, diff(roc.data$TPR))
  roc.data <- subset(roc.data, shift != 0)
  return(roc.data)
}


# code to add partial AUC area to ROC plot
fbroc.plot.add.partial.auc <- function(x, roc.plot, show.metric, show.area, perf, fill.col) {
  
  tpr.area <- grepl("TPR", perf$metric)
  tpr.m <- matrix(x$roc$TPR, nrow = 1)
  fpr.m <- matrix(x$roc$FPR, nrow = 1)
  # differentiate between TPR area and FPR area
  if (tpr.area) {
    tpr <- perf$params
    fpr <- c(get_cached_perf(tpr.m, fpr.m, tpr[1], 2),
             get_cached_perf(tpr.m, fpr.m, tpr[2], 2))
    frame.line1 <- data.frame(TPR = tpr[1], FPR = c(1, fpr[1]))
    frame.line2 <- data.frame(TPR = tpr[2], FPR = c(1, fpr[2]))
    rel.roc <- subset(x$roc, (TPR >= tpr[1]) & (TPR <= tpr[2]))
  } else {
    fpr <- perf$params
    tpr <- c(get_cached_perf(tpr.m, fpr.m, fpr[1], 1),
             get_cached_perf(tpr.m, fpr.m, fpr[2], 1))
    frame.line1 <- data.frame(FPR = fpr[1], TPR = c(0, tpr[1]))
    frame.line2 <- data.frame(FPR = fpr[2], TPR = c(0, tpr[2]))
    rel.roc <- subset(x$roc, (FPR >= fpr[1]) & (FPR <= fpr[2]))
  }
  if (show.area) { # Show partial AUC
    first.row <- data.frame(TPR = tpr[2], FPR = fpr[2], threshold = as.numeric(NA))
    last.row <- data.frame(TPR = tpr[1], FPR = fpr[1], threshold = as.numeric(NA))
    rel.roc <- rbind(first.row, rel.roc, last.row)
    if (x$tie.strategy == 2) {
      
      expand.roc <- add_roc_points(rel.roc$TPR, rel.roc$FPR)
      rel.roc <- data.frame(TPR = expand.roc[[1]],
                            FPR = expand.roc[[2]],
                            Segment = expand.roc[[3]])
    } 
    rel.roc <- unique(rel.roc[order(rel.roc$FPR, rel.roc$TPR), ])
    if (tpr.area) {
      rel.roc$TPR.MIN <- tpr[1]
      rel.roc <- rbind(rel.roc, rel.roc[nrow(rel.roc), ])
      rel.roc[nrow(rel.roc), "FPR"] <- 1
      roc.plot <- roc.plot + geom_ribbon(data = rel.roc, fill = fill.col, alpha = 0.5, 
                                         aes(ymin = TPR.MIN, ymax = TPR, y = 0))
    }
    else {
      roc.plot <- 
        roc.plot + geom_ribbon(data = rel.roc, fill = fill.col, alpha = 0.5, 
                               aes(ymin = 0, ymax = TPR, y = 0))
    }
  }
  
  roc.plot <- roc.plot + geom_line(data = frame.line1, linetype = 3, size = 1.1) + 
    geom_line(data = frame.line2, linetype = 3, size = 1.1)
  return(roc.plot)
}

fbroc.plot.add.auc <- function(x, roc.plot, fill.col) {
  rel.roc <- x$roc[nrow(x$roc):1, ]
  roc.plot <- roc.plot + geom_ribbon(data = rel.roc, fill = fill.col, alpha = 0.5, 
                                     aes(ymin = 0, ymax = TPR, y = 0))
  return(roc.plot)
}

# add performance metric visualization to roc plot (roc curve)
fbroc.plot.add.metric <- function(x, roc.plot, show.metric, show.area, perf, fill.col) {
  
  if (show.metric == "tpr") {
    extra.frame <- data.frame(FPR = perf$params, TPR = perf$Observed.Performance, Segment = 1,
                              lower = perf$CI.Performance[1], upper = perf$CI.Performance[2])
    roc.plot <- roc.plot + geom_errorbar(data = extra.frame, width = 0.02, size = 1.25,
                                         aes(ymin = lower, ymax = upper)) + 
      geom_point(data = extra.frame, size = 4)
  }
  
  if (show.metric == "fpr") {
    extra.frame <- data.frame(TPR = perf$params, FPR = perf$Observed.Performance, Segment = 1,
                              lower = perf$CI.Performance[1], upper = perf$CI.Performance[2])
    roc.plot <- roc.plot + geom_errorbarh(data = extra.frame, height = 0.02, size = 1.25,
                                          aes(xmin = lower, xmax = upper)) +
      geom_point(data = extra.frame, size = 4)
  }
  
  if (show.metric == "auc" & show.area) {
    roc.plot <- fbroc.plot.add.auc(x, roc.plot, fill.col)
  }
  
  if (show.metric == "partial.auc") { # This function does the main work
    roc.plot <- fbroc.plot.add.partial.auc(x, roc.plot, show.metric, show.area, perf, fill.col)
  }
  return(roc.plot)
}

