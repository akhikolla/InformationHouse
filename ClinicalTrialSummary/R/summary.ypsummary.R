summary.ypsummary <- function(object, ...) {
  x <- object
  if (class(x)!="ypsummary") stop("Please use the object from the 'ypsummary' function")

  digit <- max(3, getOption("digits") - 3)
  v_tau <- round(x$tau, digits = digit)
  v_repnum <- x$repnum

  #--- The average hazard ratio ---#
  x_ahrypt <- lapply(x$values$ahrypt, round, digits = digit)
  est_ahrypt <- round(x_ahrypt$estimate, digits = digit)
  lower_ahrypt <- round(x_ahrypt$lower, digits = digit)
  upper_ahrypt <- round(x_ahrypt$upper, digits = digit)
  z_ahrypt <- round(x_ahrypt$z, digits = digit)
  pv_ahrypt <- round(2 * (1 - pnorm(abs(x_ahrypt$z))), digits = digit)
  pv_ahrypt <- fun_less(pv_ahrypt)

  #--- The weighted average hazard ratio ---#
  x_ahrf <- lapply(x$values$ahrf, round, digits = digit)
  est_ahrf <- round(x_ahrf$estimate, digits = digit)
  lower_ahrf <- round(x_ahrf$lower, digits = digit)
  upper_ahrf <- round(x_ahrf$upper, digits = digit)
  z_ahrf <- round(x_ahrf$z, digits = digit)
  pv_ahrf <- round(2 * (1 - pnorm(abs(x_ahrf$z))), digits = digit)
  pv_ahrf <- fun_less(pv_ahrf)

  #--- The restricted superiority probability ratio ---#
  x_rspypt <- lapply(x$values$rspypt, round, digits = digit)
  est_rspypt <- round(x_rspypt$estimate, digits = digit)
  lower_rspypt <- round(x_rspypt$lower, digits = digit)
  upper_rspypt <- round(x_rspypt$upper, digits = digit)
  z_rspypt <- round(x_rspypt$z, digits = digit)
  pv_rspypt <- round(2 * (1 - pnorm(abs(x_rspypt$z))), digits = digit)
  pv_rspypt <- fun_less(pv_rspypt)

  #--- The restricted mean survival difference ---#
  x_mdypt <- lapply(x$values$mdypt, round, digits = digit)
  est_mdypt <- round(x_mdypt$estimate, digits = digit)
  lower_mdypt <- round(x_mdypt$lower, digits = digit)
  upper_mdypt <- round(x_mdypt$upper, digits = digit)
  z_mdypt <- round(x_mdypt$z, digits = digit)
  pv_mdypt <- round(2 * (1 - pnorm(abs(x_mdypt$z))), digits = digit)
  pv_mdypt <- fun_less(pv_mdypt)

  #--- The ratio of restricted mean times lost ---#
  x_mrypt <- lapply(x$values$mrypt, round, digits = digit)
  est_mrypt <- round(x_mrypt$estimate, digits = digit)
  lower_mrypt <- round(x_mrypt$lower, digits = digit)
  upper_mrypt <- round(x_mrypt$upper, digits = digit)
  z_mrypt <- round(x_mrypt$z, digits = digit)
  pv_mrypt <- round(2 * (1 - pnorm(abs(x_mrypt$z))), digits = digit)
  pv_mrypt <- fun_less(pv_mrypt)

  row_ahrypt <- c(v_repnum, v_tau, est_ahrypt, lower_ahrypt, upper_ahrypt,
                  z_ahrypt, pv_ahrypt)
  row_ahrf <- c(v_repnum, v_tau, est_ahrf, lower_ahrf, upper_ahrf, z_ahrf,
                pv_ahrf)
  row_rspypt <- c(v_repnum, v_tau, est_rspypt, lower_rspypt, upper_rspypt,
                  z_rspypt, pv_rspypt)
  row_mdypt <- c(v_repnum, v_tau, est_mdypt, lower_mdypt, upper_mdypt, z_mdypt,
                 pv_mdypt)
  row_mrypt <- c(v_repnum, v_tau, est_mrypt, lower_mrypt, upper_mrypt, z_mrypt,
                 pv_mrypt)

  result <- rbind(row_ahrypt, row_ahrf, row_rspypt, row_mdypt, row_mrypt)
  rownames(result) <- c("AHR", "WAHR", "RSPR", "RMSD", "RRMTL")
  colnames(result) <- c("repnum", "tau", "estimate", "lower.CI", "upper.CI", "z-value", "p-value")

  as.data.frame(result)
}
