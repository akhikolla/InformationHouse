assign_parm = function(est_nugget, est_par3, est_angle,
                       est_ratio, mod, nugget) {
  # extract relevant parameters
  r = mod$r
  psill = mod$psill
  lambda = nugget/psill
  angle = mod$angle
  ratio = mod$ratio
  par3 = mod$par3
  
  if (!est_nugget & !est_par3 & !est_angle & !est_ratio) {
    parm = c(r, psill)
    names(parm) = c("r", "psill")
  } else if (!est_nugget & est_par3 & !est_angle & !est_ratio) {
    parm = c(r, psill, par3)
    names(parm) = c("r", "psill", "par3")
  } else if (!est_nugget & !est_par3 & est_angle & !est_ratio) {
    parm = c(r, psill, angle)
    names(parm) = c("r", "psill", "angle")
  } else if (!est_nugget & !est_par3 & !est_angle & est_ratio) {
    parm = c(r, psill, ratio)
    names(parm) = c("r", "psill", "ratio")
  } else if (!est_nugget & est_par3 & est_angle & !est_ratio) {
    parm = c(r, psill, angle, par3)
    names(parm) = c("r", "psill", "angle", "par3")
  } else if (!est_nugget & est_par3 & !est_angle & est_ratio) {
    parm = c(r, psill, ratio, par3)
    names(parm) = c("r", "psill", "ratio", "par3")
  } else if (!est_nugget & !est_par3 & est_angle & est_ratio) {
    parm = c(r, psill, angle, ratio)
    names(parm) = c("r", "psill", "angle", "ratio")
  } else if (!est_nugget & est_par3 & est_angle & est_ratio) {
    parm = c(r, psill, angle, ratio, par3)
    names(parm) = c("r", "psill", "angle", "ratio", "par3")
  } else if (est_nugget & !est_par3 & !est_angle & !est_ratio) {
    parm = c(r, lambda)
    names(parm) = c("r", "lambda")
  } else if (est_nugget & est_par3 & !est_angle & !est_ratio) {
    parm = c(r, lambda, par3)
    names(parm) = c("r", "lambda", "par3")
  } else if (est_nugget & !est_par3 & est_angle & !est_ratio) {
    parm = c(r, lambda, angle)
    names(parm) = c("r", "lambda", "angle")
  } else if (est_nugget & !est_par3 & !est_angle & est_ratio) {
    parm = c(r, lambda, ratio)
    names(parm) = c("r", "lambda", "ratio")
  } else if (est_nugget & est_par3 & est_angle & !est_ratio) {
    parm = c(r, lambda, angle, par3)
    names(parm) = c("r", "lambda", "angle", "par3")
  } else if (est_nugget & est_par3 & !est_angle & est_ratio) {
    parm = c(r, lambda, ratio, par3)
    names(parm) = c("r", "lambda", "ratio", "par3")
  } else if (est_nugget & !est_par3 & est_angle & est_ratio) {
    parm = c(r, lambda, angle, ratio)
    names(parm) = c("r", "lambda", "angle", "ratio")
  } else if (est_nugget & est_par3 & est_angle & est_ratio) {
    parm = c(r, lambda, angle, ratio, par3)
    names(parm) = c("r", "lambda", "angle", "ratio", "par3")
  } 
  return(parm)
}
