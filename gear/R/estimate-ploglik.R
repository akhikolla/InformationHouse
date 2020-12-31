# all of these functions compute the log-likelihood
# based on different parameter values
# used as objective functions for the optimx::optimx
# function within the estimate.geolm_cmodStd function.

ploglik_cmodStd_r_psill = function(parm, x, y, d, weights,
                                   scmod, nugget, mu = NULL,
                                   reml = FALSE,
                                   return_ll = TRUE) {
  scmod$r = parm[1]
  scmod$psill = parm[2]
  v = evaluate(scmod, d) + nugget * diag(1/weights)
  cholv = chol(v)
  ll_xycholv(x, y, cholv, mu = mu, reml = reml,
             return_ll = return_ll)
}

ploglik_cmodStd_r_psill_par3 = function(parm, x, y, d, weights,
                                        scmod, nugget,
                                        mu = NULL, 
                                        reml = FALSE,
                                        return_ll = TRUE) {
  scmod$r = parm[1]
  scmod$psill = parm[2]
  scmod$par3 = parm[3]
  v = evaluate(scmod, d) + nugget * diag(1/weights)
  cholv = chol(v)
  ll_xycholv(x, y, cholv, mu = mu, reml = reml,
             return_ll = return_ll)
}

ploglik_cmodStd_r_psill_angle = function(parm, x, y, d, weights,
                                         scmod, nugget,
                                         mu = NULL, 
                                         reml = FALSE,
                                         return_ll = TRUE) {
  scmod$r = parm[1]
  scmod$psill = parm[2]
  scmod$angle = parm[3]
  v = evaluate(scmod, d) + nugget * diag(1/weights)
  cholv = chol(v)
  ll_xycholv(x, y, cholv, mu = mu, reml = reml, 
             return_ll = return_ll)
}

ploglik_cmodStd_r_psill_ratio = function(parm, x, y, d, weights,
                                         scmod, nugget,
                                         mu = NULL, 
                                         reml = FALSE,
                                         return_ll = TRUE) {
  scmod$r = parm[1]
  scmod$psill = parm[2]
  scmod$ratio = parm[3]
  v = evaluate(scmod, d) + nugget * diag(1/weights)
  cholv = chol(v)
  ll_xycholv(x, y, cholv, mu = mu, reml = reml, 
             return_ll = return_ll)
}

ploglik_cmodStd_r_psill_angle_ratio = function(parm, x, y, d, weights,
                                               scmod, nugget,
                                               mu = NULL, 
                                               reml = FALSE,
                                               return_ll = TRUE) {
  scmod$r = parm[1]
  scmod$psill = parm[2]
  scmod$angle = parm[3]
  scmod$ratio = parm[4]
  v = evaluate(scmod, d) + nugget * diag(1/weights)
  cholv = chol(v)
  ll_xycholv(x, y, cholv, mu = mu, reml = reml, 
             return_ll = return_ll)
}

ploglik_cmodStd_r_psill_angle_par3 = function(parm, x, y, d, weights,
                                         scmod, nugget,
                                         mu = NULL, 
                                         reml = FALSE,
                                         return_ll = TRUE) {
  scmod$r = parm[1]
  scmod$psill = parm[2]
  scmod$angle = parm[3]
  scmod$par3 = parm[4]
  v = evaluate(scmod, d) + nugget * diag(1/weights)
  cholv = chol(v)
  ll_xycholv(x, y, cholv, mu = mu, reml = reml, 
             return_ll = return_ll)
}

ploglik_cmodStd_r_psill_ratio_par3 = function(parm, x, y, d, weights,
                                         scmod, nugget,
                                         mu = NULL, 
                                         reml = FALSE,
                                         return_ll = TRUE) {
  scmod$r = parm[1]
  scmod$psill = parm[2]
  scmod$ratio = parm[3]
  scmod$par3 = parm[4]
  v = evaluate(scmod, d) + nugget * diag(1/weights)
  cholv = chol(v)
  ll_xycholv(x, y, cholv, mu = mu, reml = reml, 
             return_ll = return_ll)
}

ploglik_cmodStd_r_psill_angle_ratio_par3 = function(parm, x, y, d, weights,
                                               scmod, nugget,
                                               mu = NULL, 
                                               reml = FALSE,
                                               return_ll = TRUE) {
  scmod$r = parm[1]
  scmod$psill = parm[2]
  scmod$angle = parm[3]
  scmod$ratio = parm[4]
  scmod$par3 = parm[5]
  v = evaluate(scmod, d) + nugget * diag(1/weights)
  cholv = chol(v)
  ll_xycholv(x, y, cholv, mu = mu, reml = reml, 
             return_ll = return_ll)
}

ploglik_cmodStd_r_lambda = function(parm, x, y, d, weights,
                                    scmod, nugget, mu = NULL,
                                    reml = FALSE,
                                    return_ll = TRUE) {
  if (missing(nugget)) nugget = 0
  scmod$r = parm[1]
  v = evaluate(scmod, d) + parm[2] * diag(1/weights)
  cholv = chol(v)
  ploglik_xycholv(x, y, cholv, mu = mu, reml = reml,
                  return_ll = return_ll)
}

ploglik_cmodStd_r_lambda_par3 = function(parm, x, y, d, weights,
                                         scmod, nugget,
                                         mu = NULL, 
                                         reml = FALSE,
                                         return_ll = TRUE) {
  if (missing(nugget)) nugget = 0
  scmod$r = parm[1]
  scmod$par3 = parm[3]
  v = evaluate(scmod, d) + parm[2] * diag(1/weights)
  cholv = chol(v)
  ploglik_xycholv(x, y, cholv, mu = mu, reml = reml,
                  return_ll = return_ll)
}

ploglik_cmodStd_r_lambda_angle_ratio_par3 = function(parm, x, y, d, weights,
                                            scmod, nugget,
                                            mu = NULL, 
                                            reml = FALSE,
                                            return_ll = TRUE) {
  if (missing(nugget)) nugget = 0
  scmod$r = parm[1]
  scmod$angle = parm[3]
  scmod$ratio = parm[4]
  scmod$par3 = parm[5]
  v = evaluate(scmod, d) + parm[2] * diag(1/weights)
  cholv = chol(v)
  ploglik_xycholv(x, y, cholv, mu = mu, reml = reml,
                  return_ll = return_ll)
}

ploglik_cmodStd_r_lambda_angle_par3 = function(parm, x, y, d, weights,
                                               scmod, nugget,
                                               mu = NULL, 
                                               reml = FALSE,
                                               return_ll = TRUE) {
  if (missing(nugget)) nugget = 0
  scmod$r = parm[1]
  scmod$angle = parm[3]
  scmod$par3 = parm[4]
  v = evaluate(scmod, d) + parm[2] * diag(1/weights)
  cholv = chol(v)
  ploglik_xycholv(x, y, cholv, mu = mu, reml = reml,
                  return_ll = return_ll)
}

ploglik_cmodStd_r_lambda_ratio_par3 = function(parm, x, y, d, weights,
                                               scmod, nugget,
                                               mu = NULL, 
                                               reml = FALSE,
                                               return_ll = TRUE) {
  if (missing(nugget)) nugget = 0
  scmod$r = parm[1]
  scmod$ratio = parm[3]
  scmod$par3 = parm[4]
  v = evaluate(scmod, d) + parm[2] * diag(1/weights)
  cholv = chol(v)
  ploglik_xycholv(x, y, cholv, mu = mu, reml = reml,
                  return_ll = return_ll)
}

ploglik_cmodStd_r_lambda_angle_ratio = function(parm, x, y, d, weights,
                                                scmod, nugget,
                                                mu = NULL, 
                                                reml = FALSE,
                                                return_ll = TRUE) {
  if (missing(nugget)) nugget = 0
  scmod$r = parm[1]
  scmod$angle = parm[3]
  scmod$ratio = parm[4]
  v = evaluate(scmod, d) + parm[2] * diag(1/weights)
  cholv = chol(v)
  ploglik_xycholv(x, y, cholv, mu = mu, reml = reml,
                  return_ll = return_ll)
}

ploglik_cmodStd_r_lambda_angle = function(parm, x, y, d, weights,
                                                scmod, nugget,
                                                mu = NULL, 
                                                reml = FALSE,
                                                return_ll = TRUE) {
  if (missing(nugget)) nugget = 0
  scmod$r = parm[1]
  scmod$angle = parm[3]
  v = evaluate(scmod, d) + parm[2] * diag(1/weights)
  cholv = chol(v)
  ploglik_xycholv(x, y, cholv, mu = mu, reml = reml,
                  return_ll = return_ll)
}

ploglik_cmodStd_r_lambda_ratio = function(parm, x, y, d, weights,
                                                scmod, nugget,
                                                mu = NULL, 
                                                reml = FALSE,
                                                return_ll = TRUE) {
  if (missing(nugget)) nugget = 0
  scmod$r = parm[1]
  scmod$ratio = parm[3]
  v = evaluate(scmod, d) + parm[2] * diag(1/weights)
  cholv = chol(v)
  ploglik_xycholv(x, y, cholv, mu = mu, reml = reml,
                  return_ll = return_ll)
}
