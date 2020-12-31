## ---- eval=T, echo=FALSE------------------------------------------------------
require(psrwe)
set.seed(1000)

## ---- eval=T, echo=TRUE-------------------------------------------------------
data(ex_dta)
dta_ps <- rwe_ps(ex_dta,
                 v_covs = paste("V", 1:7, sep = ""),
                 v_grp = "Group",
                 cur_grp_level = "current",
                 nstrata = 5)

## ---- echo=TRUE, fig.width=6, fig.height=5------------------------------------
plot(dta_ps, "balance")

## ---- echo=TRUE, fig.width=6, fig.height=5------------------------------------
plot(dta_ps, "ps")

## ---- eval=T, echo=TRUE-------------------------------------------------------
ps_dist   <- rwe_ps_dist(dta_ps)
post_smps <- rwe_ps_powerp(dta_ps,
                           total_borrow = 40,
                           v_distance   = ps_dist$Dist[1:dta_ps$nstrata],
                           outcome_type = "binary",
                           v_outcome    = "Y")


## ---- echo=TRUE, fig.width=6, fig.height=5------------------------------------
traceplot(post_smps$stan_rst, pars = c("theta", "thetas"))

## ---- eval=T, echo=TRUE-------------------------------------------------------
summary(post_smps)

## ---- eval=T, echo=TRUE-------------------------------------------------------
ps_borrow <- rwe_ps_borrow(total_borrow = 40, ps_dist)
rst_cl    <- rwe_ps_cl(dta_ps, v_borrow = ps_borrow, v_outcome = "Y")
summary(rst_cl)

## ---- eval=T, echo=TRUE-------------------------------------------------------
data(ex_dta_rct)
dta_ps_2arm <- rwe_ps(ex_dta_rct,
                      v_covs = paste("V", 1:7, sep = ""),
                      v_grp = "Group",
                      cur_grp_level = "current",
                      nstrata = 5)

rst_2arm <- rwe_ps_cl2arm(dta_ps_2arm,
                          v_arm = "Arm",
                          trt_arm_level = 1,
                          outcome_type = "continuous",
                          v_outcome = "Y",
                          total_borrow = 40)

print(rst_2arm)

