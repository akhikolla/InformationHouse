### R code from vignette source 'ptsuite_vignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: ptsuite_vignette.Rnw:90-95
###################################################
library("ptsuite")
ps.options(pointsize = 12)
options(width = 72)
options(SweaveHooks = list(fig = function() par(las = 1)))
options(prompt = "R> ", continue = "+")


###################################################
### code chunk number 2: ptsuite_vignette.Rnw:404-407
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
head(d)


###################################################
### code chunk number 3: ptsuite_vignette.Rnw:426-429 (eval = FALSE)
###################################################
## set.seed(1234)
## d <- generate_pareto(100000, 1.2, 3)
## pareto_qq_test(d)


###################################################
### code chunk number 4: ptsuite_vignette.Rnw:433-436 (eval = FALSE)
###################################################
## set.seed(1234)
## exp_data <- rexp(100000, 5)
## pareto_qq_test(exp_data)


###################################################
### code chunk number 5: ptsuite_vignette.Rnw:490-493
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
pareto_test(d)


###################################################
### code chunk number 6: ptsuite_vignette.Rnw:496-499
###################################################
set.seed(1234)
exp_data <- rexp(100000, 5)
pareto_test(exp_data)


###################################################
### code chunk number 7: ptsuite_vignette.Rnw:519-521
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)


###################################################
### code chunk number 8: ptsuite_vignette.Rnw:525-526
###################################################
alpha_mle(d)


###################################################
### code chunk number 9: ptsuite_vignette.Rnw:530-531
###################################################
alpha_mle(d, FALSE)


###################################################
### code chunk number 10: ptsuite_vignette.Rnw:535-536
###################################################
alpha_mle(d, TRUE, 0.05)


###################################################
### code chunk number 11: ptsuite_vignette.Rnw:550-553
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
alpha_hills(d, 8000, FALSE)


###################################################
### code chunk number 12: ptsuite_vignette.Rnw:556-557
###################################################
alpha_hills(d, 5000, TRUE)


###################################################
### code chunk number 13: ptsuite_vignette.Rnw:560-561
###################################################
alpha_hills(d, 100000, FALSE)


###################################################
### code chunk number 14: ptsuite_vignette.Rnw:585-588
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
alpha_ls(d)


###################################################
### code chunk number 15: ptsuite_vignette.Rnw:602-605
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
alpha_percentile(d)


###################################################
### code chunk number 16: ptsuite_vignette.Rnw:620-623
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
alpha_modified_percentile(d)


###################################################
### code chunk number 17: ptsuite_vignette.Rnw:637-640
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
alpha_geometric_percentile(d)


###################################################
### code chunk number 18: ptsuite_vignette.Rnw:654-657
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
alpha_wls(d)


###################################################
### code chunk number 19: ptsuite_vignette.Rnw:671-674
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
alpha_moment(d)


###################################################
### code chunk number 20: ptsuite_vignette.Rnw:688-691
###################################################
set.seed(1234)
d <- generate_pareto(100000, 1.2, 3)
generate_all_estimates(d)


