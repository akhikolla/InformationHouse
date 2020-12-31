set.seed(10)
rl_data = geoR::grf(50, cov.pars = c(1, 1), nugget = 0.1,
                          messages = FALSE)
y = rl_data$data
x = matrix(1, nrow = length(y))
d = as.matrix(dist(rl_data$coords))
coords = rl_data$coords
geoR_ml = geoR::likfit(rl_data, ini = c(1, 1), 
                       fix.nug = FALSE, messages = FALSE,
                       lik.method = "ML")
geoR_reml = geoR::likfit(rl_data, ini = c(1, 1), 
                         fix.nug = FALSE, messages = FALSE,
                         lik.method = "REML")
# save output
fpath = system.file("testdata",  package = "gear")
fname = paste(fpath, "/ploglik_cmodStd_r_lambda.rda", sep = "")
save(rl_data, x, y, d, coords, geoR_ml, geoR_reml,
     compress = "bzip2",
     file = fname,
     version = 2)
