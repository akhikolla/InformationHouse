data(kattegat, package = "geoR")
ml = geoR::likfit(kattegat, ini = c(0.5, 0.5), fix.nug = TRUE,
                  messages = FALSE)
v = geoR::varcov.spatial(coords = kattegat$coords,
                   cov.model = ml$cov.model,
                   cov.pars = ml$cov.pars,
                   nugget = ml$nugget,
                   messages = FALSE)$varcov
y = kattegat$data
x = matrix(1, nrow = length(y))
cholv = chol(v)
geoR_ll = geoR::loglik.GRF(kattegat, obj.model = ml)
geoR_ll_reml = geoR::loglik.GRF(kattegat, obj.model = ml,
                                method = "REML")
# save output
fpath = system.file("testdata",  package = "gear")
fname = paste(fpath, "/ll_xycholv_geoR.rda", sep = "")
save(x, y, cholv, geoR_ll, geoR_ll_reml,
     compress = "bzip2",
     file = fname,
     version = 2)
