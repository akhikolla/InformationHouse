data(kattegat, package = "geoR")
y = kattegat$data
x = matrix(1, nrow = length(y))
coords = kattegat$coords
d = ganiso_d(kattegat$coords, coords2 = kattegat$coords, 
             invert = TRUE)

geoR_ml = geoR::likfit(kattegat, ini = c(1, 1), 
                       nugget = 0.1, psiA = pi/8,
                       psiR = 1.5,
                       fix.nug = TRUE, messages = FALSE)
geoR_reml = geoR::likfit(kattegat, ini = c(1, 1),
                         nugget = 0.1, psiA = pi/8,
                         psiR = 1.5,
                         fix.nug = TRUE, messages = FALSE,
                         lik.method = "REML")

geoR_ml_loglik = geoR::loglik.GRF(kattegat, obj.model = geoR_ml)
geoR_reml_loglik = geoR::loglik.GRF(kattegat, obj.model = geoR_reml)

aniso_coords = geoR::coords.aniso(coords, aniso.pars = c(pi/8, 1.5))
kattegat_noaniso = kattegat
kattegat_noaniso$coords = aniso_coords
geoR_reml_noaniso = geoR::likfit(kattegat_noaniso, ini = c(1, 1),
                                 nugget = 0.1, fix.nug = TRUE,
                                 messages = FALSE,
                                 lik.method = "REML")
geoR_reml_noaniso_loglik = geoR::loglik.GRF(kattegat_noaniso, obj.model = geoR_reml_noaniso)

# save output
fpath = system.file("testdata",  package = "gear")
fname = paste(fpath, "/ploglik_cmodStd_r_psill.rda", sep = "")
save(x, y, d, coords, geoR_ml, geoR_reml, geoR_ml_loglik,
     geoR_reml_loglik, geoR_reml_noaniso_loglik,
     compress = "bzip2",
     file = fname,
     version = 2)
