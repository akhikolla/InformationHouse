set.seed(55)
rlp_data = geoR::grf(25, cov.pars = c(1, 1),
                     nugget = 0.1,
              messages = FALSE)
y = rlp_data$data
x = matrix(1, nrow = length(y))
d = ganiso_d(rlp_data$coords, coords2 = rlp_data$coords, 
             invert = TRUE)
coords = rlp_data$coords
geoR_ml = geoR::likfit(rlp_data, ini = c(1, 1), 
                       nugget = 0.1, 
                       cov.model = "matern",
                       psiA = pi/8, psiR = 1.5,
                       kappa = 1, fix.nug = FALSE,
                       fix.kappa = FALSE, 
                       messages = FALSE)
geoR_reml = geoR::likfit(rlp_data, ini = c(1, 1), 
                         nugget = 0.1, 
                         cov.model = "matern",
                         psiA = pi/8, psiR = 1.5,
                         kappa = 1,
                         fix.nug = FALSE, messages = FALSE,
                         fix.kappa = FALSE,
                         lik.method = "REML")

aniso_coords = geoR::coords.aniso(coords, aniso.pars = geoR_reml$aniso.pars)
rlp_data_noaniso = rlp_data
rlp_data_noaniso$coords = aniso_coords
geoR_reml_noaniso = geoR::likfit(rlp_data_noaniso, ini = c(1, 1),
                                 nugget = 0.1, fix.nug = TRUE,
                                 fix.psiA = TRUE,
                                 messages = FALSE,
                                 lik.method = "REML")

geoR_ml_loglik = geoR::loglik.GRF(rlp_data, obj.model = geoR_ml)
geoR_reml_loglik = geoR::loglik.GRF(rlp_data, obj.model = geoR_reml)
geoR_reml_noaniso_loglik = geoR::loglik.GRF(rlp_data_noaniso, obj.model = geoR_reml_noaniso)

kc_ml = geoR::krige.control(obj.model = geoR_ml)
geoR_out = geoR::output.control(messages = FALSE)
geoR_ml_beta = geoR::krige.conv(rlp_data, krige = kc_ml,
                                locations = cbind(1, 1),
                                output = geoR_out)$beta.est
kc_reml = geoR::krige.control(obj.model = geoR_reml)
geoR_reml_beta = geoR::krige.conv(rlp_data, krige = kc_reml,
                                  locations = cbind(1, 1),
                                  output = geoR_out)$beta.est

# save output
fpath = system.file("testdata",  package = "gear")
fname = paste(fpath, "/ploglik_cmodStd_r_lambda_par3.rda", sep = "")
save(x, y, d, coords, geoR_ml, geoR_reml,
     geoR_ml_loglik, geoR_reml_loglik,
     geoR_reml_noaniso_loglik,
     geoR_ml_beta, geoR_reml_beta,
     compress = "bzip2",
     file = fname,
     version = 2)

