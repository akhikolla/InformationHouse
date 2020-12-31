set.seed(120)
rla_data = geoR::grf(50, cov.pars = c(1, 1), nugget = 0.1,
              messages = FALSE)
y = rla_data$data
x = matrix(1, nrow = length(y))
coords = rla_data$coords
d = ganiso_d(rla_data$coords, coords2 = rla_data$coords, 
             invert = TRUE)
geoR_ml = geoR::likfit(rla_data, ini = c(1, 1), 
                       cov.model = "matern",
                       nugget = 0.1, kappa = 1,
                       psiA = pi/8, psiR = 1.5, 
                       fix.nug = FALSE, 
                       fix.psiA = FALSE,
                       fix.psiR = TRUE,
                       fix.kappa = TRUE,
                       messages = FALSE)
geoR_reml = geoR::likfit(rla_data, ini = c(1, 1), 
                         cov.model = "matern",
                         nugget = 0.1, kappa = 1,
                         psiA = pi/8, psiR = 1.5, 
                         fix.nug = FALSE, 
                         fix.psiA = FALSE, 
                         fix.psiR = TRUE,
                         fix.kappa = TRUE, 
                         messages = FALSE,
                         lik.method = "REML")
geoR_ml$info.minimisation.function$par
geoR_reml$info.minimisation.function$par

aniso_coords = geoR::coords.aniso(coords, aniso.pars = geoR_reml$aniso.pars)
rla_data_noaniso = rla_data
rla_data_noaniso$coords = aniso_coords
geoR_reml_noaniso = geoR::likfit(rla_data_noaniso, ini = c(1, 1),
                                 nugget = 0.1, fix.nug = TRUE,
                                 fix.psiA = TRUE,
                                 messages = FALSE,
                                 lik.method = "REML")

geoR_ml_loglik = geoR::loglik.GRF(rla_data, obj.model = geoR_ml)
geoR_reml_loglik = geoR::loglik.GRF(rla_data, obj.model = geoR_reml)
geoR_reml_noaniso_loglik = geoR::loglik.GRF(rla_data_noaniso, obj.model = geoR_reml_noaniso)

kc_ml = geoR::krige.control(obj.model = geoR_ml)
geoR_out = geoR::output.control(messages = FALSE)
geoR_ml_beta = geoR::krige.conv(rla_data, krige = kc_ml,
                                locations = cbind(1, 1),
                                output = geoR_out)$beta.est
kc_reml = geoR::krige.control(obj.model = geoR_reml)
geoR_reml_beta = geoR::krige.conv(rla_data, krige = kc_reml,
                                  locations = cbind(1, 1),
                                  output = geoR_out)$beta.est

# save output
fpath = system.file("testdata",  package = "gear")
fname = paste(fpath, "/ploglik_cmodStd_r_lambda_angle.rda", sep = "")
save(x, y, d, coords, geoR_ml, geoR_reml,
     geoR_ml_loglik, geoR_reml_loglik,
     geoR_reml_noaniso_loglik,
     geoR_ml_beta, geoR_reml_beta,
     compress = "bzip2",
     file = fname,
     version = 2)
