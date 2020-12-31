set.seed(10)

# generate possible parameter values for covariance models
sigmasq = 1/rgamma(8, 2)
phi = 1/rgamma(8, 2)
error.var = 1/rgamma(8, 3)
micro = c(rep(0, 4), runif(4))
kappa = runif(8, .25, 3)
cm = sample(c("exponential", "spherical", "matern"), 8, replace = TRUE)
n = rpois(8, 150)
np = rpois(8, 300)
nsim = 10000

uk_y = uk_coords = uk_pcoords = uk_pred = uk_mspe = uk_coeff = vector("list", length(sigmasq))
ok_y = ok_coords = ok_pcoords = ok_pred = ok_mspe = ok_coeff = vector("list", length(sigmasq))
sk_y = sk_coords = sk_pcoords = sk_pred = sk_mspe = vector("list", length(sigmasq))
georout_uk_simmeans = georout_ok_simmeans = georout_sk_simmeans = vector("list", length(sigmasq))
georout_uk_simsd = georout_ok_simsd = georout_sk_simsd = vector("list", length(sigmasq))

for (i in 1:8) {
        # generate data
        fields = geoR::grf(100, cov.pars = c(sigmasq[i], phi[i]), cov.model = cm[[i]],
                           kappa = kappa[i], messages = FALSE)

        # extract response
        uk_y[[i]] = fields$data
        # get observed coordinates
        uk_coords[[i]] = fields$coords
        # generate unsampled coordinates
        uk_pcoords[[i]] = matrix(runif(np[i] * 2), ncol = 2)
        acoords = rbind(uk_coords[[i]], uk_pcoords[[i]])

        # create geoR kriging model
        modeli = geoR::krige.control(type = "ok", trend.d = "1st", trend.l = "1st",
                                     cov.model = cm[i], cov.pars = c(sigmasq[i], phi[i]), kappa = kappa[i],
                                     nugget = (error.var[i] + micro[i]), micro.scale = micro[i])

        # decide whether filtered or unfiltered
        output = geoR::output.control(signal = ifelse(i %% 2 == 0, TRUE, FALSE),
                                      messages = FALSE, n.predictive = nsim)
        georout = geoR::krige.conv(fields, loc = acoords, krige = modeli, output = output)
        uk_pred[[i]] = georout$predict
        uk_mspe[[i]] = georout$krige.var
        uk_coeff[[i]] = georout$beta.est
        georout_uk_simmeans[[i]] = georout$mean.simulations
        georout_uk_simsd[[i]] = sqrt(georout$variance.simulations)
}

for (i in 1:8) {
        fields = geoR::grf(100, cov.pars = c(sigmasq[i], phi[i]), cov.model = cm[[i]],
                           kappa = kappa[i], messages = FALSE)
        y = fields$data
        # extract response
        ok_y[[i]] = fields$data
        # get observed coordinates
        ok_coords[[i]] = fields$coords
        # generate unsampled coordinates
        ok_pcoords[[i]] = matrix(runif(np[i] * 2), ncol = 2)
        acoords = rbind(ok_coords[[i]], ok_pcoords[[i]])

        modeli = geoR::krige.control(type = "ok", trend.d = "cte", trend.l = "cte",
                                     cov.model = cm[i], cov.pars = c(sigmasq[i], phi[i]), kappa = kappa[i],
                                     nugget = (error.var[i] + micro[i]), micro.scale = micro[i])

        # decide whether signal model
        output = geoR::output.control(signal = ifelse(i %% 2 == 0, TRUE, FALSE),
                                      messages = FALSE, n.predictive = nsim)

        georout = geoR::krige.conv(fields, loc = acoords, krige = modeli, output = output)
        ok_pred[[i]] = georout$predict
        ok_mspe[[i]] = georout$krige.var
        ok_coeff[[i]] = georout$beta.est
        georout_ok_simmeans[[i]] = georout$mean.simulations
        georout_ok_simsd[[i]] = sqrt(georout$variance.simulations)
}

mus = rnorm(8, 0, sd = 25)
for (i in 1:8) {
        fields = geoR::grf(100, cov.pars = c(sigmasq[i], phi[i]), cov.model = cm[[i]],
                           kappa = kappa[i], messages = FALSE)
        # extract response
        sk_y[[i]] = fields$data
        sk_coords[[i]] = fields$coords
        sk_pcoords[[i]] = matrix(runif(np[i] * 2), ncol = 2)
        acoords = rbind(sk_coords[[i]], sk_pcoords[[i]])

        modeli = geoR::krige.control(type = "sk", trend.d = "cte", trend.l = "cte",
                                     cov.model = cm[i], cov.pars = c(sigmasq[i], phi[i]), kappa = kappa[i],
                                     nugget = (error.var[i] + micro[i]), micro.scale = micro[i], beta = mus[i])

        # decide whether signal model
        output = geoR::output.control(signal = ifelse(i %% 2 == 0, TRUE, FALSE),
                                      messages = FALSE, n.predictive = nsim)
        georout = geoR::krige.conv(fields, loc = acoords, krige = modeli, output = output)
        sk_pred[[i]] = georout$predict
        sk_mspe[[i]] = georout$krige.var
        georout_sk_simmeans[[i]] = georout$mean.simulations
        georout_sk_simsd[[i]] = sqrt(georout$variance.simulations)
}

# save output
fpath = system.file("testdata",  package = "gear")
fname = paste(fpath, "/predict_geolm_cmodStd_data_geoR.rda", sep = "")
save(sigmasq = sigmasq, phi = phi, error.var = error.var,
     micro = micro, kappa = kappa, cm = cm, n = n, np = np,
     angle = angle, ratio = ratio,
     mus = mus,
     uk_y = uk_y, uk_pcoords = uk_pcoords, uk_coords = uk_coords,
     ok_y = ok_y, ok_pcoords = ok_pcoords, ok_coords = ok_coords,
     sk_y = sk_y, sk_pcoords = sk_pcoords, sk_coords = sk_coords,
     uk_pred = uk_pred, uk_mspe = uk_mspe, uk_coeff = uk_coeff,
     ok_pred = ok_pred, ok_mspe = ok_mspe, ok_coeff = ok_coeff,
     sk_pred = sk_pred, sk_mspe = sk_mspe,
     georout_uk_simmeans, georout_uk_simsd,
     georout_ok_simmeans, georout_ok_simsd,
     georout_sk_simmeans, georout_sk_simsd,
     compress = "bzip2",
     file = fname,
     version = 2)
