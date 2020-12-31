fpath = system.file("testdata", "predict_geolm_cmodStd_data_geoR.rda", package = "gear")
load(fpath)

############ test whether predict.geolm_cmodStd results are correct
# check accuracy of method (compare to geoR results)
results.diff.uk = matrix(99, ncol = 2, nrow = 16)
# store universal kriging results
uk.results = matrix(99, nrow = 8, ncol = 3)
nsim = 0 #10000
gearout_uk_simmeans = gearout_ok_simmeans = gearout_sk_simmeans = vector("list", length(sigmasq))
gearout_uk_simsd = gearout_ok_simsd = gearout_sk_simsd = vector("list", length(sigmasq))

for (i in 1:8) {
  y = uk_y[[i]]
  coords = uk_coords[[i]]
	acoords = rbind(coords, uk_pcoords[[i]])

	# create data frames needed for gear
	data = data.frame(x1 = coords[,1], x2 = coords[,2], y = y)
	newdata = data.frame(x1 = acoords[,1], x2 = acoords[,2])

	# decide whether signal/filtered or unfiltered model
	cmod = cmod_std(model = cm[i], psill = sigmasq[i], r = phi[i],
	                par3 = kappa[i],
	                evar = ifelse(i %% 2, 0, error.var[i]),
	                fvar = ifelse(i %% 2, micro[i] + error.var[i], micro[i]))
	# gear geostatistic model
	gearmod = geolm(y ~ x1 + x2, data = data, mod = cmod,
	                coordnames = c("x1", "x2"))
	gearout = predict(gearmod,  newdata = newdata,
	                  return_type = "data.frame", nsim = nsim,
	                  dmethod = "svd")

	uk.results[i, 1] = max(abs(range(uk_pred[[i]] - gearout$pred)))
	uk.results[i, 2] = max(abs(range(uk_mspe[[i]] - gearout$mspe)))
	uk.results[i, 3] = max(abs(range(uk_coeff[[i]] - gearmod$coeff)))
	if (nsim > 0) {
	  gearout_uk_simmeans[[i]] = rowMeans(gearout$sim)
	  gearout_uk_simsd[[i]] = apply(gearout$sim, 1, sd)
	}
}

test_that("all predict.geolm_cmodStd anisotropic uk calculations are correct", {
  expect_true(max(uk.results) < 1e-10)
})

ok.results = matrix(99, nrow = 8, ncol = 3)

for (i in 1:8) {
  y = ok_y[[i]]
  coords = ok_coords[[i]]
  acoords = rbind(coords, ok_pcoords[[i]])
	# create df for geolm
	data = data.frame(x1 = coords[,1], x2 = coords[,2], y = y)
	newdata = data.frame(x1 = acoords[,1], x2 = acoords[,2])

	# decide whether signal/filtered or unfiltered model
	cmod = cmod_std(model = cm[i], psill = sigmasq[i], r = phi[i],
	                par3 = kappa[i],
	                evar = ifelse(i %% 2, 0, error.var[i]),
	                fvar = ifelse(i %% 2, micro[i] + error.var[i], micro[i]))

	# create geolm for gear package
	gearmod = geolm(y ~ 1, data = data, mod = cmod,
	                coordnames = c("x1", "x2"))

	gearout = predict(gearmod,  newdata = newdata,
	                  return_type = "data.frame", nsim = nsim,
	                  dmethod = "svd")

	ok.results[i, 1] = max(abs(range(ok_pred[[i]] - gearout$pred)))
	ok.results[i, 2] = max(abs(range(ok_mspe[[i]] - gearout$mspe)))
	ok.results[i, 3] = max(abs(range(ok_coeff[[i]] - gearmod$coeff)))
	if (nsim > 0) {
	  gearout_ok_simmeans[[i]] = rowMeans(gearout$sim)
	  gearout_ok_simsd[[i]] = apply(gearout$sim, 1, sd)
	}
}

test_that("all predict.geolm_cmodStd anisotropic ok calculations are correct", {
  expect_true(max(ok.results) < 1e-10)
})

sk.results = matrix(99, nrow = 8, ncol = 2)

for (i in 1:8) {
  y = sk_y[[i]]
  coords = sk_coords[[i]]
  acoords = rbind(coords, sk_pcoords[[i]])
  # create df for geolm
  data = data.frame(x1 = coords[,1], x2 = coords[,2], y = y)
  newdata = data.frame(x1 = acoords[,1], x2 = acoords[,2])

  # decide whether signal/filtered or unfiltered model
  cmod = cmod_std(model = cm[i], psill = sigmasq[i], r = phi[i],
                  par3 = kappa[i],
                  evar = ifelse(i %% 2, 0, error.var[i]),
                  fvar = ifelse(i %% 2, micro[i] + error.var[i], micro[i]))

	# create geolm for gear package
	gearmod = geolm(y ~ 0, data = data,
	                 coordnames = c("x1", "x2"),
	                 mod = cmod, mu = mus[i])

	gearout = predict(gearmod,  newdata = newdata,
	                  return_type = "data.frame", nsim = nsim,
	                  dmethod = "svd")

	sk.results[i, 1] = max(abs(range(sk_pred[[i]] - gearout$pred)))
	sk.results[i, 2] = max(abs(range(sk_mspe[[i]] - gearout$mspe)))
	if (nsim > 0) {
	  gearout_sk_simmeans[[i]] = rowMeans(gearout$sim)
	  gearout_sk_simsd[[i]] = apply(gearout$sim, 1, sd)
	}
}

test_that("all predict.geolm_cmodStd anisotropic sk calculations are correct", {
  expect_true(max(sk.results) < 1e-10)
})

if (nsim > 0) {
  pdf("~/Dropbox/predict_geolm_cmodStd_condsim.pdf")
  par(mfrow = c(3, 2))
  for (i in 1:8) {
    plot(georout_uk_simmeans[[i]],
         gearout_uk_simmeans[[i]])
    title(paste("mean uk", i))
    abline(0, 1)
    plot(georout_uk_simsd[[i]],
         gearout_uk_simsd[[i]])
    title(paste("sd uk", i))
    abline(0, 1)

    plot(georout_ok_simmeans[[i]],
         georout_ok_simmeans[[i]])
    title(paste("mean ok", i))
    abline(0, 1)
    plot(georout_ok_simsd[[i]],
         gearout_ok_simsd[[i]])
    title(paste("sd ok", i))
    abline(0, 1)

    plot(georout_sk_simmeans[[i]],
         georout_sk_simmeans[[i]])
    title(paste("mean sk", i))
    abline(0, 1)
    plot(georout_sk_simsd[[i]],
         georout_sk_simsd[[i]])
    title(paste("sd sk", i))
    abline(0, 1)
  }
  dev.off()
}
