fpath = system.file("testdata",  "cmod_std_aniso_data.rda", package = "gear")
load(fpath)

# check with geoR
mod = geoR_cov_aniso_params$mod
evar = geoR_cov_aniso_params$evar
fvar = geoR_cov_aniso_params$fvar
r = geoR_cov_aniso_params$r
par3 = geoR_cov_aniso_params$par3
psill = geoR_cov_aniso_params$psill
gangle = geoR_cov_aniso_params$gangle
gratio = geoR_cov_aniso_params$gratio

count = 0
gear_cov_aniso = vector("list", length(evar) * length(mod))
for (i in 1:length(mod)) {
  for (j in 1:length(evar))   {
    count = count + 1
    cmod = cmod_std(model = mod[i],
                    par3 = par3[j],
                    psill = psill[j],
                    r = r[j],
                    evar = evar[j],
                    fvar = fvar[j],
                    angle = gangle[j],
                    ratio = gratio[j],
                    invert = TRUE)
    d = ganiso_d(coords, coords, invert = TRUE)
    gear_cov_aniso[[count]] = evaluate(cmod, d)
  }
}

test_that("evaluate.cmodStd aniso is accurate (geoR)", {
  expect_equivalent(geoR_cov_aniso, gear_cov_aniso)
})


count = 0
gear_crosscov_aniso = vector("list", length(gratio))
geoR_crosscov_aniso = vector("list", length(gratio))
for (j in seq_along(gratio))   {
  cmod = cmod_std(model = "exponential",
                  par3 = par3[j],
                  psill = psill[j],
                  r = r[j],
                  evar = 0,
                  fvar = 0,
                  angle = gangle[j],
                  ratio = gratio[j],
                  invert = TRUE)
  d = ganiso_d(coords, pcoords, invert = TRUE)
  gear_crosscov_aniso[[j]] = evaluate(cmod, d)
  geoR_d = geodist(geoR_coords_aniso[[j]], geoR_pcoords_aniso[[j]])
  geoR_crosscov_aniso[[j]] = psill[j] * exp(-geoR_d/(r[j] * gratio[j]))
}

test_that("evaluate.cmodStd aniso cross-covariance is accurate (geoR)", {
  expect_equivalent(geoR_crosscov_aniso, gear_crosscov_aniso)
})
