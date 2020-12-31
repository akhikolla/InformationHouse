fpath = system.file("testdata",  "cmod_std_data.rda", package = "gear")
load(fpath)
     
# check with geoR
mod = geoR_cov_params$mod
evar = geoR_cov_params$evar
fvar = geoR_cov_params$fvar
r = geoR_cov_params$r
par3 = geoR_cov_params$par3
psill = geoR_cov_params$psill

count = 0
gear_cmod_stda = vector("list", length(evar) * length(mod))
for (i in 1:length(mod)) {
  for (j in 1:length(evar))   {
    count = count + 1
    cmod = cmod_std(model = mod[i], 
                   par3 = par3[j],
                   psill = psill[j], 
                   r = r[j], 
                   evar = evar[j],
                   fvar = fvar[j])
    gear_cmod_stda[[count]] = evaluate(cmod, d)
  }
}

test_that("evaluate.cmodStd is accurate (geoR)", {
  expect_equivalent(geoR_cov, gear_cmod_stda)
})

# comparison with spam
mod = spam_cov_params$mod
evar = spam_cov_params$evar
fvar = spam_cov_params$fvar
r = spam_cov_params$r
par3 = spam_cov_params$par3
psill = spam_cov_params$psill

count = 0
gear_cmod_stdb = vector("list", length(evar) * length(mod))
for (i in 1:length(mod)) {
  for (j in 1:length(evar)) {
    count = count + 1
    cmod = cmod_std(model = mod[i], 
                   par3 = par3[j],
                   psill = psill[j], 
                   r = r[j], 
                   evar = evar[j],
                   fvar = fvar[j])
    gear_cmod_stdb[[count]] = evaluate(cmod, d)
  }
}

test_that("evaluate.cmodStd is accurate (spam)", {
  expect_equivalent(spam_cov, gear_cmod_stdb)
})

# comparison with spam (matern)
# comparison with spam
mod = spam_matern_params$mod
evar = spam_matern_params$evar
fvar = spam_matern_params$fvar
r = spam_matern_params$r
par3 = spam_matern_params$par3
psill = spam_matern_params$psill

count = 0
gear_cmod_stdc = vector("list", length(evar) * length(mod))
for (i in 1:length(mod)) {
  for (j in 1:length(evar)) {
    count = count + 1
    cmod = cmod_std(model = mod[i], 
                   par3 = par3[j],
                   psill = psill[j], 
                   r = r[j], 
                   evar = evar[j],
                   fvar = fvar[j])
    gear_cmod_stdc[[count]] = evaluate(cmod, d)
  }
}

test_that("evaluate.cmodStd model = 'matern' is accurate (spam)", {
  expect_equivalent(spam_matern, gear_cmod_stdc)
})
