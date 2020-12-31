# save output from gstat and geoR package
set.seed(88)
n = 10
np = 5
coords = matrix(runif(n * 2), nrow = n)
pcoords = matrix(runif(np * 2), nrow = np)

# generate different combinations of covariance model
mod = c("exponential", "matern", "spherical", "gaussian")
evar = c(0, .5, 0, .5)
fvar = c(0, 0, .5, .5)
r = runif(4)
par3 = runif(4, 0, 2)
psill = rgamma(4, 2, 2)
gangle = 22 + 0:3 * 45
gratio = runif(4)
geoR_cov_aniso_params = list(mod = mod,
                         evar = evar,
                         fvar = fvar,
                         r = r,
                         par3 = par3,
                         psill = psill,
                         gangle = gangle,
                         gratio = gratio)

count = 0
geoR_cov_aniso = vector("list", length(evar) * length(mod))
# evaluate covariance for different combinations
for (i in seq_along(mod)) {
        for (j in seq_along(evar)) {
                count = count + 1
                # shrink and rotate coordinates for anisotropy
                coords_aniso = geoR::coords.aniso(coords, aniso.pars = c(gangle[j] * pi/180, 1/gratio[j]))
                # evaluate covariance
                geoR_cov_aniso[[count]] = 
                        geoR::varcov.spatial(coords_aniso, cov.model = mod[i], 
                                             kappa = par3[j],
                                             nugget = (evar[j] + fvar[j]), 
                                             cov.pars = c(psill[j], r[j] * gratio[j]))$varcov
        }
}

geoR_coords_aniso = vector("list", length(gratio))
geoR_pcoords_aniso = vector("list", length(gratio))
for (j in seq_along(gratio)) {
        # shrink and rotate coordinates for anisotropy
        geoR_coords_aniso[[j]] = geoR::coords.aniso(coords, aniso.pars = c(gangle[j] * pi/180, 1/gratio[j]))
        geoR_pcoords_aniso[[j]] = geoR::coords.aniso(pcoords, aniso.pars = c(gangle[j] * pi/180, 1/gratio[j]))
}

# save output
fpath = system.file("testdata",  package = "gear")
fname = paste(fpath, "/cmod_std_aniso_data.rda", sep = "")
save(n, coords, geoR_coords_aniso, geoR_cov_aniso, geoR_cov_aniso_params,
     np, pcoords, geoR_pcoords_aniso,
     compress = "bzip2",
     file = fname,
     version = 2)
