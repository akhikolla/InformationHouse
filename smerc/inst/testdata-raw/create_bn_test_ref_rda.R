# save output
fpath = system.file("testdata",  package = "smerc")

data(nydf)
coords = as.matrix(nydf[, c("x", "y")])
cases = nydf$cases
pop = nydf$population

set.seed(1)
bnse23 = SpatialEpi::besag_newell(geo = coords,
                                  cases = cases,
                                  population = pop,
                                  k = 23,
                                  alpha = 0.1)
fname = paste(fpath, "/bn_test_ref.rda", sep = "")
save(bnse23,
     compress = "bzip2",
     file = fname,
     version = 2)
