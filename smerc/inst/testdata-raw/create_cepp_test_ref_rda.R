# save output
fpath = system.file("testdata",  package = "smerc")

# test updated methods
data(nydf)
data(nyw)
coords = with(nydf, cbind(longitude, latitude))
set.seed(1)
cepp_test_ref = cepp.test(coords = coords,
                          cases = floor(nydf$cases),
                          pop = nydf$pop,
                          nstar = 1000, alpha = 0.1)
fname = paste(fpath, "/cepp_test_ref.rda", sep = "")
save(cepp_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)
