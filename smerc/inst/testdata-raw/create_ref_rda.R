# save output
fpath = system.file("testdata",  package = "smerc")

# test updated methods
data(nydf)
data(nyw)
coords = with(nydf, cbind(longitude, latitude))
set.seed(1)
uls_test_ref = uls.test(coords = coords, cases = floor(nydf$cases),
               pop = nydf$pop, w = nyw,
               alpha = 0.05, longlat = TRUE,
               nsim = 19, ubpop = 0.5)
fname = paste(fpath, "/uls_test_ref.rda", sep = "")
save(uls_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)

set.seed(1)
rflex_test_ref = rflex.test(coords = coords, cases = floor(nydf$cases),
                  w = nyw, k = 15,
                  pop = nydf$pop, nsim = 19,
                  alpha = 0.1, longlat = TRUE)
fname = paste(fpath, "/rflex_test_ref.rda", sep = "")
save(rflex_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)

set.seed(1)
mlink_test_ref = mlink.test(coords = coords, cases = floor(nydf$cases),
                  pop = nydf$pop, w = nyw,
                  alpha = 0.12, longlat = TRUE,
                  nsim = 19, ubpop = 0.1, ubd = 0.2)
fname = paste(fpath, "/mlink_test_ref.rda", sep = "")
save(mlink_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)

set.seed(1)
flex_test_ref = flex.test(coords = coords, cases = floor(nydf$cases),
                 w = nyw, k = 10,
                 pop = nydf$pop, nsim = 19,
                 alpha = 0.12, longlat = TRUE)
fname = paste(fpath, "/flex_test_ref.rda", sep = "")
save(flex_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)

set.seed(1)
fast_test_ref = fast.test(coords = coords, cases = floor(nydf$cases),
                 pop = nydf$pop,
                 alpha = 0.1, longlat = TRUE,
                 nsim = 19, ubpop = 0.5)
fname = paste(fpath, "/fast_test_ref.rda", sep = "")
save(fast_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)

set.seed(1)
mlf_test_ref = mlf.test(coords = coords, cases = floor(nydf$cases),
                 pop = nydf$pop, w = nyw,
                 alpha = 0.12, longlat = TRUE,
                 nsim = 19, ubpop = 0.1, ubd = 0.5)
fname = paste(fpath, "/mlf_test_ref.rda", sep = "")
save(mlf_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)

set.seed(1)
elliptic_test_ref = elliptic.test(coords = coords,
                      cases = floor(nydf$cases),
                      pop = nydf$pop, ubpop = 0.1,
                      nsim = 19,
                      alpha = 0.12)
fname = paste(fpath, "/elliptic_test_ref.rda", sep = "")
save(elliptic_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)

set.seed(1)
edmst_test_ref = edmst.test(coords = coords, cases = floor(nydf$cases),
                            pop = nydf$pop, w = nyw,
                            alpha = 0.12, longlat = TRUE,
                            nsim = 19, ubpop = 0.1, ubd = 0.2)
fname = paste(fpath, "/edmst_test_ref.rda", sep = "")
save(edmst_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)

set.seed(1)
dmst_test_ref = dmst.test(coords = coords, cases = floor(nydf$cases),
                          pop = nydf$pop, w = nyw,
                          alpha = 0.12, longlat = TRUE,
                          nsim = 19, ubpop = 0.1, ubd = 0.2)
fname = paste(fpath, "/dmst_test_ref.rda", sep = "")
save(dmst_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)

set.seed(1)
dc_test_ref = dc.test(coords = coords, cases = floor(nydf$cases),
                      pop = nydf$pop, w = nyw,
                      alpha = 0.12, longlat = TRUE,
                      nsim = 19, ubpop = 0.1, ubd = 0.2)
fname = paste(fpath, "/dc_test_ref.rda", sep = "")
save(dc_test_ref,
     compress = "bzip2",
     file = fname,
     version = 2)
