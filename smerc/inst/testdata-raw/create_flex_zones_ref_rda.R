# save output
fpath = system.file("testdata",  package = "smerc")

data(nydf)
coords = as.matrix(nydf[, c("x", "y")])
data(nyw)
flex5_zones_ref = flex.zones(coords, nyw, k = 5)

fname = paste(fpath, "/flex5_zones_ref.rda", sep = "")
save(flex5_zones_ref,
     compress = "bzip2",
     file = fname,
     version = 2)
