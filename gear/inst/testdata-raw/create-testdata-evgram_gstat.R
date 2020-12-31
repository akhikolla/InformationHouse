# save output from gstat package

data(meuse, package = "sp")
sp::coordinates(meuse) = ~ x + y
meuse_df = as.data.frame(meuse)[,1:3]
maxd = max(dist(sp::coordinates(meuse)))/2

# create gstat and geodata objects
gmeuse = gstat::gstat(id = "cadmium", formula = cadmium ~ 1, data = meuse)
# geomeuse = geoR::as.geodata(data.frame(meuse$x, meuse$y, meuse$cadmium))

# for omnidirectional standard semivarigoram
gstat_v1 = gstat::variogram(gmeuse, cutoff = maxd, width = maxd/10)

# for omnidirectional cressie semivarigoram
gstat_v2 = gstat::variogram(gmeuse, cutoff = maxd, width = maxd/10, cressie = TRUE)

# directional semivariograms
gstat_v3 = gstat::variogram(gmeuse, cutoff = maxd, width = maxd/10, alpha = c(22.5 + 0:3*45))
gstat_v4 = gstat::variogram(gmeuse, cutoff = maxd, width = maxd/10, alpha = c(70 + 0:3*45))
gstat_v5 = gstat::variogram(gmeuse, cutoff = maxd, width = maxd/10, alpha = c(0 + 0:3*45))
gstat_v6 = gstat::variogram(gmeuse, cutoff = maxd, width = maxd/10, alpha = c(115 + 0:3*45))
gstat_v7 = gstat::variogram(gmeuse, cutoff = maxd, width = maxd/10, alpha = c(160 + 0:3*45))

# test with trend
gmeuse2 = gstat::gstat(id = "cadmium", formula = cadmium ~ x + y, data = meuse)
gstat_v8 = gstat::variogram(gmeuse2, cutoff = maxd, width = maxd/10)

# cloud semivariograms
# geo_cloud = geoR::variog(geomeuse, option = "cloud", messages = FALSE)
gstat_cloud = gstat::variogram(gmeuse, cutoff = 5000, cloud = TRUE)

# geo_cloud_maxd2000 = geoR::variog(geomeuse, max.dist = 2000,
#                                   option = "cloud", messages = FALSE)
gstat_cloud_maxd2000 = gstat::variogram(gmeuse,
                                        cutoff = 2000,
                                        cloud = TRUE)
# save output
fpath = system.file("testdata",  package = "gear")
fname = paste(fpath, "/evgram_data.rda", sep = "")
save(meuse, meuse_df, maxd,
     gstat_v1, gstat_v2, gstat_v3, gstat_v4,
     gstat_v5, gstat_v6, gstat_v7, gstat_v8,
     # geo_cloud,
     gstat_cloud,
     # geo_cloud_maxd2000,
     gstat_cloud_maxd2000,
     compress = "bzip2",
     file = fname,
     version = 2)
