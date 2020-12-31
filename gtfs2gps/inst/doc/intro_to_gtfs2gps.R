## ---- eval = FALSE------------------------------------------------------------
#  install.packages("gtfs2gps")

## -----------------------------------------------------------------------------
library("gtfs2gps")
sao <- read_gtfs(system.file("extdata/saopaulo.zip", package ="gtfs2gps"))
names(sao)
sao$trips

## -----------------------------------------------------------------------------
library(magrittr)
object.size(sao) %>% format(units = "Kb")
sao_small <- gtfs2gps::filter_by_shape_id(sao, c(51338, 51956, 51657))
object.size(sao_small) %>% format(units = "Kb")

## ----sao_small_shapes_sf, message = FALSE-------------------------------------
sao_small_shapes_sf <- gtfs2gps::gtfs_shapes_as_sf(sao_small)
sao_small_stops_sf <- gtfs2gps::gtfs_stops_as_sf(sao_small)
plot(sf::st_geometry(sao_small_shapes_sf))
plot(sf::st_geometry(sao_small_stops_sf), pch = 20, col = "red", add = TRUE)
box()

## ---- message = FALSE---------------------------------------------------------
write_gtfs(sao_small, "sao_small.zip")

## ---- message = FALSE---------------------------------------------------------
  sao_gps <- gtfs2gps("sao_small.zip", parallel = FALSE, spatial_resolution = 50)
  head(sao_gps)

## ---- message = FALSE---------------------------------------------------------
  sao_gps60 <- sao_gps[1:100, ]
  
  # points
  sao_gps60_sfpoints <- gps_as_sfpoints(sao_gps60)
  
  # linestring
  sao_gps60_sflinestring <- gps_as_sflinestring(sao_gps60)

  # plot
  plot(sf::st_geometry(sao_gps60_sfpoints), pch = 20)
  plot(sf::st_geometry(sao_gps60_sflinestring), col = "blue", add = TRUE)
  box()

## ---- message = FALSE---------------------------------------------------------
poa <- system.file("extdata/poa.zip", package ="gtfs2gps")

poa_gps <- gtfs2gps(poa, parallel = FALSE, spatial_resolution = 50)

poa_gps_sflinestrig <- gps_as_sfpoints(poa_gps)

plot(sf::st_geometry(poa_gps_sflinestrig[1:200,]))

box()

## ----speed, echo = FALSE, message = FALSE-------------------------------------
knitr::include_graphics("https://github.com/ipeaGIT/gtfs2gps/blob/master/man/figures/speed.PNG")

