## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----sfg_z--------------------------------------------------------------------
library(flexpolyline)
library(sf)

coords <- matrix(
    c(8.69821, 50.10228, 10.11111,
      8.69567, 50.10201, 20.22222,
      8.69150, 50.10063, 30.33333,
      8.68752, 50.09878, 40.44444),
    ncol = 3, byrow = TRUE
)

(sfg_z <- st_linestring(coords, dim = "XYZ"))

(sfg_enc_z <- encode_sf(sfg_z))

decode_sf(sfg_enc_z)

## ----sfg_m--------------------------------------------------------------------
(sfg_m <- st_linestring(coords, dim = "XYM"))

(sfg_enc_m <- encode_sf(sfg_m))

decode_sf(sfg_enc_m)

## ----sfc----------------------------------------------------------------------
(sfc <- st_as_sfc(
  lapply(seq(1, 5), function(x) {
    st_linestring(coords[, 1:2] + runif(1, -1, 1), dim = "XY")
  }),
  crs = 4326
))

(sfc_enc <- encode_sf(sfc))

decode_sf(sfc_enc, crs = 4326)

## ----sf-----------------------------------------------------------------------
(sf <- st_as_sf(
  data.frame(
    name = c("A", "B", "C", "D", "E"),
    color = sample(c("red", "green", "blue"), 5, replace = TRUE),
    geometry = sfc
  )
))

(sf_enc <- encode_sf(sf))

decode_sf(sf_enc, crs = 4326)

