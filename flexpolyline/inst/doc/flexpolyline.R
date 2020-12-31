## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("flexpolyline")

## ---- eval=FALSE--------------------------------------------------------------
#  remotes::install_github("munterfinger/flexpolyline")

## ----cpp_encode_decode--------------------------------------------------------
library(flexpolyline)

(line <- matrix(
  c(8.69821, 50.10228, 10,
    8.69567, 50.10201, 20,
    8.69150, 50.10063, 30,
    8.68752, 50.09878, 40),
  ncol = 3, byrow = TRUE
))

(encoded <- encode(line))

(decoded <- decode(encoded))

## ----sf_encode_decode---------------------------------------------------------
(sfg <- sf::st_linestring(line, dim = "XYZ"))

(encoded <- encode_sf(sfg))

(decoded <- decode_sf(encoded))

