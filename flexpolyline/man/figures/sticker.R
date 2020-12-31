#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Name          :sticker.R
# Description   :Creates a hex sticker for the flexpolyline package.
# Author        :Merlin Unterfinger <info@munterfinger.ch>
# Date          :2020-06-10
# Version       :0.1.0
# Usage         :./sticker.R
# Notes         :
# Bash          :5.0.17
# =============================================================================

library(flexpolyline)
library(hexSticker)
library(ggplot2)
library(magrittr)
library(ggecho)

# Path
outfile <- paste0(getwd(), "/man/figures/logo.svg")

# Colors
color1 <- "#393e47"
color2 <- "#65c1c2"
color3 <- "#FFFFFF"

# Data
create_line <- function(n) {
  cbind(
    lng = runif(1, 5, 10) + cumsum(runif(n, -1, 1)),
    lat = runif(1, 45, 48) + cumsum(runif(n, -1, 1)),
    dim = runif(1, 100, 500) + cumsum(runif(n, 0, 1))
  )
}

l = 7
n = 30
enc <-
  data.frame(
    x = seq(0,n-1),
    size = seq(1,n),
    y = seq(0,n-1),
    text = sapply(seq(1,n), function(x) {create_line(l) %>% encode()})
  )

line <- data.frame(
  x = c(0.2, 0.4, 0.5, 0.7, 0.8) * 27 - 0.5,
  y = c(0.35, 0.2, 0.6, 0.3, 0.8) * 27 + 1.5,
  size = c(1.5,2,3,4,5)
)

# Plot
p <-
  ggplot(enc, aes(x = x, y = y, label = text, size = size)) +
  geom_text(color = "#000000", stat = 'echo', size_increment = 0.1, angle = -30) +
  geom_line(data = line, mapping = aes(x = x, y = y, lwd = size*10),
            stat = 'echo', size_increment = 1, color = color2) +
  geom_point(data = line, mapping = aes(x = x, y = y, size = size*10),
             stat = 'echo', size_increment = 1, color = color2) +
  geom_line(data = line, mapping = aes(x = x, y = y, lwd = size*5),
            color = color1) +
  geom_point(data = line, mapping = aes(x = x, y = y, size = size*5),
             color = color1) +
  theme_void() + theme_transparent() +
  scale_size(range = c(0.1, 3)) +
  guides(shape = FALSE, fill = FALSE, color = FALSE, size = FALSE)

# Sticker
sticker(

  # Plot
  p,
  p_size = 12,
  s_x = 1, s_y = 1, s_width = 1.7, s_height = 1.9,

  # Package name
  package = "",
  url = "flexpolyline", u_color = color3,
  u_x = 0.95, u_y = 0.2, u_size = 4, u_angle = 30,

  # Background
  h_fill = color1, h_color = color2,
  white_around_sticker = TRUE,

  # File
  filename = outfile,

)

