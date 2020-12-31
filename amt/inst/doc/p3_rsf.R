## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, message = FALSE)

## -----------------------------------------------------------------------------
library(amt)
data("deer")
deer

## -----------------------------------------------------------------------------
data("sh_forest")
sh_forest

## ---- fig.width=4, fig.height=4-----------------------------------------------
r1 <- random_points(deer)
plot(r1)

## ---- fig.width=4, fig.height=4-----------------------------------------------
r1 <- random_points(deer, n = 100)
plot(r1)

## ---- fig.width=4, fig.height=4-----------------------------------------------
hr <- hr_mcp(deer)
r1 <- random_points(hr, n = 500)
plot(r1)

## ---- fig.width=4, fig.height=4-----------------------------------------------
hr <- hr_mcp(deer)
r1 <- random_points(hr, n = 500, presence = deer)
plot(r1)

## ---- fig.width=4, fig.height=4-----------------------------------------------
hr <- hr_mcp(deer) %>% hr_isopleths() %>% 
  sf::st_buffer(dist =3e4) # add a 30km buffer
r1 <- random_points(hr, n = 500)
plot(r1)

## ---- fig.width=4, fig.height=4-----------------------------------------------
hr <- hr_mcp(deer) %>% hr_isopleths() %>% 
  sf::st_buffer(dist =3e4) # add a 30km buffer
r1 <- random_points(hr, n = 500, presence = deer)
plot(r1)

## -----------------------------------------------------------------------------
rsf1 <- deer %>% random_points() %>% 
  extract_covariates(sh_forest) %>% 
  mutate(forest = sh.forest == 1)

## -----------------------------------------------------------------------------
rsf1 %>% fit_rsf(case_ ~ forest) %>% 
  summary()


## -----------------------------------------------------------------------------
sessioninfo::session_info()

