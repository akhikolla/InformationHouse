## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(surveysd)
population <- demo.eusilc(1, prettyNames = TRUE)
population[, pWeight := 1]
pop_sample <- population[sample(1:.N, floor(.N*0.10)), ]
pop_sample[, pWeight := 10]

## -----------------------------------------------------------------------------
(gender_distribution <- xtabs(pWeight ~ gender, population))
xtabs(pWeight ~ gender, pop_sample)

## -----------------------------------------------------------------------------
pop_sample_c <- ipf(pop_sample, conP = list(gender_distribution), w = "pWeight")

## -----------------------------------------------------------------------------
dim(pop_sample)
dim(pop_sample_c)
setdiff(names(pop_sample_c), names(pop_sample))

## -----------------------------------------------------------------------------
xtabs(calibWeight ~ gender, pop_sample_c)
xtabs(pWeight ~ gender, population)

## ---- fig.align="center", out.width="100%"------------------------------------
xtabs(~ calibWeight + gender, pop_sample_c)

## ---- include = FALSE---------------------------------------------------------
overrepresented_gender <- pop_sample_c[calibWeight < 10, ][1, gender]

## -----------------------------------------------------------------------------
(con_ga <- xtabs(pWeight ~ gender + age, population))
xtabs(pWeight ~ gender + age, pop_sample)

## -----------------------------------------------------------------------------
pop_sample_c2 <- ipf(pop_sample, conP = list(con_ga), w = "pWeight")
xtabs(pWeight ~ gender + age, population)
xtabs(calibWeight ~ gender + age, pop_sample_c2)

## -----------------------------------------------------------------------------
registry_table <- xtabs(pWeight ~ region, population)

## -----------------------------------------------------------------------------
pop_sample_c2 <- ipf(pop_sample, conP = list(con_ga, registry_table), w = "pWeight")
xtabs(pWeight ~ gender + age, population)
xtabs(calibWeight ~ gender + age, pop_sample_c2)
xtabs(pWeight ~ region, population)
xtabs(calibWeight ~ region, pop_sample_c2)

## -----------------------------------------------------------------------------
(conH1 <- xtabs(pWeight ~ hsize + region, data = population[!duplicated(hid)]))
pop_sample_hh <- ipf(pop_sample, hid = "hid", conH = list(conH1), w = "pWeight",
                     bound = 10)
xtabs(calibWeight ~ hsize + region, data = pop_sample_hh[!duplicated(hid)])

## -----------------------------------------------------------------------------
ipf(pop_sample, conP = list(con_ga, registry_table), w = "pWeight",
    verbose = TRUE, epsP = 0.01)
ipf(pop_sample, conP = list(con_ga, registry_table), w = "pWeight",
    verbose = TRUE, epsP = 0.0001)

