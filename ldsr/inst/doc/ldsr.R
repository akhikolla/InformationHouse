## ----setup, include=FALSE, warning=FALSE, message=FALSE-----------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
lds <- ldsr:::NPlds
cv <- ldsr:::NPcv
u <- v <- t(ldsr::NPpc)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(ldsr)       # This package
library(data.table) # Data wrangling
library(ggplot2)    # Plotting
library(patchwork)  # Arranging multiple plots

## -----------------------------------------------------------------------------
head(NPannual)

## -----------------------------------------------------------------------------
NPpc

## ---- eval=FALSE--------------------------------------------------------------
#  # Setup doFuture as the parallel computing backend
#  doFuture::registerDoFuture()
#  future::plan(future::multiprocess)
#  # Learn LDS
#  u <- v <- t(NPpc)
#  lds <- LDS_reconstruction(NPannual, u, v, start.year = 1200, num.restarts = 20)

## -----------------------------------------------------------------------------
lds$theta

## ---- fig.width=8, fig.height=4.5---------------------------------------------
p1 <- ggplot(lds$rec[year %in% NPannual$year]) +
  geom_ribbon(aes(year, ymin = Ql, ymax = Qu), fill = 'gray90') +
  geom_line(aes(year, Q, colour = 'LDS')) +
  geom_line(aes(year, Qa, colour = 'Observation'), data = NPannual) +
  scale_colour_manual(name = NULL, values = c('black', 'darkorange')) +
  labs(x = NULL, y = 'Mean annual flow [m\u00b3/s]') +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank())

p2 <- ggplot(lds$rec[year %in% NPannual$year]) +
  geom_ribbon(aes(year, ymin = Xl, ymax = Xu), fill = 'gray90') +
  geom_line(aes(year, X)) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(x = 'Year', y = 'Catchment state [-]')

p1 / p2 + plot_layout(heights = c(1, 0.6))

## ---- fig.width=8, fig.height=4.5---------------------------------------------
p1 <- ggplot(lds$rec) +
  geom_ribbon(aes(year, ymin = Ql, ymax = Qu), fill = 'gray90') +
  geom_hline(aes(yintercept = mean(Q)), colour = 'salmon') +
  geom_line(aes(year, Q)) +
  labs(x = NULL, y = 'Mean annual flow [m\u00b3/s]') +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank())

p2 <- ggplot(lds$rec) +
  geom_ribbon(aes(year, ymin = Xl, ymax = Xu), fill = 'gray90') +
  geom_hline(yintercept = 0, colour = 'salmon') +
  geom_line(aes(year, X)) +
  theme_classic() +
  labs(x = 'Year', y = 'Catchment state [-]')

p1 / p2 + plot_layout(heights = c(1, 0.6))

## -----------------------------------------------------------------------------
set.seed(100)
Z <- make_Z(NPannual$Qa, nRuns = 30, frac = 0.25, contiguous = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  cv <- cvLDS(NPannual, u, v, start.year = 1600, num.restarts = 20, Z = Z)

## -----------------------------------------------------------------------------
cv$metrics

## -----------------------------------------------------------------------------
# Build PCR model
pcr <- PCR_reconstruction(NPannual, NPpc, start.year = 1200)
# Cross validate with the same folds as before
cvpcr <- cvPCR(NPannual, NPpc, start.year = 1200, Z = Z, metric.space = 'original')

## -----------------------------------------------------------------------------
rbind(lds = cv$metrics, pcr = cvpcr$metrics)

## ---- fig.width=8.5, fig.height=2.5-------------------------------------------
dt1 <- as.data.table(cvpcr$metrics.dist)
dt1[, model := 'PCR']
dt2 <- as.data.table(cv$metrics.dist)
dt2[, model := 'LDS']
dt <- rbind(dt1, dt2)
dt <- melt(dt, id.vars = 'model',  variable.name = 'metric')

ggplot(dt, aes(model, value)) +
  geom_boxplot() +
  stat_summary(geom = 'point', fun = mean, colour = 'red') +
  facet_wrap(vars(metric), scales = 'free', nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(strip.background = element_rect(fill = 'gray90', colour = NA))

## ---- fig.width=8, fig.height=4.5---------------------------------------------
p1 <- ggplot(lds$rec[year %in% NPannual$year]) +
  geom_ribbon(aes(year, ymin = Ql, ymax = Qu), fill = 'gray90') +
  geom_line(aes(year, Q, colour = 'LDS', linetype = 'LDS')) +
  geom_line(aes(year, Q, colour = 'PCR', linetype = 'PCR'), data = pcr$rec[year %in% NPannual$year]) +
  geom_line(aes(year, Qa, colour = 'Observation', linetype = 'Observation'), data = NPannual) +
  scale_colour_manual(name = NULL, values = c('black', 'darkorange', 'black')) +
  scale_linetype_manual(name = NULL, values = c(1, 1, 2)) +
  labs(x = NULL, y = 'Mean annual flow [m\u00b3/s]') +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank())

p2 <- ggplot(lds$rec[year %in% NPannual$year]) +
  geom_ribbon(aes(year, ymin = Xl, ymax = Xu), fill = 'gray90') +
  geom_line(aes(year, X)) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(x = 'Year', y = 'Catchment state [-]')

p1 / p2 + plot_layout(heights = c(1, 0.6))

## ---- fig.width=8, fig.height=4.5---------------------------------------------
p1 <- ggplot(lds$rec) +
  geom_ribbon(aes(year, ymin = Ql, ymax = Qu), fill = 'gray90') +
  geom_line(aes(year, Q, colour = 'LDS', linetype = 'LDS')) +
  geom_line(aes(year, Q, colour = 'PCR', linetype = 'PCR'), data = pcr$rec) +
  scale_colour_manual(name = NULL, values = c('black', 'steelblue')) +
  scale_linetype_manual(name = NULL, values = c(1, 2)) +
  labs(x = NULL, y = 'Mean annual flow [m\u00b3/s]') +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank())

p2 <- ggplot(lds$rec) +
  geom_ribbon(aes(year, ymin = Xl, ymax = Xu), fill = 'gray90') +
  geom_line(aes(year, X)) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(x = 'Year', y = 'Catchment state [-]')

p1 / p2 + plot_layout(heights = c(1, 0.6))

## -----------------------------------------------------------------------------
set.seed(100)
reps <- LDS_rep(lds$theta, u, v, years = lds$rec$year, mu = mean(log(NPannual$Qa)))

## ---- fig.width=7, fig.height=4.5---------------------------------------------
# Plot streamflow
p <- ggplot(reps) +
  geom_line(aes(year, simQ, group = rep), colour = 'gray80') +
  geom_line(aes(year, Q), data = lds$rec, colour = 'black') +
  labs(x = 'Year',
       y = 'Q [m\u00b3/s]') +
  theme_classic()
# Plot catchment state
q <- ggplot(reps) +
  geom_line(aes(year, simX, group = rep), colour = 'gray80') +
  geom_line(aes(year, X), data = lds$rec, colour = 'black') +
  labs(x = 'Year',
       y = 'Catchment state [-]') +
  theme_classic()

p / q + plot_layout(heights = c(1, 0.6))

