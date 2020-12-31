# R demo script
#
#
# Copyright (C) 2015. Javier Royuela del Val
#                     Federico Simmross Wattenberg
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 3 of the License.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.
#
#
#  Javier Royuela del Val.
#  E.T.S.I. Telecomunicacion
#  Universidad de Valladolid
#  Paseo de Belen 15, 47002 Valladolid, Spain.
#  jroyval@lpi.tel.uva.es

# Set parameter vector:
alpha <- 1.25;  # stability index
beta  <- 0.95;  # skew index
sigma <- 2.0;   # scale parameter
mu    <- 0.0;   # location parameter

pars <- c(alpha, beta, sigma, mu)

# abscissae vector
x    <- seq(from=-5, to=100, by=0.1)

# evaluation of pdf and cdf
pdf  <- stable_pdf(x, pars, parametrization=0, tol=1e-12)
cdf  <- stable_cdf(x, pars, parametrization=0, tol=1e-12)

# probabilities vector and quantile function evaluation
p  <- seq(from=0.01, to=0.99, by=0.01)
q  <- stable_q(p, pars, parametrization=0, tol=1e-12)

# Random numbers generation
N   <- 300;
rnd <- stable_rnd(N, pars, seed=1234)

# Parameter estimation with Koutrouvelis method.
# McCulloch estimator will be used for initialization.
pars_init <- stable_fit_init(rnd)
pars_K    <- stable_fit_koutrouvelis(rnd, pars_init)
# pars_ML   <- stable_fit_mle(rnd, pars_init)

# Get the pdf given by the estimated parameters
pdf_init <- stable_pdf(x, pars_init)
pdf_K    <- stable_pdf(x, pars_K)
# pdf_ML   <- stable_pdf(x, pars_ML)

print(pars)
print(pars_init)
print(pars_K)
# print(pars_ML)

# Plot histogram, true pdf and pdf from the estimated parameters
xlim <- c(-5, 30)
ylim <- c(0, 0.15)
breaks <- seq(min(rnd)-1, max(rnd)+1, 0.5);
hist(rnd, breaks=breaks, freq=FALSE, xlim=xlim, ylim=ylim)

lines(x, pdf, col = "red", lwd=2)
lines(x, pdf_init, col = "green", lwd=2)
lines(x, pdf_K, col = "blue", lwd=2)
# lines(x, pdf_ML, col = "yellow", lwd=2)

legend(20, 0.15,
       c("Histogram", "True pdf", "McCulloch", "Koutrouveils"),#,"Max. likelihood"),
       lty = c(1,1,1,1,1),
       lwd = c(2,2,2,2,2),
       col=c("black","red","green","blue"))#,"yellow"))
