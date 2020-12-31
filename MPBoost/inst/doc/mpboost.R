## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, fig.width = 4.8, fig.height = 4.8)

## ------------------------------------------------------------------------
library(MPBoost)
mpboost(N1 = 6, N2 = 6)

## ------------------------------------------------------------------------
mpboost(N1 = 6, N2 = 6, MTI = 3)

## ------------------------------------------------------------------------
mpboost(N1 = 6, N2 = 12)

## ------------------------------------------------------------------------
set.seed(1) ## Only needed for reproductibility
x1 <- mpboost(N1 = 20, N2 = 40, MTI = 4)
x1

## ------------------------------------------------------------------------
lx <- sum(x1 == 1)
ly <- sum(x1 == 2)
ratio <- lx/ly
mti <- 4
plot(cumsum(x1 == 1), cumsum(x1 == 2), xlim = c(0, lx), ylim = c(0, ly),
   xlab = expression(n[1]), ylab = expression(n[2]), lab = c(lx, ly, 7),
   type = "b", pch = 16, panel.first = grid(), col = "red", bty = "n",
   xaxs = "i", yaxs = "i", xpd = TRUE, cex.axis = 0.8)
abline(-mti/ratio, 1/ratio, lty = 3)
abline(mti/ratio, 1/ratio, lty = 3)
abline(0, 1/ratio, lty = 2)

## ------------------------------------------------------------------------
set.seed(11) ## Only needed for reproductibility
x2 <- mpboost(N1 = 20, N2 = 40, MTI = 4)
plot(cumsum(x1 == 1), cumsum(x1 == 2), xlim = c(0, lx), ylim = c(0, ly),
   xlab = expression(n[1]), ylab = expression(n[2]), lab = c(lx, ly, 7),
   type = "b", pch = 16, panel.first = grid(), col = "red", bty = "n",
   xaxs = "i", yaxs = "i", xpd = TRUE, cex.axis = 0.8)
abline(-mti/ratio, 1/ratio, lty = 3)
abline(mti/ratio, 1/ratio, lty = 3)
abline(0, 1/ratio, lty = 2)
lines(cumsum(x2 == 1), cumsum(x2 == 2), type = "b", pch = 16,
   col = rgb(0, 0, 1, alpha = 0.5), xpd = TRUE)

