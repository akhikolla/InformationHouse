### R code from vignette source 'dbscan.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: dbscan.Rnw:132-134
###################################################
options(useFancyQuotes = FALSE)
citation("dbscan")


###################################################
### code chunk number 2: dbscan.Rnw:592-593
###################################################
options(width = 75)


###################################################
### code chunk number 3: dbscan.Rnw:596-607
###################################################
library("dbscan")

set.seed(2)
n <- 400

x <- cbind(
  x = runif(4, 0, 1) + rnorm(n, sd = 0.1),
  y = runif(4, 0, 1) + rnorm(n, sd = 0.1)
  )

true_clusters <- rep(1:4, time = 100)


###################################################
### code chunk number 4: sampleData
###################################################
plot(x, col = true_clusters, pch = true_clusters)


###################################################
### code chunk number 5: kNNdistplot
###################################################
kNNdistplot(x, k = 3)
abline(h=.05, col = "red", lty=2)


###################################################
### code chunk number 6: dbscan.Rnw:641-643
###################################################
res <- dbscan(x, eps = 0.05, minPts = 3)
res


###################################################
### code chunk number 7: dbscanPlot
###################################################
plot(x, col = res$cluster + 1L, pch = res$cluster + 1L)


###################################################
### code chunk number 8: dbscanHullPlot
###################################################
hullplot(x, res)


###################################################
### code chunk number 9: dbscan.Rnw:690-691
###################################################
predict(res, x[1:25,], data = x)


###################################################
### code chunk number 10: dbscan.Rnw:700-702
###################################################
res <- optics(x, eps = 10, minPts = 10)
res


###################################################
### code chunk number 11: dbscan.Rnw:707-708
###################################################
head(res$order, n = 15)


###################################################
### code chunk number 12: opticsReachPlot
###################################################
plot(res)


###################################################
### code chunk number 13: opticsOrder
###################################################
plot(x, col = "grey")
polygon(x[res$order,], )


###################################################
### code chunk number 14: extractDBSCANReachPlot2
###################################################
res <- extractDBSCAN(res, eps_cl = .065)
plot(res)


###################################################
### code chunk number 15: extractDBSCANHullPlot2
###################################################
hullplot(x, res)


###################################################
### code chunk number 16: dbscan.Rnw:774-776
###################################################
res <- extractXi(res, xi = 0.05)
res


###################################################
### code chunk number 17: dbscan.Rnw:781-782
###################################################
res$clusters_xi


###################################################
### code chunk number 18: extractXiReachPlot
###################################################
plot(res)


###################################################
### code chunk number 19: extractXiHullPlot
###################################################
hullplot(x, res)


###################################################
### code chunk number 20: dbscan.Rnw:854-856
###################################################
dend <- as.dendrogram(res)
dend


###################################################
### code chunk number 21: opticsDendrogram
###################################################
plot(dend, ylab = "Reachability dist.", leaflab = "none")


