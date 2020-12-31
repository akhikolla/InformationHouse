## ---- eval = FALSE-------------------------------------------------------
#  library(propr)
#  data(caneToad.counts)
#  data(caneToad.groups)

## ---- echo = FALSE-------------------------------------------------------
library(propr)
data(caneToad.groups)
data(top.counts)
data(top.lr)
best <- new("propr")
best@counts <- as.data.frame(top.counts)
best@logratio <- as.data.frame(top.lr)
best@matrix <- propr:::lr2rho(top.lr)
best@metric <- "rho"
best <- best[">", .995]

## ---- eval = FALSE-------------------------------------------------------
#  keep <- apply(caneToad.counts, 2, function(x) sum(x >= 10) >= 10)
#  rho <- propr(caneToad.counts, metric = "rho", select = keep)

## ---- eval = FALSE-------------------------------------------------------
#  best <- rho[">", .995]

## ---- dpi = 66, fig.width = 8, fig.height = 8, results = "hide", fig.show = "hold", fig.keep = "last"----
plot(best)

## ---- eval = FALSE-------------------------------------------------------
#  dendrogram(best)

## ---- eval = FALSE-------------------------------------------------------
#  best <- simplify(best)

## ---- eval = FALSE-------------------------------------------------------
#  pca(best, group = caneToad.groups)

## ---- eval = FALSE-------------------------------------------------------
#  snapshot(best)

## ---- dpi = 66, fig.width = 8, fig.height = 8, results = "hide", message = FALSE----
clusts <- prism(best, k = 5)

## ---- dpi = 66, fig.width = 8, fig.height = 8, results = "hide", message = FALSE----
clusts <- bokeh(best, k = 5)

## ---- results = "hide"---------------------------------------------------
sub <- subset(best, select = (clusts == 2))

## ---- dpi = 66, fig.width = 8, fig.height = 8, results = "hide", fig.keep = "last"----
pca(sub, group = caneToad.groups)

## ------------------------------------------------------------------------
transcripts <- colnames(sub@logratio)

