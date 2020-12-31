## ----Load package, include=FALSE----------------------------------------------
require('Quartet')
require('ape') # for plot.phylo

## ----Three-four-taxon-trees, echo=FALSE, cache=TRUE, fig.width=6, fig.asp=1/3, out.width='66%', fig.align='center'----
par(mfrow = c(1, 3), mar = c(0.5, 1, 0.5, 1), cex = 1)
plot(ape::read.tree(text = '((A, B), (C, D));'),
     tip.color = Ternary::cbPalette8[c(1, 4, 7, 5)], font = 2)
plot(ape::read.tree(text = '((A, C), (B, D));'),
     tip.color = Ternary::cbPalette8[c(1, 7, 4, 5)], font = 2)
plot(ape::read.tree(text = '((A, D), (C, B));'), 
     tip.color = Ternary::cbPalette8[c(1, 5, 7, 4)], font = 2)

## ----Plot-a-quartet, echo=FALSE, cache=TRUE, fig.asp=1.3/3, fig.width=6, out.width='80%', fig.align='center'----
par(mfrow = c(1, 3))
suppressWarnings(RNGversion("3.5.0")) # Stopgap until R 3.6.0 is widely available 
set.seed(7)
trees7 <- lapply(logical(3), function (X) {
    tr <- ape::rtree(7, br = NULL)
    tr$edge.length <- rep(1, 12)
    tr$tip.label <- LETTERS[1:7]
    tr
  })
Quartet::PlotQuartet(trees7, LETTERS[1:4], cex = 1.4, font = 2)

## ----mean-proportion----------------------------------------------------------
round(vapply(4:20, function (nTip) {
 trees <- lapply(rep(nTip, 10), TreeTools::RandomTree)
 s <- ManyToManyQuartetAgreement(trees)[, , 's']
 results <- s[lower.tri(s)] / choose(nTip, 4)
 c(mean(results), sd(results))
}, c(mean = 0, sd = 0)), 3)

## ----example-tree-------------------------------------------------------------
tree_a <- ape::read.tree(text = "((1, 2), (3, (4, 5)));")

## ----plot-trees, fig.height=1.6, fig.width=2, echo=FALSE----------------------
par(mar = rep(0.3, 4))
plot(tree_a)

## ----none-in-common-----------------------------------------------------------
tree_b <- ape::read.tree(text = "((1, 5), (3, (2, 4)));")

## ----plot-tree-b, fig.height=1.6, fig.width=2, echo=FALSE---------------------
par(mar = rep(0.3, 4))
plot(tree_b)

## ----tree-c-------------------------------------------------------------------
tree_c <- ape::read.tree(text="((1, 2), ((3, 6), (4, 5)));")

## ----Add-tip-6-to-Tree-C, fig.height=1.6, fig.width=2, echo=FALSE-------------
par(mar = rep(0.3, 4))
plot(tree_c, tip.color = c(1,1,1,2,1,1))

## ----Adding-tip-6-to-Tree-B-duplicates-a-quartet, fig.height=5, fig.width=2.5, echo=FALSE----
PlotApeTree <- function (text, quartet) {
  orig <- TreeTools::RenumberTips(tree_c, as.character(1:6))
  tree <- ape::read.tree(text = text)
  PlotQuartet(list(orig, TreeTools::RenumberTips(tree, as.character(1:6))), quartet, overwritePar = FALSE, cex = 0.9)
}

par(mfrow = c(7, 2), mar = rep(0.4, 4), cex = 0.9)
PlotApeTree("(((1, 6), 5), (3, (2, 4)));", c(1, 6, 4, 5))
PlotApeTree("((1, 5), (3, ((2, 6), 4)));", c(2, 6, 4, 5))
PlotApeTree("((1, 5), ((3, 6), (2, 4)));", c(3, 6, 4, 5))
PlotApeTree("((1, 5), (3, (2, (4, 6))));", c(4, 6, 1, 2))
PlotApeTree("((1, (5, 6)), (3, (2, 4)));", c(5, 6, 1, 2))
PlotApeTree("(((1, 5), 6), (3, (2, 4)));", c(1, 5, 3, 6))
PlotApeTree("((1, 5), (3, ((2, 4), 6)));", c(4, 2, 3, 6))

