## ----set-up-------------------------------------------------------------------
tree1 <- ape::read.tree(text = '(A, ((B, (C, (D, E))), ((F, G), (H, I))));')
tree2 <- ape::read.tree(text = '(A, ((B, (C, (D, (H, I)))), ((F, G), E)));')

## ----load-package, message=FALSE----------------------------------------------
library('Quartet')

## ----quartet-status-----------------------------------------------------------
statuses <- QuartetStatus(tree1, tree2)

## ----measure-distance---------------------------------------------------------
QuartetDivergence(statuses, similarity = FALSE)

## ----all-metrics--------------------------------------------------------------
SimilarityMetrics(statuses, similarity = TRUE)

## ----visualize-quartets-------------------------------------------------------
VisualizeQuartets(tree1, tree2)

## ----partitions---------------------------------------------------------------
SimilarityMetrics(SplitStatus(tree1, tree2))

## ----multi-trees--------------------------------------------------------------
library('TreeTools', quietly = TRUE, warn.conflicts = FALSE)
oneTree <- CollapseNode(as.phylo(0, 11), 14)
twoTrees <- structure(list(bal = BalancedTree(11), pec = PectinateTree(11)),
                      class = 'multiPhylo')

status <- SharedQuartetStatus(twoTrees, cf = oneTree)
QuartetDivergence(status)

## ----one-to-many--------------------------------------------------------------
forest <- as.phylo(0:5, 11)
names(forest) <- letters[1:6]
status <- SharedQuartetStatus(forest)
QuartetDivergence(status)

## ----many-to-many-------------------------------------------------------------
status <- ManyToManyQuartetAgreement(forest)
QuartetDivergence(status, similarity = FALSE)

## ----pairwise-----------------------------------------------------------------
status <- TwoListQuartetAgreement(forest[1:4], forest[5:6])
QuartetDivergence(status, similarity = FALSE)

## ----in-one-only--------------------------------------------------------------
interestingTree <- as.phylo(42, 7)
referenceTrees <- list(BalancedTree(7), PectinateTree(7))
status <- CompareQuartetsMulti(interestingTree, referenceTrees)

