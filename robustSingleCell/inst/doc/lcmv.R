## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  results = "hide",
  message = F,
  warning = F,
  eval = F
)
library(robustSingleCell)

## ----download_data-------------------------------------------------------
#  library(robustSingleCell)
#  download_LCMV()

## ----initialize----------------------------------------------------------
#  LCMV1 <- initialize.project(datasets = "LCMV1",
#                            origins = "CD44+ cells",
#                            experiments = "Rep1",
#                            data.path = file.path(tempdir(), "LCMV"),
#                            work.path = file.path(tempdir(), "LCMV/LCMV_analysis"))

## ----read_LCMV1----------------------------------------------------------
#  LCMV1 <- read.data(LCMV1, subsample = 500)

## ----preprocess----------------------------------------------------------
#  LCMV1 <- get.variable.genes(LCMV1)
#  exhaustion_markers <- c('Pdcd1', 'Cd244', 'Havcr2', 'Ctla4', 'Cd160', 'Lag3', 'Tigit', 'Cd96')
#  LCMV1 <- add.confounder.variables(LCMV1,
#  	ribosomal.score = ribosomal.score(LCMV1),
#  	mitochondrial.score = mitochondrial.score(LCMV1),
#  	cell.cycle.score = cell.cycle.score(LCMV1),
#  	Exhaustion = controlled.mean.score(LCMV1, exhaustion_markers))

## ----PCA-----------------------------------------------------------------
#  LCMV1 <- PCA(LCMV1, local = T)

## ----cluster-------------------------------------------------------------
#  LCMV1 <- cluster.analysis(LCMV1, local = T)

## ----annotation----------------------------------------------------------
#  types = rbind(
#                  data.frame(type='Tfh',gene=c('Tcf7','Cxcr5','Bcl6')),
#                  data.frame(type='Th1',gene=c('Cxcr6','Ifng','Tbx21')),
#                  data.frame(type='Tcmp',gene=c('Ccr7','Bcl2','Tcf7')),
#                  data.frame(type='Treg',gene=c('Foxp3','Il2ra')),
#                  data.frame(type='Tmem',gene=c('Il7r','Ccr7')),
#                  data.frame(type='CD8',gene=c('Cd8a')),
#                  data.frame(type='CD4', gene = c("Cd4")),
#                  data.frame(type='Cycle',gene=c('Mki67','Top2a','Birc5'))
#  )
#  summarize(LCMV1, local = T)
#  LCMV1_cluster_names <- get.cluster.names(LCMV1, types, min.fold = 1.0, max.Qval = 0.01)
#  LCMV1 <- set.cluster.names(LCMV1, names = LCMV1_cluster_names)
#  summarize(LCMV1, local = T)

## ----plotLCMV1-----------------------------------------------------------
#  canonical_genes <- c("Cd8a", "Cd4", "Mki67", "Foxp3", "Il2ra", "Bcl6",
#                       "Cxcr5", "Cxcr6", "Ifng", "Tbx21", "Id2", "Rora",
#                       "Cxcr3", "Tcf7", "Ccr7", "Cxcr4", "Pdcd1", "Ctla4")
#  plot_simple_heatmap(LCMV1, name = "canonical", markers = canonical_genes, main = "Expression of marker genes")

## ----LCMV_2--------------------------------------------------------------
#  LCMV2 <- initialize.project(datasets = "LCMV2",
#                            origins = "CD44+ cells",
#                            experiments = "Rep2",
#                            data.path = file.path(tempdir(), "LCMV"),
#                            work.path = file.path(tempdir(), "LCMV/LCMV_analysis"))
#  LCMV2 <- read.data(LCMV2, subsample = 500)
#  LCMV2 <- get.variable.genes(LCMV2)
#  LCMV2 <- add.confounder.variables(
#    LCMV2,
#    ribosomal.score = ribosomal.score(LCMV2),
#    mitochondrial.score = mitochondrial.score(LCMV2),
#    cell.cycle.score = cell.cycle.score(LCMV2),
#    Exhaustion = controlled.mean.score(LCMV2, exhaustion_markers))
#  
#  LCMV2 <- PCA(LCMV2, local = T)
#  LCMV2 <- cluster.analysis(LCMV2, local = T) # 0.05 for KNN ratio
#  summarize(LCMV2, local = T)
#  LCMV2_cluster_names <- get.cluster.names(LCMV2, types, min.fold = 1.0, max.Qval = 0.01)
#  LCMV2 <- set.cluster.names(LCMV2, names = LCMV2_cluster_names)
#  summarize(LCMV2, local = T)
#  plot_simple_heatmap(LCMV2, name = "canonical", markers = canonical_genes, main = "Expression of marker genes")

## ----initialize_pooled---------------------------------------------------
#  pooled_env <- initialize.project(datasets = c("LCMV1", "LCMV2"),
#                            origins = c("CD44+ cells", "CD44+ cells"),
#                            experiments = c("Rep1", "Rep2"),
#                            data.path = file.path(tempdir(), "LCMV"),
#                            work.path = file.path(tempdir(), "LCMV/LCMV_analysis"))
#  pooled_env <- read.preclustered.datasets(pooled_env)
#  pooled_env <- add.confounder.variables(
#    pooled_env,
#    ribosomal.score = ribosomal.score(pooled_env),
#    mitochondrial.score = mitochondrial.score(pooled_env),
#    cell.cycle.score = cell.cycle.score(pooled_env),
#    Exhaustion = controlled.mean.score(pooled_env, exhaustion_markers))
#  pooled_env <- PCA(pooled_env, clear.previously.calculated.clustering = F, local = T)
#  summarize(pooled_env, contrast = "datasets", local = T)

## ----pooled--------------------------------------------------------------
#  cluster.similarity <- assess.cluster.similarity(pooled_env)
#  similarity <- cluster.similarity$similarity
#  map <- cluster.similarity$map
#  filtered.similarity <- get.robust.cluster.similarity(
#    pooled_env, similarity, min.sd = qnorm(.9), max.q.val = 0.01, rerun = F
#    )
#  robust.clusters <- sort(unique(c(filtered.similarity$cluster1,
#                                   filtered.similarity$cluster2)))
#  visualize.cluster.cors.heatmaps(pooled_env, pooled_env$work.path,
#                                  filtered.similarity)

## ----summary-------------------------------------------------------------
#  similarity <- filtered.similarity
#  visualize.cluster.similarity.stats(pooled_env, similarity)

## ----robust_markers------------------------------------------------------
#  differential.expression.statistics = get.robust.markers(
#     pooled_env, cluster_group1 = c('LCMV2_Tfh_CD4', 'LCMV2_Tfh_Tcmp_CD4'),
#     cluster_group2 = c('LCMV2_CD8_1', 'LCMV2_CD8_2'),
#     group1_label = 'CD4 T Cells', group2_label = 'CD8 T Cells')

## ----tSNE_overlay--------------------------------------------------------
#  plot_contour_overlay_tSNE(pooled_env, genes = c('Cd4','Cd8a'))

## ------------------------------------------------------------------------
#  plot_pair_scatter(pooled_env, gene1 = 'Cd4', gene2 = 'Cd8a',
#     cluster_group1 = c('LCMV2_Tfh_CD4', 'LCMV2_Tfh_Tcmp_CD4'),
#     cluster_group2 = c('LCMV2_CD8_1','LCMV2_CD8_2'),
#     group1_label = 'CD4 T Cells', group2_label = 'CD8 T Cells')

