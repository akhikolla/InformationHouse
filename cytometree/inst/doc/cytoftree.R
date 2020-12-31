## ----knitrsetup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(tidy = TRUE)
knitr::knit_hooks$set(small.mar = function(before, options, envir) {
    if (before) par(mar = c(0, 0, 0, 0))  # no margin
})

## ---- message = FALSE, warning = FALSE----------------------------------------
library(cytometree)
data(IMdata)
dim(IMdata)
colnames(IMdata)

## -----------------------------------------------------------------------------
zero_proportion <- apply(IMdata[,-c(1,2)], 
                         MARGIN = 2, 
                         FUN = function(x){round(prop.table(table(x==0))["TRUE"]*100,2)})
zero_proportion

## -----------------------------------------------------------------------------
num_col <- c(3:ncol(IMdata))

tree <- CytofTree(M = IMdata,
                  minleaf = 1,
                  t = 0.1,
                  verbose = FALSE,
                  force_first_markers = c("(Ir191)Dd_DNA1",
                                          "(Ir193)Dd_DNA2",
                                          "Cell_length",
                                          "(Ce140)Dd_Bead",
                                          "(In115)Dd_Dead"),
                  transformation = "asinh",
                  num_col = num_col)

max(tree$labels)

## -----------------------------------------------------------------------------
annot <- Annotation(tree, plot = FALSE, K2markers = colnames(IMdata))
annot$combinations[1:5,]

## -----------------------------------------------------------------------------
phenotypes <- list()
phenotypes[["CD4+"]] <- rbind(c("(Ir191)Dd_DNA1", 1), c("(Ir193)Dd_DNA2", 1), 
                              c("Cell_length", 0), c("(Ce140)Dd_Bead", 0), 
                              c("(In115)Dd_Dead", 0), c("(Sm154)Dd_CD14", 0), 
                              c("(Er166)Dd_CD33", 0), c("(Nd150)Dd_CD3", 1), 
                              c("(Nd143)Dd_CD4", 1))

phenotypes[["CD8+"]] <- rbind(c("(Ir191)Dd_DNA1", 1), c("(Ir193)Dd_DNA2", 1), 
                              c("Cell_length", 0), c("(Ce140)Dd_Bead", 0), 
                              c("(In115)Dd_Dead", 0), c("(Sm154)Dd_CD14", 0), 
                              c("(Er166)Dd_CD33", 0), c("(Nd150)Dd_CD3", 1), 
                              c("(Nd144)Dd_CD8", 1))

pheno_result <- RetrievePops(annot, phenotypes = phenotypes)

# CD4+
pheno_result$phenotypesinfo[[1]]

# CD8+
pheno_result$phenotypesinfo[[2]]

## ---- echo = FALSE------------------------------------------------------------
automating <- c(pheno_result$phenotypesinfo[[1]]$proportion,
                pheno_result$phenotypesinfo[[2]]$proportion)
manual <- c(0.1824389, 0.06523925)
resu <- rbind(manual, automating)
rownames(resu) <- c("Manual Gating", "Automating Gating")
colnames(resu) <- c("CD4+", "CD8+")
knitr::kable(resu, digits = 3)

