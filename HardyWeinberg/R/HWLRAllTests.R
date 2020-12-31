HWLRAllTests <- function(x,y) {
  lr.AF <- HWLRtest(x,y,scene.null="S1",scene.alt="S6")
  lr.CF <- HWLRtest(x,y,scene.null="S3",scene.alt="S6")
  lr.DF <- HWLRtest(x,y,scene.null="S4",scene.alt="S6")
  lr.BC <- HWLRtest(x,y,scene.null="S2",scene.alt="S3")
  lr.AB <- HWLRtest(x,y,scene.null="S1",scene.alt="S2")
  lr.EF <- HWLRtest(x,y,scene.null="S5",scene.alt="S6")
  lr.DE <- HWLRtest(x,y,scene.null="S4",scene.alt="S5")
  pvals <- c(lr.AF$pval,lr.CF$pval,lr.DF$pval,lr.BC$pval,lr.AB$pval,
             lr.EF$pval,lr.DE$pval)
  names(pvals) <- c("AF","CF","DF","BC","AB","EF","DE")
  return(pvals)
}
