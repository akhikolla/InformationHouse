## ------------------------------------------------------------------------
library(ProFound)
library(FITSio)

## ---- eval=FALSE---------------------------------------------------------
#  image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

## ---- fig.width=8, fig.height=8, eval=FALSE------------------------------
#  profound=profoundProFound(image, magzero=30, plot=TRUE)
#  profound=profoundProFound(image, magzero=30, plot=TRUE, skycut = 2)
#  profound=profoundProFound(image, magzero=30, plot=TRUE, tolerance = 10)
#  profound=profoundProFound(image, magzero=30, plot=TRUE, box = 50)
#  profound=profoundProFound(image, magzero=30, plot=TRUE, SBdilate = 1)
#  profound=profoundProFound(image, magzero=30, plot=TRUE, roughpedestal = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  GALEX_NUV=readFITS(system.file("extdata", 'GALEX_NUV.fits', package="magicaxis"))
#  VST_r=readFITS(system.file("extdata", 'VST_r.fits', package="magicaxis"))
#  VISTA_K=readFITS(system.file("extdata", 'VISTA_K.fits', package="magicaxis"))
#  
#  # Warp to common WCS:
#  GALEX_NUV_VST=magwarp(GALEX_NUV, VST_r$hdr)$image
#  VISTA_K_VST=magwarp(VISTA_K, VST_r$hdr)$image
#  
#  multi=profoundMultiBand(inputlist=list(GALEX_NUV_VST, VST_r$imDat, VISTA_K_VST),
#  magzero=c(20.08,0,30), detectbands=c('r','K'), multibands=c('NUV','r','K'))

## ---- fig.width=8, fig.height=8, eval=FALSE------------------------------
#  profound=profoundProFound(VST_r, roughpedestal=TRUE, SBdilate=1, plot=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  fixedRGB=profoundSegimFix(list(R=VISTA_K_VST, G=VST_r$imDat, B=GALEX_NUV_VST), segim=profound$segim)

## ---- eval=FALSE, fig.width=8, fig.height=8------------------------------
#  profoundSegimPlot(image=VST_r$imDat, segim=fixedRGB$segim)

## ---- eval=FALSE, fig.width=8, fig.height=8------------------------------
#  profound2=profoundProFound(VST_r, segim=fixedRGB$segim, plot=TRUE)

