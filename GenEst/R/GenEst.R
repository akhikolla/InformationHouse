#' @useDynLib GenEst, .registration = TRUE
#' @importFrom corpus print.corpus_frame
#' @importFrom DT dataTableOutput renderDataTable datatable
#' @importFrom graphics abline axis box hist legend lines mtext par plot
#'   plot.new points polygon rect text title
#' @importFrom grDevices dev.off devAskNewPage png rgb colors
#' @importFrom htmltools a br code div em h3 h4 HTML img p tags
#' @importFrom lubridate is.Date
#' @importFrom MASS fitdistr
#' @importFrom Rcpp sourceCpp
#' @importFrom shiny actionButton checkboxInput checkboxGroupInput column 
#'   conditionalPanel downloadButton downloadHandler fileInput fluidRow 
#'   h3 h4 h5 htmlOutput isolate mainPanel navbarPage numericInput 
#'   observeEvent outputOptions p plotOutput radioButtons reactiveValues 
#'   removeNotification renderPlot renderText renderUI runApp 
#'   selectizeInput shinyApp shinyAppDir showNotification sidebarLayout
#'   sidebarPanel tabPanel tabsetPanel textOutput updateNumericInput
#'   updateSelectizeInput updateTabsetPanel  updateNavbarPage
#' @importFrom shinyjs inlineCSS reset useShinyjs 
#' @importFrom stats .getXlevels approxfun as.formula delete.response density
#'   formula median model.matrix na.omit optim pexp pgamma plnorm pnorm pweibull
#'   qnorm quantile reformulate rnorm runif terms update.formula weighted.mean
#'   rbinom dnorm
#' @importFrom survival strata
#' @importFrom utils combn packageDescription read.csv read.csv2 write.csv
#'  write.table
#' @importFrom hellno data.frame as.data.frame
#'

#' @title Generalized estimation of mortality
#'
#' @description This package is designed to analyze searcher efficiency,
#'   carcass persistence, search schedule, and carcass observation data for
#'   the estimation of bird and bat mortality at wind and solar power 
#'   facilities.
#'
#' @name GenEst
#'
#' @section Information:
#' \code{browseVignettes("GenEst")}\cr
#' \code{packageDescription("GenEst")}\cr
#' \code{disclaimerUSGS()}\cr
#' \code{disclaimerWEST()}\cr
#'
#' @section Data sets:
#' \code{\link{mock}}\cr
#' \code{\link{wind_cleared}}\cr
#' \code{\link{wind_RP}}\cr
#' \code{\link{wind_RPbat}}\cr
#' \code{\link{solar_powerTower}}\cr
#' \code{\link{solar_PV}}\cr
#' \code{\link{solar_trough}}\cr
#'
#' @section Main command-line functions:
#'  \describe{
#'   \item{\code{\link{pkm}, \link{cpm}, \link{dwpm}}}{estimate searcher efficiency
#'    (\code{pk}), carcass persistence (\code{cp}), and (\code{dwp}) parameters}
#'   \item{\code{\link{estM}}}{estimate mortality given \code{pkm}, \code{cpm}
#'    and data}
#'   \item{\code{\link{calcSplits}}}{split mortality estimates by subcategories}
#'   \item{\code{plot}}{S3 function for \code{\link[=plot.pkm]{pkm}},
#'    \code{\link[=plot.pkmSet]{pkmSet}}, \code{\link[=plot.cpm]{cpm}},
#'    \code{\link[=plot.cpmSet]{cpmSet}}, \code{\link[=plot.estM]{estM}},
#'    \code{\link[=plot.splitFull]{splitFull}},
#'    \code{\link[=plot.splitSummary]{splitSummary}},
#'    \code{\link[=plot.gGeneric]{gGeneric}}, and
#'    \code{\link[=plot.gGenericSize]{gGenericSize}} objects}
#'   \item{\code{\link{transposeSplits}}}{transpose 2-d splits}
#'   \item{\code{summary}}{S3 function for \code{\link[=summary.estM]{estM}},
#'    \code{\link[=summary.splitFull]{splitFull}},
#'    \code{\link[=summary.gGeneric]{gGeneric}},
#'    \code{\link[=summary.gGenericSize]{gGenericSize}} objects}
#'   \item{\code{\link{aicc}}}{S3 function for extracting models' AICc values
#'     from \code{\link{pkm}}, \code{\link[=pkm]{pkmSet}},
#'    \code{\link[=pkm]{pkmSize}}, \code{\link[=pkm]{pkmSetSize}}, 
#'    \code{\link{cpm}}, \code{\link[=cpm]{cpmSet}},
#'    \code{\link[=cpm]{cpmSize}}, and \code{\link[=cpm]{cpmSetSize}} objects}
#'   \item{\code{\link{desc}}}{Calculate descriptive statistics for a fitted
#'     CP model}
#'   \item{\code{\link{estgGeneric}}}{estimate
#'    detection probability (g) for given searcher efficiency and carcass
#'    persistence model}
#'   \item{\code{runGenEst()}}{start the GUI}
#' }
#' @section Potentially useful calculation functions:
#' \code{\link{rpk}}, \code{\link{qpk}}, \code{\link{rcp}}, \code{\link{rdwp}}\cr
#' \code{\link{estg}}, \code{\link{calcg}}\cr
#' \code{\link{ppersist}}, \code{\link{SEsi}}\cr
#' \code{\link{alogit}}, \code{\link{logit}}\cr
#' \code{\link{pkLogLik}}, \code{\link{cpLogLik}}\cr
#' \code{\link{calcRate}}, \code{\link{calcTsplit}}, 
#' \code{\link{ltranspose}}\cr
#' \code{\link{refMod}}\cr
#' \code{\link{countCarcs}}\cr
#' \code{\link{simpleMplot}}\cr
#'
#' @section Potentially useful editing functions:
#' \code{\link{estgGenericSize}}\cr
#' \code{\link{prepSS}}\cr
#' \code{\link{averageSS}}\cr
#' \code{\link{tidyModelSetCP}}\cr
#' \code{\link{tidyModelSetSE}}\cr
#' \code{\link{checkDate}}\cr
#' \code{\link{dateCols}}\cr
#' \code{\link{dateToDay}}\cr
#' \code{\link{defineUnitCol}}\cr
#' \code{\link{dlModTabSE}}\cr
#' \code{\link{prettyModTabCP}}\cr
#' \code{\link{prettyModTabSE}}\cr
#' \code{\link{prettySplitTab}}\cr
#'
#' @section Other functions:
#' \code{\link{trimSetSize}}\cr
#' \code{\link{combinePreds}}\cr
#' \code{\link{combinePredsAcrossModels}}\cr
#' \code{\link{pkmSetSizeFailRemove}}\cr
#' \code{\link{pkmSetFailRemove}}\cr
#' \code{\link{cpmSetSizeFailRemove}}\cr
#' \code{\link{cpmSetSizeFail}}\cr
#' \code{\link{cpmSetFailRemove}}\cr
#' \code{\link{CO_DWP}}\cr
#' \code{\link{CPcols}}\cr
#' \code{\link{cpmCPCellPlot}}\cr
#' \code{\link{cpmFail}}\cr
#' \code{\link{cpmSetFail}}\cr
#' \code{\link{cpmSetSpecCPCellPlot}}\cr
#' \code{\link{DWPCols}}\cr
#' \code{\link{expandModelSetCP}}\cr
#' \code{\link{obsCols_SE}}\cr
#' \code{\link{pkmFail}}\cr
#' \code{\link{pkmSetAllFail}}\cr
#' \code{\link{pkmSetFail}}\cr
#' \code{\link{pkmSetSizeFail}}\cr
#' \code{\link{plotCPCells}}\cr
#' \code{\link{plotCPFigure}}\cr
#' \code{\link{plotCPHeader}}\cr
#' \code{\link{predsCols}}\cr
#' \code{\link{SEsi_left}}\cr
#' \code{\link{SEsi_right}}\cr
#' \code{\link{SEsi0}}\cr
#' \code{\link{sizeCols}}\cr
#' \code{\link{matchCells}}\cr

#' \code{checkComponents}\cr
#' \code{checkSpecificModelCP}\cr
#' \code{checkSpecificModelSE}\cr
#' \code{combinePredsAcrossModels}\cr
#' \code{CPdistOptions}\cr
#' \code{obsCols_fta}\cr
#' \code{obsCols_ltp}\cr
#' \code{pkmParamPlot}\cr
#' \code{pkmSECellPlot}\cr
#' \code{pkmSet}\cr
#' \code{pkmSetSpecParamPlot}\cr
#' \code{pkmSetSpecSECellPlot}\cr
#' \code{plotSEBoxPlots}\cr
#' \code{plotSEBoxTemplate}\cr
#' \code{plotSECells}\cr
#' \code{plotSEFigure}\cr
#' \code{plotSEHeader}\cr
#' \code{prepPredictors}\cr
#' \code{readCSV}\cr
#' \code{removeCols}\cr

#' @section Internal functions (not exported):
#' \code{_GenEst_calcRateC}\cr
#' \code{_GenEst_calcTsplitC}\cr
#' \code{calcRateC}\cr
#' \code{calcTsplitC}\cr
#' \code{aboutContent}\cr
#' \code{aboutPanel}\cr
#' \code{analysisPanel}\cr
#' \code{b}\cr
#' \code{big}\cr
#' \code{cButtonStyle}\cr
#' \code{center}\cr
#' \code{classText}\cr
#' \code{clearNotifications}\cr
#' \code{CPMainPanel}\cr
#' \code{CPPanel}\cr
#' \code{downloadCPFig}\cr
#' \code{estText}\cr
#' \code{CPSidebar}\cr
#' \code{initialReactiveValues}\cr
#' \code{createvtext}\cr
#' \code{dataDownloadWidget}\cr
#' \code{dataInputPanel}\cr
#' \code{dataInputSidebar}\cr
#' \code{dataInputWidget}\cr
#' \code{dataTabPanel}\cr
#' \code{disclaimersContent}\cr
#' \code{disclaimersPanel}\cr
#' \code{downloadCPFig}\cr
#' \code{downloadCPMod}\cr
#' \code{downloadCPMres}\cr
#' \code{downloadSEmod}\cr
#' \code{downloadgFig}\cr
#' \code{downloadMFig}\cr
#' \code{downloadSEFig}\cr
#' \code{downloadsPanel}\cr
#' \code{downloadTable}\cr
#' \code{eventReaction}\cr
#' \code{GeneralInputSidebar}\cr
#' \code{GeneralInputsPanel}\cr
#' \code{GenEstAcknowledgements}\cr
#' \code{GenEstAuthors}\cr
#' \code{GenEstGUIauthors}\cr
#' \code{GenEstInlineCSS}\cr
#' \code{GenEstLicense}\cr
#' \code{GenEstLogos}\cr
#' \code{GenEstShinyJS}\cr
#' \code{GenEstUI}\cr
#' \code{GenEstServer}\cr
#' \code{gettingStartedContent}\cr
#' \code{gettingStartedPanel}\cr
#' \code{gMainPanel}\cr
#' \code{gPanel}\cr
#' \code{gSidebar}\cr
#' \code{helpPanel}\cr
#' \code{initialOutput}\cr
#' \code{kFixedWidget}\cr
#' \code{kFixedWidgetHeader}\cr
#' \code{kFixedWidgetRow}\cr
#' \code{li}\cr
#' \code{loadedDataPanel}\cr
#' \code{MMainPanel}\cr
#' \code{modelInputWidget}\cr
#' \code{modelOutputPanel}\cr
#' \code{modelOutputWidget}\cr
#' \code{modelRunWidget}\cr
#' \code{modelSelectionWidget}\cr
#' \code{modelSelectionWidgetHeader}\cr
#' \code{modelSelectionWidgetRow}\cr
#' \code{modelSetCells}\cr
#' \code{modelSetModelCells}\cr
#' \code{modelSetModelPredictors}\cr
#' \code{modelSetPredictors}\cr
#' \code{modNamePaste}\cr
#' \code{modNameSplit}\cr
#' \code{MPanel}\cr
#' \code{msgFracNote}\cr
#' \code{msgList}\cr
#' \code{msgModDone}\cr
#' \code{msgModFail}\cr
#' \code{msgModPartialFail}\cr
#' \code{msgModRun}\cr
#' \code{msgModSENobs}\cr
#' \code{msgModWarning}\cr
#' \code{msgSampleSize}\cr
#' \code{msgSplitFail}\cr
#' \code{msgSSavgFail}\cr
#' \code{msgSSinputFail}\cr
#' \code{MSidebar}\cr
#' \code{navbar}\cr
#' \code{ol}\cr
#' \code{pickSizeclass}\cr
#' \code{plotNA}\cr
#' \code{prepSizeclassText}\cr
#' \code{preTextMaker}\cr
#' \code{reaction}\cr
#' \code{reactionMessageDone}\cr
#' \code{reactionMessageRun}\cr
#' \code{reNULL}\cr
#' \code{reVal}\cr
#' \code{SEboxes}\cr
#' \code{SEcols}\cr
#' \code{selectData}\cr
#' \code{selectedDataPanel}\cr
#' \code{SEMainPanel}\cr
#' \code{SEPanel}\cr
#' \code{SEpanel}\cr
#' \code{SESidebar}\cr
#' \code{setkNeed}\cr
#' \code{setNotSuspending}\cr
#' \code{small}\cr
#' \code{splitButtonWidget}\cr
#' \code{style}\cr
#' \code{trimSetSize}\cr
#' \code{u}\cr
#' \code{ul}\cr
#' \code{update_input}\cr
#' \code{update_output}\cr
#' \code{update_rv}\cr
#' \code{updateColNames_size}\cr
#' \code{updateSizeclasses}\cr
#' \code{updatesizeCol}\cr
#' \code{widgetMaker}\cr
#'
#' @docType package
#'
#' @keywords package
#'
NULL
