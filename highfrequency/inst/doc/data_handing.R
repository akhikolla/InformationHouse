## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(highfrequency)
head(sampleTDataRawMicroseconds)

## ----eval = TRUE--------------------------------------------------------------
summary(sampleTDataRawMicroseconds[, c("DT", "SIZE", "PRICE")])

## ----eval = TRUE--------------------------------------------------------------
tDataCleaned <- tradesCleanup(tDataRaw = sampleTDataRawMicroseconds, exchange = "N")

## ----eval = TRUE--------------------------------------------------------------
tDataCleaned$report

summary(tDataCleaned$tData[, c("DT", "SIZE", "PRICE")])

## ----eval = TRUE--------------------------------------------------------------
qDataCleaned <- quotesCleanup(qDataRaw = sampleQDataRawMicroseconds, exchange = "N")

## ----eval = TRUE--------------------------------------------------------------
qDataCleaned$report

summary(qDataCleaned$qData[, c("DT", "OFR", "OFRSIZ", "BID", "BIDSIZ", "MIDQUOTE")])

## ----eval = TRUE--------------------------------------------------------------
tqDataCleaned <- tradesCleanupUsingQuotes(tData = tDataCleaned$tData[as.Date(DT) == "2018-01-02"], 
                                           qData = qDataCleaned$qData[as.Date(DT) == "2018-01-02"])
tqDataCleaned

