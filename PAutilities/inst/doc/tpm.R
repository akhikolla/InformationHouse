## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(PAutilities)

## ----vectors-------------------------------------------------------------
set.seed(100)
algorithm <- (sample(1:100)%%2)
criterion <- (sample(1:100)%%2)

## ----TPM7----------------------------------------------------------------
  transitions <- get_transition_info(
    predictions = algorithm,
    references = criterion,
    window_size = 7
  )

## ----TPM7plot, fig.width=7, fig.height=5---------------------------------
  plot(transitions)

## ----TPM7summarize-------------------------------------------------------
  summarized1 <- summary(transitions)

## ----slots---------------------------------------------------------------
  summarized1@result


  # or:
  # slot(summarized1, "result")

## ----pipes, message=FALSE, warning=FALSE---------------------------------
  suppressPackageStartupMessages(
    library(magrittr, quietly = TRUE, verbose = FALSE)
  )

  summarized <-
    get_transition_info(algorithm, criterion, 7) %>%
    summary(.) %>%
    slot("result")

## ----occasions-----------------------------------------------------------
  # Here I'm exploiting seed changes to get different values from the same code
  # I used previously

  algorithm2 <- (sample(1:100)%%2)
  criterion2 <- (sample(1:100)%%2)

## ----add-----------------------------------------------------------------
  summarized2 <-
    get_transition_info(algorithm2, criterion2, 7) %>%
    summary(.)

  # Store the result of addition (another S4 summaryTransition object)
  
    added <- summarized1 + summarized2
  
  # Now view the result
    
    added@result
    

    # or:
    # slot(added, "result")

## ----subtract------------------------------------------------------------
  subtracted <- summarized1 - summarized2

  subtracted$differences

## ----spurious------------------------------------------------------------
  curve <- spurious_curve(trans = transitions, thresholds = 5:10)
  class(curve)
  sapply(curve, class)

## ----spur_plot, fig.width=7, fig.height=5--------------------------------
  par(mar=rep(3,4))
  plot(curve)

