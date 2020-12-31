#' @title GenEst app download funtions
#'
#' @description Handle downloading of app data and figures
#'
#' @param rv the reactive values list
#'
#' @param input list of shiny input parameters
#'
#' @param sc size class
#'
#' @param split logical indicator to use the split or not
#'
#' @param filename the name for the file writing out
#'
#' @param tablename the name of the table in the rv list
#'
#' @param csvformat format for .csv file: "" or NULL for comma-separated, 2
#'  for semi-colon separated
#'
#' @return download handler functions
#'
#' @name app_download_functions
NULL

#' @rdname app_download_functions
#'
downloadCPFig <- function(rv){
  downloadHandler(filename = "CP_fig.png",
      content = function(file){
        png(file, height = rv$figH_CP, width = rv$figW_CP, units = "px")
        plot(rv$modSet_CP, specificModel = rv$outCPdlsfig)
        dev.off()
      }
  )
}

#' @rdname app_download_functions
#'
downloadCPmod <- function(rv, input){
  downloadHandler(filename = paste0(rv$filename_CP, "_CPmod.txt"),
    content = function(file){
      cat(paste0("Data: ", rv$filename_CP),
        c("\nLast present: ", rv$ltp, "\nFirst absent: ", rv$fta, "\n"),
        file = file, sep = " ")
      cat("Confidence level: ", 100 * rv$CL, "%\n\n", sep = "",
        file = file, append = TRUE)
      selected_mods <- list()
      modChoices <- list()
      for (sci in 1:length(rv$sizeclasses)){
        selected_mods[[rv$sizeclasses[sci]]] <- paste0("dist: ",
          gsub("constant", "1", input[[paste0("modelChoices_CP", sci)]]))
      if (grepl("exponential", selected_mods[[rv$sizeclasses[sci]]]))
          selected_mods[[sci]] <- paste0(selected_mods[[sci]], "; NULL")
        modChoices[[rv$sizeclasses[sci]]] <- input[[paste0("modelChoices_CP", sci)]]
      }
      for (sci in rv$sizeclasses){
        nm <- ifelse(is.null(modChoices[[sci]]),
          "none selected", modChoices[[sci]])
        if (length(rv$sizeclasses) == 1){
          cat("Selected model: ", nm, "\n", file = file, append = TRUE)
        } else {
          cat("Selected model for", sci, ":", nm, "\n",
            file = file, append = TRUE)
        }
        if (!is.null(modChoices[[sci]])){
          suppressWarnings(
            write.table(rv$mods_CP[[sci]][[selected_mods[[sci]]]]$cell_ls,
              row.names = TRUE, file = file, quote = FALSE, append = TRUE))
        }
        cat("\n", file = file, append = TRUE)
      }
      cat("\nAIC Tables", file = file, append = TRUE)
      aicTable <- aicc(rv$mods_CP)
      for (sci in rv$sizeclasses){
        if (!(sci %in% names(aicTable))){
          cat("\n", sci, ": no model selected\n", file = file, append = TRUE)
        } else {
            cat("\n", sci, "\n", file = file, append = TRUE)
            cat("l_formula s_formula AICc deltaAICc\n",
              file = file, append = TRUE)
            suppressWarnings(write.table(aicTable[[sci]],
              col.names = FALSE, row.names = TRUE, quote = FALSE,
              file = file, append = TRUE))
            cat("\n", file = file, append = TRUE)
        }
      }
    }
  )
}

#' @rdname app_download_functions
#'
downloadSEFig <- function(rv){
  downloadHandler(filename = paste0(rv$filename_SE, "_SEfig.png"),
      content = function(file){
        png(file, height = rv$figH_SE, width = rv$figW_SE,
          pointsize = .pointsize, res = .res)
        tryCatch(
          plot(rv$modSet_SE, specificModel = rv$outSEpk),
          error = function(x){plotNA()}
        )
        dev.off()
      }
  )
}

#' @rdname app_download_functions
#'
downloadSEmod <- function(rv, input){
  downloadHandler(filename = paste0(rv$filename_SE, "_SEmod.txt"),
    content = function(file){
      cat("Data: ", rv$filename_SE, "\n",
        "Observation columns: ", rv$obsCols_SE, "\n",
        file = file, sep = " ")
      cat("Confidence level: ", 100 * rv$CL, "%\n\n", sep = "",
        file = file, append = TRUE)
      selected_mods <- list()
      modChoices <- list()
      for (sci in 1:length(rv$sizeclasses)){
        selected_mods[[rv$sizeclasses[sci]]] <-
          gsub("constant", "1", input[[paste0("modelChoices_SE", sci)]])
          modChoices[[rv$sizeclasses[sci]]] <- input[[paste0("modelChoices_SE", sci)]]
      }
      for (sci in rv$sizeclasses){
        nm <- ifelse(is.null(modChoices[[sci]]), "none selected", modChoices[[sci]])
        if (length(rv$sizeclasses) == 1){
          cat("Selected model: ", nm, "\n", file = file, append = TRUE)
        } else {
          cat("Selected model for", sci, ":", nm, "\n",
            file = file, append = TRUE)
        }
        if (!is.null(modChoices[[sci]])){
          suppressWarnings(
            write.table(rv$mods_SE[[sci]][[selected_mods[[sci]]]]$cell_pk,
              row.names = FALSE, file = file, quote = FALSE, append = TRUE))
        }
        cat("\n", file = file, append = TRUE)
      }
      cat("\nAIC Tables", file = file, append = TRUE)
      aicTable <- aicc(rv$mods_SE)
      for (sci in rv$sizeclasses){
        if (!(sci %in% names(aicTable))){
          cat("\n", sci, ": no model selected\n", file = file, append = TRUE)
        } else {
            cat("\n", sci, "\n", file = file, append = TRUE)
            cat("p_formula k_formula AICc deltaAICc\n",
              file = file, append = TRUE)
            suppressWarnings(write.table(aicTable[[sci]],
              col.names = FALSE, row.names = TRUE, quote = FALSE,
              file = file, append = TRUE))
            cat("\n", file = file, append = TRUE)
        }
      }
    }
  )
}
#' @rdname app_download_functions
#'
downloadgFig <- function(rv, sc){
  downloadHandler(filename = "g_fig.png",
      content = function(file){
        png(file, height = rv$figH_g, width = rv$figW_g, units = "px")
        plot(rv$gGeneric[[sc]], CL = rv$CL)
        dev.off()
      }
  )
}

#' @rdname app_download_functions
#'
downloadMres <- function(rv, input){
  downloadHandler(filename = paste0(rv$filename_CO, "_Mres.txt"),
    content = function(file){
      cat(
        "SE Data: ", rv$filename_SE, "\n",
        "CP Data: ", rv$filename_CP, "\n",
        "SS Data: ", rv$filename_SS, "\n",
        "DWP Data: ", rv$filename_DWP, "\n",
        "CO Data: ", rv$filename_CO, "\n",
        file = file, sep = ""
      )
      cat("Confidence level: ", 100 * rv$CL, "%\n\n", sep = "",
        file = file, append = TRUE)
      sumry <- summary(rv$Msplit, CL = rv$CL)
      if (length(attr(rv$Msplit, "vars")) > 0)
        cat("Estimated mortality by", paste(attr(rv$Msplit, "vars"),
          collapse = " and "), "\n", file = file, append = TRUE)
      else
         cat("Estimated mortality \n", file = file, append = TRUE)
      suppressWarnings(write.table(prettySplitTab(sumry), row.names = FALSE,
          file = file, append = TRUE, quote = FALSE))
      cat("\n\nSearcher Efficiency\n", file = file, append = TRUE)
      selected_mods <- list()
      modChoices <- list()
      for (sci in 1:length(rv$sizeclasses)){
        selected_mods[[rv$sizeclasses[sci]]] <-
          gsub("constant", "1", input[[paste0("modelChoices_SE", sci)]])
          modChoices[[rv$sizeclasses[sci]]] <- input[[paste0("modelChoices_SE", sci)]]
      }
      for (sci in rv$sizeclasses){
        nm <- ifelse(is.null(modChoices[[sci]]),
          "none selected", modChoices[[sci]])
        if (length(rv$sizeclasses) == 1){
          cat("Selected model: ", nm, "\n", file = file, append = TRUE)
        } else {
          cat("Selected model for", sci, ":", nm, "\n",
            file = file, append = TRUE)
        }
        if (!is.null(modChoices[[sci]])){
          suppressWarnings(
            write.table(rv$mods_SE[[sci]][[selected_mods[[sci]]]]$cell_pk,
              row.names = FALSE, file = file, quote = FALSE, append = TRUE))
        }
        cat("\n", file = file, append = TRUE)
      }
      cat("\nAIC Tables", file = file, append = TRUE)
      aicTable <- aicc(rv$mods_SE)
      for (sci in rv$sizeclasses){
        if (!(sci %in% names(aicTable))){
          cat("\n", sci, ": no model selected\n", file = file, append = TRUE)
        } else {
            cat("\n", sci, "\n", file = file, append = TRUE)
            cat("p_formula k_formula AICc deltaAICc\n",
              file = file, append = TRUE)
            suppressWarnings(write.table(aicTable[[sci]],
              col.names = FALSE, row.names = TRUE, quote = FALSE,
                file = file, append = TRUE))
            cat("\n", file = file, append = TRUE)
        }
      }

      cat(c("\nCarcass persistence\nLast present: ", rv$ltp, "\nFirst absent: ",
        rv$fta, "\n"), file = file, sep = " ", append = TRUE)
      selected_mods <- list()
      modChoices <- list()
      for (sci in 1:length(rv$sizeclasses)){
        selected_mods[[rv$sizeclasses[sci]]] <- paste0("dist: ",
          gsub("constant", "1", input[[paste0("modelChoices_CP", sci)]]))
      if (grepl("exponential", selected_mods[[rv$sizeclasses[sci]]]))
          selected_mods[[sci]] <- paste0(selected_mods[[sci]], "; NULL")
        modChoices[[rv$sizeclasses[sci]]] <- input[[paste0("modelChoices_CP", sci)]]
      }
      for (sci in rv$sizeclasses){
        nm <- ifelse(is.null(modChoices[[sci]]), "none selected", modChoices[[sci]])
        if (length(rv$sizeclasses) == 1){
          cat("Selected model: ", nm, "\n", file = file, append = TRUE)
        } else {
          cat("Selected model for", sci, ":", nm, "\n",
            file = file, append = TRUE)
        }
        if (!is.null(modChoices[[sci]])){
          suppressWarnings(
            write.table(rv$mods_CP[[sci]][[selected_mods[[sci]]]]$cell_ls,
              row.names = FALSE, file = file, quote = FALSE, append = TRUE))
        }
        cat("\n", file = file, append = TRUE)
      }
      cat("\nAIC Tables", file = file, append = TRUE)
      aicTable <- aicc(rv$mods_CP)
      for (sci in rv$sizeclasses){
        if (!(sci %in% names(aicTable))){
          cat("\n", sci, ": no model selected\n", file = file, append = TRUE)
        } else {
            cat("\n", sci, "\n", file = file, append = TRUE)
            cat("l_formula s_formula AICc deltaAICc\n",
              file = file, append = TRUE)
            suppressWarnings(write.table(aicTable[[sci]], quote = FALSE,
              col.names = FALSE, row.names = TRUE, file = file, append = TRUE))
            cat("\n", file = file, append = TRUE)
        }
      }
    }
  )
}

#' @rdname app_download_functions
#'
downloadMFig <- function(rv, split = TRUE){
  if (split){
      downloadHandler(filename = "M_fig.png",
          content = function(file){
            png(file, height = rv$figH_M, width = rv$figW_M, units = "px")
            tryCatch(
              plot(rv$Msplit, CL = rv$CL),
              error = function(x){plotNA()}
            )
            dev.off()
          }
      )
  } else {
    downloadHandler(filename = "M_fig.png",
        content = function(file){
          png(file, height = rv$figH_M, width = rv$figW_M, units = "px")
          tryCatch(
            plot(rv$M, CL = rv$CL),
            error = function(x){plotNA()}
          )
          dev.off()
        }
    )
  }
}

#' @rdname app_download_functions
#'
downloadTable <- function(filename, tablename, csvformat){
  if(is.null(tablename)) return()
  if (Sys.info()['sysname'] == "Windows"){
    colnames(tablename) <- gsub("\u0394", "delta", colnames(tablename))
  }
  downloadHandler(filename = filename, content = function(file){
    fcn <- get(paste0("write.csv", csvformat))
    fcn(x = tablename, file = file, row.names = TRUE)
  })
}
 #' @rdname app_download_functions
#'
downloadgres <- function(rv, input){
  downloadHandler(filename = "_gres.txt",
    content = function(file){
      cat(
        "SE Data: ", rv$filename_SE, "\n",
        "CP Data: ", rv$filename_CP, "\n",
        file = file, sep = ""
      )
      cat("I = ", rv$SS[["I"]], ", span = ", rv$SS[["span"]], "\n",
        file = file, append = TRUE)
      cat("Confidence level: ", 100 * rv$CL, "%\n\n", sep = "",
        file = file, append = TRUE)
      sumry <- summary(rv$gGeneric[[rv$sizeclass_g]], CL = rv$CL)
      cat("Estimated generic detection probability for ", rv$sizeclass_g, "\n",
        file = file, append = TRUE)
#     suppressWarnings(write.table(prettySplitTab(sumry), row.names = FALSE,
     suppressWarnings(write.table(sumry, row.names = FALSE,
        file = file, append = TRUE, quote = FALSE))
      cat("\n\nSearcher Efficiency\n", file = file, append = TRUE)
      selected_mods <- list()
      modChoices <- list()
      for (sci in 1:length(rv$sizeclasses)){
        selected_mods[[rv$sizeclasses[sci]]] <-
          gsub("constant", "1", input[[paste0("modelChoices_SE", sci)]])
          modChoices[[rv$sizeclasses[sci]]] <- input[[paste0("modelChoices_SE", sci)]]
      }
      for (sci in rv$sizeclasses){
        nm <- ifelse(is.null(modChoices[[sci]]),
          "none selected", modChoices[[sci]])
        if (length(rv$sizeclasses) == 1){
          cat("Selected model: ", nm, "\n", file = file, append = TRUE)
        } else {
          cat("Selected model for", sci, ":", nm, "\n",
            file = file, append = TRUE)
        }
        if (!is.null(modChoices[[sci]])){
          suppressWarnings(
            write.table(rv$mods_SE[[sci]][[selected_mods[[sci]]]]$cell_pk,
              row.names = FALSE, file = file, quote = FALSE, append = TRUE))
        }
        cat("\n", file = file, append = TRUE)
      }
      cat("\nAIC Tables", file = file, append = TRUE)
      aicTable <- aicc(rv$mods_SE)
      for (sci in rv$sizeclasses){
        if (!(sci %in% names(aicTable))){
          cat("\n", sci, ": no model selected\n", file = file, append = TRUE)
        } else {
            cat("\n", sci, "\n", file = file, append = TRUE)
            cat("p_formula k_formula AICc deltaAICc\n",
              file = file, append = TRUE)
            suppressWarnings(write.table(aicTable[[sci]],
              col.names = FALSE, row.names = TRUE, quote = FALSE,
                file = file, append = TRUE))
            cat("\n", file = file, append = TRUE)
        }
      }

      cat(c("\nCarcass persistence\nLast present: ", rv$ltp, "\nFirst absent: ",
        rv$fta, "\n"), file = file, sep = " ", append = TRUE)
      selected_mods <- list()
      modChoices <- list()
      for (sci in 1:length(rv$sizeclasses)){
        selected_mods[[rv$sizeclasses[sci]]] <- paste0("dist: ",
          gsub("constant", "1", input[[paste0("modelChoices_CP", sci)]]))
      if (grepl("exponential", selected_mods[[rv$sizeclasses[sci]]]))
          selected_mods[[sci]] <- paste0(selected_mods[[sci]], "; NULL")
        modChoices[[rv$sizeclasses[sci]]] <- input[[paste0("modelChoices_CP", sci)]]
      }
      for (sci in rv$sizeclasses){
        nm <- ifelse(is.null(modChoices[[sci]]), "none selected", modChoices[[sci]])
        if (length(rv$sizeclasses) == 1){
          cat("Selected model: ", nm, "\n", file = file, append = TRUE)
        } else {
          cat("Selected model for", sci, ":", nm, "\n",
            file = file, append = TRUE)
        }
        if (!is.null(modChoices[[sci]])){
          suppressWarnings(
            write.table(rv$mods_CP[[sci]][[selected_mods[[sci]]]]$cell_ls,
              row.names = FALSE, file = file, quote = FALSE, append = TRUE))
        }
        cat("\n", file = file, append = TRUE)
      }
      cat("\nAIC Tables", file = file, append = TRUE)
      aicTable <- aicc(rv$mods_CP)
      for (sci in rv$sizeclasses){
        if (!(sci %in% names(aicTable))){
          cat("\n", sci, ": no model selected\n", file = file, append = TRUE)
        } else {
            cat("\n", sci, "\n", file = file, append = TRUE)
            cat("l_formula s_formula AICc deltaAICc\n",
              file = file, append = TRUE)
            suppressWarnings(write.table(aicTable[[sci]], quote = FALSE,
              col.names = FALSE, row.names = TRUE, file = file, append = TRUE))
            cat("\n", file = file, append = TRUE)
        }
      }
    }
  )
}
