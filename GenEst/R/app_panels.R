#' @title app panel utility functions
#'
#' @description functions for formatting and displaying app panels
#'
#' @param dataType Toggle control for the model type of the panel. One of
#'   "SE", "CP", "SS", "DWP", or "CO".
#'
#' @param modType Toggle control for the model type of the panel. One of
#'   "SE", "CP", or "g".
#'
#' @param outType Toggle control for the model type of the panel. One of
#'   "SEFigures", "SEEstimates", "SEModComparison", "SEModSelection",
#'   "CPFigures", "CPEstimates", "CPModComparison", "CPModSelection",
#'   "gFigures", or "gSummary".
#'
#' @name app_panels
NULL

#' @name app_panels
#'
dataTabPanel <- function(dataType){

  if (!dataType %in% c("SE", "CP", "SS", "DWP", "CO")){
    stop(paste0("input dataType (", dataType, ") not supported"))
  }

  Label <- switch(dataType, "SE" = "Searcher Efficiency",
                            "CP" = "Carcass Persistence",
                            "SS" = "Search Schedule",
                            "DWP" = "Density Weighted Proportion",
                            "CO" = "Carcass Observations")
  TableName <- switch(dataType, "SE" = "data_SE",
                                "CP" = "data_CP",
                                "SS" = "data_SS",
                                "DWP" = "data_DWP",
                                "CO" = "data_CO")
  tabPanel(Label, br(), dataTableOutput(TableName))
}

#' @name app_panels
#'
selectedDataPanel <- function(modType){

  if (!modType %in% c("SE", "CP", "g")){
    stop(paste0("input modType (", modType, ") not supported"))
  }

  tName <- switch(modType,
                  "SE" = "Selected Data",
                  "CP" = "Selected Data")#,
#                  "g" = "Summary")
  Condition1 <- switch(modType,
                  "SE" = "input.obsSE == null",
                  "CP" = "input.ltp == null | input.fta == null")#,
#                  "g" = NULL)
  Condition2 <- switch(modType,
                  "SE" = "output.filename_SE != null &
                          input.obsSE != null",
                  "CP" = "output.filename_CP != null & input.ltp != null &
                          input.fta != null")#,
#                  "g" = NULL)
  Text1 <- switch(modType,
            "SE" = em("Select observation columns to view data"),
            "CP" = em("Select observation columns to view data"))#,
#            "g" = list(br(), br(), textOutput("Run 'Estimate' to view results")))
  Text2 <- switch(modType,
            "SE" = em(textOutput("filename_SE")),
            "CP" = em(textOutput("filename_CP")))#,
#            "g" = NULL)
  Data <- switch(modType,
            "SE" = dataTableOutput("selected_SE"),
            "CP" = dataTableOutput("selected_CP"))#,
#            "g" = NULL)

  tabPanel(tName, br(),
    conditionalPanel(condition = Condition1, Text1), 
    conditionalPanel(condition = Condition2, list(Text2, br(), Data))
  ) 

}

#' @name app_panels
#'
modelOutputPanel <- function(outType){

  if (!outType %in% c("SEFigures", "SEEstimates", "SEModComparison", 
                      "SEModSelection", "CPFigures", "CPEstimates", 
                      "CPModComparison", "CPModSelection", "MFigures", 
                      "MSummary", "gFigures", "gSummary")){
    stop(paste0("input outType (", outType, ") not supported"))
  }

  tName <- switch(outType,
             "SEModComparison" = "Model Comparison",
             "SEFigures" = "Figures",
             "SEEstimates" = "Estimates", 
             "SEModSelection" = "Model Selection",
             "CPModComparison" = "Model Comparison",
             "CPFigures" = "Figures",
             "CPEstimates" = "Estimates", 
             "CPModSelection" = "Model Selection",
             "MFigures" = "Figures", 
             "MSummary" = "Summary", 
             "gFigures" = "Figures", 
             "gSummary" = "Summary")

  Condition <- switch(outType,
                 "SEFigures" = 
                   c("output.SEModDone == null",
                     "output.SEModDone == 'OK'", 
                     "1 == 1"),
                 "SEEstimates" = 
                   c("output.SEModDone == null", 
                     "output.SEModDone == 'OK'", 
                     "1 == 1"),
                 "SEModComparison" =
                   c("output.SEModDone == null", 
                     "output.SEModDone == 'OK'", 
                     "1 == 1"),
                 "SEModSelection" = 
                   c("output.SEModDone == null",
                     "output.SEModDone == 'OK'", 
                     "1 == 1"), 
                 "CPFigures" = 
                   c("output.CPModDone == null", 
                     "output.CPModDone == 'OK'", 
                     "1 == 1"),
                 "CPEstimates" = 
                   c("output.CPModDone == null", 
                     "output.CPModDone == 'OK'", 
                     "1 == 1"),
                 "CPModComparison" = 
                   c("output.CPModDone == null", 
                     "output.CPModDone == 'OK'", 
                     "1 == 1"),
                 "CPModSelection" = 
                   c("output.CPModDone == null",
                     "output.CPModDone == 'OK'", 
                     "1 == 1"),
                 "MFigures" = 
                   c("input.modelChoices_SE1 == null |
                       input.modelChoices_SE1 == null |
                       output.sizeclasses_SE != output.sizeclasses_CP", 
                     "output.MModDone == null & 
                       input.modelChoices_SE1 != null &
                       input.modelChoices_CP1 != null &
                       output.sizeclasses_SE == output.sizeclasses_CP",
                     "output.MModDone == 'OK'"),
                 "MSummary" = 
                   c("input.modelChoices_SE1 == null |
                       input.modelChoices_CP1 == null |
                       output.sizeclasses_SE != output.sizeclasses_CP", 
                     "output.MModDone == null & 
                       input.modelChoices_SE1 != null &
                       input.modelChoices_CP1 != null &
                       output.sizeclasses_SE == output.sizeclasses_CP",
                     "output.MModDone == 'OK'"),
                 "gFigures" = 
                   c("output.gModDone == null", 
                     "output.gModDone == 'OK'", 
                     "1 == 1"),
                 "gSummary" = 
                   c("output.gModDone == null",
                     "output.gModDone == 'OK'", 
                     "1 == 1")
               )

  Content <- switch(outType,
               "SEFigures" = 
                 list(
                   em("Run model to view figures"),
                   list(
                     textOutput("sizeclass_SE1"), br(), 
                     plotOutput("fig_SE", inline = TRUE), br(), br(),
                     downloadButton("dlSEfig", "Download")
                   ),
                   NULL
                 ),
               "SEEstimates" =  
                 list(
                   em("Run model to view model estimates"), 
                   list(textOutput("sizeclass_SE2"), br(), 
                     textOutput("text_SE_est"), br(),
                     dataTableOutput("modTab_SE"), br(),
                     downloadButton("dlSEest", "Download")
                   ),
                   NULL
                 ),
               "SEModComparison" = 
                 list(
                   em("Run models to view model comparison"), 
                   list(textOutput("sizeclass_SE3"), br(), 
                     dataTableOutput("AICcTab_SE"), br(),
                     downloadButton("dlSEAICc", "Download")
                   ),
                   NULL
                 ),
               "SEModSelection" = 
                 list(
                   em("Run models to select models"), 
                   list(htmlOutput("modelMenu_SE"),  br(),
                     downloadButton("dlSEmod", "Download")
                   ),
                   NULL
                 ),
               "CPFigures" = 
                 list(
                   em("Run model to view figures"),
                   list(textOutput("sizeclass_CP1"), br(), 
                     plotOutput("fig_CP", inline = TRUE), br(), br(),
                     downloadButton("dlCPfig", "Download")
                   ),
                   NULL
                 ), 
               "CPEstimates" = 
                 list(
                   em("Run model to view model estimates"), 
                   list(textOutput("sizeclass_CP2"), br(), 
                     textOutput("text_CP_est"), br(),
                     dataTableOutput("modTab_CP"), br(),
                     downloadButton("dlCPest", "Download")
                   ),
                   NULL
                 ),
               "CPModComparison" = 
                 list(
                   em("Run models to view model comparison"),
                   list(textOutput("sizeclass_CP3"), br(), 
                     dataTableOutput("AICcTab_CP"), br(),
                     downloadButton("dlCPAICc", "Download")
                   ),
                   NULL
                 ), 
               "CPModSelection" = 
                  list(
                   em("Run models to select models"),
                   list(htmlOutput("modelMenu_CP"),  br(),
                     downloadButton("dlCPmod", "Download")),
                   NULL
                 ),
               "MFigures" = 
                 list(
                   em("Fitted SE and CP models must be selected before
                     detection probability can be estimated. Return to the 
                     Model Selection tabs under Searcher Efficiency and 
                     Carcass Persistence."),
                   em("Run estimate to view figure"),
                   list(plotOutput("fig_M", inline = TRUE), br(), br(),
                      downloadButton("dlMfig", "Download"),
                      conditionalPanel(
                        condition = "output.MSplitDone == 'OK' & output.nMSplits > 1",
                        br(), br(),
                        radioButtons("cscale", "Common scale for y-axes?",
                          choices = c("Yes", "No"),
                          selected = "No")
                      )
                   )
                 ),
               "MSummary" = 
                 list(
                   em("Select SE and CP models fit to matching carcass
                         classes to run model"), 
                   em("Run estimate to view summary"),
                   list(br(), dataTableOutput("table_M"), br(),
                     downloadButton("dlMres", "Download")
                   )
                 ),
               "gFigures" = 
                 list(
                   em("Run estimate to view figure"),
                   list(textOutput("sizeclass_g1"), br(), 
                     plotOutput("fig_g", inline = TRUE), br(), br(),
                     downloadButton("dlgfig", "Download")
                   ),
                   NULL
                 ), 
               "gSummary" = 
                 list(
                   em("Run estimate to view summary"),
                   list(textOutput("sizeclass_g2"), br(), 
                     br(), dataTableOutput("table_g"), br(),
                    downloadButton("dlgtab", "Download")
                   ),
                   NULL
                 )
           )


  tabPanel(tName, br(), 
    conditionalPanel(condition = Condition[1], Content[1]),
    conditionalPanel(condition = Condition[2], Content[2]),
    conditionalPanel(condition = Condition[3], Content[3])
  )

}
