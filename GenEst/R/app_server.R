#' @title The GenEst server definition function
#'
#' @description This suite of functions defines the server-side program for
#'   the GenEst user interface (UI). See the "GenEst Graphic User Interface"
#'   vignette for a more complete detailing of the codebase underlying
#'   the GenEst UI. \cr \cr \code{GenEstServer}: main server function
#'   expressed within the application. 
#'
#' @details \code{GenEstServer} is used as the main server function, and is
#'   therefore included in the \code{server.R} script of the app. This 
#'   function is not used in a standard R function sense, in that it does
#'   not return a value and is not used on its own to have side effects.
#'   The code of the function has two parts: 
#'   \enumerate{
#'     \item preamble that defines all the necessary variables and options
#'     \item \code{\link[shiny]{observeEvent}} calls, one for each event in
#'       the application. Each call to \code{\link[shiny]{observeEvent}} 
#'       includes the \code{eventExpr} (event expression) as the first 
#'       argument and the \code{handlerExpr} (handler expression) as the 
#'       second argument, which is an evaluated (via \code{\link[base]{eval}})
#'       block of code returned from \code{reaction} for the specific 
#'       event, as well as any other control switch arguments needed (such as
#'       \code{ignoreNULL}). 
#'   }
#'
#' @param input \code{input} list for the GenEst GUI.
#'
#' @param output \code{output} list for the GenEst GUI.
#'
#' @param session Environment for the GenEst GUI.
#'
#' @param rv reactive variable
#'
#' @param eventName Character name of the event. One of "clear_all",
#'   "file_SE", "file_SE_clear", "file_CP", "file_CP_clear", "file_SS",
#'   "file_SS_clear", "file_DWP", "file_DWP_clear", "file_CO",
#'   "file_CO_clear", "class", "obsSE", "predsSE", "run_SE", "run_SE_clear",
#'   "outSEclass", "outSEp", "outSEk", "ltp", "fta", "predsCP", "run_CP",
#'   "run_CP_clear", "outCPclass", "outCPdist", "outCPl", "outCPs",
#'   "run_M", "run_M_clear", "split_M", "split_M_clear", "transpose_split",
#'   "run_g", "run_g_clear", or "outgclass".
#'
#' @name app_server
#' @export
GenEstServer <- function(input, output, session){
  rv <- initialReactiveValues()
  output <- initialOutput(rv, output)
  msgs <- msgList()
  options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))
  options(DT.options = list(pageLength = 25))
  options(DT.fillContainer = FALSE)
  options(DT.autoHideNavigation = FALSE)
  options(stringsAsFactors = FALSE)
  observeEvent(input$clear_all,  eval(reaction("clear_all")))
  observeEvent(input$file_SE, eval(reaction("file_SE")))
  observeEvent(input$file_SE_clear, eval(reaction("file_SE_clear")))
  observeEvent(input$file_CP, eval(reaction("file_CP")))
  observeEvent(input$file_CP_clear, eval(reaction("file_CP_clear")))
  observeEvent(input$file_SS, eval(reaction("file_SS")))
  observeEvent(input$file_SS_clear, eval(reaction("file_SS_clear")))
  observeEvent(input$file_DWP, eval(reaction("file_DWP")))
  observeEvent(input$file_DWP_clear, eval(reaction("file_DWP_clear")))
  observeEvent(input$file_CO, eval(reaction("file_CO")))
  observeEvent(input$file_CO_clear, eval(reaction("file_CO_clear")))

  observeEvent(input$class, eval(reaction("class")), ignoreNULL = FALSE)

  observeEvent(input$obsSE, eval(reaction("obsSE")), ignoreNULL = FALSE)
  observeEvent(input$predsSE, eval(reaction("predsSE")), ignoreNULL = FALSE)
  observeEvent(input$run_SE, eval(reaction("run_SE")))
  observeEvent(input$run_SE_clear, eval(reaction("run_SE_clear")))
  observeEvent(input$outSEclass, eval(reaction("outSEclass")))
  observeEvent(input$outSEp, eval(reaction("outSEp")))
  observeEvent(input$outSEk, eval(reaction("outSEk")))

  observeEvent(input$ltp, eval(reaction("ltp")), ignoreNULL = FALSE)
  observeEvent(input$fta, eval(reaction("fta")), ignoreNULL = FALSE)
  observeEvent(input$predsCP, eval(reaction("predsCP")), ignoreNULL = FALSE)
  observeEvent(input$run_CP, eval(reaction("run_CP")))
  observeEvent(input$run_CP_clear, eval(reaction("run_CP_clear")))
  observeEvent(input$outCPclass, eval(reaction("outCPclass")))
  observeEvent(input$outCPdist, eval(reaction("outCPdist")))
  observeEvent(input$outCPl, eval(reaction("outCPl")))
  observeEvent(input$outCPs, eval(reaction("outCPs")))

  observeEvent(input$run_M, eval(reaction("run_M")))
  observeEvent(input$run_M_clear, eval(reaction("run_M_clear")))
  observeEvent(input$split_M, eval(reaction("split_M")))
  observeEvent(input$split_M_clear, eval(reaction("split_M_clear")))
  observeEvent(input$transpose_split, eval(reaction("transpose_split")))
  observeEvent(input$cscale, eval(reaction("cscale")))

  observeEvent(input$run_g, eval(reaction("run_g")))
  observeEvent(input$run_g_clear, eval(reaction("run_g_clear")))
  observeEvent(input$outgclass, eval(reaction("outgclass")))

  observeEvent(input$load_RP, eval(reaction("load_RP")))
  observeEvent(input$load_RPbat, eval(reaction("load_RPbat")))
  observeEvent(input$load_cleared, eval(reaction("load_cleared")))
  observeEvent(input$load_PV, eval(reaction("load_PV")))
  observeEvent(input$load_trough, eval(reaction("load_trough")))
  observeEvent(input$load_powerTower, eval(reaction("load_powerTower")))
  observeEvent(input$load_mock, eval(reaction("load_mock")))
}
#' @rdname app_server
#'
reaction <- function(eventName){
#   creates a handler expression to be
#   used by \code{\link[shiny]{observeEvent}} within \code{GenEstServer},
#   which includes the call to \code{eventReaction} (the function that
#   manages the reaction once the code is maevaluated), any message generation
#   or handling, and the enclosing curly braces. Calls
#   \code{reactionMessageRun} and \code{reactionMessageDone} to create the
#   event-specific reaction expression message components.
  eventOptions <- c("clear_all", "file_SE", "file_SE_clear", "file_CP",
                    "file_CP_clear", "file_SS", "file_SS_clear", "file_DWP",
                    "file_DWP_clear", "file_CO", "file_CO_clear", "class",
                    "obsSE", "predsSE", "run_SE", "run_SE_clear",
                    "outSEclass", "outSEp", "outSEk", "ltp", "fta", "predsCP",
                    "run_CP", "run_CP_clear", "outCPclass", "outCPdist",
                    "outCPl", "outCPs", "run_M", "run_M_clear", "split_M",
                    "split_M_clear", "transpose_split",
                    "run_g", "run_g_clear", "outgclass",
                    "load_RP", "load_RPbat", "load_cleared", "load_PV",
                    "load_trough", "load_powerTower", "load_mock", "cscale")

  if (missing(eventName) || (eventName %in% eventOptions) == FALSE){
    stop("eventName missing or not in list of available eventNames")
  }

  reactFun <-  'eventReaction'
  reactArgs <- paste0('"', eventName, '", rv, input, output, session')
  reactText <- paste0(reactFun, '(', reactArgs, ')')
  reactMsgRun <- reactionMessageRun(eventName)
  reactMsgDone <- reactionMessageDone(eventName)

  reactextFull <- c("{", reactMsgRun, reactText, reactMsgDone, "}")

  return(parse(text = reactextFull))
}

#' @rdname app_server
#'
reactionMessageRun <- function(eventName){
  clearEvents <- c("clear_all", "file_SE_clear", "file_CP_clear", 
                   "file_SS_clear", "file_DWP_clear", "file_CO_clear",
                   "run_SE_clear", "run_CP_clear", "run_g_clear", 
                   "run_M_clear", "split_M_clear")
  reactMsg <- NULL 
  if (eventName == "run_SE"){
    reactMsg <- 'msgs$ModSE <<- msgModRun(msgs, "SE")'
  }
  if (eventName == "run_CP"){
    reactMsg <- 'msgs$ModCP <<- msgModRun(msgs, "CP")'
  }
  if (eventName == "run_g"){
    reactMsg <- 'msgs$Modg <<- msgModRun(msgs, "g")'
  }
  if (eventName == "run_M"){
    reactMsg <- 'msgs$ModM <<- msgModRun(msgs, "M")'
  }
  if (eventName %in% clearEvents){
    reactMsg <- 'clearNotifications(msgs)'
  }
  reactMsg
}

#' @rdname app_server
#'
reactionMessageDone <- function(eventName){
  reactMsg <- NULL
  if (eventName == "run_SE"){
    reactMsg <- 'msgs$ModSE <<- msgModDone(msgs, rv, "SE")'
  }
  if (eventName == "run_CP"){
    reactMsg <- 'msgs$ModCP <<- msgModDone(msgs, rv, "CP")'
  }
  if (eventName == "run_g"){
    reactMsg <- 'msgs$Modg <<- msgModDone(msgs, rv, "g")'
  }
  if (eventName == "run_M"){
    reactMsg <- 'msgs$ModM <<- msgModDone(msgs, rv, "M")'
  }
  if (eventName == "split_M"){
    reactMsg <- 'msgs$ModM <<- msgModDone(msgs, rv, "split")'
  }

  reactMsg
}

#' @rdname app_server
#'
eventReaction <- function(eventName, rv, input, output, session){
  if (eventName == "class"){
    rv <- update_rv("run_SE_clear", rv, input)
    output <- update_output("run_SE_clear", rv, output, input)
    update_input("run_SE_clear", rv, input, session)
    rv <- update_rv("run_CP_clear", rv, input)
    output <- update_output("run_CP_clear", rv, output, input)
    update_input("run_CP_clear", rv, input, session)
    rv <- update_rv("run_M_clear", rv, input)
    output <- update_output("run_M_clear", rv, output, input)
    update_input("run_M_clear", rv, input, session)
    rv <- update_rv("split_M_clear", rv, input)
    output <- update_output("split_M_clear", rv, output, input)
    update_input("split_M_clear", rv, input, session)
    rv <- update_rv("run_g_clear", rv, input)
    output <- update_output("run_g_clear", rv, output, input)
    update_input("run_g_clear", rv, input, session)
    rv <- update_rv(eventName, rv, input)
    output <- update_output(eventName, rv, output, input)
    update_input(eventName, rv, input, session)
  } else {
    rv <- update_rv(eventName, rv, input)
    output <- update_output(eventName, rv, output, input)
    update_input(eventName, rv, input, session)
  }
}

