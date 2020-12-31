#' @title HTML parameters
#'
#' @description utility functions for simple HTML tasks
#'
#' @param ... attributes and children of the element
#'
#' @param text text to wrap in the tag
#'
#' @param buttonType "single" (for clearing a single component) or "all" (for
#'   clearing everything).
#'
#' @name app_ui_utilities
#'
NULL
#' @rdname app_ui_utilities
navbar <- function(){
  div(
    div(
      img(src = "GenEst.png", style = "margin-top: -8px;", alt = "GenEst",
        height = 40
       ), 
      small(createvtext("Short"))
    )
  )
}

#' @rdname app_ui_utilities
#'
style <- function(...){
  tags$style(...)
}

#' @rdname app_ui_utilities
#'
ol <- function(...){
  tags$ol(...)
}

#' @rdname app_ui_utilities
#'
ul <- function(...){
  tags$ul(...)
}

#' @rdname app_ui_utilities
#'
li <- function(...){
  tags$li(...)
}

#' @rdname app_ui_utilities
#'
b <- function(...){
  tags$b(...)
}

#' @rdname app_ui_utilities
#'
u <- function(...){
  tags$u(...)
}

#' @rdname app_ui_utilities
#'
small <- function(...){
  tags$small(...)
}

#' @rdname app_ui_utilities
#'
big <- function(text = NULL){
  HTML(paste0("<big>", text, "</big>"))
}

#' @rdname app_ui_utilities
#'
center <- function(text = NULL){
  HTML(paste0("<center>", text, "</center>"))
}

#' @rdname app_ui_utilities
#'
GenEstInlineCSS <- function(...){
  inlineCSS(
    list(".shiny-input-container" = "margin-bottom: 0px", 
         "#file_SE_progress" = "margin-bottom: 2px", 
         "#file_CP_progress" = "margin-bottom: 2px",
         "#file_SS_progress" = "margin-bottom: 2px",
         "#file_DWP_progress" = "margin-bottom: 2px",
         "#file_CO_progress" = "margin-bottom: 2px",
         "#file_SE_clear" = "margin-bottom: 20px", 
         "#file_CP_clear" = "margin-bottom: 20px",
         "#file_SS_clear" = "margin-bottom: 20px",
         "#file_DWP_clear" = "margin-bottom: 20px",
         "#file_CO_clear" = "margin-bottom: 20px",
         "#nsim" = "margin-bottom: 15px",
         "#CL" = "margin-bottom: 15px",
         "#sizeCol" = "margin-bottom: 15px",
         "#obsSE" = "margin-bottom: 15px",
         "#predsSE" = "margin-bottom: 15px",
         "#ltp" = "margin-bottom: 15px",
         "#fta" = "margin-bottom: 15px",
         "#predsCP" = "margin-bottom: 15px",
         "#frac" = "margin-bottom: 15px",
         "#gSearchInterval" = "margin-bottom: 15px",
         "#gSearchMax" = "margin-bottom: 20px",
         "#run_SE" = "margin-bottom: 10px",
         "#run_CP" = "margin-bottom: 10px",
         "#run_M" = "margin-bottom: 10px",
         "#split_M" = "margin-bottom: 10px",
         "#run_g" = "margin-bottom: 10px",
         "#run_SE_clear" = "margin-bottom: 20px",
         "#run_CP_clear" = "margin-bottom: 20px",
         "#run_M_clear" = "margin-bottom: 20px",
         "#run_g_clear" = "margin-bottom: 20px",
         "#split_CO" = "margin-bottom: 15px"
    ), ...
  )
}

#' @rdname app_ui_utilities
#'
GenEstShinyJS <- function(...){
  useShinyjs(...)
}

#' @rdname app_ui_utilities
#'
cButtonStyle <- function(buttonType = "single"){

  if (!buttonType %in% c("single", "all")){
    stop(paste0("button Type ", buttonType, " not supported."))
  }
  if (buttonType == "single"){
    "padding:4px; font-size:80%; background-color: #FFCD72"
  } else if (buttonType == "all"){
    "padding:6px; font-size:90%; background-color: #FF8484;
     align: center; margin-top: 20px"
  }
}



