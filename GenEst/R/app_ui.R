#' @title Create the GenEst User Interface HTML
#'
#' @description This suite of functions in \code{app_ui.R} create the HTML code
#'  underlying the GenEst user interface (UI). See the "GenEst Graphic User
#'  Interface" vignette for a more complete detailing of the codebase underlying
#'   the GenEst UI. \cr \cr \code{GenEstUI}: whole application. Calls 
#'   \code{dataInputPanel}, \code{analysisPanel}, and \code{helpPanel}.
#'
#' @details Currently there are few differences between the local and deployed
#'   versions of GenEst, and the \code{appType} toggle is only included as an 
#'   argument for functions that can produce different versions of the HTML.
#'   At this point, the only content that is different is the disclaimer text
#'   on the Help panel. 
#'
#' @param appType Toggle control for the app, \code{"base"} for local versions
#'   or \code{"deploy"} for hosted version. Currently only differentiates the
#'   disclaimer text. 
#'
#' @return Each function returns a string of HTML code, either as a
#'   \code{"shiny.tag.list"} object (in the case of \code{GenEstUI}) or a
#'   \code{"shiny.tag"} object (in the case of the other functions). \cr \cr 
#'   \code{GenEstUI}: Full GenEst user interface.
#'
#' @name app_ui
#' @export
GenEstUI <- function(appType = "base"){
  navbarPage(navbar(), collapsible = TRUE, windowTitle = createvtext("Name"),
    selected = "Help", id = "GenEstApp", 
    dataInputPanel(),
    analysisPanel(),
    helpPanel(appType),
    GenEstShinyJS(), GenEstInlineCSS()
  )
}

#' @rdname app_ui
#'
dataInputPanel <- function(){
  tabPanel("Data Input", 
    sidebarLayout(
      dataInputSidebar(), 
      loadedDataPanel())
  )
}

#' @rdname app_ui
#'
dataInputSidebar <- function(){
  sidebarPanel(width = 3, 
    h4(b(u("Select Data Files:")), style = "margin-bottom: 20px"),
    dataInputWidget("SE"), 
    dataInputWidget("CP"), 
    dataInputWidget("SS"), 
    dataInputWidget("DWP"), 
    dataInputWidget("CO"),
    br(),
    actionButton("clear_all", "Clear All", style = cButtonStyle("all")),
    br()
  )
}

#' @rdname app_ui
#'
loadedDataPanel <- function(){
  mainPanel(
    tabsetPanel(id = "LoadedDataViz",
      dataTabPanel("SE"), 
      dataTabPanel("CP"), 
      dataTabPanel("SS"), 
      dataTabPanel("DWP"), 
      dataTabPanel("CO")
    )
  )
}

#' @rdname app_ui
analysisPanel <- function(){
  tabPanel("Analyses", 
    tabsetPanel(
      GeneralInputsPanel(),
      SEPanel(),
      CPPanel(),
      MPanel(),
      gPanel()
    )
  )
}

#' @rdname app_ui
GeneralInputsPanel <- function(){
  tabPanel("General Inputs", br(), br(), 
    GeneralInputSidebar()
  )
}

#' @rdname app_ui
GeneralInputSidebar <- function(){
  sidebarPanel(width = 3,
    modelInputWidget("nsim"),
    modelInputWidget("CL"),
    modelInputWidget("class")
  )
}

#' @rdname app_ui
SEPanel <- function(){
  tabPanel("Searcher Efficiency", br(), br(), 
    SESidebar(), 
    SEMainPanel()
  )
}

#' @rdname app_ui
SESidebar <- function(){
  sidebarPanel(width = 3,
    b(u(big("Model Inputs:"))),
    br(), br(),
    modelInputWidget("obsSE"),
    modelInputWidget("predsSE"),
    modelInputWidget("kFixedInput"),
    modelRunWidget("SE"),
    modelOutputWidget("SE")
  )
}

#' @rdname app_ui
SEMainPanel <- function(){
  mainPanel(
    tabsetPanel(id = "analyses_SE",
      selectedDataPanel("SE"),
      modelOutputPanel("SEModComparison"),
      modelOutputPanel("SEFigures"),
      modelOutputPanel("SEEstimates"),
      modelOutputPanel("SEModSelection")
    )
  )
}

#' @rdname app_ui
CPPanel <- function(){
  tabPanel("Carcass Persistence", br(), br(), 
    CPSidebar(), 
    CPMainPanel()
  )
}

#' @rdname app_ui
CPSidebar <- function(){
  sidebarPanel(width = 3,
    b(u(big("Model Inputs:"))),
    br(), br(),
    modelInputWidget("ltp"),
    modelInputWidget("fta"),
    modelInputWidget("predsCP"),
    modelInputWidget("dist"),
    modelRunWidget("CP"),
    modelOutputWidget("CP")
  )
}

#' @rdname app_ui
CPMainPanel <- function(){
  mainPanel(
    tabsetPanel(id = "analyses_CP",
      selectedDataPanel("CP"),
      modelOutputPanel("CPModComparison"),
      modelOutputPanel("CPFigures"),
      modelOutputPanel("CPEstimates"),
      modelOutputPanel("CPModSelection")
    )
  )
}

#' @rdname app_ui
MPanel <- function(){
  tabPanel("Mortality Estimation", br(), br(), 
    MSidebar(),
    MMainPanel()
  )
}

#' @rdname app_ui
MSidebar <- function(){
  sidebarPanel(width = 3,
    b(u(big("Model Inputs:"))),
    br(), br(),
    modelInputWidget("xID"),
    modelInputWidget("frac"),
    modelInputWidget("DWPCol"),
    modelInputWidget("COdate"),
    modelRunWidget("M"),
    modelOutputWidget("M")
  )
}

#' @rdname app_ui
MMainPanel <- function(){
  mainPanel(
    tabsetPanel(id = "analyses_M",
      modelOutputPanel("MFigures"),
      modelOutputPanel("MSummary")
    )
  )
}

#' @rdname app_ui
gPanel <- function(){
  tabPanel("Detection Probability", br(), br(),
    gSidebar(),
    gMainPanel()
  )
}

#' @rdname app_ui
gSidebar <- function(){
  sidebarPanel(width = 3,
    b(u(big("Model Inputs:"))),
    br(), br(),
    modelInputWidget("gSearchInterval"),
    modelInputWidget("gSearchMax"),
    modelRunWidget("g"),
    modelOutputWidget("g")
  )
}

#' @rdname app_ui
gMainPanel <- function(){
  mainPanel(
    tabsetPanel(id = "analyses_g",
#      selectedDataPanel("g"),
      modelOutputPanel("gFigures"),
      modelOutputPanel("gSummary")
    )
  )
}


#' @rdname app_ui
helpPanel <- function(appType = "base"){
  tabPanel("Help",
    h5(br(), "For help, see: ",
      a("GenEst User Guide",
        href = "https://doi.org/10.3133/tm7C19", target = "_blank"), " and ",
      a("GenEst Statistical Models",
        href = " https://doi.org/10.3133/tm7A2", target = "_blank")),
    br(),
    tabsetPanel(
      gettingStartedPanel(),
      downloadsPanel(),
      aboutPanel(),
      disclaimersPanel(appType)
    )
  )
}

#' @rdname app_ui
gettingStartedPanel <- function(){
  tabPanel("Getting Started", 
    gettingStartedContent()
  )
}

#' @rdname app_ui
downloadsPanel<- function(){
  tabPanel("Example Data",
    mainPanel(
      column(10, offset = 0,
        br(),
        h3("Example data sets"),
        br(),
        dataDownloadWidget("cleared"),
        dataDownloadWidget("RP"),
        dataDownloadWidget("RPbat"),
        dataDownloadWidget("powerTower"),
        dataDownloadWidget("PV"),
        dataDownloadWidget("trough"),
        dataDownloadWidget("mock"),
        br(), br(),
        h5("NOTE: .csv versions of the data can be downloaded from the GenEst
          page at ",
          a("code.usgs.gov",
            href = "https://code.usgs.gov/ecosystems/GenEst/-/releases",
            target = "_blank"
          )
        )
      )
    )
  )
}

#' @rdname app_ui
aboutPanel <- function(){
  tabPanel("About", 
    aboutContent()
  )
}

#' @rdname app_ui
disclaimersPanel <- function(appType = "base"){
  tabPanel("Disclaimers", 
    disclaimersContent(appType)
  )
}
