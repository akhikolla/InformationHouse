library(shiny)
library(DT)

shinyUI(fluidPage(
  title = "signalHsmm",
  theme = shinythemes::shinytheme("united"),
  
  sidebarLayout(
    sidebarPanel(
      #style = "background-color: #e0e0e0;border-color: #E95420;border-width: .25rem",
      includeMarkdown("readme.md"),
      uiOutput("dynamic_ui")
    ),
    
    mainPanel(
      uiOutput("dynamic_tabset")    
    )
  )))
