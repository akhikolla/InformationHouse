#' @title GenEst Information
#'
#' @description HTML generators for app information and content
#'
#' @param type "Full" or "Short" or "Name" or "NameDate"
#'
#' @param appType "base" (for local version) or "deploy" (for hosted version)
#'
#' @return Panels and text for displaying general information about GenEst
#'
#' @name app_content
NULL

#' @rdname app_content
#'
createvtext <- function(type = "Full"){
  vnumber <- packageDescription("GenEst", fields = "Version")
  vdate <- packageDescription("GenEst", fields = "Date")
  if (type == "Full"){
    vtext <- paste0("This is version ", vnumber, " (", vdate, ")")
  }
  if (type == "Short"){
    vtext <- paste0("v", vnumber)
  }
  if (type == "Name"){
    vtext <- paste0("GenEst ", "v", vnumber)
  }
  if (type == "NameDate"){
    vtext <- paste0("GenEst ", "v", vnumber, " (", vdate, ")")
  }
  return(vtext)
}

#' @rdname app_content
#'
gettingStartedContent <- function(){
  mainPanel(
    column(10, offset = 0,
      br(), br(),
      p("GenEst is an R software package for estimating bird and bat 
        fatalities at wind and solar power facilities. Mortality estimation 
        requires five data files:"
      ),
      ol(li("searcher efficiency field trial results (SE),"),
         li("carcass persistence field trial results (CP),"),
         li("schedule for periodic carcass surveys (SS),"),
         li("fraction of total carcasses falling in the searched area at each
            unit searched or", em("density-weighted proportion"), 
            "(DWP), and"
         ),
         li("summary data from the carcass surveys (CO), including unit and
            search date of each carcass discovery along with values of other,
            optional covariates or carcass characteristics")
      ),
      p("Data sets can be uploaded under the ", code("Data Input"), " tab. ",
        "Formats are explained in the user guide. Example data sets are available for loading into GenEst under the ",
         code("Example Data"), " tab. "),
      br(), 
      p("Analysis involves several steps:"),
      ul(li("uploading data---click the ", code("Data Input"), "tab,"),
         li("entering ", code("General Input"), " parameters---click the ",  
            code("Analyses"), " tab,"),
         li("fitting searcher efficiency and carcass persistence models---
            click the ", code("Searcher Efficiency"), " and ", 
            code("Carcass Persistence"), " tabs, and "),
         li("estimating total mortality and splitting mortality estimate by
            various subcategories (such as species or sector or season) as 
            desired---click the ", code("Mortality Estimation"), " tab")
      )
    )
  )
}


#' @rdname app_content
#'
aboutContent <- function(){
  mainPanel(
    column(10, offset = 0,
      br(), 
      h3(createvtext("NameDate")),
      GenEstAuthors(),
      GenEstGUIauthors(),
      GenEstLicense(),
      GenEstAcknowledgements(),
      GenEstLogos()
    )
  )
}

#' @rdname app_content
#'
GenEstAuthors <- function(){
  HTML(
    paste0(br(), 
      b("Authors: "),
      "Daniel Dalthorp ", 
      a("(USGS)", href = "https://www.USGS.gov", target = "_blank"),
      ", Juniper Simonis ",
      a("(DAPPER Stats)", href = "https://www.dapperstats.com", 
        target = "_blank"),
      ", Lisa Madsen ",
      a("(OSU)", href = "https://oregonstate.edu", target = "_blank"),
      ", Manuela Huso ",
      a("(USGS)", href = "https://www.USGS.gov", target = "_blank"),
      ", Paul Rabie ",
      a("(WEST)", href = "https://www.west-inc.com", target = "_blank"),
      ", Jeffrey Mintz ",
      a("(USGS)", href = "https://www.USGS.gov", target = "_blank"),
      ", Robert Wolpert ",
      a("(Duke)", href = "http://www2.stat.duke.edu/~rlw", target = "_blank"),
      ", Jared Studyvin ",
      a("(WEST)", href = "https://www.west-inc.com", target = "_blank"),
      ", and Franzi Korner-Nievergelt ",
      a("(oikostat)", href = "http://www.oikostat.ch", target = "_blank"), 
      "."
    )
  )
}

#' @rdname app_content
#'
GenEstGUIauthors <- function(){
  HTML(
    paste0(br(), br(),
      b("Web Design and Graphics User Interface Programming:"), 
      " Juniper Simonis ",
      a("(DAPPER Stats)", href = "http://www.dapperstats.com", target = "_blank"),
      ", Daniel Dalthorp ",
      a("(USGS)", href = "https://www.USGS.gov", target = "_blank"),
      "."
    )
  )
}

#' @rdname app_content
GenEstLicense <- function(){
  HTML(
    paste0(br(), br(),
      "GenEst is provided under universal public domain, license ",
      a("CC0 1.0",
        href = "https://creativecommons.org/publicdomain/zero/1.0/legalcode",
        target = "_blank"
      ), "."
    )
  )
}

#' @rdname app_content
#'
GenEstAcknowledgements <- function(){
  HTML(
    paste0(br(), br(),
      "The development of GenEst was supported by ",
      a("The US Bureau of Land Management",
        href = "https://www.blm.gov", target = "_blank"
      ), 
      ", ",    
      a("The US Geological Survey",
        href = "https://www.usgs.gov", target = "_blank"
      ), 
      ", ",
      a("The National Renewable Energy",
        href = "https://www.nrel.gov", target = "_blank"
      ), 
      ", ",
      a("WEST Inc.",
        href = "https://www.west-inc.com", target = "_blank"
      ), 
      ", ",
      a("Bat Conservation International",
        href = "https://www.batcon.org", target = "_blank"
      ), 
      ", ",
      a("Avangrid Renewables",
        href = "http://www.avangridrenewables.us", target = "_blank"
      ), 
      ", ",
      a("American Wind Wildlife Institute",
        href = "https://www.awwi.org", target = "_blank"
      ), 
      ", and ",
      a("Oregon State University",
        href = "https://oregonstate.edu", target = "_blank"
      ), 
      "."
    )
  )
}

#' @rdname app_content
#'
GenEstLogos <- function(){
  HTML(
    paste0(br(), br(),
      a(img(src = "blm.jpg", height = "60", alt = "BLM logo"),
        href = "https://www.blm.gov", target = "_blank"
      ),
      a(img(src = "usgs.png", height = "60", alt = "USGS logo"),
        href = "https://www.usgs.gov", target = "_blank"
      ),
      a(img(src = "nrel.jpg", height = "60", alt = "NREL logo"),
        href = "https://www.nrel.gov", target = "_blank"
      ),
      a(img(src = "west.png", height = "60", alt = "WEST logo"),
        href = "https://www.west-inc.com", target = "_blank"
      ),
      a(img(src = "bci.jpg", height = "60", alt = "BCI logo"),
        href = "https://www.batcon.org", target = "_blank"
      ),
      a(img(src = "awwi.png", height = "60", alt = "AWWI logo"),
        href = "https://www.awwi.org", target = "_blank"
      ),
      a(img(src = "avangrid.png", height = "60", alt = "Avangrid logo"),
        href = "http://www.avangridrenewables.us", target = "_blank"
      ),
      a(img(src = "dapper.png", height = "60", alt = "DAPPER stats logo"),
        href = "https://www.dapperstats.com", target = "_blank"
      ),
      a(img(src = "oikostat.jpg", height = "60", alt = "oikostat logo"),
        href = "http://www.oikostat.ch", target = "_blank"
      ),
      a(img(src = "osu.jpg", height = "60", alt = "Oregon State logo"),
        href = "https://oregonstate.edu/", target = "_blank"
      ),
      a(img(src = "duke.png", height = "60", alt = "Duke logo"),
        href = "https://www.duke.edu", target = "_blank"
      )
    )
  )
}

#' @rdname app_content
#'
disclaimersContent <- function(appType = "base"){
  if (!appType %in% c("base", "deploy")){
    stop(paste0("input appType (", appType, ") not supported"))
  }
  mainPanel(
    column(10, 
      br(),
      h3("US Geological Survey (USGS)"),
      disclaimerUSGS(), 
      br(), 
      h3("Western EcoSystems Technology, Inc. (WEST)"),
      disclaimerWEST(appType)
    )
  )
}


#' @rdname app_content
#'
#' @description \code{disclaimerUSGS} creates the text for the USGS
#'   disclaimer.
#'
#' @export
#'
disclaimerUSGS <- function(){
  "This is an incremental update to GenEst v1.0.0 which was approved for release
  by the U.S. Geological Survey (USGS) after rigorous review in compliance with
  US Geological Survey publishing and associated with IP-101457. This update is
  preliminary or provisional and is subject to revision. It is being provided to
  meet the need for timely best science. The update has not received final
  approval by the U.S. Geological Survey (USGS). No warranty, expressed or
  implied, is made by the USGS or the U.S. Government as to the functionality of
  the software and related material nor shall the fact of release constitute any
  such warranty. The software is provided on the condition that neither the USGS
  nor the U.S. Government shall be held liable for any damages resulting from
  the authorized or unauthorized use of the software."
#  "This software has been approved for release by the U.S. Geological
#  Survey (USGS). Although the software has been subjected to rigorous
#  review, the USGS reserves the right to update the software as needed
#  pursuant to further analysis and review. No warranty, expressed or
#  implied, is made by the USGS or the U.S. Government as to the
#  functionality of the software and related material nor shall the fact of
#  release constitute any such warranty. Furthermore, the software is
#  released on condition that neither the USGS nor the U.S. Government shall
#  be held liable for any damages resulting from its authorized or
#  unauthorized use."
}

#' @rdname app_content
#'
#' @export
#'
disclaimerWEST <- function(appType){
  if (!appType %in% c("base", "deploy")){
    stop(paste0("input appType (", appType, ") not supported"))
  }
  extraBR <- switch(appType, "base" = NULL, "deploy" = c(br(), br()))
  out <- NULL

  if (appType == "deploy"){
    out <- c(out, "Western EcoSystems Technology, Inc. does not host nor 
                   maintain the Shinyapp.io website. It is advised that users 
                   not upload sensitive data containing personally 
                   identifiable information (SSN, birthdates, medical 
                   information, etc.). Western EcoSystems Technology, Inc. is 
                   not liable for any damages, including but not limited to 
                   general, compensatory, special or punitive damages, 
                   sustained by user arising out of another party or entity
                   using said sensitive data or for the use of any data by 
                   another party or entity which is obtained from viruses, 
                   Trojans or other malware. Shinyapp.io is actively 
                   maintained by the RStudio Company on Amazon Web Services.")
  }
  out <- c(out, "This program is an 'AS IS' without warranty of any kind, 
                 either expressed or implied, including but not limited to, 
                 the implied warranties of merchantability and fitness for a
                 particular purpose. The entire risk as to the quality and 
                 performance of the program is with you. Should the program 
                 prove defective, you assume all cost of all necessary 
                 servicing, repair or correction. If this program is modified 
                 and/or redistributed, Western EcoSystems Technology, Inc. is 
                 not liable for any damages, including any general, special, 
                 incidental or consequential damages arising out of the use or 
                 inability to use this program (including but not limited to 
                 loss of data or data being rendered inaccurate or losses 
                 sustained by you or third parties or a failure of the program
                 to operate with any other programs), even if such holder or 
                 other party has been advised of the possibility of such 
                 damages.")
  out
}
