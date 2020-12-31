
# function to create paths and logfile as per IOmode
#-----------------------------------------------------------------------------

create_paths <- function(outputPath) {

    if(!isTRUE(dir.exists(outputPath))){
      dir.create(outputPath)
    }

    pathScenario <- file.path(outputPath, "Scenarios")
    pathMetadata <- file.path(outputPath,"Metadata")
    pathCSV <- file.path(pathScenario,"CSV")
    pathRData <- file.path(pathScenario,"Rdata")
    
    logfilename <- filename(IOmode = "xx", arrayID = NULL)
    logfile <- paste0(outputPath, "/", logfilename)

    paths <- list(Scenario = pathScenario,
                  Metadata = pathMetadata,
                  CSV = pathCSV,
                  RData = pathRData,
                  Logfile = logfile)
    return(paths)
}


# Collating functions that act on options for scenarioGeneration
#---------------------------------------------------------------------------------
















































