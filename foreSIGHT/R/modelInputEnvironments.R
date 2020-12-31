# Model Environment Definitions
# Used to store user inputs and default parameters for the models
#--------------------------------------------------------------------------

# Is passing in ... argument to the write function a good idea? maybe not - spelling mistakes in arguments would not be detected
#     Ideal to decide the fields inside the enviornment apriori

foreSIGHT_modelEnv <- new.env(parent = emptyenv())
foreSIGHT_modelEnv$P_modelEnv <- new.env(parent = emptyenv())
foreSIGHT_modelEnv$Temp_modelEnv <- new.env(parent = emptyenv())
foreSIGHT_modelEnv$PET_modelEnv <- new.env(parent = emptyenv())
foreSIGHT_modelEnv$Radn_modelEnv <- new.env(parent = emptyenv())

write_model_env <- function(envir,            # environment to write into
                            modelInfo,        # the values of the fields
                            modelTag = NULL,  
                            datInd = NULL) {
  
  switch(modelInfo$simVar,
         "P" = {subenvir <- foreSIGHT_modelEnv$P_modelEnv},
         "Temp" = {subenvir <- foreSIGHT_modelEnv$Temp_modelEnv},
         "PET" = {subenvir <- foreSIGHT_modelEnv$PET_modelEnv},
         "Radn" = {subenvir <- foreSIGHT_modelEnv$Radn_modelEnv}
         )
  
  assign("modelTag", modelTag, envir = subenvir)
  assign("modelInfo", modelInfo, envir = subenvir)
  assign("datInd", datInd, envir = subenvir)
  invisible()
}

