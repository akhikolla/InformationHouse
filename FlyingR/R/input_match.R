# Identify and rename columns in Data
# @author Brian Masinde
# @name .colnames.match
# @param data
# @return data returns columns as a list data

.colnames.match <- function(data) {
  # data should be a dataframe
  if (is.data.frame(data) == FALSE) {
    stop("Expects data as a dataframe", call. = FALSE)
  }

  colNames <- colnames(data)

  variables <- list(
    id = c(),
    name = c(),
    wingSpan = c(),
    wingArea = c(),
    taxon = c(),
    allMass = c(),
    fatMass = c(),
    muscleMass = c()
  )


  variables$id <- grepl("id", ignore.case = TRUE, colNames)

  variables$name <- grepl("species.name|name|species|species name|species_name|scientific.name",
                        ignore.case = TRUE, colNames)

  variables$wingSpan <- grepl("ws|wing.span|wing_span|wing span|wingspan",
                   ignore.case = TRUE, colNames)

  variables$wingArea <- grepl("wa|wing.area|wing area|wingarea|wing_area",
                   ignore.case = TRUE, colNames)

  variables$taxon <- grepl("order|ordo|ord|taxon", ignore.case = TRUE, colNames)

  variables$allMass <- grepl("body.mass|empty.mass|all-up_mass|allupmass|body_mass|bodymass|allmass",
                ignore.case = TRUE, colNames)

  variables$fatMass <- grepl("fat.mass|fatmass|fat_mass",
                ignore.case = TRUE, colNames)

  variables$muscleMass <- grepl("muscle.mass|musclemass|muscle_mass",
                ignore.case = TRUE, colNames)


  for (i in 1:length(variables)) {
    # only one column for each variable
    if (sum(variables[[i]]) > 1) {
      stop(paste0("multiple columns found matching: ", names(variables)[i], "\n",
                  sep = "  "),
           call. = FALSE)
    }

    #  missing column
    if (sum(variables[[i]]) == 0) {
      if (i != 1 && i != 2) {
        stop(paste0("missing column: ", names(variables)[i], "\n",
                    sep = " "), call. = FALSE)
      }# else if(i == 8) {
      #  warning("No column matching muscle mass \n", call. = FALSE)
      #}
    }
  }

  # missing name and ID auto-gen
  if (sum(variables$name) == 0 & sum(variables$id) == 0) {
    message("Identifier column not found. Auto-gen \n")
    data$ID <- 1:nrow(data)
  }

  # non unique species names in data add a suffix of ID
  if (sum(duplicated(data$name)) != 0) {
    dups <- duplicated(data$name)
    for (i in 1:length(dups)) {
      if (dups[i] == TRUE) {
        data$name[i] <- paste(data$name[i], i, sep = "_")
      }
    }
  }

  for (i in 1:length(variables)) {
    # rename columns to standard
    if (i != 1 & sum(variables[[i]]) != 0) {
      colIndex <- which(variables[[i]] == TRUE)
      colnames(data)[colIndex] <- names(variables)[[i]]
    }
  }

  # check factors in taxon
  if (any(data$taxon != 1 & data$taxon != 2)) {
    stop("taxon column values 1 or 2", call. = FALSE)
  }

  return(data)
}
