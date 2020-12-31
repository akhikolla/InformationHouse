#' Check pedigree information.
#' 
#' \code{check.pedinfo} checks that the pedigree information provided is consistent. 
#' 
#' Member ID must be unique. Parents must precede offsprings. Sex information must match parental status, and are coded 1 and 2 for male and female respectively. An error message will be produced only if inconsistencies are found.
#' 
#' @param pedinfo dataframe.
#' @export
check.pedinfo = function(pedinfo){
  member = father = mother = sex = NULL
  pedinfo = transform(pedinfo, member = as.character(member), father = as.character(father), mother = as.character(mother, sex = as.integer(sex)))
  
  if(length(unique(pedinfo$member)) != nrow(pedinfo)){
    stop("Check for duplicated member ID.")
  }
  
  for(i in 1:nrow(pedinfo)){
    if(!is.na(pedinfo$father[i])){
      fatherID = which(pedinfo$member == pedinfo$father[i])
      if(length(fatherID) == 0){
        stop(paste("Individual ", pedinfo$member[i], "'s father not present in pedigree.", sep = ""))
      }
      
      if(fatherID >= i){
        stop(paste("Individual ", pedinfo$member[i], " appears before his/her father.", sep = ""))
      }

      if(pedinfo$sex[fatherID] != 1){
        stop(paste("Check individual ", pedinfo$member[i], "'s father sex.", sep = ""))
      }
    }
    
    if(!is.na(pedinfo$mother[i])){
      motherID = which(pedinfo$member == pedinfo$mother[i])
      if(length(motherID) == 0){
        stop(paste("Individual ", pedinfo$member[i], "'s mother not present in pedigree.", sep = ""))
      }
      
      if(motherID >= i){
        stop(paste("Individual ", pedinfo$member[i], " appears before his/her mother.", sep = ""))
      }
      
      if(pedinfo$sex[motherID] != 2){
        stop(paste("Check individual ", pedinfo$member[i], "'s mother sex.", sep = ""))
      }
    }
  }
}
