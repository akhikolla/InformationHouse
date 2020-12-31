#' Get pedigree index.
#' 
#' \code{get.pedindex} returns indices of individuals in the pedigree.
#' 
#' \code{member.set} contains member IDs of individuals of interest.
#' 
#' @param pedinfo dataframe.
#' @param member.set character vector.
#' @return An integer vector of indices for each individual of interest found in \code{pedinfo}.
#' @examples 
#' # a simple pedigree with sibling marriage
#' pedigree = as.character(rep(1, 5))
#' member = as.character(c(11, 12, 21, 22, 31))
#' sex = as.numeric(c(1, 2, 1, 2, 1))
#' father = as.character(c(NA, NA, 11, 11, 21))
#' mother = as.character(c(NA, NA, 12, 12, 22))
#' pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
#'
#' get.pedindex(pedinfo, c("22", "31"))
#' @export
get.pedindex = function(pedinfo, member.set){
  member.set = as.character(member.set)
  member = NULL
  pedinfo = transform(pedinfo, member = as.character(member))
  indices = rep(0, length(member.set))
  for(i in 1:length(member.set)){
    indices[i] = which(pedinfo$member == member.set[i])
  }
  return(indices)
}
