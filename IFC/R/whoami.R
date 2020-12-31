################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2020 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2020                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,             #
#                     Jean-Pierre Marolleau, Loïc Garçon,                      #
#                     INSERM, UPD, CHU Amiens                                  #
#                                                                              #
# DISCLAIMER:                                                                  #
# -You are using this package on your own risk!                                #
# -We do not guarantee privacy nor confidentiality.                            #
# -This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or #
# contributors be liable for any direct, indirect, incidental, special,        #
# exemplary, or consequential damages (including, but not limited to,          #
# procurement of substitute goods or services; loss of use, data, or profits;  #
# or business interruption) however caused and on any theory of liability,     #
# whether in contract, strict liability, or tort (including negligence or      #
# otherwise) arising in any way out of the use of this software, even if       #
# advised of the possibility of such damage.                                   #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with IFC. If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################

#' @title Input Identification
#' @description
#' Helper that identifies input arguments thanks to their IFC classes even if they are not or mis named.
#' @param entries arguments from the function \code{\link{whoami}} is called.
#' /!\ \code{\link{whoami}} MUST be called explicitly this way: whoami(entries = as.list(match.call())).
#' @param search a non duplicated named list of classes to search for entries.
#' @param reinit whether to reinitialize arguments to their default values in called environment. Default is TRUE.
#' @details if two argument of the same 'search' class are found an error will be thrown.
#' 'fileName' will be searched every time.\cr
#' -at first, as an argument (named or not) of the class designated in 'search' to be the "fileName",\cr
#' -otherwise, as an argument (named or not) of class `fileName`,\cr
#' -otherwise, as a named argument of name "fileName" that was not found using 'search',\cr
#' -and finally, if still not found as the first not named argument not found in 'search' of type string.
#' @return a list whose members are 'fileName': value of fileName if provided as a named argument 
#' in entries and all classes defined in 'search'
#' @keywords internal
whoami = function(entries = as.list(match.call()),
                  search = list(info = "IFC_info", 
                                param = "IFC_param",
                                offsets = "IFC_offset"),
                  reinit = TRUE) {
  eval_from = parent.frame(2)
  entry1 = entries[[1]]
  if(typeof(entry1) == "closure") {
    from = entry1
  } else {
    from = as.character(entry1)
  }
  args = entries[-1]
  L = length(args)
  
  # check for fileName class
  has_filename = names(search) %in% "fileName"
  if(any(has_filename)) { # reorder classes to place fileName 1st
    classes = c(search[has_filename], search[!has_filename])
  } else { # add "fileName" in classes
    classes = c(list("fileName" = "fileName"), search)
  }
  if(anyDuplicated(names(classes))) stop("'search' should be a named list of non duplicated elements")
  
  
  # empty 
  if(L == 0) {
    new = sapply(classes, simplify = FALSE, FUN = function(x) return(NULL))
    attr(new, "was") <- rep(0, times = length(classes))
    attr(new, "from") <- from
    return(new)
  }
  
  # retrieve arguments values
  val = lapply(1:L, FUN=function(i) {
    switch(typeof(args[[i]]),
           "symbol" = return(get(as.character(args[[i]]), eval_from)),
           "language" = return(eval(args[[i]], eval_from))
    )
    return(args[[i]])
  })
  
  # identify arguments of classes we look for
  found = sapply(1:L, FUN = function(i) {
    sapply(classes, FUN = function(k) {
      if(length(k) == 0) return(FALSE)
      return(inherits(x = val[[i]], what = k, which = FALSE))
    })
  })
  
  # count how many times a classes is present
  times = apply(found, 1, sum)
  
  # error if at least 2 arguments are with same class
  foo = times > 1
  if(any(foo)) {
    N = names(args)
    NN = sapply(1:L, FUN=function(i) {
      return(ifelse(N[i] == "", "<unk>", N[i]))
    })
    dup = which(foo)
    msg = sapply(1:length(dup), FUN = function(i) {
      idx = which(found[dup[i], ]==1)
      nam = NN[idx]
      paste0("too many elements of class `", classes[dup[i]], "`: [", paste0(paste(idx, nam, sep = "="), collapse = ",") ,"]")
    })
    stop(paste0("\nIn ", entry1, ":\n", msg, collapse = "\n"))
  }
  
  # if not error identify who is who and where it was
  iam = lapply(1:length(classes), FUN=function(i) {
    return(which(found[i, ] == 1))
  })
  was = rep(0, times = length(classes))
  was[times == 1] <- unlist(iam)
  new = lapply(iam, FUN = function(i) {
    if(length(i) != 0) return(val[[i]])
  })
  names(new) = names(classes)
  
  if(length(new$fileName) == 0) { # fileName was not found in classes
    # search fileName in named arguments that were not identified in search
    fil = (names(args) %in% "fileName")
    if(any(was > 0)) fil = fil[-was[was > 0]]
    if(any(fil)) {
      new$fileName = val[[which(fil)[1]]]
      was = c(which(fil)[1], was)
    } else {
      no_name = names(args) %in% ""
      # search for a character string in non named argument that were not identified in search
      if(any(no_name)) {
        fil = unlist(lapply((1:L)[no_name], FUN = function(i) {
          if(i %in% was) return(NULL)
          if(typeof(val[[i]]) == "character") return(i)
          return(NULL)
        }))
        if(length(fil) != 0) {
          new$fileName = val[[fil[1]]]
          was = c(fil[1], was)
        }  else {
          was = c(0, was)
        }
      } else {
        was = c(0, was)
      }
    }
  }
  attr(new, "was") <- was 
  attr(new, "from") <- from
  
  # reinit to arguments of searched classes found to their default value
  if(reinit) {
    call_from = parent.frame(1)
    form = formals(fun = from)
    mism = setdiff(na.omit(names(args)[was]), "")
    sapply(mism, FUN = function(x) assign(x = x, value = form[[x]], inherits = FALSE, envir = call_from, immediate = TRUE))
  }
  return(new)
}
