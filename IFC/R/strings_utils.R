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

#' @title String Truncation
#' @description Truncs character strings
#' @param x a string
#' @param n desired length
#' @details x will be truncated according to 'n' parameter. If x is longer than n '...' are appended.
#' @keywords internal
trunc_string = function(x, n=22) {
  x = as.character(x)
  L = nchar(x)
  if(L > n) return(paste0(substr(x, 1, n),"..."))
  return(x)
}

#' @title File Extenstion Removal
#' @description Remove file extension from file path
#' @param x a file path
#' @details file extension will be removed
#' @keywords internal
remove_ext <- function(x) {
  x = as.character(x)
  ext = getFileExt(x)
  gsub(paste0("\\.",ext,"$"), "", x)
}

#' @title Special Character Replacement
#' @description
#' Helper to replace special character.
#' @param string string where specials will be replaced if found.
#' @param replacement string replacement. Default is "_".
#' @param specials  Default is '\\\\|\\/|\\:|\\*|\\?|\\"|\\<|\\>|\\|'.
#' @keywords internal
specialr <- function(string = "", replacement = "_", specials = '\\\\|\\/|\\:|\\*|\\?|\\"|\\<|\\>|\\|') {
  assert(replacement, len = 1, typ = "character")
  assert(specials, len = 1, typ = "character")
  if(grepl(pattern = specials, x = replacement, perl = TRUE)) stop("'replacement' can't contain 'specials'")
  return(gsub(pattern = specials, replacement = replacement, x = string, perl =TRUE))
}

#' @title Name Protection
#' @description
#' Helper to protect population/region name.
#' @param name population/region names
#' @keywords internal
protectn <- function(name) {
  assert(name, typ="character")
  foo = gsub("([[:punct:]])", "\\\\\\1", name, perl=TRUE)
  paste0("(",paste0(sapply(foo, FUN = function(i) {
    return(paste0("[", i, "]"))
  }), collapse = "|"),")")
}

#' @title String Decomposition with Operators
#' @description
#' Helper that will split population definition into chunks of names and operators.
#' @param definition population definition to be splitted
#' @param all_names the names of all allowed populations
#' @param operators operators used. Default is c("And", "Or", "Not", "(", ")").
#' @param splitter the splitter that will be used to help splitting
#' @keywords internal
splitn <- function(definition, all_names, operators = c("And", "Or", "Not", "(", ")"), splitter = "[`~splitter~`]") {
  assert(definition, len=1, typ="character")
  if(substr(definition, 1, 1) == "|") definition = substr(definition, 2, nchar(definition)) # got a file where graph order start with "|"
  assert(all_names, typ="character")
  assert(operators, typ="character")
  assert(splitter, len=1, typ="character")
  return(gsub(splitter, "", strsplit(gsub(protectn(c(all_names, operators)), paste0(splitter, "\\1", splitter), definition, perl=TRUE), split=paste0(splitter, "|", splitter), fixed = TRUE)[[1]], fixed=TRUE))
}

#' @title String Decomposition with Placeholders
#' @description
#' Helper aiming to detect placeholder pattern
#' @param write_to string. Default is "\%d/\%s_fromR.\%e"
#' @details 
#' -\%s: shortname (i.e. basename without extension)\cr
#' -\%p: first parent directory\cr
#' -\%d: full path directory\cr
#' -\%e: file extension\cr
#' -\%o: object id\cr
#' -\%c: channel
#' @keywords internal
splitp = function(write_to = "%d/%s_fromR.%e") {
  assert(write_to, len = 1, typ = "character")
  foo = strsplit(write_to, split = "%", fixed = TRUE)[[1]]
  if(length(foo) > 1) {
    pre = foo[1]
    foo = gsub("^(s|p|d|e|o|c)(.*)$", "\\1%\\2", foo[-1])
    foo = unlist(strsplit(foo, split = "%", fixed = TRUE))
    out = lapply(c("d","p","s","e","o","c"), FUN=function(char) {
      which(foo == char)
    })
  } else {
    pre = ""
    out = list(0,0,0,0,0,0)
  }
  names(out) = c("dir", "parent", "short", "ext", "object", "channel")
  names(foo) <- rep("", length(foo))
  names(foo)[out[["dir"]]] <- "dir"
  names(foo)[out[["parent"]]] <- "parent"
  names(foo)[out[["short"]]] <- "short"
  names(foo)[out[["ext"]]] <- "ext"
  names(foo)[out[["object"]]] <- "object"
  names(foo)[out[["channel"]]] <- "channel"
  foo = c(pre, foo)
  out = c(out, decomp = list(foo))
  attr(out, "class") <- "splitp_obj"
  return(out)
}

#' @title File Path Decomposition
#' @description
#' Helper that will split file name into chunks
#' @param file path to file
#' @return a named vector with chunks of 'file'\cr
#' dir: full path directory of 'file'\cr
#' parent: first parent directory of 'file'\cr
#' ext: 'file' extension without leading dot\cr
#' short: 'file' with no extension nor dir\cr
#' input: 'file' path as it was provided.
#' @keywords internal
splitf <- function(file = NULL) {
  b_name = basename(file)
  dir = dirname(file)
  if(dir == "") {
    dir = suppressWarnings(normalizePath(file, mustWork = FALSE, winslash = "/"))
  } else {
    dir = suppressWarnings(normalizePath(dir, mustWork = FALSE, winslash = "/"))
  }
  ext = getFileExt(file)
  short = gsub(paste0("\\.", ext, "$"), "", b_name, ignore.case = TRUE)
  out = c("dir" = dir, "parent" = basename(dir), "ext" = ext, "short" = short, "input" = file)
  class(out) <- "splitf_obj"
  return(out)
}
# splitf <- function(file = NULL) {
#   f = normalizePath(file, mustWork = FALSE, winslash = "/")
#   dir = dirname(f)
#   b_name = basename(gsub(dir, "", f))
#   if(dir == "") {
#     dir = f
#   } else {
#     dir = suppressWarnings(normalizePath(dir, mustWork = FALSE, winslash = "/"))
#   }
#   ext = getFileExt(file)
#   short = gsub(paste0("\\.", ext, "$"), "", b_name, ignore.case = TRUE)
#   out = c("dir" = dir, "parent" = basename(dir), "ext" = ext, "short" = short, "input" = file)
#   class(out) <- "splitf_obj"
#   return(out)
# }

#' @title File Path Placeholders Formatting
#' @description
#' Helper to format splitp_obj using splitf_obj, channel and object information.
#' @param splitp_obj object returned by \code{\link{splitp}}. 
#' @param splitf_obj object returned by \code{\link{splitf}}. It will be used to substitute \%d, \%p, \%s and \%e.
#' @param channel string to be used to substitute \%c
#' @param object string to be used to substitute \%o
#' @keywords internal
formatn <- function(splitp_obj, splitf_obj, channel = "", object = "") {
  if(missing(splitf_obj)) {
    splitf_obj = list(dir = "", parent = "", file = "", ext = "")
    class(splitf_obj) <- c("splitf_obj", oldClass(splitf_obj))
  }
  # internal function all these checks are useless
  # if(missing(splitp_obj)) stop("'splitp_obj' can't be missing")
  #   assert(splitp_obj, cla = "splitp_obj")
  #   assert(splitf_obj, cla = c("splitf_obj"))
  #   channel = as.character(channel); assert(channel, len = 1, typ = "character")
  #   object = as.character(object); assert(object, len = 1, typ = "character")
  # }
  N = names(splitp_obj$decomp)
  splitp_obj$decomp[N == "dir"] <- splitf_obj["dir"]
  splitp_obj$decomp[N == "parent"] <- splitf_obj["parent"]
  splitp_obj$decomp[N == "short"] <- splitf_obj["short"]
  splitp_obj$decomp[N == "ext"] <- splitf_obj["ext"]
  splitp_obj$decomp[N == "object"] <- object
  splitp_obj$decomp[N == "channel"] <- channel
  return(paste0(splitp_obj$decomp, collapse=""))
}

#' @title Numeric to String Formatting
#' @name num_to_string
#' @description
#' Formats numeric to string used for features, images, ... values conversion when exporting to xml.
#' @param x a numeric vector.
#' @param precision number of significant decimal digits to keep when abs(x) < 1. Default is 15.
#' @return a string vector.
#' @keywords internal
num_to_string <- function(x, precision = 16) {
  return(cpp_num_to_string(x, precision))
}
