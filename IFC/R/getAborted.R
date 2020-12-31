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

#' @title Aborted Batch Files Retrieval
#' @description
#' Try to retrieve files whose processing failed during batch. This is a very beta version
#' @param aborted path to file containing aborted information.\cr
#' If missing, the default, a dialog box will be displayed to choose this file.
#' Note, that if provided 'default_batch_dir' and 'config_file' will not be used.
#' @param default_batch_dir directory where batches are stored.\cr
#' It can be found in IDEAS(R) software, under Options -> Application Defaults -> Directories -> Default Batch Report Files Directory.
#' If missing, the default, it will be deduced from IDEAS(R) config file, However, if it can't be deduced then current working directory will be used.\cr
#' This argument takes precedence over 'config_file' and filling 'default_batch_dir' prevents the use of 'config_file' argument.
#' @param config_file path to IDEAS(R) config file.\cr
#' It may depends on IDEAS(R) software installation but one may use "C:/Users/\%USER\%/AppData/Roaming/Amnis Corporation/userconfig.xml".
#' @return a list of 3 elements:\cr
#' -not_existing: a list of files paths that caused failure because they were not found during batch,\cr
#' -failed_found: a list of failed files and their unique corresponding paths,\cr
#' -failed_match: a list of failed files and their all paths that could match.
#' @export
getAborted <- function(aborted, default_batch_dir, config_file) {
  is_valid = FALSE
  if(!missing(aborted)) { 
    if(file.exists(aborted)) {
      if(grepl("^.*Aborted.txt$", aborted, ignore.case = TRUE)) {
        tmp_aborted = try(read_xml(aborted), silent = TRUE)
        if(!("try-error" %in% class(tmp_aborted))) {
          if(length(xml_find_first(tmp_aborted, "//AbortBatch")) > 0) is_valid = TRUE
        } 
      }
    }
    if(!is_valid) stop("when provided 'aborted' should be a non-case sensitive valid 'Aborted.txt' file")
  }
  if(!is_valid) { # checks for batch_dir
    batch_dir = NULL
    if(missing(default_batch_dir)) {
      if(!missing(config_file)) {
        config_file = na.omit(as.character(config_file));
        config_file = normalizePath(config_file, winslash = "/", mustWork = FALSE)
        if(length(config_file) == 1) if(file.exists(config_file)) {
          tmp_conf = read_xml(config_file, options = c("HUGE", "RECOVER", "NOENT", "NOBLANKS", "NSCLEAN"))
          batch_dir = xml_text(xml_find_first(tmp_conf, "//BatchDirectory"))
        }
      }
    } else {
      default_batch_dir = na.omit(as.character(default_batch_dir));
      if(length(default_batch_dir)==1) {
        default_batch_dir = normalizePath(default_batch_dir, winslash = "/", mustWork = FALSE)
        if(dir.exists(default_batch_dir)) {
          batch_dir = default_batch_dir
        }
      }
    }
    if(length(batch_dir) == 0) batch_dir = getwd()
    aborted_batch = list.files(path = batch_dir, full.names = TRUE, include.dirs = FALSE, recursive = TRUE, pattern = "Aborted.txt")
    aborted_batch = aborted_batch[order(sapply(aborted_batch, file.mtime))]
    old_wd = getwd()
    on.exit(setwd(old_wd), add= TRUE)
    setwd(batch_dir)
    if(.Platform$OS.type == "windows") {
      aborted = choose.files(caption = paste0("Looking for: Aborted.txt"),
                             default = ifelse(length(aborted_batch) == 0, "", aborted_batch[length(aborted_batch)]),
                             multi = FALSE, 
                             filters = cbind("Text file (*.txt)", "*.txt"))
    } else {
      aborted = file.choose()
    }
    if(length(aborted) == 0) stop("you did not choose any file")
    tmp_aborted = try(read_xml(aborted), silent = TRUE)
    if("try-error" %in% class(tmp_aborted)) stop("choosen file is not valid")
    if(length(xml_find_first(tmp_aborted, "//AbortBatch")) < 1) stop("choosen file is not valid")
  }
  tmp_aborted_attrs = xml_attrs(xml_find_first(tmp_aborted, "//AbortBatch"))
  tmp_aborted_errors = strsplit(tmp_aborted_attrs["error"], split = "\n", fixed = TRUE)[[1]][-1]
  aborted_batch_name = tmp_aborted_attrs["name"]
  pos <- regexpr(":.*$", tmp_aborted_errors)
  
  # retrieve aborted files
  aborted_batch_files = ifelse(pos > -1, substring(tmp_aborted_errors, 1, pos - 1), "")
  
  # retrieve failure message, it can be
  # -"Could not find file", in such case, ideas provides fullname of the missing file
  # it shows the basename of the input file
  # -"Object reference not set to an instance of an object", when file is corrupted
  # it shows the basename of the input file
  # -"File processing was unexpectantly stopped", it seems to happen when creating cif
  # it shows the basename of the future cif file i.e. appended with suffix
  # Others ?
  aborted_batch_reason = ifelse(pos > -1, substring(tmp_aborted_errors, pos + 2), "")
  aborted_batch_reason = gsub("\\'.*", "",  aborted_batch_reason)
  aborted_batch_reason = substring(aborted_batch_reason, 1, nchar(aborted_batch_reason) - 1)

  # retrieve files that were not found during batch
  # aborted_batch_reason == "Could not find file"
  tmp = aborted_batch_reason %in% "Could not find file"
  if(any(tmp)) {
    not_found = lapply(1:sum(tmp), FUN=function(i) {
      normalizePath(gsub("^.*\'(.*)\'.*$", "\\1", tmp_aborted_errors[which(tmp)[i]]), winslash = "/", mustWork = FALSE)
    })
    names(not_found) = aborted_batch_files[tmp]
    aborted_batch_files = aborted_batch_files[!tmp]
    aborted_batch_reason = aborted_batch_reason[!tmp]
  } else {
    not_found = list()
  }
  
  # retrieved failed files
  # aborted_batch_reason == "File processing was unexpectantly stopped" or "Object reference not set to an instance of an object"
  L = length(aborted_batch_files)
  if(L > 0) {
    batch_file = file.path(dirname(aborted), "batch.xml")
    if(file.exists(batch_file)) {
      # read submitted batch file
      tmp_batch = read_xml(file.path(dirname(aborted), "batch.xml"))
      if(xml_attr(xml_find_first(tmp_batch, "//batch"), "folder") != aborted_batch_name) {
        stop(paste0("batch name found does not correspond to expected one: '",aborted_batch_name,"'")) # should not happen
      }
      # get suffix
      batch_suffix = xml_attr(xml_find_first(tmp_batch, "//batch"), "suffix")
      # retrieve input files
      batch_files = normalizePath(sapply(xml_find_all(tmp_batch, "//file"), xml_attr, "name"), winslash = "/", mustWork = FALSE)
      # remove files we already identified thanks to "Could not find files"
      batch_files = setdiff(batch_files, unlist(not_found))
      batch_files_no_ext = gsub("\\.[[:alnum:]]+$", batch_suffix, batch_files)
      base_batch_files_no_ext = basename(batch_files_no_ext)
      
      # retrieve names of files that should have been created
      tocreate_files = paste0(batch_files_no_ext, ".daf") # daf only ?
      
      # retrieve indices of potential matches between input / export
      potential = lapply(1:L, FUN=function(i) {
        foo = gsub("\\.[[:alnum:]]+$", 
                   ifelse(aborted_batch_reason[i] == "File processing was unexpectantly stopped", "", batch_suffix),
                   aborted_batch_files[i]) == base_batch_files_no_ext
        return(which(foo))
      })
      names(potential) <- aborted_batch_files
      
      # identify failed files that match with only one possible final file
      sure = sapply(potential, length) == 1
      found_failed = potential[sure]
      if(!all(sure)) potential = potential[!sure]

      # identify failed file that were not created
      file_not_created = sapply(potential, FUN = function(i) {
        i[which(!file.exists(tocreate_files[i]))]
      })
      
      # secondly identify those that match only with one file
      sure = sapply(file_not_created, length) == 1
      found_failed = c(found_failed, file_not_created[sure])
      if(!all(sure)) potential = potential[!sure]
      
      # identify all files for which we are sure that batch has worked
      idx_files_created = setdiff(1:length(batch_files), unlist(potential))
      
      escape_loop = 100 # to be sure to escape while
      # check all remaining potential
      if(length(potential) != 0) {
        while(length(potential) != 0 && (escape_loop > 0)) {
          sure = c()
          escape_loop = escape_loop - 1
          for(i in 1:length(potential)) {
            idx = min(idx_files_created, unlist(found_failed))
            if(potential[[i]][1] > idx) {
              diff_time = diff(sapply(tocreate_files[c(min(idx, potential[[i]][1]), potential[[i]])], file.mtime))
              foo = setdiff(potential[[i]][which(diff_time < 0 | is.na(diff_time))], unlist(found_failed))
            } else { 
              idx = max(idx_files_created, unlist(found_failed))
              rev_pot = rev(potential[[i]])
              if(rev_pot[1] > idx) next # escape_loop is here to be sure that we won't have potential[[i]][1] > idx and rev_pot[1] > idx forever
              diff_time = diff(sapply(tocreate_files[c(max(idx, rev_pot[1]), rev_pot)], file.mtime))
              foo = rev(setdiff(rev_pot[which(diff_time < 0 | is.na(diff_time))], unlist(found_failed)))
            }
            if(length(foo) != 0) {
              potential[[i]] <- foo[1]
              found_failed = c(found_failed, potential[i])
              sure = c(sure, i)
            }
          }
          if(length(sure) ==0) {
            break
          } else{
            potential = potential[-sure]
          }
        }
      }
      return(list(not_existing = not_found,
                  failed_found = sapply(found_failed, simplify = FALSE, FUN=function(i) batch_files[i]),
                  failed_match = sapply(potential, simplify = FALSE, FUN=function(i) batch_files[i])))
    } else {
      stop("can't find 'batch.xml' file that correspond to 'Aborted.txt'")
    }
  }
  return(list(NULL))
}
