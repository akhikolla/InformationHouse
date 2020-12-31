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

#' @title RIF/CIF Image Values Extraction
#' @name getImagesValues
#' @description
#' Extracts the image values from RIF or CIF as what can be found in DAF files
#' @param fileName path to file.
#' @param offsets Object of class `IFC_offset`. If missing, the default, 'offsets' will be extracted from 'fileName'.\cr
#' This param is not mandatory but it may allow to save time when exporting repeated image value on same file.
#' @param objects integer vector, IDEAS objects ids numbers to extract.\cr
#' If missing, the default, images values from all objects will be extracted.
#' @param display_progress whether to display a progress bar. Default is FALSE.
#' @param fast when no 'offsets' are provided whether to fast extract 'offsets' or not. Default is TRUE.\cr
#' Meaning that 'objects' will be extracting expecting that 'objects' are stored in ascending order.\cr
#' Only apply when 'offsets' are not provided.\cr
#' Note that a warning will be sent if an object is found at an unexpected order.
#' @param ... other arguments to be passed.
#' @return A data.frame is returned.
#' @keywords internal
getImagesValues <- function(fileName, offsets, objects, display_progress = FALSE, fast = TRUE, ...) {
  dots = list(...)
  if(missing(fileName)) stop("'fileName' can't be missing")
  tmp = duplicated(fileName)
  if(any(tmp)) {
    warning(paste0("duplicated files have been removed from 'fileName': ","\n-", paste0(fileName[tmp],collapse="\n-")))
    fileName = fileName[!tmp]
  }
  if(length(fileName) != 1) stop("'fileName' should be of length 1")
  fileName = normalizePath(fileName, winslash = "/", mustWork = TRUE)
  title_progress = basename(fileName)
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  IFD = getIFD(fileName = fileName, offsets = "first", trunc_bytes = 8, force_trunc = FALSE, verbose = FALSE, verbosity = 1, bypass = FALSE, ...)
  bits = IFD[[1]]$tags$`258`$map
  tmp = read_xml(getFullTag(IFD = IFD, which = 1, tag = "33027"), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  in_use = as.logical(as.numeric(strsplit(xml_text(xml_find_first(tmp, "//Imaging//ChannelInUseIndicators_0_11")), " ", useBytes = TRUE, fixed=TRUE)[[1]]))
  rm(tmp)
  nobj = IFD[[1]]$tags$`33018`$map
  chan_number = sum(in_use)
  
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(length(tokeep) == 0) {
      warning("getImagesValues: No objects to extract, check the objects you provided.", immediate. = TRUE, call. = FALSE)
      return(data.frame())
    }
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }
  
  compute_offsets = TRUE
  if(!missing(offsets)) {
    if(!("IFC_offset" %in% class(offsets))) {
      warning("provided offsets do not match with expected ones, offsets will be recomputed", immediate. = TRUE, call. = FALSE)
    } else {
      if(attr(offsets, "checksum") != checksumXIF(fileName)) {
        warning("provided offsets do not match with expected ones, offsets will be recomputed", immediate. = TRUE, call. = FALSE)
      } else {
        compute_offsets = FALSE
      }
    }
  }
  if(compute_offsets) {
    offsets = suppressMessages(getOffsets(fileName = fileName, fast = fast, display_progress = display_progress))
  }
  sel = split(objects, ceiling(seq_along(objects)/20))
  L = length(sel)
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
    ans = lapply(1:L, FUN=function(i) {
      setPB(pb, value = i, title = title_progress, label = "extracting images values (binary)")
      t(sapply(getIFD(fileName = fileName, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], image_type = "img"), trunc_bytes = 12,
                      force_trunc = FALSE, verbose = FALSE, verbosity = 1, bypass = TRUE, ...), FUN = function(ifd) {
                        c(ifd$infos$OBJECT_ID, # id
                          ifd$curr_IFD_offset, # imgIFD
                          ifd$next_IFD_offset, # mskIFD
                          bits,                # spIFD
                          ifd$tags$`256`$map,  # w
                          ifd$tags$`257`$map,  # l
                          ifd$tags$`33012`$map,# fs
                          ifd$tags$`33016`$map,# cl
                          ifd$tags$`33017`$map,# ct
                          ifd$tags$`33071`$map,# objCenterX
                          ifd$tags$`33072`$map,# objCenterY
                          ifd$tags$`33053`$map[1:chan_number],# bgstd
                          ifd$tags$`33052`$map[1:chan_number],# bgmean
                          ifd$tags$`33054`$map[1:chan_number],# satcount
                          ifd$tags$`33055`$map[1:chan_number])# satpercent
                      }))
    })
  } else {
    ans = lapply(1:L, FUN=function(i) {
      t(sapply(getIFD(fileName = fileName, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], image_type = "img"), trunc_bytes = 12,
                      force_trunc = FALSE, verbose = FALSE, verbosity = 1, bypass = TRUE, ...), FUN = function(ifd) {
                        c(ifd$infos$OBJECT_ID, # id
                          ifd$curr_IFD_offset, # imgIFD
                          ifd$next_IFD_offset, # mskIFD
                          bits,                # spIFD
                          ifd$tags$`256`$map,  # w
                          ifd$tags$`257`$map,  # l
                          ifd$tags$`33012`$map,# fs
                          ifd$tags$`33016`$map,# cl
                          ifd$tags$`33017`$map,# ct
                          ifd$tags$`33071`$map,# objCenterX
                          ifd$tags$`33072`$map,# objCenterY
                          ifd$tags$`33053`$map[1:chan_number],# bgstd
                          ifd$tags$`33052`$map[1:chan_number],# bgmean
                          ifd$tags$`33054`$map[1:chan_number],# satcount
                          ifd$tags$`33055`$map[1:chan_number])# satpercent
                      }))
    })
  }
  if(L>1) {
    ans = do.call(what="rbind", args=ans)
  } else {
    ans = ans[[1]]
  }
  images=as.data.frame(ans, stringsAsFactors = FALSE)
  N = c("id","imgIFD","mskIFD","spIFD","w","l","fs","cl","ct","objCenterX","objCenterY",
        paste0("bgstd",(1:chan_number)),
        paste0("bgmean",(1:chan_number)),
        paste0("satcount",(1:chan_number)),
        paste0("satpercent",(1:chan_number)))
  colnames(images) <- N
  rownames(images) <- 1:nrow(images)
  if(!all(objects == images$id)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  class(images) <- c("data.frame", "IFC_images")
  return(images)
}
