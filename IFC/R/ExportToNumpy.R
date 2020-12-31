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

#' @title Numpy Export
#' @description
#' Exports IFC objects to numpy array [objects,height,width,channels]
#' @param ... arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'bypass'(=TRUE).\cr
#' If 'param' is provided the above parameters will be overwritten.\cr
#' If 'offsets' are not provided extra arguments can also be passed with ... \code{\link{getOffsets}}.\cr
#' /!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.
#' @param objects integer vector, IDEAS objects ids numbers to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`. 
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.
#' @param image_type image_type of desired offsets. Either "img" or "msk". Default is "img".
#' @param size a length 2 integer vector of final dimensions of the image, height 1st and width 2nd. Default is c(64,64).
#' @param force_width whether to use information in 'info' to fill size. Default is FALSE.
#' When set to TRUE, width of 'size' argument will be overwritten.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param python path to python. Default is Sys.getenv("RETICULATE_PYTHON").\cr
#' Note that this numpy should be available in this python to be able to export to numpy array file, otherwise 'export' will be forced to "matrix".
#' @param dtype desired array?s data-type. Default is "double". Allowed are "uint8", "int16", "uint16" or "double". If 'mode' is "raw", this parameter will be forced to "int16".
#' @param mode (\code{\link{objectParam}} argument) color mode export. Either "raw", "gray" . Default is "raw".
#' @param export export format. Either "file", "matrix". Default is "matrix".\cr
#' Note that you will need 'reticulate' package installed to be able to export to numpy array file, otherwise 'export' will be forced to "matrix".
#' @param write_to used when 'export' is "file" to compute respectively filename.
#' Exported type will be deduced from this pattern. Allowed export are '.npy'.\cr
#' Placeholders, if found, will be substituted:\cr
#' -\%d: with full path directory\cr
#' -\%p: with first parent directory\cr
#' -\%e: with extension of (without leading .)\cr
#' -\%s: with shortname (i.e. basename without extension)\cr
#' -\%o: with objects (at most 10, will be collapse with "_", if more than one).\cr
#' -\%c: with channel_id (will be collapse with "_", if more than one, composite in any will be bracketed).
#' A good trick is to use:\cr
#' -"\%d/\%s_Obj[\%o]_Ch[\%c].npy", when 'export' is "file"\cr
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' @details arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExportToNumpy}} input arguments.\cr
#' Please note that size parameter has to be supplied and could not be set to (0,) when 'object' length is not equal to one\cr
#' \code{\link{ExportToNumpy}} requires reticulate package, python and numpy installed to create npy file.\cr
#' If one of these is missing, 'export' will be set to "matrix".
#' @return Depending on 'export':\cr
#' -"matrix", an array whose dimensions are [object, height, width, channel].\cr
#' -"file", it invisibly returns path of .npy exported file. 
#' @export
ExportToNumpy <- function(...,
                          objects,
                          offsets, 
                          image_type = "img", 
                          size = c(64,64),
                          force_width = FALSE,
                          display_progress = TRUE,
                          python = Sys.getenv("RETICULATE_PYTHON"),
                          dtype = c("uint8", "int16", "uint16", "double")[3],
                          mode = c("raw", "gray")[1], 
                          export = c("file", "matrix")[2],
                          write_to, 
                          overwrite = FALSE) {
  dots=list(...)
  
  # check input
  input = whoami(entries = as.list(match.call()))
  if(!any(sapply(input, FUN = function(i) length(i) != 0))) {
    stop("can't determine what to extract with provided parameters.\n try to input at least one of: 'fileName', 'info', 'param' or 'offsets'")
  }
  
  # reattribute needed param
  offsets = input[["offsets"]]
  param = input[["param"]]
  if(length(offsets) == 0) {
    fileName = input[["fileName"]]
  } else {
    fileName = attr(offsets, "fileName_image")
  }
  
  # check mandatory param
  assert(image_type, len = 1, alw = c("img", "msk"))
  assert(mode, len = 1, alw = c("raw", "gray"))
  if(mode == "raw") dtype = "int16"
  assert(dtype, len = 1, alw = c("uint8", "int16", "uint16", "double"))
  python_back = Sys.getenv("RETICULATE_PYTHON")
  on.exit(Sys.setenv("RETICULATE_PYTHON" = python_back), add = TRUE)
  Sys.setenv("RETICULATE_PYTHON" = python)
  assert(export, len = 1, alw = c("file", "matrix"))
  if(export == "file") {
    if(!requireNamespace("reticulate")) {
      warning("ExportToNumpy: Please install 'reticulate' to export to numpy array file. 'export' has been forced to \"matrix\"")
      export = "matrix"
    } else {
      if(reticulate::py_numpy_available(initialize = TRUE)) {
        np <- reticulate::import("numpy", convert = FALSE)
      } else {
        warning("ExportToNumpy: Can't find numpy in your python installation. 'export' has been forced to \"matrix\"")
        export = "matrix"
      }
    }
  }
  size = na.omit(as.integer(size[1:2]))
  assert(size, len=2, typ="integer")
  force_width = as.logical(force_width); assert(force_width, len = 1, alw = c(TRUE,FALSE)) 
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE,FALSE))
  overwrite = as.logical(overwrite); assert(overwrite, len = 1, alw = c(TRUE,FALSE))
  assert(python, len = 1, typ = "character")
  
  # process extra parameters
  if(length(dots[["verbose"]]) == 0) { 
    verbose = FALSE
  } else {
    verbose = dots[["verbose"]]
  }
  if(length(dots[["verbosity"]]) == 0) { 
    verbosity = 1
  } else {
    verbosity = dots[["verbosity"]]
  }
  if(length(dots[["fast"]]) == 0) { 
    fast = TRUE
  } else {
    fast = dots[["fast"]]
  }
  fast = as.logical(fast); assert(fast, len = 1, alw = c(TRUE, FALSE))
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  verbosity = as.integer(verbosity); assert(verbosity, len = 1, alw = c(1, 2))
  
  param_extra = names(dots) %in% c("ifd","param","mode","export","size","force_width","bypass","verbose")
  dots = dots[!param_extra] # remove not allowed param
  param_param = names(dots) %in% c("write_to","base64_id","base64_att","overwrite",
                                   "composite","selection","random_seed",
                                   "removal","add_noise","full_range","force_range")
  dots_param = dots[param_param] # keep param_param for objectParam
  dots = dots[!param_param]
  
  # compute object param
  # 1: prefer using 'param' if found,
  # 2: otherwise use 'info' if found,
  # 3: finally look at fileName
  if(length(param) == 0) {  
    if(length(input$info) == 0) { 
      param = do.call(what = "objectParam",
                      args = c(list(fileName = fileName,
                                    export = "matrix",
                                    mode = mode,
                                    size = size, 
                                    force_width = force_width), dots_param))
    } else {
      param = do.call(what = "objectParam",
                      args = c(list(info = input$info,
                                    export = "matrix",
                                    mode = mode,
                                    size = size, 
                                    force_width = force_width), dots_param))
    }
  } else {
    param = input$param
    param$export = "matrix"
    param$mode = mode
    if(force_width) {
      param$size[1] <- size[1]
      param$size[2] <- param$channelwidth
    } else {
      param$size <- size
    }
  }
  
  fileName = param$fileName_image
  title_progress = basename(fileName)
  file_extension = getFileExt(fileName)
  
  # check objects to extract
  nobj = as.integer(param$objcount)
  N = nchar(sprintf("%1.f",abs(nobj-1)))
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(length(tokeep) == 0) {
      if(export == "file") {
        warning(paste0("ExportToNumpy: No objects to export, check the objects you provided.\n",
                       "Can't create 'write_to' =", write_to, " from file.\n", param$fileName_image),
                immediate. = TRUE, call. = FALSE)
        return(invisible(NULL))
      } else {
        warning("ExportToNumpy: No objects to export, check the objects you provided.\n", immediate. = TRUE, call. = FALSE)
        return(NULL)
      }
    }
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }

  if(length(objects)!=1) if(param$size[2] == 0) stop("'size' width should be provided when 'object' length not equal to one")
  if(length(objects)!=1) if(param$size[1] == 0) stop("'size' height should be provided when 'object' length not equal to one")
  
  # check input offsets if any
  compute_offsets = TRUE
  if(length(offsets) != 0) {
    if(!("IFC_offset" %in% class(offsets))) {
      warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
    } else {
      if(attr(offsets, "checksum") != param$checksum) {
        warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
      } else {
        compute_offsets = FALSE
      }
    }
  }
  if(compute_offsets) {
    offsets = suppressMessages(getOffsets(fileName = param$fileName_image, fast = fast, display_progress = display_progress, verbose = verbose))
  }
  
  # check export/write_to
  overwritten = FALSE
  if(export != "matrix") {
    if(missing(write_to)) {
      if(export == "file") stop("'write_to' can't be missing")
      write_to = "%s_numpy.npy"
    }
    assert(write_to, len = 1, typ = "character")
    type = getFileExt(write_to)
    assert(type, len = 1, alw = "npy")
    splitf_obj = splitf(param$fileName_image)
    splitp_obj = splitp(write_to)
    if(length(objects) > 10) {
      obj_text = paste0(sprintf(paste0("%0",N,".f"), objects[1:10]), collapse="_")
      obj_text = paste0(obj_text, "_...")
    } else {
      obj_text = paste0(sprintf(paste0("%0",N,".f"), objects), collapse="_")
    }
    chan_text = paste0(sprintf(paste0("%0",2,".f"), as.integer(param$chan_to_keep)), collapse="_")
    if(length(param$composite) != 0) {
      comp_text = paste0("_(",gsub("/",",",param$composite),")", collapse="_")
    } else {
      comp_text = ""
    }
    write_to = formatn(splitp_obj, 
                        splitf_obj, 
                        object = obj_text,
                        channel = paste0(chan_text,comp_text))
    if(export == "file") {
      dir_name = dirname(write_to)
      if(!dir.exists(dir_name)) if(!dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dir_name))
      if(file.exists(write_to)) {
        if(!overwrite) stop(paste0("file ", write_to, " already exists"))
        write_to = normalizePath(write_to, winslash = "/")
        overwritten = TRUE
      }
      message(paste0("file will be exported in :\n", normalizePath(dirname(write_to), winslash = "/")))
    }
  }
  
  # extract objects
  sel = split(objects, ceiling(seq_along(objects)/20))
  L = length(sel)
  tryCatch({
    if(display_progress) {
      pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
      ans = lapply(1:L, FUN = function(i) {
        setPB(pb, value = i, title = title_progress, label = "exporting objects to numpy")
        do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = param$fileName_image,
                                                                   offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], image_type = image_type),
                                                                   trunc_bytes = 8, 
                                                                   force_trunc = FALSE, 
                                                                   verbose = verbose, 
                                                                   verbosity = verbosity,
                                                                   bypass = TRUE), 
                                                      param = param,
                                                      verbose = verbose,
                                                      bypass = TRUE),
                                                 dots))
        
      })
    } else {
      ans = lapply(1:L, FUN = function(i) {
        do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = param$fileName_image,
                                                                   offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], image_type = image_type),
                                                                   trunc_bytes = 8, 
                                                                   force_trunc = FALSE, 
                                                                   verbose = verbose, 
                                                                   verbosity = verbosity,
                                                                   bypass = TRUE), 
                                                      param = param,
                                                      verbose = verbose,
                                                      bypass = TRUE),
                                                 dots))
      })
    }
  }, error = function(e) {
    stop(e$message, call. = FALSE)
  }, finally = {
    if(display_progress) endPB(pb)
  })
  channel_id = attr(ans[[1]][[1]], "channel_id")
  if (L > 1) {
    ans = do.call(what = "c", args = ans)
  } else {
    ans = ans[[1]]
  }
  channel_names = names(ans[[1]])
  
  # check object_ids
  if(image_type == "img") { 
    ids = sapply(ans, attr, which = "object_id")
  } else {
    ids = as.integer(gsub("^.*_(.*)$", "\\1", sapply(ans, attr, which = "offset_id")))
  }
  if(!all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  
  # create array [obj,height,width,channel]
  switch(dtype,
         "uint8" = {
           ret = aperm(array(as.integer(unlist(ans) * 255), dim = c(nrow(ans[[1]][[1]]), ncol(ans[[1]][[1]]), length(ans[[1]]), length(objects))), perm = c(4,1,2,3))
         },
         "int16" = {
           if(mode == "raw") {
             ret = aperm(array(unlist(ans), dim = c(nrow(ans[[1]][[1]]), ncol(ans[[1]][[1]]), length(ans[[1]]), length(objects))), perm = c(4,1,2,3))
           } else {
             ret = aperm(array(as.integer(unlist(ans) * 32767), dim = c(nrow(ans[[1]][[1]]), ncol(ans[[1]][[1]]), length(ans[[1]]), length(objects))), perm = c(4,1,2,3))
           }
         },
         "uint16" = {
           ret = aperm(array(as.integer(unlist(ans) * 65535), dim = c(nrow(ans[[1]][[1]]), ncol(ans[[1]][[1]]), length(ans[[1]]), length(objects))), perm = c(4,1,2,3))
         },
         {
           ret = aperm(array(unlist(ans), dim = c(nrow(ans[[1]][[1]]), ncol(ans[[1]][[1]]), length(ans[[1]]), length(objects))), perm = c(4,1,2,3))
         })
  dimnames(ret) = list("object" = ids,
                       "height" = NULL,
                       "width" = NULL,
                       "channel" = channel_id)
  attr(ret, "fileName_image") <- param$fileName_image
  attr(ret, "object_id") <- ids
  attr(ret, "offset_id") <- sapply(ans, attr, which = "offset_id")
  attr(ret, "channel_id") <- channel_id
  attr(ret, "channel_names") <- channel_names
  tryCatch({
    if(export == "file") {
      np$save(file = write_to, arr = np$array(ret, dtype=np[[dtype]], order='C'))
      message(paste0("\n######################\n", normalizePath(write_to, winslash = "/", mustWork = FALSE), "\nhas been successfully ", ifelse(overwritten, "overwritten", "exported"), "\n"))
      attr(write_to, "fileName_image") <- param$fileName_image
      attr(write_to, "object_id") <- ids
      attr(write_to, "offset_id") <- sapply(ans, attr, which = "offset_id")
      attr(write_to, "channel_id") <- channel_id
      attr(write_to, "channel_names") <- channel_names
      return(invisible(write_to))
    }
    return(ret)
  }, error = function(e) {
    stop(ifelse(export == "file", paste0(normalizePath(write_to, winslash = "/", mustWork = FALSE), " has been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n", e$message), e$message), call. = FALSE)
  })
}
