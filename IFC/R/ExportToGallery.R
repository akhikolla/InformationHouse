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

#' @title Gallery Export
#' @description
#' Exports gallery of `IFC_img` / `IFC_msk` objects
#' @param ... arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'bypass'(=TRUE).\cr
#' If 'param' is provided 'mode'(="rgb") and the above parameters will be overwritten.\cr
#' If 'offsets' are not provided extra arguments can also be passed with ... \code{\link{getOffsets}}.\cr
#' /!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.
#' @param objects integer vector, IDEAS objects ids numbers to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`. 
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.
#' @param image_type image_type of desired offsets. Either "img" or "msk". Default is "img".
#' @param layout a character vector of [acquired channels + 'composite' images] members to export. Default is missing to export everything.
#' Note that members can be missing to be removed from final gallery export.
#' Note that members not found will be automatically removed and a warning will be thrown.
#' @param export export format. Either "file", "matrix", "base64". Default is "matrix".
#' @param write_to used when 'export' is "file" or "base64" to compute respectively filename or base64 id attribute.
#' Exported type will be deduced from this pattern. Allowed export are '.bmp', '.jpg', '.jpeg', '.png', '.tif', '.tiff'.
#' Note that '.bmp' are faster but not compressed producing bigger data.\cr
#' Placeholders, if found, will be substituted:\cr
#' -\%d: with full path directory\cr
#' -\%p: with first parent directory\cr
#' -\%e: with extension (without leading .)\cr
#' -\%s: with shortname (i.e. basename without extension)\cr
#' -\%o: with objects (at most 10, will be collapse with "_", if more than one).\cr
#' -\%c: with channel_id (will be collapse with "_", if more than one, composite in any will be bracketed).
#' A good trick is to use:\cr
#' -"\%d/\%s_gallery_Obj[\%o]_Ch[\%c].tiff", when 'export' is "file"\cr
#' -"\%s_gallery.bmp", when 'export' is "base64".\cr
#' Note that if missing and 'export' is not "file", 'write_to' will be set to "\%s_gallery.bmp".
#' @param base64_id whether to add id attribute to base64 exported object. Default is TRUE.\cr
#' Only applied when 'export' is "base64".
#' @param base64_att attributes to add to base64 exported object. Default is "".\cr
#' Only applied when export is "base64". For example, use "class=draggable".\cr
#' Note that id (if base64_id is TRUE) and width and height are already used.
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' @param main main title that will be displayed on top center of the image.
#' If too large it will be clipped.
#' @param add_channels whether to add channels names. Default is TRUE.
#' @param add_ids integer, indice of column to mark objects ids number. Default is 1.
#' If add_ids < 1, no ids are added.
#' @param add_lines integer, size of separating lines between objects. Default is 1.
#' If add_lines < 1, no separating lines are added.
#' @param bg_color background color for main, channels and separating lines. Default is "grey20".
#' @param dpi integer, the resolution of the image in DPI (dots per inch). Default is 300.\cr
#' Please note that whetever this parameter is final resolution will be 96 dpi.\cr
#' However image will be scaled according this parameter and magnification factor will be equal to this parameter divided by 96.
#' @param scale a named list whose members are 'size', 'style', 'color', 'xoff', 'yoff'. Default is list() to draw no scale. Otherwise,\cr
#' -'size' positive integer. Scale's bar size in micro-meter. Default is '7'.\cr
#' This parameter can't be lesser than 6px nor higher than image width + scale text.\cr
#' -'style' a character string. Scale's bar style, either "dash" or "line". Default is "dash".\cr
#' -'color' a character string. color of the scale. Default is "white".\cr
#' -'xoff' positive integer. x offset in image to draw scale, starting from bottom left corner.\cr
#' -'yoff' positive integer. y offset in image to draw scale, starting from bottom left corner.
#' @param extract_max maximum number of objects to extract. Default is 10. Use +Inf to extract all.
#' @param sampling whether to sample objects or not. Default is FALSE.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @details arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExportToGallery}} input arguments.
#' TRICK: for exporting only ONE 'objects', set 'add_channels' = FALSE, 'add_ids' >= 1, 'force_width' = FALSE, 'dpi' = 96; this allows generating image with its original size incrusted with its id number.
#' @return Depending on 'export':\cr
#' -"matrix", a rgb array,\cr
#' -"base64", a data-uri string,\cr
#' -"file", an invisible vector of ids corresponding to the objects exported. 
#' @export
ExportToGallery <- function(...,
                            objects,
                            offsets,
                            image_type = "img", 
                            layout, 
                            export = c("file", "matrix", "base64")[2],
                            write_to, 
                            base64_id = FALSE, 
                            base64_att = "", 
                            overwrite = FALSE,
                            main = "", 
                            add_channels = TRUE, 
                            add_ids = 1, 
                            add_lines = 2,
                            bg_color = "grey20", 
                            dpi = 300, 
                            scale = list(),
                            extract_max = 10, 
                            sampling = FALSE, 
                            display_progress = TRUE) {
  dots = list(...)
  # backup last state of device ask newpage and set to FALSE
  old_ask <- devAskNewPage(ask = FALSE)
  on.exit(devAskNewPage(ask = old_ask), add = TRUE)
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
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
  sampling = as.logical(sampling); assert(sampling, len = 1, alw = c(TRUE,FALSE))
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE,FALSE))
  base64_id = as.logical(base64_id); assert(base64_id, len = 1, alw = c(TRUE,FALSE))
  if(!missing(base64_att)) {
    base64_att = na.omit(as.character(base64_att))
    assert(base64_att, len = 1, typ = "character")
  }
  overwrite = as.logical(overwrite); assert(overwrite, len = 1, alw = c(TRUE,FALSE))
  extract_max = na.omit(as.integer(extract_max)); extract_max=extract_max[extract_max>=0]
  assert(extract_max, len = 1, typ = "integer")
  dpi = na.omit(as.integer(dpi)); dpi=dpi[dpi>=0]
  assert(dpi, len = 1, typ = "integer")
  zoom = dpi / 96
  main = as.character(main); assert(main, len = 1, typ = "character")
  add_lines = na.omit(as.integer(add_lines)); add_lines=add_lines[add_lines>=0]
  assert(add_lines, len = 1, typ = "integer")
  add_ids = na.omit(as.integer(add_ids)); assert(add_ids, len = 1, typ = "integer")
  if(any(c(add_channels, add_lines > 0, main != ""))) {
    bg_color = na.omit(as.character(bg_color)); assert(bg_color, len = 1, typ = "character")
    lum = getLuminance(bg_color)
  }

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
  
  # should be checked before being passed to objectParam/objectExtract
  if(length(dots[["size"]]) == 0) {
    size = c(0,0)
  } else {
    size = dots[["size"]]
  }
  # should be checked before being passed to objectParam/objectExtract
  if(length(dots[["force_width"]]) == 0) {
    force_width = TRUE
  } else {
    force_width = dots[["force_width"]]
  }
  
  param_extra = names(dots) %in% c("ifd","param","export","write_to","mode","size","force_width","overwrite","bypass","verbose")
  dots = dots[!param_extra] # remove not allowed param
  param_param = names(dots) %in% c("base64_id","base64_att",
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
                                    mode = "rgb",
                                    size = size, 
                                    force_width = force_width), dots_param))
    } else {
      param = do.call(what = "objectParam",
                      args = c(list(info = input$info,
                                    export = "matrix",
                                    mode = "rgb",
                                    size = size, 
                                    force_width = force_width), dots_param))
    }
  } else {
    param = input$param
    param$export = "matrix"
    param$mode = "rgb"
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
        warning(paste0("ExportToGallery: No objects to export, check the objects you provided.\n",
                       "Can't create 'write_to' =", write_to, " from file.\n", param$fileName_image),
                immediate. = TRUE, call. = FALSE)
        return(invisible(NULL))
      } else {
        warning("ExportToGallery: No objects to export, check the objects you provided.\n", immediate. = TRUE, call. = FALSE)
        return(NULL)
      }
    }
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }
  extract_max = as.integer(min(extract_max, length(objects)))
  if(sampling) {objects=sample(objects,extract_max)} else {objects=objects[1:extract_max]}
  if(length(objects)!=1) if(param$size[2] == 0) stop("'size' width should be provided when 'object' length not equal to one")
  
  # check input offsets if any
  compute_offsets = TRUE
  if(!missing(offsets)) {
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
      write_to = "%s_gallery.bmp"
    }
    assert(write_to, len = 1, typ = "character")
    type = getFileExt(write_to)
    if (type == "jpg") type = "jpeg"
    if (type == "tif") type = "tiff"
    assert(type, len = 1, alw = c("bmp", "jpeg", "png", "tiff"))
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
    tryCatch({
      if(display_progress) {
        pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
        ans = lapply(1:L, FUN = function(i) {
          setPB(pb, value = i, title = title_progress, label = "exporting objects")
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
    
    # change layout
    if(missing(layout)) layout = channel_id
    layout = as.character(layout)
    mess = assert(layout, alw = channel_id, fun = "return")
    if(length(mess) != 0) warning(paste0(mess, "\n - and has been automatically removed from 'layout'"), call. = FALSE, immediate. = TRUE)
    layout = unlist(lapply(layout, FUN = function(x) {
      which(channel_id %in% x)
    }))
    if(length(layout) == 0) stop("'layout' is of length 0 which is not allowed")
    
    # check/add object_ids
    if(add_ids > 0) if(!(add_ids %in% (1:length(layout)))) {
        warning("can't find 'add_ids' in 'layout'")
        add_ids = 0
    }
    col = "white"
    if(image_type == "img") { 
      if(add_ids > 0) if(any(param$brightfield$channel)) if(channel_id[layout[add_ids]] %in% as.character(which(param$brightfield$channel))) col = "black"
      ids = sapply(ans, attr, which = "object_id")
    } else {
      ids = as.integer(gsub("^.*_(.*)$", "\\1", sapply(ans, attr, which = "offset_id")))
    }
    ans = lapply(1:length(objects), FUN = function(i) {
      img = ans[[i]][layout]
      id = ids[i]
      d = dim(img[[1]])
      if(objects[i] != id) warning(paste0("Extracted object_ids:", id, " differs from expected one:", objects[i]), call. = FALSE, immediate. = TRUE)
      if(add_ids > 0) img[[add_ids]] = addText(img[[add_ids]], text = id, col, xoff = 1, yoff = 1)
      array(do.call(rbind, args = lapply(img, FUN = function(i) { apply(i, 3, rbind)})), dim = c(d[1], d[2] * length(img), 3))
    })
    d = dim(ans[[1]])
    ret = ans[[1]]
    BG = col2rgb(bg_color)/255
    sep = NULL
    if(add_lines > 0) {
      sep = array(BG, dim = c(add_lines, dim(ret)[2], 3))
    }
    if(length(ans) > 1) for (i in 2:length(ans)) {
      ret = aperm(array(do.call(c, args = lapply(list(ret, sep, ans[[i]]), FUN = function(k) { aperm(k, perm = c(2,3,1)) })),
                  dim = c(d[2], 3, dim(ret)[1] + add_lines + dim(ans[[i]])[1])), perm = c(3,1,2))
    }
    
    # check/add title + channel names
    chan_names = channel_names[layout]
    B = c(add_channels, main != "")
    if(any(B)) {
      banner = array(BG, dim = c(sum(B) * 20, dim(ret)[2], 3))
      ret = aperm(array(do.call(c, args = lapply(list(banner, ret), FUN = function(k) { aperm(k, perm = c(2,3,1)) })),
                        dim = c(d[2], 3, dim(ret)[1] + sum(B) * 20)), perm = c(3,1,2))
    }
    d = dim(ret)
    yoff = 5
    if(main != "") {
      main_msk = texttomatrix(main)
      while(ncol(main_msk) > d[2]) {
        main = substr(main, 1, nchar(main) - 1)
        main_msk = texttomatrix(paste0(main, "..."))
      }
      main_img = objectColorize(main_msk, ifelse(lum < 128, "white", "black"))
      ret = array(sapply(1:d[3], FUN = function(x) cpp_mark(A = ret[,, x], B = main_img[, , x], mask = main_msk,
                                                            xoff = ceiling((d[2] - ncol(main_msk))/2),
                                                            yoff = yoff, invert = ifelse(lum < 128, FALSE, TRUE))), dim = d)
      yoff = yoff + 20
    }
    if(add_channels) {
      k = d[2]/length(chan_names)
      for(i in 1:length(chan_names)) {
        chan_msk = texttomatrix(chan_names[i])
        while(ncol(chan_msk) > k) {
          chan_names[i] = substr(chan_names[i], 1, nchar(chan_names[i]) -  1)
          chan_msk = texttomatrix(paste0(chan_names[i], "..."))
        }
        chan_img = objectColorize(chan_msk, ifelse(lum < 128, "white", "black"))
        ret = array(sapply(1:d[3], FUN = function(x) cpp_mark(A = ret[, , x], B = chan_img[, , x], mask = chan_msk, 
                                                              xoff = ceiling((i - 1) * k + (k - ncol(chan_msk))/2), 
                                                              yoff = yoff, invert = ifelse(lum < 128, FALSE, TRUE))), dim = d)
      }
      yoff = yoff + 20
    }
    
    # add scale bar
    if(length(scale) != 0) ret = do.call(what = "addScaleBar", args = c(list(image = ret), scale))
    if(export == "file") {
      args = list(filename = write_to, width = zoom*d[2], height = zoom*d[1], units = "px", res = 96)
      if(type == "tiff") args = c(args, list(compression = "none"))
      do.call(what = type, args = args)
      grid.raster(image = ret, interpolate = TRUE)
      dev.off(dev.cur())
      message(paste0("\n######################\n", normalizePath(write_to, winslash = "/", mustWork = FALSE), "\nhas been successfully ", ifelse(overwritten, "overwritten", "exported"), "\n"))
      return(invisible(ids))
    }
    if(export == "base64") {
      if(base64_id) {
        return(sprintf("<img id=%s %s width='%s' height='%s' src='data:image/%s;base64,%s'>",
                       write_to,
                       base64_att,
                       ncol(ret),
                       nrow(ret),
                       type,
                       base64_encode(objectWrite(x = ret, type = type, raw()))))
      } else {
        return(sprintf("<img %s width='%s' height='%s' src='data:image/%s;base64,%s'>",
                       base64_att,
                       ncol(ret),
                       nrow(ret),
                       type,
                       base64_encode(objectWrite(x = ret, type = type, raw()))))
      }
    }
    return(ret)
  }, error = function(e) {
    stop(ifelse(export == "file", paste0(normalizePath(write_to, winslash = "/", mustWork = FALSE), " has been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n", e$message), e$message), call. = FALSE)
  })
}
