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

#' @title Gallery Display
#' @description
#' Displays gallery of `IFC_img` / `IFC_msk` objects
#' @param ... arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'bypass'(=TRUE).\cr
#' If 'param' is provided 'export'(="base64") and the above parameters will be overwritten.\cr
#' If 'offsets' are not provided extra arguments can also be passed with ... to \code{\link{getOffsets}}.\cr
#' /!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.
#' @param objects integer vector, IDEAS objects ids numbers to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`. 
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.
#' @param image_type image_type of desired offsets. Either "img" or "msk". Default is "img".
#' @param layout a character vector of [acquired channels + 'composite' images] members to export. Default is missing to export everything.\cr
#' Note that members can be missing to be removed from final display.\cr
#' Note that members not found will be automatically removed and a warning will be thrown.
#' @param name id of the datatable container. Default is DisplayGallery.
#' @param caption whether to display caption name or not. Default is FALSE.
#' @param pageLength integer, number of objects to display per page. Default is 10.
#' @param pdf_pageSize string, page dimension when exporting to pdf. Default is "A2".
#' @param pdf_pageOrientation string, page orientation when exporting to pdf. Default is "landscape". Allowed are "landscape" or "portrait".
#' @param pdf_image_dpi integer, desired image resolution. Default is 96, for full resolution.
#' @param extract_max maximum number of objects to extract. Default is 10. Use +Inf to extract all.
#' @param sampling whether to sample objects or not. Default is FALSE.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param mode (\code{\link{objectParam}} argument) color mode export. Either "rgb" or "gray". Default is "rgb".
#' @details arguments of \code{\link{objectExtract}} will be deduced from \code{\link{DisplayGallery}} input arguments.\cr
#' Please note that PDF export link will be available if 'write_to' will not result in a "bmp".\cr
#' Please note that exporting to "tiff" may depend on browser capabilities.\cr
#' Please note that a warning may be sent if gallery to display contains large amount of data. This is due to use of datatable() from \pkg{DT}.\cr
#' Warning message:\cr
#' In instance$preRenderHook(instance) :\cr
#' It seems your data is too big for client-side DataTables. You may consider server-side processing: http://rstudio.github.io/DT/server.html\cr
#' For these reasons, it may be better to use "png" extension to display images.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   cif <- ExtractFromXIF(fileName = file_cif)
#'   info <- getInfo(fileName = file_cif, from = "analysis")
#'   ## randomly show at most 10 "img" objects from file
#'   DisplayGallery(info = info, image_type = "img", extract_max = 10,
#'                  sampling = TRUE, write_to = "example.png")
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return it invisibly returns a list whose members are:\cr
#' -data, data for DT::datatable(),\cr
#' -args, associated arguments to pass to DT::datatable().
#' @export
DisplayGallery <- function(..., 
                           objects,
                           offsets,
                           image_type = "img", 
                           layout, 
                           name = "DisplayGallery", 
                           caption = FALSE, 
                           pageLength = 10L, 
                           pdf_pageSize = "A2", 
                           pdf_pageOrientation = "landscape", 
                           pdf_image_dpi = 96,
                           extract_max = 10, 
                           sampling = FALSE, 
                           display_progress = TRUE, 
                           mode = c("rgb", "gray")[1]) {
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
  name = as.character(name); assert(name, len = 1, typ = "character")
  assert(image_type, len = 1, alw = c("img", "msk"))
  caption = as.logical(caption); assert(caption, len = 1, alw=c(TRUE,FALSE))
  assert(mode, len = 1, alw = c("rgb", "gray"))
  sampling = as.logical(sampling); assert(sampling, len = 1, alw = c(TRUE,FALSE))
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE,FALSE))
  extract_max = as.numeric(extract_max); extract_max = extract_max[extract_max>=0]
  assert(extract_max, len = 1, typ = "numeric")
  pageLength = na.omit(as.integer(pageLength)); pageLength = pageLength[pageLength>=0]
  assert(pageLength, len = 1, typ = "integer")
  name = as.character(name); assert(name, len = 1, typ = "character")
  cap = NULL; if(caption) cap = name
  assert(pdf_pageSize, len = 1, alw = c("4A0", "2A0", "A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10",
                                        "B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10",
                                        "C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
                                        "RA0", "RA1", "RA2", "RA3", "RA4",
                                        "SRA0", "SRA1", "SRA2", "SRA3", "SRA4",
                                        "EXECUTIVE", "FOLIO", "LEGAL", "LETTER", "TABLOID"))
  assert(pdf_pageOrientation, len = 1, alw = c("landscape", "portrait"))
  pdf_image_dpi = na.omit(as.integer(pdf_image_dpi)); pdf_image_dpi = pdf_image_dpi[pdf_image_dpi>=0]
  assert(pdf_image_dpi, len = 1, typ = "integer")

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
                                    export = "base64",
                                    mode = mode,
                                    size = size, 
                                    force_width = force_width), dots_param))
    } else {
      param = do.call(what = "objectParam",
                      args = c(list(info = input$info,
                                    export = "base64",
                                    mode = mode,
                                    size = size, 
                                    force_width = force_width), dots_param))
    }
  } else {
    param = input$param
    param$export = "base64"
    param$mode = mode
 }
  fileName = param$fileName_image
  title_progress = basename(fileName)
  
  # check objects to extract
  nobj = as.numeric(param$objcount)
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(length(tokeep) == 0) {
      warning("DisplayGallery: No objects to display, check the objects you provided.", immediate. = TRUE, call. = FALSE)
      return(invisible(NULL))
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
  
  # extract objects
  sel = split(objects, ceiling(seq_along(objects)/20))
  L=length(sel)
  tryCatch({
    if(display_progress) {
      pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
      ans = lapply(1:L, FUN=function(i) {
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
      ans = lapply(1:L, FUN=function(i) {
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
  if(L>1) {
    ans = do.call(what="c", args=ans)
  } else {
    ans = ans[[1]]
  }

  # change layout
  if(missing(layout)) layout = channel_id
  layout = as.character(layout)
  mess = assert(layout, alw = channel_id, fun = "return")
  if(length(mess) != 0) warning(paste0(mess, "\n - and has been automatically removed from 'layout'"), call. = FALSE, immediate. = TRUE)
  layout = unlist(lapply(layout, FUN = function(x) {
    which(channel_id %in% x)
  }))
  if(length(layout) == 0) stop("'layout' is of length 0 which is not allowed")
  
  # check object_ids
  if(image_type == "img") {
    ids = sapply(ans, attr, which="object_id")
    if(!all(objects == ids)) {
      warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
      dat = cbind(true_ids = names(ans), ids = ids, 
                  do.call(what = "rbind", args = lapply(ans, FUN=function(i) i[layout])))
      txt_col = 2
    } else {
      dat = cbind(ids = ids, 
                  do.call(what = "rbind", args = lapply(ans, FUN=function(i) i[layout])))
      txt_col = 1
    }
  } else {
    dat = cbind(ids = as.integer(gsub("^.*_(.*)$", "\\1", sapply(ans, attr, which="offset_id"))), 
                do.call(what = "rbind", args = lapply(ans, FUN=function(i) i[layout])))
    txt_col = 1
  }

  # disable pdf export if write_to is not tiff or png
  dt_dom = ifelse("bmp" %in% getFileExt(param$write_to), "tpr", "Btpr")
  
  # TODO, find a way to remove warning
  oop = options("DT.warn.size" = FALSE); on.exit(options(oop), add = TRUE)
  # (object.size(ans) > 1.5e6 && getOption('DT.warn.size', TRUE)) is FALSE but I am still getting warning
  # create datatable
  datatable(escape = FALSE, rownames = FALSE, extensions = "Buttons",
            selection = list(mode = 'none'), style = "bootstrap",
            caption = cap,
            elementId = name,
            data = dat,
            autoHideNavigation = TRUE,
            options = list(pageLength = pageLength,
                           dom = dt_dom,
                           autoWidth = FALSE,
                           columnDefs = list(list(orderable = FALSE, targets = "_all"),
                                             list(className = "dt-center",targets = "_all")),
                           buttons = list(list(
                             extend = 'pdf',
                             title = ifelse(caption, name, ' '),
                             text = 'PDF',
                             pageSize = pdf_pageSize,
                             extension = '.pdf',
                             header = TRUE,
                             footer = FALSE,
                             orientation = pdf_pageOrientation,
                             customize = JS("function (doc) {",
                                            "if (doc) {",
                                            "doc.margin = [0,0,0,12];",
                                            "for (var i = 1; i < doc.content[1].table.body.length; i++) {",
                                            "for (var j = 1; j < doc.content[1].table.body[i].length; j++) {",
                                            "var foo = doc.content[1].table.body[i][j].text;",
                                            "var w = foo.indexOf('width=');",
                                            "var h = foo.indexOf('height=');",
                                            "var s = foo.indexOf('src=');",
                                            "var e = foo.length;",
                                            sprintf("var wid = parseInt(foo.substring(w + 7, h - 2))*%s;", pdf_image_dpi/96*.75),
                                            sprintf("var hei = parseInt(foo.substring(h + 8, s - 2))*%s;", pdf_image_dpi/96*.75),
                                            "foo = foo.substring(s + 5, e - 2);",
                                            "doc.content[1].table.body[i][j] = { image: foo, alignment: 'center', width: wid, height: hei };",
                                            "}",
                                            "}",
                                            sprintf("doc.header = { text: '%s', alignment: 'center', fontSize: 15};", param$fileName),
                                            "doc.footer = function(currentPage, pageCount) { return [ { text: currentPage.toString() + ' - ' + pageCount, alignment: 'center' } ] };",
                                            "}",
                                            "}"),
                             exportOptions = list(
                               columns = (txt_col-1):length(layout),
                               stripHtml = FALSE)
                           )))
  )
}
