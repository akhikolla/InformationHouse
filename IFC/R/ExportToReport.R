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

#' @title Graphical and Statistic Report Generation
#' @description
#' Generates report from `IFC_data` object.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param selection when provided, indices of desired graphs.\cr
#' In such case onepage parameter is set to FALSE.\cr
#' Note that indices are read from left to right, from top to bottom. 
#' @param write_to pattern used to export file(s).
#' Placeholders, like c("\%d/\%s_fromR.pdf", "\%d/\%s_fromR.csv"), will be substituted:\cr
#' -\%d: with full path directory of 'obj$fileName'\cr
#' -\%p: with first parent directory of 'obj$fileName'\cr
#' -\%e: with extension of 'obj$fileName' (without leading .)\cr
#' -\%s: with shortname from 'obj$fileName' (i.e. basename without extension).\cr
#' Exported file(s) extension(s) will be deduced from this pattern. Note that has to be a .pdf and/or .csv.
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' Note that if TRUE, it will overwrite file. In addition a warning message will be sent.
#' @param onepage whether to generate a pdf with all graphs on one page or not. Default is TRUE.
#' @param color_mode Whether to extract colors from obj in white or black mode. Default is 'white'.
#' @param add_key whether to draw a 'global' key under title or in the first 'panel' or 'both'. Default is 'panel'.\cr
#' Accepted values are either: FALSE, 'panel', 'global', 'both' or c('panel', 'global').\cr
#' Note that it only applies when display is seen as overlaying populations.
#' @param precision when graphs is a 2D scatter with population overlay, this argument controls amount of information displayed. Default is "light".\cr
#' -"light", the default, will only display points of same coordinates that are amoung the other layers.\cr
#' -"full" will display all the layers.
#' @param trunc_labels maximum number of characters to display for labels. Default is 38.
#' @param trans transformation function for density graphs. Default is asinh.
#' @param bin default number of bin used for histogram. Default is missing.
#' @param viewport Either "ideas", "data" or "max" defining limits used for the graph. Default is "ideas".\cr
#' -"ideas" will use same limits as the one defined in ideas.\cr
#' -"data" will use data to define limits.\cr
#' -"max" will use data and regions drawn to define limits.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other parameters to be passed.
#' @details depending on 'write_to', function will create .pdf and/or .csv file(s) report with according to graphs found in 'obj'.\cr
#' - csv file if created will contain "Min.","1st Qu.","Median","Mean","3rd Qu.","Max." for each graph found for x and y (if not histogram) for drawn populations and regions.\cr
#' - pdf file if created will contain graphs and to a certain extent some stats "Min.", "Median", "Mean", "Max." (no more than 7 rows).\cr
#' Note that only graphs will be exported (no images, features values, population stats, ...) in the same layout they were created and without sizing.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   tmp <- tempdir(check = TRUE)
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf, extract_images = FALSE,
#'                         extract_offsets = FALSE, display_progress = FALSE)
#'   L = length(daf$graphs)
#'   if(L > 0) { 
#'     ## randomly export at most 5 graphs from daf
#'     sel = sample(1:L, min(5, L))
#'     ExportToReport(obj = daf, selection = sel,
#'                    write_to = paste0(tmp, "\\test.pdf"), overwrite = TRUE)
#'   }
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return It invisibly returns full path of exported .pdf and/or .csv file(s).
#' @export
ExportToReport = function(obj, selection, write_to, overwrite=FALSE, onepage=TRUE,
                          color_mode=c("white","black")[1], add_key="panel", precision=c("light","full")[1],    # parameters to pass to plotGraph
                          trunc_labels=38, trans=asinh, bin, viewport="ideas",                                  # parameters to pass to plotGraph
                          display_progress=TRUE, ...) {
  dots = list(...)
  # backup last state of graphic device
  dv <- dev.cur()
  tryCatch({
  # old_ask <- devAskNewPage(ask = FALSE)
  # on.exit(devAskNewPage(ask = old_ask), add = TRUE)

  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
  # check mandatory param
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  if(missing(write_to)) stop("'write_to' can't be missing")
  assert(write_to, typ = "character")
  file_extension = getFileExt(write_to)
  assert(file_extension, alw = c("pdf", "csv"))
  if(any(duplicated(file_extension))) stop("'write_to' has to be only one .pdf and / or only one .csv")
  title_progress = basename(obj$fileName)
  splitf_obj = splitf(obj$fileName)
  create_pdf = FALSE
  create_csv = FALSE
  export_to_pdf = NULL
  export_to_csv = NULL
  if(any(file_extension%in%"pdf")) {
    create_pdf = TRUE
    splitp_obj_pdf = splitp(write_to[file_extension=="pdf"])
    if(any(splitp_obj_pdf$channel > 0)) message("'write_to' (pdf part) has %c argument but channel information can't be retrieved with ExportToReport()")
    if(any(splitp_obj_pdf$object > 0)) message("'write_to' (pdf part) has %o argument but channel information can't be retrieved with ExportToReport()")
    export_to_pdf = formatn(splitp_obj_pdf, splitf_obj)
  }
  if(any(file_extension%in%"csv")) {
    create_csv = TRUE
    splitp_obj_csv = splitp(write_to[file_extension=="csv"])
    if(any(splitp_obj_csv$channel > 0)) message("'write_to', (csv part) has %c argument but channel information can't be retrieved with ExportToReport()")
    if(any(splitp_obj_csv$object > 0)) message("'write_to' (csv part) has %o argument but channel information can't be retrieved with ExportToReport()")
    export_to_csv = formatn(splitp_obj_csv, splitf_obj)
  }
  write_to = c(export_to_pdf, export_to_csv)
  
  assert(color_mode, len=1, alw=c("white","black"))
  assert(precision, len=1, alw=c("light","full"))
  if(!all(add_key%in%c("panel","global","both",FALSE))) stop("Accepted values for add_key are either: FALSE, 'panel', 'global', 'both' or c('panel', 'global')")
  assert(onepage, len=1, alw=c(TRUE,FALSE))
  assert(overwrite, len=1, alw=c(TRUE,FALSE))
  trunc_labels=na.omit(as.integer(trunc_labels));trunc_labels=trunc_labels[trunc_labels>=0]
  assert(trunc_labels, len=1, typ="integer")
  if(missing(bin)) {
    bin = NULL
  } else {
    bin=na.omit(as.integer(bin)); assert(bin, len=1, typ="integer")
  }
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
  overwritten = FALSE
  tmp = file.exists(write_to)
  if(any(tmp)) {
    if(!overwrite) stop(paste0("file ",paste0(write_to[tmp]," already exists"), collapse="\n"))
    if(create_pdf) if(file.exists(export_to_pdf)) export_to_pdf = normalizePath(export_to_pdf, winslash = "/")
    if(create_csv) if(file.exists(export_to_csv)) export_to_csv = normalizePath(export_to_csv, winslash = "/")
    overwritten = TRUE
  }
  if(create_pdf) {
    tryCatch({
      if(!dir.exists(dirname(export_to_pdf))) {
        if(!dir.create(dirname(export_to_pdf), showWarnings = FALSE, recursive = TRUE)) {
          stop(paste0(export_to_pdf,"\ncan't be created: check dirname"), call. = FALSE)
        }
      }
      pdf(file=export_to_pdf)
    }, error = function(e) stop(paste0(export_to_pdf,"\ncan't be created: check name / currently opened ?"), call. = FALSE))
    dev.off(dev.cur())
    export_to_pdf = normalizePath(export_to_pdf, winslash = "/")
  }
  if(create_csv) {
    tryCatch({
      if(!dir.exists(dirname(export_to_csv))) {
        if(!dir.create(dirname(export_to_csv), showWarnings = FALSE, recursive = TRUE)) {
          stop(paste0(export_to_csv,"\ncan't be created: check dirname"), call. = FALSE)
        }
      }
      write.table(x=rbind(c("pop","count","perc","x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.",
                                   "y-Min.","y-1st Qu.","y-Median","y-Mean","y-3rd Qu.","y-Max.")),
                         sep=",", row.names = FALSE, col.names = FALSE, file=export_to_csv)
      }, error = function(e) stop(paste0(export_to_csv,"\ncan't be created: check name / currently opened ?"), call. = FALSE))
    export_to_csv = normalizePath(export_to_csv, winslash = "/")
  }
  write_to = c(export_to_pdf,export_to_csv)
  message(paste0(ifelse(length(write_to)==2, "files", "file")," will be exported in :\n",paste0(normalizePath(dirname(write_to), winslash = "/"),collapse="\n")))
  tryCatch({
    # shortcuts
    G = obj$graphs
    P = obj$pops
    R = obj$regions
    if(length(G)==0) stop("there is no graph defined in 'obj'")
    # defines layout
    lay=lapply(G, FUN=function(g) with(g, c(x=xlocation,y=ylocation)))
    lay=as.data.frame(do.call(rbind, lay))
    row.names(lay)=1:length(G)
    lay=by(lay, lay$y, FUN=function(d) {
      d=d[order(d$x), ]
      d$x=seq_along(d$x)
      cbind("N"=as.integer(row.names(d)), d, stringsAsFactors=FALSE)
    })
    lay=lapply(1:length(lay), FUN=function(i) {
      d=lay[[i]]
      d$y=i
      return(d)
    })
    lay=do.call("rbind", c(lay, make.row.names=FALSE))
    lay_mat=ftable(by(lay$N, lay[,c("y","x")], FUN=function(x) x))
    if(missing(selection)) {
      selection = 1:length(G)
    } else {
      if(!all(selection%in%(1:length(G)))) stop("'selection' refers to graph absent from 'obj'")
      onepage = FALSE
    }
    if(onepage) { 
      G=G[order(lay$N)]
    } else {
      G=G[lay$N[selection]]
    }
    gl = length(G)
    if(display_progress) {
      pb_gr = newPB(session = dots$session, min = 0, max = gl, initial = 0, style = 3)
      on.exit(endPB(pb_gr), add = TRUE)
    }
    suppressWarnings({
      graphs = lapply(1:gl, FUN=function(i) {
        if(display_progress) {
          setPB(pb = pb_gr, value = i, title = title_progress, label = paste0("computing ",ifelse(create_pdf,"graphs and ",""),"stats"))
        }
        if(length(bin) == 0) {
          g = try(plotGraph(obj, G[[i]], draw=FALSE, color_mode=color_mode, add_key=add_key,
                        precision=precision, trunc_labels=trunc_labels, trans=trans, viewport=viewport), silent = TRUE)
        } else {
          g = try(plotGraph(obj, G[[i]], draw=FALSE, color_mode=color_mode, add_key=add_key,
                        precision=precision, trunc_labels=trunc_labels, trans=trans, viewport=viewport, bin=bin), silent = TRUE)
        }
        if(inherits(x = g, what = "try-error")) {
          foo = arrangeGrob(grid.text(label = paste0("Error: ", attr(x = g, which = "condition")$message), gp=gpar(col="red"), draw = FALSE),
                            top = textGrob(paste0("\n",G[[i]]$title), gp = gpar(fontsize = 8, font=2, lineheight=0.5)))
          if(G[[i]]$type=="histogram") {
            stats = matrix(NA, nrow = 1, ncol = 8)
            colnames(stats) = c("count","perc",
                                "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.")
          } else {
            stats = matrix(NA, nrow = 1, ncol = 14)
            colnames(stats) = c("count","perc",
                                "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.",
                                "y-Min.","y-1st Qu.","y-Median","y-Mean","y-3rd Qu.","y-Max.")
          }
          rownames(stats) = paste0("Error: ", G[[i]]$title)
        } else {
          foo = g$plot
          stats = g$stats
        }
        if(create_csv) {
          write.table(x=rbind(c(G[[i]]$f1,"x")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          if(G[[i]]$type!="histogram") write.table(x=rbind(c(G[[i]]$f2,"y")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          write.table(x=stats, file=export_to_csv, append=TRUE, sep=",", col.names = FALSE)
        }
        stats = stats[,!grepl("Qu", colnames(stats)),drop=FALSE]
        foo_lay = matrix(c(rep(1,9), c(NA,2,NA)), nrow=4, ncol=3, byrow = TRUE)
        if(onepage) {
          if(nrow(stats)> 7) {
            stats = stats[1:8,]
            stats[8, ] = "..."
          }
        }
        foo$vp <- viewport(x=0.5, y=0.5)
        tab = tableGrob(format(stats, scientific=FALSE, digits=5), theme = ttheme_default(base_size=4, base_family = "serif"))
        tab$vp <- viewport(x=0.5, y=unit(1,"npc") - 0.5*sum(tab$heights))
        foo = arrangeGrob(foo, tab, layout_matrix = foo_lay, respect = TRUE)
        return(foo)
      })
    })
    if(create_pdf) {
      if(display_progress) {
        pb_pdf = newPB(session = dots$session, title = title_progress, label = "writing to pdf (no update but file is being processed)", min = 0, max = gl, initial = 0, style = 3)
        on.exit(endPB(pb_pdf), add = TRUE)
      }
      if(onepage) {
        pdf(file=export_to_pdf, width = 3*max(lay$x)*2.54, height = 3*max(lay$y)*2.54, 
            family = "serif", onefile = TRUE, pagecentre = TRUE, useDingbats = FALSE)
        # on.exit(dev.off(which = dev.cur()), add = TRUE)
        # TODO add a progress bar
        # for(i in 1:gl) {
        #   pos = matrix(NA, ncol = max(lay$x), nrow = max(lay$y))
        #   pos[lay_mat == lay$N[i]] <- i
        #   print(pos)
        #   grid.arrange(grobs = graphs[lay$N[1]], top = title_progress, newpage = TRUE, layout_matrix = lay_mat, as.table = FALSE)
        #   if(i == 1) {
        #     plot(graphs[lay$N[i]])
        #     grid.arrange(grobs = graphs[lay$N[i]], top = title_progress, newpage = TRUE, layout_matrix = pos, as.table = FALSE)
        #   } else {
        #     grid.arrange(grobs = graphs[lay$N[i]], newpage = FALSE, layout_matrix = pos, as.table = FALSE)
        #   }
        #   if(display_progress) {
        #     setPB(pb_pdf, value = i, title = title_progress, label = "writing to pdf")
        #   }
        # }
        grid.arrange(grobs = graphs[lay$N], top = title_progress, newpage = TRUE, layout_matrix = lay_mat, as.table = FALSE)
      } else {
        pdf(file=export_to_pdf, paper = "a4", onefile = TRUE, pagecentre = TRUE, useDingbats = FALSE, family = "serif")
        # on.exit(dev.off(which = dev.cur()), add = TRUE)
        for(i in 1:gl) {
          grid.arrange(graphs[[i]], top = title_progress, newpage = TRUE) #, respect = TRUE)
          if(display_progress) {
            setPB(pb_pdf, value = i, title = title_progress, label = "writing to pdf")
          }
        }
      }
    }
  }, error = function(e) {
    message(paste0(ifelse(length(write_to)==2, "files have", "file has"), " been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n"))
    stop(e$message, call. = FALSE)
  })
  message(paste0("\n######################\n",ifelse(length(write_to)==2, "files have", "file has"), " been successfully ", ifelse(overwritten, "overwritten", "exported"),"\n"))
  return(invisible(write_to))
  },
  finally = {
    while(!all(dv == dev.cur())) {
      dev.off(which = rev(dev.cur())[1])
    }
  })
}
