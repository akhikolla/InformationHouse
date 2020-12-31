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

#' @title RIF/CIF File Reader
#' @description
#' Extracts data from RIF or CIF Files.
#' @param fileName path to file.
#' @param extract_features whether to extract features from file. Default is TRUE.\cr
#' If TRUE, \code{\link{ExtractFromXIF}} will try to export features. It it fails a message will be sent.\cr
#' Otherwise, graphs, pops and regions will be also extracted.
#' @param extract_images whether to extract images information from file. Default is FALSE.
#' @param extract_offsets whether to extract IFDs offsets from corresponding. Default is FALSE.\cr
#' See \code{\link{getOffsets}} for further details.
#' @param extract_stats whether to extract population statistics. Default is TRUE.
#' @param pnt_in_poly_algorithm algorithm used to determine if object belongs to a polygon region or not. Default is 1.\cr
#' Note that for the moment only 1(Trigonometry) is available.
#' @param pnt_in_poly_epsilon epsilon to determine if object belongs to a polygon region or not. It only applies when algorithm is 1. Default is 1e-12.
#' @param force_default when display information can't be retrieved whether to use default values. Default is TRUE.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @param verbosity quantity of information displayed when verbose is TRUE; 1: normal, 2: rich. Default is 1.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param fast whether to fast 'extract_offsets' or not. Default is TRUE.\cr
#' Meaning that offsets will be extracting expecting that raw object are stored in ascending order.
#' if extract_images is FALSE, a message will be thrown since extraction method does not ensure correct mapping between objects and offsets.\cr
#' if extract_images is TRUE, a warning will be sent if an object is found at an unexpected order.
#' @param recursive whether to recursively apply \code{\link{ExtractFromXIF}} on files defining input fileName when it is a merged. Default is FALSE.
#' @param ... Other arguments to be passed.
#' @source For pnt_in_poly_algorithm, Trigonometry, is an adaptation of Jeremy VanDerWal's code \url{https://github.com/jjvanderwal/SDMTools}
#' @details If extract_stats is TRUE, extract_features will be automatically forced to TRUE.\cr
#' If extract_images is TRUE, extract_offsets will be automatically forced to TRUE.\cr
#' If extract_offsets is TRUE, offsets of images and masks IFDs will be extracted.\cr
#' If extract_images is TRUE, information about images will be extracted.\cr
#' If the input fileName is a merged of several files and recursive is set to TRUE, then ExtractFromXIF will be applied recursively on these files.\cr
#' /!\ Note that features extraction is mandatory to correctly extract graphs, pops, regions and statistics values.\cr
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file, but you can also read rif
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   cif <- ExtractFromXIF(fileName = file_cif)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return A named list of class `IFC_data`, whose members are:\cr
#' -description, a list of descriptive information,\cr
#' -fileName, path of fileName input,\cr
#' -fileName_image, same as fileName,\cr
#' -features, a data.frame of features,\cr
#' -features_def, a describing how features are defined,\cr
#' -graphs, a list of graphical elements found,\cr
#' -pops, a list describing populations found,\cr
#' -regions, a list describing how regions are defined,\cr
#' -images, a data.frame describing information about images,\cr
#' -offsets, an integer vector of images and masks IFDs offsets,\cr
#' -stats, a data.frame describing populations count and percentage to parent and total population,\cr
#' -checksum, current file checksum.\cr
#' If fileName is a merged of several files returned object will be of class `IFC_data` and `Merged`.
#' If recursive is set to "TRUE", ExtractFromXIF will be applied recursively on files defining the merged.
#' and the returned object will be a list of the above-mentionned list for each of these files.
#' @export
ExtractFromXIF <- function(fileName, extract_features = TRUE, extract_images = FALSE, extract_offsets = FALSE, extract_stats = TRUE, 
                           pnt_in_poly_algorithm = 1, pnt_in_poly_epsilon = 1e-12, 
                           force_default = TRUE, verbose = FALSE, verbosity = 1, display_progress = TRUE,
                           fast = TRUE, recursive = FALSE, ...) {
  dots=list(...)
  if(missing(fileName)) stop("'fileName' can't be missing")
  tmp = duplicated(fileName)
  if(any(tmp)) {
    warning(paste0("duplicated files have been removed from 'fileName': ","\n-", paste0(fileName[tmp],collapse="\n-")))
    fileName = fileName[!tmp]
  }
  if(length(fileName) != 1) stop("'fileName' should be of length 1")
  extract_features = as.logical(extract_features); assert(extract_features, len = 1, alw = c(TRUE, FALSE))
  extract_images = as.logical(extract_images); assert(extract_images, len = 1, alw = c(TRUE, FALSE))
  extract_offsets = as.logical(extract_offsets); assert(extract_offsets, len = 1, alw = c(TRUE, FALSE))
  extract_stats = as.logical(extract_stats); assert(extract_stats, len = 1, alw = c(TRUE, FALSE))
  pnt_in_poly_algorithm = as.integer(pnt_in_poly_algorithm); assert(pnt_in_poly_algorithm, len = 1, alw = 1)
  pnt_in_poly_epsilon = as.numeric(pnt_in_poly_epsilon); pnt_in_poly_epsilon = pnt_in_poly_epsilon[pnt_in_poly_epsilon>0]; pnt_in_poly_epsilon = pnt_in_poly_epsilon[is.finite(pnt_in_poly_epsilon)]
  assert(pnt_in_poly_epsilon, len = 1, typ = "numeric")
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  force_default = as.logical(force_default); assert(force_default, len = 1, alw = c(TRUE, FALSE))
  recursive = as.logical(recursive); assert(recursive, len = 1, alw = c(TRUE, FALSE))
  endianness = cpp_checkTIFF(fileName)
  IFD = getIFD(fileName = fileName, offsets = "first", trunc_bytes = 8, verbose = verbose, verbosity = verbosity, force_trunc = FALSE, bypass = FALSE, ...)
  fileName = normalizePath(fileName, winslash = "/", mustWork = FALSE)
  title_progress = basename(fileName)
  
  ##### Initializes values
  merged = FALSE
  Files = list()
  features_def = list()
  features = data.frame()
  pops = list()
  plots = list()
  regions = list()
  stats = data.frame()
  onefile = FALSE
  V = NULL
  
  # TODO ask AMNIS how merged are defined / can be checked
  # Merged CIF Files
  if(!is.null(IFD[[1]]$tags[["33029"]])) {
    if(IFD[[1]]$tags[["33029"]]$byt != 0) V = strsplit(as.character(getFullTag(IFD = IFD, which = 1, tag="33029")), split = "|", fixed = TRUE)[[1]]
    LV = length(V)
    if(LV > 1) merged = TRUE
    if(merged & recursive) {
      Files = lapply(1:LV, FUN = function(i) {
        f = normalizePath(paste(dirname(fileName),basename(V[i]),sep="/"), winslash = "/", mustWork = FALSE)
        if(file.exists(f)) {
          ExtractFromXIF(fileName = f, pnt_in_poly_algorithm = pnt_in_poly_algorithm, pnt_in_poly_epsilon = pnt_in_poly_epsilon,
                         force_default = force_default, verbose = verbose, verbosity = verbosity, ...)
        } else {
          warning(paste0("Can't find sub-file defining merged:\n", f), call. = FALSE, immediate. = TRUE)
          out = list("description"=list(), "fileName"=V[i], "fileName_image"=V[i], "features"=features, "features_def"=features_def, "graphs"=plots, "pops"=pops, "regions"=regions, "images"=data.frame(), "offsets"=c(), "stats"=stats)
          attr(out, "class") <- c("IFC_data", "Merged")
          return(out)
        }
      })
    } else {
      onefile = TRUE
    }
  }
  # Merged RIF Files
  if(!is.null(IFD[[1]]$tags[["33030"]])) {
    if(IFD[[1]]$tags[["33030"]]$byt != 0) V = strsplit(as.character(getFullTag(IFD = IFD, which = 1, tag="33030")), split = "|", fixed = TRUE)[[1]]
    LV = length(V)
    if(LV > 1) merged = TRUE
    if(merged & recursive) {
      Files = lapply(1:LV, FUN = function(i) {
        f = normalizePath(paste(dirname(fileName),basename(V[i]),sep="/"), winslash = "/", mustWork = FALSE)
        if(file.exists(f)) {
          ExtractFromXIF(fileName = f, pnt_in_poly_algorithm = pnt_in_poly_algorithm, pnt_in_poly_epsilon = pnt_in_poly_epsilon,
                         force_default = force_default, verbose = verbose, verbosity = verbosity, ...)
        } else {
          warning(paste0("Can't find sub-file defining merged:\n", f), call. = FALSE, immediate. = TRUE)
          out = list("description"=list(), "fileName"=V[i], "fileName_image"=V[i], "features"=features, "features_def"=features_def, "graphs"=plots, "pops"=pops, "regions"=regions, "images"=data.frame(), "offsets"=c(), "stats"=stats)
          attr(out, "class") <- c("IFC_data", "Merged")
          return(out)
        }
      })
    } else {
      onefile = TRUE
    }
  }
  tmp = read_xml(getFullTag(IFD = IFD, which = 1, "33027"), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  acquisition = list("Illumination"=lapply(as_list(xml_find_first(tmp, "//Illumination")), unlist),
                     "Imaging"=lapply(as_list(xml_find_first(tmp, "//Imaging")), unlist),
                     "Display"=lapply(as_list(xml_find_first(tmp, "//Display")), unlist))
  
  infos=list("in_use"=as.logical(as.numeric(unlist(strsplit(acquisition$Imaging[["ChannelInUseIndicators_0_11"]], " ", useBytes = TRUE, fixed=TRUE)))),
             "brightfield"=list("channel"=as.logical(as.numeric(unlist(strsplit(acquisition$Illumination[["BfLedIndicators_0_11"]], " ", useBytes = TRUE, fixed=TRUE)))),
                                "power"=as.logical(as.numeric(acquisition$Illumination[["BFOnOff"]])),
                                "intensity"=as.numeric(acquisition$Illumination[["BFIntensity"]])))
  infos$volume = as.numeric(IFD[[1]]$tags[["33073"]]$map)
  
  is.binary = IFD[[1]]$tags[["33082"]]$map != 0
  if(length(is.binary)==0) {is.binary=FALSE}
  infos$collectionmode = as.numeric(acquisition$Illumination[["CollectionMode"]])
  
  if(length(acquisition$Imaging[["DafFile"]])!=0) {
    if(acquisition$Imaging[["DafFile"]]!="") {
      tmp = read_xml(acquisition$Imaging[["DafFile"]], options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
      force_default = FALSE
    } else {
      if(!force_default) stop("can't determine acquisition information")
    }
  } else {
    if(!force_default) stop("can't determine acquisition information")
  }
  if(force_default) {
    col_tmp = c("DarkOrchid", "Lime", "Yellow", "DarkOrange", "Red", "DeepPink")
    col_tmp = rep(col_tmp, 2)
    node = lapply(1:12, FUN=function(i) {
      if(infos$brightfield$channel[i]) {
        if(infos$collectionmode == 1) {
          sprintf('<image name="Ch%s" color="White" physicalChannel="%s" xmin="450" xmax="1000" xmid="725" ymid="127" scalemin="445" scalemax="1005" tokens="" baseimage="" function="" saturation="Cyan"/>', sprintf("%02.0f", i), i-1)
        } else {
          sprintf('<image name="Ch%s" color="White" physicalChannel="%s" xmin="100" xmax="300" xmid="200" ymid="127" scalemin="95" scalemax="305" tokens="" baseimage="" function="" saturation="Cyan"/>', sprintf("%02.0f", i), i-1)
        }
      } else {
        if(infos$collectionmode == 1) {
          sprintf('<image name="Ch%s" color="%s" physicalChannel="%s" xmin="0" xmax="4095" xmid="2047" ymid="127" scalemin="0" scalemax="4095" tokens="" baseimage="" function="" saturation="Cyan"/>', sprintf("%02.0f", i), col_tmp[i], i-1)
        } else {
          sprintf('<image name="Ch%s" color="%s" physicalChannel="%s" xmin="0" xmax="1023" xmid="511" ymid="127" scalemin="0" scalemax="1023" tokens="" baseimage="" function="" saturation="Cyan"/>', sprintf("%02.0f", i), col_tmp[i], i-1)
        }
      }
    })
    tmp = read_xml(paste0("<Images>",paste0(node, collapse=""),"</Images>"), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  } 
  ##### extracts description
  description=list("Assay"=xml_attrs(xml_find_all(tmp, "//Assay")),
                   "ID"=list(c("file"=fileName, "objcount"=IFD[[1]]$tags[["33018"]]$map)),
                   "Images"=xml_attrs(xml_find_all(tmp, "//image")),
                   "masks"=xml_attrs(xml_find_all(tmp, "//mask")))
  description=lapply(description, FUN=function(x) {as.data.frame(do.call(what="rbind", x), stringsAsFactors=FALSE)})
  description$Images = description$Images[order(description$Images$physicalChannel),]
  class(description$masks) <- c(class(description$masks), "IFC_masks")
  chan_number = sum(infos$in_use)
  obj_number = as.integer(description$ID$objcount)
  description$ID$objcount = obj_number
  checksum = checksumXIF(fileName)
  
  fileName_image = paste(dirname(fileName),description$ID$file,sep="/")
  if(file.exists(fileName_image)) {
    fileName_image = normalizePath(fileName_image, winslash = "/")
  } else {
    fileName_image = description$ID$file
  }
  
  for(i in c("physicalChannel","xmin","xmax","xmid","ymid","scalemin","scalemax")) description$Images[, i] = as.numeric(description$Images[, i])
  description$Images$physicalChannel = description$Images$physicalChannel + 1
  col = description$Images[,"color"]
  col[col=="Teal"] <- "Cyan4"
  col[col=="Green"] <- "Green4"
  col[col=="Lime"] <- "Chartreuse"
  description$Images[,"color"] <- col
  if("saturation"%in%names(description$Images)) {
    col = description$Images[,"saturation"]
    col[col=="Teal"] <- "Cyan4"
    col[col=="Green"] <- "Green4"
    col[col=="Lime"] <- "Chartreuse"
    description$Images[,"saturation"] <- col
  }
  
  if(extract_stats & !extract_features) {
    extract_features = TRUE
    message("'extract_features' has been forced to TRUE to extract statistics.")
  }
  
  if(extract_features) {
    ##### extracts features definition
    features_def=lapply(xml_attrs(xml_find_all(tmp, "//UDF")), FUN=function(x) as.list(x))
    feat_number = length(features_def)
    
    toread=file(description = fileName, open = "rb")
    on.exit(close(toread), add = TRUE)
    ##### extracts features values
    
    title_progress = basename(fileName)
    tryCatch({
      features = list()
      if(is.binary) {
        seek(toread, ifelse(merged | onefile, 
                            ifelse(length(IFD[[1]]$tags[["33083"]]$map)==0, stop("can't find pointer '33083' to extract features"), IFD[[1]]$tags[["33083"]]$val), 
                            ifelse(length(IFD[[1]]$tags[["33080"]]$map)==0, stop("can't find pointer '33080' to extract features"), IFD[[1]]$tags[["33080"]]$val)))
        obj_number_r = readBin(toread, what = "double", size = 4, n = 1, endian = endianness)
        feat_number_r = readBin(toread, what = "double", size = 4, n = 1, endian = endianness)
        if((length(obj_number_r) == 0) || (length(feat_number_r) == 0)) stop(fileName, "\nBinary features is of length 0")
        if(!(merged | onefile)) if(IFD[[1]]$tags[["33018"]]$map != obj_number_r) stop(fileName, "\nMismatch in object number")
        if(display_progress) {
          pb = newPB(session = dots$session, min = 0, max = obj_number_r, initial = 0, style = 3)
          tryCatch({
          features=lapply(1:obj_number_r, FUN=function(i_obj) {
            setPB(pb, value = i_obj, title = title_progress, label = "extracting features values (binary)")
            # fid=readBin(toread, "integer", n = 1, endian = endianness) # no fid found
            fv=readBin(toread, "double", size = 4, n = feat_number_r, endian = endianness)
            return(fv)
          })
        }, error = function(e) {
          stop(e$message)
        }, finally = endPB(pb))
        } else{
          features=lapply(1:obj_number_r, FUN=function(i_obj) {
            # fid=readBin(toread, "integer", n = 1, endian = endianness) # no fid found
            fv=readBin(toread, "double", size = 4, n = feat_number_r, endian = endianness)
            return(fv)
          }) 
        }
      } else { 
        # TODO
        stop("\nCan't deal with non-binary features")
        # feat_number=length(features)
        if(display_progress) {
          pb = newPB(session = dots$session, min = 0, max = feat_number, initial = 0, style = 3)
          tryCatch({
          features=lapply(1:feat_number,FUN=function(i) {
            setPB(pb, value = i, title = title_progress, label = "extracting features values (non-binary)")
            val = suppressWarnings(as.numeric(strsplit(features[i],"|", useBytes = TRUE, fixed=TRUE)[[1]]))
            val[is.na(val)] <- NaN
            val
          })
        }, error = function(e) {
          stop(e$message)
        }, finally = endPB(pb))
        } else {
          features=lapply(1:feat_number,FUN=function(i) {
            val = suppressWarnings(as.numeric(strsplit(features[i],"|", useBytes = TRUE, fixed=TRUE)[[1]]))
            val[is.na(val)] <- NaN
            val
          })
        }
      }
    }, error = function(e) {
      message(paste0(e$message, ". Features values were not exported"))
    })
    if(length(features) != 0) { # means features were extracted
      features = as.data.frame(do.call(what = "rbind", args = features), stringsAsFactors = FALSE)
      features_names = sapply(features_def, FUN=function(x) x$name)
      def_def = sapply(features_def, FUN=function(x) x$def)
      names(features_def) = features_names
      names(features) = features_names
      if(!("Object Number"%in%features_names)) {
        features_names = c(features_names, "Object Number")
        features$`Object Number` = 0:(nrow(features)-1)
        features_def = c(features_def, "Object Number" = list(name = "Object Number", type = "single", userfeaturetype = "No Parameters", def = "Object Number"))
      } else { # try to define unique object id number based on "Object Number","Camera Timer","Camera Line Number" if present
        if(all(c("Object Number","Camera Timer","Camera Line Number") %in% def_def)) {
          ids = rle(apply(sapply(c("Object Number","Camera Timer","Camera Line Number"),
                                 FUN=function(col) {
                                   foo = rle(features[,which(def_def == col)[1]])
                                   bar = lapply(1:length(foo$lengths), FUN = function(i) rep(i-1, times = foo$lengths[i]))
                                   unlist(bar)
                                 }), 1, sum))
          unique_id = unlist(sapply(1:length(ids$lengths), FUN=function(i) rep(i-1, times = ids$lengths[i])))
          features[, "Raw Number"] = features[, "Object Number"]
          features_def = c(features_def, "Raw Number" = list(name = "Raw Number", type = "single", userfeaturetype = "No Parameters", def = "Raw Number"))
          features[, "Object Number"] = unique_id
        }
      }
      if(any(duplicated(features$`Object Number`))) {
        features$`Object Number` = 0:(nrow(features)-1)
        warning(paste0("found duplicated objects when reading file: ", fileName))
      }
      rownames(features) = 0:(nrow(features)-1)
      class(features) <- c(class(features),"IFC_features")
      class(features_def) <- c(class(features),"IFC_features_def")
      
      ##### extracts graphs information
      plots=lapply(xml_attrs(xml_find_all(tmp, "//Graph")), FUN=function(x) as.list(x))
      if(length(plots)!=0) {
        plots_tmp=lapply(plots, FUN=function(plot) {
          pat=paste0("//Graph[@xlocation='",plot$xlocation,"'][@ylocation='",plot$ylocation,"']")
          sapply(c("Legend","BasePop","GraphRegion","ShownPop"), FUN=function(i_subnode){
            lapply(xml_attrs(xml_find_all(tmp, paste(pat,i_subnode,sep="//"))), FUN=function(x) as.list(x))
          })
        })
        plots=mapply(plots, plots_tmp, FUN = append, SIMPLIFY = FALSE)
        plots_tmp=c("xlocation","ylocation","scaletype","xmin","xmax","ymin","ymax","axislabelsfontsize","axistickmarklabelsfontsize",
                    "graphtitlefontsize","regionlabelsfontsize","bincount","histogramsmoothingfactor","xsize","ysize","splitterdistance")
        plots=lapply(plots, FUN=function(x) {replace(x, plots_tmp, lapply(x[plots_tmp], as.numeric))})
        plot_order=sapply(plots, FUN=function(i_plot) as.numeric(i_plot[c("xlocation", "ylocation")]))
        plots=plots[order(plot_order[1,],plot_order[2,])]
        plots=plots[order(plot_order[2,])]
        rm(list=c("plots_tmp", "plot_order"))
      }
      
      ##### TODO, add something for ChannelImage, ObjectFeatureControl, StatisticsControl
      
      ##### extracts regions information
      regions=lapply(xml_attrs(xml_find_all(tmp, "//Region")), FUN=function(x) as.list(x))
      if(length(regions) != 0) {
        names(regions)=lapply(regions, FUN=function(x) x$label)
        regions_tmp=c("cx","cy")
        regions=lapply(regions, FUN=function(x) {replace(x, regions_tmp, lapply(x[regions_tmp], as.numeric))})
        regions_tmp=lapply(regions, FUN=function(i_region) {
          pat=paste0("//Region[@label='",i_region$label,"']//axy")
          axy=do.call(cbind, args = xml_attrs(xml_find_all(tmp, pat)))
          list(x=as.numeric(axy["x",]), y=as.numeric(axy["y",]))
        })
        regions=mapply(FUN = append, regions, regions_tmp, SIMPLIFY = FALSE)
        rm(regions_tmp)
        ##### changes unknown color names in regions
        for(i in 1:length(regions)) {
          if(regions[[i]]$color=="Teal") {regions[[i]]$color="Cyan4"}
          if(regions[[i]]$color=="Green") {regions[[i]]$color="Green4"}
          if(regions[[i]]$color=="Lime") {regions[[i]]$color="Chartreuse"}
          if(regions[[i]]$lightcolor=="Teal") {regions[[i]]$lightcolor="Cyan4"}
          if(regions[[i]]$lightcolor=="Green") {regions[[i]]$lightcolor="Green4"}
          if(regions[[i]]$lightcolor=="Lime") {regions[[i]]$lightcolor="Chartreuse"}
        }
      }
      class(regions) <- "IFC_regions"
      
      ##### extracts populations information
      pops=lapply(xml_attrs(xml_find_all(tmp, "//Pop")), FUN=function(x) as.list(x))
      if(length(pops)>0) {
        names(pops)=lapply(pops, FUN=function(x) x$name)
        pops_=lapply(pops, FUN=function(i_pop) {
          pat=paste0("//Pop[@name='",i_pop$name,"']//ob")
          list(obj=as.numeric(unlist(xml_attrs(xml_find_all(tmp, pat)))))
        })
        pops=mapply(FUN = append, pops, pops_, SIMPLIFY = FALSE)
        rm(pops_)
      }
      class(pops) <- "IFC_pops"
      
      #####  retrieve name(s) of graphical population created by region applied in graph
      if(length(plots) > 0) {
        plots = lapply(plots, FUN = function(g) {
          if(length(g$GraphRegion) != 0) {
            N = sapply(g$GraphRegion, FUN = function(r) {
              foo = sapply(pops,
                           FUN = function(p) {
                             bar = (p$type == "G") && 
                               (p$region == r$name) && 
                               (p$base %in% unique(unlist(lapply(g$BasePop, FUN = function(b) b$name)))) &&
                               (g$f1 == p$fx)
                             if(regions[[r$name]]$type != "line") bar = bar && (g$f2 == p$fy)
                             return(bar)
                           })
              return(names(which(foo)))
            })
            g$GraphRegion$def = N
          }
          return(g)
        })
      }
      class(plots) <- "IFC_graphs"
    } else {
      features = data.frame()
    }
    
    l = length(pops)
    if(l>0) {
      ###### scrambles pops (for testing)
      # pops = pops[sample.int(length(pops))]
      
      ##### extracts populations dependencies/affiliations.
      ##### reorders pops
      pops = popsOrderNodes(popsGetAffiliation(pops))
      ##### determines which object belongs to each population and changes styles and colors
      pops = popsWithin(pops = pops,
                        regions = regions,
                        features = features,
                        pnt_in_poly_algorithm = pnt_in_poly_algorithm,
                        pnt_in_poly_epsilon = pnt_in_poly_epsilon,
                        display_progress = display_progress,
                        title_progress = title_progress, ...)
      
      if(extract_stats) {
        stats = data.frame(stringsAsFactors = FALSE, check.rows = FALSE, check.names = FALSE, t(sapply(names(pops), FUN=function(p) {
          count = sum(pops[[p]]$obj)
          base = pops[[p]]$base
          type = pops[[p]]$type
          if(base=="") base = "All"
          parent = sum(pops[[base]]$obj)
          c("type" = type, "parent" = base, "count" = count, "perc_parent" = count/parent*100, "perc_tot" = count/obj_number*100)
        })))
        stats[,3] = as.numeric(stats[,3])
        stats[,4] = as.numeric(stats[,4])
        stats[,5] = as.numeric(stats[,5])
      }
    }
  } else {
    features = data.frame()
  }
  # Initializes and extracts offsets if needed
  offsets = NULL
  if(extract_images) {
    if(!extract_offsets)
      message("'extract_offsets' has been forced to TRUE to extract images")
    extract_offsets = TRUE
  }
  if(extract_offsets) {
    if(extract_images) {
      offsets = suppressMessages(getOffsets(fileName, fast = fast, display_progress = display_progress))
    } else {
      offsets = getOffsets(fileName, fast = fast, display_progress = display_progress)
    }
  }
  
  images = data.frame()
  if(extract_images) {
    images = getImagesValues(fileName = fileName, offsets = offsets, fast = fast, display_progress = display_progress, ...)
    if(fast) {
      N = nchar(sprintf("%1.f",abs(obj_number-1)))
      tmp = c(paste0("img_", sprintf(paste0("%0",N,".f"), images$id)), paste0("msk_", sprintf(paste0("%0",N,".f"), images$id)))
      if(!all(offsets[tmp] == c(images$imgIFD, images$mskIFD))) {
        warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
      }
    }
  }
  
  ans = list("description"=description, "fileName"=fileName, "fileName_image"=fileName, "features"=features, "features_def"=features_def, "graphs"=plots, "pops"=pops, "regions"=regions, "images"=images, "offsets"=offsets, "stats"=stats, "checksum" = checksum)
  attr(ans, "class") <- c("IFC_data")
  if(merged) {
    out = c("Merged"=Files, ans)
    attr(out, "class") <- c("IFC_data", "Merged")
    return(out)
  }
  return(ans)
}
