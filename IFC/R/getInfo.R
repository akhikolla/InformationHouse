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

#' @title IFC File Information Extraction
#' @description
#' Retrieves rich information from RIF, CIF and DAF files.
#' @param fileName path to file..
#' @param from whether to extract information from 'acquisition' or 'analysis'. Default is 'analysis'.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @param verbosity quantity of information print to console when verbose is TRUE; 1: normal, 2: rich. Default is 1.
#' @param warn whether to send warning message when trying to read 'analysis' information from a 'rif' file. Default is TRUE.
#' @param force_default when display information can't be retrieved whether to use default values. Default is TRUE.
#' @param cifdir the path of the directory to initially look to cif file. Default is dirname(fileName). Only apply when 'fileName' is a .daf file.
#' @param ntry number of times \code{\link{getInfo}} will be allowed to find corresponding cif file. Default is +Inf. Only apply when 'fileName' is a .daf file.
#' If cif can't be found, but 'ntry' is reached, then an error will be thrown.
#' @param ... other arguments to be passed.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   info <- getInfo(fileName = file_daf, from = "analysis")
#'   ## show some information
#'   print(info$Images)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return a list of information (open .daf file in an text editor for more details) about input fileName of class `IFC_info` and `acquistion` or `analysis`, whose members are:\cr
#' -objcount, number of object in file,\cr
#' -date, date of file creation,\cr
#' -instrument, instrument identification,\cr
#' -sw_raw, version of software for raw data,\cr
#' -sw_processed, version of software for processed data,\cr
#' -channelwidth, default channel width in pixel,\cr
#' -in_use, channels used,\cr
#' -brightfield, whether brightfield is applied on channels and its intensity,\cr
#' -illumination, laser illumination parameters,\cr
#' -collectionmode, the collection mode,\cr
#' -magnification, magnification used,\cr
#' -coremode, the core mode,\cr
#' -CrossTalkMatrix. compensation matrix applied,\cr
#' -ChannelPresets, channel preset,\cr
#' -ImageDisplaySettings, image display settings,\cr
#' -Images, information about colors, range and channels,\cr
#' -masks, masks defined,\cr
#' -ViewingModes, modes of visualization,\cr
#' -checksum, checksum computed,\cr
#' -Merged_rif, character vector of path of files used to create rif, if input file was a merged,\cr
#' -Merged_cif, character vector of path of files used to create cif, if input file was a merged,\cr
#' -fileName, path of fileName input,\cr
#' -fileName_image, path of fileName_image.
#' @export
getInfo <- function(fileName,
                    from = c("acquisition","analysis")[2],
                    verbose = FALSE, 
                    verbosity = 1, 
                    warn = TRUE, 
                    force_default = TRUE,
                    cifdir = dirname(fileName), 
                    ntry = +Inf,
                    ...) {
  dots = list(...)
  if(missing(fileName)) stop("'fileName' can't be missing")
  if(missing(fileName)) stop("'fileName' can't be missing")
  tmp = duplicated(fileName)
  if(any(tmp)) {
    warning(paste0("duplicated files have been removed from 'fileName': ","\n-", paste0(fileName[tmp],collapse="\n-")))
    fileName = fileName[!tmp]
  }
  if(length(fileName) != 1) stop("'fileName' should be of length 1")
  if(!file.exists(fileName)) stop(paste0("can't find ",fileName))
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf", "cif", "rif"))
  fileName = normalizePath(fileName, winslash = "/", mustWork = FALSE)
  assert(from, len = 1, alw = c("acquisition","analysis"))
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  verbosity = as.integer(verbosity); assert(verbosity, len = 1, alw = c(1, 2))
  warn = as.logical(warn); assert(warn, len = 1, alw = c(TRUE, FALSE))
  force_default = as.logical(force_default); assert(force_default, len = 1, alw = c(TRUE, FALSE))
  cifdir = na.omit(as.character(cifdir)); assert(cifdir, len = 1, typ = "character")
  ntry = na.omit(as.numeric(ntry)); assert(ntry, len = 1, typ = "numeric")
  if(ntry < 0) ntry = 0
  if(warn & file_extension == "rif" & from == "analysis") warning("Only information from 'acquisition' can be retrieved from 'rif' file", call. = FALSE, immediate. = TRUE)
  if(file_extension == "daf") {
    toskip = cpp_scanFirst(fname = fileName, target = "</Assay>", start = 0, end = 0)
    if(toskip == 0) stop(paste0(fileName, "\ndoes not seem to be well formatted: </Assay> not found")) 
    toskip = toskip + nchar("</Assay>") - 1
    tmp_daf = read_xml(readBin(con = fileName, what = "raw", n = toskip), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
    cname = xml_attr(xml_find_first(tmp_daf, "//SOD"), attr = "file")
    fcs = xml_attr(xml_find_first(tmp_daf, "//FCS"), attr = "file")
    if(!is.na(fcs)) stop("can't extract information from .daf generated from .fcs")
    found = FALSE
    checksum = checksumDAF(fileName)
    fileName_image = file.path(cifdir, paste0(splitf(cname)[c("short","ext")], collapse = ".")) # look in cifdir 1st
    if(file.exists(fileName_image)) {
      if(checksumXIF(fileName_image) == checksum) found = TRUE
    } else {
      fileName_image = cname
    }
    if((!found) && file.exists(fileName_image)) {
      if(checksumXIF(fileName_image) == checksum) found = TRUE
    }
    
    while((interactive() && (ntry > 0) && (!found))) {
      message(paste0("daf file does not refer to: ", fileName_image))
      old_wd = getwd()
      on.exit(setwd(old_wd), add= TRUE)
      setwd(dirname(fileName))
      if(.Platform$OS.type == "windows") {
        fileName_image = choose.files(caption = paste0("Looking for: ", basename(cname)), multi = FALSE, filters = cbind("Compensated Image File (*.cif)", "*.cif"))
      } else {
        fileName_image = file.choose()
      }
      if(file.exists(fileName_image)) if(getFileExt(fileName_image)=="cif") if(checksumXIF(fileName_image) == checksum) {
        found = TRUE
        break;
      } 
      ntry = ntry - 1
    }
    if(!found) stop("can't extract information")
    fileName_image = normalizePath(fileName_image, winslash = "/")
    IFD = getIFD(fileName = fileName_image, offsets = "first", trunc_bytes = 8, force_trunc = FALSE, verbose = verbose, verbosity = verbosity, bypass = TRUE, ...)
  }
  if(file_extension == "cif" | file_extension == "rif") {
    fileName_image = fileName
    IFD = getIFD(fileName = fileName_image, offsets = "first", trunc_bytes = 8, force_trunc = FALSE, verbose = verbose, verbosity = verbosity, bypass = TRUE, ...)
  }
  Merged_cif = character()
  Merged_rif = character()
  if(!is.null(IFD[[1]]$tags[["33029"]])) {
    if(IFD[[1]]$tags[["33029"]]$byt != 0) Merged_cif = strsplit(as.character(getFullTag(IFD = IFD, which = 1, tag="33029")), "|", fixed = TRUE)[[1]]
  }
  if(!is.null(IFD[[1]]$tags[["33030"]])) {
    if(IFD[[1]]$tags[["33030"]]$byt != 0) Merged_rif = strsplit(as.character(getFullTag(IFD = IFD, which = 1, tag="33030")), "|", fixed = TRUE)[[1]]
  }
  if(file_extension == "daf" & from == "acquisition") file_extension = "cif"
  tmp_acq = read_xml(getFullTag(IFD = IFD, which = 1, "33027"), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  
  acquisition = list("Illumination"=lapply(as_list(xml_find_first(tmp_acq, "//Illumination")), unlist),
                     "Fluidics"=lapply(as_list(xml_find_first(tmp_acq, "//Fluidics")), unlist),
                     "Imaging"=lapply(as_list(xml_find_first(tmp_acq, "//Imaging")), unlist),
                     "Display"=lapply(as_list(xml_find_first(tmp_acq, "//Display")), unlist))
  
  infos = list("objcount" = IFD[[1]]$tags[["33018"]]$map, # should not exceed 4 bytes
               "date"=getFullTag(IFD = IFD, which = 1, "33004"),
               "instrument"=getFullTag(IFD = IFD, which = 1, "33006"),
               "sw_raw"=getFullTag(IFD = IFD, which = 1, "33069"),
               "sw_process"=getFullTag(IFD = IFD, which = 1, "33066")) 
  # determines channelwidth, very important for objectExtract() when force_width = TRUE
  # prefer using channelwidth extracted from ifd dedicated tag (=tag 33009) rather than the one from parsing ASSISTdb (=tag 33064)
  # TODO ask AMNIS the rules for extracting channelwidth
  channelwidth1 = IFD[[1]]$tags[["33009"]]$map # should not exceed 4 bytes
  tmp_ins = read_xml(getFullTag(IFD = IFD, which = 1, tag ="33064"), options = c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  channelwidth2 = as.numeric(xml_text(xml_find_first(tmp_ins, xpath = "//ChannelWidth")))
  
  tryCatch({
    # information in acquisition db
    lasers_nodes = grep("Filter", grep("ExLaser", names(acquisition[["Illumination"]]), value = TRUE), value = TRUE, invert = TRUE)
    lasers_on = acquisition$Illumination[grep("PowerOn", lasers_nodes, value = TRUE)]
    lasers_on = as.logical(as.integer(unlist(lasers_on[order(as.integer(gsub("\\D+", "", names(lasers_on))))])))
    lasers_power = acquisition$Illumination[grep("IntensityWatts", lasers_nodes, value = TRUE)]
    lasers_power = as.numeric(sprintf("%.2f", as.numeric(unlist(lasers_power[order(as.integer(gsub("\\D+", "", names(lasers_power))))]))))
    if(length(lasers_power) != length(lasers_on)) lasers_power = rep(NA, length(lasers_on))
    # information in instrument db
    ins_lasers = lapply(as_list(xml_find_first(tmp_ins, "Illumination")), unlist)
    lasers_nodes = grep("ExLaser", names(ins_lasers), value = TRUE)
    
    lasers_present = ins_lasers[grep("Filter|LAF", grep("Present", lasers_nodes, value = TRUE), value = TRUE, invert = TRUE)]
    lasers_present = as.logical(as.integer(unlist(lasers_present[order(as.integer(gsub("\\D+", "", names(lasers_present))))])))
    if(length(lasers_present) != length(lasers_on)) lasers_present = rep(NA, length(lasers_on))
    
    lasers_wavelength = ins_lasers[grep("Wavelength", lasers_nodes, value = TRUE)]
    lasers_wavelength = as.integer(unlist(lasers_wavelength[order(as.integer(gsub("\\D+", "", names(lasers_wavelength))))]))
    if(length(lasers_wavelength) != length(lasers_on)) lasers_wavelength = rep(NA, length(lasers_on))
    
    lasers_minpow = ins_lasers[grep("MinPower", lasers_nodes, value = TRUE)]
    lasers_minpow = as.numeric(unlist(lasers_minpow[order(as.integer(gsub("\\D+", "", names(lasers_minpow))))]))
    if(length(lasers_minpow) != length(lasers_on)) lasers_minpow = rep(NA, length(lasers_on))
    
    lasers_maxpow = ins_lasers[grep("MaxPower", lasers_nodes, value = TRUE)]
    lasers_maxpow = as.numeric(unlist(lasers_maxpow[order(as.integer(gsub("\\D+", "", names(lasers_maxpow))))]))
    if(length(lasers_maxpow) != length(lasers_on)) lasers_maxpow = rep(NA, length(lasers_on))
    
    illumination = data.frame("installed" = lasers_present,
                              "wavelength" = lasers_wavelength, 
                              "powered" = lasers_on, 
                              "power" = lasers_power, 
                              "min" = lasers_minpow, 
                              "max" = lasers_maxpow, stringsAsFactors = FALSE)
  }, error = function(e) {
    illumination = data.frame(matrix(NA, ncol = 6, nrow = 0))
    colnames(illumination) = c("installed","wavelength", "powered" ,"power" ,"min","max")
  })

  infos$channelwidth = channelwidth1
  if(length(channelwidth1)==0) infos$channelwidth = channelwidth2
  if(length(channelwidth1)!=0) if(is.na(channelwidth1)) infos$channelwidth = channelwidth2
  if(length(channelwidth1)!=0) if(channelwidth1==0) infos$channelwidth = channelwidth2
  infos$in_use = as.logical(as.numeric(unlist(strsplit(acquisition$Imaging[["ChannelInUseIndicators_0_11"]], " ", useBytes = TRUE, fixed=TRUE))))
  infos$brightfield = list("channel"=as.logical(as.numeric(unlist(strsplit(acquisition$Illumination[["BfLedIndicators_0_11"]], " ", useBytes = TRUE, fixed=TRUE)))),
                           "power"= as.logical(as.numeric(acquisition$Illumination[["BFOnOff"]])),
                           "intensity" = as.numeric(acquisition$Illumination[["BFIntensity"]]))
  infos$illumination = illumination
  # TODO ask AMNIS, if collectionmode is the good variable that determines default information
  infos$collectionmode = as.numeric(acquisition$Illumination[["CollectionMode"]])
  infos$magnification = as.numeric(acquisition$Imaging[["Magnification"]])
  infos$coremode = as.numeric(acquisition$Fluidics[["CoreMode"]])
  if(from == "analysis" & file_extension != "rif") {
    if(length(IFD[[1]]$tags[["33020"]]$map)!=0) {
      cross = getFullTag(IFD = IFD, which = 1, tag = "33020")
      infos$CrossTalkMatrix = matrix(cross, nrow = sqrt(length(cross)), byrow = TRUE)
    } else {
      infos$CrossTalkMatrix = NULL
    }
  } else {
    if(length(acquisition$Imaging$InspireCrossTalkMatrix) == 0) {
      if(length(IFD[[1]]$tags[["33020"]]$map)!=0) {
        cross = getFullTag(IFD = IFD, which = 1, tag = "33020")
        infos$CrossTalkMatrix = matrix(cross, nrow = sqrt(length(cross)), byrow = TRUE)
      } else {
        infos$CrossTalkMatrix = NULL
      }
    } else {
      infos$CrossTalkMatrix = matrix(as.numeric(strsplit(x = acquisition$Imaging$InspireCrossTalkMatrix, split=" ", fixed = TRUE)[[1]]), nrow = length(infos$in_use), byrow = TRUE)
    }
  }
  
  if(file_extension == "daf") { 
    tmp_last = tmp_daf
  } else {
    if(getFileExt(fileName) == "daf") rm(tmp_daf)
    if(length(acquisition$Imaging[["DafFile"]])!=0) {
      if(acquisition$Imaging[["DafFile"]]!="") {
        tmp_last = read_xml(acquisition$Imaging[["DafFile"]], options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
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
      tmp_last = read_xml(paste0("<Images>",paste0(node, collapse=""),"</Images>"), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
    } 
  }
  infos = c(infos, list("ChannelPresets" = to_list_node(xml_find_all(tmp_last, "//ChannelPresets")),
                        "ImageDisplaySettings" = to_list_node(xml_find_all(tmp_last, "//ImageDisplaySettings")),
                        "Images" = as.data.frame(do.call(what = "rbind", args = xml_attrs(xml_find_all(tmp_last, "//image"))), stringsAsFactors = FALSE),
                        "masks" = lapply(xml_attrs(xml_find_all(tmp_last, "//mask")), FUN=strsplit, split="|", fixed=TRUE)),
            "ViewingModes" = to_list_node(xml_find_all(tmp_last, "//ViewingModes")),
            "Merged_rif" = list(Merged_rif),
            "Merged_cif" = list(Merged_cif),
            "checksum" = checksumXIF(fileName_image),
            "fileName" = fileName,
            "fileName_image" = normalizePath(fileName_image, winslash = "/"))
  
  infos$Images = infos$Images[order(infos$Images$physicalChannel),]
  names(infos$masks) = sapply(infos$masks, FUN=function(x) x$name)
  class(infos$masks) <- c(class(infos$masks), "IFC_masks")
  if(length(infos$ViewingModes) != 0) names(infos$ViewingModes) = sapply(infos$ViewingModes, FUN=function(x) x$name)
  
  for(i in c("physicalChannel","xmin","xmax","xmid","ymid","scalemin","scalemax")) infos$Images[, i] = as.numeric(infos$Images[, i])
  infos$Images$physicalChannel = infos$Images$physicalChannel + 1
  infos$Images = infos$Images[order(infos$Images$physicalChannel), ]
  infos$Images$gamma = apply(infos$Images[,c("xmin", "xmax", "xmid", "ymid")], 1, cpp_computeGamma)
  col = infos$Images[,"color"]
  col[col=="Teal"] <- "Cyan4"
  col[col=="Green"] <- "Green4"
  col[col=="Lime"] <- "Chartreuse"
  infos$Images[,"color"] <- col
  if("saturation"%in%names(infos$Images)) {
    col = infos$Images[,"saturation"]
    col[col=="Teal"] <- "Cyan4"
    col[col=="Green"] <- "Green4"
    col[col=="Lime"] <- "Chartreuse"
    infos$Images[,"saturation"] <- col
  }
  attr(infos, "class") <- c("IFC_info",from)
  return(infos)
}
