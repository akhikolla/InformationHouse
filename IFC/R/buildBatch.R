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

#' @title Batch Builder
#' @description
#' Prepares XML node for \code{\link{ExportToBATCH}}.
#' @param files path of files to batch.
#' @param compensation path to compensation file.
#' @param analysis path to analysis file.
#' @param default_batch_dir directory where batches are stored.\cr
#' It can be found in IDEAS(R) software, under Options -> Application Defaults -> Directories -> Default Batch Report Files Directory.\cr
#' If missing, the default, it will be deduced from IDEAS(R) config file, However, if it can't be deduced then tempdir(check = TRUE) from \pkg{base} will be used.\cr
#' This argument takes precedence over 'config_file' and filling 'default_batch_dir' prevents the use of 'config_file' argument.
#' @param config_file path to IDEAS(R) config file.\cr
#' It may depends on IDEAS(R) software installation but one may use "C:/Users/\%USER\%/AppData/Roaming/Amnis Corporation/userconfig.xml".
#' @param name name of batch. Default is "Batch1".
#' @param use_acquisition whether to use acquisition as analysis template. Default is FALSE.
#' @param suffix suffix to add to files when batched. Default is "".
#' @param allow_channels_dissimilarity whether to allow building batch when all files were not acquired with same channels. Default is FALSE.
#' @param overwrite whether to overwrite files or not. Default is TRUE.
#' @param segment_rif size of file segmentation. Default is "None", for no segmentation.\cr
#' Allowed are "None", "100", "1K", "5K", "10K", "50K", "100K".
#' @param options A list of arguments to be passed.\cr
#' If missing, the default, options will be set to:\cr
#' -"Brightfield compensation"=TRUE,\cr
#' -"EDF deconvolution"=TRUE,\cr
#' -"Camera background"=TRUE,\cr
#' -"Spatial alignment"=TRUE.\cr
#' Allowed are TRUE or FALSE for all, excepted for 'Spatial aligment' which can also be path to .rif file.
#' @return a list containing batch information:\cr
#' -xml, the xml object to be written,\cr
#' -batch_dir, the directory where xml file is desired to be saved according to 'default_batch_dir' and 'config_file'.
#' @export
buildBatch <- function(files, compensation, analysis, default_batch_dir, config_file, 
                       name="Batch1", use_acquisition=FALSE, suffix="",
                       allow_channels_dissimilarity=FALSE, overwrite=TRUE, segment_rif="None",
                       options) {
  fileName_comp = ""
  fileName_align = ""
  # various checks
  if(missing(files)) stop("'files can't be missing")
  name = na.omit(as.character(name)); assert(name, len = 1, typ = "character")
  use_acquisition = as.logical(use_acquisition); assert(use_acquisition, len = 1, alw = c(TRUE,FALSE))
  suffix = na.omit(as.character(suffix)); assert(suffix, len = 1, typ = "character")
  allow_channels_dissimilarity = as.logical(allow_channels_dissimilarity); assert(allow_channels_dissimilarity, len = 1, alw = c(TRUE,FALSE))
  overwrite = as.logical(overwrite); assert(overwrite, len = 1, alw = c(TRUE,FALSE))
  allowed_chunk = c("None","100","1K","5K","10K","50K","100K")
  segment_rif = na.omit(as.character(segment_rif)); assert(segment_rif, len = 1, alw = allowed_chunk)
  # checks options
  if(missing(options)) {
    options=list("Brightfield compensation"=TRUE,
                 "EDF deconvolution"=TRUE,
                 "Camera background"=TRUE,
                 "Spatial alignment"=TRUE)
  }
  if(length(options[["Brightfield compensation"]])==0) {options[["Brightfield compensation"]]=TRUE}
  if(length(options[["EDF deconvolution"]])==0) {options[["EDF deconvolution"]]=TRUE}
  if(length(options[["Camera background"]])==0) {options[["Camera background"]]=TRUE}
  if(length(options[["Spatial alignment"]])==0) {options[["Spatial alignment"]]=TRUE}
  opt = lapply(names(options), FUN=function(x) {
    ##### can't be modified:
    # "Brightfield gains"
    # "flow speed normalization"
    # "Apply cell classifiers"
    # "Erase non-framed objects"
    # "Separate single objects"
    # "Remove clipped objects"
    # "Allow post processing"
    foo = as.logical(options[[x]])
    if(length(foo)!=1) stop(paste0("options$",x," should be of length 1"))
    if(x=="Spatial alignment") {
      if(is.na(foo)) {
        if(!file.exists(options[[x]])) {
          stop(paste0("options$`Spatial alignment` can't find file:",options[[x]]))
        } else {
          if(getFileExt(options[[x]])!="rif") stop("when provided options$`Spatial alignment` should be a rif file")
          cpp_checkTIFF(options[[x]])
          return("Y")
        }
      }
    } else {
      if(is.na(foo)) stop(paste0("options$",x," should be of class `logical`"))
    }
    return(ifelse(foo,"Y","N"))
  })
  names(opt) = names(options)
  if(opt[["Spatial alignment"]]=="Y") if(is.na(as.logical(options[["Spatial alignment"]]))) fileName_align = options[["Spatial alignment"]]
  
  # checks for batch_dir
  batch_dir = NULL
  is_tmp_dir = TRUE
  if(missing(default_batch_dir)) {
    if(!missing(config_file)) {
      config_file = na.omit(as.character(config_file));
      config_file = normalizePath(config_file, winslash = "/", mustWork = FALSE)
      if(length(config_file)!=1) {
        warning("can't find config file and default_batch_dir is not provided: tempdir(check = TRUE) will be used")
      } else {
        if(!file.exists(config_file)) {
          warning("can't find config file and default_batch_dir is not provided: tempdir(check = TRUE) will be used")
        } else {
          tmp_conf = read_xml(config_file, options = c("HUGE", "RECOVER", "NOENT", "NOBLANKS", "NSCLEAN"))
          batch_dir = xml_text(xml_find_first(tmp_conf, "//BatchDirectory"))
          is_tmp_dir = FALSE
        }
      }
    }
  } else {
    default_batch_dir = na.omit(as.character(default_batch_dir));
    if(length(default_batch_dir)!=1) {
      warning("when provided 'default_batch_dir' should be of length 1: tempdir(check = TRUE) will be used")
    } else{
      default_batch_dir = normalizePath(default_batch_dir, winslash = "/", mustWork = FALSE)
      if(!dir.exists(default_batch_dir)) {
        warning("'default_batch_dir' is invalid: tempdir(check = TRUE) will be used")
      } else {
        batch_dir = default_batch_dir
        is_tmp_dir = FALSE
      }
    }
  }
  if(length(batch_dir) == 0) batch_dir = tempdir(check = TRUE)
  attr(x = batch_dir, which = "tempdir") <- is_tmp_dir

  # checks analysis
  if(!missing(analysis)) {
    if(!file.exists(analysis)) stop(paste0("can't find analysis file:\n",analysis))
    a_Ext = getFileExt(analysis)
    if(!(a_Ext%in%c("ast","daf"))) stop("can't deal with analysis file, only .ast and .daf are supported")
  }
  # checks compensation
  if(!missing(compensation)) {
    if(!file.exists(compensation)) stop(paste0("can't find compensation file:\n",compensation))
    c_Ext = getFileExt(compensation)
    if(!(c_Ext%in%c("cif","ctm","daf"))) stop("can't deal with compensation file, only .cif, .ctm and .daf are supported for the moment")
  }
  # checks input files
  dup = duplicated(sapply(files, FUN=function(x) {
    x=tolower(normalizePath(x, winslash = "/"))
    gsub(getFileExt(x),"", x)
  }))
  if(any(dup)) {
    warning(paste0(paste0(paste0(files[dup],"\n"),collapse=""),"have been removed because batch would overwrite data being processed"))
    files = files[!dup]
  }
  if(length(files) == 0) stop("no file to batch using 'files' argument")
  info = lapply(files, getInfo, from = "acquisition")
  if(!allow_channels_dissimilarity) {
    num_channels = unique(sapply(info, FUN=function(x) nrow(x$Images)))
    if(length(num_channels) !=1 ) stop("'files' have been acquired on different amount of channels")
    ids_channels = sapply(info, FUN=function(x) c(sum(x$in_use), x$Images$physicalChannel[x$in_use])) # maybe to remove to allow
    if(typeof(ids_channels) == "list") stop("'files' have been acquired on different physical channels") # maybe to remove
    if(nlevels(as.factor(ids_channels[-1,]))!=ids_channels[1]) stop("'files' have been acquired on different physical channels") # maybe to remove
  }
  mag = unique(unlist(sapply(info, FUN=function(x) x$magnification)))
  if(length(mag)!=1) stop("'files' have been acquired with different magnifications")
  mag = paste0("Offsets",ifelse(mag==40, "",paste0(mag,"x")),"_Gen2_0_11")
  # Extracts information to build Nodes
  if(fileName_align=="") {
    offsets = c(X="",Y="")
  } else {
    IFD = getIFD(fileName = fileName_align, offsets = "first", trunc_bytes = 8, force_trunc = FALSE, bypass = FALSE)
    tmp_off = read_xml(getFullTag(IFD = IFD, which = 1, "33064"), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
    offsets = sapply(c("X","Y"), USE.NAMES = TRUE, FUN = function(off) {
      paste0(round(as.numeric(strsplit(xml_text(xml_find_first(tmp_off, xpath = paste0("//",off,mag)))," ",fixed=TRUE)[[1]]),2),collapse="|")
    })
  }
  if(!missing(compensation)) {
    if(c_Ext == "daf") {
      toskip = cpp_scanFirst(fname = compensation, target = "</Assay>", start = 0, end = 0)
      if(toskip == 0) stop(paste0(compensation, "\ndoes not seem to be well formatted: </Assay> not found")) 
      toskip = toskip + nchar("</Assay>") - 1
      tmp_comp = read_xml(readBin(compensation, what = "raw", n = toskip), options = c("HUGE", "RECOVER", "NOENT", "NOBLANKS", "NSCLEAN"))
      SOD = xml_attrs(xml_find_first(tmp_comp, xpath="//SOD"))
      
      fileName_comp = SOD["file"]
      if(!file.exists(fileName_comp)) {
        fileName_comp = paste(dirname(compensation), basename(fileName_comp), sep = "\\")
      }
      if(!file.exists(fileName_comp)) stop(paste0("can't find compensation file\n", fileName_comp),call. = FALSE)
    }
    if(c_Ext == "cif") fileName_comp = compensation
  }
  if(fileName_comp=="") {
    coeff = paste0(as.numeric(sapply(1:num_channels, FUN=function(x) 1:num_channels==x)), collapse="|")
  } else {
    IFD = getIFD(fileName = fileName_comp, offsets = "first", trunc_bytes = 8, force_trunc = FALSE, bypass = FALSE)
    coeff = paste0(getFullTag(IFD = IFD, which = 1, tag = "33020"),collapse="|")
  }
  chunk = c(1,100,1000,5000,10000,50000,100000)
  names(chunk) = allowed_chunk
  chunk = chunk[segment_rif]
  nodes = list(list(name="alignment",
                    attrs = list(name="",
                                 xoffsets=offsets["X"],
                                 yoffsets=offsets["Y"],
                                 ReferenceChannelCamera1="-1",
                                 ReferenceChannelCamera2="-1")),
               list(name="darkcurrent",
                    attrs = list(name="")),
               list(name="spectral",
                    attrs = list(name=ifelse(missing(compensation), "", normalizePath(compensation, winslash = "/")),
                                 coefficients=coeff)),
               list(name="template",
                    attrs=list(name=ifelse(missing(analysis), "", normalizePath(analysis, winslash = "/")),
                               useInspireAnalysis=ifelse(use_acquisition,"True","False"))),
               list(name="files",
                    .children=lapply(info, FUN=function(x) do.call(what=xml_new_node, args=list(name="file", attrs = list(name=x$fileName, objectCount=ifelse(getFileExt(x$fileName)!="rif","0",format(x$objcount,scientific=FALSE))))))),
               list(name="statfiles",
                    .children=lapply(unlist(lapply(info, FUN=function(x) {
                      if(getFileExt(x$fileName) == "daf") return(paste0(suffix,basename(x$fileName)))
                      if(getFileExt(x$fileName) == "cif") return(paste0(suffix,gsub("cif$","daf",basename(x$fileName_image))))
                      # segment_rif only applies to rif
                      foo = 1
                      if(segment_rif!="None") foo = ceiling(x$objcount/chunk)
                      if(foo == 1) {
                        foo = ""
                      } else {
                        foo = paste(paste0("_S", 1:foo),segment_rif,sep="_")
                      }
                      return(paste0(gsub("^(.*)\\.rif$","\\1",basename(x$fileName_image)),foo,suffix,".daf"))
                    })), FUN=function(n) xml_new_node(name="statfile", attrs = list(name=n)))))
  if(!missing(analysis)) {
    target = ifelse(a_Ext=="daf", "</Assay>", "</AssayTemplate>")
    toskip = cpp_scanFirst(fname = analysis, target = target, start = 0, end = 0)
    if(toskip == 0) {
      if(a_Ext=="daf") {
        stop(paste0(analysis, "\ndoes not seem to be well formatted: </Assay> not found")) 
      } else {
        stop(paste0(analysis, "\ndoes not seem to be well formatted: </Assay> not found")) 
      }
    }
    toskip = toskip + nchar(target) - 1
    tmp_ana = read_xml(readBin(analysis, what = "raw", n = toskip), options = c("HUGE", "RECOVER", "NOENT", "NOBLANKS", "NSCLEAN"))
    statoutputfile = list(name="statoutputfile",
                          attrs = list(name=xml_attr(xml_find_first(tmp_ana, xpath="//StatisticsReport"), attr = "title")))
    StatisticsReports = xml_find_first(tmp_ana, xpath="//StatisticsReports")
    if(length(xml_children(StatisticsReports)) != 0) {
      StatisticsReports = list(name="StatisticsReports",
                               .children = lapply(xml_children(StatisticsReports), FUN=function(x) {
                                 N = xml_name(x)
                                 A = xml_attrs(x)
                                 if(N=="StatisticsReport") A = c(list("filename"=A["title"]),A)
                                 xml_new_node(name=N, attrs = A)
                               }))
    } else {
      StatisticsReports = list(name="StatisticsReports", text = "")
    }
    nodes = c(nodes, list(statoutputfile), list(StatisticsReports))
  } else {
    statoutputfile = list(name="statoutputfile", attrs = "")
    nodes = c(nodes, list(statoutputfile))
  }
  xml = do.call(what = xml_new_node, args = list(name="batch",
                                            attrs=list(folder=name,
                                                       align=opt[["Spatial alignment"]],
                                                       dc=opt[["Camera background"]],
                                                       spectral=opt[["Brightfield compensation"]],
                                                       edf=opt[["EDF deconvolution"]],
                                                       output="current",
                                                       overwrite=ifelse(overwrite,"Y","N"),
                                                       suffix=suffix,
                                                       overridesuffix="",
                                                       alldone="Y"),
                                            .children=lapply(nodes, do.call, what=xml_new_node)))
  return(list(xml=xml, batch_dir=batch_dir))
}
