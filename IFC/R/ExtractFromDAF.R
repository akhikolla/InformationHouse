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

#' @title DAF File Reader
#' @description
#' Extracts data from DAF Files.
#' @param fileName path to file.
#' @param extract_features whether to extract features (and graphs, pops and regions) from file. Default is TRUE.
#' @param extract_images whether to extract images information from file. Default is TRUE.
#' @param extract_offsets whether to extract IFDs offsets from corresponding. Default is TRUE.\cr
#' See \code{\link{getOffsets}} for further details.
#' @param extract_stats whether to extract population statistics. Default is TRUE.
#' @param endianness The endian-ness ("big" or "little") of the target system for the file. Default is .Platform$endian.\cr
#' Endianness describes the bytes order of data stored within the files. This parameter may not be modified.
#' @param pnt_in_poly_algorithm algorithm used to determine if object belongs to a polygon region or not. Default is 1.\cr
#' Note that for the moment only 1(Trigonometry) is available.
#' @param pnt_in_poly_epsilon epsilon to determine if object belongs to a polygon region or not. It only applies when algorithm is 1. Default is 1e-12.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... Other arguments to be passed.
#' @details When extract_features is TRUE it allows eatures, graphs, pops, regions to be extracted.\cr
#' If extract_features is TRUE, extract_stats will be automatically forced to TRUE.\cr
#' If extract_stats is TRUE, extract_features will be automatically forced to TRUE.\cr
#' If extract_offsets is TRUE, extract_images will be automatically forced to TRUE.\cr
#' If extract_images is TRUE, information about images will be extracted.
#' @source For pnt_in_poly_algorithm, Trigonometry, is an adaptation of Jeremy VanDerWal's code \url{https://github.com/jjvanderwal/SDMTools}
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return A named list of class `IFC_data`, whose members are:\cr
#' -description, a list of descriptive information,\cr
#' -fileName, path of fileName input,\cr
#' -fileName_image, path of .cif image fileName is refering to,\cr
#' -features, a data.frame of features,\cr
#' -features_def, a describing how features are defined,\cr
#' -graphs, a list of graphical elements found,\cr
#' -pops, a list describing populations found,\cr
#' -regions, a list describing how regions are defined,\cr
#' -images, a data.frame describing information about images,\cr
#' -offsets, an integer vector of images and masks IFDs offsets,\cr
#' -stats, a data.frame describing populations count and percentage to parent and total population,\cr
#' -checksum, checksum of .cif image fileName is refering to computed from images values found in current daf.
#' @export
ExtractFromDAF <- function(fileName, extract_features = TRUE, extract_images = TRUE, extract_offsets = TRUE, extract_stats = TRUE,
                           endianness = .Platform$endian, pnt_in_poly_algorithm = 1, pnt_in_poly_epsilon = 1e-12, display_progress = TRUE, ...) {
  dots=list(...)
  if(missing(fileName)) stop("'fileName' can't be missing")
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = "daf")
  if(!file.exists(fileName)) stop(paste("can't find",fileName,sep=" "))
  title_progress = basename(fileName)
  extract_features = as.logical(extract_features); assert(extract_features, len = 1, alw = c(TRUE, FALSE))
  extract_images = as.logical(extract_images); assert(extract_images, len = 1, alw = c(TRUE, FALSE))
  extract_offsets = as.logical(extract_offsets); assert(extract_offsets, len = 1, alw = c(TRUE, FALSE))
  extract_stats = as.logical(extract_stats); assert(extract_stats, len = 1, alw = c(TRUE, FALSE))
  assert(endianness, len = 1, alw= c("big", "little"))
  pnt_in_poly_algorithm = as.integer(pnt_in_poly_algorithm); assert(pnt_in_poly_algorithm, len = 1, alw = 1)
  pnt_in_poly_epsilon = as.numeric(pnt_in_poly_epsilon); pnt_in_poly_epsilon = pnt_in_poly_epsilon[pnt_in_poly_epsilon>0]; pnt_in_poly_epsilon = pnt_in_poly_epsilon[is.finite(pnt_in_poly_epsilon)]
  assert(pnt_in_poly_epsilon, len = 1, typ = "numeric")
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  fileName = normalizePath(fileName, winslash = "/", mustWork = FALSE)
  toskip=cpp_scanFirst(fname = fileName, target = "</Assay>", start = 0, end = 0)
  if(toskip == 0) stop(paste0(fileName, "\ndoes not seem to be well formatted: </Assay> not found")) 
  toskip = toskip + nchar("</Assay>") - 1
  tmp=read_xml(readBin(con = fileName, what = "raw", n = toskip), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  checksum=NULL

  ##### extracts description
  description=list("Assay"=xml_attrs(xml_find_all(tmp, "//Assay")),
                   "FCS"=xml_attrs(xml_find_all(tmp, "//FCS")),
                   "SOD"=xml_attrs(xml_find_all(tmp, "//SOD")),
                   "Images"=xml_attrs(xml_find_all(tmp, "//image")),
                   "masks"=xml_attrs(xml_find_all(tmp, "//mask")))
  description=lapply(description, FUN=function(x) {as.data.frame(do.call(what="rbind", x), stringsAsFactors=FALSE)})
  if(length(description$FCS)==0) {
    description$ID = description$SOD
    is_fcs = FALSE
  } else {
    description$ID = description$FCS
    is_fcs = TRUE
  }
  obj_count = as.integer(description$ID$objcount)
  class(description$masks) <- c(class(description$masks), "IFC_masks")
  description$ID$objcount = obj_count
  chan_number = as.integer(xml_attr(xml_find_first(tmp, "//ChannelPresets"), attr = "count"))
  
  if(!is_fcs) {
    description$Images = description$Images[order(description$Images$physicalChannel),]
    checksum = checksumDAF(fileName)
    # chan_number = nrow(description$Images) # when from daf only available channels are imported

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
    
    fileName_image = paste(dirname(fileName),description$ID$file,sep="/")
    if(file.exists(fileName_image) && (checksum == checksumXIF(fileName_image))) {
      fileName_image = normalizePath(fileName_image, winslash = "/")
    } else {
      fileName_image = description$ID$file
    }
  } else {
    description$Images = as.data.frame(matrix(NA, nrow=0, ncol=12, 
                                              dimnames = list(c(), c("name","color","physicalChannel","xmin","xmax","xmid","ymid","scalemin","scalemax","tokens","baseimage","function"))))
    if(extract_offsets) {
      extract_offsets = FALSE
      message("'extract_offsets' has been forced to FALSE because no offsets can be found in .daf originated from .fcs.")
    }
    fileName_image = description$ID$file
  }

  offsets = NULL
  if(extract_offsets & !extract_images) { 
    extract_images = TRUE
    message("'extract_images' has been forced to TRUE to extract offsets.")
  }
  
  toread=file(description = fileName, open = "rb")
  on.exit(close(toread), add = TRUE)
  
  ##### checks how features and images are stored
  if(extract_features | extract_images) {
    is_binary=as.logical(na.omit(xml_attr(xml_find_first(tmp, "//Assay"), attr = "binaryfeatures")))
    if(length(is_binary)==0) {is_binary=FALSE}
    if(is_fcs) is_binary=FALSE
    if(is_binary) {
      seek(toread,toskip+3)
      ##### extracts important values
      feat_version=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      feat_number=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      obj_number=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      if(obj_count != obj_number) stop("mismatch between expected object count and features values stored")
    }
  }
  images = data.frame()
  ##### extracts images
  if(extract_images) {
    if(is_binary) {
      seek(toread, toskip + feat_number*(obj_number*8 + 4) + 15)
      SO_number=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness)) # number of SO
      if(SO_number != obj_number) stop("mismatch between expected object count and images numbers stored")
      if(display_progress) {
        pb_im = newPB(session = dots$session, min = 0, max = SO_number, initial = 0, style = 3)
        tryCatch({
        images=lapply(1:SO_number, FUN=function(i_image) {
          setPB(pb_im, value = i_image, title = title_progress, label = "extracting images values (binary)")
          id = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          imgIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used, img offsets are uint32
          mskIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used, msk offsets are uint32
          spIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used ?
          w = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          l = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          fs = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          cl = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          ct = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          objCenterX = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          objCenterY = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          n_s = cpp_int32_to_uint32(readBin(toread, "integer", size = 4,  n = 1, endian = endianness))
          bgstd = readBin(toread, "double", size = 8, n = n_s)
          n_m = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          bgmean = readBin(toread, "double", size = 8, n = n_m)
          n_c = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          satcount = readBin(toread, "double", size = 8, n = n_c)
          n_p = cpp_int32_to_uint32(readBin(toread, "integer",size = 4,  n = 1, endian = endianness))
          satpercent = readBin(toread, "double", size = 8, n = n_p)
          if(any(c(n_s, n_m, n_c, n_p) != chan_number)) stop("mismatch between expected channel numbers and channel values stored")
          c("id"=id,"imgIFD"=imgIFD,"mskIFD"=mskIFD,"spIFD"=spIFD,
            "w"=w,"l"=l,"fs"=fs,
            "cl"=cl,"ct"=ct,
            "objCenterX"=objCenterX,"objCenterY"=objCenterY,
            "bgstd"=bgstd,"bgmean"=bgmean,"satcount"=satcount,"satpercent"=satpercent)
        })
        }, error = function(e) {
          stop(e$message)
        }, finally = endPB(pb_im))
      } else{
        images=lapply(1:SO_number, FUN=function(i_image) {
          id = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          imgIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used img offsets are uint32
          mskIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used msk offsets are uint32
          spIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used ?
          w = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          l = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          fs = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          cl = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          ct = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          objCenterX = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          objCenterY = readBin(toread, "double", size = 8, n = 1, endian = endianness)
          n_s = cpp_int32_to_uint32(readBin(toread, "integer", size = 4,  n = 1, endian = endianness))
          bgstd = readBin(toread, "double", size = 8, n = n_s)
          n_m = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          bgmean = readBin(toread, "double", size = 8, n = n_m)
          n_c = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          satcount = readBin(toread, "double", size = 8, n = n_c)
          n_p = cpp_int32_to_uint32(readBin(toread, "integer",size = 4,  n = 1, endian = endianness))
          satpercent = readBin(toread, "double", size = 8, n = n_p)
          if(any(c(n_s, n_m, n_c, n_p) != chan_number)) stop("mismatch between expected channel count and channel values stored")
          c("id"=id,"imgIFD"=imgIFD,"mskIFD"=mskIFD,"spIFD"=spIFD,
            "w"=w,"l"=l,"fs"=fs,
            "cl"=cl,"ct"=ct,
            "objCenterX"=objCenterX,"objCenterY"=objCenterY,
            "bgstd"=bgstd,"bgmean"=bgmean,"satcount"=satcount,"satpercent"=satpercent)
        })
      }
      images=as.data.frame(do.call(what = "rbind", args = images), stringsAsFactors = FALSE)
    } else {
      images=do.call(what = cbind, args = xml_attrs(xml_find_all(tmp, "//SO")))
      if(ncol(images) != obj_count) stop("mismatch between expected object count and images numbers stored")
      img_tmp_tomodify=c("bgstd","bgmean","satcount","satpercent")
      img_tmp_tokeep=grep(paste0(img_tmp_tomodify, collapse="|"), dimnames(images)[[1]], invert = TRUE)
      img_tmp_new=list()
      if(!is_fcs) img_tmp_new=lapply(img_tmp_tomodify, FUN=function(k) do.call("cbind", strsplit(images[k,],"|", useBytes = TRUE, fixed=TRUE)))
      tryCatch({
        images=rbind(images[img_tmp_tokeep,],do.call("rbind",img_tmp_new))
      }, error = function(e) {
        stop("mismatch between expected channel count and channel values stored")
      })
      images=apply(images,1,as.numeric)
      rm(list=ls()[grep("img_tmp_",ls())])
      images=as.data.frame(images, stringsAsFactors = FALSE)
      if(!is_fcs) {
        names(images)=c("id","imgIFD","mskIFD","spIFD","w","l","fs","cl","ct","objCenterX","objCenterY",
                        paste0("bgstd",(1:chan_number)),
                        paste0("bgmean",(1:chan_number)),
                        paste0("satcount",(1:chan_number)),
                        paste0("satpercent",(1:chan_number)))
      }
    }
    class(images) <- c(class(images), "IFC_images")
    
    ##### extracts offsets from images in DAF
    if(extract_offsets) {
      if(nrow(images) > 0) {
        offsets = as.integer(unlist(lapply(1:nrow(images), FUN=function(i) {
          c(images$imgIFD[i], images$mskIFD[i])
        }))) 
      } 
    }
  }
  if(length(offsets) != 0) {
    N = nchar(sprintf("%1.f",abs(obj_count-1)))
    names(offsets) = paste0(c("img_","msk_"), rep(sprintf(paste0("%0",N,".f"), images$id), each = 2))
    attr(offsets, "all") = offsets
    attr(offsets, "fileName_image") = fileName_image
    # if offsets are extracted from features and fileName_image can't be found, sum of 10 first images & masks offsets is used
    # /!\ since retrieving 1st offset depends on fileName_image it is not used in this sum
    # /!\ note that using offsets extracted from daf without fileName_image may lead to different results
    # TODO ask AMNIS how checksum is computed to use same alg and place description$ID$checksum instead
    attr(offsets, "checksum") = checksumDAF(fileName)
    attr(offsets, "class") = "IFC_offset"
    attr(offsets, "first") = NULL
  }
  
  features = data.frame()
  features_def = list()
  pops = list()
  plots = list()
  regions = list()
  stats = data.frame()

  if(extract_features) {
    ##### extracts features definition
    features_def=lapply(xml_attrs(xml_find_all(tmp, "//UDF")), FUN=function(x) as.list(x))
    ##### extracts features values
    if(is_binary) {
      seek(toread, toskip+15)
      if(display_progress) {
        pb_fe = newPB(session = dots$session, min = 0, max = feat_number, initial = 0, style = 3)
        tryCatch({
        features=lapply(1:feat_number, FUN=function(i_feat) {
          setPB(pb_fe, value = i_feat, title = title_progress, label = "extracting features values (binary)")
          fid=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          fv=readBin(toread, "double", size = 8, n = obj_number, endian = endianness)
          return(c(fid,fv))
        })
      }, error = function(e) {
        stop(e$message)
      }, finally = endPB(pb_fe))
      } else {
        features=lapply(1:feat_number, FUN=function(i_feat) {
          fid=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
          fv=readBin(toread, "double", size = 8, n = obj_number, endian = endianness)
          return(c(fid,fv))
        })
      }
      if(length(features)==1) {
        features=as.data.frame(matrix(features[[1]][-1], ncol=1), stringsAsFactors = FALSE)
      } else {
        features=as.data.frame(do.call(what = cbind, args = features), stringsAsFactors = FALSE)
        features_def = features_def[order(features[1,])]
        features = features[-1,order(features[1,])]
      }
    } else {
      features=xml_attr(xml_find_all(tmp, "//UDFValues"), attr = "fv")
      feat_number=length(features)
      if(display_progress) {
        pb_fe = newPB(session = dots$session, min = 0, max = feat_number, initial = 0, style = 3)
        tryCatch({
        features=lapply(1:feat_number,FUN=function(i_feat) {
          setPB(pb_fe, value = i_feat, title = title_progress, label = "extracting features values (xml)")
          val = suppressWarnings(as.numeric(strsplit(features[i_feat],"|", useBytes = TRUE, fixed=TRUE)[[1]]))
          val[is.na(val)] <- NaN
          val
        })
      }, error = function(e) {
        stop(e$message)
      }, finally = endPB(pb_fe))
      } else {
        features=lapply(1:feat_number,FUN=function(i_feat) {
          val = suppressWarnings(as.numeric(strsplit(features[i_feat],"|", useBytes = TRUE, fixed=TRUE)[[1]]))
          val[is.na(val)] <- NaN
          val
        }) 
      }
      fid = as.numeric(xml_attr(xml_find_all(tmp, "//UDFValues"), attr = "fid"))
      features_def = features_def[order(fid)]
      if(length(features)==1) {
        features=as.data.frame(matrix(features[[1]][-1], ncol=1), stringsAsFactors = FALSE)
      } else {
        features=as.data.frame(do.call(what="cbind", args = features[order(fid)]), stringsAsFactors = FALSE)
      }
      rm(list=c("fid"))
    }
    features_names = sapply(features_def, FUN=function(x) x$name)
    def_def = sapply(features_def, FUN=function(x) x$def)
    names(features_def) = features_names
    names(features) = features_names
    
    modify_feat = FALSE
    tmp_logical = features_names == "Object Number"
    if(any(tmp_logical)) { # Object Number is found
      if(any(def_def[tmp_logical] != "Object Number")) { # Object Number is not well defined # worst case because we are forced to remove it
        # copy bad Object Number features
        names(features)[tmp_logical] <- paste0(names(features)[tmp_logical], "_copied_by_IFC") # there should be only one
        features_def[tmp_logical] <- lapply(features_def[tmp_logical], FUN = function(x) {
          x$name = paste0(x$name, "_copied_by_IFC")
          return(x)
        })
        modify_feat = TRUE
        # add new Object Number feature
        features_names = c(features_names, "Object Number")
        features$`Object Number` = 0:(nrow(features)-1)
        features_def = c(features_def, "Object Number" = list(name = "Object Number", type = "single", userfeaturetype = "No Parameters", def = "Object Number"))
      } # otherwise it is ok i.e. Object Number exists and is well def
    } else { # Object Number is not found
      if(any(def_def == "Object Number")) { # Object Number is defined but not named Object Number
        # copy it
        features_names = c(features_names, "Object Number")
        features$`Object Number` = features[, which(def_def == "Object Number")[1]] # there could be several
        features_def = c(features_def, "Object Number" = list(name = "Object Number", type = "single", userfeaturetype = "No Parameters", def = "Object Number"))
      } else {
        # create it
        features_names = c(features_names, "Object Number")
        features$`Object Number` = 0:(nrow(features)-1)
        features_def = c(features_def, "Object Number" = list(name = "Object Number", type = "single", userfeaturetype = "No Parameters", def = "Object Number"))
      }
    }
    if(any(duplicated(features$`Object Number`))) {
      features$`Object Number` = 0:(nrow(features)-1)
      warning(paste0("found duplicated objects when reading file: ", fileName))
    }
    rownames(features) = 0:(nrow(features)-1)
    class(features) <- c(class(features),"IFC_features")
    class(features_def) <- c(class(features_def),"IFC_features_def")
    ##### extracts graphs information
    plots=lapply(xml_attrs(xml_find_all(tmp, "//Graph")), FUN=function(x) as.list(x))
    if(length(plots)!=0) {
      # plots=mapply(plots, FUN=c, SIMPLIFY = FALSE)
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
      if(modify_feat) {
        plots = lapply(plots, FUN = function(g) {
          if((length(g$f1)!=0) && (g$f1=="Object Number")) g$f1 = "Object Number_copied_by_IFC"
          if((length(g$xlabel)!=0) && (g$xlabel=="Object Number")) g$xlabel = "Object Number_copied_by_IFC"
          if((length(g$f2)!=0) && (g$f2=="Object Number")) g$f2 = "Object Number_copied_by_IFC"
          if((length(g$ylabel)!=0) && (g$ylabel=="Object Number")) g$ylabel = "Object Number_copied_by_IFC"
          return(g)
        })
      }
    }
    
    ##### TODO, add something for ChannelImage, ObjectFeatureControl, StatisticsControl
    
    ##### extracts regions information
    regions=lapply(xml_attrs(xml_find_all(tmp, "//Region")), FUN=function(x) as.list(x))
    if(length(regions) != 0) {
      names(regions)=lapply(regions, FUN=function(x) x$label)
      # regions=mapply(regions, FUN=c, SIMPLIFY = FALSE)
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
      # pops=mapply(pops, fromIDEAS=TRUE, FUN=c, SIMPLIFY = FALSE)
      # pops=mapply(pops, FUN=c, SIMPLIFY = FALSE)
      pops_=lapply(pops, FUN=function(i_pop) {
        pat=paste0("//Pop[@name='",i_pop$name,"']//ob")
        list(obj=as.numeric(unlist(xml_attrs(xml_find_all(tmp, pat)))))
      })
      pops=mapply(FUN = append, pops, pops_, SIMPLIFY = FALSE)
      rm(pops_)
    }
    class(pops) <- "IFC_pops"
    
    operators = c("And","Or","Not","(",")")
    l = length(pops)
    if(l>0) {
      if(modify_feat) {
        pops = lapply(pops, FUN = function(p) {
          if((length(p$fx)!=0) && (p$fx=="Object Number")) p$fx = "Object Number_copied"
          if((length(p$fy)!=0) && (p$fy=="Object Number")) p$fy = "Object Number_copied"
          return(p)
        })
      }
      ###### scrambles pops (for testing) 
      # pops = pops[sample.int(length(pops))]; str(names(pops))
      
      ##### extracts populations dependencies/affiliations.
      ##### reorders pops
      ##### determines which object belongs to each population and changes styles and colors
      pops = popsCompute(pops = pops,
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
          c("type" = type, "parent" = base, "count" = count, "perc_parent" = count/parent*100, "perc_tot" = count/obj_count*100)
        })))
        stats[,3] = as.numeric(stats[,3])
        stats[,4] = as.numeric(stats[,4])
        stats[,5] = as.numeric(stats[,5])
      }
    }
    
    #####  retrieve name(s) of graphical population created by region applied in graph
    if(length(plots) > 0) {
      plots = lapply(plots, FUN = function(g) {
        if(length(g$GraphRegion) != 0) {
          g$GraphRegion = sapply(g$GraphRegion, FUN = function(r) {
            foo = sapply(pops,
                   FUN = function(p) {
                     bar = (p$type == "G") && 
                       (p$region == r$name) && 
                       (p$base %in% unique(unlist(lapply(g$BasePop, FUN = function(b) b$name)))) &&
                       (g$f1 == p$fx)
                     if(regions[[r$name]]$type != "line") bar = bar && (g$f2 == p$fy)
                     return(bar)
                   })
            return(list(c(r, list(def = names(which(foo))))))
          })
        }
        return(g)
      })
    }
    class(plots) <- "IFC_graphs"
  }
  ans = list("description"=description, "fileName"=fileName, "fileName_image"=fileName_image, "features"=features, "features_def"=features_def, "graphs"=plots, "pops"=pops, "regions"=regions, "images"=images, "offsets"=offsets, "stats"=stats, "checksum" = checksum)
  attr(ans, "class") <- c("IFC_data") #,"analysis")
  return(ans)
}
