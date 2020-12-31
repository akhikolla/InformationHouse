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

#' @title Add Feature to IFC_data Object
#' @description
#' Adds features to an already existing `IFC_data` object.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param features a list of features to add to obj. Each element of this list will be coerced by \code{\link{buildFeature}}.
#' @details A warning will be thrown if a provided feature is already existing in obj.\cr
#' In such a case this feature will not be added to obj.\cr
#' If any input feature is not well defined and can't be created then an error will occur.
#' @param ... Other arguments to be passed.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#'   ## copy 1st feature found in daf
#'   feat_def <- daf$features_def[[1]]
#'   if(length(feat_def) != 0) {
#'     feat_def_copy <- feat_def
#'     ## modify name and value of copied features
#'     feat_def_copy$name <- "copied_feature"
#'     feat <- daf$features[, feat_def$name]
#'     feat_copy <- feat
#'     feat_copy <- feat_copy * 10
#'     ## create new object with this new feature
#'     dafnew <- data_add_features(obj = daf, features = list(c(feat_def_copy, list(val = feat_copy))))
#'   }
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return an IFC_data object with features added.
#' @export
data_add_features <- function(obj, features, ...) {
  assert(obj, cla = "IFC_data")
  
  # try to coerce features
  features = lapply(features, FUN=function(x) do.call(what=buildFeature, args=x))
  names(features) = sapply(features, FUN=function(x) x$name)
  
  # removes duplicated inputs
  tmp = duplicated(names(features))
  if(any(tmp)) {
    warning(paste0("duplicated features automatically removed: ", names(features)[tmp]), immediate. = TRUE, call. = FALSE)
    features = features[!tmp]
  }
  
  # defines available parameters
  operators_pop = c("And","Or","Not","(",")")
  operators_daf = c(operators_pop, "/","+","*","-","ABS","COS","SIN","SQR","SQRT","False","True","false","true")
  # masks_avl = c("AdaptiveErode","Component","Dilate","Erode","Fill","Inspire","Intensity","Interface","LevelSet","Morphology",
  #               "Object","Peak","Range","Skeleton","Spot","System","Threshold","Valley","Watershed")
  # masks_other = c("Combined","Dim","Middle","Bright","Tight","Dark","Thin","Thick")
  userfeatures_avl = list("Mask Only"=c("Area", "Aspect Ratio", "Length", "Width", "Height", "Angle", "Centroid X", "Centroid Y", "Circularity", "Diameter",
                                        "Elongatedness", "Major Axis", "Minor Axis", "Perimeter", "Shape Ratio", 'Spot Area Min', "Spot Distance Min",
                                        "Thickness Max", "Thickness Min"),
                          "Mask and Image"=c("Aspect Ratio Intensity", "Modulation", "Contrast", "Gradient RMS", "Intensity", "Mean Pixel", "Median Pixel",
                                             "Max Pixel", "Raw Max Pixel", "Raw Min Pixel", "Saturation Count", "Saturation Percent", "Bright Detail Intensity R3",
                                             "Bright Detail Intensity R7", "Angle Intensity", "Centroid X Intensity", "Centroid Y Intensity", "Compactness","Gradient Max", "Internalization", 
                                             "Lobe Count","Major Axis Intensity", "Max Contour Position", "Min Pixel", "Minor Axis Intensity", "Raw Intensity", "Raw Mean Pixel",
                                             "Raw Median Pixel", "Spot Intensity Max", "Spot Intensity Min", "Std Dev", "Symmetry 2", "Symmetry 3", "Symmetry 4", 
                                             "Uncompensated Intensity", "Valley X", "Valley Y"),
                          "Image Only"=paste0("Bkgd ", c("Mean", "StdDev")),
                          "No Parameters"=c("Time", "Object Number", "Raw Centroid X", "Raw Centroid Y", "Flow Speed", "Camera Line Number", "Camera Timer", "Objects per mL", "Objects per sec"),
                          "Mask, Image and Scalar"=paste0(paste0("H ", rep(c("Contrast ","Correlation ","Energy ", "Entropy ", "Homogeneity ", "Variance "), each=2)), c("Mean", "Std")),
                          "Mask and Scalar"=c("Spot Count"),
                          "Mask and Three Images"=c("Bright Detail Colocalization 3"),
                          "Similarity"=c("Bright Detail Similarity R3", "Shift X", "Shift Y", "Similarity", "XCorr"),
                          "Delta Centroid"=paste0("Delta Centroid ", c("X","Y","XY")),
                          "Two Masks and Image"=c("Intensity Concentration Ratio"),
                          "Image and Scalar"=c("Ensquared Energy", "Diameter:"))
  exported_feats = sapply(features, FUN=function(feat) {
    if(feat$name%in%names(obj$features)) {
      warning(paste0(feat$name, "\nnot exported: trying to export an already defined feature"), immediate. = TRUE, call. = FALSE)
      return(FALSE)
    }
    def = feat$def
    if(grepl("^H ", def)) def = gsub("|Granularity:|1|20", "", def, fixed = TRUE) # removes granularity from H features
    def = strsplit(def, split = "|", fixed = TRUE)[[1]]
    def = def[!(def%in%userfeatures_avl[[feat$userfeaturetype]])] # removes possible features from definition
    if(grepl("Mask",feat$userfeaturetype)) {
      def = def[!(def%in%c(obj$description$masks$name, 
                           unlist(strsplit(obj$description$masks$def[obj$description$masks$name=="MC"], "|Or|", useBytes = TRUE, fixed=TRUE))))] # removes masks from definition
    }
    if(grepl("Image",feat$userfeaturetype)) {
      def = def[!(def%in%obj$description$Images$name)] # removes channels names from definition
    }
    if(grepl("Combined",feat$userfeaturetype)) {
      def = def[!(def%in%c(names(obj$features), names(features)))] # removes features names from definition
    }
    def = def[!(def%in%operators_daf)] # removes operators from definition
    suppressWarnings({def = as.numeric(def)}) # converts remaining to numeric
    if(length(def) != 0) {
      if(!grepl("Scalar",feat$userfeaturetype) || is.na(def)) stop(paste0(feat$name, "\nbad feature definition: ", feat$def)) # if something remains which coercion to numeric produces NA, it means that features is not well defined
    }
    if(length(feat$val) != as.integer(obj$description$ID$objcount)) stop(paste0(feat$name, "\nbad feature value length, expected: ",  obj$description$ID$objcount, ", but is: ", length(feat$val))) # TODO add some lines to allow function to automatically compute feat$val when missing
    return(TRUE)
  })
  exported_feats = features[exported_feats]
  if(length(exported_feats) == 0) return(obj)
  names(exported_feats) = sapply(exported_feats, FUN=function(x) x$name)
  
  K = class(obj$features)
  obj$features = cbind(obj$features, sapply(exported_feats, FUN=function(x) x$val))
  class(obj$features) = c(setdiff(K, "IFC_features"), "IFC_features")

  K = class(obj$features_def)
  obj$features_def = c(obj$features_def, lapply(exported_feats, FUN=function(x) x[c("name", "type", "userfeaturetype", "def")]))
  class(obj$features_def) = c(setdiff(K, "IFC_features_def"), "IFC_features_def")
  return(obj)
}
