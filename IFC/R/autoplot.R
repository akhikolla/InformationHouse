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

#' @title Automatic Parameters Detection for IFC Graphs
#' @description
#' Function intended to generate IFC graphs with minimal inputs from users.\cr
#' It is essentially based on automatic detection of graphical parameters thanks to 'shown_pops' argument.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param shown_pops one or several populations present in 'obj'. Default is NULL.\cr
#' If provided, \code{\link{autoplot}} will try to display these populations. See details when not provided.\cr
#' \cr
#' \code{\link{autoplot}} will try to determine x and y and their transformations based on 'shown_pops' parameter.
#' If all populations provided in 'shown_pops' are siblings, region(s) from which 'shown_pops' were defined 
#' will be displayed.\cr
#' In case 'shown_pops' are not siblings, they will be treated as populations and a graph will be generating with an overlay of these populations.
#' Order of this overlay is given by order of 'shown_pops'.\cr
#' Finally, changing any of the following arguments (x, x_trans, y, y_trans, type)
#' to something else than the one detected from 'shown_pops' will prevent from displaying region(s) and
#' 'shown_pops' populations will be displayed as overlay.\cr
#' However, please consider that if original type is 'histogram' changing x_trans transformation will have no impact on this.
#' @param subset a population present in 'obj'. Default is NULL.
#' Background population that will be used to generate graph.
#' This argument will not be used when graph is an histogram.
#' If this argument is filled with a different population than what can be determined thanks to 'shown_pops',
#' Then 'shown_pops' will be treated as overlay.
#' However, 'shown_pops' argument can still be used to determine x, y axis and their transformation
#' @param x feature for x-axis. Default is NULL.
#' When empty, \code{\link{autoplot}} will try to determine if automatically from 'shown_pops' argument.
#' If provided, x feature has to be a name from 'obj' features.
#' Note that providing x feature :
#' - takes precedence on automatic x-axis detection.
#' - will reset x-axis transformation to "P" except if 'x_trans' is filled.
#' @param x_trans parameter for x-axis transformation. Default is NULL.
#' If not provided, transformation will be determined thanks to 'shown_pops'.
#' It takes precedence when provided and If provided it has to be be either 'P' or coercible to a positive numeric.
#' "P' will leave x-axis as is but a positive numeric will be passed has hyper argument of \code{\link{smoothLinLog}} to transform x-axis.
#' @param y feature for y-axis. Default is NULL.
#' When empty, \code{\link{autoplot}} will try to determine it automatically from 'shown_pops' argument.
#' If provided, y feature has to be a name from obj features.
#' Note that providing y feature 
#' - takes precedence on automatic y-axis detection.
#' - will reset y-axis transformation to "P" except if 'y_trans' is filled.
#' @param y_trans parameter for y-axis transformation. Default is NULL.
#' If not provided, transformation will be determined thanks to 'shown_pops'.
#' It takes precedence when provided and has to be be either 'P' or coercible to a positive numeric.
#' "P' will leave y-axis as is but a positive numeric will be passed has hyper argument of \code{\link{smoothLinLog}} to transform y-axis.
#' Note that it is irrelevant for "histogram".
#' @param type type of plot. Default is NULL to allow \code{\link{autoplot}} to detemine 'type' automatically.
#' If provided it has to be either "histogram", "scatter", "density".
#' Note that when "histogram" is choosen, 'subset' parameter will not be used.
#' Note that "density" will be possible only when 'subset' will be automatically determined or filled with only one population.
#' Note that when \code{\link{autoplot}} has determined, thanks to 'shown_pops' that original plot is an "histogram", 
#' "Object Number" will be used as y-axis by default when 'type' is forced to "scatter" or "density".
#' @param smoothingfactor when type of graph is "histogram", whether to smooth it or not. Default is NULL. Should be an integer [0:20]
#' Note that 0 means no smoothing and other values will produce smoothing
#' @param normalize when type of graph is "histogram", whether to normalize it or not. Default is NULL. Should be a logical.
#' @param bin number of bins when graph's type is "histogram" / number of equally spaced grid points for density.
#' Default is missing to allow \code{\link{autoplot}} to determine it by itself.
#' @param viewport Either "ideas", "data" or "max" defining limits used for the graph. Default is "ideas".\cr
#' -"ideas" will use same limits as the one defined in ideas.\cr
#' -"data" will use data to define limits.\cr
#' -"max" will use data and regions drawn to define limits.
#' @param precision when graphs is a 2D scatter with population overlay, this argument controls amount of information displayed. Default is "light".\cr
#' -"light", the default, will only display points of same coordinates that are amoung the other layers.\cr
#' -"full" will display all the layers.
#' @param color_mode Whether to extract colors from obj in white or black mode. Default is "white".
#' @param draw whether to draw plot. Default is TRUE.
# #' @param sub_col,sub_pch,sub_alpha,sub_lwd,sub_lty graphical parameters for subset. Default is NULL, allowing these parameters to be determined automatically.
# #' Providing these parameters takes precedence on automatic detection.
# #' If common graphical parameters are filled subset graphical parameters will be fed with them.
# #' Otherwise, 'obj' styles will be used.
# #' @param col,pch,alpha,lwd,lty common graphical parameters for 'shown_pops'. Default is NULL, allowing these parameters to be determined automatically.
# #' Providing these parameters takes precedence on automatic detection.
# #' Provided values will be repeated as necessary.
# #' If not provided, obj styles will be used.
#' @param ... Other arguments to be passed.
#' @details when 'shown_pops' are not provided, \code{\link{autoplot}} can't determine anything.\cr
#' So, if not provided default values will be used:\cr
#' -'subset' = "All"\cr
#' -'x' = "Object Number"\cr
#' -'x_trans' = "P"\cr
#' -'y' = "Object Number"\cr
#' -'y_trans' = "P"\cr
#' -'type' = "histogram"
#' @return an \pkg{lattice} trellis object
#' @export
autoplot = function(obj, shown_pops = NULL, subset = NULL, 
                    x = NULL, x_trans = NULL,
                    y = NULL, y_trans = NULL,
                    type = NULL, smoothingfactor = NULL, normalize = NULL, bin, viewport = "ideas", # plot params
                    precision=c("light","full")[1], color_mode=c("white","black")[1], draw = TRUE,  # plot params
                    ...) {
                    # sub_col=NULL, sub_pch=NULL, sub_alpha=NULL, sub_lwd=NULL, sub_lty=NULL,      # TODO: add graphical params
                    # col=NULL, alpha=NULL, pch=NULL, border=NULL, lwd=NULL, lty=NULL, fill=TRUE) {# TODO: add graphical params
  dots = list(...)

  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  # corrects variables
  subset = unique(subset)
  shown_pops = unique(shown_pops)
  trans_x = "P"
  trans_y = "P"
  if(length(x_trans)!=0) trans_x = x_trans
  if(length(y_trans)!=0) trans_y = y_trans
  if(missing(bin)) bin = NULL

  # trick for character(0), logical(0), ...
  if(!is.null(x)) if(length(x)==0) x = NULL
  if(!is.null(y)) if(length(y)==0) y = NULL
  if(!is.null(subset)) if(length(subset)==0) subset = NULL
  if(!is.null(shown_pops)) if(length(shown_pops)==0) shown_pops = NULL
  if(!is.null(x_trans)) if(length(x_trans)==0) x_trans = NULL
  if(!is.null(y_trans)) if(length(y_trans)==0) y_trans = NULL
  if(!is.null(smoothingfactor)) if(length(smoothingfactor)==0) smoothingfactor = NULL
  if(!is.null(normalize)) if(length(normalize)==0) normalize = NULL
  if(!is.null(bin)) if(length(bin)==0) bin = NULL
  if(!is.null(type)) if(length(type)==0) type = NULL
  
  # trick when user has input "" instead of null
  if(length(x)==1) if(x=="") x = NULL 
  if(length(y)==1) if(y=="") y = NULL
  if(length(subset)==1) if(subset=="") subset = NULL
  if(length(shown_pops)==1) if(shown_pops=="") shown_pops = NULL
  if(length(x_trans)==1) if(x_trans=="") x_trans = NULL
  if(length(y_trans)==1) if(y_trans=="") y_trans = NULL
  if(length(smoothingfactor)==1) if(smoothingfactor=="") smoothingfactor = NULL
  if(length(normalize)==1) if(normalize=="") normalize = NULL
  if(length(bin)==1) if(bin=="") bin = NULL
  if(length(type)==1) if(type=="") type = NULL
  
  # various checks
  if(!all(subset%in%names(obj$pops))) stop("when provided, 'subset', has to be a population from 'obj'")
  if(length(x)>1) stop("when provided 'x' should be of length 1")
  if(length(y)>1) stop("when provided 'y' should be of length 1")
  if(length(x_trans)>1) stop("when provided 'x_trans' should be of length 1")
  if(!is.null(trans_x)) if(trans_x!="P") if(!is.numeric(as.numeric(trans_x))) stop("'x_trans', should be either \"P\" or coercible to a positive numeric")
  if(length(y_trans)>1) stop("when provided 'y_trans' should be of length 1")
  if(!is.null(trans_y)) if(trans_y!="P") if(!is.numeric(as.numeric(trans_y))) stop("y_trans, should be either \"P\" or coercible to a positive numeric")
  if(!is.null(smoothingfactor)) {
    smoothingfactor = na.omit(as.integer(smoothingfactor)); assert(smoothingfactor, len=1, alw=as.integer(0:20))
  }
  if(!is.null(normalize)) {
    normalize = as.logical(normalize); assert(normalize, len=1, alw=c(TRUE,FALSE))
  }
  if(!is.null(bin)) {
    bin = na.omit(as.integer(bin)); assert(bin, len=1, typ="integer")
  }
  if(!is.null(type)) {
    type = as.character(type); assert(type, len=1, alw = c("histogram","scatter","density"))
  }
  
  # TODO add graphical param
  # validColors = paletteIFC("palette_R")
  # validPchs = c(20, 4, 3, 1, 5, 0, 2, 18, 15, 17)  
  # if(length(col)!=0) if(!all(col%in%validColors)) stop(paste0("col argument is not valid. Available are:\n",paste0(validColors, collapse=", ")))
  # if(length(sub_col)!=0) if(!all(sub_col%in%validColors)) stop(paste0("sub_col argument is not valid. Available are:\n",paste0(validColors, collapse=", ")))
  # if(length(pch)!=0) if(!all(pch%in%validPchs)) stop(paste0("pch argument is not valid. Available are:\n",paste0(validPchs, collapse=", ")))
  # if(length(sub_pch)!=0) if(!all(sub_pch%in%validPchs)) stop(paste0("sub_pch argument is not valid. Available are:\n",paste0(validPchs, collapse=", ")))
  # if(length(alpha)!=0) if(any(c(alpha>1, alpha<0))) stop("alpha argument is not valid. alpha should be in [0,1] range")
  # if(length(sub_alpha)!=0) if(any(c(sub_alpha>1, sub_alpha<0))) stop("sub_alpha argument is not valid. sub_alpha should be in [0,1] range")
  # if(length(lty)!=0) if(!all(lty%in%(0:6))) stop("lty argument is not valid. lty should be 0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash")
  # if(length(sub_lty)!=0) if(!all(sub_lty%in%(0:6))) stop("sub_lty argument is not valid. sub_lty should be 0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash")
  # 
  
  # initializes variables
  args = list("draw" = draw, "precision" = precision, "color_mode" = color_mode, "viewport" = viewport)
  if(!is.null(bin)) args = c(args, list("bin" = bin))
  foo = list()
  original = FALSE
  keep_region = FALSE
  siblings = c()
  if(length(shown_pops)!=0) siblings = popsGetSiblings(obj, shown_pops)

  # if 'shown_pops' is given and is/are siblings then
  # if 1st shown_pops is of type graphical original graph can be retrieved
  if(length(siblings) != 0) if(length(shown_pops)!=0) if(obj$pops[[shown_pops[1]]]$type == "G") {
    foo = popsRetrieveGraph(obj, shown_pops, all_siblings = FALSE)
    original = TRUE
  }

  # if x was filled overwrites what was determined from shown_pops
  if(!is.null(x)) if(!(x %in% foo$f1)) {
    if(!(x %in% names(obj$features))) stop("can't find 'x' in 'obj' features")
    foo$f1 = x
    original = FALSE
  }
  
  # if y was filled overwrites what was determined from shown_pops
  if(!is.null(y)) if(!(y %in% foo$f2)) { 
    if(!(y %in% names(obj$features))) stop("can't find 'y' in 'obj' features")
    foo$f2 = y
    original = FALSE
  }
  
  # force type
  if(!is.null(type)) {
    if(type=="histogram") {
      # force histogram when type is set to histogram
      if(any(grepl("f2", names(foo)))) {
        foo = foo[!grepl("f2", names(foo))]
        original = FALSE
      }
    } else {
      # force scatter or density when type is not histogram and no f2 was found
      if(!any(grepl("f2", names(foo)))) {
        foo$f2 = "Object Number"
        original = FALSE
      }
    }
  }

  # register type into foo
  if(length(foo$f2)==0) {
    foo$type="histogram"
  } else {
    if(length(type)!=0) {
      foo$type = type
    } else {
      foo$type = "density"
    }
  }
  
  # if x_logrange was filled overwrites what was determined from shown_pops
  if(!is.null(x_trans)) if(!(trans_x %in% foo$xlogrange)) {
    foo$xlogrange = trans_x
    if(foo$type!="histogram") original = FALSE # regions are not modified by transformation when type is histogram
  }

  # if y_logrange was filled overwrites what was determined from shown_pops
  if(!is.null(y_trans)) if(!(trans_y %in% foo$ylogrange)) {
    foo$ylogrange = trans_y
    original = FALSE
  }
  
  # if subset was filled overwrites what was determined from shown_pops
  if(!is.null(subset)) {
    if(!all(c(subset %in% sapply(foo$BasePop, FUN=function(p) p$name), sapply(foo$BasePop, FUN=function(p) p$name) %in% subset))) {
      ss = subset
      # at this step if still original, it means that axes and their transformations were retrieved 
      # thanks to shown pop, so we can still draw regions even if user changes subset
      if(original) keep_region = TRUE
      original = FALSE
    } else {
      ss = sapply(foo$BasePop, FUN=function(p) p$name)
    }
  } else {
    ss = sapply(foo$BasePop, FUN=function(p) p$name)
  }
  # at this step if no ss is found "All" population is used
  if(length(ss)==0) ss = "All"
  ss = unique(ss)
  
  # Since density plot can contain only one BasePop, changes foo$type
  if(foo$type!="histogram") if(length(ss) > 1) {
    foo$type="scatter"
    if(length(type)!=0) if(type=="density") warning("can't plot as \"density\", 'type' has been automatically set to \"scatter\"")
  }

  if(original) {
    if(foo$type=="histogram") {
      if(!is.null(normalize)) foo$freq = ifelse(normalize, "T", "F")
      if(!is.null(smoothingfactor)) foo$histogramsmoothingfactor = smoothingfactor
    }
    if(foo$type=="density") foo$BasePop = list(list("name"=ss)) # change BasePop name
    return(do.call(what = plotGraph, args = c(args, list("obj" = obj, "graph" = foo))))
  } else {
    foo = foo[!grepl("title", names(foo))]
    if(length(foo$f1)==0) foo$f1 = "Object Number"
    if(foo$type=="histogram") {
      if(length(shown_pops) == 0) {
        shown_pops = ss
      } else {
        if(length(subset)!=0) warning("autoplot has detected graph is an histogram with overlay, 'subset' argument is ignored")
      }
      foo$BasePop = lapply(shown_pops, FUN=function(p) list("name"=p, "linestyle"="Solid", "fill" = "true"))
      foo$GraphRegion = list()
      ss = shown_pops
      if(!is.null(normalize)) foo$freq = ifelse(normalize, "T", "F")
      if(!is.null(smoothingfactor)) foo$histogramsmoothingfactor = smoothingfactor
    } else {
      if(keep_region) {
        # creates new populations depending from regions
        new_pops = list()
        for(b in ss) {
          for(g in foo$GraphRegion) {
            p = list(name = paste(g$name, b, sep = " & "), type = "G", base = b, 
                     color = obj$regions[[g$name]]$color,
                     lightModeColor = obj$regions[[g$name]]$lightcolor,
                     fx = foo$f1, fy = foo$f2,
                     style = 20,
                     region = g$name,
                     obj = obj$pops[[g$name]]$obj & obj$pops[[b]]$obj,
                     names = "")
            
            new_pops = c(new_pops, list(p))
          }
        }
        names(new_pops) = lapply(new_pops, FUN=function(p) p$name)
        obj$pops = c(obj$pops, new_pops)
      } else {
        if(length(unique(c(ss, shown_pops)))>1) {
          foo$type = "scatter"
          if(length(type)!=0) if(type=="density") warning("can't plot as \"density\", 'type' has been automatically set to \"scatter\"")
        }
      }
      if(foo$type=="scatter") {
        if(keep_region) {
          foo$BasePop = lapply(ss, FUN=function(p) list("name"=p, "linestyle"="Solid", "fill" = "true"))
        } else {
          foo$BasePop = lapply(ss, FUN=function(p) list("name"=p, "linestyle"="Solid", "fill" = "true"))
          foo$ShownPop = lapply(shown_pops, FUN=function(p) list("name"=p, "linestyle"="Solid", "fill" = "true"))
          foo$GraphRegion = list() # removes grahical regions
        }
      } else {
        foo$BasePop = lapply(ss, FUN=function(p) list("name"=p, "linestyle"="Solid", "fill" = "true"))
        if(!keep_region) foo$GraphRegion = list() # removes grahical regions
        foo$ShownPop = list()
      }
    }
    if(!keep_region) {
      P = obj$pops[unique(c(ss,shown_pops))]
      SUB = as.data.frame(do.call(what=cbind, args=lapply(P, FUN=function(p) p$obj)), stringsAsFactors = FALSE)
      SUB = apply(SUB, 1, any)
      
      xran = range(obj$features[SUB, foo$f1], na.rm = TRUE)
      if(length(foo$xlogrange)==0) foo$xlogrange = trans_x
      if(foo$xlogrange == "P") {
        xran = xran + diff(xran) * c(-0.07,0.07)
      } else {
        xran = smoothLinLog(xran, hyper = as.numeric(foo$xlogrange))
        xran = xran + diff(xran) * c(-0.07,0.07)
        xran = inv_smoothLinLog(xran, hyper = as.numeric(foo$xlogrange))
      }
      if(xran[1] == xran[2]) xran = xran[1] + c(-0.07,0.07)
      foo$xmin = xran[1]
      foo$xmax = xran[2]
      if(foo$type!="histogram") {
        yran = range(obj$features[SUB, foo$f2], na.rm = TRUE)
        if(length(foo$ylogrange)==0) foo$ylogrange = trans_y
        if(foo$ylogrange == "P") {
          yran = yran + diff(yran) * c(-0.07,0.07)
        } else {
          yran = smoothLinLog(yran, hyper = as.numeric(foo$ylogrange))
          yran = yran + diff(yran) * c(-0.07,0.07)
          yran = inv_smoothLinLog(yran, hyper = as.numeric(foo$ylogrange))
        }
        if(yran[1] == yran[2]) yran = yran[1] + c(-0.07,0.07)
        foo$ymin = yran[1]
        foo$ymax = yran[2]
      }
    }
    return(do.call(what = plotGraph, args = c(args, list("obj" = obj, "graph" = foo))))
  }
}
