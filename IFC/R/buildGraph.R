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

#' @title IFC Graph Coercion
#' @description
#' Helper to build a list to allow graph export.
#' @param type Graph's type. Either "histogram", "scatter" or "density". Default is "density".
#' @param xlocation Integer. Graph's x location. Default is 0.
#' @param ylocation Integer. Graph's x location. Default is 0.
#' @param f1 Character. Graph x axis parameter. Default is "Object Number".
#' @param f2 Character. Graph y axis parameter. Default is "Object Number". Only used when 'type' is not "histogram".
#' @param scaletype Integer. Graph scale. Either 0 (auto), 1 (manual). Default is 1.
#' @param xmin Double. Graph's xmin. Default -1.
#' @param xmax Double. Graph's xmax. Default 1.
#' @param ymin Double. Graph's xmin. Default 0.
#' @param ymax Double. Graph's xmax. Default 1.
#' @param title Character. Graph title label. Default will use names of BasePop collapse with ', '.
#' @param xlabel Character. Graph x axis label.
#' @param ylabel Character. Graph y axis label.
#' @param axislabelsfontsize Integer. Axis label font size. Default is 10. Allowed are: 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 26, 28.\cr
#' Checked but not yet implemented.
#' @param axistickmarklabelsfontsize Integer. Axis tick font size. Default is 10. Allowed are: 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 26, 28.\cr
#' Checked but not yet implemented.
#' @param graphtitlefontsize Integer. Axis title font size. Default is 12. Allowed are: 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 26, 28.\cr
#' Checked but not yet implemented.
#' @param regionlabelsfontsize Integer. Axis region font size. Default is 10. Allowed are: 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 26, 28.\cr
#' Checked but not yet implemented.
#' @param bincount Integer. Histogram bin count. Default is 0. Allowed are: 0, 8, 16, 32, 64, 128, 256, 512, 1024.
#' @param freq Character. Histogram with frequency normalization of not. Default is "T", allowed are "T" or "F".
#' @param histogramsmoothingfactor Integer. Histogram smoothing factor. Allowed are [0-20]. Only partly implemented, default is 0 for no smoothing other values will produce same smoothing. 
#' @param xlogrange determines hyper parameter of smoothLinLog transformation for x-axis. Default is "P" for no transformation.
#' @param ylogrange determines hyper parameter of smoothLinLog transformation for y-axis. Default is "P" for no transformation.
#' @param stats Character. Either "true" or "false" to display stats. Default is "false".
#' @param xsize Integer. Graph's x size. Default is 320 for small. Regular are: 320 (small), 480 (medium), 640 (big).
#' Checked but not yet implemented.
#' @param ysize Integer. Graph's y size. Default is 'ysize' + 'splitterdistance' when 'stats' is set to "true".
#' Checked but not yet implemented.
#' @param splitterdistance Integer. Default is 120.
#' Checked but not yet implemented.
#' @param xstats Character. x stats to be computed. Default is 'Count|\%Gated|Mean'.
#' It has to be a filled with the concatenation of 'Count', '\%Total', '\%Gated', 
#' '\%Plotted', 'Objects/mL', 'Mean', 'Median', 'Std. Dev.', 'MAD', 'CV',
#' 'Minimum', 'Maximum', 'Geo. Mean', 'Mode', 'Variance' and /or 'NaN', collapse with '|'.
#' Checked but not yet implemented.
#' @param ystats Character. y stats to be computed. Should be identical to 'xstats'. Default is xstats.
#' Checked but not yet implemented.
#' @param order Character. Order to display populations. 
#' When 'type' is "density" it will be BasePop[[1]]$name.
#' When 'type' is "histogram" or "density" 'ShownPop' are not allowed
#' Otherwise, it will use each of 'GraphRegion', 'BasePop' and 'ShownPop' names, collapse with '|'.
#' @param xstatsorder Character. Order of stat rows.
#' It will use each of 'GraphRegion' names & each of 'BasePop' names, reverted and collapse with '|'.
#' @param Legend Default is list(list(onoff='false',x='0',y='0',witdh='96',height='128')).
#' Not yet implemented.
#' @param BasePop Default is list(list()). See details.
#' @param GraphRegion Default is list(list()). Only allowed member are sub-list(s) with only one character component named 'name'.
#' @param ShownPop  Default is list(list()). Only allowed member are sub-list(s) with only one character component named 'name'.
#' @details Many parameters are not used or are only partly implemented, but most are checked in order to be compatible for further export.\cr
#' For 'BasePop', if left as is "All" will be used as default.\cr
#' This parameter will be built / checked according to 'type' argument.\cr
#' 'BasePop' has to be a list of list(s) and each sub-list should can contain several elements, but only "name" is mandatory.\cr
#' The sublist mebers ar:\cr
#' -"name", "linestyle", "fill",\cr
#' and only when 'type' is "density"\cr
#' -"densitybincount", "densitymin", "densitymax",\cr
#' -"densitycolors", "densitycolorslightmode", "densitycolorsdarkmode".\cr
#' Each sub-list will be created automatically with the following default values (except if explicitly provided):\cr
#' -linestyle='Solid',\cr
#' -fill='true',\cr
#' -densitybincount='128',densitymin='0',densitymax='0',\cr
#' -densitycolors='-16776961|-13447886|-256|-23296|-65536|',\cr
#' -densitycolorslightmode='-16776961|-13447886|-256|-23296|-65536|',\cr
#' -densitycolorsdarkmode='-16776961|-13447886|-256|-23296|-65536|'\cr
#' Note that when 'type' is "density", 'BasePop' should be of length one.\cr
#' and fill will be overwritten to 'true'.
#' @param ... Other arguments to be passed.
#' @return a list containing all graph information.
#' @export
buildGraph <- function(type=c("histogram","scatter","density")[3], xlocation=0, ylocation=0,
                      f1="Object Number", f2="Object Number", scaletype=1, 
                      xmin=-1, xmax=1, ymin=0, ymax=1,
                      title=paste0(unlist(lapply(BasePop, FUN=function(x) x$name)),collapse=", "),
                      xlabel=f1, ylabel=f2, 
                      axislabelsfontsize=10, axistickmarklabelsfontsize=10, graphtitlefontsize=12, regionlabelsfontsize=10,
                      bincount=0, freq=c("T","F")[1], histogramsmoothingfactor=0,
                      xlogrange="P", ylogrange="P", splitterdistance=120,
                      stats=c("true","false")[2], xsize=c(320,480,640)[1], ysize=xsize+ifelse(stats=="true",splitterdistance,0),
                      xstats="Count|%Gated|Mean", ystats=xstats,
                      order, xstatsorder, Legend,
                      BasePop=list(list()),
                      GraphRegion=list(list()),
                      ShownPop=list(list()), ...) {
  dots = list(...)
  assert(type, len=1, alw=c("histogram","scatter","density"))
  xlocation = na.omit(as.integer(xlocation));# xlocation = xlocation[xlocation>=0]
  assert(xlocation, typ="integer", len=1)
  ylocation = na.omit(as.integer(ylocation));# ylocation = ylocation[ylocation>=0] 
  assert(ylocation, typ="integer", len=1)
  assert(f1, len=1, typ="character")
  assert(stats, len=1, alw=c("true","false"))
  if(missing(Legend)) Legend=list(list(onoff="false",x="0",y="0",witdh="96",height="128"))
  ###### Removed since xsize and ysize can be freely defined
  # assert(xsize, len=1, alw=c(320,480,640))
  # assert(splitterdistance, len=1, alw=120)
  # assert(ysize, len=1, alw=xsize+ifelse(stats=="true",splitterdistance,0))
  ######
  xsize=na.omit(as.integer(xsize)); xsize = xsize[xsize>=0]
  assert(xsize, len=1, typ="integer")
  ysize=na.omit(as.integer(ysize)); ysize = ysize[ysize>=0]
  assert(ysize, len=1, typ="integer")
  splitterdistance=na.omit(as.integer(splitterdistance)); splitterdistance = splitterdistance[splitterdistance>=0]
  assert(splitterdistance, len=1, typ="integer")
  assert(xlabel, len=1, typ="character")
  assert(freq, len=1, alw=c("T","F"))
  font_size_avl = as.integer(c(8:11,(6:14)*2))
  
  if(length(axislabelsfontsize)==0) axislabelsfontsize = font_size_avl[1]; axislabelsfontsize = as.integer(axislabelsfontsize); assert(axislabelsfontsize, len=1, alw=font_size_avl)
  if(length(axistickmarklabelsfontsize)==0) axistickmarklabelsfontsize = font_size_avl[1]; axistickmarklabelsfontsize = as.integer(axistickmarklabelsfontsize); assert(axistickmarklabelsfontsize, len=1, alw=font_size_avl)
  if(length(graphtitlefontsize)==0) graphtitlefontsize = font_size_avl[1]; graphtitlefontsize = as.integer(graphtitlefontsize); assert(graphtitlefontsize, len=1, alw=font_size_avl)
  if(length(regionlabelsfontsize)==0) regionlabelsfontsize = font_size_avl[1]; regionlabelsfontsize = as.integer(regionlabelsfontsize); assert(regionlabelsfontsize, len=1, alw=font_size_avl)
  if(length(histogramsmoothingfactor)==0) histogramsmoothingfactor = as.integer(0); histogramsmoothingfactor = as.integer(histogramsmoothingfactor); assert(histogramsmoothingfactor, len=1, alw=as.integer(0:20))
  
  bincount = as.integer(bincount); assert(bincount, len=1, alw=as.integer(c(0,2^(3:10))))
  xlogrange = as.character(xlogrange)
  ylogrange = as.character(ylogrange)
  if(xlogrange != "P") {
    tmp = as.numeric(xlogrange)
    if(is.na(tmp)) stop("'xlogrange' should be coercible to positive numeric")
    if(tmp<0) stop("'xlogrange' should be coercible to positive numeric")
    xlogrange = as.character(tmp)
  }
  if(ylogrange != "P") {
    tmp = as.numeric(ylogrange)
    if(is.na(tmp)) stop("'ylogrange' should be coercible to positive numeric")
    if(tmp<0) stop("'ylogrange' should be coercible to positive numeric")
    ylogrange = as.character(tmp)
  }
  stats_alw = c("Count","%Total","%Gated","%Plotted","Objects/mL","Mean","Median","Std. Dev.","MAD","CV","Minimum","Maximum","Geo. Mean","Mode","Variance","NaN")
  assert(xstats, len=1, typ="character"); stopifnot(strsplit(xstats, split="|", fixed=TRUE)[[1]] %in% stats_alw)
  assert(ystats, len=1, typ="character"); stopifnot(xstats == ystats, strsplit(ystats, split="|", fixed=TRUE)[[1]] %in% stats_alw)
  BasePop_name_alw = c("name", "linestyle", "fill", "densitybincount", "densitymin", "densitymax", "densitycolors", "densitycolorslightmode", "densitycolorsdarkmode")
  
  # starts building args
  args=list(type=type, xlocation=xlocation, ylocation=ylocation, f1=f1)
  # clean-up BasePop, GraphRegion, ShownPop
  BasePop = lapply(BasePop, FUN=function(x) {  # removes BasePop where name is missing
    if("name"%in%names(x)) {
      return(x)
    }
    return(NULL)
  })
  if(length(BasePop)==0) BasePop = list(list("name"="All"))
  GraphRegion = lapply(GraphRegion, FUN=function(x) {  # removes GraphRegion where name is missing
    if("name"%in%names(x)) {
      return(x)
    }
    return(NULL)
  })
  ShownPop = lapply(ShownPop, FUN=function(x) {  # removes ShownPop where name is missing
    if("name"%in%names(x)) {
      return(x["name"])
    }
    return(NULL)
  })
  BasePop_style_alw = c("Dash","DashDot", "DashDotDot", "Dot", "Solid")
  BasePop_fill_alw = c("true", "false")
  if(type=="histogram") {
    BasePop_name_alw = BasePop_name_alw[1:3]
    if(freq=="T") ylabel = "Normalized Frequency"
    if(freq=="F") ylabel = "Frequency"
  } else {
    BasePop_style_alw = "Solid" # forced for non histogram
    BasePop_fill_alw = "true" # forced for non histogram
    if(type=="density") {
      if(length(BasePop)>1) {
        BasePop = BasePop[1]
        order = paste0(rep(BasePop[[1]]$name, 5), collapse = "|")
        warning("Density graphs can only display one BasePop population", call.=FALSE, immediate.=TRUE)
        # stop("Density graphs can only display one BasePop population", call.=FALSE)
      }
      if(length(ShownPop)!=0) if(length(ShownPop[[1]])>0) {
        ShownPop = list(list())
        order = paste0(rep(BasePop[[1]]$name, 5), collapse = "|")
        warning("Density graphs can't display ShownPop population", call.=FALSE, immediate.=TRUE)
        # stop("Density graphs can't display ShownPop population", call.=FALSE)
      }
    } else {
      if(type=="histogram") if(length(ShownPop)!=0) if(length(ShownPop[[1]])>0) stop("Histogram graphs can't display ShownPop population", call.=FALSE)
      BasePop_name_alw = BasePop_name_alw[1:3]
    }
    assert(f2, len=1, typ="character")
    assert(ylabel, len=1, typ="character")
    args = c(args, list(f2=f2)) # adds f2 to args for non histogram
  }
  BasePop_default = list(name="All", linestyle="Solid", fill="true",
                         densitybincount="128", densitymin="0", densitymax="0",
                         densitycolors="-16776961|-13447886|-256|-23296|-65536|",
                         densitycolorslightmode="-16776961|-13447886|-256|-23296|-65536|",
                         densitycolorsdarkmode="-16776961|-13447886|-256|-23296|-65536|")
  if(type!="density") BasePop_default = BasePop_default[1:3]
  BasePop = lapply(BasePop, FUN=function(x) {
    tmp = BasePop_name_alw %in% names(x)
    BasePop_default[tmp] = x[BasePop_name_alw[tmp]]
    if(type=="density") {
      # checks that densitybincount is ok
      densitybincount = na.omit(as.integer(BasePop_default$densitybincount)); assert(densitybincount, len=1, alw=c(128,256,512,1024))
      # checks that densitymin is ok
      densitymin = na.omit(as.integer(BasePop_default$densitymin)); assert(densitymin, len=1, alw = 0:10)
      # checks that densitymax is ok
      densitymax = na.omit(as.integer(BasePop_default$densitymax)); assert(densitymax, len=1, alw = 0:10)
      # checks that density colors are correct
      lapply(BasePop_name_alw[7:9], FUN=function(i) {
        denscols = colConv(BasePop_default[[i]])
        assert(denscols, len=5, typ="character")
      })
      # overwrites fill
      BasePop_default$fill = "true"
    }
    # checks that linestyle is ok
    assert(BasePop_default$linestyle, len=1, alw=BasePop_style_alw)
    # checks that fill is ok
    assert(BasePop_default$fill, len=1, alw=BasePop_fill_alw)
    return(BasePop_default)
  })
  # defines default order and xstatsorder
  b_names = unlist(lapply(BasePop, FUN=function(x) x$name))
  g_names = unlist(lapply(GraphRegion, FUN=function(x) x$def))
  s_names = unlist(lapply(ShownPop, FUN=function(x) x$name))
  
  # xstatsorder_tmp = gsub(" & All","", unlist(lapply(b_names, FUN=function(n) {
  #   if(length(g_names)!=0) {
  #     return(c(n,paste(g_names, n, sep=" & ")))
  #   }
  #   return(n)
  # })))
  xstatsorder_tmp = c(g_names, b_names)
  
  if(type=="histogram") order_tmp = b_names
  if(type=="density") order_tmp = rep(b_names,5)
  if(type=="scatter") {
    # order_tmp = gsub(" & All","", unlist(lapply(rev(g_names), FUN=function(n) {
    #   paste(n, rev(b_names), sep=" & ")
    # })))
    # order_tmp = c(s_names, order_tmp, b_names)
    order_tmp = c(s_names, g_names, b_names)
  }
  # checks order is possible
  # note that xstatsorder is not deeply checked ... TODO ???
  if(missing(order)) {
    order = paste0(order_tmp, collapse = "|")
  } else {
    assert(order, len=1, typ="character")
  }
  
  operators = c("And","Or","Not","(",")")
  displayed_n = splitn(definition = order, all_names = c(b_names, g_names, s_names, "Selected Bin"), operators = operators)
  displayed_n = setdiff(displayed_n, "Selected Bin")
  tmp = displayed_n %in% c(order_tmp)
  if(!all(tmp)) stop(paste0("trying to display a population not found in supplied ShownPop, BasePop, GraphRegion names: ",  paste0(displayed_n[!tmp], collapse=", ")))
  
  if(missing(xstatsorder)) {
    xstatsorder = paste0(xstatsorder_tmp, collapse = "|")
  } else {
    assert(xstatsorder, len=1, typ="character")
  }

  args = c(args, list(scaletype=scaletype, 
                      xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                      title=title, xlabel=xlabel, ylabel=ylabel, 
                      axislabelsfontsize=axislabelsfontsize, axistickmarklabelsfontsize=axistickmarklabelsfontsize, graphtitlefontsize=graphtitlefontsize,
                      regionlabelsfontsize=regionlabelsfontsize, bincount=bincount,
                      freq=freq, histogramsmoothingfactor=histogramsmoothingfactor, xlogrange=xlogrange, ylogrange=ylogrange, 
                      stats=stats, xsize=xsize, ysize=ysize, splitterdistance=splitterdistance,
                      xstats=xstats, ystats=ystats, order=order, xstatsorder=xstatsorder,
                      Legend=Legend, BasePop=BasePop, GraphRegion=GraphRegion, ShownPop=ShownPop))
  return(args)
}
