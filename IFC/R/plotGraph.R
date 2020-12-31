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

#' @title Plot and Stats Computation for IFC Graph
#' @description
#' Computes plot and stats from a IFC graph
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param graph a graph from 'obj' or a list that can be coerced by \code{\link{buildGraph}}.
#' @param draw whether to draw plot or not. Default is FALSE.
#' @param stats_print whether to print stats or not. Default is given by 'draw' argument.
#' @param color_mode whether to extract colors from 'obj' in white or black mode. Default is "white".
#' @param add_key whether to draw a "global" key under title or in the first "panel" or "both". Default is "panel".\cr
#' Accepted values are either: FALSE, "panel", "global", "both" or c("panel", "global").\cr
#' Note that it only applies when display is seen as overlaying populations.
#' @param precision when graphs is a 2D scatter with population overlay, this argument controls amount of information displayed. Default is "light".\cr
#' -"light", the default, will only display points of same coordinates that are amoung the other layers.\cr
#' -"full" will display all the layers.
#' @param trunc_labels maximum number of characters to display for labels. Default is 38.
#' @param trans transformation function for density graphs. Default is asinh.
#' @param bin number of bin used for histogram / density. Default is missing.
#' @param viewport either "ideas", "data" or "max" defining limits used for the graph. Default is "ideas".\cr
#' -"ideas" will use same limits as the one defined in ideas.\cr
#' -"data" will use data to define limits.\cr
#' -"max" will use data and regions drawn to define limits.
#' @param ... other arguments to be passed.
#' @return it invisibly returns a list whose members are:\cr
#' -plot, "trellis" object that can be displayed using plot,\cr
#' -stats, a table of satistics computed for the graph,\cr
#' -input, a list with input parameters.
#' @export
plotGraph = function(obj, graph, draw = FALSE, stats_print = draw,
                     color_mode = c("white","black")[1], add_key = "panel", precision = c("light","full")[1],
                     trunc_labels = 38, trans = asinh, bin, viewport = "ideas", ...) {
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
  assert(color_mode, len=1, alw=c("white","black"))
  color_mode=c(2,1)[c("white","black")==color_mode]
  assert(precision, len=1, alw=c("light","full"))
  if(!all(add_key%in%c("panel","global","both",FALSE))) stop("Accepted values for add_key are either: FALSE, 'panel', 'global', 'both' or c('panel', 'global')")
  trunc_labels = na.omit(as.integer(trunc_labels)); trunc_labels=trunc_labels[trunc_labels > 0]
  assert(trunc_labels, len=1, typ="integer")
  draw = as.logical(draw); assert(draw, len=1, alw=c(TRUE,FALSE))
  stats_print = as.logical(stats_print); assert(stats_print, len=1, alw=c(TRUE,FALSE))
  assert(viewport, len = 1, alw = c("ideas","data","max"))
  if(!all(typeof(trans) == "builtin") && (class(trans) == "function")) stop("'trans' should be a function")
  
  # shortcuts
  normalize = FALSE
  P = obj$pops
  R = obj$regions
  g = do.call(what=buildGraph, args=graph)
  tmp = c(g$f1, g$f2) %in% names(obj$features)
  if(!all(tmp)) stop(paste0("trying to plot a features not found in obj$features: ",  paste0(c(g$f1, g$f2)[!tmp], collapse=", ")))

  # defines binning (for histogram and density)
  if(missing(bin)) {
    if(g$type=="histogram") {
      nbin = g$bincount
    } else {
      nbin = g$BasePop[[1]]$densitybincount
    }
    nbin = na.omit(as.integer(nbin))
    if(length(nbin)==0) nbin=ifelse(g$type=="histogram", 520, 128)
    if(nbin==0) nbin=ifelse(g$type=="histogram", 520, 128)
  } else {
    nbin = na.omit(as.integer(bin)); assert(nbin, len=1, typ="integer")
  }

  # extracts graph information
  if(g$type == "histogram") {
    D = obj$features[,c("Object Number",g$f1)]
    names(D) = c("Object Number","x1")
  } else {
    D = obj$features[,c("Object Number",g$f1,g$f2)]
    names(D) = c("Object Number","x1","y1")
  }
  Xlim = c(g$xmin, g$xmax)
  Ylim = c(g$ymin, g$ymax)
  trans_x = g$xlogrange
  trans_y = g$ylogrange
  if(trans_x!="P") {
    trans_x = as.numeric(g$xlogrange)
    D[,"x2"] = smoothLinLog(D[,"x1"], hyper=trans_x, base=10)
    Xlim = smoothLinLog(Xlim, hyper=trans_x, base=10)
    Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
  } else {
    D[,"x2"] = D[,"x1"]
  }
  
  base_n = unlist(lapply(g$BasePop, FUN=function(x) x$name))
  reg_n = unlist(lapply(g$GraphRegion, FUN=function(x) x$name))
  shown_n = unlist(lapply(g$ShownPop, FUN=function(x) x$name))
  graph_n = unlist(lapply(g$GraphRegion, FUN=function(x) x$def))
  
  operators = c("And","Or","Not","(",")")
  displayed_n = unique(splitn(definition = g$order, all_names = c(base_n, graph_n, shown_n, "Selected Bin"), operators = operators))
  displayed_n = setdiff(displayed_n, "Selected Bin")
  displayed_r = rev(displayed_n)
  tmp = displayed_n %in% names(P)
  if(!all(tmp)) stop(paste0("trying to display a population not found in obj$pops: ",  paste0(displayed_n[!tmp], collapse=", ")))
  displayed_d = sapply(displayed_n, FUN=function(x) P[[x]]$obj)
  L = length(displayed_n)

  base_o = sapply(base_n, FUN=function(x) which(displayed_n%in%x))
  base_n = names(base_o)[order(base_o)]
  
  if(length(shown_n) == 0) {
    shown_o = NULL
    shown_n = NULL
  } else {
    shown_o = unlist(sapply(shown_n, FUN=function(x) which(displayed_n%in%x)))
    shown_n = names(shown_o)[order(shown_o)]
  }
  
  # subset data
  base = as.data.frame(sapply(base_n, FUN=function(x) P[[x]]$obj), stringsAsFactors = FALSE)
  sub = apply(base, 1, any)
  displayed_o = c(base_o, shown_o)
  D = cbind(D, displayed_d)
  D = D[sub, ]

  if(g$type=="histogram") {
    KEY = list("text"=list(displayed_n),
               "lines"=list(col = sapply(displayed_n, FUN=function(p) P[[p]][c("color","lightModeColor")][[color_mode]]),
                            lty = sapply(displayed_n, FUN=function(r) c(1,2,3,4,6)[match(g$BasePop[[displayed_o[r]]]$linestyle,c("Solid","Dash","Dot","DashDot","DashDotDot"))])))
    
    base_s = lapply(base_n, FUN=function(d) {
      np = sum(D[,d])
      if(np == 0) return(structure(rep(NA, 8), names = c("count","perc",
                                                          "Min.","1st Qu.","Median","Mean", "3rd Qu.","Max.")))
      p = c("count"=np, "perc"=100, summary(na.omit(D[D[,d],"x1"])))
    })
    kids_s = lapply(shown_n, FUN=function(s) {
      do.call(what = "rbind", args = lapply(base_n, FUN=function(d) {
        np = sum(D[,d])
        if(np == 0) return(structure(rep(NA, 8), names = c("count","perc",
                                                           "Min.","1st Qu.","Median","Mean", "3rd Qu.","Max.")))
        isin = D[,d] & D[,s]
        n = sum(isin)
        c("count"=n, "perc"=n/np*100, summary(na.omit(D[isin,"x1"])))
      }))
    })
    kids_r = lapply(reg_n, FUN=function(r) {
      do.call(what = "rbind", args = lapply(base_n, FUN=function(d) {
        alg = 3
        reg = R[[r]]
        coords = reg["x"]
        if(trans_x!="P") coords$x = smoothLinLog(coords$x, hyper=trans_x, base=10)
        np = sum(D[,d])
        if(np == 0) return(structure(rep(NA, 8), names = c("count","perc",
                                                           "Min.","1st Qu.","Median","Mean", "3rd Qu.","Max.")))
        isin = D[D[,d],"x2"]
        isin = (isin >= min(coords$x)) & (isin <= max(coords$x))
        n = sum(isin)
        c("count"=n, "perc"=n/np*100, summary(na.omit(D[isin,"x1"])))
      }))
    })
    stats = do.call(what=rbind, args=c(base_s, kids_s, kids_r))
    rnames = base_n
    if(length(reg_n) > 0) rnames = c(rnames, unlist(t(sapply(base_n, FUN = function(b) {if(b == "All") {graph_n} else {paste(reg_n, b, sep = " & ") }}))))
    rownames(stats) = rnames
    colnames(stats) = c(colnames(stats)[1:2], paste0("x-",colnames(stats)[3:8]))
    
    smooth = (g$histogramsmoothingfactor != 0)
    # br = do.breaks(range(D[,"x2"], na.rm = TRUE, finite = TRUE), nbin)
    if(viewport == "data") {
      Xlim = suppressWarnings(range(D[,"x1"], na.rm = TRUE, finite = TRUE))
      if(trans_x != "P") Xlim = smoothLinLog(Xlim, hyper=trans_x, base=10)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
    }
    if(viewport == "max") {
      regx = sapply(reg_n, FUN=function(r) {
        reg = R[[r]] 
        coords = reg[["x"]]
        return(c(reg$cx, coords))
      })
      Xlim = suppressWarnings(range(c(D[,"x1"], regx), na.rm = TRUE, finite = TRUE))
      if(trans_x != "P") Xlim = smoothLinLog(Xlim, hyper=trans_x, base=10)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
    }
    if(viewport == "ideas") {
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
      D[D[,"x2"] < Xlim[1], "x2"] <- Xlim[1] # D = D[(D[,"x2"] >= Xlim[1]) & (D[,"x2"] <= Xlim[2]), ]
      D[D[,"x2"] > Xlim[2], "x2"] <- Xlim[2] #
    }
    if(!all(is.finite(Xlim))) {
      Xlim = c(g$xmin, g$xmax)
      if(trans_x!="P") {
        Xlim = smoothLinLog(Xlim, hyper=trans_x, base=10)
        Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      }
      if(!all(is.finite(Xlim))) Xlim = c(-1, 1)
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
      D[D[,"x2"] < Xlim[1], "x2"] <- Xlim[1]
      D[D[,"x2"] > Xlim[2], "x2"] <- Xlim[2]
    }
    br = do.breaks(Xlim, nbin)
    
    type = "count"
    if(as.logical(g$freq)) type = "percent"
    if(nrow(D) > 0) {
      if(viewport == "ideas") {
        Ylim = c(g$ymin, g$ymax)
        if(Ylim[1] == Ylim[2]) Ylim = Ylim[1] + c(0,0.07)
      } else {
        Ylim = c(0,max(sapply(displayed_n, FUN=function(x) get_ylim(x=D[D[,x],"x2"], type=type, br=br))*1.07))
        if(Ylim[1] == Ylim[2]) Ylim = Ylim[1] + c(0,0.07)
      } 
      foo = histogram(~ D[,"x2"], auto.key=FALSE, xlim = Xlim, ylim = Ylim, main = trunc_string(g$title, trunc_labels),
                      scales =  myScales(x=list("hyper"=trans_x)), border = "transparent",
                      xlab =  trunc_string(g$xlabel, trunc_labels), ylab = g$ylabel,
                      nint = nbin, type = type, breaks = br, normalize = normalize,
                      panel = function(x, ...) { })
      for(l in L:1) {
        disp = displayed_n[l]
        if(any(D[,disp])) { # adds layer only if there is at least one point
          tmp = histogram(~ D[,"x2"], auto.key=FALSE, subset = D[,disp], alpha = 0.8,
                          col = P[[disp]][c("color","lightModeColor")][[color_mode]], border="transparent",
                          fill = as.logical(g$BasePop[[displayed_o[disp]]]$fill=="true"),
                          lty = c(1,2,3,4,6)[match(g$BasePop[[displayed_o[disp]]]$linestyle,c("Solid","Dash","Dot","DashDot","DashDotDot"))],
                          nint = nbin, type = type, breaks = br, normalize = normalize, Ylim = Ylim,
                          panel = function(x, type, breaks, normalize, fill, nint, border, col, alpha, lty, Ylim = Ylim, ...) {
                            if(smooth) {
                              pan_smooth(x=x, type=type, br=breaks, normalize=normalize, fill=fill, lwd=1, lty=lty, col=col, alpha=alpha, ylim=Ylim, bin=nint, border=border, factor=g$histogramsmoothingfactor)
                            } else {
                              pan_hist(x=x, type=type, br=breaks, normalize=normalize, fill=fill, lwd=1, lty=1, col=col, alpha=alpha, ylim=Ylim, bin=nint, border=border)
                            }
                            if(l == 1) {
                              if(any(c("panel","both")%in%add_key)) pan_key(key=c(KEY,"background"="lightgrey","alpha.background"=0.8), x = 0.02)
                              lapply(reg_n, FUN=function(r) {
                                reg = R[[r]] 
                                col = reg[c("color","lightcolor")][[color_mode]]
                                coords = reg[c("x","y")]
                                if(trans_x!="P") {
                                  coords$x = smoothLinLog(coords$x, hyper=trans_x, base=10)
                                  reg$cx = smoothLinLog(reg$cx, hyper=trans_x, base=10)
                                }
                                lab =  trunc_string(reg$label, trunc_labels)
                                panel.text(x=reg$cx, y=reg$cy*diff(Ylim), col=col, labels=lab, pos=4)
                                panel.lines(x=coords$x, y=coords$y*diff(Ylim),col=col)
                              })
                            }
                          })
          foo = foo + as.layer(tmp)
        }
      }
    } else {
      Ylim = c(g$ymin, g$ymax)
      foo = histogram(0 ~ 0, auto.key=FALSE, xlim = Xlim, ylim = Ylim, main =  trunc_string(g$title, trunc_labels), 
                      scales =  myScales(x=list("hyper"=trans_x)), border = "transparent",
                      xlab =  trunc_string(g$xlabel, trunc_labels), ylab = g$ylabel,
                      nint = nbin, type = type, normalize = normalize, Ylim = Ylim,
                      panel = function(x, Ylim = Ylim, ...) {
                        if(any(c("panel","both")%in%add_key)) pan_key(key=c(KEY,"background"="lightgrey","alpha.background"=0.8), x = 0.02)
                        lapply(reg_n, FUN=function(r) {
                          reg = R[[r]] 
                          col = reg[c("color","lightcolor")][[color_mode]]
                          coords = reg[c("x","y")]
                          if(trans_x!="P") {
                            coords$x = smoothLinLog(coords$x, hyper=trans_x, base=10)
                            reg$cx = smoothLinLog(reg$cx, hyper=trans_x, base=10)
                          }
                          lab = trunc_string(reg$label, trunc_labels)
                          panel.text(x=reg$cx, y=reg$cy*diff(Ylim), col=col, labels=lab, pos=4)
                          panel.lines(x=coords$x, y=coords$y*diff(Ylim), col=col)
                        })
                      })
    }
  } else {
    KEY = list("text"=list(displayed_r),
               "points"=list(col = sapply(P[displayed_r], FUN=function(p) p[c("color","lightModeColor")][[color_mode]]),
                             pch = sapply(P[displayed_r], FUN=function(p) p$style)))
    if(trans_y!="P") {
      trans_y = as.numeric(g$ylogrange)
      D[,"y2"] = smoothLinLog(D[,"y1"], hyper=trans_y, base=10)
      Ylim = smoothLinLog(Ylim, hyper=trans_y, base=10)
    } else {
      D[,"y2"] = D[,"y1"]
    }
    base_s = lapply(base_n, FUN=function(d) {
      np = sum(D[,d])
      if(np == 0) return(structure(rep(NA, 14), names = c("count","perc",
                                                         "Min.","1st Qu.","Median","Mean", "3rd Qu.","Max.",
                                                         "Min.","1st Qu.","Median","Mean", "3rd Qu.","Max.")))
      p = c("count"=np, "perc"=100, summary(na.omit(D[D[,d],"x1"])), summary(na.omit(D[D[,d],"y1"])))
    })
    kids_s = lapply(shown_n, FUN=function(s) {
      do.call(what = "rbind", args = lapply(base_n, FUN=function(d) {
        np = sum(D[,d])
        if(np == 0) return(structure(rep(NA, 14), names = c("count","perc",
                                                           "Min.","1st Qu.","Median","Mean","3rd Qu.","Max.",
                                                           "Min.","1st Qu.","Median","Mean", "3rd Qu.","Max.")))
        isin = D[,d] & D[,s]
        n = sum(isin)
        c("count"=n, "perc"=n/np*100, summary(na.omit(D[isin,"x1"])), summary(na.omit(D[isin,"y1"])))
      }))
    })
    kids_r = lapply(reg_n, FUN=function(r) {
      alg = 1
      reg = R[[r]]
      coords = reg[c("x","y")]
      if(trans_x!="P") coords$x = smoothLinLog(coords$x, hyper=trans_x, base=10)
      if(trans_y!="P") coords$y = smoothLinLog(coords$y, hyper=trans_y, base=10)
      if(reg$type=="oval") alg = 3
      if(reg$type=="rect") alg = 2
      do.call(what = "rbind", args = lapply(base_n, FUN=function(d) {
        np = sum(D[,d])
        if(np == 0) return(structure(rep(NA, 14), names = c("count","perc",
                                                           "Min.","1st Qu.","Median","Mean","3rd Qu.","Max.",
                                                           "Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")))
        isin = cpp_pnt_in_gate(pnts = cbind(D[D[,d],"x2"],D[D[,d],"y2"]), gate = cbind(coords$x,coords$y), algorithm = alg)
        n = sum(isin)
        c("count"=n, "perc"=n/np*100, summary(na.omit(D[isin,"x1"])), summary(na.omit(D[isin,"y1"])))
      }))
    })
    stats = do.call(what=rbind, args=c(base_s, kids_r, kids_s))
    rnames = base_n
    if(length(reg_n) > 0) rnames = c(rnames, unlist(t(sapply(base_n, FUN = function(b) {if(b == "All") {reg_n} else {paste(reg_n, b, sep = " & ") }}))))
    if(length(shown_n) > 0) rnames = c(rnames, unlist(sapply(shown_n, FUN = function(s) paste(base_n, s, sep = " & "))))
    rownames(stats) = rnames
    colnames(stats) = c(colnames(stats)[1:2], paste0("x-",colnames(stats)[3:8]), paste0("y-",colnames(stats)[9:14]))
    groups = NULL

    if(nrow(D)>0) if(g$type == "scatter") if(precision=="light") groups=apply(as.data.frame(D[,displayed_n]), 1, FUN=function(x) {
      tmp = which(x)[1]
      if(is.na(tmp)) return(NA)
      return(displayed_n[which(x)[1]])
    })
    
    if(viewport == "data") {
      Xlim = suppressWarnings(range(D[,"x1"], na.rm = TRUE, finite = TRUE))
      if(trans_x != "P") Xlim = smoothLinLog(Xlim, hyper=trans_x, base=10)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      Ylim = suppressWarnings(range(D[,"y1"], na.rm = TRUE, finite = TRUE))
      if(trans_y != "P") Ylim = smoothLinLog(Ylim, hyper=trans_y, base=10)
      Ylim = Ylim + c(-0.07,0.07)*diff(Ylim)
    }
    if(viewport == "max") {
      regx = sapply(reg_n, FUN=function(r) {
        reg = R[[r]] 
        coords = reg[["x"]]
        return(c(reg$cx, coords))
      })
      Xlim = suppressWarnings(range(c(D[,"x1"], regx), na.rm = TRUE, finite = TRUE))
      if(trans_x != "P") Xlim = smoothLinLog(Xlim, hyper=trans_x, base=10)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      regy = sapply(reg_n, FUN=function(r) {
        reg = R[[r]] 
        coords = reg[["y"]]
        return(c(reg$cy, coords))
      })
      Ylim = suppressWarnings(range(c(D[,"y1"], regy), na.rm = TRUE, finite = TRUE))
      if(trans_y != "P") Ylim = smoothLinLog(Ylim, hyper=trans_y, base=10)
      Ylim = Ylim + c(-0.07,0.07)*diff(Ylim)
    }
    if(!all(is.finite(Xlim))) {
      Xlim = c(g$xmin, g$xmax)
      if(trans_x!="P") {
        Xlim = smoothLinLog(Xlim, hyper=trans_x, base=10)
        Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      }
      if(!all(is.finite(Xlim))) Xlim = c(-1, 1)
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
    }
    if(!all(is.finite(Ylim))) {
      Ylim = c(g$ymin, g$ymax)
      if(trans_y!="P") {
        Ylim = smoothLinLog(Ylim, hyper=trans_y, base=10)
        Ylim = Ylim + c(-0.07,0.07)*diff(Ylim)
      }
      if(!all(is.finite(Ylim))) Ylim = c(-1, 1)
      if(Ylim[1] == Ylim[2]) Ylim = Ylim[1] + c(-0.07,0.07)
    }
    
    foo = xyplot(D[,"y2"] ~ D[,"x2"], auto.key=FALSE, xlim = Xlim, ylim = Ylim, main = trunc_string(g$title, trunc_labels), groups=groups,
                 scales =  myScales(x=list(hyper=trans_x), y=list(hyper=trans_y)),
                 xlab =  trunc_string(g$xlabel, trunc_labels), ylab = trunc_string(g$ylabel, trunc_labels),
                 panel = function(x, y, groups=NULL, ...) {
                   if(any(c("panel","both")%in%add_key)) if(g$type=="scatter") pan_key(key=c(KEY,"background"="lightgrey","alpha.background"=0.8), x = 0.02)
                   if(g$type == "density") panel.xyplot(x=x, y=y, pch=".", col=densCols(x=x,y=y,colramp=colorRampPalette(colConv(g$BasePop[[base_o[1]]][c("densitycolorsdarkmode","densitycolorslightmode")][[color_mode]])),nbin=nbin, transformation=trans))
                   if(g$type == "scatter") {
                     if(is.null(groups)) {
                       panel.xyplot(x=x[1], y=y[1], pch="", alpha=0)
                     } else {
                       by(data.frame("x"=x,"y"=y,"g"=groups, stringsAsFactors=FALSE), groups, FUN=function(d) {
                         disp = unique(d$g)
                         panel.xyplot(x=d$x, y=d$y, pch=P[[disp]]$style, col = P[[disp]][c("color","lightModeColor")][[color_mode]])
                       })
                     }
                   }
                   lapply(reg_n, FUN=function(r) {
                     reg = R[[r]]
                     k = reg[c("color","lightcolor")][[color_mode]]
                     coords = reg[c("x","y")]
                     if(trans_x!="P") {
                       coords$x = smoothLinLog(coords$x, hyper=trans_x, base=10)
                       reg$cx = smoothLinLog(reg$cx, hyper=trans_x, base=10)
                     }
                     if(trans_y!="P") {
                       coords$y = smoothLinLog(coords$y, hyper=trans_y, base=10)
                       reg$cy = smoothLinLog(reg$cy, hyper=trans_y, base=10)
                     }
                     if(reg$type=="rect") {
                       coords$x=c(coords$x[1],coords$x[1],coords$x[2],coords$x[2])
                       coords$y=c(coords$y[1],coords$y[2],coords$y[2],coords$y[1])
                     }
                     if(reg$type=="oval") {
                       coords = toEllipse(coords)
                     }
                     lab =  trunc_string(reg$label, trunc_labels)
                     panel.text(x=reg$cx, y=reg$cy, col=k, labels=lab, pos=4) 
                     panel.polygon(x=coords$x, y=coords$y, border=k, col="transparent", lwd=1, lty=1)
                   })
                 })
    if(precision=="full") if(g$type == "scatter") for(l in L:1) {
      disp = displayed_n[l]
      if(any(D[,disp])) { # adds layer only if there is at least one point
        tmp = xyplot(D[,"y2"] ~ D[,"x2"], pch = P[[disp]]$style, col = P[[disp]][c("color","lightModeColor")][[color_mode]], subset = D[,disp])
        foo = foo + as.layer(tmp)
      }
    }
  }
  lt <- custom.theme(bg=c("black","white")[color_mode], fg=c("white","black")[color_mode])
  lapply(names(lt), FUN=function(i) {
    if(i %in% c("add.text","axis.text","par.xlab.text","par.ylab.text","par.zlab.text","par.main.text","par.sub.text")) {
      lt[[i]]$font <- 2
      lt[[i]]$cex <- 1
      lt[[i]]$fontfamily <- "serif"
    }
    return(lt[[i]])
  })
  lt$fontsize$text <- 6
  lt$fontsize$points <- 4
  lt$grid.pars <- get.gpar()
  lt$grid.pars$fontfamily <- "serif"
  foo = update(foo, par.settings = lt)
  if(any(c("global","both")%in%add_key)) foo = update(foo, key=KEY)
  if(draw) {
    plot(foo)
    dv = dev.cur()
  }
  if(stats_print) print(stats)
  ret_order = names(D) %in% c("Object Number", "x1", "x2", "y1", "y2")
  displayed = lapply(obj$pops[displayed_n], FUN = function(p) {
    return(p[!(names(p) %in% "obj")])
  })
  ret = list("plot" = foo,
             "stats" = as.table(stats),
             "input" = list("data" = D[ ,c(which(ret_order), which(!ret_order))], 
                            "trunc_labels" = trunc_labels,
                            "title" = g$title,
                            "xlab" = g$xlab, "ylab" = g$ylab,
                            "xlim" = Xlim, "ylim" = Ylim, 
                            "trans_x" = trans_x, "trans_y" = trans_y,
                            "trans" = trans,
                            "order" = displayed_o,
                            "base" = g$BasePop,
                            "displayed" = displayed,
                            "regions" = obj$regions[reg_n],
                            "viewport" = viewport,
                            "bin" = nbin,
                            "type" = ifelse(g$type=="histogram", type, g$type),
                            "histogramsmoothingfactor" = g$histogramsmoothingfactor,
                            "normalize" = normalize,
                            "precision" = precision,
                            "mode" = color_mode))
  class(ret) <- "IFC_plot"
  return(invisible(ret))
  },
  finally = {
    while(!all(dv == dev.cur())) {
      dev.off(which = rev(dev.cur())[1])
    }
  })
}
