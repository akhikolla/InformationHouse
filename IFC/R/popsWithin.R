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

#' @title IFC_pops Checker for Object Membership to Populations
#' @description
#' Helper that will be used by popsCompute to determine which objects are within populations or not.
#' @param pops list of populations.
#' @param regions list of regions.
#' @param features dataframe of features.
#' @param pnt_in_poly_algorithm algorithm used to determine if object belongs to a polygon region or not. Default is 1.\cr
#' Note that for the moment only 1(Trigonometry) is available.
#' @param pnt_in_poly_epsilon epsilon to determine if object belongs to a polygon region or not. It only applies when algorithm is 1. Default is 1e-12.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @param ... other arguments to be passed.
#' @source For pnt_in_poly_algorithm, Trigonometry, is an adaptation of Jeremy VanDerWal's code \url{https://github.com/jjvanderwal/SDMTools}
#' @keywords internal
popsWithin <- function(pops, regions, features, pnt_in_poly_algorithm = 1, pnt_in_poly_epsilon = 1e-12,
                       display_progress = TRUE, title_progress = "", ...) {
  dots = list(...)
  # several checks
  assert(pops, cla = c("IFC_pops","Affiliated","Ordered"))
  assert(regions, cla = "IFC_regions")
  assert(features, cla = "IFC_features")
  pnt_in_poly_algorithm = as.integer(pnt_in_poly_algorithm); assert(pnt_in_poly_algorithm, len = 1, alw = 1)
  pnt_in_poly_epsilon = as.numeric(pnt_in_poly_epsilon); pnt_in_poly_epsilon = pnt_in_poly_epsilon[pnt_in_poly_epsilon>0]; pnt_in_poly_epsilon = pnt_in_poly_epsilon[is.finite(pnt_in_poly_epsilon)]
  assert(pnt_in_poly_epsilon, len = 1, typ = "numeric")
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  assert(title_progress, len = 1, typ = "character")
  
  K = class(pops)
  L = length(pops)
  obj_number = nrow(features)
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
  }
  for(i in 1:L) {
    fx_pos = NULL
    fy_pos = NULL
    pop=pops[[i]]
    # changes styles to R compatible
    style_tmp = pop$style==c("Simple Dot","Cross","Plus","Empty Circle","Empty Diamond","Empty Square","Empty Triangle","Solid Diamond","Solid Square","Solid Triangle")
    if(any(style_tmp)) pops[[i]]$style=c(20, 4, 3, 1, 5, 0, 2, 18, 15, 17)[style_tmp]
    # changes colors to R compatible
    if(pop$color=="Teal") {pops[[i]]$color="Cyan4"}
    if(pop$color=="Green") {pops[[i]]$color="Green4"}
    if(pop$color=="Lime") {pops[[i]]$color="Chartreuse"}
    if(pop$lightModeColor=="Teal") {pops[[i]]$lightModeColor="Cyan4"}
    if(pop$lightModeColor=="Green") {pops[[i]]$lightModeColor="Green4"}
    if(pop$lightModeColor=="Lime") {pops[[i]]$lightModeColor="Chartreuse"}
    
    switch(pop$type,
           "B" = { 
             pops[[i]]$obj=rep(TRUE,obj_number)
           }, 
           "G" = {
             pop_pos=which(names(regions)==pop$region)
             fx_pos=which(names(features)==pop$fx)
             x=features[,fx_pos]
             xlim=as.numeric(regions[[pop_pos]]$x)
             if(regions[[pop_pos]]$type == "line") {
               xlim=range(xlim)
               pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & x>=xlim[1] & x<=xlim[2]
             } else {
               fy_pos=which(names(features)==pop$fy)
               y=features[,fy_pos]
               ylim=as.numeric(regions[[pop_pos]]$y)
               if(regions[[pop_pos]]$xlogrange != "P") {
                 x = smoothLinLog(x, hyper = as.numeric(regions[[pop_pos]]$xlogrange), base = 10)
                 xlim = smoothLinLog(xlim, hyper = as.numeric(regions[[pop_pos]]$xlogrange), base = 10)
               }
               if(regions[[pop_pos]]$ylogrange != "P") {
                 y = smoothLinLog(y, hyper = as.numeric(regions[[pop_pos]]$ylogrange), base = 10)
                 ylim = smoothLinLog(ylim, hyper = as.numeric(regions[[pop_pos]]$ylogrange), base = 10)
               }
               switch(regions[[pop_pos]]$type, 
                      "oval" = {
                        pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & cpp_pnt_in_gate(pnts=cbind(x,y), gate = cbind(xlim,ylim), algorithm = 3)
                      },
                      "poly" = {
                        pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & cpp_pnt_in_gate(pnts=cbind(x,y), gate = cbind(xlim,ylim), algorithm = pnt_in_poly_algorithm, epsilon = pnt_in_poly_epsilon)
                      },
                      "rect" = {
                        pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & cpp_pnt_in_gate(pnts=cbind(x,y), gate = cbind(xlim,ylim), algorithm = 2)
                      })
             }
           }, 
           "C" = {
             pop_def_tmp=gsub("^And$","&",pop$split)
             pop_def_tmp=gsub("^Or$","|",pop_def_tmp)
             pop_def_tmp=gsub("^Not$","!", pop_def_tmp)
             for(i_popn in which(pop_def_tmp%in%pop$names)) {pop_def_tmp[i_popn]=paste0("`",pop_def_tmp[i_popn],"`")}
             comb_tmp=sapply(pops[pop$names], FUN=function(i_pop) i_pop$obj)
             if(obj_number == 1) {
               comb_tmp=all(comb_tmp)
             } else {
               comb_tmp=eval(parse(text=paste0(pop_def_tmp,collapse=" ")), as.data.frame(comb_tmp, stringsAsFactors = FALSE))
             }
             pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & comb_tmp
           }, 
           "T" = {
             if(length(pop$obj) != obj_number) {
               Kp = class(pop$obj)
               if(Kp%in%"numeric" | Kp%in%"integer") {
                 if((obj_number <= max(pop$obj)) | (min(pop$obj) < 0) | any(duplicated(pop$obj))) stop(paste0("trying to export a tagged population with element(s) outside of objects acquired: ", pop$name))
                 pops[[i]]$obj=rep(FALSE,obj_number)
                 pops[[i]]$obj[pop$obj+1]=TRUE
               } else {
                 if(!Kp%in%"logical") stop(paste0("trying to export a tagged population with element(s) outside of objects acquired: ", pop$name))
               }
             }
             if(sum(pops[[i]]$obj)==0) stop(paste0("trying to export a tagged population with element(s) outside of objects acquired: ", pop$name))
             if(obj_number != length(pops[[i]]$obj)) stop(paste0("trying to export a tagged population with element(s) outside of objects acquired: ", pop$name))
           })
    if(display_progress) {
      setPB(pb, value = i, title = title_progress, label = "extacting populations")
    }
  }
  class(pops) = c(setdiff(K, "IFC_pops"), "IFC_pops", "Processed")
  return(pops)
}
