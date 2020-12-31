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

#' @title IFC_pops Computation
#' @description
#' Function used to compute `IFC_pops`\cr
#' It requires pops, regions and features.
#' @param pops list of populations.
#' @param regions list of regions.
#' @param features data.frame of features.
#' @param pnt_in_poly_algorithm algorithm used to determine if object belongs to a polygon region or not. Default is 1.\cr
#' Note that for the moment only 1(Trigonometry) is available.
#' @param pnt_in_poly_epsilon epsilon to determine if object belongs to a polygon region or not. It only applies when algorithm is 1. Default is 1e-12.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @param ... other arguments to be passed.
#' @source For pnt_in_poly_algorithm, Trigonometry, is an adaptation of Jeremy VanDerWal's code \url{https://github.com/jjvanderwal/SDMTools}
#' @return an object of class `IFC_pops`.
#' @keywords internal
popsCompute <- function(pops, regions, features, pnt_in_poly_algorithm = 1, 
                        pnt_in_poly_epsilon = 1e-12, display_progress = TRUE, 
                        title_progress = "", ...) {
  dots = list(...)
  # coerce pops to buildPopulation() checking
  pops = lapply(pops, FUN=function(p) do.call(what = "buildPopulation", args = p))
  class(pops) = "IFC_pops"
  
  # compute pops
  pops = popsGetAffiliation(pops)
  pops = popsOrderNodes(pops)
  return(popsWithin(pops = pops, regions = regions, features = features, 
                    pnt_in_poly_algorithm = pnt_in_poly_algorithm, 
                    pnt_in_poly_epsilon = pnt_in_poly_epsilon, 
                    display_progress = display_progress, title_progress = title_progress, ...))
}
