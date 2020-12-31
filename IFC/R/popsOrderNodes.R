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

#' @title IFC_pops Population Level Ordering
#' @description
#' Helper that sort populations so that populations that depend from other are placed after.
#' @param pops list of populations.
#' @keywords internal
popsOrderNodes <- function(pops) {
  assert(pops, cla = c("IFC_pops","Affiliated"))
  K = class(pops)
  i=1; l=length(pops)
  while(i<l) {
    pop=pops[[i]]
    index=setdiff(c(pop$base,pop$names),"")
    index=unlist(lapply(index, function(x) which(x==names(pops))))
    index=index[index>i]
    if(length(index)!=0) {
      pops = c(pops[index],pops[setdiff(1:l,index)])
      i=1
    } else {
      i=i+1
    }
  }
  class(pops) = c(K, "Ordered")
  return(pops)
}
