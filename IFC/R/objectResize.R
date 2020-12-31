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

#' @title Object Resizing
#' @description
#' Resizes mat to new dimensions.
#' @param mat a numeric matrix.
#' @param size a length 2 integer vector of final dimensions of the image, height 1st and width 2nd. Default is c(0,0) for no change.
#' @param add_noise if TRUE adds normal noise when size is larger than mat dimensions using rnorm(), from \pkg{Rcpp}. Default is TRUE.
#' @param random_seed a single value, interpreted as an integer, or NULL to be used with set.seed() from \pkg{base} when 'add_noise' is set to TRUE. Default is NULL.
#' @param bg mean value of the background added if add_noise is TRUE. Default is 0.
#' @param sd standard deviation of the background added if add_noise is TRUE. Default is 0.
#' @return a resized matrix with padding background if desired size is larger than original mat dimensions.
#' @keywords internal
objectResize <- function(mat, size = c(0,0), add_noise = TRUE, random_seed = NULL, bg = 0, sd = 0) {
  if(length(size) != 2) return(mat)
  if(add_noise) {
    set.seed(random_seed)
    on.exit(set.seed(NULL))
  }
  return(cpp_resize(mat = mat, new_height = size[1], new_width = size[2], add_noise = add_noise, bg = bg, sd= sd))
}
