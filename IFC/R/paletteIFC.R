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

#' @title R/IDEAS Color Palette Mapping
#' @description
#' Maps colors between IDEAS and R.
#' @param x either "", "palette","palette_R", to_light, to_dark. Default is "".
#' @param col a compatible color to transform to color or lightModeColor. Default is "White".\cr
#' if 'x' == to_light, function will convert 'col' to lightModeColor.\cr
#' if 'x' == to_dark, function will convert 'col' to color.\cr
#' if 'col' is not found or 'x' is anything else then a data.frame of compatible colors is returned.
#' @return IFC palette of available colors.
#' @export
paletteIFC <- function(x = c("","palette","palette_R","to_light","to_dark")[1], col = "White") {
  assert(x, len = 1, alw = c("","palette","palette_R","to_light","to_dark"))
  Software_Colors=data.frame(matrix(c("White","LightSkyBlue","CornflowerBlue","MediumSlateBlue","Blue","Aquamarine","MediumSpringGreen","Cyan","DarkTurquoise",   
                                      "Teal","Yellow","Gold","DarkKhaki","Lime","Green","Lime","Wheat","SandyBrown","Orange",     
                                      "Tomato","Red","Pink","HotPink","Plum","Magenta","DarkOrchid","LightCoral","IndianRed",       
                                      "LightGray","Gray","Black",
                                      "Black","CornflowerBlue","CornflowerBlue","MediumSlateBlue","Blue","Aquamarine","MediumSpringGreen","Teal","DarkTurquoise",    
                                      "Teal","Gold","Gold","IndianRed","Green","Green", "Lime","DarkOrange","Tomato","DarkOrange",       
                                      "Tomato","Red","DeepPink","DeepPink","DarkOrchid","Magenta","DarkOrchid","IndianRed","IndianRed",        
                                      "Black","Gray","White"), ncol=2, byrow = FALSE), stringsAsFactors = FALSE)
  names(Software_Colors) = c("color", "lightModeColor")
  # Software_Colors = Software_Colors[!(duplicated(Software_Colors[,1]) | duplicated(Software_Colors[,2])),]
  Software_Colors$color_R = gsub("^Teal","Cyan4", Software_Colors$color)
  Software_Colors$color_R = gsub("^Green","Green4", Software_Colors$color_R)
  # Software_Colors$color_R = gsub("^Gray","Darkgray", Software_Colors$color_R)
  Software_Colors$color_R = gsub("^Lime","Chartreuse", Software_Colors$color_R)
  Software_Colors$lightModeColor_R = gsub("^Teal","Cyan4", Software_Colors$lightModeColor)
  Software_Colors$lightModeColor_R = gsub("^Green","Green4", Software_Colors$lightModeColor_R)
  # Software_Colors$lightModeColor_R = gsub("^Gray","Darkgray", Software_Colors$lightModeColor_R)
  Software_Colors$lightModeColor_R = gsub("^Lime","Chartreuse", Software_Colors$lightModeColor_R)
  if(x == "palette") return(as.character(unique(c(Software_Colors[,1], Software_Colors[,2]))))
  if(x == "palette_R") return(tolower(as.character(unique(c(Software_Colors[,3], Software_Colors[,4])))))
  columns = c(0,0)
  M = FALSE
  if(x == "to_light") { columns = c(1,3,2,4) }
  if(x == "to_dark") { columns = c(2,4,1,3) }
  if(!identical(columns,c(0,0))) M = Software_Colors[,columns[1]]==col | Software_Colors[,columns[2]]==col
  if(any(M)) return(Software_Colors[which(M), columns[3:4]])
  return(Software_Colors)
}
