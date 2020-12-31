# size.R Size of Array or Vector
#
#     Copyright (C) 2016 MeteoSwiss
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#' Size of Array or Vector
#' 
#' Return dimension of array or length of vector.
#' 
#' @param x array or vector 
#' 
#' @examples
#' tm <- toymodel()
#' 
#' sapply(tm, size)
#' 
#' @seealso \code{\link{veriApply}}
#'   
#' @keywords utilities
#' @export
size <- function(x){
 if (is.vector(x)) length(x) else dim(x) 
}