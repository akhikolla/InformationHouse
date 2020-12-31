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

.pkgenv <- new.env(emptyenv())

.onLoad <- function(libname, pkgname) {
  data_avl <- requireNamespace("IFCdata", quietly = TRUE)
  .pkgenv[["data_avl"]] <- data_avl
  invisible()
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- c(paste("Package 'IFC' version", packageVersion("IFC")),
           "\nType `citation(\"IFC\")` for citing this R package in publications.")
  if(!.pkgenv[["data_avl"]]) {
    msg <- c(msg,
             "\nTo use examples in this package, you should install 'IFCdata' package.",
             "\nTo install 'IFCdata' package, run `install.packages('IFCdata', repos = 'https://gitdemont.github.io/IFCdata/', type = 'source')`.")
  }
  packageStartupMessage(paste(strwrap(msg), collapse = "\n"))      
  invisible()
}
