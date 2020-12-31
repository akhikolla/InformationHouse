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

#' @title Batch File Writer
#' @description
#' Writes an XML file to batch files
#' @param batch list of batch nodes as created by \code{\link{buildBatch}}.
#' @return It invisibly returns full path of xml batch file.
#' @export
ExportToBATCH = function(batch) {
  path = unique(unlist(lapply(batch, FUN=function(x) x$batch_dir)))
  if(length(path)!=1) stop("you ar trying to save batch file in multiple places which is not allowed")
  if(attr(x = batch[[1]]$batch_dir, which = "tempdir")) { # batch_dir is tempdir(), allowed by CRAN
    write_to = normalizePath(paste(path,"batch.xml",sep="/"), winslash = "/", mustWork = FALSE)
    write_xml(file = write_to, xml_new_node(name="batches",.children=lapply(batch, FUN=function(x) x$xml)), encoding = "utf-8")
    message(paste0("\n######################\n", write_to, "\nhas been successfully exported\n"))
  } else { # batch_dir is not tempdir(), so file will be saved to tempdir() and user will receive a message to MANUALLY rename the file to the desired location
    write_to = normalizePath(paste(tempdir(check = TRUE),"batch.xml",sep="/"), winslash = "/", mustWork = FALSE)
    write_xml(file = write_to, xml_new_node(name="batches",.children=lapply(batch, FUN=function(x) x$xml)), encoding = "utf-8")
    message(paste0("\n######################\n",
                   write_to,
                   "\nhas been successfully exported",
                   sprintf("\nTo complete the process, please run:\nfile.copy(from = '%s', to = '%s', overwrite = TRUE)\n", write_to, normalizePath(paste(path,"batch.xml",sep="/"), winslash = "/", mustWork = FALSE))))
  }
  return(invisible(write_to))
}
