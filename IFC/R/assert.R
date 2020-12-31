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

#' @title Assert that Certain Conditions are Met
#' @description
#' Ensures that a variable respects several parameters
#' @param x variable to test
#' @param len integer vector of allowed length for x. Default is NULL, for not checking this parameter.
#' @param cla character vector of allowed classes of x. Default is NULL, for not checking this parameter.
#' @param typ character vector of allowed types of x. Default is NULL, for not checking this parameter.
#' @param alw allowed values for x. Default is NULL, for not checking this parameter.
#' @param fun function to execute when mandatory parameters are not met. Default is "stop". Allowed are "stop","warning","message","return".
#' @details /!\ alw parameter when used should be coercible to a logical, integer, numeric, complex or character vector. Otherwise, an error will be thrown.
#' @keywords internal
assert = function(x, len=NULL, cla=NULL, typ=NULL, alw=NULL, fun="stop") {
  foo = cpp_assert(x, len=len, cla=cla, typ=typ, alw=alw, fun=fun)
  if(!any(foo)) return(invisible(NULL))
  # since some conditions are not met everything is recomputed to allow user to check what unmet parameter(s).
  ele = c(len = NULL, cla = NULL, typ = NULL, alw = NULL)
  if(foo[1]) {
    tmp = len%in%length(x)
    rejected = x[!tmp]
    ele["len"] = paste0("len: is of length [",paste0(length(x), collapse=","),"]. Allowed length",ifelse(length(len)>1,"s are"," is"),": ", paste0(len,collapse=","))
  }
  if(foo[2]) {
    tmp = cla%in%class(x)
    ele["cla"] = paste0("cla: ['",paste0(class(x), collapse="','"),"'] not of required class",ifelse(length(cla)>1,"es",""),": ", paste0(paste0("'",cla[!tmp],"'"), collapse=" & "))
  }
  if(foo[3]) {
    tmp = typ%in%typeof(x)
    ele["typ"] = paste0("typ: [",paste0(typeof(x), collapse=","),"] not of required type",ifelse(length(typ)>1,"s",""),": ", paste0(paste0("'",typ[!tmp],"'"), collapse=" & "))
  }
  if(foo[4]) {
    tmp = x%in%alw
    rejected = x[!tmp]
    ele["alw"] = paste0("alw: [",paste0(head(rejected,5), collapse=","),ifelse(length(rejected)>5,", ...",""),"] ",ifelse(length(rejected)>1,"are","is")," not allowed. Allowed values are: ", paste0(alw,collapse=","))
  }
  args = list(paste0(paste0("'", paste0(as.character(substitute(x)),collapse="$"),"':\n"), paste(" -", ele, collapse = "\n")))
  if(fun == "warning") args = c(args, call.=FALSE, immediate.=TRUE)
  do.call(what = fun, args = args)
}
