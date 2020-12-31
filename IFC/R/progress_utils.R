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

#' @title Progress Bar Initializer
#' @description
#' Initializes a progress bar.
#' @param session the shiny session object, as provided by shinyServer to the server function. Default is missing, to use "txtProgressBar" or "winProgressBar" (on Windows).
#' @param title,label character strings, giving the 'title'(='message' for shiny progress bar) and the 'label'(='detail' for shiny progress bar).
#' @param min,max (finite) numeric values for the extremes of the progress bar. Must have 'min' < 'max'.
#' @param initial initial value for the progress bar.
#' @param steps (finite) numeric value for the number of individual chunk of the progress bar. Default is 21.
#' @param width only apply when 'session' is missing,the width of the progress bar. If missing, the default, will be NA for "txtProgressBar" and 300 for "winProgressBar".
#' @param style does not apply for "winProgressBar", the style of the bar. If missing, the default, will be 3 "txtProgressBar" and getShinyOption("progress.style", default = "notification") for shiny progress bar
#' @param char only apply for "txtProgressBar", the character (or character string) to form the progress bar.
#' @param file only apply for "txtProgressBar", an open connection object or "" which indicates the console: stderr() might be useful here. Default is "".
#' @details shiny progress bar will be available only if shiny package is found. 
#' @return pb an object of class `IFC_progress` containing a progress bar of class `txtProgressBar`, `winProgressBar` or `Progress`.
#' @keywords internal
newPB <- function(session,
                  title, label, 
                  min = 0, max = 1,
                  initial = 0,
                  steps = 21,
                  width,
                  style,
                  char = "=",
                  file = "") {
  fun = stop
  args = list("newPB: can't create progress bar")
  if(requireNamespace("shiny", quietly = TRUE) && length(session) != 0) {
    args = list(session = session,
                min = 0,
                max = steps- 1)
    if(missing(style)) {
      args = c(args, list(style = shiny::getShinyOption("progress.style", default = "notification")))
    } else {
      if(style[1] %in% c("old", "notification")) {
        args = c(args, list(style = style))
      } else {
        args = c(args, list(style = shiny::getShinyOption("progress.style", default = "notification")))
      }
    }
    fun = shiny::Progress$new
    bar = do.call(what = fun, args = args)
    if(is.finite(initial) && initial < max) {
      bar$set(value = initial)
    } else {
      bar$set(value = min)
    }
    typ = 3
  } else {
    if(.Platform$OS.type == "windows") {
      args = list(min = 0,
                  max = steps - 1,
                  initial = initial)
      if(missing(title)) {
        args = c(args, list(title = "R progress bar"))
      } else {
        args = c(args, list(title = title))
      }
      if(missing(label)) {
        args = c(args, list(label = " ")) # change "" to " " to be allow to pass label 
      } else {                            # with setPB even if not created at first
        args = c(args, list(label = label))
      }
      if(missing(width)) {
        args = c(args, list(width = 300))
      } else {
        args = c(args, list(width = width))
      }
      fun = winProgressBar
      bar = do.call(what = fun, args = args)
      if(is.finite(initial) && initial < max) {
        setWinProgressBar(bar, value = initial)
      } else {
        setWinProgressBar(bar, value = min)
      }
      typ = 2
    } else {
      args = list(min = 0,
                  max = steps - 1, 
                  initial = initial, 
                  char = char,
                  file = file)
      mess = c()
      if(!missing(title)) mess = title
      if(!missing(label)) mess = c(mess, label)
      if(length(mess) != 0 )cat("\n", mess, "\n", sep = "\n", file = file)
      if(missing(style)) {
        args = c(args, list(style = 3))
      } else {
        args = c(args, list(style = style))
      }
      if(missing(width)) {
        args = c(args, list(width = NA))
      } else {
        args = c(args, list(width = width))
      }
      fun = txtProgressBar
      bar = do.call(what = fun, args = args)
      if(is.finite(initial) && initial < max) {
        setTxtProgressBar(bar, value = initial)
      } else {
        setTxtProgressBar(bar, value = min)
      }
      typ = 1
    }
  } 
  ans = list(bar = bar,
             seq = seq(min, max, length.out = steps),
             steps = steps,
             typ = typ)
  class(ans) = "IFC_progress"
  return(ans)
}

#' @title Progress Bar Updater
#' @description
#' Updates a progress bar.
#' @param pb an object of class `IFC_progress` containing a progress bar of class `txtProgressBar`, `winProgressBar` or `Progress`.
#' @param value new value for the progress bar.
#' @param title,label character strings, giving the 'title'(='message' for shiny progress bar) and the 'label'(='detail' for shiny progress bar).
#' @keywords internal
setPB <- function(pb, value = NULL, title = NULL, label = NULL) {
  val = sum(value > pb$seq)
  switch(pb$typ,
         { # typ 1
           if(pb$bar$getVal() == val) return(NULL)
           return(setTxtProgressBar(pb = pb$bar, value = val, title = title, label = label))
         },
         { # typ 2
           if(getWinProgressBar(pb$bar) == val) return(NULL)
           return(setWinProgressBar(pb = pb$bar, value = val, title = title, 
                                    label = sprintf("%3.f%% - %s", 100 * val/(pb$steps-1), label)))
         },
         { # typ 3
           if(pb$bar$getValue() == val) return(NULL)
           return(pb$bar$set(value = val, message = title,
                             detail = sprintf("%3.f%% - %s", 100 * val/(pb$steps-1),label)))
         })
}

#' @title Progress Bar Terminator
#' @description
#' Terminates a progress bar.
#' @param pb an object of class `IFC_progress` containing a progress bar of class `txtProgressBar`, `winProgressBar` or `Progress`.
#' @keywords internal
endPB <- function(pb) {
  if(pb$typ == 3) {
    pb$bar$close()
  } else {
    close(con = pb$bar)
  }
}
