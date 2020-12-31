#' Time Stamps for File Names
#' 
#' Generating file name ready iso time stamps.
#' 
#' @param ts one or more POSIX time stamp
#' @param sep separators to be used for formatting
#'
#' @return Returns timestamp string in format yyyy-mm-dd_HH_MM_SS ready to be used 
#' safely in file names on various operating systems. 
#' 
#' 
#' @export
#'
#'
#' @examples
#' 
#' time_stamp()
#' time_stamp( Sys.time() - 10000 )
#' 
#' 
time_stamp <- 
  function( ts = Sys.time(), sep = c("-", "_", "_")){
    
    # process parameter
    if ( length(sep) < 3){
      sep <- rep(sep, 3)
    }
    
    # putting together format
    fmt <- 
      paste0(
        "%Y", sep[1], "%m", sep[1], "%d", 
        sep[2], 
        "%H", sep[3], "%M", sep[3], "%S"
      )
    
    # execute formatting
    format(
      x      = ts, 
      format = fmt
    )
  }