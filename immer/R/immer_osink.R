## File Name: immer_osink.R
## File Version: 0.04


immer_osink <- function(file)
{
    CDM::osink( file=file, suffix=paste0(  "__SUMMARY.Rout") )
}
