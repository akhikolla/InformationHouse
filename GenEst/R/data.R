# Mock ------------------------------
#' A mock example data set
#'
#' A template dataset used for testing purposes.  Dataset containing SE, CP, SS,
#'  DWP, and CO data. Data are mostly random without patterns.
#'
#' @format A list with 5 items:
#' \describe{
#'   \item{SE}{Searcher efficiency trial data}
#'   \item{CP}{Carcass persistence trial data}
#'   \item{SS}{Search schedule data}
#'   \item{DWP}{Density weighted proportion of area searched data}
#'   \item{CO}{Carcass observations}   
#' }
#' @source \code{mock}
"mock"

# Cleared ------------------------------
#' Wind cleared plot (60m) Search Example
#' 
#' A complete example data set for estimating fatalities from 60 m cleared plots
#'  at 23 out of 100 searches at a wind power facility.  Data on carcass 
#'  observations (CO) from a search of all terrain out to 60m from each of 100 
#'  turbines at a theoretical site, field trials for estimating carcass 
#'  persistence (CP) and searcher efficiency (SE), search schedule (SS) 
#'  parameters (for example, which turbines were searched on which days), and 
#'  density weighted proportion (DWP) of area searched at each turbine (which is 
#'  an area adjustment factor to account for incomplete  search coverage).
#'
#' @format \code{wind_cleared} is a list with 5 elements:
#' \describe{
#'   \item{\code{SE}}{Searcher efficiency trial data}
#'   \item{\code{CP}}{Carcass persistence trial data}
#'   \item{\code{SS}}{Search schedule parameters}
#'   \item{\code{DWP}}{Density weighted proportion of area searched}
#'   \item{\code{CO}}{Carcass observations}
#' }
#' 
#' @section Searcher Efficiency (\code{SE}):
#'  \code{$SE} is a data frame with each row representing the fate of a single
#'  carcass in the searcher efficiency trials. There are columns for:
#' \describe{
#'   \item{\code{pkID}}{unique ID for each carcass}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or
#'     \code{"sml"} bird}
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"}}
#'   \item{\code{Visibility}}{indicator for visibility class of the ground, with
#'      \code{"RP"} for carcasses placed on a road or turbine pad, \code{"M"}
#'      for moderate visibility (e.g., plowed field; short, sparse vegetation),
#'      or \code{"D"} for difficult visibility}
#'   \item{\code{"s1",...,"s5"}}{fate of carcass on the 1st, 2nd, 3rd, 4th, and
#'      5th search after placement. A value of 1 implies that a carcass was
#'      discovered by searchers, 0 implies the carcass was present but not
#'      discovered, and any other value is interpreted as "no search" or
#'      "carcass not present" and ignored in the model. In this data set,
#'      \code{NA} indicates that a carcass had been previously discovered and
#'      removed from the field. A user may use a variety of values to
#'      differentiate different reasons no search was conducted or the carcass
#'      was not present. For example, "SN" could be used to indicate that the
#'      turbine was not searched because of snow, or "NS" to indicate the search
#'      was not scheduled in that location at that time, or "SC" to indicate the
#'      carcass had been removed by scavengers prior to the search.}
#' }
#' 
#' @section Carcass Persistence (\code{CP}):
#'  \code{$CP} is a data frame with each row representing the fate of a single
#'  carcass in the carcass persistence trials. There are columns for:
#' \describe{
#'   \item{\code{cpID}}{unique ID for each carcass}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or 
#'   \code{"sml"} bird}
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"}}
#'   \item{\code{Visibility}}{indicator for visibility class of the ground, with
#'    \code{"RP"} for carcasses placed on a road or turbine pad, \code{"M"} for 
#'    moderate visibility (e.g., plowed field; short, sparse vegetation), or 
#'    \code{"D"} for difficult visibility.} \item{\code{LastPresent},
#'    \code{FirstAbsent}}{endpoints of the interval bracketing the time the carcass
#'    was scavenged or otherwise removed from the field. For example, 
#'    \code{LastPresent = 2.04}, \code{FirstAbsent = 3.21} indicates that the carcass was
#'    last observed 2.04 days after being placed in the field and was noted 
#'    missing 3.21 days after being placed. If the precise time of carcass 
#'    removal is known (e.g., recorded by camera), then \code{LastPresent} and
#'    \code{FirstAbsent} should be set equal to each other. If a carcass persists
#'    beyond the last day of the field trial, \code{LastPresent} is the last time it
#'    was observed and \code{FirstAbsent} is entered as \code{Inf} or \code{NA}.}
#' }
#' 
#' @section Search Schedule (\code{SS}):
#'  \code{$SS} is a data frame with a row for each date a turbine at the site
#'  was searched, a column of \code{SearchDate}s, and a column for each turbine.
#'  In addition, there is a column to indicate the \code{Season}. A column with
#'  search dates and columns for each turbine searched are required. Other
#'  columns are optional.
#' \describe{
#'   \item{\code{SearchDate}}{columns of dates on which at least one turbine was
#'    searched. Format in this data is \code{"\%Y-\%m-\%d CDT"}, but time zone 
#'    (\code{CDT}) is optional. A time stamp may be included if desired (e.g., 
#'    \code{2018-03-20 02:15:41}). Alternatively, \code{\\} can be used in place
#'     of \code{-}.}
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"} to
#'    indicate which season the search was conducted in. \code{Season} is
#'    optional but may be used as a temporal covariate for fatality estimates.}
#'   \item{\code{t1}, etc.}{unique ID for all turbines that were searched on at
#'    least one search date. Values are either 1 or 0, indicating whether the
#'    given turbine (column) was searched or not on the given date (row).}
#' }
#'  
#' @section Density Weighted Proportion (\code{DWP}):
#'  \code{$DWP} is a data frame with a row for each turbine and columns for
#'   each carcass size class. Values represent the density-weighted proportion
#'   of the searched area for each size (or the fraction of carcasses that fall
#'   in the searched area).
#' \describe{
#'   \item{\code{Turbine}}{unique ID for each turbine. IDs match those used in 
#'   the \code{$CO} data frame and the column names in the \code{$SS} data.}
#'   \item{\code{Size}}{\code{bat}, \code{sml}, \code{med}, \code{lrg}}
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"} to 
#'   indicate which season the search was conducted in. \code{Season} is 
#'   optional but may be used as a temporal covariate for fatality estimates.}}
#' 
#' @section Carcass Observations (\code{CO}):
#' \code{$CO} is a data frame with a row for carcass observed in the carcass
#'   searches and a number of columns giving information about the given carcass
#'   (date found, size, species, etc.)
#' \describe{
#'   \item{\code{carcID}}{unique identifier for each carcass: \code{"x30"},
#'   \code{"x46"}, etc.}
#'   \item{\code{Turbine}}{identifier for which turbine the given carcass was
#'     found at: \code{"t19"}, \code{"t65"}, \code{"t49"}, etc.}
#'   \item{\code{TurbineType}}{the type of turbine: \code{"X"}, \code{"Y"} or
#'     \code{"Z"}. }
#'   \item{\code{DateFound}}{dates entered in the same format as in
#'     \code{$SS$SearchDate}. Every date entered here is (and must be) included
#'      in the search schedule (\code{$SS$SearchDate})}
#'   \item{\code{Visibility}}{visibility class: \code{"RP"}, \code{"M"}, or
#'     \code{"D"}, as described in \code{$CP} and \code{$SE}}
#'   \item{\code{Species}}{species of the carcass: \code{"BA"}, \code{"BB"},
#'     \code{"BC"}, \code{"BD"}, \code{"BE"}, \code{"LA"}, \code{"LB"},
#'     \code{"LD"}, \code{"LE"}, \code{"MA"}, \code{"MB"}, \code{"SA"},
#'     \code{"SB"}, \code{"SC"}, \code{"SD"}, \code{"SE"}, \code{"SF"},
#'     \code{"SG"}}
#'   \item{\code{SpeciesGroup}}{species group: \code{"bat0"}, \code{"bat1"},
#'    \code{"brd1"}, \code{"brd2"}, \code{"brd3"}}
#'   \item{\code{Size}}{size: \code{"bat"}, \code{"lrg"}, \code{"med"},
#'    \code{"sml"}}
#'   \item{\code{Distance}}{distance from the turbine}
#' }
#'   
#' @source \code{wind_cleared}
"wind_cleared"

# RP ------------------------------
#'  Wind Road and Pad (120m) Example
#'
#'  This example dataset is based on 120 m radius road and pad searches of all 
#'  100 turbines at a theoretical site.  The simulated site consists of 100 
#'  turbines, searched on roads and pads only, out to 120 meters.  Search 
#'  schedule differs by turbine and season, with more frequent searches in the 
#'  fall, and a subset of twenty turbines searched at every scheduled search.
#'  
#' Data on carcass observations (CO) from searches, field trials for estimating 
#' carcass persistence (CP) and searcher efficiency (SE), search schedule (SS) 
#' parameters (for example, which turbines were searched on which days), and 
#' density weighted proportion (DWP) of area searched at each turbine (which is 
#' an area adjustment factor to account for incomplete search coverage).
#'
#' @format \code{wind_RP} is a list with 5 elements:
#' \describe{
#'   \item{\code{SE}}{Searcher efficiency trial data}
#'   \item{\code{CP}}{Carcass persistence trial data}
#'   \item{\code{SS}}{Search schedule parameters}
#'   \item{\code{DWP}}{Density weighted proportion of area searched}
#'   \item{\code{CO}}{Carcass observations}
#' }
#' 
#' @section Searcher Efficiency (\code{SE}):
#'  \code{$SE} is a data frame with each row representing the fate of a single
#'  carcass in the searcher efficiency trials. There are columns for:
#' \describe{
#'   \item{\code{pkID}}{unique ID for each carcass}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or
#'     \code{"sml"} bird}
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"}}
#'   \item{\code{"s1",...,"s5"}}{fate of carcass on the 1st, 2nd, 3rd, 4th, and
#'      5th search after placement. A value of 1 implies that a carcass was
#'      discovered by searchers, 0 implies the carcass was present but not
#'      discovered, and any other value is interpreted as "no search" or
#'      "carcass not present" and ignored in the model. In this data set,
#'      \code{NA} indicates that a carcass had been previously discovered and
#'      removed from the field. A user may use a variety of values to
#'      differentiate different reasons no search was conducted or the carcass
#'      was not present. For example, "SN" could be used to indicate that the
#'      turbine was not searched because of snow, or "NS" to indicate the search
#'      was not scheduled in that location at that time, or "SC" to indicate the
#'      carcass had been removed by scavengers prior to the search.}
#' }
#' @section Carcass Persistence (\code{CP}):
#'  \code{$CP} is a data frame with each row representing the fate of a single
#'  carcass in the carcass persistence trials. There are columns for:
#' \describe{
#'   \item{\code{cpID}}{unique ID for each carcass}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or 
#'   \code{"sml"} bird.} 
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"}}
#'   \item{\code{LastPresent}, \code{FirstAbsent}}{endpoints of the interval bracketing the
#'   time the carcass was scavenged or otherwise removed from the field. For 
#'   example, \code{LastPresent = 2.04}, \code{FirstAbsent = 3.21} indicates that the carcass
#'    was last observed 2.04 days after being placed in the field and was noted 
#'    missing 3.21 days after being placed. If the precise time of carcass
#'    removal is known (e.g., recorded by camera), then \code{LastPresent} and
#'    \code{FirstAbsent} should be set equal to each other. If a carcass persists
#'    beyond the last day of the field trial, \code{LastPresent} is the last time it
#'    was observed and \code{FirstAbsent} is entered as \code{Inf} or \code{NA}.}
#' }
#'
#' @section Search Schedule (\code{SS}):
#'  \code{$SS} is a data frame with a row for each date a turbine at the site
#'  was searched, a column of \code{SearchDate}s, and a column for each turbine.
#'  In addition, there is a column to indicate the \code{Season}. A column with
#'  search dates and columns for each turbine searched are required. Other
#'  columns are optional.
#' \describe{
#'   \item{\code{SearchDate}}{columns of dates on which at least one turbine was
#'    searched. Format in this data is \code{"\%Y-\%m-\%d CDT"}, but time zone
#'    (\code{CDT}) is optional. A time stamp may be included if desired
#'    (e.g., \code{2018-03-20 02:15:41}). Alternatively, \code{\\} can be used 
#'    in place of \code{-}.}
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"} to
#'    indicate which season the search was conducted in. \code{Season} is
#'    optional but may be used as a temporal covariate for fatality estimates.}
#'   \item{\code{t1}, etc.}{unique ID for all turbines that were searched on at
#'    least one search date. Values are either 1 or 0, indicating whether the
#'    given turbine (column) was searched or not on the given date (row).}
#' }
#'
#' @section Density Weighted Proportion (\code{DWP}):
#'  \code{$DWP} is a data frame with a row for each turbine and columns for
#'   each carcass size class. Values represent the density-weighted proportion
#'   of the searched area for each size (or the fraction of carcasses that fall
#'   in the searched area).
#' \describe{
#'   \item{\code{Turbine}}{unique ID for each turbine. IDs match those used in
#'    the \code{$CO} data frame and the column names in the \code{$SS} data.}
#'   \item{\code{bat}}{DWP associated with size class Bat.}
#'   \item{\code{sml}}{DWP associated with size class Small.}
#'   \item{\code{med}}{DWP associated with size class Medium.}
#'   \item{\code{lrg}}{DWP associated with size class Large.}
#' }
#'
#' @section Carcass Observations (\code{CO}):
#' \code{$CO} is a data frame with a row for carcass observed in the carcass
#'   searches and a number of columns giving information about the given
#'   carcass (date found, size, species, etc.)
#' \describe{
#'   \item{\code{carcID}}{unique identifier for each carcass: \code{"x30"},
#'   \code{"x46"}, etc.}
#'   \item{\code{Turbine}}{identifier for which turbine the given carcass was
#'     found at: \code{"t19"}, \code{"t65"}, \code{"t49"}, etc.}
#'   \item{\code{TurbineType}}{the type of turbine: \code{"X"}, \code{"Y"} or
#'     \code{"Z"}. }
#'   \item{\code{DateFound}}{dates entered in the same format as in
#'     \code{$SS$SearchDate}. Every date entered here is (and must be) included
#'      in the search schedule (\code{$SS$SearchDate}}
#'   \item{\code{Species}}{species of the carcass: \code{"BA"}, \code{"BB"},
#'     \code{"BC"}, \code{"BD"}, \code{"BE"}, \code{"LA"}, \code{"LB"},
#'     \code{"LD"}, \code{"LE"}, \code{"MA"}, \code{"MB"}, \code{"SA"},
#'     \code{"SB"}, \code{"SC"}, \code{"SD"}, \code{"SE"}, \code{"SF"},
#'     \code{"SG"}}
#'   \item{\code{SpeciesGroup}}{species group: \code{"bat0"}, \code{"bat1"},
#'    \code{"brd1"}, \code{"brd2"}, \code{"brd3"}}
#'   \item{\code{Size}}{size: \code{"bat"}, \code{"lrg"}, \code{"med"},
#'    \code{"sml"}}
#'   \item{\code{Distance}}{distance from the turbine}
#' }
#'
#' @source \code{wind_RP}
"wind_RP"

# RPbat ------------------------------
#'  Wind Bat-Only Road and Pad (120m) Example
#'
#'  This example dataset considers only bats found on 120 m radius road and pad 
#'  searches of all 100 turbines at a theoretical site.  The simulated site 
#'  consists of 100 turbines, searched on roads and pads only, out to 120 
#'  meters.  Search schedule differs by turbine and season, with more frequent 
#'  searches in the fall, and a subset of twenty turbines searched at every 
#'  scheduled search.
#'  
#' Data on carcass observations (CO) from searches, field trials for estimating 
#' carcass persistence (CP) and searcher efficiency (SE), search schedule (SS) 
#' parameters (for example, which turbines  were searched on which days), and 
#' density weighted proportion (DWP) of area searched at each  turbine (which is
#' an area adjustment factor to account for incomplete search coverage).
#'
#' @format \code{wind_RPbat} is a list with 5 elements:
#' \describe{
#'   \item{\code{SE}}{Searcher efficiency trial data}
#'   \item{\code{CP}}{Carcass persistence trial data}
#'   \item{\code{SS}}{Search schedule parameters}
#'   \item{\code{DWP}}{Density weighted proportion of area searched}
#'   \item{\code{CO}}{Carcass observations}
#' }
#' 
#' @section Searcher Efficiency (\code{SE}):
#'  \code{$SE} is a data frame with each row representing the fate of a single
#'  carcass in the searcher efficiency trials. There are columns for:
#' \describe{
#'   \item{\code{pkID}}{unique ID for each carcass}
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"}}
#'   \item{\code{"s1",...,"s5"}}{fate of carcass on the 1st, 2nd, 3rd, 4th, and
#'      5th search after placement. A value of 1 implies that a carcass was
#'      discovered by searchers, 0 implies the carcass was present but not
#'      discovered, and any other value is interpreted as "no search" or
#'      "carcass not present" and ignored in the model. In this data set,
#'      \code{NA} indicates that a carcass had been previously discovered and
#'      removed from the field. A user may use a variety of values to
#'      differentiate different reasons no search was conducted or the carcass
#'      was not present. For example, "SN" could be used to indicate that the
#'      turbine was not searched because of snow, or "NS" to indicate the search
#'      was not scheduled in that location at that time, or "SC" to indicate the
#'      carcass had been removed by scavengers prior to the search.}
#' }
#' @section Carcass Persistence (\code{CP}):
#'  \code{$CP} is a data frame with each row representing the fate of a single
#'  carcass in the carcass persistence trials. There are columns for:
#' \describe{
#'   \item{\code{cpID}}{unique ID for each carcass}
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"}}
#'   \item{\code{LastPresent}, \code{FirstAbsent}}{endpoints of the interval bracketing
#'    the time the carcass was scavenged or otherwise removed from the field.
#'    For example, \code{LastPresent = 2.04}, \code{FirstAbsent = 3.21} indicates that the
#'    carcass was last observed 2.04 days after being placed in the field and
#'    was noted missing 3.21 days after being placed. If the precise time of
#'    carcass removal is known (e.g., recorded by camera), then \code{LastPresent} and
#'    \code{FirstAbsent} should be set equal to each other. If a carcass persists
#'    beyond the last day of the field trial, \code{LastPresent} is the last time it
#'    was observed and \code{FirstAbsent} is entered as \code{Inf} or \code{NA}.}
#' }
#' @section Search Schedule (\code{SS}):
#'  \code{$SS} is a data frame with a row for each date a turbine at the site
#'  was searched, a column of \code{SearchDate}s, and a column for each turbine.
#'  In addition, there is a column to indicate the \code{Season}. A column with
#'  search dates and columns for each turbine searched are required. Other
#'  columns are optional.
#' \describe{
#'   \item{\code{SearchDate}}{columns of dates on which at least one turbine was
#'    searched. Format in this data is \code{"\%Y-\%m-\%d CDT"}, but time zone
#'    (\code{CDT}) is optional. A time stamp may be included if desired
#'    (e.g., \code{2018-03-20 02:15:41}). Alternatively, \code{\\} can be used 
#'    in place of \code{-}.}
#'   \item{\code{Season}}{\code{"spring"}, \code{"summer"}, or \code{"fall"} to
#'    indicate which season the search was conducted in. \code{Season} is
#'    optional but may be used as a temporal covariate for fatality estimates.}
#'   \item{\code{t1}, etc.}{unique ID for all turbines that were searched on at
#'    least one search date. Values are either 1 or 0, indicating whether the
#'    given turbine (column) was searched or not on the given date (row).}
#' }
#' @section Density Weighted Proportion (\code{DWP}):
#'  \code{$DWP} is a data frame with a row for each turbine and columns for
#'   each carcass size class. Values represent the density-weighted proportion
#'   of the searched area for each size (or the fraction of carcasses that fall
#'   in the searched area).
#' \describe{
#'   \item{\code{Turbine}}{unique ID for each turbine. IDs match those used in
#'    the \code{$CO} data frame and the column names in the \code{$SS} data.}
#'   \item{\code{bat}}{Contains the DWP for each turbine, with respect to size 
#'   class (in this case, bats only.}
#' }
#' 
#' @section Carcass Observations (\code{CO}):
#' \code{$CO} is a data frame with a row for carcass observed in the carcass
#'   searches and a number of columns giving information about the given
#'   carcass (date found, size, species, etc.)
#' \describe{
#'   \item{\code{carcID}}{unique identifier for each carcass: \code{"x30"},
#'   \code{"x46"}, etc.}
#'   \item{\code{Turbine}}{identifier for which turbine the given carcass was
#'     found at: \code{"t19"}, \code{"t65"}, \code{"t49"}, etc.}
#'   \item{\code{TurbineType}}{the type of turbine: \code{"X"}, \code{"Y"} or
#'     \code{"Z"}. }
#'   \item{\code{DateFound}}{dates entered in the same format as in
#'     \code{$SS$SearchDate}. Every date entered here is (and must be) included
#'      in the search schedule (\code{$SS$SearchDate}}
#'   \item{\code{Species}}{species of the carcass: \code{"BA"}, \code{"BB"},
#'     \code{"BC"}, \code{"BD"}, \code{"BE"}, \code{"LA"}, \code{"LB"},
#'     \code{"LD"}, \code{"LE"}, \code{"MA"}, \code{"MB"}, \code{"SA"},
#'     \code{"SB"}, \code{"SC"}, \code{"SD"}, \code{"SE"}, \code{"SF"},
#'     \code{"SG"}}
#'   \item{\code{SpeciesGroup}}{species group: \code{"bat0"}, \code{"bat1"},
#'    \code{"brd1"}, \code{"brd2"}, \code{"brd3"}}
#'   \item{\code{Distance}}{Distance from the turbine.}
#' }
#' @source \code{wind_RPbat}
"wind_RPbat"



# Trough ------------------------------
#' Trough-based solar thermal power simulated example
#' 
#' An example data set for estimating fatalities from a trough-based solar 
#' thermal electric power generation facility.  The simulated site is inspected 
#' daily along ten 2000 meter long transects,  which run  north-south. Observers
#' look up to 150 meters away down the rows created by troughs  (east-west).  
#' One sided distance sampling will be used, with observers looking consistently
#' in  one cardinal direction as they travel through the facility. A sitewide 
#' clearout search is implemented before the first scheduled winter search.
#'  
#' The dataset consists of five parts: Data on carcass observations (CO) from 
#' daily searches, field trials for estimating carcass persistence (CP) and 
#' searcher efficiency (SE), search schedule (SS), and density weighted 
#' proportion (DWP) of area searched for the rows within each transect (which is
#' an area adjustment factor to account for incomplete search coverage).
#'
#' @format \code{solar_trough} is a list with 5 elements:
#' \describe{
#'   \item{\code{SE}}{Searcher efficiency trial data}
#'   \item{\code{CP}}{Carcass persistence trial data}
#'   \item{\code{SS}}{Search schedule parameters}
#'   \item{\code{DWP}}{Density weighted proportion of area searched}
#'   \item{\code{CO}}{Carcass observations}
#' }
#' @section Searcher Efficiency (\code{SE}):
#'  \code{$SE} is a data frame with each row representing the fate of a single
#'  carcass in the searcher efficiency trials. There are columns for:
#' \describe{
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"}, or 
#'   \code{"fall"}}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or 
#'   \code{"sml"} bird}
#'   \item{\code{"Search1",...,"Search5"}}{fate of carcass on the 1st, 2nd, 3rd,
#'    4th, and 5th search after placement. A value of 1 implies that a carcass 
#'    was discovered by searchers, 0 implies the carcass was present but not
#'    discovered, and any other value is interpreted as "no search" or
#'    "carcass not present" and ignored in the model. In this data set,
#'    \code{NA} indicates that a carcass had been previously discovered and
#'    removed from the field. A user may use a variety of values to
#'    differentiate different reasons no search was conducted or the carcass
#'    was not present. For example, "NS" to indicate the search
#'    was not scheduled in that location at that time, or "SC" to indicate the
#'    carcass had been removed by scavengers prior to the search.}
#'   \item{\code{Distance}}{the distance a carcass was placed from the 
#'    observer's transect.}
#' }
#' @section Carcass Persistence (\code{CP}):
#'  \code{$CP} is a data frame with each row representing the fate of a single
#'  carcass in the carcass persistence trials. There are columns for:
#' \describe{
#'   \item{\code{Index}}{unique ID for each carcass}
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"}, or 
#'   \code{"fall"}}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or 
#'   \code{"sml"} bird}
#'   \item{\code{LastPresent}, \code{FirstAbsent}}{endpoints of the interval bracketing the
#'   time the carcass was scavenged or otherwise removed from the field. For 
#'   example, \code{LastPresent = 2.04}, \code{FirstAbsent = 3.21} indicates that the carcass
#'   was last observed 2.04 days after being placed in the field and was noted 
#'   missing 3.21 days after being placed. If the precise time of carcass 
#'   removal is known (e.g., recorded by camera), then \code{LastPresent} and
#'   \code{FirstAbsent} should be set equal to each other. If a carcass persists
#'   beyond the last day of the field trial, \code{LastPresent} is the last time it was
#'    observed and \code{FirstAbsent} is entered as \code{Inf} or \code{NA}.}
#' }
#'
#' @section Search Schedule (\code{SS}):
#'  \code{$SS} is a data frame with a row for each date a transect at the site 
#'  was searched, a column of \code{SearchDate}s, and a column for each 
#'  transect. In addition, there is an optional column to indicate the 
#'  \code{Season}. The columns for distinct area (array) and the date column
#'  are required, and the names of the columns for search areas must match the 
#'  names of areas used in the DWP and CO files.
#' \describe{
#'   \item{\code{SearchDate}}{columns of dates when a transect was searched. 
#'   Format in this data is \code{"\%Y-\%m-\%d CDT"}, but time zone (\code{CDT})
#'    is optional. A time stamp may be included if desired (e.g., 
#'    \code{2018-03-20 02:15:41}). Alternatively, \code{\\} can be used in
#'    place of \code{-}.}
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"}, or 
#'    \code{"fall"} to indicate which season the search was conducted in. 
#'    \code{Season} is optional but may be used as a temporal covariate for 
#'    fatality estimates.}
#' }
#' 
#' @section Density Weighted Proportion (\code{DWP}):
#'  \code{$DWP} is a data frame with a row for each transect and columns for 
#'  each carcass size class (labels must match those of the class factors in the
#'   carcass observation file). Values represent the density-weighted proportion
#'   of the searched area for each size (or the fractionof carcasses that fall
#'   in the searched area). Since the whole site was searched, DWP is uniformly
#'   set equal to 1.
#' \describe{
#'   \item{\code{Unit}}{unique ID for each transect.  IDs match those used in 
#'   the \code{$CO} data frame and the column names in the \code{$SS} data.}
#'   \item{\code{bat}}{DWP associated with size class Bat}
#'   \item{\code{sml}}{DWP associated with size class Small}
#'   \item{\code{med}}{DWP associated with size class Medium}
#'   \item{\code{lrg}}{DWP associated with size class Large}
#' }
#'
#' @section Carcass Observations (\code{CO}):
#' \code{$CO} is a data frame with a row for carcass observed in the carcass 
#' searches and a number of columns giving information about the given carcass 
#' (date found, size, species, etc.)
#' \describe{
#'   \item{\code{Index}}{unique identifier for each carcass.}
#'   \item{\code{Unit}}{identifier for which transect the given carcass was 
#'   found at.  Values must match with DWP Transect values Search Schedule 
#'   column names.}
#'   \item{\code{Species}}{species of the carcass: \code{"BA"}, \code{"BB"},
#'     \code{"BC"}, \code{"BD"}, \code{"BE"}, \code{"LA"}, \code{"LB"},
#'     \code{"LD"}, \code{"LE"}, \code{"MA"}, \code{"MB"}, \code{"SA"},
#'     \code{"SB"}, \code{"SC"}, \code{"SD"}, \code{"SE"}, \code{"SF"},
#'     \code{"SG"}}
#'   \item{\code{Size}}{size: \code{"bat"}, \code{"lrg"}, \code{"med"},
#'    \code{"sml"}}
#'   \item{\code{Row}}{Optional indicator of which row within an array a carcass
#'    was found at.}
#'   \item{\code{Distance}}{The perpendicular distance from the searcher's 
#'    transect at which the carcass was discovered at.}
#'   \item{\code{DateFound}}{dates entered in the same format as in 
#'    \code{$SS$SearchDate}. 
#'   Every date entered here is (and must be) included in the search schedule 
#'   (\code{$SS$SearchDate})}
#'   \item{\code{X}}{UTM Easting of carcass.}
#'   \item{\code{Y}}{UTM Northing of carcass.}
#'  }
#'
#' @source \code{solar_trough}
"solar_trough"


# PV ------------------------------
#' Photovoltaic Example Dataset
#' 
#' An example data set for estimating fatalities from a large photovoltaic solar
#'  generation facility.
#'  
#' The simulated site is organized into 300 arrays of panels.  As observers walk
#'  north-south along paths between arrays, they look east or west down rows 
#'  between solar panels 150 meters long, with 38 searchable rows per array.  
#'  Observers consistently look for animals down one cardinal direction, making 
#'  this a one-sided distance sample.  Searches are scheduled on a seven day
#'  rotation, with 60 arrays searched per weekday.  A sitewide clearout search 
#'  is implemented before the first scheduled winter search.
#'  
#' The dataset consists of five parts: Data on carcass observations (CO) from 
#' array searches, field trials for estimating carcass persistence (CP) and 
#' searcher efficiency (SE), search schedule (SS), and density weighted 
#' proportion (DWP) of area searched at each array (which is an area adjustment 
#' factor to account for incomplete search coverage).
#'
#' @format \code{solar_PV} is a list with 5 elements:
#' \describe{
#'   \item{\code{SE}}{Searcher efficiency trial data}
#'   \item{\code{CP}}{Carcass persistence trial data}
#'   \item{\code{SS}}{Search schedule parameters}
#'   \item{\code{DWP}}{Density weighted proportion of area searched}
#'   \item{\code{CO}}{Carcass observations}
#' }
#' @section Searcher Efficiency (\code{SE}):
#'  \code{$SE} is a data frame with each row representing the fate of a single
#'  carcass in the searcher efficiency trials. There are columns for:
#' \describe{
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"}, or 
#'   \code{"fall"}}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or
#'     \code{"sml"} bird}
#'   \item{\code{"Search1",...,"Search5"}}{fate of carcass on the 1st, 2nd, 3rd,
#'    4th, and 5th search after placement. A value of 1 implies that a carcass 
#'    was discovered by searchers, 0 implies the carcass was present but not
#'      discovered, and any other value is interpreted as "no search" or
#'      "carcass not present" and ignored in the model. In this data set,
#'      \code{NA} indicates that a carcass had been previously discovered and
#'      removed from the field. A user may use a variety of values to
#'      differentiate different reasons no search was conducted or the carcass
#'      was not present. For example, "NS" to indicate the search
#'      was not scheduled in that location at that time, or "SC" to indicate the
#'      carcass had been removed by scavengers prior to the search.}
#'   \item{\code{Distance}}{the distance a carcass was placed from the 
#'   observer's transect. Used in determining probability to detect with 
#'   distance sampling.}
#' }
#' @section Carcass Persistence (\code{CP}):
#'  \code{$CP} is a data frame with each row representing the fate of a single
#'  carcass in the carcass persistence trials. There are columns for:
#' \describe{
#'   \item{\code{Index}}{unique ID for each carcass}
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"}, or 
#'   \code{"fall"}}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or
#'    \code{"sml"} bird}
#'   \item{\code{LastPresent}, \code{FirstAbsent}}{endpoints of the interval bracketing
#'    the time the carcass was scavenged or otherwise removed from the field.
#'    For example, \code{LastPresent = 2.04}, \code{FirstAbsent = 3.21} indicates that the
#'    carcass was last observed 2.04 days after being placed in the field and
#'    was noted missing 3.21 days after being placed. If the precise time of
#'    carcass removal is known (e.g., recorded by camera), then \code{LastPresent} and
#'    \code{FirstAbsent} should be set equal to each other. If a carcass persists
#'    beyond the last day of the field trial, \code{LastPresent} is the last time it
#'    was observed and \code{FirstAbsent} is entered as \code{Inf} or \code{NA}.}
#' }
#'
#'
#' @section Search Schedule (\code{SS}):
#'  \code{$SS} is a data frame with a row for each date an array at the site was
#'   searched, a column of \code{SearchDate}s, and a column for each array. In 
#'   addition, there is an optional column to indicate the \code{Season}. The 
#'   columns for distinct area (array) and the date column are required, and the
#'   names of the columns for search areas must match the names of areas used in 
#'   the DWP and CO files.
#' \describe{
#'   \item{\code{SearchDate}}{columns of dates when arrays were searched. Format
#'    in this data is \code{"\%Y-\%m-\%d CDT"}, but time zone (\code{CDT}) is 
#'    optional. A time stamp may be included if desired (e.g., 
#'    \code{2018-03-20 02:15:41}). Alternatively, \code{\\} can be used in
#'    place of \code{-}.}
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"},
#'    or \code{"fall"} to indicate which season the search was conducted in. 
#'    \code{Season} is optional but may be used as a temporal covariate for 
#'    fatality estimates.}
#' }
#' 
#' @section Density Weighted Proportion (\code{DWP}):
#'  \code{$DWP} is a data frame with a row for each array and columns for each 
#'  carcass size class (labels must match those of the class factors in the 
#'  carcass observation file). Values represent the density-weighted proportion 
#'  of the searched area for each size (or the fraction of carcasses that fall 
#'  in the searched area).  In this example, observers walk along transects 
#'  separated by 150 meters, and search coverage is assumed to be 100%, i.e.,
#'  DWP = 1 for each unit. This requires that carcasses be placed at random
#'  locations in the field, even at distances from the transects that would make
#'  it unlikely to observe small carcasses.
#' \describe{
#'   \item{\code{Unit}}{unique ID for each array.  IDs match those used in the 
#'   \code{$CO} data frame and the column names in the \code{$SS} data.}
#'   \item{\code{bat}}{DWP associated with size class Bat}
#'   \item{\code{sml}}{DWP associated with size class Small}
#'   \item{\code{med}}{DWP associated with size class Medium}
#'   \item{\code{lrg}}{DWP associated with size class Large}
#' }
#'
#' @section Carcass Observations (\code{CO}):
#' \code{$CO} is a data frame with a row for carcass observed in the carcass
#'   searches and a number of columns giving information about the given
#'   carcass (date found, size, species, etc.)
#' \describe{
#'   \item{\code{Index}}{unique identifier for each carcass.}
#'   \item{\code{Unit}}{identifier for which unit the given carcass was found 
#'   at: \code{"arc19"}, \code{"arc65"}, etc, for arcs in the outer heliostat 
#'   field, or \code{"center"}, indicating the inner heliostat field.}
#'   \item{\code{Species}}{species of the carcass: \code{"BA"}, \code{"BB"},
#'     \code{"BC"}, \code{"BD"}, \code{"BE"}, \code{"LA"}, \code{"LB"},
#'     \code{"LD"}, \code{"LE"}, \code{"MA"}, \code{"MB"}, \code{"SA"},
#'     \code{"SB"}, \code{"SC"}, \code{"SD"}, \code{"SE"}, \code{"SF"},
#'     \code{"SG"}}
#'   \item{\code{Size}}{size: \code{"bat"}, \code{"lrg"}, \code{"med"},
#'    \code{"sml"}}
#'   \item{\code{Row}}{Optional indicator of which row within an array a carcass
#'    was found at.}
#'   \item{\code{Distance}}{The perpendicular distance from the searcher's 
#'   transect at which the carcass was discovered at.}
#'   \item{\code{DateFound}}{dates entered in the same format as in 
#'   \code{$SS$SearchDate}. Every date entered here is (and must be) included in
#'    the search schedule 
#'   (\code{$SS$SearchDate}}
#'   \item{\code{X}}{UTM Easting of carcass.}
#'   \item{\code{Y}}{UTM Northing of carcass.}
#' }
#'
#' @source \code{solar_PV}
"solar_PV"


# Power Tower ------------------------------
#' Power Tower Example Dataset
#' 
#' An example data set for estimating fatalities from a concentrating
#' solar-thermal (power tower) generation facility.  
#'  
#' The simulated site consists of a single tower generating approximately 130 
#' MW.  The tower is surrounded by a 250 meter radius circular inner field of 
#' heliostats, searched on a weekly schedule.  From the inner circle, 18 
#' concentric rings of heliostats 50 meters deep extend to the boundaries of the
#' simulated site.  Rings are subdivided into 8 arcs each, with arcs 1-8 
#' immediately adjacent to the central circle.  Arcs are search using distance 
#' sampling techniques on a weekly schedule, with 29 arcs searched per weekday.
#'  
#' There are two sources of mortality simulated: flux and non-flux (collision or
#' unknown cause).Flux carcasses are generated (weibull) about the tower, with 
#' 5\% to be found in the outer field. Non-flux mortality is assumed uniform 
#' across the site.
#'   
#' The dataset consists of five parts: Data on carcass observations (CO) from 
#' inner and outer heliostat searches, field trials for estimating carcass 
#' persistence (CP) and searcher efficiency (SE), search schedule (SS), and 
#' density weighted proportion (DWP) of area searched at each turbine (which is 
#' an area adjustment factor to account for incomplete search coverage).
#'
#' @format \code{solar_powerTower} is a list with 5 elements:
#' \describe{
#'   \item{\code{SE}}{Searcher efficiency trial data}
#'   \item{\code{CP}}{Carcass persistence trial data}
#'   \item{\code{SS}}{Search schedule parameters}
#'   \item{\code{DWP}}{Density weighted proportion of area searched}
#'   \item{\code{CO}}{Carcass observations}
#' }
#' @section Searcher Efficiency (\code{SE}):
#'  \code{$SE} is a data frame with each row representing the fate of a single
#'  carcass in the searcher efficiency trials. There are columns for:
#' \describe{
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"}, or 
#'   \code{"fall"}}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or
#'     \code{"sml"} bird}
#'   \item{\code{Field}}{indicates carcass placed in inner or outer heliostat 
#'   field, with levels \code{"inner"} or \code{outer}.}
#'   \item{\code{"Search1",...,"Search5"}}{fate of carcass on the 1st, 2nd, 3rd,
#'    4th, and 5th search after placement. A value of 1 implies that a carcass 
#'    was discovered by searchers, 0 implies the carcass was present but not
#'    discovered, and any other value is interpreted as "no search" or
#'    "carcass not present" and ignored in the model. In this data set,
#'    \code{NA} indicates that a carcass had been previously discovered and
#'    removed from the field. A user may use a variety of values to
#'    differentiate different reasons no search was conducted or the carcass
#'    was not present. For example, "NS" to indicate the search
#'    was not scheduled in that location at that time, or "SC" to indicate the
#'    carcass had been removed by scavengers prior to the search.}
#' }
#' @section Carcass Persistence (\code{CP}):
#'  \code{$CP} is a data frame with each row representing the fate of a single
#'  carcass in the carcass persistence trials. There are columns for:
#' \describe{
#'   \item{\code{cpID}}{unique ID for each carcass}
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"}, or 
#'   \code{"fall"}}
#'   \item{\code{Size}}{\code{"bat"}; or \code{"lrg"}, \code{"med"}, or
#'    \code{"sml"} bird}
#'   \item{\code{LastPresent}, \code{FirstAbsent}}{endpoints of the interval bracketing
#'    the time the carcass was scavenged or otherwise removed from the field.
#'    For example, \code{LastPresent = 2.04}, \code{FirstAbsent = 3.21} indicates that the
#'    carcass was last observed 2.04 days after being placed in the field and
#'    was noted missing 3.21 days after being placed. If the precise time of
#'    carcass removal is known (e.g., recorded by camera), then \code{LastPresent} and
#'    \code{FirstAbsent} should be set equal to each other. If a carcass persists
#'    beyond the last day of the field trial, \code{LastPresent} is the last time it
#'    was observed and \code{FirstAbsent} is entered as \code{Inf} or \code{NA}.}
#' }
#'
#' @section Search Schedule (\code{SS}):
#'  \code{$SS} is a data frame with a row for each date an arc at the site
#'  was searched, a column of \code{SearchDate}s, and a column for each arc, and
#'  one column at the end for the inner heliostat field, labeled \code{center}.
#'  In addition, there is a column to indicate the \code{Season}. A column with 
#'  search dates and columns for each distinct area (arcs and center) searched 
#'  are required. Other columns are optional.
#' \describe{
#'   \item{\code{SearchDate}}{columns of dates on which an arc was searched. 
#'   Format in this data is \code{"\%Y-\%m-\%d CDT"}, but time zone (\code{CDT})
#'    is optional. A time stamp may be included if desired (e.g., 
#'    \code{2018-03-20 02:15:41}). Alternatively, \code{\\} can be used in place
#'     of \code{-}.}
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"}, or 
#'   \code{"fall"} to indicate which season the search was conducted in. 
#'   \code{Season} is optional but may be used as a temporal covariate for 
#'   fatality estimates.}
#' }
#' 
#' @section Density Weighted Proportion (\code{DWP}):
#'  \code{$DWP} is a data frame with a row for each arc and columns for each
#'  carcass size class (labels must match those of the class factors in the
#'  carcass observation file). Values represent the density-weighted proportion
#'  of the searched area for each size (or the fraction of carcasses that fall
#'  in the searched area).  In this example, within the inner field (center)
#'  observers are unobstructed in ability to discover carcasses, for a DWP of 1.
#'  In the outer heliostat field observers walk along transects separated by 50 
#'  meters, but the entire area is surveyed, so DWP = 1.
#' \describe{
#'   \item{\code{Unit}}{unique ID for each arc, plus one labeled \code{center} 
#'   for the inner heliostat field.  IDs match those used in the \code{$CO} data
#'    frame and the column names in the \code{$SS} data.}
#'   \item{\code{bat}}{DWP associated with size class Bat}
#'   \item{\code{sml}}{DWP associated with size  class Small}
#'   \item{\code{med}}{DWP associated with size  class Medium}
#'   \item{\code{lrg}}{DWP associated with size  class Large}
#' }
#'
#' @section Carcass Observations (\code{CO}):
#' \code{$CO} is a data frame with a row for carcass observed in the carcass
#'   searches and a number of columns giving information about the given
#'   carcass (date found, size, species, etc.)
#' \describe{
#'   \item{\code{carcID}}{unique identifier for each carcass.}
#'   \item{\code{Unit}}{identifier for which unit the given carcass was found 
#'   at: \code{"arc19"}, \code{"arc65"}, etc, for arcs in the outer heliostat 
#'   field, or \code{"center"}, indicating the inner heliostat field.}
#'   \item{\code{Species}}{species of the carcass: \code{"BA"}, \code{"BB"},
#'     \code{"BC"}, \code{"BD"}, \code{"BE"}, \code{"LA"}, \code{"LB"},
#'     \code{"LD"}, \code{"LE"}, \code{"MA"}, \code{"MB"}, \code{"SA"},
#'     \code{"SB"}, \code{"SC"}, \code{"SD"}, \code{"SE"}, \code{"SF"},
#'     \code{"SG"}}
#'   \item{\code{Size}}{size: \code{"bat"}, \code{"lrg"}, \code{"med"},
#'    \code{"sml"}}
#'   \item{\code{Season}}{\code{"winter"}, \code{"spring"}, \code{"summer"}, or 
#'   \code{"fall"}}
#'   \item{\code{Flux}}{An optional field indicating whether there Was evidence 
#'   the animal was killed by flux. \code{"TRUE"}, or \code{"False"}.}
#'   \item{\code{Field}}{Optional indicator of whether the animal found in the 
#'   \code{"inner"} or \code{"outer"} heliostat field?}
#'   \item{\code{Ring}}{Optional note animals found in the outer heliostat field
#'    indicating which concentric ring the carcass was found in.}
#'   \item{\code{Distance}}{Optional note animals found in the outer heliostat 
#'   field representing the perpendicular distance from the searcher the carcass
#'    was discovered at.}
#'   \item{\code{DateFound}}{dates entered in the same format as in 
#'   \code{$SS$SearchDate}. Every date entered here is (and must be) included
#'      in the search schedule (\code{$SS$SearchDate}}
#'   \item{\code{X}}{Distance in meters from the Western edge of the facility.}
#'   \item{\code{Y}}{Distance in meters from the Southern edge of the facility.}
#' }
#'
#' @source \code{solar_powerTower}
"solar_powerTower"
