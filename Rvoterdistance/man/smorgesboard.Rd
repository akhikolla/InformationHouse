\name{smorgesboard}
\alias{smorgesboard}
\title{Calculates nearest drop box or polling location, Smorgesboard back}
\description{Given a set of lat-long coordinates for each voter, and a set of coordinates for all drop boxes or polling locations; this function returns a dataframe length data1 (usually voter file), including haversine, mile, and kilometer distance output, as well as any other data2 variables (perhaps address).}
\usage{
smorgesboard(data1, data2 , lat_long1_char, lat_long2_char)
}
\arguments{
  \item{data1}{Dataset, probably a voter file, including lat/long coordinates}

  \item{data2}{Dataset, probably of drop box locations/polling locations, including lat/long coordinates}

  \item{lat_long1_char}{Character vector, latitude/longitude column names found in data1}

  \item{lat_long2_char}{Character vector, latitude/longitude column names found in data2}
  
}
\value{
A data frame of length data1, with all columns from data2 and distance_haversine, distance_mile, and distance_km appended.
}
\references{
Haversine: CC Robusto, 1957
}
\author{
Loren Collingwood <loren.collingwood@ucr.edu>
}
\seealso{dist_km, dist_mile, nearest_dbox}
\examples{

data(meck_ev)
str(voter_meck)
str(early_meck)

# Grab smorgesboard of distance information from polling location
vote_distance <- smorgesboard(voter_meck, early_meck[,-1], c("lat", "long"), c("lat", "long"))
head(vote_distance)

}