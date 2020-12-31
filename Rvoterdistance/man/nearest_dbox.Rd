\name{nearest_dbox}
\alias{nearest_dbox}
\title{Calculates nearest drop box or polling location}
\description{Given a set of lat-long coordinates for each voter, and a set of coordinates for all drop boxes or polling locations, nearest_dbox() calculates the nearest drop box or polling location for each voter, in haversines. The function ports to C++, which greatly expedites speed.}
\usage{
nearest_dbox(lat1d_vec, lon1d_vec, lat2d_vec, lon2d_vec)
}
\arguments{
  \item{lat1d_vec}{Numeric vector, latitude coordinate of voter}

  \item{lon1d_vec}{Numeric vector, longitude coordinate of voter}

  \item{lat2d_vec}{Numeric vector, latitude coordinate of drop box, polling location}

  \item{lon2d_vec}{Numeric vector, longitude coordinate of drop box, polling location}

}
\value{
A numeric vector of minimum distances for each voter to their nearest drop box or polling location
}
\references{
Haversine: CC Robusto, 1957
}
\author{
Loren Collingwood <loren.collingwood@ucr.edu>
}
\examples{
data(king_dbox)
# Haversine distance between voter and drop boxes, King County
hav_calc <- nearest_dbox (king_geo$Residence_Addresses_Latitude, 
king_geo$Residence_Addresses_Longitude, 
dbox$lat, dbox$long)

summary(hav_calc)

data(meck_ev)
# Voter and early vote location, Mecklenburg County
hav_meck <- nearest_dbox (voter_meck$lat, voter_meck$long,
				early_meck$lat, early_meck$long)
summary(hav_meck)

}