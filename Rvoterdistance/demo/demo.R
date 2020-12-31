# Distance from voters to drop boxes
data(king_dbox)

# Haversine distance between voter and drop boxes, King County
hav_calc <- nearest_dbox (king_geo$Residence_Addresses_Latitude, king_geo$Residence_Addresses_Longitude, dbox$lat, dbox$long)

summary(hav_calc)

# Read in Mecklenburg County data
data(meck_ev)

# Voter and early vote location, Mecklenburg County
hav_meck <- nearest_dbox (voter_meck$lat, voter_meck$long,
				early_meck$lat, early_meck$long)

summary(hav_meck)
head(hav_meck)

# Convert distance to miles
hav_mile <- dist_mile(num_vec=hav_meck, vec_only=TRUE)
head(hav_mile)



