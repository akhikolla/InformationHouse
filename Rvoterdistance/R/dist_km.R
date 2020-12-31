dist_km <- function(lat1d_vec, lon1d_vec, lat2d_vec, lon2d_vec, num_vec=NULL, vec_only=FALSE) {

	if (vec_only) {
		km <- num_vec / 1000
	} else {
		hav <- nearest_dbox(lat1d_vec, lon1d_vec, lat2d_vec, lon2d_vec)
		km <- hav / 1000
	}
	return(km)
}
