dist_mile <- function(lat1d_vec, lon1d_vec, lat2d_vec, lon2d_vec, num_vec=NULL, vec_only=FALSE) {
	if (vec_only) {
		mile <- num_vec / 1609.34
	} else {
		hav <- nearest_dbox(lat1d_vec, lon1d_vec, lat2d_vec, lon2d_vec)
		mile <- hav / 1609.34
	}
	return(mile)
}
