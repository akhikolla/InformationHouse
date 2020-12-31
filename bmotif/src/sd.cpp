#include <Rcpp.h>
#include <math.h>       /* sqrt */

using namespace Rcpp;

//[[Rcpp::export]]
double sd_m2 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int x = 0; x < NZ; x++) {
		for (int p = 0; p < NP; p++) {
			for (int q = 0; q < NP; q++) {
				if ((W(x,p) * W(x,q) != 0) && (p != q)) {
					d = (W(x,p) + W(x,q)) / 2.0 - m;
					num +=  d * d;
				}
			}
		}
	}
	return(0.5 * num);
}

//[[Rcpp::export]]
double sd_m3 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int w = 0; w < NZ; w++) {
			for (int x = 0; x < NZ; x++) {
				if ((W(w,p) * W(x,p) != 0) && (w != x)) {
					d = (W(w,p) + W(x, p)) / 2.0 - m;
					num +=  d * d;
				}
			}
		}
	}
	return(0.5 * num);
}


//[[Rcpp::export]]
double sd_m4 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int w = 0; w < NZ; w++) {
			for (int x = 0; x < NZ; x++) {
				if (w == x) {
					continue;
				}
				for (int y = 0; y < NZ; y++) {
					if ((W(w,p) * W(x,p) * W(y,p) != 0) && (w != y) && (x != y)) {
						d = (W(w,p) + W(x, p) + W(y,p)) / 3.0 - m;
						num +=  d * d;
					}
				}
			}
		}
	}
	return(num / 6.0);
}


//[[Rcpp::export]]
double sd_m5 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
			if (p == q) {
				continue;
			}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if ((W(x,p) * W(x,q) * W(y,q) != 0) && (W(y,p) == 0) && (x!=y)) {
						d = (W(x, p) + W(x,q) + W(y,q)) / 3.0 - m;
						num +=  d * d;
					}
				}
			}
		}
	}
	return(num);
}

//[[Rcpp::export]]
double sd_m6 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
			if (p == q) {
				continue;
			}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if ((W(x,p) * W(x,q) * W(y,p) * W(y,q) != 0) && (x != y)) {
						d = (W(x, p) + W(x,q) + W(y,p) + W(y,q)) / 4.0 - m;
						num +=  d * d;
					}
				}
			}
		}
	}
	return(0.25 * num);
}


//[[Rcpp::export]]
double sd_m7 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int x = 0; x < NZ; x++) {
		for (int p = 0; p < NP; p++) {
			for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
				for (int r = 0; r < NP; r++) {
					if ((W(x,p) * W(x,q) * W(x,r) != 0) && (p != r) && (q != r)) {
						d = (W(x,p) + W(x,q) + W(x,r)) / 3.0 - m;
						num +=  d * d;
					}
				}
			}
		}
	}
	return(num / 6.0);
}


//[[Rcpp::export]]
double sd_m8 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int w = 0; w < NZ; w++) {
			for (int x = 0; x < NZ; x++) {
				if (w == x) {
					continue;
				}
				for (int y = 0; y < NZ; y++) {
					if (w == y || x == y) {
						continue;
					}
					for (int z = 0; z < NZ; z++) {
						if ((W(w,p) * W(x,p) * W(y,p) * W(z,p) != 0) && (w != z) && (y != z) && (x != z)) {
							d = 0.25 * (W(w,p) + W(x, p) + W(y,p) + W(z,p)) - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(num / 24.0);
}


//[[Rcpp::export]]
double sd_m9 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if (x == y) {
						continue;
					}
					for (int z = 0; z < NZ; z++) {
						if ((W(x,p) * W(y,p) * W(z,p) * W(z,q) != 0 && W(x,q) == 0) && (W(y,q) == 0) && (x != z) && (y != z)) {
							d = 0.25 * (W(x,p) + W(y,p) + W(z,p) + W(z,q)) - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(0.5 * num);
}

//[[Rcpp::export]]
double sd_m10 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if (x == y) {
						continue;
					}
					for (int z = 0; z < NZ; z++) {
						if ((W(x,p) * W(y,p) * W(y,q) * W(z,q) != 0) && (W(x,q) == 0) && (W(z,p) == 0) && (x != z) && (y != z)) {
							d = 0.25 * (W(x,p) + W(y,p) + W(y,q) + W(z,q)) - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(0.5 * num);
}

//[[Rcpp::export]]
double sd_m11 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if (x == y) {
						continue;
					}
					for (int z = 0; z < NZ; z++) {
						if ((W(x,p) * W(y,p) * W(y,q) * W(z,p) * W(z,q) != 0) && (W(x,q) == 0) && (x != z) && (y != z)) {
							d = (W(x,p) + W(y,p) + W(y,q) + W(z,p) + W(z,q)) / 5.0 - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(0.5 * num);
}

//[[Rcpp::export]]
double sd_m12 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if (x == y) {
						continue;
					}
					for (int z = 0; z < NZ; z++) {
						if ((W(x,p) * W(x,q) * W(y,p) * W(y,q) * W(z,p) * W(z,q) != 0) && (x != z) && (y != z)) {
							d = (W(x,p) + W(x,q) + W(y,p) + W(y,q) + W(z,p) + W(z,q)) / 6.0 - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(num / 12.0);
}

//[[Rcpp::export]]
double sd_m13 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
			for (int r = 0; r < NP; r++) {
				if (p == r || q == r) {
					continue;
				}				
				for (int x = 0; x < NZ; x++) {
					for (int y = 0; y < NZ; y++) {
						if ((W(x,p) * W(y,p) * W(y,q) * W(y,r) != 0) && (W(x,q) == 0) & (W(x,r) == 0) && (x != y)) {
							d = 0.25 * (W(x,p) + W(y,p) + W(y,q) + W(y,r)) - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(0.5 * num);
}

//[[Rcpp::export]]
double sd_m14 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
			for (int r = 0; r < NP; r++) {
				if (p == r || q == r) {
					continue;
				}				
				for (int x = 0; x < NZ; x++) {
					for (int y = 0; y < NZ; y++) {
						if ((W(x,p) * W(x,q) * W(y,q) * W(y,r) != 0) && (W(x,r) == 0) && (W(y,p) == 0) && (x != y)) {
							d = 0.25 * (W(x,p) + W(x,q) + W(y,q) + W(y,r)) - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(0.5 * num);
}

//[[Rcpp::export]]
double sd_m15 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
			for (int r = 0; r < NP; r++) {
				if (p == r || q == r) {
					continue;
				}				
				for (int x = 0; x < NZ; x++) {
					for (int y = 0; y < NZ; y++) {
						if ((W(x,p) * W(x,q) * W(y,p) * W(y,q) * W(y,r) != 0) && (W(x,r) == 0) && (x != y)) {
							d = (W(x,p) + W(x,q) + W(y,p) + W(y,q) + W(y,r)) / 5.0 - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(0.5 * num);
}

//[[Rcpp::export]]
double sd_m16 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
			for (int r = 0; r < NP; r++) {
				if (p == r || q == r) {
					continue;
				}				
				for (int x = 0; x < NZ; x++) {
					for (int y = 0; y < NZ; y++) {
						if ((W(x,p) * W(x,q) * W(x,r) * W(y,p) * W(y,q) * W(y,r) != 0) && (x != y)) {
							d = (W(x,p) + W(x,q) + W(x,r) + W(y,p) + W(y,q) + W(y,r)) / 6.0 - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(num / 12.0);
}


//[[Rcpp::export]]
double sd_m17 (int NZ, int NP, NumericMatrix W, double m) {
	double num = 0;
	double d = 0;
	for (int x = 0; x < NZ; x++) {
		for (int p = 0; p < NP; p++) {
			for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
				for (int r = 0; r < NP; r++) {
					if (p == r || q == r) {
						continue;
					}
					for (int s = 0; s < NP; s++) {
						if ((W(x,p) * W(x,q) * W(x,r) * W(x,s) != 0) && (p != s) && (q != s) && (r != s)) {
							d = 0.25 * (W(x,p) + W(x,q) + W(x,r) + W(x,s)) - m;
							num +=  d * d;
						}
					}
				}
			}
		}
	}
	return(num / 24.0);
}
