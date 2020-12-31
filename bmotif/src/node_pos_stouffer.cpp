#include <Rcpp.h>
#include <math.h>       /* sqrt */

using namespace Rcpp;

/* A few comments:
NumericVectors are automatically filled with zeros
mw is just summing the mean weights of the entire motif
nw is summing the mean weights of the links that the node is in
con is summing the quotient of (sum of node's link weights) / (sum of all motif link weights)
*/

//[[Rcpp::export]]
List np_m2_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v3_py(NP);
	NumericVector v4_py(NZ);
	double d = 0;
	double py_q = 0;
	for (int x = 0; x < NZ; x++) {
		for (int p = 0; p < NP; p++) {
			for (int q = 0; q < NP; q++) {
				if ((W(x,p) * W(x,q) != 0) && (p != q)) {
					d = 0.5 * (W(x,p) + W(x,q));
					py_q = d / (W(x,p) + d + W(x,q));

					v4_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;
					v3_py[p] += W(x,p) * py_q;
				}
			}
		}
	}

	v4_py = 0.5 * v4_py;
	List L = List::create(v3_py, v4_py);
	return L;
}

//[[Rcpp::export]]
List np_m3_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v5_py(NP);
	NumericVector v6_py(NZ);
	double d = 0;
	double py_q = 0;
	for (int p = 0; p < NP; p++) {
		for (int w = 0; w < NZ; w++) {
			for (int x = 0; x < NZ; x++) {
				if ((W(w,p) * W(x,p) != 0) && (w != x)) {
					d = 0.5 * (W(w,p) + W(x, p));
					py_q = d / (W(x,p) + d + W(w,p));

					v5_py[p] += d * py_q;
					v6_py[x] += W(x,p) * py_q;
				}
			}
		}
	}
	v5_py = 0.5 * v5_py;
	List L = List::create(v5_py, v6_py);
	return L;
}


//[[Rcpp::export]]
List np_m4_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v7_py(NP);
	NumericVector v8_py(NZ);
	double d = 0;
	double py_q = 0;
	for (int p = 0; p < NP; p++) {
		for (int w = 0; w < NZ; w++) {
			for (int x = 0; x < NZ; x++) {
				if (w == x) {
					continue;
				}
				for (int y = 0; y < NZ; y++) {
					if ((W(w,p) * W(x,p) * W(y,p) != 0) && (w != y) && (x != y)) {
						d = (W(w,p) + W(x, p) + W(y,p)) / 3.0;
						py_q = d / (d + W(w,p) + W(x,p) + W(y,p));

						v7_py[p] += d * py_q;
						v8_py[w] += W(w,p) * py_q;
					}
				}
			}
		}
	}

	v7_py = v7_py / 6.0;
	v8_py = 0.5 * v8_py;
	List L = List::create(v7_py, v8_py);
	return L;
}

//[[Rcpp::export]]
List np_m5_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v12_py(NZ);
	NumericVector v11_py(NZ);
	NumericVector v9_py(NP);
	NumericVector v10_py(NP);
	double d = 0;
	double py_q = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
			if (p == q) {
				continue;
			}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if ((W(x,p) * W(x,q) * W(y,q) != 0) && (W(y,p) == 0) && (x!=y)) {
						d = (W(x, p) + W(x,q) + W(y,q)) / 3.0;
						py_q = d / (0.5 * (W(x,p) + W(x,q)) + W(y,q) + 0.5 * (W(x,q) + W(y,q)) + W(x,p));
						
						v12_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;
						v11_py[y] += W(y,q) * py_q;
						v10_py[q] += 0.5 * (W(x,q) + W(y,q)) * py_q;
						v9_py[p] += W(x,p) * py_q;						
					}
				}
			}
		}
	}
	// now want to return a list with all vectors 
	List L = List::create(v9_py, v10_py, v11_py, v12_py);
	return L;
}


//[[Rcpp::export]]
List np_m6_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v13_py(NP);
	NumericVector v14_py(NZ);
	double py_q = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
			if (p == q) {
				continue;
			}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if ((W(x,p) * W(x,q) * W(y,p) * W(y,q) != 0) && (x != y)) {
						py_q = 0.25;
						v14_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;
						v13_py[p] += 0.5 * (W(x,p) + W(y,p)) * py_q;
					}
				}
			}
		}
	}
	v13_py = 0.5 * v13_py;
	v14_py = 0.5 * v14_py;
	
	List L = List::create(v13_py, v14_py);
	return L;
}


//[[Rcpp::export]]
List np_m7_py (int NZ, int NP, NumericMatrix W) {

	NumericVector v15_py(NP);
	NumericVector v16_py(NZ);
	double d = 0;
	double py_q = 0;
	for (int x = 0; x < NZ; x++) {
		for (int p = 0; p < NP; p++) {
			for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
				for (int r = 0; r < NP; r++) {
					if ((W(x,p) * W(x,q) * W(x,r) != 0) && (p != r) && (q != r)) {
						d = (W(x,p) + W(x,q) + W(x,r)) / 3.0;
						py_q = d / (W(x,p) + d + W(x,q) + W(x,r));
						v15_py[p] +=  W(x,p) * py_q;
						v16_py[x] += d * py_q;
					}
				}
			}
		}
	}

	v16_py = v16_py / 6.0;
	v15_py = 0.5 * v15_py;	
	List L = List::create(v15_py, v16_py);
	return L;
}


//[[Rcpp::export]]
List np_m8_py (int NZ, int NP, NumericMatrix W) {

	NumericVector v17_py(NP);
	NumericVector v18_py(NZ);
	double d = 0;
	double py_q = 0;
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
							d = 0.25 * (W(w,p) + W(x, p) + W(y,p) + W(z,p));
							py_q = d / (W(x,p) + d + W(w,p) + W(y,p) + W(z,p)); //quotient in pymfinder

							v17_py[p] += d * py_q;
							v18_py[x] += W(x,p) * py_q;
						}
					}
				}
			}
		}
	}

	v17_py = v17_py / 24.0;
	v18_py = v18_py / 6.0;
	List L = List::create(v17_py, v18_py);
	return L;
}



//[[Rcpp::export]]
List np_m9_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v19_py(NP);
	NumericVector v20_py(NP);
	NumericVector v21_py(NZ);
	NumericVector v22_py(NZ);	
	double d = 0;
	double py_q = 0;
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
							d = 0.25 * (W(x,p) + W(y,p) + W(z,p) + W(z,q));
							py_q = d / ((W(x,p) + W(y,p) + W(z,p)) / 3.0 + W(z,q) + W(x,p) + W(y,p) + 0.5 * (W(z,p) + W(z,q)));

							v20_py[p] += (W(x,p) + W(y,p) + W(z,p)) * py_q / 3.0;
							v19_py[q] += W(z,q) * py_q;
							v21_py[x] += W(x,p) * py_q;
							v22_py[z] += 0.5 * (W(z,p) + W(z,q)) * py_q;							
						}
					}
				}
			}
		}
	}
	v19_py = 0.5 * v19_py;
	v20_py = 0.5 * v20_py;
	v22_py = 0.5 * v22_py;
	List L = List::create(v19_py, v20_py, v21_py, v22_py);
	return L;
}


//[[Rcpp::export]]
List np_m10_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v23_py(NP);
	NumericVector v24_py(NZ);
	NumericVector v25_py(NZ);
	
	double d = 0;
	double py_q = 0;
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
							d = 0.25 * (W(x,p) + W(y,p) + W(y,q) + W(z,q));
							py_q = d / (W(x,p) + 0.5 * (W(y,p) + W(y,q)) + W(z,q) + 0.5 * (W(x,p) + W(y,p)) + 0.5 * (W(y,q) + W(z,q)));
							
							v23_py[p] += 0.5 * (W(x,p) + W(y,p)) * py_q;
							v24_py[x] += W(x,p) * py_q;
							v25_py[y] += 0.5 * (W(y,p) + W(y,q)) * py_q;							
						}
					}
				}
			}
		}
	}

	v25_py = 0.5 * v25_py;
	List L = List::create(v23_py, v24_py, v25_py);
	return L;
}


//[[Rcpp::export]]
List np_m11_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v26_mw(NP);
	NumericVector v27_mw(NP);
	NumericVector v28_mw(NZ);
	NumericVector v29_mw(NZ);

	NumericVector v26_nw(NP);
	NumericVector v27_nw(NP);
	NumericVector v28_nw(NZ);
	NumericVector v29_nw(NZ);

	NumericVector v26_con(NP);
	NumericVector v27_con(NP);
	NumericVector v28_con(NZ);
	NumericVector v29_con(NZ);

	NumericVector v26_py(NP);
	NumericVector v27_py(NP);
	NumericVector v28_py(NZ);
	NumericVector v29_py(NZ);
	
	double d = 0;
	double py_q = 0;
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
							d = (W(x,p) + W(y,p) + W(y,q) + W(z,p) + W(z,q)) / 5.0;
							py_q = d / (W(x,p) + 0.5 * (W(y,p) + W(y,q)) + 0.5 * (W(z,p) + W(z,q)) + (W(x,p) + W(y,p) + W(z,p)) / 3.0 + 0.5 * (W(y,q) + W(z,q)));				
							
							v28_py[x] += W(x,p) * py_q;
							v29_py[y] += 0.5 * (W(y,p) + W(y,q)) * py_q;
							v27_py[p] += (W(x,p) + W(y,p) + W(z,p)) * py_q / 3.0;
							v26_py[q] += 0.5 * (W(y,q) + W(z,q)) * py_q;
						}
					}
				}
			}
		}
	}
	v26_py = 0.5 * v26_py;
	v27_py = 0.5 * v27_py;
	v28_py = 0.5 * v28_py;	
	List L = List::create(v26_py, v27_py, v28_py, v29_py);
	return L;
}


//[[Rcpp::export]]
List np_m12_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v30_py(NP);
	NumericVector v31_py(NZ);	
	double d = 0;
	double py_q = 0;
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
							d = (W(x,p) + W(x,q) + W(y,p) + W(y,q) + W(z,p) + W(z,q)) / 6.0;
							py_q = d /( (W(x,p) + W(y,p) + W(z,p)) / 3.0 + (W(x,q) + W(y,q) + W(z,q)) / 3.0 + 0.5 * (W(x,p) + W(x,q)) + 0.5 * (W(y,p) + W(y,q)) + 0.5 * (W(z,p) + W(z,q)));

							v30_py[p] += (W(x,p) + W(y,p) + W(z,p)) * py_q / 3.0;
							v31_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;							
						}
					}
				}
			}
		}
	}
	v31_py = 0.25 * v31_py;
	v30_py = v30_py / 6.0;
	List L = List::create(v30_py, v31_py);
	return L;
}


//[[Rcpp::export]]
List np_m13_py (int NZ, int NP, NumericMatrix W) {

	NumericVector v32_py(NP);
	NumericVector v33_py(NP);
	NumericVector v34_py(NZ);
	NumericVector v35_py(NZ);	
	double d = 0;
	double py_q = 0;
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
							d = 0.25 * (W(x,p) + W(y,p) + W(y,q) + W(y,r));
							py_q = d / (W(x,p) + (W(y,p) + W(y,q) + W(y,r)) / 3.0 + 0.5 * (W(x,p) + W(y,p)) + W(y,q) + W(y,r));
							v34_py[x] += W(x,p) * py_q;
							v35_py[y] += (W(y,p) + W(y,q) + W(y,r)) / 3.0 * py_q;
							v33_py[p] += 0.5 * (W(x,p) + W(y,p)) * py_q;
							v32_py[q] += W(y,q) * py_q;
							
						}
					}
				}
			}
		}
	}	
	v34_py = 0.5 * v34_py;
	v35_py = 0.5 * v35_py;
	v33_py = 0.5 * v33_py;	
	List L = List::create(v32_py, v33_py, v34_py, v35_py);
	return L;
}


//[[Rcpp::export]]
List np_m14_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v36_py(NP);
	NumericVector v37_py(NP);
	NumericVector v38_py(NZ);
	double d = 0;
	double py_q = 0;
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
							d = 0.25 * (W(x,p) + W(x,q) + W(y,q) + W(y,r));
							py_q = d / (0.5 * (W(x,p) + W(x,q)) + 0.5 * (W(y,q) + W(y,r)) + W(x,p) + W(y,r) + 0.5 * (W(x,q) + W(y,q)));
							v38_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;
							v36_py[p] += W(x,p) * py_q;
							v37_py[q] += 0.5 * (W(x,q) + W(y,q)) * py_q;							
						}
					}
				}
			}
		}
	}
	v37_py = 0.5 * v37_py;
	List L = List::create(v36_py, v37_py, v38_py);
	return L;
}



//[[Rcpp::export]]
List np_m15_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v39_py(NP);
	NumericVector v40_py(NP);
	NumericVector v41_py(NZ);
	NumericVector v42_py(NZ);	
	double d = 0;
	double py_q = 0;
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
							d = (W(x,p) + W(x,q) + W(y,p) + W(y,q) + W(y,r)) / 5.0;
							py_q = d / (W(y,r) + 0.5 * (W(x,p) + W(y,p)) + 0.5 * (W(x,q) + W(y,q)) + 0.5 * (W(x,p) + W(x,q)) + (W(y,p) + W(y,q) + W(y,r)) / 3.0);

							v39_py[r] += W(y,r) * py_q;
							v40_py[p] += 0.5 * (W(x,p) + W(y,p)) * py_q;
							v41_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;
							v42_py[y] += (W(y,p) + W(y,q) + W(y,r)) * py_q / 3.0;							
						}
					}
				}
			}
		}
	}
	v39_py = 0.5 * v39_py;
	v41_py = 0.5 * v41_py;
	v42_py = 0.5 * v42_py;	
	List L = List::create(v39_py, v40_py, v41_py, v42_py);
	return L;
}


//[[Rcpp::export]]
List np_m16_py (int NZ, int NP, NumericMatrix W) {

	NumericVector v43_py(NP);
	NumericVector v44_py(NZ);
	
	double d = 0;
	double py_q = 0;
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
							d = (W(x,p) + W(x,q) + W(x,r) + W(y,p) + W(y,q) + W(y,r)) / 6.0;
							py_q = d / (0.5 * (W(x,p) + W(y,p)) + 0.5 * (W(x,q) + W(y,q)) + 0.5 * (W(x,r) + W(y,r)) + (W(x,p) + W(x,q) + W(x,r)) / 3.0 + (W(y,p) + W(y,q) + W(y,r)) / 3.0);

							v43_py[p] += 0.5 * (W(x,p) + W(y,p)) * py_q;
							v44_py[x] += (W(x,p) + W(x,q) + W(x,r)) * py_q / 3.0; 
						}
					}
				}
			}
		}
	}

	v43_py = 0.25 * v43_py;
	v44_py = v44_py / 6.0;	
	List L = List::create(v43_py, v44_py);
	return L;
}


//[[Rcpp::export]]
List np_m17_py (int NZ, int NP, NumericMatrix W) {
	NumericVector v45_py(NP);
	NumericVector v46_py(NZ);	
	double d = 0;
	double py_q = 0;
	
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
							d = 0.25 * (W(x,p) + W(x,q) + W(x,r) + W(x,s));
							py_q = d / (d + W(x,p) + W(x,q) + W(x,r) + W(x,s));

							v46_py[x] += d* py_q;
							v45_py[p] += W(x,p) * py_q;							
						}
					}
				}
			}
		}
	}
	v46_py = v46_py / 24.0;
	v45_py = v45_py / 6.0;	
	List L = List::create(v45_py, v46_py);
	return L;
}
