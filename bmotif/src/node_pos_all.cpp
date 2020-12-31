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
List np_m2_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v3_mw(NP);
	NumericVector v4_mw(NZ);
	NumericVector v3_nw(NP);
	NumericVector v4_nw(NZ);
	NumericVector v3_con(NP);
	NumericVector v4_con(NZ);
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
					v4_mw[x] += d;
					v3_mw[p] += d;
					v4_nw[x] += 0.5 * (W(x,p) + W(x,q));
					v3_nw[p] += W(x,p);
					v4_con[x] += 1;
					v3_con[p] += 0.5 * W(x,p) / d;
					v4_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;
					v3_py[p] += W(x,p) * py_q;
				}
			}
		}
	}
	v4_mw = 0.5 * v4_mw;
	v4_nw = 0.5 * v4_nw;
	v4_con = 0.5 * v4_con;
	v4_py = 0.5 * v4_py;
	List L = List::create(v3_mw, v3_nw, v3_con, v3_py, v4_mw, v4_nw, v4_con, v4_py);
	return L;
}

//[[Rcpp::export]]
List np_m3_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v5_mw(NP);
	NumericVector v6_mw(NZ);
	NumericVector v5_nw(NP);
	NumericVector v6_nw(NZ);
	NumericVector v5_con(NP);
	NumericVector v6_con(NZ);
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
					v5_mw[p] += d;
					v6_mw[x] += d;
					v5_nw[p] += d;
					v6_nw[x] += W(x,p);
					v5_con[p] += 1;
					v6_con[x] += 0.5 * W(x,p) / d;
					v5_py[p] += d * py_q;
					v6_py[x] += W(x,p) * py_q;
				}
			}
		}
	}
	v5_mw = 0.5 * v5_mw;
	v5_nw = 0.5 * v5_nw;
	v5_con = 0.5 * v5_con;
	v5_py = 0.5 * v5_py;
	List L = List::create(v5_mw, v5_nw, v5_con, v5_py, v6_mw, v6_nw, v6_con, v6_py);
	return L;
}


//[[Rcpp::export]]
List np_m4_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v7_mw(NP);
	NumericVector v8_mw(NZ);
	NumericVector v7_nw(NP);
	NumericVector v8_nw(NZ);
	NumericVector v7_con(NP);
	NumericVector v8_con(NZ);
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
						v7_mw[p] += d;
						v8_mw[w] += d;
						v7_nw[p] += d;
						v8_nw[w] += W(w,p);
						v7_con[p] += 1;
						v8_con[w] += W(w,p) / d / 3.0;
						v7_py[p] += d * py_q;
						v8_py[w] += W(w,p) * py_q;
					}
				}
			}
		}
	}
	v7_mw = v7_mw / 6.0;
	v8_mw = 0.5 * v8_mw;
	v7_nw = v7_nw / 6.0;
	v8_nw = 0.5 * v8_nw;
	v7_con = v7_con / 6.0;
	v8_con = 0.5 * v8_con;
	v7_py = v7_py / 6.0;
	v8_py = 0.5 * v8_py;
	List L = List::create(v7_mw, v7_nw, v7_con, v7_py, v8_mw, v8_nw, v8_con, v8_py);
	return L;
}

//[[Rcpp::export]]
List np_m5_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v12_mw(NZ);
	NumericVector v11_mw(NZ);
	NumericVector v9_mw(NP);
	NumericVector v10_mw(NP);
	
	NumericVector v12_nw(NZ);
	NumericVector v11_nw(NZ);
	NumericVector v9_nw(NP);
	NumericVector v10_nw(NP);
	
	NumericVector v12_con(NZ);
	NumericVector v11_con(NZ);
	NumericVector v9_con(NP);
	NumericVector v10_con(NP);
	
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
						v12_mw[x] += d;
						v11_mw[y] += d;
						v10_mw[q] += d;
						v9_mw[p] += d;
						
						v12_nw[x] += 0.5 * (W(x,p) + W(x,q));
						v11_nw[y] += W(y,q);
						v10_nw[q] += 0.5 * (W(x,q) + W(y,q));
						v9_nw[p] += W(x,p);
						
						v12_con[x] += (W(x,p) + W(x,q)) / d / 3.0;
						v11_con[y] += W(y,q) / d / 3.0;
						v10_con[q] += (W(x,q) + W(y,q)) / d / 3.0;
						v9_con[p] += W(x,p) / d / 3.0;
						
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
	List L = List::create(v9_mw, v9_nw, v9_con, v9_py, v10_mw, v10_nw, v10_con, v10_py, v11_mw, v11_nw, v11_con, v11_py, v12_mw, v12_nw, v12_con, v12_py);
	return L;
}


//[[Rcpp::export]]
List np_m6_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v13_mw(NP);
	NumericVector v14_mw(NZ);
	NumericVector v13_nw(NP);
	NumericVector v14_nw(NZ);
	NumericVector v13_con(NP);
	NumericVector v14_con(NZ);
	NumericVector v13_py(NP);
	NumericVector v14_py(NZ);
	double d = 0;
	double py_q = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
			if (p == q) {
				continue;
			}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if ((W(x,p) * W(x,q) * W(y,p) * W(y,q) != 0) && (x != y)) {
						d = (W(x, p) + W(x,q) + W(y,p) + W(y,q)) / 4.0;
						py_q = 0.25;
						v14_mw[x] += d;
						v13_mw[p] += d;
						v14_nw[x] += 0.5 * (W(x,p) + W(x,q));
						v13_nw[p] += 0.5 * (W(x,p) + W(y,p));
						v14_con[x] += 0.25 * (W(x,p) + W(x,q)) / d;
						v13_con[p] += 0.25 * (W(x,p) + W(y,p)) / d;
						v14_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;
						v13_py[p] += 0.5 * (W(x,p) + W(y,p)) * py_q;
					}
				}
			}
		}
	}
	v13_mw = 0.5 * v13_mw;
	v14_mw = 0.5 * v14_mw;
	v13_nw = 0.5 * v13_nw;
	v14_nw = 0.5 * v14_nw;
	v13_con = 0.5 * v13_con;
	v14_con = 0.5 * v14_con;
	v13_py = 0.5 * v13_py;
	v14_py = 0.5 * v14_py;
	
	List L = List::create(v13_mw, v13_nw, v13_con, v13_py, v14_mw, v14_nw, v14_con, v14_py);
	return L;
}


//[[Rcpp::export]]
List np_m7_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v15_mw(NP);
	NumericVector v16_mw(NZ);
	NumericVector v15_nw(NP);
	NumericVector v16_nw(NZ);
	NumericVector v15_con(NP);
	NumericVector v16_con(NZ);
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
						v15_mw[p] += d;
						v16_mw[x] += d;
						v15_nw[p] += W(x,p);
						v16_nw[x] += d;
						v15_con[p] +=  W(x,p) / d / 3.0;
						v16_con[x] += 1;
						v15_py[p] +=  W(x,p) * py_q;
						v16_py[x] += d * py_q;
					}
				}
			}
		}
	}
	v16_mw = v16_mw / 6.0;
	v15_mw = 0.5 * v15_mw;
	v16_nw = v16_nw / 6.0;
	v15_nw = 0.5 * v15_nw;
	v16_con = v16_con / 6.0;
	v15_con = 0.5 * v15_con;
	v16_py = v16_py / 6.0;
	v15_py = 0.5 * v15_py;	
	List L = List::create(v15_mw, v15_nw, v15_con, v15_py, v16_mw, v16_nw, v16_con, v16_py);
	return L;
}


//[[Rcpp::export]]
List np_m8_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v17_mw(NP);
	NumericVector v18_mw(NZ);
	NumericVector v17_nw(NP);
	NumericVector v18_nw(NZ);
	NumericVector v17_con(NP);
	NumericVector v18_con(NZ);
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
							v17_mw[p] += d;
							v18_mw[x] += d;
							v17_nw[p] += d;
							v18_nw[x] += W(x,p);
							v17_con[p] += 1;
							v18_con[x] += 0.25 * W(x,p) / d;
							v17_py[p] += d * py_q;
							v18_py[x] += W(x,p) * py_q;
						}
					}
				}
			}
		}
	}
	v17_mw = v17_mw / 24.0;
	v18_mw = v18_mw / 6.0;
	v17_nw = v17_nw / 24.0;
	v18_nw = v18_nw / 6.0;
	v17_con = v17_con / 24.0;
	v18_con = v18_con / 6.0;
	v17_py = v17_py / 24.0;
	v18_py = v18_py / 6.0;
	List L = List::create(v17_mw, v17_nw, v17_con, v17_py, v18_mw, v18_nw, v18_con, v18_py);
	return L;
}



//[[Rcpp::export]]
List np_m9_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v19_mw(NP);
	NumericVector v20_mw(NP);
	NumericVector v21_mw(NZ);
	NumericVector v22_mw(NZ);
	
	NumericVector v19_nw(NP);
	NumericVector v20_nw(NP);
	NumericVector v21_nw(NZ);
	NumericVector v22_nw(NZ);

	NumericVector v19_con(NP);
	NumericVector v20_con(NP);
	NumericVector v21_con(NZ);
	NumericVector v22_con(NZ);
	
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
							v20_mw[p] += d;
							v19_mw[q] += d;
							v21_mw[x] += d;
							v22_mw[z] += d;

							v20_nw[p] += (W(x,p) + W(y,p) + W(z,p)) / 3.0;
							v19_nw[q] += W(z,q);
							v21_nw[x] += W(x,p);
							v22_nw[z] += 0.5 * (W(z,p) + W(z,q));

							v20_con[p] += 0.25 * (W(x,p) + W(y,p) + W(z,p)) / d;
							v19_con[q] += 0.25 * W(z,q) / d;
							v21_con[x] += 0.25 * W(x,p) / d;
							v22_con[z] += 0.25 * (W(z,p) + W(z,q)) / d;

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
	v19_mw = 0.5 * v19_mw;
	v20_mw = 0.5 * v20_mw;
	v22_mw = 0.5 * v22_mw;
	
	v19_nw = 0.5 * v19_nw;
	v20_nw = 0.5 * v20_nw;
	v22_nw = 0.5 * v22_nw;
	
	v19_con = 0.5 * v19_con;
	v20_con = 0.5 * v20_con;
	v22_con = 0.5 * v22_con;
	
	v19_py = 0.5 * v19_py;
	v20_py = 0.5 * v20_py;
	v22_py = 0.5 * v22_py;
	List L = List::create(v19_mw, v19_nw, v19_con, v19_py, v20_mw, v20_nw, v20_con, v20_py, v21_mw, v21_nw, v21_con, v21_py, v22_mw, v22_nw, v22_con, v22_py);
	return L;
}


//[[Rcpp::export]]
List np_m10_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v23_mw(NP);
	NumericVector v24_mw(NZ);
	NumericVector v25_mw(NZ);
	
	NumericVector v23_nw(NP);
	NumericVector v24_nw(NZ);
	NumericVector v25_nw(NZ);	

	NumericVector v23_con(NP);
	NumericVector v24_con(NZ);
	NumericVector v25_con(NZ);

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
							v23_mw[p] += d;
							v24_mw[x] += d;
							v25_mw[y] += d;
							
							v23_nw[p] += 0.5 * (W(x,p) + W(y,p));
							v24_nw[x] += W(x,p);
							v25_nw[y] += 0.5 * (W(y,p) + W(y,q));							

							v23_con[p] += 0.25 * (W(x,p) + W(y,p)) / d;
							v24_con[x] += 0.25 * W(x,p) / d;
							v25_con[y] += 0.25 * (W(y,p) + W(y,q)) / d;

							v23_py[p] += 0.5 * (W(x,p) + W(y,p)) * py_q;
							v24_py[x] += W(x,p) * py_q;
							v25_py[y] += 0.5 * (W(y,p) + W(y,q)) * py_q;							
						}
					}
				}
			}
		}
	}
	v25_mw = 0.5 * v25_mw;
	v25_nw = 0.5 * v25_nw;
	v25_con = 0.5 * v25_con;
	v25_py = 0.5 * v25_py;
	List L = List::create(v23_mw, v23_nw, v23_con, v23_py, v24_mw, v24_nw, v24_con, v24_py, v25_mw, v25_nw, v25_con, v25_py);
	return L;
}


//[[Rcpp::export]]
List np_m11_c (int NZ, int NP, NumericMatrix W) {
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
							v28_mw[x] += d;
							v29_mw[y] += d;
							v27_mw[p] += d;
							v26_mw[q] += d;
							
							v28_nw[x] += W(x,p);
							v29_nw[y] += 0.5 * (W(y,p) + W(y,q));
							v27_nw[p] += (W(x,p) + W(y,p) + W(z,p)) / 3.0;
							v26_nw[q] += 0.5 * (W(y,q) + W(z,q));
							
							v28_con[x] += 0.2 * W(x,p) / d;
							v29_con[y] += 0.2 * (W(y,p) + W(y,q)) / d;
							v27_con[p] += 0.2 * (W(x,p) + W(y,p) + W(z,p)) / d;
							v26_con[q] += 0.2 * (W(y,q) + W(z,q)) / d;							
							
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
	v26_mw = 0.5 * v26_mw;
	v27_mw = 0.5 * v27_mw;
	v28_mw = 0.5 * v28_mw;
	
	v26_nw = 0.5 * v26_nw;
	v27_nw = 0.5 * v27_nw;
	v28_nw = 0.5 * v28_nw;

	v26_con = 0.5 * v26_con;
	v27_con = 0.5 * v27_con;
	v28_con = 0.5 * v28_con;

	v26_py = 0.5 * v26_py;
	v27_py = 0.5 * v27_py;
	v28_py = 0.5 * v28_py;	
	List L = List::create(v26_mw, v26_nw, v26_con, v26_py, v27_mw, v27_nw, v27_con, v27_py, v28_mw, v28_nw, v28_con, v28_py, v29_mw, v29_nw, v29_con, v29_py);
	return L;
}


//[[Rcpp::export]]
List np_m12_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v30_mw(NP);
	NumericVector v31_mw(NZ);
	NumericVector v30_nw(NP);
	NumericVector v31_nw(NZ);
	NumericVector v30_con(NP);
	NumericVector v31_con(NZ);
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
							v30_mw[p] += d;
							v31_mw[x] += d;
							
							v30_nw[p] += (W(x,p) + W(y,p) + W(z,p)) / 3.0;
							v31_nw[x] += 0.5 * (W(x,p) + W(x,q));
							
							v30_con[p] += (W(x,p) + W(y,p) + W(z,p)) / d / 6.0;
							v31_con[x] += (W(x,p) + W(x,q)) / d / 6.0;

							v30_py[p] += (W(x,p) + W(y,p) + W(z,p)) * py_q / 3.0;
							v31_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;							
						}
					}
				}
			}
		}
	}
	v31_mw = 0.25 * v31_mw;
	v30_mw = v30_mw / 6.0;
	v31_nw = 0.25 * v31_nw;
	v30_nw = v30_nw / 6.0;
	v31_con = 0.25 * v31_con;
	v30_con = v30_con / 6.0;
	v31_py = 0.25 * v31_py;
	v30_py = v30_py / 6.0;
	List L = List::create(v30_mw, v30_nw, v30_con, v30_py, v31_mw, v31_nw, v31_con, v31_py);
	return L;
}


//[[Rcpp::export]]
List np_m13_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v32_mw(NP);
	NumericVector v33_mw(NP);
	NumericVector v34_mw(NZ);
	NumericVector v35_mw(NZ);
	
	NumericVector v32_nw(NP);
	NumericVector v33_nw(NP);
	NumericVector v34_nw(NZ);
	NumericVector v35_nw(NZ);
	
	NumericVector v32_con(NP);
	NumericVector v33_con(NP);
	NumericVector v34_con(NZ);
	NumericVector v35_con(NZ);
	
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
							v34_mw[x] += d;
							v35_mw[y] += d;
							v33_mw[p] += d;
							v32_mw[q] += d;
							
							v34_nw[x] += W(x,p);
							v35_nw[y] += (W(y,p) + W(y,q) + W(y,r)) / 3.0;
							v33_nw[p] += 0.5 * (W(x,p) + W(y,p));
							v32_nw[q] += W(y,q);							
							
							v34_con[x] += 0.25 * W(x,p) / d;
							v35_con[y] += 0.25 * (W(y,p) + W(y,q) + W(y,r)) / d;
							v33_con[p] += 0.25 * (W(x,p) + W(y,p)) / d;
							v32_con[q] += 0.25 * W(y,q) / d;
							
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
	v34_mw = 0.5 * v34_mw;
	v35_mw = 0.5 * v35_mw;
	v33_mw = 0.5 * v33_mw;
	
	v34_nw = 0.5 * v34_nw;
	v35_nw = 0.5 * v35_nw;
	v33_nw = 0.5 * v33_nw;
	
	v34_con = 0.5 * v34_con;
	v35_con = 0.5 * v35_con;
	v33_con = 0.5 * v33_con;
	
	v34_py = 0.5 * v34_py;
	v35_py = 0.5 * v35_py;
	v33_py = 0.5 * v33_py;	
	List L = List::create(v32_mw, v32_nw, v32_con, v32_py, v33_mw, v33_nw, v33_con, v33_py, v34_mw, v34_nw, v34_con, v34_py, v35_mw, v35_nw, v35_con, v35_py);
	return L;
}


//[[Rcpp::export]]
List np_m14_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v36_mw(NP);
	NumericVector v37_mw(NP);
	NumericVector v38_mw(NZ);

	NumericVector v36_nw(NP);
	NumericVector v37_nw(NP);
	NumericVector v38_nw(NZ);

	NumericVector v36_con(NP);
	NumericVector v37_con(NP);
	NumericVector v38_con(NZ);

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
							v38_mw[x] += d;
							v36_mw[p] += d;
							v37_mw[q] += d;
							
							v38_nw[x] += 0.5 * (W(x,p) + W(x,q));
							v36_nw[p] += W(x,p);
							v37_nw[q] += 0.5 * (W(x,q) + W(y,q));
							
							v38_con[x] += 0.25 * (W(x,p) + W(x,q)) / d;
							v36_con[p] += 0.25 * W(x,p) / d;
							v37_con[q] += 0.25 * (W(x,q) + W(y,q)) / d;

							v38_py[x] += 0.5 * (W(x,p) + W(x,q)) * py_q;
							v36_py[p] += W(x,p) * py_q;
							v37_py[q] += 0.5 * (W(x,q) + W(y,q)) * py_q;							
						}
					}
				}
			}
		}
	}
	v37_mw = 0.5 * v37_mw;
	v37_nw = 0.5 * v37_nw;
	v37_con = 0.5 * v37_con;
	v37_py = 0.5 * v37_py;
	List L = List::create(v36_mw, v36_nw, v36_con, v36_py, v37_mw, v37_nw, v37_con, v37_py, v38_mw, v38_nw, v38_con, v38_py);
	return L;
}



//[[Rcpp::export]]
List np_m15_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v39_mw(NP);
	NumericVector v40_mw(NP);
	NumericVector v41_mw(NZ);
	NumericVector v42_mw(NZ);
	
	NumericVector v39_nw(NP);
	NumericVector v40_nw(NP);
	NumericVector v41_nw(NZ);
	NumericVector v42_nw(NZ);

	NumericVector v39_con(NP);
	NumericVector v40_con(NP);
	NumericVector v41_con(NZ);
	NumericVector v42_con(NZ);

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
							v39_mw[r] += d;
							v40_mw[p] += d;
							v41_mw[x] += d;
							v42_mw[y] += d;
							
							v39_nw[r] += W(y,r);
							v40_nw[p] += 0.5 * (W(x,p) + W(y,p));
							v41_nw[x] += 0.5 * (W(x,p) + W(x,q));
							v42_nw[y] += (W(y,p) + W(y,q) + W(y,r)) / 3.0;
							
							v39_con[r] += W(y,r) / d / 5.0;
							v40_con[p] += (W(x,p) + W(y,p)) / d / 5.0;
							v41_con[x] += (W(x,p) + W(x,q)) / d / 5.0;
							v42_con[y] += (W(y,p) + W(y,q) + W(y,r)) / d / 5.0;

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
	v39_mw = 0.5 * v39_mw;
	v41_mw = 0.5 * v41_mw;
	v42_mw = 0.5 * v42_mw;
	
	v39_nw = 0.5 * v39_nw;
	v41_nw = 0.5 * v41_nw;
	v42_nw = 0.5 * v42_nw;

	v39_con = 0.5 * v39_con;
	v41_con = 0.5 * v41_con;
	v42_con = 0.5 * v42_con;

	v39_py = 0.5 * v39_py;
	v41_py = 0.5 * v41_py;
	v42_py = 0.5 * v42_py;	
	List L = List::create(v39_mw, v39_nw, v39_con, v39_py, v40_mw, v40_nw, v40_con, v40_py, v41_mw, v41_nw, v41_con, v41_py, v42_mw, v42_nw, v42_con, v42_py);
	return L;
}


//[[Rcpp::export]]
List np_m16_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v43_mw(NP);
	NumericVector v44_mw(NZ);
	
	NumericVector v43_nw(NP);
	NumericVector v44_nw(NZ);

	NumericVector v43_con(NP);
	NumericVector v44_con(NZ);	
	
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
							v43_mw[p] += d;
							v44_mw[x] += d;
							
							v43_nw[p] += 0.5 * (W(x,p) + W(y,p));
							v44_nw[x] += (W(x,p) + W(x,q) + W(x,r)) / 3.0;
							
							v43_con[p] += (W(x,p) + W(y,p)) / d / 6.0;
							v44_con[x] += (W(x,p) + W(x,q) + W(x,r)) / d / 6.0;
							
							v43_py[p] += 0.5 * (W(x,p) + W(y,p)) * py_q;
							v44_py[x] += (W(x,p) + W(x,q) + W(x,r)) * py_q / 3.0; 
						}
					}
				}
			}
		}
	}
	v43_mw = 0.25 * v43_mw;
	v44_mw = v44_mw / 6.0;
	v43_nw = 0.25 * v43_nw;
	v44_nw = v44_nw / 6.0;
	v43_con = 0.25 * v43_con;
	v44_con = v44_con / 6.0;
	v43_py = 0.25 * v43_py;
	v44_py = v44_py / 6.0;	
	List L = List::create(v43_mw, v43_nw, v43_con, v43_py, v44_mw, v44_nw, v44_con, v44_py);
	return L;
}


//[[Rcpp::export]]
List np_m17_c (int NZ, int NP, NumericMatrix W) {
	NumericVector v45_mw(NP);
	NumericVector v46_mw(NZ);
	
	NumericVector v45_nw(NP);
	NumericVector v46_nw(NZ);
	
	NumericVector v45_con(NP);
	NumericVector v46_con(NZ);

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
							v46_mw[x] += d;
							v45_mw[p] += d;
							
							v46_nw[x] += 0.25 * (W(x,p) + W(x,q) + W(x,r) + W(x,s));
							v45_nw[p] += W(x,p);

							v46_con[x] += 1;
							v45_con[p] += 0.25 * W(x,p) / d;

							v46_py[x] += d* py_q;
							v45_py[p] += W(x,p) * py_q;							
						}
					}
				}
			}
		}
	}
	v46_mw = v46_mw / 24.0;
	v45_mw = v45_mw / 6.0;
	
	v46_nw = v46_nw / 24.0;
	v45_nw = v45_nw / 6.0;

	v46_con = v46_con / 24.0;
	v45_con = v45_con / 6.0;

	v46_py = v46_py / 24.0;
	v45_py = v45_py / 6.0;	
	List L = List::create(v45_mw, v45_nw, v45_con, v45_py, v46_mw, v46_nw, v46_con, v46_py);
	return L;
}
