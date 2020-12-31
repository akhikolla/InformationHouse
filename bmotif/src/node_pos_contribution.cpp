#include <Rcpp.h>
#include <math.h>       /* sqrt */

using namespace Rcpp;

//[[Rcpp::export]]
List np_m2_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v3(NP);
	NumericVector v4(NZ);
	double d = 0;
	for (int x = 0; x < NZ; x++) {
		for (int p = 0; p < NP; p++) {
			for (int q = 0; q < NP; q++) {
				if ((W(x,p) * W(x,q) != 0) && (p != q)) {
					d = W(x,p) + W(x,q);
					v4[x] += 1;
					v3[p] += W(x,p) / d;
				}
			}
		}
	}
	v4 = 0.5 * v4;
	List L = List::create(v3, v4);
	return L;
}

//[[Rcpp::export]]
List np_m3_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v5(NP);
	NumericVector v6(NZ);
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int w = 0; w < NZ; w++) {
			for (int x = 0; x < NZ; x++) {
				if ((W(w,p) * W(x,p) != 0) && (w != x)) {
					d = W(w,p) + W(x, p);
					v5[p] += 1;
					v6[w] += W(w,p) / d;
				}
			}
		}
	}
	v5 = 0.5 * v5;
	List L = List::create(v5, v6);
	return L;
}

//[[Rcpp::export]]
List np_m4_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v7(NP);
	NumericVector v8(NZ);
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int w = 0; w < NZ; w++) {
			for (int x = 0; x < NZ; x++) {
				if (w == x) {
					continue;
				}
				for (int y = 0; y < NZ; y++) {
					if ((W(w,p) * W(x,p) * W(y,p) != 0) && (w != y) && (x != y)) {
						d = (W(w,p) + W(x, p) + W(y,p));
						v7[p] += 1;
						v8[w] += W(w,p) / d;
					}
				}
			}
		}
	}
	v7 = v7 / 6.0;
	v8 = 0.5 * v8;
	List L = List::create(v7, v8);
	return L;
}

//[[Rcpp::export]]
List np_m5_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v12(NZ); //this is automatically filled with zeros
	NumericVector v11(NZ);
	NumericVector v9(NP);
	NumericVector v10(NP);
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
			if (p == q) {
				continue;
			}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if ((W(x,p) * W(x,q) * W(y,q) != 0) && (W(y,p) == 0) && (x!=y)) {
						d = (W(x, p) + W(x,q) + W(y,q));
						v12[x] += (W(x,p) + W(x,q)) / d;
						v11[y] += W(y,q) / d;
						v10[q] += (W(x,q) + W(y,q)) / d;
						v9[p] += W(x,p) / d;
					}
				}
			}
		}
	}
	// now want to return a list with all vectors 
	List L = List::create(v9, v10, v11, v12);
	return L;
}

//[[Rcpp::export]]
List np_m6_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v13(NP);
	NumericVector v14(NZ);
	double d = 0;
	for (int p = 0; p < NP; p++) {
		for (int q = 0; q < NP; q++) {
			if (p == q) {
				continue;
			}
			for (int x = 0; x < NZ; x++) {
				for (int y = 0; y < NZ; y++) {
					if ((W(x,p) * W(x,q) * W(y,p) * W(y,q) != 0) && (x != y)) {
						d = (W(x, p) + W(x,q) + W(y,p) + W(y,q));
						v14[x] += (W(x,p) + W(x,q)) / d;
						v13[p] += (W(x,p) + W(y,p)) / d;
					}
				}
			}
		}
	}
	v13 = 0.5 * v13;
	v14 = 0.5 * v14;
	List L = List::create(v13, v14);
	return L;
}

//[[Rcpp::export]]
List np_m7_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v15(NP);
	NumericVector v16(NZ);
	double d = 0;	
	for (int x = 0; x < NZ; x++) {
		for (int p = 0; p < NP; p++) {
			for (int q = 0; q < NP; q++) {
				if (p == q) {
					continue;
				}
				for (int r = 0; r < NP; r++) {
					if ((W(x,p) * W(x,q) * W(x,r) != 0) && (p != r) && (q != r)) {
						d = (W(x,p) + W(x,q) + W(x,r));
						v15[p] += W(x,p) / d;
						v16[x] += 1;
					}
				}
			}
		}
	}
	v16 = v16 / 6.0;
	v15 = 0.5 * v15;
	List L = List::create(v15, v16);
	return L;
}

//[[Rcpp::export]]
List np_m8_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v17(NP);
	NumericVector v18(NZ);
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
							d = (W(w,p) + W(x, p) + W(y,p) + W(z,p));
							v17[p] += 1;
							v18[x] += W(x,p) / d;
						}
					}
				}
			}
		}
	}
	v17 = v17 / 24.0;
	v18 = v18 / 6.0;
	List L = List::create(v17, v18);
	return L;
}

//[[Rcpp::export]]
List np_m9_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v19(NP);
	NumericVector v20(NP);
	NumericVector v21(NZ);
	NumericVector v22(NZ);
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
							d = (W(x,p) + W(y,p) + W(z,p) + W(z,q));
							v20[p] += (W(x,p) + W(y,p) + W(z,p)) / d;
							v19[q] += W(z,q) / d;
							v21[x] += W(x,p) / d;
							v22[z] += (W(z,p) + W(z,q)) / d;
						}
					}
				}
			}
		}
	}
	v19 = 0.5 * v19;
	v20 = 0.5 * v20;
	v22 = 0.5 * v22;
	List L = List::create(v19, v20, v21, v22);
	return L;
}

//[[Rcpp::export]]
List np_m10_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v23(NP);
	NumericVector v24(NZ);
	NumericVector v25(NZ);
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
							d = (W(x,p) + W(y,p) + W(y,q) + W(z,q));
							v23[p] += (W(x,p) + W(y,p)) / d;
							v24[x] += W(x,p) / d;
							v25[y] += (W(y,p) + W(y,q)) / d;
						}
					}
				}
			}
		}
	}
	v25 = 0.5 * v25;	
	List L = List::create(v23, v24, v25);
	return L;
}

//[[Rcpp::export]]
List np_m11_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v26(NP);
	NumericVector v27(NP);
	NumericVector v28(NZ);
	NumericVector v29(NZ);
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
							d = (W(x,p) + W(y,p) + W(y,q) + W(z,p) + W(z,q));
							v28[x] += W(x,p) / d;
							v29[y] += (W(y,p) + W(y,q)) / d;
							v27[p] += (W(x,p) + W(y,p) + W(z,p)) / d;
							v26[q] += (W(y,q) + W(z,q)) / d;
						}
					}
				}
			}
		}
	}
	v26 = 0.5 * v26;
	v27 = 0.5 * v27;
	v28 = 0.5 * v28;
	List L = List::create(v26, v27, v28, v29);
	return L;
}

//[[Rcpp::export]]
List np_m12_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v30(NP);
	NumericVector v31(NZ);
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
							d = (W(x,p) + W(x,q) + W(y,p) + W(y,q) + W(z,p) + W(z,q));
							v30[p] += (W(x,p) + W(y,p) + W(z,p)) / d;
							v31[x] += (W(x,p) + W(x,q)) / d;
						}
					}
				}
			}
		}
	}
	v31 = 0.25 * v31;
	v30 = v30 / 6.0;
	List L = List::create(v30, v31);
	return L;
}

//[[Rcpp::export]]
List np_m13_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v32(NP);
	NumericVector v33(NP);
	NumericVector v34(NZ);
	NumericVector v35(NZ);
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
							d = (W(x,p) + W(y,p) + W(y,q) + W(y,r));
							v34[x] += W(x,p) / d;
							v35[y] += (W(y,p) + W(y,q) + W(y,r)) / d;
							v33[p] += (W(x,p) + W(y,p)) / d;
							v32[q] += W(y,q) / d;
						}
					}
				}
			}
		}
	}
	v34 = 0.5 * v34;
	v35 = 0.5 * v35;
	v33 = 0.5 * v33;
	List L = List::create(v32, v33, v34, v35);
	return L;
}

//[[Rcpp::export]]
List np_m14_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v36(NP);
	NumericVector v37(NP);
	NumericVector v38(NZ);
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
							d = (W(x,p) + W(x,q) + W(y,q) + W(y,r));
							v38[x] += (W(x,p) + W(x,q)) / d;
							v36[p] += W(x,p) / d;
							v37[q] += (W(x,q) + W(y,q)) / d;
						}
					}
				}
			}
		}
	}
	v37 = 0.5 * v37;
	List L = List::create(v36, v37, v38);
	return L;
}

//[[Rcpp::export]]
List np_m15_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v39(NP);
	NumericVector v40(NP);
	NumericVector v41(NZ);
	NumericVector v42(NZ);
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
							d = (W(x,p) + W(x,q) + W(y,p) + W(y,q) + W(y,r));
							v39[r] += W(y,r) / d;
							v40[p] += (W(x,p) + W(y,p)) / d;
							v41[x] += (W(x,p) + W(x,q)) / d;
							v42[y] += (W(y,p) + W(y,q) + W(y,r)) / d;
						}
					}
				}
			}
		}
	}
	v39 = 0.5 * v39;
	v41 = 0.5 * v41;
	v42 = 0.5 * v42;
	List L = List::create(v39, v40, v41, v42);
	return L;
}

//[[Rcpp::export]]
List np_m16_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v43(NP);
	NumericVector v44(NZ);
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
							d = (W(x,p) + W(x,q) + W(x,r) + W(y,p) + W(y,q) + W(y,r));
							v43[p] += (W(x,p) + W(y,p)) / d;
							v44[x] += (W(x,p) + W(x,q) + W(x,r)) / d;
						}
					}
				}
			}
		}
	}
	v43 = 0.25 * v43;
	v44 = v44 / 6.0;
	List L = List::create(v43, v44);
	return L;
}


//[[Rcpp::export]]
List np_m17_con (int NZ, int NP, NumericMatrix W) {
	NumericVector v45(NP);
	NumericVector v46(NZ);
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
							d = W(x,p) + W(x,q) + W(x,r) + W(x,s);
							v46[x] += 1;
							v45[p] += W(x,p) / d;
						}
					}
				}
			}
		}
	}
	v46 = v46 / 24.0;
	v45 = v45 / 6.0;
	List L = List::create(v45, v46);
	return L;
}
