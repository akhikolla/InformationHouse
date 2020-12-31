#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double c_s_h(double t, double h){
  return (1 / (1 + exp(-t / h)));
}

// [[Rcpp::export]]
/* theta_sh_0 */
NumericVector c_dscrt(NumericMatrix y, NumericVector z, NumericVector l){
  IntegerVector ydim(y.attr("dim"));
  int yr = ydim[0];
  int yc = ydim[1];

  int n = z.size(); /* find the number of subjects */

  int i, j;
  double temp;
  NumericVector temp_y(n);

  for(i = 0; i < yr; i++){
    temp = 0;
    for(j = 0; j < yc; j++){
      temp += y[i + yr * j] * l[j];
    }
    temp_y[i] = temp;
  }

  double psij, temp_psij;
  NumericVector psii(n);

  for(i = 0; i < n; i++){
    psij = 0;
    for(j = 0; j < n; j++){
      if(i != j){
        if(((temp_y[i] - temp_y[j]) * (z[i] - z[j])) > 0){
          temp_psij = 1;
        }
        else if((temp_y[i] - temp_y[j]) == 0 || (z[i] - z[j]) == 0){
          temp_psij = 0.5;
        }
        else{
          temp_psij = 0;
        }
      psij += temp_psij;
      }
    }
    psii[i] = psij;
  }

  NumericVector theta_h_p(2);
  /* formula 9 */
  NumericVector s_psii(n);
  for(i = 0; i < n; i++){
    s_psii[i] = psii[i];
  }

  theta_h_p[0] = sum(s_psii) / (n * (n - 1));

  /*formula 10 & 11 */
  NumericVector temp_psii(n);
  for(i = 0; i < n; i++){
    temp_psii[i] = pow((psii[i] / (n - 1)) - theta_h_p[0], 2);
  }
  NumericVector s_temp_psii(n);
  for(i = 0; i < n; i++){
    s_temp_psii[i] = temp_psii[i];
  }
  theta_h_p[1] = sum(s_temp_psii) / ((n / 2) * ((n / 2) - 1));

return theta_h_p;
}  

// [[Rcpp::export]]
/*  theta_sh */
NumericVector c_cntin(NumericMatrix y, NumericVector z, NumericVector l, double h){
  IntegerVector ydim(y.attr("dim"));
//  int yr = ydim[0];
  int yc = ydim[1];
  
  int n = z.size(); /* find the number of subjects */

  int i, j, k;

  double psij, temp_y;
  NumericVector psii(n);
  double temp_psij, a, b;

  for(i = 0; i < n; i++){
    psij = 0;
    for(j = 0; j < n; j++){
      if(i != j){
        temp_y = 0;
        for(k = 0; k < yc; k++){
          temp_y += (y[i + k * n] - y[j + k * n]) * l[k];
        }
        a = c_s_h(temp_y, h);
        b = c_s_h(z[i] - z[j], h);
        /* formula 2.5 */
          temp_psij = a * b + (1 - a) * (1 - b);
        /* s.h(t(y[i, ] - y[j, ]) %*% l, h) * (2 * s.h((z[i] - z[j]), h) - 1) + 1 - s.h((z[i] - z[j]), h) */
          psij += temp_psij;
      }
    }
    psii[i] = psij;
  }

  NumericVector theta_sh_h_p(2);
  /* formula 9 */
  NumericVector s_psii(n);
  for(i = 0; i < n; i++){
    s_psii[i] = psii[i];
  }
  theta_sh_h_p[0] = sum(s_psii) / (n * (n - 1));

  /*formula 10 & 11 */
  NumericVector temp_psii(n);
  for(i = 0; i < n; ++i){
    temp_psii[i] = pow((psii[i] / (n - 1)) - theta_sh_h_p[0], 2);
  }
  NumericVector s_temp_psii(n);
  for(i = 0; i < n; i++){
    s_temp_psii[i] = temp_psii[i];
  }
  theta_sh_h_p[1] = sum(s_temp_psii) / ((n / 2) * ((n / 2) - 1));

return theta_sh_h_p;
}

// [[Rcpp::export]]
/* step 1 */
NumericVector c_d_theta_sh_h_p(NumericMatrix y, NumericVector z, NumericVector l, double h){
  IntegerVector ydim(y.attr("dim"));
//  int yr = ydim[0];
  int yc = ydim[1];

  int n = z.size(); /* the number of subjects */

  int i, j, k;
  double temp_y, a, b;
  NumericVector temp_sij(yc), sij(yc);
  /* formula 2.6 */
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(i != j){
				temp_y = 0;
				for(k = 0; k < yc; k++){
					temp_y += (y[i + k * n] - y[j + k * n]) * l[k];
				}
				a = c_s_h(temp_y, h);
				b = c_s_h(z[i] - z[j], h);
				for(k = 0; k < yc; k++){
					temp_sij[k] = (a * (1 - a) * ((2 * b) - 1) * (y[i + k * n] - y[j + k * n])) / h;
					sij[k] += temp_sij[k];
				}
			}
		}
	}
  
  NumericVector d_theta_sh_h_p(yc);
	for(i = 0; i < yc; i++){
		d_theta_sh_h_p[i] = sij[i] / (n * (n - 1));
	}
return d_theta_sh_h_p;
}

//double c_optimal_delta(NumericMatrix y, NumericVector z, NumericVector l, double h, NumericVector ind_d_l){
//  int i, j;
//  double l_i[50 * l.size()];
//  for(i = 0; i < l.size(); i++){
//    for(j = 0; j < 50; j++){
//      l_i[i * 50 + j] = l[i];
//    }
//  }
//  	
//	double delta[50];
//  for(i = 0; i < 50; i++){
//    delta[i] = (5 / 49) * i;
//  }
//
//  double m[50 * ind_d_l.size()];
//  for(i = 0; i < ind_d_l.size(); i++){
//    for(j = 0; j < 50; j++){
//      m[i * 50 + j] = delta[j] * ind_d_l[i];
//    }
//  }
//
//  for(i = 0; i < (50 * l.size()); i++){
//    l_i[i] = l_i[i] + m[i];
//  }
//  
//  double l_i_max[50], temp;
//  for(i = 0; i < 50; i++){
//    temp = l_i[i];
//    for(j = 0; j < ind_d_l.size(); j++){
//      if(temp < l_i[i + j * 50]){
//        temp = l_i[i + j * 50];
//      }
//    }
//    l_i_max[i] = temp;
//  }
//
//  for(i = 0; i < ind_d_l.size(); i++){
//    for(j = 0; j < 50; j++){
//      l_i[i * 50 + j] = l_i[i * 50 + j] / l_i_max[i];
//    }
//  }
//
//  double max_theta = 0; /* AUC = theta */
//  NumericVector temp_cntin;
//  NumericVector temp_l_i(l.size());
//  double delta_star = 3;
//  for(i = 1; i < 50; i++){
//    for(j = 0; j < l.size(); j++){
//      temp_l_i[i] = l_i[i + j * 50];
//    }
//    temp_cntin = c_cntin(y, z, temp_l_i, h);
//    if(max_theta < temp_cntin[0]){
//      max_theta = temp_cntin[0];
//      delta_star = delta[i];
//    }
//  }
//
//return delta_star;
//}
