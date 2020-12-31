#include <Rcpp.h>
using namespace Rcpp;

//////// Gaussian kernel

NumericVector erf_cpp(NumericVector x){
  return(2. * pnorm(x * sqrt(2.)) - 1.);
}


// [[Rcpp::export]]
NumericVector mi_gauss_cpp(NumericMatrix Mu, NumericVector sigma){
  NumericVector mis(Mu.nrow(), 1.);

  for(int i = 0; i < Mu.nrow(); i++){
    for(int j = 0; j < Mu.ncol(); j++){
      mis(i) *= 0.5 * sqrt(M_PI) * sigma(j) * (erf((1 - Mu(i,j))/sigma(j)) + erf(Mu(i,j)/sigma(j)));
    }
  }
  return(mis);
}

// -1/2 * sqrt(pi/2.) * sigma * exp(-(mu1 - mu2)^2 / (2*sigma^2)) * (erf((mu1 + mu2 - 2) / (sqrt(2.) * sigma)) - erf((mu1 + mu2) / (sqrt(2.) * sigma))) 

// // [[Rcpp::export]]
// NumericMatrix Mijs_gauss_cpp(NumericMatrix Mu1, NumericMatrix Mu2, NumericVector sigma){
//   NumericMatrix Mijs(Mu1.nrow(), Mu2.nrow());
//   Mijs.fill(1.);
//   
//   for(int i = 0; i < Mu1.nrow(); i++){
//     for(int j = 0; j < Mu2.nrow(); j++){
//       for(int k = 0; k < Mu1.ncol(); k++){
//         Mijs(i,j) *= -1./2. * sqrt(M_PI/2.) * sigma(k) * exp(-(Mu1(i,k) - Mu2(j,k)) * (Mu1(i,k) - Mu2(j,k)) / (2. * sigma(k) * sigma(k))) * (erf((Mu1(i, k) + Mu2(j, k) - 2.) / (sqrt(2.) * sigma(k))) - erf((Mu1(i, k) + Mu2(j, k)) / (sqrt(2.) * sigma(k)))) ;
//       }
//     }
//   }
//   
//   NumericVector mis1 = mi_gauss_cpp(Mu1, sigma);
//   NumericVector mis2 = mi_gauss_cpp(Mu2, sigma);
//   
//   for(int i = 0; i < Mu1.nrow(); i++){
//     for(int j = 0; j < Mu2.nrow(); j++){
//       Mijs(i,j) -= mis1(i) * mis2(j);
//     }
//   }
//   
//   
//   return(Mijs);
// }

// [[Rcpp::export]]
NumericMatrix Wijs_gauss_cpp(NumericMatrix Mu1, NumericMatrix Mu2, NumericVector sigma){
  int m1c = Mu1.ncol();
  int m2r = Mu2.nrow();
  NumericMatrix Wijs(Mu1.nrow(), m2r);
  Wijs.fill(1.);
  double a,b;
  
  for(int i = 0; i < Mu1.nrow(); i++){
    for(int j = 0; j < m2r; j++){
      const double* ptr_s = (const double*) &sigma(0);
      for(int k = 0; k < m1c; k++, ptr_s++){
        a = Mu1(i, k);
        b = Mu2(j, k);
        Wijs(i,j) *= -1./2. * sqrt(M_PI/2.) * *ptr_s * exp(-(a - b) * (a - b) / (2. * *ptr_s * *ptr_s)) * (erf((a + b - 2.) / (sqrt(2.) * *ptr_s)) - erf((a + b) / (sqrt(2.) * *ptr_s))) ;
      }
    }
  }
  
  return(Wijs);
}


// // [[Rcpp::export]]
// NumericMatrix Mijs_gauss_sym_cpp(NumericMatrix Mu, NumericVector sigma){
//   NumericMatrix Mijs(Mu.nrow(), Mu.nrow());
//   Mijs.fill(1.);
//   
//   for(int i = 0; i < Mu.nrow(); i++){
//     for(int j = 0; j <= i; j++){
//       for(int k = 0; k < Mu.ncol(); k++){
//         if(i == j){
//           Mijs(i, i) *= -1./2. * sqrt(M_PI/2.) * sigma(k) * (erf((2 * Mu(i, k)  - 2.) / (sqrt(2.) * sigma(k))) - erf((2 * Mu(i, k)) / (sqrt(2.) * sigma(k))));
//         }else{
//           Mijs(i,j) = Mijs(j,i) *= -1./2. * sqrt(M_PI/2.) * sigma(k) * exp(-(Mu(i,k) - Mu(j,k)) * (Mu(i,k) - Mu(j,k)) / (2. * sigma(k) * sigma(k))) * (erf((Mu(i, k) + Mu(j, k) - 2.) / (sqrt(2.) * sigma(k))) - erf((Mu(i, k) + Mu(j, k)) / (sqrt(2.) * sigma(k))));
//         }
//       }
//     }
//   }
//   
//   NumericVector mis = mi_gauss_cpp(Mu, sigma);
//   
//   for(int i = 0; i < Mu.nrow(); i++){
//     for(int j = 0; j <= i; j++){
//       if(i == j){
//         Mijs(i,i) -= mis(i) * mis(i);
//       }else{
//         Mijs(i,j) = Mijs(j,i) -= mis(i) * mis(j);
//       }
//       
//     }
//   }
//   
//   
//   return(Mijs);
// }


// [[Rcpp::export]]
NumericMatrix Wijs_gauss_sym_cpp(NumericMatrix Mu, NumericVector sigma){
  int mc = Mu.ncol();
  NumericMatrix Wijs(Mu.nrow(), Mu.nrow());
  Wijs.fill(1.);
  
  double a,b;
  
  for(int i = 0; i < Mu.nrow(); i++){
    for(int j = 0; j <= i; j++){
      const double* ptr_s = (const double*) &sigma(0);
      for(int k = 0; k < mc; k++, ptr_s++){
        if(i == j){
          a = Mu(i, k);
          Wijs(i, i) *= -1./2. * sqrt(M_PI/2.) * *ptr_s * (erf((2. * a  - 2.) / (sqrt(2.) * *ptr_s)) - erf((2. * a) / (sqrt(2.) * *ptr_s)));
        }else{
          a = Mu(i, k);
          b = Mu(j, k);
          Wijs(i,j) = Wijs(j,i) *= -1./2. * sqrt(M_PI/2.) * *ptr_s * exp(-(a - b) * (a - b) / (2. * *ptr_s * *ptr_s)) * (erf((a + b - 2.) / (sqrt(2.) * *ptr_s)) - erf((a + b) / (sqrt(2.) * *ptr_s)));
        }
      }
    }
  }
  
  return(Wijs);
}

// // // Jiangeng function
// // 
// // [[Rcpp::export]]
// NumericVector d_gauss_cpp(NumericMatrix X, NumericVector x, int m, NumericVector theta){
//   NumericVector dis(X.nrow());
//   NumericVector pro(X.nrow(), 1.);
// 
//   for(int i = 0; i < X.nrow(); i++){
//     for(int j = 0; j < X.ncol(); j++){
//       pro(i) *= exp(-(X(i, j) - x(j)) * (X(i, j) - x(j)) / theta(j));
//     }
//     dis(i) = 2 / theta(m) * (X(i, m) - x(m)) * pro(i);
//   }
//   return(dis);
// }


// Alternative: not recomputing the covariance

// [[Rcpp::export]]
NumericVector d_gauss_cpp(NumericVector X, double x, double sigma){
  NumericVector dis(X.length());  
  
  for(int i = 0; i < X.length(); i++){
    dis(i) = 2. / sigma * (X(i) - x);
  }
  return(dis);
}



double c1i_gauss(double x1, double X, double sigma){
  double tmp = -1./2. * sqrt(M_PI/2.) * sigma * exp(-(X - x1) * (X - x1) / (2. * sigma * sigma)) * (erf((X + x1 - 2.) / (sqrt(2.) * sigma)) - erf((X + x1) / (sqrt(2.) * sigma)));
  if(tmp == 0.) return(0.);
  
  return((0.5*exp(-(x1 - X)*(x1 - X) / (2.*sigma*sigma))* (exp(-(x1 + X)*(x1 + X)/(2.*sigma*sigma)) -
         exp(-(2.-(x1 + X))*(2.-(x1 + X))/(2.*sigma * sigma))) - sqrt(2.*M_PI)/4./sigma*(x1 - X) * exp(-(x1 - X)*(x1 - X)/(2.*sigma*sigma)) * 
         (erf((x1 + X)/(sqrt(2.)*sigma)) - erf((x1 + X - 2.)/(sqrt(2.)*sigma))))/tmp);
}

// [[Rcpp::export]]
double c2_gauss_cpp(double x, double t, double w){
  if(w == 0.) return(0.);
  double tmp = -1./2. * sqrt(M_PI/2.) * t * (erf((2. * x  - 2.) / (sqrt(2.) * t)) - erf((2. * x) / (sqrt(2.) * t)));
  if(tmp == 0.) return(0.);
  return((exp(-2. * x * x / (t * t)) - exp(-2.*(1. - x) * (1. - x) / (t * t)))*w/tmp);
}


// [[Rcpp::export]]
NumericVector c1_gauss_cpp(NumericVector X, double x, double sigma, NumericVector W){
  NumericVector cis(X.length());  
  
  for(int i = 0; i < X.length(); i++){
    cis(i) = c1i_gauss(x, X(i), sigma) * W(i);
  }
  
  return(cis);
}

//////// Matern 5/2 kernel

double A_2_cpp(double x){
  return((8. + 5. * x + x * x) * exp(-x));
}

// [[Rcpp::export]]
NumericVector mi_mat52_cpp(NumericMatrix Mu, NumericVector sigma){
  NumericVector mis(Mu.nrow(), 1.);

  for(int i = 0; i < Mu.nrow(); i++){
    for(int j = 0; j < Mu.ncol(); j++){
      mis(i) *= sigma(j)/(3.*sqrt(5.)) * (16. - A_2_cpp(sqrt(5.) * Mu(i,j)/sigma(j)) - A_2_cpp(sqrt(5.) * (1. - Mu(i,j))/sigma(j)));
    }
  }
  return(mis);
}


// Mij_mat52 <- function(a, b, t){

//* // [[Rcpp::export]]
// NumericMatrix Mijs_mat52_cpp(NumericMatrix Mu1, NumericMatrix Mu2, NumericVector sigma){
//   
//   NumericMatrix Mijs(Mu1.nrow(), Mu2.nrow());
//   Mijs.fill(1.);
//   
//   double a,b,t;
//   double p11, p12, p21, p22, p3, p41, p42;
//   
//   for(int i = 0; i < Mu1.nrow(); i++){
//     for(int j = 0; j < Mu2.nrow(); j++){
//       for(int k = 0; k < Mu1.ncol(); k++){
//         
//         if(Mu2(j,k) < Mu1(i,k)){
//           a = Mu2(j,k);
//           b = Mu1(i,k);
//         }else{
//           a = Mu1(i,k);
//           b = Mu2(j,k);
//         }
//         t = sigma(k);
//         p11 = -(-7*sqrt(5.)*pow(t,3) - 2*(10*b+13*a)*pow(t,2)-2*sqrt(5.)*a*(8*b+3*a)*t-20*a*a*b) * exp(sqrt(5.)*(-b-a)/t)/(24*t*t);
//         p12 = -(20*pow(a,3)-2*(11*sqrt(5.)*t+10*b+20*a)*a*a+2*(23*t*t+8*sqrt(5.)*b*t+14*sqrt(5.)*a*t+20*a*b+10*a*a)*a-7*sqrt(5.)*t*t*t-2*(10*b+13*a)*t*t-2*sqrt(5.)*a*(8*b+3*a)*t-20*a*a*b)*exp(sqrt(5.)*(a-b)/t)/(24*t*t);
//         p21 = -(525*pow(t,4)+(6*pow(5, 7/2.)*b+12*pow(5,5/2.)*a)*pow(t,3)+(2500*b*b+2500*a*b+250*a*a)*t*t+(16*pow(5,7/2.)*a*b*b+4*pow(5, 7/2.)*a*a*b)*t+2500*a*a*b*b)*exp(-sqrt(5.)*b/t-sqrt(5.)*a/t)/(72*pow(5,5/2.)*pow(t,3));
//         p22 = (2500*pow(a,4)+(-4*pow(5, 9/2.)*t-5000*b-5000*a)*pow(a,3)+(5250*t*t+(36*pow(5, 7/2.)*b+24*pow(5,7/2.)*a)*t+2500*b*b+10000*a*b+2500*a*a)*a*a+(-42*pow(5, 5/2.)*pow(t,3)+(-7500*b-3000*a)*t*t+(-16*pow(5, 7/2.)*b*b-8*pow(5,9/2.)*a*b-4*pow(5,7/2.)*a*a)*t-5000*a*b*b-5000*a*a*b)*a+525*pow(t,4)+(6*pow(5,7/2.)*b+12*pow(5,5/2.)*a)*pow(t,3)+(2500*b*b+2500*a*b+250*a*a)*t*t+(16*pow(5,7/2.)*a*b*b+4*pow(5,7/2.)*a*a*b)*t+2500*a*a*b*b)*exp(sqrt(5.)*a/t-sqrt(5.)*b/t)/(72*pow(5,5/2.)*pow(t,3));
//         p3 = (b-a)*(54*pow(t,4)+(54*sqrt(5.)*b-54*sqrt(5.)*a)*pow(t,3)+(105*b*b-210*a*b+105*a*a)*t*t+(3*pow(5,3/2.)*pow(b,3)-9*pow(5,3/2.)*a*b*b+9*pow(5,3/2.)*a*a*b-3*pow(5,3/2.)*pow(a,3))*t+5*pow(b,4)-20*a*pow(b,3)+30*a*a*b*b-20*pow(a,3)*b+5*pow(a,4))*exp(sqrt(5.)*(a-b)/t)/(54*pow(t,4));
//         p41 = -(t*(t*(9*t*(7*t-pow(5,3/2.)*(b+a-2))+10*b*(5*b+17*a-27)+10*(5*a*a-27*a+27))-8*pow(5,3/2.)*(a-1)*(b-1)*(b+a-2))+50*(a-1)*(a-1)*(b-2)*b+50*(a-1)*(a-1))*exp(-sqrt(5.)*(-b-a+2)/t)/(36*sqrt(5.)*pow(t,3));
//         p42 = t*t*(63*t*t+9*pow(5,3/2.)*b*t-9*pow(5,3/2.)*a*t+50*b*b-100*a*b+50*a*a)*exp(-sqrt(5.)*(b-a)/t)/(36*sqrt(5.)*pow(t,3));
//         
//         Mijs(i,j) *= -p11 + p12 + p21 + p22 + p3 + p41 + p42;
//         
//       }
//     }
//   }
//   
//   NumericVector mis1 = mi_mat52_cpp(Mu1, sigma);
//   NumericVector mis2 = mi_mat52_cpp(Mu2, sigma);
//   
//   for(int i = 0; i < Mu1.nrow(); i++){
//     for(int j = 0; j < Mu2.nrow(); j++){
//       Mijs(i,j) -= mis1(i) * mis2(j);
//     }
//   }
//   
//   return(Mijs);   
// }

// [[Rcpp::export]]
NumericMatrix Wijs_mat52_cpp(NumericMatrix Mu1, NumericMatrix Mu2, NumericVector sigma){
  
  int m1c = Mu1.ncol();
  int m2r = Mu2.nrow();
  double tmp;
  NumericMatrix Wijs(Mu1.nrow(), m2r);
  Wijs.fill(1.);
  
  double a,a2,b,b2,t,t2;
  double p1, p3, p4;
  
  for(int i = 0; i < Mu1.nrow(); i++){
    for(int j = 0; j < m2r; j++){
      const double* ptr_s = (const double*) &sigma(0);
      for(int k = 0; k < m1c; k++, ptr_s++){
        a = Mu1(i,k);
        b = Mu2(j,k);
        
        if(b < a){
          tmp = b; 
          b = a;
          a = tmp;
        }
        
        t = *ptr_s;
        t2 = t*t;
        a2 = a*a;
        b2 = b*b;
        p1 = (2*t2*(63*t2 + 9*5.*sqrt(5.)*b*t-9*5.*sqrt(5.)*a*t+50*b2-100*a*b+50*a2)*exp(2*sqrt(5.)*a/t)-63*t2*t2-9*5.*sqrt(5.)*(b+a)*t*t2-10*(5*b2+17*a*b+5*a2)*t2-8*5.*sqrt(5.)*a*b*(b+a)*t-50*a2*b2)*exp(-sqrt(5.)*(b+a)/t)/(36*sqrt(5.)*t*t2);
        p3 = (b-a)*(54*t2*t2+(54*sqrt(5.)*b-54*sqrt(5.)*a)*t*t2+(105*b2-210*a*b+105*a2)*t2+(3*5.*sqrt(5.)*b2*b-9*5.*sqrt(5.)*a*b2+9*5.*sqrt(5.)*a2*b-3*5.*sqrt(5.)*a2*a)*t+5*b2*b2-20*a*b2*b+30*a2*b2-20*a2*a*b+5*a2*a2)*exp(sqrt(5.)*(a-b)/t)/(54*t2*t2);
        p4 = -((t*(t*(9*t*(7*t-5.*sqrt(5.)*(b+a-2))+10*b*(5*b+17*a-27)+10*(5*a2-27*a+27))-8*5.*sqrt(5.)*(a-1)*(b-1)*(b+a-2))+50*(a-1)*(a-1)*(b-2)*b+50*(a-1)*(a-1))*exp(2*sqrt(5.)*b/t))*exp(-sqrt(5.)*(b-a+2)/t)/(36*sqrt(5.)*t*t2);
        
        Wijs(i,j) *= p1 + p3 + p4;
        
      }
    }
  } 
  
  return(Wijs);   
}


//* // [[Rcpp::export]]
// NumericMatrix Mijs_mat52_sym_cpp(NumericMatrix Mu, NumericVector sigma){
//   NumericMatrix Mijs(Mu.nrow(), Mu.nrow());
//   Mijs.fill(1.);
//   
//   double a,b,t;
//   double p11, p12, p21, p22, p3, p41, p42;
//   
//   for(int i = 0; i < Mu.nrow(); i++){
//     for(int j = 0; j <= i; j++){
//       for(int k = 0; k < Mu.ncol(); k++){
//         
//         if(Mu(j,k) < Mu(i,k)){
//           a = Mu(j,k);
//           b = Mu(i,k);
//         }else{
//           a = Mu(i,k);
//           b = Mu(j,k);
//         }
//         t = sigma(k);
//         p11 = -(-7*sqrt(5.)*pow(t,3) - 2*(10*b+13*a)*pow(t,2)-2*sqrt(5.)*a*(8*b+3*a)*t-20*a*a*b) * exp(sqrt(5.)*(-b-a)/t)/(24*t*t);
//         p12 = -(20*pow(a,3)-2*(11*sqrt(5.)*t+10*b+20*a)*a*a+2*(23*t*t+8*sqrt(5.)*b*t+14*sqrt(5.)*a*t+20*a*b+10*a*a)*a-7*sqrt(5.)*t*t*t-2*(10*b+13*a)*t*t-2*sqrt(5.)*a*(8*b+3*a)*t-20*a*a*b)*exp(sqrt(5.)*(a-b)/t)/(24*t*t);
//         p21 = -(525*pow(t,4)+(6*pow(5, 7/2.)*b+12*pow(5,5/2.)*a)*pow(t,3)+(2500*b*b+2500*a*b+250*a*a)*t*t+(16*pow(5,7/2.)*a*b*b+4*pow(5, 7/2.)*a*a*b)*t+2500*a*a*b*b)*exp(-sqrt(5.)*b/t-sqrt(5.)*a/t)/(72*pow(5,5/2.)*pow(t,3));
//         p22 = (2500*pow(a,4)+(-4*pow(5, 9/2.)*t-5000*b-5000*a)*pow(a,3)+(5250*t*t+(36*pow(5, 7/2.)*b+24*pow(5,7/2.)*a)*t+2500*b*b+10000*a*b+2500*a*a)*a*a+(-42*pow(5, 5/2.)*pow(t,3)+(-7500*b-3000*a)*t*t+(-16*pow(5, 7/2.)*b*b-8*pow(5,9/2.)*a*b-4*pow(5,7/2.)*a*a)*t-5000*a*b*b-5000*a*a*b)*a+525*pow(t,4)+(6*pow(5,7/2.)*b+12*pow(5,5/2.)*a)*pow(t,3)+(2500*b*b+2500*a*b+250*a*a)*t*t+(16*pow(5,7/2.)*a*b*b+4*pow(5,7/2.)*a*a*b)*t+2500*a*a*b*b)*exp(sqrt(5.)*a/t-sqrt(5.)*b/t)/(72*pow(5,5/2.)*pow(t,3));
//         p3 = (b-a)*(54*pow(t,4)+(54*sqrt(5.)*b-54*sqrt(5.)*a)*pow(t,3)+(105*b*b-210*a*b+105*a*a)*t*t+(3*pow(5,3/2.)*pow(b,3)-9*pow(5,3/2.)*a*b*b+9*pow(5,3/2.)*a*a*b-3*pow(5,3/2.)*pow(a,3))*t+5*pow(b,4)-20*a*pow(b,3)+30*a*a*b*b-20*pow(a,3)*b+5*pow(a,4))*exp(sqrt(5.)*(a-b)/t)/(54*pow(t,4));
//         p41 = -(t*(t*(9*t*(7*t-pow(5,3/2.)*(b+a-2))+10*b*(5*b+17*a-27)+10*(5*a*a-27*a+27))-8*pow(5,3/2.)*(a-1)*(b-1)*(b+a-2))+50*(a-1)*(a-1)*(b-2)*b+50*(a-1)*(a-1))*exp(-sqrt(5.)*(-b-a+2)/t)/(36*sqrt(5.)*pow(t,3));
//         p42 = t*t*(63*t*t+9*pow(5,3/2.)*b*t-9*pow(5,3/2.)*a*t+50*b*b-100*a*b+50*a*a)*exp(-sqrt(5.)*(b-a)/t)/(36*sqrt(5.)*pow(t,3));
//         
//         if(i == j){
//           Mijs(i,j) *= -p11 + p12 + p21 + p22 + p3 + p41 + p42;
//         }else{
//           Mijs(j,i) = Mijs(i,j) *= -p11 + p12 + p21 + p22 + p3 + p41 + p42;
//         }
//         
//       }
//     }
//   }
//   
//   NumericVector mis = mi_mat52_cpp(Mu, sigma);
//   
//   
//   for(int i = 0; i < Mu.nrow(); i++){
//     for(int j = 0; j <= i; j++){
//       if(i == j){
//         Mijs(j, i) -= mis(i) * mis(j);
//       }else{
//         Mijs(i,j) = Mijs(j, i) -= mis(i) * mis(j);
//       }
//     }
//   }
//   
//   return(Mijs);   
// }

// [[Rcpp::export]]
NumericMatrix Wijs_mat52_sym_cpp(NumericMatrix Mu, NumericVector sigma){
  int m1c = Mu.ncol();
  NumericMatrix Wijs(Mu.nrow(), Mu.nrow());
  Wijs.fill(1.);
  
  double a,a2,b,b2,t,t2;
  double p1, p3, p4;
  double tmp;
  
  for(int i = 0; i < Mu.nrow(); i++){
    for(int j = 0; j <= i; j++){
      const double* ptr_s = (const double*) &sigma(0);
      for(int k = 0; k < m1c; k++, ptr_s++){
        a = Mu(j,k);
        b = Mu(i,k);
        if(b < a){
          tmp = b;
          b = a;
          a = tmp;
        }
        
        t = *ptr_s;
        t2 = t*t;
        a2 = a*a;
        b2 = b*b;
        
        if(i == j){
          Wijs(i,j) *= (exp(-2*sqrt(5.)*a/t)*(63*t2*t2*exp(2*sqrt(5.)*a/t)-50*a2*a2-16*5*sqrt(5.)*t*a2*a-270*t2*a2-18*5*sqrt(5.)*t2*t*a-63*t2*t2)-exp(-2*sqrt(5.)/t)*((t*(t*(10*(5*a2-27*a+27)+9*t*(7*t-5*sqrt(5.)*(2*a-2))+10*a*(22*a-27))-8*5*sqrt(5.)*(a-1)*(a-1)*(2*a-2))+50*(a-2)*(a-1)*(a-1)*a+50*(a-1)*(a-1))*exp(2*sqrt(5.)*a/t)-63*t2*t2*exp(2*sqrt(5.)/t)))/(36*sqrt(5.)*t2*t);;
        }else{
          p1 = (2*t2*(63*t2 + 9*5.*sqrt(5.)*b*t-9*5.*sqrt(5.)*a*t+50*b2-100*a*b+50*a2)*exp(2*sqrt(5.)*a/t)-63*t2*t2-9*5.*sqrt(5.)*(b+a)*t*t2-10*(5*b2+17*a*b+5*a2)*t2-8*5.*sqrt(5.)*a*b*(b+a)*t-50*a2*b2)*exp(-sqrt(5.)*(b+a)/t)/(36*sqrt(5.)*t*t2);
          p3 = (b-a)*(54*t2*t2+(54*sqrt(5.)*b-54*sqrt(5.)*a)*t*t2+(105*b2-210*a*b+105*a2)*t2+(3*5.*sqrt(5.)*b2*b-9*5.*sqrt(5.)*a*b2+9*5.*sqrt(5.)*a2*b-3*5.*sqrt(5.)*a2*a)*t+5*b2*b2-20*a*b2*b+30*a2*b2-20*a2*a*b+5*a2*a2)*exp(sqrt(5.)*(a-b)/t)/(54*t2*t2);
          p4 = -((t*(t*(9*t*(7*t-5.*sqrt(5.)*(b+a-2))+10*b*(5*b+17*a-27)+10*(5*a2-27*a+27))-8*5.*sqrt(5.)*(a-1)*(b-1)*(b+a-2))+50*(a-1)*(a-1)*(b-2)*b+50*(a-1)*(a-1))*exp(2*sqrt(5.)*b/t))*exp(-sqrt(5.)*(b-a+2)/t)/(36*sqrt(5.)*t*t2);
          
          Wijs(j,i) = Wijs(i,j) *= p1 + p3 + p4;
        }
        
      }
    }
  }
  
  return(Wijs);   
}


double c1i_mat52(double a, double b, double t){
  double p1,p3,p4,t2,a2,b2,dw;
  
  bool boo = true;
  
  if(b < a){
    double c = b;
    b = a;
    a = c;
    boo = false;
  }
  
  
  t2 = t*t;
  a2 = a*a;
  b2 = b*b;
  
  p1 = (2*t2*(63*t2 + 9*5.*sqrt(5.)*b*t-9*5.*sqrt(5.)*a*t+50*b2-100*a*b+50*a2)*exp(2*sqrt(5.)*a/t)-63*t2*t2-9*5.*sqrt(5.)*(b+a)*t*t2-10*(5*b2+17*a*b+5*a2)*t2-8*5.*sqrt(5.)*a*b*(b+a)*t-50*a2*b2)*exp(-sqrt(5.)*(b+a)/t)/(36*sqrt(5.)*t*t2);
  p3 = (b-a)*(54*t2*t2+(54*sqrt(5.)*b-54*sqrt(5.)*a)*t*t2+(105*b2-210*a*b+105*a2)*t2+(3*5.*sqrt(5.)*b2*b-9*5.*sqrt(5.)*a*b2+9*5.*sqrt(5.)*a2*b-3*5.*sqrt(5.)*a2*a)*t+5*b2*b2-20*a*b2*b+30*a2*b2-20*a2*a*b+5*a2*a2)*exp(sqrt(5.)*(a-b)/t)/(54*t2*t2);
  p4 = -((t*(t*(9*t*(7*t-5.*sqrt(5.)*(b+a-2))+10*b*(5*b+17*a-27)+10*(5*a2-27*a+27))-8*5.*sqrt(5.)*(a-1)*(b-1)*(b+a-2))+50*(a-1)*(a-1)*(b-2)*b+50*(a-1)*(a-1))*exp(2*sqrt(5.)*b/t))*exp(-sqrt(5.)*(b-a+2)/t)/(36*sqrt(5.)*t*t2);
  
  if((p1 + p3 + p4) == 0.) return(0.); 
  
  if(!boo){
    dw = -((2*5*sqrt(5.)*exp(2*sqrt(5.)/t)*a2*a2*a+(-100*t-2*5*5*sqrt(5.)*b)*exp(2*sqrt(5.)/t)*a2*a2+(18*5*sqrt(5.)*t2+400*b*t+4*5*5*sqrt(5.)*b2)*exp(2*sqrt(5.)/t)*a2*a+((150*t2*t+(24*5*sqrt(5.)-24*5*sqrt(5.)*b)*t2+(150*b2-300*b+150)*t)*exp(2*sqrt(5.)*b/t)+
      (-210*t2*t-54*5*sqrt(5.)*b*t2-600*b2*t-4*5*5*sqrt(5.)*b2*b)*exp(2*sqrt(5.)/t))*a2+((-3*5*5*sqrt(5.)*t2*t2+(270*b-570)*t2*t+(-12*5*sqrt(5.)*b2+72*5*sqrt(5.)*b-12*5*5*sqrt(5.))*t2+(-300*b2+600*b-300)*t)*exp(2*sqrt(5.)*b/t)+(42*sqrt(5.)*t2*t2+420*b*t2*t+
      54*5*sqrt(5.)*b2*t2+400*b2*b*t+2*5*5*sqrt(5.)*b2*b2)*exp(2*sqrt(5.)/t))*a+(54*t2*t2*t+(108*sqrt(5.)-33*sqrt(5.)*b)*t2*t2+(30*b2-330*b+450)*t2*t+(12*5*sqrt(5.)*b2-48*5*sqrt(5.)*b+36*5*sqrt(5.))*t2+(150*b2-300*b+150)*t)*exp(2*sqrt(5.)*b/t)+(-42*sqrt(5.)*b*t2*t2-
      210*b2*t2*t-18*5*sqrt(5.)*b2*b*t2-100*b2*b2*t-2*5*sqrt(5.)*b2*b2*b)*exp(2*sqrt(5.)/t))*exp(2*sqrt(5.)*a/t)+(-150*t2*t-24*5*sqrt(5.)*b*t2-150*b2*t)*exp(2*sqrt(5.)/t)*a2+(-3*5*5*sqrt(5.)*t2*t2-270*b*t2*t-12*5*sqrt(5.)*b2*t2)*exp(2*sqrt(5.)/t)*a+(-54*t2*t2*t-33*sqrt(5.)*b*t2*t2-
      30*b2*t2*t)*exp(2*sqrt(5.)/t))*exp(-sqrt(5.)*a/t-sqrt(5.)*b/t-2*sqrt(5.)/t)/(108*t2*t2*t);
  }else{
    dw = -(((150*t2*t+(24*5*sqrt(5.)-24*5*sqrt(5.)*a)*t2+(150*a2-300*a+150)*t)*exp(2*sqrt(5.)*a/t)*b2+(-3*5*5*sqrt(5.)*t2*t2+(270*a-570)*t2*t+(-12*5*sqrt(5.)*a2+72*5*sqrt(5.)*a-12*5*5*sqrt(5.))*t2+(-300*a2+600*a-300)*t)*exp(2*sqrt(5.)*a/t)*b+
      (54*t2*t2*t+(108*sqrt(5.)-33*sqrt(5.)*a)*t2*t2+(30*a2-330*a+450)*t2*t+(12*5*sqrt(5.)*a2-48*5*sqrt(5.)*a+36*5*sqrt(5.))*t2+(150*a2-300*a+150)*t)*exp(2*sqrt(5.)*a/t))*exp(2*sqrt(5.)*b/t)+2*5*sqrt(5.)*exp(2*sqrt(5.)*a/t+2*sqrt(5.)/t)*b2*b2*b+
      (100*t-2*5*5*sqrt(5.)*a)*exp(2*sqrt(5.)*a/t+2*sqrt(5.)/t)*b2*b2+(18*5*sqrt(5.)*t2-400*a*t+4*5*5*sqrt(5.)*a2)*exp(2*sqrt(5.)*a/t+2*sqrt(5.)/t)*b2*b+((210*t2*t-54*5*sqrt(5.)*a*t2+600*a2*t-4*5*5*sqrt(5.)*a2*a)*exp(2*sqrt(5.)*a/t+
      2*sqrt(5.)/t)+(-150*t2*t-24*5*sqrt(5.)*a*t2-150*a2*t)*exp(2*sqrt(5.)/t))*b2+((42*sqrt(5.)*t2*t2-420*a*t2*t+54*5*sqrt(5.)*a2*t2-400*a2*a*t+2*5*5*sqrt(5.)*a2*a2)*exp(2*sqrt(5.)*a/t+2*sqrt(5.)/t)+
      (-3*5*5*sqrt(5.)*t2*t2-270*a*t2*t-12*5*sqrt(5.)*a2*t2)*exp(2*sqrt(5.)/t))*b+(-42*sqrt(5.)*a*t2*t2+210*a2*t2*t-18*5*sqrt(5.)*a2*a*t2+100*a2*a2*t-2*5*sqrt(5.)*a2*a2*a)*exp(2*sqrt(5.)*a/t+2*sqrt(5.)/t)+
      (-54*t2*t2*t-33*sqrt(5.)*a*t2*t2-30*a2*t2*t)*exp(2*sqrt(5.)/t))*exp(-sqrt(5.)*b/t-sqrt(5.)*a/t-2*sqrt(5.)/t)/(108*t2*t2*t);
  }
  
  return(dw/(p1+p3+p4));
}

// [[Rcpp::export]]
double c2_mat52_cpp(double x, double t, double w){
  double x2 = x*x;
  double t2 = t*t;
  if(w == 0.) return(0.);
  
  double tmp = (exp(-2*sqrt(5.)*x/t)*(63*t2*t2*exp(2*sqrt(5.)*x/t)-50*x2*x2-16*5*sqrt(5.)*t*x2*x-270*t2*x2-18*5*sqrt(5.)*t2*t*x-63*t2*t2)-exp(-2*sqrt(5.)/t)*((t*(t*(10*(5*x2-27*x+27)+9*t*(7*t-5*sqrt(5.)*(2*x-2))+10*x*(22*x-27))-8*5*sqrt(5.)*(x-1)*(x-1)*(2*x-2))+50*(x-2)*(x-1)*(x-1)*x+50*(x-1)*(x-1))*exp(2*sqrt(5.)*x/t)-63*t2*t2*exp(2*sqrt(5.)/t)))/(36*sqrt(5.)*t2*t);
  if(tmp == 0.) return(0.);
  double dw = -((25*x2*x2-2*(3*5*sqrt(5.)*t+50)*x2*x+3*(t*(25*t+6*5*sqrt(5.))+50)*x2-2*(3*t*(t*(3*sqrt(5.)*t+25)+3*5*sqrt(5.))+50)*x+9*t2*t2+18*sqrt(5.)*t2*t+75*t2+6*5*sqrt(5.)*t+25)*exp(4*sqrt(5.)*x/t)-25*exp(2*sqrt(5.)/t)*x2*x2-6*5*sqrt(5.)*t*exp(2*sqrt(5.)/t)*x2*x-75*t2*exp(2*sqrt(5.)/t)*x2-18*sqrt(5.)*t2*t*exp(2*sqrt(5.)/t)*x-9*t2*t2*exp(2*sqrt(5.)/t))*exp(-2*sqrt(5.)*(x+1)/t)/(9*t2*t2);
  return(dw*w/tmp);
}


// [[Rcpp::export]]
NumericVector c1_mat52_cpp(NumericVector X, double x, double sigma, NumericVector W){
  NumericVector cis(X.length());  
  
  for(int i = 0; i < X.length(); i++){
    cis(i) += c1i_mat52(X(i), x, sigma) * W(i);
  }
  
  return(cis);
}

// [[Rcpp::export]]
NumericVector d_mat52_cpp(NumericVector X, double x, double sigma){
  NumericVector s(X.length());
  double tmp;
  
  for(int i = 0; i < X.length(); i++){
    tmp = (x - X(i))/sigma;
    if(tmp > 0){
      s(i) = ((10./3. - 5.) * tmp - 5 * sqrt(5.)/3. * tmp * tmp) / (1. + sqrt(5.) * tmp + 5./3. * tmp * tmp);
    }else{
      if(tmp == 0){
        s(i) = 0;
      }else{
        tmp = std::abs(tmp);
        s(i) = -((10./3. - 5.) * tmp - 5. * sqrt(5.)/3. * tmp * tmp) / ((1. + sqrt(5.) * tmp + 5./3. * tmp * tmp));
      }
    }
  }
  return(s/sigma);
}

//////// Matern 3/2 kernel

double A_1_cpp(double x){
  return((2. + x) * exp(-x));
}

// [[Rcpp::export]]
NumericVector mi_mat32_cpp(NumericMatrix Mu, NumericVector sigma){
  NumericVector mis(Mu.nrow(), 1.);
  
  for(int i = 0; i < Mu.nrow(); i++){
    for(int j = 0; j < Mu.ncol(); j++){
      mis(i) *= sigma(j)/(sqrt(3.)) * (4. - A_1_cpp(sqrt(3.) * Mu(i,j)/sigma(j)) - A_1_cpp(sqrt(3.) * (1. - Mu(i,j))/sigma(j)));
    }
  }
  return(mis);
}

// [[Rcpp::export]]
NumericMatrix Wijs_mat32_cpp(NumericMatrix Mu1, NumericMatrix Mu2, NumericVector sigma){
  
  int m1c = Mu1.ncol();
  int m2r = Mu2.nrow();
  double tmp;
  NumericMatrix Wijs(Mu1.nrow(), m2r);
  Wijs.fill(1.);
  
  double a,b,t,t2;
  double p;
  
  for(int i = 0; i < Mu1.nrow(); i++){
    for(int j = 0; j < m2r; j++){
      const double* ptr_s = (const double*) &sigma(0);
      for(int k = 0; k < m1c; k++, ptr_s++){
        a = Mu1(i,k);
        b = Mu2(j,k);
        
        if(b < a){
          tmp = b; 
          b = a;
          a = tmp;
        }
        
        t = *ptr_s;
        t2 = t*t;
        p = ((t*(5*sqrt(3.)*t+9*b-9*a)*exp((2*sqrt(3.)*a)/t)-5*sqrt(3.)*t2-9*(b+a)*t-2*3*sqrt(3.)*a*b)*exp(-(sqrt(3.)*(b+a))/t))/(12*t)+((b-a)*(2*t2+2*sqrt(3.)*(b-a)*t+b*b-2*a*b+a*a)*exp(-(sqrt(3.)*(b-a))/t))/(2*t2)-(((t*(5*t-3*sqrt(3.)*(b+a-2))+6*(a-1)*b-6*a+6)*exp((2*sqrt(3.)*b)/t)-t*(5*t+3*sqrt(3.)*(b-a))*exp((2*sqrt(3.))/t))*exp(-(sqrt(3.)*(b-a+2))/t))/(4*sqrt(3.)*t);
        
        Wijs(i,j) *= p;
        
      }
    }
  } 
  
  return(Wijs);   
}

// [[Rcpp::export]]
NumericMatrix Wijs_mat32_sym_cpp(NumericMatrix Mu, NumericVector sigma){
  int m1c = Mu.ncol();
  NumericMatrix Wijs(Mu.nrow(), Mu.nrow());
  Wijs.fill(1.);
  
  double a,b,t,t2;
  double p;
  double tmp;
  
  for(int i = 0; i < Mu.nrow(); i++){
    for(int j = 0; j <= i; j++){
      const double* ptr_s = (const double*) &sigma(0);
      for(int k = 0; k < m1c; k++, ptr_s++){
        a = Mu(j,k);
        b = Mu(i,k);
        if(b < a){
          tmp = b;
          b = a;
          a = tmp;
        }
        
        t = *ptr_s;
        t2 = t*t;
        
        if(i == j){
          Wijs(i,j) *= (15*t2-(t*(15*t-2*9*sqrt(3.)*(a-1))+18*(a-1)*(a-1))*exp((2*sqrt(3.)*a)/t-(2*sqrt(3.))/t))/(4*3*sqrt(3.)*t) -((5*t2+2*3*sqrt(3.)*a*t+6*a*a)*exp(-(2*sqrt(3.)*a)/t)-5*t2)/(4*sqrt(3.)*t);
        }else{
          p = ((t*(5*sqrt(3.)*t+9*b-9*a)*exp((2*sqrt(3.)*a)/t)-5*sqrt(3.)*t2-9*(b+a)*t-2*3*sqrt(3.)*a*b)*exp(-(sqrt(3.)*(b+a))/t))/(12*t)+((b-a)*(2*t2+2*sqrt(3.)*(b-a)*t+b*b-2*a*b+a*a)*exp(-(sqrt(3.)*(b-a))/t))/(2*t2)-(((t*(5*t-3*sqrt(3.)*(b+a-2))+6*(a-1)*b-6*a+6)*exp((2*sqrt(3.)*b)/t)-t*(5*t+3*sqrt(3.)*(b-a))*exp((2*sqrt(3.))/t))*exp(-(sqrt(3.)*(b-a+2))/t))/(4*sqrt(3.)*t);
          
          Wijs(j,i) = Wijs(i,j) *= p;
        }
        
      }
    }
  }
  
  return(Wijs);   
}

// [[Rcpp::export]]
NumericVector d_mat32_cpp(NumericVector X, double x, double sigma){
  NumericVector s(X.length());
  double tmp;
  
  for(int i = 0; i < X.length(); i++){
    tmp = (x - X(i))/sigma;
    if(tmp > 0){
      s(i) = -3*tmp / (1 + sqrt(3.) * tmp);
    }else{
      if(tmp == 0){
        s(i) = 0;
      }else{
        tmp = -tmp;
        s(i) = 3*tmp / (1 + sqrt(3.) * tmp);
      }
    }
  }
  return(s/sigma);
}

double c1i_mat32(double a, double b, double t){
  double a2, b2, p,t2,dw;
  
  bool boo = true;
  
  if(b < a){
    double c = b;
    b = a;
    a = c;
    boo = false;
  }
  
  a2 = a*a;
  b2 = b*b;
  
  t2 = t*t;
  
  p = ((t*(5*sqrt(3.)*t+9*b-9*a)*exp((2*sqrt(3.)*a)/t)-5*sqrt(3.)*t2-9*(b+a)*t-2*3*sqrt(3.)*a*b)*exp(-(sqrt(3.)*(b+a))/t))/(12*t)+((b-a)*(2*t2+2*sqrt(3.)*(b-a)*t+b*b-2*a*b+a*a)*exp(-(sqrt(3.)*(b-a))/t))/(2*t2)-(((t*(5*t-3*sqrt(3.)*(b+a-2))+6*(a-1)*b-6*a+6)*exp((2*sqrt(3.)*b)/t)-t*(5*t+3*sqrt(3.)*(b-a))*exp((2*sqrt(3.))/t))*exp(-(sqrt(3.)*(b-a+2))/t))/(4*sqrt(3.)*t);

  if(p == 0.) return(0.);
  
  if(!boo){
    dw =-(((2*sqrt(3.)*exp((2*sqrt(3.))/t)*a2*a+(-6*t-2*3*sqrt(3.)*b)*exp((2*sqrt(3.))/t)*a2+(((6*b-6)*t-3*sqrt(3.)*t2)*exp((2*sqrt(3.)*b)/t)+(2*sqrt(3.)*t2+12*b*t+2*3*sqrt(3.)*b2)*exp((2*sqrt(3.))/t))*a+(2*t2*t+(4*sqrt(3.)-sqrt(3.)*b)*t2+(6-6*b)*t)*exp((2*sqrt(3.)*b)/t)+(-2*sqrt(3.)*b*t2-6*b2*t-2*sqrt(3.)*b2*b)*exp((2*sqrt(3.))/t))*exp((2*sqrt(3.)*a)/t)+(-3*sqrt(3.)*t2-6*b*t)*exp((2*sqrt(3.))/t)*a+(-2*t2*t-sqrt(3.)*b*t2)*exp((2*sqrt(3.))/t))*exp((-sqrt(3.)*a-sqrt(3.)*b-2*sqrt(3.))/t))/(4*t2*t);
  }else{
    dw = ((t*((3*sqrt(3.)*t-6*a+6)*b-t*(2*t-sqrt(3.)*(a-4))+6*a-6)*exp((2*sqrt(3.)*(b+a))/t)-2*sqrt(3.)*exp((2*sqrt(3.)*(a+1))/t)*b2*b-2*(3*t-3*sqrt(3.)*a)*exp((2*sqrt(3.)*(a+1))/t)*b2-exp((2*sqrt(3.))/t)*(2*sqrt(3.)*t2*exp((2*sqrt(3.)*a)/t)-12*a*t*exp((2*sqrt(3.)*a)/t)+2*3*sqrt(3.)*a2*exp((2*sqrt(3.)*a)/t)-3*sqrt(3.)*t2-6*a*t)*b+2*a*(sqrt(3.)*t2-3*a*t+sqrt(3.)*a2)*exp((2*sqrt(3.)*(a+1))/t)+t2*(2*t+sqrt(3.)*a)*exp((2*sqrt(3.))/t))*exp(-(sqrt(3.)*(b+a+2))/t))/(4*t2*t);
  }
  
  return(dw/p);
}

// [[Rcpp::export]]
double c2_mat32_cpp(double x, double t, double w){
  double x2 = x*x;
  double t2 = t*t;
  
  if(w == 0.) return(0.);
  
  double tmp = (15*t2-(t*(15*t-2*9*sqrt(3.)*(x-1))+18*(x-1)*(x-1))*exp((2*sqrt(3.)*x)/t-(2*sqrt(3.))/t))/(4*3*sqrt(3.)*t) -((5*t2+2*3*sqrt(3.)*x*t+6*x*x)*exp(-(2*sqrt(3.)*x)/t)-5*t2)/(4*sqrt(3.)*t);
  
  if(tmp == 0.) return(0.);
  
  double dw = -(((3*x*x-2*(sqrt(3.)*t+3)*x+t2+2*sqrt(3.)*t+3)*exp((4*sqrt(3.)*x)/t)-3*exp((2*sqrt(3.))/t)*x2-2*sqrt(3.)*t*exp((2*sqrt(3.))/t)*x-t2*exp((2*sqrt(3.))/t))*exp(-(2*sqrt(3.)*(x+1))/t))/t2;
  return(dw*w/tmp);
}

// [[Rcpp::export]]
NumericVector c1_mat32_cpp(NumericVector X, double x, double sigma, NumericVector W){
  NumericVector cis(X.length());
  
  for(int i = 0; i < X.length(); i++){
    cis(i) += c1i_mat32(X(i), x, sigma) * W(i);
  }
  
  return(cis);
}
