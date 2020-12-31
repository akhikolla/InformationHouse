/* Routines for GeDS R-package
 *
 * Routines are based on the algorithm in Kaishev et alt. (2016).
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 *
 * These functions are distributed WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.	 See the GNU General Public License for more details.
 */


#include <Rcpp.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;


// [[Rcpp::export]]
int whmx(NumericVector vettore) {
  int i, k=0, lun = vettore.size();
  double massimo = vettore(0);
  for(i=0; i<lun; i++){
    if(vettore(i) > massimo) {
      massimo = vettore(i);
      k = i;
    }
  }
  return k;
}

NumericVector Knotnewbis(NumericVector wht, NumericVector restmp, NumericVector x,
                  NumericVector dcm, NumericVector oldknots, NumericVector oldintknots){
  int indice, dcl, j, kk, u=dcm.size();
  int dcu, nold = oldknots.size(), noldint = oldintknots.size(), i, xdim = x.size();
  double sup, inf, newknot;
  bool temp, temp2, temp3;
  for (kk=0; kk< u; kk++) {
    indice = whmx(wht);
    if (indice==0) {
      dcl = 0;
      } else {
        dcl = dcm[indice-1];
        }

    dcu = dcm[indice]-1;
    sup = x(dcu);
    inf = x(dcl);
      newknot = std::inner_product(restmp.begin() + dcl,restmp.begin() + dcu + 1, x.begin() + dcl, 0.0);
      newknot = newknot /  std::accumulate(restmp.begin() + dcl, restmp.begin() + dcu + 1, 0.0);  //  / std::accumulate()



    NumericVector sortedknots(nold+1);
    sortedknots(nold) = newknot;
    for (i=0; i< nold; i++){sortedknots(i) = oldknots(i);
    }
    std::sort(sortedknots.begin(),sortedknots.end());

    temp2 = true;
    for(i=0;i<nold-3;i++){ // dipende dal nold
      temp = false;
      for(j=0; j < xdim; j++){
        temp = temp || ((sortedknots(i) < x(j)) && (sortedknots(i+3) > x(j)));
      }
      temp2 = temp2 && temp;
    }

    temp3 = false;
    if (inf != sup){
      for(j=0; j<noldint; j++){
        temp3 = temp3 || ((inf<=oldintknots(j))  && (sup>=oldintknots(j))) ; //se i nodi vecchi stanno tra le x
      }
    } else {
      for(j=0; j<noldint; j++){
        temp3 = temp3 || inf == oldintknots(j) ; //se i nodi vecchi stanno tra le x
      }
    }

//    break;
    if ((!temp2) || temp3) {
      wht(indice) = 0;
    } else {
      break;
    }
  }
  NumericVector ris(2);
  ris(0) = newknot;
  ris(1) = indice+1;
  return ris;
}


// [[Rcpp::export]]
NumericVector Knotnewtest(NumericVector wht, NumericVector restmp, NumericVector x,
                      NumericVector dcm, NumericVector oldknots, double tol){
  int indice, dcl, jj, kk, u=dcm.size();
  int dcu, nold = oldknots.size();
  int noldint = nold-6, i, xdim = x.size();
  long double sup, inf, newknot;
  bool temp, temp2, temp3;
  for (kk=0; kk< u; kk++) {
    indice = whmx(wht);
    if (indice==0) {
      dcl = 0;
    } else {
      dcl = dcm[indice-1];
    }

    dcu = dcm[indice]-1;
    sup = x(dcu);
    inf = x(dcl);



    newknot = std::inner_product(restmp.begin() + dcl,restmp.begin() + dcu + 1, x.begin() + dcl, 0.0);
    newknot = newknot /  std::accumulate(restmp.begin() + dcl, restmp.begin() + dcu + 1, 0.0);  //  / std::accumulate()



    NumericVector sortedknots(nold+1);
    sortedknots(nold) = newknot;
    for (i=0; i< nold; i++){sortedknots(i) = oldknots(i);
    }
    std::sort(sortedknots.begin(),sortedknots.end());

    temp2 = true;
    for(i=0;i<nold-2;i++){ // dipende dal nold
      temp = false;
      for(jj=0; jj < xdim; jj++){
        temp = temp || ((sortedknots(i) < x(jj)) && (sortedknots(i+3) > x(jj)));
        if(temp) {break;}
      }
      temp2 = temp2 && temp;
      if(!temp2){break;}
    }

    temp3 = false;
    if(noldint>0){
      NumericVector oldintknots(noldint);
      for(i=0; i<noldint; i++){
        oldintknots(i) = oldknots(6+i);
      }
      if (inf != sup){
        for(jj=0; jj<noldint; jj++){
          temp3 = temp3 || ((inf<=oldintknots(jj))  && (sup>=oldintknots(jj))) ; //se i nodi vecchi stanno tra le x
        }
      } else {
        for(jj=0; jj<noldint; jj++){
          temp3 = temp3 || std::abs(inf - oldintknots(jj)) < tol ; //se i nodi vecchi stanno tra le x
        }
      }
    }
    //    break;
    if ((!temp2) || temp3) {
      wht(indice) = 0;
    } else {
      break;
    }
  }
  NumericVector ris(2);
  ris(0) = newknot;
  ris(1) = indice+1;
  return ris;
}

// [[Rcpp::export]]
NumericVector Knotnew(NumericVector wht, NumericVector restmp, NumericVector x,
                      NumericVector dcm, NumericVector oldknots, double tol){
  int indice, dcl, jj, kk, u=dcm.size();
  int dcu, nold = oldknots.size();
  int noldint = nold-6, i, xdim = x.size();
  long double sup, inf, newknot;
  bool temp, temp2, temp3;
  for (kk=0; kk< u; kk++) {
    indice = whmx(wht);
    if (indice==0) {
      dcl = 0;
    } else {
      dcl = dcm[indice-1];
    }

    dcu = dcm[indice]-1;
    sup = x(dcu);
    inf = x(dcl);

    temp3 = false;
    if(noldint>0){
      NumericVector oldintknots(noldint);
      for(i=0; i<noldint; i++){
        oldintknots(i) = oldknots(6+i);
      }
      if (inf != sup){
        for(jj=0; jj<noldint; jj++){
          temp3 = temp3 || ((inf<=oldintknots(jj))  && (sup>=oldintknots(jj))) ; //se i nodi vecchi stanno tra le x
        }
      } else {
        for(jj=0; jj<noldint; jj++){
          temp3 = temp3 || std::abs(inf - oldintknots(jj)) < tol ; //se i nodi vecchi stanno tra le x
        }
      }
    }

    if (temp3) {
      wht(indice) = 0;
    } else {


    newknot = std::inner_product(restmp.begin() + dcl,restmp.begin() + dcu + 1, x.begin() + dcl, 0.0);
    newknot = newknot /  std::accumulate(restmp.begin() + dcl, restmp.begin() + dcu + 1, 0.0);  //  / std::accumulate()



    NumericVector sortedknots(nold+1);
    sortedknots(nold) = newknot;
    for (i=0; i< nold; i++){sortedknots(i) = oldknots(i);
    }
    std::sort(sortedknots.begin(),sortedknots.end());

    temp2 = true;
    for(i=0;i<nold-2;i++){ // dipende dal nold
      temp = false;
      for(jj=0; jj < xdim; jj++){
        temp = temp || ((sortedknots(i) < x(jj)) && (sortedknots(i+3) > x(jj)));
        if(temp) {break;}
      }
      temp2 = temp2 && temp;
      if(!temp2){break;}
    }

    //    break;
    if ((!temp2)) {
      wht(indice) = 0;
    } else {
      break;
    }
    }
    }
  NumericVector ris(2);
  ris(0) = newknot;
  ris(1) = indice+1;
  return ris;
}



/*

  tmp2 <- cbind(tmp[1:(j+3)],tmp[4:(j+6)])
  tmp3 <- logical(j+3)
  for(kk in 1:(j+3)){
    tmp3[kk]<-any(((tmp2[kk,1]<X) * (tmp2[kk,2]>X)))
  }
#check both conditions - no knots between Xs and Xs between knots
  if(all(tmp3) &&       (((d[indice]!=1) && (!any((nodi>=inf) * (nodi<=sup))))  ||
     ((d[indice]==1) && (!any(nodi==inf))) )
  ) {break } else {
    w[indice] <- NA
  }
*/

// [[Rcpp::export]]
NumericVector makenewknots(NumericVector knots, int degree) {
  /* 13Ys 13nodi 100data */
  // assumendo che siano tutti, 16
  int i, j, lun;  // n=512
  double cumtemp;
  lun = knots.size();
  NumericVector newknots(lun-(degree-2));// 16 - 3 = 13
  for(i=0; i< (lun-(degree-2)); ++i){
    for(j=0, cumtemp = 0.0; j<(degree-1); ++j) cumtemp += knots[i+j];
    newknots(i) = cumtemp / (degree - 1);
  }
  return newknots;
}

// [[Rcpp::export]]
NumericVector makeEpsilonsb(NumericVector data, NumericVector Xs, NumericVector Ys, int degree) {
  int lun = Xs.size();  /* 13Ys 13nodi 100data */
      // assumendo che siano tutti, 16
      int i, j, k, p, n = data.size();  // n=512
      double cumtemp;
      p = lun - degree;
      NumericVector epsilon(p);// 16 - 3 = 13
      for(i=0; i<p ; ++i){
        for(j=1, cumtemp = 0.0; j<degree; ++j) cumtemp += Xs[i+j];
        /*assumo che kk siano i nodi completi, non quelli interni */
        epsilon(i) = cumtemp / (degree - 1);
      }
      return epsilon;
}




NumericVector ctrlpolyfun(NumericVector data, NumericVector Xs, NumericVector Ys, int degree) {
  NumericVector epsilons;
  int lun = Xs.size();
  // assumendo che siano tutti, 16
  int i, j, k, p, n = data.size();  // n=512
  epsilons = makeEpsilonsb(data, Xs, Ys, degree);
  NumericVector vettore(n);
  p = lun - degree;
    for (i=0; i<n; i++){
      for (j = 1; j < p; j++) {
        if (epsilons[j] >= data[i]){k =  j; break;}
        }
      vettore[i] =Ys(k - 1) + (Ys(k) - Ys(k - 1)) * (data(i) - epsilons(k - 1)) / (epsilons(k) - epsilons(k - 1));
          }
    return vettore;
}

// vettore(i) = (Ys(k)-Ys(k-1))*(data(i)-epsilons(k-1))/(epsilons(k)-epsilons(k-1));


// [[Rcpp::export]]
NumericMatrix makeRatSplines(NumericMatrix matrice, NumericVector h) {
  int i, j, k, nr=matrice.nrow(), nc=matrice.ncol();
  NumericVector norm(nc);
  NumericMatrix matriceb(nr, nc);
  double s;
  for (j=k=0; j<nc; j++) {
    for (i=0, s=0.0; i<nr; i++,k++) {
      matriceb[k] = matrice[k]*h[i];
      s += matriceb[k];
      }
    norm[j] = s;
  }
  for (j=k=0; j<nc; j++) {
    for (i=0; i<nr; i++,k++) {
      matriceb[k] = matriceb[k]/norm[j];
    }
  }
  return matriceb;
}


// [[Rcpp::export]]
NumericVector makeWeights(NumericMatrix x){
  int i, nr = x.nrow();
  double r, s, t;
  NumericVector coseni(nr);
  for(i=0;i<nr-2;i++){
    r = (x(i,1)-x(i+1,1))*(x(i+2,1)-x(i+1,1))+(x(i,2)-x(i+1,2))*(x(i+2,2)-x(i+1,2));
    s = pow((x(i,1)-x(i+1,1)),2)+pow((x(i,2)-x(i+1,2)),2);
    t = pow((x(i+2,1)-x(i+1,1)),2)+pow((x(i+2,1)-x(i+1,2)),2);
    coseni(i) = r/(s*t);
    }
  return coseni;
}

// [[Rcpp::export]]
NumericMatrix tensorProd(NumericMatrix Xmat, NumericMatrix Ymat){
  int ii, jjx, jjy, kk, ib = Xmat.nrow(), jx = Xmat.ncol(), jy = Ymat.ncol();
  NumericMatrix ris(ib,jx*jy);
  for(jjx=0, kk=0; jjx<jx; jjx++){
    for(jjy=0;jjy<jy; jjy++){
      for(ii=0; ii<ib; ii++){
        ris(ii,kk) = Xmat(ii,jjx)*Ymat(ii,jjy);
      }
      kk++;
    }
  }
  return ris;
}



