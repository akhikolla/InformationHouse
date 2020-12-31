#include <Rcpp.h>
using namespace Rcpp;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List splitt2(NumericMatrix X, NumericMatrix Y, int m_feature, NumericVector Index, NumericMatrix Inv_Cov_Y, int Command, NumericVector ff){
  int Index_size = Index.size();
  NumericVector Index1(Index_size);
  for (int ii = 0; ii < Index_size; ii++) {
    Index1[ii] = Index[ii] - 1;
  }
  int nrowx = Index_size, ncolx = X.ncol();
  NumericMatrix x(nrowx, ncolx);
  for (int ii = 0; ii < nrowx; ii++) {
    for (int jj = 0; jj < ncolx; jj++) {
      x(ii, jj) = X(Index1[ii], jj);           // change possible
    }
  }
  int nrowy = Index_size, ncoly = Y.ncol();
  NumericMatrix y(nrowy, ncoly);
  for (int ii = 0; ii < nrowy; ii++) {
    for (int jj = 0; jj < ncoly; jj++) {
      y(ii, jj) = Y(Index1[ii], jj);
    }
  }
  NumericVector f(m_feature);
  for (int ii = 0; ii < m_feature; ii++) {
    f[ii] = ff[ii];
  }
  double min_score = 0;
  double D;
  NumericVector Threshold_value(2);
  int Feature_number;
  IntegerVector Idx_left, Idx_right;
  NumericVector xj(nrowx), xjj(nrowx);
  NumericMatrix Result(m_feature, (nrowx - 1)), DL((nrowx - 1), m_feature), DR((nrowx - 1), m_feature);
  //NumericVector Idx(nrowx);
  for (int mm = 0; mm<m_feature; mm++) {
    NumericVector xj(nrowx), xjj(nrowx);
    for (int ll = 0; ll<nrowx; ll++) {
      xj(ll) = x(ll, (f[mm] - 1));
      xjj(ll) = xj(ll);
    }
    std::sort(xj.begin(), xj.end());
    NumericVector Idx(nrowx);
    for (int ii = 0; ii<nrowx; ii++) {
      for (int jj = 0; jj<nrowx; jj++) {
        if (xj[ii] == xjj[jj]) {
          Idx[ii] = jj;
          xjj[jj] = 0.00000001;
          break;
        }
      }
    }

    for (int nn = 0; nn<(nrowx - 1); nn++) {
      NumericVector IIdx(nn + 1), IIdx2(nrowx - (nn + 1));
      //NumericVector Idx_left_temp(IIdx.size()), Idx_right_temp(IIdx2.size());
      NumericMatrix yk_left(IIdx.size(), ncoly), yk_right(IIdx2.size(), ncoly);
      for (int uu = 0; uu<ncoly; uu++) {
        for (int pp = 0; pp<IIdx.size(); pp++) {
          yk_left(pp, uu) = y(Idx(pp), uu);
        }
        for (int qq = 0; qq<IIdx2.size(); qq++) {
          int Ikdx = qq + IIdx.size();
          yk_right(qq, uu) = y(Idx(Ikdx), uu);
        }
      }
      /* Calculating node cost for left node */
      double D_left = 0;
      if (Command == 1) {
        int N = yk_left.nrow();
        double sum_of_elems = 0;
        for (int i = 0; i<N; ++i) {
          sum_of_elems += yk_left[i];
        }
        double ybar = sum_of_elems / N;
        for (int i = 0; i<N; i++) {
          D_left += pow((yk_left[i] - ybar), 2);
        }
      }
      else if (Command == 2) {
        int nrow = yk_left.nrow(), ncol = yk_left.ncol();
        NumericVector ybar(ncol);
        NumericMatrix yhat(nrow,ncol);
        for (int i = 0; i<ncol; i++) {
          double sum_of_elements = 0;
          for (int j = 0; j<nrow; j++) {
            sum_of_elements += yk_left(j, i);
          }
          ybar[i] = sum_of_elements / nrow;
          for (int j = 0; j<nrow; j++) {
            yhat(j,i) = yk_left(j, i) - ybar[i];
          }
        }

        NumericMatrix mult(nrow,ncol);
        /* Initializing elements of matrix mult to 0.*/
        int yhat_nrow = yhat.nrow();//sizeof(yhat) / sizeof(yhat[0]);
        for (int i = 0; i<yhat_nrow; ++i)
          for (int j = 0; j<Inv_Cov_Y.ncol(); ++j)
          {
            mult(i,j) = 0;
          }
          /* Multiplying matrix a and b and storing in array mult. */
          for (int i = 0; i<yhat_nrow; ++i)
            for (int j = 0; j<Inv_Cov_Y.ncol(); ++j)
              for (int k = 0; k<Inv_Cov_Y.nrow(); ++k)
              {
                mult(i,j) += yhat(i,k) * Inv_Cov_Y(k, j);
              }

              NumericMatrix mult2(nrow,nrow);
        int r1 = mult.nrow();//sizeof(mult) / sizeof(mult[0]);
        int c1 = mult.ncol();//sizeof(mult[0]) / sizeof(mult[0][0]);
        int c2 = yhat.nrow();//sizeof(yhat) / sizeof(yhat[0]);
        //int r2=sizeof(yhat[0]) / sizeof(yhat[0][0]);
        for (int i = 0; i<r1; ++i)
          for (int j = 0; j<c2; ++j)
          {
            mult2(i,j) = 0;
          }

          for (int i = 0; i<r1; ++i)
            for (int j = 0; j<c2; ++j)
              for (int k = 0; k<c1; ++k)
              {
                mult2(i,j) += mult(i,k) * yhat(j,k);
              }

              for (int kk = 0; kk<nrow; kk++) {
                D_left += mult2(kk,kk);
              }
      }
      DL(nn, mm) = D_left;
      /* Calculating node cost for right node */
      double D_right = 0;
      if (Command == 1) {
        int N = yk_right.nrow();
        double sum_of_elems = 0;
        for (int i = 0; i<N; ++i) {
          sum_of_elems += yk_right[i];
        }
        double ybar = sum_of_elems / N;
        for (int i = 0; i<N; i++) {
          D_right += pow((yk_right[i] - ybar), 2);
        }
      }
      else if (Command == 2) {
        int nrow = yk_right.nrow(), ncol = yk_right.ncol();
        NumericVector ybar(ncol);
        NumericMatrix yhat(nrow,ncol);
        for (int i = 0; i<ncol; i++) {
          double sum_of_elements = 0;
          for (int j = 0; j<nrow; j++) {
            sum_of_elements += yk_right(j, i);
          }
          ybar[i] = sum_of_elements / nrow;
          for (int j = 0; j<nrow; j++) {
            yhat(j,i) = yk_right(j, i) - ybar[i];
          }
        }

        NumericMatrix mult(nrow,ncol);
        /* Initializing elements of matrix mult to 0.*/
        int yhat_nrow = yhat.nrow();//sizeof(yhat) / sizeof(yhat[0]);
        for (int i = 0; i<yhat_nrow; ++i)
          for (int j = 0; j<Inv_Cov_Y.ncol(); ++j)
          {
            mult(i,j) = 0;
          }
          /* Multiplying matrix a and b and storing in array mult. */
          for (int i = 0; i<yhat_nrow; ++i)
            for (int j = 0; j<Inv_Cov_Y.ncol(); ++j)
              for (int k = 0; k<Inv_Cov_Y.nrow(); ++k)
              {
                mult(i,j) += yhat(i,k) * Inv_Cov_Y(k, j);
              }

              NumericMatrix mult2(nrow,nrow);
        int r1 = mult.nrow();//sizeof(mult) / sizeof(mult[0]);
        int c1 = mult.ncol();//sizeof(mult[0]) / sizeof(mult[0][0]);
        int c2 = yhat.nrow();//sizeof(yhat) / sizeof(yhat[0]);
        //int r2=sizeof(yhat[0]) / sizeof(yhat[0][0]);
        for (int i = 0; i<r1; ++i)
          for (int j = 0; j<c2; ++j)
          {
            mult2(i,j) = 0;
          }

          for (int i = 0; i<r1; ++i)
            for (int j = 0; j<c2; ++j)
              for (int k = 0; k<c1; ++k)
              {
                mult2(i,j) += mult(i,k) * yhat(j,k);
              }

              for (int kk = 0; kk<nrow; kk++) {
                D_right += mult2(kk,kk);
              }
      }
      DR(nn, mm) = D_right;
      D = D_left + D_right;
      //Result(mm,nn)=D;
      int temp1, temp2, temp3, temp4;
      if (mm == 0 && nn == 0) {
        min_score = D;
        int Feature_number_temp;
        Feature_number_temp = f[mm];
        Threshold_value = (xj[nn]+xj[nn + 1])/2 ;
        //Threshold_value[1] = xj[nn + 1];
        NumericVector Idx_left_temp(IIdx.size()), Idx_right_temp(IIdx2.size());
        for (int tt = 0; tt<IIdx.size(); tt++) {
          //Idx_left_temp[tt] = Idx[tt];
          temp1=Idx[tt];
          Idx_left_temp[tt]=Index1[temp1];
        }
        for (int rr = 0; rr<IIdx2.size(); rr++) {
          //Idx_right_temp[rr] = Idx[rr + IIdx.size()];
          temp2=Idx[rr + IIdx.size()];
          Idx_right_temp[rr]=Index1[temp2];
        }
        Feature_number = Feature_number_temp;
        Idx_left = Idx_left_temp;
        Idx_right = Idx_right_temp;
      }
      if (D< min_score) {
        min_score = D;
        int Feature_number_temp;
        Feature_number_temp = f[mm];
        Threshold_value = (xj[nn]+xj[nn + 1])/2 ;
        //Threshold_value[1] = xj[nn + 1];
        NumericVector Idx_left_temp(IIdx.size()), Idx_right_temp(IIdx2.size());
        for (int tt = 0; tt<IIdx.size(); tt++) {
          //Idx_left_temp[tt] = Idx[tt];
          temp3=Idx[tt];
          Idx_left_temp[tt]=Index1[temp3];
        }
        for (int rr = 0; rr<IIdx2.size(); rr++) {
          //Idx_right_temp[rr] = Idx[rr + IIdx.size()];
          temp4=Idx[rr + IIdx.size()];
          Idx_right_temp[rr]=Index1[temp4];
        }
        Feature_number = Feature_number_temp;
        Idx_left = Idx_left_temp;
        Idx_right = Idx_right_temp;
      }
    }
  }
  for (int ii = 0; ii < Idx_left.size(); ii++) {
    Idx_left[ii] = Idx_left[ii] + 1;
  }
  for (int ii = 0; ii < Idx_right.size(); ii++) {
    Idx_right[ii] = Idx_right[ii] + 1;
  }
  return Rcpp::List::create(Rcpp::Named("Idx_left") = Idx_left,
                            Rcpp::Named("Idx_right") = Idx_right,
                            Rcpp::Named("Feature_number") = Feature_number,
                            Rcpp::Named("Threshold_value") = Threshold_value);
}

// [[Rcpp::export]]
double Node_cost(NumericMatrix y, NumericMatrix Inv_Cov_Y, int Command){
  double D=0;
  if (Command==1){
    int N=y.nrow();
    double sum_of_elems=0;
    for (int i=0;i<N; ++i){
      sum_of_elems +=y[i];
    }
    double ybar=sum_of_elems/N ;
    for (int i=0;i<N;i++){
      D +=pow((y[i]-ybar),2);
    }
  } else if (Command==2){
    int nrow=y.nrow(), ncol=y.ncol();
    NumericVector ybar(ncol);
    NumericMatrix yhat(nrow,ncol);
    for (int i=0;i<ncol;i++){
      double sum_of_elements = 0;
      for (int j=0;j<nrow;j++){
        sum_of_elements +=y(j,i);
      }
      ybar[i]=sum_of_elements/nrow ;
      for (int j=0; j<nrow;j++){
        yhat(j,i) =y(j,i)-ybar[i];
      }
    }

    NumericMatrix mult(nrow,ncol);
    /* Initializing elements of matrix mult to 0.*/
    int yhat_nrow=yhat.nrow();//sizeof(yhat) / sizeof(yhat[0]);
    for(int i=0; i<yhat_nrow; ++i)
      for(int j=0; j<Inv_Cov_Y.ncol(); ++j)
      {
        mult(i,j)=0;
      }
      /* Multiplying matrix a and b and storing in array mult. */
      for(int i=0; i<yhat_nrow; ++i)
        for(int j=0; j<Inv_Cov_Y.ncol(); ++j)
          for(int k=0; k<Inv_Cov_Y.nrow(); ++k)
          {
            mult(i,j)+=yhat(i,k)*Inv_Cov_Y(k,j);
          }

          NumericMatrix mult2(nrow,nrow);
    int r1=mult.nrow();//sizeof(mult) / sizeof(mult[0]);
    int c1=mult.ncol();//sizeof(mult[0]) / sizeof(mult[0][0]);
    int c2=yhat.nrow();//sizeof(yhat) / sizeof(yhat[0]);
    //int r2=sizeof(yhat[0]) / sizeof(yhat[0][0]);
    for(int i=0; i<r1; ++i)
      for(int j=0; j<c2; ++j)
      {
        mult2(i,j)=0;
      }

      for(int i=0; i<r1; ++i)
        for(int j=0; j<c2; ++j)
          for(int k=0; k<c1; ++k)
          {
            mult2(i,j)+=mult(i,k)*yhat(j,k);
          }

          for(int kk=0; kk<nrow; kk++){
            D +=mult2(kk,kk);
          }
  }

  return D;
}
