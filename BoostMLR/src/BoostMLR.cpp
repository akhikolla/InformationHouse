#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double Sum_C(NumericVector x){
  int n = x.size();
  double Sum = 0;
  for(int i = 0;i < n;++i){
    Sum = Sum + x[i];
  }
  return Sum;
}

// [[Rcpp::export]]
double Sum_C_NA(NumericVector x){
  double Sum;
  if(all(is_na(x) )){
    Sum = NA_REAL;
  }else
  {
    int n = x.size();
    double x_temp;
    Sum = 0;
    for(int i = 0;i < n;++i){
      x_temp = x[i];
      if( (NumericVector::is_na( x_temp ) )) {
        x_temp = 0;
      }
      Sum = Sum + x_temp;
    }
  }
  return Sum;
}

// [[Rcpp::export]]
int length_C_NA(NumericVector x){
  int n = x.size();
  if(any(is_na(x))){
    n = n - sum(is_na(x));  
  }
  return n;
}

// [[Rcpp::export]]
double Mean_C(NumericVector x){
  int n = x.size();
  return Sum_C(x)/n;
}

// [[Rcpp::export]]
double Mean_C_NA(NumericVector x){
  return Sum_C_NA(x)/length_C_NA(x);
}

// [[Rcpp::export]]
IntegerVector Which_C(double x, NumericVector x_set){
  int N = x_set.size();
  IntegerVector out_Temp(N);
  int count = 0;
  for(int j = 0; j < N; ++j){
    if(x == x_set[j]){
      count = count + 1;
      out_Temp[count - 1] = j;
    }
  }
  if(count == 0){
    return(0);
  }else
  {
    IntegerVector out(count);
    for(int i = 0; i < count; ++i){
      out[i] = out_Temp[i];
    }
    return out;
  }  
}


// [[Rcpp::export]]
IntegerVector Which_C_NA(double x, NumericVector x_set){
  if(NumericVector::is_na(x)){
    IntegerVector out(1);
    out = NA_INTEGER;
    return out;
  }else
  {
    if(all(is_na(x_set)) ){
      IntegerVector out(1);
      out = NA_INTEGER;
      return out;
    }else
    {
      int N = x_set.size();
      IntegerVector out_Temp(N);
      int count = 0;
      for(int j = 0; j < N; ++j){
        if(x == x_set[j]){
          count = count + 1;
          out_Temp[count - 1] = j;
        }
      }
      if(count == 0){
        return(0);
      }else
      {
        IntegerVector out(count);
        for(int i = 0; i < count; ++i){
          out[i] = out_Temp[i];
        }
        return out;
      } 
    }
  }
}


// [[Rcpp::export]]
int Which_Min_C(NumericVector x){
  int Which_Min;
  int n = x.size();
  if(n > 1){
    if(is_true(all(x == x[0]))){
      Which_Min = -1;
    }else
    {
      int n = x.size();
      double Min_Value;
      int count = 0;
      for(int i = 0; i < n; ++i){
        if(NumericVector::is_na(x[i])){
          continue;
        }
        else
        {
          if(count == 0){
            Min_Value = x[i];
            count = count + 1;
          }
          if( x[i] <= Min_Value){
            Min_Value = x[i];
            Which_Min = i;
          }
        }
      }
    }
  }else
  {
    Which_Min = 0;
  }
  return Which_Min;
}


// [[Rcpp::export]]
int Which_Min_C_NA(NumericVector x){
  int Which_Min;
  if(all(is_na(x))){
    Which_Min = NA_INTEGER;
  }else
  {
    int n = x.size();
    if(n > 1){
      if(is_true(all(x == x[0]))){
        Which_Min = -1;
      }else
      {
        int n = x.size();
        double Min_Value;
        int count = 0;
        for(int i = 0; i < n; ++i){
          if(NumericVector::is_na(x[i])){
            continue;
          }
          else
          {
            if(count == 0){
              Min_Value = x[i];
              count = count + 1;
            }
            if( x[i] <= Min_Value){
              Min_Value = x[i];
              Which_Min = i;
            }
          }
        }
      }
    }else
    {
      Which_Min = 0;
    }
  }
  return Which_Min;
}

// [[Rcpp::export]]
int Which_Max_C(NumericVector x){
  int Which_Max;
  int n = x.size();
  if(n > 1){
    if(is_true(all(x == x[0]))){
      Which_Max = -1;
    }else
    {
      double Max_Value;
      int count = 0;
      for(int i = 0; i < n; ++i){
        if(NumericVector::is_na(x[i])){
          continue;
        }
        else
        {
          if(count == 0){
            Max_Value = x[i];
            count = count + 1;
          }
          if( x[i] >= Max_Value){
            Max_Value = x[i];
            Which_Max = i;
          }
        }
      }    
    }
  }else
  {
    Which_Max = 0;
  }
  return Which_Max;
}


// [[Rcpp::export]]
int Which_Max_C_NA(NumericVector x){
  int Which_Max;
  if(all(is_na(x))){
    Which_Max = NA_INTEGER;
  }else
  {
    int n = x.size();
    if(n > 1){
      if(is_true(all(x == x[0]))){
        Which_Max = -1;
      }else
      {
        double Max_Value;
        int count = 0;
        for(int i = 0; i < n; ++i){
          if(NumericVector::is_na(x[i])){
            continue;
          }
          else
          {
            if(count == 0){
              Max_Value = x[i];
              count = count + 1;
            }
            if( x[i] >= Max_Value){
              Max_Value = x[i];
              Which_Max = i;
            }
          }
        }    
      }
    }else
    {
      Which_Max = 0;
    } 
  }
  return Which_Max;
}


// [[Rcpp::export]]
List StdVar_C(NumericMatrix MyMat){
  int n = MyMat.nrow();
  int K = MyMat.ncol();
  NumericMatrix NewMyMat(n,K);
  NumericVector x_Mean(K);
  double x_Std_Error_Temp;
  NumericVector x_Std_Error(K);
  for(int k = 0; k < K; ++k){
    NumericVector Column_k(n);
    for(int i = 0; i < n; ++i){
      Column_k[i] = MyMat(i,k);
    }
    x_Mean[k] = Mean_C(Column_k);
    x_Std_Error_Temp = sqrt(Sum_C( pow( Column_k - x_Mean[k],2.0  ) ));
    if(x_Std_Error_Temp == 0){
      x_Std_Error[k] = 1;
    }else
    {
      x_Std_Error[k] = x_Std_Error_Temp;
    }
    for(int i = 0; i < n; ++i){
      NewMyMat(i,k) = (Column_k[i] - x_Mean[k])/x_Std_Error[k];
    }
  }
  return List::create(
    _["Std_Matrix"] = NewMyMat, 
    _["Std_Mean"] = x_Mean, 
    _["Std_Error"] = x_Std_Error
  );
}

// [[Rcpp::export]]
List StdVar_C_NA(NumericMatrix MyMat){
  int n = MyMat.nrow();
  int K = MyMat.ncol();
  NumericMatrix NewMyMat(n,K);
  NumericVector x_Mean(K);
  double x_Std_Error_Temp;
  NumericVector x_Std_Error(K);
  for(int k = 0; k < K; ++k){
    NumericVector Column_k(n);
    for(int i = 0; i < n; ++i){
      Column_k[i] = MyMat(i,k);
    }
    x_Mean[k] = Mean_C_NA(Column_k);
    x_Std_Error_Temp = sqrt(Sum_C_NA( pow( Column_k - x_Mean[k],2.0  ) ));
    if(x_Std_Error_Temp == 0){
      x_Std_Error[k] = 1;
    }else
    {
      x_Std_Error[k] = x_Std_Error_Temp;
    }
    for(int i = 0; i < n; ++i){
      NewMyMat(i,k) = (Column_k[i] - x_Mean[k])/x_Std_Error[k];
    }
  }
  return List::create(
    _["Std_Matrix"] = NewMyMat, 
    _["Std_Mean"] = x_Mean, 
    _["Std_Error"] = x_Std_Error
  );
}

// [[Rcpp::export]]
IntegerVector Match_C(NumericVector x_subset, NumericVector x_set){
  int n = x_subset.size();
  int N = x_set.size();
  IntegerVector out(n);
  for(int ii = 0; ii < n; ++ii){
    out[ii] = NA_INTEGER;
  }
  for(int i = 0; i < n; ++i){
    double x_subset_Temp = x_subset[i];
    for(int j = 0; j < N; ++j){
      if(x_subset_Temp == x_set[j]){
        out[i] = j;
        break;
      }
    }
  }
  return out;
}


// [[Rcpp::export]]
IntegerVector Match_C_NA(NumericVector x_subset, NumericVector x_set){
  int n = x_subset.size();
  int N = x_set.size();
  IntegerVector out(n);
  for(int ii = 0; ii < n; ++ii){
    out[ii] = NA_INTEGER;
  }
  if(all(is_na(x_subset))){
    return out;
  }else
  {
    if(all(is_na(x_set))){
      return out;
    }else
    {
      for(int i = 0; i < n; ++i){
        double x_subset_Temp = x_subset[i];
        if(NumericVector::is_na(x_subset_Temp) ){
          out[i] = NA_INTEGER;
        }else
        {
          for(int j = 0; j < N; ++j){
            double x_set_Temp = x_set[j];
            if(NumericVector::is_na(x_set_Temp)){
              continue;
            }else{
              if(x_subset_Temp == x_set_Temp){
                out[i] = j;
                break;
              }               
            }
          }
        }
      }
      return out;
    }
  }
}


// [[Rcpp::export]]
IntegerVector Approx_Match_C(NumericVector x, NumericVector y){
  int n = x.size();
  IntegerVector out(n);
  for(int i = 0; i < n; ++i){
    out[i] = Which_Min_C( abs(x[i] - y) );
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector Approx_Match_C_NA(NumericVector x, NumericVector y){
  int n = x.size();
  IntegerVector out(n);
  for(int i = 0; i < n; ++i){
    out[i] = Which_Min_C_NA( abs(x[i] - y) );
  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix Diag_Matrix_C(NumericVector x){
  int n = x.size();
  NumericMatrix x_Mat(n,n);
  for(int i = 0;i < n;++i){
    for(int j = 0; j < n; ++j){
      if(i == j){
        x_Mat(i,j) = x[i];
      }else
      {
        x_Mat(i,j) = 0; 
      } 
    }
  }
  return x_Mat;
}


// [[Rcpp::export]]
NumericMatrix Matrix_Sum_C(NumericMatrix x, NumericMatrix y){
  int Nr = x.nrow();
  int Nc = x.ncol();
  if(Nr != y.nrow() || Nc != y.ncol() ){
    stop("Dimensions do not match");
  }
  NumericMatrix Matrix_Sum(Nr,Nc);
  for(int i = 0; i < Nr; ++i){
    for(int j = 0; j < Nc; ++j){
      Matrix_Sum(i,j) = x(i,j) + y(i,j); 
    }
  }
  return Matrix_Sum;
}

// [[Rcpp::export]]
NumericMatrix Matrix_Sum_C_NA(NumericMatrix x, NumericMatrix y){
  int Nr = x.nrow();
  int Nc = x.ncol();
  if(Nr != y.nrow() || Nc != y.ncol() ){
    stop("Dimensions do not match");
  }
  NumericMatrix Matrix_Sum(Nr,Nc);
  for(int i = 0; i < Nr; ++i){
    for(int j = 0; j < Nc; ++j){
      double x_Temp = x(i,j);
      double y_Temp = y(i,j);
      if( NumericVector::is_na(x_Temp) ){
        x_Temp = 0;
      }
      if( NumericVector::is_na(y_Temp) ){
        y_Temp = 0;
      }
      Matrix_Sum(i,j) = x_Temp + y_Temp; 
    }
  }
  return Matrix_Sum;
}


// [[Rcpp::export]]
double l2Dist_Vector_C(NumericVector x1, NumericVector x2, List ID){
  int n = ID.size();
  NumericVector x_Mean(n);
  for(int i = 0; i < n;++i){
    IntegerVector ID_n = ID[i];
    int ni = ID_n.size();
    NumericVector x1_Temp(ni);
    NumericVector x2_Temp(ni);
    for(int j = 0; j < ni; ++j){
      x1_Temp[j] = x1[ ID_n[j]  ];
      x2_Temp[j] = x2[ ID_n[j]  ];
    }
    x_Mean[i] = Mean_C(pow(x1_Temp - x2_Temp,2.0));
  }
  return sqrt(Mean_C(x_Mean));
}


// [[Rcpp::export]]
double l2Dist_Vector_C_NA(NumericVector x1, NumericVector x2, List ID){
  double out;
  if( all(is_na(x1)) ){
    out = NA_REAL;
  }else
  {
    if(all(is_na(x2))){
      out = NA_REAL;
    }else
    {
      int n = ID.size();
      NumericVector x_Mean(n);
      for(int i = 0; i < n;++i){
        IntegerVector ID_n = ID[i];
        int ni = ID_n.size();
        NumericVector x1_Temp(ni);
        NumericVector x2_Temp(ni);
        for(int j = 0; j < ni; ++j){
          x1_Temp[j] = x1[ ID_n[j]  ];
          x2_Temp[j] = x2[ ID_n[j]  ];
        }
        x_Mean[i] = Mean_C_NA(pow(x1_Temp - x2_Temp,2.0));
      }
      out = sqrt(Mean_C_NA(x_Mean));
    }
  }
  return out;
}

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a) {
  Rcpp::NumericVector b = Rcpp::clone(a);
  std::random_shuffle(b.begin(), b.end(), randWrapper);
  return b;
}

// [[Rcpp::export]]
Rcpp::IntegerVector int_randomShuffle(Rcpp::IntegerVector a) {
  Rcpp::IntegerVector b = Rcpp::clone(a);
  std::random_shuffle(b.begin(), b.end(), randWrapper);
  return b;
}

// [[Rcpp::export]]
NumericVector RemoveNA(NumericVector x){
  int n = x.size();
  if(any(is_na(x))){
    int n_new = length_C_NA(x);
    NumericVector x_new(n_new);
    int count = 0;
    for(int i = 0; i < n;++i){
      if(NumericVector::is_na(x[i]) ){
        continue; 
      }else
      {
        x_new[count] = x[i];
        count = count + 1;
      }
    }
    return x_new;
  }else
  {
    return x;
  }
}


// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
NumericVector stl_sort_NA(NumericVector x) {
  NumericVector x_new = RemoveNA(x);
  NumericVector y = clone(x_new);
  std::sort(y.begin(), y.end());
  return y;
}


// [[Rcpp::export]]
Rcpp::NumericVector stl_sort_reverse(Rcpp::NumericVector x) {
  Rcpp::NumericVector y = Rcpp::clone(x);
  y.sort(true);
  return y;
}

// [[Rcpp::export]]
Rcpp::NumericVector stl_sort_reverse_NA(Rcpp::NumericVector x) {
  NumericVector x_new = RemoveNA(x);
  Rcpp::NumericVector y = Rcpp::clone(x_new);
  y.sort(true);
  return y;
}


// [[Rcpp::export]]
NumericVector sort_unique_C(NumericVector x) {
 NumericVector out = sort_unique(x);
  return out;
}

// [[Rcpp::export]]
NumericVector sort_unique_C_NA(NumericVector x) {
  NumericVector x_new = RemoveNA(x);
  NumericVector out = sort_unique(x_new);
  return out;
}


// [[Rcpp::export]]
IntegerVector Which_Max_Matrix(NumericMatrix x) {
  IntegerVector Output(2); 
  int N_R = x.nrow();
  int N_C = x.ncol();
  NumericVector Max_Col(N_C);
  for(int j = 0; j < N_C; ++j){
    NumericVector x_Col(N_R); 
    for(int i = 0; i < N_R; ++i){
      x_Col[i] = x(i,j);
    }
    Max_Col[j] = Which_Max_C(x_Col);
  }
  if(is_true(all(Max_Col == -1) )){
    Output[0] = -1;
    Output[1] = -1;
  }else
  {
    NumericVector x_Col_Max(N_C); 
    for(int jj = 0; jj < N_C; ++jj){
      x_Col_Max[jj] = x( Max_Col[jj] ,jj);
    }
    int Max_Col_update = Which_Max_C(x_Col_Max);
    Output[0] = Max_Col[Max_Col_update];
    Output[1] = Max_Col_update;
  }
  return Output;
}


// [[Rcpp::export]]
IntegerVector Which_Max_Matrix_NA(NumericMatrix x) {
  IntegerVector Output(2); 
  int N_R = x.nrow();
  int N_C = x.ncol();
  IntegerVector Max_Col(N_C);
  for(int j = 0; j < N_C; ++j){
    NumericVector x_Col(N_R); 
    for(int i = 0; i < N_R; ++i){
      x_Col[i] = x(i,j);
    }
    Max_Col[j] = Which_Max_C_NA(x_Col);
  }
  if(is_true(all(Max_Col == -1) )){
    Output[0] = -1;
    Output[1] = -1;
  }else
  {
    NumericVector x_Col_Max(N_C); 
    for(int jj = 0; jj < N_C; ++jj){
      if(NumericVector::is_na(Max_Col[jj])){
        x_Col_Max[jj] = NA_REAL;
      }else
      {
        x_Col_Max[jj] = x( Max_Col[jj] ,jj);  
      }
    }
    int Max_Col_update = Which_Max_C_NA(x_Col_Max);
    Output[0] = Max_Col[Max_Col_update];
    Output[1] = Max_Col_update;
  }
  return Output;
}

// [[Rcpp::export]]
NumericVector rowSums_C(NumericMatrix x){
  int N_R = x.nrow();
  int N_C = x.ncol();
  NumericVector output(N_R);
  for(int i = 0; i < N_R; ++i){
    NumericVector Col_i(N_C);
    for(int j = 0; j < N_C; ++j){
      Col_i[j] = x(i,j);
    }
    output[i] = Sum_C(Col_i);
  }
  return output;
}

// [[Rcpp::export]]
NumericVector rowSums_C_NA(NumericMatrix x){
  int N_R = x.nrow();
  int N_C = x.ncol();
  NumericVector output(N_R);
  for(int i = 0; i < N_R; ++i){
    NumericVector Col_i(N_C);
    for(int j = 0; j < N_C; ++j){
      Col_i[j] = x(i,j);
    }
    output[i] = Sum_C_NA(Col_i);
  }
  return output;
}

// [[Rcpp::export]]
LogicalVector isNA(IntegerVector x) {
  int n = x.size();
  LogicalVector out(n);
  
  for (int i = 0; i < n; ++i) {
    out[i] = IntegerVector::is_na(x[i]);
  }
  return out;
}

// [[Rcpp::export]]
double Rho_Inv_C(double Rho_Value,double N_Value){
  double eps = 0.01;
  double N_Rho = (N_Value - 1.0);
  double out = 0.0;
  if(N_Rho == 0){
    out = 0.0;
  }else
  {
    double rho_inv = 0.0;
    rho_inv = Rho_Value + (1/N_Rho);
    if(Rho_Value < 0 && std::abs(rho_inv) < eps){
      out = ((-1/N_Rho) + eps)/(N_Rho*eps);
    }else
    {
      out = Rho_Value/(1 + (N_Rho*Rho_Value));
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat MatrixInversion_Equicorrelation_C(int N_Value, double phi, double rho){
  arma::mat Output;
  arma::mat Ii(N_Value,N_Value);
  Ii.fill(0);
  for(int j = 0; j < N_Value; ++j){
    Ii(j,j) = 1;
  }
  arma::mat Ji(N_Value,N_Value);
  Ji.fill(1);
  double N_Value_double = 1.0*N_Value;
  
  double rho_inv = Rho_Inv_C(rho,N_Value_double);
  arma::mat rho_inv_matrix(N_Value,N_Value);
  rho_inv_matrix.fill(rho_inv);
  
  double phi_rho_inv = 1/((1-rho)*phi);
  arma::mat phi_rho_inv_matrix(N_Value,N_Value);
  phi_rho_inv_matrix.fill(phi_rho_inv);
  
  Output = phi_rho_inv_matrix%((Ii - (rho_inv_matrix%Ji) ));
  return Output;
}

// [[Rcpp::export]]
NumericVector Matrix_Vector_Multiplication_C(NumericVector x,arma::mat y){
  int n = x.size();
  int y_row = y.n_rows;
  int y_col = y.n_cols;
  if(y_row != y_col){
    stop("Matrix must be a symmetric square");
  }
  if(n != y_col){
    stop("Number of columns differ from a length of the vector");
  }
  NumericVector Output(n);
  for(int i = 0; i < n; ++i){
    NumericVector Output_Temp(n);
    for(int j = 0; j < n; ++j){
      Output_Temp[j] = x[j]*y(j,i);
    }
    Output[i] = Sum_C_NA(Output_Temp);
  }
  return Output;
}

// [[Rcpp::export]]
List DataProcessing_C(NumericMatrix Org_x,
                      NumericMatrix Org_y,
                      NumericVector id, 
                      NumericVector tm,
                      bool x_miss)
{
  NumericVector unq_id = sort_unique_C_NA(id);
  int n = unq_id.size();
  
  List id_index(n);
  for(int i = 0; i < n; ++i){
    id_index[i] = Which_C_NA(unq_id[i],id);
  }
  
  IntegerVector ni(n);
  for(int i = 0; i < n; ++i){
    IntegerVector id_index_Temp = id_index[i];
    ni[i] = id_index_Temp.size();
  }
  
  int N = Org_y.nrow();
  int L = Org_y.ncol();
  int K = Org_x.ncol();
  
  IntegerVector unlist_id_index(N);
  int count = 0;
  for(int i = 0; i < n; ++i){
    NumericVector id_index_Temp = id_index[i];
    for(int j = 0; j < ni[i]; ++j){
      unlist_id_index[count] = id_index_Temp[j];
      count = count + 1;
    }
  }
  
  NumericVector id_New(N);
  for(int ii = 0; ii < N; ++ii){
    id_New[ii] = id[ unlist_id_index[ii] ];
  }
  
  LogicalVector id_Match(N);
  for(int ii = 0; ii < N; ++ii){
    id_Match[ii] = (id[ii] == id_New[ii]);
  }
  
  if(is_false(all( id_Match == TRUE) )){
    NumericVector tm_New(N);
    NumericMatrix Org_x_New(N,K);
    NumericMatrix Org_y_New(N,L);
    for(int ii = 0; ii < N; ++ii){
      tm_New[ii] = tm[ unlist_id_index[ii]  ];
      id_New[ii] = id[ unlist_id_index[ii]  ];
      for(int k = 0; k < K; ++k){
        Org_x_New(ii,k) = Org_x( unlist_id_index[ii], k  );
      }
      for(int l = 0; l < L; ++l){
        Org_y_New(ii,l) = Org_y( unlist_id_index[ii], l  );
      }
    }
    tm = tm_New;
    id = id_New;
    Org_x = Org_x_New;
    Org_y = Org_y_New;  
  }
  
  NumericMatrix x(N,K);
  NumericVector x_Mean(K);
  NumericVector x_Std_Error(K);
  if(x_miss){
    x = Org_x;
    x_Mean = 0.0;
    x_Std_Error = 1.0;
  }
  else
  {
    List Std_x = StdVar_C_NA(Org_x);
    x = as<NumericMatrix>(Std_x["Std_Matrix"]);
    x_Mean = as<NumericVector>(Std_x["Std_Mean"]);
    x_Std_Error = as<NumericVector>(Std_x["Std_Error"]);
  }
  
  List Std_y = StdVar_C_NA(Org_y);
  NumericMatrix y = as<NumericMatrix>(Std_y["Std_Matrix"]);
  NumericVector y_Mean = as<NumericVector>(Std_y["Std_Mean"]);
  NumericVector y_Std_Error = as<NumericVector>(Std_y["Std_Error"]);
  
  List Data = List::create(
    _["Org_x"] = Org_x, 
    _["Org_y"] = Org_y,
    _["id"] = id,
    _["tm"] = tm,
    _["x"] = x, 
    _["y"] = y,
    _["x_Mean"] = x_Mean,
    _["x_Std_Error"] = x_Std_Error,
    _["y_Mean"] = y_Mean,
    _["y_Std_Error"] = y_Std_Error
  );
  List Dimensions = List::create(
    _["n"] = n,
    _["K"] = K,
    _["L"] = L,
    _["ni"] = ni,
    _["N"] = N
  );
  
  List Index = List::create(
    _["unq_id"] = unq_id,
    _["id_index"] = id_index
  );
  
  return List::create(
    _["Data"] = Data,
    _["Dimensions"] = Dimensions,
    _["Index"] = Index
  );
}
  
// [[Rcpp::export]]
List BoostMLR_C(NumericMatrix Org_x,
                NumericMatrix Org_y,
                NumericVector id,
                NumericVector tm,
                NumericMatrix x,
                NumericMatrix y,
                NumericVector x_Mean,
                NumericVector x_Std_Error,
                NumericVector y_Mean,
                NumericVector y_Std_Error,
                int n,
                int K,
                int L,
                int H,
                IntegerVector Dk,
                IntegerVector ni,
                int N,
                NumericVector unq_id,
                NumericVector unq_tm,
                List unq_x,
                List id_index,
                NumericMatrix Bt,
                List Bx,
                List Bx_Scale,
                NumericMatrix Time_Add_New,
                LogicalVector Time_Unmatch,
                double nu,
                int M,
                bool Mod_Grad,
                LogicalVector UseRaw,
                NumericVector Lambda_Ridge_Vec,
                bool Ridge_Penalty,
                bool Shrink,
                double lower_perc,
                double upper_perc,
                double Lambda_Scale,
                int NLambda,
                bool VarFlag,
                NumericVector rho,
                NumericVector phi,
                bool Verbose)
  
{
  List tm_index(n);
  for(int i = 0; i < n; ++i){
    NumericVector tm_Temp(ni[i]);
    IntegerVector id_index_Temp = id_index[i];
    for(int j = 0; j < ni[i]; ++j){
      tm_Temp[j] = tm[id_index_Temp[j]];
    }
    tm_index[i] = Match_C_NA(tm_Temp,unq_tm);
  }

  List Bt_H(H);
  for(int h = 0; h < H; ++h){
    List Bt_n(n);
    for(int i = 0; i < n; ++i){
      IntegerVector tm_index_Temp = tm_index[i];
      NumericVector Bt_ni(ni[i]);
      for(int j = 0; j < ni[i]; ++j){
        if(IntegerVector::is_na( tm_index_Temp[j])){
          Bt_ni[j] = NA_REAL;
        }else
        {
          Bt_ni[j] = Bt(tm_index_Temp[j], h);  
        }
      }
      Bt_n[i] = Bt_ni;
    }
    Bt_H[h] = Bt_n;
  }

  List x_index(K);
  for(int k = 0; k < K; ++k){
    if(UseRaw[k]){
      x_index[k] = NA_INTEGER;
    }else
    {
      NumericVector unq_x_Temp = unq_x[k];
      List x_index_Subject(n);
      for(int i = 0; i < n; ++i){
        NumericVector x_Temp(ni[i]);
        IntegerVector id_index_Temp = id_index[i];
        for(int j = 0; j < ni[i]; ++j){
          x_Temp[j] = x(id_index_Temp[j] , k);
        }
        x_index_Subject[i] = Match_C_NA(x_Temp,unq_x_Temp);
      }
      x_index[k] = x_index_Subject;
    }
  }
  
  int count = -1;
  List Time_Add_index(K);
  for(int k = 0; k < K; ++k){
    if(Time_Unmatch[k]){
      count = count + 1;
      List Time_Add_index_Subject(n);
      for(int i = 0; i < n; ++i){
        NumericVector Time_Add_Temp(ni[i]);
        IntegerVector id_index_Temp = id_index[i];
        for(int j = 0; j < ni[i]; ++j){
          Time_Add_Temp[j] = Time_Add_New(id_index_Temp[j],count); 
        }
        Time_Add_index_Subject[i] = Approx_Match_C_NA(Time_Add_Temp,unq_tm);
      }
      Time_Add_index[k] = Time_Add_index_Subject;
     }
      else {
      Time_Add_index[k] = NA_INTEGER;
    }
  }
  
  IntegerVector DkT(K);
  for(int k = 0; k < K; ++k){
    if(UseRaw[k]){
      DkT[k] = Dk[k];
    }else
    {
      if(Time_Unmatch[k]){
        DkT[k] = Dk[k]*H;
      }else
      {
        DkT[k] = Dk[k];  
      }
    }
  }
  
  List Bx_K(K);
  for(int k = 0; k < K; ++k){
    NumericMatrix Bx_K_Temp = Bx[k];
    List Bx_Dk(DkT[k]);
    if(UseRaw[k]){
      for(int d = 0; d < DkT[k]; ++d){
        List Bx_n(n);
        for(int i = 0; i < n; ++i){
          IntegerVector id_index_Temp = id_index[i];
          NumericVector Bx_ni(ni[i]);
          for(int j = 0; j < ni[i]; ++j){
            Bx_ni[j] = Bx_K_Temp( id_index_Temp[j],d) ;
          }
          Bx_n[i] = Bx_ni;
        }
        Bx_Dk[d] = Bx_n;
      }
    }
    else{
      count = -1;
      if(Time_Unmatch[k]){
        List x_index_Temp = x_index[k];
        List Time_Add_index_Temp = Time_Add_index[k];
          for(int d = 0; d < Dk[k]; ++d){
            for(int h = 0; h < H; ++h){
              count = count + 1;
              List Bx_n(n);
              for(int i = 0; i < n; ++i){
                IntegerVector x_index_Temp_n = x_index_Temp[i];
                IntegerVector Time_Add_index_Temp_n = Time_Add_index_Temp[i];
                NumericVector Bx_ni(ni[i]);
                for(int j = 0; j < ni[i]; ++j){
                  if(IntegerVector::is_na( x_index_Temp_n[j])){
                    Bx_ni[j] = NA_REAL;
                  }else
                  {
                    if(IntegerVector::is_na( Time_Add_index_Temp_n[j])){
                      Bx_ni[j] = NA_REAL;
                    }else
                    {
                      Bx_ni[j] = Bx_K_Temp(x_index_Temp_n[j],d)*Bt(Time_Add_index_Temp_n[j],h);
                    }
                  }
                }
                Bx_n[i] = Bx_ni;
              }
              Bx_Dk[count] = Bx_n;
            }
          } 
      }
      else {
        List x_index_Temp = x_index[k];
        for(int d = 0; d < DkT[k]; ++d){
          List Bx_n(n);
          for(int i = 0; i < n; ++i){
            IntegerVector x_index_Temp_n = x_index_Temp[i];
            NumericVector Bx_ni(ni[i]);
            for(int j = 0; j < ni[i]; ++j){
              if(IntegerVector::is_na( x_index_Temp_n[j])){
                Bx_ni[j] = NA_REAL;
              }else
              {
                Bx_ni[j] = Bx_K_Temp(x_index_Temp_n[j],d);
              }
            }
            Bx_n[i] = Bx_ni;
          }
          Bx_Dk[d] = Bx_n;
        }
      }
    }
    Bx_K[k] = Bx_Dk;
  }
  
  Dk = DkT;
  
  List Bxt(K);
  for(int k = 0; k < K; ++k){
    List Bx_Temp_K = Bx_K[k];
    List Bxt_Dk(Dk[k]);
    for(int d = 0; d < Dk[k]; ++d){
      List Bx_Temp_Dk = Bx_Temp_K[d];
      List Bxt_H(H);
      for(int h = 0; h < H; ++h){
        List Bt_Temp_H = Bt_H[h];
        List Bxt_n(n);
        for(int i = 0; i < n; ++i){
          NumericVector Bx_Temp_n = Bx_Temp_Dk[i];
          NumericVector Bt_Temp_n = Bt_Temp_H[i];
          NumericVector Bxt_ni(ni[i]);
          for(int j = 0; j < ni[i]; ++j){
            Bxt_ni[j] = Bx_Temp_n[j]*Bt_Temp_n[j];
          }
          Bxt_n[i] = Bxt_ni;
        }
        Bxt_H[h] = Bxt_n;
      }
      Bxt_Dk[d] = Bxt_H;
    }
    Bxt[k] = Bxt_Dk;
  }
  
  NumericMatrix mu_zero(N,L);
  for(int l = 0; l < L; ++l){
    for(int i = 0; i < N; ++i){
      mu_zero(i,l) = 0;
    }
  }
  
  NumericMatrix mu(N,L);
  mu = mu_zero;
  
  NumericVector Vec_zero(N);
  for(int i = 0; i < N; ++i){
    Vec_zero[i] = 0.0;
  }
  
  List Beta(K);
  for(int k = 0; k < K; ++k){
    List Beta_K(Dk[k]);
    for(int d = 0; d < Dk[k]; ++d){
      List Beta_Dk(H);
      for(int h = 0; h < H; ++h){
        NumericVector Beta_H(L);
        for(int l = 0; l < L; ++l){
          Beta_H[l] = 0.0;
        }
        Beta_Dk[h] = Beta_H;
      }
      Beta_K[d] = Beta_Dk;
    }
    Beta[k] = Beta_K;
  }

  NumericMatrix Error_Rate(M,L);
  IntegerMatrix Variable_Select(M,H);
  IntegerMatrix Response_Select(M,H);
  List Beta_Hat_List(M);
  List Sum_Beta_Hat_List(M);
  List Beta_Hat_List_Iter(M);
  List mu_List(M);
  List Lambda_List(M);
  NumericMatrix Phi(M,L);
  NumericMatrix Rho(M,L);
  
  for(int m = 0; m < M; ++m){
    for(int h = 0; h < H; ++h){
      Variable_Select(m,h) = NA_INTEGER;
      Response_Select(m,h) = NA_INTEGER;
    }
    Beta_Hat_List[m] = NA_REAL;
    Sum_Beta_Hat_List[m] = NA_REAL;
    Beta_Hat_List_Iter[m] = NA_REAL;
    mu_List[m] = NA_REAL;
    Lambda_List[m] = NA_REAL;
    for(int l = 0; l < L; ++l){
      Error_Rate(m,l) = NA_REAL;
      Phi(m,l) = NA_REAL;
      Rho(m,l) = NA_REAL;
    }
  }
  List List_Trace_Bxt_gm(M);
  List lower_Beta_Hat_Noise(K);
  List upper_Beta_Hat_Noise(K);
  List Beta_Hat_Old(K);
  Beta_Hat_Old = Beta;
  List V_inv(L);
  
  // Boosting iteration starts here...  
  
  for(int m = 0; m <  M; ++m){
    
    if(Verbose){
      if( (m+1)%(M/10) == 0.0){
        double FracBoosting;
        FracBoosting = ( (m+1)*100)/M;
        int Perc_boosting = FracBoosting;
         Rcout << Perc_boosting << "%" << " ";
      }
    }
    
    if(m == 0 || VarFlag == TRUE){
      for(int l = 0; l < L; ++l){
        List Vi_inv(n);
        for(int i = 0; i < n;++i){
          arma::mat Vi_inv_Temp;
          Vi_inv_Temp = MatrixInversion_Equicorrelation_C(ni[i],phi[l],rho[l]);
          Vi_inv[i] = Vi_inv_Temp;
        }
        V_inv[l] = Vi_inv;
      }
    }
    
    List gm(L);
    for(int l = 0; l < L; ++l){
      List gm_n(n);
      List V_inv_l = V_inv[l];
      for(int i = 0; i < n; ++i){
        arma::mat V_inv_i = V_inv_l[i];
        IntegerVector id_index_Temp = id_index[i];
        NumericVector gm_ni(ni[i]);
        for(int j = 0; j < ni[i]; ++j){
          gm_ni[j] = y(id_index_Temp[j],l) - mu(id_index_Temp[j],l);
        }
        if(Mod_Grad){
          gm_n[i] = gm_ni;  
        }
        else {
          gm_n[i] = Matrix_Vector_Multiplication_C(gm_ni,V_inv_i);
        }
      }
      gm[l] = gm_n;
    }

    List Temp_Trace_Bxt_gm(K);
    List Lambda(K);
    List Beta_Hat(K);
    if(m > 0){
      for(int k = 0; k < K; ++k){
        List Beta_Hat_Old_Temp_Dk = Beta_Hat_Old[k];
        List Bxt_Temp_K = Bxt[k];
        NumericVector Bx_Scale_K = Bx_Scale[k];
        List Lambda_Dk(Dk[k]);
        List Beta_Hat_Dk(Dk[k]);
        List Temp_Trace_Bxt_gm_Dk(Dk[k]);
        for(int d = 0; d < Dk[k]; ++d){
          List Beta_Hat_Old_Temp_H = Beta_Hat_Old_Temp_Dk[d];
          List Bxt_Temp_Dk = Bxt_Temp_K[d];
          List Lambda_H(H);
          List Beta_Hat_H(H);
          List Temp_Trace_Bxt_gm_H(H);
          for(int h = 0; h < H; ++h){
            NumericVector Beta_Hat_Old_Temp_L = Beta_Hat_Old_Temp_H[h];
            List Bxt_Temp_H = Bxt_Temp_Dk[h];
            NumericVector Lambda_L(L);
            NumericVector Beta_Hat_L(L);
            NumericVector Temp_Trace_Bxt_gm_L(L);
            for(int l = 0; l < L; ++l){
              double Beta_Hat_Old_Temp = Beta_Hat_Old_Temp_L[l];
              if(Beta_Hat_Old_Temp == 0.0){
                Beta_Hat_L[l] = 0.0;
              } else
              {
                List V_inv_l = V_inv[l];
                List gm_Temp_L = gm[l];
                NumericVector Bxt_n_2(n);
                NumericVector Bxt_gm_n(n);
                NumericVector gm_n_2(n);
                for(int i = 0; i < n; ++i){
                  arma::mat V_inv_i = V_inv_l[i];
                  NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                  NumericVector gm_Temp_n = gm_Temp_L[i];
                  NumericVector Bxt_V = Matrix_Vector_Multiplication_C(Bxt_Temp_n,V_inv_i);
                  NumericVector Bxt_ni_2(ni[i]);
                  NumericVector Bxt_gm_ni(ni[i]);
                  NumericVector gm_ni_2(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    Bxt_ni_2[j] = Bxt_V[j] * Bxt_Temp_n[j];
                    Bxt_gm_ni[j] = Bxt_V[j] * gm_Temp_n[j];
                    gm_ni_2[j] = gm_Temp_n[j] * gm_Temp_n[j];
                  }
                  Bxt_n_2[i] = Sum_C_NA(Bxt_ni_2);
                  Bxt_gm_n[i] = Sum_C_NA(Bxt_gm_ni);
                  gm_n_2[i] = Sum_C_NA(gm_ni_2);
                }
                double Trace_Bxt    = Sum_C_NA(Bxt_n_2);
                double Trace_Bxt_gm = Sum_C_NA(Bxt_gm_n);
                Temp_Trace_Bxt_gm_L[l] = Trace_Bxt_gm;
                if(Trace_Bxt < 0.001 ){
                  Beta_Hat_L[l] = 0.0;
                }else
                {
                  if(Ridge_Penalty){
                 /* 
                 Come back to this later since I don't know how to calculate
                  (Vi)^{-1/2} require for calculating lambda for Ridge penalty.
                  Therefore, default setting Ridge_Penalty = FALSE should not be changed.
                  Trace_gm specified below will be incorrect because it does not include Vi_inv matrix.
                 */
                    double Trace_gm       = Sum_C_NA(gm_n_2);
                    double Beta_Hat_NL    = Beta_Hat_Old_Temp;
                    double Trace_eps      = Trace_gm + (Beta_Hat_NL * Beta_Hat_NL * Trace_Bxt) - (2 * Beta_Hat_NL * Trace_Bxt_gm);
                    double Lambda_L_Temp  = (Trace_Bxt/(Trace_gm - Trace_eps))*Lambda_Scale;
                    if(Lambda_L_Temp > 5000){
                      Lambda_L[l] = 5000;
                    }else
                    {
                      if(Lambda_L_Temp < 0){
                        Lambda_L[l] = Lambda_Ridge_Vec[k];
                      }else
                      {
                        Lambda_L[l] = Lambda_L_Temp;  
                      }
                    }
                  }else
                  {
                    Lambda_L[l]         = Lambda_Ridge_Vec[k];   
                  }
                  Beta_Hat_L[l]         = (Trace_Bxt_gm/(Trace_Bxt + Lambda_L[l]))/Bx_Scale_K[d];
                }
              }
            }
            Lambda_H[h] = Lambda_L;
            Beta_Hat_H[h] = Beta_Hat_L; 
            Temp_Trace_Bxt_gm_H[h] = Temp_Trace_Bxt_gm_L;
          }
          Lambda_Dk[d] = Lambda_H;
          Beta_Hat_Dk[d] = Beta_Hat_H;
          Temp_Trace_Bxt_gm_Dk[d] = Temp_Trace_Bxt_gm_H;
        }
        Lambda[k] = Lambda_Dk;
        Beta_Hat[k] = Beta_Hat_Dk;
        Temp_Trace_Bxt_gm[k] = Temp_Trace_Bxt_gm_Dk;
      }
    }else
    {
      for(int k = 0; k < K; ++k){
        List Bxt_Temp_K = Bxt[k];
        NumericVector Bx_Scale_K = Bx_Scale[k];
        List Lambda_Dk(Dk[k]);
        List Beta_Hat_Dk(Dk[k]);
        List Temp_Trace_Bxt_gm_Dk(Dk[k]);
        for(int d = 0; d < Dk[k]; ++d){
          List Bxt_Temp_Dk = Bxt_Temp_K[d];
          List Lambda_H(H);
          List Beta_Hat_H(H);
          List Temp_Trace_Bxt_gm_H(H);
          for(int h = 0; h < H; ++h){
            List Bxt_Temp_H = Bxt_Temp_Dk[h];
            NumericVector Lambda_L(L);
            NumericVector Beta_Hat_L(L);
            NumericVector Temp_Trace_Bxt_gm_L(L);
            for(int l = 0; l < L; ++l){
              List V_inv_l = V_inv[l];
              List gm_Temp_L = gm[l];
              NumericVector Bxt_n_2(n);
              NumericVector Bxt_gm_n(n);
              NumericVector gm_n_2(n);
              for(int i = 0; i < n; ++i){
                arma::mat V_inv_i = V_inv_l[i];
                NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                NumericVector gm_Temp_n = gm_Temp_L[i];
                NumericVector Bxt_V = Matrix_Vector_Multiplication_C(Bxt_Temp_n,V_inv_i);
                NumericVector Bxt_ni_2(ni[i]);
                NumericVector Bxt_gm_ni(ni[i]);
                NumericVector gm_ni_2(ni[i]);
                for(int j = 0; j < ni[i]; ++j){
                  Bxt_ni_2[j] = Bxt_V[j] * Bxt_Temp_n[j];
                  Bxt_gm_ni[j] = Bxt_V[j] * gm_Temp_n[j];
                  gm_ni_2[j] = gm_Temp_n[j] * gm_Temp_n[j];
                }
                Bxt_n_2[i] = Sum_C_NA(Bxt_ni_2);
                Bxt_gm_n[i] = Sum_C_NA(Bxt_gm_ni);
                gm_n_2[i] = Sum_C_NA(gm_ni_2);
              }
              double Trace_Bxt    = Sum_C_NA(Bxt_n_2);
              double Trace_Bxt_gm = Sum_C_NA(Bxt_gm_n);
              Temp_Trace_Bxt_gm_L[l] = Trace_Bxt_gm;
              if(Trace_Bxt < 0.001 ){
                Beta_Hat_L[l] = 0.0;
              }else
              {
                if(Ridge_Penalty){
                  /* 
                   Come back to this later since I don't know how to calculate
                   (Vi)^{-1/2} require for calculating lambda for Ridge penalty.
                   Therefore, default setting Ridge_Penalty = FALSE should not be changed.
                   Trace_gm specified below will be incorrect because it does not include Vi_inv matrix.
                   */
                  double Trace_gm       = Sum_C_NA(gm_n_2);
                  double Beta_Hat_NL    = Trace_Bxt_gm/(Trace_Bxt);
                  double Trace_eps      = Trace_gm + (Beta_Hat_NL * Beta_Hat_NL * Trace_Bxt) - (2 * Beta_Hat_NL * Trace_Bxt_gm);
                  double Lambda_L_Temp  = (Trace_Bxt/(Trace_gm - Trace_eps))*Lambda_Scale;
                  if(Lambda_L_Temp > 5000){
                    Lambda_L[l] = 5000;
                  }else
                  {
                    if(Lambda_L_Temp < 0){
                      Lambda_L[l] = Lambda_Ridge_Vec[k];
                    }else
                    {
                      Lambda_L[l] = Lambda_L_Temp;  
                    }
                  }
                }else
                {
                  Lambda_L[l]         = Lambda_Ridge_Vec[k];
                }
                Beta_Hat_L[l]         = (Trace_Bxt_gm/(Trace_Bxt + Lambda_L[l]))/Bx_Scale_K[d]; 
              }
            }
            Lambda_H[h] = Lambda_L;
            Beta_Hat_H[h] = Beta_Hat_L; 
            Temp_Trace_Bxt_gm_H[h] = Temp_Trace_Bxt_gm_L;
          }
          Lambda_Dk[d] = Lambda_H;
          Beta_Hat_Dk[d] = Beta_Hat_H;
          Temp_Trace_Bxt_gm_Dk[d] = Temp_Trace_Bxt_gm_H;
        }
        Lambda[k] = Lambda_Dk;
        Beta_Hat[k] = Beta_Hat_Dk;
        Temp_Trace_Bxt_gm[k] = Temp_Trace_Bxt_gm_Dk;
      } 
    }
    
    List_Trace_Bxt_gm[m] = Temp_Trace_Bxt_gm;
    Lambda_List[m] = Lambda;
    
    if(Shrink){
      if(m == 0){
        for(int k = 0; k < K; ++k){
          List Bxt_Temp_K = Bxt[k];
          NumericVector Bx_Scale_K = Bx_Scale[k];
          List lower_Beta_Hat_Dk(Dk[k]);
          List upper_Beta_Hat_Dk(Dk[k]);
          for(int d = 0; d < Dk[k]; ++d){
            List Bxt_Temp_Dk = Bxt_Temp_K[d];
            List lower_Beta_Hat_H(H);
            List upper_Beta_Hat_H(H);
            for(int h = 0; h < H; ++h){
              List Bxt_Temp_H = Bxt_Temp_Dk[h];
              NumericVector lower_Beta_Hat_L(L);
              NumericVector upper_Beta_Hat_L(L);
              for(int l = 0; l < L; ++l){
                List V_inv_l = V_inv[l];
                NumericVector Beta_Hat_Noise_Temp(NLambda);
                for(int l_m = 0; l_m < NLambda; ++l_m){
                List gm_Temp_L = gm[l];
                NumericVector gm_unlist(N);
                int count = 0;
                  for(int i = 0; i < n; ++i){
                    NumericVector gm_Temp_n = gm_Temp_L[i];
                    for(int j = 0; j < ni[i]; ++j){
                      gm_unlist[count] = gm_Temp_n[j];
                      count = count + 1;
                    }
                  }
                NumericVector gm_unlist_noise = randomShuffle(gm_unlist);
                  count = 0;
                    List gm_noise_n(n);
                    for(int i = 0; i < n; ++i){
                      NumericVector gm_noise_ni(ni[i]);
                      for(int j = 0; j < ni[i]; ++j){
                        gm_noise_ni[j] = gm_unlist_noise[count];
                        count = count + 1;
                      }
                      gm_noise_n[i] = gm_noise_ni;
                    }
                NumericVector Bxt_n_2(n);
                NumericVector Bxt_gm_n(n);
                NumericVector gm_noise_n_2(n);
                for(int i = 0; i < n; ++i){
                  arma::mat V_inv_i = V_inv_l[i];
                  NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                  NumericVector gm_noise_Temp_n = gm_noise_n[i];
                  NumericVector Bxt_V = Matrix_Vector_Multiplication_C(Bxt_Temp_n,V_inv_i);
                  NumericVector Bxt_ni_2(ni[i]);
                  NumericVector Bxt_gm_ni(ni[i]);
                  NumericVector gm_noise_ni_2(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    Bxt_ni_2[j] = Bxt_V[j] * Bxt_Temp_n[j];
                    Bxt_gm_ni[j] = Bxt_V[j] * gm_noise_Temp_n[j];
                    gm_noise_ni_2[j] = gm_noise_Temp_n[j] * gm_noise_Temp_n[j];
                  }
                  Bxt_n_2[i] = Sum_C_NA(Bxt_ni_2);
                  Bxt_gm_n[i] = Sum_C_NA(Bxt_gm_ni);
                  gm_noise_n_2[i] = Sum_C_NA(gm_noise_ni_2);
                }
                double Trace_Bxt    = Sum_C_NA(Bxt_n_2);
                double Trace_Bxt_gm = Sum_C_NA(Bxt_gm_n);
                if(Trace_Bxt == 0.0){
                  Beta_Hat_Noise_Temp[l_m]  = 0.0;
                }else
                {
                  if(Ridge_Penalty){
                    /* 
                     Come back to this later since I don't know how to calculate
                     (Vi)^{-1/2} require for calculating lambda for Ridge penalty.
                     Therefore, default setting Ridge_Penalty = FALSE should not be changed.
                     Trace_gm specified below will be incorrect because it does not include Vi_inv matrix.
                     */
                    double Trace_gm       = Sum_C_NA(gm_noise_n_2);
                    double Beta_Hat_NL    = Trace_Bxt_gm/(Trace_Bxt);
                    double Trace_eps      = Trace_gm + (Beta_Hat_NL * Beta_Hat_NL * Trace_Bxt) - (2 * Beta_Hat_NL * Trace_Bxt_gm);
                    double Lambda_L_Temp  = (Trace_Bxt/(Trace_gm - Trace_eps))*Lambda_Scale;
                    if(Lambda_L_Temp > 5000){
                      Lambda_L_Temp = 5000;
                    }else
                    {
                      if(Lambda_L_Temp < 0){
                        Lambda_L_Temp = Lambda_Ridge_Vec[k];
                      } //else
                        //{
                        //  Lambda_L_Temp = Lambda_L_Temp;  
                        //}
                    }
                    Beta_Hat_Noise_Temp[l_m]         = (Trace_Bxt_gm/(Trace_Bxt + Lambda_L_Temp))/Bx_Scale_K[d];
                  }else
                  {
                    Beta_Hat_Noise_Temp[l_m]         = (Trace_Bxt_gm/(Trace_Bxt + Lambda_Ridge_Vec[k]))/Bx_Scale_K[d]; 
                  }
                }
               }
                NumericVector sort_Beta_Hat_Noise = stl_sort_NA(Beta_Hat_Noise_Temp);
                int lower_range = (NLambda*lower_perc);
                int upper_range = (NLambda*upper_perc);
                lower_Beta_Hat_L[l] = sort_Beta_Hat_Noise[lower_range];
                upper_Beta_Hat_L[l] = sort_Beta_Hat_Noise[upper_range];
              }
              lower_Beta_Hat_H[h] = lower_Beta_Hat_L;
              upper_Beta_Hat_H[h] = upper_Beta_Hat_L;
            }
            lower_Beta_Hat_Dk[d] = lower_Beta_Hat_H;
            upper_Beta_Hat_Dk[d] = upper_Beta_Hat_H;
          }
          lower_Beta_Hat_Noise[k] = lower_Beta_Hat_Dk;
          upper_Beta_Hat_Noise[k] = upper_Beta_Hat_Dk;
        }
      }

      List Beta_Hat_New(K);
      for(int k = 0; k < K; ++k){
        List Beta_Hat_K = Beta_Hat[k];
        List lower_Beta_Hat_Noise_K = lower_Beta_Hat_Noise[k];
        List upper_Beta_Hat_Noise_K = upper_Beta_Hat_Noise[k];
        List Beta_Hat_New_Dk(Dk[k]);
        for(int d = 0; d < Dk[k]; ++d){
          List Beta_Hat_Dk = Beta_Hat_K[d];
          List lower_Beta_Hat_Noise_Dk = lower_Beta_Hat_Noise_K[d];
          List upper_Beta_Hat_Noise_Dk = upper_Beta_Hat_Noise_K[d];
          List Beta_Hat_New_H(H);
          for(int h = 0; h < H; ++h){
            NumericVector Beta_Hat_H = Beta_Hat_Dk[h];
            NumericVector lower_Beta_Hat_Noise_H = lower_Beta_Hat_Noise_Dk[h];
            NumericVector upper_Beta_Hat_Noise_H = upper_Beta_Hat_Noise_Dk[h];
            NumericVector Beta_Hat_New_L(L);
            for(int l = 0; l < L; ++l){
              double Beta_Hat_L = Beta_Hat_H[l];
              double lower_Beta_Hat_Noise_L = lower_Beta_Hat_Noise_H[l];
              double upper_Beta_Hat_Noise_L = upper_Beta_Hat_Noise_H[l];
              if( (Beta_Hat_L >=  lower_Beta_Hat_Noise_L) && (Beta_Hat_L <= upper_Beta_Hat_Noise_L) ){
                Beta_Hat_New_L[l] = 0.0;
              }else
              {
                Beta_Hat_New_L[l] = Beta_Hat_L;
              }
            }
            Beta_Hat_New_H[h] = Beta_Hat_New_L;
          }
          Beta_Hat_New_Dk[d] = Beta_Hat_New_H;
        }
        Beta_Hat_New[k] = Beta_Hat_New_Dk;
      }
      Beta_Hat = Beta_Hat_New;
    }

    Beta_Hat_List_Iter[m] = Beta_Hat;
    Beta_Hat_Old = Beta_Hat;
    
    List List_Mat_Sum_Beta_Hat(H);
    for(int h = 0; h < H; ++h){
      NumericMatrix Mat_Sum_Beta_Hat(K,L);
      for(int k = 0; k < K; ++k){
        List Beta_Hat_Temp_K = Beta_Hat[k]; //Temp_Trace_Bxt_gm[k];
        for(int l = 0; l < L; ++l){
          for(int d = 0; d < Dk[k]; ++d){
            double Mult_Factor = 1.0; //(Dk[k]*H);
            List Beta_Hat_Temp_Dk = Beta_Hat_Temp_K[d];
              NumericVector Beta_Hat_Temp_H = Beta_Hat_Temp_Dk[h];
              Mat_Sum_Beta_Hat(k,l) = Mat_Sum_Beta_Hat(k,l) + ((Beta_Hat_Temp_H[l]*Beta_Hat_Temp_H[l])/Mult_Factor);
          }
        }
      }
      List_Mat_Sum_Beta_Hat[h] = Mat_Sum_Beta_Hat;
    }
  
  LogicalVector Sum_Beta_Zero(H);
  for(int h = 0; h < H; ++h){
    NumericMatrix Mat_Sum_Beta_Hat = List_Mat_Sum_Beta_Hat[h];
    NumericVector Vec_Sum_Beta_Hat_K(K);
    for(int k = 0; k < K; ++k){
      NumericVector Vec_Sum_Beta_Hat_L(L);
      for(int l = 0; l < L; ++l){
        Vec_Sum_Beta_Hat_L[l] = Mat_Sum_Beta_Hat(k,l);
      }
      Vec_Sum_Beta_Hat_K[k] = Sum_C_NA(Vec_Sum_Beta_Hat_L);
    }
    Sum_Beta_Zero[h] = is_true(all(Vec_Sum_Beta_Hat_K == 0.0));
  }
  
  if(is_true(all(Sum_Beta_Zero)) ){
    M = (m-1);
    break;
  }
  
  Sum_Beta_Hat_List[m] = List_Mat_Sum_Beta_Hat;
  IntegerVector km_H(H);
  IntegerVector lm_H(H);
  for(int h = 0; h < H; ++h){
    NumericMatrix Mat_Sum_Beta_Hat = List_Mat_Sum_Beta_Hat[h];
    IntegerVector km_lm = Which_Max_Matrix_NA(Mat_Sum_Beta_Hat);
    km_H[h] = km_lm[0];
    lm_H[h] = km_lm[1];
    Variable_Select(m,h) = km_lm[0];
    Response_Select(m,h) = km_lm[1];
  }
    
    List Sum_Beta(K);
    for(int k = 0; k < K; ++k){
      List Beta_Hat_K = Beta_Hat[k];
      List Beta_K = Beta[k];
      List Sum_Beta_Dk(Dk[k]);
      for(int d = 0; d < Dk[k]; ++d){
        List Beta_Hat_Dk = Beta_Hat_K[d];
        List Beta_Dk = Beta_K[d];
        List Sum_Beta_H(H);
        for(int h = 0; h < H; ++h){
          NumericVector Beta_Hat_H = Beta_Hat_Dk[h];
          NumericVector Beta_H = Beta_Dk[h];
          NumericVector Sum_Beta_L(L);
          IntegerVector km_lm_Temp(2);
          km_lm_Temp[0] = km_H[h];
          km_lm_Temp[1] = lm_H[h];
          for(int l = 0; l < L; ++l){
            if(is_true(any( km_lm_Temp == -1 ))){
              Sum_Beta_L[l] = Beta_H[l];
            }else
            {
              if(k == km_H[h] && l == lm_H[h]){
                Sum_Beta_L[l] = Beta_H[l] + (nu*Beta_Hat_H[l]);  
              }else
              {
                Sum_Beta_L[l] = Beta_H[l];
              } 
            }
          }
          Sum_Beta_H[h] = Sum_Beta_L;
        }
        Sum_Beta_Dk[d] = Sum_Beta_H;
      }
      Sum_Beta[k] = Sum_Beta_Dk;
    }
    
    Beta = Sum_Beta;
    Beta_Hat_List[m] = Beta;
    
    mu = mu_zero;
    for(int k = 0; k < K; ++k){
      List Bxt_Temp_K = Bxt[k];
      List Beta_Temp_K = Beta[k];
      for(int d = 0; d < Dk[k]; ++d){
        List Bxt_Temp_Dk = Bxt_Temp_K[d];
        List Beta_Temp_Dk = Beta_Temp_K[d];
        for(int h = 0; h < H; ++h){
          List Bxt_Temp_H = Bxt_Temp_Dk[h];
          NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
          NumericMatrix mu_Hat(N,L);
          for(int l = 0; l < L; ++l){
            for(int i = 0; i < n; ++i){
              IntegerVector id_index_Temp = id_index[i];
              NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
              for(int j = 0; j < ni[i]; ++j){
                mu_Hat( id_index_Temp[j], l) = Bxt_Temp_n[j]*Beta_Temp_H[l];
              }
            }
          }
          mu = Matrix_Sum_C_NA(mu, mu_Hat);
        }
      }
    }
    
    NumericMatrix Org_mu(N,L);
    for(int l = 0; l < L; ++l){
      for(int i = 0; i < N; ++i){
        Org_mu(i,l) = ( mu(i,l) * y_Std_Error[l] ) + y_Mean[l];
      }
    }
    
    mu_List[m] = Org_mu;
    
    
    if(VarFlag){
      for(int l = 0; l < L; ++l){
        NumericVector Resid(N);
        for(int i = 0; i < N; ++i){
         Resid[i] = Org_y(i,l) - Org_mu(i,l); 
        }
        double phi_Temp;
        double phi_Init = 0.0;
        int phi_count = 0;
        for(int i = 0; i < n; ++i){
          IntegerVector id_index_Temp = id_index[i];
          for(int j = 0; j < ni[i]; ++j){
           phi_Init = phi_Init + (Resid[id_index_Temp[j]]*Resid[id_index_Temp[j]]);
           phi_count = phi_count + 1;  
          }
        }
        if(phi_count <= K){
          phi_Temp = phi_Init/phi_count;
        } else
        {
          phi_Temp = phi_Init/(phi_count - K);  
        }
        
        if(phi_Temp <= 0.0){
          phi[l] = 1.0;
        } else 
          {
            phi[l] = phi_Temp;
          }
        Phi(m,l) = phi[l];
        
        double rho_Temp;
        double rho_Init = 0.0;
        int rho_count = 0;
        for(int i = 0; i < n; ++i){
          IntegerVector id_index_Temp = id_index[i];
          if(ni[i] > 1){
            for(int j = 0; j < (ni[i]-1); ++j){
              for(int jj = (j+1); jj < ni[i]; ++jj){
                rho_Init = rho_Init + (Resid[id_index_Temp[j]]*Resid[id_index_Temp[jj]]);
                rho_count = rho_count + 1;  
              }
            }            
          }
        }   
        if(rho_count <= K){
          rho_Temp = rho_Init/(rho_count*phi[l]);
        } else
        {
          rho_Temp = rho_Init/((rho_count - K)*phi[l]);  
        }
        if(rho_Temp < -1 || rho_Temp > 1){
          rho[l] = 0.0;
        } else
        {
          rho[l] = rho_Temp;
        }
        Rho(m,l) = rho[l];
      }
    }
    
    
    for(int l = 0; l < L; ++l){
      NumericVector Org_y_Temp(N);
      NumericVector Org_mu_Temp(N);
      for(int i = 0; i < N; ++i){
        Org_y_Temp[i] = Org_y(i,l);
        Org_mu_Temp[i] = Org_mu(i,l);
      }
      Error_Rate(m,l) = l2Dist_Vector_C_NA(Org_y_Temp, Org_mu_Temp, id_index)/l2Dist_Vector_C_NA(Org_y_Temp, Vec_zero, id_index); 
    }
    
  } 
  
  List Tm_Beta_C(K);
  for(int k = 0; k < K; ++k){
    if(UseRaw[k]){
      List Beta_K = Beta[k];
      List Tm_Beta_Dk(Dk[k]);
      for(int d = 0; d < Dk[k]; ++d){
        List Beta_Dk = Beta_K[d];
        List Tm_Beta_H(H);
        for(int h = 0; h < H; ++h){
          NumericVector Beta_H = Beta_Dk[h];
          List Bt_n = Bt_H[h];
          List Tm_Beta_L(L);
          for(int l = 0; l < L; ++l){
            List Tm_Beta_n(n);
            for(int i = 0; i < n; ++i){
              NumericVector Bt_i = Bt_n[i];
              NumericVector Tm_Beta_ni(ni[i]);
              for(int j = 0; j < ni[i]; ++j){
                Tm_Beta_ni[j] = (Bt_i[j]*Beta_H[l]*y_Std_Error[l])/x_Std_Error[k];
              }
              Tm_Beta_n[i] = Tm_Beta_ni;
            }
            Tm_Beta_L[l] = Tm_Beta_n;
          }
          Tm_Beta_H[h] = Tm_Beta_L;
        }
        Tm_Beta_Dk[d] = Tm_Beta_H;
      }
      Tm_Beta_C[k] = Tm_Beta_Dk;
    }else
    {
      Tm_Beta_C[k] = NA_REAL;
    }
  }
  
  List Data = List::create(
    _["Org_x"] = Org_x, 
    _["Org_y"] = Org_y,
    _["id"] = id,
    _["tm"] = tm,
    _["x"] = x, 
    _["y"] = y,
    _["x_Mean"] = x_Mean,
    _["x_Std_Error"] = x_Std_Error,
    _["y_Mean"] = y_Mean,
    _["y_Std_Error"] = y_Std_Error
  );
  List Dimensions = List::create(
    _["n"] = n,
    _["K"] = K,
    _["L"] = L,
    _["H"] = H,
    _["Dk"] = Dk,
    _["ni"] = ni,
    _["N"] = N
  );
  
  List Index = List::create(
    _["unq_id"] = unq_id,
    _["unq_tm"] = unq_tm,
    _["unq_x"] = unq_x,
    _["id_index"] = id_index,
    _["tm_index"] = tm_index,
    _["x_index"] = x_index
  );
  
  List BS = List::create(
    _["Bt"] = Bt,
    _["Bx"] = Bx,
    _["Bt_H"] = Bt_H,
    _["Bx_K"] = Bx_K,
    _["Bxt"] = Bxt,
    _["Bx_Scale"] = Bx_Scale
  );
  
  List Regulate = List::create(
    _["M"] = M,
    _["nu"] = nu,
    _["Lambda_Ridge_Vec"] = Lambda_Ridge_Vec,
    _["Shrink"] = Shrink,
    _["Ridge_Penalty"] = Ridge_Penalty,
    _["Lambda_Scale"] = Lambda_Scale,
    _["NLambda"] = NLambda,
    _["lower_perc"] = lower_perc,
    _["upper_perc"] = upper_perc
  );
  
  List Beta_Estimate = List::create(
    _["Beta"] = Beta,
    _["Beta_Hat_List"] = Beta_Hat_List,
    _["Beta_Hat_List_Iter"] = Beta_Hat_List_Iter,
    _["Sum_Beta_Hat_List"] = Sum_Beta_Hat_List,
    _["Tm_Beta_C"] = Tm_Beta_C,
    _["lower_Beta_Hat_Noise"] = lower_Beta_Hat_Noise,
    _["upper_Beta_Hat_Noise"] = upper_Beta_Hat_Noise,
    _["List_Trace_Bxt_gm"] = List_Trace_Bxt_gm
  );
  
  return List::create(
    _["Data"] = Data,
    _["Dimensions"] = Dimensions,
    _["Index"] = Index,
    _["BS"] = BS,
    _["Regulate"] = Regulate,
    _["Beta_Estimate"] = Beta_Estimate,
    _["Error_Rate"] = Error_Rate,
    _["Variable_Select"] = Variable_Select,
    _["Response_Select"] = Response_Select,
    _["mu"] = mu,
    _["mu_List"] = mu_List,
    _["Phi"] = Phi,
    _["Rho"] = Rho,
    _["Lambda_List"] = Lambda_List,
    _["mu_zero"] = mu_zero,
    _["Vec_zero"] = Vec_zero
  );
}




// [[Rcpp::export]]
List update_BoostMLR_C(NumericMatrix Org_x,
                       NumericMatrix Org_y,
                       NumericVector id,
                       NumericVector tm,
                       NumericMatrix x,
                       NumericMatrix y,
                       NumericVector x_Mean,
                       NumericVector x_Std_Error,
                       NumericVector y_Mean,
                       NumericVector y_Std_Error,
                       int n,
                       int K,
                       int L,
                       int H,
                       IntegerVector Dk,
                       IntegerVector ni,
                       int N,
                       NumericVector unq_id,
                       NumericVector unq_tm,
                       List unq_x,
                       List id_index,
                       List tm_index,
                       List x_index,
                       NumericMatrix Bt,
                       List Bx,
                       List Bt_H,
                       List Bx_K,
                       List Bxt,
                       List Bx_Scale,
                       double nu,
                       int M,
                       int M_New,
                       LogicalVector UseRaw,
                       bool Shrink,
                       bool Ridge_Penalty,
                       NumericVector Lambda_Ridge_Vec,
                       double Lambda_Scale,
                       int NLambda,
                       double lower_perc,
                       double upper_perc,
                       List Lambda_List,
                       NumericMatrix mu,
                       List mu_List,
                       NumericMatrix mu_zero,
                       NumericVector Vec_zero,
                       NumericMatrix Error_Rate,
                       IntegerMatrix Variable_Select,
                       IntegerMatrix Response_Select,
                       List Beta_Hat_List,
                       List Sum_Beta_Hat_List,
                       List Beta,
                       List Beta_Hat_List_Iter,
                       List lower_Beta_Hat_Noise,
                       List upper_Beta_Hat_Noise,
                       List List_Trace_Bxt_gm,
                       bool Mod_Grad,
                       bool VarFlag,
                       NumericVector phi,
                       NumericVector rho,
                       NumericMatrix Phi,
                       NumericMatrix Rho,
                       bool Verbose)
{
  
  NumericMatrix UP_Error_Rate(M_New,L);
  NumericMatrix UP_Phi(M_New,L);
  NumericMatrix UP_Rho(M_New,L);
  IntegerMatrix UP_Variable_Select(M_New,H);
  IntegerMatrix UP_Response_Select(M_New,H);
  List UP_Beta_Hat_List(M_New);
  List UP_Sum_Beta_Hat_List(M_New);
  List UP_Beta_Hat_List_Iter(M_New);
  List UP_mu_List(M_New);
  List UP_Lambda_List(M_New);
  List UP_List_Trace_Bxt_gm(M_New);
  
  for(int m = 0; m < M_New; ++m){
    if(m < M){
      for(int h = 0; h < H; ++h){
        int Int_Variable_Select;
        Int_Variable_Select = Variable_Select(m,h);
        int Int_Response_Select;
        Int_Response_Select = Response_Select(m,h);
        UP_Variable_Select(m,h) = Int_Variable_Select;
        UP_Response_Select(m,h) = Int_Response_Select;
      }
    }else
    {
      for(int h = 0; h < H; ++h){
        UP_Variable_Select(m,h) = NA_INTEGER;
        UP_Response_Select(m,h) = NA_INTEGER;
      }
    }
  }
  
  for(int m = 0; m < M_New; ++m){
    
    if(m < M){
      UP_Beta_Hat_List[m] = Beta_Hat_List[m];
      UP_Sum_Beta_Hat_List[m] = Sum_Beta_Hat_List[m];
      UP_Beta_Hat_List_Iter[m] = Beta_Hat_List_Iter[m];
      UP_mu_List[m] = mu_List[m];
      UP_Lambda_List[m] = Lambda_List[m];
      UP_List_Trace_Bxt_gm[m] = List_Trace_Bxt_gm[m];
      for(int l = 0; l < L; ++l){
        UP_Error_Rate(m,l) = Error_Rate(m,l);
        UP_Phi(m,l) = Phi(m,l);
        UP_Rho(m,l) = Rho(m,l);
      }
    }else
    {
      UP_Beta_Hat_List[m] = NA_REAL;
      UP_Sum_Beta_Hat_List[m] = NA_REAL;
      UP_Beta_Hat_List_Iter[m] = NA_REAL;
      UP_mu_List[m] = NA_REAL;
      UP_Lambda_List[m] = NA_REAL;
      UP_List_Trace_Bxt_gm[m] = List_Trace_Bxt_gm;
      for(int l = 0; l < L; ++l){
        UP_Error_Rate(m,l) = NA_REAL;
        UP_Phi(m,l) = NA_REAL;
        UP_Rho(m,l) = NA_REAL;
      }
    }
  }
  
  Variable_Select = UP_Variable_Select;
  Response_Select = UP_Response_Select;
  Beta_Hat_List = UP_Beta_Hat_List;
  Sum_Beta_Hat_List = UP_Sum_Beta_Hat_List;
  Beta_Hat_List_Iter = UP_Beta_Hat_List_Iter;
  mu_List = UP_mu_List;
  Lambda_List = UP_Lambda_List;
  List_Trace_Bxt_gm = UP_List_Trace_Bxt_gm;
  Error_Rate = UP_Error_Rate;
  Phi = UP_Phi;
  Rho = UP_Rho;
  List Beta_Hat_Old(K);
  Beta_Hat_Old = Beta;
  List V_inv(L);
  
  for(int m = M; m <  M_New; ++m){
    
    if(Verbose){
      if( (m + 1 - M)%((M_New - M)/10) == 0.0){
        double FracBoosting;
        FracBoosting = ((m + 1 - M)*100)/(M_New - M);
        int Perc_boosting = FracBoosting;
         Rcout << Perc_boosting << "%" << " ";
      }
    }
    
    if(m == M || VarFlag == TRUE){
      for(int l = 0; l < L; ++l){
        List Vi_inv(n);
        for(int i = 0; i < n;++i){
          arma::mat Vi_inv_Temp;
          Vi_inv_Temp = MatrixInversion_Equicorrelation_C(ni[i],phi[l],rho[l]);
          Vi_inv[i] = Vi_inv_Temp;
        }
        V_inv[l] = Vi_inv;
      }
    }
    
    List gm(L);
    for(int l = 0; l < L; ++l){
      List gm_n(n);
      List V_inv_l = V_inv[l];
      for(int i = 0; i < n; ++i){
        arma::mat V_inv_i = V_inv_l[i];
        IntegerVector id_index_Temp = id_index[i];
        NumericVector gm_ni(ni[i]);
        for(int j = 0; j < ni[i]; ++j){
          gm_ni[j] = y(id_index_Temp[j],l) - mu(id_index_Temp[j],l);
        }
        if(Mod_Grad){
          gm_n[i] = gm_ni;  
        }
        else {
          gm_n[i] = Matrix_Vector_Multiplication_C(gm_ni,V_inv_i);
        }
      }
      gm[l] = gm_n;
    }
    
    List Temp_Trace_Bxt_gm(K);
    List Lambda(K);
    List Beta_Hat(K);
    if(m > 0){
      for(int k = 0; k < K; ++k){
        List Beta_Hat_Old_Temp_Dk = Beta_Hat_Old[k];
        List Bxt_Temp_K = Bxt[k];
        NumericVector Bx_Scale_K = Bx_Scale[k];
        List Lambda_Dk(Dk[k]);
        List Beta_Hat_Dk(Dk[k]);
        List Temp_Trace_Bxt_gm_Dk(Dk[k]);
        for(int d = 0; d < Dk[k]; ++d){
          List Beta_Hat_Old_Temp_H = Beta_Hat_Old_Temp_Dk[d];
          List Bxt_Temp_Dk = Bxt_Temp_K[d];
          List Lambda_H(H);
          List Beta_Hat_H(H);
          List Temp_Trace_Bxt_gm_H(H);
          for(int h = 0; h < H; ++h){
            NumericVector Beta_Hat_Old_Temp_L = Beta_Hat_Old_Temp_H[h];
            List Bxt_Temp_H = Bxt_Temp_Dk[h];
            NumericVector Lambda_L(L);
            NumericVector Beta_Hat_L(L);
            NumericVector Temp_Trace_Bxt_gm_L(L);
            for(int l = 0; l < L; ++l){
              double Beta_Hat_Old_Temp = Beta_Hat_Old_Temp_L[l];
              if(Beta_Hat_Old_Temp == 0.0){
                Beta_Hat_L[l] = 0.0;
              } else
              {
                List V_inv_l = V_inv[l];
                List gm_Temp_L = gm[l];
                NumericVector Bxt_n_2(n);
                NumericVector Bxt_gm_n(n);
                NumericVector gm_n_2(n);
                for(int i = 0; i < n; ++i){
                  arma::mat V_inv_i = V_inv_l[i];
                  NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                  NumericVector gm_Temp_n = gm_Temp_L[i];
                  NumericVector Bxt_V = Matrix_Vector_Multiplication_C(Bxt_Temp_n,V_inv_i);
                  NumericVector Bxt_ni_2(ni[i]);
                  NumericVector Bxt_gm_ni(ni[i]);
                  NumericVector gm_ni_2(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    Bxt_ni_2[j] = Bxt_V[j] * Bxt_Temp_n[j];
                    Bxt_gm_ni[j] = Bxt_V[j] * gm_Temp_n[j];
                    gm_ni_2[j] = gm_Temp_n[j] * gm_Temp_n[j];
                  }
                  Bxt_n_2[i] = Sum_C_NA(Bxt_ni_2);
                  Bxt_gm_n[i] = Sum_C_NA(Bxt_gm_ni);
                  gm_n_2[i] = Sum_C_NA(gm_ni_2);
                }
                double Trace_Bxt    = Sum_C_NA(Bxt_n_2);
                double Trace_Bxt_gm = Sum_C_NA(Bxt_gm_n);
                Temp_Trace_Bxt_gm_L[l] = Trace_Bxt_gm;
                if(Trace_Bxt < 0.001 ){
                  Beta_Hat_L[l] = 0.0;
                }else
                {
                  if(Ridge_Penalty){
                    /* 
                     Come back to this later since I don't know how to calculate
                     (Vi)^{-1/2} require for calculating lambda for Ridge penalty.
                     Therefore, default setting Ridge_Penalty = FALSE should not be changed.
                     Trace_gm specified below will be incorrect because it does not include Vi_inv matrix.
                     */
                    double Trace_gm       = Sum_C_NA(gm_n_2);
                    double Beta_Hat_NL    = Beta_Hat_Old_Temp;
                    double Trace_eps      = Trace_gm + (Beta_Hat_NL * Beta_Hat_NL * Trace_Bxt) - (2 * Beta_Hat_NL * Trace_Bxt_gm);
                    double Lambda_L_Temp  = (Trace_Bxt/(Trace_gm - Trace_eps))*Lambda_Scale;
                    if(Lambda_L_Temp > 5000){
                      Lambda_L[l] = 5000;
                    }else
                    {
                      if(Lambda_L_Temp < 0){
                        Lambda_L[l] = Lambda_Ridge_Vec[k];
                      }else
                      {
                        Lambda_L[l] = Lambda_L_Temp;  
                      }
                    }
                  }else
                  {
                    Lambda_L[l]         = Lambda_Ridge_Vec[k];   
                  }
                  Beta_Hat_L[l]         = (Trace_Bxt_gm/(Trace_Bxt + Lambda_L[l]))/Bx_Scale_K[d];
                }
              }
            }
            Lambda_H[h] = Lambda_L;
            Beta_Hat_H[h] = Beta_Hat_L; 
            Temp_Trace_Bxt_gm_H[h] = Temp_Trace_Bxt_gm_L;
          }
          Lambda_Dk[d] = Lambda_H;
          Beta_Hat_Dk[d] = Beta_Hat_H;
          Temp_Trace_Bxt_gm_Dk[d] = Temp_Trace_Bxt_gm_H;
        }
        Lambda[k] = Lambda_Dk;
        Beta_Hat[k] = Beta_Hat_Dk;
        Temp_Trace_Bxt_gm[k] = Temp_Trace_Bxt_gm_Dk;
      }
    }else
    {
      for(int k = 0; k < K; ++k){
        List Bxt_Temp_K = Bxt[k];
        NumericVector Bx_Scale_K = Bx_Scale[k];
        List Lambda_Dk(Dk[k]);
        List Beta_Hat_Dk(Dk[k]);
        List Temp_Trace_Bxt_gm_Dk(Dk[k]);
        for(int d = 0; d < Dk[k]; ++d){
          List Bxt_Temp_Dk = Bxt_Temp_K[d];
          List Lambda_H(H);
          List Beta_Hat_H(H);
          List Temp_Trace_Bxt_gm_H(H);
          for(int h = 0; h < H; ++h){
            List Bxt_Temp_H = Bxt_Temp_Dk[h];
            NumericVector Lambda_L(L);
            NumericVector Beta_Hat_L(L);
            NumericVector Temp_Trace_Bxt_gm_L(L);
            for(int l = 0; l < L; ++l){
              List V_inv_l = V_inv[l];
              List gm_Temp_L = gm[l];
              NumericVector Bxt_n_2(n);
              NumericVector Bxt_gm_n(n);
              NumericVector gm_n_2(n);
              for(int i = 0; i < n; ++i){
                arma::mat V_inv_i = V_inv_l[i];
                NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                NumericVector gm_Temp_n = gm_Temp_L[i];
                NumericVector Bxt_V = Matrix_Vector_Multiplication_C(Bxt_Temp_n,V_inv_i);
                NumericVector Bxt_ni_2(ni[i]);
                NumericVector Bxt_gm_ni(ni[i]);
                NumericVector gm_ni_2(ni[i]);
                for(int j = 0; j < ni[i]; ++j){
                  Bxt_ni_2[j] = Bxt_V[j] * Bxt_Temp_n[j];
                  Bxt_gm_ni[j] = Bxt_V[j] * gm_Temp_n[j];
                  gm_ni_2[j] = gm_Temp_n[j] * gm_Temp_n[j];
                }
                Bxt_n_2[i] = Sum_C_NA(Bxt_ni_2);
                Bxt_gm_n[i] = Sum_C_NA(Bxt_gm_ni);
                gm_n_2[i] = Sum_C_NA(gm_ni_2);
              }
              double Trace_Bxt    = Sum_C_NA(Bxt_n_2);
              double Trace_Bxt_gm = Sum_C_NA(Bxt_gm_n);
              Temp_Trace_Bxt_gm_L[l] = Trace_Bxt_gm;
              if(Trace_Bxt < 0.001 ){
                Beta_Hat_L[l] = 0.0;
              }else
              {
                if(Ridge_Penalty){
                  /* 
                   Come back to this later since I don't know how to calculate
                   (Vi)^{-1/2} require for calculating lambda for Ridge penalty.
                   Therefore, default setting Ridge_Penalty = FALSE should not be changed.
                   Trace_gm specified below will be incorrect because it does not include Vi_inv matrix.
                   */
                  double Trace_gm       = Sum_C_NA(gm_n_2);
                  double Beta_Hat_NL    = Trace_Bxt_gm/(Trace_Bxt);
                  double Trace_eps      = Trace_gm + (Beta_Hat_NL * Beta_Hat_NL * Trace_Bxt) - (2 * Beta_Hat_NL * Trace_Bxt_gm);
                  double Lambda_L_Temp  = (Trace_Bxt/(Trace_gm - Trace_eps))*Lambda_Scale;
                  if(Lambda_L_Temp > 5000){
                    Lambda_L[l] = 5000;
                  }else
                  {
                    if(Lambda_L_Temp < 0){
                      Lambda_L[l] = Lambda_Ridge_Vec[k];
                    }else
                    {
                      Lambda_L[l] = Lambda_L_Temp;  
                    }
                  }
                }else
                {
                  Lambda_L[l]         = Lambda_Ridge_Vec[k];
                }
                Beta_Hat_L[l]         = (Trace_Bxt_gm/(Trace_Bxt + Lambda_L[l]))/Bx_Scale_K[d]; 
              }
            }
            Lambda_H[h] = Lambda_L;
            Beta_Hat_H[h] = Beta_Hat_L; 
            Temp_Trace_Bxt_gm_H[h] = Temp_Trace_Bxt_gm_L;
          }
          Lambda_Dk[d] = Lambda_H;
          Beta_Hat_Dk[d] = Beta_Hat_H;
          Temp_Trace_Bxt_gm_Dk[d] = Temp_Trace_Bxt_gm_H;
        }
        Lambda[k] = Lambda_Dk;
        Beta_Hat[k] = Beta_Hat_Dk;
        Temp_Trace_Bxt_gm[k] = Temp_Trace_Bxt_gm_Dk;
      } 
    }
    
    List_Trace_Bxt_gm[m] = Temp_Trace_Bxt_gm;
    Lambda_List[m] = Lambda;
    
    if(Shrink){
      if(m == 0){
        for(int k = 0; k < K; ++k){
          List Bxt_Temp_K = Bxt[k];
          NumericVector Bx_Scale_K = Bx_Scale[k];
          List lower_Beta_Hat_Dk(Dk[k]);
          List upper_Beta_Hat_Dk(Dk[k]);
          for(int d = 0; d < Dk[k]; ++d){
            List Bxt_Temp_Dk = Bxt_Temp_K[d];
            List lower_Beta_Hat_H(H);
            List upper_Beta_Hat_H(H);
            for(int h = 0; h < H; ++h){
              List Bxt_Temp_H = Bxt_Temp_Dk[h];
              NumericVector lower_Beta_Hat_L(L);
              NumericVector upper_Beta_Hat_L(L);
              for(int l = 0; l < L; ++l){
                List V_inv_l = V_inv[l];
                NumericVector Beta_Hat_Noise_Temp(NLambda);
                for(int l_m = 0; l_m < NLambda; ++l_m){
                  List gm_Temp_L = gm[l];
                  NumericVector gm_unlist(N);
                  int count = 0;
                  for(int i = 0; i < n; ++i){
                    NumericVector gm_Temp_n = gm_Temp_L[i];
                    for(int j = 0; j < ni[i]; ++j){
                      gm_unlist[count] = gm_Temp_n[j];
                      count = count + 1;
                    }
                  }
                  NumericVector gm_unlist_noise = randomShuffle(gm_unlist);
                  count = 0;
                  List gm_noise_n(n);
                  for(int i = 0; i < n; ++i){
                    NumericVector gm_noise_ni(ni[i]);
                    for(int j = 0; j < ni[i]; ++j){
                      gm_noise_ni[j] = gm_unlist_noise[count];
                      count = count + 1;
                    }
                    gm_noise_n[i] = gm_noise_ni;
                  }
                  NumericVector Bxt_n_2(n);
                  NumericVector Bxt_gm_n(n);
                  NumericVector gm_noise_n_2(n);
                  for(int i = 0; i < n; ++i){
                    arma::mat V_inv_i = V_inv_l[i];
                    NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                    NumericVector gm_noise_Temp_n = gm_noise_n[i];
                    NumericVector Bxt_V = Matrix_Vector_Multiplication_C(Bxt_Temp_n,V_inv_i);
                    NumericVector Bxt_ni_2(ni[i]);
                    NumericVector Bxt_gm_ni(ni[i]);
                    NumericVector gm_noise_ni_2(ni[i]);
                    for(int j = 0; j < ni[i]; ++j){
                      Bxt_ni_2[j] = Bxt_V[j] * Bxt_Temp_n[j];
                      Bxt_gm_ni[j] = Bxt_V[j] * gm_noise_Temp_n[j];
                      gm_noise_ni_2[j] = gm_noise_Temp_n[j] * gm_noise_Temp_n[j];
                    }
                    Bxt_n_2[i] = Sum_C_NA(Bxt_ni_2);
                    Bxt_gm_n[i] = Sum_C_NA(Bxt_gm_ni);
                    gm_noise_n_2[i] = Sum_C_NA(gm_noise_ni_2);
                  }
                  double Trace_Bxt    = Sum_C_NA(Bxt_n_2);
                  double Trace_Bxt_gm = Sum_C_NA(Bxt_gm_n);
                  if(Trace_Bxt == 0.0){
                    Beta_Hat_Noise_Temp[l_m]  = 0.0;
                  }else
                  {
                    if(Ridge_Penalty){
                      /* 
                       Come back to this later since I don't know how to calculate
                       (Vi)^{-1/2} require for calculating lambda for Ridge penalty.
                       Therefore, default setting Ridge_Penalty = FALSE should not be changed.
                       Trace_gm specified below will be incorrect because it does not include Vi_inv matrix.
                       */
                      double Trace_gm       = Sum_C_NA(gm_noise_n_2);
                      double Beta_Hat_NL    = Trace_Bxt_gm/(Trace_Bxt);
                      double Trace_eps      = Trace_gm + (Beta_Hat_NL * Beta_Hat_NL * Trace_Bxt) - (2 * Beta_Hat_NL * Trace_Bxt_gm);
                      double Lambda_L_Temp  = (Trace_Bxt/(Trace_gm - Trace_eps))*Lambda_Scale;
                      if(Lambda_L_Temp > 5000){
                        Lambda_L_Temp = 5000;
                      }else
                      {
                        if(Lambda_L_Temp < 0){
                          Lambda_L_Temp = Lambda_Ridge_Vec[k];
                        } //else
                          // {
                          // Lambda_L_Temp = Lambda_L_Temp;  
                          // }
                      }
                      Beta_Hat_Noise_Temp[l_m]         = (Trace_Bxt_gm/(Trace_Bxt + Lambda_L_Temp))/Bx_Scale_K[d];
                    }else
                    {
                      Beta_Hat_Noise_Temp[l_m]         = (Trace_Bxt_gm/(Trace_Bxt + Lambda_Ridge_Vec[k]))/Bx_Scale_K[d]; 
                    }
                  }
                }
                NumericVector sort_Beta_Hat_Noise = stl_sort_NA(Beta_Hat_Noise_Temp);
                int lower_range = (NLambda*lower_perc);
                int upper_range = (NLambda*upper_perc);
                lower_Beta_Hat_L[l] = sort_Beta_Hat_Noise[lower_range];
                upper_Beta_Hat_L[l] = sort_Beta_Hat_Noise[upper_range];
              }
              lower_Beta_Hat_H[h] = lower_Beta_Hat_L;
              upper_Beta_Hat_H[h] = upper_Beta_Hat_L;
            }
            lower_Beta_Hat_Dk[d] = lower_Beta_Hat_H;
            upper_Beta_Hat_Dk[d] = upper_Beta_Hat_H;
          }
          lower_Beta_Hat_Noise[k] = lower_Beta_Hat_Dk;
          upper_Beta_Hat_Noise[k] = upper_Beta_Hat_Dk;
        }
      }
      
      List Beta_Hat_New(K);
      for(int k = 0; k < K; ++k){
        List Beta_Hat_K = Beta_Hat[k];
        List lower_Beta_Hat_Noise_K = lower_Beta_Hat_Noise[k];
        List upper_Beta_Hat_Noise_K = upper_Beta_Hat_Noise[k];
        List Beta_Hat_New_Dk(Dk[k]);
        for(int d = 0; d < Dk[k]; ++d){
          List Beta_Hat_Dk = Beta_Hat_K[d];
          List lower_Beta_Hat_Noise_Dk = lower_Beta_Hat_Noise_K[d];
          List upper_Beta_Hat_Noise_Dk = upper_Beta_Hat_Noise_K[d];
          List Beta_Hat_New_H(H);
          for(int h = 0; h < H; ++h){
            NumericVector Beta_Hat_H = Beta_Hat_Dk[h];
            NumericVector lower_Beta_Hat_Noise_H = lower_Beta_Hat_Noise_Dk[h];
            NumericVector upper_Beta_Hat_Noise_H = upper_Beta_Hat_Noise_Dk[h];
            NumericVector Beta_Hat_New_L(L);
            for(int l = 0; l < L; ++l){
              double Beta_Hat_L = Beta_Hat_H[l];
              double lower_Beta_Hat_Noise_L = lower_Beta_Hat_Noise_H[l];
              double upper_Beta_Hat_Noise_L = upper_Beta_Hat_Noise_H[l];
              if( (Beta_Hat_L >=  lower_Beta_Hat_Noise_L) && (Beta_Hat_L <= upper_Beta_Hat_Noise_L) ){
                Beta_Hat_New_L[l] = 0.0;
              }else
              {
                Beta_Hat_New_L[l] = Beta_Hat_L;
              }
            }
            Beta_Hat_New_H[h] = Beta_Hat_New_L;
          }
          Beta_Hat_New_Dk[d] = Beta_Hat_New_H;
        }
        Beta_Hat_New[k] = Beta_Hat_New_Dk;
      }
      Beta_Hat = Beta_Hat_New;
    }
    
    Beta_Hat_List_Iter[m] = Beta_Hat;
    Beta_Hat_Old = Beta_Hat;
    
    List List_Mat_Sum_Beta_Hat(H);
    for(int h = 0; h < H; ++h){
      NumericMatrix Mat_Sum_Beta_Hat(K,L);
      for(int k = 0; k < K; ++k){
        List Beta_Hat_Temp_K = Beta_Hat[k]; //Temp_Trace_Bxt_gm[k];
        for(int l = 0; l < L; ++l){
          for(int d = 0; d < Dk[k]; ++d){
            double Mult_Factor = 1.0; //(Dk[k]*H);
            List Beta_Hat_Temp_Dk = Beta_Hat_Temp_K[d];
            NumericVector Beta_Hat_Temp_H = Beta_Hat_Temp_Dk[h];
            Mat_Sum_Beta_Hat(k,l) = Mat_Sum_Beta_Hat(k,l) + ((Beta_Hat_Temp_H[l]*Beta_Hat_Temp_H[l])/Mult_Factor);
          }
        }
      }
      List_Mat_Sum_Beta_Hat[h] = Mat_Sum_Beta_Hat;
    }
    
    LogicalVector Sum_Beta_Zero(H);
    for(int h = 0; h < H; ++h){
      NumericMatrix Mat_Sum_Beta_Hat = List_Mat_Sum_Beta_Hat[h];
      NumericVector Vec_Sum_Beta_Hat_K(K);
      for(int k = 0; k < K; ++k){
        NumericVector Vec_Sum_Beta_Hat_L(L);
        for(int l = 0; l < L; ++l){
          Vec_Sum_Beta_Hat_L[l] = Mat_Sum_Beta_Hat(k,l);
        }
        Vec_Sum_Beta_Hat_K[k] = Sum_C_NA(Vec_Sum_Beta_Hat_L);
      }
      Sum_Beta_Zero[h] = is_true(all(Vec_Sum_Beta_Hat_K == 0.0));
    }
    
    if(is_true(all(Sum_Beta_Zero)) ){
      M_New = (m-1);
      break;
    }
    
    Sum_Beta_Hat_List[m] = List_Mat_Sum_Beta_Hat;
    IntegerVector km_H(H);
    IntegerVector lm_H(H);
    for(int h = 0; h < H; ++h){
      NumericMatrix Mat_Sum_Beta_Hat = List_Mat_Sum_Beta_Hat[h];
      IntegerVector km_lm = Which_Max_Matrix_NA(Mat_Sum_Beta_Hat);
      km_H[h] = km_lm[0];
      lm_H[h] = km_lm[1];
      Variable_Select(m,h) = km_lm[0];
      Response_Select(m,h) = km_lm[1];
    }
    
    List Sum_Beta(K);
    for(int k = 0; k < K; ++k){
      List Beta_Hat_K = Beta_Hat[k];
      List Beta_K = Beta[k];
      List Sum_Beta_Dk(Dk[k]);
      for(int d = 0; d < Dk[k]; ++d){
        List Beta_Hat_Dk = Beta_Hat_K[d];
        List Beta_Dk = Beta_K[d];
        List Sum_Beta_H(H);
        for(int h = 0; h < H; ++h){
          NumericVector Beta_Hat_H = Beta_Hat_Dk[h];
          NumericVector Beta_H = Beta_Dk[h];
          NumericVector Sum_Beta_L(L);
          IntegerVector km_lm_Temp(2);
          km_lm_Temp[0] = km_H[h];
          km_lm_Temp[1] = lm_H[h];
          for(int l = 0; l < L; ++l){
            if(is_true(any( km_lm_Temp == -1 ))){
              Sum_Beta_L[l] = Beta_H[l];
            }else
            {
              if(k == km_H[h] && l == lm_H[h]){
                Sum_Beta_L[l] = Beta_H[l] + (nu*Beta_Hat_H[l]);  
              }else
              {
                Sum_Beta_L[l] = Beta_H[l];
              } 
            }
          }
          Sum_Beta_H[h] = Sum_Beta_L;
        }
        Sum_Beta_Dk[d] = Sum_Beta_H;
      }
      Sum_Beta[k] = Sum_Beta_Dk;
    }
    
    Beta = Sum_Beta;
    Beta_Hat_List[m] = Beta;
    
    mu = mu_zero;
    for(int k = 0; k < K; ++k){
      List Bxt_Temp_K = Bxt[k];
      List Beta_Temp_K = Beta[k];
      for(int d = 0; d < Dk[k]; ++d){
        List Bxt_Temp_Dk = Bxt_Temp_K[d];
        List Beta_Temp_Dk = Beta_Temp_K[d];
        for(int h = 0; h < H; ++h){
          List Bxt_Temp_H = Bxt_Temp_Dk[h];
          NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
          NumericMatrix mu_Hat(N,L);
          for(int l = 0; l < L; ++l){
            for(int i = 0; i < n; ++i){
              IntegerVector id_index_Temp = id_index[i];
              NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
              for(int j = 0; j < ni[i]; ++j){
                mu_Hat( id_index_Temp[j], l) = Bxt_Temp_n[j]*Beta_Temp_H[l];
              }
            }
          }
          mu = Matrix_Sum_C_NA(mu, mu_Hat);
        }
      }
    }
    
    NumericMatrix Org_mu(N,L);
    for(int l = 0; l < L; ++l){
      for(int i = 0; i < N; ++i){
        Org_mu(i,l) = ( mu(i,l) * y_Std_Error[l] ) + y_Mean[l];
      }
    }
    
    mu_List[m] = Org_mu;
    
    
    if(VarFlag){
      for(int l = 0; l < L; ++l){
        NumericVector Resid(N);
        for(int i = 0; i < N; ++i){
          Resid[i] = Org_y(i,l) - Org_mu(i,l); 
        }
        double phi_Temp;
        double phi_Init = 0.0;
        int phi_count = 0;
        for(int i = 0; i < n; ++i){
          IntegerVector id_index_Temp = id_index[i];
          for(int j = 0; j < ni[i]; ++j){
            phi_Init = phi_Init + (Resid[id_index_Temp[j]]*Resid[id_index_Temp[j]]);
            phi_count = phi_count + 1;  
          }
        }
        if(phi_count <= K){
          phi_Temp = phi_Init/phi_count;
        } else
        {
          phi_Temp = phi_Init/(phi_count - K);  
        }
        
        if(phi_Temp <= 0.0){
          phi[l] = 1.0;
        } else 
        {
          phi[l] = phi_Temp;
        }
        Phi(m,l) = phi[l];
        
        double rho_Temp;
        double rho_Init = 0.0;
        int rho_count = 0;
        for(int i = 0; i < n; ++i){
          IntegerVector id_index_Temp = id_index[i];
          if(ni[i] > 1){
            for(int j = 0; j < (ni[i]-1); ++j){
              for(int jj = (j+1); jj < ni[i]; ++jj){
                rho_Init = rho_Init + (Resid[id_index_Temp[j]]*Resid[id_index_Temp[jj]]);
                rho_count = rho_count + 1;  
              }
            }            
          }
        }   
        if(rho_count <= K){
          rho_Temp = rho_Init/(rho_count*phi[l]);
        } else
        {
          rho_Temp = rho_Init/((rho_count - K)*phi[l]);  
        }
        if(rho_Temp < -1 || rho_Temp > 1){
          rho[l] = 0.0;
        } else
        {
          rho[l] = rho_Temp;
        }
        Rho(m,l) = rho[l];
      }
    }
    
    
    for(int l = 0; l < L; ++l){
      NumericVector Org_y_Temp(N);
      NumericVector Org_mu_Temp(N);
      for(int i = 0; i < N; ++i){
        Org_y_Temp[i] = Org_y(i,l);
        Org_mu_Temp[i] = Org_mu(i,l);
      }
      Error_Rate(m,l) = l2Dist_Vector_C_NA(Org_y_Temp, Org_mu_Temp, id_index)/l2Dist_Vector_C_NA(Org_y_Temp, Vec_zero, id_index); 
    }
    
  } 
  
  List Tm_Beta_C(K);
  for(int k = 0; k < K; ++k){
    if(UseRaw[k]){
      List Beta_K = Beta[k];
      List Tm_Beta_Dk(Dk[k]);
      for(int d = 0; d < Dk[k]; ++d){
        List Beta_Dk = Beta_K[d];
        List Tm_Beta_H(H);
        for(int h = 0; h < H; ++h){
          NumericVector Beta_H = Beta_Dk[h];
          List Bt_n = Bt_H[h];
          List Tm_Beta_L(L);
          for(int l = 0; l < L; ++l){
            List Tm_Beta_n(n);
            for(int i = 0; i < n; ++i){
              NumericVector Bt_i = Bt_n[i];
              NumericVector Tm_Beta_ni(ni[i]);
              for(int j = 0; j < ni[i]; ++j){
                Tm_Beta_ni[j] = (Bt_i[j]*Beta_H[l]*y_Std_Error[l])/x_Std_Error[k];
              }
              Tm_Beta_n[i] = Tm_Beta_ni;
            }
            Tm_Beta_L[l] = Tm_Beta_n;
          }
          Tm_Beta_H[h] = Tm_Beta_L;
        }
        Tm_Beta_Dk[d] = Tm_Beta_H;
      }
      Tm_Beta_C[k] = Tm_Beta_Dk;
    }else
    {
      Tm_Beta_C[k] = NA_REAL;
    }
  }
  
  List Data = List::create(
    _["Org_x"] = Org_x, 
    _["Org_y"] = Org_y,
    _["id"] = id,
    _["tm"] = tm,
    _["x"] = x, 
    _["y"] = y,
    _["x_Mean"] = x_Mean,
    _["x_Std_Error"] = x_Std_Error,
    _["y_Mean"] = y_Mean,
    _["y_Std_Error"] = y_Std_Error
  );
  List Dimensions = List::create(
    _["n"] = n,
    _["K"] = K,
    _["L"] = L,
    _["H"] = H,
    _["Dk"] = Dk,
    _["ni"] = ni,
    _["N"] = N
  );
  
  List Index = List::create(
    _["unq_id"] = unq_id,
    _["unq_tm"] = unq_tm,
    _["unq_x"] = unq_x,
    _["id_index"] = id_index,
    _["tm_index"] = tm_index,
    _["x_index"] = x_index
  );
  
  List BS = List::create(
    _["Bt"] = Bt,
    _["Bx"] = Bx,
    _["Bt_H"] = Bt_H,
    _["Bx_K"] = Bx_K,
    _["Bxt"] = Bxt,
    _["Bx_Scale"] = Bx_Scale
  );
  
  List Regulate = List::create(
    _["M"] = M_New,
    _["nu"] = nu,
    _["Lambda_Ridge_Vec"] = Lambda_Ridge_Vec,
    _["Shrink"] = Shrink,
    _["Ridge_Penalty"] = Ridge_Penalty,
    _["Lambda_Scale"] = Lambda_Scale,
    _["NLambda"] = NLambda,
    _["lower_perc"] = lower_perc,
    _["upper_perc"] = upper_perc
  );
  
  List Beta_Estimate = List::create(
    _["Beta"] = Beta,
    _["Beta_Hat_List"] = Beta_Hat_List,
    _["Beta_Hat_List_Iter"] = Beta_Hat_List_Iter,
    _["Sum_Beta_Hat_List"] = Sum_Beta_Hat_List,
    _["Tm_Beta_C"] = Tm_Beta_C,
    _["lower_Beta_Hat_Noise"] = lower_Beta_Hat_Noise,
    _["upper_Beta_Hat_Noise"] = upper_Beta_Hat_Noise,
    _["List_Trace_Bxt_gm"] = List_Trace_Bxt_gm
  );
  
  return List::create(
    _["Data"] = Data,
    _["Dimensions"] = Dimensions,
    _["Index"] = Index,
    _["BS"] = BS,
    _["Regulate"] = Regulate,
    _["Beta_Estimate"] = Beta_Estimate,
    _["Error_Rate"] = Error_Rate,
    _["Variable_Select"] = Variable_Select,
    _["Response_Select"] = Response_Select,
    _["mu"] = mu,
    _["mu_List"] = mu_List,
    _["Phi"] = Phi,
    _["Rho"] = Rho,
    _["Lambda_List"] = Lambda_List,
    _["mu_zero"] = mu_zero,
    _["Vec_zero"] = Vec_zero
  );
}

// [[Rcpp::export]]
List predict_BoostMLR_C(NumericMatrix Org_x,
                        NumericVector tm,
                        NumericVector id,
                        NumericMatrix Org_y,
                        NumericVector x_Mean,
                        NumericVector x_Std_Error,
                        NumericVector y_Mean,
                        NumericVector y_Std_Error,
                        int K,
                        int L,
                        int H,
                        IntegerVector Dk,
                        NumericVector unq_tm,
                        List unq_x,
                        NumericMatrix Bt,
                        List Bx,
                        LogicalVector UseRaw,
                        NumericMatrix Time_Add_New,
                        LogicalVector Time_Unmatch,
                        List Beta,
                        List Beta_Hat_List,
                        bool testFlag,
                        int M,
                        double nu,
                        bool Time_Varying,
                        bool vimpFlag,
                        bool vimpFlag_Coef,
                        double eps)
{                           
  
  NumericVector unq_id = sort_unique_C_NA(id);
  int n = unq_id.size();
  
  List id_index(n);
  for(int i = 0; i < n; ++i){
    id_index[i] = Which_C_NA(unq_id[i],id);
  }
  
  IntegerVector ni(n);
  for(int i = 0; i < n; ++i){
    IntegerVector id_index_Temp = id_index[i];
    ni[i] = id_index_Temp.size();
  }
  
  int N = Org_y.nrow();
  
  IntegerVector unlist_id_index(N);
  int count = 0;
  for(int i = 0; i < n; ++i){
    NumericVector id_index_Temp = id_index[i];
    for(int j = 0; j < ni[i]; ++j){
      unlist_id_index[count] = id_index_Temp[j];
      count = count + 1;
    }
  }
  
  NumericVector id_New(N);
  for(int ii = 0; ii < N; ++ii){
    id_New[ii] = id[ unlist_id_index[ii] ];
  }
  
  LogicalVector id_Match(N);
  for(int ii = 0; ii < N; ++ii){
    id_Match[ii] = (id[ii] == id_New[ii]);
  }
  
  if(is_false(all( id_Match == TRUE) )){
    NumericVector tm_New(N);
    NumericMatrix Org_x_New(N,K);
    NumericMatrix Org_y_New(N,L);
    for(int ii = 0; ii < N; ++ii){
      tm_New[ii] = tm[ unlist_id_index[ii]  ];
      id_New[ii] = id[ unlist_id_index[ii]  ];
      for(int k = 0; k < K; ++k){
        Org_x_New(ii,k) = Org_x( unlist_id_index[ii], k  );
      }
      for(int l = 0; l < L; ++l){
        Org_y_New(ii,l) = Org_y( unlist_id_index[ii], l  );
      }
    }
    tm = tm_New;
    id = id_New;
    Org_x = Org_x_New;
    Org_y = Org_y_New;  
  }
  
  List tm_index(n);
  for(int i = 0; i < n; ++i){
    NumericVector tm_Temp(ni[i]);
    IntegerVector id_index_Temp = id_index[i];
    for(int j = 0; j < ni[i]; ++j){
      tm_Temp[j] = tm[ id_index_Temp[j] ];
    }
    tm_index[i] = Approx_Match_C_NA(tm_Temp,unq_tm);
  }
  
  List Bt_H(H);
  for(int h = 0; h < H; ++h){
    List Bt_n(n);
    for(int i = 0; i < n; ++i){
      IntegerVector tm_index_Temp = tm_index[i];
      NumericVector Bt_ni(ni[i]);
      for(int j = 0; j < ni[i]; ++j){
        if(IntegerVector::is_na(tm_index_Temp[j])){
          Bt_ni[j] = NA_REAL;
        }else
        {
          Bt_ni[j] = Bt(tm_index_Temp[j], h);          
        }
      }
      Bt_n[i] = Bt_ni;
    }
    Bt_H[h] = Bt_n;
  }
  
  int n_unq_tm = Bt.nrow(); 
  IntegerVector Index_Bt(n_unq_tm);
  for(int i = 0; i < n_unq_tm; ++i){
    Index_Bt[i] = i;
  }
  
  List unq_x_New(K);
  for(int k = 0; k < K; ++k){
    NumericVector unq_x_Temp = unq_x[k];
    int unq_x_Temp_size = unq_x_Temp.size();
    NumericVector unq_x_New_Temp(unq_x_Temp_size);
    for(int i = 0; i < unq_x_Temp_size; ++i){
      unq_x_New_Temp[i] = (unq_x_Temp[i]*x_Std_Error[k]) + x_Mean[k];
    }
    unq_x_New[k] = unq_x_New_Temp;
  }
  
  List x_index(K);
  for(int k = 0; k < K; ++k){
    if(UseRaw[k]){
      x_index[k] = NA_INTEGER;
    }
    else
    {
      NumericVector unq_x_Temp = unq_x_New[k];
      List x_index_Subject(n);
      for(int i = 0; i < n; ++i){
        NumericVector x_Temp(ni[i]);
        IntegerVector id_index_Temp = id_index[i];
        for(int j = 0; j < ni[i]; ++j){
          x_Temp[j] = Org_x(id_index_Temp[j] , k);
        }
        x_index_Subject[i] = Approx_Match_C_NA(x_Temp,unq_x_Temp);
      }
      x_index[k] = x_index_Subject;
    }
  }
  
  List Bx_K(K);
  for(int k = 0; k < K; ++k){
    NumericMatrix Bx_K_Temp = Bx[k];
    List Bx_Dk(Dk[k]);
    if(UseRaw[k]){
      for(int d = 0; d < Dk[k]; ++d){
        List Bx_n(n);
        for(int i = 0; i < n; ++i){
          IntegerVector id_index_Temp = id_index[i];
          NumericVector Bx_ni(ni[i]);
          for(int j = 0; j < ni[i]; ++j){
            Bx_ni[j] = (Org_x(id_index_Temp[j],k) - x_Mean[k])/x_Std_Error[k];
          }
          Bx_n[i] = Bx_ni;
        }
        Bx_Dk[d] = Bx_n;
      }
    }
    else{
      List x_index_Temp = x_index[k];
      for(int d = 0; d < Dk[k]; ++d){
        List Bx_n(n);
        for(int i = 0; i < n; ++i){
          IntegerVector x_index_Temp_n = x_index_Temp[i];
          NumericVector Bx_ni(ni[i]);
          for(int j = 0; j < ni[i]; ++j){
            if(IntegerVector::is_na(x_index_Temp_n[j])){
              Bx_ni[j] = NA_REAL;
            }else
            {
              Bx_ni[j] = Bx_K_Temp(x_index_Temp_n[j],d);  
            }
          }
          Bx_n[i] = Bx_ni;
        }
        Bx_Dk[d] = Bx_n;
      }
    }
    Bx_K[k] = Bx_Dk;
  }
  
  List Bxt(K);
  for(int k = 0; k < K; ++k){
    List Bx_Temp_K = Bx_K[k];
    List Bxt_Dk(Dk[k]);
    for(int d = 0; d < Dk[k]; ++d){
      List Bx_Temp_Dk = Bx_Temp_K[d];
      List Bxt_H(H);
      for(int h = 0; h < H; ++h){
        List Bt_Temp_H = Bt_H[h];
        List Bxt_n(n);
        for(int i = 0; i < n; ++i){
          NumericVector Bx_Temp_n = Bx_Temp_Dk[i];
          NumericVector Bt_Temp_n = Bt_Temp_H[i];
          NumericVector Bxt_ni(ni[i]);
          for(int j = 0; j < ni[i]; ++j){
            Bxt_ni[j] = Bx_Temp_n[j]*Bt_Temp_n[j];
          }
          Bxt_n[i] = Bxt_ni;
        }
        Bxt_H[h] = Bxt_n;
      }
      Bxt_Dk[d] = Bxt_H;
    }
    Bxt[k] = Bxt_Dk;
  }
  
  NumericMatrix mu_zero(N,L);
  for(int l = 0; l < L; ++l){
    for(int i = 0; i < N; ++i){
      mu_zero(i,l) = 0;
    }
  }
  
  NumericMatrix mu_zero_vec(N,1);
  for(int i = 0; i < N; ++i){
    mu_zero_vec(i,0) = 0;
  }

  
  NumericMatrix Org_mu(N,L);
  List mu_List(M);
  for(int m = 0; m < M; ++m){
    mu_List[m] = NA_REAL;
  }
  NumericMatrix Org_mu_Mopt(N,L);
  NumericVector Vec_zero(N);
  for(int i = 0; i < N; ++i){
    Vec_zero[i] = 0.0;
  }
  
  NumericMatrix Error_Rate(M,L);
  for(int m = 0; m < M; ++m){
    for(int l = 0; l < L; ++l){
      Error_Rate(m,l) = NA_REAL;
    }
  }
  
  IntegerVector Mopt(L);
  NumericVector rmse(L);
  for(int l = 0; l < L; ++l){
    Mopt[l] = NA_INTEGER;
    rmse[l] = NA_REAL;
  }
  int Mopt_Max;
  
  List vimp(L);
  List vimp_Coef(L);
  NumericMatrix mu_vimp(N,1);
  NumericMatrix Org_mu_vimp(N,1);
  double Error_Rate_vimp;
  NumericMatrix Org_x_New(N,K);
    
  if(testFlag){
    
    for(int m = 0;m < M; ++m){
      
      List Beta_Hat = Beta_Hat_List[m];
      
      NumericMatrix mu = mu_zero;
      for(int k = 0; k < K; ++k){
        List Bxt_Temp_K = Bxt[k];
        List Beta_Temp_K = Beta_Hat[k];
        for(int d = 0; d < Dk[k]; ++d){
          List Bxt_Temp_Dk = Bxt_Temp_K[d];
          List Beta_Temp_Dk = Beta_Temp_K[d];
          for(int h = 0; h < H; ++h){
            List Bxt_Temp_H = Bxt_Temp_Dk[h];
            NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
            NumericMatrix mu_Hat(N,L);
            for(int l = 0; l < L; ++l){
              for(int i = 0; i < n; ++i){
                IntegerVector id_index_Temp = id_index[i];
                NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                for(int j = 0; j < ni[i]; ++j){
                  mu_Hat( id_index_Temp[j], l) = Bxt_Temp_n[j]*Beta_Temp_H[l];
                }
              }
            }
            mu = Matrix_Sum_C_NA(mu, mu_Hat);
          }
        }
      }
      
      for(int l = 0; l < L; ++l){
        for(int i = 0; i < N; ++i){
          Org_mu(i,l) = ( mu(i,l) * y_Std_Error[l] ) + y_Mean[l];
        }
      }
      
      /* 
       Work on code: mu_List[m] = Org_mu; in future. mu_List save Org_mu 
       at every iteration, however when you come out of the loop, all elements of list
       will have Org_mu corresponding to the last boosting iteration.
       For now, I have different way of calculating mu_Mopt rather
       than using mu_List. mu_List will not be save for output until
       I fixed the issue.
       */
      mu_List[m] = Org_mu;
      
      for(int l = 0; l < L; ++l){
        NumericVector Org_y_Temp(N);
        NumericVector Org_mu_Temp(N);
        for(int i = 0; i < N; ++i){
          Org_y_Temp[i] = Org_y(i,l);
          Org_mu_Temp[i] = Org_mu(i,l);
        }
        Error_Rate(m,l) = l2Dist_Vector_C_NA(Org_y_Temp, Org_mu_Temp, id_index)/l2Dist_Vector_C_NA(Org_y_Temp, Vec_zero, id_index); 
      }
      
    }
    
  } else
  {
    NumericMatrix mu = mu_zero;
    for(int k = 0; k < K; ++k){
      List Bxt_Temp_K = Bxt[k];
      List Beta_Temp_K = Beta[k];
      for(int d = 0; d < Dk[k]; ++d){
        List Bxt_Temp_Dk = Bxt_Temp_K[d];
        List Beta_Temp_Dk = Beta_Temp_K[d];
        for(int h = 0; h < H; ++h){
          List Bxt_Temp_H = Bxt_Temp_Dk[h];
          NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
          NumericMatrix mu_Hat(N,L);
          for(int l = 0; l < L; ++l){
            for(int i = 0; i < n; ++i){
              IntegerVector id_index_Temp = id_index[i];
              NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
              for(int j = 0; j < ni[i]; ++j){
                mu_Hat( id_index_Temp[j], l) = Bxt_Temp_n[j]*Beta_Temp_H[l];
              }
            }
          }
          mu = Matrix_Sum_C_NA(mu, mu_Hat);
        }
      }
    }
    
    for(int l = 0; l < L; ++l){
      for(int i = 0; i < N; ++i){
        Org_mu(i,l) = ( mu(i,l) * y_Std_Error[l] ) + y_Mean[l];
      }
    }
  }
  
  if(testFlag){
    for(int l = 0; l < L; ++l){
      NumericVector Error_Rate_L(M);
      for(int m = 0; m < M; ++m){
        Error_Rate_L[m] = Error_Rate(m,l);
      }
      double Min_Error_Rate = min(Error_Rate_L);
      NumericVector Diff_Error(M);
      for(int m = 0; m < M; ++m){
        Diff_Error[m] = ( Error_Rate_L[m] - Min_Error_Rate );
      }
      Diff_Error = abs(Diff_Error);
      for(int m = 0; m < M; ++m){
        if(Diff_Error[m] < eps){
          Mopt[l] = m;
          rmse[l] = Error_Rate(Mopt[l],l);
          break;
        }
      }
    }
    Mopt_Max = max(Mopt);
  }
  
  if(testFlag){
    NumericMatrix mu_Mopt = mu_zero;
    for(int l = 0; l < L; ++l){
    int Mopt_Temp = Mopt[l];
    List Beta_Hat = Beta_Hat_List[Mopt_Temp];
    for(int k = 0; k < K; ++k){
      List Bxt_Temp_K = Bxt[k];
      List Beta_Temp_K = Beta_Hat[k];
      for(int d = 0; d < Dk[k]; ++d){
        List Bxt_Temp_Dk = Bxt_Temp_K[d];
        List Beta_Temp_Dk = Beta_Temp_K[d];
        for(int h = 0; h < H; ++h){
          List Bxt_Temp_H = Bxt_Temp_Dk[h];
          NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
          NumericMatrix mu_Hat(N,L);
            for(int i = 0; i < n; ++i){
              IntegerVector id_index_Temp = id_index[i];
              NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
              for(int j = 0; j < ni[i]; ++j){
                mu_Hat( id_index_Temp[j], l) = Bxt_Temp_n[j]*Beta_Temp_H[l];
              }
            }
          mu_Mopt = Matrix_Sum_C_NA(mu_Mopt, mu_Hat);
        }
      }
    }
  }
    for(int l = 0; l < L; ++l){
      for(int i = 0; i < N; ++i){
        Org_mu_Mopt(i,l) = ( mu_Mopt(i,l) * y_Std_Error[l] ) + y_Mean[l];
      }
    }
}  
  
  //List Beta_Mopt = Beta_Hat_List[Mopt_Max];
  
  if(vimpFlag){
    
    for(int l = 0; l < L; ++l){
    List Beta_Mopt = Beta_Hat_List[ Mopt[l]  ];
    int H_vimp;
    if(Time_Varying){
      H_vimp = H+1;
    }
    else {
      H_vimp = H;
    }
    NumericMatrix vimp_Temp(K,H_vimp);
    
    for(int kk = 0; kk < K; ++kk){
     
     List vimp_Bxt(K);
     for(int k = 0; k < K; ++k){
       if(k == kk){
         NumericVector Org_x_kk(N);
         for(int i = 0; i < N; ++i){
           Org_x_kk[i] = Org_x(i,kk);
         }
         Org_x_kk = randomShuffle(Org_x_kk);
         NumericMatrix Bx_K_Temp = Bx[k];
         List Bx_Dk(Dk[k]);
         if(UseRaw[k]){
           for(int d = 0; d < Dk[k]; ++d){
             List Bx_n(n);
             for(int i = 0; i < n; ++i){
               IntegerVector id_index_Temp = id_index[i];
               NumericVector Bx_ni(ni[i]);
               for(int j = 0; j < ni[i]; ++j){
                 Bx_ni[j] = ( Org_x_kk[ id_index_Temp[j] ] - x_Mean[k])/x_Std_Error[k];
               }
               Bx_n[i] = Bx_ni;
             }
             Bx_Dk[d] = Bx_n;
           }
         }else
         {
           NumericVector unq_x_Temp = unq_x_New[k];
           List x_index_Subject(n);
           for(int i = 0; i < n; ++i){
             NumericVector x_Temp(ni[i]);
             IntegerVector id_index_Temp = id_index[i];
             for(int j = 0; j < ni[i]; ++j){
               x_Temp[j] = Org_x_kk[id_index_Temp[j]];
             }
             x_index_Subject[i] = Approx_Match_C_NA(x_Temp,unq_x_Temp);
           }
           for(int d = 0; d < Dk[k]; ++d){
             List Bx_n(n);
             for(int i = 0; i < n; ++i){
               IntegerVector x_index_Temp_n = x_index_Subject[i];
               NumericVector Bx_ni(ni[i]);
               for(int j = 0; j < ni[i]; ++j){
                 if(IntegerVector::is_na(x_index_Temp_n[j])){
                   Bx_ni[j] = NA_REAL;
                 }else
                 {
                   Bx_ni[j] = Bx_K_Temp(x_index_Temp_n[j],d);  
                 }
               }
               Bx_n[i] = Bx_ni;
             }
             Bx_Dk[d] = Bx_n;
           }
         }
         
         List Bxt_Dk(Dk[k]);
         for(int d = 0; d < Dk[k]; ++d){
           List Bx_Temp_Dk = Bx_Dk[d];
           List Bxt_H(H);
           for(int h = 0; h < H; ++h){
             List Bt_Temp_H = Bt_H[h];
             List Bxt_n(n);
             for(int i = 0; i < n; ++i){
               NumericVector Bx_Temp_n = Bx_Temp_Dk[i];
               NumericVector Bt_Temp_n = Bt_Temp_H[i];
               NumericVector Bxt_ni(ni[i]);
               for(int j = 0; j < ni[i]; ++j){
                 Bxt_ni[j] = Bx_Temp_n[j]*Bt_Temp_n[j];
               }
               Bxt_n[i] = Bxt_ni;
             }
             Bxt_H[h] = Bxt_n;
           }
           Bxt_Dk[d] = Bxt_H;
         }
         vimp_Bxt[k] = Bxt_Dk;
         
       }else
       {
         vimp_Bxt[k] = Bxt[k];
       }
     }
     
     mu_vimp = mu_zero_vec;
     for(int k = 0; k < K; ++k){
       List Bxt_Temp_K = vimp_Bxt[k];
       List Beta_Temp_K = Beta_Mopt[k];
       for(int d = 0; d < Dk[k]; ++d){
         List Bxt_Temp_Dk = Bxt_Temp_K[d];
         List Beta_Temp_Dk = Beta_Temp_K[d];
         for(int h = 0; h < H; ++h){
           List Bxt_Temp_H = Bxt_Temp_Dk[h];
           NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
           NumericMatrix mu_Hat(N,1);
             for(int i = 0; i < n; ++i){
               IntegerVector id_index_Temp = id_index[i];
               NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
               for(int j = 0; j < ni[i]; ++j){
                 mu_Hat( id_index_Temp[j], 0) = Bxt_Temp_n[j]*Beta_Temp_H[l];
               }
             }
           mu_vimp = Matrix_Sum_C_NA(mu_vimp, mu_Hat);
         }
       }
     }
     
       for(int i = 0; i < N; ++i){
         Org_mu_vimp(i,0) = ( mu_vimp(i,0) * y_Std_Error[l] ) + y_Mean[l];
       }
     
       NumericVector Org_y_Temp(N);
       NumericVector Org_mu_Temp(N);
       for(int i = 0; i < N; ++i){
         Org_y_Temp[i] = Org_y(i,l);
         Org_mu_Temp[i] = Org_mu_vimp(i,0);
       }
       Error_Rate_vimp = l2Dist_Vector_C_NA(Org_y_Temp, Org_mu_Temp, id_index)/l2Dist_Vector_C_NA(Org_y_Temp, Vec_zero, id_index);
       vimp_Temp(kk,0) = (Error_Rate_vimp - rmse[l])/rmse[l];
    }
    
    if(Time_Varying){
    for(int kk = 0; kk < K; ++kk){
      for(int hh = 0; hh < H; ++hh){
      List vimp_Bxt(K);
      for(int k = 0; k < K; ++k){
        List vimp_Bt_H(H);
        if(k == kk){
          
          IntegerVector Index_Bt_noise = int_randomShuffle(Index_Bt);
          
          NumericMatrix vimp_Bt(n_unq_tm,H);
          for(int i = 0; i < n_unq_tm; ++i){
            for(int j = 0; j < H; ++j){
              if(j == hh){
                vimp_Bt(i,j) = Bt(Index_Bt_noise[i],j);  
              }else
              {
                vimp_Bt(i,j) = Bt(Index_Bt[i],j); 
              }
            }
          }
          for(int h = 0; h < H; ++h){
            List Bt_n(n);
            for(int i = 0; i < n; ++i){
              IntegerVector tm_index_Temp = tm_index[i];
              NumericVector Bt_ni(ni[i]);
              for(int j = 0; j < ni[i]; ++j){
                if(IntegerVector::is_na(tm_index_Temp[j])){
                  Bt_ni[j] = NA_REAL;
                }else
                {
                  Bt_ni[j] = vimp_Bt(tm_index_Temp[j], h);  
                }
              }
              Bt_n[i] = Bt_ni;
            }
            vimp_Bt_H[h] = Bt_n;
          }
          
          List Bx_Temp_K = Bx_K[k];
          List Bxt_Dk(Dk[k]);
          for(int d = 0; d < Dk[k]; ++d){
          List Bx_Temp_Dk = Bx_Temp_K[d]; 
          List Bxt_H(H);
          for(int h = 0; h < H; ++h){
          List Bt_Temp_H = vimp_Bt_H[h];
            List Bxt_n(n);
            for(int i = 0; i < n; ++i){
              NumericVector Bx_Temp_n = Bx_Temp_Dk[i];
              NumericVector Bt_Temp_n = Bt_Temp_H[i];
              NumericVector Bxt_ni(ni[i]);
              for(int j = 0; j < ni[i]; ++j){
              Bxt_ni[j] = Bx_Temp_n[j]*Bt_Temp_n[j];  
              }
              Bxt_n[i] = Bxt_ni;
            }
            Bxt_H[h] = Bxt_n;
          }
          Bxt_Dk[d] = Bxt_H;
          }
          vimp_Bxt[k] = Bxt_Dk;
        }else
        {
          vimp_Bxt[k] = Bxt[k];
        }
      }
      
      mu_vimp = mu_zero_vec;
      for(int k = 0; k < K; ++k){
        List Bxt_Temp_K = vimp_Bxt[k];
        List Beta_Temp_K = Beta_Mopt[k];
        for(int d = 0; d < Dk[k]; ++d){
          List Bxt_Temp_Dk = Bxt_Temp_K[d];
          List Beta_Temp_Dk = Beta_Temp_K[d];
          for(int h = 0; h < H; ++h){
            List Bxt_Temp_H = Bxt_Temp_Dk[h];
            NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
            NumericMatrix mu_Hat(N,1);
              for(int i = 0; i < n; ++i){
                IntegerVector id_index_Temp = id_index[i];
                NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                for(int j = 0; j < ni[i]; ++j){
                  mu_Hat( id_index_Temp[j], 0) = Bxt_Temp_n[j]*Beta_Temp_H[l];
                }
              }
            mu_vimp = Matrix_Sum_C_NA(mu_vimp, mu_Hat);
          }
        }
      }
      
        for(int i = 0; i < N; ++i){
          Org_mu_vimp(i,0) = ( mu_vimp(i,0) * y_Std_Error[l] ) + y_Mean[l];
        }
      
        NumericVector Org_y_Temp(N);
        NumericVector Org_mu_Temp(N);
        
        for(int i = 0; i < N; ++i){
          Org_y_Temp[i] = Org_y(i,l);
          Org_mu_Temp[i] = Org_mu_vimp(i,0);
        }
        
        Error_Rate_vimp = l2Dist_Vector_C_NA(Org_y_Temp, Org_mu_Temp, id_index)/l2Dist_Vector_C_NA(Org_y_Temp, Vec_zero, id_index);
        vimp_Temp(kk, hh+1 ) = (Error_Rate_vimp - rmse[l])/rmse[l];
      }
    }
  }
    vimp[l] = vimp_Temp;
  }
}
  

  
  if(vimpFlag_Coef){
    
    for(int l = 0; l < L; ++l){
      List Beta_Mopt = Beta_Hat_List[ Mopt[l]  ];
      int H_vimp;
      if(Time_Varying){
        H_vimp = H+1;
      }
      else {
        H_vimp = H;
      }
      NumericMatrix vimp_Temp(K,H_vimp);
      
      for(int kk = 0; kk < K; ++kk){
        
        mu_vimp = mu_zero_vec;
        for(int k = 0; k < K; ++k){
          List Bxt_Temp_K = Bxt[k];
          List Beta_Temp_K = Beta_Mopt[k];
          for(int d = 0; d < Dk[k]; ++d){
            List Bxt_Temp_Dk = Bxt_Temp_K[d];
            List Beta_Temp_Dk = Beta_Temp_K[d];
            for(int h = 0; h < H; ++h){
              List Bxt_Temp_H = Bxt_Temp_Dk[h];
              NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
              double Beta_Temp_L;
              if(k == kk){
                Beta_Temp_L = 0;  
              } else
              {
                Beta_Temp_L = Beta_Temp_H[l];
              }
              NumericMatrix mu_Hat(N,1);
              for(int i = 0; i < n; ++i){
                IntegerVector id_index_Temp = id_index[i];
                NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                for(int j = 0; j < ni[i]; ++j){
                  mu_Hat( id_index_Temp[j], 0) = Bxt_Temp_n[j]*Beta_Temp_L;
                }
              }
              mu_vimp = Matrix_Sum_C_NA(mu_vimp, mu_Hat);
            }
          }
        }
        
        for(int i = 0; i < N; ++i){
          Org_mu_vimp(i,0) = ( mu_vimp(i,0) * y_Std_Error[l] ) + y_Mean[l];
        }
        
        NumericVector Org_y_Temp(N);
        NumericVector Org_mu_Temp(N);
        for(int i = 0; i < N; ++i){
          Org_y_Temp[i] = Org_y(i,l);
          Org_mu_Temp[i] = Org_mu_vimp(i,0);
        }
        Error_Rate_vimp = l2Dist_Vector_C_NA(Org_y_Temp, Org_mu_Temp, id_index)/l2Dist_Vector_C_NA(Org_y_Temp, Vec_zero, id_index);
        vimp_Temp(kk,0) = (Error_Rate_vimp - rmse[l])/rmse[l];
      }
      
      if(Time_Varying){
        for(int kk = 0; kk < K; ++kk){
          for(int hh = 0; hh < H; ++hh){
            
            mu_vimp = mu_zero_vec;
            for(int k = 0; k < K; ++k){
              List Bxt_Temp_K = Bxt[k];
              List Beta_Temp_K = Beta_Mopt[k];
              for(int d = 0; d < Dk[k]; ++d){
                List Bxt_Temp_Dk = Bxt_Temp_K[d];
                List Beta_Temp_Dk = Beta_Temp_K[d];
                for(int h = 0; h < H; ++h){
                  List Bxt_Temp_H = Bxt_Temp_Dk[h];
                  NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
                  double Beta_Temp_L;
                  if(k == kk && h == hh){
                    Beta_Temp_L = 0;  
                  } else
                  {
                    Beta_Temp_L = Beta_Temp_H[l];
                  }
                  NumericMatrix mu_Hat(N,1);
                  for(int i = 0; i < n; ++i){
                    IntegerVector id_index_Temp = id_index[i];
                    NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                    for(int j = 0; j < ni[i]; ++j){
                      mu_Hat( id_index_Temp[j], 0) = Bxt_Temp_n[j]*Beta_Temp_L;
                    }
                  }
                  mu_vimp = Matrix_Sum_C_NA(mu_vimp, mu_Hat);
                }
              }
            }
            
            for(int i = 0; i < N; ++i){
              Org_mu_vimp(i,0) = ( mu_vimp(i,0) * y_Std_Error[l] ) + y_Mean[l];
            }
            
            NumericVector Org_y_Temp(N);
            NumericVector Org_mu_Temp(N);
            
            for(int i = 0; i < N; ++i){
              Org_y_Temp[i] = Org_y(i,l);
              Org_mu_Temp[i] = Org_mu_vimp(i,0);
            }
            
            Error_Rate_vimp = l2Dist_Vector_C_NA(Org_y_Temp, Org_mu_Temp, id_index)/l2Dist_Vector_C_NA(Org_y_Temp, Vec_zero, id_index);
            vimp_Temp(kk, hh+1 ) = (Error_Rate_vimp - rmse[l])/rmse[l];
          }
        }
      }
      vimp_Coef[l] = vimp_Temp;
    }
  }  
  
  
  
  List Data = List::create(
    _["Org_x"] = Org_x, 
    _["Org_y"] = Org_y,
    _["id"] = id,
    _["tm"] = tm,
    _["x_Mean"] = x_Mean,
    _["x_Std_Error"] = x_Std_Error,
    _["y_Mean"] = y_Mean,
    _["y_Std_Error"] = y_Std_Error
  );
  
  List Dimensions = List::create(
    _["n"] = n,
    _["K"] = K,
    _["L"] = L,
    _["H"] = H,
    _["Dk"] = Dk,
    _["ni"] = ni,
    _["N"] = N,
    _["n_unq_tm"] = n_unq_tm
  );
  
  List Index = List::create(
    _["unq_id"] = unq_id,
    _["unq_tm"] = unq_tm,
    _["unq_x"] = unq_x,
    _["id_index"] = id_index,
    _["tm_index"] = tm_index,
    _["x_index"] = x_index,
    _["unq_x_New"] = unq_x_New,
    _["Index_Bt"] = Index_Bt
  );
  
  List BS = List::create(
    _["Bt"] = Bt,
    _["Bxt"] = Bxt,
    _["Bx_K"] = Bx_K,
    _["Bt_H"] = Bt_H,
    _["Bx"] = Bx
  );
  
  List Pred_Object = List::create(
    _["Vec_zero"] = Vec_zero,
    _["mu_zero_vec"] = mu_zero_vec
  );
  
  return List::create(
    _["Data"] = Data,
    _["UseRaw"] = UseRaw,
    _["Dimensions"] = Dimensions,
    _["Index"] = Index,
    _["BS"] = BS,
    _["Beta_Hat_List"] = Beta_Hat_List,
    _["Org_mu"] = Org_mu,
    _["Org_mu_Mopt"] = Org_mu_Mopt,
    _["Error_Rate"] = Error_Rate,
    _["Mopt"] = Mopt,
    _["rmse"] = rmse,
    _["vimp"] = vimp,
    _["vimp_Coef"] = vimp_Coef,
    _["Pred_Object"] = Pred_Object
  );
}


// [[Rcpp::export]]
List vimp_BoostMLR_C(NumericMatrix Org_x,
                     NumericMatrix Org_y,
                     NumericVector tm,
                     NumericVector id,
                     NumericVector x_Mean,
                     NumericVector x_Std_Error,
                     NumericVector y_Mean,
                     NumericVector y_Std_Error,
                     int n,
                     IntegerVector ni,
                     int N,
                     int L,
                     int K,
                     int p,
                     int H,
                     IntegerVector Dk,
                     int n_unq_tm,
                     LogicalVector UseRaw,
                     List id_index,
                     List tm_index,
                     List unq_x_New,
                     IntegerVector Index_Bt,
                     IntegerVector vimp_set,
                     bool joint,
                     NumericMatrix Bt,
                     List Bt_H,
                     List Bx,
                     List Bxt,
                     List Bx_K,
                     List Beta_Hat_List,
                     IntegerVector Mopt,
                     double nu,
                     NumericVector rmse,
                     bool Time_Varying,
                     NumericVector Vec_zero,
                     NumericMatrix mu_zero_vec)
{
  List vimp(L);
  NumericMatrix mu_vimp(N,1);
  NumericMatrix Org_mu_vimp(N,1);
  NumericMatrix Org_x_Noise(N,K);
  double Error_Rate_vimp;
  IntegerVector Index_N(N);
  for(int i = 0;i < N; ++i){
    Index_N[i] = i;
  }
  for(int l = 0; l < L; ++l){
    List Beta_Mopt = Beta_Hat_List[ Mopt[l]  ];
    int H_vimp;
    if(Time_Varying){
      H_vimp = H+1;
    }
    else {
      H_vimp = H;
    }
    NumericMatrix vimp_Temp(p,H_vimp);
    
    for(int kk = 0; kk < p; ++kk){
      
      IntegerVector Index_Permute(N);
      Index_Permute = int_randomShuffle(Index_N);
      Org_x_Noise = Org_x;
      int n_vimp_set = vimp_set.size();
      for(int i = 0; i < N; ++i){
        if(joint){
          for(int k_vimp = 0; k_vimp < n_vimp_set; ++k_vimp){
            Org_x_Noise(i, vimp_set[k_vimp]) = Org_x(Index_Permute[i], vimp_set[k_vimp]);
          }
        }
        else {
          Org_x_Noise(i,vimp_set[kk]) = Org_x(Index_Permute[i], vimp_set[kk]);
        }
      }
      List vimp_Bxt(K);
      if(joint){
        for(int k = 0; k < K; ++k){
          if(is_true(any( vimp_set == k ))){
            NumericMatrix Bx_K_Temp = Bx[k];
            List Bx_Dk(Dk[k]);
            if(UseRaw[k]){
              for(int d = 0; d < Dk[k]; ++d){
                List Bx_n(n);
                for(int i = 0; i < n; ++i){
                  IntegerVector id_index_Temp = id_index[i];
                  NumericVector Bx_ni(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    Bx_ni[j] = ( Org_x_Noise(id_index_Temp[j],k) - x_Mean[k])/x_Std_Error[k];
                  }
                  Bx_n[i] = Bx_ni;
                }
                Bx_Dk[d] = Bx_n;
              }
            } else
            {
              NumericVector unq_x_Temp = unq_x_New[k];
              List x_index_Subject(n);
              for(int i = 0; i < n; ++i){
                NumericVector x_Temp(ni[i]);
                IntegerVector id_index_Temp = id_index[i];
                for(int j = 0; j < ni[i]; ++j){
                  x_Temp[j] = Org_x_Noise(id_index_Temp[j],k);
                }
                x_index_Subject[i] = Approx_Match_C_NA(x_Temp,unq_x_Temp);
              }
              for(int d = 0; d < Dk[k]; ++d){
                List Bx_n(n);
                for(int i = 0; i < n; ++i){
                  IntegerVector x_index_Temp_n = x_index_Subject[i];
                  NumericVector Bx_ni(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    if(IntegerVector::is_na(x_index_Temp_n[j])){
                      Bx_ni[j] = NA_REAL;
                    }else
                    {
                      Bx_ni[j] = Bx_K_Temp(x_index_Temp_n[j],d);
                    }
                  }
                  Bx_n[i] = Bx_ni;
                }
                Bx_Dk[d] = Bx_n;
              }
            }
            
            List Bxt_Dk(Dk[k]);
            for(int d = 0; d < Dk[k]; ++d){
              List Bx_Temp_Dk = Bx_Dk[d];
              List Bxt_H(H);
              for(int h = 0; h < H; ++h){
                List Bt_Temp_H = Bt_H[h];
                List Bxt_n(n);
                for(int i = 0; i < n; ++i){
                  NumericVector Bx_Temp_n = Bx_Temp_Dk[i];
                  NumericVector Bt_Temp_n = Bt_Temp_H[i];
                  NumericVector Bxt_ni(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    Bxt_ni[j] = Bx_Temp_n[j]*Bt_Temp_n[j];
                  }
                  Bxt_n[i] = Bxt_ni;
                }
                Bxt_H[h] = Bxt_n;
              }
              Bxt_Dk[d] = Bxt_H;
            }
            vimp_Bxt[k] = Bxt_Dk;
            
          }else
          {
            vimp_Bxt[k] = Bxt[k];
          }
        }
      }
      else {
        for(int k = 0; k < K; ++k){
          if(k == vimp_set[kk]){
            NumericMatrix Bx_K_Temp = Bx[k];
            List Bx_Dk(Dk[k]);
            if(UseRaw[k]){
              for(int d = 0; d < Dk[k]; ++d){
                List Bx_n(n);
                for(int i = 0; i < n; ++i){
                  IntegerVector id_index_Temp = id_index[i];
                  NumericVector Bx_ni(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    Bx_ni[j] = ( Org_x_Noise(id_index_Temp[j],k) - x_Mean[k])/x_Std_Error[k];
                  }
                  Bx_n[i] = Bx_ni;
                }
                Bx_Dk[d] = Bx_n;
              }
            } else
            {
              NumericVector unq_x_Temp = unq_x_New[k];
              List x_index_Subject(n);
              for(int i = 0; i < n; ++i){
                NumericVector x_Temp(ni[i]);
                IntegerVector id_index_Temp = id_index[i];
                for(int j = 0; j < ni[i]; ++j){
                  x_Temp[j] = Org_x_Noise(id_index_Temp[j],k);
                }
                x_index_Subject[i] = Approx_Match_C_NA(x_Temp,unq_x_Temp);
              }
              for(int d = 0; d < Dk[k]; ++d){
                List Bx_n(n);
                for(int i = 0; i < n; ++i){
                  IntegerVector x_index_Temp_n = x_index_Subject[i];
                  NumericVector Bx_ni(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    if(IntegerVector::is_na(x_index_Temp_n[j])){
                      Bx_ni[j] = NA_REAL;
                    }else
                    {
                      Bx_ni[j] = Bx_K_Temp(x_index_Temp_n[j],d);
                    }
                  }
                  Bx_n[i] = Bx_ni;
                }
                Bx_Dk[d] = Bx_n;
              }
            }
            
            List Bxt_Dk(Dk[k]);
            for(int d = 0; d < Dk[k]; ++d){
              List Bx_Temp_Dk = Bx_Dk[d];
              List Bxt_H(H);
              for(int h = 0; h < H; ++h){
                List Bt_Temp_H = Bt_H[h];
                List Bxt_n(n);
                for(int i = 0; i < n; ++i){
                  NumericVector Bx_Temp_n = Bx_Temp_Dk[i];
                  NumericVector Bt_Temp_n = Bt_Temp_H[i];
                  NumericVector Bxt_ni(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    Bxt_ni[j] = Bx_Temp_n[j]*Bt_Temp_n[j];
                  }
                  Bxt_n[i] = Bxt_ni;
                }
                Bxt_H[h] = Bxt_n;
              }
              Bxt_Dk[d] = Bxt_H;
            }
            vimp_Bxt[k] = Bxt_Dk;
            
          }else
          {
            vimp_Bxt[k] = Bxt[k];
          }
        }        
      }

      mu_vimp = mu_zero_vec;
      for(int k = 0; k < K; ++k){
        List Bxt_Temp_K = vimp_Bxt[k];
        List Beta_Temp_K = Beta_Mopt[k];
        for(int d = 0; d < Dk[k]; ++d){
          List Bxt_Temp_Dk = Bxt_Temp_K[d];
          List Beta_Temp_Dk = Beta_Temp_K[d];
          for(int h = 0; h < H; ++h){
            List Bxt_Temp_H = Bxt_Temp_Dk[h];
            NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
            NumericMatrix mu_Hat(N,1);
            for(int i = 0; i < n; ++i){
              IntegerVector id_index_Temp = id_index[i];
              NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
              for(int j = 0; j < ni[i]; ++j){
                mu_Hat( id_index_Temp[j], 0) = Bxt_Temp_n[j]*Beta_Temp_H[l];
              }
            }
            mu_vimp = Matrix_Sum_C_NA(mu_vimp, mu_Hat);
          }
        }
      }
      
      for(int i = 0; i < N; ++i){
        Org_mu_vimp(i,0) = ( mu_vimp(i,0) * y_Std_Error[l] ) + y_Mean[l];
      }
      
      NumericVector Org_y_Temp(N);
      NumericVector Org_mu_Temp(N);
      for(int i = 0; i < N; ++i){
        Org_y_Temp[i] = Org_y(i,l);
        Org_mu_Temp[i] = Org_mu_vimp(i,0);
      }
      Error_Rate_vimp = l2Dist_Vector_C_NA(Org_y_Temp, Org_mu_Temp, id_index)/l2Dist_Vector_C_NA(Org_y_Temp, Vec_zero, id_index);
      vimp_Temp(kk,0) = (Error_Rate_vimp - rmse[l])/rmse[l];
    }
    
    if(Time_Varying){
    for(int kk = 0; kk < p; ++kk){
      for(int hh = 0; hh < H; ++hh){
        List vimp_Bxt(K);
        if(joint){
          for(int k = 0; k < K; ++k){
            List vimp_Bt_H(H);
            if(is_true(any( vimp_set == k ))){
              
              IntegerVector Index_Bt_noise = int_randomShuffle(Index_Bt);
              
              NumericMatrix vimp_Bt(n_unq_tm,H);
              for(int i = 0; i < n_unq_tm; ++i){
                for(int j = 0; j < H; ++j){
                  if(j == hh){
                    vimp_Bt(i,j) = Bt(Index_Bt_noise[i],j);
                  }else
                  {
                    vimp_Bt(i,j) = Bt(Index_Bt[i],j);
                  }
                }
              }
              for(int h = 0; h < H; ++h){
                List Bt_n(n);
                for(int i = 0; i < n; ++i){
                  IntegerVector tm_index_Temp = tm_index[i];
                  NumericVector Bt_ni(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    if(IntegerVector::is_na(tm_index_Temp[j])){
                      Bt_ni[j] = NA_REAL;
                    }else
                    {
                      Bt_ni[j] = vimp_Bt(tm_index_Temp[j], h);
                    }
                  }
                  Bt_n[i] = Bt_ni;
                }
                vimp_Bt_H[h] = Bt_n;
              }
              
              List Bx_Temp_K = Bx_K[k];
              List Bxt_Dk(Dk[k]);
              for(int d = 0; d < Dk[k]; ++d){
                List Bx_Temp_Dk = Bx_Temp_K[d];
                List Bxt_H(H);
                for(int h = 0; h < H; ++h){
                  List Bt_Temp_H = vimp_Bt_H[h];
                  List Bxt_n(n);
                  for(int i = 0; i < n; ++i){
                    NumericVector Bx_Temp_n = Bx_Temp_Dk[i];
                    NumericVector Bt_Temp_n = Bt_Temp_H[i];
                    NumericVector Bxt_ni(ni[i]);
                    for(int j = 0; j < ni[i]; ++j){
                      Bxt_ni[j] = Bx_Temp_n[j]*Bt_Temp_n[j];
                    }
                    Bxt_n[i] = Bxt_ni;
                  }
                  Bxt_H[h] = Bxt_n;
                }
                Bxt_Dk[d] = Bxt_H;
              }
              vimp_Bxt[k] = Bxt_Dk;
            }else
            {
              vimp_Bxt[k] = Bxt[k];
            }
          }
        }
        else {
          for(int k = 0; k < K; ++k){
            List vimp_Bt_H(H);
            if(k == vimp_set[kk]){
              
              IntegerVector Index_Bt_noise = int_randomShuffle(Index_Bt);
              
              NumericMatrix vimp_Bt(n_unq_tm,H);
              for(int i = 0; i < n_unq_tm; ++i){
                for(int j = 0; j < H; ++j){
                  if(j == hh){
                    vimp_Bt(i,j) = Bt(Index_Bt_noise[i],j);
                  }else
                  {
                    vimp_Bt(i,j) = Bt(Index_Bt[i],j);
                  }
                }
              }
              for(int h = 0; h < H; ++h){
                List Bt_n(n);
                for(int i = 0; i < n; ++i){
                  IntegerVector tm_index_Temp = tm_index[i];
                  NumericVector Bt_ni(ni[i]);
                  for(int j = 0; j < ni[i]; ++j){
                    if(IntegerVector::is_na(tm_index_Temp[j])){
                      Bt_ni[j] = NA_REAL;
                    }else
                    {
                      Bt_ni[j] = vimp_Bt(tm_index_Temp[j], h);
                    }
                  }
                  Bt_n[i] = Bt_ni;
                }
                vimp_Bt_H[h] = Bt_n;
              }
              
              List Bx_Temp_K = Bx_K[k];
              List Bxt_Dk(Dk[k]);
              for(int d = 0; d < Dk[k]; ++d){
                List Bx_Temp_Dk = Bx_Temp_K[d];
                List Bxt_H(H);
                for(int h = 0; h < H; ++h){
                  List Bt_Temp_H = vimp_Bt_H[h];
                  List Bxt_n(n);
                  for(int i = 0; i < n; ++i){
                    NumericVector Bx_Temp_n = Bx_Temp_Dk[i];
                    NumericVector Bt_Temp_n = Bt_Temp_H[i];
                    NumericVector Bxt_ni(ni[i]);
                    for(int j = 0; j < ni[i]; ++j){
                      Bxt_ni[j] = Bx_Temp_n[j]*Bt_Temp_n[j];
                    }
                    Bxt_n[i] = Bxt_ni;
                  }
                  Bxt_H[h] = Bxt_n;
                }
                Bxt_Dk[d] = Bxt_H;
              }
              vimp_Bxt[k] = Bxt_Dk;
            }else
            {
              vimp_Bxt[k] = Bxt[k];
            }
          }
        }
        
        mu_vimp = mu_zero_vec;
        for(int k = 0; k < K; ++k){
          List Bxt_Temp_K = vimp_Bxt[k];
          List Beta_Temp_K = Beta_Mopt[k];
          for(int d = 0; d < Dk[k]; ++d){
            List Bxt_Temp_Dk = Bxt_Temp_K[d];
            List Beta_Temp_Dk = Beta_Temp_K[d];
            for(int h = 0; h < H; ++h){
              List Bxt_Temp_H = Bxt_Temp_Dk[h];
              NumericVector Beta_Temp_H = Beta_Temp_Dk[h];
              NumericMatrix mu_Hat(N,1);
              for(int i = 0; i < n; ++i){
                IntegerVector id_index_Temp = id_index[i];
                NumericVector Bxt_Temp_n = Bxt_Temp_H[i];
                for(int j = 0; j < ni[i]; ++j){
                  mu_Hat( id_index_Temp[j], 0) = Bxt_Temp_n[j]*Beta_Temp_H[l];
                }
              }
              mu_vimp = Matrix_Sum_C_NA(mu_vimp, mu_Hat);
            }
          }
        }
        
        for(int i = 0; i < N; ++i){
          Org_mu_vimp(i,0) = ( mu_vimp(i,0) * y_Std_Error[l] ) + y_Mean[l];
        }
        
        NumericVector Org_y_Temp(N);
        NumericVector Org_mu_Temp(N);
        
        for(int i = 0; i < N; ++i){
          Org_y_Temp[i] = Org_y(i,l);
          Org_mu_Temp[i] = Org_mu_vimp(i,0);
        }
        
        Error_Rate_vimp = l2Dist_Vector_C_NA(Org_y_Temp, Org_mu_Temp, id_index)/l2Dist_Vector_C_NA(Org_y_Temp, Vec_zero, id_index);
        vimp_Temp(kk, hh+1 ) = (Error_Rate_vimp - rmse[l])/rmse[l];
      }
    }
  }
    vimp[l] = vimp_Temp;
  }
  
  return List::create(
    _["vimp"] = vimp
  );
}

