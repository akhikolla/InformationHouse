
#include <Rcpp.h>
using namespace Rcpp;



// Binary search to find the location of the largest time that is smaller of the pattern
// 
// @param array A vector of sorted numbers to search the pattern in
// @param pattern A number to be searched in the array
// @return The location of the largest time that is smaller of the pattern
// 
// @details This is an inner function inside SG package, not for the users.
// [[Rcpp::export]]
int binary_search_km_ge (NumericVector array,  double pattern) { //ge stands for greater
  //this function goal is to return the location of the largest time that is smaller of the pattern
  int array_length=array.size();
  int i = 0, j = array_length - 1;
  while (i <= j) {
    int k = (i + j) / 2;
    if (array[k] == pattern) {
      return (k-1);
      // Note that K cannot be equal to the earliest time time since the earliest time is zero
    }
    else if (array[k] < pattern) {
      i = k + 1;
    }
    else {
      j = k - 1;
    }
  }
  if (j == -1){
    return 0;
  }
  else{
    return j;
  }
}

// Binary search to find the location of the largest time that is smaller or equal of the pattern
// 
// @param array A vector of sorted numbers to search the pattern in
// @param pattern A number to be searched in the array
// @return The location of the largest time that is smaller or equal of the pattern
// 
// @details This is an inner function inside SG package, not for the users.
// [[Rcpp::export]]
int binary_search_km_g (NumericVector array,  double pattern) { //g stands for greater 
  //this function goal is to return tthe location of the largest time that is smaller or equal to the pattern
  int array_length=array.size();
  int i = 0, j = array_length - 1;
  while (i <= j) {
    int k = (i + j) / 2;
    if (array[k] == pattern) {
      return (k);
      // Note that K cannot be equal to the earliest time time since the earliest time is zero
    }
    else if (array[k] < pattern) {
      i = k + 1;
    }
    else {
      j = k - 1;
    }
  }
  if (j == -1){
    return 0;
  }
  else{
    return j;
  }
}



// A function to return the test statistic of the SG test
// 
// @param s0 A vector of the survival KM estimates in group 0 for all sorted unique times in group 0
// @param s1 A vector of the survival KM estimates in group 1 for all sorted unique times in group 1
// @param time0 all sorted unique times in group 0
// @param time1 all sorted unique times in group 1
// @param time The follow up time for all the data
// @param delta A binary status vector, where 0 stands for censored observations and 1 stands for events (corresponds with time)
// @param trt Group vector, must contain 0 and 1 values only (corresponds with time)
// @param tau The maximum time in which we can estimate the Kaplan Meier in the two groups
// @return A list with the two test stats and table usage
// 
// @details This is an inner function inside SG package, not for the users.
// [[Rcpp::export]]
List hhgsurv_test_stat(NumericVector s0,NumericVector s1, NumericVector time0, NumericVector time1,
                        NumericVector time,IntegerVector delta,IntegerVector trt,double tau){
  //Note that the vectors for km must contain a value for time=0
  int n = time.size();
  int n1 = sum(trt);
  int n0 = n-n1;
  int table_possible = n*(n-1);
  
  double A11;
  double A12;
  double A21;
  double A22;
  double A1_;
  double A_1;
  double A2_;
  double A_2;
  int n_i;
  int n_j;
  
  double time_low;
  double time_up;
  
  double obs_drop;
  
  int sum_tables = 0;
  double chisq_stat = 0;
  double lr_stat = 0;
  
  double temp1;
  double temp2;
  double temp3;
  double temp4;
  
  
  //defining s and time vectors that will be initialized for each different i 
  NumericVector si;
  NumericVector sj;
  NumericVector ti;
  NumericVector tj;
  
  
  //This will contain the 2 KM's values for the two relevant time points 
  double km_i_time_low;
  double km_i_time_up;
  double km_j_time_low;
  double km_j_time_up;
  
  
  IntegerVector ind = seq(0,n-1);
  ind = ind[(delta==1) & (time<=tau)];
  for(IntegerVector::iterator i = ind.begin(); i != ind.end(); ++i) {
    
    if (trt[*i] == 0) {
      ti = time0;
      tj = time1;
      si = s0;
      sj = s1;
      n_i = n0;
      n_j = n1;
      
    }
    else {
      ti = time1;
      tj = time0;
      si = s1;
      sj = s0;
      n_i = n1;
      n_j = n0;
    }
    
    
    for(IntegerVector::iterator j = ind.begin(); j != ind.end(); ++j) {
      if(((2*time[*i]-time[*j])<=tau) && (*j!=*i)){
        sum_tables += 1;
        if (trt[*i]==trt[*j]){
          obs_drop = 1;
        }
        else{
          obs_drop = 0;
        }
        
        if (time[*i] > time[*j]) {
          time_low = time[*j];
          time_up = 2*time[*i]-time[*j];
        }
        else {
          time_low = 2*time[*i]-time[*j];
          time_up = time[*j];
        }
        km_i_time_low = si[binary_search_km_ge(ti,time_low)];
        km_i_time_up = si[binary_search_km_g(ti,time_up)];
        km_j_time_low = sj[binary_search_km_ge(tj,time_low)];
        km_j_time_up = sj[binary_search_km_g(tj,time_up)];
        
        A11 = (km_i_time_low-km_i_time_up)*n_i - 1 -obs_drop;
        A12 = (km_j_time_low-km_j_time_up)*n_j -(1-obs_drop);
        A21 = (1-km_i_time_low+km_i_time_up)*n_i;
        A22 = (1-km_j_time_low+km_j_time_up)*n_j;
        
        A1_ = A11+A12;
        A2_ = A21+A22;
        A_1 = A11+A21;
        A_2 = A12+A22;
        
        if (A1_>0 && A2_>0 && A_1>0 && A_2>0){
          chisq_stat += ((n-2)*pow(A12*A21-A11*A22,2.0))/(A1_*A2_*A_1*A_2);
          if (A11>0){
            temp1 = A11*log((n-2)*A11/(A1_*A_1));
          }
          else{
            temp1 = 0;
          }
          if (A12>0){
            temp2 = A12*log((n-2)*A12/(A1_*A_2));
          }
          else{
            temp2 = 0;
          }
          if (A21>0){
            temp3 = A21*log((n-2)*A21/(A2_*A_1));
          }
          else{
            temp3 = 0;
          }
          if (A22>0){
            temp4 = A22*log((n-2)*A22/(A2_*A_2));
          }
          else{
            temp4 = 0;
          }
          
          lr_stat += 2*(temp1 + temp2 + temp3 + temp4);
        }
        
      }
    }
  }
  List lst; lst["chisq_stat"] = chisq_stat/sum_tables;
  lst["lr_stat"] = lr_stat/sum_tables; lst["tab_usage"] = (double)sum_tables/table_possible;
  return lst;
}




// An implementation of Kaplan Meier calculation in C++
// 
// @param time The follow up time 
// @param delta A binary status vector, where 0 stands for censored observations and 1 stands for events (corresponds with time)
// @return The function returns two vectors inside a list \cr
// 
// \code{time} - returns the sorted unique time of the data \cr 
// \code{s} - returns the corresponding KM survival estimates to time
// 
// @details This is an inner function inside SG package, not for the users.
// [[Rcpp::export]]
List KM_C(NumericVector time,IntegerVector delta){
  //Note that time should have no ties, ties should be solved in r before
  NumericVector time_sort_unique = sort_unique(time);
  int n = time_sort_unique.size();
  int n_event;
  int risk_to_drop;
  //First we sort the time and the censorship vector
  // Now he have finished sorting time and delta by time
  //Now we create the Kaplan Meier estimator
  int n_risk = delta.size();
  NumericVector event_vec; 
  NumericVector s_estimator(n+1);
  s_estimator[0] = 1;
  for (int i = 0; i < n; i++) {
    event_vec = delta[time==time_sort_unique[i]];
    n_event = sum(event_vec);
    s_estimator[i+1] = (s_estimator[i])*(1 - ((double)n_event/n_risk));
    // the next two rows are to replace risk_to_drop=event_vec.size, but for some reason it doesn't work
    event_vec[event_vec==0] = 1;
    risk_to_drop = sum(event_vec);
    n_risk -= risk_to_drop;
  }
  time_sort_unique.insert(time_sort_unique.begin(),0);
  List lst; lst["time"] = time_sort_unique; lst["s"] = s_estimator;
  return(lst);
}




// A function to calculate the permutation test statistics
// 
// @param trt The original data Group vector, must contain 0 and 1 values only
// @param ptrt_mat The permuted treatment matrix for all permutations
// @param time_original The original data follow up time 
// @param delta_orginial The original data binary status vector, where 0 stands for censored observations and 1 stands for events 
// @param imputed_altern_time_vec The time vector for observations that switched groups in permuted sample
// @param imputed_altern_delta_vec The delta vector for observations that switched groups in permuted sample
// @param n_perm The number of permutations
// @return A list with the two permutations test stats vectors and a permutation table usage vector 
// 
// @details This is an inner function inside SG package, not for the users.
// [[Rcpp::export]]
List get_perm_stats(IntegerVector trt,IntegerMatrix ptrt_mat,NumericVector time_original,
                    IntegerVector delta_orginial,NumericVector imputed_altern_time_vec,
                    IntegerVector imputed_altern_delta_vec,int n_perm){
  int n = trt.size();

  
  NumericVector tab_usage_perm(n_perm);
  
  IntegerVector ptrt(n);
  
  //All of these vairables will evantually enter hhgsurv function
  NumericVector s0;
  NumericVector s1;
  NumericVector time0;
  NumericVector time1;
  NumericVector ptime(n);
  IntegerVector pdelta(n);
  
  //These variables are in order to create tau (and only for that)
  NumericVector time_event0;
  NumericVector time_event1;
  //M represents a really large number
  double M = 2*max(time_original) - min(time_original) +5;
  double max_ev0;
  double max_ev1;
  double max_obs0;
  double max_obs1;
  double max0;
  double max1;
  double tau;
  
  //These four variables are in order to hold the values of the test stat
  NumericVector chisq_stat(n_perm);
  NumericVector lr_stat(n_perm);
  List temp0;
  List temp1;
  
  List test_stat_temp;
  
  
  for (int i = 0; i < n_perm; i++) {
    ptrt = ptrt_mat(_,i);
    ptime[ptrt==trt] = time_original[ptrt==trt];
    ptime[ptrt!=trt] = imputed_altern_time_vec[ptrt!=trt];
    pdelta[ptrt==trt] = delta_orginial[ptrt==trt];
    pdelta[ptrt!=trt] = imputed_altern_delta_vec[ptrt!=trt];
    temp0 = KM_C(ptime[ptrt==0],pdelta[ptrt==0]);
    temp1 = KM_C(ptime[ptrt==1],pdelta[ptrt==1]);
    s0 = temp0["s"];
    time0 = temp0["time"];
    s1 = temp1["s"];
    time1 = temp1["time"];
    
    //Creating tau
    time_event0 = ptime[(ptrt==0) & (pdelta==1)];
    time_event1 = ptime[(ptrt==1) & (pdelta==1)];
    max_ev0 = max(time_event0);
    max_ev1 = max(time_event1);
    max_obs0 = max(time0);
    max_obs1 = max (time1);
    
    if (max_ev0 == max_obs0){
      max0 = M;
    }
    else{
      max0 = max_ev0;
    }
    if (max_ev1 == max_obs1){
      max1 = M;
    }
    else{
      max1 = max_ev1;
    }
    
    tau = min(NumericVector::create(max0,max1));
    test_stat_temp = hhgsurv_test_stat(s0,s1,time0,time1,ptime,pdelta,ptrt,tau);
    chisq_stat[i] = test_stat_temp["chisq_stat"];
    lr_stat[i] = test_stat_temp["lr_stat"];
    tab_usage_perm[i] = test_stat_temp["tab_usage"];
    
  }
  List lst; lst["chisq_stat"]=chisq_stat; lst["lr_stat"]=lr_stat; lst["tab_usage_perm"]=tab_usage_perm;
  return(lst);
}



// A function to return the test statistic of the KONP K-sample test
// 
// @param s_group A list with the survival KM estimates for all sorted unique times in each group, each group is an element in the list
// @param time_group  A list with all sorted unique in each group, each group is an element in the list
// @param n_vec A vector with the sample sizes in each group
// @param time The follow up time for all the data
// @param delta A binary status vector, where 0 stands for censored observations and 1 stands for events (corresponds with time)
// @param trt Group vector, contains the values 1,..,K only, which correspond with time
// @param tau_k A vector contains maximum times in which we can estimate the Kaplan Meier in each group
// @param tau The maximum time in which we can estimate the Kaplan Meier in at least two groups
// @return A list with the two test statistic
// 
// @details This is an inner function inside SG package, not for the users.
// [[Rcpp::export]]
List hhgsurv_test_stat_K_sample(List s_group,List time_group, IntegerVector n_vec,
                       NumericVector time,IntegerVector delta,IntegerVector trt,NumericVector tau_k, double tau){
  //Note that the vectors for km must contain a value for time=0
  int n = time.size();
  int K = s_group.size();
  IntegerVector group = seq(1,K);
  IntegerVector ind_group;
  int k;
  
  double A11;
  double A12;
  double A21;
  double A22;
  double A1_;
  double A_1;
  double A2_;
  double A_2;

  
  double time_low;
  double time_up;
  
  double obs_drop;
  
  int sum_tables = 0;
  double chisq_stat = 0;
  double lr_stat = 0;
  
  double temp1;
  double temp2;
  double temp3;
  double temp4;
  
  
  NumericVector time_k;
  NumericVector s_k;
  int n_k;
  int n_current;

  
  //This will contain the 2 KM's values for the two relevant time points 
  double km_k_time_low;
  double km_k_time_up;
  
  
  IntegerVector ind = seq(0,n-1);
  ind = ind[(delta==1) & (time<=tau)];
  for(IntegerVector::iterator i = ind.begin(); i != ind.end(); ++i) {
    
    k = trt[*i];
    time_k = time_group[k-1];
    s_k = s_group[k-1];
    n_k = n_vec[k-1];

    
    
    for(IntegerVector::iterator j = ind.begin(); j != ind.end(); ++j) {
      if(((2*time[*i]-time[*j])<=tau) && (*j!=*i) && ((2*time[*i]-time[*j])) <= tau_k[k-1] && (time[*j] <= tau_k[k-1])){
        sum_tables += 1;
        if (trt[*i]==trt[*j]){
          obs_drop = 1;
        }
        else{
          obs_drop = 0;
        }
        
        if (time[*i] > time[*j]) {
          time_low = time[*j];
          time_up = 2*time[*i]-time[*j];
        }
        else {
          time_low = 2*time[*i]-time[*j];
          time_up = time[*j];
        }
        NumericVector s_current;
        ind_group = group[(group!=k) & (tau_k>=time_up)]; // this are the names of all groups taken into account in the A2 calculation
        double s_minus_k_low_sum = 0; //contains the number of observations not grom group k that are larger or equal to time_low
        double s_minus_k_up_sum = 0; //contains the number of observations not grom group k that are smaller or equal to time_up
        double n_minus_k = 0;
        for(IntegerVector::iterator g = ind_group.begin(); g != ind_group.end(); ++g){
           s_current = s_group[*g-1];
          s_minus_k_low_sum += n_vec[*g-1]*s_current[binary_search_km_ge(time_group[*g-1],time_low)];
          s_minus_k_up_sum += n_vec[*g-1]*s_current[binary_search_km_g(time_group[*g-1],time_up)];
          n_minus_k +=  n_vec[*g-1];
        }
        
        km_k_time_low = s_k[binary_search_km_ge(time_k,time_low)]; // This is done as in the two sample case
        km_k_time_up = s_k[binary_search_km_g(time_k,time_up)];   // This is done as in the two sample case
        
        A11 = (km_k_time_low-km_k_time_up)*n_k - 1 -obs_drop;
        A12 = s_minus_k_low_sum-s_minus_k_up_sum -(1-obs_drop);
        A21 = (1-km_k_time_low+km_k_time_up)*n_k;
        A22 = n_minus_k - s_minus_k_low_sum + s_minus_k_up_sum;

        
        A1_ = A11+A12;
        A2_ = A21+A22;
        A_1 = A11+A21;
        A_2 = A12+A22;
        
        n_current = n_minus_k+n_k;
        
        if (A1_>0 && A2_>0 && A_1>0 && A_2>0){
          chisq_stat += ((n_current-2)*pow(A12*A21-A11*A22,2.0))/(A1_*A2_*A_1*A_2);
          if (A11>0){
            temp1 = A11*log((n_current-2)*A11/(A1_*A_1));
          }
          else{
            temp1 = 0;
          }
          if (A12>0){
            temp2 = A12*log((n_current-2)*A12/(A1_*A_2));
          }
          else{
            temp2 = 0;
          }
          if (A21>0){
            temp3 = A21*log((n_current-2)*A21/(A2_*A_1));
          }
          else{
            temp3 = 0;
          }
          if (A22>0){
            temp4 = A22*log((n_current-2)*A22/(A2_*A_2));
          }
          else{
            temp4 = 0;
          }
          
          lr_stat += 2*(temp1 + temp2 + temp3 + temp4);
        }
        
      }
    }
  }
  List lst; lst["chisq_stat"] = chisq_stat/sum_tables;
  lst["lr_stat"] = lr_stat/sum_tables; lst["sum_tables"] = sum_tables;
  return lst;
}



// A function to calculate the permutation test statistics for the KONP K-sample test
// 
// @param ptrt_mat The permuted treatment matrix for all permutations
// @param imputed_time_matrix A matrix contianing the time each obsevations will recieve for each group it belongs
// @param imputed_delta_matrix A matrix contianing the delta each obsevations will recieve for each group it belongs
// @param n_perm The number of permutations
// @param n_vec A vector with the sample sizes in each group
// @return A list with the two permutations test stats vectors 
// 
// @details This is an inner function inside SG package, not for the users.
// [[Rcpp::export]]
List get_perm_stats_K_sample(IntegerMatrix ptrt_mat,NumericMatrix imputed_time_matrix,
                             IntegerMatrix imputed_delta_matrix,int n_perm, IntegerVector n_vec){
  
  int n = imputed_time_matrix.nrow();
  
  int K = imputed_time_matrix.ncol();
  
  IntegerVector ptrt(n);
  NumericVector ptime(n);
  IntegerVector pdelta(n);

  //These variables are in order to create tau (and only for that)
  NumericVector time_event_k;
  //M represents a really large number
  double M = std::numeric_limits<double>::infinity();
  double max_ev_k;
  double max_obs_k;
  double max_k;
  double tau;
  
  //variables to enter to the test stat function
  List sk(K);
  List time_k(K);
  NumericVector tau_k(K);
  NumericVector tau_k_temp(K);
  
  
  //These four variables are in order to hold the values of the test stat
  NumericVector chisq_stat(n_perm);
  NumericVector lr_stat(n_perm);
  List temp; // This will contain the output to KM-C function

  List test_stat;
  
  
  for (int p = 0; p < n_perm; p++) {
    ptrt = ptrt_mat(_,p);
    for (int i = 0; i < n; i++)
    {
      ptime[i] = imputed_time_matrix(i,ptrt[i]-1);
      pdelta[i] = imputed_delta_matrix(i,ptrt[i]-1);
    }
    for (int k = 0; k < K; k++)
    {
      //we use ptrt==(k+1) and not k since indexing in c++ start from 0
      temp = KM_C(ptime[ptrt==(k+1)],pdelta[ptrt==(k+1)]);
      sk[k] = temp["s"];
      time_k[k] = temp["time"];
      
      //Creating tau
      time_event_k = ptime[(ptrt==(k+1)) & (pdelta==1)];
      max_ev_k = max(time_event_k);
      max_obs_k = max((NumericVector)time_k[k]);

      if (max_ev_k == max_obs_k){
        max_k = M;
      }
      else{
        max_k = max_ev_k;
      }
      tau_k[k] = max_k;
    }
    tau_k_temp = tau_k;
    std::sort(tau_k_temp.begin(), tau_k_temp.end(), std::greater<double>()); // sorting in a decreasing order
    tau = tau_k_temp[1]; // Second largest value of tau_k
    test_stat = hhgsurv_test_stat_K_sample(sk,time_k,n_vec,ptime,pdelta,ptrt,tau_k,tau);

      
    chisq_stat[p] = test_stat["chisq_stat"];
    lr_stat[p] = test_stat["lr_stat"];
  }
  List lst; lst["chisq_stat"]=chisq_stat; lst["lr_stat"]=lr_stat;;
  return(lst);
}

