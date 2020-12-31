#include <Rcpp.h>
using namespace Rcpp;

// Global variables used by the recursion
int len;

double * lprob_term;
double * chisq_term;
double * llr_term;

double * pr;
int * ipt;

bool * ipt0;
bool * bounds;

double stat_lprob;
double stat_prob;
double stat_chisq;
double stat_llr;

double newProbs_prob;
double newProbs_chisq;
double newProbs_llr;

// Recursive function called from the main routine:
// The idea here is to start with the initial point ipt
// and add resp. subtract a total of rp resp. rm to/from distinct entries.
// The recursion is over the categories given by n_cat.
// Here some recursive calls have been replaced by for-loops,
// which results in a considerable speedup. This complicates the recursion quite a bit.
// An earlier version of the recursion, which may help to understand the basic idea
// can be found below (commented out).
void recurse(int rp,int rm,int n_cat,double lprob,double chisq,double llr){
  if(n_cat > 2){
    int lenN = n_cat*len;
    int c = ipt[n_cat] + rm;
    int c1 = 0;
    if(c < 0){
      c1 += c;
      c = 0;
    }
    c += lenN;
    for(int newrm = c1; newrm >= rm; newrm--){
      recurse(rp,newrm,n_cat-1,
              lprob + lprob_term[c],
              chisq + chisq_term[c],
              llr + llr_term[c]);
      c++;
    }
    for(int newrp = rp-1; newrp >= 0; newrp--){
      recurse(newrp,rm,n_cat-1,
              lprob + lprob_term[c],
              chisq + chisq_term[c],
              llr + llr_term[c]);
      c++;
    }
  }
  if(n_cat == 2){
    int len2 = 2*len;
    int c = len2 + ipt[2];
    int c1 = len + ipt[1];
    int c0 = ipt[0];
    
    if(rp == 0){
      if(rm == 0){
        double prob = std::exp(lprob + lprob_term[c] + lprob_term[c1] + lprob_term[c0]);
        
        if(prob > stat_prob) newProbs_prob += prob;
        if(chisq + chisq_term[c] + chisq_term[c1] + chisq_term[c0] < stat_chisq) newProbs_chisq += prob;
        if(llr + llr_term[c] + llr_term[c1] + llr_term[c0] < stat_llr) newProbs_llr += prob;
      }
      else{
        c += rm;
        
        for(int newrm = 0;newrm >= rm;newrm--){
          if(c >= len2) recurse(0,newrm,1,lprob + lprob_term[c],chisq + chisq_term[c],llr + llr_term[c]);
          c++;
        }
      }
    }
    else if(rm == 0){
      for(int newrp = rp;newrp >= 0;newrp--){
        recurse(newrp,0,1,lprob + lprob_term[c],chisq + chisq_term[c],llr + llr_term[c]);
        c++;
      }
    }
    else{
      int c_pm = c + rm;
      int c1_pm = c1 + rp;
      int c0_pm = c0;
      
      int c_mp;
      int c1_mp = c1;
      int c0_mp = c0 + rp;
      
      if(c_pm >= len2){
        recurse(rp,0,1,lprob + lprob_term[c_pm],chisq + chisq_term[c_pm],llr + llr_term[c_pm]);
        c_pm++;
        c0_pm--;
        c1_mp--;
      }
      else{
        c_pm -= len2;
        if(c_pm < 0){
          c0_pm += c_pm;
          c1_mp += c_pm;
          c_pm = 0;
        }
        c_pm += len2;
      }
      c_mp = c_pm;
      
      //pm: Add to cat 1, subtract from cat 0
      double prob_pm = std::exp(lprob + lprob_term[c_pm] + lprob_term[c1_pm] + lprob_term[c0_pm]);
      
      while(c_pm < c && c0_pm >= 0){
        if(prob_pm > stat_prob) newProbs_prob += prob_pm;
        if(chisq + chisq_term[c_pm] + chisq_term[c1_pm] + chisq_term[c0_pm] < stat_chisq) newProbs_chisq += prob_pm;
        if(llr + llr_term[c_pm] + llr_term[c1_pm] + llr_term[c0_pm] < stat_llr) newProbs_llr += prob_pm;
        
        c_pm++;
        prob_pm *= ((c0_pm)/pr[0]*pr[2]/(c_pm-len2));
        c0_pm--;
      }
      if(ipt[0] >= -rm){
        while(c_pm < c + rp){
          if(prob_pm > stat_prob) newProbs_prob += prob_pm;
          if(chisq + chisq_term[c_pm] + chisq_term[c1_pm] + chisq_term[c0_pm] < stat_chisq) newProbs_chisq += prob_pm;
          if(llr + llr_term[c_pm] + llr_term[c1_pm] + llr_term[c0_pm] < stat_llr) newProbs_llr += prob_pm;
          
          c_pm++;
          prob_pm *= ((c1_pm-len)/pr[1]*pr[2]/(c_pm-len2));
          c1_pm--;
        }
      }
      
      //mp: Subtract from cat 1, add to cat 0
      double prob_mp = std::exp(lprob + lprob_term[c_mp] + lprob_term[c1_mp] + lprob_term[c0_mp]);
      
      while(c_mp < c && c1_mp >= len){
        if(prob_mp > stat_prob) newProbs_prob += prob_mp;
        if(chisq + chisq_term[c_mp] + chisq_term[c1_mp] + chisq_term[c0_mp] < stat_chisq) newProbs_chisq += prob_mp;
        if(llr + llr_term[c_mp] + llr_term[c1_mp] + llr_term[c0_mp] < stat_llr) newProbs_llr += prob_mp;
        
        c_mp++;
        prob_mp *= ((c1_mp-len)/pr[1]*pr[2]/(c_mp-len2));
        c1_mp--;
      }
      if(ipt[1] >= -rm){
        while(c_mp < c + rp){
          
          if(prob_mp > stat_prob) newProbs_prob += prob_mp;
          if(chisq + chisq_term[c_mp] + chisq_term[c1_mp] + chisq_term[c0_mp] < stat_chisq) newProbs_chisq += prob_mp;
          if(llr + llr_term[c_mp] + llr_term[c1_mp] + llr_term[c0_mp] < stat_llr) newProbs_llr += prob_mp;
          
          c_mp++;
          prob_mp *= ((c0_mp)/pr[0]*pr[2]/(c_mp-len2));
          c0_mp--;
        }
      }
      
      if(ipt[1]+ipt[0] >= -rm) recurse(0,rm,1,lprob + lprob_term[c+rp],chisq + chisq_term[c+rp],llr + llr_term[c+rp]);
      
    }
    
  }
  if(n_cat == 1){
    if(rp == 0){
      int c1 = ipt[1] + len;
      int c0 = ipt[0] + rm;
      if(c0 < 0){
        c1 += c0;
        c0 = 0;
      }
      
      double prob = std::exp(lprob + lprob_term[c1] + lprob_term[c0]);
      
      while(c1 >= len && c0 <= ipt[0]){
        if(prob > stat_prob) newProbs_prob += prob;
        if(chisq + chisq_term[c1] + chisq_term[c0] < stat_chisq) newProbs_chisq += prob;
        if(llr + llr_term[c1] + llr_term[c0] < stat_llr) newProbs_llr += prob;
        
        c0++;
        prob *= ((c1-len)/pr[1]*pr[0]/c0);
        c1--;
      }
    }
    else if(rm == 0){
      int c1 = ipt[1] + len;
      int c0 = ipt[0] + rp;
      
      double prob = std::exp(lprob + lprob_term[c1] + lprob_term[c0]);
      
      while(c0 >= ipt[0]){
        if(prob > stat_prob) newProbs_prob += prob;
        if(chisq + chisq_term[c1] + chisq_term[c0] < stat_chisq) newProbs_chisq += prob;
        if(llr + llr_term[c1] + llr_term[c0] < stat_llr) newProbs_llr += prob;
        
        c1++;
        prob *= ((c0)/pr[0]*pr[1]/(c1-len));
        c0--;
      }
    }
  }
}

// Adjusted recursive function which allows to avoid enumerating
// areas with negligible probability, only replaces the part
// of recurse for n_cat > 2. Do not call with n_cat < 3!
void recurse0(int rp,int rm,int n_cat,double lprob,double chisq,double llr){
  int lenN = n_cat*len;
  int c = lenN;

  for(int newrp = rp; newrp >= 0; newrp--){
    if(bounds[c]){
      if(ipt0[n_cat-1] && n_cat > 3){
        recurse0(newrp,rm,n_cat-1,
                 lprob + lprob_term[c],
                 chisq + chisq_term[c],
                 llr + llr_term[c]);
      }
      else recurse(newrp,rm,n_cat-1,
                    lprob + lprob_term[c],
                    chisq + chisq_term[c],
                    llr + llr_term[c]);
    }
    c++;
  }
}


// Old simpler recursive function:
// void recurse(int rp,int rm,int n_cat,double lprob,double chisq,double llr){
//   if(n_cat > 1 || (n_cat == 1 && (rp == 0 || rm == 0))){
//     for(int newrm = 0; newrm >= rm; newrm--){
//       int c = ipt[n_cat] + rm - newrm;
//       if(c >= 0){
//         recurse(rp,newrm,n_cat-1,
//                 lprob + lprob_term[n_cat*len + c],
//                 chisq + chisq_term[n_cat*len + c],
//                 llr + llr_term[n_cat*len + c]);
//       }
//     }
//     for(int newrp = 0; newrp < rp; newrp++){
//       int c = ipt[n_cat] + rp - newrp;
//       recurse(newrp,rm,n_cat-1,
//               lprob + lprob_term[n_cat*len + c],
//               chisq + chisq_term[n_cat*len + c],
//               llr + llr_term[n_cat*len + c]);
//     }
//   }
//   if(n_cat == 1 && rp != 0 && rm != 0){
//     int c = ipt[n_cat] + rm;
//     if(c >= 0){
//       recurse(rp,0,n_cat-1,
//               lprob + lprob_term[n_cat*len + c],
//               chisq + chisq_term[n_cat*len + c],
//               llr + llr_term[n_cat*len + c]);
//     }
//     c = ipt[n_cat] + rp;
//     recurse(0,rm,n_cat-1,
//             lprob + lprob_term[n_cat*len + c],
//             chisq + chisq_term[n_cat*len + c],
//             llr + llr_term[n_cat*len + c]);
//   }
//   if(n_cat == 0){
//     int c = ipt[n_cat] + rm + rp;
//     if(c >= 0){
//       lprob += lprob_term[n_cat*len + c];
//       chisq += chisq_term[n_cat*len + c];
//       llr += llr_term[n_cat*len + c];
//       double prob = std::exp(lprob);
//       if(lprob > stat_lprob){
//         newProbs_prob += prob;
//       }
//       if(chisq < stat_chisq){
//         newProbs_chisq += prob;
//       }
//       if(2*llr < stat_llr){
//         newProbs_llr += prob;
//       }
//     }
//   }
// }


//' @title C++ Function for Exact Multinomial Test
//' @export
//' 
//' @description C++ function computing exact multinomial p-values. Does not perform any safety checks. Incorrect input may result in unwanted behavior.
//' Use only through \code{\link{multinom.test}} with method = "exact" is recommended.
//' 
//' @param x Vector of non-negative integers - number of times each outcome was observed.
//' @param p A vector of positive numbers - the hypothesized probabilities for each outcome. Need to sum to 1!
//' @param theta Parameter - p-values less than theta will not be determined precisely.
//' 
//' @details The outcomes should be ordered by the hypothesized probabilities from largest to smallest for optimal performance.
//' 
//' @return Returns a vector containing three values which are the p-values computed from the probability mass, Pearson's chi-square and the log-likelihood ratio statistic.
//' Values below the threshold theta are upper bounds only and not exact p-values!
//' 
//' @seealso \code{\link{multinom.test}}
//' 
// [[Rcpp::export]]
NumericVector multinom_test_cpp(IntegerVector x, NumericVector p,double theta = 0.0001){ // assumes x and p to be ordered in decreasing order of p
  int n = sum(x);
  int m = p.length();
  len = n+1;

  // precompute values which are used frequently
  double *lfac = new double[n+1];
  lfac[0] = 0;
  for(int i = 0;i<n;i++){
    lfac[i+1] = lfac[i] + std::log(i+1);
  }
  double *lint = new double[n+1];
  lint[0] = 0;
  for(int i = 1;i<n+1;i++){
    lint[i] = std::log(i);
  }
  double *lp = new double[m];
  for(int j = 0;j<m;j++){
    lp[j] = std::log(p[j]);
  }
  double *expect = new double[m];
  for(int j = 0;j<m;j++){
    expect[j] = p[j]*n;
  }
  double *expect_inv = new double[m];
  for(int j = 0;j<m;j++){
    expect_inv[j] = 1/expect[j];
  }
  double *lexpect = new double[m];
  for(int j = 0;j<m;j++){
    lexpect[j] = std::log(expect[j]);
  }
  
  // precompute all possible summands of the test statistics
  lprob_term = new double[m*len];
  for(int j = 0;j<m;j++){
    for(int i = 0;i<n+1;i++){
      lprob_term[len*j + i] = i*lp[j] - lfac[i];
    }
  }
  chisq_term = new double[m*len];
  for(int j = 0;j<m;j++){
    for(int i = 0;i<n+1;i++){
      chisq_term[len*j + i] = i*expect_inv[j]*i;
    }
  }
  llr_term = new double[m*len];
  for(int j = 0;j<m;j++){
    for(int i = 0;i<n+1;i++){
      llr_term[len*j + i] = 2*i*(lint[i] - lexpect[j]);
    }
  }
  
  // compute Chisq and LLR values for observation
  double obs_chisq = -n;
  for(int j = 0;j<m;j++){
    obs_chisq += chisq_term[len*j + x[j]];
  }
  double obs_llr = 0;
  for(int j = 0;j<m;j++){
    obs_llr += llr_term[len*j + x[j]];
  }
  
  // compute center (point with integer values close to the expectation)
  NumericVector expected = n*p;
  NumericVector init_point = floor(expected);
  NumericVector decimalplaces = expected - init_point;
  double missing = n - sum(init_point);
  for(int i = 0;i < missing;i++){
    init_point[i]++;
  } 
  double init_lprob = lfac[n];
  double obs_lprob = lfac[n];
  for(int i = 0;i < m;i++){
    init_lprob += init_point[i]*std::log(p[i]) - lfac[int(init_point[i])];
    obs_lprob += x[i]*std::log(p[i]) - lfac[x[i]];
  }
  if(obs_lprob >= init_lprob){
    init_point = x;
    init_lprob = obs_lprob;
  }
  
  // write p and initial point to global variables accessible by the recursive function
  pr = new double[m];
  ipt = new int[m];
  for(int j = 0;j<m;j++){
    lp[j] = std::log(p[j]);
    pr[j] = p[j];
    ipt[j] = init_point[j];
  }
  
  // set boundaries for the test statistics
  stat_lprob = obs_lprob + 0.00000001*fabs(obs_lprob);
  stat_prob = std::exp(stat_lprob);
  stat_chisq = obs_chisq - 0.00000001*fabs(obs_chisq);
  stat_llr = obs_llr - 0.00000001*fabs(obs_llr);
  
  // bounds for small expectations to avoid enumerating many points with negligible probability
  ipt0 = new bool[m];
  bounds = new bool[m*len];
  for(int j = 0;j<m;j++){
    ipt0[j] = (ipt[j] == 0);
    if(ipt0[j]){
      double pbinom = 0; // null probability to observe at least i in category j
      double lp_comp = std::log(1-pr[j]); // complementary probability
      for(int i = n;i > -1;i--){
        pbinom += std::exp(lfac[n] - lfac[i] - lfac[n-i] + i*lp[j] + (n-i)*lp_comp);
        bounds[len*j + i] = pbinom > theta * 0.00000001;
      }
    }
  }
  
  // initial values
  double sumProbs_prob = 0;
  double sumProbs_chisq = 0;
  double sumProbs_llr = 0;
  
  double lprob;
  double chisq;
  double llr;
  
  int r = 0;
  
  do{ 
    newProbs_prob = 0;
    newProbs_chisq = 0;
    newProbs_llr = 0;
    
    lprob = lfac[n];
    chisq = -n;
    llr = 0;
    
    if(m == 2){ // binomial case
      int c1 = ipt[1] + len - r;
      int c0 = ipt[0] + r;
      if(c0 < len){
        double prob = std::exp(lprob + lprob_term[c1] + lprob_term[c0]);
        if(prob > stat_prob) newProbs_prob += prob;
        if(chisq + chisq_term[c1] + chisq_term[c0] < stat_chisq) newProbs_chisq += prob;
        if(llr + llr_term[c1] + llr_term[c0] < stat_llr) newProbs_llr += prob;
      }
      c1 += 2*r;
      c0 -= 2*r;
      if(c0 >= 0 && r > 0){
        double prob = std::exp(lprob + lprob_term[c1] + lprob_term[c0]);
        if(prob > stat_prob) newProbs_prob += prob;
        if(chisq + chisq_term[c1] + chisq_term[c0] < stat_chisq) newProbs_chisq += prob;
        if(llr + llr_term[c1] + llr_term[c0] < stat_llr) newProbs_llr += prob;
      }
    }
    else if(ipt0[m-1] && m > 3) recurse0(r,-r,m-1,lprob,chisq,llr);
    else recurse(r,-r,m-1,lprob,chisq,llr);
    
    sumProbs_prob += newProbs_prob;
    sumProbs_chisq += newProbs_chisq;
    sumProbs_llr += newProbs_llr;
    
    try{
      Rcpp::checkUserInterrupt();
    }
    catch(Rcpp::internal::InterruptedException e) 
    {
      delete[] lprob_term;
      delete[] chisq_term;
      delete[] llr_term;
      delete[] ipt;
      delete[] pr;
      delete[] ipt0;
      delete[] bounds;
      delete[] lfac;
      delete[] lint;
      delete[] lp;
      delete[] expect;
      delete[] expect_inv;
      delete[] lexpect;
      Rcpp::stop("User interrupted or timed out.");
    }
    
    r++;
  } while ((newProbs_prob > 0 || newProbs_chisq > 0 || newProbs_llr > 0 || r < 2)
             && (sumProbs_prob <= 1-theta || sumProbs_chisq <= 1-theta || sumProbs_llr <= 1-theta)  && r < n);
  
  delete[] lprob_term;
  delete[] chisq_term;
  delete[] llr_term;
  delete[] ipt;
  delete[] pr;
  delete[] ipt0;
  delete[] bounds;
  delete[] lfac;
  delete[] lint;
  delete[] lp;
  delete[] expect;
  delete[] expect_inv;
  delete[] lexpect;
  
  // return {1-sumProbs_prob,1-sumProbs_chisq,1-sumProbs_llr};
  return NumericVector::create(1-sumProbs_prob,1-sumProbs_chisq,1-sumProbs_llr);
}


