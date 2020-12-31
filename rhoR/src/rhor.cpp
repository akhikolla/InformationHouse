#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// NumericMatrix toNumericMatrix(arma::imat x) {
//   int nRows=x.n_rows;
//   
//   NumericMatrix y(nRows, x.n_cols);
//   for (int j=0; j<x.n_rows; j++) {
//     for (int i=0; i<x.n_cols; i++) {
//       y(j,i) = x(j, i);
//     }
//   }
//   
//   return y;
// }

// // @name sample_weighted
// // @title sample_weighted
// // @param x vector
// // @param n integer
// // @export
// //' @export
// // [[Rcpp::export]]
// arma::vec sample_weighted(arma::vec xx, arma::vec probs, int n) {
//   // NumericVector ret(n);
//   // NumericVector xx(x);
//   arma::vec ret(n);
//   
//   // IntegerVector rng = seq(0, xx.size()-1);
//   arma::ivec rng = regspace<ivec>(0, xx.size() - 1);
//   for(int w = 0; w < n; w++) {
//     // int tot = sum(xx);
//     // NumericVector probs = xx / tot;
//     // arma::vec probs2 = toArmaVec(probs);
//   //   IntegerVector s = Rcpp::sugar::sample(rng, 1, false,  probs);
//     arma::ivec s = Rcpp::RcppArmadillo::sample(rng, 1, false, probs);
//     xx[s[0]] -= 1;
//     ret[w] = s[0];
//   }
//   // return(xx);
//   return(ret + 1);
// }

// // ' @export
// // [[Rcpp::export]]
// arma::ivec sample_weighted2(arma::imat xx, int n, bool forR = true) {
//   arma::ivec ret = arma::ivec(n);
//   NumericMatrix x = toNumericMatrix(xx);
//   IntegerVector rng = seq(0, x.size() - 1);
//   
//   for(int w = 0; w < n; w++) {
//     double tot = accu(xx);
//     NumericVector probs = x / tot;
//     IntegerVector s = sample(rng, 1, false,  probs);
//     x[s[0]] -= 1;
//     ret[w] = s[0];
//   }
//   
//   if(forR == true) {
//     return(ret + 1);
//   } else {
//     return(ret);
//   }
// }

//' @title sample_contingency_table
//' @name sample_contingency_table
//' 
//' @param xx contingency table matrix
//' @param n int size of the contingency table
//' @param forR bool if true, add 1 to the results accounting for R indices starting at 1
//' 
//' @export
// [[Rcpp::export]]
arma::ivec sample_contingency_table(arma::imat xx, int n, bool forR = true) {
  arma::ivec ret = arma::ivec(n);
  
  int maxInt;
  int s;
  for(int w = 0; w < n; w++) {
    maxInt = accu(xx);
    s = randi<int>(distr_param(0, maxInt));
    ret[w] = (s <= xx(0) ? 0 : (s <= (xx(0)+xx(1))) ? 1 : (s <= (xx(0)+xx(1)+xx(2))) ? 2 : 3);
    xx(ret[w])--;
  }
  
  if(forR == true) return(ret + 1);
  else return(ret);
}

//' @name getBootPvalue_c
//' @title getBootPvalue_c
//' 
//' @param distribution vector of calculated kappas
//' @param result double calculated kappa to compare against
//' 
//' @description 
//' returns the percentage of the time that the distribution was greater or equal to the observed kappa
//' if the result is less than the mean of the distribution, than the p value is 1
//' else return the number of times that the distribution is greater than the result as a percentage of the total number 
//' of items in the distribution
//'
//' @return double calculated p-value
//'
//' @export
// [[Rcpp::export]]
double getBootPvalue_c(arma::vec distribution, double result) {
  if(result < mean(distribution)) {
    return(1.0);
  } 

  else {
    arma::uvec matched = find(distribution >= result);
    
    double rho = (matched.size() * 1.0) / distribution.size();
    return( rho );
  }
}

// Validate a Baserate/Kappa combo against the supplied Precision. 
// Return true if Precision is greater than the calculated right value
// [[Rcpp::export]]
bool check_BRK_combo(double BR, double P, double K) {
  double right = ((2 * BR * K) - (2 * BR) - K) / (K - 2);
  
  return(P > right);
}

//' @name recall
//' @title recall
//' 
//' @param kappa double
//' @param BR double
//' @param P double
//' 
//' @return Recall calculated from provided kappa, BR, and P
//' 
//' @export
// [[Rcpp::export]]
double recall(double kappa, double BR, double P) {
  double top = kappa * P;
  double R = top / (2*P-2*BR-kappa+2*BR*kappa);
  return(R);
}

// Find a valid precision and kappa combo from a distribution of Kappas and Precisions given a baserate
//
// [[Rcpp::export]]
NumericVector find_valid_pk(
  arma::colvec kappaDistribution, arma::vec kappaProbability,
  arma::colvec precisionDistribution, arma::vec precisionProbability,
  double baserate
) {
  arma::ivec kappaRng = regspace<arma::ivec>(0, kappaDistribution.size() - 1);

  arma::ivec whKappa;
  if (kappaProbability.n_elem > 0) {
    whKappa = Rcpp::RcppArmadillo::sample(kappaRng, 1, false, kappaProbability);
  } else {
    whKappa = Rcpp::RcppArmadillo::sample(kappaRng, 1, false);
  }
  double currKappa = kappaDistribution[whKappa.at(0)];
  
  arma::ivec precRng = regspace<arma::ivec>(0, precisionDistribution.size() - 1);
  arma::ivec whPrec = Rcpp::RcppArmadillo::sample(precRng, 1, false);
  double currPrec = precisionDistribution[whPrec.at(0)];
  
  if (!check_BRK_combo(baserate, currPrec, currKappa)) {
    double precisionMin = (2 * baserate * currKappa - 2 * baserate - currKappa) / (currKappa - 2);
    arma::uvec indicies = find(precisionDistribution > precisionMin);

    if (indicies.size() == 0) {
      return(find_valid_pk(kappaDistribution, kappaProbability, precisionDistribution, precisionProbability, baserate));
    }

    precisionDistribution = precisionDistribution.elem(indicies);
    
    if (precisionProbability.size() > 0) {
      precisionProbability = precisionProbability.elem(indicies);
    }
    
    precRng = regspace<ivec>(0, precisionDistribution.size() - 1);
    whPrec = Rcpp::RcppArmadillo::sample(precRng, 1, false);
    currPrec = precisionDistribution[whPrec.at(0)];
  }
  
  return(Rcpp::NumericVector::create(currPrec, currKappa));
}


//' @name generateKPs_c
//' @title generate_kp_list
//' 
//' @param numNeeded int
//' @param baserate double
//' @param kappaMin double
//' @param kappaMax double
//' @param precisionMin double
//' @param precisionMax double
//' @param distributionType int 0 - normal (default), 1 - bell
//' @param distributionLength long
//' 
//' @return matrix of kappa and precision values (column 1 as precision)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix generate_kp_list(
    int numNeeded, double baserate, 
    double kappaMin, double kappaMax, 
    double precisionMin, double precisionMax, 
    int distributionType = 0, long distributionLength = 10000
) {
  double kappaStep = ((kappaMax - kappaMin) / (distributionLength - 1));

  colvec kappaDistribution = regspace(kappaMin, kappaStep, kappaMax);

  arma::vec kappaProbability;
  if(distributionType == 1) {
    kappaProbability = normpdf(kappaDistribution, 0.9, 0.1);
  }

  double precStep = (precisionMax - precisionMin) / (10000 - 1);
  colvec precisionDistribution = regspace(precisionMin, precStep, precisionMax);

  arma::vec precisionProbability;
  Rcpp::NumericMatrix KPs(numNeeded,2);

  for(int i=0; i < numNeeded; i++) {
    KPs(i, _) = find_valid_pk(
      kappaDistribution,
      kappaProbability,
      precisionDistribution,
      precisionProbability, 
      baserate
    );
  }

  return (KPs);
}

//' @name contingency_table
//' @title contingency_table
//' @description Create a contingency table using the provied precision, recall, baserate, and length.
//' 
//' @param precision double
//' @param rec double
//' @param length int
//' @param baserate double
//' 
//' @export
// [[Rcpp::export]]
arma::imat contingency_table(double precision, double rec, int length, double baserate) {
  int gold1s = max(NumericVector::create(round(baserate * length), 1));
  int gold0s = length - gold1s;
  int TP = max(NumericVector::create(round(gold1s * rec), 1));
  int FP = min(NumericVector::create(gold0s, max(NumericVector::create(round(TP * (1 - precision) / precision),  1))));
  
  arma::imat ct = arma::imat(2,2);
  ct(0,0) = TP;
  ct(0,1) = gold1s - TP;
  ct(1,0) = FP;
  ct(1,1) = gold0s - FP;
  
  return(ct);
}

//' @name random_contingency_table
//' @title random_contingency_table
//' @param setLength [TBD]
//' @param baserate [TBD]
//' @param kappaMin [TBD]
//' @param kappaMax [TBD]
//' @param minPrecision [TBD]
//' @param maxPrecision [TBD]
//' 
//' @export
// [[Rcpp::export]]
arma::imat random_contingency_table(
    int setLength,double baserate, 
    double kappaMin, double kappaMax, 
    double minPrecision = 0, double maxPrecision = 1
) {
  NumericMatrix KP = generate_kp_list(1, baserate, kappaMin, kappaMax, minPrecision, maxPrecision);
  double kappa = KP(0,1);
  double precision = KP(0,0);
  double rec = recall(kappa, baserate, precision);
  
  arma::imat ct = contingency_table(precision, rec, setLength, baserate);
  
  return(ct);
}

//' @title kappa_ct
//' @description Calculate kappa from a contingency table
//' @param ct [TBD]
//' 
//' @export
// [[Rcpp::export]]
double kappa_ct(arma::imat ct) {
  double a = ct(0,0); //gold 1 & silver 1
  double c = ct(0,1); //gold 1 & silver 0
  double b = ct(1,0); //gold 0 & silver 1
  double d = ct(1,1); //gold 0 & silver 0
  double size = accu(ct);
  
  double pZero = (a+d) / size;
  double pYes = ((a + b) / size) * ((a + c) / size);
  double pNo =  ((c + d) / size) * ((b + d) / size);
  double pE = (pYes + pNo);
  double k = (pZero - pE) / (1 - pE);
  
  return(k);
}

// [[Rcpp::export]]
arma::imat getHand_ct(arma::imat ct, int handSetLength, double handSetBaserate) {
  int positives = ceil(handSetLength * handSetBaserate);
  arma::ivec positiveIndices; 
  arma::imat otherCT(ct);
  
  if(positives > 0) {
    arma::irowvec gold1s = ct.row(0);
    positiveIndices = sample_contingency_table(gold1s, positives, false);
  
    double sumOnes = sum(positiveIndices == 0);
    double sumTwos = sum(positiveIndices == 1);
    positiveIndices *= 2;
    otherCT(0,0) = otherCT(0,0) - sumOnes;
    otherCT(0,1) = otherCT(0,1) - sumTwos;
  }
  
  arma::ivec otherAsVector = otherCT.as_col();
  arma::ivec otherIndices = sample_contingency_table(otherCT, handSetLength - positives, false);
  arma::ivec allIndices = join_cols<ivec>(positiveIndices, otherIndices);
  arma::imat newCT = arma::imat(2,2);
  
  newCT(0,0) = sum(allIndices == 0);
  newCT(1,0) = sum(allIndices == 1);
  newCT(0,1) = sum(allIndices == 2);
  newCT(1,1) = sum(allIndices == 3);
  
  return(newCT);
}

//' @name getHand_kappa
//' @title getHand_kappa
//' 
//' @description This function returns kappa calculated from a Handset taken from  a larger Contingency Table
//' 
//' @param ct KPs matrix of kappa (column 1) and precision (column 2) values
//' @param handSetLength The length of the \code{\link[=getTestSet]{testSet}} (ignored unless \emph{data} is an observed kappa value)
//' @param handSetBaserate baserate to inflate the sampled contingency table to
//' 
//' @return Kappa as double
//' 
//' @export
// [[Rcpp::export]]
double getHand_kappa(arma::imat ct, int handSetLength, double handSetBaserate) {
  arma::imat newCT = getHand_ct(ct, handSetLength, handSetBaserate);
  
  double k = kappa_ct(newCT);
  return(k);
}

// This function calculates rho given the parameters of the initial set and 
// the kappa of the observed handset and is called by rhoK()
//
// [[Rcpp::export]]
double calcRho_c(
  double x, 
  double OcSBaserate, 
  int testSetLength, 
  double testSetBaserateInflation = 0,
  int OcSLength = 10000, 
  int replicates = 800, 
  double ScSKappaThreshold = 0.9,
  double ScSKappaMin = 0.40,
  double ScSPrecisionMin = 0.6, 
  double ScSPrecisionMax = 1.0,
  NumericMatrix KPs = NumericMatrix(0)
) {
  if (KPs.size() == 1 && KPs[0] == 0) {
    KPs = generate_kp_list(replicates, OcSBaserate, ScSKappaMin, ScSKappaThreshold, ScSPrecisionMin, ScSPrecisionMax);
  }

  if(KPs.nrow() < replicates) {
    replicates = KPs.nrow();
  }

  arma::vec savedKappas = arma::vec(replicates);
  for (int KProw = 0; KProw < replicates; KProw++) {
    NumericVector KP = KPs(KProw, _);
    double precision = KP[0];
    double kappa = KP[1];
    double rec = recall(kappa, OcSBaserate, precision);

    arma::imat fullCT = contingency_table(precision, rec, OcSLength, OcSBaserate);
    double kk = getHand_kappa(fullCT, testSetLength, testSetBaserateInflation);
    savedKappas[KProw] = kk;
  };

  return(getBootPvalue_c(savedKappas, x));
}

// // @name rhoCT_c
// // @title rhoCT_ct
// // @param contingencyTable [TBD]
// // @param testSetBaserateInflation [TBD]
// // @param replicates [TBD]
// // @param OcSLength [TBD]
// // @param OcSBaserate [TBD]
// // @param ScSKappaThreshold [TBD]
// // @param ScSKappaMin [TBD]
// // @param ScSPrecisionMin [TBD]
// // @param ScSPrecisionMax [TBD]
// //' @export
// // [[Rcpp::export]]
// arma::vec rhoCT_c(
//   arma::imat contingencyTable,
//   double testSetBaserateInflation = 0.2,
//   int replicates = 800,
//   int OcSLength = 10000,
//   double OcSBaserate = -1.0,
//   double ScSKappaThreshold = 0.65,
//   double ScSKappaMin = 0.40,
//   double ScSPrecisionMin = 0.6,
//   double ScSPrecisionMax = 1
// ){
//   // Observed kappa
//   double kappa = kappa_ct(contingencyTable);
//   double size = accu(contingencyTable);
//   double gold1;
//   
//   if(OcSBaserate == -1.0) {
//     gold1 = accu(contingencyTable.row(0));
//     OcSBaserate = gold1 / size;
//   }
//   
//   double rho = calcRho_c(
//     kappa, OcSBaserate,
//     size, testSetBaserateInflation,
//     OcSLength, replicates,
//     ScSKappaThreshold, ScSKappaMin, ScSPrecisionMin, ScSPrecisionMax
//   );
//   
//   return(NumericVector::create(kappa, rho));
// }


// //' @export
// // [[Rcpp::export]]
// bool run_rho_in_c(
//     double threshold,  double infl,
//     int handsetLength = 20,
//     double handsetBaserate = 0.2,
//     int replicates = 10000
// ) {
//   arma::imat ct = random_contingency_table(10000, infl, 0.4, threshold, 0.2, 0.8);
//   // handset --> testset
//   arma::imat handset = getHand_ct(ct, handsetLength, handsetBaserate);
//   arma::vec res = rhoCT_c(handset, handsetBaserate, replicates, 10000, -1, 0.65, 0.4, 0.6, 1);
//   
//   return(res[0]>threshold);
// }

