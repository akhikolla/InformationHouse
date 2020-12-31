#include <Rcpp.h>
#include <math.h>       /* pow */
using namespace Rcpp;
//' @title C++ Function for Weighted K-Means
//' @name KMeansW
//' @description This function does a weighted K-means clustering.
//' @param nclust The number of clusters.
//' @param start The current cluster membership vector.
//' @param weight The vector of length \code{nrows(data)} with weights with nonnegative elements.
//' @param data The concatenated data, with N rows and M columns. Currently, the columns are clustered.
//' @param eps Numerical absolute convergence criteria for the K-means.
//' @param IterMax Integer giving the maximum number of iterations allowed for the K-means.
//' @param cm Numeric vector of class indicators.
//' @param M Matrix of cluster means.
//' @return A list with the folowing values.
//' \item{centers}{the \code{nclust} by M matrix \code{centers} of cluster means.} 
//' \item{cluster}{vector of length N with cluster memberships.} 
//' \item{loss}{vector of length \code{IterMax} with the first entries containing the loss.} 
//' \item{iterations}{the number of iterations used (corresponding to the number 
//' of nonzero entries in \code{loss})} 
//' @examples 
//' set.seed(1)
//' clustmem <- sample.int(n = 10, size = 100, replace = TRUE)
//' mat <- rbind(matrix(rnorm(30*4, mean = 3), nrow = 30), 
//'              matrix(rnorm(30*4, mean = -2), nrow = 30), 
//'              matrix(rnorm(40*4, mean = 0), nrow = 40))
//' wt <- runif(100)
//' testMeans <- lsbclust:::ComputeMeans(cm = clustmem, data = mat, weight = wt, nclust = 3)
//' testK <- lsbclust:::KMeansW(start = clustmem, data = mat, weight = wt, nclust = 3)
// [[Rcpp::export]]
NumericMatrix ComputeMeans(IntegerVector cm, NumericMatrix data, NumericVector weight, int nclust){
  // Compute the matrix of weighted means given a cluster membership vector cm and 
  // a data matrix data
  int n = data.nrow();
  int m = data.ncol();
  
  NumericMatrix M(nclust, m); // Means
  NumericMatrix sumwdat(nclust, m); // 
  NumericMatrix sumw(nclust, m);
  // Initialize sumw and sumwdat to zero
  // k over clusters, j over variables
  for (int k = 0; k < nclust; k++) {
    for (int j = 0; j < m; j++) {
      sumw(k, j)    = 0;
      sumwdat(k, j) = 0;
    }
  }
  // For each variable j (column j of data ) and cluster k compute the sum of the weights of 
  // individuals in that cluster and the weighted sum of the observed data.
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < nclust; k++) {
      if (cm[i] == k + 1){                // The + 1 is needed as c++ start indexes at 0
      for (int j = 0; j < m; j++) {
        sumw(k, j)    += weight[i];
        sumwdat(k, j) += weight[i]*data(i, j);
      }
      } 
    }
  }
  // For each variable j (column j of data ) and cluster k compute the weighted mean   
  for (int k = 0; k < nclust; k++) {
    for (int j = 0; j < m; j++) {
      if (sumw(k, j) > 0){
        M(k, j) = sumwdat(k, j) / sumw(k, j);       
      } else {
        M(k, j) = 0;
      }
    }
  }
  return M;
}

//' @rdname KMeansW
// [[Rcpp::export]]
List AssignCluster(NumericMatrix data, NumericVector weight, NumericMatrix M, int nclust){
  // Assign objects to nearest cluster mean. 
  int n = data.nrow();
  int m = data.ncol();
  double minLoss;
  int    minLossInd; 
  double maxLoss;
  int    maxLossInd;
  int    EmptyClusterID = -1;
  double sumLoss = 0;
  
  
  IntegerVector cm(n);                 // cluster memberships
  NumericVector loss(n);               // loss per object
  NumericVector lossClust(nclust);     // loss per cluster (needed for cluster assignment)
  NumericVector sumw(nclust);          // sum of weights per clusters (needed too discover empty clusters)
  
  // Initialize sumw to zero
  for (int k = 0; k < nclust; k++) {
    sumw[k] = 0;
  }
  
  // Loop over objects i.   
  for (int i = 0; i < n; i++) {
    // Initialize lossClust to 0
    for (int k = 0; k < nclust; k++) {
      lossClust[k] = 0;
    }
    minLoss = 0;
    // Compute for the squared Euclidean distance to each cluster center 
    for (int k = 0; k < nclust; k++) {
      for (int j = 0; j < m; j++) {
        lossClust[k] += pow(data(i, j) - M(k, j),2); 
      }
      // Check if distance to cluster k is smaller than best cluster so far
      if (k == 0 || lossClust[k] < minLoss){
        minLoss    = weight[i] * lossClust[k];
        minLossInd = k;
      }
    }
    cm[i]   = minLossInd + 1;  // the + 1 is needed for use in R (as c++ indexes arrays from 0)
    loss[i] = weight[i] * minLoss;
    sumw[minLossInd] += weight[i];
    if (i == 0 || loss[i] > maxLoss){
      maxLoss    = loss[i];
      maxLossInd = i;
    }
  }
  // Check for empty clusters
  for (int k = 0; k < nclust; k++) {
    if (sumw[k] == 0){
      EmptyClusterID = k;
    }
  }    
  // 
  while (EmptyClusterID >= 0){
    // Place object maxLossInd into cluster EmptyClusterID
    sumw[EmptyClusterID] = weight[maxLossInd];
    cm[maxLossInd] = EmptyClusterID + 1; 
    loss[maxLossInd] = 0;
    // Check if empty clusters remain
    EmptyClusterID = -1;
    for (int k = 0; k < nclust; k++) {
      if (sumw[k] == 0){
        EmptyClusterID = k;
      }
    }
    // Check if there is an empty cluster left
    if (EmptyClusterID >= 0){
      // Find the loss and index of worst fitting cluster
      for (int i = 0; i < n; i++) {
        if (i == 0 || loss[i] > maxLoss){
          maxLoss    = loss[i];
          maxLossInd = i;
        }
      }
    }
  }
  // Compute overall loss
  for (int i = 0; i < n; i++) {
    sumLoss += loss[i];
  }
  
  return (Rcpp::List::create(Rcpp::Named("cm") = cm, Rcpp::Named("loss") = loss, 
  Rcpp::Named("sumLoss") = sumLoss));
}


// Main function
//' @rdname KMeansW
// [[Rcpp::export]]
List KMeansW(int nclust, IntegerVector start, NumericMatrix data, NumericVector weight, 
              double eps = 1e-8, int IterMax = 100) {
  int n = data.nrow();
  int m = data.ncol();
  double sumLoss = 0;
  double sumLossPrev = 0;
  int iter = 0;
  IntegerVector cm(n);
  NumericMatrix M(nclust, m);
  NumericMatrix sumwdat(nclust, m);
  NumericMatrix sumw(nclust, m);
  NumericVector AllLoss(IterMax);
  List out;  
  
  //Compute (weighted) means for current cluster memberships
  cm = start; // cluster membership vector
  // Do k-means iterations
  while (iter < 2 || sumLossPrev - sumLoss > eps){ // && iter < IterMax
    iter += 1;
    // printf ("Iteration: %i;\t", iter);
    sumLossPrev = sumLoss;
    M  = ComputeMeans(cm, data, weight, nclust);
    out = AssignCluster(data, weight, M, nclust);
    cm = as<IntegerVector>(out["cm"]);
    sumLoss = as<double>(out["sumLoss"]);
    // printf ("Loss: %f\n", sumLoss);
    AllLoss[iter - 1] = sumLoss;
  }  
  return (Rcpp::List::create(Rcpp::Named("centers")= M, Rcpp::Named("cluster") = cm, 
          Rcpp::Named("loss") = AllLoss, Rcpp::Named("iterations") = iter));
  
}
