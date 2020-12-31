#include <Rcpp.h>
using namespace Rcpp;

bool Pdom(const double* ptr_mat, int i, int j, int nobj, int nr);
std::vector<int> kung(int i_begin, int i_end, const double* ptr_mat, int nobj, int nr);

//' Non dominated indices of a matrix
//' @param mat matrix of objective values, of size n x nobj
//' @return indices of non-dominated points
//' @details Use Kung non-domination sorting
//' @noRd
//' @references
//' Kung, H. T., Luccio, F., & Preparata, F. P. (1975). On finding the maxima of a set of vectors. Journal of the ACM (JACM), 22(4), 469-476.
//' @examples
//' d <- 6
//' n <- 1000
//' test <- matrix(runif(d * n), n)
//' library(emoa)
//' test <- test[order(test[,1]),]
//' indPF_ref <- which(!is_dominated(t(test)))
//' indPF <- nonDomInd(test)
//' all(indPF == indPF_ref)
//'
//' library(microbenchmark)
//' microbenchmark(is_dominated(t(test)), nonDomInd(test))
// [[Rcpp::export]]
std::vector<int> nonDomInd_cpp(NumericMatrix mat){
  const double* ptr_mat = (const double*) &mat(0,0); // pointer to the matrix elements

  return(kung(1, mat.nrow(), ptr_mat, mat.ncol(), mat.nrow()));
}

// Kung sorting
//' @param i_begin, i_end first index in the matrix to consider
//' @param ptr_mat pointer to matrix with dimensions:
//' @param nobj,nr number of columns and rows
//' @noRd
std::vector<int> kung(int i_begin, int i_end, const double* ptr_mat, int nobj, int nr){

  std::vector<int> tmp1;
  if(i_begin == i_end){
    tmp1.push_back(i_begin);
    return(tmp1);
  }

  tmp1 = kung(i_begin, (i_begin + i_end)/2, ptr_mat, nobj, nr);
  std::vector<int> tmp2 = kung((i_begin + i_end)/2 + 1, i_end, ptr_mat, nobj, nr);

  std::size_t j, i = 0;
  std::size_t sizetmp1 = tmp1.size();
  std::size_t sizetmp2 = tmp2.size();

  while(i < sizetmp2){
    j = 0;
    while(j < sizetmp1){
      if(!Pdom(ptr_mat, tmp1[j]-1, tmp2[i]-1, nobj, nr)){ // elements starts from 0 in C++
        j++;
      }else break;
    }
    if(j == sizetmp1) tmp1.push_back(tmp2[i]);
    i++;
  }

  return(tmp1);
}

// true if x_i dominates x_j
bool Pdom(const double* ptr_mat, int i, int j, int nobj, int nr){
  for(int k = 0; k < nobj; k++){
    if(*(ptr_mat + i + k * nr) > *(ptr_mat + j + k * nr)) return(false);
  }
  return(true);
}


//' Determines which elements in a set are dominated by reference points
//' @title Non-dominated points with respect to a reference
//' @param points matrix (one point per row) that are compared to a reference \code{ref} (i.e., not between themselves)
//' @param ref matrix (one point per row) of reference (faster if they are already Pareto optimal)
//' @noRd
//' @examples
//' \dontrun{
//' d <- 6
//' n <- 1000
//' n2 <- 1000
//'
//' test <- matrix(runif(d * n), n)
//' ref <- matrix(runif(d * n), n)
//' indPF <- nonDom(ref)
//'
//' system.time(res <- nonDomSet(test, ref[indPF,,drop = FALSE]))
//'
//' res2 <- rep(NA, n2)
//' library(emoa)
//' t0 <- Sys.time()
//' for(i in 1:n2){
//'   res2[i] <- !is_dominated(t(rbind(test[i,, drop = FALSE], ref[indPF,])))[1]
//' }
//' print(Sys.time() - t0)
//'
//' all(res == res2)
//' }
// [[Rcpp::export]]
LogicalVector nonDomSet(NumericMatrix points, NumericMatrix ref){
  LogicalVector res(points.nrow());
  int j, k;

  for(int i = 0; i < points.nrow(); i++){
    j = 0;
    while(j < ref.nrow()){
      for(k = 0; k < points.ncol(); k++){
        if(points(i,k) <= ref(j, k)){
          break;
        }
      }
      if(k == points.ncol()) break;
      j++;
    }
    if(j == ref.nrow()) res(i) = true;
  }

  return(res);

}




