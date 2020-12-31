#include <Rcpp.h>
using namespace Rcpp;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))


//' Fast estimation of allele and genotype frequencies under Hardy-Weinberg equilibrium
//' 
//' Alleles are assumed to be numerated from 1 and up with no missing label. Thus if the largest value in either allele1 or allele2 is K then we assume that there can be at least K possible alleles.
//' Genotypes are sorted such the the smallest allele comes first, i.e., 2x1 -> 1x2, and 2x3 -> 2x3
//' 
//' @param allele1 An integer vector (starting with values 1 upwards) of first alleles
//' @param allele2 An integer vector (starting with values 1 upwards) of second alleles
//' @param min_alleles A minimum number of unique alleles available
//' @return A list with three variables: allele_freq for estimated allele frequencies, genotype_freq for estimated genotype_frequencies (under HWE assumption), obs_genotype is the frequency of the genotypes, available_genotypes is the number of available genotypes used for the estimation, and unique_alleles is the number of unique alleles (matches the length of allele_freq)
//' @author Claus Ekstrom <claus@@rprimer.dk>
//' @examples
//' al1 <- sample(1:5, size=1000, replace=TRUE, prob=c(.4, .2, .2, .1, .1))
//' al2 <- sample(1:5, size=1000, replace=TRUE, prob=c(.4, .2, .2, .1, .1))
//' hwe_frequencies(al1, al2)
//'
//' @export
// [[Rcpp::export]]
List hwe_frequencies(IntegerVector allele1, IntegerVector allele2, int min_alleles=0) {
  int m1 = max(allele1);
  int m2 = max(allele2);
  int nunique = (m1 > m2) ? m1 : m2;
  nunique = (nunique > min_alleles) ? nunique : min_alleles;
  int N = allele1.size();

  NumericVector allelefreq(nunique);
  NumericVector genotypefreq(nunique * (nunique+1) / 2);
  IntegerMatrix obs(nunique, nunique);
  IntegerVector obsgenotype(nunique * (nunique+1) / 2);


  if (allele1.size() != allele2.size()) {
    stop("The lengths of allele1 and allele2 must match");
  }

  // Compute allele frequencies
  for(int i = 0; i < N; i++){
    allelefreq[allele1[i]-1]++;
    allelefreq[allele2[i]-1]++;
    obs(MIN(allele1[i], allele2[i])-1, MAX(allele1[i], allele2[i])-1)++;
  }
  allelefreq = allelefreq / (2.0*N);

  // Now the genotype frequencies from lower.tri form
  // There are K * (K+1) / 2 of those
  int index = 0;
  for(int i = 0; i < nunique; i++){
    for(int j = i; j < nunique; j++){
      genotypefreq[index] = allelefreq[i]*allelefreq[j];
      if (i != j) 
        genotypefreq[index] *= 2;

      obsgenotype[index] = obs(i,j);

      index++;
    }
  }

  Rcpp::List RVAL =  Rcpp::List::create(Rcpp::Named("allele_freq") = allelefreq,
					Rcpp::Named("genotype_freq") = genotypefreq,
                                        Rcpp::Named("obs_genotype") = obsgenotype,
                                        Rcpp::Named("available_genotypes") = N,
					Rcpp::Named("unique_alleles") = nunique);

  //  RVAL.attr("class") = "htest";
  
  return(RVAL);



}

/*

List sim_hwegt_cpp(NumericVector allelefreq, int N=0) {

  // Margin: 0 - N, 1 - rows, 2 - both
  arma::colvec af(allelefreq.begin(), allelefreq.nrow(), allelefreq.ncol(), false);

  int K = af.size();

  arma::colvec allele_frequency(K);
  arma::colvec gt_frequency(K*(K+1)/2);
  arma::colvec obs(K*(K+1)/2);

  // Sanity checks
  if (N<1)
    Rcpp::stop("The number of observations to simulate must be greater than 0");

  // Check sum to one


  arma::uvec sequence = arma::linspace<arma::uvec>(0, K-1, K);
  arma::uvec allele1 = Rcpp::RcppArmadillo::sample(sequence, N, true, allelefreq);
  arma::uvec allele2 = Rcpp::RcppArmadillo::sample(sequence, N, true, allelefreq);


  for (int i=0; i<N; i++) {
    allele_frequency[allele1[i]]++;
    allele_frequency[allele2[i]]++;
    
  }
  
  


  NumericVector stat = NumericVector::create(_["X-squared"] = originaltt) ;

}

*/
