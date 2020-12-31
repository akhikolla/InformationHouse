#include <Rcpp.h>
using namespace Rcpp;


/// INDIVIDUAL SCORING FUNCTIONS
// only scorecls_c is actually used in algorithm, others for testing

// define function pointer to switch scoring function
typedef double (*scorefun)(int k, int j, int c, NumericVector seq,
			   int M, NumericMatrix csim);

// s(i, j, C) = -M + sum_{k=i}^j( \sigma(x_k, C)  )
// [[Rcpp::export]]
double icor(int k, int j, int c, NumericVector seq,
	    int M, NumericMatrix csim) {
  
  // sum of similarities of positions i:j to cluster c
  // where k is the running index from k=i to k<=j
  double scr = -M;
  for ( ; k <= j; k++ ) 
    scr += csim( k, c ); 
  return scr;
}

// s(i, j, C) = -M + sum_{k=i}^j( \Delta(C_k, C)  )
// note that the simplest scoring function `ccls' is reflected by
// a similarity matrix `csim' where only the diagonal is 1 and the
// rest is 0
// [[Rcpp::export]]
double ccor(int k, int j, int c, NumericVector seq,
	    int M, NumericMatrix csim) {

  // sum of similarities of clusters at positions i:j to cluster c
  // where k is the running index from k=i to k<=j
  double scr = -M;
  for ( ; k <= j; k++ ) 
    scr += csim( seq( k )-1, c );  // TODO: in segmentClusters, don't decrease seqr
  return scr;
}

// handles the switch between scoring functions, based on the
// string passed to \code{\link{calculateScoring}}.
XPtr<scorefun> getScorefun(std::string fstr) {
    if (fstr == "icor")
        return(XPtr<scorefun>(new scorefun(&icor)));
    else if (fstr == "ccor")
        return(XPtr<scorefun>(new scorefun(&ccor)));
    else
        return XPtr<scorefun>(R_NilValue); // runtime error as NULL no XPtr
}




// TODO: rm dependency on seq, since on seq.length is required here!
// Note that \code{seq} is defined differently here then
// then in the wrapper interfaces and MUST be a sequence of positive
// integers

//' segmenTier's core dynamic programming routine in Rcpp
//' 
//' @details This is \code{\link{segmenTier}}'s core dynamic programming
//' routine. It constructs the total score matrix S(i,c), based on
//' the passed scoring function ("icor" or "ccor"), and length penalty
//' \code{M}. "Nuisance" cluster "0" can have a smaller penalty \code{Mn}
//' to allow for shorter distances between "real" segments.
//'
//' Scoring function "icor" calculates the sum of similarities of
//' data at positions k:i to cluster centers c over all k and i.
//' The similarities are calculated e.g., as a (Pearson) correlation between
//' the data at individual positions and the tested cluster c center.
//'
//' Scoring function "ccor" calculates the sum of similarities
//' between the clusters at positions k:i to cluster c over all k and i.
//'
//' Scoring function "ccls" is a special case of "ccor" and is NOT handled
//' here, but is reflected in the cluster similarity matrix \code{csim}. It
//' is handled and automatically constructed in the R wrapper 
//' \code{\link{segmentClusters}}, and merely counts the 
//' number of clusters in sequence k:i, over all k and i, that are identical
//' to the tested cluster \code{c}, and sub-tracts 
//' a penalty for the count of non-identical clusters.
//' @param seq the cluster sequence (where clusters at positions k:i are
//' considered). Note, that unlike the R wrapper, clustering numbers
//' here are 0-based, where 0 is the nuisance cluster.
//' @param C the list of clusters, including nuisance cluster '0', see 
//' \code{seq}
//' @param score the scoring function to be used, one of "ccor" or "icor",
//' an apt similarity matrix must be supplied via option \code{csim}
//' @param M minimal sequence length; Note, that this is not a strict
//' cut-off but defined as an accumulating penalty that must be
//' "overcome" by good score
//' @param Mn minimal sequence length for nuisance cluster, Mn<M will allow
//' shorter distances between segments
//' @param csim a matrix, providing either the cluster-cluster (scoring 
//' function "ccor") or the position-cluster similarity function
//' (scoring function "icor")
//' @param multi if multiple \code{k} are found which return the same maximal
//' score, should the "max" (shorter segment) or "min" (longer segment) be used?
//' This has little effect on real-life large data sets, since the situation
//' will rarely occur. Default is "max".
//' @return Returns the total score matrix \code{S(i,c)} and the matrix 
//' \code{K(i,c)} which stores the position \code{k} which delivered
//' the maximal score at position \code{i}. This is used in the back-tracing
//' phase.
//' @references Machne, Murray & Stadler (2017)
//'     <doi:10.1038/s41598-017-12401-8>
//' @export
// [[Rcpp::export]]
List calculateScore(NumericVector seq, NumericVector C, 
		    std::string score, NumericMatrix csim, int M, int Mn,
		    String multi="max") {

  // result matrices S(i,c) and K(i,c)
  int N = seq.length();  // TODO: get N and M from SM
  int L = C.length();
  NumericMatrix S(N,L); // total score matrix
  NumericMatrix K(N,L); // for backtracing
  NumericMatrix S1(N,L); // S(k,i,c) dynamically constructed
  
  // initialize matrix to 0 and first seq cluster to 1
  // S(0,-1) = 0  wins over S(0,c) = -Inf; 
  std::fill( S.begin(), S.end(), 0.0 );
  std::fill( K.begin(), K.end(), 1 ) ;

  // scoring function "ccls" is a special case of "ccor"
  // and should be handled externally, i.e., an appropriate
  // csim should have been defined! See function 
  // clusterSegments in segment.R
  if ( score=="ccls" ) 
    score = "ccor";

  // get scoring function
  XPtr<scorefun> xpfun = getScorefun(score);
  scorefun scoref = *xpfun;
  
  // calculate score function s(1, i, C) from 1 to i
  // TODO: can this be done dynamically below?
  //       -> would avoid L loops over N
  for ( int c=0; c<L; c++ )  {
    
    int m;
    if ( c==0 ) m = Mn; // c==0 is nuisance cluster
    else m = M;

    // Rcpp::Rcout << "HALLO: " << m << std::endl;

    //TODO: is nuisance used? it seems that Mn has
    // no effect and the passed cluster sequence starts at 1!?
    for ( int i=0; i<N; i++ ) 
      S1(i,c) = scoref(0, i, c, seq, m, csim); // s(1, i, C)

    // initialization (basis case) 
    // S(1,C) = s(1, 1, C) = −M + ∆(x1 , C) for C = C and S(1,C0) = 0. 
    S(0,c) = S1(0,c); //note: this was 0 in prev. implementation !?
    // S(2,C) = s(1, 2, C); since no max_{D!=C} S(j−1, D) is defined yet
    S(1,c) = S1(1,c); // TODO: can this be solved below within the loop?
  }

  //S(i,C) = max_{j<=i} max_{D!=C} ( S(j−1, D) + s(j, i, C) )
  // go through sequence of clusters
  // start at index 2, since 0&1 were initialized already
  for ( int i=2; i<N; i++ ) {
    for ( int c=0; c<L; c++ ) {
      
      int jmax = i; // j<=i
      NumericVector scr(jmax); // store values from j=0 to j=i-1

      // max_D ( S(j-1,D) + score(j,i,c) )

      scr[0] = S1(i, c); // s(1,i,c)  for j=0
      for ( int j=1; j<jmax; j++ ) {
	
	// s(j,i,C) = -M + s(1,i,C) - s(1,j-1,C)
	scr[j] = -M + S1(i, c) - S1(j-1, c);
	
	// + max_D S(j-1,D) : add this sub-max on the fly 
	double mxsj = - std::numeric_limits<double>::infinity();
	for ( int cp=0; cp<L; cp++ ) 
	  if ( cp!=c ) 
	    if ( S(j-1,cp) > mxsj ) mxsj = S(j-1,cp);
	scr[j] += mxsj;
      }
      // max_j ( max_D ( S(j-1,D) + score(j,i,c) ) )
      float mxsc = max( scr );
      S(i,c) = mxsc; // TODO: why not directly max(scr) but via mxsc?

      // store which j was used for back-tracing
      if ( multi=="max" ) { 
	// "max" means the closest j from the left - shorter
	// i.e. min(i-j)
	NumericVector rcs = clone<NumericVector>(scr);
	std::reverse(rcs.begin(), rcs.end());
	K(i,c) = jmax - which_max( rcs );
      } else {
	K(i,c) = which_max( scr ) + 1;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("S1") = S1, /* scoring fnc s(0,i,c)*/
			    Rcpp::Named("S") = S,   /* scoring mat S(i,c) */
			    Rcpp::Named("K") = K);  /* back-tracing matrix */
}

