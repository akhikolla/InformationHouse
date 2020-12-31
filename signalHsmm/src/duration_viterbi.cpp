#include <Rcpp.h>
using namespace Rcpp;

//' Compute most probable path with extended Viterbi algorithm.
//' 
//' Viterbi algorithm for Hidden Markov Model with duration
//' @param aa_sample \code{character} vector representing single aminoacid sequence.
//' @param pipar probabilities of initial state in Markov Model.
//' @param tpmpar matrix of transition probabilities between states.
//' @param od matrix of response probabilities. Eg. od[1,2] is a probability of signal 2 in state 1.
//' @param params matrix of probability distribution for duration. Eg. params[10,2] is probability of duration of time 10 in state 2.
//' @export
//' @return A list of length four:
//' \itemize{
//'  \item{path}{ a vector of most probable path}
//'  \item{viterbi}{ values of probability in all intermediate points,}
//'  \item{psi}{ matrix that gives for every signal and state the previous state in viterbi path,}
//'  \item{duration}{ matrix that gives for every signal and state gives the duration in that state on viterbi path.}
//'  }
//' @note All computations are on logarithms of probabilities.
// [[Rcpp::export]]
List duration_viterbi(NumericVector aa_sample, NumericVector pipar, NumericMatrix tpmpar, NumericMatrix od, NumericMatrix params) {
  int maxDuration = params.nrow();
  int nstates = pipar.length();
  
//  for(int i=0; i<nstates; i++){
//      Rprintf("%f ", pipar(i));
//  }
//  Rprintf("\n");
//  
//  for(int i=0; i<tpmpar.nrow(); i++){
//    for(int j=0; j<tpmpar.ncol(); j++){
//      Rprintf("%f ", tpmpar(i, j));
//    }
//    Rprintf("\n");
//  }
//  
//  for(int i=0; i<params.nrow(); i++){
//    for(int j=0; j<params.ncol(); j++){
//      Rprintf("%f ", params[i, j]);
//    }
//    Rprintf("\n");
//  }
  
  NumericMatrix viterbi(aa_sample.length(), nstates); //probabiliy values for viterbi path that ends in specific state and signal
  NumericMatrix psi(aa_sample.length(), nstates); //the previous state for viterbi path that ends in specific state and signal
  NumericMatrix dura(aa_sample.length(), nstates); //the current duration for viterbi path that ends in specific state and signal
  
  //first signal is treated seperately
  for(unsigned int j=0; j<nstates; j++) {
    viterbi(0,j) = log(pipar(j)) + log( od(j,aa_sample(0)) );
    psi(0,j) = 0;
    dura(0,j) = 0;
  }  
  
  double previous;
  double transition = 0.0;
  
  for(int i=1; i<aa_sample.length(); i++){
    //For each state we will compute the probability of the viterbi path that ends in that state
    for(int j=0; j<nstates; j++){
      double max = R_NegInf;
      int maxInd = 0; //previous state that is maximising probability
      int maximDura = 0;  //duration that is maximising probability
      for(int k=0; k<nstates; k++){  
        int dMax=std::min(maxDuration-1, i);
        for(int d=0; d<=dMax; d++){
          // for every possible previous state, and for every possible duration in current state
          if (i-d==0){ //if duration is as long as number of signal considered
            if(j == 0){ //only first state is accepted
              previous = 1.0;
              transition = 1.0;
            } else { //other states are discarded with log-probability -Inf
              previous = R_NegInf;
              transition = 1.0;
            }
          } else{ //there is some previous state
            previous = viterbi(i-d-1, k); //previous probability on this viterbi path
            transition = log(tpmpar(k, j)); //probability of transition to current state from prevois state
          }
          double duration = log(params(d, j)); //probability of duration that lasts time d
          double responses = 0; //probability of generating d signal in current state
          for(int l=i-d; l<=i; l++){
            responses += log(od(j, aa_sample[l]));
          }
//          if(j==1){
//            Rprintf("i %d k %d d %d \n", i,k,d);
//            Rprintf("The value of previous is %f\n", previous);
//            Rprintf("The value of transition is %f\n", transition );
//            Rprintf("The value of duration is %f\n", duration );
//            Rprintf("The value of responses is %f\n", responses );
//          }
          if(previous + transition + duration + responses > max){ 
            //if that path is better than the best yet found store it
            max = previous + transition + duration + responses;
            maxInd = k;
            maximDura = d;
          }
        }
        
      }
      //assign information about the best path that ends on signal i and state j
      viterbi(i, j) = max;
      psi(i, j) = maxInd;
      dura(i, j) = maximDura;
    }  
  }
  
  //now we extract information about the path. We look for the most probable path
  IntegerVector path(aa_sample.size());
  //the last state
  int seqLength = aa_sample.length();
  NumericVector prob_path = viterbi(seqLength-1,_);
  
//  for(int i=0; i<prob_path.length(); i++){
//    Rprintf("%f ", prob_path(i));
//  }
//  Rprintf("\n %d \n", which_max(viterbi(seqLength-1,_)));
  
  path(seqLength-1) = which_max(viterbi(seqLength-1,_));
  int i = seqLength-2;
  int last = seqLength-1;
  while(i > 0){
    if(last - i <= dura(last, path[last])){ //we are still in the same state
      path(i) = path(last);
    }
    else{ //time to change the state for the previous, stored in matrix psi
      path(i) = psi(last, path(last));
      last = i;
    }
    i -= 1;
  }
  path(0) = 0;
  
  return List::create(
    _["path"] = path,
    _["viterbi"] = viterbi,
    _["psi"] = psi,
    _["duration"] = dura
    );
}
