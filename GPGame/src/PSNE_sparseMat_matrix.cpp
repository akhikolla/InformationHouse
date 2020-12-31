#include <Rcpp.h>
using namespace Rcpp;

// Computes the pure strategy Nash Equilibria in a N-player finite game.
// ---------------------------------------------------------------------
// input:
//
// NS(i) : total number of pure strategies of player (i)
//
// Poff : sparse matrix representation of Payoffs
// Poff(i,_) (i^th row): columns 1:N are payoffs for the N objectives for strategies in columns N+1:2N
//
// output : isNash, vector giving for each line of Poffs if it is a Nash equilibrium (0) or not (1)
//
// ---------------------------------------------------------------------
// A. Habbal Inria 2016-05
// V. Picheny INRA 2016-06
// M. Binois Chicago Booth 2016-07
// @details
// WARNING: Poffs are supposed to correspond to [Z11, ..., Z1p, Z21, ..., Z2p, ..., Znp] (simulations par objectives)
// [[Rcpp::export]]
LogicalMatrix PSNE_sparseMat(NumericVector NS, NumericMatrix Poffs, IntegerMatrix expindices){
  int nplay = NS.length();
  int nalt = Poffs.nrow();
  int nsim = Poffs.ncol()/nplay;
  bool isStrat_i;
  double bestPoffi_j; // Store best payoff for player i and strategies -k
  LogicalMatrix isNash(nalt, nsim); // 0 if possible Nash equilibrium
  //IntegerVector strat_i(nplay); // store the strategy -i considered
  std::vector<int> ind_old; // Store addresses of isNash of current best(s) payoffs

  // Initialization
  for(int i = 0; i < nalt; i++)
    for(int j = 0; j < nsim; j++)
      isNash(i, j) = true;

  for(int s = 0; s < nsim; s++){

    for(int i = 0; i < nplay; i++){
      for(int j = 0; j < (nalt - 1); j++){

        // Continue if not previously screened out
        if(isNash(j, s)){
          bestPoffi_j = Poffs(j, i * nsim + s); // Store best payoff for player i and strategies -k

          ind_old.push_back(j);

          // for(int p = 0; p < nplay; p++){
          //   strat_i(p) = Poffs(j, nplay + p);
          // }

          for(int k = 0; k < nalt; k++){
            if(k == j)
              continue;

            isStrat_i = true;
            for(int p = 0; p < nplay; p++){
              // if(s != i && Poffs(k, nplay + s) != strat_i(s)){
              if(p != i && expindices(k, p) != expindices(j, p)){
                isStrat_i = false;
                break;
              }
            }

            if(isStrat_i){
              if(Poffs(k, i * nsim + s) < bestPoffi_j){
                for(std::size_t p = 0, max = ind_old.size(); p != max; p++){
                  isNash(ind_old[p], s) = false;
                }
                ind_old.clear();
                ind_old.push_back(k);
                bestPoffi_j = Poffs(k, i * nsim + s);
              }else{
                if(Poffs(k, i * nsim + s) == bestPoffi_j){
                  ind_old.push_back(k);
                }else{
                  isNash(k, s) = false;
                }
              }
            }

          }
          ind_old.clear();
        }
      }
    }
  }

  return isNash;
}

// Computes the pure strategy Nash Equilibria in a N-player finite game.
// ---------------------------------------------------------------------
// input:
//
// NS(i) : total number of pure strategies of player (i)
//
// Poff : sparse matrix representation of Payoffs
// Poff(i,_) (i^th row): columns 1:N are payoffs for the N objectives for strategies in columns N+1:2N
//
// output : isNash, vector giving for each line of Poffs if it is a Nash equilibrium (0) or not (1)
//
// ---------------------------------------------------------------------
// A. Habbal Inria 2016-05
// V. Picheny INRA 2016-06
// M. Binois Chicago Booth 2016-07
// @details WARNING It is assumed that the strategies of the last player are sorted in increasing order
// [[Rcpp::export]]
LogicalMatrix PSNE_sparseMat_sorted(NumericVector NS, NumericMatrix Poffs, IntegerMatrix expindices){
  int nplay = NS.length();
  int nalt = Poffs.nrow();
  int nsim = Poffs.ncol()/nplay;
  bool isStrat_i;
  double bestPoffi_j; // Store best payoff for player i and strategies -k
  LogicalMatrix isNash(nalt, nsim); // true if possible Nash equilibrium
  // IntegerVector strat_i(nplay); // store the strategy -i considered
  std::vector<int> ind_old; // Store indices of isNash of current best(s) payoffs

  // Initialization
  for(int i = 0; i < nalt; i++)
    for(int j = 0; j < nsim; j++)
      isNash(i, j) = true;

  IntegerVector start_index(NS(nplay - 1));
  int tmp = 0;
  // Compute index of first element for strategy of the last player
  start_index(tmp) = expindices(tmp, nplay - 1);
  for(int i = 1; i < nalt; i++){
    if(expindices(i, nplay - 1) > expindices(i - 1, nplay - 1)){
      tmp++;
      start_index(tmp) = i;
    }
  }

  for(int s = 0; s < nsim; s++){
    for(int i = 0; i < nplay; i++){
      for(int j = 0; j < (nalt - 1); j++){

        // Continue if not previously screened out
        if(isNash(j, s)){
          bestPoffi_j = Poffs(j, i * nsim + s); // Store best payoff for player i and strategies -k

          ind_old.push_back(j);

          // for(int p = 0; p < nplay; p++){
          //   strat_i(p) = expindices(j, p);
          // }

          if(i < (nplay - 1)){
            tmp = start_index(expindices(j, nplay - 1));
          }else{
            tmp = 0;
          }

          for(int k = tmp; k < nalt; k++){

            if(k == j)
              continue;

            // Since the strategies of the last player are sorted, no need to check them all
            if(expindices(j, nplay - 1) < expindices(k, nplay - 1) && i < nplay - 1)
              break;

            isStrat_i = true;
            for(int p = 0; p < nplay; p++){
              if(p != i && expindices(k, p) != expindices(j, p)){
                isStrat_i = false;
                break;
              }
            }

            if(isStrat_i){
              if(Poffs(k, i * nsim + s) < bestPoffi_j){
                for(std::size_t p = 0, max = ind_old.size(); p != max; p++){
                  isNash(ind_old[p], s) = false;
                }
                ind_old.clear();
                ind_old.push_back(k);
                bestPoffi_j = Poffs(k, i * nsim + s);
              }else{
                if(Poffs(k, i * nsim + s) == bestPoffi_j){
                  ind_old.push_back(k);
                }else{
                  isNash(k, s) = false;
                }
              }
            }

          }
          ind_old.clear();
        }
      }
    }
  }

  return isNash;
}


// Get Nash Poffs in the crossed case
// [[Rcpp::export]]
NumericMatrix getPoffs(LogicalMatrix isNash, NumericMatrix Poffs, int nsim, int nobj){
  // Count number of Nash equilibriums
  int tmp = 0;
  for(int i = 0; i < isNash.nrow(); i++){
    for(int j = 0; j < isNash.ncol(); j++){
      if(isNash(i,j))
        tmp++;
    }
  }

  NumericMatrix NashPoffs(std::max(tmp, 1), nobj);

  if(tmp > 0){
    tmp = 0;
    for(int j = 0; j < isNash.ncol(); j++){
      for(int i = 0; i < isNash.nrow(); i++){

        if(isNash(i, j)){
          for(int k = 0; k < nobj; k++){
            NashPoffs(tmp, k) = Poffs(i, k * nsim + j);
          }
          tmp++;
        }

      }
    }
  }else{
    for(int i = 0; i < isNash.ncol(); i++){
      NashPoffs(0, i) = NA_REAL;
    }
  }


  return(NashPoffs);
}




