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
// @title NE equilibrium cross case
// @param combisim matrix containing the indices of combinations of simulations
// @param ncross number of times each simulation is used
// @details WARNING It is assumed that the strategies of the last player are sorted in increasing order
// WARNING2: Poffs are supposed to correspond to [Z11, ..., Z1p, Z21, ..., Z2p, ..., Znp] (simulations par objectives)
// [[Rcpp::export]]
LogicalMatrix PSNE_sparseMat_cross(NumericVector NS, NumericMatrix Poffs, IntegerMatrix expindices, IntegerMatrix combisim, int ncross){
  int nplay = NS.length();
  int nalt = Poffs.nrow();
  int nsim = Poffs.ncol()/nplay;
  bool isStrat_i;
  double bestPoffi_j; // Store best payoff for player i and strategies -k
  //IntegerVector isNash(nalt); // 0 if possible Nash equilibrium
  // IntegerVector strat_i(nplay); // store the strategy -i considered
  std::vector<int> ind_old; // Store indices of isNash of current best(s) payoffs

  LogicalMatrix isNash(nalt, combisim.nrow());

  // Initialization
  for(int i = 0; i < nalt; i++){
    for(int j = 0; j < combisim.nrow(); j++){
      isNash(i, j) = true;
    }
  }

  IntegerVector ind_sims(nplay); // Store indices of Poffs considered
  IntegerVector indNash(ncross); // Which columns of isNash should be updated simulatneously

  int tmp = 0;

  IntegerVector start_index(NS(nplay - 1));

  // Compute index of first element for strategy of the last player
  start_index(tmp) = expindices(tmp, nplay - 1);
  for(int i = 1; i < nalt; i++){
    if(expindices(i, nplay - 1) > expindices(i - 1, nplay - 1)){
      tmp++;
      start_index(tmp) = i;
    }
  }


  for(int cb = 0; cb < combisim.nrow(); cb++){

    // Determine indices of Poffs subcase considered
    for(int i = 0; i < nplay; i++){
      ind_sims(i) = i * nsim + combisim(cb, i);
    }

    for(int i = 0; i < nplay; i++){

      // Determine which columns of isNash are turned off simultaneously
      tmp = 0;
      for(int j = 0; j < combisim.nrow(); j++){
        if(combisim(j, i) == combisim(cb, i)){
          indNash(tmp) = j;
          tmp++;
        }
      }

      for(int j = 0; j < (nalt - 1); j++){

        // Continue if not previously screened out
        if(isNash(j, cb)){
          bestPoffi_j = Poffs(j, ind_sims(i)); // Store best payoff for player i and strategies -k

          ind_old.push_back(j);

          // for(int p = 0; p < nplay; p++){
          //   strat_i(p) = expindices(j, p);
          // }

          for(int k = start_index(expindices(j, nplay - 1)); k < nalt; k++){

            if(k == j || (expindices(j, nplay - 1) > expindices(k, nplay - 1) && i < nplay - 1))
              continue;

            // Since the strategies of the last player are sorted, no need to check them all
            if(expindices(j, nplay - 1) < expindices(k, nplay - 1) && i < nplay - 1)
              break;

            isStrat_i = true;
            for(int s = 0; s < nplay; s++){
              if(s != i && expindices(k, s) != expindices(j, s)){
                isStrat_i = false;
                break;
              }
            }

            if(isStrat_i){
              if(Poffs(k, ind_sims(i)) < bestPoffi_j){
                for(std::size_t p = 0, max = ind_old.size(); p != max; p++){
                  for(int q = 0; q < ncross; q++)
                    isNash(ind_old[p], indNash(q)) = false;
                }
                ind_old.clear();
                ind_old.push_back(k);
                bestPoffi_j = Poffs(k, ind_sims(i));
              }else{
                if(Poffs(k, ind_sims(i)) == bestPoffi_j){
                  ind_old.push_back(k);
                }else{
                  for(int q = 0; q < ncross; q++)
                    isNash(k, indNash(q)) = false;
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
NumericMatrix getPoffsCross(LogicalMatrix isNash, NumericMatrix Poffs, IntegerMatrix combisim, int nsim){
  // Count number of Nash equilibriums
  int tmp = 0;
  for(int i = 0; i < isNash.nrow(); i++){
    for(int j = 0; j < isNash.ncol(); j++){
      if(isNash(i,j))
        tmp++;
    }
  }

  NumericMatrix NashPoffs(std::max(tmp, 1), combisim.ncol());

  if(tmp > 0){
    tmp = 0;
    for(int j = 0; j < isNash.ncol(); j++){
      for(int i = 0; i < isNash.nrow(); i++){

        if(isNash(i, j)){
          for(int k = 0; k < combisim.ncol(); k++){
            NashPoffs(tmp, k) = Poffs(i, k * nsim + combisim(j, k));
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
