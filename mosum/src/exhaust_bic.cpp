#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

//' Get integer vector of changepoint indices,
//' based on bool-field representation of combinations.
//' E.g. for combination index 11 [=1011]: get_comb_ind(c(T,F,T,T))=11
//' @keywords internal
// [[Rcpp::export]]
unsigned get_comb_ind(const std::vector<bool> &active) {
  const unsigned m = active.size();
  unsigned res = 0;
  for (unsigned j=0; j<m; ++j) {
    res += active[j] * (1 << j);
  }
  return res;
}

//' Helping function for algorithm 2: Pre-compute the partial sums
//' S_i = sum{j=k_i+1}^{k_{i+1}}x_i and the partial sums of squared
//' T_i = sum{j=k_i+1}^{k_{i+1}}x_i^2
//' between the (sorted) candidates k_i and k_{i+1} in cand.
//' Output: data frame with 4 columns k_i | k_{i+1} | S_i | T_i
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix extract_sub(const IntegerVector &cand, const NumericVector &x) {
  const unsigned m = cand.length();
  NumericMatrix res(m-1, 4);
  unsigned i=0; // position in cand vector
  int j=cand[i]; // position in x vector
  double sum=0.0;
  double sum_sq=0.0;
  while (i+1<m) {
    sum += x[j];
    sum_sq += x[j]*x[j];
    if (j+1 == cand[i+1]) {
      res(i, 0) = cand[i]+1;
      res(i, 1) = cand[i+1];
      res(i, 2) = sum;
      res(i, 3) = sum_sq;
      sum = 0.0;
      sum_sq = 0.0;
      ++i;
    }
    ++j;
  }
  return res;
}

//' Starting value to iterate (in lexicographical order) 
//' over all bit permutaions having l bits set to 1.
//' E.g.: start_bit_permutations(2) = 3 [=0..011].
//' @keywords internal
// [[Rcpp::export]]
unsigned start_bit_permutations(unsigned l) {
  return (1 << l) - 1;
}

//' Next value to iterate (in lexicographical order) over all bit 
//' permutaions having l bits set to 1.
//' Example sequence (2 bits): {0011, 0101, 0110, 1001, 1010, 1100}.
//' Source: https://stackoverflow.com/questions/1851134/generate-all-binary-strings-of-length-n-with-k-bits-set
//' @keywords internal
// [[Rcpp::export]]
unsigned next_bit_permutation(unsigned v) {
  // Note: The following (commented) alternative is faster on Unix machines:
  // unsigned int t = v | (v - 1);
  // unsigned int w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));
  unsigned int t = (v | (v - 1)) + 1;
  unsigned int w = t | ((((t & -t) / (v & -v)) >> 1) - 1);
  return w;
}

//' Is index i_child a child of index i_parent?
//' ASSERT: i_child is of the form (i_parent XOR i_help),
//'         with i_help having exactly one non-zero bit
//' @keywords internal
// [[Rcpp::export]]
bool is_child(unsigned i_child, unsigned i_parent) {
  return (i_child < i_parent);
}

//' Get number of non-zero bits of a 32bit integer
//' Source: https://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer
//' @keywords internal
// [[Rcpp::export]]
unsigned numberOfSetBits(uint32_t i) {
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//' Does the combination comb of changepoints contain the changepoint k_ind?
//' @keywords internal
// [[Rcpp::export]]
bool comb_contains_cpt(unsigned comb, unsigned k_ind) {
  return comb & (1 << k_ind);
}

//' Compute the Local cost terms of combination icomb (for RSS resp. sBIC).
//' Use pre-computed partial sum matrix sub_sums (see extract_sub) for speedup
//' @keywords internal
// [[Rcpp::export]]
double get_local_costs(unsigned icomb, const NumericMatrix &sub_sums) {
  const unsigned m = sub_sums.nrow() - 1;
  double res = 0.0;
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  for (unsigned j=0; j<=m; ++j) {
    A += sub_sums(j,3);
    B += sub_sums(j,2);
    C += sub_sums(j,1) - sub_sums(j,0) + 1.0;
    if (j==m || comb_contains_cpt(icomb, j)) {
      res += A-B*B/C;
      A = 0.0;
      B = 0.0;
      C = 0.0;
    }
  }
  return res;
}

//' where is leftmost one?
//' https://www.geeksforgeeks.org/find-significant-set-bit-number/
//' @keywords internal
// [[Rcpp::export]]
int setBitNumber(int n) 
{ 
  // Below steps set bits after 
  // MSB (including MSB) 
  
  // Suppose n is 273 (binary 
  // is 100010001). It does following 
  // 100010001 | 010001000 = 110011001 
  n |= n >> 1; 
  
  // This makes sure 4 bits 
  // (From MSB and including MSB) 
  // are set. It does following 
  // 110011001 | 001100110 = 111111111 
  n |= n >> 2; 
  
  n |= n >> 4; 
  n |= n >> 8; 
  n |= n >> 16; 
  
  // Increment n by 1 so that 
  // there is only one set bit 
  // which is just before original 
  // MSB. n now becomes 1000000000 
  n = n + 1; 
  
  // Return original MSB after shifting. 
  // n now becomes 100000000 
  return (n >> 1); 
} 

//' Algorithm II (Local change-point search with SC)
//' 
//' Input cand: =mathcal D, conflicting changepoints candidate set
//' Input sub_sums: Pre-computed partial sums, as obtained by extract_sub
//' Input strength: Exponent for penalty
//' Input log_penalty: log (or polynomial) penalty term?
//' Input n: Overall length of data
//' Input auc: =|mathcal C|, total number of currently active changepoints 
//'       (+candidates)
//' Input min_cost: Minimal RSS with all the candidates
//' 
//' Output sc: (Mx2) matrix (M=2^m with m=|cand|) containing RSS/cost and SC
//'         terms for all combinations within cand. Combinations are indexed
//'         by their implicit integer representation, i.e. sc[0,] corresponds
//'         to the empty set, sc[3,] to {k_1,k_2} [0..011], etc.
//'         Note: Row May be Inf, if combination was not visited in algorithm.
//' Output est_cpts: Integer Vector of estimated changepoints
//' Output final: Bool Vector indicating if combinations are final states
//' Output num_cpts: For debugging purposes
//' @keywords internal
// [[Rcpp::export]]
List exhaust_sc(const IntegerVector &cand,
                const NumericMatrix &sub_sums,
                double strength,
                bool log_penalty,
                unsigned n,
                unsigned auc,
                double min_cost) {
  const unsigned m = cand.length();
  const unsigned M = (1 << m);
  const double n_half = (double)n / 2.0;
  
  const double sc_penalty = (log_penalty ? 
                               std::pow(std::log((double)n), strength) :
                               std::pow((double)n, strength));
  
  const double INF = std::numeric_limits<double>::infinity();
  
  std::vector<bool> flag(M, true);
  NumericVector sc_vals(M, INF);
  NumericVector cost_vals(M, INF); // local costs
  IntegerVector num_cpts(M, NA_INTEGER);
  std::vector<int> final(0);
  int m_star = m;
  
  const double min_cost_local = sum(sub_sums(_,3) -
                                    sub_sums(_,2)*sub_sums(_,2) /
                                      (sub_sums(_,1)-sub_sums(_,0)+1.0));
  // M-1 [=1...1] represents ALL changes
  // 0 [=0...0] represents NO change
  cost_vals[M-1] = min_cost_local;
  sc_vals[M-1] = n_half * log(min_cost / double(n)) + auc*sc_penalty;
  num_cpts[M-1] = m;
  
  // iterate over all combination lengths
  int l = m;
  while (l > 0) {
    
    // iterate over all combinations of length l
    // step 1: pruning: inherit FALSE flags in next generation
    unsigned i_parent = start_bit_permutations(l);
    int count = 0;
    while (i_parent<M) {
      if (!flag[i_parent]) {
        for (unsigned i_child_help=0; i_child_help<m; ++i_child_help) {
          const unsigned i_child = i_parent^(1 << i_child_help);
          if(flag[i_child]){
            if (is_child(i_child, i_parent)) flag[i_child] = false;
          }
        }
      } else{
        count += 1;
        final.push_back(i_parent);
        if(m_star > l) m_star = l;
      }
      i_parent = next_bit_permutation(i_parent);
    }
    if(count == 0) break;
    
    // iterate over all combinations of length l
    // step 2: Compute SCs and update flags
    i_parent = start_bit_permutations(l);
    while (i_parent < M) {
      if(flag[i_parent]){
        for (unsigned i_child_help=0; i_child_help<m; ++i_child_help) {
          const unsigned i_child = i_parent^(1 << i_child_help);
          if (is_child(i_child, i_parent) && flag[i_child]) {
            
            // step 2.1: compute SC
            if (cost_vals[i_child]==INF) {
              cost_vals[i_child] = get_local_costs(i_child, sub_sums);
              const double child_cost = min_cost - min_cost_local +
                cost_vals[i_child];
              num_cpts[i_child] = l-1;
              const double child_auc = auc - m + num_cpts[i_child];
              sc_vals[i_child] = n_half * std::log(child_cost / double(n)) +
                child_auc*sc_penalty;
            }
            
            // step 2.2: pruning
            if (sc_vals[i_parent] < sc_vals[i_child]) flag[i_child] = false;
          }
        }
      }
      i_parent = next_bit_permutation(i_parent);
    }
    l -= 1;
  }
  
  if((cost_vals[0]!=INF) & flag[0]) {
    final.push_back(0);
    m_star = 0;
  }
  
  // index of final combination (minimize sc, if several)
  unsigned final_ind_star = 0;
  double min_sc = INF;
  //unsigned j;
  //for (int i = 0; i < final.size(); i++) {
  // j = final[i];
  //  if (num_cpts[j] >= m_star & num_cpts[j] <= m_star + 2){
  //    if(sc_vals[j] < min_sc){
  //      min_sc = sc_vals[j];
  //      final_ind_star = j;    
  //    }
  //  }
  //}
  unsigned left, right, j, jj;
  for (unsigned i = 0; i < final.size(); i++) {
    j = final[i];
    if ((num_cpts[j] >= m_star) & (num_cpts[j] <= m_star + 2)){
      if(sc_vals[j] < min_sc){
        min_sc = sc_vals[j];
        final_ind_star = j;    
      }
      if(num_cpts[j] >= 1){
        left = setBitNumber(j);
        if(num_cpts[j] >= 2){
          right = j^(j & (j - 1));  
        } else{
          right = left; 
        }
        jj = j - left;
        if (cost_vals[jj]==INF) {
          cost_vals[jj] = get_local_costs(jj, sub_sums);
          const double child_cost = min_cost - min_cost_local + cost_vals[jj];
          num_cpts[jj] = num_cpts[j] - 1;
          const double child_auc = auc - m + num_cpts[jj];
          sc_vals[jj] = n_half * std::log(child_cost / double(n)) + child_auc*sc_penalty;
        }
        if(sc_vals[jj] < min_sc){
          min_sc = sc_vals[jj];
          final_ind_star = jj;    
        }
        if(num_cpts[j] >= 2){
          jj = j - right;
          if (cost_vals[jj]==INF) {
            cost_vals[jj] = get_local_costs(jj, sub_sums);
            const double child_cost = min_cost - min_cost_local + cost_vals[jj];
            num_cpts[jj] = num_cpts[j] - 1;
            const double child_auc = auc - m + num_cpts[jj];
            sc_vals[jj] = n_half * std::log(child_cost / double(n)) + child_auc*sc_penalty;
          }
          if(sc_vals[jj] < min_sc){
            min_sc = sc_vals[jj];
            final_ind_star = jj;    
          }
          jj = j - left - right;
          if (cost_vals[jj]==INF) {
            cost_vals[jj] = get_local_costs(jj, sub_sums);
            const double child_cost = min_cost - min_cost_local + cost_vals[jj];
            num_cpts[jj] = num_cpts[j] - 2;
            const double child_auc = auc - m + num_cpts[jj];
            sc_vals[jj] = n_half * std::log(child_cost / double(n)) + child_auc*sc_penalty;
          }
          if(sc_vals[jj] < min_sc){
            min_sc = sc_vals[jj];
            final_ind_star = jj;    
          }
        }
      }
    }
  }
  
  // get estimated changepoints, accoring to final combination
  IntegerVector est_cpts(numberOfSetBits(final_ind_star));
  unsigned est_cpts_ind = 0;
  for (unsigned j=0; j<m; ++j) {
    if (comb_contains_cpt(final_ind_star, j)) {
      est_cpts[est_cpts_ind] = cand[j];
      ++est_cpts_ind;
    }
  }
  
  NumericMatrix sc(M, 2);
  sc(_,0) = cost_vals;
  sc(_,1) = sc_vals;
  
  List res;
  res["sc"] = sc;
  res["est_cpts"] = est_cpts;
  res["final"] = final_ind_star;
  res["num_cpts"] = num_cpts;
  // remove
  res["finals"] = final;
  

  return res;
}

