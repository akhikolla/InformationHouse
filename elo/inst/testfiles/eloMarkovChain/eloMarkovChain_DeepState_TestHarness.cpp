// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// eloMarkovChain_DeepState_TestHarness_generation.cpp and eloMarkovChain_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

List eloMarkovChain(NumericMatrix teamA, NumericMatrix teamB, NumericVector winsA, NumericVector weightsA, NumericVector weightsB, NumericVector weights, NumericMatrix k, int nTeams);

TEST(elo_deepstate_test,eloMarkovChain_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix teamA  = RcppDeepState_NumericMatrix();
  qs::c_qsave(teamA,"/home/akhila/fuzzer_packages/fuzzedpackages/elo/inst/testfiles/eloMarkovChain/inputs/teamA.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "teamA values: "<< teamA << std::endl;
  NumericMatrix teamB  = RcppDeepState_NumericMatrix();
  qs::c_qsave(teamB,"/home/akhila/fuzzer_packages/fuzzedpackages/elo/inst/testfiles/eloMarkovChain/inputs/teamB.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "teamB values: "<< teamB << std::endl;
  NumericVector winsA  = RcppDeepState_NumericVector();
  qs::c_qsave(winsA,"/home/akhila/fuzzer_packages/fuzzedpackages/elo/inst/testfiles/eloMarkovChain/inputs/winsA.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "winsA values: "<< winsA << std::endl;
  NumericVector weightsA  = RcppDeepState_NumericVector();
  qs::c_qsave(weightsA,"/home/akhila/fuzzer_packages/fuzzedpackages/elo/inst/testfiles/eloMarkovChain/inputs/weightsA.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "weightsA values: "<< weightsA << std::endl;
  NumericVector weightsB  = RcppDeepState_NumericVector();
  qs::c_qsave(weightsB,"/home/akhila/fuzzer_packages/fuzzedpackages/elo/inst/testfiles/eloMarkovChain/inputs/weightsB.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "weightsB values: "<< weightsB << std::endl;
  NumericVector weights  = RcppDeepState_NumericVector();
  qs::c_qsave(weights,"/home/akhila/fuzzer_packages/fuzzedpackages/elo/inst/testfiles/eloMarkovChain/inputs/weights.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "weights values: "<< weights << std::endl;
  NumericMatrix k  = RcppDeepState_NumericMatrix();
  qs::c_qsave(k,"/home/akhila/fuzzer_packages/fuzzedpackages/elo/inst/testfiles/eloMarkovChain/inputs/k.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "k values: "<< k << std::endl;
  IntegerVector nTeams(1);
  nTeams[0]  = RcppDeepState_int();
  qs::c_qsave(nTeams,"/home/akhila/fuzzer_packages/fuzzedpackages/elo/inst/testfiles/eloMarkovChain/inputs/nTeams.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "nTeams values: "<< nTeams << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    eloMarkovChain(teamA,teamB,winsA,weightsA,weightsB,weights,k,nTeams[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}