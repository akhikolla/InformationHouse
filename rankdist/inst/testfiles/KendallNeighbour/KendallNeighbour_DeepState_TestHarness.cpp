#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix KendallNeighbour(NumericVector rank);

TEST(rankdist_deepstate_test,KendallNeighbour_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector rank  = RcppDeepState_NumericVector();
  qs::c_qsave(rank,"/home/akhila/fuzzer_packages/fuzzedpackages/rankdist/inst/testfiles/KendallNeighbour/inputs/rank.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "rank values: "<< rank << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    KendallNeighbour(rank);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
