#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

CharacterVector shorten_olc(CharacterVector olcs, NumericVector lats, NumericVector longs);

TEST(olctools_deepstate_test,shorten_olc_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  CharacterVector olcs  = RcppDeepState_CharacterVector();
  qs::c_qsave(olcs,"/home/akhila/fuzzer_packages/fuzzedpackages/olctools/inst/testfiles/shorten_olc/inputs/olcs.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "olcs values: "<< olcs << std::endl;
  NumericVector lats  = RcppDeepState_NumericVector();
  qs::c_qsave(lats,"/home/akhila/fuzzer_packages/fuzzedpackages/olctools/inst/testfiles/shorten_olc/inputs/lats.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "lats values: "<< lats << std::endl;
  NumericVector longs  = RcppDeepState_NumericVector();
  qs::c_qsave(longs,"/home/akhila/fuzzer_packages/fuzzedpackages/olctools/inst/testfiles/shorten_olc/inputs/longs.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "longs values: "<< longs << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    shorten_olc(olcs,lats,longs);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
