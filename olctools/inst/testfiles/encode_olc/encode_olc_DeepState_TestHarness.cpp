#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

CharacterVector encode_olc(NumericVector lats, NumericVector longs, IntegerVector length);

TEST(olctools_deepstate_test,encode_olc_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector lats  = RcppDeepState_NumericVector();
  qs::c_qsave(lats,"/home/akhila/fuzzer_packages/fuzzedpackages/olctools/inst/testfiles/encode_olc/inputs/lats.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "lats values: "<< lats << std::endl;
  NumericVector longs  = RcppDeepState_NumericVector();
  qs::c_qsave(longs,"/home/akhila/fuzzer_packages/fuzzedpackages/olctools/inst/testfiles/encode_olc/inputs/longs.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "longs values: "<< longs << std::endl;
  IntegerVector length  = RcppDeepState_IntegerVector();
  qs::c_qsave(length,"/home/akhila/fuzzer_packages/fuzzedpackages/olctools/inst/testfiles/encode_olc/inputs/length.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "length values: "<< length << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    encode_olc(lats,longs,length);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
