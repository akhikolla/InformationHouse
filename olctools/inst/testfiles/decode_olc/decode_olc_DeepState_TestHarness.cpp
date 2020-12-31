#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

DataFrame decode_olc(CharacterVector olcs);

TEST(olctools_deepstate_test,decode_olc_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  CharacterVector olcs  = RcppDeepState_CharacterVector();
  qs::c_qsave(olcs,"/home/akhila/fuzzer_packages/fuzzedpackages/olctools/inst/testfiles/decode_olc/inputs/olcs.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "olcs values: "<< olcs << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    decode_olc(olcs);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
