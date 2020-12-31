#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

LogicalVector validate_olc(CharacterVector codes);

TEST(olctools_deepstate_test,validate_olc_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  CharacterVector codes  = RcppDeepState_CharacterVector();
  qs::c_qsave(codes,"/home/akhila/fuzzer_packages/fuzzedpackages/olctools/inst/testfiles/validate_olc/inputs/codes.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "codes values: "<< codes << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    validate_olc(codes);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
