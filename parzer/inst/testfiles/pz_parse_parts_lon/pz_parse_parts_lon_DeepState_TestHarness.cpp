#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

DataFrame pz_parse_parts_lon(CharacterVector x);

TEST(parzer_deepstate_test,pz_parse_parts_lon_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  CharacterVector x  = RcppDeepState_CharacterVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/parzer/inst/testfiles/pz_parse_parts_lon/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    pz_parse_parts_lon(x);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
