// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// fixstar2_mag_DeepState_TestHarness_generation.cpp and fixstar2_mag_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::List fixstar2_mag(Rcpp::CharacterVector starname);

TEST(swephR_deepstate_test,fixstar2_mag_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  CharacterVector starname  = RcppDeepState_CharacterVector();
  qs::c_qsave(starname,"/home/akhila/fuzzer_packages/fuzzedpackages/swephR/inst/testfiles/fixstar2_mag/inputs/starname.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "starname values: "<< starname << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    fixstar2_mag(starname);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}