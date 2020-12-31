#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

bool isLabelName(Rcpp::CharacterVector lblToCheck, Rcpp::CharacterVector lbl);

TEST(phylobase_deepstate_test,isLabelName_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  CharacterVector lblToCheck  = RcppDeepState_CharacterVector();
  qs::c_qsave(lblToCheck,"/home/akhila/fuzzer_packages/fuzzedpackages/phylobase/inst/testfiles/isLabelName/inputs/lblToCheck.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "lblToCheck values: "<< lblToCheck << std::endl;
  CharacterVector lbl  = RcppDeepState_CharacterVector();
  qs::c_qsave(lbl,"/home/akhila/fuzzer_packages/fuzzedpackages/phylobase/inst/testfiles/isLabelName/inputs/lbl.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "lbl values: "<< lbl << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    isLabelName(lblToCheck,lbl);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
