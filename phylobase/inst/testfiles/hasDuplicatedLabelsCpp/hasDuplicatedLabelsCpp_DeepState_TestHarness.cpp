#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

bool hasDuplicatedLabelsCpp(Rcpp::CharacterVector label);

TEST(phylobase_deepstate_test,hasDuplicatedLabelsCpp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  CharacterVector label  = RcppDeepState_CharacterVector();
  qs::c_qsave(label,"/home/akhila/fuzzer_packages/fuzzedpackages/phylobase/inst/testfiles/hasDuplicatedLabelsCpp/inputs/label.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "label values: "<< label << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    hasDuplicatedLabelsCpp(label);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
