#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double eho(NumericVector param, NumericVector numbuys, NumericVector numsells);

TEST(pinbasic_deepstate_test,eho_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector param  = RcppDeepState_NumericVector();
  qs::c_qsave(param,"/home/akhila/fuzzer_packages/fuzzedpackages/pinbasic/inst/testfiles/eho/inputs/param.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "param values: "<< param << std::endl;
  NumericVector numbuys  = RcppDeepState_NumericVector();
  qs::c_qsave(numbuys,"/home/akhila/fuzzer_packages/fuzzedpackages/pinbasic/inst/testfiles/eho/inputs/numbuys.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "numbuys values: "<< numbuys << std::endl;
  NumericVector numsells  = RcppDeepState_NumericVector();
  qs::c_qsave(numsells,"/home/akhila/fuzzer_packages/fuzzedpackages/pinbasic/inst/testfiles/eho/inputs/numsells.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "numsells values: "<< numsells << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    eho(param,numbuys,numsells);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
