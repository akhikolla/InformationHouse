#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

List duration_viterbi(NumericVector aa_sample, NumericVector pipar, NumericMatrix tpmpar, NumericMatrix od, NumericMatrix params);

TEST(signalHsmm_deepstate_test,duration_viterbi_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector aa_sample  = RcppDeepState_NumericVector();
  qs::c_qsave(aa_sample,"/home/akhila/fuzzer_packages/fuzzedpackages/signalHsmm/inst/testfiles/duration_viterbi/inputs/aa_sample.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "aa_sample values: "<< aa_sample << std::endl;
  NumericVector pipar  = RcppDeepState_NumericVector();
  qs::c_qsave(pipar,"/home/akhila/fuzzer_packages/fuzzedpackages/signalHsmm/inst/testfiles/duration_viterbi/inputs/pipar.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "pipar values: "<< pipar << std::endl;
  NumericMatrix tpmpar  = RcppDeepState_NumericMatrix();
  qs::c_qsave(tpmpar,"/home/akhila/fuzzer_packages/fuzzedpackages/signalHsmm/inst/testfiles/duration_viterbi/inputs/tpmpar.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "tpmpar values: "<< tpmpar << std::endl;
  NumericMatrix od  = RcppDeepState_NumericMatrix();
  qs::c_qsave(od,"/home/akhila/fuzzer_packages/fuzzedpackages/signalHsmm/inst/testfiles/duration_viterbi/inputs/od.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "od values: "<< od << std::endl;
  NumericMatrix params  = RcppDeepState_NumericMatrix();
  qs::c_qsave(params,"/home/akhila/fuzzer_packages/fuzzedpackages/signalHsmm/inst/testfiles/duration_viterbi/inputs/params.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "params values: "<< params << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    duration_viterbi(aa_sample,pipar,tpmpar,od,params);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
