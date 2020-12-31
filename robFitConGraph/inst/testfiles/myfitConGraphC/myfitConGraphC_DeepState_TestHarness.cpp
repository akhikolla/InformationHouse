#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

List myfitConGraphC(NumericMatrix amat, NumericMatrix S, double n, double tol);

TEST(robFitConGraph_deepstate_test,myfitConGraphC_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix amat  = RcppDeepState_NumericMatrix();
  qs::c_qsave(amat,"/home/akhila/fuzzer_packages/fuzzedpackages/robFitConGraph/inst/testfiles/myfitConGraphC/inputs/amat.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "amat values: "<< amat << std::endl;
  NumericMatrix S  = RcppDeepState_NumericMatrix();
  qs::c_qsave(S,"/home/akhila/fuzzer_packages/fuzzedpackages/robFitConGraph/inst/testfiles/myfitConGraphC/inputs/S.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "S values: "<< S << std::endl;
  std::ofstream n_stream;
  double n  = RcppDeepState_double();
  n_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/robFitConGraph/inst/testfiles/myfitConGraphC/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::ofstream tol_stream;
  double tol  = RcppDeepState_double();
  tol_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/robFitConGraph/inst/testfiles/myfitConGraphC/inputs/tol");
  tol_stream << tol;
  std::cout << "tol values: "<< tol << std::endl;
  tol_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    myfitConGraphC(amat,S,n,tol);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
