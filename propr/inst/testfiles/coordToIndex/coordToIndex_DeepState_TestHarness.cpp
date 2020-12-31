#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

IntegerVector coordToIndex(IntegerVector row, IntegerVector col, const int N);

TEST(propr_deepstate_test,coordToIndex_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector row  = RcppDeepState_IntegerVector();
  qs::c_qsave(row,"/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/coordToIndex/inputs/row.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "row values: "<< row << std::endl;
  IntegerVector col  = RcppDeepState_IntegerVector();
  qs::c_qsave(col,"/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/coordToIndex/inputs/col.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "col values: "<< col << std::endl;
  std::ofstream N_stream;
  int N  = RcppDeepState_int();
  N_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/coordToIndex/inputs/N");
  N_stream << N;
  std::cout << "N values: "<< N << std::endl;
  N_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    coordToIndex(row,col,N);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
