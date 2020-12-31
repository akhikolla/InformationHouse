// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// do_par_in_hash_int_DeepState_TestHarness_generation.cpp and do_par_in_hash_int_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

LogicalVector do_par_in_hash_int(IntegerVector x, IntegerVector table, int nThread);

TEST(hutilscpp_deepstate_test,do_par_in_hash_int_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector x  = RcppDeepState_IntegerVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/hutilscpp/inst/testfiles/do_par_in_hash_int/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  IntegerVector table  = RcppDeepState_IntegerVector();
  qs::c_qsave(table,"/home/akhila/fuzzer_packages/fuzzedpackages/hutilscpp/inst/testfiles/do_par_in_hash_int/inputs/table.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "table values: "<< table << std::endl;
  IntegerVector nThread(1);
  nThread[0]  = RcppDeepState_int();
  qs::c_qsave(nThread,"/home/akhila/fuzzer_packages/fuzzedpackages/hutilscpp/inst/testfiles/do_par_in_hash_int/inputs/nThread.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "nThread values: "<< nThread << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    do_par_in_hash_int(x,table,nThread[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}