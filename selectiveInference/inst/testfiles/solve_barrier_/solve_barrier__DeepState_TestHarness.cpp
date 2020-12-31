// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// solve_barrier__DeepState_TestHarness_generation.cpp and solve_barrier__DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::List solve_barrier_(Rcpp::NumericVector conjugate_arg, Rcpp::NumericMatrix precision, Rcpp::NumericVector feasible_point, int max_iter, int min_iter, double value_tol, double initial_step);

TEST(selectiveInference_deepstate_test,solve_barrier__test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector conjugate_arg  = RcppDeepState_NumericVector();
  qs::c_qsave(conjugate_arg,"/home/akhila/fuzzer_packages/fuzzedpackages/selectiveInference/inst/testfiles/solve_barrier_/inputs/conjugate_arg.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "conjugate_arg values: "<< conjugate_arg << std::endl;
  NumericMatrix precision  = RcppDeepState_NumericMatrix();
  qs::c_qsave(precision,"/home/akhila/fuzzer_packages/fuzzedpackages/selectiveInference/inst/testfiles/solve_barrier_/inputs/precision.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "precision values: "<< precision << std::endl;
  NumericVector feasible_point  = RcppDeepState_NumericVector();
  qs::c_qsave(feasible_point,"/home/akhila/fuzzer_packages/fuzzedpackages/selectiveInference/inst/testfiles/solve_barrier_/inputs/feasible_point.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "feasible_point values: "<< feasible_point << std::endl;
  IntegerVector max_iter(1);
  max_iter[0]  = RcppDeepState_int();
  qs::c_qsave(max_iter,"/home/akhila/fuzzer_packages/fuzzedpackages/selectiveInference/inst/testfiles/solve_barrier_/inputs/max_iter.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "max_iter values: "<< max_iter << std::endl;
  IntegerVector min_iter(1);
  min_iter[0]  = RcppDeepState_int();
  qs::c_qsave(min_iter,"/home/akhila/fuzzer_packages/fuzzedpackages/selectiveInference/inst/testfiles/solve_barrier_/inputs/min_iter.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "min_iter values: "<< min_iter << std::endl;
  NumericVector value_tol(1);
  value_tol[0]  = RcppDeepState_double();
  qs::c_qsave(value_tol,"/home/akhila/fuzzer_packages/fuzzedpackages/selectiveInference/inst/testfiles/solve_barrier_/inputs/value_tol.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "value_tol values: "<< value_tol << std::endl;
  NumericVector initial_step(1);
  initial_step[0]  = RcppDeepState_double();
  qs::c_qsave(initial_step,"/home/akhila/fuzzer_packages/fuzzedpackages/selectiveInference/inst/testfiles/solve_barrier_/inputs/initial_step.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "initial_step values: "<< initial_step << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    solve_barrier_(conjugate_arg,precision,feasible_point,max_iter[0],min_iter[0],value_tol[0],initial_step[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}