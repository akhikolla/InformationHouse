// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// sequence_hit_ids_cpp_DeepState_TestHarness_generation.cpp and sequence_hit_ids_cpp_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector sequence_hit_ids_cpp(NumericVector con, NumericVector subcon, NumericVector pos, NumericVector term_i, double length);

TEST(corpustools_deepstate_test,sequence_hit_ids_cpp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector con  = RcppDeepState_NumericVector();
  qs::c_qsave(con,"/home/akhila/fuzzer_packages/fuzzedpackages/corpustools/inst/testfiles/sequence_hit_ids_cpp/inputs/con.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "con values: "<< con << std::endl;
  NumericVector subcon  = RcppDeepState_NumericVector();
  qs::c_qsave(subcon,"/home/akhila/fuzzer_packages/fuzzedpackages/corpustools/inst/testfiles/sequence_hit_ids_cpp/inputs/subcon.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "subcon values: "<< subcon << std::endl;
  NumericVector pos  = RcppDeepState_NumericVector();
  qs::c_qsave(pos,"/home/akhila/fuzzer_packages/fuzzedpackages/corpustools/inst/testfiles/sequence_hit_ids_cpp/inputs/pos.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "pos values: "<< pos << std::endl;
  NumericVector term_i  = RcppDeepState_NumericVector();
  qs::c_qsave(term_i,"/home/akhila/fuzzer_packages/fuzzedpackages/corpustools/inst/testfiles/sequence_hit_ids_cpp/inputs/term_i.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "term_i values: "<< term_i << std::endl;
  NumericVector length(1);
  length[0]  = RcppDeepState_double();
  qs::c_qsave(length,"/home/akhila/fuzzer_packages/fuzzedpackages/corpustools/inst/testfiles/sequence_hit_ids_cpp/inputs/length.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "length values: "<< length << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    sequence_hit_ids_cpp(con,subcon,pos,term_i,length[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}