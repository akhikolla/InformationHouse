#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

int get_k_star(NumericVector x_star, int k_hat, int G_l, int G_r);

TEST(mosum_deepstate_test,get_k_star_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector x_star  = RcppDeepState_NumericVector();
  qs::c_qsave(x_star,"/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/get_k_star/inputs/x_star.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x_star values: "<< x_star << std::endl;
  std::ofstream k_hat_stream;
  int k_hat  = RcppDeepState_int();
  k_hat_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/get_k_star/inputs/k_hat");
  k_hat_stream << k_hat;
  std::cout << "k_hat values: "<< k_hat << std::endl;
  k_hat_stream.close();
  std::ofstream G_l_stream;
  int G_l  = RcppDeepState_int();
  G_l_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/get_k_star/inputs/G_l");
  G_l_stream << G_l;
  std::cout << "G_l values: "<< G_l << std::endl;
  G_l_stream.close();
  std::ofstream G_r_stream;
  int G_r  = RcppDeepState_int();
  G_r_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/get_k_star/inputs/G_r");
  G_r_stream << G_r;
  std::cout << "G_r values: "<< G_r << std::endl;
  G_r_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    get_k_star(x_star,k_hat,G_l,G_r);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
