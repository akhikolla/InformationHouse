#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector bilinear(NumericVector x_breaks, NumericVector y_breaks, NumericMatrix grid, NumericVector x, NumericVector y);

TEST(particles_deepstate_test,bilinear_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector x_breaks  = RcppDeepState_NumericVector();
  qs::c_qsave(x_breaks,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/bilinear/inputs/x_breaks.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x_breaks values: "<< x_breaks << std::endl;
  NumericVector y_breaks  = RcppDeepState_NumericVector();
  qs::c_qsave(y_breaks,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/bilinear/inputs/y_breaks.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "y_breaks values: "<< y_breaks << std::endl;
  NumericMatrix grid  = RcppDeepState_NumericMatrix();
  qs::c_qsave(grid,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/bilinear/inputs/grid.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "grid values: "<< grid << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/bilinear/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  NumericVector y  = RcppDeepState_NumericVector();
  qs::c_qsave(y,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/bilinear/inputs/y.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "y values: "<< y << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    bilinear(x_breaks,y_breaks,grid,x,y);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
