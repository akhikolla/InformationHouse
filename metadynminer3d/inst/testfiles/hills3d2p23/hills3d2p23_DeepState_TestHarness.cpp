#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector hills3d2p23(NumericVector cv1, NumericVector cv2, NumericVector cv3, NumericVector width1, NumericVector width2, NumericVector width3, NumericVector heights, int n, int tmin, int tmax);

TEST(metadynminer3d_deepstate_test,hills3d2p23_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector cv1  = RcppDeepState_NumericVector();
  qs::c_qsave(cv1,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/cv1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cv1 values: "<< cv1 << std::endl;
  NumericVector cv2  = RcppDeepState_NumericVector();
  qs::c_qsave(cv2,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/cv2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cv2 values: "<< cv2 << std::endl;
  NumericVector cv3  = RcppDeepState_NumericVector();
  qs::c_qsave(cv3,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/cv3.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cv3 values: "<< cv3 << std::endl;
  NumericVector width1  = RcppDeepState_NumericVector();
  qs::c_qsave(width1,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/width1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "width1 values: "<< width1 << std::endl;
  NumericVector width2  = RcppDeepState_NumericVector();
  qs::c_qsave(width2,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/width2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "width2 values: "<< width2 << std::endl;
  NumericVector width3  = RcppDeepState_NumericVector();
  qs::c_qsave(width3,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/width3.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "width3 values: "<< width3 << std::endl;
  NumericVector heights  = RcppDeepState_NumericVector();
  qs::c_qsave(heights,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/heights.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "heights values: "<< heights << std::endl;
  std::ofstream n_stream;
  int n  = RcppDeepState_int();
  n_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::ofstream tmin_stream;
  int tmin  = RcppDeepState_int();
  tmin_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/tmin");
  tmin_stream << tmin;
  std::cout << "tmin values: "<< tmin << std::endl;
  tmin_stream.close();
  std::ofstream tmax_stream;
  int tmax  = RcppDeepState_int();
  tmax_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d2p23/inputs/tmax");
  tmax_stream << tmax;
  std::cout << "tmax values: "<< tmax << std::endl;
  tmax_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    hills3d2p23(cv1,cv2,cv3,width1,width2,width3,heights,n,tmin,tmax);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
