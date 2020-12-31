#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector hills3d1p12(NumericVector cv1, NumericVector cv2, NumericVector cv3, double width1, double width2, double width3, NumericVector heights, int n, int tmin, int tmax);

TEST(metadynminer3d_deepstate_test,hills3d1p12_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector cv1  = RcppDeepState_NumericVector();
  qs::c_qsave(cv1,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/cv1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cv1 values: "<< cv1 << std::endl;
  NumericVector cv2  = RcppDeepState_NumericVector();
  qs::c_qsave(cv2,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/cv2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cv2 values: "<< cv2 << std::endl;
  NumericVector cv3  = RcppDeepState_NumericVector();
  qs::c_qsave(cv3,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/cv3.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cv3 values: "<< cv3 << std::endl;
  std::ofstream width1_stream;
  double width1  = RcppDeepState_double();
  width1_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/width1");
  width1_stream << width1;
  std::cout << "width1 values: "<< width1 << std::endl;
  width1_stream.close();
  std::ofstream width2_stream;
  double width2  = RcppDeepState_double();
  width2_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/width2");
  width2_stream << width2;
  std::cout << "width2 values: "<< width2 << std::endl;
  width2_stream.close();
  std::ofstream width3_stream;
  double width3  = RcppDeepState_double();
  width3_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/width3");
  width3_stream << width3;
  std::cout << "width3 values: "<< width3 << std::endl;
  width3_stream.close();
  NumericVector heights  = RcppDeepState_NumericVector();
  qs::c_qsave(heights,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/heights.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "heights values: "<< heights << std::endl;
  std::ofstream n_stream;
  int n  = RcppDeepState_int();
  n_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::ofstream tmin_stream;
  int tmin  = RcppDeepState_int();
  tmin_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/tmin");
  tmin_stream << tmin;
  std::cout << "tmin values: "<< tmin << std::endl;
  tmin_stream.close();
  std::ofstream tmax_stream;
  int tmax  = RcppDeepState_int();
  tmax_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/hills3d1p12/inputs/tmax");
  tmax_stream << tmax;
  std::cout << "tmax values: "<< tmax << std::endl;
  tmax_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    hills3d1p12(cv1,cv2,cv3,width1,width2,width3,heights,n,tmin,tmax);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
