#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector fe3dp3(NumericVector cv1, NumericVector cv2, NumericVector cv3, NumericVector width1, NumericVector width2, NumericVector width3, NumericVector heights, double x, double y, double z, double p3, int tmin, int tmax);

TEST(metadynminer3d_deepstate_test,fe3dp3_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector cv1  = RcppDeepState_NumericVector();
  qs::c_qsave(cv1,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/cv1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cv1 values: "<< cv1 << std::endl;
  NumericVector cv2  = RcppDeepState_NumericVector();
  qs::c_qsave(cv2,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/cv2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cv2 values: "<< cv2 << std::endl;
  NumericVector cv3  = RcppDeepState_NumericVector();
  qs::c_qsave(cv3,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/cv3.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cv3 values: "<< cv3 << std::endl;
  NumericVector width1  = RcppDeepState_NumericVector();
  qs::c_qsave(width1,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/width1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "width1 values: "<< width1 << std::endl;
  NumericVector width2  = RcppDeepState_NumericVector();
  qs::c_qsave(width2,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/width2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "width2 values: "<< width2 << std::endl;
  NumericVector width3  = RcppDeepState_NumericVector();
  qs::c_qsave(width3,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/width3.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "width3 values: "<< width3 << std::endl;
  NumericVector heights  = RcppDeepState_NumericVector();
  qs::c_qsave(heights,"/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/heights.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "heights values: "<< heights << std::endl;
  std::ofstream x_stream;
  double x  = RcppDeepState_double();
  x_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  std::ofstream y_stream;
  double y  = RcppDeepState_double();
  y_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  std::ofstream z_stream;
  double z  = RcppDeepState_double();
  z_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/z");
  z_stream << z;
  std::cout << "z values: "<< z << std::endl;
  z_stream.close();
  std::ofstream p3_stream;
  double p3  = RcppDeepState_double();
  p3_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/p3");
  p3_stream << p3;
  std::cout << "p3 values: "<< p3 << std::endl;
  p3_stream.close();
  std::ofstream tmin_stream;
  int tmin  = RcppDeepState_int();
  tmin_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/tmin");
  tmin_stream << tmin;
  std::cout << "tmin values: "<< tmin << std::endl;
  tmin_stream.close();
  std::ofstream tmax_stream;
  int tmax  = RcppDeepState_int();
  tmax_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/metadynminer3d/inst/testfiles/fe3dp3/inputs/tmax");
  tmax_stream << tmax;
  std::cout << "tmax values: "<< tmax << std::endl;
  tmax_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    fe3dp3(cv1,cv2,cv3,width1,width2,width3,heights,x,y,z,p3,tmin,tmax);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
