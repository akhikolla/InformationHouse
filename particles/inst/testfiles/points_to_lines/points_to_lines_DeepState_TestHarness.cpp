#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

List points_to_lines(NumericMatrix line1, NumericMatrix line2, NumericMatrix point);

TEST(particles_deepstate_test,points_to_lines_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix line1  = RcppDeepState_NumericMatrix();
  qs::c_qsave(line1,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/points_to_lines/inputs/line1.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "line1 values: "<< line1 << std::endl;
  NumericMatrix line2  = RcppDeepState_NumericMatrix();
  qs::c_qsave(line2,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/points_to_lines/inputs/line2.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "line2 values: "<< line2 << std::endl;
  NumericMatrix point  = RcppDeepState_NumericMatrix();
  qs::c_qsave(point,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/points_to_lines/inputs/point.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "point values: "<< point << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    points_to_lines(line1,line2,point);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
