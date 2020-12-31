#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix collision(NumericMatrix pos, NumericMatrix vel, NumericVector radii, double strength);

TEST(particles_deepstate_test,collision_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix pos  = RcppDeepState_NumericMatrix();
  qs::c_qsave(pos,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/collision/inputs/pos.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "pos values: "<< pos << std::endl;
  NumericMatrix vel  = RcppDeepState_NumericMatrix();
  qs::c_qsave(vel,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/collision/inputs/vel.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "vel values: "<< vel << std::endl;
  NumericVector radii  = RcppDeepState_NumericVector();
  qs::c_qsave(radii,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/collision/inputs/radii.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "radii values: "<< radii << std::endl;
  std::ofstream strength_stream;
  double strength  = RcppDeepState_double();
  strength_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/collision/inputs/strength");
  strength_stream << strength;
  std::cout << "strength values: "<< strength << std::endl;
  strength_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    collision(pos,vel,radii,strength);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
