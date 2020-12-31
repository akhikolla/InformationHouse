#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix nbody(NumericMatrix pos, NumericVector strength, double theta, double min_dist, double max_dist, double alpha);

TEST(particles_deepstate_test,nbody_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix pos  = RcppDeepState_NumericMatrix();
  qs::c_qsave(pos,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/nbody/inputs/pos.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "pos values: "<< pos << std::endl;
  NumericVector strength  = RcppDeepState_NumericVector();
  qs::c_qsave(strength,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/nbody/inputs/strength.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "strength values: "<< strength << std::endl;
  std::ofstream theta_stream;
  double theta  = RcppDeepState_double();
  theta_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/nbody/inputs/theta");
  theta_stream << theta;
  std::cout << "theta values: "<< theta << std::endl;
  theta_stream.close();
  std::ofstream min_dist_stream;
  double min_dist  = RcppDeepState_double();
  min_dist_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/nbody/inputs/min_dist");
  min_dist_stream << min_dist;
  std::cout << "min_dist values: "<< min_dist << std::endl;
  min_dist_stream.close();
  std::ofstream max_dist_stream;
  double max_dist  = RcppDeepState_double();
  max_dist_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/nbody/inputs/max_dist");
  max_dist_stream << max_dist;
  std::cout << "max_dist values: "<< max_dist << std::endl;
  max_dist_stream.close();
  std::ofstream alpha_stream;
  double alpha  = RcppDeepState_double();
  alpha_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/nbody/inputs/alpha");
  alpha_stream << alpha;
  std::cout << "alpha values: "<< alpha << std::endl;
  alpha_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    nbody(pos,strength,theta,min_dist,max_dist,alpha);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
