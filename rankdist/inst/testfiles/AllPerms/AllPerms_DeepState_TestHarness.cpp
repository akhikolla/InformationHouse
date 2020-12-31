#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix AllPerms(int nobj);

TEST(rankdist_deepstate_test,AllPerms_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  std::ofstream nobj_stream;
  int nobj  = RcppDeepState_int();
  nobj_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rankdist/inst/testfiles/AllPerms/inputs/nobj");
  nobj_stream << nobj;
  std::cout << "nobj values: "<< nobj << std::endl;
  nobj_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    AllPerms(nobj);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
