#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

arma::imat random_contingency_table(int setLength, double baserate, double kappaMin, double kappaMax, double minPrecision, double maxPrecision);

TEST(rhoR_deepstate_test,random_contingency_table_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  std::ofstream setLength_stream;
  int setLength  = RcppDeepState_int();
  setLength_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/random_contingency_table/inputs/setLength");
  setLength_stream << setLength;
  std::cout << "setLength values: "<< setLength << std::endl;
  setLength_stream.close();
  std::ofstream baserate_stream;
  double baserate  = RcppDeepState_double();
  baserate_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/random_contingency_table/inputs/baserate");
  baserate_stream << baserate;
  std::cout << "baserate values: "<< baserate << std::endl;
  baserate_stream.close();
  std::ofstream kappaMin_stream;
  double kappaMin  = RcppDeepState_double();
  kappaMin_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/random_contingency_table/inputs/kappaMin");
  kappaMin_stream << kappaMin;
  std::cout << "kappaMin values: "<< kappaMin << std::endl;
  kappaMin_stream.close();
  std::ofstream kappaMax_stream;
  double kappaMax  = RcppDeepState_double();
  kappaMax_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/random_contingency_table/inputs/kappaMax");
  kappaMax_stream << kappaMax;
  std::cout << "kappaMax values: "<< kappaMax << std::endl;
  kappaMax_stream.close();
  std::ofstream minPrecision_stream;
  double minPrecision  = RcppDeepState_double();
  minPrecision_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/random_contingency_table/inputs/minPrecision");
  minPrecision_stream << minPrecision;
  std::cout << "minPrecision values: "<< minPrecision << std::endl;
  minPrecision_stream.close();
  std::ofstream maxPrecision_stream;
  double maxPrecision  = RcppDeepState_double();
  maxPrecision_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/random_contingency_table/inputs/maxPrecision");
  maxPrecision_stream << maxPrecision;
  std::cout << "maxPrecision values: "<< maxPrecision << std::endl;
  maxPrecision_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    random_contingency_table(setLength,baserate,kappaMin,kappaMax,minPrecision,maxPrecision);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
