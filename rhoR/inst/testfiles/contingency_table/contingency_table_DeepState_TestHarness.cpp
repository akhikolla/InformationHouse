#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

arma::imat contingency_table(double precision, double rec, int length, double baserate);

TEST(rhoR_deepstate_test,contingency_table_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  std::ofstream precision_stream;
  double precision  = RcppDeepState_double();
  precision_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/contingency_table/inputs/precision");
  precision_stream << precision;
  std::cout << "precision values: "<< precision << std::endl;
  precision_stream.close();
  std::ofstream rec_stream;
  double rec  = RcppDeepState_double();
  rec_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/contingency_table/inputs/rec");
  rec_stream << rec;
  std::cout << "rec values: "<< rec << std::endl;
  rec_stream.close();
  std::ofstream length_stream;
  int length  = RcppDeepState_int();
  length_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/contingency_table/inputs/length");
  length_stream << length;
  std::cout << "length values: "<< length << std::endl;
  length_stream.close();
  std::ofstream baserate_stream;
  double baserate  = RcppDeepState_double();
  baserate_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/contingency_table/inputs/baserate");
  baserate_stream << baserate;
  std::cout << "baserate values: "<< baserate << std::endl;
  baserate_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    contingency_table(precision,rec,length,baserate);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
