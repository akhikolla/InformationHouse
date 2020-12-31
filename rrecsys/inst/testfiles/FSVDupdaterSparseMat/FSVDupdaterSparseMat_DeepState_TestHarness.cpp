#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

List FSVDupdaterSparseMat(NumericMatrix sparseRatingMat, double learningRate, double regCoef, int nrfeat, int steps, int nr_users, int nr_items);

TEST(rrecsys_deepstate_test,FSVDupdaterSparseMat_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix sparseRatingMat  = RcppDeepState_NumericMatrix();
  qs::c_qsave(sparseRatingMat,"/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/FSVDupdaterSparseMat/inputs/sparseRatingMat.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "sparseRatingMat values: "<< sparseRatingMat << std::endl;
  std::ofstream learningRate_stream;
  double learningRate  = RcppDeepState_double();
  learningRate_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/FSVDupdaterSparseMat/inputs/learningRate");
  learningRate_stream << learningRate;
  std::cout << "learningRate values: "<< learningRate << std::endl;
  learningRate_stream.close();
  std::ofstream regCoef_stream;
  double regCoef  = RcppDeepState_double();
  regCoef_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/FSVDupdaterSparseMat/inputs/regCoef");
  regCoef_stream << regCoef;
  std::cout << "regCoef values: "<< regCoef << std::endl;
  regCoef_stream.close();
  std::ofstream nrfeat_stream;
  int nrfeat  = RcppDeepState_int();
  nrfeat_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/FSVDupdaterSparseMat/inputs/nrfeat");
  nrfeat_stream << nrfeat;
  std::cout << "nrfeat values: "<< nrfeat << std::endl;
  nrfeat_stream.close();
  std::ofstream steps_stream;
  int steps  = RcppDeepState_int();
  steps_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/FSVDupdaterSparseMat/inputs/steps");
  steps_stream << steps;
  std::cout << "steps values: "<< steps << std::endl;
  steps_stream.close();
  std::ofstream nr_users_stream;
  int nr_users  = RcppDeepState_int();
  nr_users_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/FSVDupdaterSparseMat/inputs/nr_users");
  nr_users_stream << nr_users;
  std::cout << "nr_users values: "<< nr_users << std::endl;
  nr_users_stream.close();
  std::ofstream nr_items_stream;
  int nr_items  = RcppDeepState_int();
  nr_items_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/FSVDupdaterSparseMat/inputs/nr_items");
  nr_items_stream << nr_items;
  std::cout << "nr_items values: "<< nr_items << std::endl;
  nr_items_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    FSVDupdaterSparseMat(sparseRatingMat,learningRate,regCoef,nrfeat,steps,nr_users,nr_items);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
