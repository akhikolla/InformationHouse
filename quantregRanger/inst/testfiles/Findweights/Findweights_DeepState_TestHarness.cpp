#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::NumericVector Findweights(NumericVector ONv, NumericVector NNv, NumericVector WEv, int nobs, int nnew, int ntree, double thres, NumericVector counti, int normalise);

TEST(quantregRanger_deepstate_test,Findweights_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector ONv  = RcppDeepState_NumericVector();
  qs::c_qsave(ONv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweights/inputs/ONv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ONv values: "<< ONv << std::endl;
  NumericVector NNv  = RcppDeepState_NumericVector();
  qs::c_qsave(NNv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweights/inputs/NNv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "NNv values: "<< NNv << std::endl;
  NumericVector WEv  = RcppDeepState_NumericVector();
  qs::c_qsave(WEv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweights/inputs/WEv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "WEv values: "<< WEv << std::endl;
  std::ofstream nobs_stream;
  int nobs  = RcppDeepState_int();
  nobs_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweights/inputs/nobs");
  nobs_stream << nobs;
  std::cout << "nobs values: "<< nobs << std::endl;
  nobs_stream.close();
  std::ofstream nnew_stream;
  int nnew  = RcppDeepState_int();
  nnew_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweights/inputs/nnew");
  nnew_stream << nnew;
  std::cout << "nnew values: "<< nnew << std::endl;
  nnew_stream.close();
  std::ofstream ntree_stream;
  int ntree  = RcppDeepState_int();
  ntree_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweights/inputs/ntree");
  ntree_stream << ntree;
  std::cout << "ntree values: "<< ntree << std::endl;
  ntree_stream.close();
  std::ofstream thres_stream;
  double thres  = RcppDeepState_double();
  thres_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweights/inputs/thres");
  thres_stream << thres;
  std::cout << "thres values: "<< thres << std::endl;
  thres_stream.close();
  NumericVector counti  = RcppDeepState_NumericVector();
  qs::c_qsave(counti,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweights/inputs/counti.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "counti values: "<< counti << std::endl;
  std::ofstream normalise_stream;
  int normalise  = RcppDeepState_int();
  normalise_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweights/inputs/normalise");
  normalise_stream << normalise;
  std::cout << "normalise values: "<< normalise << std::endl;
  normalise_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    Findweights(ONv,NNv,WEv,nobs,nnew,ntree,thres,counti,normalise);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
