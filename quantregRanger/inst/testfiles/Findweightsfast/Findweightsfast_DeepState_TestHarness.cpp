#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::NumericVector Findweightsfast(NumericVector OrdNv, NumericVector NNv, NumericVector filterednodes, IntegerVector index, IntegerVector newindex, NumericVector WEv, int nobs, int nnew, int ntree, double thres, int l);

TEST(quantregRanger_deepstate_test,Findweightsfast_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector OrdNv  = RcppDeepState_NumericVector();
  qs::c_qsave(OrdNv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/OrdNv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "OrdNv values: "<< OrdNv << std::endl;
  NumericVector NNv  = RcppDeepState_NumericVector();
  qs::c_qsave(NNv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/NNv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "NNv values: "<< NNv << std::endl;
  NumericVector filterednodes  = RcppDeepState_NumericVector();
  qs::c_qsave(filterednodes,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/filterednodes.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "filterednodes values: "<< filterednodes << std::endl;
  IntegerVector index  = RcppDeepState_IntegerVector();
  qs::c_qsave(index,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/index.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "index values: "<< index << std::endl;
  IntegerVector newindex  = RcppDeepState_IntegerVector();
  qs::c_qsave(newindex,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/newindex.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "newindex values: "<< newindex << std::endl;
  NumericVector WEv  = RcppDeepState_NumericVector();
  qs::c_qsave(WEv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/WEv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "WEv values: "<< WEv << std::endl;
  std::ofstream nobs_stream;
  int nobs  = RcppDeepState_int();
  nobs_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/nobs");
  nobs_stream << nobs;
  std::cout << "nobs values: "<< nobs << std::endl;
  nobs_stream.close();
  std::ofstream nnew_stream;
  int nnew  = RcppDeepState_int();
  nnew_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/nnew");
  nnew_stream << nnew;
  std::cout << "nnew values: "<< nnew << std::endl;
  nnew_stream.close();
  std::ofstream ntree_stream;
  int ntree  = RcppDeepState_int();
  ntree_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/ntree");
  ntree_stream << ntree;
  std::cout << "ntree values: "<< ntree << std::endl;
  ntree_stream.close();
  std::ofstream thres_stream;
  double thres  = RcppDeepState_double();
  thres_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/thres");
  thres_stream << thres;
  std::cout << "thres values: "<< thres << std::endl;
  thres_stream.close();
  std::ofstream l_stream;
  int l  = RcppDeepState_int();
  l_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsfast/inputs/l");
  l_stream << l;
  std::cout << "l values: "<< l << std::endl;
  l_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    Findweightsfast(OrdNv,NNv,filterednodes,index,newindex,WEv,nobs,nnew,ntree,thres,l);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
