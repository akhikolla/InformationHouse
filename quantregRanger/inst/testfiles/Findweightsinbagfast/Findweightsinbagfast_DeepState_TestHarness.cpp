#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::NumericVector Findweightsinbagfast(NumericVector ONv, NumericVector OrdNv, NumericVector filterednodes, IntegerVector index, IntegerVector newindex, IntegerVector inbag, NumericVector WEv, int nobs, int ntree, double thres, int l);

TEST(quantregRanger_deepstate_test,Findweightsinbagfast_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector ONv  = RcppDeepState_NumericVector();
  qs::c_qsave(ONv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/ONv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ONv values: "<< ONv << std::endl;
  NumericVector OrdNv  = RcppDeepState_NumericVector();
  qs::c_qsave(OrdNv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/OrdNv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "OrdNv values: "<< OrdNv << std::endl;
  NumericVector filterednodes  = RcppDeepState_NumericVector();
  qs::c_qsave(filterednodes,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/filterednodes.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "filterednodes values: "<< filterednodes << std::endl;
  IntegerVector index  = RcppDeepState_IntegerVector();
  qs::c_qsave(index,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/index.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "index values: "<< index << std::endl;
  IntegerVector newindex  = RcppDeepState_IntegerVector();
  qs::c_qsave(newindex,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/newindex.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "newindex values: "<< newindex << std::endl;
  IntegerVector inbag  = RcppDeepState_IntegerVector();
  qs::c_qsave(inbag,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/inbag.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "inbag values: "<< inbag << std::endl;
  NumericVector WEv  = RcppDeepState_NumericVector();
  qs::c_qsave(WEv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/WEv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "WEv values: "<< WEv << std::endl;
  std::ofstream nobs_stream;
  int nobs  = RcppDeepState_int();
  nobs_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/nobs");
  nobs_stream << nobs;
  std::cout << "nobs values: "<< nobs << std::endl;
  nobs_stream.close();
  std::ofstream ntree_stream;
  int ntree  = RcppDeepState_int();
  ntree_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/ntree");
  ntree_stream << ntree;
  std::cout << "ntree values: "<< ntree << std::endl;
  ntree_stream.close();
  std::ofstream thres_stream;
  double thres  = RcppDeepState_double();
  thres_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/thres");
  thres_stream << thres;
  std::cout << "thres values: "<< thres << std::endl;
  thres_stream.close();
  std::ofstream l_stream;
  int l  = RcppDeepState_int();
  l_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfast/inputs/l");
  l_stream << l;
  std::cout << "l values: "<< l << std::endl;
  l_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    Findweightsinbagfast(ONv,OrdNv,filterednodes,index,newindex,inbag,WEv,nobs,ntree,thres,l);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
