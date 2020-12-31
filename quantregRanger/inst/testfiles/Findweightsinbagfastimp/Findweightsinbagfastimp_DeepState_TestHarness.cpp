#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::List Findweightsinbagfastimp(NumericVector ONv, NumericVector ONvp, NumericVector OrdNv, NumericVector filterednodes, IntegerVector index, IntegerVector newindex, IntegerVector inbag, NumericVector WEv, NumericVector WEvp, int npred, int nobs, int ntree, double thres, int l, IntegerVector countbreak);

TEST(quantregRanger_deepstate_test,Findweightsinbagfastimp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector ONv  = RcppDeepState_NumericVector();
  qs::c_qsave(ONv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/ONv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ONv values: "<< ONv << std::endl;
  NumericVector ONvp  = RcppDeepState_NumericVector();
  qs::c_qsave(ONvp,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/ONvp.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ONvp values: "<< ONvp << std::endl;
  NumericVector OrdNv  = RcppDeepState_NumericVector();
  qs::c_qsave(OrdNv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/OrdNv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "OrdNv values: "<< OrdNv << std::endl;
  NumericVector filterednodes  = RcppDeepState_NumericVector();
  qs::c_qsave(filterednodes,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/filterednodes.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "filterednodes values: "<< filterednodes << std::endl;
  IntegerVector index  = RcppDeepState_IntegerVector();
  qs::c_qsave(index,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/index.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "index values: "<< index << std::endl;
  IntegerVector newindex  = RcppDeepState_IntegerVector();
  qs::c_qsave(newindex,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/newindex.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "newindex values: "<< newindex << std::endl;
  IntegerVector inbag  = RcppDeepState_IntegerVector();
  qs::c_qsave(inbag,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/inbag.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "inbag values: "<< inbag << std::endl;
  NumericVector WEv  = RcppDeepState_NumericVector();
  qs::c_qsave(WEv,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/WEv.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "WEv values: "<< WEv << std::endl;
  NumericVector WEvp  = RcppDeepState_NumericVector();
  qs::c_qsave(WEvp,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/WEvp.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "WEvp values: "<< WEvp << std::endl;
  std::ofstream npred_stream;
  int npred  = RcppDeepState_int();
  npred_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/npred");
  npred_stream << npred;
  std::cout << "npred values: "<< npred << std::endl;
  npred_stream.close();
  std::ofstream nobs_stream;
  int nobs  = RcppDeepState_int();
  nobs_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/nobs");
  nobs_stream << nobs;
  std::cout << "nobs values: "<< nobs << std::endl;
  nobs_stream.close();
  std::ofstream ntree_stream;
  int ntree  = RcppDeepState_int();
  ntree_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/ntree");
  ntree_stream << ntree;
  std::cout << "ntree values: "<< ntree << std::endl;
  ntree_stream.close();
  std::ofstream thres_stream;
  double thres  = RcppDeepState_double();
  thres_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/thres");
  thres_stream << thres;
  std::cout << "thres values: "<< thres << std::endl;
  thres_stream.close();
  std::ofstream l_stream;
  int l  = RcppDeepState_int();
  l_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/l");
  l_stream << l;
  std::cout << "l values: "<< l << std::endl;
  l_stream.close();
  IntegerVector countbreak  = RcppDeepState_IntegerVector();
  qs::c_qsave(countbreak,"/home/akhila/fuzzer_packages/fuzzedpackages/quantregRanger/inst/testfiles/Findweightsinbagfastimp/inputs/countbreak.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "countbreak values: "<< countbreak << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    Findweightsinbagfastimp(ONv,ONvp,OrdNv,filterednodes,index,newindex,inbag,WEv,WEvp,npred,nobs,ntree,thres,l,countbreak);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
