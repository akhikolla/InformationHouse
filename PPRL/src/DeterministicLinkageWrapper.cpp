// Enable C++11 via this plugin
// [[Rcpp::plugins("cpp11")]]
#include "DeterministicLinkageWrapper.h"


//'Wrapper function for Deterministic Linkage

//' @param IDA_
//' @param dataA_

//' @return

//[[Rcpp::export(".DeterministicLinkagec")]]
DataFrame DeterministicLinkagec(SEXP IDA_, SEXP dataA_, SEXP blockingdataA_,
                                SEXP IDB_, SEXP dataB_, SEXP blockingdataB_, SEXP method_,
                                SEXP blocking_, SEXP threshold_, SEXP lenNgram_,  LogicalVector ind_c0_ =
                                  LogicalVector::create(0), LogicalVector ind_c1_ =
                                  LogicalVector::create(0), int counterSim = 1 ) {
  //transform data from input
  vector < string > method = Rcpp::as < std::vector<string> > (method_);
  string blocking = Rcpp::as < string > (blocking_);
  //in_c only used for JWMcLWL
  vector<bool> ind_c0 = Rcpp::as < std::vector<bool> > (ind_c0_); // only used for JWMcLWL
  vector<bool> ind_c1 = Rcpp::as < std::vector<bool> > (ind_c1_);
  vector<int> lenNgram = Rcpp::as < std::vector<int> > (lenNgram_);
  float threshold = Rcpp::as<float>(threshold_);
  // Handle IDs, Ids can either be of type int or character
  vector<string> ID1, ID2;
  vector<int> ID1Int, ID2Int;
  int ID1size, ID2size;

  if (TYPEOF(IDA_) == STRSXP) { // case vector of ID1 is of type character
    ID1 = Rcpp::as < std::vector<string> > (IDA_);
    ID1size = ID1.size();
  } else if (TYPEOF(IDA_) == INTSXP) { // case vector of ID1 is of type int
    ID1Int = Rcpp::as < std::vector<int> > (IDA_);
    ID1size = ID1Int.size();
    Rcpp::Rcout << "Info: ID1 is of type int. It will be transformed to character."
         << endl;
    if (ID1size > 0) {
      for (int j = 0; j < ID1size; j++) {
        ID1.push_back(to_string(ID1Int[j])); //Copy the vector to the string
      }
    }
  } else {
    Rcpp::Rcerr << "Input ID1 must be a vector of Type character or int.\n";
    return 0;
  }
  if (TYPEOF(IDB_) == STRSXP) { // case vector of ID2 is of type character
    ID2 = Rcpp::as < std::vector<string> > (IDB_);
    ID2size = ID2.size();
  } else if (TYPEOF(IDB_) == INTSXP) { // case vector of ID2 is of type int
    ID2Int = Rcpp::as < std::vector<int> > (IDB_);
    ID2size = ID2Int.size();
    Rcpp::Rcout << "Info: ID2 is of type int. It will be transformed to character."
         << endl;
    if (ID2size > 0) {
      for (int j = 0; j < ID2size; j++) {
        ID2.push_back(to_string(ID2Int[j])); //Copy the vector to the string
      }
    }
  } else {
    Rcpp::Rcerr << "Input ID2 must be a vector of Type character or int.\n";
    return 0;
  }
  //Handle variables, vars can either be of type int, float or character
  vector < string > varsV1 = prepareData(dataA_, "Linking vars1", false);
  vector < string > varsV2 = prepareData(dataB_, "Linking vars2", false);

  // Handle blocking variables,  blocking vars can either be of type int, float or character
  vector<string> blockingvarsV1, blockingvarsV2;
  if (!Rf_isNull(blockingdataA_)) {
    blockingvarsV1 = prepareData(blockingdataA_, "Blocking vars1", false);
  } else {
    for (unsigned i = 0; i < varsV1.size(); i++)
      blockingvarsV1.push_back("0");
  }

  if (!Rf_isNull(blockingdataB_)) {
    blockingvarsV2 = prepareData(blockingdataB_, "Blocking vars2", false);
  } else {
    for (unsigned i = 0; i < varsV2.size(); i++)
      blockingvarsV2.push_back("0");
  }

  //check sizes
  if (varsV1.size() % ID1size != 0) {
    Rcpp::Rcerr << " ID-Vector A and Input-Vector A must have the same size. "
         << varsV1.size() << " " << ID1size << endl;
    return 0;
  }
  if (varsV2.size() % ID2size != 0) {
    Rcpp::Rcerr << " ID-Vector B and Input-Vector B must have the same size. "
         << endl;
    return 0;
  }
  //set StringVectorData
  MTB_StringVectorData d1 = MTB_StringVectorData("var1", ID1, varsV1,
                                                 blockingvarsV1);
  MTB_StringVectorData d2 = MTB_StringVectorData("var2", ID2, varsV2,
                                                 blockingvarsV2);
  //set merge configuration

  MergingConfiguration merging;
  vector<MergingConfiguration> mc;
  for (unsigned i = 0; i < method.size(); i++) {
    merging.setAlgorithm(method[i]);
    merging.setBlocking(blocking);
    merging.setInd_c0(ind_c0[i]);
    merging.setInd_c1(ind_c1[i]);
    merging.setJaroWeightFactor(1.0);
    merging.setThreshold(threshold);
    if (lenNgram[i] > 4||lenNgram[i] < 1) {
      lenNgram[i] = 2;
      Rcpp::Rcerr
        << "Length of ngrams must be between 1 and 4. Length of the qgram will be set to 2."
        << endl;

    }
    merging.setLenNgram(lenNgram[i]);
    mc.push_back(merging);

  }

  MTB_Result res;
  res = distanceWrapper(d1, d2, mc, counterSim);
  if (counterSim == 1) {
    return DataFrame::create(Named("ID1") = wrap(res.getID1()),
                             Named("ID2") = wrap(res.getID2()),
                             Named(mc[counterSim - 1].getAlgorithmName()) = wrap(
                             res.getRes()),
                             Named("stringsAsFactors") = false);
  } else if (counterSim > 1) {
    return DataFrame::create(Named(mc[counterSim - 1].getAlgorithmName()) =
                             wrap(res.getRes()),
                             Named("stringsAsFactors") = false);
  }

  return DataFrame::create(Named("stringsAsFactors") = false);
}



/**
 * Calculates distance without em
 */
MTB_Result distanceWrapper(MTB_StringVectorData d1, MTB_StringVectorData d2,
                           vector<MergingConfiguration> mc, int counterSim) {
  MTB_Result res;
  double distance;
  vector<MTB_MergeData> md = initMergeData(d1, d2, mc); // blocking takes place
  MTB_InternalProbabilisticComparisonsObservation ipco;
  MTB_StringVectorData m1, m2;

  for (MTB_MergeData m : md) {
    m1 = m.getData1();
    m2 = m.getData2();
    for (unsigned i = 0; i < m1.getData().size(); i++) {
      for (unsigned j = 0; j < m2.getData().size(); j++) {

        distance = ipco.calculateDistance(mc[counterSim - 1],
                                          m1.getData()[i], m2.getData()[j]);
        res.addResult(m1.getData()[i], m1.getID()[i],
                      m1.getBlockingData()[i], m2.getData()[j], m2.getID()[j],
                                                                          m2.getBlockingData()[j], distance);
      }
    }
  }
  return res;
} // end DistanceWrapper

