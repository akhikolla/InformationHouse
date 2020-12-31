// Enable C++11 via this plugin
// [[Rcpp::plugins("cpp11")]]
#include "ProbabilisticLinkageWrapper.h"


//'Wrapper function for Probabilistic Linkage

//' @param IDA_ is the input string
//' @param dataA_ is the size of the qgram

//' @return the input string transformed into a vector of qgrams

//[[Rcpp::export(".ProbabilisticLinkagec")]]
DataFrame ProbabilisticLinkagec(SEXP IDA_, SEXP dataA_, SEXP blockingdataA_,
                                SEXP IDB_, SEXP dataB_, SEXP blockingdataB_, SEXP method_,
                                SEXP blocking_, SEXP threshold_, SEXP lenNgram_,
                                LogicalVector ind_c0_ = LogicalVector::create(0),
                                LogicalVector ind_c1_ = LogicalVector::create(0),
                                NumericVector m_ = NumericVector::create(0.9),
                                NumericVector u_ = NumericVector::create(0.1),
                                NumericVector p_ = NumericVector::create(0.05),
                                NumericVector e = NumericVector::create(0.0004),
                                double upper = 0.0, double lower = 0.0,
                                double jaroWeightFactor = 1.0) {
  //transform data from input
  vector < string > method = Rcpp::as < std::vector<string> > (method_);
  string blocking = Rcpp::as < string > (blocking_);
  vector<float> m = Rcpp::as < std::vector<float> > (m_);
  vector<float> u = Rcpp::as < std::vector<float> > (u_);
  vector<float> epsilon = Rcpp::as < std::vector<float> > (e);
  vector<float> p = Rcpp::as < std::vector<float> > (p_);
  //in_c only used for JWMcLWL
  vector<bool> ind_c0 = Rcpp::as < std::vector<bool> > (ind_c0_); // only used for JWMcLWL
  vector<bool> ind_c1 = Rcpp::as < std::vector<bool> > (ind_c1_);
  vector<int> lenNgram = Rcpp::as < std::vector<int> > (lenNgram_);
  float threshold = Rcpp::as<float>(threshold_);
  //  // Handle IDs, Ids can either be of type int or character
  vector<string> ID1, ID2;
  vector<int> ID1Int, ID2Int;
  int ID1size, ID2size;
  if (TYPEOF(IDA_) == STRSXP) { // case vector of IDA_ is of type character
    ID1 = Rcpp::as < std::vector<string> > (IDA_);
    ID1size = ID1.size();
  } else if (TYPEOF(IDA_) == INTSXP) { // case vector of IDA_ is of type int
    ID1Int = Rcpp::as < std::vector<int> > (IDA_);
    ID1size = ID1Int.size();
    Rcpp::Rcout << "Info: IDA is of type int. It will be transformed to character."
                << endl;
    if (ID1size > 0) {
      for (int j = 0; j < ID1size; j++) {
        ID1.push_back(to_string(ID1Int[j])); //Copy the vector to the string
      }
    }
  } else {
    Rcpp::Rcerr << "Input IDA must be a vector of Type character or int.\n ";
    return 0;
  }

  if (TYPEOF(IDB_) == STRSXP ) { // case vector of ID2 is of type character
    ID2 = Rcpp::as < std::vector<string> > (IDB_);
    ID2size = ID2.size();
  } else if (TYPEOF(IDB_) == INTSXP) { // case vector of IDB_ is of type int
    ID2Int = Rcpp::as < std::vector<int> > (IDB_);
    ID2size = ID2Int.size();
    Rcpp::Rcout << "Info: IDB is of type int. It will be transformed to character."
                << endl;
    if (ID2size > 0) {
      for (int j = 0; j < ID2size; j++) {
        ID2.push_back(to_string(ID2Int[j])); //Copy the vector to the string
      }
    }
  } else {
    Rcpp::Rcerr << "Input IDB must be a vector of Type character or int.\n";
    return 0;
  }
  if (ID1.size() <= 1|| ID2.size() <=1){
    Rcpp::Rcerr << "Length of the ID vector has to be bigger than one.\n" << endl;
    return DataFrame::create(Named("stringsAsFactors") = false);
  }
  //Handle variables, vars can either be of type int, float or character
  vector < vector<string> > varsV1;
  vector < vector<string> > varsV2;
  vector < MTB_StringVectorData>  d1;
  vector < MTB_StringVectorData > d2;
  vector<string> blockingvarsV1, blockingvarsV2;
  int dataSize = 1;
  //if input is a list
  if (TYPEOF(dataA_) == 19) {
    Rcpp::List dataAlist(dataA_);
    Rcpp::List dataBlist(dataB_);
    dataSize= dataAlist.size();
    for(int i = 0; i < dataSize; i++ ){

      varsV1.push_back(prepareData(dataAlist[i], "Linking vars1", true));
      varsV2.push_back(prepareData(dataBlist[i], "Linking vars2", true));
    }
    if (!Rf_isNull(blockingdataA_) ) {
      blockingvarsV1 = prepareData(blockingdataA_, "Blocking vars1", false);
    } else {
      for (unsigned i = 0; i < varsV1[0].size(); i++)
        blockingvarsV1.push_back("0");
    }

    if (!Rf_isNull(blockingdataB_) ) {
      blockingvarsV2 = prepareData(blockingdataB_, "Blocking vars2", false);
    } else {
      for (unsigned i = 0; i < varsV2[0].size(); i++)
        blockingvarsV2.push_back("0");
    }
    for(int i = 0; i < dataSize; i++ ){
      d1.push_back(MTB_StringVectorData("var1", ID1, varsV1[i], blockingvarsV1));
      d2.push_back(MTB_StringVectorData("var2", ID2, varsV2[i], blockingvarsV2));
    }
  } else { //if input is a vector
    varsV1.push_back(prepareData(dataA_, "Linking vars1", true));
    varsV2.push_back(prepareData(dataB_, "Linking vars2", true));

    // // Handle blocking variables,  blocking vars can either be of type int, float or character

    if (!Rf_isNull(blockingdataA_) ) {
      blockingvarsV1 = prepareData(blockingdataA_, "Blocking vars1", false);
    } else {
      for (unsigned i = 0; i < varsV1[0].size(); i++)
        blockingvarsV1.push_back("0");
    }

    if (!Rf_isNull(blockingdataB_) ) {
      blockingvarsV2 = prepareData(blockingdataB_, "Blocking vars2", false);
    } else {
      for (unsigned i = 0; i < varsV2[0].size(); i++)
        blockingvarsV2.push_back("0");
    }
    //set StringVectorData
    d1.push_back(MTB_StringVectorData("var1", ID1, varsV1[0], blockingvarsV1));
    d2.push_back(MTB_StringVectorData("var2", ID2, varsV2[0], blockingvarsV2));
  }

  //set merge configuration

  MergingConfiguration merging;
  vector<MergingConfiguration> mc;
  for (int i = 0; i < dataSize; i++) {
    merging.setAlgorithm(method[i]);
    merging.setBlocking(blocking);
    merging.setM(m[i]);
    merging.setU(u[i]);
    merging.setEpsilon(epsilon[i]);
    merging.setP(p[i]);
    merging.setInd_c0(ind_c0[i]);
    merging.setInd_c1(ind_c1[i]);
    merging.setJaroWeightFactor(jaroWeightFactor);
    merging.setLower(lower);
    merging.setUpper(upper);
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
  vector<float> tani; //vector holding resulting tanimoto values
  vector <string > match;//classification
  //results are calculated and added
  res = emWrapper(d1[0], d2[0], mc, 0);
  tani = res.getRes();
  if (dataSize ==1 ){
    for (unsigned i = 0; i < tani.size(); i++){
      if (tani[i] >= threshold){
        match.push_back("M");
      } else {
        match.push_back("NM");
      }
    }
    return DataFrame::create(Named("ID1") = wrap(res.getID1()),
                             Named("ID2") = wrap(res.getID2()),
                             Named("sim") = wrap(res.getRes()),
                             Named("classification") =  wrap(match),
                             Named("stringsAsFactors") = false);
  }
  else if (dataSize > 1){

    for (int i = 1 ; i < dataSize; i++){ // res = emWrapper(d1, d2, mc);
      res = emWrapper(d1[i], d2[i], mc, i);
      std::transform(tani.begin(), tani.end(), res.getRes().begin(), tani.begin(), std::plus<float>()); // iterates over tani and res.getRes(), applies plus that adds floats, and writes the sum back to tani
    }
    for (unsigned i = 0; i < tani.size(); i++){
      if (tani[i] >= threshold){
        match.push_back("M");
      } else {
        match.push_back("NM");
      }
    }
    return DataFrame::create(Named("ID1") = wrap(res.getID1()),
                             Named("ID2") = wrap(res.getID2()),
                             Named("sim") = wrap(tani),
                             Named("classification") = wrap(match),
                             Named("stringsAsFactors") = false);
  }
  return DataFrame::create(Named("stringsAsFactors") = false);
}

/**
* Calculates distance with em
*/
MTB_Result emWrapper(MTB_StringVectorData d1, MTB_StringVectorData d2,
                     vector<MergingConfiguration> mc, int counterSim) {
  MTB_Result res;
  string classification;
  // set thread to calculate either absolute values or use em algorithm
  MTB_ProbabilityCalculation pcth;
  vector<double> mv;
  mv.push_back(mc[counterSim].getM());
  vector<double> uv;
  uv.push_back(mc[counterSim].getU());
  double p = mc[counterSim].getP();
  double epsilon = mc[counterSim].getEpsilon();
  vector<int> cf = pcth.countFrequencies(d1.getData(), d2.getData());
  MTB_EMAlgorithm emAlg(cf, mv, uv, p, epsilon);
  mc = pcth.run(mc, d1.getData(), d1.getID(), d2.getData(), d2.getID(),
                emAlg);
  vector<MTB_MergeData> mergeData = initMergeData(d1, d2, mc); // blocking takes place
  MTB_InternalProbabilisticComparisonsObservation ipco;
  MTB_StringVectorData m1, m2;
  for (MTB_MergeData m : mergeData) {

    m1 = m.getData1();
    m2 = m.getData2();
    for (unsigned i = 0; i < m1.getData().size(); i++) {
      for (unsigned j = 0; j < m2.getData().size(); j++) {

        ipco.calculateQuality(mc[counterSim], m1.getData()[i],m1.getID()[i],
                              m2.getData()[j], m2.getID()[j],
                                                         mc[counterSim].getJaroWeightFactor(), 0, &res);
      }
    }
  }
  return res;
} // end emWrapper

