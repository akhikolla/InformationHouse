// Enable C++11 via this plugin
// [[Rcpp::plugins("cpp11")]]


#include "CreateCLKWrapper.h"

/**
* Wrapper function for CreateCLK
*/
// [[Rcpp::export]]
DataFrame CreateCLK(SEXP ID, DataFrame data, SEXP password, int k =20, IntegerVector padding = IntegerVector::create(0), IntegerVector qgram = IntegerVector::create(2), unsigned lenBloom =1000 ) {
  DataFrame res;
  if (data.size() >1){
    res = CreateCLKc(ID, data, password, k, padding, qgram, lenBloom );
  } else if (data.size() == 1) {
    Rcpp::Rcerr << "Please use CreateBF\n";
    //res = CreateBFc(ID, data[0], password, k, padding[0], qgram[0], lenBloom );
  }
  return res;
}

DataFrame CreateCLKc(SEXP ID, DataFrame data, SEXP password, int k , IntegerVector padding , IntegerVector qgram, unsigned lenBloom ) {
  int IDsize;
  if(TYPEOF(ID)==STRSXP){
    vector<string> IDtemp = Rcpp::as<std::vector<string> > (ID);
    IDsize = IDtemp.size();
  }
  else if(TYPEOF(ID)==INTSXP) {
    vector<int> IDtemp = Rcpp::as<std::vector<int> > (ID);
    IDsize = IDtemp.size();
  }
  else {
    Rcpp::Rcerr << "Input ID must be a vector of Type character or int.\n";
    return NULL;
  }
  if (data.nrows() != IDsize){
    Rcpp::Rcerr << " ID-Vector and Input-Dataframe must have the same size. " << endl;
    return NULL;
  }

  vector<string>  varsV(data.nrows());
  vector<float>  varsVFloat(data.nrows());
  vector<int> varsVInt(data.nrows());
  vector<string> inputVector(data.size());
  vector<vector<string> > varsM ;
  CharacterVector CLKs(data.nrows());
  vector<int> padding_ =Rcpp::as<std::vector<int> > (padding);
  vector<int> qgram_ =Rcpp::as<std::vector<int> > (qgram);

  //check password vector
  vector<string> password_;
  if(TYPEOF(password)!=STRSXP){
    Rcpp::Rcerr << "Please select a password for each input variable in a vector of class character.";
  }
  else{
    password_ =  Rcpp::as<std::vector<string> > (password);
    if (data.size() != (unsigned)password_.size()){
      Rcpp::Rcerr << "vector of passwords must have the same size as the input data.frame."  ;
      return NULL;}
  }
  //Sizes of Vectors are checked and adapted
  if (data.size() > (unsigned)padding_.size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be fill with zeros."<< endl;
    for (int i = padding_.size() ; i< data.size(); i++){
      padding_.push_back(0);
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
  if (data.size() < (unsigned)padding_.size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be cut."<< endl;
    while((unsigned)padding_.size() > data.size()){
      padding_.pop_back();
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }

  if (data.size() > (unsigned)qgram_.size()){
    Rcpp::Rcerr << "Vector qgrams must have the same size as the input data.frame. Qgrams will be fill with 2s." << endl;
    for (int i = qgram_.size() ; i< data.size(); i++){
      qgram_.push_back(2);
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
  if (data.size() < (unsigned)qgram_.size()){
    Rcpp::Rcerr << "Vector qgram must have the same size as the input data.frame. Qgram will be cut."<< endl;
    while((unsigned) qgram_.size() > data.size()){
      qgram_.pop_back();
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }

  //Matrix is filled with input, Types are set to string
  for (int i = 0 ; i< data.size(); i++){
    if(TYPEOF(data[i])==STRSXP) {

      varsV=Rcpp::as<std::vector<string> > (data[i]);
      varsM.push_back(varsV);
    }
    if(TYPEOF(data[i])==REALSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains floats. Data will be transformed to characters." << endl;
      varsVFloat = Rcpp::as<std::vector<float> > (data[i]);

      if(varsVFloat.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          varsV[j] = to_string(varsVFloat[j]);//Copy the vector to the string
        }
      }
      varsM.push_back(varsV);
    }
    if(TYPEOF(data[i])==INTSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains integers or factors. Make sure the input data are not factors when you want to use characters. Data will be transformed to characters." << endl;
      varsVInt = Rcpp::as<std::vector<int> > (data[i]);
      if(varsVInt.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          varsV[j]=to_string(varsVInt[j]);//Copy the vector to the string
        }
      }
      varsM.push_back(varsV);
    }
  }

  //for each row
  for (int i = 0 ; i<data.nrows(); i++){
    for (int j = 0 ; j<data.size(); j++){
      inputVector[j]=varsM[j][i];
    }
    CLKs[i]=CreateCLKBigramSeed(inputVector, k, padding_, qgram_, lenBloom, password_);
  }

  return DataFrame::create(Named("ID") = ID ,
                           Named("CLKs") = CLKs,
                           Named("stringsAsFactors") = false);
}


/**
* Wrapper function for CreateMarkovCLK
*/
// [[Rcpp::export]]
DataFrame CreateMarkovCLK(SEXP ID, DataFrame data, SEXP password, NumericMatrix markovTable, int k1 =20, int k2 =4,
                          IntegerVector padding =IntegerVector::create(0), IntegerVector qgram= IntegerVector::create(2), unsigned lenBloom =1000, bool includeOriginalBigram = true, bool v = false) {
  List dimnames = markovTable.attr("dimnames");
  vector<string> rownames = dimnames[0];
  vector<string> colnames = dimnames[1];
  NumericVector zz1;
  vector<float> zz = Rcpp::as<std::vector<float> > (zz1);
  vector<vector<float>> markovTableC;
  for (unsigned i = 0; i < colnames.size(); i++ ){
    zz1 = markovTable( i, _); // rows
    zz = Rcpp::as<std::vector<float> > (zz1);
    markovTableC.push_back(zz);
  }

  int IDsize;
  if(TYPEOF(ID)==STRSXP){
    vector<string> IDtemp = Rcpp::as<std::vector<string> > (ID);
    IDsize = IDtemp.size();
  }
  else if(TYPEOF(ID)==INTSXP) {
    vector<int> IDtemp = Rcpp::as<std::vector<int> > (ID);
    IDsize = IDtemp.size();
  }
  else {
    Rcpp::Rcerr << "Input ID must be a vector of Type character or int.\n";
    return NULL;
  }
  if (data.nrows() != IDsize){
    Rcpp::Rcerr << " ID-Vector and Input-Dataframe must have the same size. " << endl;
    return NULL;
  }

  vector<string>  varsV(data.nrows());
  vector<float>  varsVFloat(data.nrows());
  vector<int> varsVInt(data.nrows());
  vector<string> inputVector(data.size());
  vector<vector<string> > varsM ;
  CharacterVector CLKs(data.nrows());
  vector<int> padding_ =Rcpp::as<std::vector<int> > (padding);
  vector<int> qgram_ =Rcpp::as<std::vector<int> > (qgram);
  //check password vector
  vector<string> password_;
  if(TYPEOF(password)!=STRSXP){
    Rcpp::Rcerr << "Please select a password for each input variable in a vector of class character.";
  }
  else{
    password_ =  Rcpp::as<std::vector<string> > (password);
    if (data.size() != (unsigned)password_.size()){
      Rcpp::Rcerr << "vector of password must have the same size as the input data.frame." <<data.nrows() << " "<< password_.size() ;
      return NULL;}
  }
  //Sizes of Vectors are checked and adapted
  if (data.size() > (unsigned)padding_.size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be fill with zeros."<< endl;
    for (int i = padding_.size() ; i < data.size(); i++){
      padding_.push_back(0);
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
  if (data.size() < (unsigned)padding_.size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be cut."<< endl;
    while((unsigned)padding_.size() > data.size()){
      padding_.pop_back();
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }

  if (data.size() > (unsigned)qgram_.size()){
    Rcpp::Rcerr << "Vector qgrams must have the same size as the input data.frame. Qgrams will be fill with 2s." << endl;
    for (int i = qgram_.size() ; i< data.size(); i++){
      qgram_.push_back(2);
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
  if (data.size() < (unsigned)qgram_.size()){
    Rcpp::Rcerr << "Vector qgram must have the same size as the input data.frame. Qgram will be cut."<< endl;
    while((unsigned)qgram_.size() > data.size()){
      qgram_.pop_back();
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
  //Matrix is filled with input, Types are set to string
  for (int i = 0 ; i< data.size(); i++){
    if(TYPEOF(data[i])==STRSXP) {

      varsV=Rcpp::as<std::vector<string> > (data[i]);
      varsM.push_back(varsV);
    }
    if(TYPEOF(data[i])==REALSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains floats. Data will be transformed to characters." << endl;
      varsVFloat = Rcpp::as<std::vector<float> > (data[i]);

      if(varsVFloat.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          varsV[j] = to_string(varsVFloat[j]);//Copy the vector to the string
        }
      }
      varsM.push_back(varsV);
    }
    if(TYPEOF(data[i])==INTSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains integers or factors. Make sure the input data are not factors when you want to use characters. Data will be transformed to characters." << endl;
      varsVInt = Rcpp::as<std::vector<int> > (data[i]);
      if(varsVInt.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          varsV[j]=to_string(varsVInt[j]);//Copy the vector to the string
        }
      }
      varsM.push_back(varsV);
    }
  }

  //for each row
  for (int i = 0 ; i<data.nrows(); i++){
    for (int j = 0 ; j<data.size(); j++){
      inputVector[j]=varsM[j][i];
    }
    CLKs[i]=CreateMarkovCLKc(inputVector, k1, k2, padding_, qgram_, lenBloom, password_,  markovTableC, rownames, colnames, includeOriginalBigram, v);
  }


  return DataFrame::create(Named("ID") = ID ,
                           Named("CLKs") = CLKs,
                           Named("stringsAsFactors") = false);
}


/**
* Wrapper function for CreateEnsembleCLK
*/
// [[Rcpp::export]]
DataFrame CreateEnsembleCLK(SEXP ID, DataFrame data, SEXP password, int NumberOfCLK =1 , int k =20, IntegerVector padding =IntegerVector::create(0), IntegerVector qgram= IntegerVector::create(2), unsigned lenBloom =1000) {
  int IDsize;
  if(TYPEOF(ID)==STRSXP){
    vector<string> IDtemp = Rcpp::as<std::vector<string> > (ID);
    IDsize = IDtemp.size();
  }
  else if(TYPEOF(ID)==INTSXP) {
    vector<int> IDtemp = Rcpp::as<std::vector<int> > (ID);
    IDsize = IDtemp.size();
  }
  else {
    Rcpp::Rcerr << "Input ID must be a vector of Type character or int.\n";
    return NULL;
  }
  //data = DataFrame::create(data,_["stringsAsFactors"] = false );
  if (data.nrows() != IDsize){
    Rcpp::Rcerr << " ID-Vector and Input-Dataframe must have the same size. " << endl;
    return NULL;
  }

  vector<string>  varsV(data.nrows());
  vector<float>  varsVFloat(data.nrows());
  vector<int> varsVInt(data.nrows());
  vector<string> inputVector(data.size());
  vector<vector<string> > varsM ;
  CharacterVector CLKs(data.nrows());
  vector<int> padding_ =Rcpp::as<std::vector<int> > (padding);
  vector<int> qgram_ =Rcpp::as<std::vector<int> > (qgram);

  //check password vector
  vector<string> password_;
  if(TYPEOF(password)!=STRSXP){
    Rcpp::Rcerr << "Please select a password for each input variable in a vector of class character.";
  }
  else{
    password_ =  Rcpp::as<std::vector<string> > (password);
    if (data.size() != (unsigned)password_.size()){
      Rcpp::Rcerr << "vector of password must have the same size as the input data.frame." <<data.nrows() << " "<< password_.size() ;
      return NULL;}
  }
  //Sizes of Vectors are checked and adapted
  if (data.size() > (unsigned)padding_.size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be fill with zeros."<< endl;
    for (int i = padding_.size() ; i< data.size(); i++){
      padding_.push_back(0);
    }
  }
  if (data.size() < (unsigned)padding_.size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be cut."<< endl;
    while((unsigned)padding_.size() > data.size()){
      padding_.pop_back();
    }
  }

  if (data.size() > (unsigned)qgram_.size()){
    Rcpp::Rcerr << "Vector qgrams must have the same size as the input data.frame. Qgrams will be fill with 2s." << endl;
    for (int i = qgram_.size() ; i< data.size(); i++){
      qgram_.push_back(2);
    }
  }
  if (data.size() < (unsigned)qgram_.size()){
    Rcpp::Rcerr << "Vector qgram must have the same size as the input data.frame. Qgram will be cut."<< endl;
    while((unsigned)qgram_.size() > data.size()){
      qgram_.pop_back();
    }
  }
  //Matrix is filled with input, Types are set to string
  for (int i = 0 ; i < data.size(); i++){
    if(TYPEOF(data[i])==STRSXP) {

      varsV=Rcpp::as<std::vector<string> > (data[i]);
      varsM.push_back(varsV);
    }
    if(TYPEOF(data[i])==REALSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains floats. Data will be transformed to characters." << endl;
      varsVFloat = Rcpp::as<std::vector<float> > (data[i]);

      if(varsVFloat.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          varsV[j] = to_string(varsVFloat[j]);//Copy the vector to the string
        }
      }
      varsM.push_back(varsV);
    }
    if(TYPEOF(data[i])==INTSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains integers or factors. Make sure the input data are not factors when you want to use characters. Data will be transformed to characters." << endl;
      varsVInt = Rcpp::as<std::vector<int> > (data[i]);
      if(varsVInt.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          varsV[j]=to_string(varsVInt[j]);//Copy the vector to the string
        }
      }
      varsM.push_back(varsV);
    }
  }

  //for each row
  for (int i = 0 ; i<data.nrows(); i++){
    for (int j = 0 ; j<data.size(); j++){
      inputVector[j]=varsM[j][i];
    }
    CLKs[i]=CreateEnsembleCLKc(inputVector, k, NumberOfCLK, padding_, qgram_, lenBloom, password_ );
  }

  return DataFrame::create(Named("ID") = ID ,
                           Named("CLK") = CLKs,
                           Named("stringsAsFactors") = false);
}


/**
* Wrapper function for CreateBF with just one vector as input
*/
//[[Rcpp::export]]
DataFrame CreateBF(SEXP ID, SEXP data,SEXP password, int k =20, int padding =1, int qgram= 2, unsigned lenBloom =1000 ) {
  DataFrame res = CreateBFc(ID, data, password, k, padding, qgram, lenBloom );
  return res;
}

DataFrame CreateBFc(SEXP ID, SEXP data,SEXP password, int k =20, int padding =1, int qgram= 2, unsigned lenBloom =1000 ) {
  if(TYPEOF(data)!=STRSXP){
    Rcpp::Rcerr << "Input data must be a vector of Type character.\n" << TYPEOF(data);
    return 0;
  }
  int IDsize;
  if(TYPEOF(ID)==STRSXP){
    vector<string> IDtemp = Rcpp::as<std::vector<string> > (ID);
    IDsize = IDtemp.size();
  }
  else if(TYPEOF(ID)==INTSXP) {
    vector<int> IDtemp = Rcpp::as<std::vector<int> > (ID);
    IDsize = IDtemp.size();

  }
  else {
    Rcpp::Rcerr << "Input ID must be a vector of Type character or int.\n";
    return 0;
  }

  vector<string> inputVector = Rcpp::as<std::vector<string> > (data);
  CharacterVector CLKs(inputVector.size()) ;
  string password_ =  as<string>(password);

  int qgram_ = qgram;
  int padding_ = padding;
  if (inputVector.size() != (unsigned)IDsize){
    Rcpp::Rcerr << " ID-Vector and Input-vector must have the same size. " <<inputVector.size() << " " << IDsize << endl;
    return 0;
  }


  //for each row

  for (unsigned i = 0 ; i < inputVector.size(); i++){
    CLKs[i]= CreateBFBigramSeed(inputVector[i], k, padding_, qgram_, lenBloom, password_);

  }
  return DataFrame::create(Named("ID") = ID,
                           Named("CLKs") = CLKs,
                           Named("stringsAsFactors") = false);
}


/**
* Wrapper function for CreateRecordLevelBF
*/
//[[Rcpp::export]]
DataFrame CreateRecordLevelBF(SEXP ID ,DataFrame data, SEXP password, int lenRLBF =1000,  int k =20,
                              IntegerVector padding = IntegerVector::create(0), IntegerVector qgram= IntegerVector::create(2),
                              IntegerVector lenBloom = IntegerVector::create(500) , std::string method = "StaticUniform", NumericVector weigths =NumericVector::create(1)) {

  int IDsize, n = 0;
  if(TYPEOF(ID)==STRSXP){
    vector<string> IDtemp = Rcpp::as<std::vector<string> > (ID);
    IDsize = IDtemp.size();
  }
  else if(TYPEOF(ID)==INTSXP) {
    vector<int> IDtemp = Rcpp::as<std::vector<int> > (ID);
    IDsize = IDtemp.size();
  }
  else {
    Rcpp::Rcerr << "Input ID must be a vector of Type character or int.\n";
    return 0;
  }
  if (data.nrows()!=IDsize){
    Rcpp::Rcerr << " ID-Vector and Input-Dataframe must have the same size. " << endl;
    return 0;
  }
  vector<string>  varsV(data.nrows());
  vector<float>  varsVFloat(data.nrows());
  vector<int> varsVInt(data.nrows());
  vector<string> inputVector(data.size());
  vector<vector<string> > varsM ;

  vector<int> padding_ =Rcpp::as<std::vector<int> > (padding);
  vector<int> qgram_ =Rcpp::as<std::vector<int> > (qgram);
  vector<int> lenBloomV =Rcpp::as<std::vector<int> > (lenBloom);
  unsigned lenBloom_ = lenBloomV[0];

  //check password vector
  vector<string> password_;
  if(TYPEOF(password)!=STRSXP){
    Rcpp::Rcerr << "Please select a password for each input variable in a vector of class character.";
  }
  else{
    password_ =  Rcpp::as<std::vector<string> > (password);
    if (data.size() != (unsigned)password_.size()){
      Rcpp::Rcerr << "vector of passwords must have the same size as the input data.frame."  ;
      return NULL;}
  }
  vector<float> weigths_ = Rcpp::as<std::vector<float> > (weigths);
  CLK* rlbf = new CLK(lenRLBF);
  DataFrame RLBFout;
  vector<string> RLBF;
  //LogicalMatrix RLBFres(data.nrows(), lenRLBF);
  char *CLKout = new char[lenRLBF+1] ;

  if ((unsigned int)data.size()>padding_.size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be fill with zeros."<< endl;
    for (int i = padding_.size() ; i< data.size(); i++){
      padding_.push_back(0);
    }
  }
  if ((unsigned int)data.size()<padding_.size()){
    Rcpp::Rcerr << "Vector padding must have the same size as the input data.frame. Padding will be cut."<< endl;
    while((unsigned int)padding_.size()> (unsigned int)data.size()){
      padding_.pop_back();
    }
  }
  if ((unsigned int)data.size()>qgram_.size()){
    Rcpp::Rcerr << "Vector qgrams must have the same size as the input data.frame. Qgrams will be fill with 2s." << endl;
    for (int i = qgram_.size() ; i< data.size(); i++){
      qgram_.push_back(2);
    }
    //return NULL; Alternative zu Auffüllen ist Programm abzubrechen.
  }
  if ((unsigned int)data.size()<qgram_.size()){
    Rcpp::Rcerr << "Vector qgram must have the same size as the input data.frame. Qgram will be cut."<< endl;
    while((unsigned int)qgram_.size()> (unsigned int)data.size()){
      qgram_.pop_back();
    }
  }

  //Matrix is filled with input, Types are set to string
  for (int i = 0 ; i< data.size(); i++){
    if(TYPEOF(data[i])==STRSXP) {
      varsV=Rcpp::as<std::vector<string> > (data[i]);
      varsM.push_back(varsV);
    }
    if(TYPEOF(data[i])==REALSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains floats. Data will be transformed to characters." << endl;
      varsVFloat = Rcpp::as<std::vector<float> > (data[i]);

      if(varsVFloat.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          varsV[j] = to_string(varsVFloat[j]);//Copy the vector to the string
        }
      }
      varsM.push_back(varsV);
    }
    if(TYPEOF(data[i])==INTSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains integers or factors. Make sure the input data are not factors when you want to use characters. Data will be transformed to characters." << endl;
      varsVInt = Rcpp::as<std::vector<int> > (data[i]);
      if(varsVInt.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          varsV[j]=to_string(varsVInt[j]);//Copy the vector to the string
        }
      }
      varsM.push_back(varsV);
    }
  }
  if (method != "StaticUniform" &&method != "StaticWeigthed" && method != "DynamicUniform" && method != "DynamicWeigthed"){
    Rcpp::Rcerr << "'method' can't be '" <<method <<"', it must be eihter 'StaticUniform', 'StaticWeigthed', 'DynamicUniform' or 'DynamicWeigthed'\n";
    delete rlbf;
    delete[] CLKout;
    return RLBFout;
  }
  if (method == "StaticWeigthed"||method == "DynamicWeigthed"){ // set length of the parts of the rlbf if it is weigthed
    if ((unsigned)weigths_.size() != data.size()){
      Rcpp::Rcerr << "Vector weights must have the same size as the input data.frame.";
      delete rlbf;
      delete[] CLKout;
      return RLBFout;
    }
    float sum = accumulate(weigths_.begin(), weigths_.end(), 0.0f);
    if (accumulate(weigths_.begin(), weigths_.end(), 0.0f) !=1.0){ // sum of the has to be 1
      Rcpp::Rcerr << "Sum of weights has to be 1.0. but is " << sum << "\n";
      delete rlbf;
      delete[] CLKout;
      return RLBFout;
    }
    //set sizes of the parts
    rlbf ->setPartsize(data.size(), &weigths_[0]);

  } else{ // set length of the parts of the rlbf if it is not weigthed
    weigths_.clear();
    for (int i  = 0 ; i<data.size(); i++){
      if (data.size()>0){
        weigths_.push_back((float)1/data.size());
      }
    }
    rlbf ->setPartsize(data.size(), &weigths_[0]);
  }

  if (method == "StaticUniform"||method == "StaticWeigthed"){ // set length of Bloom Filter if "Static"
    if (lenBloom.size()>1){
      Rcpp::Rcerr << "Length of the vector of length of Bloom filters is bigger than 1, all BF have to have the same size, the first value will be taken.";
    }
    for (int i  = 0 ; i<data.size(); i++){
      n = i;
    }
  } else { // set length of Bloom Filter if "Dynamic"
    if (data.size() == (unsigned)lenBloomV.size()){
      for (int i  = 0 ; i<data.size(); i++){
        n = i;
      }
    } else {
      Rcpp::Rcerr << "Length of the vector of length of Bloom filters has to be the same as the number of variables.";
      delete rlbf;

      delete[] CLKout;
      return RLBFout;
    }
  }
  for (int i = 0 ; i<data.nrows(); i++){
    for (int j = 0 ; j<data.size(); j++){
      inputVector[j]=varsM[j][i];
    }
    rlbf->clear();
    CreateRBF(rlbf, inputVector, k, padding_, qgram_, lenBloom_, lenRLBF, password_);
    rlbf->copyToString(CLKout, lenRLBF);
    RLBF.push_back(CLKout);
  }

  delete rlbf;
  delete[] CLKout;
  return DataFrame::create(Named("ID") = ID ,
                           Named("RLBF") = RLBF,
                           Named("stringsAsFactors") = false);
}


