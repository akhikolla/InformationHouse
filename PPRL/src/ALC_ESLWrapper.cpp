// Enable C++11 via this plugin
// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include "ALC_ESLWrapper.h"

/**
 * Wrapper function for CreateALC
 */
// [[Rcpp::export]]
DataFrame CreateALC(CharacterVector ID ,DataFrame data, LogicalVector soundex, SEXP password) {
  vector<string>  dataV(data.nrows());
  vector<float>  dataVFloat(data.nrows());
  vector<int> dataVInt(data.nrows());
  vector<string> inputVector(data.size());
  vector<vector<string> > dataM ;
  CharacterVector ALC(data.nrows());
  vector<bool> soundex_ = Rcpp::as<std::vector<bool> > (soundex);
  string password_ =  as<string>(password);

  //Sizes of Vectors are checked and adapted
  if (data.size()!=soundex.size()){
    Rcpp::Rcerr << "Vector soundex must have the same size as the input data.frame. Please check both!"<< endl;
   return NULL;
  }


  //Matrix is filled with input, Types are set to string
  for (int i = 0 ; i< data.size(); i++){
    if(TYPEOF(data[i])==STRSXP) {
      dataV=Rcpp::as<std::vector<string> > (data[i]);
      dataM.push_back(dataV);
    }
    if(TYPEOF(data[i])==REALSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains floats. Data will be transformed to characters." << endl;
      dataVFloat = Rcpp::as<std::vector<float> > (data[i]);

      if(dataVFloat.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          dataV[j] = to_string(dataVFloat[j]);//Copy the vector to the string
        }
      }
      dataM.push_back(dataV);
    }
    if(TYPEOF(data[i])==INTSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains integers or factors. Make sure the input data are not factors when you want to use characters. Data will be transformed to characters." << endl;
      dataVInt = Rcpp::as<std::vector<int> > (data[i]);
      if(dataVInt.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          dataV[j]=to_string(dataVInt[j]);//Copy the vector to the string
        }
      }
      dataM.push_back(dataV);
    }
  }

  for (int i = 0 ; i<data.nrows(); i++){
    for (int j = 0 ; j<data.size(); j++){
      inputVector[j]=dataM[j][i];
    }
    ALC[i]=createALC(inputVector, soundex_, password_);
  }

  return DataFrame::create(Named("ID") = ID ,
                           Named("ALC") = ALC,
                           Named("stringsAsFactors") = false);
}

/**
 * Wrapper function for Create581
 */
// [[Rcpp::export]]
DataFrame Create581(CharacterVector ID ,DataFrame data, List code,  SEXP password) {
  vector<string>  dataV(data.nrows());
  vector<float>  dataVFloat(data.nrows());
  vector<int> dataVInt(data.nrows());
  vector<string> inputVector(data.size());
  vector<vector<string> > dataM ;
  CharacterVector ESL(data.nrows());
  vector<vector<int> > code_(code.size()) ;
  //Sizes of Vectors are checked and adapted
  if (data.size()!=code.size()){
    Rcpp::Rcerr << "Matrix must have the same size as the input data.frame. Please check both!"<< endl;
    return NULL;
  }


  Rcpp::List xlist(code);
  int n = xlist.size();
  for(int i=0; i<n; i++) {
    IntegerVector y(xlist[i]);
    code_[i] = Rcpp::as<std::vector<int> > (y);
   }
  string password_ =  as<string>(password);



  //Matrix is filled with input, Types are set to string
  for (int i = 0 ; i< data.size(); i++){
    if(TYPEOF(data[i])==STRSXP) {
      dataV=Rcpp::as<std::vector<string> > (data[i]);
      dataM.push_back(dataV);
    }
    if(TYPEOF(data[i])==REALSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains floats. Data will be transformed to characters." << endl;
      dataVFloat = Rcpp::as<std::vector<float> > (data[i]);

      if(dataVFloat.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          dataV[j] = to_string(dataVFloat[j]);//Copy the vector to the string
        }
      }
      dataM.push_back(dataV);
    }
    if(TYPEOF(data[i])==INTSXP) {
      Rcpp::Rcerr << "Warning: Column " << i << " contains integers or factors. Make sure the input data are not factors when you want to use characters. Data will be transformed to characters." << endl;
      dataVInt = Rcpp::as<std::vector<int> > (data[i]);
      if(dataVInt.size()>0)
      {
        for(int j=0; j<data.nrows(); j++){
          dataV[j]=to_string(dataVInt[j]);//Copy the vector to the string
        }
      }
      dataM.push_back(dataV);
    }
  }

  for (int i = 0 ; i<data.nrows(); i++){
    for (int j = 0 ; j<data.size(); j++){
      inputVector[j]=dataM[j][i];
    }
    ESL[i]=createESL(inputVector, code_, password_);
  }

  return DataFrame::create(Named("ID") = ID ,
                           Named("581") = ESL,
                           Named("stringsAsFactors") = false);
}

// /**
//  * Wrapper function for CreateESL
//  * Allgemeiner als Creaste581
//  */
// // [[Rcpp::export]]
// DataFrame CreateESL(CharacterVector ID ,DataFrame data, List code,  SEXP password) {
//   vector<string>  dataV(data.nrows());
//   vector<float>  dataVFloat(data.nrows());
//   vector<int> dataVInt(data.nrows());
//   vector<string> inputVector(data.size());
//   vector<vector<string> > dataM ;
//   CharacterVector ESL(data.nrows());
//   vector<vector<int> > code_(code.size()) ;
//   //Sizes of Vectors are checked and adapted
//   if (data.size()!=code.size()){
//     Rcpp::Rcerr << "Matrix must have the same size as the input data.frame. Please check both!"<< endl;
//     return NULL;
//   }
//
//
//   Rcpp::List xlist(code);
//   int n = xlist.size();
//   for(int i=0; i<n; i++) {
//     IntegerVector y(xlist[i]);
//     code_[i] = Rcpp::as<std::vector<int> > (y);
//   }
//   string password_ =  as<string>(password);
//
//
//
//   //Matrix is filled with input, Types are set to string
//   for (int i = 0 ; i< data.size(); i++){
//     if(TYPEOF(data[i])==STRSXP) {
//       dataV=Rcpp::as<std::vector<string> > (data[i]);
//       dataM.push_back(dataV);
//     }
//     if(TYPEOF(data[i])==REALSXP) {
//       Rcpp::Rcerr << "Warning: Column " << i << " contains floats. Data will be transformed to characters." << endl;
//       dataVFloat = Rcpp::as<std::vector<float> > (data[i]);
//
//       if(dataVFloat.size()>0)
//       {
//         for(int j=0; j<data.nrows(); j++){
//           dataV[j] = to_string(dataVFloat[j]);//Copy the vector to the string
//         }
//       }
//       dataM.push_back(dataV);
//     }
//     if(TYPEOF(data[i])==INTSXP) {
//       Rcpp::Rcerr << "Warning: Column " << i << " contains integers or factors. Make sure the input data are not factors when you want to use characters. Data will be transformed to characters." << endl;
//       dataVInt = Rcpp::as<std::vector<int> > (data[i]);
//       if(dataVInt.size()>0)
//       {
//         for(int j=0; j<data.nrows(); j++){
//           dataV[j]=to_string(dataVInt[j]);//Copy the vector to the string
//         }
//       }
//       dataM.push_back(dataV);
//     }
//   }
//
//   for (int i = 0 ; i<data.nrows(); i++){
//     for (int j = 0 ; j<data.size(); j++){
//       inputVector[j]=dataM[j][i];
//     }
//     ESL[i]=createESL(inputVector, code_, password_);
//   }
//
//   return DataFrame::create(Named("ID") = ID ,
//                            Named("ESL") = 581,
//                            Named("stringsAsFactors") = false);
// }
//
//
