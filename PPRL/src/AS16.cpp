#include <Rcpp.h>
#include <random>
#include <string>

using namespace Rcpp;
using namespace std;

int computeHWNew(string s);
vector<int> computeHWs(vector<string>);
string CreateAS16cNew(string r, string password);
float CompareArmknechtCLKNew(string d1, string d2s);


/**
 * Function for Create
 *
 * This method is not necessary for compare, it is just to see the output of the method create
 *
 * @param id1 IDs of CLK1
 * @param CLK1 vector of CLKs to be encoded
 * @param id2 IDs of CLK2
 * @param CLK2 vector of CLKs to be encoded
 * @param password codeword to create vectors
 * @param t Tanomito similartiy to be tested
 */

// [[Rcpp::export]]
DataFrame CreateAS16(CharacterVector ID, CharacterVector data,  SEXP password) {
    //Algorithm 1 Procedure Create
    vector<int> hw;
    vector<string> d;
    vector<string> CLKs = Rcpp::as < vector < string > >(data);
    string passwordc =  as<string>(password);
    hw = computeHWs(CLKs);

    for (int i = 0; i < data.length(); i ++){
      d.push_back(CreateAS16cNew(CLKs[i], passwordc));
    }

    return DataFrame::create(
      Named("ID") = Rcpp::wrap(ID),
      Named("d") = Rcpp::wrap(d),
      Named("hw") = Rcpp::wrap(hw),
      Named("stringsAsFactors") = false);
    }


/**
 * Creates two vectors from CLKs according to Armknechts method create
 * and compares them to each other
 *
 * @param d1_ encoded record 1
 * @param h1 Hamming weights of the (unknown) record 1
 * @param d2_ encoded record 2
 * @param h2 Hamming weights of the (unknown) record 2
 * @param password codeword to create vectors
 * @param t Tanomito similartiy to be tested
 */
//[[Rcpp::export]]
DataFrame CompareAS16( CharacterVector IDA, CharacterVector dataA,  CharacterVector IDB, CharacterVector dataB, SEXP password, float t =0.85) {
  vector<string> CLKsA = Rcpp::as < vector < string > >(dataA);
  vector<string> CLKsB = Rcpp::as < vector < string > >(dataB);
  string passwordc =  as<string>(password);
  float x;
  string rA;
  CharacterVector ID1, ID2;
  NumericVector value1, value2, value3, value4 ;
  for(unsigned i = 0; i< CLKsA.size(); i++){
    rA = CreateAS16cNew(CLKsA[i], passwordc);
    for(unsigned j = 0; j< CLKsB.size(); j++){
      x = CompareArmknechtCLKNew(rA, CreateAS16cNew(CLKsB[j], passwordc));
      if (x<=((1-t)/(t+1))){
        ID1.push_back(IDA[i]);
        ID2.push_back(IDB[j]);
        value1.push_back(x);
        value2.push_back(((1-t)/(t+1)));
        value3.push_back(x <= ((1-t)/(t+1)));
        value4.push_back(((1-x)/(x+1)));
      }
    }
  }

   return DataFrame::create(Named("ID1") = Rcpp::wrap(ID1) ,
                          Named("ID2") = Rcpp::wrap(ID2),
                          //Named("ctr/h1+h2") = value1,
                          // Named("1-t/t+1") = value2,
                          Named("similarity") = Rcpp::wrap(value4),
                          // Named("Res") = (value3 == 1),
                          Named("stringsAsFactors") = false);
}


  //convert input from R to cpp
/**
 * Algorithm 1 Procedure Create
 *
 * @param clkR Record clk of length n
 */
string CreateAS16cNew(string r, string password){
  int n = r.size();
  //string w(n, '0');
  string d(n, '0');
  bool b;
  //Generate lenBloom/2 randoms bit b dependened on password as seed
  seed_seq seed(password.begin(), password.end());
  default_random_engine engine(seed);
  vector<int> randomBits(n/2);
  uniform_int_distribution<int> distr(0, n - 1);
  generate(begin(randomBits), end(randomBits), [&]() {return distr(engine);});
    for (int i = 1; i < n/2; i++){
      b = randomBits[i]%2;
       if (b){
         if (r[2*i-1] == '0'){
           d[2*i-1] = '1';
         }
         if (r[2*i] == '0'){
           d[2*i] = '1';
         }
          //w[2*i+1] = '1';
          //w[2*i] = '1';
       }else {
         if (r[2*i-1] == '1'){
           d[2*i-1] = '1';
         }
         if (r[2*i] == '1'){
           d[2*i] = '1';
         }
       }
  }
  return d;
}

float CompareArmknechtCLKNew(string d1, string d2) {
  float ctr=0; //counter, this value will track the estimate for HD(r1, r2).
  int h1 = computeHWNew(d1);
  int h2 = computeHWNew(d2);
  // Compute v = (v[1], . . . , v[n]) := d1 ⊕ d2. his yields the XOR of the two records plus a code word.
  //this is done within the if brackets

   for (int i =0; i < d1.size() ; i++){
     i++; //another i++, because two consecutive bit are compared as couples/pairs
     if(((int)d1[i] - '0' + (int)d2[i] - '0') % 2 != ((int)d1[i+1]-'0' + (int)d2[i+1]-'0') % 2){ // A difference between the two records has been found for sure.
       ctr++;
       }
   }
  return (ctr/(h1+h2));
}

int computeHWNew(string s){
  int hw = count(s.begin(), s.end(), '1');
  return hw;
}

vector<int> computeHWs(vector<string> data){
  vector<int> hw;
  for (int i = 0; i < data.size(); i++){
    hw.push_back(computeHWNew(data[i]));
  }
  return hw;
}
