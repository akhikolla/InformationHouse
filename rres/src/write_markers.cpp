#include <Rcpp.h>
#include <fstream>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
void write_markers_cpp(NumericVector marker, NumericVector freq, NumericMatrix genotype, CharacterVector memberID, std::string outfile) {
  std::ofstream markerfile;
  markerfile.open(outfile.c_str());
  markerfile << "map marker positions\n";
  for(int i = 0; i < marker.length(); i++){
    markerfile << marker[i] << "\n";
  }
  
  for(int i = 0; i < freq.length(); i++){
    markerfile << "set marker " << i+1 << " allele frequencies  " << freq[i] << " " << 1-freq[i] << "\n";
  }
  markerfile << "set markers " << freq.length() << " data\n";
  

  for(int i = 0; i < memberID.length(); i++){
    markerfile << memberID[i];
    for(int j = 0; j < genotype.ncol(); j++){
      markerfile << " " << genotype(i, j);
    }
    markerfile << "\n";
  }
  markerfile.close();
}
