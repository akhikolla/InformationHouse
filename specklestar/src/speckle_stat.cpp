#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "image_helper.h"
using namespace Rcpp;

//' Statistics of speckles
//'
//' Calculate statistics of speckles in the series of 512 x 512
//' speckle images and filter "bad" frames
//'
//' @param filename a character string with the path name to a file.
//' @param threshold an integer (default is 50000).
//' @return The list with 2 elements 'badFrames' and 'hist': \cr
//' 1 number of bad frames, \cr
//' 2 double vector of speckle statistics.
//' @examples
//' obj_filename <- system.file("extdata", "ads15182_550_2_frames.dat", package = "specklestar")
//' spec_stat <- speckle_stat(obj_filename)
//' @export
// [[Rcpp::export]]
List speckle_stat(String filename, int threshold = 50000) {
  std::ifstream file(filename.get_cstring(), std::ios::binary);
  file.seekg(0, std::ios::end);
  long file_length = file.tellg();
  int N_frame = file_length / (sizeof(unsigned short) * IMAGE_SIZE);
  file.seekg(0, std::ios::beg);

  unsigned short piData[IMAGE_SIZE];
  NumericVector badFrames;
  NumericVector histData(65535);

  for(int f = 0; f < N_frame; f++) {
    file.read((char*)piData, IMAGE_SIZE * sizeof(unsigned short));
    for(int i = 0; i < IMAGE_SIZE; i++) {
      histData[piData[i]]++;
    }

    if (IsOverThresholdFrame(piData, threshold))
        badFrames.push_back(f + 1);
  }

  file.close();

  return List::create(Named("badFrames") = badFrames,
                      Named("hist") = histData);
}
