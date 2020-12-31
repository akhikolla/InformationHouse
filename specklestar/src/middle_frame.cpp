#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include "image_helper.h"
using namespace Rcpp;

//' Middle frame
//'
//' Average image of the series of 512 x 512 px images
//'
//' @param filename A string.
//' @param subtrahend 512 x 512 matrix to subtract.
//' @param threshold An integer (default 50000).
//' @return The 512 x 512 matrix of middle speckle image.
//' @examples
//' obj_filename <- system.file("extdata", "ads15182_550_2_frames.dat", package = "specklestar")
//' zero_matrix <- matrix(0, 512, 512)
//' mf <- middle_frame(obj_filename, subtrahend = zero_matrix)
//' @export
// [[Rcpp::export]]
NumericVector middle_frame(String filename, NumericMatrix subtrahend, int threshold = 50000) {
  std::ifstream file(filename.get_cstring(), std::ios::binary);

  file.seekg(0, std::ios::end);
  long file_length = file.tellg();
  int N_frame = file_length / (sizeof(unsigned short) * IMAGE_SIZE);

  unsigned short piData[IMAGE_SIZE];
  memset(piData, 0, IMAGE_SIZE * sizeof(unsigned short));

  NumericMatrix meanData(512, 512);
  int n_good_frames = 0;

  file.seekg(0, std::ios::beg);
  for(int f = 0; f < N_frame; f++) {
    file.read((char*)piData, IMAGE_SIZE * sizeof(unsigned short));
    if (IsOverThresholdFrame(piData, threshold)) continue;
    n_good_frames++;

    for(int i = 0; i < IMAGE_SIZE; i++) {
      meanData[i] += piData[i] - subtrahend[i];
    }
  }

  for(int i = 0; i < IMAGE_SIZE; i++) {
    meanData[i] = meanData[i] / n_good_frames;
  }

  Rcout << n_good_frames << " averaged frames\n";

  file.close();
  return meanData;
}
