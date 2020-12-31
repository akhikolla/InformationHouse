#include <Rcpp.h>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <math.h>
#include "image_helper.h"
using namespace Rcpp;

//' Power spectrum calculation
//'
//' Power spectrum of the difference of neighboring frames
//' in the series of speckle images
//'
//' @param filename a character string with the path name to a file.
//' @param threshold an integer (default is 50000).
//' @return The 513 x 1024 double matrix of power spectrum.
//' @examples
//' obj_filename <- system.file("extdata", "ads15182_550_2_frames.dat", package = "specklestar")
//' pow_spec_diff <- speckle_ps_diff(obj_filename)
//' @export
// [[Rcpp::export]]
NumericVector speckle_ps_diff(String filename, int threshold = 50000) {
  long frameSize = IMAGE_SIZE * sizeof(unsigned short);

  std::ifstream file(filename.get_cstring(), std::ios::binary);
  file.seekg(0, std::ios::end);
  long file_length = file.tellg();
  file.seekg(0, std::ios::beg);

  NumericMatrix outData(513, 1024);
  NumericMatrix big_dData(1024, 1024);
  int good_frames = 0;

  int N_frame = file_length / frameSize;
  unsigned short piData1[IMAGE_SIZE]
               , piData2[IMAGE_SIZE]
               ;
  memset(piData1, 0, frameSize);
  memset(piData2, 0, frameSize);

  fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 1024 * (1024 / 2 + 1));

  // int state = 0;
  bool state = true;
  int j = 0;
  while( file.read((char*)piData1, frameSize) ) {
    if( !file ) break;
    j++;
    if (IsOverThresholdFrame(piData1, threshold)) continue;
    good_frames += 1;
    break;
  }
  for(; file && j < N_frame; j++) {
    file.read((char*)(state ? piData1 : piData2), frameSize);
    if( !file ) break;
    if (IsOverThresholdFrame((state ? piData1 : piData2), threshold)) continue;

    for(int i = 0; i < 512; i++) {
      for(int j = 0; j < 512; j++) {
        big_dData[j * 1024 + i] = (double)(state ? piData1 : piData2)[i * 512 + j]
                                - (double)(state ? piData2 : piData1)[i * 512 + j];
      }
    }

    //fftw_plan p = fftw_plan_dft_r2c_2d(1024, 1024, big_dData.begin(), out, FFTW_ESTIMATE);
    fftw_plan p = fftw_plan_dft_r2c_2d(1024, 1024, big_dData.begin(), out, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    for (int i = 0; i < 1024 * 513; i++) outData[i] += out[i][0] * out[i][0] + out[i][1] * out[i][1];
    // state = 0 == state ? 1 : 0;
    state = !state;

    good_frames += 1;
  }
  fftw_free(out);
  file.close();

  Rcout << good_frames << " processed frames\n";

  return outData;
}
