#include <Rcpp.h>
using namespace Rcpp;

//' Speckle Generator
//'
//' Generate model 512 x 512 x 2 (bytes) speckle image of binary star
//'
//' @param rho a separation (an arcsec).
//' @param theta a positional angle.
//' @param dm a magnitude difference.
//' @param seeing a number.
//' @param speckle_sigma a number.
//' @param wind a wind speed.
//' @return The vector of model speckle image.
//' @examples
//' speckle_vector <- speckle_generator(rho = 0.5, theta = 70,
//' dm = 0.3, seeing = 20, speckle_sigma = 1, wind = 0)
//' speckle_matrix <- matrix(speckle_vector, nrow = 512, ncol = 512)
//' @export
// [[Rcpp::export]]
NumericVector speckle_generator(double rho, double theta, double dm, double seeing, double speckle_sigma, double wind) {
  int N_speckle = 300;
  int n_x = 512;
  int n_y = n_x;

  double m1 = 1000;
  double m2 = m1 / pow(2.512, dm);

  double rho_x = rho * (512 / 4.4) * cos((theta + 180) * M_PI / 180);
  double rho_y = rho * (512 / 4.4) * sin((theta + 180) * M_PI / 180);

  double stellar_center_x = R::rnorm(n_x / 2, wind);
  double stellar_center_y = R::rnorm(n_y / 2, wind);

  NumericVector speckle_field(262144);
  NumericVector gaussian_2d(262144);

  for (int i = 0; i < N_speckle; i++) {
    double x0 = R::rnorm(stellar_center_x, seeing) - rho_x / 2;
    double y0 = R::rnorm(stellar_center_y, seeing) - rho_y / 2;
    for (int x = 0;  x < n_x; x++) {
      for (int y = 0; y < n_y; y++) {
        gaussian_2d[x * n_y + y] = m1 * exp(-(pow((x - x0), 2)/(2 * pow(speckle_sigma, 2)) + pow((y - y0), 2)/(2 * pow(speckle_sigma, 2)))) + \
          m2 * exp(-(pow((x - x0 - rho_x), 2)/(2 * pow(speckle_sigma, 2)) + pow((y - y0 - rho_y), 2)/(2 * pow(speckle_sigma, 2))));
      }
    }
    for (int c = 0 ; c < 262144 ; c++) {
      speckle_field[c] = speckle_field[c] + gaussian_2d[c];
    }
  }
    return speckle_field;
}
