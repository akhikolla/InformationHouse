#include "cpphutils.h"


// #nocov start
// [[Rcpp::export]]
void showValue(const char* what, double x) {
  Rcout << " " << what << " \t " << x << std::endl;
  // return 0;
}
// #nocov end


double sinhalfsq (double x) {
  const double o = sin(x / 2);
  return o * o;
  // return 0.5 - 0.5 * cos(x); // slowest!
  // return (sin(x/2) * sin(x/2));
}

// [[Rcpp::export]]
double haversine_distance (double olat1, double olon1, double olat2, double olon2,
                           bool unitless = false) {
  // double pi = 3.1415926535897;
  const double lat1 = olat1 * (M_PI / 180) ;
  const double lat2 = olat2 * (M_PI / 180) ;
  const double lon1 = olon1 * (M_PI / 180) ;
  const double lon2 = olon2 * (M_PI / 180) ;
  // std::abs doesn't work with the precision req


  // const double delta_lat = (lat1 > lat2) ? (lat1 - lat2) : (lat2 - lat1) ;
  const double delta_lat = std::fabs(lat1 - lat2);
  // const double delta_lon = (lon1 > lon2) ? (lon1 - lon2) : (lon2 - lon1) ;
  const double delta_lon = std::fabs(lon1 - lon2);

  // 6371 * 2 * asin(sqrt(sin(d_lat / 2)^2 + cos(lat1) * cos(lat2) * sin(d_lon / 2)^2))
  double out = 0;
  if (unitless) {
    // just assume sin(x) = x
    out = delta_lat * delta_lat + delta_lon * delta_lon * cos(lat1) * cos(lat2);
    return out;
  }

  double den = cos(lat1) * cos(lat2) * sinhalfsq(delta_lon);
  out = sinhalfsq(delta_lat);
  out += den;
  out = sqrt(out);
  out = asin(out);
  out *= 6371;
  out *= 2;
  return out;
}

double do_euclid_dist (double x1,
                       double x2,
                       double y1,
                       double y2,
                       bool unitless = false) {
  double out = 0;
  double x = x2 - x1;
  double y = y2 - y1;
  out = x * x + y * y;
  if (unitless) {
    return out;
  }
  return sqrt(out);
}

// unitless =>

// [[Rcpp::export]]
NumericVector haversineDistance(NumericVector lat1, NumericVector lon1, NumericVector lat2, NumericVector lon2, bool unitless = false) {
  R_xlen_t N = lat1.length();
  if (N != lon1.length() || N != lat2.length() || N != lon2.length()) {
    stop("Lengths of input vectors differ.");
  }
  NumericVector out(N);
  if (unitless) {
    for (R_xlen_t i = 0; i < N; ++i) {
      out[i] = haversine_distance(lat1[i], lon1[i], lat2[i], lon2[i], true);
    }
  } else {
    for (R_xlen_t i = 0; i < N; ++i) {
      out[i] = haversine_distance(lat1[i], lon1[i], lat2[i], lon2[i]);
    }
  }
  return out;
}





// [[Rcpp::export]]
NumericVector theEuclidDistance (NumericVector x1,
                                 NumericVector x2,
                                 NumericVector y1,
                                 NumericVector y2,
                                 bool unitless = false) {
  R_xlen_t N = x1.size();
  R_xlen_t y1N = y1.size();
  R_xlen_t x2N = x2.size();
  R_xlen_t y2N = y2.size();
  if (N != y1N || N != x2N || N != y2N) {
    stop("x and y lengths differ.");
  }
  NumericVector out(N);
  double x1i = 0;
  double x2i = 0;
  double y1i = 0;
  double y2i = 0;
  for (R_xlen_t i = 0; i < N; ++i) {
    x1i = x1[i];
    x2i = x2[i];
    y1i = y1[i];
    y2i = y2[i];
    out[i] = do_euclid_dist(x1i, x2i, y1i, y2i, unitless);
  }
  return out;
}

// [[Rcpp::export]]
double hausdorffEuclid (NumericVector x,
                        NumericVector y) {
  R_xlen_t N = x.size();
  double maxmin_dist = 0;
  for (R_xlen_t i = 0; i < N; ++i) {
    double xi = x[i];
    double yi = y[i];
    double min_dist_i = 0;
    for (int j = 0; j < N; ++j) {
      if (i == j) {
        continue;
      }
      double xj = x[j];
      double yj = y[j];
      double dij = 0;
      dij = do_euclid_dist(xi, xj, yi, yj);
      if (min_dist_i == 0 || dij < min_dist_i) {
        min_dist_i = dij;
      }
    }
    if (min_dist_i > maxmin_dist) {
      maxmin_dist = min_dist_i;
    }
  }
  return maxmin_dist;
}

// 00 -> bottom left
// 01 -> top left
// 10 -> bottom right
// 11 -> top right
// [[Rcpp::export]]
IntegerVector EmptiestQuarter (NumericVector x,
                               NumericVector y,
                               double minx = 1,
                               double maxx = -1,
                               double miny = 1,
                               double maxy = -1) {
  NumericVector x_range(4);
  x_range[0] = minx, x_range[1] = maxx;
  if (minx > maxx) {
    x_range = NumericVector(do_range_dbl(x));
    minx = x_range[0];
    maxx = x_range[1];
  }
  NumericVector y_range(4);
  y_range[0] = miny, y_range[1] = maxy;
  if (miny > maxy) {
    y_range = NumericVector(do_range_dbl(y));
    miny = y_range[0];
    maxy = y_range[1];
  }


  double xcentre = minx + (maxx - minx) / 2;
  double ycentre = miny + (maxy - miny) / 2;
  R_xlen_t N = x.size();
  int q = -1;
  int o = -1;

  int min_points = N;
  double x0, x1, y0, y1 = 0;
  for (q = 0; q < 4; ++q) {

    int q_points = 0;
    switch (q) {
    case 0:
      x0 = minx;
      x1 = xcentre;
      y0 = miny;
      y1 = ycentre;
      break;
    case 1:
      x0 = minx;
      x1 = xcentre;
      y0 = ycentre;
      y1 = maxy;
      break;
    case 2:
      x0 = xcentre;
      x1 = maxx;
      y0 = miny;
      y1 = ycentre;
      break;
    case 3:
      x0 = xcentre;
      x1 = maxx;
      y0 = ycentre;
      y1 = maxy;
      break;
    }
    for (R_xlen_t i = 0; i < N; ++i) {
      double xi = x[i];
      double yi = y[i];
      if (xi > x0 && xi < x1 && yi > y0 && yi < y1) {
        ++q_points;
      }
    }
    if (q_points < min_points) {
      min_points = q_points;
      o = q;
    }
  }
  IntegerVector out(2);
  out[0] = o;
  out[1] = min_points;
  return out;
}

// [[Rcpp::export]]
IntegerVector theEmptiestQuarters (NumericVector x,
                                   NumericVector y,
                                   double minx = 1,
                                   double maxx = -1,
                                   double miny = 1,
                                   double maxy = -1,
                                   int depth = 4) {
  NumericVector x_range(4);
  x_range[0] = minx, x_range[1] = maxx;
  if (minx > maxx) {
    x_range = NumericVector(do_range_dbl(x));
    minx = x_range[0];
    maxx = x_range[1];
  }
  NumericVector y_range(4);
  y_range[0] = miny, y_range[1] = maxy;
  if (miny > maxy) {
    y_range = NumericVector(do_range_dbl(y));
    miny = y_range[0];
    maxy = y_range[1];
  }
  IntegerVector out(depth);
  double x0 =0;
  double x1 = 0;
  double y0 = 0;
  double y1 = 0;
  x0 += minx;
  x1 += maxx;
  y0 += miny;
  y1 += maxy;
  double dx, dy = 0;
  int o = EmptiestQuarter(x, y, minx, maxx, miny, maxy)[0];
  for (R_xlen_t i = 0; i < depth; ++i) {
    out[i] = o;
    dx = x1 - x0;
    dx /= 2;
    dy = y1 - y0;
    dy /= 2;
    switch (o) {
    case 0:
      x1 -= dx;
      y1 -= dy;
      break;
    case 1:
      x1 -= dx;
      y0 += dy;
      break;
    case 2:
      x0 += dx;
      y1 -= dy;
      break;
    case 3:
      x0 += dx;
      y0 += dy;
      break;
    }
    IntegerVector eq = EmptiestQuarter(x, y, x0, x1, y0, y1);
    o = eq[0];
    int nPoints = eq[1];
    if (nPoints == 0) {
      for (int j = i + 1; j < depth; ++j) {
        out[j] = -1;
      }
      break;
    }
  }

  return out;
}

// IntegerVector FirstEmptyQtr (NumericVector x,
//                              NumericVector y,
//                              double minx = 1,
//                              double maxx = -1,
//                              double miny = 1,
//                              double maxy = -1,
//                              int depth = 4) {
//   for (int xp = 1; xp <= depth; ++xp) {
//     int two_xp = pow(2, xp);
//     for ()
//   }
// }

// [[Rcpp::export]]
int which_min_HaversineDistance (NumericVector lat1,
                                 NumericVector lon1,
                                 double lat2,
                                 double lon2,
                                 double upperBound = 10) {
  R_xlen_t N = lat1.length();
  if (N != lon1.length()) {
    stop("length(lat1) != length(lat2).");
  }

  // Find a way to detect upper bound at runtime
  double max_delta_lat = 0;
  double max_delta_lon = 0;
  double new_dist = 0;

  if (upperBound > 0) {
    while (max_delta_lat < 2 && new_dist < upperBound) {
      max_delta_lat += 0.001;
      new_dist = haversine_distance(lat2, lon2, lat2 + max_delta_lat, lon2);
    }
    new_dist = 0;
    while (max_delta_lon < 2 && new_dist < upperBound) {
      max_delta_lon += 0.001;
      new_dist = haversine_distance(lat2, lon2, lat2, lon2 + max_delta_lon);
    }
  }

  int out = 1;
  double cur_dist = 0;
  double min_dist = 50000; // bigger than the circumference of Earth
  min_dist = haversine_distance(lat1[0], lon1[0], lat2, lon2);
  double delta_lat = 0;
  double delta_lon = 0;
  bool skip = false;

  for (R_xlen_t i = 1; i < N; ++i) {
    // If the delta lat is greater than the upper bound, don't calculate haversine distance
    if (upperBound > 0) {
      delta_lat = lat1[i] - lat2;
      skip = (delta_lat > 0) ? (delta_lat > max_delta_lat) : (delta_lat < -max_delta_lat);
      if (skip) {
        continue ;
      }
      delta_lon = lon1[i] - lon2;
      skip = (delta_lon > 0) ? (delta_lon > max_delta_lon) : (delta_lon < -max_delta_lon);
      if (skip) {
        continue ;
      }
    }

    cur_dist = haversine_distance(lat1[i], lon1[i], lat2, lon2);
    if (cur_dist < min_dist) {
      min_dist = cur_dist;
      out = i + 1;
    }
  }
  return out;
}

// [[Rcpp::export]]
List match_min_Haversine (NumericVector lat1,
                          NumericVector lon1,
                          NumericVector lat2,
                          NumericVector lon2,
                          IntegerVector tabl,
                          double cartR = -1,
                          double dist0_km = 0.01,
                          bool verify_cartR = false,
                          bool do_verify_box = false,
                          bool excl_self = false,
                          int ncores = 1) {
  int N1 = lat1.length();
  if (N1 != lon1.length()) {
    stop("length(lat1) != length(lon1).");
  }
  int N2 = lat2.length();
  if (N2 != lon2.length()) {
    stop("length(lat2) != length(lon2)");
  }
  if (excl_self && N2 != N1) {
    warning("`excl_self = true`, yet lengths of `lat1` and `lat2` differ. This implies matching positions in lat1,lon1 and lat2,lon2 do not reflect the same points.");
  }

  int N3 = tabl.length();
  bool use_tbl = N3 == N2;

  double lati = 0;
  double loni = 0;
  double latj = 0;
  double lonj = 0;

  IntegerVector out(N1);
  NumericVector out2(N1);

  // half-equatorial circumference: used as an 'infinity' for
  // min_dist while also available to check we have actually
  // achieved a minimum distance. (Should be 1 and around 20,000.)
  double BIGDIST = haversine_distance(0, 0, 0, 179.99, true);
  double BIGDISTKM = haversine_distance(0, 0, 0, 179.99, false);

  bool do_verify_cartR = verify_cartR;
  bool do_check_cartR = cartR > 0;

  int k = 0;

  for (R_xlen_t i = 0; i < N1; ++i) {
    if (ncores == 1 && (i % 16) == 0) {
      Rcpp::checkUserInterrupt();
    }
    lati = lat1[i];
    loni = lon1[i];

    // Unitless distances to monitor for the minimum
    double min_dist = BIGDIST;
    double cur_dist = 0;

    // Use this to check the 'near enough' distance after a minimum candidate
    double min_dist_km = BIGDISTKM;
    double max_lati = lati + cartR;
    double min_lati = lati - cartR;
    double max_loni = loni + cartR;
    double min_loni = loni - cartR;

    k = 0;
    for (int j = 0; j < N2; ++j) {
      if (excl_self && j == i) {
        continue ;
      }
      if (excl_self && j > i) {
        // Since excl_self is only used when
        // the two sets are the same,
        // we can get more performance
        // by just replaying the values
        // we've already found to be closest.

        // Perhaps not? Consider x = {1, 2, 4, 8, 16}
      }

      latj = lat2[j];
      lonj = lon2[j];

      // Don't calculate haversine distance when the current
      // point is outside the box defined by cartR
      if (do_check_cartR) {
        if (latj < min_lati ||
            latj > max_lati ||
            lonj < min_loni ||
            lonj > max_loni) {
          continue;
        }
      }

      // unitless if we just need to compare to min_dist
      cur_dist = haversine_distance(lati, loni, latj, lonj, true);
      if (cur_dist < min_dist) {
        min_dist = cur_dist;

        // There was an issue with undefined behaviour;
        // k being out of range of j appeared to be the
        // cause so this was included in case k was somehow
        // not being reset within every j-loop.
        k = 0;
        k += j;

        // Calculate the real distance
        min_dist_km = haversine_distance(lati, loni, latj, lonj);
        // If within the "close-enough" distance, break to the next item to geocode.
        if (min_dist_km < dist0_km) {
          break;
        }
      }
    }

    if (do_verify_box) {
      // The half-length of the square to check within
      double box_r = do_euclid_dist(loni, lon2[k], lati, lat2[k]);
      double cur_dist_km_new = 0;

      // Box is around the *target*
      double box_max_lat = lati + box_r;
      double box_min_lat = lati - box_r;
      double box_max_lon = loni + box_r;
      double box_min_lon = loni - box_r;
      for (int j = 0; j < N2; ++j) {
        double lat2j = lat2[j], lon2j = lon2[j];
        if (lat2j > box_min_lat &&
            lat2j < box_max_lat &&
            lon2j > box_min_lon &&
            lon2j < box_max_lon) {
          cur_dist_km_new = haversine_distance(lati, loni, lat2j, lon2j);
          if (cur_dist_km_new < min_dist_km) {
            k = 0;
            k += j;
            min_dist_km = cur_dist_km_new;
          }
        }
      }
    }

    if (do_verify_cartR && min_dist_km == BIGDISTKM) {
      // We have failed to identify a small distance
      // Likely reason: too ambitious cartR
      do_check_cartR = false; // set from here onwards -- if it happens once, likely to happen again

      for (int j = 0; j < N2; ++j) {
        double cur_dist_km = 0;
        cur_dist_km = haversine_distance(lati, loni, latj, lonj, false);
        if (cur_dist_km < min_dist_km) {
          k = 0;
          k += j;
          min_dist_km = cur_dist_km;
          if (min_dist_km < dist0_km) {
            break;
          }
        }
      }
    }

    if (use_tbl) {
      if (k >= tabl.length()) {
        showValue("k", k); // #nocov
        showValue("i", i); // #nocov
        stop("k >= tabl.length"); // #nocov
      }
      out[i] = tabl[k];
      out2[i] = min_dist_km;
    } else {
      if (k >= lon2.length()) {
        showValue("k", k); // #nocov
        showValue("i", i); // #nocov
        stop("k >= tabl.length"); // #nocov
      }
      ++k;  // for R indexing
      out[i] = k;
      out2[i] = min_dist_km;
    }
  }

  List mat = List::create(Named("pos") = out,
                          Named("dist") = out2);
  return mat;
}





