#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
LogicalVector do_and3_na(LogicalVector x, LogicalVector y, LogicalVector z,
                         int nThread = 1,
                         int na_value = 0,
                         int maxCall = 3) {
  if (maxCall < 0) {
    stop("Internal error: maxCall exceeded.");
  }

  // length-0 = ignore
  // length-1 = treat as if it were recycled: TRUE => ignore; FALSE => all FALSE; NA => NA.
  R_xlen_t nx = x.length();
  R_xlen_t ny = y.length();
  R_xlen_t nz = z.length();

  R_xlen_t n = (nx >= ny && nx >= nz) ? nx : ((ny >= nz) ? ny : nz);
  if (n == 0) {
    return x;
  }
  if (n > nx) {
    if (ny >= nz) {
      return do_and3_na(y, x, z, nThread, na_value, --maxCall);
    } else {
      return do_and3_na(z, y, x, nThread, na_value, --maxCall);
    }
  }
  if (nz > ny) {
    return do_and3_na(x, z, y, nThread, na_value, --maxCall);
  }

  if ((nx > 1 && nx < n) ||
      (ny > 1 && ny < n) ||
      (nz > 1 && ny < n)) {
    stop("Internal error nx, ny, nz have incompatible/wrong lengths.");
  }

  LogicalVector out = no_init(n);

  bool na_becomes_true  = na_value > 0;
  bool na_becomes_false = na_value < 0;

  for (R_xlen_t i = 0; i < n; ++i) {
    R_xlen_t jx = (nx == n) ? i : 0;
    R_xlen_t jy = (ny == n) ? i : 0;
    R_xlen_t jz = (nz == n) ? i : 0;

    if (x[jx] == FALSE || y[jy] == FALSE || z[jz] == FALSE) {
      out[i] = FALSE;
      continue;
    }
    if (x[jx] == TRUE && y[jy] == TRUE && z[jz] == TRUE) {
      out[i] = TRUE;
      continue;
    }

    out[i] = (na_becomes_false) ? FALSE : (na_becomes_true ? TRUE : NA_LOGICAL);
  }
  return out;
}

inline int opn(bool eq, bool gt, bool lt) {
  // != == >= <=  >  <
  //  0  1  2  3  4  5
  return !(eq || gt || lt) ? 0 : (eq ? (gt ? 2 : (lt ? 3 : 1)) : (gt ? 4 : 5));
}


// [[Rcpp::export(rng = false)]]
IntegerVector do_which2_yr(IntegerVector Year,
                           int yr,
                           bool consider_yr,
                           IntegerVector x,
                           int xa,
                           bool eqx,
                           bool gtx,
                           bool ltx,
                           IntegerVector y,
                           int ya,
                           bool eqy,
                           bool gty,
                           bool lty) {
  int nx = x.size();
  if (y.size() != nx) {
    stop("Internal error: lengths differ.");
  }
  if ((gtx && ltx) ||
      (gty && lty)) {
    stop("Internal error: both > and < requested.");
  }
  std::vector<int> o;
  o.reserve(nx);
  const int opx = opn(eqx, gtx, ltx);
  const int opy = opn(eqy, gty, lty);
  const int opxy = opx + 10*opy;
  if (consider_yr && Year.size() != nx) {
    stop("Internal error: lengths differ (Year).");
  }

  for (int i = 0; i < nx; i++) {
    bool yr_bool = true;
    if (consider_yr) {
      yr_bool = Year[i] == yr;
    }
    switch (opxy) {
    case 0:
      if (yr_bool && x[i] != xa && y[i] != ya) o.push_back(i + 1);
      continue;
    case 1:
      if (yr_bool && x[i] == xa && y[i] != ya) o.push_back(i + 1);
      continue;
    case 2:
      if (yr_bool && x[i] >= xa && y[i] != ya) o.push_back(i + 1);
      continue;
    case 3:
      if (yr_bool && x[i] <= xa && y[i] != ya) o.push_back(i + 1);
      continue;
    case 4:
      if (yr_bool && x[i] >  xa && y[i] != ya) o.push_back(i + 1);
      continue;
    case 5:
      if (yr_bool && x[i] <  xa && y[i] != ya) o.push_back(i + 1);
      continue;

    case 10:
      if (yr_bool && x[i] != xa && y[i] == ya) o.push_back(i + 1);
      continue;
    case 11:
      if (yr_bool && x[i] == xa && y[i] == ya) o.push_back(i + 1);
      continue;
    case 12:
      if (yr_bool && x[i] >= xa && y[i] == ya) o.push_back(i + 1);
      continue;
    case 13:
      if (yr_bool && x[i] <= xa && y[i] == ya) o.push_back(i + 1);
      continue;
    case 14:
      if (yr_bool && x[i] >  xa && y[i] == ya) o.push_back(i + 1);
      continue;
    case 15:
      if (yr_bool && x[i] <  xa && y[i] == ya) o.push_back(i + 1);
      continue;

    case 20:
      if (yr_bool && x[i] != xa && y[i] >= ya) o.push_back(i + 1);
      continue;
    case 21:
      if (yr_bool && x[i] == xa && y[i] >= ya) o.push_back(i + 1);
      continue;
    case 22:
      if (yr_bool && x[i] >= xa && y[i] >= ya) o.push_back(i + 1);
      continue;
    case 23:
      if (yr_bool && x[i] <= xa && y[i] >= ya) o.push_back(i + 1);
      continue;
    case 24:
      if (yr_bool && x[i] >  xa && y[i] >= ya) o.push_back(i + 1);
      continue;
    case 25:
      if (yr_bool && x[i] <  xa && y[i] >= ya) o.push_back(i + 1);
      continue;

    case 30:
      if (yr_bool && x[i] != xa && y[i] <= ya) o.push_back(i + 1);
      continue;
    case 31:
      if (yr_bool && x[i] == xa && y[i] <= ya) o.push_back(i + 1);
      continue;
    case 32:
      if (yr_bool && x[i] >= xa && y[i] <= ya) o.push_back(i + 1);
      continue;
    case 33:
      if (yr_bool && x[i] <= xa && y[i] <= ya) o.push_back(i + 1);
      continue;
    case 34:
      if (yr_bool && x[i] >  xa && y[i] <= ya) o.push_back(i + 1);
      continue;
    case 35:
      if (yr_bool && x[i] <  xa && y[i] <= ya) o.push_back(i + 1);
      continue;

    case 40:
      if (yr_bool && x[i] != xa && y[i] > ya) o.push_back(i + 1);
      continue;
    case 41:
      if (yr_bool && x[i] == xa && y[i] > ya) o.push_back(i + 1);
      continue;
    case 42:
      if (yr_bool && x[i] >= xa && y[i] > ya) o.push_back(i + 1);
      continue;
    case 43:
      if (yr_bool && x[i] <= xa && y[i] > ya) o.push_back(i + 1);
      continue;
    case 44:
      if (yr_bool && x[i] >  xa && y[i] > ya) o.push_back(i + 1);
      continue;
    case 45:
      if (yr_bool && x[i] <  xa && y[i] > ya) o.push_back(i + 1);
      continue;

    case 50:
      if (yr_bool && x[i] != xa && y[i] < ya) o.push_back(i + 1);
      continue;
    case 51:
      if (yr_bool && x[i] == xa && y[i] < ya) o.push_back(i + 1);
      continue;
    case 52:
      if (yr_bool && x[i] >= xa && y[i] < ya) o.push_back(i + 1);
      continue;
    case 53:
      if (yr_bool && x[i] <= xa && y[i] < ya) o.push_back(i + 1);
      continue;
    case 54:
      if (yr_bool && x[i] >  xa && y[i] < ya) o.push_back(i + 1);
      continue;
    case 55:
      if (yr_bool && x[i] <  xa && y[i] < ya) o.push_back(i + 1);
      continue;
    }
  }
  return wrap(o);
}


// [[Rcpp::export(rng = false)]]
IntegerVector do_which_in(IntegerVector x, IntegerVector tbl) {
  R_xlen_t n = x.size();
  if (n >= INT_MAX) {
    stop("Length exceeds integer max.");  // # nocov
  }

  int m = (int)n;
  int tbln = tbl.size();
  std::vector<int> o;
  o.reserve(m);
  for (int i = 0; i < m; ++i) {
    int xi = x[i];
    if (xi == NA_INTEGER) {
      continue;
    }
    for (int j = 0; j < tbln; ++j) {
      if (x[i] == tbl[j]) {
        o.push_back(i + 1);
        break;
      }
    }
  }
  return wrap(o);
}

// # nocov start
std::vector<int> do_intersect3_stdint(std::vector<int> x,
                                      std::vector<int> y,
                                      std::vector<int> z) {
  int nx = x.size();
  int ny = y.size();
  int nz = z.size();

  std::vector<int> o;
  o.reserve(nx);


  int j = 0;
  int k = 0;
  for (int i = 0; i < nx; ++i) {
    int xi = x[i];
    while (j < ny && y[j] < xi) {
      ++j;
    }
    if (j >= ny) {
      break; // not possible to have intersections
    }
    int yj = y[j];
    if (xi > yj) {
      continue;
    }
    while (k < nz && z[k] < xi) {
      ++k;
    }
    if (k >= nz) {
      break;
    }
    int zk = z[k];
    if (xi == yj && xi == zk) {
      o.push_back(xi);
    }
  }
  return o;
}
// # nocov end

// [[Rcpp::export(rng = false)]]
bool test_intersect3_stdint(int M = 0) {
  std::vector<int> a = {1, 2, 3};
  std::vector<int> b = {1, 2, 3};
  std::vector<int> c = {1, 2, 3};
  if (M == 4) {
    a.push_back(M);
  }
  if (M == 5) {
    a.push_back(5);
    a.push_back(7);
    b.push_back(4);
    c.push_back(6);
  }

  std::vector<int> d = do_intersect3_stdint(a, b, c);
  // # nocov start
  for (int i = 0; i < 3; ++i) {
    if (d[i] != a[i]) {
      return false;
    }
  }
  // # nocov end
  return true;
}


// [[Rcpp::export(rng = false)]]
IntegerVector do_count_logical(LogicalVector x, int nThread = 1) {
  int n = x.length();
  int trues = 0;
  int nas = 0;
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+:trues,nas)
#endif
  for (int i = 0; i < n; ++i) {
    if (x[i] == NA_LOGICAL) {
      nas += 1;
    } else if (x[i]) {
      trues += 1;
    }
  }
  int falses = n - trues - nas;

  IntegerVector o(3);
  o[0] = falses;
  o[1] = trues;
  o[2] = nas;
  return o;
}

// [[Rcpp::export(rng = false)]]
DoubleVector do_count_logical_long(LogicalVector x, int nThread = 1) {
  R_xlen_t n = x.length();
  R_xlen_t trues = 0;
  R_xlen_t nas = 0;
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+:trues,nas)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    if (x[i] == NA_LOGICAL) {
      nas += 1;
    } else if (x[i]) {
      trues += 1;
    }
  }
  R_xlen_t falses = n - trues - nas;

  DoubleVector o(3);
  o[0] = falses;
  o[1] = trues;
  o[2] = nas;
  return o;
}


bool single_ox_x1_x2(int x, int oix, int x1, int x2) {
  switch(oix) {
  case OP_NE:
    return x != x1;

  case OP_EQ:
    return x == x1;

  case OP_GE:
    return x >= x1;

  case OP_LE:
    return x <= x1;

  case OP_GT:
    return x >  x1;

  case OP_LT:
    return x <  x1;

  case OP_BW:
    return x >= x1 && x <= x2;

  case OP_BO:
    return x >  x1 && x <  x2;

  case OP_BC:
    return x <= x1 || x >= x2;
  }
  return false; // # nocov
}

bool single_ox_x1_x2(double x, int oix, double x1, double x2) {
  switch(oix) {
  case OP_NE:
    return x != x1;

  case OP_EQ:
    return x == x1;

  case OP_GE:
    return x >= x1;

  case OP_LE:
    return x <= x1;

  case OP_GT:
    return x >  x1;

  case OP_LT:
    return x <  x1;

  case OP_BW:
    return x >= x1 && x <= x2;

  case OP_BO:
    return x >  x1 && x <  x2;

  case OP_BC:
    return x <= x1 || x >= x2;
  }
  return false; // # nocov
}

// [[Rcpp::export(rng = false)]]
bool do_in_int(int x, IntegerVector table) {
  int tn = table.length();
  for (int i = 0; i < tn; ++i) {
    if (x == table[i]) {
      return true;
    }
  }
  return false;
}

// [[Rcpp::export(rng = false)]]
LogicalVector do_par_in(IntegerVector x, IntegerVector table, int nThread = 1) {
  R_xlen_t n = x.length();
  int tn = table.length();
  std::vector<int> tango;
  tango.reserve(tn);
  for (int ti = 0; ti < tn; ++ti) {
    tango.push_back(table[ti]);
  }
  sort(tango.begin(), tango.end());
  int a = tango[0];
  int b = tango[tn - 1];

  LogicalVector out = no_init(n);

#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) shared(out, tango)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    bool oi = false;
    int xi = x[i];
    if (xi >= a && xi <= b) {
      for (int j = 0; j < tn; ++j) {
        int tj = tango[j];
        if (x[i] == tj) {
          oi = true;
          break;
        }
      }
    }
    out[i] = oi;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
LogicalVector do_par_in_hash_int(IntegerVector x, IntegerVector table, int nThread = 1) {
  std::unordered_set<int> H;
  int tn = table.length();
  for (int t = 0; t < tn; ++t) {
    H.insert(table[t]);
  }

  R_xlen_t N = x.length();
  LogicalVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    int xi = x[i];
    bool oi = H.count(xi);
    out[i] = oi;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
LogicalVector do_par_in_hash_dbl(DoubleVector x, DoubleVector table, int nThread = 1) {
  std::unordered_set<double> H;
  int tn = table.length();
  for (int t = 0; t < tn; ++t) {
    H.insert(table[t]);
  }

  R_xlen_t N = x.length();
  LogicalVector out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (R_xlen_t i = 0; i < N; ++i) {
    double xi = x[i];
    bool oi = H.count(xi);
    out[i] = oi;
  }
  return out;
}








