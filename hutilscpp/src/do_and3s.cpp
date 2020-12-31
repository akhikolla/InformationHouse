#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
LogicalVector do_and3_par(IntegerVector x, int ox, int x1, int x2,
                          IntegerVector y, int oy, int y1, int y2,
                          IntegerVector z, int oz, int z1, int z2,
                          LogicalVector A,
                          LogicalVector B,
                          LogicalVector C,
                          CharacterVector nom,
                          int nThread = 1) {
  int nx = x.length();
  int nA = A.length();
  int n = (nx > nA) ? nx : nA;
  bool useX = x.length() == n;
  bool useY = y.length() == n;
  bool useZ = z.length() == n;

  // Which variables are bare logicals
  bool A_lgl = A.length() == n;
  bool B_lgl = B.length() == n;
  bool C_lgl = C.length() == n;

  // Is the 1st, 2nd, 3rd expression usable or should we just set it to false?
  bool e1 = useX || A_lgl;
  bool e2 = useY || B_lgl;
  bool e3 = useZ || C_lgl;

  // Are the expressions preceded by `!` -- i.e the opposite
  bool A_opposite = A_lgl && ox == 1;
  bool B_opposite = B_lgl && oy == 1;
  bool C_opposite = C_lgl && oz == 1;

  if (useX && A_lgl) {
    stop("Internal error: useX && A_lgl"); // # nocov
  }
  if (useY && B_lgl) {
    stop("Internal error: useY && B_lgl"); // # nocov
  }
  if (useZ && C_lgl) {
    stop("Internal error: useZ && C_lgl"); // # nocov
  }

  if (!e1 && !e2 && !e3) {
    return LogicalVector(n);  // # nocov
  }

  LogicalVector out = no_init(n);

  if (useX && useY && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        single_ox_x1_x2(y[i], oy, y1, y2) &&
        single_ox_x1_x2(z[i], oz, z1, z2);
    }
    return out;
  }

  if (useX && useY && C_opposite) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        single_ox_x1_x2(y[i], oy, y1, y2) &&
        !C[i];
    }
    return out;
  }

  if (useX && useY && C_lgl) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        single_ox_x1_x2(y[i], oy, y1, y2) &&
        C[i];
    }
    return out;
  }

  // e3 = false
  if (useX && useY && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        single_ox_x1_x2(y[i], oy, y1, y2);
    }
    return out;
  }

  // B_lgl (but B_opposite must always precede!)
  if (useX && B_opposite && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        !B[i] &&
        single_ox_x1_x2(z[i], oz, z1, z2);
    }
    return out;
  }

  if (useX && B_lgl && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        B[i] &&
        single_ox_x1_x2(z[i], oz, z1, z2);
    }
    return out;
  }

  if (useX && B_opposite && C_opposite) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        !B[i] &&
        !C[i];
    }
    return out;
  }

  if (useX && B_lgl && C_opposite) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        B[i] &&
        !C[i];
    }
    return out;
  }

  if (useX && B_opposite && C_lgl) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        !B[i] &&
        C[i];
    }
    return out;
  }

  if (useX && B_lgl && C_lgl) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        B[i] &&
        C[i];
    }
    return out;
  }

  if (useX && B_opposite && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        !B[i];
    }
    return out;
  }

  if (useX && B_lgl && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) &&
        B[i];
    }
    return out;
  }

  if (useX && !e2 && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2);
    }
    return out;
  }
  // // A_opposite

  if (A_opposite && useY && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        single_ox_x1_x2(y[i], oy, y1, y2) &&
        single_ox_x1_x2(z[i], oz, z1, z2);
    }
    return out;
  }

  if (A_opposite && useY && C_opposite) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        single_ox_x1_x2(y[i], oy, y1, y2) &&
        !C[i];
    }
    return out;
  }

  if (A_opposite && useY && C_lgl) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        single_ox_x1_x2(y[i], oy, y1, y2) &&
        C[i];
    }
    return out;
  }

  // e3 = false
  if (A_opposite && useY && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        single_ox_x1_x2(y[i], oy, y1, y2);
    }
    return out;
  }

  // B_lgl
  if (A_opposite && B_opposite && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        !B[i] &&
        single_ox_x1_x2(z[i], oz, z1, z2);
    }
    return out;
  }

  if (A_opposite && B_lgl && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        B[i] &&
        single_ox_x1_x2(z[i], oz, z1, z2);
    }
    return out;
  }

  if (A_opposite && B_opposite && C_opposite) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        !B[i] &&
        !C[i];
    }
    return out;
  }

  if (A_opposite && B_lgl && C_opposite) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        B[i] &&
        !C[i];
    }
    return out;
  }

  if (A_opposite && B_opposite && C_lgl) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        !B[i] &&
        C[i];
    }
    return out;
  }

  if (A_opposite && B_lgl && C_lgl) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        B[i] &&
        C[i];
    }
    return out;
  }

  if (A_opposite && B_opposite && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        !B[i];
    }
    return out;
  }

  if (A_opposite && B_lgl && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i] &&
        B[i];
    }
    return out;
  }

  if (A_opposite && !e2 && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        !A[i];
    }
    return out;
  }

  // // A_lgl

  if (A_lgl && useY && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        single_ox_x1_x2(y[i], oy, y1, y2) &&
        single_ox_x1_x2(z[i], oz, z1, z2);
    }
    return out;
  }

  if (A_lgl && useY && C_opposite) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        single_ox_x1_x2(y[i], oy, y1, y2) &&
        !C[i];
    }
    return out;
  }

  if (A_lgl && useY && C_lgl) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        single_ox_x1_x2(y[i], oy, y1, y2) &&
        C[i];
    }
    return out;
  }

  // e3 = false
  if (A_lgl && useY && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        single_ox_x1_x2(y[i], oy, y1, y2);
    }
    return out;
  }

  // B_lgl
  if (A_lgl && B_opposite && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        !B[i] &&
        single_ox_x1_x2(z[i], oz, z1, z2);
    }
    return out;
  }

  if (A_lgl && B_lgl && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        B[i] &&
        single_ox_x1_x2(z[i], oz, z1, z2);
    }
    return out;
  }

  if (A_lgl && B_opposite && C_opposite) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        !B[i] &&
        !C[i];
    }
    return out;
  }

  if (A_lgl && B_lgl && C_opposite) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        B[i] &&
        !C[i];
    }
    return out;
  }

  if (A_lgl && B_opposite && C_lgl) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        !B[i] &&
        C[i];
    }
    return out;
  }

  if (A_lgl && B_lgl && C_lgl) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        B[i] &&
        C[i];
    }
    return out;
  }

  if (A_lgl && B_opposite && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        !B[i];
    }
    return out;
  }

  if (A_lgl && B_lgl && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i] &&
        B[i];
    }
    return out;
  }

  if (A_lgl && !e2 && !e3) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        A[i];
    }
    return out;
  }

#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
  for (int i = 0; i < n; ++i) {
    bool oi = false;

    // 1st expression
    if (e1) {
      oi = useX ? single_ox_x1_x2(x[i], ox, x1, x2) : (A_opposite ? !A[i] : A[i]);
      if (!oi) {
        out[i] = FALSE;
        continue;
      }
    }

    // 2nd expression
    // #nocov start
    if (e2) {
      oi = useY ? single_ox_x1_x2(y[i], oy, y1, y2) : (B_opposite ? !B[i] : B[i]);
      if (!oi) {
        out[i] = FALSE;
        continue;
      }
    }
    // #nocov end
    // 3rd expression
    if (e3) {
      oi = useZ ? single_ox_x1_x2(z[i], oz, z1, z2) : (C_opposite ? !C[i] : C[i]);
      if (!oi) {
        out[i] = FALSE;
        continue;
      }
    }

    out[i] = TRUE;

  }
  return out;
}

// [[Rcpp::export]]
R_xlen_t do_sum3s_par(IntegerVector x, int ox, int x1, int x2,
                      DoubleVector xd, double xd1, double xd2,
                      IntegerVector y, int oy, int y1, int y2,
                      DoubleVector yd, double yd1, double yd2,
                      IntegerVector z, int oz, int z1, int z2,
                      DoubleVector zd, double zd1, double zd2,
                      LogicalVector A,
                      LogicalVector B,
                      LogicalVector C,
                      bool ampersand = true,
                      int nThread = 1) {
  // ampersand TRUE => sum_and3,  FALSE  => sum_or3

  R_xlen_t lengths[9] = {x.length(), xd.length(),
                         y.length(), yd.length(),
                         z.length(), zd.length(),
                         A.length(), B.length(), C.length()};

  R_xlen_t n = *std::max_element(lengths, lengths + 9);

  // 0 use none, 1 use int, 2 use double, 3 use logical, 4 use opposite
  const int Xcase = (x.length() == n) + 2 * (xd.length() == n) + 3 * (A.length() == n && ox != 1) + 4 * (A.length() == n && ox == 1);
  const int Ycase = (y.length() == n) + 2 * (yd.length() == n) + 3 * (B.length() == n && oy != 1) + 4 * (B.length() == n && oy == 1);
  const int Zcase = (z.length() == n) + 2 * (zd.length() == n) + 3 * (C.length() == n && oz != 1) + 4 * (C.length() == n && oz == 1);


  R_xlen_t out = 0;

#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread) reduction(+ : out)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    bool oi = Xcase == 0;

    // 1st expression
    switch(Xcase) {
    case 1:
      oi = single_ox_x1_x2(x[i], ox, x1, x2);
      break;
    case 2:
      oi = single_ox_x1_x2(xd[i], ox, xd1, xd2);
      break;
    case 3:
      oi = A[i];
      break;
    case 4:
      oi = !A[i];
      break;
    }
    if (ampersand) {
      if (!oi) {
        continue;
      }
    } else {
      if (oi) {
        out += 1;
        continue;
      }
    }


    // 2nd expression

    switch(Ycase) {
    case 1:
      oi = single_ox_x1_x2(y[i], oy, y1, y2);
      break;
    case 2:
      oi = single_ox_x1_x2(yd[i], oy, yd1, yd2);
      break;
    case 3:
      oi = B[i];
      break;
    case 4:
      oi = !B[i];
      break;
    }
    if (ampersand) {
      if (!oi) {
        continue;
      }
    } else {
      if (oi) {
        out += 1;
        continue;
      }
    }

    // 3rd expression
    switch(Zcase) {
    case 1:
      oi = single_ox_x1_x2(z[i], oz, z1, z2);
      break;
    case 2:
      oi = single_ox_x1_x2(zd[i], oz, zd1, zd2);
      break;
    case 3:
      oi = C[i];
      break;
    case 4:
      oi = !C[i];
      break;
    }
    if (ampersand) {
      if (!oi) {
        continue;
      }
      out += 1;
    } else {
      if (oi) {
        out += 1;
        continue;
      }
    }

  }
  return out;
}


