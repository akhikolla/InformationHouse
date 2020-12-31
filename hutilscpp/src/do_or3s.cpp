#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
LogicalVector do_or3_par(IntegerVector x, int ox, int x1, int x2,
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
    return LogicalVector(n); // # nocov
  }

  LogicalVector out = no_init(n);

  if (useX && useY && useZ) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
    for (int i = 0; i < n; ++i) {
      out[i] =
        single_ox_x1_x2(x[i], ox, x1, x2) ||
        single_ox_x1_x2(y[i], oy, y1, y2) ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
        single_ox_x1_x2(y[i], oy, y1, y2) ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
        single_ox_x1_x2(y[i], oy, y1, y2) ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
        !B[i] ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
        B[i] ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
        !B[i] ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
        B[i] ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
        !B[i] ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
        B[i] ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
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
        single_ox_x1_x2(x[i], ox, x1, x2) ||
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
        !A[i] ||
        single_ox_x1_x2(y[i], oy, y1, y2) ||
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
        !A[i] ||
        single_ox_x1_x2(y[i], oy, y1, y2) ||
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
        !A[i] ||
        single_ox_x1_x2(y[i], oy, y1, y2) ||
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
        !A[i] ||
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
        !A[i] ||
        !B[i] ||
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
        !A[i] ||
        B[i] ||
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
        !A[i] ||
        !B[i] ||
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
        !A[i] ||
        B[i] ||
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
        !A[i] ||
        !B[i] ||
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
        !A[i] ||
        B[i] ||
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
        !A[i] ||
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
        !A[i] ||
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
        A[i] ||
        single_ox_x1_x2(y[i], oy, y1, y2) ||
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
        A[i] ||
        single_ox_x1_x2(y[i], oy, y1, y2) ||
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
        A[i] ||
        single_ox_x1_x2(y[i], oy, y1, y2) ||
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
        A[i] ||
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
        A[i] ||
        !B[i] ||
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
        A[i] ||
        B[i] ||
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
        A[i] ||
        !B[i] ||
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
        A[i] ||
        B[i] ||
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
        A[i] ||
        !B[i] ||
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
        A[i] ||
        B[i] ||
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
        A[i] ||
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
        A[i] ||
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
      if (oi) {
        out[i] = oi;
        continue;
      }
    }

    // 2nd expression
    // #nocov start
    if (e2) {
      oi = useY ? single_ox_x1_x2(y[i], oy, y1, y2) : (B_opposite ? !B[i] : B[i]);
      if (oi) {
        out[i] = TRUE;
        continue;
      }
    }
    // #nocov end
    // 3rd expression
    if (e3) {
      oi = useZ ? single_ox_x1_x2(z[i], oz, z1, z2) : (C_opposite ? !C[i] : C[i]);
      if (oi) {
        out[i] = TRUE;
        continue;
      }
    }

    out[i] = FALSE;

  }
  return out;
}

/*
namespace h3a {

template <int RTYPE>
LogicalVector and1s(const Vector<RTYPE>& x, int opx,
                    const Vector<RTYPE>& x1,
                    const Vector<RTYPE>& x2,
                    int nThread) {
  R_xlen_t N = x.length();
  R_xlen_t nx1 = x1.length();

  const bool r1 = nx1 == 1;
  const bool rn = nx1 == N;


  LogicalVector out = no_init(N);

  switch(opx) {
  case 1:
    switch(TYPEOF(x)) {
    case INTSXP:
      switch(TYPEOF(x1)) {
      case INTSXP:
        if (r1) {
          const int rhs1 = x1[0];
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
          for (R_xlen_t i = 0; i < N; ++i) {
            out[i] = x[i] != rhs1;
          }
          return out;
        }
        if (rn) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
          for (R_xlen_t i = 0; i < N; ++i) {
            out[i] = x[i] != x1[i];
          }
          return out;
        }
      case REALSXP:
        if (r1) {
          const double rhs1 = x1[0];
          const bool isnt_integerish =
            !R_finite(rhs1) || rhs1 > 2147483647 || rhs < -2147483647 ||
            (rhs1 != ((int)rhs1));

#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
          for (R_xlen_t i = 0; i < N; ++i) {
            out[i] = isnt_integerish ? TRUE : (x[i] != ((int)(rhs1)));
          }
        }
        if (rn) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
          for (R_xlen_t i = 0; i < N; ++i) {
            const double rhs1 = x1[i];
            const bool isnt_integerish =
              !R_finite(rhs1) || rhs1 > 2147483647 || rhs < -2147483647 ||
              (rhs1 != ((int)rhs1));
            out[i] = isnt_integerish ? TRUE : (x[i] != ((int)(rhs1)));
          }
        }
      }



    case REALSXP:
      switch(TYPEOF(x1)) {
      case INTSXP:
        if (r1) {
          const int rhs1 = x1[0];
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
          for (R_xlen_t i = 0; i < N; ++i) {
            out[i] = x[i] != rhs1;
          }
          return out;
        }
        if (rn) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
          for (R_xlen_t i = 0; i < N; ++i) {
            out[i] = x[i] != x1[i];
          }
          return out;
        }
      case REALSXP:
        // int != dbl
        // if it's not integerish then all true
        if (r1) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
          for (R_xlen_t i = 0; i < N; ++i) {
            out[i] = x[i] != x1[0];
          }
        }
        if (rn) {
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(nThread)
#endif
          for (R_xlen_t i = 0; i < N; ++i) {
            const double rhs1 = x1[i];
            out[i] = x[i] != rhs1
          }
        }
      }
    }
  }
  return Rcpp::LogicalVector();
}
// [[Rcpp::export(rng = false)]]
LogicalVector do_and1(SEXP x, int opx,
                      SEXP x1,
                      SEXP x2,
                      int nThread = 1) {
  R_xlen_t N = x.length();
  Vector<RTYPE> out = no_init(N);
#if defined _OPENMP && _OPENMP >= 201511
#pragma omp parallel for num_threads(n)
#endif
  for (int i = 0; i < N; ++i) {
    out[i] = (x[i] <= 0) ? 0 : x[i];
  }
  return out;
}
}
*/








