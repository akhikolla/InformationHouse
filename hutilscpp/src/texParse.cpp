#include "cpphutils.h"

// #nocov start
int showValuea(const char* what, double x) {
  // Rcout << "The value " << what << " is " << x << std::endl;
  return 0;
}
// #nocov end

// [[Rcpp::export]]
int validate_nchar1(CharacterVector x, bool return_size = false) {
  R_xlen_t N = x.length();

  for (R_xlen_t i = 0; i < N; ++i) {
    if (x[i].size() > 1) {
      return return_size ? x[i].size() : i + 1;
    }
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
int max_charsize(CharacterVector x) {
  R_xlen_t N = x.length();
  int o = 0;
  for (R_xlen_t i = 0; i < N; ++i) {
    int xis = x[i].size();
    o = xis > o ? xis : o;
  }
  return o;
}

// [[Rcpp::export]]
bool is_space(CharacterVector x) {
  char x0 = x[0][0];
  return x0 == SPACE;
}

//' @name where_square_bracket_opens
//' @param x Character vector of characters.
//' @param i position of closing bracket.
//'
//' @return
//' -1 if x[i] does not closing bracket
//'  0 if bracket never closes
//'  j the location of the closing brace
//'
//' @noRd

// [[Rcpp::export(rng = false)]]
R_xlen_t where_square_bracket_opens(CharacterVector x, R_xlen_t i = 0) {
  R_xlen_t N = x.length();

  if (i < 0 || i >= N) {
    return -1;
  }
  char xi = x[i][0];
  if (xi != STOP_SQBRK) {
    return -1;
  }
  int depth = 0;
  for (R_xlen_t k = i; k >= 0; --k) {
    char xk = x[k][0];
    depth += (xk == STOP_SQBRK) - (xk == OPEN_SQBRK);
    if (depth == 0) {
      return k;
    }
  }

  return 0;
}
#if false
/// [[Rcpp::export(rng = false)]]
R_xlen_t where_brace_opens(CharacterVector x, R_xlen_t i = 0) {
  R_xlen_t N = x.length();

  if (i < 0 || i >= N) {
    return -1;
  }
  char xi = x[i][0];
  if (xi != STOP_BRACE) {
    return -1;
  }
  int depth = 0;
  for (R_xlen_t k = i; k >= 0; --k) {
    char xk = x[k][0];
    depth += (xk == STOP_BRACE) - (xk == OPEN_BRACE);
    if (depth == 0) {
      return k;
    }
  }

  return 0;
}

/// [[Rcpp::export(rng = false)]]
R_xlen_t locate_prev_brace(CharacterVector x, R_xlen_t i = 0) {
  R_xlen_t N = x.length();
  if (i >= N) {
    return -1;
  }
  while (--i >= 0) {
    char xi = x[i][0];
    if (xi == OPEN_BRACE) {
      return i;
    }
  }
  return -1;
}

/// [[Rcpp::export(rng = false)]]
R_xlen_t locate_next_brace(CharacterVector x, R_xlen_t i = 0) {
  R_xlen_t N = x.length();
  if (i < 0) {
    return -1;
  }
  while (++i <= N) {
    char xi = x[i][0];
    if (xi == STOP_BRACE) {
      return i;
    }
  }
  return -1;
}


/// [[Rcpp::export(rng = false)]]
IntegerVector tex_group(CharacterVector x) {
  int n = x.length();
  if (n > 1e7) {
    stop("Too large");
  }
  int N = 100 * n; // proxy for number of chars
  std::vector<char> o;
  o.reserve(N);

  std::vector<int> tg;
  tg.reserve(N);

  int current_tex_group = 0;
  for (int li = 0; li < n; ++li) {
    int ncharli = x[li].size();
    for (int ci = 0; ci < ncharli; ++ci) {
      char xci = x[li][ci];
      o.push_back(xci);
      current_tex_group += (xci == OPEN_BRACE) - (xci == STOP_BRACE);
      tg.push_back(current_tex_group);
    }
  }

  return wrap(tg);
}





/// [[Rcpp::export(rng = false)]]
IntegerVector tmp_mark_work(CharacterVector x, CharacterVector Command) {
  R_xlen_t N = x.length();
  R_xlen_t CN = Command.length();

  IntegerVector out(N);
  char c0 = Command[0][0];
  R_xlen_t i = 0;
  while (i < (N - CN)) {
    char xi = x[i][0];
    if (xi != c0) {
      ++i;
      continue;
    }
    // cj is the index throughout Command
    R_xlen_t c_j = 1;
    bool i_starts_command = true;
    while (c_j < CN) {
      char Cj = Command[c_j][0];
      char xj = x[i + c_j][0];
      if (Cj != xj) {
        i_starts_command = false;
        break;
      }
      ++c_j;
    }
    if (!i_starts_command) {
      ++i; // can't be tempted to i += c_j because \text{\textbf}
      continue;
    }
    out[i] += 1;
    for (R_xlen_t c_k = 1; c_k < CN; ++c_k) {
      ++i;
      out[i] += 1;
    }
  }
  return out;
}

/// [[Rcpp::export(rng = false)]]
IntegerVector MatchArgCont(CharacterVector x, CharacterVector Command) {
  R_xlen_t N = x.length();
  R_xlen_t CN = Command.length();

  // Different to max(nchar()); this fn is for purpose of array size
  int max_nchar_x = max_charsize(x);
  int max_nchar_c = max_charsize(Command);
  if (max_nchar_c > max_nchar_x) {
    // Implies Command can never occur in x
    return IntegerVector(N);
  }
  if (max_nchar_c != 1) {
    stop("Unexpected strings in Command");
  }
  IntegerVector out = no_init(N);
  char c0 = Command[0][0];

  constexpr char SPACE = 32;
  constexpr char OPEN_BRACE = 123;
  constexpr char STOP_BRACE = 125;
  constexpr char OPEN_SQBRK = 91;
  constexpr char STOP_SQBRK = 93;


  for (R_xlen_t i = CN; i < N; ++i) {
    char xi = x[i][0];
    if (xi == 123) {
      // i = "{"
      R_xlen_t j = i;
      R_xlen_t k = CN;
      bool maybe_command = true;
      while (--j >= 0 && maybe_command) {
        char xj = x[j][0];
        if (xj == SPACE) {
          continue;
        }
        if (xj == STOP_SQBRK) {

        }
      }
    }
  }


  return out;


}
#endif

// [[Rcpp::export]]
List extractMandatory (CharacterVector x, CharacterVector command, int nCommands) {
  R_xlen_t N = x.length();
  int command_len = command.length();
  // const char *xp = x[0];
  // const char *cp = command[0];

  bool finish_extract = false;
  bool within_brace = false;
  // bool within_optional = false;
  int opt_group = 0;
  int cj = 0;


  int command_no = 0;

  CharacterVector support(N);
  IntegerVector commandNo(nCommands);
  IntegerVector commandOpeners(nCommands);
  IntegerVector commandClosers(nCommands);
  for (R_xlen_t i = 0; i < N; ++i) {
    int k = 0;
    if (x[i] == command[0]) {
      showValuea("k", k);
      for (int ci = 1; ci < command_len; ++ci) {
        // Consider command = \abc and document ends with ab
        if (i + ci >= N) {
          break;
        }
        const char *cii = command[ci];
        const char *xic = x[i + ci];
        // if (x[i + ci] == command[ci]) {
        showValuea("i + ci", i + ci);
        showValuea("cj", cj);
        if (*xic == *cii) {
          ++cj;
        } else {
          cj = 0;
          break;
        }
      }
      // Zero if at the same group level as the command name
      int rel_group = 0;
      if (cj == command_len - 1) {
        // Now wait until we see an opening brace
        within_brace = false;
        finish_extract = false;
        k = i + cj;
        while (!within_brace && k < N - 1) {
          showValuea("++k1", k);
          ++k;
          if (x[k] == " " || x[k] == "") {
            // \abc {def}
            continue;
          }
          if (x[k] == "[") {
            ++opt_group;
            int rel_opt_group = 1;
            while (rel_opt_group && k < N - 1) {
              // just keep moving forward until we get out of the current
              // optional group.
              ++k;
              showValuea("++k1_rel_group", k);

              if (x[k] == "]") {
                --rel_opt_group;
                --opt_group;
              } else {
                if (x[k] == "{") {
                  ++rel_group;
                  while (rel_group && k < N - 1) {
                    ++k;
                    if (x[k] == "}") {
                      --rel_group;
                    } else {
                      if (x[k] == "{") {
                        ++rel_group;
                      }
                    }
                  }
                }
                if (x[k] == "[") {
                  ++rel_opt_group;
                  ++opt_group;
                }
              }
            }
            showValuea("++k2", k);
            // Current position is on the closing ']'
            ++k; // move to next character once we're done with optional group
          }
          if (k >= N) { // in case the document is not well-formed.
            break;
          }
          within_brace = (rel_group == 0) && x[k] == "{";

          // abc{xyz} but not abcd{xyz}
          if (x[k] != "" && x[k] != " " && x[k] != "{") {
            showValuea("x[k]", k);
            k = N + 1;  // effectively break both while loops
          } else {
            commandOpeners[command_no] = k + 1;
          }
        }

        showValuea("at while(k < N), k", k);
        showValuea("at while(k < N), N", N);
        while (k < N) {
          support[k] = x[k];
          // #nocov start
          if (command_no >= nCommands) {
            showValuea("command_no = ", command_no);
            showValuea("k = ", k);
            stop("command_no overflow");
          }
          // #nocov end
          // R indexing
          commandNo[command_no] = command_no + 1;
          if (x[k] == "{") {
            ++rel_group;
          } else {
            if (x[k] == "}") {
              --rel_group;
            }
          }


          finish_extract = (rel_group == 0) && x[k] == "}";
          if (finish_extract) {
            commandClosers[command_no] = k + 1;
            ++command_no;
            break;
          }
          ++k;
        }
      }
    }
  }

  List out = List::create(Named("support") = support,
                          Named("Openers") = commandOpeners,
                          Named("Closers") = commandClosers);


  return out;
}



