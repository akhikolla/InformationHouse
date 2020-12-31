#include "cpphutils.h"

// [[Rcpp::export(rng = false)]]
int do_op2M(CharacterVector x) {
  int nchar = x[0].size();
  char x00 = x[0][0];
  char x01 = (nchar >= 2) ? x[0][1] : x00;
  switch(x00) {
  case '!':
    return OP_NE;
    break;
  case '=':
    return OP_EQ;
    break;
  case '>':
    return (nchar == 1) ? OP_GT : OP_GE;
    break;
  case '<':
    return (nchar == 1) ? OP_LT : OP_LE;
    break;
  case '%':
    switch (x01) {
    case 'i':
      return OP_IN;
      break;
    case 'b':
      return OP_BW;
      break;
    case '(':
      return OP_BO;
      break;
    case ']':
      return OP_BC;
      break;
    }
  }
  return 0;
}
