#include <Rcpp.h>
#include "CLK.h"

using namespace Rcpp;

#ifndef WOLFRAMRULE30_H
#define WOLFRAMRULE30_H

CharacterVector WolframRule30(CharacterVector CLK, int lenBloom, int t);
void WolframRule30c(CLK* clkin, CLK* clkout, int t);

#endif
