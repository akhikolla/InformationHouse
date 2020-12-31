#include <Rcpp.h>
#include "CLK.h"

using namespace Rcpp;

#ifndef WOLFRAMRULE90_H
#define WOLFRAMRULE90_H

CharacterVector WolframRule90(CharacterVector CLK, int lenBloom, int t);
void WolframRule90c(CLK* clkin, CLK* clkout, int t);

#endif
