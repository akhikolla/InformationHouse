//
//  Function.h
//  
//
//  Created by sedki on 09/04/2014.
//
//

#ifndef ____Function__
#define ____Function__

#include <iostream>

long double Quad_Form(colvec x, colvec mu, mat S_Inv);

long double ldcppmvt(colvec x, colvec mu, mat SInv, double SLogDet);
#endif /* defined(____Function__) */
