/***************************************************************************
                              selectRegGen.hpp  
                             ----------------------------
    copyright         : INRIA, INRA, Université Paris-Sud 11
    years             : 2008-2009 all rights reserved
    authors           : Fatou DIA, Marie-Laure MARTIN-MAGNIETTE, Cathy MAUGIS
    email             : cathy.maugis@math.u-psud.fr 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *  
     Permission to use, copy, modify, and distribute this software and
     its documentation for any purpose and without fee is hereby
     granted, provided that the above copyright notice appear in all
     copies and that both that the copyright notice and this
     permission notice and warranty disclaimer appear in supporting
     documentation, and that the name of Lucent or any of its entities
     not be used in advertising or publicity pertaining to
     distribution of the software without specific, written prior
     permission.

 ***************************************************************************/

#include "vect.hpp"

#ifndef SELECTREGGEN_H 
#define SELECTREGGEN_H 


class SelectRegGen{
  
  Vect v;    

public:
  /**** Constructors ****/
  SelectRegGen();  
  SelectRegGen(Vect v); 

  /**** Methods ****/
  
  //Exclusion step
  void exclusion_reggen(vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& stopreg, int& nummodel, int& InitialProjectsNb);

  //Inclusion step
  void inclusion_reggen(vector<int> varSelect, vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& stopreg, int& nummodel, int& InitialProjectsNb);

  //backward stepwise algorithm for the multidimensional multivariate regression
  vector<int> selectReggen(vector<int> varSelect, vector<int>& varNonSig, int nummodel, int& InitialProjectsNb);

};

#endif
