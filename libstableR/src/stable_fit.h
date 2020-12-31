/* stable/stable_fit.h
 * 
 * Functions employed by different methods of estimation implemented
 * in Libstable.
 *
 * Copyright (C) 2013. Javier Royuela del Val
 *                     Federico Simmross Wattenberg
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  Javier Royuela del Val.
 *  E.T.S.I. Telecomunicación
 *  Universidad de Valladolid
 *  Paseo de Belén 15, 47002 Valladolid, Spain.
 *  jroyval@lpi.tel.uva.es    
 */

#ifndef _STABLE_FIT_H_
#define _STABLE_FIT_H_

#include "stable.h"


gsl_complex stable_samplecharfunc_point(const double* x,
            const unsigned int N, double t);

void stable_samplecharfunc(const double* x, const unsigned int Nx,
             const double* t, const unsigned int Nt, gsl_complex *z);

#endif
