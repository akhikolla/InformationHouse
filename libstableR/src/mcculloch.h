/* stable/mcculloch.h
 * 
 * Stable Distribution Parameter Estimation by Quantile Method
 * (McCulloch 96)
 * Based on original GAUSS implementation available in
 * <http://www.econ.ohio-state.edu/jhm/jhm.html>.
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
#ifndef _MCCULLOCH_H_
#define _MCCULLOCH_H_

int stab(const double *x, const unsigned int n, unsigned int symm, double *alpha,
         double *beta, double *c, double *zeta);

void cztab(double *x, unsigned int n, double *cn, double *zn);

void czab(double alpha, double beta, double cn, double q50, double *c, double *zeta);

#endif
