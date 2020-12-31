/* stable/stable_integration.h
 * 
 * Functions to perform numerical integration used when calculating
 * the PDF and CDF of alpha-stable distributions. Based on GSL
 * numerical quadrature methods.
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
 * 
 */
#ifndef _STABLE_INTEGRATION_H_
#define _STABLE_INTEGRATION_H_

#include "stable.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_math.h>


int stable_integration_METHODNAME(unsigned short method, char *s);

void
stable_integration(StableDist *dist,double(function)(double, void*),
                   double a, double b,
                   double epsabs, double epsrel, unsigned short limit,
                   double *result, double *abserr, unsigned short method);

#endif
