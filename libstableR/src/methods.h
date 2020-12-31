/* stable/methods.h
 * 
 * Numerical methods employed by Libstable. Some methods extracted from:
 *  Press, W.H. et al. Numerical Recipes in C: the art of scientific
 *    computing. Cambridge University Press, 1994. Cambridge, UK.
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
#ifndef _METHODS_H_
#define _METHODS_H_

#define MAX(a,b,c,d) ( a > b ? c : d)

#include "mcculloch.h"
#include <stddef.h>

double gammaln(const double x);

double zbrent(double (*func)(double x, void *args), void * args, double x1, double x2, double value, const double tol, int *warn);

double dfridr(double (*func)(double, void * args),void * args, double x,  double h, double *err);

double quadstep(double (*func)(double x, void *args),void * args, double a, double b,double fa, double fc,double fb, const double epsabs, const double epsrel, double *abserr, int *warn, size_t *fcnt);

double qromb(double (*func)(double, void *),void *args, double a, double b, double epsabs, double epsrel, int K, int JMAX, int method, int *warn, size_t *fcnt, double *err);

double trapzd(double (*func)(double, void *), void * args, double a, double b, int n, double s);

void polint(const double xa[], const double ya[], const int n, double x, double *y, double *dy);

void ratint(const double xa[], const double ya[], const int n, double x, double *y, double *dy);

void medfilt(const double xa[], const double ya[], int n, double N);

void vector_step(double **x, double min, double max, double step, int *n);

void vector_npoints(double **x, double min, double max, int n, double * step);
#endif
