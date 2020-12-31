/* stable/mcculloch.c
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

#include "mcculloch.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

//ENA IS THE INDEX NU SUB ALPHA DESCRIBED IN MCCULLOCH (1986).
double ena[16][5] = {
        {2.4388, 2.4388, 2.4388, 2.4388, 2.4388},
        {2.5120, 2.5117, 2.5125, 2.5129, 2.5148},
        {2.6080, 2.6093, 2.6101, 2.6131, 2.6174},
        {2.7369, 2.7376, 2.7387, 2.7420, 2.7464},
        {2.9115, 2.9090, 2.9037, 2.8998, 2.9016},
        {3.1480, 3.1363, 3.1119, 3.0919, 3.0888},
        {3.4635, 3.4361, 3.3778, 3.3306, 3.3161},
        {3.8824, 3.8337, 3.7199, 3.6257, 3.5997},
        {4.4468, 4.3651, 4.1713, 4.0052, 3.9635},
        {5.2172, 5.0840, 4.7778, 4.5122, 4.4506},
        {6.3140, 6.0978, 5.6241, 5.2195, 5.1256},
        {7.9098, 7.5900, 6.8606, 6.2598, 6.1239},
        {10.4480, 9.9336, 8.7790, 7.9005, 7.6874},
        {14.8378, 13.9540, 12.0419, 10.7219, 10.3704},
        {23.4831, 21.7682, 18.3320, 16.2163, 15.5841},
        {44.2813, 40.1367, 33.0018, 29.1399, 27.7822}
};

// ENB IS THE INDEX NU SUB BETA.
double enb[16][5] = {
        {0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {0.0000, 0.0179, 0.0357, 0.0533, 0.0710},
        {0.0000, 0.0389, 0.0765, 0.1133, 0.1480},
        {0.0000, 0.0626, 0.1226, 0.1784, 0.2281},
        {0.0000, 0.0895, 0.1736, 0.2478, 0.3090},
        {0.0000, 0.1183, 0.2282, 0.3199, 0.3895},
        {0.0000, 0.1478, 0.2849, 0.3942, 0.4686},
        {0.0000, 0.1769, 0.3422, 0.4703, 0.5458},
        {0.0000, 0.2062, 0.3993, 0.5473, 0.6210},
        {0.0000, 0.2362, 0.4561, 0.6240, 0.6934},
        {0.0000, 0.2681, 0.5134, 0.6993, 0.7616},
        {0.0000, 0.3026, 0.5726, 0.7700, 0.8248},
        {0.0000, 0.3415, 0.6343, 0.8339, 0.8805},
        {0.0000, 0.3865, 0.6994, 0.8900, 0.9269},
        {0.0000, 0.4408, 0.7678, 0.9362, 0.9620},
        {0.0000, 0.5095, 0.8381, 0.9700, 0.9847}
};

// ENC IS THE INDEX NU SUB C
double enc[16][5] = {
        {1.9078, 1.9078, 1.9078, 1.9078, 1.9078},
        {1.9140, 1.9150, 1.9160, 1.9185, 1.9210},
        {1.9210, 1.9220, 1.9275, 1.9360, 1.9470},
        {1.9270, 1.9305, 1.9425, 1.9610, 1.9870},
        {1.9330, 1.9405, 1.9620, 1.9970, 2.0430},
        {1.9390, 1.9520, 1.9885, 2.0450, 2.1160},
        {1.9460, 1.9665, 2.0220, 2.1065, 2.2110},
        {1.9550, 1.9845, 2.0670, 2.1880, 2.3330},
        {1.9650, 2.0075, 2.1255, 2.2945, 2.4910},
        {1.9800, 2.0405, 2.2050, 2.4345, 2.6965},
        {2.0000, 2.0850, 2.3115, 2.6240, 2.9735},
        {2.0400, 2.1490, 2.4610, 2.8865, 3.3565},
        {2.0980, 2.2445, 2.6765, 3.2650, 3.9125},
        {2.1890, 2.3920, 3.0040, 3.8440, 4.7755},
        {2.3370, 2.6355, 3.5425, 4.8085, 6.2465},
        {2.5880, 3.0735, 4.5340, 6.6365, 9.1440}
};

// ZA IS THE INDEX NU SUB ZETA.
// DELTA IS ESTIMATED FROM ZETA, SO NU SUB DELTA IS NOT USED.
double za[16][5] = {
        {0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {0.0000,-0.0166,-0.0322,-0.0488,-0.0644},
        {0.0000,-0.0302,-0.0615,-0.0917,-0.1229},
        {0.0000,-0.0434,-0.0878,-0.1321,-0.1785},
        {0.0000,-0.0556,-0.1113,-0.1699,-0.2315},
        {0.0000,-0.0660,-0.1340,-0.2060,-0.2830},
        {0.0000,-0.0751,-0.1542,-0.2413,-0.3354},
        {0.0000,-0.0837,-0.1733,-0.2760,-0.3896},
        {0.0000,-0.0904,-0.1919,-0.3103,-0.4467},
        {0.0000,-0.0955,-0.2080,-0.3465,-0.5080},
        {0.0000,-0.0980,-0.2230,-0.3830,-0.5760},
        {0.0000,-0.0986,-0.2372,-0.4239,-0.6525},
        {0.0000,-0.0956,-0.2502,-0.4688,-0.7424},
        {0.0000,-0.0894,-0.2617,-0.5201,-0.8534},
        {0.0000,-0.0779,-0.2718,-0.5807,-0.9966},
        {0.0000,-0.0610,-0.2790,-0.6590,-1.1980}
};
// The above are used in transposed form ([16,5])



double frctl (const double *xx, double p, unsigned int n)
{
// This proc obtains selected quantiles off the ordered vector x
// La he modificado ligeramente respecto a la original para contemplar
// los casos 0<p<0.5/n y 1>p>(n-0.5)/n
double zi, th;
int i;
zi = p*n - .5; //modo R-5: h=Np+1/2
i = floor(zi);
if (zi<0) return xx[0];
else if (zi>n-1) return xx[n-1];
th = zi - i;
return ((1-th)*xx[i] + th * xx[i+1]);
}

int stab(const double *x, const unsigned int n, unsigned int symm,
         double *alpha, double *beta, double *c, double *zeta)
{
/* Main function. x must be ordered to get the percentiles */

  double ah[5], bh[16], q05, q25, q50, q75, q95;
  double d, cn, cnu, an, bn, sign;//, delta=0.0;
  double aa, b, bb, t1, t2, s1, s2, dt, ds, s, t;
  int i, j, ii, jj;

  q05 = frctl(x,.05,n);
  q25 = frctl(x,.25,n);
  q50 = frctl(x,.50,n);
  q75 = frctl(x,.75,n);
  q95 = frctl(x,.95,n);
  
  d = q95 - q05;
  cn = q75 - q25;
  if (cn == 0)
    {
      perror("Warning -- interquartile range is 0\n");
      return(-1);
    }

  an = d/cn;
  if (an < 2.4388)
    {
      *alpha=1.95;
      *beta=0;
      *c=cn/1.9078;
      //delta=q50;
      *zeta=q50;

 /*     printf("alpha: %lf beta: %lf c: %lf d: %lf z: %lf",
             *alpha,*beta,*c,delta,*zeta);*/
      return 2;
    }


  switch (symm)
    {
      case 0:
        bn = (q05 + q95 - 2 * q50)/d;
        //The sign of beta is saved and restored at the end.
        sign = 1.0;
        if (bn < 0) sign = -1;
        else if (bn == 0) sign = 0;
        bn = fabs(bn);

        //printf("an= %lf bn= %lf\n",an, bn);

        for (j=0;j<5;j++)
          {
            for (i=1;i<15;i++)
              {
                if (an <= ena[i][j]) break;
              }
            t = (an - ena[i-1][j])/(ena[i][j]-ena[i-1][j]);
            if (t>1) ah[j] = 0.5;
            else ah[j] = 2 - (i-1+t)*.1;
            //printf("ah[%d]=%lf, i:%d %.4f\n",j,ah[j],i,t);
          }
        for (i=1;i<16;i++)
          {
            for (j=1;j<4;j++)
              {
                if (bn <= enb[i][j]) break;
              }
            t = (bn - enb[i][j-1])/(enb[i][j]-enb[i][j-1]);
            if (t>1) bh[i] = 1.0;
            else bh[i] = (j-1+t)*.25;
            //printf("bh[%d]=%lf, j:%d %.4f\n",i,bh[i],j,t);
          }
        bh[0] = 2 * bh[1] - bh[2]; //Linear extrapolation

        for (j=1;j<5;j++)
          {
            jj=j;
            i=floor((2-ah[j])*10+1);
            if (i<1) i=1;
            else if (i>15) i=15;
            aa = 2 - .1*(i-1);
            s2 = -(ah[j] - aa) * 10;
            b = (1 - s2) * bh[i-1] + s2 * bh[i];
            //printf("%d  %d,  ",i,j);
            if (b < j/4.0) break;
          }
        j=jj;

        bb = .25 * (j-1);
        t1 = (bh[i-1] - bb) / .25;
        t2 = (bh[i] - bb) / .25;
        s1 = - (ah[j-1] - aa) * 10;
        dt = t2 - t1;
        ds = s2 - s1;
        s = (s1 + t1 * ds) / (1 - ds * dt);
        t = t1 + s * dt;
        *alpha = aa - s * .1;
        if (*alpha < .5) *alpha = .5;
        *beta = bb + t * .25;
        //printf("beta: %lf %lf %lf %lf",*beta,bb,t,s);
        if (*beta > 1) *beta = 1.0;
        *beta = *beta * sign;

        *c = enc[i-1][j-1] * (1-s) * (1-t) + enc[i][j-1] * s * (1-t) +
             enc[i-1][j] * t * (1-s) + enc[i][j] * t * s;
        *c = cn/(*c);
        *zeta = za[i-1][j-1] * (1-s) * (1-t) + za[i][j-1] * s * (1-t) +
                za[i-1][j] * t * (1-s) + za[i][j] * t * s;
        *zeta = q50 + *c * sign * (*zeta);

        break;
      case 1:
        for(i = 1;i<15;i++)
          {
            ii = i;
            if (an <= ena[i][0]) break;
          }
        i = ii;

        t = (an - ena[ i-1][0]) / (ena[i][0] - ena[i-1][0]);
        *alpha = (21 - i - t) / 10;
        if (*alpha < .5) *alpha = .5;
        cnu = enc[ i-1][0] * (1-t) + enc[i][0] * t;
        *c = cn/cnu;
        *beta=0;
        //delta=q50;
        *zeta=q50;
        break;
    }

  //printf("alpha: %lf beta: %lf c: %lf d: %lf z: %lf\n",
    //       *alpha,*beta,*c,delta,*zeta);

  return 0;

}

void cztab(double *x, unsigned int n, double *cn, double *zn)
{
/* Return McCulloch's interquartile range and median */
  if (cn!=NULL) *cn = frctl(x,.75,n) - frctl(x,.25,n);
  if (zn!=NULL) *zn = frctl(x,.5 ,n);
  return;
}

void czab(double alpha, double beta, double cn, double q50,
          double *c, double *zeta)
{
/* Estimate c and zeta from current alpha and beta estimations, the interquartile range
 * and the median */

  double s,t;
  int i,j,sign;

  sign=1;
  if (beta<0) sign=-1;
  else if(beta==0) sign = 0;

  i=floor((2-alpha)*10+1);
  if (i<1) i=1;
  else if (i>15) i=15;
  j=floor(beta/0.25+1);
  if (j<1) j=1;
  else if (j>4) j=4;

  t=beta/0.25-j+1;
  s=(2-alpha)/0.1-i+1;

    *c = enc[i-1][j-1] * (1-s) * (1-t) + enc[i][j-1] * s * (1-t) +
         enc[i-1][j] * t * (1-s) + enc[i][j] * t * s;
    *c = cn/(*c);
    *zeta = za[i-1][j-1] * (1-s) * (1-t) + za[i][j-1] * s * (1-t) +
              za[i-1][j] * t * (1-s) + za[i][j] * t * s;
    *zeta = q50 + *c * sign * (*zeta);


  return;
}
