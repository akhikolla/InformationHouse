/* stable/methods.c
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

#include "methods.h"
#include "stable.h"
//#include "stable_common.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FUNC(x) ((*func)(x,args))

//extern FILE * velogout;
/*----------------------------------------------------------------------------*/
/*                             INTERPOLATION                                  */
/*----------------------------------------------------------------------------*/

void polint(const double xa[], const double ya[], const int n, double x,
            double *y, double *dy)
{
	int i,m,ns=0;
	double den,dif,dift,ho,hp,w;
	double c[n],d[n];
	dif=fabs(x-xa[0]);
						// Here we ﬁnd the index ns of the closest table entry,
	for (i=0;i<n;i++) {
	    if ( (dift=fabs(x-xa[i])) < dif) {
	         ns=i;
	         dif=dift;
	    }
	                                //and initialize the tableau of c’s and d’s.
	    c[i]=ya[i];
	    d[i]=ya[i];
	}
	                                   //This is the initial approximation to y.
	*y=ya[ns--];
	                                   //For each column of the tableau,
	for (m=1;m<n;m++) {
                               //we loop over the current c’s and d’s and update
		for (i=0;i<n-m;i++) {
                                          //them.
	         ho=xa[i]-x;
	         hp=xa[i+m]-x;
	         w=c[i+1]-d[i];
	         if ( (den=ho-hp) == 0.0) return;
								//This error can occur only if two input xa’s
								//are (to within roundoff) identical.
	         den=w/den;
	                            //Here the c’s and d’s are updated.
	         d[i]=hp*den;
	         c[i]=ho*den;
	    }
	    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));

   /* After each column in the tableau is completed, we decide which correction,
   c or d, we want to add to our accumulating value of y, i.e., which path to
   take through the tableau—forking up or down. We do this in such a way as to
   take the most “straight line” route through the tableau to its apex, updating
   ns accordingly to keep track of where we are. This route keeps the partial
   approximations centered (insofar as possible) on the target x. The last dy
   added is thus the error indication. */
	}
}

void ratint(const double xa[],const double ya[], const int n, double x,
            double * y, double * dy)
/* Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine
returns a value of y and an accuracy estimate dy. The value returned is that of
the diagonal rational function, evaluated at x, which passes through the
n points (xa_i, ya_i), i = 1, 2, ... n*/

{
	int m,i,ns=0;
	double w,t,hh,h,dd;
	double c[n];
	double d[n];
	hh=fabs(x-xa[0]);
	
	for (i=0;i<n;i++) {	
		h=fabs(x-xa[i]);
		if (h == 0.0) {
			*y=ya[i];
			*dy=(double)0.0;
			return;
		} else if (h < hh) {
			ns=i;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]+TINY;
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<n-m;i++) {
			w=c[i+1]-d[i];
			h=xa[i+m]-x;
			t=(xa[i]-x)*d[i]/h;
			dd=t-c[i+1];
			if (dd == 0.0) return;
			dd=w/dd;
			d[i]=c[i+1]*dd;
			c[i]=t*dd;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	return;
}


/*----------------------------------------------------------------------------*/
/*                               INTEGRATION                                  */
/*----------------------------------------------------------------------------*/

double trapzd(double (*func)(double, void *), void * args,
              double a, double b, int n, double s)
/* This routine computes the nth stage of refinement of an extended trapezoidal
rule. Func is input as a pointer to the function to be integrated between limits
a and b, also input. When called with n=1, the routine returns the crudest
estimate of int(a,b,f(x)dx). Subsequent calls with n=2,3,... (in that sequential
 order) will improve the accuracy by adding 2n-2 additional interior points. */
{
     double x,tnm,sum,del;
   //static double s;          //No puedo usar variables estaticas en multihilo
     int it,j;
     if (n == 1) {
          return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
     } else {
          it = 1<<(n-2);
          tnm=it;
          del=(b-a)/tnm;     //This is the spacing of the points to be added.
          x=a+0.5*del;
          //for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
          for (sum=0.0,j=1;j<=it;x=a+0.5*del+del*j++) {
            sum += FUNC(x);
          }
          s=0.5*(s+(b-a)*sum/tnm);  //This replaces s by its refined value.
          return s;
     }
}

/*Este esta extraido del MATLAB. Se trata de un metodo adaptativo de biseccion*/
double quadstep(double (*func)(double,void *),void *args,
                double a, double b,double fa, double fc, double fb,
                const double epsabs, const double epsrel,
                double *abserr, int *warn, size_t *fcnt)
{
    double h=b-a;
    double c=0.5*(a+b);
    int aux;
    size_t aux2;

    if (warn==NULL) warn=&aux;
    if (fcnt==NULL) fcnt=&aux2;

    if (fabs(h) < EPS || c==a || c==b) {*warn=1;return (h*fc);}
	
    double d=0.5*(a+c);
    double e=0.5*(c+b);
    double fd=FUNC(d);
    double fe=FUNC(e);
    *fcnt+=2;

  #ifdef DEBUG
    if(isnan(fa)) {Rprintf("a");}
    if(isnan(fd)) {Rprintf("d");}
    if(isnan(fc)) {Rprintf("c");}
    if(isnan(fe)) {Rprintf("e");}
    if(isnan(fb)) {Rprintf("b");}
  #endif

    double Q1=(h/6.0)*(fa+4.0*fc+fb);
    double Q2=(h/12.0)*(fa+4.0*fd+2.0*fc+4.0*fe+fb);
    double Q = Q2+(Q2-Q1)/15.0;
    *abserr = fabs(Q2-Q);

    if (isnan(Q)) {*warn=3;return (h*fc);}
    if ( *abserr <= fabs(Q)*epsrel || *abserr <= epsabs ) {*warn=0;return(Q);}
    if (*fcnt > 10000) {*warn=2;return (Q);}
			
    int warnac=0, warncb=0;
    double abserrac=0.0, abserrcb=0.0;
    double Qac = quadstep (func,args,a,c,fa,fd,fc,epsabs,epsrel,&abserrac,&warnac,fcnt);
    double Qcb = quadstep (func,args,c,b,fc,fe,fb,epsabs,epsrel,&abserrcb,&warncb,fcnt);
    *warn=((warnac > warncb) ? (warnac) : (warncb) );
    *abserr=sqrt(abserrac*abserrac+abserrcb*abserrcb);
    //fprintf(logout,"\tfcnt %5d, Qac = %1.6e, Qcb = %1.6e, Qac+Qcb = %1.6e.\n",
    //        *fcnt,Qac,Qcb,Qac+Qcb);

    return (Qac+Qcb);
}

#define JMAXP (JMAX+1)
/* Here EPS is the fractional accuracy desired, as determined by the
extrapolation error estimate; JMAX limits the total number of steps; K is the
number of points used in the extrapolation.*/
double qromb(double (*func)(double, void *),void *args,
             double a, double b, double epsabs, double epsrel, int K, int JMAX,
             int method, int *warn, size_t *fcnt, double *err)
/* Returns the integral of the function func from a to b. Integration is
performed by Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule.*/
{
  double ss,dss;
                             //These store the successive trapezoidal approxi-
  double s[JMAXP],h[JMAXP+1];//mations and their relative stepsizes.
  int j;
  void (*extrapol)(const double *, const double *, const int, double,
                   double *, double *);

  int aux;
  size_t aux2;

  if (warn==NULL) warn=&aux;
  if (fcnt==NULL) fcnt=&aux2;

  if (method==1) {extrapol = &polint;}
  else if (method==2) {extrapol = &ratint;}
  else {perror("\nERROR\n"); return(-1);}

  h[1]=1.0;
  s[0]=0.0;
  *fcnt += 2; //La primera vez evaluará en los extremos del intervalo.
  for (j=1;j<=JMAX;j++)
    {
	  s[j]=trapzd(func,args,a,b,j,s[j-1]);
	  //printf("s[%2d]=%1.6e\t",j,s[j]);
	  if (j>1) *fcnt += (1 << (j-2));
	  if (j >= K)
        {
          (*extrapol)(&h[j-K+1],&s[j-K+1],K,0.0,&ss,&dss);
          if (fabs(dss) <= epsrel*fabs(ss) || fabs(dss) <= epsabs) {//
            *warn=0;
            *err=fabs(dss);
            return ss;
	        }
        }
      h[j+1]=0.25*h[j];
      /*This is a key step: The factor is 0.25 even though the stepsize is
      decreased by only 0.5. This makes the extrapolation a polynomial in h2 as
      allowed by equation (4.2.1), not just a polynomial in h.*/
	}
	*warn=2; // si llego aqui me he pasao de iters
	*err=fabs(dss);
	return ss;
}

/*----------------------------------------------------------------------------*/
/*                              ROOT FINDING                                  */
/*----------------------------------------------------------------------------*/
#define ITMAX 200					// Maximum allowed number of iterations.
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double zbrent(double (*func)(double,void *),void * args,
              double x1, double x2, double value, double tol, int *warn)
/*
Using Brent’s method, find the root of a function func known to lie between x1
and x2. The root, returned as zbrent, will be refined until its accuracy is tol.
*/
{
  int iter;
  double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
  double fa=((*func)(a,args))-value,fb=((*func)(b,args))-value;
  double fc,p,q,r,s,tol1,xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    {
      //nrerror("Root must be bracketed in zbrent");
    #ifdef DEBUG
      Rprintf("f(a=%1.16lf)=%e, f(b=%1.16lf)=%1.16e, value=%1.16e, Root must be bracketed in zbrent\n",
              x1,fa+value,x2,fb+value,value);
    #endif
    if( fabs(fa) < fabs(fb))
      {
        *warn=-2;
        return a;
      }
    else
      {
        *warn=-1;
        return b;
      }
    }

  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
             c=a;         //Rename a, b, c and adjust bounding interval d.
             fc=fa;                                   
             e=d=b-a;
         }
         if (fabs(fc) < fabs(fb)) {
             a=b;
             b=c;
             c=a;
             fa=fb;
             fb=fc;
             fc=fa;
         }
                // Convergence check.
         tol1=2.0*EPS*fabs(b)+0.5*tol;
         xm=0.5*(c-b);
         if (fabs(xm) <= tol1 || fb == 0.0)
           {
             *warn=0;
             return b;
           }
         if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
                  // Attempt inverse quadratic interpolation.
             s=fb/fa;
             if (a == c) {
                  p=2.0*xm*s;
                  q=1.0-s;
             } else {
                  q=fa/fc;
                  r=fb/fc;
                  p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                  q=(q-1.0)*(r-1.0)*(s-1.0);
             }
                                                //  Check whether in bounds.
             if (p > 0.0) q = -q;
             p=fabs(p);
             min1=3.0*xm*q-fabs(tol1*q);
             min2=fabs(e*q);
             if (2.0*p < (min1 < min2 ? min1 : min2)) {
                                                 // Accept interpolation.
                  e=d;
                  d=p/q;
             } else {                  // Interpolation failed, use bisection.
                  d=xm;
                  e=d;
             }
         } else {				// Bounds decreasing too slowly, use bisection.
             d=xm;
             e=d;
      }
      a=b;		// Move last best guess to a.
      fa=fb;
                                    
      	if (fabs(d) > tol1)
          	b += d;  // Evaluate new trial root.
	else
		b += SIGN(tol1,xm);
      	fb=((*func)(b,args))-value;
  }
  #ifdef DEBUG
  Rprintf("Maximum number of iterations exceeded in zbrent\n");
  #endif
  *warn=-3;
	return 1.0/0.0;//Never get here.
  return 0.0;
}

/* ------------------------------------------
             FUNCTION LN(GAMMA(X))
   ------------------------------------------ */

double gammaln(const double xx)
{
	int j;
	double x, y, tmp, ser;
	static const double cof[6]= {76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
		-0.5395239384953e-5};
	
	y=x=xx;
	tmp=x+5.5;
	tmp-=(x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++) ser+= cof[j]/++y;

	return -tmp+log(2.5066282746310005*ser/x);

}

#define CON 1.4                //Stepsize is decreased by CON at each iteration.
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10        // Sets maximum size of tableau.
#define SAFE 2.0       // Return when error is SAFE worse than the best so far.

#define FMAX(a,b) (a>b ? a : b)



double dfridr(double (*func)(double, void * args), void * args,
              double x, double h, double *err)
/* Returns the derivative of a function func at a point x by Ridders’ method of
polynomial extrapolation. The value h is input as an estimated initial stepsize;
it need not be small, but rather should be an increment in x over which func
changes substantially. An estimate of the error in the derivative
is returned as err.*/
{
 int i,j;
 double errt,fac,hh,*a,ans=0;

 if (h == 0.0) {perror("h must be nonzero in dfridr.");return HUGE_VAL;}
 a=(double*)malloc(sizeof(double)*NTAB*NTAB);
 hh=h;
 a[0+NTAB*0]=((*func)(x+hh,args)-(*func)(x-hh,args))/(2.0*hh);
 *err=BIG;
 for (i=1;i<NTAB;i++) {
// Successive columns in the Neville tableau will go to smaller stepsizes and
// higher orders of extrapolation.
  hh /= CON;   // Try new, smaller stepsize
  a[0+NTAB*i]=((*func)(x+hh,args)-(*func)(x-hh,args))/(2.0*hh);
  fac=CON2;
//Compute extrapolations of various orders,
//requiring no new function evaluations.
  for (j=1;j<=i;j++) {
   a[j+NTAB*i]=(a[j-1+NTAB*i]*fac-a[j-1+NTAB*(i-1)])/(fac-1.0);
   fac=CON2*fac;
   errt=FMAX(fabs(a[j+NTAB*i]-a[j-1+NTAB*i]),
        fabs(a[j+NTAB*i]-a[j-1+NTAB*(i-1)]));
//The error strategy is to compare each new extrapolation to one order lower,
//both at the present stepsize and the previous one.
//If error is decreased, save the improved answer.
   if (errt <= *err) {
    *err=errt;
    ans=a[j+NTAB*i];
   }
  }
  if (fabs(a[i+NTAB*i]-a[i-1+NTAB*(i-1)]) >= SAFE*(*err)) break;
  // If higher order is worse by a significant factor SAFE, then quit early.
 }

 free(a);
 return ans;
}



void vector_step(double **x, double min, double max, double step, int * n)
{
  int i,m;
  double aux;

  aux=(max-min)/step;

  if (aux<0 || isnan(aux) || isinf(aux))
    {
      *n=0;
      (*x)=NULL;
      perror("Warning: Empty vector");
      return;
    }

  m=(int)(aux)+1;
  (*x)=(double*)malloc(m*sizeof(double));

  if ((*x) == NULL)
    {
      perror("Error while creating x array");
      return;
    }

  for(i=0;i<m;i++)
    {
      (*x)[i]=min+i*step;
    }

  *n=m;

  return;
}


void vector_npoints(double **x, double min, double max, int n, double * step)
{
	int i;

	*step=(max-min)/((double)(n)-1.0);

	(*x)=(double*)malloc(n*sizeof(double));
	if ((*x) == NULL)
	{
		perror("Error while creating x array");
		return;
	}

	for(i=0;i<n;i++)
	{
		(*x)[i]=min+i*(*step);
	}

	return;
}
