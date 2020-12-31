/* stable/stable_pdf.c
 * 
 * Code for computing the PDF of an alpha-estable distribution.
 * Expresions presented in [1] are employed.
 *
 * [1] Nolan, J. P. Numerical Calculation of Stable Densities and
 *     Distribution Functions Stochastic Models, 1997, 13, 759-774
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
#include "stable.h"
#include "stable_integration.h"

#include "methods.h"
#include <pthread.h>

double stable_pdf_g1(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;
  double g, V, aux;

//  g   = dist->beta_;
//  aux = theta+dist->theta0_;
//  V   = M_PI_2-theta;

//  if ((g==1 && aux < THETA_TH*1.1) || (g==-1 && V<THETA_TH*1.1)) {
//    V = dist->Vbeta1;// printf("");
//  }
//  else {
    aux = (dist->beta_*theta+M_PI_2)/cos(theta);
    V = sin(theta)*aux/dist->beta_ + log(aux) + dist->k1;
//  }
  #ifdef DEBUG
  integ_eval++;
  #endif

  g = V + dist->xxipow;
  //Actually we compute log(g)
  //Taylor expansion: exp(-x) ~ 1-x for x ~ 0 
  //If g<1.52e-8 -> exp(-g)=(1-g) -> g·exp(-g) = g·(1-g) with double precision
  if (isnan(g)) return 0.0;
  if((g=exp(g)) < 1.522e-8 ) return (1.0-g)*g;
  g = exp(-g)*g;
  if (isnan(g) || g<0) return 0.0;

  return g;
}

double stable_pdf_g2(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;
  double g, cos_theta,aux,V;

  cos_theta = cos(theta);
  aux = (dist->theta0_+theta)*dist->alpha;
  V = log(cos_theta/sin(aux))*dist->alphainvalpha1 +
       + log(cos(aux-theta)/cos_theta) + dist->k1;

  #ifdef DEBUG
  integ_eval++;
  #endif

  g = V + dist->xxipow;
  // g>6.55 -> exp(g-exp(g)) < 2.1E-301
  if(g>6.55 || g<-700) return 0.0;
  else  g=exp(g);
  g = exp(-g)*g;
  if (isnan(g) || isinf(g) || g<0) {return 0.0;}
  return g;
}

double stable_pdf_g(double theta, void *args)
{
  StableDist * dist = (StableDist *)args;
  if (dist->ZONE == ALPHA_1) {
    return stable_pdf_g1(theta, args); }
  else if (dist->ZONE == CAUCHY) {
    return -1.0; }
  else {
    return stable_pdf_g2(theta, args); }
}

double stable_g_aux1(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;
  double g, V,aux;

  aux = (dist->beta_*theta+M_PI_2)/cos(theta);
  V = sin(theta)*aux/dist->beta_ + log(aux) + dist->k1;
  g = V + dist->xxipow;

  #ifdef DEBUG
  integ_eval++;
  #endif

  return g;
}

double stable_g_aux2(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;
  double g, cos_theta,aux,V;

  cos_theta = cos(theta);
  aux = (dist->theta0_+theta)*dist->alpha;
  V = log(cos_theta/sin(aux))*dist->alphainvalpha1 +
       + log(cos(aux-theta)/cos_theta) + dist->k1;

  g = V + dist->xxipow;

  #ifdef DEBUG
  integ_eval++;
  #endif

  return g;
}

double stable_g_aux(double theta, void *args)
{
  StableDist *dist = (StableDist *)args;

  if (dist->ZONE==ALPHA_1)
    return stable_g_aux1(theta,args);
  else
    return stable_g_aux2(theta,args);
}

void * thread_init_pdf(void *ptr_args)
{
  StableArgsPdf *args = (StableArgsPdf *)ptr_args;
  int counter_ = 0;

  while (counter_ < args->Nx)
    {
      args->pdf[counter_]=(*(args->ptr_funcion))(args->dist,args->x[counter_],
                                                 &(args->err[counter_]));
      counter_++;
    }
    
  pthread_exit(NULL);
  
  return NULL;
}

void stable_pdf(StableDist *dist, const double* x, const unsigned int Nx,
                double *pdf, double *err)
{
  int Nx_thread[THREADS],
      initpoint[THREADS],
      k,flag=0;
  void *status;
  pthread_t threads[THREADS];
  StableArgsPdf args[THREADS];
  
  if (err==NULL) {flag=1;err=(double*)malloc(Nx*sizeof(double));}

  // Evaluation points are distributed between the threads
  Nx_thread[0] = Nx/THREADS;
  if (0 < Nx%THREADS) Nx_thread[0]++;

  initpoint[0] = 0;
  for(k=1;k<THREADS;k++)
    {
      Nx_thread[k] = Nx/THREADS;
      if (k < Nx%THREADS) Nx_thread[k]++;
      initpoint[k] = initpoint[k-1] + Nx_thread[k-1];
    }

  // Threads are created with a copy of the distribution
  for(k=0; k<THREADS; k++)
    {
      args[k].ptr_funcion = dist->stable_pdf_point;

      args[k].dist = stable_copy(dist);
      args[k].pdf  = pdf+initpoint[k];
      args[k].x    = x+initpoint[k];
      args[k].Nx   = Nx_thread[k];
      args[k].err  = err+initpoint[k];

      if(pthread_create(&threads[k], NULL, thread_init_pdf, (void *)&args[k]))
        {
          perror("Error en la creacion de hilo");
          if (flag==1) free(err);
          return;
        }
    }

  // Wait until every thread finishes
  for(k=0; k<THREADS; k++)
    {
      pthread_join(threads[k], &status);
    }

  // Free distribution copies
  for(k=0; k<THREADS; k++)
    {
      stable_free(args[k].dist);
    }
    
  if (flag==1) free(err);
}

/******************************************************************************/
/*   PDF computation quadrature method                                        */
/******************************************************************************/

double
stable_integration_pdf_low(StableDist *dist, double(*integrand)(double,void*),
                  double(*integ_aux)(double,void*), double *err)
/* low precision method: two integration intervals --one around the maximum and
   another one for the rest */
{
  int warnz[5],k;
  double pdf=0,
         pdf_aux=0,pdf1=0,/*pdf2=0.0,pdf3=0.0,*/
         err_aux=0;
  double theta[5];
  //int method_;

  #ifdef DEBUG
  int aux_eval=0;
  #endif

  theta[0] = -dist->theta0_+THETA_TH; warnz[0]=0;
  theta[4] = M_PI_2 - THETA_TH;
  theta[2] = zbrent(integ_aux,(void*)dist,theta[0],theta[4],
                         0.0,1e-6*(theta[4]-theta[0]),&k);

  switch (k)
    {
      case 0:   //Max en el interior del intervalo de integracion.
    // Create integrand sub-interval around the maximum
        if (theta[2]-theta[0] < theta[4]-theta[2]) {
          theta[2]=theta[2]*2.0-theta[0];
        }
        else
         {
           pdf_aux=theta[0];
           theta[0]=theta[4];
           theta[4]=pdf_aux;
           theta[2]=theta[2]*2.0-theta[0];
         }
        break;

      case -2: // Max found in the left border
        // This is possible if beta=+-1 and alpha<1
        pdf1=(integrand)(theta[0],(void*)dist);

        theta[2]=zbrent(integrand,(void*)dist,theta[0],theta[4],
                         pdf1*1e-6,1e-6*(theta[4]-theta[0]),&warnz[2]);
        break;

      case  -1: //Max found in the right corner
        // This is possible if beta=+-1 and alpha<1
        theta[1]=theta[4];
        theta[4]=theta[0];
        theta[0]=theta[1];

        pdf1=(integrand)(theta[0],(void*)dist);

        theta[2]=zbrent(integrand,(void*)dist,theta[4],theta[0],
                         pdf1*1e-6,1e-6*(theta[0]-theta[4]),&warnz[2]);
        break;

      default: // Never get here
        theta[1] = 0.5*(theta[4]-theta[2]);
        theta[3] = 0.5*(theta[2]+theta[0]);
        break;
    }

  #ifdef DEBUG
  aux_eval=integ_eval;
  integ_eval=0;
  #endif

  stable_integration(dist,integrand,theta[0],theta[2],
                             absTOL,relTOL,IT_MAX,
                             &pdf_aux,&err_aux,STABLE_QAG2);
  pdf=fabs(pdf_aux);
  *err=err_aux*err_aux;

  #ifdef DEBUG
  warnz[0]=integ_eval;
  integ_eval=0;
  #endif

  stable_integration(dist,integrand,theta[2],theta[4],
                             max(pdf*relTOL,absTOL)*0.5,relTOL,IT_MAX,
                             &pdf_aux,&err_aux,STABLE_QAG2);
  pdf+=fabs(pdf_aux);
  *err+=err_aux*err_aux;
  #ifdef DEBUG
  warnz[1]=integ_eval;
  integ_eval=0;
  #endif

  *err=sqrt(*err)/pdf;

  #ifdef DEBUG
  fprintf(FINTEG,"%+1.3e % 1.3e % 1.3e",x,pdf,*err);
  fprintf(FINTEG," %+1.3e %+1.3e %+1.3e %+1.3e %+1.3e",theta[0],theta[1],theta[2],theta[3],theta[4]);
  fprintf(FINTEG," % 1.3e % 1.3e % 1.3e % 1.3e", pdf1,pdf2,pdf3,fabs(pdf_aux));
  fprintf(FINTEG," %d %d %d %d %d %d\n",
          warnz[0],warnz[1],warnz[2],integ_eval,aux_eval,
          warnz[0]+warnz[1]+warnz[2]+integ_eval+aux_eval);
  Rprintf("abstols % 1.3e % 1.3e % 1.3e % 1.3e \n",absTOL,max(pdf1*relTOL,absTOL)*0.5,max((pdf2+pdf1)*relTOL,absTOL)*0.25,max((pdf3+pdf2+pdf1)*relTOL,absTOL)*0.25);
  #endif

  return pdf;
}

double
stable_integration_pdf(StableDist *dist, double(*integrand)(double,void*),
                  double(*integ_aux)(double,void*), double *err)
{
/* Get here if:
       x >> xi for alpha > 1
       x ~  xi for alpha < 1
       x >> 0  for alpha = 1 and beta < 0
       x << 0  for alpha = 1 and beta > 0 */

/* Strategy:
     - Find the maximun of the integrand: theta[2]
     - Find the point where the integrand goes bellow a threshold
     - 1 - Integrate in a symmetric interval around the maximum
     - 2 - Integrate the remaining portion above the threshold
     - 3 y 4 - Integrate below the threshold
     - Add 1 to 4*/

  int warnz[5],k;
  double pdf=0,
         pdf_aux=0,pdf1=0,pdf2=0.0,pdf3=0.0,
         err_aux=0;
  double theta[5];

  #ifdef DEBUG
  int aux_eval=0;
  #endif

  theta[0] = -dist->theta0_+THETA_TH; warnz[0]=0;
  theta[4] = M_PI_2 - THETA_TH;
  theta[2] = zbrent(integ_aux,(void*)dist,theta[0],theta[4],
                         0.0,1e-9*(theta[4]-theta[0]),&k);

  switch (k)
    {
      case 0:   // Maximum found in the interior of the integration interval
        // Find points where integrand goes below threshold
        pdf1 = (integ_aux)(theta[0],(void*)dist);
        pdf2 = (integ_aux)(theta[4],(void*)dist);

        if (fabs(AUX1)>fabs(pdf1)) {
          theta[1]=theta[0]+1e-2*(theta[2]-theta[0]);
        }
        else {
          theta[1] = zbrent(integ_aux,(void*)dist,theta[0],theta[2],
                        AUX1, 1e-9*(theta[2]-theta[0]),&warnz[1]);
		}

        if (fabs(AUX2)>fabs(pdf2)) {
          theta[3]=theta[4]-1e-2*(theta[4]-theta[2]);
        }
        else {
	 	  theta[3] = zbrent(integ_aux,(void*)dist,theta[2],theta[4],
                        AUX2, 1e-9*(theta[4]-theta[2]),&warnz[3]);
        }

    // Symmetric interval around the maximum
        if (theta[2]-theta[1] < theta[3]-theta[2]) {
          theta[2]=theta[2]*2.0-theta[1];
        }
        else
         {
           pdf_aux=theta[0];                 //            .
           theta[0]=theta[4];                //           /|\    .
           theta[4]=pdf_aux;                 //         /  | \   .
           pdf_aux=theta[3];                 //_______/ |  |  \___
           theta[3]=theta[1];                //4      3 2  |   1 0
           theta[1]=pdf_aux;
           theta[2]=theta[2]*2.0-theta[1];
         }
        break;

      case -2: // Max in left border
        theta[1]=theta[0];
        pdf1=(integrand)(theta[1],(void*)dist);

        theta[2]=zbrent(integrand,(void*)dist,theta[1],theta[4],
                         pdf1*1e-6,1e-9*(theta[4]-theta[1]),&warnz[2]);
        pdf1=stable_pdf_g(theta[2],(void*)dist);
        theta[3]=zbrent(integrand,(void*)dist,theta[2],theta[4],
                         pdf1*1e-6,1e-9*(theta[4]-theta[2]),&warnz[2]);
        pdf1=stable_pdf_g(theta[3],(void*)dist);

        break;

      case  -1: // Max in right border
        // puede pasar si beta=+-1 y alpha<1
        theta[1]=theta[4];
        theta[4]=theta[0];
        pdf1=(integrand)(theta[1],(void*)dist);

        theta[2]=zbrent(integrand,(void*)dist,theta[4],theta[1],
                         pdf1*1e-6,1e-9*(theta[1]-theta[4]),&warnz[2]);
        pdf1=stable_pdf_g(theta[2],(void*)dist);
        theta[3]=zbrent(integrand,(void*)dist,theta[4],theta[2],
                         pdf1*1e-6,1e-9*(theta[2]-theta[4]),&warnz[3]);
 
        theta[0]=theta[1];
        break;

      default: // Nunca llegara aqui
        theta[1] = 0.5*(theta[4]-theta[2]);
        theta[3] = 0.5*(theta[2]+theta[0]);
        break;
    }

  #ifdef DEBUG
  aux_eval=integ_eval;
  integ_eval=0;
  #endif

  stable_integration(dist,integrand,theta[1],theta[2],
                             absTOL,relTOL,IT_MAX,
                             &pdf_aux,&err_aux,METHOD1);
  pdf1=fabs(pdf_aux);
  *err=err_aux*err_aux;

  #ifdef DEBUG
  warnz[0]=integ_eval;
  integ_eval=0;
  #endif

  stable_integration(dist,integrand,theta[2],theta[3],
                             max(pdf1*relTOL,absTOL)*0.25,relTOL,IT_MAX,
                             &pdf_aux,&err_aux,METHOD2);
  pdf2=fabs(pdf_aux);
  *err+=err_aux*err_aux;
  #ifdef DEBUG
  warnz[1]=integ_eval;
  integ_eval=0;
  #endif

  stable_integration(dist,integrand,theta[3],theta[4],
                            max((pdf2+pdf1)*relTOL,absTOL)*0.25,relTOL,IT_MAX,
                            &pdf_aux,&err_aux,METHOD3);
  pdf3=fabs(pdf_aux);
  *err+=err_aux*err_aux;
  #ifdef DEBUG
  warnz[2]=integ_eval;
  integ_eval=0;
  #endif

  stable_integration(dist,integrand,theta[0],theta[1],
                             max((pdf3+pdf2+pdf1)*relTOL,absTOL)*0.25,relTOL,IT_MAX,
                             &pdf_aux,&err_aux,STABLE_QAG2);
  *err+=err_aux*err_aux;


  // Add from lowest to highest contribution (minimize numerical error)
  pdf=fabs(pdf_aux)+pdf3+pdf2+pdf1;
  *err=sqrt(*err)/pdf;

  #ifdef DEBUG
  fprintf(FINTEG,"%+1.3e % 1.3e % 1.3e",x,pdf,*err);
  fprintf(FINTEG," %+1.3e %+1.3e %+1.3e %+1.3e %+1.3e",theta[0],theta[1],theta[2],theta[3],theta[4]);
  fprintf(FINTEG," % 1.3e % 1.3e % 1.3e % 1.3e", pdf1,pdf2,pdf3,fabs(pdf_aux));
  fprintf(FINTEG," %d %d %d %d %d %d\n",
          warnz[0],warnz[1],warnz[2],integ_eval,aux_eval,
          warnz[0]+warnz[1]+warnz[2]+integ_eval+aux_eval);
  Rprintf("abstols % 1.3e % 1.3e % 1.3e % 1.3e \n",absTOL,max(pdf1*relTOL,absTOL)*0.5,max((pdf2+pdf1)*relTOL,absTOL)*0.25,max((pdf3+pdf2+pdf1)*relTOL,absTOL)*0.25);

  #endif

  return pdf;
}

/******************************************************************************/
/*   PDF of special cases                                               */
/******************************************************************************/

double
stable_pdf_point_GAUSS(StableDist *dist, const double x, double *err)
{
  double x_=(x-dist->mu_0)/dist->sigma;
  *err = 0.0;

  return 0.5*sqrt(M_1_PI)/dist->sigma*exp(-x_*x_*0.25);
}

double
stable_pdf_point_CAUCHY(StableDist *dist, const double x, double *err)
{
  double x_=(x-dist->mu_0)/dist->sigma;
  *err = 0.0;

  return M_1_PI/(1+x_*x_)/dist->sigma;
}   

double
stable_pdf_point_LEVY(StableDist *dist, const double x, double *err)
{
  double xxi=(x-dist->mu_0)/dist->sigma-dist->xi;
  *err=0.0;

  if (xxi>0 && dist->beta>0)
    return sqrt(dist->sigma*0.5*M_1_PI) *
           exp(-dist->sigma*0.5/(xxi*dist->sigma)) /
           pow(xxi*dist->sigma,1.5);
  else if (xxi<0 && dist->beta<0)
    return sqrt(dist->sigma*0.5*M_1_PI) *
           exp(-dist->sigma*0.5/(fabs(xxi)*dist->sigma)) /
           pow(fabs(xxi)*dist->sigma,1.5);
  else return 0.0;
}

/******************************************************************************/
/*   General PDF                                                     */
/******************************************************************************/

double
stable_pdf_point_PEQXXIP(StableDist *dist, const double x, double *err)
{
/* Get here if:
       x ~ xi  con alpha > 1
       x >> xi con alpha < 1
       x << 0  con alpha = 1 y beta < 0
       x >> 0  con alpha = 1 y beta > 0 */

/* Removed in final version */
    return 0.0;
}

double
stable_pdf_point_MEDXXIP(StableDist *dist, const double x, double *err)
{
  return 0.0;
}

double
stable_pdf_point_ALPHA_1(StableDist *dist, const double x, double *err)
{
  double pdf=0;
  double x_;//, xxi;

  #ifdef DEBUG
  integ_eval=0;
  #endif

  x_=(x-dist->mu_0)/dist->sigma;

  dist->beta_ = dist->beta;

  if (dist->beta < 0.0)
    {
      x_ = -x_;
      dist->beta_ = -dist->beta;
    }

  else {dist->beta_ = dist->beta;}

  dist->xxipow = (-M_PI*x_*dist->c2_part);

  pdf = stable_integration_pdf(dist,&stable_pdf_g1,&stable_g_aux1,err);
  pdf = dist->c2_part*pdf;
  return pdf/dist->sigma;
}

double stable_pdf_point_STABLE(StableDist *dist, const double x, double *err)
{
  double pdf=0;
  double x_, xxi;

  #ifdef DEBUG
  int aux_eval=0;
  integ_eval=0;
  #endif

  x_=(x-dist->mu_0)/dist->sigma;
  xxi=x_-dist->xi;

  if (fabs(xxi) <= XXI_TH)
    {
      *err=0;

      pdf = exp(gammaln(1.0+1.0/dist->alpha)) *
            cos(dist->theta0)/(M_PI*dist->S);

      #ifdef DEBUG
      Rprintf("Aproximando x a zeta para alpha = %f, beta = %f, zeta = %f : pdf = %f\n",
               dist->alpha,dist->beta,dist->xi,pdf);
            
      fprintf(FINTEG,"%1.3e\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%d\t%d\t%d\t%d\n",
              x,pdf,*err,0.0,0.0,0.0,aux_eval,0,aux_eval,aux_eval<<1);
      #endif
      return pdf/dist->sigma;
    }
  else if (xxi<0)   /* pdf(x<xi,a,b) = pdf(-x,a,-b)*/
    { 
      xxi=-xxi;
      dist->theta0_=-dist->theta0; /*theta0(a,-b)=-theta0(a,b)*/
      dist->beta_=-dist->beta;
    }
  else
    {
      dist->theta0_=dist->theta0;
      dist->beta_=dist->beta;
    }
  dist->xxipow=dist->alphainvalpha1*log(fabs(xxi));

  /* If theta0~=-PI/2 pdf = 0.0 */
  /* This is possible if beta=+-1, alpha < 1 and x < xi*/
  if (fabs(dist->theta0_+M_PI_2)<2*THETA_TH)
    {
     #ifdef DEBUG
      Rprintf("Intervalo de integracion nulo\n");
     #endif

     return 0.0;
    }
  
  pdf = stable_integration_pdf(dist,&stable_pdf_g2,&stable_g_aux2,err);

  pdf = dist->c2_part/xxi*pdf;

  return pdf/dist->sigma;
}

/******************************************************************************/
/*   General PDF point                                                      */
/******************************************************************************/

double
stable_pdf_point(StableDist *dist, const double x, double *err)
{
  double temp;

  if (err==NULL) err=&temp;

  return (dist->stable_pdf_point)(dist,x,err);
}
