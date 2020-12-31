#include <R.h>
#include<Rinternals.h>
#include <R_ext/Applic.h>
#include<Rmath.h>
typedef struct int_struct{
  SEXP f; /* function  */
  SEXP env; /* environment to evaluate the function */
} int_struct, *IntStruct;
static void Rintfn(double *x, int n, void *ex){
  SEXP args,resultsxp,tmp;
  int i;
  IntStruct IS=(IntStruct)ex;

  PROTECT(args=allocVector(REALSXP,n));
  for(i=0; i<n; i++)REAL(args)[i]=x[i];
  
  PROTECT(tmp=lang2(IS->f,args));
  PROTECT(resultsxp=eval(tmp,IS->env));

  if(length(resultsxp)!=n)
    error("evaluation of function gave a result of wrong length");
  if(TYPEOF(resultsxp)==INTSXP){
    resultsxp=coerceVector(resultsxp,REALSXP);
  }else if(TYPEOF(resultsxp)!=REALSXP)
    error("evaluation of function gave a result of wrong type");

  for(i=0; i<n; i++){
    x[i]=REAL(resultsxp)[i];
    if(!R_FINITE(x[i]))
      error("non-finite function value");
  }
  UNPROTECT(3);
  return;
}
SEXP mloglik0(SEXP args){
  int_struct is_nu,is_g;
  SEXP res,jtms,T,gcall,nucall,tmp_in,tmp_gout,tmp_nuout;/* ,tmp; */
  double lo, up, epsabs, epsrel, intgrv, abserr, *work,tmp;
  int neval,ier,limit,lenw,last,*iwork;
  int i,j,njmp;

  args=CDR(args);jtms=CAR(args);
  args=CDR(args);njmp=asInteger(CAR(args));
  args=CDR(args);T=CAR(args);
  args=CDR(args);is_nu.f=CAR(args);
  args=CDR(args);is_g.f=CAR(args);
  args=CDR(args);is_nu.env=CAR(args);is_g.env=CAR(args);
  args=CDR(args);epsabs = asReal(CAR(args));
  args=CDR(args);epsrel = asReal(CAR(args));
  args=CDR(args);limit = asInteger(CAR(args));
  lenw = 4 * limit;
  iwork = (int *) R_alloc((size_t) limit, sizeof(int));
  work = (double *) R_alloc((size_t) lenw, sizeof(double));
  PROTECT(res = allocVector(REALSXP, 1));

  /* njmp=length(jtms); */
  if(njmp>1){
    PROTECT(tmp_in=allocVector(REALSXP,njmp*(njmp-1)/2));
    for(i=1; i<=njmp-1; i++){
      for(j=0; j<i; j++){
	REAL(tmp_in)[i*(i-1)/2+j] = REAL(jtms)[i]-REAL(jtms)[j];
      }
    }
    PROTECT(gcall=lang2(is_g.f,tmp_in));
    PROTECT(tmp_gout=eval(gcall,is_g.env));
  
    PROTECT(nucall=lang2(is_nu.f,jtms));
    PROTECT(tmp_nuout=eval(nucall,is_nu.env));
    REAL(res)[0] = -log(REAL(tmp_nuout)[0]);
    for(i=1; i<=njmp-1; i++){
      tmp = REAL(tmp_nuout)[i];
      for(j=0;j<i;j++)tmp += REAL(tmp_gout)[i*(i-1)/2+j];
      if(tmp<=0){REAL(res)[0]=exp(400.0);UNPROTECT(6);return(res);}
      REAL(res)[0] -= log(tmp);
    }
  } else REAL(res)[0] -= 0.0;

  lo=0.0; up=REAL(T)[0];
  Rdqags(Rintfn, (void*)&is_nu,
  	 &lo, &up, &epsabs, &epsrel, &intgrv,
  	 &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  REAL(res)[0] += intgrv;
  if(ier>0)warning("warning:ier=%d>0 in integrating nu!\n",ier);

  if(njmp>0){
    for(i=0; i<njmp; i++){
      up=REAL(T)[0]-REAL(jtms)[i];
      Rdqags(Rintfn, (void*)&is_g,
	     &lo, &up, &epsabs, &epsrel, &intgrv,
	     &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
      REAL(res)[0] += intgrv;
      if(ier>0)warning("warning:ier=%d>0 in integrating g from 0 to 1-t_%d!\n",ier,i+1);
    }
  }

  UNPROTECT(njmp>0?6:1);
  return(res);
}
SEXP mloglik1(SEXP args){
  int_struct is_nu,is_g,is_Ig;
  SEXP res,jtms,T,gcall,nucall,Igcall,tmp_in,tmp_gout,tmp_nuout,Igout,Igin;/* ,tmp; */
  double lo, up, epsabs, epsrel, intgrv, abserr, *work,tmp;
  int neval,ier,limit,lenw,last,*iwork;
  int i,j,njmp;

  args=CDR(args);jtms=CAR(args);
  args=CDR(args);T=CAR(args);
  args=CDR(args);is_nu.f=CAR(args);
  args=CDR(args);is_g.f=CAR(args);
  args=CDR(args);is_Ig.f=CAR(args);
  args=CDR(args);is_nu.env=is_g.env=is_Ig.env=CAR(args);
  args=CDR(args);epsabs = asReal(CAR(args));
  args=CDR(args);epsrel = asReal(CAR(args));
  args=CDR(args);limit = asInteger(CAR(args));
  lenw = 4 * limit;
  iwork = (int *) R_alloc((size_t) limit, sizeof(int));
  work = (double *) R_alloc((size_t) lenw, sizeof(double));
  PROTECT(res = allocVector(REALSXP, 1));

  njmp=length(jtms);
  PROTECT(tmp_in=allocVector(REALSXP,njmp*(njmp-1)/2));
  for(i=1; i<=njmp-1; i++){
    for(j=0; j<i; j++){
      REAL(tmp_in)[i*(i-1)/2+j] = REAL(jtms)[i]-REAL(jtms)[j];
    }
  }
  PROTECT(gcall=lang2(is_g.f,tmp_in));
  PROTECT(tmp_gout=eval(gcall,is_g.env));
  
  PROTECT(nucall=lang2(is_nu.f,jtms));
  PROTECT(tmp_nuout=eval(nucall,is_nu.env));
  REAL(res)[0] = -log(REAL(tmp_nuout)[0]);
  for(i=1; i<=njmp-1; i++){
    tmp = REAL(tmp_nuout)[i];
    for(j=0;j<i;j++)tmp += REAL(tmp_gout)[i*(i-1)/2+j];
    REAL(res)[0] -= log(tmp);
  }

  lo=0.0; up=REAL(T)[0];
  Rdqags(Rintfn, (void*)&is_nu,
  	 &lo, &up, &epsabs, &epsrel, &intgrv,
  	 &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  REAL(res)[0] += intgrv;
  if(ier>0)warning("warning:ier=%d>0 in integrating nu!\n",ier);

  PROTECT(Igin=allocVector(REALSXP,njmp));
  for(i=0;i<njmp;i++)REAL(Igin)[i]=REAL(T)[0]-REAL(jtms)[i];
  PROTECT(Igcall=lang2(is_Ig.f,Igin));
  PROTECT(Igout=eval(Igcall,is_Ig.env));  
  for(i=0; i<njmp; i++)REAL(res)[0] += REAL(Igout)[i];

  UNPROTECT(9);
  return(res);
}

SEXP mloglik1a(SEXP args){
  int_struct is_nu,is_g,is_Ig;
  SEXP res,jtms,T,gcall,nucall,Igcall,tmp_in,tmp_gout,tmp_nuout,Igout,Igin;/* ,tmp; */
  double lo, up, epsabs, epsrel, intgrv, abserr, *work,tmp;
  int neval,ier,limit,lenw,last,*iwork;
  int i,j,njmp;

  args=CDR(args);jtms=CAR(args);
  args=CDR(args);njmp=asInteger(CAR(args));
  args=CDR(args);T=CAR(args);
  args=CDR(args);is_nu.f=CAR(args);
  args=CDR(args);is_g.f=CAR(args);
  args=CDR(args);is_Ig.f=CAR(args);
  args=CDR(args);is_nu.env=is_g.env=is_Ig.env=CAR(args);
  args=CDR(args);epsabs = asReal(CAR(args));
  args=CDR(args);epsrel = asReal(CAR(args));
  args=CDR(args);limit = asInteger(CAR(args));
  lenw = 4 * limit;
  iwork = (int *) R_alloc((size_t) limit, sizeof(int));
  work = (double *) R_alloc((size_t) lenw, sizeof(double));
  PROTECT(res = allocVector(REALSXP, 1));

  /* njmp=length(jtms); */
  if(njmp>0){
    PROTECT(tmp_in=allocVector(REALSXP,njmp*(njmp-1)/2));
    for(i=1; i<=njmp-1; i++){
      for(j=0; j<i; j++){
	REAL(tmp_in)[i*(i-1)/2+j] = REAL(jtms)[i]-REAL(jtms)[j];
      }
    }
    PROTECT(gcall=lang2(is_g.f,tmp_in));
    PROTECT(tmp_gout=eval(gcall,is_g.env));
    PROTECT(nucall=lang2(is_nu.f,jtms));
    PROTECT(tmp_nuout=eval(nucall,is_nu.env));
    REAL(res)[0] = -log(REAL(tmp_nuout)[0]);
    for(i=1; i<=njmp-1; i++){
      tmp = REAL(tmp_nuout)[i];
      for(j=0;j<i;j++)tmp += REAL(tmp_gout)[i*(i-1)/2+j];
      if(tmp<=0){REAL(res)[0]=exp(400.0);UNPROTECT(6);return(res);}
      REAL(res)[0] -= log(tmp);
    }
  }
  else{ REAL(res)[0] = 0.0; }

  lo=0.0; up=REAL(T)[0];
  Rdqags(Rintfn, (void*)&is_nu,
  	 &lo, &up, &epsabs, &epsrel, &intgrv,
  	 &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
  REAL(res)[0] += intgrv;
  if(ier>0)warning("warning:ier=%d>0 in integrating nu!\n",ier);
  if(njmp>0){
    PROTECT(Igin=allocVector(REALSXP,njmp));
    for(i=0;i<njmp;i++)REAL(Igin)[i]=REAL(T)[0]-REAL(jtms)[i];
    PROTECT(Igcall=lang2(is_Ig.f,Igin));
    PROTECT(Igout=eval(Igcall,is_Ig.env));  
    for(i=0; i<njmp; i++)REAL(res)[0] += REAL(Igout)[i];
  }
  UNPROTECT(njmp>0?9:1);
  return(res);
}

SEXP mloglik1b(SEXP args){
  int_struct is_nu,is_g,is_Ig,is_Inu;
  SEXP res,jtms,T,gcall,nucall,Igcall,Inucall,tmp_in,tmp_gout,tmp_nuout,Igout,Igin,Inu_out;/* ,tmp; */
  /* double lo, up, epsabs, epsrel, intgrv, abserr, *work,tmp; */
  double tmp;  
  /* int neval,ier,limit,lenw,last,*iwork;  */
  int i,j,njmp;

  args=CDR(args);jtms=CAR(args);
  args=CDR(args);njmp=asInteger(CAR(args));
  args=CDR(args);T=CAR(args);
  args=CDR(args);is_nu.f=CAR(args);
  args=CDR(args);is_g.f=CAR(args);
  args=CDR(args);is_Ig.f=CAR(args);
  args=CDR(args);is_Inu.f=CAR(args);  
  args=CDR(args);is_nu.env=is_g.env=is_Ig.env=is_Inu.env=CAR(args);

  PROTECT(res = allocVector(REALSXP, 1));

  /* njmp=length(jtms); */
  if(njmp>0){
    PROTECT(tmp_in=allocVector(REALSXP,njmp*(njmp-1)/2));
    for(i=1; i<=njmp-1; i++){
      for(j=0; j<i; j++){
	REAL(tmp_in)[i*(i-1)/2+j] = REAL(jtms)[i]-REAL(jtms)[j];
      }
    }
    PROTECT(gcall=lang2(is_g.f,tmp_in));
    PROTECT(tmp_gout=eval(gcall,is_g.env));
    PROTECT(nucall=lang2(is_nu.f,jtms));
    PROTECT(tmp_nuout=eval(nucall,is_nu.env));
    REAL(res)[0] = -log(REAL(tmp_nuout)[0]);
    for(i=1; i<=njmp-1; i++){
      tmp = REAL(tmp_nuout)[i];
      for(j=0;j<i;j++)tmp += REAL(tmp_gout)[i*(i-1)/2+j];
      if(tmp<=0){REAL(res)[0]=exp(400.0);UNPROTECT(6);return(res);}
      REAL(res)[0] -= log(tmp);
    }
  }
  else{ REAL(res)[0] = 0.0; }

  PROTECT(Inucall=lang2(is_Inu.f,T));
  PROTECT(Inu_out=eval(Inucall,is_Inu.env));

  REAL(res)[0] += REAL(Inu_out)[0];

  if(njmp>0){
    PROTECT(Igin=allocVector(REALSXP,njmp));
    for(i=0;i<njmp;i++)REAL(Igin)[i]=REAL(T)[0]-REAL(jtms)[i];
    PROTECT(Igcall=lang2(is_Ig.f,Igin));
    PROTECT(Igout=eval(Igcall,is_Ig.env));  
    for(i=0; i<njmp; i++)REAL(res)[0] += REAL(Igout)[i];
  }
  UNPROTECT(njmp>0?11:3);
  return(res);
}
