#include <Rcpp.h>
#include <R_ext/BLAS.h>
using namespace Rcpp;

#define ILENGTH(p) (7+(p))
#define DLENGTH(p) (5*(p)+(p)*(p))

namespace {
    // Interface
    enum {lar, lasso, forward, larsols, larsxscale, larswscale, larsnoscale};
    // Initialization of ilars and dlars opt=(n,p,scale,meth,est) 
    // ilars=(p,act,drop,ign,scale,meth,est,order) length=7+p 
    // dlars=(b,cov,d,r,u,a) length=3*p+p*p+2*p 
    void gglarsinit(int n, int p, double *x, double *y, int scale, double *wx,
		    int meth, int est, int *ilars, double *dlars);
    void ggsetilars(int p, int meth, int est, int *ilars);
    void ggsetxty(int n, int p, double *x, double *y, double *dlars);
    void ggsetr(int n, int p, double *x, double *dlars);
    void ggsetscale(int scale, double *wx, int *ilars, double *dlars);
    // Next beta in lars path is returned in beta; number of non zero elements in nbeta
    // if the return is false we have reached the end of the lars path
    bool gglarsnext(int *ilars, double *dlars, double *beta, int *nbeta);

    enum {p_,act_,drop_,ign_,scale_,meth_,est_,order_};

    // Constants
    const double EPS2=std::numeric_limits<double>::epsilon(), EPS=sqrt(EPS2);

    // Utilities    
    void ggzeros(int n, double *y);
    void ggidamax(int p, double *x, int *id, double *xmax);
    void gggivens(int n, double *x, int incx, double *y, int incy);
    void ggrswap(int which, int p, double *r);
    void gglarsmove(int from, int to, int *ilars, double *dlars);
    void gglarsadd(int *ilars, double *dlars);
    void gglarsdrop(int *ilars, double *dlars);
    void gglarsdir(int *ilars, double *dlars);
    double gglarsstep(int *ilars, double *dlars);
}



// [[Rcpp::export]]
LogicalVector alasso(NumericMatrix x, NumericVector y, double gamma, double P) {
    bool go=TRUE;
    int i, n=x.nrow(), p=x.ncol(), n2=n/2, nb, nbopt=0, ione=1;
    double bi, crit, critopt, k=log(static_cast<double>(n));
    IntegerVector ifor(ILENGTH(p));
    NumericVector b(p), bopt(p), r(n), dfor(DLENGTH(p));
    gglarsinit(x.nrow(),p,x.begin(),y.begin(),larsnoscale,NULL,forward,forward,
	       ifor.begin(),dfor.begin());
    IntegerVector ilasso=clone(ifor);
    NumericVector dlasso=clone(dfor);
    while (gglarsnext(ifor.begin(),dfor.begin(),b.begin(),&nb) && (nb<n2));
    ggsetilars(p,lasso,larsols,ilasso.begin());
    ggsetscale(larswscale,b.begin(),ilasso.begin(),dlasso.begin());
    do {
	go = gglarsnext(ilasso.begin(),dlasso.begin(),b.begin(),&nb);
	std::copy(y.begin(),y.end(),r.begin());
	for (i=0; i<p; i++) {
	    if (fabs(b[i])>EPS2) {
		bi = -b[i];
		F77_CALL(daxpy)(&n,&bi,x.begin()+n*i,&ione,r.begin(),&ione);
	    }
	}
	crit = n*log(mean(r*r))+k*nb+2*gamma*R::lchoose(P,nb);
	if ((nbopt==0)||(crit<critopt)) {
	    critopt = crit;
	    nbopt = nb;
	    std::copy(b.begin(),b.end(),bopt.begin());
	}
    } while(go);
    return abs(bopt)>EPS;
}



    



namespace {
    
    void gglarsinit(int n, int p, double *x, double *y, int scale, double *wx,
		    int meth, int est, int *ilars, double *dlars) {
	ggsetilars(p,meth,est,ilars);
	ggzeros(3*p+p*p,dlars) ;
	ggsetxty(n,p,x,y,dlars);
	ggsetr(n,p,x,dlars);
	ggsetscale(scale,wx,ilars,dlars);
    }
    

    void ggsetilars(int p, int meth, int est, int *ilars) {
	int i, *order=ilars+order_;
	ilars[p_]=p ; ilars[act_]=0 ; ilars[drop_]= -1; ilars[ign_]=0;
	ilars[meth_]=meth; ilars[est_]=est ;
	for ( i=0 ; i<p ; i++) order[i] = i ;
    }

    void ggsetxty(int n, int p, double *x, double *y, double *dlars) {
	char tran='T' ;
	int ione=1;
	double *cov=dlars+p, one=1.0 ,zero=0.0;
	/* cov=X'y */
	F77_CALL(dgemv)(&tran,&n,&p,&one,x,&n,y,&ione,&zero,cov,&ione);	
    }

    void ggsetr(int n, int p, double *x, double *dlars) {
	/* computation of the R factor of the QR decomposition of X in r
	   we make, row by row,  a copy x so that x is unchanged on output */
	int i, j, ione=1;
	double *r=dlars+3*p, *u=r+p*p;
	for ( i=0 ; i<n ; i++) {
	    F77_CALL(dcopy)(&p, x+i, &n, u, &ione);
	    for (j=0 ; j<p ; j++) gggivens(p-j, r+j*p+j, p, u+j, ione);
	}	
    }

    void ggsetscale(int scale, double *wx, int *ilars, double *dlars) {
	int i, p=ilars[p_], pi, ione=1;
	double *cov=dlars+p, *d=cov+p, *r=d+p, g;
	/* d is used to rescale the column of x  */
	ilars[scale_] = scale;
	for ( i=0, pi=1 ; i<p ; i++, pi++, r+=p ) {
	    g = F77_CALL(dnrm2)(&pi,r,&ione);
	    if ((g<EPS) || ((scale==larswscale) && (fabs(wx[i])<EPS))) {
		d[i] = 1.0 ; cov[i] = 0.0 ; ggzeros(pi,r);
	    } else if (scale != larsnoscale) {
		g = (scale==larsxscale) ? (1.0/g) : wx[i] ;
		d[i] = g ; cov[i] *= g ; F77_CALL(dscal)(&pi,&g,r,&ione);
	    } else {
		d[i] = 1.0 ;
	    }
	}
    }
    
    bool gglarsnext(int *ilars, double *dlars, double *beta, int *nbeta){
	int i, ione=1, p=ilars[p_], *order=ilars+order_;
	double g, *b=dlars, *cov=b+p, *d=cov+p, *r=d+p, *u=r+p*p, *a=u+p,
	    eps = 1000*EPS2;
	/* if last iteration wasn't a drop, add a variable to the active set */
	if (ilars[drop_]<0) gglarsadd(ilars,dlars); 
	/* compute step directions; u for beta; a for cov */
	gglarsdir(ilars, dlars);
	/* compute steplength  */
	g = gglarsstep(ilars, dlars) ;
	/* do the step */
	F77_CALL(daxpy)(ilars+act_, &g, u, &ione, b, &ione) ;
	g = -g ; F77_CALL(daxpy)(&p, &g, a, &ione, cov, &ione) ;
	/* if needed, drop a variable from the active set */
	if (ilars[drop_]>=0) gglarsdrop(ilars,dlars);
	/* output */
	*nbeta = ilars[act_];
	g += 1 ; ggzeros(p, beta); 
	if ((ilars[est_]!=larsols) || (g<eps)) {
	    for (i=0 ; i<ilars[act_] ; i++) beta[order[i]]= d[i]*b[i] ;
	} else {
	    for (i=0 ; i<ilars[act_] ; i++) beta[order[i]]= d[i]*(b[i]+g*u[i]) ;
	}
	/* stopping condition */
	ggidamax(p,cov,&i,&g);
	return (((g<eps) || ((ilars[act_]+ilars[ign_])==p)))? false : true;
    }


    void ggzeros(int n, double *y) {
	int i;
	for ( i=0 ; i<n ; i++) y[i] = 0.0 ;
    }

    void ggidamax(int p, double *x, int *id, double *xmax) {
	int i, iid;
	double d, cmax;
	for (iid=0, cmax=fabs(x[iid]), i=1; i<p; i++) {
	    d = fabs(x[i]);
	    if (d>cmax) {
		iid = i; cmax=d;
	    } 
	}
	*id = iid; *xmax = cmax;
    }


    /* apply givens rotation to x and y zeroing the first elements of y */
    void gggivens(int n, double *x, int incx, double *y, int incy) {
	double t1=x[0], t2=y[0], g1, g2;
	if (fabs(t2) > 0.0) {
	    F77_CALL(drotg)(&t1, &t2, &g1, &g2);
	    F77_CALL(drot)(&n, x, &incx, y, &incy, &g1, &g2);    
	}
    }


    void ggrswap(int which, int p, double *r) {
	int len=p-which, ione=1;
	double *col=r+which*p , *x=col+which, *y=x+1;
	F77_CALL(dswap)(&p, col, &ione, col+p, &ione);
	gggivens(len, x, p, y, p) ;
	y[0] = 0.0 ;
    } 

#define SWAP(V,i,T) {T=V[i];V[i]=V[(i)+1];V[(i)+1]=T;}
    void gglarsmove(int from, int to, int *ilars, double *dlars) {
	int step, itmp, p=ilars[p_], *order=ilars+order_;
	double tmp, *b=dlars, *cov=b+p, *d=cov+p, *r=d+p;
	if (from==to) {
	    return;
	} else if (from < to ) {
	    step = 1 ;
	} else {
	    from--; to--; step = -1 ;
	}
	while (from != to ) {
	    SWAP(order,from,itmp);
	    SWAP(b,from,tmp);
	    SWAP(cov,from,tmp);
	    SWAP(d,from,tmp);
	    ggrswap(from,p,r);
	    from += step;
	}
    }
#undef SWAP


    void gglarsdrop(int *ilars, double *dlars) {
	int i, p=ilars[p_]; 
	double *b=dlars, *cov=b+p, *d=cov+p, *r=d+p , *u=r+p*p;
	ilars[ign_] = 0 ;
	ilars[act_] -= 1 ;
	gglarsmove(ilars[drop_],ilars[act_],ilars,dlars);
	b[ilars[act_]] = u[ilars[act_]] = 0.0 ;
    }

    void gglarsadd(int *ilars, double *dlars) {
	bool add=false;
	int i, p=ilars[p_] , act=ilars[act_], ign=ilars[ign_], id, ncov, ione=1;
	double cmax, *b=dlars, *cov=b+p, *d=cov+p, *r=d+p;
	while (!add && (act+ign<p) ) {
	    ggidamax(p-act-ign,cov+act+ign,&id,&cmax);
	    id += act+ign; cmax *= (1-EPS);
	    for ( i=act+ign ; i<p ; i++) {
		if ((i==id) || (fabs(cov[i])>=cmax)) {
		    gglarsmove(i,ilars[act_],ilars,dlars);
		    if (fabs(r[ilars[act_]*(p+1)])>EPS) {
			ilars[act_] += 1 ;
			add = true;
		    } else {
			ign++;
		    } 
		}
	    }
	}
	ilars[ign_] = ign ;
    }


    void gglarsdir(int *ilars, double *dlars) {
	int i, p=ilars[p_], act=ilars[act_], nact=p-act, ione=1;
	char no='N', tran='T', up='U';
	double *b=dlars, *cov=b+p, *d=cov+p, *r=d+p, *u=r+p*p, *a=u+p,
	    cmax, z=0.0, o=1.0, *r12=r+act*p, *a2=a+act ;
	F77_CALL(dcopy)(&act,cov,&ione,u,&ione);
	F77_CALL(dcopy)(&act,cov,&ione,a,&ione);
	F77_CALL(dtrsv)(&up, &tran, &no, &act, r, &p, u, &ione) ;
	if (nact) 
	    F77_CALL(dgemv)(&tran, &act, &nact, &o, r12, &p, u, &ione, &z, a2, &ione); 
	F77_CALL(dtrsv)(&up, &no, &no, &act, r, &p, u, &ione) ;
    }


    double gglarsstep(int *ilars, double *dlars) {
	int i, p=ilars[p_], act=ilars[act_], f=act+ilars[ign_], method=ilars[meth_];
	double g, v, *b=dlars, *cov=b+p, *d=cov+p, *r=d+p, 
	    *u=r+p*p, *a=u+p, eps = 1000*EPS2, cmax=fabs(cov[0]);
	if (method==forward) return 1.0 ;
	for (i=f, g=1.0 ; i<p ; i++) {
	    v = (cmax-cov[i])/(cmax-a[i]) ;
	    if (v>eps) g = R::fmin2(g,v) ;
	    v = (cmax+cov[i])/(cmax+a[i]) ;
	    if (v>eps) g = R::fmin2(g,v) ;
	}
	if (method==lasso) {
	    ilars[drop_] = -1 ;
	    for (i=0 ; i<act ; i++) {
		v = -b[i]/u[i] ;
		if ((v>eps) && (v<g) ) {
		    g = v ;
		    ilars[drop_] = i ;
		}
	    }
	}
	return g;
    }
}
