#include <Rcpp.h>
#include <R_ext/BLAS.h>
using namespace Rcpp;


namespace {
    //--------------------
    // Constants 
    //--------------------
    const double eps2=std::numeric_limits<double>::epsilon(), eps=sqrt(eps2), one=1, onem=-1;
    const int ione=1, MMedIter=100;
    //--------------------
    // Functions
    //--------------------
    // column of x are randomly permutated
    inline void ggmperm(int p, int n, double *x);
    // length(a)=n; set r=x to replace the values of vector with their ranks
    inline void ggrank(int n, double *x, double *r, int *a);
    // add one row to the Cholesky factor X'X
    inline void gggivens(int p, double *r, double *x);
    // update the mean and the cholesky factor of cov; length(w)=p
    inline void ggupdbr(int p, int n, double *xb, double *r, double *x, double *w);
    // subgroup sums
    inline void ggsubgroupsums(int p, int n, int m, double *x, double *sx);
    // x[,i] <- x[,i]-m
    inline void ggcenter(int p, int n, double *x, double *m);
    // y <- solve(t(r),x); it's safe to set y=x to transform x in place
    inline void ggorth(int p, int n, double *r, double *x, double *y);
    // spatial median;x is overwritten with the centered data;length(t)=max(n,p)
    inline void ggmmedinit(int p, int n, double *x, double *m, double *t);
    inline bool ggmmedstep(int p, int n, double *x, double *m, double *t);
    inline int ggmmed(int p, int n, double *x, double *m, double *t);
    // spatial median of the subgroup means + classic length(w)=pm+max(m,p)
    inline void ggmc(int p, int n, int m, double *x, double *l, double *r, double *w);
    // spatial median + cov(x[,i]-x[,i-1]); length(w)=max(m,p)
    inline void ggmc1(int p, int m, double *x, double *l, double *r, double *z, double *w);
    // length(iw)=n length(w)=2*n 
    inline void ggsignedrank(int p, int n, double *x, int *iw, double *w);
    inline double isolatedobj(int p, int n1, int n2, double *sf, double *si, double *sl);
    inline double stepobj(int p, int n1, int n2, double *sf, double *si, double *sl);    
    inline double maxstand(int ncp, double *stat, double *a, double *b);
    inline void scan(int p, int n, int first, int last, bool isolated, bool step, int lmin,
		     double *s, int &best, double &gain, int *cp);
    // length(iwork)=2+4*ncp+nm
    // length(work)=p*(m+1)+std::max(m,p)+2*nm
    inline void ggforward(int p, int n, int m, double *x, bool isolated, bool step, int lmin, int ncp,
			  double *l, double *r, double *sc, int *cp, double *stat,
			  int *iwork, double *work);
}

//-----------------------------------
// R Interface
//-----------------------------------

// [[Rcpp::export]]
List MPHASE1(NumericVector xx, bool isolated, bool step, int ncp, int lmin, int nperm) {
    IntegerVector dim=xx.attr("dim");
    int i, j, p=dim[0], n=dim[1], m=dim[2], nm=n*m, pnm=p*nm;
    IntegerVector stepssteps(1+2*ncp), iwork(2+4*ncp+nm);
    NumericVector scsc(pnm), ll(p), statstat(ncp), aa(ncp), bb(ncp),
	work(ncp*nperm+2*pnm+p*(m+1)+std::max(m,p)+2*nm);
    NumericMatrix rr(p,p);
    double *x=xx.begin(), *sc=scsc.begin(), *l=ll.begin(), *r=rr.begin(),
	*stat=statstat.begin(), *a=aa.begin(), *b=bb.begin(),
	*pstat=work.begin(), *xperm=pstat+ncp*nperm, *w=xperm+pnm, 
	*psi, pv, wobs, adj=nperm/(nperm-1.0);
    int *steps=stepssteps.begin(), *iw=iwork.begin();
    std::copy(x,x+pnm,xperm);
    for (i=1, psi=pstat; i<=nperm; i++, psi+=ncp) {
	ggmperm(p,nm,xperm);
	ggforward(p,n,m,xperm,isolated,step,lmin,ncp,l,r,sc,steps,psi,iw,w);
	for (j=0; j<ncp; j++) {
	    a[j] += (psi[j]-a[j])/i;
	    b[j] += (psi[j]*psi[j]-b[j])/i;
	}
    }
    for (j=0; j<ncp; j++) b[j] =std::max(eps,sqrt(adj*(b[j]-a[j]*a[j])));
    ggforward(p,n,m,x,isolated,step,lmin,ncp,l,r,sc,steps,stat,iw,w);
    wobs = maxstand(ncp,stat,a,b);
    for (i=0, psi=pstat, pv=0; i<nperm; i++, psi+=ncp) {
	if (maxstand(ncp,psi,a,b)>wobs) pv++;
    }
    pv /= nperm;
    scsc.attr("dim") = dim;
    int fcp = steps[0];
    CharacterVector type(fcp);
    IntegerVector time(fcp);
    NumericVector ss(fcp), sa(fcp), sb(fcp);
    for (i=0, steps++; i<fcp; i++, steps+=2) {
	time[i] = steps[0];
	type[i] = (steps[1]==1) ? "Step" : "Isolated";
	ss[i] = stat[i];
	sa[i] = a[i];
	sb[i] = b[i];
    }
    return List::create(_["center"]=ll, _["r"]=rr, _["signed.ranks"]=scsc,
			_["forward"] =DataFrame::create(_["type"]=type,_["time"]=time,
							_["T"]=ss,_["a"]=sa,_["b"]=sb),
			_["Wobs"]=wobs,_["p.value"]=pv);
}


namespace {

    // --------------------------------------------------
    // IMPLEMENTATION
    // --------------------------------------------------
 
    inline void ggmperm(int p, int n, double *x) {
	int i;
	while (n) {
	    i = floor(n*unif_rand()); n--;
	    F77_CALL(dswap)(&p,x+i*p,&ione,x+n*p,&ione);
	}
    }

    struct Comparator
    {
	Comparator(const double *x) : data(x) {}
	bool operator()(int left, int right) const { return data[left] < data[right]; }
	const double *data;
    };

    // length(a)=n; set r=x to replace the values of vector with their ranks
    inline void ggrank(int n, double *x, double *r, int *a) {
	int i, j, s;
	double d;
	for (i=0; i<n; i++) a[i] = i;
	std::sort(a,a+n,Comparator(x));
	for (i=0; i<n;) {
	    for (s=i, j=i+1; (j<n)&&(x[a[j]]<=x[a[i]]); j++) s += j;
	    for (d=1+((double)s)/(j-i); i<j; i++) r[a[i]] = d;
	}   
    }


    // add one row to the Cholesky factor X'X
    inline void gggivens(int p, double *r, double *x) {
	int i, p1=p+1;
	double t1, t2, g1, g2;
	for (i=p; i; i--, r+=p1, x++) {
	    t1 = r[0]; t2 = x[0];
	    if (fabs(t2)>0.0) {
		F77_CALL(drotg)(&t1,&t2,&g1,&g2);
		F77_CALL(drot)(&i, r, &p, x, &ione, &g1, &g2);
	    }
	}
    }



    // update the mean and the cholesky factor of cov; length(w)=p
    inline void ggupdbr(int p, int n, double *xb, double *r, double *x, double *w) {
	int i;
	double d, n1=n+1.0, sn=sqrt(n/n1);
	for (i=0; i<p; i++) {
	    d = x[i]-xb[i]; xb[i] += d/n1; w[i] = sn*d;
	}
	gggivens(p, r, w);
    } 

    inline void ggsubgroupsums(int p, int n, int m, double *x, double *sx) {
	int i, j;
	for (i=0; i<m; i++, sx+=p) {
	    std::copy(x,x+p,sx);
	    for (j=1, x+=p; j<n; j++, x+=p) F77_CALL(daxpy)(&p,&one,x,&ione,sx,&ione);  
	}
    }

    // x[,i] <- x[,i]-m
    inline void ggcenter(int p, int n, double *x, double *m) {
	for (int i=0; i<n; i++, x+=p) F77_CALL(daxpy)(&p,&onem,m,&ione,x,&ione);
    }

    // it's safe to set y=x to transform x in place
    inline void ggorth(int p, int n, double *r, double *x, double *y) {
	int i, j;
	double s, rii, *ri;
	for (; n ; n--, x+=p, y+=p) {
	    for (i=0, ri=r ; i<p ; i++, ri+=p) {
		rii = ri[i];
		if (fabs(rii)<eps) {
		    y[i] = 0.0;
		} else {
		    for (j=0, s=0.0; j<i; j++) s += ri[j]*y[j];
		    y[i] = (x[i]-s)/rii;
		}
	    } 
	}
    }

    inline void ggnorm2(int p, int n, double *x, double *dist) {
	int i, j;
	double s;
	for (i=0; i<n; i++, x+=p) {
	    for (j=0, s=0.0; j<p; j++) s += x[j]*x[j];
	    dist[i] = s;
	}
    }

    // length(t)=n
    void ggmmedinit(int p, int n, double *x, double *m, double *t) {
	int i, j, k, half=n/2;
	for (j=0; j<p ; j++) {
	    F77_CALL(dcopy)(&n, x+j, &p, t, &ione);
	    std::nth_element(t, t+half, t+n);
	    m[j] = t[half];
	}
	ggcenter(p,n,x,m);
    }

    // length(t)=p
    bool ggmmedstep(int p, int n, double *x, double *m, double *t) {
	bool notdone = false;
	int i, j;
	double w, s, dj, *xi;
	std::fill(t,t+p,0.0);
	for (i=0, s=0.0, xi=x ; i<n; i++, xi+=p) {
	    for (j=0, w=0.0; j<p; j++) {
		dj = xi[j]; w += dj*dj;
	    }
	    w = 1/std::max(sqrt(w),eps);
	    s += w; w /= s;
	    for (j=0; j<p; j++) t[j] += w*(xi[j]-t[j]);
	}
	ggcenter(p,n,x,t);
	for (j=0; j<p; j++) {
	    if (fabs(t[j])>eps*std::max(1.0,fabs(m[j]))) notdone = true;
	    m[j] += t[j];
	}
	return notdone;
    }

    // length(t)=max(p,n)
    int ggmmed(int p, int n, double *x, double *m, double *t) {
	int i=0;
	ggmmedinit(p,n,x,m,t);
	while(ggmmedstep(p,n,x,m,t) && (i<MMedIter)) i++;
	return i;
    }

    // length(w)=pm+max(p,m)
    void ggmc(int p, int n, int m, double *x, double *l, double *r, double *w) {
	int i, j, p2=p*p, pm=p*m;
	double *xb=w, *xbi=w, sn=1.0/sqrt(static_cast<double>(m*(n-1)));
	std::fill(r,r+p2,0.0); std::fill(xb,xb+pm,0.0); w=w+pm;
	for (i=0; i<m; i++, xbi+=p) for (j=0; j<n; j++,x+=p) ggupdbr(p,j,xbi,r,x,w);
	F77_CALL(dscal)(&p2,&sn,r,&ione);
	ggorth(p,m,r,xb,xb);
	ggmmed(p,m,xb,l,w);
    }
		 
    // length(w)=max(p,m)
    void ggmc1(int p, int m, double *x, double *l, double *r, double *z, double *w) {
	int i, j, p2=p*p;
	double *xi, *xi1, d=sqrt(0.5/(m-1));
	std::fill(r,r+p2,0.0);
	for (i=1, xi1=x, xi=x+p; i<m; i++, xi1=xi, xi+=p) {
	    for (j=0; j<p; j++) w[j] = xi[j]-xi1[j];
	    gggivens(p, r, w);
	}
	F77_CALL(dscal)(&p2, &d, r, &ione);
	ggorth(p,m,r,x,z);
	ggmmed(p,m,z,l,w);
    }

    // length(iw)=n length(w)=2*n 
    void ggsignedrank(int p, int n, double *x, int *iw, double *w) {
	int i, n1=n+1;
	double u, *xi, *d=w+n; 
	ggnorm2(p,n,x,d);
	ggrank(n,d,w,iw);
	for (i=0, xi=x; i<n; i++, xi+=p) {
	    if (d[i]>eps2) {
		u = sqrt(R::qchisq(w[i]/n1,p,1,0)/d[i]);
		F77_CALL(dscal)(&p,&u,xi,&ione);
	    } else {
		std::fill(xi,xi+p,0.0);
	    }
	}
    }


    inline double stepobj(int p, int n1, int n2, double *sf, double *si, double *sl) {
	int r;
	double g, d;
	for (r=0, g=0; r<p; r++) {
	    d = ((si[r]-sf[r])/n1)-(sl[r]-si[r])/n2;
	    g += d*d;
	}
	return n1*n2*g/(n1+n2);
    }

    inline double isolatedobj(int p, int n1, int n2, double *sf, double *si, double *sl) {
	int r;
	double g, d, di, *sii=si-p;
	for (r=0, g=0; r<p; r++) {
	    di = si[r]-sii[r];
	    d = (di/n1)-(sl[r]-sf[r]-di)/n2;
	    g += d*d;
	}
	return n1*n2*g/(n1+n2);
    }

    inline void scan(int p, int n, int first, int last, bool isolated, bool step, int lmin,
		     double *s, int &best, double &gain, int &type, int *cp) {
	int i, n1, n2, nm, nm1, nmin=n*lmin;
	double g, *sf=s+first*p, *sl=s+last*p, *si=sf+p;
	for (i=first+1, nm=0; i<=last; i++) if (cp[i]) nm++;
	gain = 0.0; nm = n*(last-first-nm); nm1 = nm-n;
	for (i=first+1, n1=0; i<=last; i++, si+=p) {
	    if (!cp[i]) {
		n1 += n; n2 = nm-n1;
		if (step && (n1>nmin) && (n2>nmin)) {
		    g = stepobj(p, n1, n2, sf, si, sl);
		    if (R_FINITE(g) && (g>(gain))) {
			best=i; gain=g; type=1;
		    }
		}
		if (isolated) {
		    g = isolatedobj(p, n, nm1, sf, si, sl);
		    if (R_FINITE(g) && (g>(gain))) {
			best=i; gain=g; type=2;
		    }
		}
	    }
	}
    }

    inline double maxstand(int ncp, double *stat, double *a, double *b) {
	double v = (stat[0]-a[0])/b[0];
	for (int i=1; i<ncp; i++) v = std::max(v,(stat[i]-a[i])/b[i]);
	return v;
    }
    
    inline void ggforward(int p, int n, int m, double *x, bool isolated, bool step, int lmin, int ncp,
			  double *l, double *r, double *sc, int *cp, double *stat,
			  int *iwork, double *work) {
	int i, j, nm=n*m, pnm=p*nm, nc1=ncp+1, iopt, ntau, *bcp,
	    *tau = iwork, *split = tau+nc1, *type=split+nc1, *icp=type+nc1;
	double *sb, *s, *si, *sii, *gi, st;
	if (n==1) {
	} else {
	}
	// computation of the signed ranks
	if (n==1) {
	    sb = sc; s = work; 
	    ggmc1(p, m, x, l, r, sc, work);
	    ggsignedrank(p, m, sc, iwork, work);
	} else {
	    sb = work; s = sb+p*m;
	    ggmc(p, n, m, x, l, r, work);
	    ggorth(p, nm, r, x, sc);
	    ggcenter(p, nm, sc, l);
	    ggsignedrank(p,nm,sc,iwork,work);
	    ggsubgroupsums(p,n,m,sc,sb);	    
	}
	// accumulation of the subgroup sums
	si = s+p; std::fill(s,si,0.0); std::copy(sb,sb+p*m,si); 
	for (i=2, sii=si+p; i<=m; i++, si=sii, sii+=p) {
	    F77_CALL(daxpy)(&p, &one, si, &ione, sii, &ione);
	}
	// main loop
	std::fill(icp,icp+m+1,0);
	*cp=0;  bcp=cp+1; ntau = 2; tau[0] = 0 ; tau[1] = -m; gi = s+(m+1)*p;
	while (*cp < ncp) {
	    for (i=iopt=1 ; i<ntau ; i++) {
		if (tau[i]<0) {
		    tau[i] = -tau[i];
		    scan(p, n, tau[i-1], tau[i], isolated, step, lmin, s, split[i], gi[i], type[i], icp);
		}
		if (gi[i]>gi[iopt]) iopt=i;
	    }
	    if (*cp==0) stat[0]=gi[iopt]; else stat[*cp]=stat[*cp-1]+gi[iopt];
	    if (gi[iopt]<eps) {
		for (i=*cp, st=(*cp==0)?0.0:stat[*cp]; i<ncp; i++) stat[i] = st;
		return;
	    }
	    bcp[0]=split[iopt]; bcp[1]=type[iopt]; *cp+=1; bcp+=2;
	    icp[split[iopt]] = 1;
	    if (type[iopt]==1) {
		memmove(tau+iopt+1,tau+iopt,(ntau-iopt)*sizeof(int));
		memmove(split+iopt+1,split+iopt,(ntau-iopt)*sizeof(int));
		memmove(type+iopt+1,type+iopt,(ntau-iopt)*sizeof(int));
		memmove(gi+iopt+1,gi+iopt,(ntau-iopt)*sizeof(double));
		ntau++; tau[iopt] = -split[iopt+1]; tau[iopt+1] = -tau[iopt+1];
	    } else {
		ggcenter(p,m+1-split[iopt],s+split[iopt]*p,sb+split[iopt]*p-p);
		tau[iopt] = -tau[iopt];
	    }
	}
	for (i=0, bcp=cp+1; i<*cp; i++, bcp+=2) if (bcp[1]==1) bcp[0]+=1;
    }
    
}


