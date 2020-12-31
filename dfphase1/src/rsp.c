#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

static void ggperm(int n, double *y) {
    int i;
    double yi;
    GetRNGstate();
    while (n>1) {
	i = floor(n*unif_rand()); n -= 1;
	yi = y[i]; y[i] = y[n]; y[n] = yi;
    }
    PutRNGstate();
}

static double lobj(int nall,double mall,int n1,double m1,int n2,double m2){
    return(n1*n2*(m1-m2)*(m1-m2)/nall);
}

static double sobj(int nall,double mall,int n1,double m1,int n2,double m2){
    return(n1*log(mall/m1)+n2*log(mall/m2));
}

static void scan(int first, int last, int lmin, int n, double *sy,
		 double (*obj)(int,double,int,double,int,double),
		 int *best, double *gain) {
    int i, n1, n2, nm=n*(last-first);
    double mu, mu1, mu2, g;
    mu = (sy[last]-sy[first])/nm;
    for (i=first+lmin, *gain=0.0; i<=last-lmin; i++) {
	n1 = n*(i-first);
	n2 = nm-n1;
	mu1 = (sy[i]-sy[first])/n1;
	mu2 = (sy[last]-sy[i])/n2;
	g = obj(nm,mu,n1,mu1,n2,mu2);
	if (R_FINITE(g) && (g>(*gain))) {
	    *best = i;
	    *gain = g; 
	}
    }
} 


static void ggsx(int m, int n, double *x, double *sx) {
    int i, j;
    for (i=0 ; i<m ; i++,sx++,x++) {
	for (j=0, *sx=0.0 ; j<n ; j++) {
	    *sx += x[j*m]; 
	}
    }
}

static void ggsx2(int m, int n, double *x, double *sx) {
    int i, j;
    double xij;
    for (i=0 ; i<m ; i++,sx++,x++) {
	for (j=0, *sx=0.0 ; j<n ; j++) {
	    xij = x[j*m];
	    *sx += xij*xij; 
	}
    }
}

static void ggcumsum(int m, double *sx, double *cs) {
    int i;
    for (i=0, cs[0]=0.0; i<m ; i++) cs[i+1] = cs[i]+sx[i];
}

static double ggstat(int m, int n, double *sx) {
    int i ;
    double me, stat;
    for (i=0, me=0.0; i<m; i++) me += sx[i];
    me /= (m*n);
    for (i=0, stat=0.0 ; i<m ; i++) stat = fmax2(stat,fabs((sx[i]/n)-me));
    return(stat);
} 


/* 
   ipar=(m,n,nsteps,lmin,type,nperm)
   length(iwork)=2*(nsteps+1) length(work)=m+1+max(m,nsteps+2)
*/

#define NINT (*steps)
#define LOW(i) steps[2*(i)-1]
#define UP(i) steps[2*(i)]
enum type_ {level_=1,scale_};
void ggfitstep1(int *ipar, double *y, int *steps, double *stat,
                int *iwork, double *work) {
    int m=ipar[0], n=ipar[1], nsteps=ipar[2], lmin=ipar[3], type=ipar[4],
	i, iopt, ntau, *tau = iwork, *split = tau+nsteps+1;
    double *sy=work, *gi=work+m+1;
    static double (*obj)(int,double,int,double,int,double);
    if (type==level_) {
	obj = lobj;
	ggsx(m,n,y,gi);
	ggcumsum(m,gi,sy);
	if (n>1) stat[nsteps] = ggstat(m,n,gi);
    } else {
	obj = sobj;
	ggsx2(m,n,y,gi);
	ggcumsum(m,gi,sy);
	if (n>1) stat[nsteps] = ggstat(m,n,gi);
    }
    NINT=1; LOW(NINT)=0 ; UP(NINT)=m;
    ntau = 2; tau[0] = 0 ; tau[1] = -m;
    while (NINT <= nsteps) {
	for (i=iopt=1 ; i<ntau ; i++) {
	    if (tau[i]<0) {
		tau[i] = -tau[i];
		scan(tau[i-1],tau[i],lmin,n,sy,obj,&split[i],&gi[i]);
	    }
	    if (gi[i]>gi[iopt]) iopt=i;
	}
	if (gi[iopt]<DOUBLE_EPS) {
	    for (i=NINT-1; i<nsteps; i++) stat[i] = stat[i-1];
	    return;
	} 
	if (NINT==1) stat[0]=gi[iopt]; else stat[NINT-1]=stat[NINT-2]+gi[iopt];
	NINT+=1; LOW(NINT)=tau[iopt-1]; UP(NINT)=split[iopt]; 
	memmove(tau+iopt+1,tau+iopt,(ntau-iopt)*sizeof(int));
	memmove(split+iopt+1,split+iopt,(ntau-iopt)*sizeof(int));
	memmove(gi+iopt+1,gi+iopt,(ntau-iopt)*sizeof(double));
	ntau++; tau[iopt] = -split[iopt+1]; tau[iopt+1] = -tau[iopt+1];
    }
}

void ggstepfactor(int *level, int *steps, int *factor) {
    int i, j;
    for (i=1 ; i<=imin2(*level,NINT) ; i++) {
	for (j=LOW(i) ; j<UP(i) ; j++) 
	    factor[j] = i;
    }
}



void ggdotrsp(int *ipar, double *y, int *steps, double *stat, double *perm) {
    int m=ipar[0], n=ipar[1], nstat=(n==1)?ipar[2]:(ipar[2]+1), nperm=ipar[5],
	nm=n*m, *iwork, *psteps, i;
    double *work;
    iwork = (int *) R_alloc(1+2*(nstat+1)+2*(nstat+1),sizeof(int));
    work = (double *) R_alloc(m+1+imax2(nstat+2,m),sizeof(double));
    psteps = iwork+2*(nstat+1);
    ggfitstep1(ipar, y, steps, stat, iwork, work);
    for (i=0; i<nperm ; i++, perm+=nstat) {
	ggperm(nm,y);
	ggfitstep1(ipar,y,psteps,perm,iwork,work);
    }
}

