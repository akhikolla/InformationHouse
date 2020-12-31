#include <Rcpp.h>
using namespace Rcpp;


/// dpseg in Rcpp, see R/dpseg.R for documentation

// initialize a matrix to NA
// [[Rcpp::export]]
NumericMatrix na_matrix(int n){
  NumericMatrix m(n,n) ;
  std::fill( m.begin(), m.end(), NumericVector::get_na() ) ;
  return m ;
}


// [[Rcpp::export]]
NumericVector backtrace_c(NumericVector imax, int jumps=0) {

  int end = imax.size();
  NumericVector ends (end);
  ends[0] = end;
  int cnt = 1; 
  while ( end>1 ) {
    end = imax[end-1] - jumps;
    ends[cnt] = end;
    cnt += 1;
  }
  // cut and reverse
  NumericVector res (cnt);
  for ( int i=0; i<cnt; i++ ) 
    res[cnt-i-1] = ends[i];
  return res;

}


/// S_j = max_{j-maxl < i <= j-minl} S_[i-jumps] + fitscore(i,j) - P
/// NOTE: argument scoref is a dummy to be consistent with generic recursion
// [[Rcpp::export]]
List recursion_linreg_c(NumericVector x, NumericVector y, int maxl,
			bool jumps=0, double P=0, int minl=3, double S0=1.0,
			String type="var", bool storem=0, bool storev=1,
			String scoref="dummy") {

  int N = x.length();
  int maxn = maxl-minl+1; // length of considered S[i] vector
  
  // main results: S_j and positions i that yielded max
  NumericVector S(N);
  NumericVector imax(N);
  // vector "S_i + fitscore(i+1,j) - P" for max_{j-maxl < i <= j-minl}
  // si as Rcpp vector only to use Rcpp::which_max
  NumericVector si(maxn); // score function i+1 to j
  
  // incremental variance calculation vectors
  // NOTE: using NumericVector causes overflow problems for large numbers
  std::vector<long double> sx(N); 
  std::vector<long double> sy(N); 
  std::vector<long double> sx2(N); 
  std::vector<long double> sy2(N); 
  std::vector<long double> sxy(N);   

  // optionally stored values
  // TODO: save memory by declaring but only optionally allocating variables?
  //NumericVector ivar, islp, icpt, irsq, icor, vari, slpi, cpti, rsqi, cori;
  //if ( storev ) {
  // returned values
  NumericVector ivar(N); // variances i+1 to j for imax
  NumericVector islp(N); // slopes i+1 to j for imax
  NumericVector icpt(N); // intercepts i+1 to j for imax
  NumericVector irsq(N); // r-squared i+1 to j for imax
  NumericVector icor(N); // correlation i+1 to j for imax
  
  // internal 
  std::vector<double> vari(maxn); // variance i+1 to j
  std::vector<double> slpi(maxn); // slope i+1 to j
  std::vector<double> cpti(maxn); // intercept i+1 to j
  std::vector<double> rsqi(maxn); // r2 i+1 to j
  std::vector<double> cori(maxn); // correlation i+1 to j
  
  if ( storev ) {
    std::fill( ivar.begin(), ivar.end(), 0.0 );
    std::fill( islp.begin(), islp.end(), 0.0 );
    std::fill( icpt.begin(), icpt.end(), 0.0 );
    std::fill( irsq.begin(), irsq.end(), 0.0 );
    std::fill( icor.begin(), icor.end(), 0.0 );
  }
  
  // score function matrix: only required for educational movie
  NumericMatrix SCR(0,0);
  if ( storem )
    SCR = na_matrix(N); 
  
  // INITIALIZATION
  // initialize S and imax
  //S0 = 1;//-P;
  std::fill( S.begin(), S.end(), -P );
  std::fill( imax.begin(), imax.end(), 0.0 );

  // initialize x, y, x2, y2, xy
  std::fill( sx.begin(), sx.end(),0.0 );
  std::fill( sy.begin(), sy.end(),0.0 );
  std::fill( sx2.begin(), sx2.end(),0.0 );
  std::fill( sy2.begin(), sy2.end(),0.0 );
  std::fill( sxy.begin(), sxy.end(),0.0 );
  sx[0] = x[0];
  sy[0] = y[0];
  sx2[0] = x[0]*x[0];
  sy2[0] = y[0]*y[0];
  sxy[0] = y[0]*x[0];

  // variances
  long double Syy, Sxx, Sxy;
  // RECURSION
  for ( int j=1; j<N; j++ ) {

    // incremental calculation of variances
    sx[j]  =  sx[j-1] + x[j];        // \sum x
    sy[j]  =  sy[j-1] + y[j];        // \sum y
    sx2[j] = sx2[j-1] + x[j]*x[j];   // \sum x^2
    sy2[j] = sy2[j-1] + y[j]*y[j];   // \sum y^2
    sxy[j] = sxy[j-1] + x[j]*y[j];   // \sum x*y


    // variance from 1 to j - to get variance for first segment
    int n = j+1; // length n
    // Syy = \sum y^2 - (\sum y)^2/n
    Syy = sy2[j] - sy[j]*sy[j]/n;
    // Sxx = \sum x^2 - (\sum x)^2/n
    Sxx = sx2[j] - sx[j]*sx[j]/n;
    // Sxy = \sum x*y - \sum(x)*\sum(y)/n
    Sxy = sxy[j] - sx[j]*sy[j]/n;

    
    // calculate all scores i->j
    std::fill( si.begin(), si.end(),
	       - std::numeric_limits<double>::infinity() );
    
    int idx = std::max(maxl-(j+1),0); // start at i=1 at starts
    for ( ; idx < maxn; idx++ ) {

      // current index i = (j+1)-maxl+1+idx - 1-based index
      int i = (j+1)-maxl+idx; // 0-based array index i
      int ij = i-1 ; // sum of squares subtraction index
      n = j-ij; // length of segment i -> j

      long double scr, var, rsq, slp, cor;

      // sum values from i to j : subtract values from 1 to i-1
      if ( ij>=0 ) {
	// Syy = \sum y^2 - (\sum y)^2/n
	Syy = (sy2[j]-sy2[ij]) - (sy[j]-sy[ij])*(sy[j]-sy[ij])/n;
	// Sxx = \sum x^2 - (\sum x)^2/n
	Sxx = (sx2[j]-sx2[ij]) - (sx[j]-sx[ij])*(sx[j]-sx[ij])/n;
	// Sxy = \sum x*y - \sum(x)*\sum(y)/n
	Sxy = (sxy[j]-sxy[ij]) - (sx[j]-sx[ij])*(sy[j]-sy[ij])/n;
      }
      
      // slope 
      slp = Sxy/Sxx;

      // goodness of fit measures for scoring function
      var = (Syy - Sxy*slp) / (n-1);
      rsq = slp*slp*Sxx/Syy;
      if ( rsq>=0.0 )
	cor =  std::sqrt(rsq);
      else cor = 0.0;

      // score
      scr =  - std::numeric_limits<double>::infinity(); // compiler warning
      if ( type=="var" ) 
	scr = -var;
      else if ( type=="cor" )
	scr = cor -1;
      else if ( type=="r2" )
	scr = rsq -1;

      // CALCULATE SCORE
      double si1;
      if ( i-jumps == -1 )
	si1 =  S0;
      else
	si1 = S[i-jumps]; 
      if ( Sxx!=0.0 ) // avoid! see https://de.wikipedia.org/wiki/Verschiebungssatz_(Statistik)
	si[idx] = si1 + scr - P;
      else warning("numeric cancellation: Sxx==0.0");
      
      // store values
      if ( storev ) {

	// slope & intercept i+1,j
	slpi[idx] = slp;
	if ( ij >= 0 )
	  cpti[idx] = ((sy[j]-sy[ij]) - slp*(sx[j]-sx[ij]))/n;
	else
	  cpti[idx] = ((sy[j]) - slp*(sx[j]))/n;
	  
	// goodness of fit measures for scoring function
	vari[idx] = var;
	rsqi[idx] = rsq;
	cori[idx] = cor;
      }
      if ( storem )
	SCR(i,j) = scr;
    }
    
    // GET MAXIMUM SCORE
    // S[j] = max S[i-shift] + score(i,j)
    S[j] = max(si);
    idx = which_max(si);

    // ... store which i yielded maximum - 1-based index
    // i = (j+1)-maxl+idx // 0-based
    imax[j] = j-maxl+idx+1+1; // 1-based
    // ... store parameters i+1 to j
    if ( storev ) {
      ivar[j] = vari[idx];
      islp[j] = slpi[idx];
      icpt[j] = cpti[idx];
      irsq[j] = rsqi[idx];
      icor[j] = cori[idx];
    }
  }

  // store recorded values
  List values;
  if ( storev ) 
    values = List::create(Named("ivar") =ivar, /* variance i+1 to j */
			  Named("islp") =islp, /* slope i+1 to j */
			  Named("icpt") =icpt, /* intercept i+1 to j*/
			  Named("irsq") =irsq, /* r2 i+1 to j */
			  Named("icor") =icor);/* cor i+1 to j */
  
  return List::create(Named("S") = S,    /* scoring vector S(j) */
		      Named("imax") = imax, /* back-tracing  */
		      Named("values") = values,
		      Named("SCR") = SCR);  
}

/// S_j = max_{j-maxl < i <= j-minl} S_[i-jumps] + scorematrix(i,j) - P
/// NOTE: arguments scoref is a dummy to be consistent with generic recursion
// [[Rcpp::export]]
List recursion_matrix(NumericVector x, NumericMatrix y, int maxl,
		      bool jumps=0, double P=0, int minl=3, double S0=1.0,
		      String type="matrix", bool storem=0, bool storev=1,
		      String scoref="dummy") {

  int N = y.nrow();
  int maxn = maxl-minl+1; // length of considered S[i] vector
  
  // main results: S_j and positions i that yielded max
  NumericVector S(N);
  NumericVector imax(N);
  // vector "S_i + fitscore(i+1,j) - P" for max_{j-maxl < i <= j-minl}
  // si as Rcpp vector only to use Rcpp::which_max
  NumericVector si(maxn); // score function i+1 to j
  
  
  // INITIALIZATION
  // initialize S and imax
  //S0 = 1;//-P;
  std::fill( S.begin(), S.end(), -P );
  std::fill( imax.begin(), imax.end(), 0.0 );


  // RECURSION
  for ( int j=1; j<N; j++ ) {

    // calculate all scores i->j
    std::fill( si.begin(), si.end(),
	       - std::numeric_limits<double>::infinity() );
    
    int idx = std::max(maxl-(j+1),0); // start at i=1 at starts
    for ( ; idx < maxn; idx++ ) {

      // current index i = (j+1)-maxl+1+idx - 1-based index
      int i = (j+1)-maxl+idx; // 0-based array index i
  

      // CALCULATE SCORE
      double si1;
      if ( i-jumps == -1 )
	si1 =  S0;
      else
	si1 = S[i-jumps]; 
      si[idx] = si1 + y(i,j) - P;
    }
    
    // GET MAXIMUM SCORE
    // S[j] = max S[i-shift] + score(i,j)
    S[j] = max(si);
    idx = which_max(si);

    // ... store which i yielded maximum - 1-based index
    // i = (j+1)-maxl+idx // 0-based
    imax[j] = j-maxl+idx+1+1; // 1-based
  }

  List values;
    
  return List::create(Named("S") = S,    /* scoring vector S(j) */
		      Named("imax") = imax); /* back-tracing  */
}

/// S_j = max_{j-maxl < i <= j-minl} S_i + fitscore(i+1,j) - P
/// NOTE: outdated, kept for demonstration purposes; S0
/// jumps and S0 have no effect, just kept for consistency with wrapper
// [[Rcpp::export]]
List scoref_c(NumericVector x, NumericVector y, 
	      int minl, int maxl, double P, double S0,
	      const String type="var", bool jumps=1,
	      bool storem=0, bool storev=1) {
  
  int N = x.length();
  int maxn = maxl+1; // length of considered S[i] vector
  
  // main results: S_j and positions i that yielded max
  NumericVector S(N);
  NumericVector imax(N);
  // vector "S_i + fitscore(i+1,j) - P" for max_{j-maxl < i <= j-minl}
  // si as Rcpp vector only to use Rcpp::which_max
  NumericVector si(maxn); // score function i+1 to j
  
  // incremental variance calculation vectors
  // NOTE: using NumericVector causes overflow problems for large numbers
  std::vector<long double> sx(N); 
  std::vector<long double> sy(N); 
  std::vector<long double> sx2(N); 
  std::vector<long double> sy2(N); 
  std::vector<long double> sxy(N);   

  // optionally stored values

  // returned values
  NumericVector ivar(N); // variances i+1 to j for imax
  NumericVector jvar(N); // variance 1 to j (i=0)
  NumericVector islp(N); // slopes i+1 to j for imax
  NumericVector jslp(N); // slope 1 to j (i=0)
  NumericVector icpt(N); // intercepts i+1 to j for imax
  NumericVector jcpt(N); // intercepts 1 to j (i=0)
  NumericVector irsq(N); // r-squared i+1 to j for imax
  NumericVector jrsq(N); // r-squared 1 to j (i=0)
  NumericVector icor(N); // correlation i+1 to j for imax
  NumericVector jcor(N); // correlation 1 to j (i=0)
  
  // internal 
  std::vector<double> vari(maxn); // variance i+1 to j
  std::vector<double> slpi(maxn); // slope i+1 to j
  std::vector<double> cpti(maxn); // intercept i+1 to j
  std::vector<double> rsqi(maxn); // r2 i+1 to j
  std::vector<double> cori(maxn); // correlation i+1 to j
  
  std::fill( ivar.begin(), ivar.end(), 0.0 );
  std::fill( jvar.begin(), jvar.end(), 0.0 );
  std::fill( islp.begin(), islp.end(), 0.0 );
  std::fill( jslp.begin(), jslp.end(), 0.0 );
  std::fill( icpt.begin(), icpt.end(), 0.0 );
  std::fill( jcpt.begin(), jcpt.end(), 0.0 );
  std::fill( irsq.begin(), irsq.end(), 0.0 );
  std::fill( jrsq.begin(), jrsq.end(), 0.0 );
  std::fill( icor.begin(), icor.end(), 0.0 );
  std::fill( jcor.begin(), jcor.end(), 0.0 );
  
  jvar[0] = 0.0;
  jslp[0] = 0.0;
  jcpt[0] = 0.0;
  jrsq[0] = 0.0;
  jcor[0] = 0.0;
    
  
  // score function matrix: only required for educational movie
  NumericMatrix SCR(0,0);
  if ( storem )
    SCR = na_matrix(N); 
  
  // INITIALIZATION
  // initialize S and imax
  std::fill( S.begin(), S.end(), -P );
  std::fill( imax.begin(), imax.end(), 0.0 );

  // initialize x, y, x2, y2, xy
  std::fill( sx.begin(), sx.end(),0.0 );
  std::fill( sy.begin(), sy.end(),0.0 );
  std::fill( sx2.begin(), sx2.end(),0.0 );
  std::fill( sy2.begin(), sy2.end(),0.0 );
  std::fill( sxy.begin(), sxy.end(),0.0 );
  sx[0] = x[0];
  sy[0] = y[0];
  sx2[0] = x[0]*x[0];
  sy2[0] = y[0]*y[0];
  sxy[0] = y[0]*x[0];

  // variances
  long double Syy, Sxx, Sxy;
  // RECURSION
  for ( int j=1; j<N; j++ ) {

    // incremental calculation of variances
    sx[j]  =  sx[j-1] + x[j];        // \sum x
    sy[j]  =  sy[j-1] + y[j];        // \sum y
    sx2[j] = sx2[j-1] + x[j]*x[j];   // \sum x^2
    sy2[j] = sy2[j-1] + y[j]*y[j];   // \sum y^2
    sxy[j] = sxy[j-1] + x[j]*y[j];   // \sum x*y


    // variance from 1 to j - to get variance for first segment
    int n = j+1; // length n
    // Syy = \sum y^2 - (\sum y)^2/n
    Syy = sy2[j] - sy[j]*sy[j]/n;
    // Sxx = \sum x^2 - (\sum x)^2/n
    Sxx = sx2[j] - sx[j]*sx[j]/n;
    // Sxy = \sum x*y - \sum(x)*\sum(y)/n
    Sxy = sxy[j] - sx[j]*sy[j]/n;

    if ( storev ) {
      // slope 1 to j
      jslp[j] = Sxy/Sxx;
      // intercept 1 to j - TODO: why off?? correlates with slope!?
      jcpt[j] = sy[j]/n - jslp[j]*(sx[j])/n;
      // variance 1 to j
      jvar[j] = (Syy - Sxy*jslp[j]) / (n-1);
      // R2= b1^2 * Sxx/Syy
      jrsq[j] = jslp[j]*jslp[j]*Sxx/Syy;
      // correlation
      jcor[j] = std::pow(jrsq[j], 0.5);
    }
    
    // calculate all scores i->j
    std::fill( si.begin(), si.end(),
	       - std::numeric_limits<double>::infinity() );
    
    int idx = std::max(maxl-j,0);
    for ( ; idx < maxl; idx++ ) {
      
      int i = j-maxl+idx; // current index i
      n = j-i; // length of segment i+1 -> j
      double scr, var, rsq, slp, cor;
      
      //if ( n < minl ) continue;	

      // sum values from i+1 to j : subtract values from 1 to i
      // Syy = \sum y^2 - (\sum y)^2/n
      Syy = (sy2[j]-sy2[i]) - (sy[j]-sy[i])*(sy[j]-sy[i])/n;
      // Sxx = \sum x^2 - (\sum x)^2/n
      Sxx = (sx2[j]-sx2[i]) - (sx[j]-sx[i])*(sx[j]-sx[i])/n;
      // Sxy = \sum x*y - \sum(x)*\sum(y)/n
      Sxy = (sxy[j]-sxy[i]) - (sx[j]-sx[i])*(sy[j]-sy[i])/n;
      
      // slope 
      slp = Sxy/Sxx;

      // goodness of fit measures for scoring function
      var = (Syy - Sxy*slp) / (n-1);
      rsq = slp*slp*Sxx/Syy;
      if ( rsq>=0.0 )
	cor =  std::sqrt(rsq); //NOTE: abs(cor)
      else cor = 0;

      // score
      scr = - std::numeric_limits<double>::infinity(); // compiler warning
      if ( type=="var" ) 
	scr = -var;
      else if ( type=="cor" )
	scr = cor;
      else if ( type=="r2" )
	scr = rsq;
      
      // CALCULATE SCORE
      if ( j-i >= minl ) 
	if ( Sxx!=0.0 )  // avoid! see https://de.wikipedia.org/wiki/Verschiebungssatz_(Statistik)
	  si[idx] = S[i] + scr - P;

      // store values
      if ( storev ) {

	// slope & intercept i+1,j
	slpi[idx] = slp;
	cpti[idx] = ((sy[j]-sy[i]) - slp*(sx[j]-sx[i]))/n;

	// goodness of fit measures for scoring function
	vari[idx] = var;
	rsqi[idx] = rsq;
	cori[idx] =  cor;
      }
      if ( storem )
	SCR(i,j) = scr;
    }
    
    // GET MAXIMUM SCORE
    S[j] = max(si);
    idx = which_max(si);

    // ... store which i yielded maximum
    imax[j] = j-maxl + idx+1;
    // ... store parameters i+1 to j
    if ( storev ) {
      ivar[j] = vari[idx];
      islp[j] = slpi[idx];
      icpt[j] = cpti[idx];
      irsq[j] = rsqi[idx];
      icor[j] = cori[idx];
    }
  }

  // recorded values
  List values;
  if ( storev ) 
    values = List::create(Named("ivar") =ivar, /* variance i+1 to j */
			  Named("jvar") =jvar, /* variance 1 to j */
			  Named("islp") =islp, /* slope i+1 to j */
			  Named("jslp") =jslp, /* slope 1 to j */
			  Named("icpt") =icpt, /* intercept i+1 to j*/
			  Named("jcpt") =jcpt, /* intercept 1 to j */
			  Named("irsq") =irsq, /* r2 i+1 to j */
			  Named("jrsq") =jrsq, /* r2 1 to j */
			  Named("icor") =icor, /* cor i+1 to j */
			  Named("jcor") =jcor);/* cor 1 to j */
  
  return List::create(Named("S") = S,    /* scoring vector S(j) */
		      Named("imax") = imax, /* back-tracing  */
		      Named("values") = values,
		      Named("SCR") = SCR);  
}

/// S_j = max_{j-maxl < i <= j-minl} S_i + fitscore(i+1,j) - P
/// NOTE: this implementation has segments <minl at start (j<minl)
/// and is equivalent to the current R implementation;
/// jumps and S0 have no effect, just kept for consistency with wrapper
// [[Rcpp::export]]
List scoref_oldc(NumericVector x, NumericVector y, 
		 int minl, int maxl, double P, double S0,
		 const String type="var", bool jumps=0,
		 bool storem=0, bool storev=1) {

  int N = x.length();

  // main results: S_j and positions i that yielded max
  NumericVector S(N);
  NumericVector imax(N);
  NumericVector ieps(N); // variances i+1 to j for imax

  // score function matrix: only required for educational movie
  NumericMatrix SCR(0,0);
  if ( storem )
    SCR = na_matrix(N); 

  // vector "S_i + fitscore(i+1,j) - P" for max_{j-maxl < i <= j-minl}
  // si as Rcpp vector only to use Rcpp::which_max
  NumericVector si(maxl-minl+1); // score function i to j
  std::vector<double> epsi(maxl-minl+1); // variance i to j
  
  // incremental variance calculation vectors
  // NOTE: using NumericVector causes overflow problems for large numbers
  std::vector<long double> sx(N); 
  std::vector<long double> sy(N); 
  std::vector<long double> sx2(N); 
  std::vector<long double> sy2(N); 
  std::vector<long double> sxy(N); 

  // initialize S and imax
  std::fill( S.begin(), S.end(), -P );
  std::fill( imax.begin(), imax.end(), 0 );
  std::fill( ieps.begin(), ieps.end(), 0 );

  // initialize x, y, x2, y2, xy
  std::fill( sx.begin(), sx.end(),0 );
  std::fill( sy.begin(), sy.end(),0 );
  std::fill( sx2.begin(), sx2.end(),0 );
  std::fill( sy2.begin(), sy2.end(),0 );
  std::fill( sxy.begin(), sxy.end(),0 );
  sx[0] = x[0];
  sy[0] = y[0];
  sx2[0] = x[0]*x[0];
  sy2[0] = y[0]*y[0];
  sxy[0] = y[0]*x[0];

  for ( int j=1; j<N; j++ ) {

    // fill up variances
    sx[j]  =  sx[j-1] + x[j];        // \sum x
    sy[j]  =  sy[j-1] + y[j];        // \sum y
    sx2[j] = sx2[j-1] + x[j]*x[j];   // \sum x^2
    sy2[j] = sy2[j-1] + y[j]*y[j];   // \sum y^2
    sxy[j] = sxy[j-1] + x[j]*y[j];   // \sum x*y

    // calculate all scores i->j
    if ( j-minl >=  std::max(0, j-maxl) ) {
  
      std::fill( si.begin(), si.end(),
		 - std::numeric_limits<double>::infinity() );

      int idx = std::max(maxl-j,0);
      for ( ; idx <= (maxl-minl); idx++ ) {

	int i = j-maxl+idx; // current index i
	int n = j-i; // length of segment i -> j
	
	// variances
	long double Syy, Sxx, Sxy;

	// Syy = \sum y^2 - (\sum y)^2/n
	Syy = (sy2[j]-sy2[i]) - (sy[j]-sy[i])*(sy[j]-sy[i])/n;
	// Sxx = \sum x^2 - (\sum x)^2/n
	Sxx = (sx2[j]-sx2[i]) - (sx[j]-sx[i])*(sx[j]-sx[i])/n;
	// Sxy = \sum x*y - sum(x)*sum(y)/n
	Sxy = (sxy[j]-sxy[i]) - (sx[j]-sx[i])*(sy[j]-sy[i])/n;

	// SCORING FUNCTION: variance, from i to j:
	epsi[idx] = (Syy - (Sxy*Sxy)/Sxx) / (n-1);

	// store full matrix
	if ( storem )
	  SCR(i,j) = -epsi[idx];

	// CALCULATE SCORE
	if ( Sxx!=0.0 ) // TODO: why can Sxx be 0?
	  si[idx] = S[i] - epsi[idx] - P;
	
      }
      
      // get maximum 
      S[j] = max(si);
      idx = which_max(si);
      // ... and store which i yielded maximum
      imax[j] = j-maxl + idx+1;
      // ... and variances i+1 to j
      ieps[j] = epsi[idx];
      
    }
  }

  List values = List::create(Named("ivar") = ieps); /* variance i+1 to j */
		      //Named("sx") = sx, Named("sy")=sy, 
		      //Named("sx2") = sx2, Named("sy2")=sy2, 
		      //Named("sxy") = sxy, 
  
  return List::create(Named("S") = S,    /* scoring vector S(j) */
		      Named("imax") = imax, /* back-tracing  */
		      Named("values") = values, /* variance */
		      Named("SCR") = SCR);  
}
