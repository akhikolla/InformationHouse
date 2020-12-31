

#include "interp.h"


// [[Rcpp::export]]
List interpDeltri(NumericVector x, NumericVector y,
                     NumericVector zD,
                     List t, // data xD and yD contained here!
                     CharacterVector input = "points",
                     CharacterVector output = "grid") {

  List T(t);
  int nT = T.size();
  int nG = x.size();
  int mG = y.size();

  NumericMatrix z;
  // initialize return matrix with NA:
  if(as<std::string>(output)=="grid"){
    NumericMatrix z = NumericMatrix(nG,mG,NumericVector (nG*mG,NumericVector::get_na()).begin());
  }
  if(as<std::string>(output)=="points"){
    NumericMatrix z = NumericMatrix(nG,1,NumericVector (nG,NumericVector::get_na()).begin());
  }

  List ret;

  // bounding box for triangles:
  IntegerVector jTsw(nT);
  IntegerVector kTsw(nT);
  IntegerVector jTne(nT);
  IntegerVector kTne(nT);

  try {

    if(as<std::string>(output)=="grid"){
      // get bounding boxes (SW <-> NE) for all triangles:
      for(int i=0; i<nT; i++) {
        SEXP Ti = T[i];

        DataFrame Triangle(Ti);

        NumericVector xT = Triangle["x"];
        NumericVector yT = Triangle["y"];
        NumericVector zT = Triangle["z"];

        // bounding box for triangle i
        double xsw=min(xT);
        double ysw=min(yT);
        double xne=max(xT);
        double yne=max(yT);

        // translate bounding box into grid indices
        jTsw[i]=0;
        kTsw[i]=0;
        jTne[i]=nG-1;
        kTne[i]=mG-1;

        for(int j=0; j<nG; j++){
          if(x[j]<xsw) jTsw[i]=j;
          if(x[nG-j-1]>xne) jTne[i]=nG-j-1;
        }
        for(int k=0; k<mG; k++){
          if(y[k]<ysw) kTsw[i]=1;
          if(y[mG-k-1]>yne) kTne[i]=mG-k-1;
        }
      }
    }

    // iterate over triangles
    for(int i=0; i<nT; i++) {
      SEXP Ti = T[i];

      DataFrame Triangle(Ti);

      NumericVector xT = Triangle["x"];
      NumericVector yT = Triangle["y"];
      NumericVector zT = Triangle["z"];

      if(as<std::string>(output)=="grid"){
        // iterate only over grid points (j,k) inside bounding box of triangle i
        for(int j=jTsw[i]; j<jTne[i]; j++) {
          for(int k=kTsw[i]; k<kTne[i]; k++) {
            // calculate barycentric coordinates:
            double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[k] - yT[2])) /
              ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
            double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[k] - yT[2])) /
              ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
            double c = 1 - a - b;
            // check if inside triangle, handle only yet untouched grid points
            if(R_IsNA(z(j,k))){
              if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
                // perform barycentric interpolation:
                z(j,k)=a*zT[0]+b*zT[1]+c*zT[2];
              }
            }
          }
        }
      } else if(as<std::string>(output)=="points"){
        // iterate over output points
        for(int j=0; j<nG; j++) {
          // calculate barycentric coordinates:
          double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[j] - yT[2])) /
            ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
          double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[j] - yT[2])) /
            ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
          double c = 1 - a - b;
          // check if inside triangle, handle only yet untouched grid points
          if(R_IsNA(z(j,0))){
            if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
              // perform barycentric interpolation:
              z(j,0)=a*zT[0]+b*zT[1]+c*zT[2];
            }
          }
        }
      } else Rf_error("invalid output specification!");
    }

    ret=List::create(_("x")=x, _("y")=y, _("z")=z);

    return ret ;

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return List::create();             // not reached

}

/*

  plan for interpShull:

  type1: Akima (iii)
  type2: my (iii)'

  part 0:
  get triangulation (works)
  get bounding boxes for triangles to speed up code later (works)


  part 1: operate on data grid / triangles:
  if not linear:
  until no estimates missing or marked for re-estimate:

  iterate over data:
  if not yet estimated or marked for re-estimate:
  estimate up to 5 partial derivatives per data point depending on degree

  if type2:
  iterate over triangles:
  if not yet estimated or marked for re-estimate:
  estimate up to 3 directional derivatives per triangle side depending on degree

  iterate over triangles:
  calculate spline parameters
  compare polynomial derivatives with estimates
  if differences too big:
  mark data point or (if type2) triangle side for re-estmatiion with lower degree

  part 2: operate on the output grid:
  (from here on it is more or less already done)

  iterate over triangles:
  iterate over output grid or points:
  if point is in bounding box (only for grid points):
  if point is in triangle (use barycentric coordinates):
  if linear:
  convex linear combination (barymetric interpolation)
  else (spline):
  calculate polynomial


*/



// [[Rcpp::export]]
List interpShull(NumericVector x, NumericVector y,
                    NumericVector xD, NumericVector yD,
                    NumericVector zD,
                    bool linear=true,
                    CharacterVector input = "points",
                    CharacterVector output = "grid"){

  int nxD=xD.size();
  int nyD=yD.size();

  if(nxD!=nyD)
    ::Rf_error("length of xD and yD dont match!");

  int nG = x.size();
  int mG = y.size();


  NumericMatrix z;

  // initialize return matrix with NA:
  if(as<std::string>(output)=="grid"){
    z = NumericMatrix(nG,mG,NumericVector (nG*mG,NumericVector::get_na()).begin());
  }
  if(as<std::string>(output)=="points"){
    z = NumericMatrix(nG,1,NumericVector (nG,NumericVector::get_na()).begin());
  }

  List ret;

  try{

    // part 0

    // do s-Hull triangulation:
    // call shDt
    Triang tXY=shDt(Rcpp::as<std::vector<double> >(xD),
		    Rcpp::as<std::vector<double> >(yD));
    // note: triangles are enumerated counterclockwise



    int nT=tXY.nT;

    // Rcout << "get bounding boxes" <<std::endl;

    // get bounding boxes (SW <-> NE) for all triangles:

    IntegerVector jTsw(nT);
    IntegerVector kTsw(nT);
    IntegerVector jTne(nT);
    IntegerVector kTne(nT);

    IntegerVector iT = IntegerVector(3);
    NumericVector xT = NumericVector(3);
    NumericVector yT = NumericVector(3);
    NumericVector zT = NumericVector(3);

    for(int i=0; i<nT; i++) {

      iT[0]=tXY.i1[i]; iT[1]=tXY.i2[i]; iT[2]=tXY.i3[i];
      xT[0]=xD[tXY.i1[i]]; xT[1]=xD[tXY.i2[i]]; xT[2]=xD[tXY.i3[i]];
      yT[0]=yD[tXY.i1[i]]; yT[1]=yD[tXY.i2[i]]; yT[2]=yD[tXY.i3[i]];



      // bounding box for triangle i
      double xsw=min(xT);
      double ysw=min(yT);
      double xne=max(xT);
      double yne=max(yT);

      // translate bounding box into grid indices, start with complete grid:
      jTsw[i]=0;
      kTsw[i]=0;
      jTne[i]=nG-1;
      kTne[i]=mG-1;

      // Rcout << "bb to grid indices ..." << std::endl;

      for(int j=0; j<nG; j++){
	if(x[j]<xsw) jTsw[i]=j;
	if(x[nG-j-1]>xne) jTne[i]=nG-j-1;
      }
      for(int k=0; k<mG; k++){
	if(y[k]<ysw) kTsw[i]=1;
	if(y[nG-k-1]>yne) kTne[i]=mG-k-1;
      }
    }

    // part 2
    // TODO for !linear

    // part 3



    // iterate over triangles
    for(int i=0;i<nT;i++){
      iT[0]=tXY.i1[i]; iT[1]=tXY.i2[i]; iT[2]=tXY.i3[i];
      xT[0]=xD[tXY.i1[i]]; xT[1]=xD[tXY.i2[i]]; xT[2]=xD[tXY.i3[i]];
      yT[0]=yD[tXY.i1[i]]; yT[1]=yD[tXY.i2[i]]; yT[2]=yD[tXY.i3[i]];
      zT[0]=zD[tXY.i1[i]]; zT[1]=zD[tXY.i2[i]]; zT[2]=zD[tXY.i3[i]];


      if(as<std::string>(output)=="grid"){
	// iterate only over grid points (j,k) inside bounding box of triangle i
	for(int j=jTsw[i]; j<jTne[i]; j++) {
	  for(int k=kTsw[i]; k<kTne[i]; k++) {
	    // calculate barycentric coordinates:
	    double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[k] - yT[2])) /
	      ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
	    double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[k] - yT[2])) /
	      ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
	    double c = 1 - a - b;
	    // check if inside triangle, handle only yet untouched grid points
	    //if(R_IsNA(z(j,k))){
	    if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
	      z(j,k)=a*zT[0]+b*zT[1]+c*zT[2];
	    }
	  }
	}
      } else if(as<std::string>(output)=="points"){
	// iterate over output points
	for(int j=0; j<nG; j++) {
	  // calculate barycentric coordinates:
	  double a = ((yT[1] - yT[2])*(x[j] - xT[2]) + (xT[2] - xT[1])*(y[j] - yT[2])) /
	    ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
	  double b = ((yT[2] - yT[0])*(x[j] - xT[2]) + (xT[0] - xT[2])*(y[j] - yT[2])) /
	    ((yT[1] - yT[2])*(xT[0] - xT[2]) + (xT[2] - xT[1])*(yT[0] - yT[2]));
	  double c = 1 - a - b;
	  // check if inside triangle, handle only yet untouched grid points
	  //if(R_IsNA(z(j,k))){
	  if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
	    z(j,0)=a*zT[0]+b*zT[1]+c*zT[2];
	  }
	}
      }


    } // triangle


    ret=List::create(_("x")=x, _("y")=y, _("z")=z);

        return ret;

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return List::create();             // not reached

}

