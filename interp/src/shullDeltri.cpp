
#include "interp.h"


// [[Rcpp::export(name="shull.deltri")]]
List shullDeltri(NumericVector x, NumericVector y) {

  std::vector<Shx> pts;
  std::vector<Triad> triads;
  std::vector<int> outx;

  int nx=x.size();
  int ny=y.size();

  List ret;

  if(nx!=ny)
    ::Rf_error("length of x and y dont match!");

  try {

    // do s-Hull triangulation:
    // call shDt

    Triang tXYZ=shDt(Rcpp::as<std::vector<double> >(x),
		     Rcpp::as<std::vector<double> >(y));


    int nT=tXYZ.nT;

    tXYZ.xc=std::vector<double>(nT);
    tXYZ.yc=std::vector<double>(nT);
    tXYZ.rc=std::vector<double>(nT);
    tXYZ.ar=std::vector<double>(nT);
    tXYZ.rt=std::vector<double>(nT);


    for(int i=0; i<nT; i++){
      CC cc;
      cc=circum(x[tXYZ.i1[i]],y[tXYZ.i1[i]],
		x[tXYZ.i2[i]],y[tXYZ.i2[i]],
		x[tXYZ.i3[i]],y[tXYZ.i3[i]]);
      // inscribed radius: area/k , k=1/2(la+lb+lc)
      // ratio=
      double ir=cc.ar/(0.5*(sqrt((x[tXYZ.i2[i]]-x[tXYZ.i1[i]])*
				 (x[tXYZ.i2[i]]-x[tXYZ.i1[i]])+
				 (y[tXYZ.i2[i]]-y[tXYZ.i1[i]])*
				 (y[tXYZ.i2[i]]-y[tXYZ.i1[i]]))+
			    sqrt((x[tXYZ.i3[i]]-x[tXYZ.i2[i]])*
				 (x[tXYZ.i3[i]]-x[tXYZ.i2[i]])+
				 (y[tXYZ.i3[i]]-y[tXYZ.i2[i]])*
				 (y[tXYZ.i3[i]]-y[tXYZ.i2[i]]))+
			    sqrt((x[tXYZ.i1[i]]-x[tXYZ.i3[i]])*
				 (x[tXYZ.i1[i]]-x[tXYZ.i3[i]])+
				 (y[tXYZ.i1[i]]-y[tXYZ.i3[i]])*
				 (y[tXYZ.i1[i]]-y[tXYZ.i3[i]]))));
      tXYZ.rt[i]=ir/cc.rc;
      
      tXYZ.i1[i]++; // start enumeration with 1 in R, ATTENTION:
      tXYZ.i2[i]++; // from now to exiting this function do not access
      tXYZ.i3[i]++; // x or y elements through  tXYZ.i1, tXYZ.i2, tXYZ.i3 !!!
      tXYZ.xc[i]=cc.xc;
      tXYZ.yc[i]=cc.yc;
      tXYZ.rc[i]=cc.rc;
      tXYZ.ar[i]=cc.ar;
    }
    // get convex hull and arcs
    std::vector<int> cp1=std::vector<int>(nx);
    std::vector<int> cp2=std::vector<int>(nx);
    // count
    int nCH=0;
    // check if neigbour triangle is not present (-1), means that
    // arc is part of the convex hull,
    for(int i=0; i<nT; i++){
      // store arcs on convexhull if triangle neighbour is indexed as -1
      if(tXYZ.j1[i]==-1){
        cp1[nCH]=tXYZ.i2[i];  
        cp2[nCH]=tXYZ.i3[i];  
        nCH++;
        }
      if(tXYZ.j2[i]==-1){
        cp1[nCH]=tXYZ.i3[i];  
        cp2[nCH]=tXYZ.i1[i];  
        nCH++;
      }
      if(tXYZ.j3[i]==-1){
        cp1[nCH]=tXYZ.i1[i]; 
        cp2[nCH]=tXYZ.i2[i]; 
        nCH++;
      }
    }
    //Rcout << "nCH=" << nCH << std::endl;

    // initialize in correct size
    tXYZ.ch=std::vector<int>(nCH);
    tXYZ.nch=nCH;
    // number of arcs: only arcs on the convex hull counted once,
    // others twice
    int nArcs=(3*nT+nCH)/2;
    tXYZ.a1=std::vector<int>(nArcs);
    tXYZ.a2=std::vector<int>(nArcs);
    tXYZ.k1=std::vector<int>(nT);
    tXYZ.k2=std::vector<int>(nT);
    tXYZ.k3=std::vector<int>(nT);
    tXYZ.na=nArcs;

    // get arcs
    int ia=0;
    for(int i=0; i<nT; i++){
      if(ia>nArcs)
        Rf_error("error counting arcs!");

      // store arcs if not already done:
      bool found=false;
      for(int j=0; j<ia; j++){
        if((tXYZ.a1[j]==tXYZ.i3[i]) & (tXYZ.a2[j]==tXYZ.i2[i])){ 
          found=true;
          tXYZ.k1[i]=j;
          break;
        }
      }
      if(!found){
        tXYZ.a1[ia]=tXYZ.i2[i];  
        tXYZ.a2[ia]=tXYZ.i3[i];  
        tXYZ.k1[i]=ia;
        ia++;
      }
      found=false;
      for(int j=0; j<ia; j++){
        if((tXYZ.a1[j]==tXYZ.i1[i]) & (tXYZ.a2[j]==tXYZ.i3[i])){ 
          found=true;
          tXYZ.k2[i]=j;
          break;
        }
      }
      if(!found){
          tXYZ.a1[ia]=tXYZ.i3[i];
          tXYZ.a2[ia]=tXYZ.i1[i];
          tXYZ.k2[i]=ia;
          ia++;
      }
      found=false;
      for(int j=0; j<ia; j++){
        if((tXYZ.a1[j]==tXYZ.i2[i]) & (tXYZ.a2[j]==tXYZ.i1[i])){
          found=true;
          tXYZ.k3[i]=j;
          break;
        }
      }
      if(!found){
          tXYZ.a1[ia]=tXYZ.i1[i];
          tXYZ.a2[ia]=tXYZ.i2[i];
          tXYZ.k3[i]=ia;
          ia++;
      }
    }


    // first arc on convex hull


    tXYZ.ch[0]=cp1[0];
    tXYZ.ch[1]=cp2[0];
    int j=1;

    while(j<nCH-1){
      // search next arc starting with last endpoint
      for(int i=1; i<nCH; i++){
        if(tXYZ.ch[j]==cp1[i]){
          tXYZ.ch[j+1]=cp2[i];
        }
      }
      j++;
    }
    /*
    for(int i=0; i<tXYZ.nch; i++){
      tXYZ.ch[i]++; // start enumeration with 1 in R
    }
    */
    IntegerMatrix trlist=IntegerMatrix(nT,9);
    for(int i=0; i<nT; i++){
      trlist(i,0)=tXYZ.i1[i];
      trlist(i,1)=tXYZ.i2[i];
      trlist(i,2)=tXYZ.i3[i];
      trlist(i,3)=tXYZ.j1[i]+1; // index shift 0 -> 1
      trlist(i,4)=tXYZ.j2[i]+1;
      trlist(i,5)=tXYZ.j3[i]+1;
      trlist(i,6)=tXYZ.k1[i]+1;
      trlist(i,7)=tXYZ.k2[i]+1;
      trlist(i,8)=tXYZ.k3[i]+1;
    }
    NumericMatrix cclist=NumericMatrix(nT,5);
    for(int i=0; i<nT; i++){
      cclist(i,0)=tXYZ.xc[i];
      cclist(i,1)=tXYZ.yc[i];
      cclist(i,2)=tXYZ.rc[i];
      cclist(i,3)=tXYZ.ar[i];
      cclist(i,4)=tXYZ.rt[i];
    }

    ret=List::create(_("n")=nx, _("x")=x, _("y")=y,
		     _("nt")=nT,
		     _("trlist")=trlist,
		     _("cclist")=cclist,
		     _("nch")=tXYZ.nch, _("ch")=tXYZ.ch,
		     _("na")=tXYZ.na, _("a1")=tXYZ.a1,
		     _("a2")=tXYZ.a2);

    return ret ;

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return List::create();             // not reached

}


Triang shDt(std::vector<double> x, std::vector<double> y){

  // Note: circumcircles and convex hull only done in shullDeltri
  //       as this is not needed for the application within Akimas
  //       spline routines.

  Triang Txy;

  std::vector<Shx> pts;
  std::vector<Triad> triads;
  std::vector<int> outx;

  int nx=x.size();
  int ny=y.size();

  if(nx!=ny)
    ::Rf_error("length of x and y dont match!");

  try {

    // triangulation

    // Rcout << "start triangulation" << std::endl;

    for(int i=0; i<nx; i++){
      Shx pt;
      pt.id=i;
      pt.r=x[i];
      pt.c=y[i];
      pts.push_back(pt);
    }

    // Note: points are already deduplicated in R!
    int n_dups = de_duplicate(pts, outx);
    if (n_dups != 0) {
      stop("shull: duplicate points found");
    }
    int ierr = s_hull_pro(pts, triads);
    if (ierr != 1) {
      if (ierr == -1) {
        stop("shull: less than 3 points, aborting");
      } else if (ierr == -2) {
        stop("shull: linear structure, aborting");
      } else if (ierr == -3) {
        stop("shull: cannot triangulate this set");
      } else if (ierr == -4) {
        stop("shull: cannot triangulate this set");
      } else if (ierr == -5) {
        stop("shull: triangle flipping error");
      } else if (ierr == -6) {
        stop("shull: triangle flipping error");
      } else {
        stop("shull: unspecified error");
      }
    }

    int nT = triads.size();

    Txy.i1=std::vector<int>(nT);
    Txy.i2=std::vector<int>(nT);
    Txy.i3=std::vector<int>(nT);
    Txy.j1=std::vector<int>(nT);
    Txy.j2=std::vector<int>(nT);
    Txy.j3=std::vector<int>(nT);
    // not used here, dummy allocation
    Txy.k1=std::vector<int>(1);
    Txy.k2=std::vector<int>(1);
    Txy.k3=std::vector<int>(1);
    Txy.ch=std::vector<int>(1);
    Txy.a1=std::vector<int>(1);
    Txy.a2=std::vector<int>(1);
    Txy.xc=std::vector<double>(1);
    Txy.yc=std::vector<double>(1);
    Txy.rc=std::vector<double>(1);
    Txy.ar=std::vector<double>(1);
    Txy.rt=std::vector<double>(1);


    for(int i=0; i<nT; i++){
      // check for counter clockwise orientation:
      double orientation=(x[triads[i].c]-x[triads[i].b])*(y[triads[i].b]-y[triads[i].a]) +
        (y[triads[i].c]-y[triads[i].b])*(x[triads[i].a]-x[triads[i].b]);
      if(orientation==0.0)
        Rf_error("triangle collapsed!");
      //Rcout << "tr " << i << " or: " << orientation << std::endl;
      if(orientation<0.0){
        // a,b,c is already counter clockwise:
        Txy.i1[i]=triads[i].a;
        Txy.i2[i]=triads[i].b;
        Txy.i3[i]=triads[i].c;
        Txy.j1[i]=triads[i].bc;
        Txy.j2[i]=triads[i].ac;
        Txy.j3[i]=triads[i].ab;
      } else {
        // force counter clockwise
        Txy.i1[i]=triads[i].a;
        Txy.i2[i]=triads[i].c;
        Txy.i3[i]=triads[i].b;
        Txy.j1[i]=triads[i].bc ;
        Txy.j2[i]=triads[i].ab;
        Txy.j3[i]=triads[i].ac;
      }
    }
    Txy.nT=nT;


    return Txy;

  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  // not reached:
  Txy.nT=0;
  return Txy;
}



CC circum(double x1,double y1, double x2,double y2, double x3,double y3){
  // https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
  CC ret;


  double D = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
  if(D==0)
    Rf_error("three points coincide or are collinear!");
  ret.xc =
(((x1 - x3) * (x1 + x3) + (y1 - y3) * (y1 + y3)) / 2 * (y2 - y3)
    -  ((x2 - x3) * (x2 + x3) + (y2 - y3) * (y2 + y3)) / 2 * (y1 - y3))
    / D;

  ret.yc =
(((x2 - x3) * (x2 + x3) + (y2 - y3) * (y2 + y3)) / 2 * (x1 - x3)
    -  ((x1 - x3) * (x1 + x3) + (y1 - y3) * (y1 + y3)) / 2 * (x2 - x3))
    / D;

  ret.rc = sqrt((x1 - ret.xc)*(x1 - ret.xc) + (y1 - ret.yc)*(y1 - ret.yc));

  // http://www.mathopenref.com/coordtrianglearea.html
  ret.ar = std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0;

  return ret;
}

// [[Rcpp::export]]
List triFind(int nT, NumericVector xT, NumericVector yT,
             IntegerVector i1, IntegerVector i2, IntegerVector i3,
             NumericVector x, NumericVector y){

  int n=x.size();
  NumericVector li1=NumericVector(n), li2=NumericVector(n),
    li3=NumericVector(n), tr=NumericVector(n);
  NumericMatrix bc=NumericMatrix(n,3);

  for(int j=0; j<n; j++){
    li1[j]=0;
    li2[j]=0;
    li3[j]=0;
  }
  for(int j=0; j<n; j++){
    for(int i=0; i<nT; i++){
      // get barycentric coordinates
      double a = ((yT[i2[i]-1] - yT[i3[i]-1])*(x[j] - xT[i3[i]-1]) + (xT[i3[i]-1] - xT[i2[i]-1])*(y[j] - yT[i3[i]-1])) /
      ((yT[i2[i]-1] - yT[i3[i]-1])*(xT[i1[i]-1] - xT[i3[i]-1]) + (xT[i3[i]-1] - xT[i2[i]-1])*(yT[i1[i]-1] - yT[i3[i]-1]));
      double b = ((yT[i3[i]-1] - yT[i1[i]-1])*(x[j] - xT[i3[i]-1]) + (xT[i1[i]-1] - xT[i3[i]-1])*(y[j] - yT[i3[i]-1])) /
        ((yT[i2[i]-1] - yT[i3[i]-1])*(xT[i1[i]-1] - xT[i3[i]-1]) + (xT[i3[i]-1] - xT[i2[i]-1])*(yT[i1[i]-1] - yT[i3[i]-1]));

      double c = 1 - a - b;
      bc(j,0)=a;
      bc(j,1)=b;
      bc(j,2)=c;

      //Rcout << "tr: " << i << " (a,b,c): (" << a << ", " << b << ", " << c << ")" << std::endl;
      if(0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1){
        // triangle found
        li1[j]=i1[i];
        li2[j]=i2[i];
        li3[j]=i3[i];
        tr[j]=i;
        break;
      }
    }
  }
  return List::create(_("i1")=li1, _("i2")=li2, _("i3")=li3, _("tr")=tr, _("bc")=bc);
}

// [[Rcpp::export]]
LogicalVector left(double x1,double y1, double x2, double y2,
		   NumericVector x0, NumericVector y0, double eps=1E-16){

  int n=x0.size();
  LogicalVector ret(n);

  for(int i=0; i<n; i++){
    ret[i]=((x2-x1)*(y0[i]-y1)-(x0[i]-x1)*(y2-y1)>=eps);
    //ret[i]=((x2-x1)*(y0[i]-y1)>=(x0[i]-x1)*(y2-y1));
   }

  return ret;
}

// [[Rcpp::export]]
LogicalVector on(double x1,double y1, double x2, double y2,
		 NumericVector x0, NumericVector y0, double eps=1E-16){

  int n=x0.size();
  LogicalVector ret(n);

  for(int i=0; i<n; i++){
    ret[i]=(std::abs(((x2-x1)*(y0[i]-y1)-(x0[i]-x1)*(y2-y1)))<eps);
   }

  return ret;
}

// [[Rcpp::export]]
LogicalVector inHull(List triObj,
		     NumericVector x, NumericVector y, double eps=1E-16){

  int n=x.size();
  LogicalVector ret(n);
  for(int i=0; i<n; i++){
    ret[i]=true;
  }
  List tO(triObj);

  int nCH=tO["nchull"];
  NumericVector xD = tO["x"];
  NumericVector yD = tO["y"];
  IntegerVector ch = tO["chull"]; // starts with 1 ! use index shift -1 !!!
  LogicalVector lft;


  for(int j=0; j<nCH; j++){
    if(j<nCH-1){
      lft=left(xD[ch[j]-1],yD[ch[j]-1],xD[ch[j+1]-1],yD[ch[j+1]-1],x,y,eps);
      for(int i=0; i<n; i++){
	ret[i]=ret[i] & lft[i];
      }
    } else {
      // last segment closing the hull:
      lft=left(xD[ch[j]-1],yD[ch[j]-1],xD[ch[0]-1],yD[ch[0]-1],x,y,eps);
      for(int i=0; i<n; i++){
	ret[i]=ret[i] & lft[i];
      }
    }
  }

  return ret;
}

// [[Rcpp::export]]
LogicalVector onHull(List triObj,
		     NumericVector x, NumericVector y, double eps=1E-16){

  int n=x.size();
  LogicalVector ret(n);
  for(int i=0; i<n; i++){
    ret[i]=false;
  }
  List tO(triObj);

  int nCH=tO["nchull"];
  NumericVector xD = tO["x"];
  NumericVector yD = tO["y"];
  IntegerVector ch = tO["chull"]; // starts with 1 ! use index shift -1 !!!
  LogicalVector onH;


  for(int j=0; j<nCH; j++){
    if(j<nCH-1){
      onH=on(xD[ch[j]-1],yD[ch[j]-1],xD[ch[j+1]-1],yD[ch[j+1]-1],x,y,eps);
      for(int i=0; i<n; i++){
	ret[i]=ret[i] | onH[i];
      }
    } else {
      onH=on(xD[ch[j]-1],yD[ch[j]-1],xD[ch[0]-1],yD[ch[0]-1],x,y,eps);
      for(int i=0; i<n; i++){
	ret[i]=ret[i] | onH[i];
      }
    }
  }

  return ret;
}
