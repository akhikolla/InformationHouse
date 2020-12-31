
#include "interp.h"

// [[Rcpp::export(name="locpoly.partderiv.grid")]]
List partDerivGrid(NumericVector x, NumericVector y,
		   NumericVector xD, NumericVector yD,
		   NumericVector zD,
		   CharacterVector kernel="gaussian",
		   NumericVector h=NumericVector::create(0.25,0.25),
		   CharacterVector solver="QR",
		   int degree=3,
		   bool smoothpde=false,
		   bool akimaweight=false,
		   int nweight=25) {
  // Estimate up to third order partial derivatives at x,y locations:
  // apply local polynomial regression of order up to 3


  List ret;
  int nD = xD.size();
  int nG = x.size();
  int mG = y.size();
  
  int p=0;
  if(degree==0)
    p=1; // local constant trend
  if(degree==1)
    p=3; // local linear trend
  else if(degree==2)
    p=6; // local quadratic trend
  else if(degree==3)
    p=10; // local cubic trend
  else if(degree>3)
    Rf_error("degree>3 !");

  // initialize return matrices
  NumericMatrix Ze = NumericMatrix(nG,mG);
  NumericMatrix Zx = NumericMatrix(nG,mG);
  NumericMatrix Zy = NumericMatrix(nG,mG);


  NumericMatrix Zxy;
  NumericMatrix Zx2;
  NumericMatrix Zy2;

  if(degree>=2){
    Zxy = NumericMatrix(nG,mG);
    Zx2 = NumericMatrix(nG,mG);
    Zy2 = NumericMatrix(nG,mG);
  }

  NumericMatrix Zx2y;
  NumericMatrix Zxy2;
  NumericMatrix Zx3;
  NumericMatrix Zy3;

  if(degree>=3){
    Zx2y = NumericMatrix(nG,mG);
    Zxy2 = NumericMatrix(nG,mG);
    Zx3  = NumericMatrix(nG,mG);
    Zy3  = NumericMatrix(nG,mG);
  }

  //Rcout << "size is " << nD << std::endl;
  PDEst pde;
  pde.est=VectorXd(p);  // do better in pd*?
  pde.se=VectorXd(p);
  
  NN nn=nN(xD,yD);

  for(int i=0; i<nG; i++){
    for(int j=0; j<mG; j++){
      if(smoothpde)
	pde=pDsmooth(xD,yD,zD,nn, x[i],y[j],kernel,h,as<std::string>(solver),degree,nweight,akimaweight);
      else
	pde=pD(xD,yD,zD,nn, x[i],y[j],kernel,h,as<std::string>(solver),degree);

      // extract partial derivatives from betahat

      //Rcout << "betahat " << std::endl;
      //Rcout <<  betahat   << std::endl;
      //Rcout << "#### END #### " << i << std::endl;

      Ze(i,j)  =  pde.est(0); // local estimate, not really needed, but
      // don't throw it away, use it for checking
      // get partial derivatives by using betahat for Taylor series:
      if(degree>=1){
        Zx(i,j)  =  pde.est[1];
        Zy(i,j)  =  pde.est[2];
      }
      if(degree>=2){
        Zxy(i,j) =  pde.est[3];
        Zx2(i,j) =  pde.est[4];
        Zy2(i,j) =  pde.est[5];
      }
      if(degree>=3){
        Zx2y(i,j) = pde.est[6];
        Zxy2(i,j) = pde.est[7];
        Zx3(i,j) =  pde.est[8];
        Zy3(i,j) =  pde.est[9];
      }
    }
  }

  if(degree==0){
    ret=List::create(_("z")=Ze);
  }

  if(degree==1){
    ret=List::create(_("z")=Ze, _("zx")=Zx, _("zy")=Zy);
  }

  if(degree==2){
    ret=List::create(_("z")=Ze, _("zx")=Zx, _("zy")=Zy,
                     _("zxy")=Zxy, _("zxx")=Zx2, _("zyy")=Zy2);
  }
  if(degree==3){
    ret=List::create(_("z")=Ze, _("zx")=Zx, _("zy")=Zy,
                     _("zxy")=Zxy, _("zxx")=Zx2, _("zyy")=Zy2,
                     _("zxxy")=Zx2y, _("zxyy")=Zxy2, _("zxxx")=Zx3, _("zyyy")=Zy3);
  }
  return ret;
}



// [[Rcpp::export(name="locpoly.partderiv.points")]]
List partDerivPoints(NumericVector x, NumericVector y,
		     NumericVector xD, NumericVector yD, NumericVector zD,
		     CharacterVector kernel="gaussian",
		     NumericVector h=NumericVector::create(0.25,0.25),
		     CharacterVector solver="QR",
		     int degree=3,
		     bool smoothpde=false,
		     bool akimaweight=false,
		     int nweight=25) {
  // Estimate up to third order partial derivatives at data grid locations:
  // apply local polynomial regression of order up to 3

  
  List ret;
  int nD = xD.size();
  int nP = x.size();


  
  int p=0;

  if(degree==0)
    p=1; // local constant trend
  if(degree==1)
    p=3; // local linear trend
  else if(degree==2)
    p=6; // local quadratic trend
  else if(degree==3)
    p=10; // local cubic trend
  else if(degree>3)
    Rf_error("degree>3 !");
  
  // initialize return vectors
  NumericVector Ze = NumericVector(nP);
  NumericVector Zx = NumericVector(nP);
  NumericVector Zy = NumericVector(nP);

  NumericVector Zxy;
  NumericVector Zx2;
  NumericVector Zy2;

  if(degree>=2){
    Zxy = NumericVector(nP);
    Zx2 = NumericVector(nP);
    Zy2 = NumericVector(nP);
  }

  NumericVector Zx2y;
  NumericVector Zxy2;
  NumericVector Zx3;
  NumericVector Zy3;

  if(degree>=3){
    Zx2y = NumericVector(nP);
    Zxy2 = NumericVector(nP);
    Zx3  = NumericVector(nP);
    Zy3  = NumericVector(nP);
  }

  PDEst pde;
  pde.est=VectorXd(p);  // do better in pd*?
  pde.se=VectorXd(p);

  NN nn=nN(xD,yD);

  for(int i=0; i<nP; i++){
    if(smoothpde)
      pde=pDsmooth(xD,yD,zD,nn, xD[i],yD[i],kernel,h,as<std::string>(solver),degree,nweight,akimaweight);
    else
      pde=pD(xD,yD,zD,nn, xD[i],yD[i],kernel,h,as<std::string>(solver),degree);

    Ze[i]  =  pde.est[0]; // local estimate, not really needed, but
    // don't throw it away, use it for checking
    // get partial derivatives by using betahat for Taylor series:

    if(degree>=1){
      Zx[i]  =  pde.est[1];
      Zy[i]  =  pde.est[2];
    }
    if(degree>=2){
      Zxy[i] =  pde.est[3];
      Zx2[i] =  pde.est[4];
      Zy2[i] =  pde.est[5];
    }
    if(degree>=3){
      Zx2y[i] = pde.est[6];
      Zxy2[i] = pde.est[7];
      Zx3[i] = pde.est[8];
      Zy3[i] = pde.est[9];
    }
  }
  
  if(degree==0){
    ret=List::create(_("z")=Ze);
  }

  if(degree==1){
    ret=List::create(_("z")=Ze, _("zx")=Zx, _("zy")=Zy);
  }

  if(degree==2){
    ret=List::create(_("z")=Ze, _("zx")=Zx, _("zy")=Zy,
                     _("zxy")=Zxy, _("zxx")=Zx2, _("zyy")=Zy2);
  }
  if(degree==3){
    ret=List::create(_("z")=Ze, _("zx")=Zx, _("zy")=Zy,
                     _("zxy")=Zxy, _("zxx")=Zx2, _("zyy")=Zy2,
                     _("zxxy")=Zx2y, _("zxyy")=Zxy2, _("zxxx")=Zx3, _("zyyy")=Zy3);
  }
  return ret;
}


PDEst pD(NumericVector xD, NumericVector yD, NumericVector zD, NN nn,
         double x, double y, CharacterVector kernel, NumericVector h,
         std::string solver, int degree){


  int nD=xD.size();
  double xRange=max(xD)-min(xD);
  double yRange=max(yD)-min(yD);
  if((h.size()!=2) & (h.size()!=1))
    Rf_error("bandwidth parameter h is not a vector of 2 or 1 elements!");
  double bwX, bwY;
  // global bandwidth:
  if(h.size()==2){
    bwX=h[0]*xRange;
    bwY=h[1]*yRange;
    //Rcout << "global bw: (" << bwX << ", " << bwY << ")" << std::endl;
  }
  // initialize nearest neigbour structure for local bandwidth:
  NN lnn;

  if(h.size()==1){

    NumericVector xtmp(1);
    xtmp[0]=x;
    NumericVector ytmp(1);
    ytmp[0]=y;
    // FIXME: for partDerivData only one call to nN is necessary,
    // outside this for loop!!! this generates the runtime difference to
    // global bandwidth!!
    lnn=extendNN(nn, xD,yD,xtmp,ytmp);
    //Rcout << "distance matrix" << std::endl;
    //Rcout << nn.ind << std::endl;
    //Rcout << nn.dist << std::endl;
  }

  //Rcout << "data point " << i << std::endl;

  // setup design matrix,
  // 3, 6 or 10 columns for 1st, 2nd or 3rd degree bivariate polynomial:
  // X=(1, (x-x0), (y-y0),
  //   (x-x0)(y-y0), (x-x0)^2, (y-y0)^2,
  //   (x-x0)^2(y-y0), (x-x0)(y-y0)^2, (x-x0)^3, (y-y0)^3)
  int p;
  if(degree==0)
    p=1; // local constant trend
  else if(degree==1)
    p=3; // local linear trend
  else if(degree==2)
    p=6; // local quadratic trend
  else if(degree==3)
    p=10; // local cubic trend

  MatrixXd X(nD,p);
  for(int j=0; j<nD; j++){
    X(j,0)=1.0L;
    if(degree>=1){
      X(j,1)=x-xD[j];
      X(j,2)=y-yD[j];
      if(degree>=2){
        X(j,3)=(x-xD[j])*(y-yD[j]);
        X(j,4)=(x-xD[j])*(x-xD[j]);
        X(j,5)=(y-yD[j])*(y-yD[j]);
      }
      if(degree>=3){
        X(j,6)=(x-xD[j])*(x-xD[j])*(y-yD[j]);
        X(j,7)=(x-xD[j])*(y-yD[j])*(y-yD[j]);
        X(j,8)=(x-xD[j])*(x-xD[j])*(x-xD[j]);
        X(j,9)=(y-yD[j])*(y-yD[j])*(y-yD[j]);
      }
    }
  }




  // build diagonal weight matrix, better use Diagonal matrix type

  Eigen::DiagonalMatrix<double,Eigen::Dynamic> W(nD);


  // local bandwidth:
  if(h.size()==1){
    if(h[0]>1)
      Rf_error("local bandwidth parameter >1 !");
    int nnX=h[0]*nD+1; // +1: the actual point has a duplicate!
    // the 2nd order local polynomial needs at least 6 data locations:
    if(nnX<=p){
      // Rf_warning("local bandwidth parameter to small, increasing");
      nnX=min(IntegerVector(nD,p));
    }
    if(nnX==nD)
      nnX=nD-1;
    bwX=lnn.dist(0,nnX); // FIXME: use lnn.dist() ??????????
    bwY=bwX;
    //Rcout << "local bw: (" << bwX << ", " << bwY << ")" << std::endl;
  }


  for(int j=0; j<nD; j++){
    // fill in sqrt(w_j), to be able to use transformed data
    // W^0.5*X and W^0.5*y later:
    W.diagonal()(j)=sqrt(kern2d(xD[j], x, bwX,
                                yD[j], y, bwY,
                                as<std::string>(kernel)));
  }


  //Rcout << "W^0.5 is " << std::endl;
  //Rcout << W << std::endl;


  // replace Xm' with X'*W^0.5 to get weighted least squares:
  const MatrixXd Xm(W*X);

  // Rcout << "design matrix is" << std::endl;
  // Rcout << Xh << std::endl;

  // Rcout << "Wm is " << std::endl;
  // Rcout << Wm << std::endl;

  // solve normal equations, use:
  // https://cran.r-project.org/web/packages/RcppEigen/
  // vignettes/RcppEigen-Introduction.pdf section 4:

  // replace y with W^0.5*y to get weighted least squares:
  const VectorXd yd(W*as<MapVecd>(zD));

  // Rcout << "ys is " << std::endl;
  // Rcout << yd << std::endl;


  const int n(Xm.rows());//, p(Xm.cols());

  PDEst pde;
  pde.betahat=VectorXd(p);
  for(int i=0;i<p;i++){
    pde.betahat[i]=0.0;
  }
  pde.est=VectorXd(p);
  pde.se=VectorXd(p);
  Eigen::JacobiSVD<MatrixXd> svdXm(Xm);
  pde.cond = svdXm.singularValues()(0)
    / svdXm.singularValues()(svdXm.singularValues().size()-1);

  if(solver=="LLt"){
    // this is the LLt Cholesky solver, section 4.1 of
    // https://cran.r-project.org/web/packages/RcppEigen/vignettes/RcppEigen-Introduction.pdf:
    const LLT<MatrixXd> llt(AtA(Xm));
    pde.betahat=llt.solve(Xm.adjoint() * yd);
    const VectorXd fitted(Xm * pde.betahat);
    const VectorXd resid(yd - fitted);
    const int df(n - p);
    const double s(resid.norm() / std::sqrt(double(df)));
    // FIXME
    //pde.se=W.inverse()*(s * llt.matrixL().solve(MatrixXd::Identity(p, p))
    //                    .colwise().norm());
  } else
    // if(solver=="LDLt"){
    // } else
    if(solver=="QR"){
      // this is the unpivoted QR decomposition solver, section 4.2
      const HouseholderQR<MatrixXd> QR(Xm);
      pde.betahat=QR.solve(yd);
      //const VectorXd fitted(Xm * pde.betahat);
      //const int df(n - p);
      // FIXME: memory errors detected by ASAN and valgrind:
      //pde.se=W.inverse()*(QR.matrixQR().topRows(p).triangularView<Upper>()
      //                    .solve(MatrixXd::Identity(p,p)).rowwise().norm());
    } else
      if(solver=="SVD"){
        // this is the SVD based solver, section 4.4:
        const Eigen::JacobiSVD<MatrixXd> UDV(Xm.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV));
        const ArrayXd Dp(Dplus(UDV.singularValues()));
        const int r((Dp > 0).count());
        const MatrixXd VDp(UDV.matrixV() * Dp.matrix().asDiagonal());
        pde.betahat=VDp * UDV.matrixU().adjoint() * yd;
        const VectorXd fitted(Xm * pde.betahat);
        const VectorXd resid(yd - fitted);
        const int df(nD - p);
        const double s(resid.norm() / std::sqrt(double(df)));
	// FIXME
        // pde.se=W.inverse()*(s * VDp.rowwise().norm());
      } else
        if(solver=="Eigen"){
          // this is the eigen decomposition based solver, section 4.5:
          const Eigen::SelfAdjointEigenSolver<MatrixXd> VLV(AtA(Xm));
          const ArrayXd Dp(Dplus(VLV.eigenvalues()).sqrt());
          const int r((Dp > 0).count());
          const MatrixXd VDp(VLV.eigenvectors() * Dp.matrix().asDiagonal());
          pde.betahat=VDp * VDp.adjoint() * Xm.adjoint() * yd;
          const VectorXd fitted(Xm * pde.betahat);
          const VectorXd resid(yd - fitted);
          const int df(nD - p);
          const double s(resid.norm() / std::sqrt(double(df)));
	  // FIXME
          //pde.se=W.inverse()*(s * VDp.rowwise().norm());
        } else
          if(solver=="CPivQR"){
            // this is the column based pivoted QR solver, section 4.6:
            const CPivQR PQR(Xm);
            const Permutation Pmat(PQR.colsPermutation());
            const int r(PQR.rank());
            VectorXd fitted, se;
            if (r == Xm.cols()) { // full rank case
              // Rcout << "pQR full rank" << std::endl;
              pde.betahat = PQR.solve(yd);
              //fitted = Xm * pde.betahat;
	      // FIXME
              //pde.se = W.inverse() * (Pmat * PQR.matrixQR().topRows(p).triangularView<Upper>()
              //                        .solve(MatrixXd::Identity(p, p)).rowwise().norm());
            } else {
              // Rcout << "pQR no full rank " << r << " < " << Xm.cols() << std::endl;
              MatrixXd Rinv(PQR.matrixQR().topLeftCorner(r, r)
                            .triangularView<Upper>().solve(MatrixXd::Identity(r, r)));
              VectorXd effects(PQR.householderQ().adjoint() * yd);
              pde.betahat.fill(::NA_REAL);
              pde.betahat.head(r) = Rinv * effects.head(r);
              pde.betahat = Pmat * pde.betahat;
	      // FIXME
              //se.fill(::NA_REAL);
              //se.head(r) = Rinv.rowwise().norm();
              //se = W.inverse() * (Pmat * se);
              // create fitted values from effects
              //effects.tail(Xm.rows() - r).setZero();
              //fitted = PQR.householderQ() * effects;
            }
          } else
            Rf_error("unknown solver");
  /* results, we use only pde.betahat and se for now:
     return List::create(Named("coefficients") = pde.betahat,
     Named("fitted.values") = fitted,
     Named("residuals") = resid,
     Named("s") = s,
     Named("df.residual") = df,
     Named("rank") = p,
     Named("Std. Error") = se);
  */
  /*
    Rcout << "pde.est:" << std::endl;
    Rcout << pde.est << std::endl;
    Rcout << "pde.se:" << std::endl;
    Rcout << pde.se << std::endl;
  */
  // prepare return vector
  pde.est[0] =  pde.betahat[0]; // local estimate, not really needed, but
  // don't throw it away, use it for checking
  // get partial derivatives by treating pde.betahat as Taylor series coefficients:
  if(degree>=1){
    pde.est[1]  = -pde.betahat[1]; // -1=-1 * 1! because we use 1/n!(x_0-x)^n in taylor
    pde.est[2]  = -pde.betahat[2]; //            series so we have: 1/n!(-1)^n(x-x_0)^n
  }
  if(degree>=2){
    pde.est[3] =    pde.betahat(3); // factor 1= 1!*1!
    pde.est[4] = 2L*pde.betahat(4); // factor 2= 2!*0!
    pde.est[5] = 2L*pde.betahat(5); // factor 2= 0!*2!
  }
  if(degree>=3){ // not used for akima splines, but keep it anyway:
    pde.est[6] = -2L*pde.betahat(6); // factor 2= 2!*1!, "-" see above
    pde.est[7] = -2L*pde.betahat(7); // factor 2= 1!*2!
    pde.est[8] = -6L*pde.betahat(8); // factor 6= 3!*0!
    pde.est[9] = -6L*pde.betahat(9); // factor 6= 0!*3!
  }

  return pde;
}

PDEst pDsmooth(NumericVector xD, NumericVector yD, NumericVector zD, NN nn,
               double x, double y, CharacterVector kernel, NumericVector h,
               std::string solver, int degree, int n, bool akimaweight){
  // estimate derivatives for up to n (or better p) nearest neighbours,
  // return average according to Akimas weigthing scheme

  int p;
  if(degree==0)
    p=1; // local constant trend,
  //FIXME: use Akimas plane with only two points, e.g. by
  //else if(degree==0.5)
  //  p=2; //
  else if(degree==1)
    p=3; // local linear trend
  else if(degree==2)
    p=6; // local quadratic trend
  else if(degree==3)
    p=10; // local cubic trend
  if(n==0){
    n=p;
  } else {
    if(n>xD.size())
      n=xD.size();
  }


  VectorXd Z(n), L(n);
  VectorXd Zx(n), Lx(n);
  VectorXd Zy(n), Ly(n);

  VectorXd Zxy(n);
  VectorXd Zxx(n);
  VectorXd Zyy(n);

  VectorXd Zxxy(n);
  VectorXd Zxyy(n);
  VectorXd Zxxx(n);
  VectorXd Zyyy(n);
  /*
    NN lnn;

    if(h.size()==1){

    NumericVector xtmp(1);
    xtmp[0]=x;
    NumericVector ytmp(1);
    ytmp[0]=y;
    // FIXME: for partDerivData only one call to nN is necessary,
    // outside this for loop!!! this generates the runtime difference to
    // global bandwidth!!
    lnn=extendNN(nn, xD,yD,xtmp,ytmp);
    Rcout << "distance matrix" << std::endl;
    Rcout << lnn.ind << std::endl;
    Rcout << lnn.dist << std::endl;
    }
  */
  // TODO
  PDEst pde=pD(xD,yD,zD,nn,x,y,kernel,h,solver,degree);

  PDEst lsfit=pD(xD,yD,zD,nn,x,y,kernel,h,solver,1);

  // Rcout << " pde: " << pde.est << std::endl;
  for(int i=0;i<n;i++){
    if(i==0){
      Z[i]=pde.est[0];
      L[i]=lsfit.est[0];
      if(degree>0){
        Zx[i]=pde.est[1];
        Zy[i]=pde.est[2];
        Lx[i]=lsfit.est[1];
        Ly[i]=lsfit.est[2];
        if(degree>1){
          Zxy[i]=pde.est[3];
          Zxx[i]=pde.est[4];
          Zyy[i]=pde.est[5];
          if(degree>2){
            Zxxy[i]=pde.est[6];
            Zxyy[i]=pde.est[7];
            Zxxx[i]=pde.est[8];
            Zyyy[i]=pde.est[9];
          }
        }
      }
    } else {
      double lx,ly;

      lx=xD[nn.ind(0,i)];
      ly=yD[nn.ind(0,i)];

      // prepare vectors for later calculation of derivatives via scalar product
      // reproduce in Maxima with
      /*

        beta:[b0,b1,b2,b3,b4,b5,b6,b7,b8,b9]$
        fvec(x,y):=[1,(x-x0),(y-y0),(x-x0)*(y-y0),(x-x0)^2,(y-y0)^2,(x-x0)^2*(y-y0),(x-x0)*(y-y0)^2,(x-x0)^3,(y-y0)^3]$
        pol(x,y):=beta.fvec(x,y)$
        diff(pol(x,y),x,1);
        diff(pol(x,y),y,1);
        diff(diff(pol(x,y),x,1),y,1);
        diff(pol(x,y),x,2);
        diff(pol(x,y),y,2);
        diff(diff(pol(x,y),x,2),y,1);
        diff(diff(pol(x,y),x,1),y,2);
        diff(pol(x,y),x,3);
        diff(pol(x,y),y,3);

      */


      VectorXd fvec(p);

      VectorXd fxvec(p);
      VectorXd fyvec(p);

      VectorXd fxyvec(p);
      VectorXd fxxvec(p);
      VectorXd fyyvec(p);

      VectorXd fxxyvec(p);
      VectorXd fxyyvec(p);
      VectorXd fxxxvec(p);
      VectorXd fyyyvec(p);


      fvec[0]=1.0;


      fxvec[0]=0.0;
      fyvec[0]=0.0;

      fxxvec[0]=0.0;
      fxyvec[0]=0.0;
      fyyvec[0]=0.0;

      fxxyvec[0]=0.0;
      fxyyvec[0]=0.0;
      fxxxvec[0]=0.0;
      fyyyvec[0]=0.0;


      if(degree>0){
	fvec[1]=(lx-x);
	fxvec[1]=1.0;
	fyvec[1]=0.0;
	fxyvec[1]=0.0;
	fxxvec[1]=0.0;
	fyyvec[1]=0.0;
	fxxyvec[1]=0.0;
	fxyyvec[1]=0.0;
	fxxxvec[1]=0.0;
	fyyyvec[1]=0.0;

	fvec[2]=(ly-y);
	fxvec[2]=0.0;
	fyvec[2]=1.0;
	fxyvec[2]=0.0;
	fxxvec[2]=0.0;
	fyyvec[2]=0.0;
	fxxyvec[2]=0.0;
	fxyyvec[2]=0.0;
	fxxxvec[2]=0.0;
	fyyyvec[2]=0.0;


	if(degree>1){
	  fvec[3]=(lx-x)*(ly-y);
	  fxvec[3]=(ly-y);
	  fyvec[3]=(lx-x);
	  fxyvec[3]=1.0;
	  fxxvec[3]=0.0;
	  fyyvec[3]=0.0;
	  fxxyvec[3]=0.0;
	  fxyyvec[3]=0.0;
	  fxxxvec[3]=0.0;
	  fyyyvec[3]=0.0;

	  fvec[4]=(lx-x)*(lx-x);
	  fxvec[4]=2.0*(lx-x);
	  fyvec[4]=0.0;
	  fxyvec[4]=0.0;
	  fxxvec[4]=2.0;
	  fyyvec[4]=0.0;
	  fxxyvec[4]=0.0;
	  fxyyvec[4]=0.0;
	  fxxxvec[4]=0.0;
	  fyyyvec[4]=0.0;

	  fvec[5]=(ly-y)*(ly-y);
	  fxvec[5]=0.0;
	  fyvec[5]=2.0*(ly-y);
	  fxyvec[5]=0.0;
	  fxxvec[5]=0.0;
	  fyyvec[5]=2.0;
	  fxxyvec[5]=0.0;
	  fxyyvec[5]=0.0;
	  fxxxvec[5]=0.0;
	  fyyyvec[5]=0.0;


	  if(degree>2){
	    fvec[6]  = (lx-x)*(lx-x)*(ly-y);
	    fxvec[6] = 2.0*(lx-x)*(ly-y);
	    fyvec[6] = (lx-x)*(lx-x);
	    fxyvec[6] = 2.0*(lx-x);
	    fxxvec[6] = 2.0*(ly-y);
	    fyyvec[6] = 0.0;
	    fxxyvec[6]=2.0;
	    fxyyvec[6]=0.0;
	    fxxxvec[6]=0.0;
	    fyyyvec[6]=0.0;

	    fvec[7]  = (lx-x)*(ly-y)*(ly-y);
	    fxvec[7] = (ly-y)*(ly-y);
	    fyvec[7] = 2.0*(lx-x)*(ly-y);
	    fxyvec[7] = 2.0*(ly-y);
	    fxxvec[7] = 0.0;
	    fyyvec[7] = 2.0*(lx-x);
	    fxxyvec[7]=0.0;
	    fxyyvec[7]=2.0;
	    fxxxvec[7]=0.0;
	    fyyyvec[7]=0.0;

	    fvec[8]  = (lx-x)*(lx-x)*(lx-x);
	    fxvec[8] = 3.0*(lx-x)*(lx-x);
	    fyvec[8] = 0.0;
	    fxyvec[8] = 0.0;
	    fxxvec[8] = 6.0*(lx-x);
	    fyyvec[8] = 0.0;
	    fxxyvec[8]=0.0;
	    fxyyvec[8]=0.0;
	    fxxxvec[8]=6.0;
	    fyyyvec[8]=0.0;

	    fvec[9]  = (ly-y)*(ly-y)*(ly-y);
	    fxvec[9] = 0.0;
	    fyvec[9] = 3.0*(ly-y)*(ly-y);
	    fxyvec[9] = 0.0;
	    fxxvec[9] = 0.0;
	    fyyvec[9] = 6.0*(ly-y);
	    fxxyvec[9]=0.0;
	    fxyyvec[9]=0.0;
	    fxxxvec[9]=0.0;
	    fyyyvec[9]=6.0;

	  }
	}
      }


      // FIXME, check coefficients to find bug

      Z[i]=1.0/20.0*pde.betahat.transpose()*fvec;
      if(degree>=1){
        Zx[i]=-1.0/10.0*pde.betahat.transpose()*fxvec;
        Zy[i]=-1.0/10.0*pde.betahat.transpose()*fyvec;
        if(degree>=2){
          Zxy[i]=1.0/3.0*pde.betahat.transpose()*fxyvec;
          Zxx[i]=1.0/3.0*pde.betahat.transpose()*fxxvec;
          Zyy[i]=1.0/3.0*pde.betahat.transpose()*fyyvec;
          if(degree>=3){
            Zxxy[i]=-pde.betahat.transpose()*fxxyvec;
            Zxyy[i]=-pde.betahat.transpose()*fxyyvec;
            Zxxx[i]=-pde.betahat.transpose()*fxxxvec;
            Zyyy[i]=-pde.betahat.transpose()*fyyyvec;
          }
        }
      }
    }
  }

  // Prepare Akimas weighting scheme
  VectorXd pdmean(n);
  VectorXd pdsd(n);
  VectorXd weight(n),vweight(n),pweight(n);

  // first part: 5- (2-) dimensional normal density values as weights,
  // estimate componentwise parameters, use product density
  if(akimaweight){
    pdmean[0]=Z.sum()/n;

    pdsd[0]=1.0/(n-1)*((Z.array()-pdmean[0]).array()*(Z.array()-pdmean[0]).array()).sum();


    if(degree>=1){
      pdmean[1]=Zx.sum()/n;
      pdsd[1]=1.0/(n-1)*((Z.array()-pdmean[1]).array()*(Z.array()-pdmean[1]).array()).sum();
      pdmean[2]=Zy.sum()/n;
      pdsd[2]=1.0/(n-1)*((Z.array()-pdmean[2]).array()*(Z.array()-pdmean[2]).array()).sum();
      if(degree>=2){
        pdmean[3]=Zxy.sum()/n;
        pdsd[3]=1.0/(n-1)*((Z.array()-pdmean[3]).array()*(Z.array()-pdmean[3]).array()).sum();
        pdmean[4]=Zxx.sum()/n;
        pdsd[4]=1.0/(n-1)*((Z.array()-pdmean[4]).array()*(Z.array()-pdmean[4]).array()).sum();
        pdmean[5]=Zxy.sum()/n;
        pdsd[5]=1.0/(n-1)*((Z.array()-pdmean[5]).array()*(Z.array()-pdmean[5]).array()).sum();
        if(degree>=3){
          pdmean[6]=Zxxy.sum()/n;
          pdsd[6]=1.0/(n-1)*((Z.array()-pdmean[6]).array()*(Z.array()-pdmean[6]).array()).sum();
          pdmean[7]=Zxyy.sum()/n;
          pdsd[7]=1.0/(n-1)*((Z.array()-pdmean[7]).array()*(Z.array()-pdmean[7]).array()).sum();
          pdmean[8]=Zxxx.sum()/n;
          pdsd[8]=1.0/(n-1)*((Z.array()-pdmean[8]).array()*(Z.array()-pdmean[8]).array()).sum();
          pdmean[9]=Zyyy.sum()/n;
          pdsd[9]=1.0/(n-1)*((Z.array()-pdmean[9]).array()*(Z.array()-pdmean[9]).array()).sum();
        }
      }
    }


    //Rcout << "pdmean: " << pdmean << std::endl;
    //Rcout << "pdsd: " << pdsd << std::endl;
    // this doesnt work, why?
    // weight=dnorm(wrap(Zx),pdmean[1],pdsd[1]);
    pweight=(myDnorm(Zx,pdmean[1],pdsd[1])).array()*
      (myDnorm(Zy,pdmean[2],pdsd[2])).array();
    vweight=(Zx.array()-Lx.array())*(Zx.array()-Lx.array())+
      (Zy.array()-Ly.array())*(Zy.array()-Ly.array());
    if(degree>=2){
      pweight=pweight.array()*(myDnorm(Zxy,pdmean[3],pdsd[3])).array()*
        (myDnorm(Zxx,pdmean[4],pdsd[4])).array()*
        (myDnorm(Zyy,pdmean[5],pdsd[5])).array();
      vweight=vweight.array()+Zxy.array()*Zxy.array()+
        Zxx.array()*Zxx.array()+
        Zyy.array()*Zyy.array();
    }
    double wsp=pweight.sum();
    double wsv=vweight.sum();
    pweight = pweight.array()/wsp;
    // vweight = vweight.array()/wsv;
    // TODO: if vweight != 0
    bool gtZ=true;
    //Rcout << "use volatility weights: " << gtZ << std::endl;
    for(int i=0; i<n; i++)
      if(vweight[i]<10e-16) gtZ=false;
    if(gtZ)
      weight=pweight.array()/vweight.array();
    else
      weight=pweight;
    // weight=pweight;
    //Rcout << "pweight: " << pweight << std::endl;
    //Rcout << "vweight: " << vweight << std::endl;
    //Rcout << "weight: "  << weight << std::endl;
  } else {
    for(int i=0;i<n;i++)
      weight[i]=1.0/n;
  }
  //Rcout << "n=" << n << std::endl;
  //Rcout << "weight: " << weight << std::endl;

  // second part: inverse distance in R^5 (R^2) to parameters of a
  // plane estimated by least squares
  // TODO

  pde.est[0]=weight.transpose()*Z;
  if(degree>=1){
    pde.est[1]=weight.transpose()*Zx;
    pde.est[2]=weight.transpose()*Zy;
    if(degree>=2){
      pde.est[3]=weight.transpose()*Zxy;
      pde.est[4]=weight.transpose()*Zxx;
      pde.est[5]=weight.transpose()*Zxy;
      if(degree>=3){
        pde.est[6]=weight.transpose()*Zxxy;
        pde.est[7]=weight.transpose()*Zxyy;
        pde.est[8]=weight.transpose()*Zxxx;
        pde.est[9]=weight.transpose()*Zyyy;
      }
    }
  }


  // note: pde.betahat is here meaningless !!
  return pde;
}


// [[Rcpp::export(name="nearest.neighbours")]]
List nearestNeighbours(NumericVector x, NumericVector y){

  NN ans=nN(x,y);

  List ret=List::create(_("index")=(ans.ind.array()+1).matrix(), _("dist")=ans.dist);

  return ret;
}

NN nN(NumericVector x, NumericVector y){
  NN ret;

  //Rcout << "x: " << x << std::endl;
  //Rcout << "y: " << y << std::endl;
  int n=x.size();
  if(y.size()!=n)
    Rf_error("sizes of x and y dont match!");

  ret.ind=MatrixXi(n,n).setZero();
  ret.dist=MatrixXd(n,n).setZero();
  // FIXME: exclude case i==j !!!, return matrix should be n x (n-1) 
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      double dij=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
      //Rcout << "dist: " << dij << std::endl;
      // simply record first neighbour
      if(j==0){
        ret.dist(i,j)=dij;
        ret.ind(i,j)=j;
      } else //if(i!=j)
        { // sort in other neighbours, skip i==j
          for(int k=0; k<j; k++){
            if(ret.dist(i,k)>dij){
              // shift right to make room for insert
              for(int l=j;l>k;l--){
                ret.dist(i,l)=ret.dist(i,l-1);
                ret.ind(i,l)=ret.ind(i,l-1);
              }
              // insert
              //Rcout << "point " << i << ", insert " << j << " at " << k <<std::endl;
              ret.dist(i,k)=dij;
              ret.ind(i,k)=j;
              break;
            } else {
              //append
              //Rcout << "point " << i << ", append " << j <<std::endl;
              ret.dist(i,j)=dij;
              ret.ind(i,j)=j;
            }
          }
        }
    }
  }
  return ret;
}

NN nN(VectorXd x, VectorXd y){
  return nN(wrap(x),wrap(y));
}

NN extendNN(NN nn, NumericVector X, NumericVector Y,
            NumericVector x, NumericVector y){
  NN ret;

  //Rcout << "x: " << x << std::endl;
  //Rcout << "y: " << y << std::endl;
  int n=x.size();
  int N=X.size();
  if(y.size()!=n)
    Rf_error("sizes of x and y dont match!");
  if(Y.size()!=N)
    Rf_error("sizes of X and Y dont match!");
  //Rcout << nn.ind.rows() << ", " << nn.ind.cols() << ", " << N << std::endl;
  if((nn.ind.rows()!=N) | (nn.ind.cols()!=N))
    Rf_error("sizes of nn and X and Y dont match!");

  ret.ind=MatrixXi(n+N,n+N).setZero();
  ret.dist=MatrixXd(n+N,n+N).setZero();
  VectorXd xtmp = VectorXd(n+N);
  xtmp << Rcpp::as<Eigen::Map<Eigen::VectorXd> >(X),
    Rcpp::as<Eigen::Map<Eigen::VectorXd> >(x);
  VectorXd ytmp  = VectorXd(n+N);
  ytmp << Rcpp::as<Eigen::Map<Eigen::VectorXd> >(Y),
    Rcpp::as<Eigen::Map<Eigen::VectorXd> >(y);
  ret.ind.block(0,0,N,N)=nn.ind;
  ret.dist.block(0,0,N,N)=nn.dist;

  for(int i=0; i<n+N; i++){
    for(int j=0; j<n+N; j++){
      if(((i<N) & (j>=N)) | (i>+N)){
        double dij=sqrt((xtmp[i]-xtmp[j])*(xtmp[i]-xtmp[j])+(ytmp[i]-ytmp[j])*(ytmp[i]-ytmp[j]));
        //Rcout << "dist: " << dij << std::endl;
        // simply record first  neighbour
        // sort in other neighbours
        for(int k=0; k<j; k++){
          if(ret.dist(i,k)>dij){
            // shift right
            for(int l=j;l>k;l--){
              ret.dist(i,l)=ret.dist(i,l-1);
              ret.ind(i,l)=ret.ind(i,l-1);
            }
            // insert
            //Rcout << "point " << i << ", insert " << j << " at " << k <<std::endl;
            ret.dist(i,k)=dij;
            ret.ind(i,k)=j;
            break;
          }
          //append
          //Rcout << "point " << i << ", append " << j <<std::endl;
          ret.dist(i,j)=dij;
          ret.ind(i,j)=j;
        }
      }
    }
  }
  return ret;
}


VectorXd myDnorm(VectorXd x, double mu, double sd){

  int n=x.size();
  VectorXd ret(n);

  for(int i=0; i<n; i++){
    ret[i]=1.0/sqrt(2.0*M_PI)/sd*exp(-1.0/2.0/sd*(x[i]-mu)*(x[i]-mu));
  }

  return ret;
}
