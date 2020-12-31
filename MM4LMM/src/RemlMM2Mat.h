// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
// #define ANY(_X_) (std::any_of(_X_.begin(), _X_.end(), [](bool x) {return x;})) 




#include <Rcpp.h>
#include "Inversion.h"
using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;




List MM_Reml2Mat(VectorXd & Y, MatrixXd & Xinit , MatrixXd & U, VectorXd & D , VectorXd & Init , int MaxIter, double CritVar , double CritLogLik){
        int NbObs(Y.size());
	MatrixXd X = U*Xinit;
	int NbCof(Xinit.cols());
	VectorXd Var(NbObs);
        Var.setZero();
	VectorXd Var_inv(NbObs);
	MatrixXd xVarx(NbCof,NbCof);
	MatrixXd xVarx_inv(NbCof,NbCof);
  	MatrixXd Varx(NbObs,NbCof);

	VectorXd I = VectorXd::Ones(NbObs);
	VectorXd Beta(X.cols());
	VectorXd SigmaAnc = Init;
        VectorXd SigmaOld = Init;
        VectorXd SigmaNew = Init;
        VectorXd SigmaNewBis = Init;
	double deltaAnc = SigmaNew(0)/SigmaNew(1);
	double deltaNew = deltaAnc;
	double deltaOld = deltaAnc;
	double deltaNewBis = deltaAnc;

        VectorXd InvFixed(Y);

	VectorXd LogLik(MaxIter+20);
        LogLik.setZero();
        double crit = CritVar+1;
  	double critLL = CritLogLik + 1;
        int iteration = 0;
        int it = 0;

	double logdetVar , logdetxVarx , detxVarx;
	double A,B,C,AB,BC,Delta;

	while((((crit>CritVar)||(critLL>CritLogLik))&&(iteration<MaxIter))||(it==0)){
		Var = deltaNew * D + I;
		Var_inv = Var.cwiseInverse();

		Varx = Var_inv.asDiagonal() * X ;

		logdetVar = Var.array().log().sum();

		xVarx = X.transpose() * Varx;
		sym_inverse(xVarx , xVarx_inv , logdetxVarx , detxVarx , 0);


        	InvFixed = Var_inv.asDiagonal() * Y - Varx * ( xVarx_inv * ( Varx.transpose() * Y ) ) ;

        	VectorXd Pdiag = Var_inv - (Varx * xVarx_inv * Varx.transpose()).diagonal();

		A = (Y.transpose() * InvFixed).sum();
		B = (InvFixed.transpose()*D.asDiagonal() * InvFixed).sum();


    MatrixXd T = xVarx_inv.selfadjointView<Lower>() * Varx.transpose() * D.asDiagonal() * Varx;

		C = D.dot(Var_inv) - T.trace();
		BC = B*C*pow(deltaNew,2);
		AB = A - B * deltaNew;
		Delta = pow(BC,2) + 4 * (NbObs - NbCof) * AB * BC;

		deltaOld = deltaNew;
		SigmaOld = SigmaNew;
		deltaNew = (- BC + sqrt( Delta )) / (2 * AB * C);
		SigmaNew(1) = A/(NbObs-NbCof);		
		SigmaNew(0) = deltaOld * SigmaNew(1);

		
		LogLik(iteration) = -((NbObs-NbCof)*log(SigmaNew(1)) + logdetxVarx + logdetVar + A / SigmaNew(1))/2;
   
		if (iteration > 0){
			if (LogLik(iteration) < LogLik(iteration - 1)){

				deltaNew = deltaNewBis;
				Var = deltaNew * D + I;
				Var_inv = Var.cwiseInverse();
				Varx = Var_inv.asDiagonal() * X ;
				logdetVar = Var.array().log().sum();

				xVarx = X.transpose() * Varx;
				sym_inverse(xVarx , xVarx_inv , logdetxVarx , detxVarx , 0);


        			InvFixed = Var_inv.asDiagonal() * Y - Varx * ( xVarx_inv * ( Varx.transpose() * Y ) ) ;

				A = (Y.transpose() * InvFixed).sum();
				B = (InvFixed.transpose()*D.asDiagonal() * InvFixed).sum();

				MatrixXd T = xVarx_inv.selfadjointView<Lower>() * Varx.transpose() * D.asDiagonal() * Varx;
				C = D.dot(Var_inv) - T.trace();
				BC = B*C*pow(deltaNew,2);
				AB = A - B * deltaNew;
				Delta = pow(BC,2) + 4 * (NbObs - NbCof) * AB * BC;

				deltaOld = deltaNew;
				SigmaOld = SigmaNew;
				deltaNew = (- BC + sqrt( Delta )) / (2 * AB * C);
				SigmaNew(1) = A/(NbObs-NbCof);		
				SigmaNew(0) = deltaOld * SigmaNew(1);

				LogLik(iteration) = -((NbObs-NbCof)*log(SigmaNew(1)) + logdetxVarx + logdetVar + (Y.transpose() * InvFixed).sum() / SigmaNew(1))/2;
                	}
		}
			 
       		it ++;
		if (it == 2 ){
			deltaNewBis = deltaNew;
                	double r = deltaOld - deltaAnc;
                        double v = deltaNew - deltaOld - r;

                        double Alpha = - sqrt(pow(r,2)/pow(v,2));

                        if (Alpha > -1){
                        	Alpha = -1;
                        }			

                        deltaNew = deltaAnc - 2*Alpha*r + pow(Alpha,2)*v;

			if (deltaNew<=0) deltaNew = v + r + deltaOld;

                        it = 0;
		}

		VectorXd DiffAbs = SigmaNew-SigmaOld;
   		crit = DiffAbs.lpNorm<Infinity>();
    		if (iteration >0) critLL = LogLik(iteration) - LogLik(iteration-1); 
		SigmaAnc = SigmaOld;

		iteration ++;
		
	}
	Var = deltaNew * D + I;
	Var_inv = Var.cwiseInverse();
	Varx = Var_inv.asDiagonal() * X ;
	
	logdetVar = Var.array().log().sum();

	xVarx = X.transpose() * Varx;
	sym_inverse(xVarx , xVarx_inv , logdetxVarx , detxVarx , 0);


        InvFixed = Var_inv.asDiagonal() * Y - Varx * ( xVarx_inv * ( Varx.transpose() * Y ) ) ;
	

	SigmaNew(1) = (Y.transpose() * InvFixed).sum() / (NbObs - NbCof);
	SigmaNew(0) = deltaNew * SigmaNew(1);	

	LogLik(iteration) = -((NbObs-NbCof)*log(SigmaNew(1)) + logdetxVarx + logdetVar + (Y.transpose() * InvFixed).sum() / SigmaNew(1))/2;	
	double LogLikRestrain = LogLik(iteration);
        
	Beta = xVarx_inv.selfadjointView<Lower>() * (Varx.transpose() * Y);

	return List::create(Rcpp::Named("Beta")=Beta,Rcpp::Named("Sigma2")=SigmaNew,Rcpp::Named("VarBeta") = xVarx_inv * SigmaNew(1),Rcpp::Named("LogLik (Reml)")=LogLikRestrain,Rcpp::Named("NbIt")=iteration,Rcpp::Named("Method")="Reml");
        
	
}                         
                          


