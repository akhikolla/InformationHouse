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


List MM_ML2Mat(VectorXd & Y, MatrixXd & Xinit , MatrixXd & U, VectorXd & D , VectorXd & Init , int MaxIter, double CritVar, double CritLogLik){
	int NbVar = 2; 
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

	VectorXd LogLik(MaxIter);
        LogLik.setZero();
        double crit = CritVar+1;
  	double critLL = CritLogLik + 1;
        int iteration = 0;
        int it = 0;

	while((((crit>CritVar)||(critLL>CritLogLik))&&(iteration<MaxIter))||(it==0)){
		Var = SigmaNew(0) * D + SigmaNew(1) * I;
		Var_inv = Var.cwiseInverse();
		Varx = Var_inv.asDiagonal() * X ;

		double logdetVar , logdetxVarx , detxVarx , traceP , tracePD;
		logdetVar = Var.array().log().sum();

		xVarx = X.transpose() * Varx;
		sym_inverse(xVarx , xVarx_inv , logdetxVarx , detxVarx , 0);
		
		Beta = xVarx_inv.selfadjointView<Lower>() * (Varx.transpose() * Y);

        	VectorXd InvFixed = Var_inv.asDiagonal() * (Y - X*Beta) ;
	
                LogLik(iteration) = -(NbObs * log(2*M_PI) + logdetVar + (Y-X*Beta).transpose() * InvFixed)/2;

		if (iteration > 0){
			if (LogLik(iteration) < LogLik(iteration-1)){
				SigmaNew = SigmaNewBis;
				Var = SigmaNew(0) * D + SigmaNew(1) * I;
				Var_inv = Var.cwiseInverse();
				Varx = Var_inv.asDiagonal() * X ;

				double logdetVar , logdetxVarx , detxVarx;
				logdetVar = Var.array().log().sum();
	
				xVarx = X.transpose() * Varx;
				sym_inverse(xVarx , xVarx_inv , logdetxVarx , detxVarx , 0);
		
				Beta = xVarx_inv.selfadjointView<Lower>() * (Varx.transpose() * Y);

        			VectorXd InvFixed = Var_inv.asDiagonal() * (Y - X*Beta) ;
	
                		LogLik(iteration) = -(NbObs * log(2*M_PI) + logdetVar + (Y-X*Beta).transpose() * InvFixed)/2;
			}
		}

                VectorXd Toto(NbVar);
		Toto.setZero();

		traceP = Var_inv.sum();
		tracePD = D.dot(Var_inv);
		
                Toto(0) = sqrt(InvFixed.dot(D.cwiseProduct(InvFixed))/tracePD);
                Toto(1) = sqrt(InvFixed.dot(InvFixed) / traceP);

		SigmaOld = SigmaNew;
		SigmaNew = Toto.cwiseProduct(SigmaOld);		
       
		it ++;
		if (it == 2 ){
			SigmaNewBis = SigmaNew;

                	VectorXd r = SigmaOld - SigmaAnc;


                        VectorXd v = SigmaNew - SigmaOld - r;

                        double Alpha = - sqrt(r.squaredNorm()/v.squaredNorm());

                        if (Alpha > -1){
                        	Alpha = -1;
                        }			

                        SigmaNew = SigmaAnc - 2*Alpha*r + pow(Alpha,2)*v;

                   	for (int i=0 ; i<NbVar ; i++){
				if (SigmaNew[i]<=0) SigmaNew = v + r + SigmaOld;
                        }
                          

                        it = 0;
		}

		VectorXd DiffAbs = SigmaNew-SigmaOld;

   		crit = DiffAbs.lpNorm<Infinity>();
    		if (iteration >0) critLL = LogLik(iteration) - LogLik(iteration-1); 
                
                SigmaAnc = SigmaOld;


                iteration ++;
	}
        
	double LogLikRestrain = LogLik(iteration-1);

	return List::create(Rcpp::Named("Beta")=Beta,Rcpp::Named("Sigma2")=SigmaNew,Rcpp::Named("VarBeta") = xVarx_inv,Rcpp::Named("LogLik (ML)")=LogLikRestrain,Rcpp::Named("NbIt")=iteration,Rcpp::Named("Method")="ML");
	
}                         
                          


