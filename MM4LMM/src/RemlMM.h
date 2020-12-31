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


List MM_RemlRcpp(VectorXd & Y, MatrixXd & X, List & VarList , VectorXd & Init , int MaxIter, double CritVar, double CritLogLik){
	int NbVar(VarList.size()); 
        int NbObs(Y.size());
	int NbCof(X.cols());
	MatrixXd Var(NbObs,NbObs);
        Var.setZero();
	MatrixXd Var_inv(NbObs,NbObs);
	MatrixXd xVarx(NbCof,NbCof);
	MatrixXd xVarx_inv(NbCof,NbCof);
  	MatrixXd Varx(NbObs,NbCof);
	MatrixXd P(NbObs,NbObs);

	VectorXd InvFixed(NbObs);
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

		Var.setZero();
		for (int i=0; i<NbVar ; i++){
			MatrixXd MatVar(as<Map<MatrixXd> >(VarList[i]));
			Var += MatVar * SigmaNew(i);
		}

		double logdetVar , detVar , logdetxVarx , detxVarx;

		sym_inverse(Var , Var_inv , logdetVar , detVar , 0);
		Varx = Var_inv * X;
		xVarx = X.transpose() * Varx;
		sym_inverse(xVarx , xVarx_inv , logdetxVarx , detxVarx , 0);
	
		P = Var_inv - Varx * xVarx_inv.selfadjointView<Lower>() * Varx.transpose() ;
        	InvFixed = P.selfadjointView<Lower>() * Y;


		LogLik(iteration) = -(logdetxVarx + logdetVar + Y.transpose() * InvFixed)/2;



		if (iteration > 0){
			if (LogLik(iteration)<LogLik(iteration-1)){

				SigmaNew = SigmaNewBis;

				Var.setZero();
				for (int i=0; i<NbVar ; i++){
					MatrixXd MatVar(as<Map<MatrixXd> >(VarList[i]));
					Var += MatVar * SigmaNew(i);
				}

				double logdetVar , detVar , logdetxVarx , detxVarx;

				sym_inverse(Var , Var_inv , logdetVar , detVar , 0);
				Varx = Var_inv * X;
				xVarx = X.transpose() * Varx;
				sym_inverse(xVarx , xVarx_inv , logdetxVarx , detxVarx , 0);
	
				P = Var_inv - Varx * xVarx_inv.selfadjointView<Lower>() * Varx.transpose() ;
        			InvFixed = P.selfadjointView<Lower>() * Y;
	


				LogLik(iteration) = -(logdetxVarx + logdetVar + Y.transpose() * InvFixed)/2;
			}
		}

                VectorXd Toto(NbVar);
		Toto.setZero();



		for (int i=0 ; i<NbVar ; i++){
                	MatrixXd MatVar(as<Map<MatrixXd> >(VarList[i]));
			VectorXd Vec = MatVar.selfadjointView<Lower>() * InvFixed;
			
                        Toto(i) = sqrt((InvFixed.transpose() * Vec).sum()/trace_of_product(P,MatVar));
                }

		SigmaOld = SigmaNew;

		SigmaNew = Toto.cwiseProduct(SigmaOld);		

		VectorXd DiffAbs = SigmaNew-SigmaOld;

   		crit = DiffAbs.lpNorm<Infinity>();
    		if (iteration >0) critLL = LogLik(iteration) - LogLik(iteration-1);
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
                          
                          
			DiffAbs = SigmaNew-SigmaOld;

                	crit = DiffAbs.lpNorm<Infinity>();  

                        it = 0;
		}
                
                SigmaAnc = SigmaOld;


                iteration ++;
	}

	double LogLikRestrain = LogLik(iteration-1);
	Beta = xVarx_inv.selfadjointView<Lower>() * (Varx.transpose() * Y);

	return List::create(Rcpp::Named("Beta")=Beta,Rcpp::Named("Sigma2")=SigmaNew,Rcpp::Named("VarBeta") = xVarx_inv,Rcpp::Named("LogLik (Reml)")=LogLikRestrain,Rcpp::Named("NbIt")=iteration , Rcpp::Named("Method")="Reml");
	
}                         
                          


