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


List MM_RemlRcppHen(VectorXd & Y, MatrixXd & X, MatrixXd & Z , List & GList , List & GinvList , MatrixXd & Rinv , VectorXd logdetV , VectorXd & Init , int MaxIter, double CritVar , double CritLogLik){
	int NbVar(GList.size());
 
        int NbObs(Y.size());
	int NbCof(X.cols());
	int NbAlea(Z.cols());

	double logdetG;


	MatrixXd Ginv(NbAlea,NbAlea);
	Ginv.setZero();
	MatrixXd GinvCompl(NbAlea+NbCof,NbAlea+NbCof);
	GinvCompl.setZero();


	MatrixXd xRx = X.transpose() * Rinv * X;
	MatrixXd xRx_inv(xRx);
	double detR,detRinv;
	sym_inverse(xRx,xRx_inv,detR,detRinv,0);
	MatrixXd R(as<Map<MatrixXd> >(GList[NbVar-1]));

	VectorXd DimG(NbVar);
	DimG.setZero();

	MatrixXd Hen(NbAlea+NbCof,NbAlea+NbCof);
	Hen.setZero();
	MatrixXd MME(NbAlea+NbCof,NbAlea+NbCof);
	MME.setZero();
	MatrixXd MME_inv(MME);
	MME_inv.setZero();

	MatrixXd L(NbAlea+NbCof,NbAlea+NbCof);
	L.setZero();
	MatrixXd L_inv(NbAlea+NbCof,NbAlea+NbCof);
	L_inv.setZero();

	MatrixXd CHen(NbAlea,NbAlea);
	MatrixXd C_inv(NbAlea,NbAlea);

	VectorXd VecMME(NbAlea+NbCof);

	MatrixXd RinvX = Rinv.selfadjointView<Lower>() * X;
	MatrixXd RinvZ = Rinv.selfadjointView<Lower>() * Z;
	MatrixXd xRinvZ = X.transpose() * RinvZ;
	MatrixXd zSz = Z.transpose()*RinvZ - xRinvZ.transpose() * xRx_inv.selfadjointView<Lower>() * xRinvZ;


	VectorXd VecMMEz = RinvZ.transpose() * Y;
	VectorXd VecMMEx = RinvX.transpose() * Y;
	
	for (int i=0 ; i < NbCof ; i++){
		VecMME[i] = VecMMEx[i];
	}
	for (int i = NbCof ; i < NbAlea+NbCof ; i++){
		VecMME[i] = VecMMEz[i-NbCof];
	}



	Hen.block(0,0,NbCof,NbCof) = X.transpose() * RinvX;
	Hen.block(0,NbCof,NbCof,NbAlea) = X.transpose() * RinvZ;
	Hen.block(NbCof,0,NbAlea,NbCof) = Z.transpose() * RinvX;
	Hen.block(NbCof,NbCof,NbAlea,NbAlea) = Z.transpose() * RinvZ;
	MatrixXd GinvC(NbAlea,NbAlea);
	VectorXd Py(NbObs);
	VectorXd ZPy(NbObs);

	VectorXd Beta(NbCof);
	VectorXd Uhat(NbAlea);
	VectorXd Pred(NbObs);
	VectorXd Residual(NbObs);
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


		int Tmp = 0;

		for (int i=0; i<(NbVar-1) ; i++){
			MatrixXd MatG_inv(as<Map<MatrixXd> >(GinvList[i]));
			DimG[i+1]=MatG_inv.cols();
			Tmp += DimG[i];
			Ginv.block(Tmp,Tmp,DimG[i+1],DimG[i+1]) = MatG_inv * SigmaNew[NbVar-1]/SigmaNew[i];
		}


		GinvCompl.block(NbCof,NbCof,NbAlea,NbAlea) = Ginv;		

		MME = Hen + GinvCompl;

		double logdetC,logdetCinv;
		sym_inverse(MME,MME_inv,logdetC,logdetCinv,0);
		Pred = MME_inv.selfadjointView<Lower>() * VecMME;


		for (int i = 0; i < NbCof ; i++){
			Beta[i] = Pred[i];
		}
		for (int i = NbCof ; i < NbAlea+NbCof ; i++){
			Uhat[i-NbCof] = Pred[i];
		}



		Residual = Y - X*Beta - Z*Uhat;

		Py = Rinv.selfadjointView<Lower>() * Residual;
		ZPy = Ginv.selfadjointView<Lower>() * Uhat;

		logdetG = 0;
		for (int i=0 ; i < NbVar-1 ; i++){
			logdetG += logdetV[i] + DimG[i+1]*log(SigmaNew[i]/SigmaNew[NbVar-1]);
		}	

		LogLik(iteration) = -(logdetC + logdetV[NbVar-1] + logdetG + (NbObs - NbCof) *log(SigmaNew[NbVar-1]) + (Y.transpose()*Py).sum() / SigmaNew[NbVar-1])/2;

		if (iteration > 0){
			if (LogLik(iteration) < LogLik(iteration-1)){

				SigmaNew = SigmaNewBis;
				int Tmp = 0;

				for (int i=0; i<(NbVar-1) ; i++){
					MatrixXd MatG_inv(as<Map<MatrixXd> >(GinvList[i]));
					DimG[i+1]=MatG_inv.cols();
					Tmp += DimG[i];
					Ginv.block(Tmp,Tmp,DimG[i+1],DimG[i+1]) = MatG_inv * SigmaNew[NbVar-1]/SigmaNew[i];
				}


				GinvCompl.block(NbCof,NbCof,NbAlea,NbAlea) = Ginv;		
		
				MME = Hen + GinvCompl;

				double logdetC,logdetCinv;
				sym_inverse(MME,MME_inv,logdetC,logdetCinv,0);
				Pred = MME_inv.selfadjointView<Lower>() * VecMME;


				for (int i = 0; i < NbCof ; i++){
					Beta[i] = Pred[i];
				}
				for (int i = NbCof ; i < NbAlea+NbCof ; i++){
					Uhat[i-NbCof] = Pred[i];
				}
		
				Residual = Y - X*Beta - Z*Uhat;

				Py = Rinv.selfadjointView<Lower>() * Residual;
				ZPy = Ginv.selfadjointView<Lower>() * Uhat;


				SigmaOld = SigmaNew;
				logdetG = 0;
				for (int i=0 ; i < NbVar-1 ; i++){
					logdetG += logdetV[i] + DimG[i+1]*log(SigmaNew[i]/SigmaNew[NbVar-1]);
				}	

                		LogLik(iteration) = -(logdetC + logdetV[NbVar-1] + logdetG + (NbObs - NbCof) *log(SigmaNew[NbVar-1]) + (Y.transpose()*Py).sum() / SigmaNew[NbVar-1])/2;

				


			}
		}		

		SigmaOld = SigmaNew;


		C_inv = MME_inv.block(NbCof,NbCof,NbAlea,NbAlea);
		VectorXd yPy = Py.transpose() * R * Py;
		SigmaNew[NbVar-1] = SigmaOld[NbVar-1] * sqrt(yPy.sum()/(SigmaOld[NbVar-1]*(NbObs - NbCof - trace_of_product(zSz,C_inv))));
		GinvC = Ginv.selfadjointView<Lower>() * C_inv;
	
		int Stop = 0;
		int Start = 0;

		for (int i = 0 ; i < NbVar-1 ; i++){
			Stop += DimG[i+1];
			Start += DimG[i];
			VectorXd ZPySig(Stop-Start);
			ZPySig.setZero();

			for (int j = Start ; j < Stop ; j++){
				ZPySig[j-Start] = ZPy[j];
			}




			MatrixXd MatVar(as<Map<MatrixXd> >(GList[i]));
			


			VectorXd Vec = MatVar.selfadjointView<Lower>() * ZPySig;



			SigmaNew[i] = SigmaOld[i] * sqrt((ZPySig.transpose() * Vec).sum()/ ((SigmaOld[NbVar-1]*SigmaOld[NbVar-1] / SigmaOld[i])*(DimG[i+1]-(GinvC.block(Start,Start,DimG[i+1],DimG[i+1]).trace())  ) ) );
		}	


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

	Hen.block(0,0,NbAlea,NbAlea) = Z.transpose() * RinvZ + Ginv;
	Hen.block(0,NbAlea,NbAlea,NbCof) = Z.transpose() * RinvX;
	Hen.block(NbAlea,0,NbCof,NbAlea) = X.transpose() * RinvZ;
	Hen.block(NbAlea,NbAlea,NbCof,NbCof) = X.transpose() * RinvX;

	MatrixXd Hen_inv(Hen);
	double logdetHen,logdetHeninv;
	sym_inverse(Hen,Hen_inv,logdetHen,logdetHeninv,0);

	MatrixXd VarBeta = SigmaOld[NbVar-1]*Hen_inv.bottomRightCorner(NbCof,NbCof);

	double LogLikRestrain = LogLik(iteration-1);

	return List::create(Rcpp::Named("Beta")=Beta,Rcpp::Named("Sigma2")=SigmaNew,Rcpp::Named("VarBeta") = VarBeta,Rcpp::Named("LogLik (Reml)")=LogLikRestrain,Rcpp::Named("NbIt")=iteration,Rcpp::Named("Method") = "Reml");
}                         
                          


