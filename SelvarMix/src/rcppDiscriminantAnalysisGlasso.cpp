#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>

using namespace std;
using namespace Rcpp;
using namespace arma;


#include "discriminant.h"
#include "Function.h"


discriminant::discriminant(NumericMatrix X_,  IntegerVector labels_, const int nbClust_, double l, double r){
    nbClust = nbClust_;
    Xd = as<mat>(X_);
    n = Xd.n_rows;
    p = Xd.n_cols;
    
    labels = as< vector<int> >(labels_);
    idx.resize(nbClust);
    for(int k = 0; k < nbClust; k++)
        idx[k].clear();
    
    for(int i = 0; i < (int)Xd.n_rows; i++)
    {
        idx[labels[i]-1].push_back(i);
    }
    
    prop = zeros<rowvec>(nbClust);
    for(int k = 0; k < nbClust; k++)
        prop(k) = (double)idx[k].size()/(double)n;
    
    lambda = l;
    rho = r;
};

void discriminant::Initialization(void)
{
    Environment glasso("package:glasso");
    Function Rglasso = glasso["glasso"];
    
    Mu = zeros<mat>(p, nbClust);
    CovarianceMatrix = cube(p, p, nbClust);
    EmpiricalCovariance = cube(p, p, nbClust);
    PrecisionMatrix = cube(p, p, nbClust);
    
    for(int k = 0; k < nbClust;  k++)
    {
        uvec aux_idx = conv_to<uvec>::from(idx[k]);
        mat X_aux = Xd.rows(aux_idx);
        Mu.col(k) = trans(mean(X_aux));
        EmpiricalCovariance.slice(k) = cov(X_aux);
        List Glasso = Rglasso(Named("s") = EmpiricalCovariance.slice(k),
                              Named("rho") = (2*rho/(double)idx[k].size()),
                              Named("penalize.diagonal") = 0,
                              Named("thr") = 0.001,
                              Named("maxit") = 1000);
        CovarianceMatrix.slice(k) = as<mat>(Glasso[0]);
        PrecisionMatrix.slice(k) =   as<mat>(Glasso[1]);
    }
    
    
    //cout<< " .... Initialization Done .... " <<endl;
    
};




double discriminant::PenLogLik(void)
{
    double SLogDet = 0., sign, val, penloglik=0;
    for(int k = 0; k < nbClust; k++)
    {
        log_det(val, sign, PrecisionMatrix.slice(k));
        SLogDet = -log(sign*exp(val));
        uvec aux_idx = conv_to<uvec>::from(idx[k]);
        mat X_aux = Xd.rows(aux_idx);
        for(int i = 0; i < (int) idx[k].size(); i++)
            penloglik += log(prop(k)) + ldcppmvt(trans(X_aux.row(i)), Mu.col(k), PrecisionMatrix.slice(k), SLogDet);
    };
    
    penloglik -= (lambda * accu(abs(Mu))) - (rho * accu(abs(PrecisionMatrix)));
    return penloglik;
};




void discriminant::UpdateMeans(void)
{
    long double T, Tabs;
    mat MusNew = zeros<mat>(p, nbClust), W;
    colvec MuPrec, MuNew;
    
    
    for(int k = 0; k < nbClust; k++){
        uvec aux_idx = conv_to<uvec>::from(idx[k]);
        mat X_aux = Xd.rows(aux_idx);
        W = PrecisionMatrix.slice(k);
        MuPrec = Mu.col(k);
        MuNew  = zeros<colvec>(p);
        
        for(int j = 0; j < p; j++){
            Tabs = 0.;
            for(int i = 0; i < (int)idx[k].size(); i++)
                Tabs +=  dot(X_aux.row(i) - trans(MuPrec), trans(W.col(j))) + (MuPrec(j)*W(j,j));
            
            
            Tabs = abs(Tabs);
            
            if(Tabs <= lambda) MuNew(j) = 0.;
            else{
                
                //le second membre de l'equation (15)  papier Zhou, Pan , Shen 2009 EJS à k et j fixés
                T = 0.;
                for(int i = 0; i < (int) idx[k].size(); i++)
                    T +=  dot(X_aux.row(i),W.row(j));
                
                T -= ((double)idx[k].size()) * (dot(W.row(j),trans(MuPrec)) - MuPrec(j) * W(j,j));
                if(T < 0.) MuNew(j) = (T + lambda)/(((double)idx[k].size()) * W(j,j));
                else
                    MuNew(j) = (T - lambda)/(((double)idx[k].size())* W(j,j));
            }
        }
        MusNew.col(k) = MuNew;
    }
    
    Mu = MusNew;
};

void discriminant::GetEmpiricalCovariance(void){
    mat S = zeros<mat>(p,p);
    for(int k = 0; k < nbClust; k++)
    {
        uvec aux_idx = conv_to<uvec>::from(idx[k]);
        mat X_aux = Xd.rows(aux_idx);
        S.zeros(p,p);
        for(int i = 0; i < (int)idx[k].size(); i++)
            S = S + ((trans(X_aux.row(i))- Mu.col(k)) * trans(trans(X_aux.row(i))- Mu.col(k)));
        
        EmpiricalCovariance.slice(k) = S/(double)idx[k].size();
    };
}


void discriminant::UpdateCovarianceMatrices(void){
    Environment glasso("package:glasso");
    Function Rglasso = glasso["glasso"];
    for(int k = 0; k < nbClust; k++)
    {
        List Glasso = Rglasso(Named("s") = EmpiricalCovariance.slice(k),
                              Named("rho") = (2*rho/(double)idx[k].size()),
                              Named("penalize.diagonal") = 0,
                              Named("thr") = 0.001,
                              Named("maxit") = 1000);
        CovarianceMatrix.slice(k) = as<mat>(Glasso[0]);
        PrecisionMatrix.slice(k) =   as<mat>(Glasso[1]);
    };
}

rowvec discriminant::VarRole(void){
    rowvec  MuSum = zeros<rowvec>(p);
    rowvec  alive = ones<rowvec>(p);
    
    MuSum = trans(sum(abs(Mu), 1));
    for(int j = 0; j < p; ++j)
        if(MuSum(j) == 0.)
            alive(j) = 0;
    
    return alive;
}

int discriminant::GetNbClust(void)
{
    return(nbClust);
}



//[[Rcpp::export]]
IntegerVector rcppDiscriminantAnalysisGlasso(NumericMatrix X_,  IntegerVector labels_, const int nbClust, double l, double r){
    /*******************************************************************/
    /************ une analyse discriminante pénalisée ******************/
    /*******************************************************************/
    //cout << " je vais commencer " << endl;
    discriminant MyParameter(X_, labels_, nbClust,l, r);
    //cout << " .... K = .... " << MyParameter.GetNbClust() <<" .... lambda = .... " << l << " ....rho = .... " << r<<endl;
    MyParameter.Initialization();
    
    
    double PenLogLik_1  =  MyParameter.PenLogLik(), PenLogLik_0 = 0.0;
    int itr = 0;
    while((abs(PenLogLik_1 - PenLogLik_0) > 0.01) && (itr < 10))
    {
        PenLogLik_0 = PenLogLik_1;
        
        //                    Maximization :
        
        //Calcul des moyennes ....mise à jour des moyennes
        MyParameter.UpdateMeans();
        //Calcul des  matrices de covariances empiriques
        MyParameter.GetEmpiricalCovariance();
        //mise à jour des matrices de convariances et de précision
        MyParameter.UpdateCovarianceMatrices();
        //Calcul de la log-vraisemblance pénalisée
        PenLogLik_1 = MyParameter.PenLogLik();
        itr++;
    };
    
    //cout <<" ... Done ... K = ... " << MyParameter.GetNbClust() <<" ... lambda = ... " << l << " ... rho = ... " << r <<" ... PenLogLik = ... " << PenLogLik_1 << endl;
    return wrap(MyParameter.VarRole());
    
    
};
