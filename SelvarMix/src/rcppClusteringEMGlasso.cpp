#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>


using namespace std;
using namespace Rcpp;
using namespace arma;

#include "Mixture.h"
#include "Function.h"

Mixture::Mixture(List InputList, double l, double r)
{
    Xd = as<mat>(InputList[0]);
    prop = as<rowvec>(InputList[1]);
    n = Xd.n_rows;
    p = Xd.n_cols;
    nbClust  = prop.size();
    
    Mu = as<mat>(InputList[2]);
    
    NumericVector vecArray(as<NumericVector>(InputList[3]));
    CovarianceMatrix = cube(vecArray.begin(), p, p, nbClust, false);
    EmpiricalCovariance = cube(vecArray.begin(), p, p, nbClust, false);// just an initialization
    
    NumericVector vecArray1(as<NumericVector>(InputList[4]));
    PrecisionMatrix = cube(vecArray1.begin(), p, p, nbClust, false);
    
    ProbCond = as<mat>(InputList[5]);
    lambda = l;
    rho = r;
};


Mixture::Mixture(List InputList)
{
    Xd = as<mat>(InputList[0]);
    prop = as<rowvec>(InputList[1]);
    n = Xd.n_rows;
    p = Xd.n_cols;
    nbClust  = prop.size();
    
    Mu = as<mat>(InputList[2]);
    
    NumericVector vecArray(as<NumericVector>(InputList[3]));
    CovarianceMatrix = cube(vecArray.begin(), p, p, nbClust, false);
    EmpiricalCovariance = cube(vecArray.begin(), p, p, nbClust, false);// just an initialization
    
    NumericVector vecArray1(as<NumericVector>(InputList[4]));
    PrecisionMatrix = cube(vecArray1.begin(), p, p, nbClust, false);
    
    ProbCond = as<mat>(InputList[5]);
    lambda = 0.;
    rho = 0.;
};


double Mixture::PenLogLik(void)
{
    mat lD = zeros<mat>(n, nbClust), D = zeros<mat>(n, nbClust), SInv = zeros<mat>(p,p), L = zeros<mat>(p,p),  T = zeros<mat>(n, nbClust);
    double SLogDet = 0., sign, val;
    colvec D_rows_sums = zeros<colvec>(D.n_rows);
    
    for(int k = 0; k < nbClust; k++)
    {
        
        log_det(val, sign, PrecisionMatrix.slice(k));
        SLogDet = -log(sign*exp(val));
        for(int i = 0; i < n; i++)
            lD(i,k) = log(prop(k)) + ldcppmvt(trans(Xd.row(i)), Mu.col(k), PrecisionMatrix.slice(k), SLogDet);
        
    };
    D = exp(lD);
    D_rows_sums = sum(D, 1);
    
    double penloglik = sum(log(D_rows_sums)) - (lambda * accu(abs(Mu))) - (rho * accu(abs(PrecisionMatrix)));
    return penloglik;
};



void Mixture::GetProbCond(void){
    mat T = zeros<mat>(n, nbClust), lD = zeros<mat>(n, nbClust), D = zeros<mat>(n, nbClust);
    colvec D_rows_sums = zeros<colvec>(n);
    double SLogDet = 0., sign, val;
    
    for(int k = 0; k < nbClust; k++)
    {
        
        log_det(val, sign, PrecisionMatrix.slice(k));
        SLogDet = -log(sign*exp(val));
        
        for(int i = 0; i < n; i++)
        {
            lD(i,k) = log(prop(k)) + ldcppmvt(trans(Xd.row(i)), Mu.col(k), PrecisionMatrix.slice(k), SLogDet);
            
        };
    };
    
    //ici une astuce nuémérique de Gérard Govaert
    colvec tilde = max(lD, 1);
    for(int k = 0; k < nbClust; k++)
        lD.col(k) = lD.col(k) - tilde;
    
    D = exp(lD);
    D_rows_sums = sum(D, 1);
    
    for(int i = 0; i < n; i++)
        T.row(i) = D.row(i)/D_rows_sums(i);
    
    ProbCond = T;
};


void Mixture::GetClassesSizes(void){
    prop = mean(ProbCond, 0);
    
}



void Mixture::UpdateMeans(void)
{
    
    long double T, Tabs;
    rowvec ProbCond_col_sums =  sum(ProbCond, 0);
    mat MusNew = zeros<mat>(p, nbClust), W;
    colvec MuPrec, MuNew;
    
    
    for(int k = 0; k < nbClust; k++){
        W = PrecisionMatrix.slice(k);
        MuPrec = Mu.col(k);
        MuNew  = zeros<colvec>(p);
        
        for(int j = 0; j < p; j++){
            Tabs = 0.;
            for(int i = 0; i < n; i++)
                Tabs += ProbCond(i,k) * (dot(Xd.row(i) - trans(MuPrec), trans(W.col(j))) + (MuPrec(j)*W(j,j)));
            
            
            Tabs = abs(Tabs);
            
            if(Tabs <= lambda) MuNew(j) = 0.;
            else{
                
                //le second membre de l'equation (15)  papier Zhou, Pan , Shen 2009 EJS à k et j fixés
                T = 0.;
                for(int i = 0; i < n; i++)
                    T +=  ProbCond(i,k) * dot(Xd.row(i),W.row(j));
                
                T -= ProbCond_col_sums(k) * (dot(W.row(j),trans(MuPrec)) - MuPrec(j) * W(j,j));
                if(T < 0.) MuNew(j) = (T + lambda)/(ProbCond_col_sums(k) * W(j,j));
                else
                    MuNew(j) = (T - lambda)/ (ProbCond_col_sums(k) * W(j,j));
            }
        }
        MusNew.col(k) = MuNew;
    }
    
    Mu = MusNew;
};

void Mixture::GetEmpiricalCovariance(void){
    rowvec ProbCond_cols_sums = zeros<rowvec>(nbClust);
    ProbCond_cols_sums = sum(ProbCond, 0);
    mat S = zeros<mat>(p,p);
    for(int k = 0; k < nbClust; k++)
    {
        S.zeros(p,p);
        for(int i = 0; i < n; i++)
            S = S + (ProbCond(i,k) * ((trans(Xd.row(i))- Mu.col(k)) * trans(trans(Xd.row(i))- Mu.col(k))));
        
        EmpiricalCovariance.slice(k) = S/ProbCond_cols_sums(k);
    };
}

void Mixture::UpdateCovarianceMatrices(void){
    Environment glasso("package:glasso");
    Function Rglasso = glasso["glasso"];
    rowvec ProbCond_cols_sums = zeros<rowvec>(nbClust);
    ProbCond_cols_sums = sum(ProbCond, 0);
    
    for(int k = 0; k < nbClust; k++)
    {
        long double rhotilde = (2 * rho)/ProbCond_cols_sums(k);
        List Glasso = Rglasso(Named("s") = EmpiricalCovariance.slice(k),
                              Named("rho") = rhotilde,
                              Named("penalize.diagonal") = 0,
                              Named("thr") = 0.001,
                              Named("maxit") = 1000);
        CovarianceMatrix.slice(k) = as<mat>(Glasso[0]);
        PrecisionMatrix.slice(k) =   as<mat>(Glasso[1]);
    };
}


rowvec Mixture::VarRole(void){
    rowvec  MuSum = zeros<rowvec>(p);
    rowvec  alive = ones<rowvec>(p);
    
    MuSum = trans(sum(abs(Mu), 1));
    for(int j = 0; j < p; ++j)
        if(MuSum(j) == 0.)
            alive(j) = 0;
    
    return alive;
}

int Mixture::GetNbClust(void)
{
    return(nbClust);
}






//[[Rcpp::export]]
IntegerVector rcppClusteringEMGlasso(List InputList, double l, double r){
    /******************************************************/
    /************ Les étapes EM pénalisé ******************/
    /******************************************************/
    Mixture MyMixture(InputList, l, r);
    //cout << " .... K = .... " << MyMixture.GetNbClust() <<" .... lambda = .... " << l << " ....rho = .... " << r<<endl;
    double PenLogLik_1  =  MyMixture.PenLogLik(), PenLogLik_0 = 0.0;
    int itr = 0;
    while((abs(PenLogLik_1 - PenLogLik_0) > 0.001) && (itr < 250))
    {
        PenLogLik_0 = PenLogLik_1;
        //                    l'étape E :
        //Calcul des proba conditionnelles à base de X, prop, Mu, V, m
        MyMixture.GetProbCond();
        //                    l'étape M :
        //Calcul des proportions : les moyennes des colonnes de ProbCond
        MyMixture.GetClassesSizes();
        //Calcul des moyennes ....mise à jour des moyennes
        MyMixture.UpdateMeans();
        //Calcul des  matrices de covariances empiriques
        MyMixture.GetEmpiricalCovariance();
        //mise à jour des matrices de convariances et de précision
        MyMixture.UpdateCovarianceMatrices();
        //Calcul de la log-vraisemblance pénalisée
        PenLogLik_1 = MyMixture.PenLogLik();
        itr++;
        
    };
    
    //cout <<" ... Done ... K = ... " << MyMixture.GetNbClust() <<" ... lambda = ... " << l << " ... rho = ... " << r <<" ... PenLogLik = ... " << PenLogLik_1 << endl;
    
    return wrap(MyMixture.VarRole());
};




