//
//  Mixture.h
//
//
//  Created by sedki on 21/11/2013.
//
//

#ifndef _Mixture_h
#define _Mixture_h



class Mixture{
 private:
    mat Xd;
    int nbClust;
    int p;
    int n;
    cube CovarianceMatrix;
    cube PrecisionMatrix;
    cube EmpiricalCovariance;
    mat Mu;
    rowvec prop;
    double lambda;
    double rho;
    mat ProbCond;
    
 public:
    Mixture(List InputList, double l, double r); // Constructor
    Mixture(List InputList); // Constructor
    ~Mixture(){} ; // Destructor
    double PenLogLik(void);
    void GetProbCond(void);
    void GetClassesSizes(void);
    void UpdateMeans(void);
    void GetEmpiricalCovariance(void);
    void UpdateCovarianceMatrices(void);
    int GetNbClust(void);
    rowvec VarRole(void);
    
};
#endif
