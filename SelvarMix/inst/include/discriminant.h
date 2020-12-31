//
//  discriminant.h
//  
//
//  Created by sedki on 08/04/2014.
//
//

#ifndef _discriminant_h
#define _discriminant_h
class discriminant{
private:
    mat Xd;
    int nbClust;
    int p;
    int n;
    vector<int> labels; 
    vector< vector<int> > idx; 
    cube CovarianceMatrix;
    cube PrecisionMatrix;
    cube EmpiricalCovariance;
    mat Mu;
    rowvec prop;
    double lambda;
    double rho;
    
    
public:
    discriminant(NumericMatrix X_,  IntegerVector labels_,  int nbClust_, double l, double r); // Constructor
    ~discriminant(){} ; // Destructor
    double PenLogLik(void);
    void Initialization(void);
    void UpdateMeans(void);
    void GetEmpiricalCovariance(void);
    void UpdateCovarianceMatrices(void);
    int GetNbClust(void);
    rowvec VarRole(void);
    
};
#endif
