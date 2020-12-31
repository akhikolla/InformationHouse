
#include "vect.hpp"

//****************************************************************************//
Vect::Vect(){}

Vect::Vect(NumericMatrix Data,vector<int> experiments) 
{
  this->Data = Data;
  this->experiments = experiments;
}

Vect::Vect(NumericMatrix Data)
{
  this->Data = Data;
  this->initExperiments();
}
//Initialization of experiments **********************************************// 
void Vect::initExperiments()
{
  for (int i=1; i <= Data.ncol(); ++i)
    experiments.push_back(i);

  
}

//construction of a matrice according to a variable vector 
mat Vect::const_matrix(vector<int> vecteur)
{
  mat Xr(Data.begin(), Data.nrow(), Data.ncol(), false); 
  mat M = zeros<mat>(Data.nrow(), vecteur.size());
  for(int i = 0; i < (int)vecteur.size(); ++i)
    M.col(i) = Xr.col(vecteur[i] - 1); 
  
  return M;
}//end Vect::const_matrix



//Calculation of BICreg 
List Vect::bicReggen(vector<int> vectH, vector<int> vectY, int numr)
{
  double reg = 0.0, sign, val;  
  
  // Ici, H est la matrice réponse. 
  mat H=Vect::const_matrix(vectH);
  int n = H.n_rows,v = H.n_cols;


  //construction of the matrix X 
  int a; 
  if (vectY.empty())
    a=0;
  else
    a=vectY.size();
  
  
  mat X(n,a+1);  
  if (vectY.empty())
    X.col(0) = ones<colvec>(n);
  else
    { 
      mat Y = Vect::const_matrix(vectY);
      Y.insert_cols(0,  ones<colvec>(n));
      X = Y;
    }
  //Parameter estimation 
  mat XtX = X.t() * X;
  mat B = inv_sympd(XtX) * X.t() *H;
  //mat B = pinv(XtX) * X.t() *H;
  double lambda;
  mat A=X*B;
  
  if (numr==3)  //(r=[LC])
    { 
      mat Omega = (1.0/n)*(H.t()*(H-A));
      Omega = 2*M_PI*Omega;
      log_det(val, sign, Omega);
      double det = log(sign*exp(val));
      lambda = ((a+1)*v) + (0.5*v*(v+1));
      reg = (-n*det)-(n*v)-(lambda*log(n)); 
    }
  
  if (numr==2) //(r=[LB])
    { 
      mat H_A = (1.0/n)*(H - A)%(H - A);
      rowvec sigma2 = sum(H_A, 0);
      sigma2 = log(sigma2);
      lambda=(v*(a+1)) +v;
      reg=-(n*v*log(2*M_PI)) - (n* sum(sigma2)) -(n*v) - (lambda*log(n));   
    }
  
  if (numr==1) //(r=[LI])
    { 
      mat Aux=H-A;
      double sigma=(1.0/(n*v))* accu(Aux % Aux);
      lambda=(v*(a+1)) + 1;
      reg=-(n*v*log(2*M_PI*sigma)) -(n*v) - (lambda*log(n));
    }
  return List::create(Named("bicvalue") = reg, 
                      Named("B") = B);
                    
}//end Vect::bicReggen

//removing a part of a vector
vector<int> Vect::enlever_var(vector<int>& vecteur, vector<int>& varenv)
{
  int taille=0,nb=0;
  vector<int>::iterator it;
  vector<int>res; 
  res = vecteur;
  taille = nb = res.size();
  for (int i=0; i < (int)varenv.size(); ++i)
     {
        for (it=res.begin();it!=res.end();++it)
           {
              if (*it == varenv[i])
                {	  
                   res.erase(it);
                   --nb; 
                   res.resize(taille);
                }
           }
     }
  vector<int>tmp = res;
  res.clear();
  for(int i=0;i<nb;++i)
    res.push_back(tmp[i]);   
  tmp.clear();
  return res;
}//end Vect::enlever_var

//concatenation of two vectors 
vector<int> Vect::ajouter_var(vector<int>& vecteur, vector<int>& varajout)
{
  vector<int> res = vecteur;
  vector<int>::iterator it;
  res.insert(res.end(),varajout.begin(),varajout.end());
  sort(res.begin(),res.end());
return res;
}//end Vect::ajouter_var


