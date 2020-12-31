#include <Rcpp.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h> //qsort
#include <limits.h>
#include <math.h>
#include <vector>

using namespace std;
using namespace Rcpp;
//// Base of this code by 
//   J. Chave and F. Jabot
///  last update 05-23-2008

void calcLogKDA(std::vector<long double>& K, 
                long double J, 
                int numspecies, 
                std::vector<int> Abund)  {           
   
    if (Abund.size() < 1) return;
   
    sort(Abund.begin(), Abund.end());

    J=0;
    long double SPP = numspecies;
    for(int s=0;s<SPP;s++) {
        J += Abund[s];
	  }   

    int MaxA = 0;
	  MaxA = Abund[ Abund.size() - 1]; //MaxA = Abund[(int)SPP-1];
	

   
    // abundance distribution
    std::vector<int> Phi(MaxA+1,0); //int *Phi = new int[MaxA+1]; for(s=0;s<=MaxA;s++) Phi[s]=0;
    
   for (int s=0;s<SPP;s++) Phi[Abund[s]]++;

   // Number of distinct abundances
   int NDA=0;
   for (int s=0;s<=MaxA;s++) if(Phi[s] > 0) {NDA++;}
    
    
    //cerr << "Start computing Stirling numbers ...\n";
    // FIRST STAGE: compute the Stirling numbers
    // I use the relationship S(n,m)=S(n-1,m-1)-(n-1)S(n-1,m)
    // where n >= m >= 1
    // The relation of interest is sum_m S(n,m)/S(n,1)*S(m,1) * x^m
    // defined to be equal to sum_m T(n,m) * x^m
    // The recurrence relation on T(n,m) is 
    // T(n,m)= T(n-1,m) + T(n-1,m-1)*(m-1)/(n-1)
   
    
	  std::vector<int> f(NDA,0); //int *f = new int[NDA];
	  std::vector<int> g(NDA,0); //int *g = new int[NDA];
    int i = 0;    
    //for(s=0;s<NDA;s++) {f[s]=0;g[s]=0;}
    for(int n=0;n<=MaxA;n++) {
	  if(Phi[n] > 0) {             
        f[i] = Phi[n];                                  
        g[i] = n;                                        
        i++;
        }
	  }
    //long double **T= new long double*[NDA];          // T(n,m) just for the n which are useful
    //T[0] = new long double[g[0]+1];
 
    
	  std::vector< long double > fill;
	  std::vector< std::vector< long double > > T(NDA,fill);
	
	  std::vector<long double> T0(g[0]+1);
	  T[0] = T0;
	  T[0][0]=0;T[0][1]=1;
	
	  if(g[0]!=1)
	  {
		    std::vector<long double> lS2(g[0]+1); //        long double *lS2 = new long double[g[0]+1]; 
        lS2[0]=0;lS2[1]=1;
        for (int n=2;n<=g[0];n++) {
            std::vector<long double> lS1(n+1); //long double *lS1 = new long double[n+1];                
            for(int im=0;im<=n-1;im++) {
                lS1[im] = lS2[im];
            }
            lS1[n]=0;
            for(int im=2;im<=n;im++) {
                lS2[im] = lS1[im]+lS1[im-1]*(im-1)/(n-1); 
            }
        }
        for(int im=2;im<=g[0];im++) {
            T[0][im]=lS2[im];
        }
	  }
	
	
    for (int in=1;in<i;in++) {
        std::vector<long double> Tin(g[in]+1); //T[in]= new long double[g[in]+1];
    		T[in] = Tin;
    		T[in][0]=0;T[in][1]=1;
		
        std::vector<long double> lS2(g[in]+1); //long double *lS2 = new long double[g[in]+1];         
        for(int im=0;im<=g[in-1];im++) {
            lS2[im] = T[in-1][im];
        }
        for (int n=g[in-1]+1;n<=g[in];n++) {
            std::vector<long double> lS1(n+1); //long double *lS1 = new long double[n+1];                
            for(int im=0;im<=n-1;im++) {
                lS1[im] = lS2[im];
            }
            lS1[n]=0;
            for(int im=2;im<=n;im++) {
                lS2[im] = lS1[im]+lS1[im-1]*(im-1)/(n-1); 
            }
        }
        for(int im=2;im<=g[in];im++) {
            T[in][im]=lS2[im];
        }
    }
    // After this stage we have stored in T[i][m] T(g[i],m)
    // with T(n,m) = S(n,m)*S(m,1)/S(n,1) for i>0

    // cerr << "Start computing ln(K(D,A)) ...\n";
    // SECOND STAGE: compute the K(D,A)
    // I follow Etienne's route. Compute the product of polynomials 
    // of length J
    K.clear();
    K.resize(J+2,0.0); //K = new long double[int(J)+1];
		
    std::vector<long double> poly2(J+1,0.0);//long double *poly2 = new long double[int(J)+1];
    K[0]=1;
    int degree = 0;
    int spe=0;
    for(int i=0;i<NDA;i++) // loop over number of distinct abundances
        for(int j=0;j<f[i];j++){ // loop over abundances per class             
            for(int nn=0;nn<=degree;nn++)
                for(int mm=1;mm<=g[i];mm++){
                    if (K[nn]>0){                   
                       poly2[nn+mm] += T[i][mm]*K[nn];
                    }
                    
                }              
            degree += g[i];
            for(int nn=0;nn<=degree;nn++){            
                K[nn] = (poly2[nn]/powl(10,(4500.0/SPP)));
                poly2[nn] = 0.0;
            }
            spe++;
        }
    
    for(int i=int(SPP);i<=J;i++){
        K[i] = logl(K[i]);                                    // now K[A]=ln(K(D,A)/10^4500) in Etienne's paper
    }
    //for(i=0;i<NDA;i++) delete[] T[i];
	//delete[] T;
    //delete[] poly2;
    //delete[] f;
    //delete[] g;

// search of "infinite" values in K[A]
    int borneinf=(SPP-1);
    int bornesup=(J+1);
    long double maxlog=11333.2;
    int infinity=0;
    for(int i= SPP;i<=J;i++){
        if ((K[i]>maxlog) || (K[i] < -maxlog)) {
            infinity=1;
            break;
        }
        borneinf++;
    }  //after that, borneinf=indice next to infinity but before
    for(int i=0;i<=J-SPP;i++){
        if ((K[J-i]>maxlog)||(K[(int)J-i]<-maxlog)) {
            infinity=1;
            break;
        }
        bornesup--;
    }    //after that, bornesup=indice next to infinity but after
  if (infinity==1) {
  //cerr << "WARNING : the sample is too large to compute an exact likelihood, the program is thus doing approximations. The following results are to be taken with caution"<<endl;
  //cerr << "Value of A above which K(D,A) is computed approximately ="<<borneinf<<endl;
  //cerr << "Value of A below which K(D,A) is computed approximately ="<<bornesup<<endl;

  //fitting of the infinite values of K[A] by a polynom of degree 3
    //computing of the derivatives at the critic points
    
    // Rcpp::Rcout << "Infinity == 1 !! You made it!\n" ;
    
    if(borneinf > (int)K.size()) return;
    if(bornesup > (int)(K.size()-1)) return;
    long double Kprimeinf = K[borneinf]-K[borneinf-1];

     
    long double Kprimesup1 = -50.0;  //out of range check
    if((bornesup+1) < (int)K.size()) Kprimesup1 = K[bornesup+1];
    long double Kprimesup2 = 0.0;   //out of range check
    if(bornesup > -0.1) Kprimesup2 = K[bornesup];

    long double Kprimesup = Kprimesup1-Kprimesup2; //K[bornesup+1]-K[bornesup];
    // definition of the parameters of the fitted polynom aX^3+bX^2+cX+d
    long double a,b,c,d;
    //inversion of the linear system of equations (with the Gauss method)
    long double borneinf2=(long double)borneinf*(long double)borneinf;
    long double borneinf3=(long double)borneinf2*(long double)borneinf;
    long double bornesup2=(long double)bornesup*(long double)bornesup;
    long double bornesup3=(long double)bornesup2*(long double)bornesup;
    d=(Kprimesup-3*bornesup2*K[borneinf]/borneinf3+(2*bornesup/(long double)borneinf-3*bornesup2/borneinf2)*(Kprimeinf-3*K[borneinf]/(long double)borneinf)-((1+3*bornesup2/borneinf2-4*bornesup/(long double)borneinf)/(bornesup-2*bornesup2/(long double)borneinf+bornesup3/borneinf2))*(K[bornesup]-bornesup3*K[borneinf]/borneinf3+(bornesup2/(long double)borneinf-bornesup3/borneinf2)*(Kprimeinf-3*K[borneinf]/(long double)borneinf)))/((6*bornesup2/borneinf3)-(6*bornesup/borneinf2)-((1+3*bornesup2/borneinf2-4*bornesup/(long double)borneinf)/(bornesup-2*bornesup2/(long double)borneinf+bornesup3/borneinf2))*(1-3*bornesup2/borneinf2+2*bornesup3/borneinf3));
    c=((K[bornesup]-bornesup3*K[borneinf]/borneinf3+(bornesup2/(long double)borneinf-bornesup3/borneinf2)*(Kprimeinf-3*K[borneinf]/(long double)borneinf))-d*(1-3*bornesup2/borneinf2+2*bornesup3/borneinf3))/(bornesup-2*bornesup2/(long double)borneinf+bornesup3/borneinf2);
    b=(Kprimeinf-3*K[borneinf]/(long double)borneinf+2*c+3*d/(long double)borneinf)/(0.0-(long double)borneinf);
    a=(K[borneinf]-b*borneinf2-c*(long double)borneinf-d)/borneinf3;
    
    //reconstruction of K[A] with the fitted polynom
    for (int i=borneinf+1;i<bornesup;i++) {
     K[i]=(a*i*i*i+b*i*i+c*i+d);
    }
   }
}




// [[Rcpp::export]]
NumericVector calcKDA(NumericVector A)
{
  
  // Rcpp::Rcout << "Hello! This is the calcKDA function!\n";
    //convert abundances from A to Species
	int numspecies = A.size();
	std::vector<int> Abund(numspecies); //Abund = new int[numspecies];       
	                   
	int J = 0;
	for(int s = 0; s < numspecies; ++s) {
	   if(s > (int)Abund.size()) break;
	   Abund[s] = A[s];
	   J += Abund[s];
	}
	//call calcLogKDA
	std::vector<long double> K;
	calcLogKDA(K,J,numspecies,Abund);
	//return K
	   
	//int sizeofK =  J + 1; //I hope!
	int sizeofK = K.size();	
	NumericVector out(sizeofK);
	for(int i = 0; i < sizeofK; ++i) {
	   out[i] = K[i] + 4500.0 * logl(10);
	}

  return out;	
}
