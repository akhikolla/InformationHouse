/*
 Purpose:
 Generate fully sequential space-filling design using FSSF-fr algorithm.
 
 Authors:
 Boyang Shang and Daniel Apley
 
 Inputs:
 d = dimension of the design space;
 nMax = maximum number of samples required or maximum number of sample points to add
 N = size of candidate set. Large N will lead to better space-filling performance, but worse computation expense. N = -1 corresponds to default value, and will be computed as N = 1000d + 2nMax;
 Preference = "minimax" or "maximin". Default is "minimax".  If 'minimax', the early-stage design points will be away from the boundary but relatively close to each other. If 'maximin', the early-stage design will be far away from each other but a little closer to the boundary. Data type is string.
 ScaleVector = Lengthparameters of different inputs. This is an array of size d. For example, if d = 2, ScaleVector = (theta_1, theta_2), the the distance between (x_1, x_2), and (z_1, z_2) will be computed as (z_1-x_1)^2/theta_1^2 + (z_2-x_2)^2/theta_2^2. Default is NA, in which case ScaleVector will be set as a unit vector.
 Dinit = an array of size n_init * d. This is the initial design provided by the user, and the fssf_f function will select additional design points taking Dinit into account.
 
 Outputs:
 Result_matrix = the produced fully-sequential space-filling design. It is a nMax * d array with double precision.
 
 */

# include <RcppArmadillo.h>
# include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

#include <iostream>
#include <fstream>
#include <list>
#include <set>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <map>
#include <math.h>
#include "sobol_dataset.h"


using namespace std;


typedef std::vector<std::vector<double> > my_vector_of_vectors_t;
typedef std::vector<std::vector<size_t> > my_vector_of_vectors_index;


// [[Rcpp::export()]]
NumericMatrix fssf_fr(int d, int nMax, int N = -1, String Preference = "minimax", Nullable<NumericVector> ScaleVector = R_NilValue, Nullable<NumericMatrix> Dinit = R_NilValue){
    
    //size of candidate set
    if (N == -1) {
        N = 1000*d + 2*nMax;
    }
    
    if (N < nMax) {
        Rcout << "Candidate set size is too small." << endl;
        return 0;
    }
    
    //squared factor
    int factor2;
    
    if (Preference == "minimax") {
        factor2 = 2 * d;
    }else{
        factor2 = d * d;
    }
    
    
    //skip is the number of samples to skip in sobol original sequence; we set is as a random number between 0 and 2000.
    int skip;
    //srand((unsigned)time(0));
    //Rcpp::set.seed((unsigned)time(0));
    //skip = (rand()%2000)+1;
    skip = R::runif(0,1) * 1999 + 1;
    
    
    //genereage sobol sequence of size N, dimension d
    //2d array data;
    double** Data = new double*[N];
    for (int i=0; i<N; ++i) {
        Data[i] = new double[d];
    }
    
    
    long long int seed = ( long long ) skip;
    
    for (int j = 0; j < N; j++ )
    {
        i8_sobol ( d, &seed, &Data[j][0] );
    }
    
    //scalevector for each dimension
    vector<double> myScaleVector(d, 1.0);
    if (ScaleVector.isNotNull()) {
        NumericVector RScaleVector(ScaleVector);
        for (int j = 0; j < d; j++) {
            myScaleVector[j] = 1.0/RScaleVector[j];
        }
    }
    
    
    vector<int> s;
    
    my_vector_of_vectors_t d_matrix;
    
    d_matrix.resize(2);
    d_matrix[0].resize(N);
    d_matrix[1].resize(N);
    
    double* ptr_1;
    double* ptr_2;
    
    
    //consider all reflected points of data to the boundry
    double dist2boundary;
    int clsidx;
    
    
    for (auto i=0; i<N; ++i) {
        
        
        for (int j=0; j<d; ++j) {
            if (Data[i][j] < 0.5) {
                dist2boundary = 2*Data[i][j];
            }else{
                dist2boundary = 2*(1-Data[i][j]);
            }
            
            dist2boundary *= myScaleVector[j];
            
            if (j==0) {
                d_matrix[0][i] = dist2boundary;
                clsidx = 0;
                
            }else{
                if (dist2boundary < d_matrix[0][i]) {
                    d_matrix[0][i] = dist2boundary;
                    clsidx = j;
                }
            }
        }
        
        
        
         if (ScaleVector.isNotNull()) {
             double scalingFactor = factor2;
             for (int tt = 0; tt < d; ++tt) {
                 //square of the scalingFactor
                 scalingFactor *=  pow(myScaleVector[tt]/myScaleVector[clsidx], 2.0/d);
             }
             
//Rcout << scalingFactor << endl;
             
             d_matrix[0][i] = scalingFactor * pow(d_matrix[0][i],2);
             
             
         }else{
             d_matrix[0][i] = factor2 * pow(d_matrix[0][i],2);
         }
        
        
        
    }
    
    for (int j = 0; j < d; ++j) {
        myScaleVector[j] = pow(myScaleVector[j], 2);
    }
    
    //handle initial design
    if (Dinit.isNotNull()) {
        NumericMatrix RDinit(Dinit);
        
        int nrow = RDinit.nrow();
        
        for (int i = 0; i < nrow; ++i) {
            for (int j = 0; j < N; ++j) {
                d_matrix[1][j] = 0.0;
                
                for (auto t = 0; t<d; ++t) {
                    d_matrix[1][j] += myScaleVector[t] * pow(Data[j][t] - RDinit(i,t), 2);
                }
                
                if (d_matrix[1][j] < d_matrix[0][j]) {
                    d_matrix[0][j] = d_matrix[1][j];
                }
            }
        }
    }
    
    
    vector<bool> Tracking(N, true);
    
    int max_idx;
    double max_dist;
    int count;
    
    
    for (auto k=0; k<nMax; ++k) {
        
        
        max_dist=0.0;
        
        
        for(int i=0; i<N; i++){
            
            if (d_matrix[0][i] >= max_dist and Tracking[i] == true) {
                max_idx = i;
                max_dist = d_matrix[0][i];
            }
            
        }
        
        s.push_back(max_idx);
        Tracking[max_idx] = false;
        
        for (auto j = 0; j<N; ++j) {
            
            if (Tracking[j]==true) {
                d_matrix[1][j] = 0.0;
                ptr_1 = &(Data[j][0]);
                ptr_2 = &(Data[s[k]][0]);
                
                for (auto t=0; t<d; ++t) {
                    d_matrix[1][j] += myScaleVector[t] * pow(*(ptr_1++) - *(ptr_2++), 2);
                }
                
                if (d_matrix[1][j] < d_matrix[0][j]) {
                    d_matrix[0][j] = d_matrix[1][j];
                }
            }
        }
        
    }
    
    
    
    
    NumericMatrix Result_matrix(nMax,d);
    
    
    count = 0;
    for (auto it:s){
        
        for (int j=0; j<d; ++j) {
            Result_matrix(count,j) = Data[it][j];
        }
        
        count++;
        
        if (count == nMax) {
            break;
        }
        
    }
    
    
    for (int i=0; i<N; ++i) {
        delete [] Data[i];
    }
    delete [] Data;
    
    
    return Result_matrix;
    
    
}
