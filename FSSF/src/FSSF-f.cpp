/*
 Purpose:
 Generate fully sequential space-filling design using FSSF-f algorithm.
 
 Authors:
 Boyang Shang and Daniel Apley
 
 Inputs:
 d = dimension of the design space;
 nMax = maximum number of (additional) samples required
 N = size of candidate set. Large N will lead to better space-filling performance, but worse computation expense. N = -1 corresponds to default value, and will be computed as N = 1000d + 2nMax;
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
NumericMatrix fssf_f(int d, int nMax, int N = -1, Nullable<NumericVector> ScaleVector = R_NilValue, Nullable<NumericMatrix> Dinit = R_NilValue){
    
    //size of candidate set
    if (N == -1) {
        N = 1000*d + 2*nMax;
    }
    
    if (N < nMax) {
        Rcout << "Candidate set size is too small." << endl;
        return 0;
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
    
    //generate sobol set of size N, dimension d;
    
    long long int seed = ( long long ) skip;
    
    
    for (int j = 0; j < N; j++ )
    {
        i8_sobol ( d, &seed, &Data[j][0] );
    }
    
    vector<double> myScaleVector(d, 1.0);
    if (ScaleVector.isNotNull()) {
        NumericVector RScaleVector(ScaleVector);
        for (int j = 0; j < d; j++) {
            myScaleVector[j] = 1.0/pow(RScaleVector[j],2);
        }
    }
    
    
    vector<int> s;
    int count;
    
    int value_in_push_back;
    value_in_push_back = R::runif(0,1) * (N-1);
    //s.push_back(rand() % N);
    
    
    my_vector_of_vectors_t d_matrix;
    
    d_matrix.resize(2);
    d_matrix[0].resize(N);
    d_matrix[1].resize(N);
    
    double* ptr_1;
    double* ptr_2;
    
    //figure out the first design point
    if (Dinit.isNull()) {
        
        //if the user did not provide an initial design set
        s.push_back(value_in_push_back);
        //initialize the array d_matrix
        for (auto j=0; j<N; ++j) {
            d_matrix[0][j] = 0.0;
            ptr_1 = &(Data[j][0]);
            ptr_2 = &(Data[s[0]][0]);
            for (auto t=0; t<d; ++t) {
 
                d_matrix[0][j] += myScaleVector[t] * pow(*(ptr_1++) - *(ptr_2++), 2);
                
            }
        }
    }else{
        
        NumericMatrix RDinit(Dinit);
        //if the user provided an initial design set, we figure out the first design point by computing the d_matrix using D_init
        int nrow = RDinit.nrow();
        
        //initialize d_matrix using the first row of Dinit
        for (auto j = 0; j<N; ++j) {
            d_matrix[0][j] = 0.0;
            
            for (auto t = 0; t<d; ++t) {
                d_matrix[0][j] += myScaleVector[t] * pow(Data[j][t] - RDinit(0,t), 2);
            }
        }
        
        if (nrow > 1) {
            for (int i = 1; i < nrow; ++i) {
                for (auto j = 0; j<N; ++j) {
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
        
        //select first design point to add
        value_in_push_back = 0;
        double max_dist = d_matrix[0][0];
        
        for (int j = 1; j < N; ++j) {
            if (d_matrix[0][j] > max_dist) {
                max_dist = d_matrix[0][j];
                value_in_push_back = j;
            }
        }
        
        s.push_back(value_in_push_back);
        
        for (auto j=0; j<N; ++j) {
            d_matrix[1][j] = 0.0;
            ptr_1 = &(Data[j][0]);
            ptr_2 = &(Data[s[0]][0]);
            for (auto t=0; t<d; ++t) {
                
                d_matrix[1][j] += myScaleVector[t] * pow(*(ptr_1++) - *(ptr_2++), 2);
                
            }
            
            if (d_matrix[1][j] < d_matrix[0][j]) {
                d_matrix[0][j] = d_matrix[1][j];
            }
            
        }
    }
    
    
    
    
   
    
    vector<bool> Tracking(N, true);
    Tracking[s[0]] = false;
    
    
    int max_idx;
    double max_dist;
    
    
    
    for (auto k=1; k<nMax; ++k) {
        
        
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
