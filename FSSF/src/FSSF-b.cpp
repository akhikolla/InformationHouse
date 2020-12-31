/*
 Purpose:
 Generate fully sequential space-filling design using FSSF-b algorithm.
 
 Authors:
 Boyang Shang and Daniel Apley
 
 Inputs:
 d = dimension of the design space;
 nMax = maximum number of samples required
 N = size of candidate set. Large N will lead to better space-filling performance, but worse computation expense. N = -1 corresponds to default value, and will be computed as N = 1000d + 2nMax;
 eps = the error bound for approximate nearest neighbor search. Larger eps will make the program run faster, but cause worse space-filling performance. Default value is 0.5.
 
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
//#include <ctime>
#include <map>
#include <new>
#include "ANN/ANN.h"
#include <iomanip>
#include "sobol_dataset.h"


using namespace std;

// [[Rcpp::export()]]
NumericMatrix fssf_b(int d, int nMax, int N = -1, double eps = 0.5){
    
    //size of candidate set
    if (N == -1) {
        N = 1000*d + 2*nMax;
    }
    
    if (N < nMax) {
        Rcout << "Candidate set size is too small." << endl;
        return 0;
    }
    
    //number of nearest neighbor to search for each time
    int K;
    if (d >= 1 && d < 4) {
        K = 10;
    }else if(d >= 4 && d < 6){
        K = 15;
    }else if(d >= 6 && d < 8){
        K = 20;
    }else if(d >= 8 && d < 10){
        K = 30;
    }else{
        K = 40;
    }
    
    //skip is the number of samples to skip in sobol original sequence; we set is as a random number between 0 and 2000.
    int skip;
    //srand((unsigned)time(0));
    //Rcpp::set.seed((unsigned)time(0));
    //skip = (rand()%2000)+1;
    skip = R::runif(0,1) * 1999 + 1;
    
    //genereage sobol sequence of size N, dimension d
    ANNpointArray Data;
    //allocate memory
    Data = annAllocPts(N, d);
    ANNcoord** newData;
    
    //generate sobol set of size N, dimension d;
    
    long long int seed = ( long long ) skip;
    
    for (int j = 0; j < N; j++ )
    {
        i8_sobol ( d, &seed, &Data[j][0] );
    }
    
    //the index matrix and the distance matrix
    //index matrix:
    ANNidx** A1 = new ANNidx*[N];
    for(int i = 0; i < N; ++i)
        A1[i] = new ANNidx[K+1];
    //distance matrix
    ANNdist** A2 = new ANNdist*[N];
    for(int i = 0; i < N; ++i)
        A2[i] = new ANNdist[K+1];
    
    //an empty doubly-linked list to hold integers
    list<int> De;
    //an empty doubly-linked list to hold integers
    list<int> D;
    //an continuously stored dynamic array for remaining sample indices
    vector<int> Re;
    //dynamic to store the number of undeleted nearest neighbors for every undeleted sample
    vector<int> Rn;
    //the full index set;
    set<int> full_index;
    //convenient variables
    int n_delete,nhat;
    nhat=N;
    n_delete=0;
    //indicator of whether or not to re-compute the nearest neighbors
    bool startover = true;
    //an dynamic array of linked lists
    vector< list<int> > Id;
    
    //a tree of a with distances as keys and indices as values
    std::multimap<double,size_t> mymap;
    std::pair <std::multimap<double,size_t>::iterator, std::multimap<double,size_t>::iterator> ret;
    
    //convenient variables
    size_t s,t,de;
    double ds,dt,Dde;
    int count;
    
    for (int n=N; n>2; --n) {
        
        if (startover == true) {
            
            //get the indices of the remaining samples
            n_delete = De.size();
            
            if (n_delete != 0) {
                Re.clear();
                for (int i=0; i<N; ++i) {
                    full_index.insert(full_index.end(), i);
                }
                
                
                for (auto ee: De){
                    full_index.erase(ee);
                }
                
                for(auto ee : full_index) {
                    Re.push_back(ee);
                }
                
                
                
                full_index.clear();
                
            }else{
                for (int i=0; i<N; ++i) {
                    Re.push_back(i);
                }
            }
            
            int nhat_before = nhat;
            nhat = Re.size();
            
            
            //if there are less than K samples left
            
            if (nhat <= K) {
                K = nhat-1;
            }
            
            
            if (Rn.size() > 0) {
                Rn.clear();
            }
            
            Rn.assign(nhat,K);
            
            //change sizes of index matrix A1 and distance matrix A2
            
            if (nhat < N) {
                for(int i = 0; i < nhat_before; ++i) {
                    delete [] A1[i];
                    delete [] A2[i];
                }
                delete [] A1;
                delete [] A2;
                
                A1 = new ANNidx*[nhat];
                for(int i = 0; i < nhat; ++i)
                    A1[i] = new ANNidx[K+1];
                //distance matrix
                A2 = new ANNdist*[nhat];
                for(int i = 0; i < nhat; ++i)
                    A2[i] = new ANNdist[K+1];
            }
            
            
            if(n_delete !=0){
                
                //only allocate memory for the location array of each row of newData
                newData = new (nothrow) ANNpoint[nhat];
                if(!newData){
                    Rcout << "memory allocation failed!" << endl;
                    return 0;
                }
                
                for (int ee=0; ee<nhat; ++ee) {
                    
                    newData[ee] = Data[Re[ee]];
                    
                }
                
                //re-construct kd-tree
                ANNkd_tree* kdTree;
                kdTree = new (nothrow) ANNkd_tree(newData, nhat, d);
                if (!kdTree) {
                    Rcout << "memory allocation for kdtree is not successful!" << endl;
                }
                
                //do a (K+1)NN search
                for (int i=0; i<nhat; ++i) {
                    kdTree->annkSearch(newData[i], K+1, A1[i], A2[i], eps);
                }
                
                delete kdTree;
                
                //delete memory allocated for newData
                delete [] newData;
                
            } else{
                ANNkd_tree* kdTree;
                kdTree = new ANNkd_tree(Data, N, d);
                
                for (int i=0; i<nhat; ++i) {
                    kdTree->annkSearch(Data[i], K+1, A1[i], A2[i], eps);
                }
                delete kdTree;
            }
            
            if (D.size() > 0) {
                D.clear();
            }
            //make sure the first column corresponds the query point itself
            for (auto i=0; i<nhat; ++i) {
                if (A1[i][0] != i) {
                    for (auto j=1; j<K+1; ++j) {
                        if (A1[i][j] == i) {
                            A1[i][j] = A1[i][0];
                            break;
                        }
                    }
                }
            }
            
            Id.clear();
            Id.resize(nhat);
            for (auto i=0; i<nhat; ++i) {
                for (auto j=1; j<K+1; ++j) {
                    Id[A1[i][j]].push_front(i);
                }
            }
            
            //first delete all previous containers
            if(mymap.size() > 0)
                mymap.clear();
            //initialize a map for fast query
            for (auto i=0; i<nhat; ++i) {
                mymap.insert(std::pair<double,size_t>(A2[i][1],i));
            }
            
            startover = false;
            
        }
        
        ret = mymap.equal_range(mymap.begin()->first);
        count = std::distance(ret.first,ret.second);
        
        if (count==2) {
            //If the candidate pair is unique
            //find the first element in the map: mymap and stores the corresponding value
            s = mymap.begin()->second;
            t = A1[s][1];
            
            if (A1[t][1] != s) {
                ds = A2[t][1];
                
            }else{
                ds = A2[t][2];
            }
            
            dt = A2[s][2];
            
            
            if (ds >= dt) {
                de = s;
                
                
            } else {
                de=t;
                
            }
            
        }else{
            
            //If the candidate pair is not unique
            de = ret.first->second;
            double MinDist = A2[de][2];
            
            for (std::multimap<double,size_t>::iterator it=ret.first; it!=ret.second; ++it){
                
                if (A2[it->second][2] < MinDist) {
                    
                    de = it->second;
                    MinDist = A2[de][2];
                }
            }
            
        }
        
        
        Dde = A2[de][1];
        
        D.push_back(de);
        De.push_back(Re[de]);
        
        
        //find the location so of Dde in mymap, and delete the node whose value equals de
        
        ret = mymap.equal_range(Dde);
        
        for (std::multimap<double,size_t>::iterator it=ret.first; it!=ret.second; ++it){
            
            if (it->second == de) {
                
                mymap.erase(it);
                break;
            }
        }
        
        Rn[de] = 0;
        
        if(Id[de].size()>0){
            
            for (auto u: Id[de]){
                
                --Rn[u];
                
                if (Rn[u]==1) {
                    startover = true;
                    break;
                    
                }
                
                int j {1};
                for (j=1; j<K+1; ++j) {
                    if (A1[u][j]==de) {
                        
                        break;
                        
                    }
                }
                
                
                if (j==1) {
                    
                    //first check whether pair (A2[u][1],u) exists in mymap or not
                    ret = mymap.equal_range(A2[u][1]);
                    count = std::distance(ret.first,ret.second);
                    
                    //if the key exists
                    if (count > 0) {
                        //check that among those nodes whose key equals A2[u][1], whether there is one whose value equals u
                        for (std::multimap<double,size_t>::iterator it=ret.first; it!=ret.second; ++it){
                            //if such a node indeed exsits
                            if (it->second == u) {
                                //remove this node from mymap
                                mymap.erase(it);
                                
                                //insert a new node to mymap
                                mymap.insert(std::pair<double,size_t>(A2[u][2],u));
                                
                                break;
                            }
                        }
                    }else{
                        continue;
                    }
                }
                
                if(j < K){
                    for (auto i=j; i<K; ++i) {
                        A1[u][i] = A1[u][i+1];
                        A2[u][i] = A2[u][i+1];
                    }
                }
                
                
            }
            
        }
        
        
    }
    
    
    set<int> R;
    for (auto i=0; i<Re.size(); ++i) {
        R.insert(Re[i]);
    }
    
    
    for(auto i: D){
        R.erase(Re[i]);
    }
    
    if (R.size()==2) {
        
        for(auto f : R) {
            De.push_back(f);
        }
    }else{
        
        Rcout << "The final vector is not of the right size! Something is wrong." << endl;
        return 0;
    }
    
    NumericMatrix Result_matrix(nMax,d);
    
    
    count = 1;
    
    
    for (auto it = De.crbegin(); it != De.crend(); ++it){
        
        for (int j=0; j<d; ++j) {
            Result_matrix(count-1,j) = Data[*it][j];
        }
        
        if (count == nMax) {
            break;
        }
        ++count;
        
    }
    
    
    
    
    for(int i = 0; i < nhat; ++i) {
        delete [] A1[i];
        delete [] A2[i];
    }
    delete [] A1;
    delete [] A2;
    
    annDeallocPts(Data, N);
    annClose();
    
    return Result_matrix;
}
