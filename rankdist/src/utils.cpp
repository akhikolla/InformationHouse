#include <Rcpp.h>
#include <math.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix KendallNeighbour(NumericVector rank){
    int nobj = rank.size();
    int nrow = nobj-1;
    int ncol = nobj;
    NumericMatrix ret(nrow,ncol);
    NumericVector::iterator rankitr;
    for (int i = 0; i<nrow; ++i){
        rankitr=rank.begin();
        for (int j=0; j<ncol; ++j){
            if (rankitr[j]!=i+1 && rankitr[j]!=i+2){
                ret(i,j) = rankitr[j];
            } else if (rankitr[j] == i+1){
                ret(i,j) = i+2;
            } else {
                ret(i,j) = i+1;
            }
            
        }
    }
    return(ret);
}

// [[Rcpp::export]]
NumericMatrix CayleyNeighbour(NumericVector rank){
    typedef NumericVector::iterator vec_iterator;
    vec_iterator rankitr;
    rankitr = rank.begin();
    int max_rank = 0;
    int nobj = rank.size();
    for (int j=0; j<nobj; ++j){
        if(max_rank<rankitr[j])
          max_rank = rankitr[j];
    }
    
    int nrow = nobj*(nobj-1)/2-(nobj-max_rank+1)*(nobj-max_rank)/2;
    int ncol = nobj;
    NumericMatrix ret(nrow,ncol);
    int ind1=0;
    int ind2=1;
    for (int i = 0; i<nrow; ++i){
        rankitr=rank.begin();
        for (int j=0; j<ncol; ++j){
            ret(i,j) = rankitr[j];
        }
        while(ret(i,ind1)==ret(i,ind2)){
          if(ret(i,ind2)==ret(i,ind1) && ind2==ncol-1){
            ++ind1;
            ind2 = ind1;
          }
          ++ind2;
        }
            
        
        
        double tmp = ret(i,ind1);
        ret(i,ind1) = ret(i,ind2);
        ret(i,ind2) = tmp;
        ++ind2;
        if (ind2 == ncol){
            ++ind1;
            ind2 = ind1+1;
        }
    }
    return(ret);
}

// [[Rcpp::export]]
double LogC(NumericVector fai){
    double acc_log = 0;
    double acc_e;
    double acc_w;
    int i_max = fai.size();
    NumericVector w(i_max);
    NumericVector::iterator w_itr=w.begin();
    NumericVector::iterator fai_itr=fai.begin();
    
    for (int i=0; i<i_max; ++i){
        w_itr[i] = 0;
        for (int j = i; j<i_max; ++j){
            w_itr[i] = w_itr[i] + fai_itr[j];
        }
    }
    
    for (int i=1; i<=i_max; ++i){
        acc_e = 0;
        for (int j=1;j<=i;++j){
            acc_w = 0;
            for (int k=j; k<=i; ++k){
                acc_w = acc_w - w_itr[i_max-k];
            }
            acc_e = acc_e + exp(acc_w);
        }
        acc_log = acc_log + log(acc_e + 1);
    }
    return (acc_log);
}


// [[Rcpp::export]]
NumericVector CWeightGivenPi(NumericVector r1, NumericVector r2){

	int nobj = r2.size();
	int nrow = r1.size()/nobj;
	typedef NumericVector::iterator vec_iterator;
	vec_iterator itr1=r1.begin(),itr2=r2.begin();
	
	double L,I;
	NumericVector w(nrow*(nobj-1));
	vec_iterator itrw=w.begin();
	for (int k=0; k<nrow; ++k){
		for (int i=0; i<nobj; ++i){
			I=0;
			L=0;
			for (int j=0 ;j<nobj; ++j){
				if( ((itr1[i*nrow+k]>itr1[j*nrow+k]) && (itr2[i]<itr2[j])) || ((itr1[i*nrow+k]<itr1[j*nrow+k]) && (itr2[i]>itr2[j]))){
					++I;
				}
			}
			L= (itr1[i*nrow+k] + itr2[i] + I)/2;
			if(itr1[i*nrow+k] <= (L-1)){
				for (int p=itr1[i*nrow+k]-1; p<L-1;++p){
					itrw[p*nrow+k] = itrw[p*nrow+k]+0.5;
				}
			}
			if (itr2[i]<=(L-1)){
				for (int p=itr2[i]-1; p<L-1;++p){
					itrw[p*nrow+k] = itrw[p*nrow+k]+0.5;
				}
			}
		}
		for (int i=1; i<nobj-1; ++i){
			itrw[i*nrow+k] = itrw[(i-1)*nrow+k] + itrw[i*nrow+k];
		}
		// handle partial rankings
		int max_rank=0;
		for (int i=0 ;i<nobj; ++i){
		  if (max_rank < itr1[i*nrow+k])
		    max_rank = itr1[i*nrow+k];
		}
		for (int i=max_rank-1; i<nobj-1;++i){
		  itrw[i*nrow+k] = 0;
		}
	}
	return w;
}


// [[Rcpp::export]]
NumericMatrix FindV(NumericMatrix obs, NumericVector pi0){
    NumericMatrix ret(obs.nrow(),obs.ncol()-1);
    NumericVector pi0_perm(pi0.size());
    NumericVector pi0_dummy(pi0.size());
    
    int index;
    for (double i=0; i<pi0.size(); ++i){
        index = int(pi0(i)-1);
        pi0_perm(index) = i; 
    }
    for (int i=0;i<obs.nrow();++i){
        pi0_dummy=seq_len(pi0.size());
        for (int j=0;j<ret.ncol();++j){
            ret(i,j) = std::fabs(obs(i,pi0_perm(j))-pi0_dummy(j));
            for (int k=j; k<ret.ncol();++k){
                if (obs(i,pi0_perm(k)) < obs(i,pi0_perm(j))){
                    --pi0_dummy(k);
                }
            }
        }
    }
    
    return(ret);
}

// [[Rcpp::export]]
double LogC_Component(NumericVector fai){
    int t = fai.size()+1;
    double acc_log=0,acc_term1=0,acc_term2=0;
    for (int i=0;i<t-1;++i){
        acc_term1 += log(1-exp(-1*(t-i)*fai(i)));
        acc_term2 += log(1-exp(-1*fai(i)));
    }
    acc_log = acc_term1 - acc_term2;
    return (acc_log);
}

// [[Rcpp::export]]
int cycle_decomp(NumericVector comp){
    int nobj = comp.size();
    int num_cycles=0, num_visited=0, item= 0;
    bool*visited = new bool[ nobj ];
    for (int i = 0 ; i < nobj; i ++ )visited[ i ] =false;
    while(num_visited < nobj ){
        item = num_cycles;
        while ( visited[ item ]) item++;
        num_cycles++;
        int max_item_inobjcycle= 0;
        do{
            if ( item > max_item_inobjcycle ) max_item_inobjcycle = item;
            visited[ item ] =true;
            num_visited++;
            item = comp(item)-1;
        }while ( !visited[ item ] );
    }
    delete [] visited;
    return (nobj - num_cycles );
}

// [[Rcpp::export]]
NumericVector FindCayley(NumericMatrix obs, NumericVector pi0){
    int nrow = obs.nrow();
    int nobj = obs.ncol();
    NumericVector ret(nrow);
    NumericVector sigma_inv(nobj), comp(nobj);
    for (int this_obs=0; this_obs<nrow; ++this_obs){
        for(int j = 0 ; j < nobj ; j++) sigma_inv(pi0(j) - 1) = j + 1;
        for(int i = 0 ; i < nobj ; i++) comp(i) = obs(this_obs, sigma_inv(i) - 1);
        ret(this_obs) = cycle_decomp(comp);
    }
    return ret;
}

// [[Rcpp::export]]
NumericMatrix Wtau(NumericMatrix obs, NumericVector pi0){
    int param_len = obs.ncol()*(obs.ncol()-1)/2;
    int nobj = obs.ncol();
    NumericMatrix ret(obs.nrow(), param_len);
    int filling_pos = 0;
    // each observation
    for (int ind_obs=0; ind_obs<obs.nrow(); ++ind_obs){
        filling_pos = 0;
        // columns
        for (int i=1; i<nobj; ++i){
            // rows
            for (int j=0; j<i; ++j){
                double equal_sign = (pi0(j) - pi0(i)) * (obs(ind_obs,j) - obs(ind_obs,i));
                ret(ind_obs, filling_pos) = (equal_sign < 0) ? 1:0;
                ++filling_pos;
            }
        }
    }
    
    return (ret);
}



//[[Rcpp::export]]
NumericMatrix AllPerms(int nobj){
    int fac_nobj = 1;
    for(int b = 1; b <= nobj; b++) {
        fac_nobj *= b;
    }
    NumericMatrix ret(fac_nobj ,nobj);
    int *ranking = new int[nobj];
    for (int i=0; i<nobj; ++i){
        ranking[i] = i+1;
    }
    int current_row = 0;
    do {
        for (int i=0; i<nobj; ++i){
            ret(current_row, i) = ranking[i];
        }
        ++current_row;
    } while ( std::next_permutation(ranking,ranking+nobj) );
    return(ret);
}




