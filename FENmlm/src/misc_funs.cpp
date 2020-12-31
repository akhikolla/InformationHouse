#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <vector>

using namespace Rcpp;
using namespace std;

// This file contains the algorithm to get the partial derivative wrt dummies

// Different types:
// int, IntegerVector
// NumericVector, NumericMatrix
// bool
// Rcpp::stop("problem here");
// R_CheckUserInterrupt(); // to allow the user to stop the algo if the process is too long

// Some mathematical expressions that might be problematic
// exp log
// Attention: pour les fonctions log ou exp: il faut absolument du numeric vector (et pas du integer ou ca cree du overloading)


// [[Rcpp::export]]
NumericVector cpp_lgamma(NumericVector x){
	// simple function to compute lgamma of a vector

	int n = x.length();
	NumericVector res(n);

	for(int i=0 ; i<n ; i++){
		res[i] = lgamma(x[i]);
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpp_log_a_exp(double a, NumericVector mu, NumericVector exp_mu){
	// faster this way

	int n = mu.length();
	NumericVector res(n);

	for(int i=0 ; i<n ; i++){
		if(mu[i] < 200){
			res[i] = log(a + exp_mu[i]);
		} else {
			res[i] = mu[i];
		}
	}

	return(res);
}

// [[Rcpp::export]]
NumericVector cpp_partialDerivative_other(int iterMax, int Q, int N, double epsDeriv, NumericVector ll_d2,	NumericVector dx_dother, NumericVector init, IntegerMatrix dumMat, IntegerVector nbCluster){
	// takes in:
	// dumMat: the matrix of dummies (n X c) each obs => cluster // must be in cpp index!!!
	// init: the initialisation of the sum of derivatives vector
	// ll_d2: the second derivative
	// dx_dother: the vector of dx_dother

	int iter;

	int i, q, c;
	int index;
	int sum_cases=0;
	bool ok;
	double new_value;
	IntegerVector start(Q), end(Q);

	for(q=0 ; q<Q ; q++){
		// the total number of clusters (eg if man/woman and 10 countries: total of 12 cases)
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	NumericVector clusterDeriv(sum_cases); // the derivatives per cluster, stacked in one vector
	NumericVector sum_lld2(sum_cases);

	// Creation of the sum_lld2
	for(i=0 ; i<N ; i++){
		for(q=0 ; q<Q ; q++){
			index = start[q] + dumMat(i, q);
			sum_lld2[index] += ll_d2(i);
		}
	}

	// the result to return
	NumericVector S(N);
	for(i=0 ; i<N ; i++){
		S[i] = init(i);
	}

	ok = true;
	iter = 0;
	while( ok & (iter<iterMax) ){
		iter++;
		ok = false;

		for(q=0 ; q<Q ; q++){
			R_CheckUserInterrupt();

			// init of the vector
			for(c=start[q] ; c<end[q] ; c++){
				clusterDeriv(c) = 0;
			}

			for(i=0 ; i<N ; i++){
				index = start[q] + dumMat(i, q);
				clusterDeriv(index) += dx_dother(i) + S(i)*ll_d2(i);
			}

			// on finit de calculer clusterDeriv + controle
			for(c=start[q] ; c<end[q] ; c++){
				new_value = -clusterDeriv(c)/sum_lld2[c];
				clusterDeriv(c) = new_value;
				if(fabs(new_value) > epsDeriv){
					ok = true;
				}
			}

			// on ajoute la derivee a S:
			for(i=0 ; i<N ; i++){
				index = start[q] + dumMat(i, q);
				S[i] += clusterDeriv(index);
			}

		}
	}

	// Rprintf("other, nb iter=%i\n", iter);
	if(iter == iterMax){
		Rprintf("[Getting cluster deriv. other] Max iterations reached (%i)\n", iterMax);
	}

	return(S);
}

// Function to get the conditional sum of a matrix
// [[Rcpp::export]]
NumericMatrix cpp_tapply_sum(int Q, NumericMatrix x, IntegerVector dum){
	// Q: nber of classes
	// N: nber of observations
	// x: a matrix
	// dum: the N vector of clusters

	int N = x.nrow();
	int K = x.ncol();

	NumericMatrix res(Q, K);
	int i, q, k;

	for(i=0 ; i<N ; i++){
		q = dum(i) - 1; // we take 1 off => different indexation in C

		for(k=0 ; k<K ; k++){
			res(q, k) += x(i, k);
		}
	}

	return(res);
}

// Function to get the conditional sum of a vector
// [[Rcpp::export]]
NumericVector cpp_tapply_vsum(int Q, NumericVector x, IntegerVector dum){
	// Q: nber of classes
	// x: a matrix
	// dum: the N vector of clusters

	int N = x.size();

	NumericVector res(Q);
	int i, q;

	for(i=0 ; i<N ; i++){
		q = dum(i) - 1; // we take 1 off => different indexation in C
		res(q) += x(i);
	}

	return(res);
}

// similar a table but faster
// [[Rcpp::export]]
NumericVector cpp_table(int Q, IntegerVector dum){
	// Q: nber of classes
	// dum: the N vector of clusters

	int N = dum.size();

	NumericVector res(Q);
	int i, q;

	for(i=0 ; i<N ; i++){
		q = dum(i) - 1; // we take 1 off => different indexation in C
		res(q)++;
	}

	return(res);
}


// [[Rcpp::export]]
IntegerVector cpp_unclassFactor(NumericVector x){
	// x: a sorted integer vector

	int N = x.size();

	IntegerVector res(N);
	int k=1;
	res[0] = 1;

	for(int i=1 ; i<N ; i++){
		if(x(i-1)!=x(i)) k++;
		res[i] = k;
	}

	return(res);
}



//
// Getting the cluster coefficients
//

// [[Rcpp::export]]
NumericMatrix cpp_get_fe_2(SEXP clusterSize,
                           SEXP i_sorted_index_j, SEXP i_sorted_sumFE,
                           SEXP j_sorted_index_i, SEXP j_sorted_sumFE,
                           SEXP r_cumtable_i, SEXP r_cumtable_j){

	int *psize = INTEGER(clusterSize);
	int n_i = psize[0], n_j = psize[1];
	int nb_coef = n_i + n_j;


	int *is_ind_j = INTEGER(i_sorted_index_j);
	double *is_sumFE = REAL(i_sorted_sumFE);
	int *cumtable_i = INTEGER(r_cumtable_i);

	int *js_ind_i = INTEGER(j_sorted_index_i);
	double *js_sumFE = REAL(j_sorted_sumFE);
	int *cumtable_j = INTEGER(r_cumtable_j);

	// vectors
	vector<bool> isRef_i(n_i, false);
	vector<bool> isRef_j(n_j, false);
	vector<double> cluster_coef_i(n_i);
	vector<double> cluster_coef_j(n_j);
	vector<bool> to_visit_i(n_i, true);
	vector<bool> to_visit_j(n_j, true);
	vector<int> cluster_pending_i(n_i);
	vector<int> cluster_pending_j(n_j);

	int n_done = 0, n_pending_i = 0, n_pending_j = 0, i, j_start = 0, j, m, u;
	while(n_done < nb_coef){

		if(n_pending_j == 0){
			// we set the first guy encountered to 0
			for(j=j_start ; to_visit_j[j] == false ; j++){}
			j_start = j + 1;
			isRef_j[j] = true;
			to_visit_j[j] = false;
			n_done++;
			cluster_coef_j[j] = 0;
			n_pending_j = 1;
			cluster_pending_j[0] = j;
		}

		// if(n_pending_i == 0){
		// 	// we set the first guy encountered to 0
		// 	for(i=i_start ; to_visit_i[i] == false ; i++){}
		// 	i_start = i + 1;
		// 	isRef_i[i] = true;
		// 	to_visit_i[i] = false;
		// 	n_done++;
		// 	cluster_coef_i[i] = 0;
		// 	n_pending_i = 1;
		// 	cluster_pending_i[0] = i;
		// }

		n_pending_i = 0;
		for(m=0 ; m<n_pending_j ; m++){
			j = cluster_pending_j[m];
			for(u = j == 0 ? 0 : cumtable_j[j-1] ; u<cumtable_j[j] ; u++){
				i = js_ind_i[u];
				if(to_visit_i[i]){
					cluster_coef_i[i] = js_sumFE[u] - cluster_coef_j[j];

					to_visit_i[i] = false;
					n_done++;
					cluster_pending_i[n_pending_i] = i;
					n_pending_i++;
				}
			}
		}

		n_pending_j = 0;
		for(m=0 ; m<n_pending_i ; m++){
			i = cluster_pending_i[m];
			for(u = i == 0 ? 0 : cumtable_i[i-1] ; u<cumtable_i[i] ; u++){
				j = is_ind_j[u];
				if(to_visit_j[j]){
					cluster_coef_j[j] = is_sumFE[u] - cluster_coef_i[i];

					to_visit_j[j] = false;
					n_done++;
					cluster_pending_j[n_pending_j] = j;
					n_pending_j++;
				}
			}
		}

	}

	NumericMatrix res(nb_coef, 4);

	for(i=0 ; i<n_i ; i++){
		res(i, 0) = 1;
		res(i, 1) = i + 1;
		res(i, 2) = cluster_coef_i[i];
		res(i, 3) = isRef_i[i];
	}

	for(j=0 ; j<n_j ; j++){
		res(j+n_i, 0) = 2;
		res(j+n_i, 1) = j + 1;
		res(j+n_i, 2) = cluster_coef_j[j];
		res(j+n_i, 3) = isRef_j[j];
	}

	return(res);
}



// [[Rcpp::export]]
List cpp_get_fe_gnl(int Q, int N, NumericVector S, IntegerMatrix dumMat, IntegerVector nbCluster, IntegerVector obsCluster){
	// This function returns a list of the cluster coefficients for each cluster
	// dumMat: the matrix of cluster ID for each observation, with cpp index style
	// Q, N: nber of clusters / obs
	// nbCluster: vector of the number of cases per cluster
	// obsCluster: the integer vector that is equal to order(dum[[g]])
	// RETURN:
	// a list for each cluster of the cluster coefficient value
	// + the last element is the number of clusters that have been set as references (nb_ref)


	int iter=0, iterMax=10000;
	int iter_loop=0, iterMax_loop=10000;


	// Creation of the indices to put all the cluster values into a single vector
	int sum_cases=0;
	IntegerVector start(Q), end(Q), nb_ref(Q); // nb_ref takes the nb of elements set as ref
	int q;

	for(q=0 ; q<Q ; q++){
		// the total number of clusters (eg if c1: man/woman and c2: 10 countries: total of 12 cases)
		sum_cases += nbCluster(q);
		if(q == 0){
			start(q) = 0;
			end(q) = nbCluster(q);
		} else {
			start(q) = start(q-1) + nbCluster(q-1);
			end(q) = end(q-1) + nbCluster(q);
		}
	}

	NumericVector cluster_values(sum_cases);

	// Now we create the vector of observations for each cluster
	// we need a strating and an end vector as well
	IntegerVector start_cluster(sum_cases), end_cluster(sum_cases);

	int i, k, index;

	for(q=0 ; q<Q ; q++){

		// table cluster: nber of elements for each cluster class
		IntegerVector tableCluster(nbCluster(q));
		for(i=0 ; i<N ; i++){
			k = dumMat(i, q);
			tableCluster(k) += 1; // the number of obs per case
		}

		// now creation of the start/end vectors
		for(k=0 ; k<nbCluster(q) ; k++){
			index = start(q) + k;
			if(k == 0){
				start_cluster[index] = 0;
				end_cluster[index] = tableCluster[k];
			} else {
				start_cluster[index] = end_cluster[index-1];
				end_cluster[index] = end_cluster[index-1] + tableCluster[k];
			}
		}
	}

	// matrix of the clusters that have been computed
	IntegerMatrix mat_done(N, Q);

	// vector of elements to loop over
	IntegerVector id2do(N);
	int nb2do = N;
	for(i=0 ; i<nb2do ; i++){
		id2do(i) = i;
	}

	// Other indices and variables
	int qui_max, obs;
	int rs, rs_max;
	int j, ind, id_cluster;
	double other_value;
	bool first, ok;

	//
	// THE MAIN LOOP
	//

	while(iter < iterMax){
		iter++;

		//
		// Finding the row where to put the 0s
		//


		if(iter == 1){
			// 1st iter, we select the first element
			qui_max = 0;
		} else {
			// we find the row that has the maximum of items done

			qui_max = 0;
			rs_max = 0;
			for(i=0 ; i<nb2do ; i++){
				obs = id2do[i];

				rs = 0;
				for(q=0 ; q<Q ; q++){
					rs += mat_done(obs, q);
				}

				if(rs == Q-2){
					// if rs == Q-2 => its the maximum possible, no need to go further
					qui_max = obs;
					break;
				} else if(rs<Q && rs>rs_max){
					// this is to handle complicated cases with more than 3+ clusters
					qui_max = obs;
					rs_max = rs;
				} else if(qui_max == 0 && rs == 0){
					// used to initialize qui_max
					qui_max = obs;
				}
			}
		}

		//
		// Putting the 0s, ie setting the references
		//

		// the first element is spared
		first = true;
		for(q=0 ; q<Q ; q++){
			if(mat_done(qui_max, q) == 0){
				if(first){
					// we spare the first element
					first = false;
				} else {
					// we set the cluster to 0
					// 1) we find the cluster
					id_cluster = dumMat(qui_max, q);
					// Rprintf("Cluster: %i\n", id_cluster + 1);
					// 2) we get the index of the cluster vector
					index = start[q] + id_cluster;
					// 3) we set the cluster value to 0
					cluster_values(index) = 0;
					// 4) we update the mat_done matrix for the elements of this cluster
					for(i=start_cluster[index] ; i<end_cluster[index] ; i++){
						obs = obsCluster(i, q);
						mat_done(obs, q) = 1;
					}
					// 5) => we save the information on which cluster was set as a reference
					nb_ref(q)++;
				}
			}
		}

		//
		// LOOP OF ALL OTHER UPDATES (CRITICAL)
		//

		iter_loop = 0;
		while(iter_loop < iterMax_loop){
			iter_loop++;

			R_CheckUserInterrupt();

			//
			// Selection of indexes (new way) to be updated
			//

			IntegerVector qui_indexes(sum_cases), qui_obs(sum_cases); // init at the max
			int nb_index = 0;

			for(i=0 ; i<nb2do ; i++){
				// we compute the rowsum of the obs that still needs to be done
				obs = id2do[i];

				rs = 0;
				for(q=0 ; q<Q ; q++){
					rs += mat_done(obs, q);
				}

				if(rs == Q-1){
					// means: needs to be updated
					for(q=0 ; mat_done(obs, q)!=0 ; q++){
						// selection of the q that is equal to 0
					}

					index = start(q) + dumMat(obs, q);

					ok = true;
					for(j=0 ; j<nb_index ; j++){
						if(index == qui_indexes[j]){
							ok = false;
							break;
						}
					}

					if(ok){
						// the index was not already there
						qui_indexes[nb_index] = index; // the 'cluster index'
						qui_obs[nb_index] = obs; // the observation to be updated
						nb_index++;
					}
				}
			}

			if(nb_index == 0) break;

			// Rprintf("nb = %i\n", nb_index);

			// => we obtain a unique list of index to be updated

			//
			// We update each index
			//

			for(ind=0 ; ind<nb_index ; ind++){

				int index_select = qui_indexes[ind];
				int q_select=0;

				other_value = 0;
				// Computing the sum of the other cluster values
				// and finding the cluster to be updated (q_select)
				for(q=0 ; q<Q ; q++){
					// we can loop over all q because cluster_values is initialized to 0
					obs = qui_obs[ind];
					index = start(q) + dumMat(obs, q);
					other_value += cluster_values(index);

					if(index == index_select){
						q_select = q;
					}

				}

				// the index to update
				cluster_values(index_select) = S(obs) - other_value;

				// Update of the mat_done and the id2do
				for(i=start_cluster[index_select] ; i<end_cluster[index_select] ; i++){
					obs = obsCluster(i, q_select);
					mat_done(obs, q_select) = 1;
				}
			}
		}

		// Check that everything is all right
		if(iter_loop == iterMax_loop){
			Rprintf("Problem getting FE, maximum iterations reached (2nd order loop).");
		}

		// now the control
		IntegerVector id2do_new(nb2do);
		int nb2do_new = 0;

		for(i=0 ; i<nb2do ; i++){
			obs = id2do[i];

			rs = 0;
			for(q=0 ; q<Q ; q++){
				rs += mat_done(obs, q);
			}

			if(rs < Q){
				id2do_new[nb2do_new] = obs;
				nb2do_new++;
			}
		}

		if(nb2do_new == 0) break;

		// update id2do
		IntegerVector id2do(nb2do_new);
		int nb2do = nb2do_new;

		for(i=0 ; i<nb2do ; i++){
			id2do[i] = id2do_new[i];
		}

	}

	if(iter == iterMax){
		Rprintf("Problem getting FE, maximum iterations reached (1st order loop).");
	}

	// final formatting and save
	List res(Q + 1);
	for(q=0 ; q<Q ; q++){
		NumericVector quoi(nbCluster(q));
		for(k=0 ; k<nbCluster(q) ; k++){
			index = start(q) + k;
			quoi(k) = cluster_values(index);
		}
		res(q) = quoi;
	}

	res(Q) = nb_ref;

	return(res);
}

