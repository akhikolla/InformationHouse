#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
#define NDEBUG 1
using namespace Rcpp ;
// [[Rcpp::export]]
IntegerVector csample_num( IntegerVector x,
                           int size,
                           bool replace,
                           NumericVector prob = NumericVector::create()
) {
  RNGScope scope;
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

//######################################################################################################################//

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix add_rows(NumericMatrix prior_tree_table_temp,int grow_node){
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
  M(grow_node-1,5)=0;
  M(grow_node-1,6)=0;
  M(grow_node-1,0)=grow_node+1;
  M(grow_node-1,1)=grow_node+2;
  M.insert_rows(grow_node,2);
  M(grow_node,4)=-1;
  M(grow_node+1,4)=-1;
  NumericMatrix t=as<NumericMatrix>(wrap(M));
  IntegerVector rname=seq_len(t.nrow());
  
  List dimnms = // two vec. with static names
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  t.attr("dimnames") = dimnms;
  
  return(t);
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix addcol(NumericMatrix prior_tree_matrix_temp,int grow_node,NumericVector ld_obs,NumericVector rd_obs){
  int ncol=prior_tree_matrix_temp.ncol();
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
  M.insert_cols(ncol,1);
  for(int i =0;i<ld_obs.size();i++){
    // try{
    //   if(ld_obs[i]>prior_tree_matrix_temp.nrow()){
    //     throw std::range_error("can't add col because ld row index is out of range");
    //   }
    // }catch(...){
    //   ::Rf_error("there is a problem adding col to mat don't know why");
    // }
    M(ld_obs[i],ncol)=grow_node+1;
  }
  for(int i =0;i<rd_obs.size();i++){
    // try{
    //   if(rd_obs[i]>prior_tree_matrix_temp.nrow()){
    //     throw std::range_error("can't add col because rd row index is out of range");
    //   }
    // }catch(...){
    //   ::Rf_error("there is a problem adding rd col to mat");
    // }    
    M(rd_obs[i],ncol)=grow_node+2;
  }
  return(wrap(M));
} 

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_daughter_to_end_tree(int grow_node,NumericMatrix prior_tree_table_temp,double left_daughter){
  int nrow=prior_tree_table_temp.nrow();
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
  // Rcout << "Line 82";
  M(grow_node-1,5)=0;
  M(grow_node-1,6)=0;
  M.insert_rows(nrow,2);
  // Rcout << "Line 86";
  
  M(grow_node-1,0)=left_daughter;
  M(grow_node-1,1)=left_daughter+1;
  M(left_daughter-1,4)=-1;
  M(left_daughter,4)=-1;
  // Rcout << "Line 92";
  
  NumericMatrix s=as<NumericMatrix>(wrap(M));
  IntegerVector rname=seq_len(s.nrow());
  
  List dimnms = // two vec. with static names
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  s.attr("dimnames") = dimnms;
  
  return(s);
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_daughter_to_end_mat(int d,NumericMatrix prior_tree_matrix_temp,double left_daughter,NumericVector ld_obs,NumericVector rd_obs){
  int ncol_mat=prior_tree_matrix_temp.ncol();
  arma::mat N=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
  arma::vec colmat=N.col(d);
  NumericVector colmat2=wrap(colmat);
  
  if(d+1==ncol_mat){
    N.insert_cols(ncol_mat,1);
    int nrow_mat=prior_tree_matrix_temp.nrow();
    NumericVector colmatzero(nrow_mat);
    colmatzero[ld_obs]=left_daughter;
    colmatzero[rd_obs]=left_daughter+1;
    //colmat=Rcpp::as<arma::vec>(colmatzero);
    //N.col(d+1)=colmat;
    N.col(d+1)=Rcpp::as<arma::vec>(colmatzero);
    
  }else{
    //else just update existing column
    colmat2[ld_obs]=left_daughter;
    colmat2[rd_obs]=left_daughter+1;
    //colmat=Rcpp::as<arma::vec>(colmat2);
    //N.col(d)=colmat; 
    N.col(d)=Rcpp::as<arma::vec>(colmat2);
    
  }
  
  return(wrap(N));
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector remove_zero(NumericVector nodes_at_depth){
  arma::vec nodes_at_depth2=Rcpp::as<arma::vec>(nodes_at_depth);
  arma::vec ret=nodes_at_depth2.elem(arma::find(nodes_at_depth2!=0));
  return(wrap(ret));
}

//######################################################################################################################//

// [[Rcpp::export]]
IntegerVector order_intvec_(IntegerVector x) {
  IntegerVector sorted = clone(x).sort();
  std::reverse(sorted.begin(), sorted.end());
  
  return match(sorted, x);
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector get_gnp(NumericVector nodes_at_depth,int grow_node){
  arma::uvec grow_node_pos=arma::find(as<arma::vec>(nodes_at_depth)==grow_node);
  
  return(wrap(arma::conv_to<arma::vec>::from(grow_node_pos)));  
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector find_term_nodes(NumericMatrix tree_table){
  arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
  arma::vec colmat=arma_tree.col(4);
  arma::uvec term_nodes=arma::find(colmat==-1);
  term_nodes=term_nodes+1;
  
  return(wrap(arma::conv_to<arma::vec>::from(term_nodes)));
} 

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::uvec find_term_obs(NumericMatrix tree_matrix_temp,double terminal_node){
  arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
  arma::uvec term_obs;
  
  for(int j=0;j<tree_matrix_temp.ncol();j++){
    arma::vec colmat=arma_tree_mat.col(j);
    term_obs=arma::find(colmat==terminal_node);
    if(term_obs.size()>0){
      break;
    }
  }
  
  return(term_obs);
}

//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double likelihood_function(NumericVector y_temp,NumericMatrix treetable_temp,NumericMatrix obs_to_nodes_temp,double a,double mu,double nu,double lambda){
  double tree_log_lik;
  NumericVector terminal_nodes=find_term_nodes(treetable_temp);
  double b=terminal_nodes.size();
  IntegerVector n(b);
  double term1=0;
  double term2=0;
  arma::uvec term_obs;
  double ni;
  arma::vec y_temp2=as<arma::vec>(y_temp);
  for(int i=0;i< b;i++){
    term_obs=find_term_obs(obs_to_nodes_temp,terminal_nodes[i]);
    arma::vec y_k=y_temp2.elem(term_obs);
    ni=term_obs.size();
    //n[i]=ni;
    double ybar=0;
    if(y_k.size()!=0){
      ybar=mean(y_k);
    }
    term1+=log(ni+a);
    arma::vec y_k_sq=pow(y_k,2);
    double sum_yksq=sum(y_k_sq);
    double b2=pow(ni*ybar +a*mu,2)/(ni+a);
    
    
    //term2+=(sum_yksq+a*pow(mu,2)-b2+nu*lambda);
    
    
    term2+=(sum_yksq+a*pow(mu,2)-b2);
  }
  //tree_log_lik=(b/2)*log(a)-0.5*term1-((y_temp.size()+nu)/2)*log(term2);
  
  
  
  tree_log_lik=(b/2)*log(a)-0.5*term1-((y_temp.size()+nu))*0.5*log(term2+nu*lambda);
  
  
  // for(int i=0;i<b;i++){
  //   if(n[i]<=5){
  //     tree_log_lik=tree_log_lik-(100000);
  //   }
  //   else{
  //     tree_log_lik=tree_log_lik;
  //   }
  // }
  return(tree_log_lik);
}  
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::uvec find_internal_nodes(NumericMatrix treetable){
  
  arma::mat arma_tree(treetable.begin(),treetable.nrow(), treetable.ncol(), false); 
  arma::vec colmat=arma_tree.col(4);
  arma::uvec term_nodes=arma::find(colmat==1);
  term_nodes=term_nodes+1;
  
  return(term_nodes);
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double find_prev_nonterm(arma::uvec find_nonterm,NumericVector prev){
  
  double ret=0;
  int z=prev.size();
  for(int j=0;j<z;j++){
    arma::uvec term_equal = arma::find(find_nonterm==prev[j]);    
    ret+=term_equal.size(); 
  }  
  
  return(ret);
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::uvec find_nodes_to_update(arma::uvec all_ld,double left_daughter){
  arma::uvec gr_ld=arma::find(all_ld>=left_daughter);  
  return(gr_ld);
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_tree_to_middle(NumericVector node_to_update,NumericMatrix prior_tree_table_temp,int grow_node,double left_daughter){
  
  // Rcout << "Line 300  set_tree_to_middle()" ;
  
  for(int i=0;i<node_to_update.size();i++){
    if(prior_tree_table_temp(node_to_update[i],0) && prior_tree_table_temp(node_to_update[i],1)!=0){
      prior_tree_table_temp(node_to_update[i],0)+=2;
      prior_tree_table_temp(node_to_update[i],1)+=2;
    }
  }
  
  prior_tree_table_temp(grow_node-1,5)=0;
  prior_tree_table_temp(grow_node-1,6)=0;
  
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
  M.insert_rows(left_daughter-1,2);
  M(left_daughter-1,4)=-1;
  M(left_daughter,4)=-1;      
  
  M(grow_node-1,0)=left_daughter;
  M(grow_node-1,1)=left_daughter+1;
  NumericMatrix t=as<NumericMatrix>(wrap(M));
  IntegerVector rname=seq_len(t.nrow());
  
  List dimnms = // two vec. with static names
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  t.attr("dimnames") = dimnms;
  
  return(t);
}

//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix update_grow_obs(NumericMatrix prior_tree_matrix_temp,double grow_node,double left_daughter,int d,NumericVector ld_obs,NumericVector rd_obs){
  arma::mat prior_tree_matrix_temp2(prior_tree_matrix_temp.begin(),prior_tree_matrix_temp.nrow(),prior_tree_matrix_temp.ncol(),false);
  arma::vec ptm2=prior_tree_matrix_temp2.col(d);
  NumericVector ptm=wrap(ptm2);
  ptm[ld_obs]=left_daughter;
  ptm[rd_obs]=left_daughter+1;
  prior_tree_matrix_temp(_,d)=ptm;  
  
  return(prior_tree_matrix_temp);
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix find_obs_to_update_grow(NumericMatrix prior_tree_matrix_temp,double left_daughter,int d,NumericVector ld_obs,NumericVector rd_obs){
  
  int rows=prior_tree_matrix_temp.nrow();
  int cols=prior_tree_matrix_temp.ncol();
  int elements=rows*cols;
  std::vector<double> rows_obs(elements);
  std::vector<double> cols_obs(elements);
  int count=0;
  for(int i=0;i<prior_tree_matrix_temp.nrow();i++){
    for(int j=0;j<prior_tree_matrix_temp.ncol();j++){
      if(prior_tree_matrix_temp(i,j)>=left_daughter){
        rows_obs[count]=i;
        cols_obs[count]=j;
        count++;
      }
    }
  }
  rows_obs.resize(count);
  cols_obs.resize(count);
  
  if(rows_obs.size()!=0){
    for(int k=0;k< count;k++){
      if(prior_tree_matrix_temp(rows_obs[k],cols_obs[k])<left_daughter){
        prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=prior_tree_matrix_temp(rows_obs[k],cols_obs[k]);
      }else if(prior_tree_matrix_temp(rows_obs[k],cols_obs[k])==0){
        //prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=0;
      }else{   
        //int temp=prior_tree_matrix_temp(rows_obs[k],cols_obs[k])+2;
        //prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=temp;
        prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=prior_tree_matrix_temp(rows_obs[k],cols_obs[k])+2;
        
      }
    }
  }
  
  
  arma::mat prior_tree_matrix_temp2(prior_tree_matrix_temp.begin(),prior_tree_matrix_temp.nrow(),prior_tree_matrix_temp.ncol(),false);	// copy as arma mat with new name
  
  //update prior_tree_matrix
  //insert extra column for the new split node 
  if(prior_tree_matrix_temp.ncol()==d+1){								// If the number of columns is greater than d+1
    //arma::mat M=Rcpp::as<arma::mat>(prior_tree_matrix_temp);		// M equals arma mat copy of the matrix
    //M.insert_cols(prior_tree_matrix_temp.ncol(),1);  				// insert one column after last column
    //NumericMatrix prior_tree_matrix=as<NumericMatrix>(wrap(M));		// convert back to NumericMatrix and save as prior_tree_matrix
    
    prior_tree_matrix_temp2.insert_cols(prior_tree_matrix_temp.ncol(),1);  				// insert one column after last column
  }
  
  
  arma::vec ptm2=prior_tree_matrix_temp2.col(d+1);
  NumericVector ptm=wrap(ptm2);
  ptm[ld_obs]=left_daughter;
  ptm[rd_obs]=left_daughter+1;
  prior_tree_matrix_temp(_,d+1)=ptm;  
  
  return(prior_tree_matrix_temp);
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_daughter_obs(arma::mat& xmat,NumericVector obs_to_update,int split_var,double split_point){
  
  //arma::uvec ld_obs_ind(obs_to_update.size());
  //arma::uvec rd_obs_ind(obs_to_update.size());
  
  // Rcout << "line 414.\n";
  arma::colvec sv_col;
  List daughter_obs(2);
  
  // Rcout << "line 419.\n";
  sv_col=xmat.col(split_var-1);
  // Rcout << "line 421.\n";
  
  arma::vec obs_to_update_arma=as<arma::vec>(obs_to_update);
  arma::uvec obs_to_update_arma_u=as<arma::uvec>(obs_to_update);
  
  arma::uvec ld_obs_ind = arma::find(sv_col.elem(obs_to_update_arma_u)<=split_point);
  arma::uvec rd_obs_ind = arma::find(sv_col.elem(obs_to_update_arma_u)>split_point);
  // Rcout << "line 426.\n";
  
  //NumericVector ld_ind2(as<NumericVector>(wrap(ld_obs_ind)));
  //NumericVector rd_ind2(as<NumericVector>(wrap(rd_obs_ind)));
  // Rcout << "ld_ind2 = "<< ld_ind2 << " .\n" ;
  // Rcout << "ld_ind2 = "<< rd_ind2 << " .\n" ;
  // Rcout << "obs_to_update = "<< obs_to_update << " .\n" ;
  
  //NumericVector ld_obs=obs_to_update[ld_ind2];
  //NumericVector rd_obs=obs_to_update[rd_ind2];
  // Rcout << "line 433.\n";
  arma::vec obs_to_update_arma_left=obs_to_update_arma.elem(ld_obs_ind);
  arma::vec obs_to_update_arma_right=obs_to_update_arma.elem(rd_obs_ind);
  
  //daughter_obs[0]=ld_obs;
  //daughter_obs[1]=rd_obs;
  //NumericVector ld_obs=wrap(obs_to_update_arma_left);
  //NumericVector rd_obs=wrap(obs_to_update_arma_right);
  
  daughter_obs[0]=wrap(obs_to_update_arma_left);
  daughter_obs[1]=wrap(obs_to_update_arma_right);
  
  // Rcout << "line 437.\n";
  
  return(daughter_obs);
  
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

int find_term_cols(NumericMatrix tree_matrix_temp,int terminal_node){
  
  arma::mat tree_matrix_temp2(tree_matrix_temp.begin(),tree_matrix_temp.nrow(),tree_matrix_temp.ncol(),false);
  //int count=0;
  //std::vector<double> term_cols(tree_matrix_temp.ncol());
  int term_col=0;
  for(int j=0;j<tree_matrix_temp.ncol();j++){
    
    arma::vec tempcol=tree_matrix_temp2.col(j);
    arma::uvec term_nodes=find(tempcol==terminal_node);
    
    if(term_nodes.size()>0){  
      //term_cols[count]=j;
      //count++;
      term_col=j;
      break;
    }
  }
  
  //term_cols.resize(count);
  
  //return(wrap(term_cols));
  return(term_col);
  
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector get_grow_obs(arma::mat& xmat,NumericVector grow_obs,int split_var){
  
  arma::vec sv_col=xmat.col(split_var-1);
  arma::uvec grow_obs2(as<arma::uvec>(grow_obs));
  arma::vec get_min=sv_col.elem(grow_obs2);
  
  return(wrap(get_min));
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List grow_tree(arma::mat& xmat,//NumericVector y,
               NumericMatrix prior_tree_matrix,int grow_node,NumericMatrix prior_tree_table,int splitvar,
               double splitpoint,//NumericVector terminal_nodes,
               NumericVector grow_obs,
               int d//,NumericVector get_min,arma::mat& data_curr_node
)
{
  // Rcout << "Line 489.\n";
  
  NumericMatrix prior_tree_matrix_temp=clone(prior_tree_matrix);
  NumericMatrix prior_tree_table_temp=clone(prior_tree_table);
  // Rcout << "Line 496.\n";
  double yy=xmat.n_cols;
  IntegerVector xx=seq_len(yy);
  // Rcout << "Line 499.\n";
  prior_tree_table_temp(grow_node-1,3)=splitpoint;
  prior_tree_table_temp(grow_node-1,2)=splitvar;
  prior_tree_table_temp(grow_node-1,4)=1;
  //get data subset for left and right daughter nodes
  // Rcout << "Line 504.\n";
  // Rcout << "grow_obs= " << grow_obs << ".\n";
  // Rcout << "splitvar= " << splitvar << ".\n";
  // Rcout << "splitpoint= " << splitpoint << ".\n";
  
  List daughter_obs=get_daughter_obs(xmat,grow_obs,splitvar,splitpoint);
  // Rcout << "Line 506.\n";
  NumericVector ld_obs=daughter_obs[0];
  NumericVector rd_obs=daughter_obs[1];
  
  // Rcout << "Line 510.\n";
  
  if(prior_tree_table_temp.nrow()==grow_node){
    // Rcout << "Line 513";
    
    prior_tree_table_temp=add_rows(prior_tree_table_temp,grow_node);
    prior_tree_matrix_temp=addcol(prior_tree_matrix_temp,grow_node,ld_obs,rd_obs);  
  }else{
    // Rcout << "Line 518.\n";
    
    //if grow node is in the middle of the tree
    NumericVector nodes_d;
    nodes_d=prior_tree_matrix_temp(_,d);
    NumericVector nodes_at_depth=sort_unique(nodes_d);
    NumericVector nodes_at_depth1=remove_zero(nodes_at_depth);
    NumericVector gn_pos=get_gnp(nodes_at_depth1, grow_node);
    arma::uvec prev_uvec= arma::find(as<arma::vec>(nodes_at_depth1)<grow_node);
    arma::vec nontermvec=Rcpp::as<arma::vec>(nodes_at_depth1);
    nontermvec=nontermvec.elem(prev_uvec);
    NumericVector prev= as<NumericVector>(wrap(nontermvec));
    double prev_nonterm=0;
    if(prev.size()!=0){
      arma::uvec find_nonterm=find_internal_nodes(prior_tree_table);
      //should only find internal nodes at the current depth
      prev_nonterm=find_prev_nonterm(find_nonterm,prev);
    }
    double left_daughter=grow_node +2*(prev_nonterm)+(nodes_at_depth1.size()-gn_pos[0]);
    NumericVector ptt=prior_tree_table(_,1);
    arma::uvec node_to_update=find_nodes_to_update(as<arma::uvec>(ptt),left_daughter);
    //increase the node number of nodes after the grow node by two (because 2 daughter nodes are now appended to grow node)
    //do this for all observations except those that already belong to a terminal node (a value of 0)
    
    // Rcout << "Line 542.\n";
    // Rcout << "node_to_update= " << node_to_update << ".\n";
    
    if(node_to_update.size()==0){
      if(prior_tree_matrix_temp.ncol()>d+1){
        // Rcout << "Line 559.\n";
        
        prior_tree_table_temp=set_daughter_to_end_tree(grow_node,prior_tree_table_temp,left_daughter);
        // Rcout << "Line 562.\n";
        
        prior_tree_matrix_temp=update_grow_obs(prior_tree_matrix_temp,grow_node,left_daughter,d+1,ld_obs,rd_obs);
      }else{
        // Rcout << "Line 552.\n";
        //if the daughter node number already exists in the tree and existing node numbers have to be updated
        //daughter nodes need to be added to the end of the table not in the center of it
        prior_tree_table_temp=set_daughter_to_end_tree(grow_node,prior_tree_table_temp,left_daughter);
        // Rcout << "Line 568.\n";
        
        prior_tree_matrix_temp=set_daughter_to_end_mat(d,prior_tree_matrix_temp,left_daughter,ld_obs,rd_obs);
      }
    }else{
      // Rcout << "Line 559.\n";
      //if the daughter node number already exists in the tree and existing node numbers have to be updated
      prior_tree_table_temp=set_tree_to_middle(wrap(arma::conv_to<arma::vec>::from(node_to_update)),prior_tree_table_temp,grow_node,left_daughter);
      // Rcout << "Line 562.\n";
      prior_tree_matrix_temp=find_obs_to_update_grow(prior_tree_matrix_temp,left_daughter,d,ld_obs,rd_obs);    
    }
  }
  
  // Rcout << "Line 546.\n";
  
  List ret(2);
  ret[0]=prior_tree_table_temp;
  ret[1]=prior_tree_matrix_temp;
  
  return(ret);
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_daughter(int left_daughter,int right_daughter,IntegerVector ld_obs,IntegerVector rd_obs,NumericMatrix tree_matrix_temp,double term_cols){
  
  arma::mat tree_matrix_temp2(tree_matrix_temp.begin(),tree_matrix_temp.nrow(),tree_matrix_temp.ncol(),false);
  arma::vec arma_col=tree_matrix_temp2.col(term_cols+1);
  NumericVector col(as<NumericVector>(wrap(arma_col)));
  col[ld_obs]=left_daughter;
  col[rd_obs]=right_daughter;  
  tree_matrix_temp(_,term_cols+1)=col;
  
  return(tree_matrix_temp);
}

//######################################################################################################################//

// [[Rcpp::export]]

IntegerVector order_(NumericVector x) { //p osition of largest value, then second largest, and so on
  NumericVector sorted = clone(x).sort();
  std::reverse(sorted.begin(), sorted.end());
  
  return match(sorted, x);
}
//######################################################################################################################//

// [[Rcpp::export]]

IntegerVector orderforOW(NumericVector x) {	// gives vector of position of smallest value, then position of second smallest value, and so on.
  NumericVector sorted = clone(x).sort();		// sorted is x in ascending order
  //std::reverse(sorted.begin(), sorted.end()); // reverse so that it is in descending order. Could use one line with std::sort(clone(x).begin(), clone(x).end(), std::greater<>())
  
  return match(sorted, x);	//  match is the Rcpp sugar version of the R function match, which returns a vector of the positions of the first matches of the first argument in the second.
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double secondKindStirlingNumber(int n, int k) {
  if(k>n)
    throw std::range_error("Sterling number undefined for k>n");
  if(k==0 && n==0)
    return 1;
  if (n == 0 || k == 0 || k > n) 
    return 0; 
  if (k == 1 || k == n) 
    return 1; 
  
  arma::mat sf=arma::zeros(n + 1,n + 1);
  for (int i = 0; i < k+1; i++) {
    sf(i,i) = 1;
  }
  for(int i=1; i< n+1 ; i++){
    sf(i,1)=1;
  }
  for (int i = 3; i < n + 1; i++) {
    for (int j = 2; j < k + 1; j++) {
      sf(i,j) = j * sf(i - 1,j) + sf(i - 1,j - 1);
    }
  }
  return sf(n,k);
}
// //######################################################################################################################//
// // Returns count of different partitions of n 
// // elements in k subsets
// 
// //' @export
// // [[Rcpp::export]]
// double countP(int n, int k) 
// { 
//   // Base cases 
//   if(k>n)
//     throw std::range_error("Sterling number undefined for k>n");
//   if(k==0 && n==0)
//     return 1;
//   if (n == 0 || k == 0 || k > n) 
//     return 0; 
//   if (k == 1 || k == n) 
//     return 1; 
//   
//   // S(n+1, k) = k*S(n, k) + S(n, k-1) 
//   return k*countP(n-1, k) + countP(n-1, k-1); 
// } 
//######################################################################################################################//
#include <math.h>       /* tgamma */
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double get_tree_prior(double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t,
                      double num_obs, double num_vars, double lambda_poisson,
                      NumericMatrix tree_table,NumericMatrix tree_matrix,double alpha,double beta){
  
  
  if(spike_tree==1){
    arma::uvec internal_nodes_prop=find_internal_nodes(tree_table);
    arma::mat tree_table2(tree_table.begin(),tree_table.nrow(),tree_table.ncol(),false);
    
    //Rcout << "line 600.internal_nodes_prop = " << internal_nodes_prop <<".\n";
    
    double k_temp=internal_nodes_prop.size()+1;
    
    //Rcout << "line 600.k_temp = " << k_temp <<".\n";
    
    arma::mat split_var_rows=tree_table2.rows(internal_nodes_prop-1);
    arma::vec split_var_vec=split_var_rows.col(2);
    arma::vec uniquesplitvars=arma::unique(split_var_vec);
    double q_temp=uniquesplitvars.n_elem;
    // Rcout << "line 606.\n";
    //Rcout << "line 600.q_temp = " << q_temp <<".\n";
    // Rcout << "line 600.lambda_poisson = " << lambda_poisson <<".\n";
    // Rcout << "line 600.num_obs = " << num_obs <<".\n";
    // Rcout << "line double(tgamma(1))= " << double(std::tgamma(1)) <<".\n";
    // Rcout << "line double(tgamma(100)) = " << double(std::tgamma(100)) <<".\n";
    // Rcout << "line double(tgamma(400)) = " << double(std::tgamma(400)) <<".\n";
    // Rcout << "line double(tgamma(400)) = " << double(std::tgamma(double(400))) <<".\n";
    // Rcout << "line double(lgamma(400)) = " << double(std::lgamma(double(400))) <<".\n";
    
    
    //Rcout << "line double(lgamma(q_temp+1)) = " << double(std::lgamma(double(q_temp+1))) <<".\n";
    
    //int qint = q_temp;
    //int num_obsint = num_obs;
    
    //FIRST CALCULATE THE log of denom and right_truncatin
    //Then take the exponential
    //then take the difference
    
    double denom=1;
    for(int i=0; i<q_temp+1;i++){
      //denom= denom-(pow(lambda_poisson,double(i))*exp(-lambda_poisson)/double(tgamma(i+1))); 
      denom = denom-exp(i*log(lambda_poisson)-lambda_poisson-std::lgamma(double(i+1)));
    }
    double right_truncation=1;
    for(int i=0; i<num_obs+1;i++){
      //right_truncation= right_truncation-(pow(lambda_poisson,double(i))*std::exp(-lambda_poisson)/double(tgamma(i+1))); 
      right_truncation= right_truncation-exp(i*log(lambda_poisson)-lambda_poisson-std::lgamma(double(i+1)));
    }
    //Rcout << " right_truncation= " << right_truncation << ".\n";
    denom=denom-right_truncation;
    //Rcout << " denom= " << denom << ".\n";
    
    //double propsplit=(tgamma(num_vars+1)/(tgamma(q_temp+1)*(tgamma(num_vars-q_temp+1)))*tgamma(q_temp+1)*tgamma(num_vars-q_temp+1)/tgamma(num_vars+2));
    //Beta(1,1) on binomial probability implies uniform prior on numbers of variables from 0 to p (i.e. number of vars has prob (1/(p+1)))
    //double propsplit=(1/double(num_vars+1))*pow(lambda_poisson,k_temp)*
    //  exp(-lambda_poisson)*(1/double(tgamma(k_temp)*denom))*
    //  (1/(tgamma(num_obs)*pow(q_temp,k_temp-1-q_temp)*
    //  tgamma(q_temp+1)/double(tgamma(num_obs-k_temp+1))));
    
    // double propsplit=(1/double(num_vars+1))*
    //   pow(lambda_poisson,k_temp)*
    //   exp(-lambda_poisson)*(1/double(tgamma(k_temp)*denom))*
    //   (1/exp(lgamma(num_obs)+(k_temp-1-q_temp)*log(q_temp)+
    //   lgamma(q_temp+1)-(lgamma(num_obs-k_temp+1))));
    // 
    // double propsplit=(1/double(num_vars+1))*
    //   exp(  k_temp*log(lambda_poisson)-
    //   -lambda_poisson-std::lgamma(k_temp)-denom  )*
    //   (1/exp(std::lgamma(num_obs)+(k_temp-1-q_temp)*log(q_temp)+
    //   std::lgamma(q_temp+1)-(std::lgamma(num_obs-k_temp+1))));
    
    if(q_temp==0){
      
      
      if(s_t_hyperprior==1){
        double propsplit=//(1/double(num_vars+1))*
          exp(std::lgamma(num_vars+1)-std::lgamma(q_temp+1)-std::lgamma(num_vars-q_temp+1)+
          q_temp*log(p_s_t)+(num_vars-q_temp)*log(1-(p_s_t))+
          k_temp*log(lambda_poisson)-
          lambda_poisson-std::lgamma(k_temp+1)-denom)  ;
        //Rcout << " propsplit= " << propsplit << ".\n";
        return(propsplit);
      }else{
        double propsplit=//(1/double(num_vars+1))*
          exp(std::lgamma(num_vars+1)-std::lgamma(q_temp+1)-std::lgamma(num_vars-q_temp+1)+
          std::lgamma(q_temp+a_s_t)+std::lgamma(num_vars-q_temp+b_s_t)-std::lgamma(num_vars+a_s_t+b_s_t)+
          k_temp*log(lambda_poisson)-
          lambda_poisson-std::lgamma(k_temp+1)-denom)  ;
        //Rcout << " propsplit= " << propsplit << ".\n";
        return(propsplit);
        
      }
      
      
      
      
      
    }else{
      if(s_t_hyperprior==1){
        double propsplit=//(1/double(num_vars+1))*
          exp(  std::lgamma(num_vars+1)-std::lgamma(q_temp+1)-std::lgamma(num_vars-q_temp+1)+
          std::lgamma(q_temp+a_s_t)+std::lgamma(num_vars-q_temp+b_s_t)-std::lgamma(num_vars+a_s_t+b_s_t)+
          k_temp*log(lambda_poisson)-
          lambda_poisson-std::lgamma(k_temp+1)-denom  -
          (std::lgamma(num_obs)+log(secondKindStirlingNumber(k_temp-1,q_temp))+//(k_temp-1-q_temp)*log(q_temp)+
          std::lgamma(q_temp+1)-(std::lgamma(num_obs-k_temp+1))));
        //Rcout << " propsplit= " << propsplit << ".\n";
        return(propsplit);
      }else{
        double propsplit=//(1/double(num_vars+1))*
          exp(  std::lgamma(num_vars+1)-std::lgamma(q_temp+1)-std::lgamma(num_vars-q_temp+1)+
          q_temp*log(p_s_t)+(num_vars-q_temp)*log(1-(p_s_t))+
          k_temp*log(lambda_poisson)-
          lambda_poisson-std::lgamma(k_temp+1)-denom  -
          (std::lgamma(num_obs)+log(secondKindStirlingNumber(k_temp-1,q_temp))+//(k_temp-1-q_temp)*log(q_temp)+
          std::lgamma(q_temp+1)-(std::lgamma(num_obs-k_temp+1))));
        //Rcout << " propsplit= " << propsplit << ".\n";
        return(propsplit);
      }
      
    }
    
    // Rcout << "line 621.\n";
    //Rcout << " propsplit= " << propsplit << ".\n";
    
    
    
    //return(propsplit);
    
  }else{
    
    double propsplit=1;
    IntegerVector d;
    IntegerVector d2;
    //int col=tree_matrix.ncol();
    //std::vector<int> int_nodes_index(100*col); // Why  100* ? 
    //std::vector<int> term_nodes_index(100*col); //
    
    //int index_count=0;  
    //int index_count2=0;
    arma::uvec internal_nodes_prop=find_internal_nodes(tree_table);
    //NumericVector terminal_nodes_prop_wrapped=find_term_nodes(tree_table);
    
    arma::mat tree_table2(tree_table.begin(),tree_table.nrow(),tree_table.ncol(),false);
    arma::vec colmat=tree_table2.col(4);
    arma::uvec terminal_nodes_prop=arma::find(colmat==-1)+1;
    
    //arma::uvec terminal_nodes_prop=as<arma::uvec>(terminal_nodes_prop_wrapped);
    arma::mat tree_matrix2(tree_matrix.begin(),tree_matrix.nrow(),tree_matrix.ncol(),false);
    int count=internal_nodes_prop.size();
    int count_term_nodes=terminal_nodes_prop.size();
    
    if(count==0){
      propsplit=1-alpha;
    }else{
      
      for(int k=0;k<count;k++){ 
        for(int j=0;j<tree_matrix.ncol();j++){
          arma::vec armacol=tree_matrix2.col(j);
          arma::uvec found=arma::find(armacol==internal_nodes_prop[k]);      
          if(found.size()>0){        
            //int_nodes_index[index_count]=j;
            //index_count++;
            propsplit*=alpha*pow((j+1),-beta) ; 
            break;
          }        
        }
        //int_nodes_index.resize(index_count);
        //if(int_nodes_index.size()!=0){      
        //  d=unique(as<IntegerVector>(wrap(int_nodes_index)));
        //  double d1=d[0];
        //  propsplit*=alpha*pow((d1+1),-beta) ;  
        //}
        //std::vector<int> temp(col);
        //int_nodes_index=temp;
        //index_count=0;
      } 
      
      for(int k=0;k<count_term_nodes;k++){ 
        for(int j=0;j<tree_matrix.ncol();j++){
          arma::vec armacol=tree_matrix2.col(j);
          arma::uvec found=arma::find(armacol==terminal_nodes_prop[k]);      
          if(found.size()>0){        
            //term_nodes_index[index_count2]=j;
            //index_count2++;
            propsplit*=1-alpha*pow((j+1),-beta) ;  
            
            break;
          }        
        }
        //term_nodes_index.resize(index_count2);
        //if(term_nodes_index.size()!=0){      
        //  d2=unique(as<IntegerVector>(wrap(term_nodes_index)));
        //  double d1=d2[0];
        //  propsplit*=1-alpha*pow((d1+1),-beta) ;  
        //}
        //std::vector<int> temp2(col);
        //term_nodes_index=temp2;
        //index_count2=0;
      }
    }
    
    return(propsplit);
  }
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericMatrix start_tree(double start_mean,double start_sd){
  
  NumericMatrix treemat(1,7);
  double rand=R::rnorm(start_mean,start_sd);
  NumericVector testrow = NumericVector::create(0,0,0,0,-1,rand,0);
  for(int k=0;k<1;k++){
    for(int j=0;j<7;j++){
      treemat(k,j)=testrow[j];
    }
  }
  List dimnms = // two vec. with static names
    List::create(CharacterVector::create("1"),
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  treemat.attr("dimnames") = dimnms;
  
  return(treemat);
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericMatrix start_tree2(){
  
  NumericMatrix treemat(1,7);
  //double rand=R::rnorm(start_mean,start_sd);
  NumericVector testrow = NumericVector::create(0,0,0,0,-1,0,0);
  //for(int k=0;k<1;k++){
  for(int j=0;j<7;j++){
    //treemat(k,j)=testrow[j];
    treemat(0,j)=testrow[j];
  }
  //}
  List dimnms = // two vec. with static names
    List::create(CharacterVector::create("1"),
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  treemat.attr("dimnames") = dimnms;
  
  return(treemat);
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericMatrix start_matrix(int n){
  NumericMatrix mat(n,1);
  std::fill(mat.begin(), mat.end(), 1);
  return(mat);
}
//######################################################################################################################//

// [[Rcpp::export]]
List evaluate_model_occams_window(NumericVector tree_lik,double lowest_BIC,double c,List tree_list,List tree_mat_list,IntegerVector tree_parent){
  IntegerVector sorted_lik_index=order_(tree_lik);
  std::vector<double> to_be_removed(tree_lik.size());
  int s=0;
  
  // check tree is in Occam's window if it isn't then delete it from list                    
  while((tree_lik[sorted_lik_index[s]-1])-(lowest_BIC) > c){
    
    // if(s==(tree_lik.size()-1)){
    //  break;
    // }
    
    //delete tree from tree list
    //if((tree_lik[sorted_lik_index[s]-1])-(lowest_BIC)>c){
    //set indicator for index of trees to be removed 
    to_be_removed[s]=sorted_lik_index[s]-1;
    
    if(s==(tree_lik.size()-1)){
      s+=1;
      break;
    }
    
    s+=1;
    //}
  }
  
  to_be_removed.resize(s);
  IntegerVector remove_order_index(to_be_removed.size());
  //delete elements from the higest index down 
  remove_order_index=order_(wrap(to_be_removed));
  
  for(int j=0;j<s;j++){    
    tree_list.erase(to_be_removed[remove_order_index[j]-1]);
    tree_mat_list.erase(to_be_removed[remove_order_index[j]-1]);
    tree_lik.erase(to_be_removed[remove_order_index[j]-1]);
    tree_parent.erase(to_be_removed[remove_order_index[j]-1]);
  }
  List ret(4);
  ret(0)=tree_lik;
  ret(1)=tree_list;
  ret(2)=tree_mat_list;
  ret(3)=tree_parent;
  
  return(ret);  
}

//######################################################################################################################//

// [[Rcpp::export]]
List evaluate_model_occams_window_exact(NumericVector tree_lik,double lowest_BIC,double c,List tree_list,List tree_mat_list,IntegerVector tree_parent, 
                                        List tree_pred_list){
  // Rcout << "Line 847.\n";
  // Rcout << "tree_lik =" << tree_lik << ".\n";
  
  IntegerVector sorted_lik_index=order_(tree_lik);
  std::vector<double> to_be_removed(tree_lik.size());
  int s=0;
  
  
  
  // Rcout << "Line 852.\n";
  // Rcout << "lowest_BIC =" << lowest_BIC << ".\n";
  // Rcout << "sorted_lik_index =" << sorted_lik_index << ".\n";
  // check tree is in Occam's window if it isn't then delete it from list                    
  while((tree_lik[sorted_lik_index[s]-1])-(lowest_BIC) > c){
    // Rcout << "Line 856.\n";
    
    // if(s==(tree_lik.size()-1)){
    //   // Rcout << "Line 858.\n";
    //   break;
    // }
    
    //delete tree from tree list
    //if((tree_lik[sorted_lik_index[s]-1])-(lowest_BIC)>c){
    //set indicator for index of trees to be removed 
    // Rcout << "Line 864.\n";
    to_be_removed[s]=sorted_lik_index[s]-1;
    
    if(s==(tree_lik.size()-1)){
      s+=1;
      break;
    }
    
    // Rcout << "Line 866.\n";
    s+=1;
    //}
  }
  // Rcout << "Line 870.\n";
  
  to_be_removed.resize(s);
  IntegerVector remove_order_index(to_be_removed.size());
  //delete elements from the higest index down 
  remove_order_index=order_(wrap(to_be_removed));
  // Rcout << "Line 872.\n";
  
  for(int j=0;j<s;j++){    
    tree_list.erase(to_be_removed[remove_order_index[j]-1]);
    tree_mat_list.erase(to_be_removed[remove_order_index[j]-1]);
    tree_lik.erase(to_be_removed[remove_order_index[j]-1]);
    tree_parent.erase(to_be_removed[remove_order_index[j]-1]);
    tree_pred_list.erase(to_be_removed[remove_order_index[j]-1]);
    
  }
  // Rcout << "Line 882.\n";
  
  List ret(5);
  ret(0)=tree_lik;
  ret(1)=tree_list;
  ret(2)=tree_mat_list;
  ret(3)=tree_parent;
  ret(4)=tree_pred_list;
  
  return(ret);  
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector get_testdata_term_obs(NumericMatrix test_data,NumericMatrix tree_data//,NumericVector term_node_means
) {
  //Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.
  
  arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
  arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);
  NumericVector terminal_nodes=find_term_nodes(tree_data);
  //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
  //NumericVector tree_predictions;
  //for each internal node find the observations that belong to the terminal nodes
  NumericVector predictions(test_data.nrow());
  if(terminal_nodes.size()==1){
    double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
    predictions=rep(nodemean,test_data.nrow());
  }
  else{
    for(int i=0;i<terminal_nodes.size();i++){
      //arma::mat subdata=testd;
      int curr_term=terminal_nodes[i];
      int row_index;
      int term_node=terminal_nodes[i];
      if(curr_term % 2==0){
        row_index=terminal_nodes[i];
      }else{
        row_index=terminal_nodes[i]-1;
      }
      //save the left and right node data into arma uvec
      
      arma::vec left_nodes=arma_tree.col(0);
      arma::vec right_nodes=arma_tree.col(1);
      arma::mat node_split_mat;    
      node_split_mat.set_size(0,3);
      
      while(row_index!=1){
        //for each terminal node work backwards and see if the parent node was a left or right node
        //append split info to a matrix 
        int rd=0;
        arma::uvec parent_node=arma::find(left_nodes == term_node);
        
        if(parent_node.size()==0){
          parent_node=arma::find(right_nodes == term_node);
          rd=1;
        }      
        node_split_mat.insert_rows(0,1);
        node_split_mat(0,0)=tree_data(parent_node[0],2);
        node_split_mat(0,1)=tree_data(parent_node[0],3);
        node_split_mat(0,2)=rd;
        row_index=parent_node[0] +1;
        term_node=parent_node[0]+1;
      }  
      //fill in the predicted value for tree
      arma::uvec pred_indices;
      int split= node_split_mat(0,0)-1;
      arma::vec tempvec = testd.col(split);
      double temp_split = node_split_mat(0,1);
      if(node_split_mat(0,2)==0){
        pred_indices = arma::find(tempvec <= temp_split);
      }else{
        pred_indices = arma::find(tempvec > temp_split);
      }
      
      arma::uvec temp_pred_indices;
      //arma::vec data_subset = testd.col(split);
      //data_subset=data_subset.elem(pred_indices);
      int n=node_split_mat.n_rows;
      
      for(int j=1;j<n;j++){
        int curr_sv=node_split_mat(j,0);
        double split_p = node_split_mat(j,1);
        
        arma::vec data_subset = testd.col(curr_sv-1);
        data_subset=data_subset.elem(pred_indices);
        
        if(node_split_mat(j,2)==0){
          //split is to the left
          temp_pred_indices=arma::find(data_subset <= split_p);
        }else{
          //split is to the right
          temp_pred_indices=arma::find(data_subset > split_p);
        }
        pred_indices=pred_indices.elem(temp_pred_indices);
        
        if(pred_indices.size()==0){
          continue;
        }
      }
      double nodemean=tree_data(terminal_nodes[i]-1,5);
      IntegerVector predind=as<IntegerVector>(wrap(arma::conv_to<arma::ivec>::from(pred_indices)));
      predictions[predind]= nodemean;
    } 
  }
  return(predictions);
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List get_initial_resids(NumericMatrix test_data,List List_of_lists_tree_tables,NumericVector ytrain
) {
  //Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.
  
  List List_of_resid_lists(List_of_lists_tree_tables.size());
  List new_pred_list(List_of_lists_tree_tables.size());
  
  for(int k=0;k<List_of_lists_tree_tables.size();k++){
    
    List One_sum_of_tree_list = List_of_lists_tree_tables[k];  
    //create list of outcomes from which to take away predictions in each round to create residuals
    NumericMatrix new_pred_mat(test_data.nrow(),One_sum_of_tree_list.size());
    
    List temp_resid_list(One_sum_of_tree_list.size()) ; 
    for(int q=0;q<One_sum_of_tree_list.size();q++){
      temp_resid_list[q]=ytrain;
    }
    
    
    for(int l=0;l<One_sum_of_tree_list.size();l++){
      
      
      NumericMatrix tree_data = One_sum_of_tree_list[l];
      
      arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
      arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);
      NumericVector terminal_nodes=find_term_nodes(tree_data);
      //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
      //NumericVector tree_predictions;
      //for each internal node find the observations that belong to the terminal nodes
      NumericVector predictions(test_data.nrow());
      if(terminal_nodes.size()==1){
        double nodemean=mean(ytrain)/One_sum_of_tree_list.size();				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
        predictions=rep(nodemean,test_data.nrow());
      }
      else{
        for(int i=0;i<terminal_nodes.size();i++){
          //arma::mat subdata=testd;
          int curr_term=terminal_nodes[i];
          int row_index;
          int term_node=terminal_nodes[i];
          if(curr_term % 2==0){
            row_index=terminal_nodes[i];
          }else{
            row_index=terminal_nodes[i]-1;
          }
          //save the left and right node data into arma uvec
          
          arma::vec left_nodes=arma_tree.col(0);
          arma::vec right_nodes=arma_tree.col(1);
          arma::mat node_split_mat;    
          node_split_mat.set_size(0,3);
          
          while(row_index!=1){
            //for each terminal node work backwards and see if the parent node was a left or right node
            //append split info to a matrix 
            int rd=0;
            arma::uvec parent_node=arma::find(left_nodes == term_node);
            
            if(parent_node.size()==0){
              parent_node=arma::find(right_nodes == term_node);
              rd=1;
            }      
            node_split_mat.insert_rows(0,1);
            node_split_mat(0,0)=tree_data(parent_node[0],2);
            node_split_mat(0,1)=tree_data(parent_node[0],3);
            node_split_mat(0,2)=rd;
            row_index=parent_node[0] +1;
            term_node=parent_node[0]+1;
          }  
          //fill in the predicted value for tree
          arma::uvec pred_indices;
          int split= node_split_mat(0,0)-1;
          arma::vec tempvec = testd.col(split);
          double temp_split = node_split_mat(0,1);
          if(node_split_mat(0,2)==0){
            pred_indices = arma::find(tempvec <= temp_split);
          }else{
            pred_indices = arma::find(tempvec > temp_split);
          }
          
          arma::uvec temp_pred_indices;
          //arma::vec data_subset = testd.col(split);
          //data_subset=data_subset.elem(pred_indices);
          int n=node_split_mat.n_rows;
          
          for(int j=1;j<n;j++){
            int curr_sv=node_split_mat(j,0);
            double split_p = node_split_mat(j,1);
            
            arma::vec data_subset = testd.col(curr_sv-1);
            data_subset=data_subset.elem(pred_indices);
            
            if(node_split_mat(j,2)==0){
              //split is to the left
              temp_pred_indices=arma::find(data_subset <= split_p);
            }else{
              //split is to the right
              temp_pred_indices=arma::find(data_subset > split_p);
            }
            pred_indices=pred_indices.elem(temp_pred_indices);
            
            if(pred_indices.size()==0){
              continue;
            }
          }
          IntegerVector predind=as<IntegerVector>(wrap(arma::conv_to<arma::ivec>::from(pred_indices)));
          NumericVector ys_in_node= ytrain[predind];
          double nodemean=mean(ys_in_node)/One_sum_of_tree_list.size();
          predictions[predind]= nodemean;
        } 
      }
      
      //Add predictions to tree predictions
      for(int q=0;q<One_sum_of_tree_list.size();q++){
        if(q==l){
          
        }else{
          NumericVector temp_for_update = temp_resid_list[q];
          temp_resid_list[q]=temp_for_update-predictions;
        }
      }
      
      new_pred_mat(_,l)=predictions;
      
    }
    List_of_resid_lists[k]=temp_resid_list;
    new_pred_list[k]=new_pred_mat;
  }
  
  List ret(2);
  ret[0]=List_of_resid_lists;
  ret[1]=new_pred_list;
  
  return(ret);
  
}
//######################################################################################################################//

// [[Rcpp::export]]
List resize(const List& x, int n ){
  List y(n) ;
  for( int i=0; i<n; i++) y[i] = x[i] ;
  
  return y ;
}
//######################################################################################################################//

// [[Rcpp::export]]
List resize_bigger( const List& x, int n ){
  int oldsize = x.size() ;
  List y(n) ;
  for( int i=0; i<oldsize; i++) y[i] = x[i] ;
  return y ;
}
//###########################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat J(NumericMatrix obs_to_nodes_temp,NumericVector tree_term_nodes){
  //this function will make a binary nxb matrix where each column assigns observations to terminal nodes
  //Rcout << "Line 1218.\n";
  
  //create J matrix with correct dimensions and fill with zeros
  arma::mat Jmat(obs_to_nodes_temp.nrow(), tree_term_nodes.size());
  Jmat.zeros();
  
  //for each terminal node get the observations associated with it and set column
  for(int i=0;i<tree_term_nodes.size();i++){
    // Rcout << "Line 1153.\n";
    
    double tn=tree_term_nodes[i];
    //Rcout << "Line 1229.\n";
    
    arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
    // Rcout << "Line 1232.\n";
    // 
    // //assign term_obs to the correct index of J
    // IntegerVector term_obs2=as<IntegerVector>(wrap(term_obs));
    //  Rcout << "Line 1236.\n";
    // 
    // NumericVector obs_col(obs_to_nodes_temp.nrow());
    // // Rcout << "Line 1166.\n";
    // obs_col[term_obs2]=1;
    // // Rcout << "Line 1168.\n";
    // arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
    // // Rcout << "Line 1170.\n";
    // Jmat.col(i)= colmat;
    
    arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
    colmat.elem(term_obs).fill(1);
    Jmat.col(i)= colmat;
  }
  //Rcout << "Line 1250.\n";
  
  return(Jmat);
}

//###########################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector mu_vector(List sum_treetable,int n){
  NumericVector mu_vec;
  
  for(int j=0;j<sum_treetable.size();j++){    
    NumericMatrix curr_tree=sum_treetable[j];
    NumericVector tree_term_nodes=find_term_nodes(curr_tree);
    NumericVector term_means1=curr_tree(_,5);
    NumericVector term_means=remove_zero(term_means1);
    
    for(int i=0;i<term_means.size();i++){
      mu_vec.push_back(term_means[i]);  
    }
  }
  return(mu_vec);
}
//###########################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat W(List sum_treetable ,List sum_obs_to_nodes,int n){
  //this will take in a list of obs to node matrices for each tree in the sum make the J matrix assigning observations to terminal nodes
  //J is an nxb_j matrix. It will then iteratively append the J matrix to itself to get the overall W matrix which has dimensions nxsumb_j
  
  //create empty matrix to which we will append individual J matrices
  arma::mat W(n,0);
  int upsilon=0;
  for(int j=0;j<sum_obs_to_nodes.size();j++){
    
    NumericMatrix curr_tree=sum_treetable[j];
    NumericMatrix curr_obs_nodes=sum_obs_to_nodes[j];
    //Rcout << "Line 1288.\n"; 
    NumericVector tree_term_nodes=find_term_nodes(curr_tree);
    //Rcout << "Line 1290.\n"; 
    
    int b_j=tree_term_nodes.size();
    //will make J as we go in BART-BMA no need to create it again here....
    arma::mat Jmat=J(curr_obs_nodes,tree_term_nodes);
    //Rcout << "Line 1295.\n"; 
    
    W.insert_cols(upsilon,Jmat);
    upsilon+=b_j;
  }
  
  return(W);
  // ret[1]=mu_vec;
}
//######################################################################################################################//
// 
// #include <RcppDist.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double likelihood_function_new(NumericVector y_temp,NumericMatrix treetable_temp,NumericMatrix obs_to_nodes_temp,double a,double mu,double nu,double lambda){
//   
//   int n=y_temp.size();
//   NumericVector tree_term_nodes=find_term_nodes(treetable_temp);
//   //int b_j=tree_term_nodes.size();
//   //will make J as we go in BART-BMA no need to create it again here....
//   arma::mat Wmat=J(obs_to_nodes_temp,tree_term_nodes);
// 
//   arma::vec yvec=Rcpp::as<arma::vec>(y_temp);
//   arma::mat y(n,1);
//   y.col(0)=yvec;
//   
// 
//   arma::mat WWt=(1/a)*(Wmat*Wmat.t());
//   arma::mat In(n,n);
//   arma::mat covar=lambda*(In.eye()+WWt);
//   arma::vec  muvec=arma::zeros(y.n_rows) ;
//   
// 
//   arma::vec rel= dmvt(y.t(), muvec,
//                       covar, nu, 
//                       true); 
//   
// 
//   double rel2=as<double>(wrap(rel));
//   
//   // double tree_log_lik;
//   // NumericVector terminal_nodes=find_term_nodes(treetable_temp);
//   // double b=terminal_nodes.size();
//   // IntegerVector n(b);
//   // double term1=0;
//   // double term2=0;
//   // arma::uvec term_obs;
//   // int ni;
//   // arma::vec y_temp2=as<arma::vec>(y_temp);
//   // for(int i=0;i< b;i++){
//   //   term_obs=find_term_obs(obs_to_nodes_temp,terminal_nodes[i]);
//   //   arma::vec y_k=y_temp2.elem(term_obs);
//   //   ni=term_obs.size();
//   //   n[i]=ni;
//   //   double ybar=0;
//   //   if(y_k.size()!=0){
//   //     ybar=mean(y_k);
//   //   }
//   //   term1+=log(ni+a);
//   //   arma::vec y_k_sq=pow(y_k,2);
//   //   double sum_yksq=sum(y_k_sq);
//   //   double b2=pow(ni*ybar +a*mu,2)/(ni+a);
//   //   term2+=(sum_yksq+a*pow(mu,2)-b2+nu*lambda);
//   // }
//   // tree_log_lik=(b/2)*log(a)-0.5*term1-((y_temp.size()+nu)/2)*log(term2);
//   // for(int i=0;i<b;i++){
//   //   if(n[i]<=5){
//   //     tree_log_lik=tree_log_lik-(100000);
//   //   }
//   //   else{
//   //     tree_log_lik=tree_log_lik;
//   //   }
//   // }
//   return(rel2);
// }
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double likelihood_function2(NumericVector y_temp,NumericMatrix treetable_temp,NumericMatrix obs_to_nodes_temp,double a,double mu,double nu,double lambda){
  
  int n=y_temp.size();
  NumericVector tree_term_nodes=find_term_nodes(treetable_temp);
  //int b_j=tree_term_nodes.size();
  //will make J as we go in BART-BMA no need to create it again here....
  arma::mat Wmat=J(obs_to_nodes_temp,tree_term_nodes);
  double b=Wmat.n_cols;
  
  arma::vec yvec=Rcpp::as<arma::vec>(y_temp);
  arma::mat y(n,1);
  y.col(0)=yvec;
  
  
  //get exponent
  double expon=(n+nu)*0.5;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;
  
  //get t(y)inv(psi)J
  arma::mat ytW=y.t()*Wmat;
  //get t(J)inv(psi)J  
  arma::mat WtW=Wmat.t()*Wmat;
  //get jpsij +aI
  arma::mat aI(b,b);
  aI=a*aI.eye();
  arma::mat sec_term=WtW+aI;
  //arma::mat sec_term_inv=sec_term.i();
  arma::mat sec_term_inv=inv_sympd(sec_term);
  //get t(J)inv(psi)y
  arma::mat third_term=Wmat.t()*y;
  //get m^TV^{-1}m
  arma::mat mvm= ytW*sec_term_inv*third_term;
  
  
  //arma::mat rel=(b/2)*log(a)-(1/2)*log(det(sec_term))-expon*log(nu*lambda - mvm +yty);
  arma::mat rel=(b*0.5)*log(a)-0.5*real(arma::log_det(sec_term))-expon*log(nu*lambda - mvm +yty);
  
  
  double rel2=as<double>(wrap(rel));
  return(rel2);
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List likelihood_function2_exact(NumericVector y_temp,NumericMatrix treetable_temp,NumericMatrix obs_to_nodes_temp,double a,double mu,double nu,double lambda){
  
  //double likelihood_function2(NumericVector y_temp,NumericMatrix treetable_temp,NumericMatrix obs_to_nodes_temp,double a,double mu,double nu,double lambda){
  
  // Rcout << "Line 1420.\n";
  
  int n=y_temp.size();
  // Rcout << "Line 1328.\n";
  // Rcout << "treetable_temp.ncol() = " << treetable_temp.ncol() << ".\n";
  NumericVector tree_term_nodes=find_term_nodes(treetable_temp);
  // Rcout << "Line 1426.\n";
  
  //int b_j=tree_term_nodes.size();
  //will make J as we go in BART-BMA no need to create it again here....
  arma::mat Wmat=J(obs_to_nodes_temp,tree_term_nodes);
  double b=Wmat.n_cols;
  // Rcout << "Line 1432.\n";
  
  arma::vec yvec=Rcpp::as<arma::vec>(y_temp);
  arma::mat y(n,1);
  y.col(0)=yvec;
  // Rcout << "Line 1325";
  
  
  //get exponent
  double expon=(n+nu)*0.5;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;
  
  //get t(y)inv(psi)J
  arma::mat ytW=y.t()*Wmat;
  //get t(J)inv(psi)J  
  arma::mat WtW=Wmat.t()*Wmat;
  //get jpsij +aI
  arma::mat aI(b,b);
  aI=a*aI.eye();
  arma::mat sec_term=WtW+aI;
  //arma::mat sec_term_inv=sec_term.i();
  arma::mat sec_term_inv=inv_sympd(sec_term);
  //get t(J)inv(psi)y
  arma::mat third_term=Wmat.t()*y;
  //get m^TV^{-1}m
  arma::mat mvm= ytW*sec_term_inv*third_term;
  
  // Rcout << "Line 1461";
  
  //arma::mat rel=(b/2)*log(a)-(1/2)*log(det(sec_term))-expon*log(nu*lambda - mvm +yty);
  arma::mat rel=(b*0.5)*log(a)-0.5*real(arma::log_det(sec_term))-expon*log(nu*lambda - mvm +yty);
  
  // Rcout << "Line 1466";
  
  double rel2=as<double>(wrap(rel));
  
  //double rel2=arma::as_scalar(rel);
  // Rcout << "Line 1471";
  
  arma::vec predsoutput=Wmat*sec_term_inv*third_term;
  
  // Rcout << "Line 1372";
  
  List ret(2);
  ret[0]=rel2;
  ret[1]=wrap(predsoutput);
  
  
  //return(rel2);
  return(ret);
  
}
//############################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sumtree_likelihood_function(NumericVector y_temp,List sum_treetable ,List sum_obs_to_nodes,int n,double a,double nu,double lambda){
  //make W and mu matrices for the sum of trees
  arma::mat Wmat=W(sum_treetable,sum_obs_to_nodes,n);
  double b=Wmat.n_cols;
  arma::vec yvec=Rcpp::as<arma::vec>(y_temp);
  arma::mat y(n,1);
  y.col(0)=yvec;
  //get exponent
  double expon=(n+nu+b)*0.5;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;
  
  //get t(y)inv(psi)J
  arma::mat ytW=y.t()*Wmat;
  //get t(J)inv(psi)J  
  arma::mat WtW=Wmat.t()*Wmat;
  //get jpsij +aI
  arma::mat aI(b,b);
  aI=a*aI.eye();
  arma::mat sec_term=WtW+aI;
  arma::mat sec_term_inv=sec_term.i();
  //get t(J)inv(psi)y
  arma::mat third_term=Wmat.t()*y;
  //get m^TV^{-1}m
  arma::mat mvm= ytW*sec_term_inv*third_term;
  arma::mat rel=-expon*log(nu*lambda - mvm +yty);
  double rel2=as<double>(wrap(rel));
  return(rel2);
}  
//############################################################################################################################//
// 
// #include <RcppDist.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double sumtree_likelihood_function_new(NumericVector y_temp,List sum_treetable ,List sum_obs_to_nodes,int n,double a,double nu,double lambda){
//   //make W and mu matrices for the sum of trees
// 
//   arma::mat Wmat=W(sum_treetable,sum_obs_to_nodes,n);
// 
//   //double b=Wmat.n_cols;
//   arma::vec yvec=Rcpp::as<arma::vec>(y_temp);
//   arma::mat y(n,1);
//   y.col(0)=yvec;
//   //get exponent
//   //double expon=(n+nu+b)/2;
//   //get y^Tpsi^{-1}y
//   // arma::mat psi_inv=psi.i();
//   //arma::mat yty=y.t()*y;
// 
//   arma::mat WWt=(1/a)*(Wmat*Wmat.t());
//   arma::mat In(n,n);
//   arma::mat covar=lambda*(In.eye()+WWt);
//   
//   arma::vec  muvec=arma::zeros(y.n_rows) ;
// 
//   arma::vec rel= dmvt(y.t(), muvec,
//                       covar, nu, 
//                       true); 
// 
// 
//   double rel2=as<double>(wrap(rel));
//   
//   //get t(y)inv(psi)J
//   //arma::mat ytW=y.t()*Wmat;
//   //get t(J)inv(psi)J  
//   //arma::mat WtW=Wmat.t()*Wmat;
//   //get jpsij +aI
//   //arma::mat aI(b,b);
//   //aI=a*aI.eye();
//   //arma::mat sec_term=WtW+aI;
//   //arma::mat sec_term_inv=sec_term.i();
//   //get t(J)inv(psi)y
//   //arma::mat third_term=Wmat.t()*y;
//   //get m^TV^{-1}m
//   //arma::mat mvm= ytW*sec_term_inv*third_term;
//   //arma::mat rel=-expon*log(nu*lambda - mvm +yty);
//   //double rel2=as<double>(wrap(rel));
//   
//   return(rel2);
// }
//############################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sumtree_likelihood_function2(NumericVector y_temp,List sum_treetable ,List sum_obs_to_nodes,int n,double a,double nu,double lambda){
  //make W and mu matrices for the sum of trees
  arma::mat Wmat=W(sum_treetable,sum_obs_to_nodes,n);
  double b=Wmat.n_cols;
  arma::vec yvec=Rcpp::as<arma::vec>(y_temp);
  arma::mat y(n,1);
  y.col(0)=yvec;
  //get exponent
  double expon=(n+nu)*0.5;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;
  
  //get t(y)inv(psi)J
  arma::mat ytW=y.t()*Wmat;
  //get t(J)inv(psi)J  
  arma::mat WtW=Wmat.t()*Wmat;
  //get jpsij +aI
  arma::mat aI(b,b);
  aI=a*aI.eye();
  arma::mat sec_term=WtW+aI;
  //arma::mat sec_term_inv=sec_term.i();
  arma::mat sec_term_inv=inv_sympd(sec_term);  
  //get t(J)inv(psi)y
  arma::mat third_term=Wmat.t()*y;
  //get m^TV^{-1}m
  arma::mat mvm= ytW*sec_term_inv*third_term;
  //arma::mat rel=(b/2)*log(a)-(1/2)*log(det(sec_term))-expon*log(nu*lambda - mvm +yty);
  arma::mat rel=(b/2)*log(a)-0.5*log(det(sec_term))-expon*log(nu*lambda - mvm +yty);
  
  double rel2=as<double>(wrap(rel));
  return(rel2);
}  
//############################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sumtree_likelihood_function2_exact(NumericVector y_temp,List sum_treetable ,List sum_obs_to_nodes,int n,double a,double nu,double lambda){
  
  //double sumtree_likelihood_function2(NumericVector y_temp,List sum_treetable ,List sum_obs_to_nodes,int n,double a,double nu,double lambda){
  //make W and mu matrices for the sum of trees
  
  // Rcout<< "Line 1617 .\n";
  
  arma::mat Wmat=W(sum_treetable,sum_obs_to_nodes,n);
  // Rcout << "Line 1620.\n";
  
  double b=Wmat.n_cols;
  arma::vec yvec=Rcpp::as<arma::vec>(y_temp);
  arma::mat y(n,1);
  y.col(0)=yvec;
  //get exponent
  double expon=(n+nu)*0.5;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;
  
  //get t(y)inv(psi)J
  arma::mat ytW=y.t()*Wmat;
  //get t(J)inv(psi)J  
  arma::mat WtW=Wmat.t()*Wmat;
  //get jpsij +aI
  arma::mat aI(b,b);
  aI=a*aI.eye();
  arma::mat sec_term=WtW+aI;
  //arma::mat sec_term_inv=sec_term.i();
  // Rcout<< "Line 1492 .\n";
  
  arma::mat sec_term_inv=inv_sympd(sec_term);  
  //get t(J)inv(psi)y
  // Rcout<< "Line 1496 .\n";
  
  arma::mat third_term=Wmat.t()*y;
  //get m^TV^{-1}m
  arma::mat mvm= ytW*sec_term_inv*third_term;
  //arma::mat rel=(b/2)*log(a)-(1/2)*log(det(sec_term))-expon*log(nu*lambda - mvm +yty);
  // Rcout<< "Line 1659 .\n";
  
  arma::mat rel=(b/2)*log(a)-0.5*real(arma::log_det(sec_term))-expon*log(nu*lambda - mvm +yty);
  // Rcout<< "Line 1652 .\n";
  
  double rel2=as<double>(wrap(rel));
  // Rcout<< "Line 1655 .\n";
  
  
  arma::vec predsoutput=Wmat*sec_term_inv*third_term;
  
  List ret(2);
  ret[0]=rel2;
  ret[1]=wrap(predsoutput);
  
  
  //return(rel2);
  return(ret);
}  
//############################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sumtree_likelihood_function3(NumericVector y_temp,List sum_treetable ,List sum_obs_to_nodes,int n,double a,double nu,double lambda){
  
  //Cholesky attempt
  
  //make W and mu matrices for the sum of trees
  arma::mat Wmat=W(sum_treetable,sum_obs_to_nodes,n);
  double b=Wmat.n_cols;
  arma::vec yvec=Rcpp::as<arma::vec>(y_temp);
  arma::mat y(n,1);
  y.col(0)=yvec;
  //get exponent
  double expon=(n+nu)*0.5;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;
  
  //get t(y)inv(psi)J
  arma::mat ytW=y.t()*Wmat;
  //get t(J)inv(psi)J  
  arma::mat WtW=Wmat.t()*Wmat;
  //get jpsij +aI
  arma::mat aI(b,b);
  aI=a*aI.eye();
  arma::mat sec_term=WtW+aI;
  
  //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
  //obtain the log of the root of the determinant
  double rootisum = arma::sum(log(rooti.diag()));
  
  
  //arma::mat sec_term_inv=sec_term.i();
  //get t(J)inv(psi)y
  //arma::mat third_term=Wmat.t()*y;
  //get m^TV^{-1}m
  //arma::mat mvm= ytW*sec_term_inv*third_term;
  
  //t(:)*W^T*Y
  arma::mat LtWtY= rooti*(Wmat.t()*y);
  
  
  arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
  double rel2=as<double>(wrap(rel));
  return(rel2);
} 
//############################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sumtree_likelihood_function4(NumericVector y_temp,List sum_treetable ,List sum_obs_to_nodes,int n,double a,double nu,double lambda){
  
  //Eigenvalue attempt
  
  //make W and mu matrices for the sum of trees
  arma::mat Wmat=W(sum_treetable,sum_obs_to_nodes,n);
  double b=Wmat.n_cols;
  arma::vec yvec=Rcpp::as<arma::vec>(y_temp);
  arma::mat y(n,1);
  y.col(0)=yvec;
  //get exponent
  double expon=(n+nu)*0.5;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;
  
  //get t(y)inv(psi)J
  arma::mat ytW=y.t()*Wmat;
  //get t(J)inv(psi)J  
  arma::mat WtW=Wmat.t()*Wmat;
  //get jpsij +aI
  arma::mat aI(b,b);
  aI=a*aI.eye();
  arma::mat sec_term=WtW+aI;
  
  
  
  arma::mat sec_term_inv=sec_term.i();
  //get t(J)inv(psi)y
  arma::mat third_term=Wmat.t()*y;
  //get m^TV^{-1}m
  arma::mat mvm= ytW*sec_term_inv*third_term;
  
  double logdet = sum(arma::log(arma::eig_sym(sec_term)));
  
  
  //arma::mat rel=(b/2)*log(a)-(1/2)*logdet-expon*log(nu*lambda - mvm +yty);
  arma::mat rel=(b/2)*log(a)-0.5*logdet-expon*log(nu*lambda - mvm +yty);
  
  double rel2=as<double>(wrap(rel));
  return(rel2);
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                    NumericVector resids,arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,
                    double a,double mu,double nu,double lambda,double c,double lowest_BIC,int parent
                      ,NumericMatrix cp_mat,double alpha,double beta,int maxOWsize, 
                      unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split//,int first_round
){
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;
  NumericMatrix treetable_c=treetable;
  NumericMatrix treemat_c=tree_mat;
  
  NumericVector terminal_nodes=find_term_nodes(treetable_c);
  //IntegerVector change_node1;
  int list_size=1000;
  std::vector<double> tree_lik(list_size);
  List proposal_tree;
  //List ret(9);
  List ret(4);
  
  bool no_tree_err=0;
  //List likeliest_tree;
  List tree_list(list_size);
  List tree_mat_list(list_size);
  int count=0;
  //std::vector<int> tree_parent(list_size);
  //int best_sv;
  //double best_sp;
  double tree_prior=0;
  //List changetree;
  double BIC;
  //int p;
  List eval_model;
  //NumericVector int_nodes;
  //arma::colvec curr_col=data.col(0);
  //arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[0]);
  //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[0]));
  //arma::mat data_curr_node=data.rows(grow_obs);
  //double d=d1[0];
  //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),cp_mat(0,0)+1);
  double lik;
  
  for(int l=0;l<terminal_nodes.size();l++){
    //loop over each terminal node
    arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[l]);
    //depth of tree at current terminal node
    //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[l]));
    arma::mat data_curr_node=data.rows(grow_obs);
    //double d=d1[0];
    int d = find_term_cols(treemat_c,terminal_nodes[l]);
    int w=cp_mat.nrow();
    if(data_curr_node.n_rows<=min_num_obs_for_split){
      throw std::range_error("not enough obs in node to grow any further");
      //continue;
    }
    for(int k=0;k<w;k++){
      split_var=cp_mat(k,0)+1;
      //arma::colvec curr_cols=data.col(split_var-1);
      //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),split_var);
      // The following lines are unnecessary because get_min.size()=data_curr_node.n_rows above
      // if(get_min.size()<=2){
      //   throw std::range_error("obs in this terminal node are too small");
      // }
      
      double split_point=cp_mat(k,1);
      arma::vec curr_cols2=data_curr_node.col(split_var-1);
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));
      
      if(ld_prop.size()<=min_num_obs_after_split || rd_prop.size()<=min_num_obs_after_split){
        continue;
      }
      proposal_tree=grow_tree(data,//resids,
                              treemat_c,terminal_nodes[l],treetable_c,
                              split_var,split_point,//terminal_nodes,
                              wrap(arma::conv_to<arma::vec>::from(grow_obs)),
                              d//,get_min,data_curr_node
      );
      
      // Test lines below have been removed
      // NumericMatrix test =proposal_tree[0];
      // NumericMatrix test1 =proposal_tree[1];
      
      // if(test1.ncol()==3){
      //   NumericVector u1=unique(test1(_,0));
      //   NumericVector u2=unique(test1(_,1));
      //   NumericVector u3=unique(test1(_,2));
      // }
      
      // get_best_split should only be used in the first_round. Removing the if condition below 
      //if(first_round==1){
      lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);
      
      // }else{
      //   //have a sum of trees
      //   lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);  
      // }
      tree_prior=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,proposal_tree[0],proposal_tree[1],alpha,beta);
      //int_nodes=find_term_nodes(proposal_tree[0]);
      //p=int_nodes.size();
      //BIC=-2*(lik+log(tree_prior))+p*log(data.n_rows); 
      BIC=-2*(lik+log(tree_prior));  
      
      if(BIC<lowest_BIC){
        lowest_BIC=BIC;
        //best_sv=split_var;
        //best_sp=split_point;
        //likeliest_tree=proposal_tree;
        tree_list[count]=proposal_tree[0];
        tree_mat_list[count]=proposal_tree[1];
        tree_lik[count]=BIC;
        //tree_parent[count]=parent;
        count++;
        if(count==(tree_list.size()-1)){
          list_size=list_size*2;
          tree_list=resize_bigger(tree_list,list_size);
          tree_mat_list=resize_bigger(tree_mat_list,list_size);
          tree_lik.resize(list_size);
          //tree_parent.resize(list_size);
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){
          // if(is<NumericMatrix>(proposal_tree[0])){
          // }else{
          //   throw std::range_error("proposal tree not a matrix");
          // }
          tree_list[count]=proposal_tree[0];
          tree_mat_list[count]=proposal_tree[1];
          tree_lik[count]=BIC;
          //tree_parent[count]=parent;
          count++;
          if(count==(tree_list.size()-1)){
            list_size=list_size*2;
            tree_list=resize_bigger(tree_list,list_size);
            tree_mat_list=resize_bigger(tree_mat_list,list_size);
            tree_lik.resize(list_size);
            //tree_parent.resize(list_size);						  
          }
        }
      }
    }  
  }
  tree_list=resize(tree_list,count);
  tree_mat_list=resize(tree_mat_list,count);
  tree_lik.resize(count);
  //tree_parent.resize(count);
  IntegerVector tree_parent(count, parent);
  
  
  if(count>0){
    
    if(less_greedy==1){
      ret[0]=tree_list;
      ret[1]=tree_lik;
      ret[2]=tree_mat_list;
      ret[3]=tree_parent;
      
      return (ret);	
    }
    
    //eval_model=evaluate_model_occams_window(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));
    eval_model=evaluate_model_occams_window(wrap(tree_lik),lowest_BIC, c,wrap(tree_list),wrap(tree_mat_list),tree_parent);
    NumericVector testlik =eval_model[0];
    List testtree =eval_model[1];    
    List testmat =eval_model[2]; 
    IntegerVector testpar =eval_model[3];
    
    if(testlik.size()>0){
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){
        IntegerVector owindices=orderforOW(testlik);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        for(int t=0;t<maxOWsize;t++){
          temp_olik[t]=testlik[owindices[t]];
          temp_otrees[t]=testtree[owindices[t]];
          temp_omat[t]=testmat[owindices[t]];
          temp_oparent[t]=testpar[owindices[t]];
        }
        testlik=temp_olik;
        testtree=temp_otrees;
        testmat=temp_omat;
        testpar=temp_oparent;
      }
      // ret[0]=lowest_BIC;
      // ret[1]=best_sv;
      // ret[2]=best_sp;
      // ret[3]=likeliest_tree;
      // ret[4]=testtree;
      // ret[5]=testlik;
      // ret[6]=testmat;
      // ret[7]=testpar;
      // ret[8]=no_tree_err;
      
      
      ret[0]=testtree;
      ret[1]=testlik;
      ret[2]=testmat;
      ret[3]=testpar;
      
      return (ret);
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;
      List gr(1);
      gr[0]=no_tree_err;
      return(gr);
    }
  }else{
    no_tree_err=1;
    List gr(1);
    gr[0]=no_tree_err;
    return(gr);
  }
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_2(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                      NumericVector resids,arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,
                      double a,double mu,double nu,double lambda,double c,double lowest_BIC,int parent
                        ,List cp_matlist,double alpha,double beta,int maxOWsize, unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split//,int first_round
){
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;
  NumericMatrix treetable_c=treetable;
  NumericMatrix treemat_c=tree_mat;
  
  NumericVector terminal_nodes=find_term_nodes(treetable_c);
  //IntegerVector change_node1;
  int list_size=1000;
  std::vector<double> tree_lik(list_size);
  List proposal_tree;
  //List ret(9);
  List ret(4);
  
  bool no_tree_err=0;
  //List likeliest_tree;
  List tree_list(list_size);
  List tree_mat_list(list_size);
  int count=0;
  //std::vector<int> tree_parent(list_size);
  //int best_sv;
  //double best_sp;
  double tree_prior=0;
  //List changetree;
  double BIC;
  //int p;
  List eval_model;
  //NumericVector int_nodes;
  //arma::colvec curr_col=data.col(0);
  //arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[0]);
  //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[0]));
  //arma::mat data_curr_node=data.rows(grow_obs);
  //double d=d1[0];
  //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),cp_mat(0,0)+1);
  double lik;
  
  for(int l=0;l<terminal_nodes.size();l++){
    NumericMatrix cp_mat=cp_matlist[l];
    //loop over each terminal node
    arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[l]);
    //depth of tree at current terminal node
    //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[l]));
    arma::mat data_curr_node=data.rows(grow_obs);
    //double d=d1[0];
    int d = find_term_cols(treemat_c,terminal_nodes[l]);
    int w=cp_mat.nrow();
    if(data_curr_node.n_rows<=min_num_obs_for_split){
      throw std::range_error("not enough obs in node to grow any further");
      //continue;
    }
    for(int k=0;k<w;k++){
      split_var=cp_mat(k,0)+1;
      //arma::colvec curr_cols=data.col(split_var-1);
      //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),split_var);
      // The following lines are unnecessary because get_min.size()=data_curr_node.n_rows above
      // if(get_min.size()<=2){
      //   throw std::range_error("obs in this terminal node are too small");
      // }
      
      double split_point=cp_mat(k,1);
      arma::vec curr_cols2=data_curr_node.col(split_var-1);
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));
      
      if(ld_prop.size()<=min_num_obs_after_split || rd_prop.size()<=min_num_obs_after_split){
        continue;
      }
      proposal_tree=grow_tree(data,//resids,
                              treemat_c,terminal_nodes[l],treetable_c,
                              split_var,split_point,//terminal_nodes,
                              wrap(arma::conv_to<arma::vec>::from(grow_obs)),
                              d//,get_min,data_curr_node
      );
      
      // Test lines below have been removed
      // NumericMatrix test =proposal_tree[0];
      // NumericMatrix test1 =proposal_tree[1];
      
      // if(test1.ncol()==3){
      //   NumericVector u1=unique(test1(_,0));
      //   NumericVector u2=unique(test1(_,1));
      //   NumericVector u3=unique(test1(_,2));
      // }
      
      // get_best_split should only be used in the first_round. Removing the if condition below 
      //if(first_round==1){
      lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);
      
      // }else{
      //   //have a sum of trees
      //   lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);  
      // }
      tree_prior=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,proposal_tree[0],proposal_tree[1],alpha,beta);
      //int_nodes=find_term_nodes(proposal_tree[0]);
      //p=int_nodes.size();
      //BIC=-2*(lik+log(tree_prior))+p*log(data.n_rows); 
      BIC=-2*(lik+log(tree_prior));  
      
      if(BIC<lowest_BIC){
        lowest_BIC=BIC;
        //best_sv=split_var;
        //best_sp=split_point;
        //likeliest_tree=proposal_tree;
        tree_list[count]=proposal_tree[0];
        tree_mat_list[count]=proposal_tree[1];
        tree_lik[count]=BIC;
        //tree_parent[count]=parent;
        count++;
        if(count==(tree_list.size()-1)){
          list_size=list_size*2;
          tree_list=resize_bigger(tree_list,list_size);
          tree_mat_list=resize_bigger(tree_mat_list,list_size);
          tree_lik.resize(list_size);
          //tree_parent.resize(list_size);
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){
          // if(is<NumericMatrix>(proposal_tree[0])){
          // }else{
          //   throw std::range_error("proposal tree not a matrix");
          // }
          tree_list[count]=proposal_tree[0];
          tree_mat_list[count]=proposal_tree[1];
          tree_lik[count]=BIC;
          //tree_parent[count]=parent;
          count++;
          if(count==(tree_list.size()-1)){
            list_size=list_size*2;
            tree_list=resize_bigger(tree_list,list_size);
            tree_mat_list=resize_bigger(tree_mat_list,list_size);
            tree_lik.resize(list_size);
            //tree_parent.resize(list_size);						  
          }
        }
      }
    }  
  }
  tree_list=resize(tree_list,count);
  tree_mat_list=resize(tree_mat_list,count);
  tree_lik.resize(count);
  //tree_parent.resize(count);
  IntegerVector tree_parent(count, parent);
  if(count>0){
    if(less_greedy==1){
      ret[0]=tree_list;
      ret[1]=tree_lik;
      ret[2]=tree_mat_list;
      ret[3]=tree_parent;
      
      return (ret);	
    }
    //eval_model=evaluate_model_occams_window(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));
    eval_model=evaluate_model_occams_window(wrap(tree_lik),lowest_BIC,c,wrap(tree_list),wrap(tree_mat_list),tree_parent);
    NumericVector testlik =eval_model[0];
    List testtree =eval_model[1];    
    List testmat =eval_model[2]; 
    IntegerVector testpar =eval_model[3];
    
    if(testlik.size()>0){
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){
        IntegerVector owindices=orderforOW(testlik);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        for(int t=0;t<maxOWsize;t++){
          temp_olik[t]=testlik[owindices[t]];
          temp_otrees[t]=testtree[owindices[t]];
          temp_omat[t]=testmat[owindices[t]];
          temp_oparent[t]=testpar[owindices[t]];
        }
        testlik=temp_olik;
        testtree=temp_otrees;
        testmat=temp_omat;
        testpar=temp_oparent;
      }
      // ret[0]=lowest_BIC;
      // ret[1]=best_sv;
      // ret[2]=best_sp;
      // ret[3]=likeliest_tree;
      // ret[4]=testtree;
      // ret[5]=testlik;
      // ret[6]=testmat;
      // ret[7]=testpar;
      // ret[8]=no_tree_err;
      
      
      ret[0]=testtree;
      ret[1]=testlik;
      ret[2]=testmat;
      ret[3]=testpar;
      
      return (ret);
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;
      List gr(1);
      gr[0]=no_tree_err;
      return(gr);
    }
  }else{
    no_tree_err=1;
    List gr(1);
    gr[0]=no_tree_err;
    return(gr);
  }
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_sum(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                        arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,
                        double a,double mu,double nu,double lambda,double c,double lowest_BIC,
                        int parent,NumericMatrix cp_mat,double alpha,double beta,int maxOWsize,//int first_round,
                        List sum_trees,List sum_trees_mat,NumericVector y_scaled,IntegerVector parent2,int i,
                        unsigned int min_num_obs_for_split,unsigned int min_num_obs_after_split){
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;
  NumericMatrix treetable_c=treetable;
  NumericMatrix treemat_c=tree_mat;
  
  NumericVector terminal_nodes=find_term_nodes(treetable_c);
  //IntegerVector change_node1;
  int list_size=1000;
  std::vector<double> tree_lik(list_size);
  List proposal_tree;
  //List ret(9);
  List ret(4);
  bool no_tree_err=0;
  //List likeliest_tree;
  List tree_list(list_size);
  List tree_mat_list(list_size);
  int count=0;
  //std::vector<int> tree_parent(list_size);
  //int best_sv;
  //double best_sp;
  double tree_prior=1;
  //List changetree;
  double BIC;
  //int p;
  //int p_other=0;
  List eval_model;
  //NumericVector int_nodes;
  //NumericVector other_int_nodes;
  //arma::colvec curr_col=data.col(0);
  //arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[0]);
  //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[0]));
  //arma::mat data_curr_node=data.rows(grow_obs);
  //double d=d1[0];
  //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),cp_mat(0,0)+1);
  double lik;
  
  for(int l=0;l<terminal_nodes.size();l++){
    //loop over each terminal node
    arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[l]);
    //depth of tree at current terminal node
    //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[l]));
    arma::mat data_curr_node=data.rows(grow_obs);
    //double d=d1[0];
    int d = find_term_cols(treemat_c,terminal_nodes[l]);
    
    int w=cp_mat.nrow();
    if(data_curr_node.n_rows<=min_num_obs_for_split){
      throw std::range_error("not enough obs in node to grow any further");
      //continue;
    }
    for(int k=0;k<w;k++){
      //p_other=0;
      split_var=cp_mat(k,0)+1;
      //arma::colvec curr_cols=data.col(split_var-1);
      //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),split_var);
      
      // The following lines are unnecessary because get_min.size()=data_curr_node.n_rows above
      // if(get_min.size()<=2){
      //   throw std::range_error("obs in this terminal node are too small");
      // }
      
      double split_point=cp_mat(k,1);
      arma::vec curr_cols2=data_curr_node.col(split_var-1);
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));
      
      if(ld_prop.size()<=min_num_obs_after_split || rd_prop.size()<=min_num_obs_after_split){
        continue;
      }
      proposal_tree=grow_tree(data,//resids,
                              treemat_c,terminal_nodes[l],treetable_c,split_var,
                              split_point,//terminal_nodes,
                              wrap(arma::conv_to<arma::vec>::from(grow_obs)),
                              d//,get_min,data_curr_node
      );
      
      //Test lines below have been removed
      //NumericMatrix test =proposal_tree[0];
      //NumericMatrix test1 =proposal_tree[1];
      
      //if(test1.ncol()==3){
      //  NumericVector u1=unique(test1(_,0));
      //  NumericVector u2=unique(test1(_,1));
      //  NumericVector u3=unique(test1(_,2));
      //}
      
      //It should not be possible for get_best_split to be used outside of the first round, therefore removing the if condition below
      //if(first_round==1){
      //  lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);
      // }else{
      SEXP s = sum_trees[parent2[i]];
      if(is<List>(s)){
        List sum_trees2=sum_trees[parent2[i]];
        List sum_trees_mat2=sum_trees_mat[parent2[i]];
        sum_trees2.push_back(proposal_tree[0]);
        sum_trees_mat2.push_back(proposal_tree[1]);
        lik=sumtree_likelihood_function2(y_scaled,sum_trees2,sum_trees_mat2,y_scaled.size(),a,nu,lambda);
        for(int t=0;t<sum_trees2.size();t++){
          NumericMatrix tree=sum_trees2[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          NumericMatrix mat=sum_trees_mat2[t];
          tree_prior*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
      }else{
        NumericMatrix sum_trees2=sum_trees[parent2[i]];
        NumericMatrix sum_trees_mat2=sum_trees_mat[parent2[i]];
        //other_int_nodes = find_term_nodes(sum_trees2);
        //p_other=other_int_nodes.size();
        List st(2);
        List st_mat(2);
        st[0]=sum_trees2;
        st[1]=proposal_tree[0];
        st_mat[0]=sum_trees_mat2;
        st_mat[1]=proposal_tree[1];
        // return(st);
        lik=sumtree_likelihood_function2(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        for(int t=0;t<st.size();t++){
          NumericMatrix tree=st[t];
          NumericMatrix mat=st_mat[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
      }
      // }It should not be possible for get_best_split to be used outside o
      
      
      // FIXED (I think) at the moment tree prior is only for current tree need to get it for entire sum of tree list.
      
      //int_nodes=find_term_nodes(proposal_tree[0]);
      //p=int_nodes.size()+p_other; //COULD ADD 1 for the variance parameter, and more numbers for other parameters. (This might influence probability-weights, but not the ranking of BICs)
      //BIC=-2*(lik+log(tree_prior))+p_other*log(data.n_rows);  
      BIC=-2*(lik+log(tree_prior));  
      if(BIC<lowest_BIC){
        lowest_BIC=BIC;
        //best_sv=split_var;
        //best_sp=split_point;
        //likeliest_tree=proposal_tree;
        tree_list[count]=proposal_tree[0];
        tree_mat_list[count]=proposal_tree[1];
        tree_lik[count]=BIC;
        //tree_parent[count]=parent;
        count++;
        if(count==(tree_list.size()-1)){
          list_size=list_size*2;
          tree_list=resize_bigger(tree_list,list_size);
          tree_mat_list=resize_bigger(tree_mat_list,list_size);
          tree_lik.resize(list_size);
          //tree_parent.resize(list_size);
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){
          if(is<NumericMatrix>(proposal_tree[0])){
            //std::cout<<"its a matrix "<<"\n";
          }else{
            throw std::range_error("proposal tree not a matrix");
          }
          tree_list[count]=proposal_tree[0];
          tree_mat_list[count]=proposal_tree[1];
          tree_lik[count]=BIC;
          //tree_parent[count]=parent;
          count++;
          if(count==(tree_list.size()-1)){
            list_size=list_size*2;
            tree_list=resize_bigger(tree_list,list_size);
            tree_mat_list=resize_bigger(tree_mat_list,list_size);
            tree_lik.resize(list_size);
            //tree_parent.resize(list_size);						  
          }
        }
      }
    }  
  }
  tree_list=resize(tree_list,count);
  tree_mat_list=resize(tree_mat_list,count);
  tree_lik.resize(count);
  //tree_parent.resize(count);
  IntegerVector tree_parent(count, parent);
  
  if(count>0){
    if(less_greedy==1){
      ret[0]=tree_list;
      ret[1]=tree_lik;
      ret[2]=tree_mat_list;
      ret[3]=tree_parent;
      
      return (ret);	
    }
    //eval_model=evaluate_model_occams_window(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));
    eval_model=evaluate_model_occams_window(wrap(tree_lik),lowest_BIC,c,wrap(tree_list),wrap(tree_mat_list),tree_parent);
    NumericVector testlik =eval_model[0];
    List testtree =eval_model[1];    
    List testmat =eval_model[2]; 
    IntegerVector testpar =eval_model[3];
    
    if(testlik.size()>0){
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){
        IntegerVector owindices=orderforOW(testlik);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        for(int t=0;t<maxOWsize;t++){
          temp_olik[t]=testlik[owindices[t]];
          temp_otrees[t]=testtree[owindices[t]];
          temp_omat[t]=testmat[owindices[t]];
          temp_oparent[t]=testpar[owindices[t]];
        }
        testlik=temp_olik;
        testtree=temp_otrees;
        testmat=temp_omat;
        testpar=temp_oparent;
      }
      // ret[0]=lowest_BIC;
      // ret[1]=best_sv;
      // ret[2]=best_sp;
      // ret[3]=likeliest_tree;
      // ret[4]=testtree;
      // ret[5]=testlik;
      // ret[6]=testmat;
      // ret[7]=testpar;
      // ret[8]=no_tree_err;
      
      ret[0]=testtree;
      ret[1]=testlik;
      ret[2]=testmat;
      ret[3]=testpar;
      
      return (ret);
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;
      List gr(1);
      gr[0]=no_tree_err;
      return(gr);
    }
  }else{
    no_tree_err=1;
    List gr(1);
    gr[0]=no_tree_err;
    return(gr);
  }
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_sum_2(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                          arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,
                          double a,double mu,double nu,double lambda,double c,double lowest_BIC,
                          int parent,List cp_matlist,double alpha,double beta,int maxOWsize,//int first_round,
                          List sum_trees,List sum_trees_mat,NumericVector y_scaled,IntegerVector parent2,int i,
                          unsigned int min_num_obs_for_split,unsigned int min_num_obs_after_split){
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;
  NumericMatrix treetable_c=treetable;
  NumericMatrix treemat_c=tree_mat;
  
  NumericVector terminal_nodes=find_term_nodes(treetable_c);
  //IntegerVector change_node1;
  int list_size=1000;
  std::vector<double> tree_lik(list_size);
  List proposal_tree;
  //List ret(9);
  List ret(4);
  bool no_tree_err=0;
  //List likeliest_tree;
  List tree_list(list_size);
  List tree_mat_list(list_size);
  int count=0;
  //std::vector<int> tree_parent(list_size);
  //int best_sv;
  //double best_sp;
  double tree_prior=1;
  //List changetree;
  double BIC;
  //int p;
  //int p_other=0;
  List eval_model;
  //NumericVector int_nodes;
  //NumericVector other_int_nodes;
  //arma::colvec curr_col=data.col(0);
  //arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[0]);
  //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[0]));
  //arma::mat data_curr_node=data.rows(grow_obs);
  //double d=d1[0];
  //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),cp_mat(0,0)+1);
  double lik;
  
  for(int l=0;l<terminal_nodes.size();l++){
    //loop over each terminal node
    arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[l]);
    //depth of tree at current terminal node
    //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[l]));
    arma::mat data_curr_node=data.rows(grow_obs);
    //double d=d1[0];
    int d = find_term_cols(treemat_c,terminal_nodes[l]);
    NumericMatrix cp_mat=cp_matlist[l];
    int w=cp_mat.nrow();
    if(data_curr_node.n_rows<=min_num_obs_for_split){
      throw std::range_error("not enough obs in node to grow any further");
      //continue;
    }
    for(int k=0;k<w;k++){
      //p_other=0;
      split_var=cp_mat(k,0)+1;
      //arma::colvec curr_cols=data.col(split_var-1);
      //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),split_var);
      
      // The following lines are unnecessary because get_min.size()=data_curr_node.n_rows above
      // if(get_min.size()<=2){
      //   throw std::range_error("obs in this terminal node are too small");
      // }
      
      double split_point=cp_mat(k,1);
      arma::vec curr_cols2=data_curr_node.col(split_var-1);
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));
      
      if(ld_prop.size()<=min_num_obs_after_split || rd_prop.size()<=min_num_obs_after_split){
        continue;
      }
      proposal_tree=grow_tree(data,//resids,
                              treemat_c,terminal_nodes[l],treetable_c,split_var,
                              split_point,//terminal_nodes,
                              wrap(arma::conv_to<arma::vec>::from(grow_obs)),
                              d//,get_min,data_curr_node
      );
      
      //Test lines below have been removed
      //NumericMatrix test =proposal_tree[0];
      //NumericMatrix test1 =proposal_tree[1];
      
      //if(test1.ncol()==3){
      //  NumericVector u1=unique(test1(_,0));
      //  NumericVector u2=unique(test1(_,1));
      //  NumericVector u3=unique(test1(_,2));
      //}
      
      //It should not be possible for get_best_split to be used outside of the first round, therefore removing the if condition below
      //if(first_round==1){
      //  lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);
      // }else{
      SEXP s = sum_trees[parent2[i]];
      if(is<List>(s)){
        List sum_trees2=sum_trees[parent2[i]];
        List sum_trees_mat2=sum_trees_mat[parent2[i]];
        sum_trees2.push_back(proposal_tree[0]);
        sum_trees_mat2.push_back(proposal_tree[1]);
        lik=sumtree_likelihood_function2(y_scaled,sum_trees2,sum_trees_mat2,y_scaled.size(),a,nu,lambda);
        for(int t=0;t<sum_trees2.size();t++){
          NumericMatrix tree=sum_trees2[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          NumericMatrix mat=sum_trees_mat2[t];
          tree_prior*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
      }else{
        NumericMatrix sum_trees2=sum_trees[parent2[i]];
        NumericMatrix sum_trees_mat2=sum_trees_mat[parent2[i]];
        //other_int_nodes = find_term_nodes(sum_trees2);
        //p_other=other_int_nodes.size();
        List st(2);
        List st_mat(2);
        st[0]=sum_trees2;
        st[1]=proposal_tree[0];
        st_mat[0]=sum_trees_mat2;
        st_mat[1]=proposal_tree[1];
        // return(st);
        lik=sumtree_likelihood_function2(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        for(int t=0;t<st.size();t++){
          NumericMatrix tree=st[t];
          NumericMatrix mat=st_mat[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
      }
      // }It should not be possible for get_best_split to be used outside o
      
      
      // FIXED (I think) at the moment tree prior is only for current tree need to get it for entire sum of tree list.
      
      //int_nodes=find_term_nodes(proposal_tree[0]);
      //p=int_nodes.size()+p_other; //COULD ADD 1 for the variance parameter, and more numbers for other parameters. (This might influence probability-weights, but not the ranking of BICs)
      //BIC=-2*(lik+log(tree_prior))+p_other*log(data.n_rows);  
      BIC=-2*(lik+log(tree_prior));  
      if(BIC<lowest_BIC){
        lowest_BIC=BIC;
        //best_sv=split_var;
        //best_sp=split_point;
        //likeliest_tree=proposal_tree;
        tree_list[count]=proposal_tree[0];
        tree_mat_list[count]=proposal_tree[1];
        tree_lik[count]=BIC;
        //tree_parent[count]=parent;
        count++;
        if(count==(tree_list.size()-1)){
          list_size=list_size*2;
          tree_list=resize_bigger(tree_list,list_size);
          tree_mat_list=resize_bigger(tree_mat_list,list_size);
          tree_lik.resize(list_size);
          //tree_parent.resize(list_size);
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){
          if(is<NumericMatrix>(proposal_tree[0])){
            //std::cout<<"its a matrix "<<"\n";
          }else{
            throw std::range_error("proposal tree not a matrix");
          }
          tree_list[count]=proposal_tree[0];
          tree_mat_list[count]=proposal_tree[1];
          tree_lik[count]=BIC;
          //tree_parent[count]=parent;
          count++;
          if(count==(tree_list.size()-1)){
            list_size=list_size*2;
            tree_list=resize_bigger(tree_list,list_size);
            tree_mat_list=resize_bigger(tree_mat_list,list_size);
            tree_lik.resize(list_size);
            //tree_parent.resize(list_size);						  
          }
        }
      }
    }  
  }
  tree_list=resize(tree_list,count);
  tree_mat_list=resize(tree_mat_list,count);
  tree_lik.resize(count);
  //tree_parent.resize(count);
  IntegerVector tree_parent(count, parent);
  
  if(count>0){
    if(less_greedy==1){
      ret[0]=tree_list;
      ret[1]=tree_lik;
      ret[2]=tree_mat_list;
      ret[3]=tree_parent;
      
      return (ret);	
    }
    //eval_model=evaluate_model_occams_window(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));
    eval_model=evaluate_model_occams_window(wrap(tree_lik),lowest_BIC,c,wrap(tree_list),wrap(tree_mat_list),tree_parent);
    NumericVector testlik =eval_model[0];
    List testtree =eval_model[1];    
    List testmat =eval_model[2]; 
    IntegerVector testpar =eval_model[3];
    
    if(testlik.size()>0){
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){
        IntegerVector owindices=orderforOW(testlik);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        for(int t=0;t<maxOWsize;t++){
          temp_olik[t]=testlik[owindices[t]];
          temp_otrees[t]=testtree[owindices[t]];
          temp_omat[t]=testmat[owindices[t]];
          temp_oparent[t]=testpar[owindices[t]];
        }
        testlik=temp_olik;
        testtree=temp_otrees;
        testmat=temp_omat;
        testpar=temp_oparent;
      }
      // ret[0]=lowest_BIC;
      // ret[1]=best_sv;
      // ret[2]=best_sp;
      // ret[3]=likeliest_tree;
      // ret[4]=testtree;
      // ret[5]=testlik;
      // ret[6]=testmat;
      // ret[7]=testpar;
      // ret[8]=no_tree_err;
      
      ret[0]=testtree;
      ret[1]=testlik;
      ret[2]=testmat;
      ret[3]=testpar;
      
      return (ret);
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;
      List gr(1);
      gr[0]=no_tree_err;
      return(gr);
    }
  }else{
    no_tree_err=1;
    List gr(1);
    gr[0]=no_tree_err;
    return(gr);
  }
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_exact(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                          NumericVector resids,arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,
                          double a,double mu,double nu,double lambda,double c,double lowest_BIC,int parent
                            ,NumericMatrix cp_mat,double alpha,double beta,int maxOWsize, unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split//,int first_round
){
  //this function will search through all predictive split points and return those within Occam's Window.
  
  // Rcout << "Line 2733 .\n";
  
  int split_var;
  NumericMatrix treetable_c=treetable;
  NumericMatrix treemat_c=tree_mat;
  
  NumericVector terminal_nodes=find_term_nodes(treetable_c);
  //IntegerVector change_node1;
  int list_size=1000;
  std::vector<double> tree_lik(list_size);
  List proposal_tree;
  //List ret(9);
  //List ret(4);
  List ret(5);
  
  bool no_tree_err=0;
  //List likeliest_tree;
  List tree_list(list_size);
  List tree_mat_list(list_size);
  
  List tree_preds(list_size);
  
  
  int count=0;
  //std::vector<int> tree_parent(list_size);
  //int best_sv;
  //double best_sp;
  double tree_prior=0;
  //List changetree;
  double BIC;
  //int p;
  List eval_model;
  //NumericVector int_nodes;
  //arma::colvec curr_col=data.col(0);
  //arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[0]);
  //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[0]));
  //arma::mat data_curr_node=data.rows(grow_obs);
  //double d=d1[0];
  //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),cp_mat(0,0)+1);
  double lik;
  List lik_listoutput;
  // Rcout << "Line 2774 .\n";
  
  for(int l=0;l<terminal_nodes.size();l++){
    //loop over each terminal node
    arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[l]);
    //depth of tree at current terminal node
    //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[l]));
    arma::mat data_curr_node=data.rows(grow_obs);
    //double d=d1[0];
    int d = find_term_cols(treemat_c,terminal_nodes[l]);
    int w=cp_mat.nrow();
    if(data_curr_node.n_rows<=min_num_obs_for_split){
      throw std::range_error("not enough obs in node to grow any further");
      //continue;
    }
    for(int k=0;k<w;k++){
      
      // Rcout << "Line 2791 in get_best_split_exact(). before get_tree_prior() .\n";
      
      split_var=cp_mat(k,0)+1;
      //arma::colvec curr_cols=data.col(split_var-1);
      //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),split_var);
      // The following lines are unnecessary because get_min.size()=data_curr_node.n_rows above
      // if(get_min.size()<=2){
      //   throw std::range_error("obs in this terminal node are too small");
      // }
      
      double split_point=cp_mat(k,1);
      // Rcout << "Line 2642 in get_best_split_exact().\n";
      
      arma::vec curr_cols2=data_curr_node.col(split_var-1);
      // Rcout << "Line 2645 in get_best_split_exact() .\n";
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));
      
      if(ld_prop.size()<=min_num_obs_after_split || rd_prop.size()<=min_num_obs_after_split){
        continue;
      }
      
      // Rcout << "Line 2814 in get_best_split_exact() .\n";
      
      proposal_tree=grow_tree(data,//resids,
                              treemat_c,terminal_nodes[l],treetable_c,
                              split_var,split_point,//terminal_nodes,
                              wrap(arma::conv_to<arma::vec>::from(grow_obs)),
                              d//,get_min,data_curr_node
      );
      
      // Test lines below have been removed
      // NumericMatrix test =proposal_tree[0];
      // NumericMatrix test1 =proposal_tree[1];
      
      // if(test1.ncol()==3){
      //   NumericVector u1=unique(test1(_,0));
      //   NumericVector u2=unique(test1(_,1));
      //   NumericVector u3=unique(test1(_,2));
      // }
      
      // get_best_split_exact should only be used in the first_round. Removing the if condition below 
      //if(first_round==1){
      //lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);
      
      // Rcout << "Line 2837 in get_best_split_exact() .\n";
      
      //NumericMatrix proptabtemp = proposal_tree[0];
      //NumericMatrix propmattemp = proposal_tree[1];
      
      // Rcout << "proposal_tree[0] = "<< proptabtemp << " .\n";
      // Rcout << "proposal_tree[1] = "<< propmattemp << " .\n";
      
      lik_listoutput=likelihood_function2_exact(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);
      // Rcout << "Line 2846 in get_best_split_exact().  .\n";
      
      lik=as<double>(lik_listoutput[0]);
      NumericVector temp_predvec=lik_listoutput[1];
      
      
      // }else{
      //   //have a sum of trees
      //   lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);  
      // }
      
      // Rcout << "Line 2857 in get_best_split_exact(). before get_tree_prior() .\n";
      
      tree_prior=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,proposal_tree[0],proposal_tree[1],alpha,beta);
      
      // Rcout << "Line 2696 in get_best_split_exact(). before get_tree_prior() .\n";
      // Rcout << "k = "<< k << " .\n";
      // Rcout << "w = "<< w << " .\n";
      // Rcout << "l = "<< l << " .\n";
      // Rcout << "terminal_nodes.size() = "<< terminal_nodes.size() << " .\n";
      
      
      //int_nodes=find_term_nodes(proposal_tree[0]);
      //p=int_nodes.size();
      //BIC=-2*(lik+log(tree_prior))+p*log(data.n_rows); 
      BIC=-2*(lik+log(tree_prior));  
      
      if(BIC<lowest_BIC){
        lowest_BIC=BIC;
        //best_sv=split_var;
        //best_sp=split_point;
        //likeliest_tree=proposal_tree;
        tree_list[count]=proposal_tree[0];
        tree_mat_list[count]=proposal_tree[1];
        tree_lik[count]=BIC;
        //tree_parent[count]=parent;
        
        tree_preds[count]=temp_predvec;
        
        
        
        count++;
        if(count==(tree_list.size()-1)){
          list_size=list_size*2;
          tree_list=resize_bigger(tree_list,list_size);
          tree_mat_list=resize_bigger(tree_mat_list,list_size);
          tree_lik.resize(list_size);
          //tree_parent.resize(list_size);
          
          tree_preds=resize_bigger(tree_preds,list_size);
          
          
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){
          // if(is<NumericMatrix>(proposal_tree[0])){
          // }else{
          //   throw std::range_error("proposal tree not a matrix");
          // }
          tree_list[count]=proposal_tree[0];
          tree_mat_list[count]=proposal_tree[1];
          tree_lik[count]=BIC;
          //tree_parent[count]=parent;
          
          tree_preds[count]=temp_predvec;
          
          count++;
          if(count==(tree_list.size()-1)){
            list_size=list_size*2;
            tree_list=resize_bigger(tree_list,list_size);
            tree_mat_list=resize_bigger(tree_mat_list,list_size);
            tree_lik.resize(list_size);
            //tree_parent.resize(list_size);
            
            tree_preds=resize_bigger(tree_preds,list_size);
            
          }
        }
      }
    }  
  }
  
  // Rcout << "Line 2928 in get_best_split_exact().\n";
  
  tree_list=resize(tree_list,count);
  tree_mat_list=resize(tree_mat_list,count);
  tree_lik.resize(count);
  //tree_parent.resize(count);
  IntegerVector tree_parent(count, parent);
  
  tree_preds=resize(tree_preds,count);
  
  // Rcout << "Line 2930 in get_best_split_exact().\n";
  
  if(count>0){
    if(less_greedy==1){
      ret[0]=tree_list;
      ret[1]=tree_lik;
      ret[2]=tree_mat_list;
      ret[3]=tree_parent;
      ret[4]=tree_preds;
      
      return (ret);	
    }
    //eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));
    //eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),tree_parent);
    eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,c,wrap(tree_list),wrap(tree_mat_list),tree_parent,
                                                  tree_preds);
    
    NumericVector testlik =eval_model[0];
    List testtree =eval_model[1];    
    List testmat =eval_model[2]; 
    IntegerVector testpar =eval_model[3];
    
    List testpredlist=eval_model[4];
    
    if(testlik.size()>0){
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){
        IntegerVector owindices=orderforOW(testlik);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        
        List temp_opreds(maxOWsize);
        
        
        for(int t=0;t<maxOWsize;t++){
          temp_olik[t]=testlik[owindices[t]];
          temp_otrees[t]=testtree[owindices[t]];
          temp_omat[t]=testmat[owindices[t]];
          temp_oparent[t]=testpar[owindices[t]];
          temp_opreds[t]=testpredlist[owindices[t]];
        }
        testlik=temp_olik;
        testtree=temp_otrees;
        testmat=temp_omat;
        testpar=temp_oparent;
        testpredlist=temp_opreds;
      }
      // ret[0]=lowest_BIC;
      // ret[1]=best_sv;
      // ret[2]=best_sp;
      // ret[3]=likeliest_tree;
      // ret[4]=testtree;
      // ret[5]=testlik;
      // ret[6]=testmat;
      // ret[7]=testpar;
      // ret[8]=no_tree_err;
      
      
      ret[0]=testtree;
      ret[1]=testlik;
      ret[2]=testmat;
      ret[3]=testpar;
      ret[4]=testpredlist;
      
      return (ret);
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;
      List gr(1);
      gr[0]=no_tree_err;
      return(gr);
    }
  }else{
    no_tree_err=1;
    List gr(1);
    gr[0]=no_tree_err;
    return(gr);
  }
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_2_exact(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                            NumericVector resids,arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,
                            double a,double mu,double nu,double lambda,double c,double lowest_BIC,int parent
                              ,List cp_matlist,double alpha,double beta,int maxOWsize, unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split//,int first_round
){
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;
  NumericMatrix treetable_c=treetable;
  NumericMatrix treemat_c=tree_mat;
  
  NumericVector terminal_nodes=find_term_nodes(treetable_c);
  //IntegerVector change_node1;
  int list_size=1000;
  std::vector<double> tree_lik(list_size);
  List proposal_tree;
  //List ret(9);
  //List ret(4);
  List ret(5);
  
  bool no_tree_err=0;
  //List likeliest_tree;
  List tree_list(list_size);
  List tree_mat_list(list_size);
  
  List tree_preds(list_size);
  
  
  int count=0;
  //std::vector<int> tree_parent(list_size);
  //int best_sv;
  //double best_sp;
  double tree_prior=0;
  //List changetree;
  double BIC;
  //int p;
  List eval_model;
  //NumericVector int_nodes;
  //arma::colvec curr_col=data.col(0);
  //arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[0]);
  //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[0]));
  //arma::mat data_curr_node=data.rows(grow_obs);
  //double d=d1[0];
  //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),cp_mat(0,0)+1);
  double lik;
  List lik_listoutput;
  
  for(int l=0;l<terminal_nodes.size();l++){
    NumericMatrix cp_mat=cp_matlist[l];
    //loop over each terminal node
    arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[l]);
    //depth of tree at current terminal node
    //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[l]));
    arma::mat data_curr_node=data.rows(grow_obs);
    //double d=d1[0];
    int d = find_term_cols(treemat_c,terminal_nodes[l]);
    int w=cp_mat.nrow();
    if(data_curr_node.n_rows<=min_num_obs_for_split){
      throw std::range_error("not enough obs in node to grow any further");
      //continue;
    }
    for(int k=0;k<w;k++){
      split_var=cp_mat(k,0)+1;
      //arma::colvec curr_cols=data.col(split_var-1);
      //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),split_var);
      // The following lines are unnecessary because get_min.size()=data_curr_node.n_rows above
      // if(get_min.size()<=2){
      //   throw std::range_error("obs in this terminal node are too small");
      // }
      
      double split_point=cp_mat(k,1);
      arma::vec curr_cols2=data_curr_node.col(split_var-1);
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));
      
      if(ld_prop.size()<=min_num_obs_after_split || rd_prop.size()<=min_num_obs_after_split){
        continue;
      }
      proposal_tree=grow_tree(data,//resids,
                              treemat_c,terminal_nodes[l],treetable_c,
                              split_var,split_point,//terminal_nodes,
                              wrap(arma::conv_to<arma::vec>::from(grow_obs)),
                              d//,get_min,data_curr_node
      );
      
      // Test lines below have been removed
      // NumericMatrix test =proposal_tree[0];
      // NumericMatrix test1 =proposal_tree[1];
      
      // if(test1.ncol()==3){
      //   NumericVector u1=unique(test1(_,0));
      //   NumericVector u2=unique(test1(_,1));
      //   NumericVector u3=unique(test1(_,2));
      // }
      
      // get_best_split_exact should only be used in the first_round. Removing the if condition below 
      //if(first_round==1){
      lik_listoutput=likelihood_function2_exact(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);
      lik=as<double>(lik_listoutput[0]);
      NumericVector temp_predvec=lik_listoutput[1];
      
      
      // }else{
      //   //have a sum of trees
      //   lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);  
      // }
      tree_prior=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,proposal_tree[0],proposal_tree[1],alpha,beta);
      //int_nodes=find_term_nodes(proposal_tree[0]);
      //p=int_nodes.size();
      //BIC=-2*(lik+log(tree_prior))+p*log(data.n_rows); 
      BIC=-2*(lik+log(tree_prior));  
      
      if(BIC<lowest_BIC){
        lowest_BIC=BIC;
        //best_sv=split_var;
        //best_sp=split_point;
        //likeliest_tree=proposal_tree;
        tree_list[count]=proposal_tree[0];
        tree_mat_list[count]=proposal_tree[1];
        tree_lik[count]=BIC;
        //tree_parent[count]=parent;
        
        tree_preds[count]=temp_predvec;
        
        
        
        count++;
        if(count==(tree_list.size()-1)){
          list_size=list_size*2;
          tree_list=resize_bigger(tree_list,list_size);
          tree_mat_list=resize_bigger(tree_mat_list,list_size);
          tree_lik.resize(list_size);
          //tree_parent.resize(list_size);
          
          tree_preds=resize_bigger(tree_preds,list_size);
          
          
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){
          // if(is<NumericMatrix>(proposal_tree[0])){
          // }else{
          //   throw std::range_error("proposal tree not a matrix");
          // }
          tree_list[count]=proposal_tree[0];
          tree_mat_list[count]=proposal_tree[1];
          tree_lik[count]=BIC;
          //tree_parent[count]=parent;
          
          tree_preds[count]=temp_predvec;
          
          count++;
          if(count==(tree_list.size()-1)){
            list_size=list_size*2;
            tree_list=resize_bigger(tree_list,list_size);
            tree_mat_list=resize_bigger(tree_mat_list,list_size);
            tree_lik.resize(list_size);
            //tree_parent.resize(list_size);
            
            tree_preds=resize_bigger(tree_preds,list_size);
            
          }
        }
      }
    }  
  }
  tree_list=resize(tree_list,count);
  tree_mat_list=resize(tree_mat_list,count);
  tree_lik.resize(count);
  //tree_parent.resize(count);
  IntegerVector tree_parent(count, parent);
  
  tree_preds=resize(tree_preds,count);
  
  
  if(count>0){
    if(less_greedy==1){
      ret[0]=tree_list;
      ret[1]=tree_lik;
      ret[2]=tree_mat_list;
      ret[3]=tree_parent;
      ret[4]=tree_preds;
      
      return (ret);	
    }
    //eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));
    //eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),tree_parent);
    eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,c,wrap(tree_list),wrap(tree_mat_list),tree_parent,
                                                  tree_preds);
    
    NumericVector testlik =eval_model[0];
    List testtree =eval_model[1];    
    List testmat =eval_model[2]; 
    IntegerVector testpar =eval_model[3];
    
    List testpredlist=eval_model[4];
    
    if(testlik.size()>0){
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){
        IntegerVector owindices=orderforOW(testlik);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        
        List temp_opreds(maxOWsize);
        
        
        for(int t=0;t<maxOWsize;t++){
          temp_olik[t]=testlik[owindices[t]];
          temp_otrees[t]=testtree[owindices[t]];
          temp_omat[t]=testmat[owindices[t]];
          temp_oparent[t]=testpar[owindices[t]];
          temp_opreds[t]=testpredlist[owindices[t]];
        }
        testlik=temp_olik;
        testtree=temp_otrees;
        testmat=temp_omat;
        testpar=temp_oparent;
        testpredlist=temp_opreds;
      }
      // ret[0]=lowest_BIC;
      // ret[1]=best_sv;
      // ret[2]=best_sp;
      // ret[3]=likeliest_tree;
      // ret[4]=testtree;
      // ret[5]=testlik;
      // ret[6]=testmat;
      // ret[7]=testpar;
      // ret[8]=no_tree_err;
      
      
      ret[0]=testtree;
      ret[1]=testlik;
      ret[2]=testmat;
      ret[3]=testpar;
      ret[4]=testpredlist;
      
      return (ret);
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;
      List gr(1);
      gr[0]=no_tree_err;
      return(gr);
    }
  }else{
    no_tree_err=1;
    List gr(1);
    gr[0]=no_tree_err;
    return(gr);
  }
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_sum_exact(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                              arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,
                              double a,double mu,double nu,double lambda,double c,double lowest_BIC,
                              int parent,NumericMatrix cp_mat,double alpha,double beta,int maxOWsize,//int first_round,
                              List sum_trees,List sum_trees_mat,NumericVector y_scaled,IntegerVector parent2,int i,
                              unsigned int min_num_obs_for_split,unsigned int min_num_obs_after_split){
  //this function will search through all predictive split points and return those within Occam's Window.
  
  // Rcout << " line 3296. no update .\n";
  
  int split_var;
  NumericMatrix treetable_c=treetable;
  NumericMatrix treemat_c=tree_mat;
  //Rcout << " line 3223. no update .\n";
  
  NumericVector terminal_nodes=find_term_nodes(treetable_c);
  //IntegerVector change_node1;
  int list_size=1000;
  std::vector<double> tree_lik(list_size);
  List proposal_tree;
  //List ret(9);
  //List ret(4);
  List ret(5);
  //Rcout << " line 3233. no update .\n";
  
  bool no_tree_err=0;
  //List likeliest_tree;
  List tree_list(list_size);
  List tree_mat_list(list_size);
  
  List tree_preds(list_size);
  //Rcout << " line 3241. no update .\n";
  
  int count=0;
  //std::vector<int> tree_parent(list_size);
  //int best_sv;
  //double best_sp;
  double tree_prior=1;
  //List changetree;
  double BIC;
  //int p;
  //int p_other=0;
  List eval_model;
  //NumericVector int_nodes;
  //NumericVector other_int_nodes;
  //arma::colvec curr_col=data.col(0);
  //arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[0]);
  //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[0]));
  //arma::mat data_curr_node=data.rows(grow_obs);
  //double d=d1[0];
  //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),cp_mat(0,0)+1);
  double lik;
  List lik_list;
  NumericVector temppredoutput;
  
  // Rcout << " line 3343. no update .\n";
  
  for(int l=0;l<terminal_nodes.size();l++){
    //loop over each terminal node
    arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[l]);
    //depth of tree at current terminal node
    //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[l]));
    arma::mat data_curr_node=data.rows(grow_obs);
    //double d=d1[0];
    int d = find_term_cols(treemat_c,terminal_nodes[l]);
    
    // Rcout << " line 3354. no update .\n";
    
    int w=cp_mat.nrow();
    if(data_curr_node.n_rows<=min_num_obs_for_split){
      throw std::range_error("not enough obs in node to grow any further");
      //continue;
    }
    for(int k=0;k<w;k++){
      // Rcout << " line 3362. no update .\n";
      
      //p_other=0;
      split_var=cp_mat(k,0)+1;
      //arma::colvec curr_cols=data.col(split_var-1);
      //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),split_var);
      
      // The following lines are unnecessary because get_min.size()=data_curr_node.n_rows above
      // if(get_min.size()<=2){
      //   throw std::range_error("obs in this terminal node are too small");
      // }
      
      double split_point=cp_mat(k,1);
      arma::vec curr_cols2=data_curr_node.col(split_var-1);
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));
      
      if(ld_prop.size()<=min_num_obs_after_split || rd_prop.size()<=min_num_obs_after_split){
        continue;
      }
      
      
      proposal_tree=grow_tree(data,//resids,
                              treemat_c,terminal_nodes[l],treetable_c,split_var,
                              split_point,//terminal_nodes,
                              wrap(arma::conv_to<arma::vec>::from(grow_obs)),
                              d//,get_min,data_curr_node
      );
      
      //Test lines below have been removed
      //NumericMatrix test =proposal_tree[0];
      //NumericMatrix test1 =proposal_tree[1];
      
      //if(test1.ncol()==3){
      //  NumericVector u1=unique(test1(_,0));
      //  NumericVector u2=unique(test1(_,1));
      //  NumericVector u3=unique(test1(_,2));
      //}
      
      //It should not be possible for get_best_split_exact to be used outside of the first round, therefore removing the if condition below
      //if(first_round==1){
      //  lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);
      // }else{
      
      //Rcout << " line 3406. no update .\n";
      
      SEXP s = sum_trees[parent2[i]];
      if(is<List>(s)){
        List sum_trees2=sum_trees[parent2[i]];
        List sum_trees_mat2=sum_trees_mat[parent2[i]];
        sum_trees2.push_back(proposal_tree[0]);
        sum_trees_mat2.push_back(proposal_tree[1]);
        //lik=sumtree_likelihood_function2_exact(y_scaled,sum_trees2,sum_trees_mat2,y_scaled.size(),a,nu,lambda);
        lik_list=sumtree_likelihood_function2_exact(y_scaled,sum_trees2,sum_trees_mat2,y_scaled.size(),a,nu,lambda);
        lik=as<double>(lik_list[0]);
        temppredoutput=lik_list[1];
        
        for(int t=0;t<sum_trees2.size();t++){
          NumericMatrix tree=sum_trees2[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          NumericMatrix mat=sum_trees_mat2[t];
          tree_prior*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
      }else{
        NumericMatrix sum_trees2=sum_trees[parent2[i]];
        NumericMatrix sum_trees_mat2=sum_trees_mat[parent2[i]];
        //other_int_nodes = find_term_nodes(sum_trees2);
        //p_other=other_int_nodes.size();
        List st(2);
        List st_mat(2);
        st[0]=sum_trees2;
        st[1]=proposal_tree[0];
        st_mat[0]=sum_trees_mat2;
        st_mat[1]=proposal_tree[1];
        // return(st);
        //lik=sumtree_likelihood_function2_exact(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        lik_list=sumtree_likelihood_function2_exact(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        lik=as<double>(lik_list[0]);
        temppredoutput=lik_list[1];
        
        for(int t=0;t<st.size();t++){
          NumericMatrix tree=st[t];
          NumericMatrix mat=st_mat[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
      }
      // }It should not be possible for get_best_split_exact to be used outside o
      
      //Rcout << " line 3454. no update .\n";
      
      // FIXED (I think) at the moment tree prior is only for current tree need to get it for entire sum of tree list.
      
      //int_nodes=find_term_nodes(proposal_tree[0]);
      //p=int_nodes.size()+p_other; //COULD ADD 1 for the variance parameter, and more numbers for other parameters. (This might influence probability-weights, but not the ranking of BICs)
      //BIC=-2*(lik+log(tree_prior))+p_other*log(data.n_rows);  
      BIC=-2*(lik+log(tree_prior));  
      if(BIC<lowest_BIC){
        lowest_BIC=BIC;
        //best_sv=split_var;
        //best_sp=split_point;
        //likeliest_tree=proposal_tree;
        tree_list[count]=proposal_tree[0];
        tree_mat_list[count]=proposal_tree[1];
        tree_lik[count]=BIC;
        //tree_parent[count]=parent;
        
        tree_preds[count]=temppredoutput;
        
        count++;
        if(count==(tree_list.size()-1)){
          list_size=list_size*2;
          tree_list=resize_bigger(tree_list,list_size);
          tree_mat_list=resize_bigger(tree_mat_list,list_size);
          tree_lik.resize(list_size);
          //tree_parent.resize(list_size);
          
          tree_preds=resize_bigger(tree_preds,list_size);
          
          
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){
          if(is<NumericMatrix>(proposal_tree[0])){
            //std::cout<<"its a matrix "<<"\n";
          }else{
            throw std::range_error("proposal tree not a matrix");
          }
          tree_list[count]=proposal_tree[0];
          tree_mat_list[count]=proposal_tree[1];
          tree_lik[count]=BIC;
          //tree_parent[count]=parent;
          
          tree_preds[count]=temppredoutput;
          
          count++;
          if(count==(tree_list.size()-1)){
            list_size=list_size*2;
            tree_list=resize_bigger(tree_list,list_size);
            tree_mat_list=resize_bigger(tree_mat_list,list_size);
            tree_lik.resize(list_size);
            //tree_parent.resize(list_size);
            
            tree_preds=resize_bigger(tree_preds,list_size);
            
          }
        }
      }
    }  
  }
  tree_list=resize(tree_list,count);
  tree_mat_list=resize(tree_mat_list,count);
  tree_lik.resize(count);
  //tree_parent.resize(count);
  
  tree_preds=resize(tree_preds,count);
  
  // Rcout << " line 3518. no update .\n";
  
  IntegerVector tree_parent(count, parent);
  
  if(count>0){
    if(less_greedy==1){
      ret[0]=tree_list;
      ret[1]=tree_lik;
      ret[2]=tree_mat_list;
      ret[3]=tree_parent;
      ret[4]=tree_preds;
      
      return (ret);	
    }
    //eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));
    //eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),tree_parent);
    eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,c,wrap(tree_list),wrap(tree_mat_list),tree_parent,
                                                  tree_preds);
    NumericVector testlik =eval_model[0];
    List testtree =eval_model[1];    
    List testmat =eval_model[2]; 
    IntegerVector testpar =eval_model[3];
    List testpredlist=eval_model[4];
    
    // Rcout << " line 3541. no update .\n";
    
    if(testlik.size()>0){
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){
        IntegerVector owindices=orderforOW(testlik);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        List temp_opreds(maxOWsize);
        
        for(int t=0;t<maxOWsize;t++){
          temp_olik[t]=testlik[owindices[t]];
          temp_otrees[t]=testtree[owindices[t]];
          temp_omat[t]=testmat[owindices[t]];
          temp_oparent[t]=testpar[owindices[t]];
          temp_opreds[t]=testpredlist[owindices[t]];
          
        }
        testlik=temp_olik;
        testtree=temp_otrees;
        testmat=temp_omat;
        testpar=temp_oparent;
        testpredlist=temp_opreds;
        
      }
      // ret[0]=lowest_BIC;
      // ret[1]=best_sv;
      // ret[2]=best_sp;
      // ret[3]=likeliest_tree;
      // ret[4]=testtree;
      // ret[5]=testlik;
      // ret[6]=testmat;
      // ret[7]=testpar;
      // ret[8]=no_tree_err;
      
      ret[0]=testtree;
      ret[1]=testlik;
      ret[2]=testmat;
      ret[3]=testpar;
      ret[4]=testpredlist;
      return (ret);
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;
      List gr(1);
      gr[0]=no_tree_err;
      return(gr);
    }
  }else{
    no_tree_err=1;
    List gr(1);
    gr[0]=no_tree_err;
    return(gr);
  }
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_sum_2_exact(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                                arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,
                                double a,double mu,double nu,double lambda,double c,double lowest_BIC,
                                int parent,List cp_matlist,double alpha,double beta,int maxOWsize,//int first_round,
                                List sum_trees,List sum_trees_mat,NumericVector y_scaled,IntegerVector parent2,int i,
                                unsigned int min_num_obs_for_split,unsigned int min_num_obs_after_split){
  //this function will search through all predictive split points and return those within Occam's Window.
  
  // Rcout << "Line 3603. \n";
  
  
  int split_var;
  NumericMatrix treetable_c=treetable;
  NumericMatrix treemat_c=tree_mat;
  
  NumericVector terminal_nodes=find_term_nodes(treetable_c);
  //IntegerVector change_node1;
  int list_size=1000;
  std::vector<double> tree_lik(list_size);
  List proposal_tree;
  //List ret(9);
  //List ret(4);
  List ret(5);
  
  
  bool no_tree_err=0;
  //List likeliest_tree;
  List tree_list(list_size);
  List tree_mat_list(list_size);
  
  List tree_preds(list_size);
  
  
  int count=0;
  //std::vector<int> tree_parent(list_size);
  //int best_sv;
  //double best_sp;
  double tree_prior=1;
  //List changetree;
  double BIC;
  //int p;
  //int p_other=0;
  List eval_model;
  //NumericVector int_nodes;
  //NumericVector other_int_nodes;
  //arma::colvec curr_col=data.col(0);
  //arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[0]);
  //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[0]));
  //arma::mat data_curr_node=data.rows(grow_obs);
  //double d=d1[0];
  //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),cp_mat(0,0)+1);
  double lik;
  List lik_list;
  NumericVector temppredoutput;
  
  // Rcout << "Line 3650. \n";
  
  for(int l=0;l<terminal_nodes.size();l++){
    
    // Rcout << "Line 3654. \n";
    
    
    //loop over each terminal node
    arma::uvec grow_obs=find_term_obs(treemat_c,terminal_nodes[l]);
    //depth of tree at current terminal node
    //NumericVector d1=unique(find_term_cols(treemat_c,terminal_nodes[l]));
    arma::mat data_curr_node=data.rows(grow_obs);
    //double d=d1[0];
    int d = find_term_cols(treemat_c,terminal_nodes[l]);
    NumericMatrix cp_mat=cp_matlist[l];
    int w=cp_mat.nrow();
    if(data_curr_node.n_rows<=min_num_obs_for_split){
      throw std::range_error("not enough obs in node to grow any further");
      //continue;
    }
    for(int k=0;k<w;k++){
      
      // Rcout << "Line 3672. \n";
      
      //p_other=0;
      split_var=cp_mat(k,0)+1;
      //arma::colvec curr_cols=data.col(split_var-1);
      //NumericVector get_min=get_grow_obs(data,wrap(grow_obs),split_var);
      
      // The following lines are unnecessary because get_min.size()=data_curr_node.n_rows above
      // if(get_min.size()<=2){
      //   throw std::range_error("obs in this terminal node are too small");
      // }
      
      double split_point=cp_mat(k,1);
      arma::vec curr_cols2=data_curr_node.col(split_var-1);
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));
      
      if(ld_prop.size()<=min_num_obs_after_split || rd_prop.size()<=min_num_obs_after_split){
        continue;
      }
      // Rcout << "Line 3693. \n";
      
      proposal_tree=grow_tree(data,//resids,
                              treemat_c,terminal_nodes[l],treetable_c,split_var,
                              split_point,//terminal_nodes,
                              wrap(arma::conv_to<arma::vec>::from(grow_obs)),
                              d//,get_min,data_curr_node
      );
      // Rcout << "Line 3701. \n";
      
      //Test lines below have been removed
      //NumericMatrix test =proposal_tree[0];
      //NumericMatrix test1 =proposal_tree[1];
      
      //if(test1.ncol()==3){
      //  NumericVector u1=unique(test1(_,0));
      //  NumericVector u2=unique(test1(_,1));
      //  NumericVector u3=unique(test1(_,2));
      //}
      
      //It should not be possible for get_best_split_exact to be used outside of the first round, therefore removing the if condition below
      //if(first_round==1){
      //  lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);
      // }else{
      SEXP s = sum_trees[parent2[i]];
      if(is<List>(s)){
        List sum_trees2=sum_trees[parent2[i]];
        List sum_trees_mat2=sum_trees_mat[parent2[i]];
        sum_trees2.push_back(proposal_tree[0]);
        sum_trees_mat2.push_back(proposal_tree[1]);
        //lik=sumtree_likelihood_function2_exact(y_scaled,sum_trees2,sum_trees_mat2,y_scaled.size(),a,nu,lambda);
        
        lik_list=sumtree_likelihood_function2_exact(y_scaled,sum_trees2,sum_trees_mat2,y_scaled.size(),a,nu,lambda);
        lik=as<double>(lik_list[0]);
        temppredoutput=lik_list[1];
        
        for(int t=0;t<sum_trees2.size();t++){
          NumericMatrix tree=sum_trees2[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          NumericMatrix mat=sum_trees_mat2[t];
          tree_prior*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
      }else{
        NumericMatrix sum_trees2=sum_trees[parent2[i]];
        NumericMatrix sum_trees_mat2=sum_trees_mat[parent2[i]];
        //other_int_nodes = find_term_nodes(sum_trees2);
        //p_other=other_int_nodes.size();
        List st(2);
        List st_mat(2);
        st[0]=sum_trees2;
        st[1]=proposal_tree[0];
        st_mat[0]=sum_trees_mat2;
        st_mat[1]=proposal_tree[1];
        // return(st);
        //lik=sumtree_likelihood_function2_exact(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        
        lik_list=sumtree_likelihood_function2_exact(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        lik=as<double>(lik_list[0]);
        temppredoutput=lik_list[1];
        
        for(int t=0;t<st.size();t++){
          NumericMatrix tree=st[t];
          NumericMatrix mat=st_mat[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
      }
      // }It should not be possible for get_best_split_exact to be used outside o
      
      
      // FIXED (I think) at the moment tree prior is only for current tree need to get it for entire sum of tree list.
      
      //int_nodes=find_term_nodes(proposal_tree[0]);
      //p=int_nodes.size()+p_other; //COULD ADD 1 for the variance parameter, and more numbers for other parameters. (This might influence probability-weights, but not the ranking of BICs)
      //BIC=-2*(lik+log(tree_prior))+p_other*log(data.n_rows);  
      BIC=-2*(lik+log(tree_prior));  
      if(BIC<lowest_BIC){
        lowest_BIC=BIC;
        //best_sv=split_var;
        //best_sp=split_point;
        //likeliest_tree=proposal_tree;
        tree_list[count]=proposal_tree[0];
        tree_mat_list[count]=proposal_tree[1];
        tree_lik[count]=BIC;
        //tree_parent[count]=parent;
        
        tree_preds[count]=temppredoutput;
        
        
        count++;
        if(count==(tree_list.size()-1)){
          list_size=list_size*2;
          tree_list=resize_bigger(tree_list,list_size);
          tree_mat_list=resize_bigger(tree_mat_list,list_size);
          tree_lik.resize(list_size);
          //tree_parent.resize(list_size);
          
          tree_preds=resize_bigger(tree_preds,list_size);
          
          
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){
          if(is<NumericMatrix>(proposal_tree[0])){
            //std::cout<<"its a matrix "<<"\n";
          }else{
            throw std::range_error("proposal tree not a matrix");
          }
          tree_list[count]=proposal_tree[0];
          tree_mat_list[count]=proposal_tree[1];
          tree_lik[count]=BIC;
          //tree_parent[count]=parent;
          
          tree_preds[count]=temppredoutput;
          
          
          count++;
          if(count==(tree_list.size()-1)){
            list_size=list_size*2;
            tree_list=resize_bigger(tree_list,list_size);
            tree_mat_list=resize_bigger(tree_mat_list,list_size);
            tree_lik.resize(list_size);
            //tree_parent.resize(list_size);
            
            tree_preds=resize_bigger(tree_preds,list_size);
            
            
          }
        }
      }
    }  
  }
  tree_list=resize(tree_list,count);
  tree_mat_list=resize(tree_mat_list,count);
  tree_lik.resize(count);
  //tree_parent.resize(count);
  IntegerVector tree_parent(count, parent);
  
  tree_preds=resize(tree_preds,count);
  
  // Rcout << "Line 3835. \n";
  
  if(count>0){
    if(less_greedy==1){
      ret[0]=tree_list;
      ret[1]=tree_lik;
      ret[2]=tree_mat_list;
      ret[3]=tree_parent;
      ret[4]=tree_preds;
      
      return (ret);	
    }
    //eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));
    //eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),tree_parent);
    eval_model=evaluate_model_occams_window_exact(wrap(tree_lik),lowest_BIC,c,wrap(tree_list),wrap(tree_mat_list),tree_parent,
                                                  tree_preds);
    NumericVector testlik =eval_model[0];
    List testtree =eval_model[1];    
    List testmat =eval_model[2]; 
    IntegerVector testpar =eval_model[3];
    List testpredlist=eval_model[4];
    // Rcout << "Line 3856. \n";
    
    if(testlik.size()>0){
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){
        IntegerVector owindices=orderforOW(testlik);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        List temp_opreds(maxOWsize);
        
        for(int t=0;t<maxOWsize;t++){
          temp_olik[t]=testlik[owindices[t]];
          temp_otrees[t]=testtree[owindices[t]];
          temp_omat[t]=testmat[owindices[t]];
          temp_oparent[t]=testpar[owindices[t]];
          temp_opreds[t]=testpredlist[owindices[t]];
          
        }
        testlik=temp_olik;
        testtree=temp_otrees;
        testmat=temp_omat;
        testpar=temp_oparent;
        testpredlist=temp_opreds;
        
      }
      // ret[0]=lowest_BIC;
      // ret[1]=best_sv;
      // ret[2]=best_sp;
      // ret[3]=likeliest_tree;
      // ret[4]=testtree;
      // ret[5]=testlik;
      // ret[6]=testmat;
      // ret[7]=testpar;
      // ret[8]=no_tree_err;
      
      ret[0]=testtree;
      ret[1]=testlik;
      ret[2]=testmat;
      ret[3]=testpar;
      ret[4]=testpredlist;
      return (ret);
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;
      List gr(1);
      gr[0]=no_tree_err;
      return(gr);
    }
  }else{
    no_tree_err=1;
    List gr(1);
    gr[0]=no_tree_err;
    return(gr);
  }
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericVector update_mean_var(NumericMatrix tree_table,NumericMatrix tree_matrix,NumericVector resids,double a){
  //// Rcout << "Get to start of update_mean_var. \n";
  //List update_params(1);
  NumericVector terminal_nodes;
  arma::uvec term_obs;
  terminal_nodes= find_term_nodes(tree_table);
  NumericVector Tj(terminal_nodes.size());
  NumericVector new_mean(terminal_nodes.size());
  //arma::vec armaresids=as<arma::vec>(resids);
  
  for(int k=0;k< terminal_nodes.size();k++){    
    term_obs=find_term_obs(tree_matrix,terminal_nodes[k]);
    //get the number of observations in node k
    Tj[k]=term_obs.size();
    NumericVector  get_mean(term_obs.size());
    for(int i=0;i<Tj[k];i++){
      get_mean[i]=resids[term_obs[i]];
    } 
    double sum_resids=std::accumulate(get_mean.begin(),get_mean.end(),0.0);
    new_mean[k]=sum_resids/(Tj[k]+a);
    arma::uvec temp;
    term_obs=temp;
  }  
  
  return(new_mean);
}
//######################################################################################################################//

// [[Rcpp::export]]
List update_predictions(NumericMatrix tree_table,NumericMatrix tree_matrix,NumericVector new_mean,int n){
  
  List updated_preds(2);
  NumericVector new_preds(n);
  NumericVector terminal_nodes;
  arma::uvec term_obs;
  terminal_nodes=find_term_nodes(tree_table);
  
  for(int k=0;k<terminal_nodes.size();k++){
    term_obs=find_term_obs(tree_matrix,terminal_nodes[k]);        
    //update the terminal node mean of the selected tree nodes:
    tree_table(terminal_nodes[k]-1,5)= new_mean[k];
    IntegerVector term_obs2=wrap(arma::conv_to<arma::ivec>::from(term_obs));
    new_preds[term_obs2]=new_mean[k];
  }
  updated_preds[0]=tree_table;
  updated_preds[1]=new_preds;
  
  return(updated_preds);
}
//######################################################################################################################//

using namespace Rcpp;
using namespace std;

const double flagval = __DBL_MIN__; 
inline double flag(double a, bool b) { return b ? a : flagval; }

// [[Rcpp::export]]
NumericVector subsetter(NumericVector a, LogicalVector b) {
  NumericVector a1=clone(a);
  transform(a1.begin(), a1.end(), b.begin(), a1.begin(), flag);
  NumericVector res = NumericVector(sum(b));  
  remove_copy(a1.begin(), a1.end(), res.begin(), flagval);
  
  return res;    
}
//######################################################################################################################//

// [[Rcpp::export]]
IntegerVector order_inc_(NumericVector x) {
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}

//######################################################################################################################//

// [[Rcpp::export]]
List min_which2(NumericVector array,int n,double minout,int whichout){
  // Function to find minimum of an array with n elements that is put in min
  minout=array[0];
  whichout=0;
  
  for(int i=1;i<n;i++){
    if(array[i]< minout){
      minout= array[i];
      whichout=i;
    }
  }
  List ret(2);
  ret[0]=minout;
  ret[1]=whichout;
  
  return(ret);
}
//######################################################################################################################//

#include <Rmath.h>
// [[Rcpp::export]]
double mll_meanvar2(double x, double x2, int n){
  double sigsq=(x2-((x*x)/n))/n;
  if(sigsq<=0){sigsq=0.00000000001;}
  
  return(n*(log(2*M_PI)+log(sigsq)+1)); /* M_PI is in Rmath.h  */
  
  //if(sigsq<=0){//sigsq=0.00000000001;
  //
  //return(n*(log(2*M_PI)+-26.32844)); /* M_PI is in Rmath.h  */
  //}else{
  //  return(n*(log(2*M_PI)+log(sigsq)+1)); /* M_PI is in Rmath.h  */
  //}
  
  //if(sigsq<=0){//sigsq=0.00000000001;
  //  
  //  return(n*(1.837877+-24.32844)); /* M_PI is in Rmath.h  */
  //}else{
  //  return(n*(1.837877+log(sigsq)+1)); /* M_PI is in Rmath.h  */
  //}
}
//######################################################################################################################//

// [[Rcpp::export]]
IntegerVector PELT_meanvar_norm2(NumericVector resp,double pen)
  // 0 by default, nonzero indicates error in code 
{
  int n=resp.size();
  NumericVector y2=cumsum(pow(resp,2));
  y2.push_front(0);
  NumericVector y=cumsum(resp);
  y.push_front(0);
  IntegerVector cptsout(n,0);    
  IntegerVector lastchangecpts(2*(n+1));
  NumericVector lastchangelike(n+1);
  IntegerVector checklist(n+1);  
  int nchecklist;
  double minout;
  NumericVector tmplike(n+1);
  IntegerVector tmpt(n+1);
  int tstar,i,whichout,nchecktmp;
  //double mll_meanvar();
  //void min_which();
  lastchangelike[0]= -pen;
  lastchangecpts[0]=0; lastchangecpts[n]=0;
  double x=y[1];
  double x2=y2[1];
  lastchangelike[1]=mll_meanvar2(x,x2,1);
  lastchangecpts[1]=0; lastchangecpts[n+1]=1;
  lastchangelike[2]=mll_meanvar2(y[2],y2[2],2);
  lastchangecpts[2]=0; lastchangecpts[n+2]=2;
  lastchangelike[3]=mll_meanvar2(y[3],y2[3],3);
  lastchangecpts[3]=0; lastchangecpts[n+3]=3;
  
  minout=lastchangelike[checklist[0]] + mll_meanvar2(x,x2,0)+pen;
  whichout=0;
  
  nchecklist=2;
  checklist[0]=0;
  checklist[1]=2;
  
  for(tstar=4;tstar<(n+1);tstar++){
    //R_CheckUserInterrupt(); // checks if user has interrupted the R session and quits if true 
    
    for(i=0;i<nchecklist;i++){
      tmplike[i]=lastchangelike[checklist[i]] + mll_meanvar2(y[tstar]- y[checklist[i]],y2[tstar]-y2[checklist[i]],tstar-checklist[i])+pen;
    }
    List mw=min_which2(tmplike,nchecklist,minout,whichout); //updates minout and whichout with min and which element 
    NumericVector tempmin=mw[0];
    minout=tempmin[0];
    lastchangelike[tstar]=minout;
    whichout=mw[1];
    lastchangecpts[tstar]=checklist[whichout]; 
    lastchangecpts[n+tstar]=tstar;
    
    // Update checklist for next iteration, first element is next tau 
    nchecktmp=0;
    for(i=0;i<nchecklist;i++){
      if(tmplike[i]<= (lastchangelike[tstar]+pen)){
        checklist[nchecktmp]=checklist[i];
        nchecktmp+=1;
      }
    }
    checklist[nchecktmp]=tstar-1;  // atleast 2 obs per seg
    nchecktmp+=1;
    nchecklist=nchecktmp;
  } // end taustar
  
  // put final set of changepoints together
  int ncpts=0;
  int last=n;
  while(last!=0){
    cptsout[ncpts]=lastchangecpts[n+last];
    last=lastchangecpts[last];
    ncpts+=1;
  }
  
  IntegerVector cptsoutret=cptsout[cptsout>0];
  std::sort(cptsoutret.begin(), cptsoutret.end());
  
  return(cptsoutret);
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double SS(arma::vec x, arma::vec y, double split){            
  double meanLeft, meanRight;
  int n = x.n_rows;
  arma::vec meanAll(n);
  arma::vec xLeft = y.elem(arma::find(x <= split));
  arma::vec xRight = y.elem(arma::find(x > split));
  meanLeft = mean(xLeft);
  meanRight = mean(xRight);
  
  for(int j=0;j<n;j++) {
    
    if(x[j]<=split){
      meanAll[j] = meanLeft;
    }else{
      meanAll[j] = meanRight;
    }  
  }
  arma::vec test_resids=y-meanAll;
  double tSS=as<double>(wrap(trans(y-meanAll)*(y-meanAll)));
  
  return tSS;
}
//######################################################################################################################//
// [[Rcpp::export]]

List gridCP(arma::vec x, arma::vec y, int gridSize = 10) {
  
  NumericVector out(gridSize-2);
  NumericVector cp_strength(gridSize-2);
  arma::vec no_split(gridSize-2);
  double xGrid, gridStep = (max(x)-min(x))/((double)gridSize-1.0);
  xGrid = min(x);
  List currSS(2);
  
  for(int i=1;i<(gridSize-1);i++) {
    xGrid += gridStep;
    arma::vec ld_size= y.elem(arma::find(x <= xGrid));
    arma::vec rd_size=y.elem(arma::find(x > xGrid));    
    
    if(ld_size.size()>2 && rd_size.size()>2)
    {
      out[i-1]=xGrid;
      double testSS=SS(x,y,xGrid);
      cp_strength[i-1] = testSS;
      no_split[i-1]=0;
    }else{
      no_split[i-1]=1;
    }
  }
  arma::uvec to_remove=find(no_split ==1);
  //IntegerVector remove_order_index=as<IntegerVector>(wrap(order_(as<NumericVector> (wrap(to_remove)))));
  //IntegerVector remove_order_index=order_(as<NumericVector> (wrap(arma::conv_to<arma::ivec>::from(to_remove))));
  //IntegerVector remove_order_index=wrap(arma::conv_to<arma::ivec>::from(to_remove));
  int num_to_remove = to_remove.size()-1;
  
  if(to_remove.size()>0){
    //for(int k=0;k<remove_order_index.size();k++){
    for(unsigned int k=0;k<to_remove.size();k++){
      out.erase(to_remove[num_to_remove-k]);
      cp_strength.erase(to_remove[num_to_remove-k]);
    }
  }
  if(out.size()>0){
    currSS[0]= out;
    currSS[1]= cp_strength;
    
    return currSS;
  }else{
    List ret(2);
    ret[0]=0;
    ret[1]=0;
    
    return ret;
  }
}
//######################################################################################################################//
// [[Rcpp::export]]

arma::field<arma::vec> gridCP_arma(arma::vec x, arma::vec y, int gridSize = 10) {
  
  arma::vec out(gridSize-2);
  arma::vec cp_strength(gridSize-2);
  arma::vec no_split(gridSize-2);
  double xGrid, gridStep = (max(x)-min(x))/((double)gridSize-1.0);
  xGrid = min(x);
  arma::field<arma::vec> currSS(2);
  
  for(int i=1;i<(gridSize-1);i++) {
    xGrid += gridStep;
    arma::vec ld_size= y.elem(arma::find(x <= xGrid));
    arma::vec rd_size=y.elem(arma::find(x > xGrid));    
    
    if(ld_size.size()>2 && rd_size.size()>2)
    {
      out[i-1]=xGrid;
      double testSS=SS(x,y,xGrid);
      cp_strength[i-1] = testSS;
      no_split[i-1]=0;
    }else{
      no_split[i-1]=1;
    }
  }
  arma::uvec to_remove=find(no_split ==1);
  //IntegerVector remove_order_index=as<IntegerVector>(wrap(order_(as<NumericVector> (wrap(to_remove)))));
  
  if(to_remove.size()>0){
    //perhaps should be ordered, but not sure why
    
    out=out(to_remove);
    cp_strength=cp_strength(to_remove);
    
    //out.shed_rows(to_remove);
    //cp_strength.shed_rows(to_remove);
    //for(int k=0;k<to_remove.size();k++){
    //out.erase(to_remove[remove_order_index[k]-1]);
    
    //cp_strength.erase(to_remove[remove_order_index[k]-1]);
    //}
  }
  if(out.size()>0){
    currSS[0]= out;
    currSS[1]= cp_strength;
    
    return currSS;
  }else{
    //List ret(2);
    arma::vec veczerotemp = {0}; ;
    currSS[0]= veczerotemp;
    currSS[1]= veczerotemp; 
    
    return currSS;
  }
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List make_gridpoint_cpmat(NumericMatrix data,NumericVector resp,int gridsize,int num_cp){
  int err1=0;
  int numvar = data.ncol();
  int numrows=gridsize*numvar-1;
  int numcols=3;
  arma::mat change_points(numrows,numcols);
  int row_count=-1;
  //IntegerVector sorted_var;
  
  arma::vec resp_ar=Rcpp::as<arma::vec>(resp);
  
  for(int i=0;i<data.ncol();i++){
    NumericVector coli=data(_,i);
    //arma::vec coli_ar=as<arma::vec>(coli);   
    
    //arma::vec t=Rcpp::as<arma::vec>(resp);
    List ans1 = gridCP(coli,resp_ar,gridsize);
    NumericVector ans2=ans1[0];
    NumericVector cp_strength=ans1[1];
    
    if(ans1.size()!=1 && ans2[0]!=0){
      for(int j=0;j<ans2.size();j++){
        row_count+=1;
        change_points(row_count,0)=i;
        change_points(row_count,1)=ans2[j];
        change_points(row_count,2)=cp_strength[j];
      }   
    }
  }
  if(row_count+1!=(int) change_points.n_rows){
    change_points.shed_rows(row_count+1,change_points.n_rows-1);
  }
  arma::vec te=change_points.col(2);
  //NumericVector col_to_order=as<NumericVector>(wrap(te));
  //IntegerVector ordered_dev=order_inc_(col_to_order);
  //ordered_dev=ordered_dev-1;
  change_points.shed_col(2);
  int cp=change_points.n_rows;
  double cp_prop=(double)num_cp/(double)100;
  int num_cp2=round(cp*(cp_prop));
  num_cp=round(num_cp2);
  
  if(num_cp==0 && cp!=0){
    num_cp=cp;
  }
  if(cp<num_cp){
    num_cp=change_points.n_rows;
  }
  if(num_cp==0){
    err1=1;
  }
  if(err1==0){
    //arma::uvec t=Rcpp::as<arma::uvec>(ordered_dev);
    arma::uvec t=arma::sort_index(te);
    
    t=t.subvec(0,num_cp-1);
    change_points=change_points.rows(t);
  }
  List ret(2);
  ret[0]=wrap(change_points);
  ret[1]=err1;
  
  return(ret);
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List make_gridpoint_cpmat_arma(arma::mat data,arma::vec resp,int gridsize,int num_cp){
  int err1=0;
  int numvar = data.n_cols;
  int numrows=gridsize*numvar-1;
  int numcols=3;
  arma::mat change_points(numrows,numcols);
  int row_count=-1;
  //IntegerVector sorted_var;
  
  for(int i=0;i<numvar;i++){
    //NumericVector coli=data(_,i);
    //arma::vec coli_ar=as<arma::vec>(coli);   
    //arma::vec t=Rcpp::as<arma::vec>(resp);
    arma::vec coli_ar=data.col(i);   
    //arma::vec t=Rcpp::as<arma::vec>(resp);
    arma::field<arma::vec> ans1 = gridCP_arma(coli_ar,resp,gridsize);
    
    //List ans1 = gridCP(coli,t,gridsize);
    arma::vec ans2=ans1[0];
    arma::vec cp_strength=ans1[1];
    
    if(ans1.size()!=1 && ans2[0]!=0){
      for(unsigned int j=0;j<ans2.size();j++){
        row_count+=1;
        change_points(row_count,0)=i;
        change_points(row_count,1)=ans2[j];
        change_points(row_count,2)=cp_strength[j];
      }   
    }
  }
  if(row_count+1!=(int) change_points.n_rows){
    change_points.shed_rows(row_count+1,change_points.n_rows-1);
  }
  arma::vec te=change_points.col(2);
  //NumericVector col_to_order=as<NumericVector>(wrap(te));
  //IntegerVector ordered_dev=order_inc_(col_to_order);
  //ordered_dev=ordered_dev-1;
  change_points.shed_col(2);
  int cp=change_points.n_rows;
  double cp_prop=(double)num_cp/(double)100;
  int num_cp2=round(cp*(cp_prop));
  num_cp=round(num_cp2);
  
  if(num_cp==0 && cp!=0){
    num_cp=cp;
  }
  if(cp<num_cp){
    num_cp=change_points.n_rows;
  }
  if(num_cp==0){
    err1=1;
  }
  if(err1==0){
    //arma::uvec t=Rcpp::as<arma::uvec>(ordered_dev);
    arma::uvec t=arma::sort_index(te);
    
    t=t.subvec(0,num_cp-1);
    change_points=change_points.rows(t);
  }
  List ret(2);
  ret[0]=wrap(change_points);
  ret[1]=err1;
  
  return(ret);
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List make_pelt_cpmat(NumericMatrix data,NumericVector resp,double pen,int num_cp){
  int err1=0;
  int n=data.nrow();
  arma::mat change_points;
  change_points.set_size(10000000,3);
  int row_count=-1;
  IntegerVector sorted_var;
  
  for(int i=0;i<data.ncol();i++){
    NumericVector coli=data(_,i);
    arma::vec coli_ar=as<arma::vec>(coli);
    sorted_var=order_inc_(coli);
    sorted_var=sorted_var-1;
    IntegerVector ans = PELT_meanvar_norm2(resp[sorted_var],pen);
    ans=ans-1;
    arma::vec resp_ar=as<arma::vec>(resp);
    //NumericVector allMean(n);
    arma::vec allMean_ar(n);
    
    if(ans.size()!=1 && ans[0]!=n){
      double meanRight;
      double meanLeft;      
      IntegerVector pelt_sp=sorted_var[ans];
      NumericVector split_values=coli[pelt_sp];
      NumericVector unique_sp=unique(split_values);
      //NumericVector all_sp=coli[pelt_sp];
      
      //if(unique_sp.size()<all_sp.size()){
      if(unique_sp.size()<split_values.size()){
        //IntegerVector index=match(unique_sp,all_sp);
        IntegerVector index=match(unique_sp,split_values);
        arma::ivec indexarma=as<arma::ivec>(index);
        IntegerVector t4=ifelse(is_na(index),1,0);
        arma::ivec t4arma=as<arma::ivec>(t4);
        arma::uvec t42=find(t4arma==0);
        IntegerVector index2(t42.size());
        
        arma::ivec index3=indexarma.elem(t42);
        
        index2=wrap(index3);
        index2=index2-1;
        ans=index2;
      }
      
      for(int j=0;j<ans.size()-1;j++){
        arma::vec xLeft;
        arma::vec xRight;
        double splitpoint;
        //if(unique_sp.size()<all_sp.size()){
        //  splitpoint = all_sp[j];
        if(unique_sp.size()<split_values.size()){
          splitpoint = split_values[j];
          xLeft = resp_ar.elem(arma::find(coli_ar<=splitpoint));
          xRight = resp_ar.elem(arma::find(coli_ar>splitpoint));
        }else{
          splitpoint = coli[sorted_var[ans[j]]];
          xLeft = resp_ar.elem(arma::find(coli_ar<=splitpoint));
          xRight = resp_ar.elem(arma::find(coli_ar>splitpoint));
        }
        if(xLeft.size()>2 && xRight.size()>2)
        {
          row_count +=1;
          meanLeft = mean(xLeft);
          meanRight = mean(xRight);
          arma::uvec left_ind=find(coli_ar<=splitpoint);
          arma::uvec right_ind=find(coli_ar>splitpoint);
          //NumericVector left_ind2=wrap(left_ind);
          //NumericVector right_ind2=wrap(right_ind);
          //allMean[left_ind2]=meanLeft;
          //allMean[right_ind2]=meanRight;
          //arma::vec allMean_ar=as<arma::vec>(allMean);
          
          allMean_ar.elem(left_ind).fill(meanLeft);
          allMean_ar.elem(right_ind).fill(meanRight);
          
          
          double cp_strength;
          cp_strength=as<double>(wrap(trans(resp_ar-allMean_ar)*(resp_ar-allMean_ar)));
          change_points(row_count,0)=i;
          change_points(row_count,1)=coli[sorted_var[ans[j]]];
          change_points(row_count,2)=cp_strength;
        }
      }          
    }
  }
  change_points.shed_rows(row_count+1,change_points.n_rows-1);
  arma::vec te=change_points.col(2);
  //NumericVector col_to_order=as<NumericVector>(wrap(te));
  //IntegerVector ordered_dev=order_inc_(col_to_order);
  //ordered_dev=ordered_dev-1;
  change_points.shed_col(2);
  int cp=change_points.n_rows;
  double cp_prop=(double)num_cp/(double)100;
  int num_cp2=round(cp*(cp_prop));
  num_cp=round(num_cp2);
  if(num_cp==0 && cp!=0){
    num_cp=cp;
  }
  if(cp<num_cp){
    num_cp=change_points.n_rows;
  }
  if(num_cp==0){
    err1=1;
  }
  if(err1==0){
    //arma::uvec t=Rcpp::as<arma::uvec>(ordered_dev);
    arma::uvec t=arma::sort_index(te);
    
    t=t.subvec(0,num_cp-1);
    change_points=change_points.rows(t);
  }
  List ret(2);
  ret[0]=wrap(change_points);
  ret[1]=err1;
  
  return(ret);
}
//###################################################################################//

// [[Rcpp::export]]

List get_best_trees(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                    arma::mat& D1,NumericMatrix resids,double a,double mu,double nu,double lambda,double c,
                    double sigma_mu,List tree_table,List tree_mat,double lowest_BIC,//int first_round,
                    IntegerVector parent,List cp_mat_list,//IntegerVector err_list,
                    NumericMatrix test_data,double alpha,double beta,bool is_test_data,double pen,int num_cp,bool split_rule_node,bool gridpoint,int maxOWsize,int num_splits,int gridsize, bool zero_split,
                    unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split
){
  List eval_model;
  NumericVector lik_list;
  List best_subset;
  int overall_size=1000;
  List overall_trees(overall_size);
  NumericVector overall_lik2;
  IntegerVector overall_parent2;
  List overall_mat(overall_size);
  std::vector<int> overall_parent(overall_size);
  std::vector<double> overall_lik(overall_size);
  int overall_count=0;  
  //Rcout << "Line 2800";
  if(zero_split==1){
    //Rcout << "Line 2802";
    
    //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
    overall_trees[0]= tree_table[0];
    overall_mat[0]= tree_mat[0];
    overall_parent[0]=-1;
    overall_parent2[0]=-1;
    double lik_temp=likelihood_function(resids,tree_table[0],tree_mat[0],a,mu,nu,lambda);
    double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
    //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
    double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp));
    
    overall_lik[0]= lowest_BIC_temp;
    //Rcout << "Next zero split tree lowest_BIC_temp = " << lowest_BIC_temp << ".\n"; 
    NumericVector templikvec(1);
    templikvec[0]=lowest_BIC_temp;
    overall_lik2=templikvec;
    
    
    overall_count=1;
    
  }
  //Rcout << "Line 1759 .\n";
  NumericVector test_preds;
  
  //for(int j=0;j<0;j++){
  
  for(int j=0;j<num_splits;j++){
    int lsize=1000;
    List table_subset_curr_round(lsize);
    std::vector<double> lik_subset_curr_round(lsize);
    List mat_subset_curr_round(lsize);
    std::vector<int> parent_curr_round(lsize);
    int count=0;
    for(int i=0;i<tree_table.size();i++){
      //NumericMatrix temp_list=cp_mat_list[0];
      //Rcout << "Line 2366. j = " << j << " i = "<<  i << ".\n";
      parent=-1;
      
      if(split_rule_node==1){
        if(j==0){
          best_subset=get_best_split(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                     lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
                                     min_num_obs_for_split,min_num_obs_after_split//,first_round
          );
        }else{
          best_subset=get_best_split(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                     lowest_BIC,parent[0],cp_mat_list[i],alpha,beta,maxOWsize,
                                     min_num_obs_for_split,min_num_obs_after_split//,first_round
          ); 
        }
      }else{
        best_subset=get_best_split(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                   lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
                                   min_num_obs_for_split,min_num_obs_after_split//,first_round
        ); 
      }
      
      
      
      
      //Rcout << "Line 2392. j = " << j << " i = "<<  i << ".\n";
      if(best_subset.size()==1){
        continue;
      }
      // List temp_trees=best_subset[4];
      // List temp_mat=best_subset[6];
      // lik_list=best_subset[5];
      // IntegerVector temp_parent=best_subset[7];
      List temp_trees=best_subset[0];
      List temp_mat=best_subset[2];
      lik_list=best_subset[1];
      IntegerVector temp_parent=best_subset[3];
      if(temp_parent.size()!= temp_trees.size()){
        throw std::range_error("there should be a parent for each tree!!!");
      }
      if(lik_list.size()==0){
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      if(min(lik_list)<lowest_BIC){
        lowest_BIC=min(lik_list);
      }
      //Rcout << "temp_parent" << temp_parent << " .\n";
      for(int k=0;k<temp_trees.size();k++){
        table_subset_curr_round[count]=temp_trees[k];
        lik_subset_curr_round[count]=lik_list[k];
        mat_subset_curr_round[count]=temp_mat[k];
        parent_curr_round[count]=temp_parent[k];
        count++;
        //Rcout << "Line 2419. j = " << j << " i = "<<  i << ".\n";
        if(count==(lsize-1)){
          lsize=lsize*2;
          table_subset_curr_round=resize_bigger(table_subset_curr_round,lsize);
          mat_subset_curr_round=resize_bigger(mat_subset_curr_round,lsize);
          lik_subset_curr_round.resize(lsize);
          parent_curr_round.resize(lsize);
        }
      }
    }
    table_subset_curr_round=resize(table_subset_curr_round,count);
    mat_subset_curr_round=resize(mat_subset_curr_round,count);
    lik_subset_curr_round.resize(count);
    parent_curr_round.resize(count);
    
    //NumericVector testparvec = wrap(parent_curr_round);
    //Rcout << "parent_curr_round" << testparvec << " .\n";
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0
      break;															// break out of the for-loop,
    }
    
    
    List  eval_modeltemp=evaluate_model_occams_window(as<NumericVector>(wrap(lik_subset_curr_round)),
                                                      lowest_BIC,
                                                      log(c),
                                                      table_subset_curr_round,
                                                      mat_subset_curr_round,
                                                      as<IntegerVector>(wrap(parent_curr_round)));
    
    
    lik_subset_curr_round=Rcpp::as<std::vector<double>>(eval_modeltemp[0]);
    table_subset_curr_round=eval_modeltemp[1];
    mat_subset_curr_round=eval_modeltemp[2];
    //overall_count=overall_trees.size();
    parent_curr_round=Rcpp::as<std::vector<int>>(eval_modeltemp[3]);
    
    
    
    
    if(table_subset_curr_round.size()==0){
      break;
    }
    
    for(int k=0;k<table_subset_curr_round.size();k++){
      overall_trees[overall_count]=table_subset_curr_round[k];
      overall_lik[overall_count]=lik_subset_curr_round[k];
      overall_mat[overall_count]=mat_subset_curr_round[k];
      overall_parent[overall_count]=parent_curr_round[k];
      overall_count++;
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
      }
    }
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    
    if(less_greedy==1){
      
    }else{
      eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
      overall_lik2=eval_model[0];
      overall_trees=eval_model[1];
      overall_mat=eval_model[2];
      overall_count=overall_trees.size();
      overall_parent2=eval_model[3];
      //add in check to see if OW accepted more than the top maxOW models...
      if(overall_lik2.size()>maxOWsize){
        IntegerVector owindices=orderforOW(overall_lik2);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        //Rcout << "Line 1849. j = " << j << ". \n";
        //now only select those elements
        for(int t=0;t<maxOWsize;t++){  
          temp_olik[t]=overall_lik2[owindices[t]];
          temp_otrees[t]=overall_trees[owindices[t]];
          temp_omat[t]= overall_mat[owindices[t]];
          temp_oparent[t]=overall_parent2[owindices[t]];
        }
        
        overall_lik2=temp_olik;
        overall_trees=temp_otrees;
        overall_mat=temp_omat;
        overall_count=overall_trees.size();
        overall_parent2=temp_oparent;
      }
      
    }
    
    tree_table=table_subset_curr_round;
    //IntegerVector temp1(table_subset_curr_round.size(),1);
    //err_list=temp1;
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }
    tree_mat=mat_subset_curr_round;
    parent=parent_curr_round;
    
    if(split_rule_node==1){
      NumericVector temp_preds;
      List updated_curr_preds;
      NumericVector new_mean;
      lowest_BIC=min(overall_lik2);
      
      
      //NumericMatrix curr_resids(resids.nrow(),resids.ncol());
      
      List temp(table_subset_curr_round.size());
      
      cp_mat_list=temp;
      
      
      //Rcout << "Line 2519. j = " << j << ". \n";
      //Rcout << "table_subset_curr_round.size() = " << table_subset_curr_round.size() << ". \n";
      //Rcout << "parent.size() = " << parent.size() << ". \n";
      //Rcout << "parent_curr_round.size() = " << parent_curr_round.size() << ". \n";
      //Rcout << "resids.ncol() = " << resids.ncol() << ". \n";
      //Rcout << "parent = " << parent << ". \n";
      
      for(int k=0;k<table_subset_curr_round.size();k++){
        //Rcout << "Line 2521. j = " << j << ". \n";
        
        //Rcout << "parent_curr_round[k] = " << parent_curr_round[k] << ". \n";
        
        
        NumericVector terminal_nodes;
        
        if(parent_curr_round[k]==-1){
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a);
        }else{    
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a);
        }  
        //Rcout << "Line 2533. j = " << j << ". \n";
        terminal_nodes=find_term_nodes(table_subset_curr_round[k]);
        updated_curr_preds=update_predictions(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,D1.n_rows);
        NumericVector test_res;
        
        
        
        
        if(parent_curr_round[k]==-1){
          test_res=resids(_,0); 
        }else{
          test_res=resids(_,parent_curr_round[k]);
        }
        
        NumericVector curr_test_res=updated_curr_preds[1];
        //Rcout << "Line 2548. j = " << j << ". \n";
        
        // if(parent_curr_round[k]==-1){
        //   curr_resids(_,0)=test_res-curr_test_res;
        // }else{
        //   curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;
        // }
        NumericVector temp_curr_resids=test_res-curr_test_res;
        
        
        List cp_mat_list1;
        if(gridpoint==0){
          cp_mat_list1=make_pelt_cpmat(wrap(D1),temp_curr_resids,pen,num_cp);
        }else{
          cp_mat_list1=make_gridpoint_cpmat(wrap(D1),temp_curr_resids,gridsize,num_cp);
        }
        
        cp_mat_list[k]=cp_mat_list1[0];
        
        //Rcout << "Line 2567. j = " << j << ". k = " << k << " . \n";
        
        
        
      }
      
      
      // List temp(0);
      // 
      // cp_mat_list=temp;
      // 
      // for(int f=0;f<curr_resids.ncol();f++){
      //   List cp_mat_list1;
      //   if(gridpoint==0){
      //     cp_mat_list1=make_pelt_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
      //   }else{
      //     cp_mat_list1=make_gridpoint_cpmat(wrap(D1),curr_resids(_,f),gridsize,num_cp);
      //   }
      //   
      //   cp_mat_list.push_back(cp_mat_list1[0]);      
      // }
      
    }  //end of if-statement split_rule_node==1
  }
  
  
  
  //Rcout << "overall_lik" << overall_lik.size() << " .\n";
  //Rcout << "overall_count" << overall_count << " .\n";
  //Rcout << "overall_trees" << overall_trees.size() << " .\n";
  
  overall_trees=resize(overall_trees,overall_count);
  overall_mat=resize(overall_mat,overall_count); 
  overall_lik.resize(overall_count);
  overall_parent.resize(overall_count);
  
  
  if(less_greedy==1){
    
    eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
    
  }
  //Rcout << "overall_lik" << overall_lik.size() << " .\n";
  //Rcout << "overall_trees" << overall_trees.size() << " .\n";
  //Rcout << "overall_count" << overall_count << " .\n";
  
  
  NumericVector temp_preds;
  List updated_preds;
  NumericVector new_mean;  
  NumericMatrix overall_test_preds(test_data.nrow(),overall_trees.size());  
  NumericMatrix overallpreds(D1.n_rows,overall_trees.size());
  lowest_BIC=min(overall_lik2);
  //Rcout << "Line 2596. \n";
  for(int k=0;k<overall_trees.size();k++){
    //NumericVector terminal_nodes;
    //Rcout << "Line 1947. k = " << k << ". \n";
    if(overall_parent2[k]==-1){
      //Rcout << "Line 1949. k = " << k << ". \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,0),a);
    }else{   
      //Rcout << "Line 1952. k = " << k << ". \n";
      //Rcout << "Line 1952. overall_parent2[k] = " << overall_parent2[k] << ". \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a);
    }
    //Rcout << "Line 1953. k = " << k << ". \n";
    //terminal_nodes=find_term_nodes(overall_trees[k]);
    updated_preds=update_predictions(overall_trees[k],overall_mat[k],new_mean,D1.n_rows);
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs(test_data,overall_trees[k]//,new_mean
    );
    temp_preds=updated_preds[1];
    overallpreds(_,k)=temp_preds;
    if(is_test_data)overall_test_preds(_,k)=test_preds;
    //Rcout << "Line 1961. k = " << k << ". \n";
  }
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
  //arma::colvec predicted_values=sum(M1,1);
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  //Rcout << "Line 3130. \n";
  
  
  //arma::colvec predicted_test_values=sum(M2,1);
  // List ret(8);
  // ret[0]=overall_lik2;
  // ret[1]=overall_trees;
  // ret[2]=overall_mat;
  // ret[3]=predicted_values;
  // ret[4]=overall_parent2;
  // ret[5]=wrap(M1);
  // ret[6]=lowest_BIC;
  // ret[7]=wrap(M2);
  List ret(7);
  ret[0]=overall_lik2;
  ret[1]=overall_trees;
  ret[2]=overall_mat;
  ret[3]=overall_parent2;
  ret[4]=wrap(M1);
  ret[5]=lowest_BIC;
  ret[6]=wrap(M2);
  return(ret);
}
//###################################################################################//

// [[Rcpp::export]]

List get_best_trees_update_splits(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                                  arma::mat& D1,NumericMatrix resids,double a,double mu,double nu,double lambda,double c,
                                  double sigma_mu,List tree_table,List tree_mat,double lowest_BIC,//int first_round,
                                  IntegerVector parent,List cp_mat_list,//IntegerVector err_list,
                                  NumericMatrix test_data,double alpha,double beta,bool is_test_data,double pen,int num_cp,bool split_rule_node,bool gridpoint,int maxOWsize,int num_splits,int gridsize, bool zero_split,
                                  unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split
){
  List eval_model;
  NumericVector lik_list;
  List best_subset;
  int overall_size=1000;
  List overall_trees(overall_size);
  NumericVector overall_lik2;
  IntegerVector overall_parent2;
  List overall_mat(overall_size);
  std::vector<int> overall_parent(overall_size);
  std::vector<double> overall_lik(overall_size);
  int overall_count=0;  
  //Rcout << "Line 1747";
  if(zero_split==1){
    //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
    overall_trees[0]= tree_table[0];
    overall_mat[0]= tree_mat[0];
    overall_parent[0]=-1;
    overall_parent2[0]=-1;
    double lik_temp=likelihood_function(resids,tree_table[0],tree_mat[0],a,mu,nu,lambda);
    double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
    //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
    double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp));
    
    overall_lik[0]= lowest_BIC_temp;
    
    
    NumericVector templikvec(1);
    templikvec[0]=lowest_BIC_temp;
    overall_lik2=templikvec;    
    
    
    overall_count=1;  
  }
  //Rcout << "Line 1759 .\n";
  NumericVector test_preds;
  //for(int j=0;j<0;j++){
  
  for(int j=0;j<num_splits;j++){
    int lsize=1000;
    List table_subset_curr_round(lsize);
    std::vector<double> lik_subset_curr_round(lsize);
    List mat_subset_curr_round(lsize);
    std::vector<int> parent_curr_round(lsize);
    int count=0;
    for(int i=0;i<tree_table.size();i++){
      //NumericMatrix temp_list=cp_mat_list[0];
      //Rcout << "Line 3136. j = " << j << " i = "<<  i << ".\n";
      parent=-1;
      
      if(split_rule_node==1){
        if(j==0){
          best_subset=get_best_split(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                     lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
                                     min_num_obs_for_split,min_num_obs_after_split//,first_round
          ); 
        }else{
          best_subset=get_best_split_2(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                       lowest_BIC,parent[0],cp_mat_list[i],alpha,beta,maxOWsize,
                                       min_num_obs_for_split,min_num_obs_after_split//,first_round
          ); 
        }
      }else{
        throw std::range_error("get_best_trees_update_splits should only apply when split_rule_node==1");
        
        best_subset=get_best_split_2(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                     lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
                                     min_num_obs_for_split,min_num_obs_after_split//,first_round
        ); 
      }
      
      
      
      
      //Rcout << "Line 2392. j = " << j << " i = "<<  i << ".\n";
      if(best_subset.size()==1){
        continue;
      }
      // List temp_trees=best_subset[4];
      // List temp_mat=best_subset[6];
      // lik_list=best_subset[5];
      // IntegerVector temp_parent=best_subset[7];
      List temp_trees=best_subset[0];
      List temp_mat=best_subset[2];
      lik_list=best_subset[1];
      IntegerVector temp_parent=best_subset[3];
      if(temp_parent.size()!= temp_trees.size()){
        throw std::range_error("there should be a parent for each tree!!!");
      }
      if(lik_list.size()==0){
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      if(min(lik_list)<lowest_BIC){
        lowest_BIC=min(lik_list);
      }
      //Rcout << "temp_parent" << temp_parent << " .\n";
      for(int k=0;k<temp_trees.size();k++){
        table_subset_curr_round[count]=temp_trees[k];
        lik_subset_curr_round[count]=lik_list[k];
        mat_subset_curr_round[count]=temp_mat[k];
        parent_curr_round[count]=temp_parent[k];
        count++;
        //Rcout << "Line 2419. j = " << j << " i = "<<  i << ".\n";
        if(count==(lsize-1)){
          lsize=lsize*2;
          table_subset_curr_round=resize_bigger(table_subset_curr_round,lsize);
          mat_subset_curr_round=resize_bigger(mat_subset_curr_round,lsize);
          lik_subset_curr_round.resize(lsize);
          parent_curr_round.resize(lsize);
        }
      }
    }
    table_subset_curr_round=resize(table_subset_curr_round,count);
    mat_subset_curr_round=resize(mat_subset_curr_round,count);
    lik_subset_curr_round.resize(count);
    parent_curr_round.resize(count);
    
    //NumericVector testparvec = wrap(parent_curr_round);
    //Rcout << "parent_curr_round" << testparvec << " .\n";
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0
      break;															// break out of the for-loop,
    }
    
    
    List  eval_modeltemp=evaluate_model_occams_window(as<NumericVector>(wrap(lik_subset_curr_round)),
                                                      lowest_BIC,
                                                      log(c),
                                                      table_subset_curr_round,
                                                      mat_subset_curr_round,
                                                      as<IntegerVector>(wrap(parent_curr_round)));
    
    
    lik_subset_curr_round=Rcpp::as<std::vector<double>>(eval_modeltemp[0]);
    table_subset_curr_round=eval_modeltemp[1];
    mat_subset_curr_round=eval_modeltemp[2];
    //overall_count=overall_trees.size();
    parent_curr_round=Rcpp::as<std::vector<int>>(eval_modeltemp[3]);
    
    
    
    
    
    if(table_subset_curr_round.size()==0){
      break;
    }
    
    
    for(int k=0;k<table_subset_curr_round.size();k++){
      overall_trees[overall_count]=table_subset_curr_round[k];
      overall_lik[overall_count]=lik_subset_curr_round[k];
      overall_mat[overall_count]=mat_subset_curr_round[k];
      overall_parent[overall_count]=parent_curr_round[k];
      overall_count++;
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
      }
    }
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    
    if(less_greedy==1){
      
    }else{
      eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
      overall_lik2=eval_model[0];
      overall_trees=eval_model[1];
      overall_mat=eval_model[2];
      overall_count=overall_trees.size();
      overall_parent2=eval_model[3];
      //add in check to see if OW accepted more than the top maxOW models...
      if(overall_lik2.size()>maxOWsize){
        IntegerVector owindices=orderforOW(overall_lik2);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        //Rcout << "Line 1849. j = " << j << ". \n";
        //now only select those elements
        for(int t=0;t<maxOWsize;t++){  
          temp_olik[t]=overall_lik2[owindices[t]];
          temp_otrees[t]=overall_trees[owindices[t]];
          temp_omat[t]= overall_mat[owindices[t]];
          temp_oparent[t]=overall_parent2[owindices[t]];
        }
        
        overall_lik2=temp_olik;
        overall_trees=temp_otrees;
        overall_mat=temp_omat;
        overall_count=overall_trees.size();
        overall_parent2=temp_oparent;
      }
      
    }
    
    
    tree_table=table_subset_curr_round;
    //IntegerVector temp1(table_subset_curr_round.size(),1);
    //err_list=temp1;
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }
    tree_mat=mat_subset_curr_round;
    parent=parent_curr_round;
    
    if(split_rule_node==1){
      NumericVector temp_preds;
      List updated_curr_preds;
      NumericVector new_mean;
      lowest_BIC=min(overall_lik2);
      
      
      //NumericMatrix curr_resids(resids.nrow(),resids.ncol());
      
      List temp(table_subset_curr_round.size());
      
      cp_mat_list=temp;
      
      
      //Rcout << "Line 2519. j = " << j << ". \n";
      //Rcout << "table_subset_curr_round.size() = " << table_subset_curr_round.size() << ". \n";
      //Rcout << "parent.size() = " << parent.size() << ". \n";
      //Rcout << "parent_curr_round.size() = " << parent_curr_round.size() << ". \n";
      //Rcout << "resids.ncol() = " << resids.ncol() << ". \n";
      //Rcout << "parent = " << parent << ". \n";
      
      for(int k=0;k<table_subset_curr_round.size();k++){
        //Rcout << "Line 2521. j = " << j << ". \n";
        
        //Rcout << "parent_curr_round[k] = " << parent_curr_round[k] << ". \n";
        
        
        //NumericVector terminal_nodes;
        
        if(parent_curr_round[k]==-1){
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a);
        }else{    
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a);
        }  
        //Rcout << "Line 2533. j = " << j << ". \n";
        NumericVector terminal_nodes=find_term_nodes(table_subset_curr_round[k]);
        updated_curr_preds=update_predictions(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,D1.n_rows);
        NumericVector test_res;
        
        
        
        
        if(parent_curr_round[k]==-1){
          test_res=resids(_,0); 
        }else{
          test_res=resids(_,parent_curr_round[k]);
        }
        
        //NumericVector curr_test_res=updated_curr_preds[1];
        //Rcout << "Line 2548. j = " << j << ". \n";
        
        // if(parent_curr_round[k]==-1){
        //   curr_resids(_,0)=test_res-curr_test_res;
        // }else{
        //   curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;
        // }
        //NumericVector temp_curr_resids=test_res-curr_test_res;
        
        //Loop over terminal nodes
        List Tempcpmatlist(terminal_nodes.size());
        for(int nodeind=0;nodeind<terminal_nodes.size();nodeind++){
          //Rcout << "Line 3369. j = " << j << ". k = " << k << " . \n";
          
          arma::uvec term_obsarma=find_term_obs(mat_subset_curr_round[k],terminal_nodes[nodeind]);
          //IntegerVector term_obs=wrap(term_obsarma);
          IntegerVector term_obs=wrap(arma::conv_to<arma::ivec>::from(term_obsarma));
          
          //Rcout << "Line 3373. j = " << j << ". k = " << k << " . \n";
          //Rcout << "term_obs = " << term_obs << " . \n";
          //Rcout << "term_obs.size() = " << term_obs.size() << " . \n";
          //Rcout << "temp_curr_resids.size() = " << temp_curr_resids.size() << " . \n";
          
          NumericVector tempsubset = test_res[term_obs];
          
          
          List cp_mat_list1;
          arma::mat tempdata_subset = D1.rows(term_obsarma);
          if(gridpoint==0){
            //cp_mat_list1=make_pelt_cpmat(wrap(tempdata_subset),temp_curr_resids[term_obs],pen,num_cp);
            cp_mat_list1=make_pelt_cpmat(wrap(tempdata_subset),tempsubset,pen,num_cp);
          }else{
            //cp_mat_list1=make_gridpoint_cpmat(wrap(tempdata_subset),temp_curr_resids[term_obs],gridsize,num_cp);
            cp_mat_list1=make_gridpoint_cpmat(wrap(tempdata_subset),tempsubset,gridsize,num_cp);
          }
          //Rcout << "Line 3386. j = " << j << ". k = " << k << " . \n";
          
          Tempcpmatlist[nodeind]=cp_mat_list1[0];
          //Rcout << "Line 3389. j = " << j << ". k = " << k << " . \n";
          
        }
        //Rcout << "Line 3384. j = " << j << ". k = " << k << " . \n";
        
        cp_mat_list[k]=Tempcpmatlist;
        //Rcout << "Line 3395. j = " << j << ". k = " << k << " . \n";
        
        
        //cp_mat_list[k]=cp_mat_list1[0];
        
        //Rcout << "Line 3387. j = " << j << ". k = " << k << " . \n";
        
        
        
      }
      
      
      // List temp(0);
      // 
      // cp_mat_list=temp;
      // 
      // for(int f=0;f<curr_resids.ncol();f++){
      //   List cp_mat_list1;
      //   if(gridpoint==0){
      //     cp_mat_list1=make_pelt_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
      //   }else{
      //     cp_mat_list1=make_gridpoint_cpmat(wrap(D1),curr_resids(_,f),gridsize,num_cp);
      //   }
      //   
      //   cp_mat_list.push_back(cp_mat_list1[0]);      
      // }
      
    }//end of if-statement split_rule_node==1
  }
  overall_trees=resize(overall_trees,overall_count);
  overall_mat=resize(overall_mat,overall_count); 
  overall_lik.resize(overall_count);
  overall_parent.resize(overall_count);
  
  if(less_greedy==1){
    
    eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
  }
  
  
  
  
  NumericVector temp_preds;
  List updated_preds;
  NumericVector new_mean;  
  NumericMatrix overall_test_preds(test_data.nrow(),overall_trees.size());  
  NumericMatrix overallpreds(D1.n_rows,overall_trees.size());
  lowest_BIC=min(overall_lik2);
  //Rcout << "Line 3434. \n";
  for(int k=0;k<overall_trees.size();k++){
    //NumericVector terminal_nodes;
    //Rcout << "Line 1947. k = " << k << ". \n";
    if(overall_parent2[k]==-1){
      //Rcout << "Line 1949. k = " << k << ". \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,0),a);
    }else{   
      //Rcout << "Line 1952. k = " << k << ". \n";
      //Rcout << "Line 1952. overall_parent2[k] = " << overall_parent2[k] << ". \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a);
    }
    //Rcout << "Line 1953. k = " << k << ". \n";
    //terminal_nodes=find_term_nodes(overall_trees[k]);
    updated_preds=update_predictions(overall_trees[k],overall_mat[k],new_mean,D1.n_rows);
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs(test_data,overall_trees[k]//,new_mean
    );
    temp_preds=updated_preds[1];
    overallpreds(_,k)=temp_preds;
    if(is_test_data)overall_test_preds(_,k)=test_preds;
    //Rcout << "Line 1961. k = " << k << ". \n";
  }
  //Rcout << "Line 3457. \n";
  
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
  //arma::colvec predicted_values=sum(M1,1);
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  //arma::colvec predicted_test_values=sum(M2,1);
  // List ret(8);
  // ret[0]=overall_lik2;
  // ret[1]=overall_trees;
  // ret[2]=overall_mat;
  // ret[3]=predicted_values;
  // ret[4]=overall_parent2;
  // ret[5]=wrap(M1);
  // ret[6]=lowest_BIC;
  // ret[7]=wrap(M2);
  List ret(7);
  ret[0]=overall_lik2;
  ret[1]=overall_trees;
  ret[2]=overall_mat;
  ret[3]=overall_parent2;
  ret[4]=wrap(M1);
  ret[5]=lowest_BIC;
  ret[6]=wrap(M2);
  return(ret);
}
//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_sum(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                        arma::mat& D1,NumericMatrix resids,double a,double mu,double nu,double lambda,
                        double c,double sigma_mu,List tree_table,List tree_mat,double lowest_BIC,//int first_round,
                        IntegerVector parent,List cp_mat_list,IntegerVector err_list,NumericMatrix test_data,double alpha,double beta,bool is_test_data,double pen,int num_cp,bool split_rule_node,bool gridpoint,int maxOWsize,List prev_sum_trees,List prev_sum_trees_mat,NumericVector y_scaled,int num_splits,int gridsize,bool zero_split,
                        unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split
){
  //Rcout << "Get to start of get_best_trees_sum. \n";
  
  List eval_model;
  NumericVector lik_list;
  List best_subset;
  int overall_size=1000;
  List overall_trees(overall_size);
  NumericVector overall_lik2;
  IntegerVector overall_parent2;
  List overall_mat(overall_size);
  int overall_count=0;  
  std::vector<int> overall_parent(overall_size);
  std::vector<double> overall_lik(overall_size);
  NumericVector test_preds;
  
  
  //////   //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  if(zero_split==1){
    for(int q=0; q<parent.size();q++){
      //Rcout << "q= "<< q << ". \n";
      //Rcout << "length of parent = "<< parent.size() << ". \n";
      //Rcout << "length of prev_sum_trees = "<< prev_sum_trees.size() << ". \n";
      
      // if(parent.size()!= prev_sum_trees.size() )throw std::range_error("length of parent !=length of prev_sum_trees ");
      
      SEXP s_temp = prev_sum_trees[parent[q]];
      //Rcout <<"line 1999.\n";
      
      if(is<List>(s_temp)){
        //Rcout << "s is a list. \n";
        List sum_trees2_temp=prev_sum_trees[parent[q]];
        List sum_trees_mat2_temp=prev_sum_trees_mat[parent[q]];
        sum_trees2_temp.push_back(tree_table[0]);
        sum_trees_mat2_temp.push_back(tree_mat[0]);
        double lik_temp=sumtree_likelihood_function2(y_scaled,sum_trees2_temp,sum_trees_mat2_temp,y_scaled.size(),a,nu,lambda);  
        double tree_prior_temp=1;
        //int p_other=0;
        //NumericVector other_int_nodes;
        //Rcout << "Get to loop over t. \n";
        
        for(int t=0;t<sum_trees2_temp.size();t++){
          NumericMatrix tree=sum_trees2_temp[t];
          NumericMatrix mat=sum_trees_mat2_temp[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior_temp*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
        //Rcout << "Finish Loop. \n";
        
        //double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other)*log(D1.n_rows);  
        double BIC=-2*(lik_temp+log(tree_prior_temp));  
        
        ////Rcout << "Get to fill in overall_. \n";
        
        overall_trees[overall_count]= tree_table[0];
        overall_mat[overall_count]= tree_mat[0];
        overall_parent[overall_count]=parent[q];
        //Rcout << "parent[q] = " << parent[q] <<  ".\n";
        //Rcout << "When q=" << q << " parent[q]=" << parent[q] << ". overall_count =" << overall_count << ".\n";
        
        //double lik_temp=likelihood_function(resids[0],tree_table[0],tree_mat[0],a,mu,nu,lambda);
        //double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
        //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
        overall_lik[overall_count]= BIC;
        //Rcout << "Get to end of adding no split trees. \n";
        
        overall_count++;
      }else{
        //Rcout << "s is not a list. \n";
        NumericMatrix sum_trees2_temp=prev_sum_trees[parent[q]];
        NumericMatrix sum_trees_mat2=prev_sum_trees_mat[parent[q]];
        List st(2);
        List st_mat(2);
        st[0]=sum_trees2_temp;
        st[1]=tree_table[0];
        st_mat[0]=sum_trees_mat2;
        st_mat[1]=tree_mat[0];
        // return(st);
        double lik_temp=sumtree_likelihood_function2(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);  
        double tree_prior_temp=1;
        //int p_other=0;
        //NumericVector other_int_nodes;
        for(int t=0;t<st.size();t++){
          NumericMatrix tree=st[t];
          NumericMatrix mat=st_mat[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior_temp*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
        
        //double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other)*log(D1.n_rows);  
        double BIC=-2*(lik_temp+log(tree_prior_temp));  
        
        
        overall_trees[overall_count]= tree_table[0];
        overall_mat[overall_count]= tree_mat[0];
        overall_parent[overall_count]=parent[q];
        //Rcout << "When q=" << q << " parent[q]=" << parent[q] << ". overall_count =" << overall_count << ".\n";
        //double lik_temp=likelihood_function(resids[0],tree_table[0],tree_mat[0],a,mu,nu,lambda);
        //double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
        //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
        overall_lik[overall_count]= BIC;
        overall_count++;
        //Rcout << "Get to end of adding no split trees. \n";
        
      }
      //Rcout <<"line 2073.\n";
      //Rcout << "overall_count = " << overall_count << ".\n";
      //Rcout << "overall_size = " << overall_size << ".\n";
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        //Rcout <<"line 2081.\n";
        
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
      }
    }
    //int overall_count=0; 
    //Rcout <<"line 2082.\n";
    //Rcout << "overall_count = " << overall_count << ".\n";
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    
    
    
    eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
    //Rcout <<"line 2094.\n";
    
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=order_(overall_lik2);
      owindices=owindices-1;
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);
      List temp_otrees(maxOWsize);
      List temp_omat(maxOWsize);
      IntegerVector temp_oparent(maxOWsize);
      //Rcout <<"line 2106.\n";
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){  
        temp_olik[t]=overall_lik2[owindices[t]];
        temp_otrees[t]=overall_trees[owindices[t]];
        temp_omat[t]= overall_mat[owindices[t]];
        temp_oparent[t]=overall_parent2[owindices[t]];
      }
      
      overall_lik2=temp_olik;
      overall_trees=temp_otrees;
      overall_mat=temp_omat;
      overall_count=overall_trees.size();
      overall_parent2=temp_oparent;
    }
    //Rcout <<"line 2122.\n";
    
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }
  }
  
  //Rcout <<"Get past zero split tree block of code.\n";
  
  
  for(int j=0;j<num_splits;j++){
    int lsize=1000;
    List table_subset_curr_round(lsize);
    std::vector<double> lik_subset_curr_round(lsize);
    List mat_subset_curr_round(lsize);
    std::vector<int> parent_curr_round(lsize);
    int count=0;
    //Rcout <<"LENGTH OF TREE TABLE LIST = " << tree_table.size() << " !!!!!!!!!.\n";
    for(int i=0;i<tree_table.size();i++){
      // if(first_round==1){
      //   parent=-1;
      //   //NumericMatrix temp_list=cp_mat_list[0];
      //   best_subset=get_best_split(spike_tree,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
      //                              lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
      //                              min_num_obs_for_split,min_num_obs_after_split//,first_round
      //                                ); 
      //   
      // }else{
      if(err_list[i]==0){
        //NumericMatrix test_tree=tree_table[i];
        //NumericMatrix test_treemat=tree_mat[i];
        //NumericMatrix test_cpmat= cp_mat_list[parent[i]];
        //need to append current tree_table[i] to its parent sum_of_trees   
        
        
        if(split_rule_node==1){
          if(j==0){
            best_subset=get_best_split_sum(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                           parent[i],cp_mat_list[i],
                                                                alpha,beta,maxOWsize,//first_round,
                                                                prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                min_num_obs_for_split,min_num_obs_after_split);   
          }else{
            best_subset=get_best_split_sum(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                           parent[i],cp_mat_list[i],
                                                                alpha,beta,maxOWsize,//first_round,
                                                                prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                min_num_obs_for_split,min_num_obs_after_split);
          }
        }else{
          //Rcout << " line 3757. no update .\n";
          
          //Rcout << " i = " << i << ".\n";
          //Rcout << " parent[i] = " << parent[i] << ".\n";
          
          best_subset=get_best_split_sum(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                         parent[i],cp_mat_list[parent[i]],
                                                              alpha,beta,maxOWsize,//first_round,
                                                              prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                              min_num_obs_for_split,min_num_obs_after_split);
        }
        
        
        
        // return(best_subset);
      }else if(err_list[i]==1){
        //Rcout << "CONTINUE. error list.  \n";
        continue;
      }else{
        List ret_list(6);
        ret_list[0]=9999;
        ret_list[1]=err_list[i];
        ret_list[2]=i;
        ret_list[3]=j;
        ret_list[4]=tree_table;
        ret_list[5]=err_list;
        return(ret_list);
        throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");
      }
      //}
      //Rcout << "Get past get_best_split. \n";
      
      if(best_subset.size()==1){
        //Rcout << "CONTINUE. \n";
        continue;
      }
      // List temp_trees=best_subset[4];
      // List temp_mat=best_subset[6];
      // lik_list=best_subset[5];
      // IntegerVector temp_parent=best_subset[7];
      List temp_trees=best_subset[0];
      List temp_mat=best_subset[2];
      lik_list=best_subset[1];
      IntegerVector temp_parent=best_subset[3];
      if(temp_parent.size()!= temp_trees.size()){
        throw std::range_error("there should be a parent for each tree!!!");
      }
      if(lik_list.size()==0){
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      
      if(min(lik_list)<lowest_BIC){
        lowest_BIC=min(lik_list);
      }
      for(int k=0;k<temp_trees.size();k++){
        table_subset_curr_round[count]=temp_trees[k];
        lik_subset_curr_round[count]=lik_list[k];
        mat_subset_curr_round[count]=temp_mat[k];
        parent_curr_round[count]=temp_parent[k];
        count++;
        
        if(count==(lsize-1)){
          lsize=lsize*2;
          table_subset_curr_round=resize_bigger(table_subset_curr_round,lsize);
          mat_subset_curr_round=resize_bigger(mat_subset_curr_round,lsize);
          lik_subset_curr_round.resize(lsize);
          parent_curr_round.resize(lsize);
        }
      }
    }
    table_subset_curr_round=resize(table_subset_curr_round,count);
    mat_subset_curr_round=resize(mat_subset_curr_round,count);
    lik_subset_curr_round.resize(count);
    parent_curr_round.resize(count);
    
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0
      break;															// break out of the for-loop,
    }
    
    
    List  eval_modeltemp=evaluate_model_occams_window(as<NumericVector>(wrap(lik_subset_curr_round)),
                                                      lowest_BIC,
                                                      log(c),
                                                      table_subset_curr_round,
                                                      mat_subset_curr_round,
                                                      as<IntegerVector>(wrap(parent_curr_round)));
    
    
    lik_subset_curr_round=Rcpp::as<std::vector<double>>(eval_modeltemp[0]);
    table_subset_curr_round=eval_modeltemp[1];
    mat_subset_curr_round=eval_modeltemp[2];
    //overall_count=overall_trees.size();
    parent_curr_round=Rcpp::as<std::vector<int>>(eval_modeltemp[3]);
    
    
    
    if(table_subset_curr_round.size()==0){
      //Rcout << "BREAK. \n";
      break;
    }
    //Rcout << "inner loop round WITHOUT CONTINUE OR BREAK. \n";
    
    for(int k=0;k<table_subset_curr_round.size();k++){
      overall_trees[overall_count]=table_subset_curr_round[k];
      overall_lik[overall_count]=lik_subset_curr_round[k];
      overall_mat[overall_count]=mat_subset_curr_round[k];
      overall_parent[overall_count]=parent_curr_round[k];
      overall_count++;
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
      }
    }
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    
    
    if(less_greedy==1){
      
    }else{
      //Rcout << "overall_parent[0] BEFORE OW EVALUATION = " << overall_parent[0] << ".\n";
      eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
      overall_lik2=eval_model[0];
      overall_trees=eval_model[1];
      overall_mat=eval_model[2];
      overall_count=overall_trees.size();
      overall_parent2=eval_model[3];
      //Rcout << "overall_parent2[0] AFTER OW EVALUATION" << overall_parent2[0] << ".\n";
      //Rcout << "overall_parent2[0] AFTER OW EVALUATION" << overall_parent2[0] << ".\n";
      
      //add in check to see if OW accepted more than the top maxOW models...
      if(overall_lik2.size()>maxOWsize){
        //Rcout << "MORE THAN MAXOWSIZE!!!!!!!!!!!" << overall_parent2[0] << ".\n";
        
        //find the maxOWsize best models and continue with those!
        IntegerVector owindices=orderforOW(overall_lik2);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        
        //now only select those elements
        for(int t=0;t<maxOWsize;t++){  
          temp_olik[t]=overall_lik2[owindices[t]];
          temp_otrees[t]=overall_trees[owindices[t]];
          temp_omat[t]= overall_mat[owindices[t]];
          temp_oparent[t]=overall_parent2[owindices[t]];
        }
        
        overall_lik2=temp_olik;
        overall_trees=temp_otrees;
        overall_mat=temp_omat;
        overall_count=overall_trees.size();
        overall_parent2=temp_oparent;
        //Rcout << "overall_parent2[0] AFTER REMOVING EXTRA MODELS AND REARRAGING = " << overall_parent2[0] << ".\n";
        
      }
      
    }
    
    
    
    tree_table=table_subset_curr_round;
    //IntegerVector temp1(table_subset_curr_round.size(),1);
    //err_list=temp1;
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }
    tree_mat=mat_subset_curr_round;
    parent=parent_curr_round;
    
    if(split_rule_node==1){
      NumericVector temp_preds;
      List updated_curr_preds;
      NumericVector new_mean;
      lowest_BIC=min(overall_lik2);
      
      
      //NumericMatrix curr_resids(resids.nrow(),resids.ncol());
      
      
      List temp(table_subset_curr_round.size());
      
      cp_mat_list=temp;
      
      
      
      for(int k=0;k<table_subset_curr_round.size();k++){
        //NumericVector terminal_nodes;
        
        if(parent_curr_round[k]==-1){
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a);
        }else{    
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a);
        }  
        
        //terminal_nodes=find_term_nodes(table_subset_curr_round[k]);
        updated_curr_preds=update_predictions(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,D1.n_rows);
        NumericVector test_res;
        
        if(parent_curr_round[k]==-1){
          test_res=resids(_,0); 
        }else{
          test_res=resids(_,parent_curr_round[k]);
        }
        
        NumericVector curr_test_res=updated_curr_preds[1];
        //Rcout << "Line 2548. j = " << j << ". \n";
        
        // if(parent_curr_round[k]==-1){
        //   curr_resids(_,0)=test_res-curr_test_res;
        // }else{
        //   curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;
        // }
        NumericVector temp_curr_resids=test_res-curr_test_res;
        
        List cp_mat_list1;
        if(gridpoint==0){
          cp_mat_list1=make_pelt_cpmat(wrap(D1),temp_curr_resids,pen,num_cp);
        }else{
          cp_mat_list1=make_gridpoint_cpmat(wrap(D1),temp_curr_resids,gridsize,num_cp);
        }
        
        cp_mat_list[k]=cp_mat_list1[0];
        
      }
      
      
      
      //List temp(0);
      
      //cp_mat_list=temp;
      
      // for(int f=0;f<curr_resids.ncol();f++){
      //   List cp_mat_list1;
      //   if(gridpoint==0){
      //     cp_mat_list1=make_pelt_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
      //   }else{
      //     cp_mat_list1=make_gridpoint_cpmat(wrap(D1),curr_resids(_,f),gridsize,num_cp);
      //   }
      //   
      //   cp_mat_list.push_back(cp_mat_list1[0]);      
      // }
      
    }   //end of if-statement split_rule_node==1
  }
  //Rcout << "Get to end of outer loop. \n";
  //Rcout << "overall_trees.size()= " << overall_trees.size() << " \n";
  //Rcout << "overall_mat.size()= " << overall_mat.size() << " \n";
  //Rcout << "overall_lik.size()= " << overall_lik.size() << " \n";
  //Rcout << "overall_parent.size()= " << overall_parent.size() << " \n";
  //Rcout << "overall_parent2.size()= " << overall_parent2.size() << " \n";
  
  overall_trees=resize(overall_trees,overall_count);
  overall_mat=resize(overall_mat,overall_count); 
  overall_lik.resize(overall_count);
  overall_parent.resize(overall_count);
  
  
  if(less_greedy==1){
    
    eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
  }
  
  
  
  NumericVector temp_preds;
  List updated_preds;
  NumericVector new_mean;  
  NumericMatrix overall_test_preds(test_data.nrow(),overall_trees.size());  
  NumericMatrix overallpreds(D1.n_rows,overall_trees.size());
  lowest_BIC=min(overall_lik2);
  //Rcout << "Get to start of update mean loop. \n";
  //Rcout << "overall_trees.size()= " << overall_trees.size() << " \n";
  //Rcout << "overall_mat.size()= " << overall_mat.size() << " \n";
  //Rcout << "overall_lik.size()= " << overall_lik.size() << " \n";
  //Rcout << "overall_parent.size()= " << overall_parent.size() << " \n";
  //Rcout << "overall_parent2.size()= " << overall_parent2.size() << " \n";
  
  for(int k=0;k<overall_trees.size();k++){
    //NumericVector terminal_nodes;
    
    if(overall_parent2[k]==-1){
      //Rcout << "within update mean var. overall_parent2[k]= " << overall_parent2[k] << " \n";
      //Rcout << "k= " << k << " \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,0),a);
    }else{
      //Rcout << "within update mean var. overall_parent2[k]= " << overall_parent2[k] << " \n";
      //Rcout << "k= " << k << " \n";
      
      //Rcout << "a= " << a << " \n";
      SEXP tree_checktype = overall_trees[k];
      if(is<NumericMatrix>(tree_checktype)){
        //Rcout << "tree_checktype is a NumericMatrix. \n";
      }else{
        //Rcout << "tree_checktype is NOT a NumericMatrix. \n";
      }
      SEXP mat_checktype = overall_mat[k];
      if(is<NumericMatrix>(mat_checktype)){
        //Rcout << "mat_checktype is a NumericMatrix. \n";
      }else{
        //Rcout << "mat_checktype is NOT a NumericMatrix. \n";
      }
      //  SEXP vec_checktype = resids(_,overall_parent2[k]);
      //  if(is<NumericVector>(vec_checktype)){
      //Rcout << "vec_checktype is a NumericMatrix. \n";
      //  }else{
      //Rcout << "vec_checktype is NOT a NumericMatrix. \n";
      //  }
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a);
    }
    //Rcout << "Get past update mean_var. k= " << k << " \n";
    
    //terminal_nodes=find_term_nodes(overall_trees[k]);
    updated_preds=update_predictions(overall_trees[k],overall_mat[k],new_mean,D1.n_rows);
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs(test_data,overall_trees[k]//,new_mean
    );
    temp_preds=updated_preds[1];
    overallpreds(_,k)=temp_preds;
    if(is_test_data)overall_test_preds(_,k)=test_preds;
  }
  //Rcout << "Get to end of update mean loop. \n";
  
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
  //arma::colvec predicted_values=sum(M1,1);
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  //arma::colvec predicted_test_values=sum(M2,1);
  //Rcout << "Get to end of get_est_trees_sum. \n";
  // List ret(8);
  // ret[0]=overall_lik2;
  // ret[1]=overall_trees;
  // ret[2]=overall_mat;
  // ret[3]=predicted_values;
  // ret[4]=overall_parent2;
  // ret[5]=wrap(M1);
  // ret[6]=lowest_BIC;
  // ret[7]=wrap(M2);
  List ret(7);
  ret[0]=overall_lik2;
  ret[1]=overall_trees;
  ret[2]=overall_mat;
  ret[3]=overall_parent2;
  ret[4]=wrap(M1);
  ret[5]=lowest_BIC;
  ret[6]=wrap(M2);
  return(ret);
}
//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_sum_update_splits(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                                      arma::mat& D1,NumericMatrix resids,double a,double mu,double nu,double lambda,
                                      double c,double sigma_mu,List tree_table,List tree_mat,double lowest_BIC,//int first_round,
                                      IntegerVector parent,List cp_mat_list,IntegerVector err_list,NumericMatrix test_data,double alpha,double beta,bool is_test_data,double pen,int num_cp,bool split_rule_node,bool gridpoint,int maxOWsize,List prev_sum_trees,List prev_sum_trees_mat,NumericVector y_scaled,int num_splits,int gridsize,bool zero_split,
                                      unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split
){
  //Rcout << "Get to start of get_best_trees_sum. \n";
  
  List eval_model;
  NumericVector lik_list;
  List best_subset;
  int overall_size=1000;
  List overall_trees(overall_size);
  NumericVector overall_lik2;
  IntegerVector overall_parent2;
  List overall_mat(overall_size);
  int overall_count=0;  
  std::vector<int> overall_parent(overall_size);
  std::vector<double> overall_lik(overall_size);
  NumericVector test_preds;
  
  
  //////   //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  if(zero_split==1){
    for(int q=0; q<parent.size();q++){
      //Rcout << "q= "<< q << ". \n";
      //Rcout << "length of parent = "<< parent.size() << ". \n";
      //Rcout << "length of prev_sum_trees = "<< prev_sum_trees.size() << ". \n";
      
      // if(parent.size()!= prev_sum_trees.size() )throw std::range_error("length of parent !=length of prev_sum_trees ");
      
      SEXP s_temp = prev_sum_trees[parent[q]];
      //Rcout <<"line 1999.\n";
      
      if(is<List>(s_temp)){
        //Rcout << "s is a list. \n";
        List sum_trees2_temp=prev_sum_trees[parent[q]];
        List sum_trees_mat2_temp=prev_sum_trees_mat[parent[q]];
        sum_trees2_temp.push_back(tree_table[0]);
        sum_trees_mat2_temp.push_back(tree_mat[0]);
        double lik_temp=sumtree_likelihood_function2(y_scaled,sum_trees2_temp,sum_trees_mat2_temp,y_scaled.size(),a,nu,lambda);  
        double tree_prior_temp=1;
        //int p_other=0;
        //NumericVector other_int_nodes;
        //Rcout << "Get to loop over t. \n";
        
        for(int t=0;t<sum_trees2_temp.size();t++){
          NumericMatrix tree=sum_trees2_temp[t];
          NumericMatrix mat=sum_trees_mat2_temp[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior_temp*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
        //Rcout << "Finish Loop. \n";
        
        //double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other)*log(D1.n_rows);  
        double BIC=-2*(lik_temp+log(tree_prior_temp));  
        
        ////Rcout << "Get to fill in overall_. \n";
        
        overall_trees[overall_count]= tree_table[0];
        overall_mat[overall_count]= tree_mat[0];
        overall_parent[overall_count]=parent[q];
        //Rcout << "parent[q] = " << parent[q] <<  ".\n";
        //Rcout << "When q=" << q << " parent[q]=" << parent[q] << ". overall_count =" << overall_count << ".\n";
        
        //double lik_temp=likelihood_function(resids[0],tree_table[0],tree_mat[0],a,mu,nu,lambda);
        //double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
        //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
        overall_lik[overall_count]= BIC;
        //Rcout << "Get to end of adding no split trees. \n";
        
        overall_count++;
      }else{
        //Rcout << "s is not a list. \n";
        NumericMatrix sum_trees2_temp=prev_sum_trees[parent[q]];
        NumericMatrix sum_trees_mat2=prev_sum_trees_mat[parent[q]];
        List st(2);
        List st_mat(2);
        st[0]=sum_trees2_temp;
        st[1]=tree_table[0];
        st_mat[0]=sum_trees_mat2;
        st_mat[1]=tree_mat[0];
        // return(st);
        double lik_temp=sumtree_likelihood_function2(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);  
        double tree_prior_temp=1;
        //int p_other=0;
        //NumericVector other_int_nodes;
        for(int t=0;t<st.size();t++){
          NumericMatrix tree=st[t];
          NumericMatrix mat=st_mat[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior_temp*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
        
        //double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other)*log(D1.n_rows);  
        double BIC=-2*(lik_temp+log(tree_prior_temp));  
        
        
        overall_trees[overall_count]= tree_table[0];
        overall_mat[overall_count]= tree_mat[0];
        overall_parent[overall_count]=parent[q];
        //Rcout << "When q=" << q << " parent[q]=" << parent[q] << ". overall_count =" << overall_count << ".\n";
        //double lik_temp=likelihood_function(resids[0],tree_table[0],tree_mat[0],a,mu,nu,lambda);
        //double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
        //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
        overall_lik[overall_count]= BIC;
        overall_count++;
        //Rcout << "Get to end of adding no split trees. \n";
        
      }
      //Rcout <<"line 2073.\n";
      //Rcout << "overall_count = " << overall_count << ".\n";
      //Rcout << "overall_size = " << overall_size << ".\n";
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        //Rcout <<"line 2081.\n";
        
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
      }
    }
    //int overall_count=0; 
    //Rcout <<"line 2082.\n";
    //Rcout << "overall_count = " << overall_count << ".\n";
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
    //Rcout <<"line 2094.\n";
    
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=order_(overall_lik2);
      owindices=owindices-1;
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);
      List temp_otrees(maxOWsize);
      List temp_omat(maxOWsize);
      IntegerVector temp_oparent(maxOWsize);
      //Rcout <<"line 2106.\n";
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){  
        temp_olik[t]=overall_lik2[owindices[t]];
        temp_otrees[t]=overall_trees[owindices[t]];
        temp_omat[t]= overall_mat[owindices[t]];
        temp_oparent[t]=overall_parent2[owindices[t]];
      }
      
      overall_lik2=temp_olik;
      overall_trees=temp_otrees;
      overall_mat=temp_omat;
      overall_count=overall_trees.size();
      overall_parent2=temp_oparent;
    }
    //Rcout <<"line 2122.\n";
    
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }
  }
  
  //Rcout <<"get sum trees update splits. Get past zero split tree block of code.\n";
  
  
  for(int j=0;j<num_splits;j++){
    int lsize=1000;
    List table_subset_curr_round(lsize);
    std::vector<double> lik_subset_curr_round(lsize);
    List mat_subset_curr_round(lsize);
    std::vector<int> parent_curr_round(lsize);
    int count=0;
    //Rcout << "begin loop over splits. j == " << j << ".\n";
    
    //Rcout <<"LENGTH OF TREE TABLE LIST = " << tree_table.size() << " !!!!!!!!!.\n";
    for(int i=0;i<tree_table.size();i++){
      //Rcout << "begin loop over models. j == " << j << ". i == "<< i << ".\n";
      
      // if(first_round==1){
      //   parent=-1;
      //   //NumericMatrix temp_list=cp_mat_list[0];
      //   best_subset=get_best_split(spike_tree,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
      //                              lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
      //                              min_num_obs_for_split,min_num_obs_after_split//,first_round
      //                                ); 
      //   
      // }else{
      if(err_list[i]==0){
        //NumericMatrix test_tree=tree_table[i];
        //NumericMatrix test_treemat=tree_mat[i];
        //NumericMatrix test_cpmat= cp_mat_list[parent[i]];
        //need to append current tree_table[i] to its parent sum_of_trees   
        
        
        if(split_rule_node==1){
          if(j==0){
            best_subset=get_best_split_sum(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                           parent[i],cp_mat_list[parent[i]],
                                                                alpha,beta,maxOWsize,//first_round,
                                                                prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                min_num_obs_for_split,min_num_obs_after_split);   
          }else{
            //Rcout << " line 4302. with update .\n";
            
            //Rcout << " i = " << i << ".\n";
            //Rcout << " parent[i] = " << parent[i] << ".\n";
            best_subset=get_best_split_sum_2(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                             parent[i],cp_mat_list[i],
                                                                  alpha,beta,maxOWsize,//first_round,
                                                                  prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                  min_num_obs_for_split,min_num_obs_after_split);
          }
        }else{
          throw std::range_error("get_best_trees_update_splits should only apply when split_rule_node==1");
          best_subset=get_best_split_sum(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                         parent[i],cp_mat_list[i],
                                                              alpha,beta,maxOWsize,//first_round,
                                                              prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                              min_num_obs_for_split,min_num_obs_after_split);
        }
        
        
        
        // return(best_subset);
      }else if(err_list[i]==1){
        //Rcout << "j == << " << j << "  \n";
        
        throw std::range_error("err_list[i] is 0");
        
        //Rcout << "CONTINUE. error list.  \n";
        continue;
      }else{
        List ret_list(6);
        ret_list[0]=9999;
        ret_list[1]=err_list[i];
        ret_list[2]=i;
        ret_list[3]=j;
        ret_list[4]=tree_table;
        ret_list[5]=err_list;
        //return(ret_list);
        throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");
      }
      //}
      //Rcout << "Get past get_best_split. j == " << j << ".\n";
      
      if(best_subset.size()==1){
        //Rcout << "CONTINUE. \n";
        continue;
      }
      // List temp_trees=best_subset[4];
      // List temp_mat=best_subset[6];
      // lik_list=best_subset[5];
      // IntegerVector temp_parent=best_subset[7];
      List temp_trees=best_subset[0];
      List temp_mat=best_subset[2];
      lik_list=best_subset[1];
      IntegerVector temp_parent=best_subset[3];
      if(temp_parent.size()!= temp_trees.size()){
        throw std::range_error("there should be a parent for each tree!!!");
      }
      if(lik_list.size()==0){
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      
      if(min(lik_list)<lowest_BIC){
        lowest_BIC=min(lik_list);
      }
      for(int k=0;k<temp_trees.size();k++){
        table_subset_curr_round[count]=temp_trees[k];
        lik_subset_curr_round[count]=lik_list[k];
        mat_subset_curr_round[count]=temp_mat[k];
        parent_curr_round[count]=temp_parent[k];
        count++;
        
        if(count==(lsize-1)){
          lsize=lsize*2;
          table_subset_curr_round=resize_bigger(table_subset_curr_round,lsize);
          mat_subset_curr_round=resize_bigger(mat_subset_curr_round,lsize);
          lik_subset_curr_round.resize(lsize);
          parent_curr_round.resize(lsize);
        }
      }
    }
    table_subset_curr_round=resize(table_subset_curr_round,count);
    mat_subset_curr_round=resize(mat_subset_curr_round,count);
    lik_subset_curr_round.resize(count);
    parent_curr_round.resize(count);
    
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0
      break;															// break out of the for-loop,
    }
    
    
    List  eval_modeltemp=evaluate_model_occams_window(as<NumericVector>(wrap(lik_subset_curr_round)),
                                                      lowest_BIC,
                                                      log(c),
                                                      table_subset_curr_round,
                                                      mat_subset_curr_round,
                                                      as<IntegerVector>(wrap(parent_curr_round)));
    
    
    lik_subset_curr_round=Rcpp::as<std::vector<double>>(eval_modeltemp[0]);
    table_subset_curr_round=eval_modeltemp[1];
    mat_subset_curr_round=eval_modeltemp[2];
    //overall_count=overall_trees.size();
    parent_curr_round=Rcpp::as<std::vector<int>>(eval_modeltemp[3]);
    
    
    
    
    if(table_subset_curr_round.size()==0){
      //Rcout << "BREAK. \n";
      break;
    }
    //Rcout << "inner loop round WITHOUT CONTINUE OR BREAK. \n";
    
    for(int k=0;k<table_subset_curr_round.size();k++){
      overall_trees[overall_count]=table_subset_curr_round[k];
      overall_lik[overall_count]=lik_subset_curr_round[k];
      overall_mat[overall_count]=mat_subset_curr_round[k];
      overall_parent[overall_count]=parent_curr_round[k];
      overall_count++;
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
      }
    }
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    //Rcout << "overall_parent[0] BEFORE OW EVALUATION = " << overall_parent[0] << ".\n";
    
    
    if(less_greedy==1){
      
    }else{
      eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
      overall_lik2=eval_model[0];
      overall_trees=eval_model[1];
      overall_mat=eval_model[2];
      overall_count=overall_trees.size();
      overall_parent2=eval_model[3];
      //Rcout << "overall_parent2[0] AFTER OW EVALUATION" << overall_parent2[0] << ".\n";
      //Rcout << "overall_parent2[0] AFTER OW EVALUATION" << overall_parent2[0] << ".\n";
      
      //add in check to see if OW accepted more than the top maxOW models...
      if(overall_lik2.size()>maxOWsize){
        //Rcout << "MORE THAN MAXOWSIZE!!!!!!!!!!!" << overall_parent2[0] << ".\n";
        
        //find the maxOWsize best models and continue with those!
        IntegerVector owindices=orderforOW(overall_lik2);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        
        //now only select those elements
        for(int t=0;t<maxOWsize;t++){  
          temp_olik[t]=overall_lik2[owindices[t]];
          temp_otrees[t]=overall_trees[owindices[t]];
          temp_omat[t]= overall_mat[owindices[t]];
          temp_oparent[t]=overall_parent2[owindices[t]];
        }
        
        overall_lik2=temp_olik;
        overall_trees=temp_otrees;
        overall_mat=temp_omat;
        overall_count=overall_trees.size();
        overall_parent2=temp_oparent;
        //Rcout << "overall_parent2[0] AFTER REMOVING EXTRA MODELS AND REARRAGING = " << overall_parent2[0] << ".\n";
        
      }
      
    }
    
    tree_table=table_subset_curr_round;
    IntegerVector temp1(table_subset_curr_round.size(),0);
    err_list=temp1;
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
    }
    tree_mat=mat_subset_curr_round;
    parent=parent_curr_round;
    
    if(split_rule_node==1){
      NumericVector temp_preds;
      List updated_curr_preds;
      NumericVector new_mean;
      lowest_BIC=min(overall_lik2);
      
      
      //NumericMatrix curr_resids(resids.nrow(),resids.ncol());
      
      
      List temp(table_subset_curr_round.size());
      
      cp_mat_list=temp;
      
      
      
      for(int k=0;k<table_subset_curr_round.size();k++){
        //NumericVector terminal_nodes;
        
        if(parent_curr_round[k]==-1){
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a);
        }else{    
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a);
        }  
        
        NumericVector terminal_nodes=find_term_nodes(table_subset_curr_round[k]);
        updated_curr_preds=update_predictions(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,D1.n_rows);
        NumericVector test_res;
        
        if(parent_curr_round[k]==-1){
          test_res=resids(_,0); 
        }else{
          test_res=resids(_,parent_curr_round[k]);
        }
        
        //NumericVector curr_test_res=updated_curr_preds[1];
        //Rcout << "Line 2548. j = " << j << ". \n";
        
        // if(parent_curr_round[k]==-1){
        //   curr_resids(_,0)=test_res-curr_test_res;
        // }else{
        //   curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;
        // }
        //NumericVector temp_curr_resids=test_res-curr_test_res;
        
        //Loop over terminal nodes
        List Tempcpmatlist(terminal_nodes.size());
        for(int nodeind=0;nodeind<terminal_nodes.size();nodeind++){
          
          //IntegerVector term_obs=wrap(find_term_obs(mat_subset_curr_round[k],terminal_nodes[nodeind]));
          arma::uvec term_obsarma=find_term_obs(mat_subset_curr_round[k],terminal_nodes[nodeind]);
          //IntegerVector term_obs=wrap(term_obsarma);
          IntegerVector term_obs=wrap(arma::conv_to<arma::ivec>::from(term_obsarma));
          
          //Rcout << "term_obs = " << term_obs << ".\n";
          //Rcout << "test_res= " << test_res << ".\n";
          
          //Rcout << "test_res[term_obs] = " << test_res[term_obs] << ".\n";
          NumericVector tempsubset = test_res[term_obs];
          //Rcout << "tempsubset = " << tempsubset << ".\n";
          
          
          List cp_mat_list1;
          arma::mat tempdata_subset = D1.rows(term_obsarma);
          if(gridpoint==0){
            //cp_mat_list1=make_pelt_cpmat(wrap(tempdata_subset),temp_curr_resids[term_obs],pen,num_cp);
            cp_mat_list1=make_pelt_cpmat(wrap(tempdata_subset),tempsubset,pen,num_cp);
          }else{
            //cp_mat_list1=make_gridpoint_cpmat(wrap(tempdata_subset),temp_curr_resids[term_obs],gridsize,num_cp);
            cp_mat_list1=make_gridpoint_cpmat(wrap(tempdata_subset),tempsubset,gridsize,num_cp);
          }
          Tempcpmatlist[nodeind]=cp_mat_list1[0];
        }
        cp_mat_list[k]=Tempcpmatlist;
        
      }
      
      
      
      //List temp(0);
      
      //cp_mat_list=temp;
      
      // for(int f=0;f<curr_resids.ncol();f++){
      //   List cp_mat_list1;
      //   if(gridpoint==0){
      //     cp_mat_list1=make_pelt_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
      //   }else{
      //     cp_mat_list1=make_gridpoint_cpmat(wrap(D1),curr_resids(_,f),gridsize,num_cp);
      //   }
      //   
      //   cp_mat_list.push_back(cp_mat_list1[0]);      
      // }
      
    }   //end of if-statement split_rule_node==1
  }
  //Rcout << "Get to end of outer loop. \n";
  //Rcout << "overall_trees.size()= " << overall_trees.size() << " \n";
  //Rcout << "overall_mat.size()= " << overall_mat.size() << " \n";
  //Rcout << "overall_lik.size()= " << overall_lik.size() << " \n";
  //Rcout << "overall_parent.size()= " << overall_parent.size() << " \n";
  //Rcout << "overall_parent2.size()= " << overall_parent2.size() << " \n";
  
  overall_trees=resize(overall_trees,overall_count);
  overall_mat=resize(overall_mat,overall_count); 
  overall_lik.resize(overall_count);
  overall_parent.resize(overall_count);
  
  
  if(less_greedy==1){
    
    eval_model=evaluate_model_occams_window(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
  }
  
  NumericVector temp_preds;
  List updated_preds;
  NumericVector new_mean;  
  NumericMatrix overall_test_preds(test_data.nrow(),overall_trees.size());  
  NumericMatrix overallpreds(D1.n_rows,overall_trees.size());
  lowest_BIC=min(overall_lik2);
  //Rcout << "Get to start of update mean loop. \n";
  //Rcout << "overall_trees.size()= " << overall_trees.size() << " \n";
  //Rcout << "overall_mat.size()= " << overall_mat.size() << " \n";
  //Rcout << "overall_lik.size()= " << overall_lik.size() << " \n";
  //Rcout << "overall_parent.size()= " << overall_parent.size() << " \n";
  //Rcout << "overall_parent2.size()= " << overall_parent2.size() << " \n";
  
  for(int k=0;k<overall_trees.size();k++){
    //NumericVector terminal_nodes;
    
    if(overall_parent2[k]==-1){
      //Rcout << "within update mean var. overall_parent2[k]= " << overall_parent2[k] << " \n";
      //Rcout << "k= " << k << " \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,0),a);
    }else{
      //Rcout << "within update mean var. overall_parent2[k]= " << overall_parent2[k] << " \n";
      //Rcout << "k= " << k << " \n";
      
      //Rcout << "a= " << a << " \n";
      SEXP tree_checktype = overall_trees[k];
      if(is<NumericMatrix>(tree_checktype)){
        //Rcout << "tree_checktype is a NumericMatrix. \n";
      }else{
        //Rcout << "tree_checktype is NOT a NumericMatrix. \n";
      }
      SEXP mat_checktype = overall_mat[k];
      if(is<NumericMatrix>(mat_checktype)){
        //Rcout << "mat_checktype is a NumericMatrix. \n";
      }else{
        //Rcout << "mat_checktype is NOT a NumericMatrix. \n";
      }
      //  SEXP vec_checktype = resids(_,overall_parent2[k]);
      //  if(is<NumericVector>(vec_checktype)){
      //Rcout << "vec_checktype is a NumericMatrix. \n";
      //  }else{
      //Rcout << "vec_checktype is NOT a NumericMatrix. \n";
      //  }
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a);
    }
    //Rcout << "Get past update mean_var. k= " << k << " \n";
    
    //terminal_nodes=find_term_nodes(overall_trees[k]);
    updated_preds=update_predictions(overall_trees[k],overall_mat[k],new_mean,D1.n_rows);
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs(test_data,overall_trees[k]//,new_mean
    );
    temp_preds=updated_preds[1];
    overallpreds(_,k)=temp_preds;
    if(is_test_data)overall_test_preds(_,k)=test_preds;
  }
  //Rcout << "Get to end of update mean loop. \n";
  
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
  //arma::colvec predicted_values=sum(M1,1);
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  //arma::colvec predicted_test_values=sum(M2,1);
  //Rcout << "Get to end of get_est_trees_sum. \n";
  // List ret(8);
  // ret[0]=overall_lik2;
  // ret[1]=overall_trees;
  // ret[2]=overall_mat;
  // ret[3]=predicted_values;
  // ret[4]=overall_parent2;
  // ret[5]=wrap(M1);
  // ret[6]=lowest_BIC;
  // ret[7]=wrap(M2);
  List ret(7);
  ret[0]=overall_lik2;
  ret[1]=overall_trees;
  ret[2]=overall_mat;
  ret[3]=overall_parent2;
  ret[4]=wrap(M1);
  ret[5]=lowest_BIC;
  ret[6]=wrap(M2);
  return(ret);
}
//###################################################################################//

// [[Rcpp::export]]

List get_best_trees_exact(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                          arma::mat& D1,NumericMatrix resids,double a,double mu,double nu,double lambda,double c,
                          double sigma_mu,List tree_table,List tree_mat,double lowest_BIC,//int first_round,
                          IntegerVector parent,List cp_mat_list,//IntegerVector err_list,
                          NumericMatrix test_data,double alpha,double beta,bool is_test_data,double pen,int num_cp,bool split_rule_node,bool gridpoint,int maxOWsize,int num_splits,int gridsize, bool zero_split,
                          unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split
){
  List eval_model;
  NumericVector lik_list;
  List best_subset;
  int overall_size=1000;
  List overall_trees(overall_size);
  NumericVector overall_lik2;
  IntegerVector overall_parent2;
  List overall_mat(overall_size);
  std::vector<int> overall_parent(overall_size);
  std::vector<double> overall_lik(overall_size);
  
  List overall_predvecs(overall_size);
  
  
  int overall_count=0;  
  // Rcout << "Line 6275";
  if(zero_split==1){
    // Rcout << "Line 3016";
    
    //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
    overall_trees[0]= tree_table[0];
    overall_mat[0]= tree_mat[0];
    overall_parent[0]=-1;
    overall_parent2[0]=-1;
    //double lik_temp=likelihood_function(resids,tree_table[0],tree_mat[0],a,mu,nu,lambda);
    //Rcout << "Line 3024";
    
    List lik_templist=likelihood_function2_exact(resids,tree_table[0],tree_mat[0],a,mu,nu,lambda);
    //Rcout << "Line 3027";
    
    double lik_temp= as<double>(lik_templist[0]);
    NumericVector temp_predvec=lik_templist[1];
    
    // Rcout << "Line 6156";
    // Rcout << "Line 6307. lik_temp= " << lik_temp << ".\n";
    NumericMatrix tabletemptest = tree_table[0];
    // Rcout << "Line 6307. num rows= " <<tabletemptest.nrow() << ".\n";
    
    
    double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
    // Rcout << "Line 6309. tree_prior_temp= " << tree_prior_temp << ".\n";
    //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
    double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp));
    // Rcout << "Line 6311. lowest_BIC_temp= " << lowest_BIC_temp << ".\n";
    overall_lik[0]= lowest_BIC_temp;
    //Rcout << "Next zero split tree lowest_BIC_temp = " << lowest_BIC_temp << ".\n"; 
    NumericVector templikvec(1);
    templikvec[0]=lowest_BIC_temp;
    overall_lik2=templikvec;
    
    overall_predvecs[0]=temp_predvec;
    
    overall_count=1;
    
  }
  
  // Rcout << "Line 6311 .\n";
  NumericVector test_preds;
  
  //for(int j=0;j<0;j++){
  
  for(int j=0;j<num_splits;j++){
    int lsize=1000;
    List table_subset_curr_round(lsize);
    std::vector<double> lik_subset_curr_round(lsize);
    List mat_subset_curr_round(lsize);
    std::vector<int> parent_curr_round(lsize);
    
    List predvecs_curr_round(lsize);
    
    
    int count=0;
    // Rcout << "Line 6327. tree_table.size() = " << tree_table.size() << ".\n";
    
    for(int i=0;i<tree_table.size();i++){
      //NumericMatrix temp_list=cp_mat_list[0];
      // Rcout << "Line 6331. j = " << j << " i = "<<  i << ".\n";
      
      parent=-1;
      
      if(split_rule_node==1){
        if(j==0){
          best_subset=get_best_split_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                           lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
                                           min_num_obs_for_split,min_num_obs_after_split//,first_round
          ); 
        }else{
          best_subset=get_best_split_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                           lowest_BIC,parent[0],cp_mat_list[i],alpha,beta,maxOWsize,
                                           min_num_obs_for_split,min_num_obs_after_split//,first_round
          ); 
        }
      }else{
        best_subset=get_best_split_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                         lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
                                         min_num_obs_for_split,min_num_obs_after_split//,first_round
        ); 
      }
      
      
      
      
      // Rcout << "Line 6357. j = " << j << " i = "<<  i << ".\n";
      
      if(best_subset.size()==1){
        continue;
      }
      // List temp_trees=best_subset[4];
      // List temp_mat=best_subset[6];
      // lik_list=best_subset[5];
      // IntegerVector temp_parent=best_subset[7];
      List temp_trees=best_subset[0];
      List temp_mat=best_subset[2];
      lik_list=best_subset[1];
      IntegerVector temp_parent=best_subset[3];
      
      List temp_predlist=best_subset[4];
      
      // Rcout << "Line 6372. j = " << j << " i = "<<  i << ".\n";
      
      if(temp_parent.size()!= temp_trees.size()){
        throw std::range_error("there should be a parent for each tree!!!");
      }
      if(lik_list.size()==0){
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      if(min(lik_list)<lowest_BIC){
        lowest_BIC=min(lik_list);
      }
      //Rcout << "temp_parent" << temp_parent << " .\n";
      for(int k=0;k<temp_trees.size();k++){
        table_subset_curr_round[count]=temp_trees[k];
        lik_subset_curr_round[count]=lik_list[k];
        mat_subset_curr_round[count]=temp_mat[k];
        parent_curr_round[count]=temp_parent[k];
        predvecs_curr_round[count]=temp_predlist[k];
        
        count++;
        
        // Rcout << "Line 6392. j = " << j << " i = "<<  i << ".\n";
        
        if(count==(lsize-1)){
          lsize=lsize*2;
          table_subset_curr_round=resize_bigger(table_subset_curr_round,lsize);
          mat_subset_curr_round=resize_bigger(mat_subset_curr_round,lsize);
          lik_subset_curr_round.resize(lsize);
          parent_curr_round.resize(lsize);
          predvecs_curr_round=resize_bigger(predvecs_curr_round,lsize);
          
        }
      }
    }
    
    // Rcout << "Line 6406. j = " << j <<  ".\n";
    
    table_subset_curr_round=resize(table_subset_curr_round,count);
    mat_subset_curr_round=resize(mat_subset_curr_round,count);
    lik_subset_curr_round.resize(count);
    parent_curr_round.resize(count);
    predvecs_curr_round=resize(predvecs_curr_round,count);
    
    //NumericVector testparvec = wrap(parent_curr_round);
    //Rcout << "parent_curr_round" << testparvec << " .\n";
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0
      break;															// break out of the for-loop,
    }
    
    List eval_modeltemp=evaluate_model_occams_window_exact(as<NumericVector>(wrap(lik_subset_curr_round)),
                                                           lowest_BIC,
                                                           log(c),
                                                           table_subset_curr_round,
                                                           mat_subset_curr_round,
                                                           as<IntegerVector>(wrap(parent_curr_round)),
                                                           predvecs_curr_round);
    
    
    lik_subset_curr_round=Rcpp::as<std::vector<double>>(eval_modeltemp[0]);
    table_subset_curr_round=eval_modeltemp[1];
    mat_subset_curr_round=eval_modeltemp[2];
    //overall_count=overall_trees.size();
    parent_curr_round=Rcpp::as<std::vector<int>>(eval_modeltemp[3]);
    
    predvecs_curr_round=eval_modeltemp[4];
    
    
    
    if(table_subset_curr_round.size()==0){
      // Rcout << "Line 6431. j = " << j <<  ".\n";
      
      break;
    }
    
    
    for(int k=0;k<table_subset_curr_round.size();k++){
      overall_trees[overall_count]=table_subset_curr_round[k];
      overall_lik[overall_count]=lik_subset_curr_round[k];
      overall_mat[overall_count]=mat_subset_curr_round[k];
      overall_parent[overall_count]=parent_curr_round[k];
      
      overall_predvecs[overall_count]=predvecs_curr_round[k];
      
      overall_count++;
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
        
        overall_predvecs=resize_bigger(overall_predvecs,overall_size);
        
      }
    }
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    
    overall_predvecs=resize(overall_predvecs,overall_count);
    
    // Rcout << "Line 6448. j = " << j <<  ".\n";
    
    if(less_greedy==1){
      
    }else{
      //eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
      eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                    overall_predvecs);
      overall_lik2=eval_model[0];
      overall_trees=eval_model[1];
      overall_mat=eval_model[2];
      overall_count=overall_trees.size();
      overall_parent2=eval_model[3];
      overall_predvecs=eval_model[4];
      
      //add in check to see if OW accepted more than the top maxOW models...
      if(overall_lik2.size()>maxOWsize){
        IntegerVector owindices=orderforOW(overall_lik2);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        
        List temp_opreds(maxOWsize);
        //Rcout << "Line 1849. j = " << j << ". \n";
        //now only select those elements
        for(int t=0;t<maxOWsize;t++){  
          temp_olik[t]=overall_lik2[owindices[t]];
          temp_otrees[t]=overall_trees[owindices[t]];
          temp_omat[t]= overall_mat[owindices[t]];
          temp_oparent[t]=overall_parent2[owindices[t]];
          
          temp_opreds[t]=overall_predvecs[owindices[t]];
        }
        
        overall_lik2=temp_olik;
        overall_trees=temp_otrees;
        overall_mat=temp_omat;
        overall_count=overall_trees.size();
        overall_parent2=temp_oparent;
        
        overall_predvecs=temp_opreds;
      }
      
      
    }
    
    tree_table=table_subset_curr_round;
    //IntegerVector temp1(table_subset_curr_round.size(),1);
    //err_list=temp1;
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }
    tree_mat=mat_subset_curr_round;
    parent=parent_curr_round;
    
    if(split_rule_node==1){
      NumericVector temp_preds;
      List updated_curr_preds;
      NumericVector new_mean;
      lowest_BIC=min(overall_lik2);
      
      
      //NumericMatrix curr_resids(resids.nrow(),resids.ncol());
      
      List temp(table_subset_curr_round.size());
      
      cp_mat_list=temp;
      
      
      // Rcout << "Line 6536. j = " << j << ". \n";
      //Rcout << "table_subset_curr_round.size() = " << table_subset_curr_round.size() << ". \n";
      //Rcout << "parent.size() = " << parent.size() << ". \n";
      //Rcout << "parent_curr_round.size() = " << parent_curr_round.size() << ". \n";
      //Rcout << "resids.ncol() = " << resids.ncol() << ". \n";
      //Rcout << "parent = " << parent << ". \n";
      
      for(int k=0;k<table_subset_curr_round.size();k++){
        //Rcout << "Line 3260. j = " << j << ". \n";
        
        //Rcout << "parent_curr_round[k] = " << parent_curr_round[k] << ". \n";
        
        
        NumericVector terminal_nodes;
        
        if(parent_curr_round[k]==-1){
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a);
        }else{    
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a);
        }  
        //Rcout << "Line 2533. j = " << j << ". \n";
        terminal_nodes=find_term_nodes(table_subset_curr_round[k]);
        updated_curr_preds=update_predictions(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,D1.n_rows);
        NumericVector test_res;
        
        
        
        
        if(parent_curr_round[k]==-1){
          test_res=resids(_,0); 
        }else{
          test_res=resids(_,parent_curr_round[k]);
        }
        
        NumericVector curr_test_res=updated_curr_preds[1];
        //Rcout << "Line 2548. j = " << j << ". \n";
        
        // if(parent_curr_round[k]==-1){
        //   curr_resids(_,0)=test_res-curr_test_res;
        // }else{
        //   curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;
        // }
        NumericVector temp_curr_resids=test_res-curr_test_res;
        
        
        List cp_mat_list1;
        if(gridpoint==0){
          cp_mat_list1=make_pelt_cpmat(wrap(D1),temp_curr_resids,pen,num_cp);
        }else{
          cp_mat_list1=make_gridpoint_cpmat(wrap(D1),temp_curr_resids,gridsize,num_cp);
        }
        
        cp_mat_list[k]=cp_mat_list1[0];
        
        //Rcout << "Line 2567. j = " << j << ". k = " << k << " . \n";
        
        
        
      }
      
      
      // List temp(0);
      // 
      // cp_mat_list=temp;
      // 
      // for(int f=0;f<curr_resids.ncol();f++){
      //   List cp_mat_list1;
      //   if(gridpoint==0){
      //     cp_mat_list1=make_pelt_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
      //   }else{
      //     cp_mat_list1=make_gridpoint_cpmat(wrap(D1),curr_resids(_,f),gridsize,num_cp);
      //   }
      //   
      //   cp_mat_list.push_back(cp_mat_list1[0]);      
      // }
      
    }//end of if-statement split_rule_node==1
  }
  
  
  // Rcout << "Line 6615.\n";
  // Rcout << "overall_lik = " << wrap(overall_lik) << " .\n";
  // Rcout << "overall_lik2 = " << overall_lik2 << " .\n";
  
  // Rcout << "overall_lik.size()" << overall_lik.size() << " .\n";
  //Rcout << "overall_count" << overall_count << " .\n";
  //Rcout << "overall_trees" << overall_trees.size() << " .\n";
  
  overall_trees=resize(overall_trees,overall_count);
  overall_mat=resize(overall_mat,overall_count); 
  overall_lik.resize(overall_count);
  overall_parent.resize(overall_count);
  
  overall_predvecs=resize(overall_predvecs,overall_count);
  // Rcout << "Line 6626.\n";
  // Rcout << "overall_lik = " << wrap(overall_lik) << " .\n";
  
  if(less_greedy==1){
    // Rcout << "Line 6629.\n";
    
    eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                  overall_predvecs);
    // Rcout << "Line 6633.\n";
    overall_lik2=eval_model[0];
    
    // Rcout << "Line 6634.\n";
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
    overall_predvecs=eval_model[4];
    // Rcout << "Line 6640.\n";
    
  }
  
  //Rcout << "overall_lik" << overall_lik.size() << " .\n";
  //Rcout << "overall_trees" << overall_trees.size() << " .\n";
  //Rcout << "overall_count" << overall_count << " .\n";
  
  // Rcout << "Line 6643.\n";
  
  
  NumericVector temp_preds;
  NumericVector temp_true_preds;
  
  List updated_preds;
  NumericVector new_mean;  
  NumericMatrix overall_test_preds(test_data.nrow(),overall_trees.size());  
  NumericMatrix overallpreds(D1.n_rows,overall_trees.size());
  lowest_BIC=min(overall_lik2);
  
  
  NumericMatrix overallpreds_total(D1.n_rows,overall_trees.size());
  
  // Rcout << "Line 6508. \n";
  for(int k=0;k<overall_trees.size();k++){
    //NumericVector terminal_nodes;
    //Rcout << "Line 1947. k = " << k << ". \n";
    if(overall_parent2[k]==-1){
      //Rcout << "Line 1949. k = " << k << ". \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,0),a);
    }else{   
      //Rcout << "Line 1952. k = " << k << ". \n";
      //Rcout << "Line 1952. overall_parent2[k] = " << overall_parent2[k] << ". \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a);
    }
    //Rcout << "Line 1953. k = " << k << ". \n";
    //terminal_nodes=find_term_nodes(overall_trees[k]);
    updated_preds=update_predictions(overall_trees[k],overall_mat[k],new_mean,D1.n_rows);
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs(test_data,overall_trees[k]//,new_mean
    );
    temp_preds=updated_preds[1];
    overallpreds(_,k)=temp_preds;
    if(is_test_data)overall_test_preds(_,k)=test_preds;
    //Rcout << "Line 1961. k = " << k << ". \n";
    
    temp_true_preds=overall_predvecs[k];
    overallpreds_total(_,k)=temp_true_preds;
    
  }
  // Rcout << "Line 6685.\n";
  
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
  //arma::colvec predicted_values=sum(M1,1);
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  //Rcout << "Line 3393. \n";
  
  
  //arma::colvec predicted_test_values=sum(M2,1);
  // List ret(8);
  // ret[0]=overall_lik2;
  // ret[1]=overall_trees;
  // ret[2]=overall_mat;
  // ret[3]=predicted_values;
  // ret[4]=overall_parent2;
  // ret[5]=wrap(M1);
  // ret[6]=lowest_BIC;
  // ret[7]=wrap(M2);
  List ret(8);
  ret[0]=overall_lik2;
  ret[1]=overall_trees;
  ret[2]=overall_mat;
  ret[3]=overall_parent2;
  ret[4]=wrap(M1);
  ret[5]=lowest_BIC;
  ret[6]=wrap(M2);
  ret[7]=overallpreds_total;
  return(ret);
}
//###################################################################################//

// [[Rcpp::export]]

List get_best_trees_update_splits_exact(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                                        arma::mat& D1,NumericMatrix resids,double a,double mu,double nu,double lambda,double c,
                                        double sigma_mu,List tree_table,List tree_mat,double lowest_BIC,//int first_round,
                                        IntegerVector parent,List cp_mat_list,//IntegerVector err_list,
                                        NumericMatrix test_data,double alpha,double beta,bool is_test_data,double pen,int num_cp,bool split_rule_node,bool gridpoint,int maxOWsize,int num_splits,int gridsize, bool zero_split,
                                        unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split
){
  List eval_model;
  NumericVector lik_list;
  List best_subset;
  int overall_size=1000;
  List overall_trees(overall_size);
  NumericVector overall_lik2;
  IntegerVector overall_parent2;
  List overall_mat(overall_size);
  std::vector<int> overall_parent(overall_size);
  std::vector<double> overall_lik(overall_size);
  
  List overall_predvecs(overall_size);
  
  
  int overall_count=0;  
  // Rcout << "Line 2800";
  if(zero_split==1){
    //Rcout << "Line 2802";
    
    //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
    overall_trees[0]= tree_table[0];
    overall_mat[0]= tree_mat[0];
    overall_parent[0]=-1;
    overall_parent2[0]=-1;
    //double lik_temp=likelihood_function(resids,tree_table[0],tree_mat[0],a,mu,nu,lambda);
    List lik_templist=likelihood_function2_exact(resids,tree_table[0],tree_mat[0],a,mu,nu,lambda);
    double lik_temp= as<double>(lik_templist[0]);
    NumericVector temp_predvec=lik_templist[1];
    
    
    double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
    //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
    double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp));
    
    overall_lik[0]= lowest_BIC_temp;
    //Rcout << "Next zero split tree lowest_BIC_temp = " << lowest_BIC_temp << ".\n"; 
    NumericVector templikvec(1);
    templikvec[0]=lowest_BIC_temp;
    overall_lik2=templikvec;
    
    overall_predvecs[0]=temp_predvec;
    
    overall_count=1;
    
  }
  //Rcout << "Line 1759 .\n";
  NumericVector test_preds;
  //for(int j=0;j<0;j++){
  
  for(int j=0;j<num_splits;j++){
    int lsize=1000;
    List table_subset_curr_round(lsize);
    std::vector<double> lik_subset_curr_round(lsize);
    List mat_subset_curr_round(lsize);
    std::vector<int> parent_curr_round(lsize);
    
    List predvecs_curr_round(lsize);
    
    int count=0;
    for(int i=0;i<tree_table.size();i++){
      //NumericMatrix temp_list=cp_mat_list[0];
      //Rcout << "Line 3136. j = " << j << " i = "<<  i << ".\n";
      parent=-1;
      
      if(split_rule_node==1){
        if(j==0){
          best_subset=get_best_split_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                           lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
                                           min_num_obs_for_split,min_num_obs_after_split//,first_round
          ); 
        }else{
          best_subset=get_best_split_2_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                             lowest_BIC,parent[0],cp_mat_list[i],alpha,beta,maxOWsize,
                                             min_num_obs_for_split,min_num_obs_after_split//,first_round
          ); 
        }
      }else{
        throw std::range_error("get_best_trees_update_splits should only apply when split_rule_node==1");
        
        best_subset=get_best_split_2_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
                                           lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
                                           min_num_obs_for_split,min_num_obs_after_split//,first_round
        ); 
      }
      
      
      
      
      //Rcout << "Line 2392. j = " << j << " i = "<<  i << ".\n";
      //Rcout << "Line 2392. j = " << j << " i = "<<  i << ".\n";
      if(best_subset.size()==1){
        continue;
      }
      // List temp_trees=best_subset[4];
      // List temp_mat=best_subset[6];
      // lik_list=best_subset[5];
      // IntegerVector temp_parent=best_subset[7];
      List temp_trees=best_subset[0];
      List temp_mat=best_subset[2];
      lik_list=best_subset[1];
      IntegerVector temp_parent=best_subset[3];
      
      List temp_predlist=best_subset[4];
      
      if(temp_parent.size()!= temp_trees.size()){
        throw std::range_error("there should be a parent for each tree!!!");
      }
      if(lik_list.size()==0){
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      if(min(lik_list)<lowest_BIC){
        lowest_BIC=min(lik_list);
      }
      //Rcout << "temp_parent" << temp_parent << " .\n";
      for(int k=0;k<temp_trees.size();k++){
        table_subset_curr_round[count]=temp_trees[k];
        lik_subset_curr_round[count]=lik_list[k];
        mat_subset_curr_round[count]=temp_mat[k];
        parent_curr_round[count]=temp_parent[k];
        predvecs_curr_round[count]=temp_predlist[k];
        
        count++;
        //Rcout << "Line 2419. j = " << j << " i = "<<  i << ".\n";
        if(count==(lsize-1)){
          lsize=lsize*2;
          table_subset_curr_round=resize_bigger(table_subset_curr_round,lsize);
          mat_subset_curr_round=resize_bigger(mat_subset_curr_round,lsize);
          lik_subset_curr_round.resize(lsize);
          parent_curr_round.resize(lsize);
          predvecs_curr_round=resize_bigger(predvecs_curr_round,lsize);
          
        }
      }
    }
    table_subset_curr_round=resize(table_subset_curr_round,count);
    mat_subset_curr_round=resize(mat_subset_curr_round,count);
    lik_subset_curr_round.resize(count);
    parent_curr_round.resize(count);
    predvecs_curr_round=resize(predvecs_curr_round,count);
    
    //NumericVector testparvec = wrap(parent_curr_round);
    //Rcout << "parent_curr_round" << testparvec << " .\n";
    
    
    if(table_subset_curr_round.size()==0){
      break;
    }
    
    List eval_modeltemp=evaluate_model_occams_window_exact(as<NumericVector>(wrap(lik_subset_curr_round)),
                                                           lowest_BIC,
                                                           log(c),
                                                           table_subset_curr_round,
                                                           mat_subset_curr_round,
                                                           as<IntegerVector>(wrap(parent_curr_round)),
                                                           predvecs_curr_round);
    
    lik_subset_curr_round=Rcpp::as<std::vector<double>>(eval_modeltemp[0]);
    table_subset_curr_round=eval_modeltemp[1];
    mat_subset_curr_round=eval_modeltemp[2];
    //overall_count=overall_trees.size();
    parent_curr_round=Rcpp::as<std::vector<int>>(eval_modeltemp[3]);
    predvecs_curr_round=eval_modeltemp[4];
    
    
    
    
    
    if(table_subset_curr_round.size()==0){
      break;
    }
    
    
    for(int k=0;k<table_subset_curr_round.size();k++){
      overall_trees[overall_count]=table_subset_curr_round[k];
      overall_lik[overall_count]=lik_subset_curr_round[k];
      overall_mat[overall_count]=mat_subset_curr_round[k];
      overall_parent[overall_count]=parent_curr_round[k];
      
      overall_predvecs[overall_count]=predvecs_curr_round[k];
      
      overall_count++;
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
        
        overall_predvecs=resize_bigger(overall_predvecs,overall_size);
        
      }
    }
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    
    overall_predvecs=resize(overall_predvecs,overall_count);
    
    
    if(less_greedy==1){
      
    }else{
      //eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
      eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                    overall_predvecs);
      overall_lik2=eval_model[0];
      overall_trees=eval_model[1];
      overall_mat=eval_model[2];
      overall_count=overall_trees.size();
      overall_parent2=eval_model[3];
      overall_predvecs=eval_model[4];
      
      //add in check to see if OW accepted more than the top maxOW models...
      if(overall_lik2.size()>maxOWsize){
        IntegerVector owindices=orderforOW(overall_lik2);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        
        List temp_opreds(maxOWsize);
        //Rcout << "Line 1849. j = " << j << ". \n";
        //now only select those elements
        for(int t=0;t<maxOWsize;t++){  
          temp_olik[t]=overall_lik2[owindices[t]];
          temp_otrees[t]=overall_trees[owindices[t]];
          temp_omat[t]= overall_mat[owindices[t]];
          temp_oparent[t]=overall_parent2[owindices[t]];
          
          temp_opreds[t]=overall_predvecs[owindices[t]];
        }
        
        overall_lik2=temp_olik;
        overall_trees=temp_otrees;
        overall_mat=temp_omat;
        overall_count=overall_trees.size();
        overall_parent2=temp_oparent;
        
        overall_predvecs=temp_opreds;
      }
      
      
    }
    
    tree_table=table_subset_curr_round;
    //IntegerVector temp1(table_subset_curr_round.size(),1);
    //err_list=temp1;
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }
    tree_mat=mat_subset_curr_round;
    parent=parent_curr_round;
    
    if(split_rule_node==1){
      NumericVector temp_preds;
      List updated_curr_preds;
      NumericVector new_mean;
      lowest_BIC=min(overall_lik2);
      
      
      //NumericMatrix curr_resids(resids.nrow(),resids.ncol());
      
      List temp(table_subset_curr_round.size());
      
      cp_mat_list=temp;
      
      
      //Rcout << "Line 2519. j = " << j << ". \n";
      //Rcout << "table_subset_curr_round.size() = " << table_subset_curr_round.size() << ". \n";
      //Rcout << "parent.size() = " << parent.size() << ". \n";
      //Rcout << "parent_curr_round.size() = " << parent_curr_round.size() << ". \n";
      //Rcout << "resids.ncol() = " << resids.ncol() << ". \n";
      //Rcout << "parent = " << parent << ". \n";
      
      for(int k=0;k<table_subset_curr_round.size();k++){
        //Rcout << "Line 2521. j = " << j << ". \n";
        
        //Rcout << "parent_curr_round[k] = " << parent_curr_round[k] << ". \n";
        
        
        //NumericVector terminal_nodes;
        
        if(parent_curr_round[k]==-1){
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a);
        }else{    
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a);
        }  
        //Rcout << "Line 2533. j = " << j << ". \n";
        NumericVector terminal_nodes=find_term_nodes(table_subset_curr_round[k]);
        updated_curr_preds=update_predictions(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,D1.n_rows);
        NumericVector test_res;
        
        
        
        
        if(parent_curr_round[k]==-1){
          test_res=resids(_,0); 
        }else{
          test_res=resids(_,parent_curr_round[k]);
        }
        
        //NumericVector curr_test_res=updated_curr_preds[1];
        //Rcout << "Line 2548. j = " << j << ". \n";
        
        // if(parent_curr_round[k]==-1){
        //   curr_resids(_,0)=test_res-curr_test_res;
        // }else{
        //   curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;
        // }
        //NumericVector temp_curr_resids=test_res-curr_test_res;
        
        //Loop over terminal nodes
        List Tempcpmatlist(terminal_nodes.size());
        for(int nodeind=0;nodeind<terminal_nodes.size();nodeind++){
          //Rcout << "Line 3369. j = " << j << ". k = " << k << " . \n";
          
          arma::uvec term_obsarma=find_term_obs(mat_subset_curr_round[k],terminal_nodes[nodeind]);
          //IntegerVector term_obs=wrap(term_obsarma);
          IntegerVector term_obs=wrap(arma::conv_to<arma::ivec>::from(term_obsarma));
          
          //Rcout << "Line 3373. j = " << j << ". k = " << k << " . \n";
          //Rcout << "term_obs = " << term_obs << " . \n";
          //Rcout << "term_obs.size() = " << term_obs.size() << " . \n";
          //Rcout << "temp_curr_resids.size() = " << temp_curr_resids.size() << " . \n";
          
          NumericVector tempsubset = test_res[term_obs];
          
          
          List cp_mat_list1;
          arma::mat tempdata_subset = D1.rows(term_obsarma);
          if(gridpoint==0){
            //cp_mat_list1=make_pelt_cpmat(wrap(tempdata_subset),temp_curr_resids[term_obs],pen,num_cp);
            cp_mat_list1=make_pelt_cpmat(wrap(tempdata_subset),tempsubset,pen,num_cp);
          }else{
            //cp_mat_list1=make_gridpoint_cpmat(wrap(tempdata_subset),temp_curr_resids[term_obs],gridsize,num_cp);
            cp_mat_list1=make_gridpoint_cpmat(wrap(tempdata_subset),tempsubset,gridsize,num_cp);
          }
          //Rcout << "Line 3386. j = " << j << ". k = " << k << " . \n";
          
          Tempcpmatlist[nodeind]=cp_mat_list1[0];
          //Rcout << "Line 3389. j = " << j << ". k = " << k << " . \n";
          
        }
        //Rcout << "Line 3384. j = " << j << ". k = " << k << " . \n";
        
        cp_mat_list[k]=Tempcpmatlist;
        //Rcout << "Line 3395. j = " << j << ". k = " << k << " . \n";
        
        
        //cp_mat_list[k]=cp_mat_list1[0];
        
        //Rcout << "Line 3387. j = " << j << ". k = " << k << " . \n";
        
        
        
      }
      
      
      // List temp(0);
      // 
      // cp_mat_list=temp;
      // 
      // for(int f=0;f<curr_resids.ncol();f++){
      //   List cp_mat_list1;
      //   if(gridpoint==0){
      //     cp_mat_list1=make_pelt_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
      //   }else{
      //     cp_mat_list1=make_gridpoint_cpmat(wrap(D1),curr_resids(_,f),gridsize,num_cp);
      //   }
      //   
      //   cp_mat_list.push_back(cp_mat_list1[0]);      
      // }
      
    }//end of if-statement split_rule_node==1
  }
  overall_trees=resize(overall_trees,overall_count);
  overall_mat=resize(overall_mat,overall_count); 
  overall_lik.resize(overall_count);
  overall_parent.resize(overall_count);
  
  overall_predvecs=resize(overall_predvecs,overall_count);
  
  
  if(less_greedy==1){
    
    eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                  overall_predvecs);    
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
    overall_predvecs=eval_model[4];
    
  }
  
  
  //Rcout << "overall_lik" << overall_lik.size() << " .\n";
  //Rcout << "overall_trees" << overall_trees.size() << " .\n";
  //Rcout << "overall_count" << overall_count << " .\n";
  
  
  NumericVector temp_preds;
  NumericVector temp_true_preds;
  
  List updated_preds;
  NumericVector new_mean;  
  NumericMatrix overall_test_preds(test_data.nrow(),overall_trees.size());  
  NumericMatrix overallpreds(D1.n_rows,overall_trees.size());
  lowest_BIC=min(overall_lik2);
  
  
  NumericMatrix overallpreds_total(D1.n_rows,overall_trees.size());
  
  //Rcout << "Line 2596. \n";
  for(int k=0;k<overall_trees.size();k++){
    //NumericVector terminal_nodes;
    //Rcout << "Line 1947. k = " << k << ". \n";
    if(overall_parent2[k]==-1){
      //Rcout << "Line 1949. k = " << k << ". \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,0),a);
    }else{   
      //Rcout << "Line 1952. k = " << k << ". \n";
      //Rcout << "Line 1952. overall_parent2[k] = " << overall_parent2[k] << ". \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a);
    }
    //Rcout << "Line 1953. k = " << k << ". \n";
    //terminal_nodes=find_term_nodes(overall_trees[k]);
    updated_preds=update_predictions(overall_trees[k],overall_mat[k],new_mean,D1.n_rows);
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs(test_data,overall_trees[k]//,new_mean
    );
    temp_preds=updated_preds[1];
    overallpreds(_,k)=temp_preds;
    if(is_test_data)overall_test_preds(_,k)=test_preds;
    //Rcout << "Line 1961. k = " << k << ". \n";
    
    temp_true_preds=overall_predvecs[k];
    overallpreds_total(_,k)=temp_true_preds;
    
  }
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
  //arma::colvec predicted_values=sum(M1,1);
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  //Rcout << "Line 3130. \n";
  
  
  //arma::colvec predicted_test_values=sum(M2,1);
  // List ret(8);
  // ret[0]=overall_lik2;
  // ret[1]=overall_trees;
  // ret[2]=overall_mat;
  // ret[3]=predicted_values;
  // ret[4]=overall_parent2;
  // ret[5]=wrap(M1);
  // ret[6]=lowest_BIC;
  // ret[7]=wrap(M2);
  List ret(8);
  ret[0]=overall_lik2;
  ret[1]=overall_trees;
  ret[2]=overall_mat;
  ret[3]=overall_parent2;
  ret[4]=wrap(M1);
  ret[5]=lowest_BIC;
  ret[6]=wrap(M2);
  ret[7]=overallpreds_total;
  return(ret);
}
//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_sum_exact(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                              arma::mat& D1,NumericMatrix resids,double a,double mu,double nu,double lambda,
                              double c,double sigma_mu,List tree_table,List tree_mat,double lowest_BIC,//int first_round,
                              IntegerVector parent,List cp_mat_list,IntegerVector err_list,NumericMatrix test_data,double alpha,double beta,bool is_test_data,double pen,int num_cp,bool split_rule_node,bool gridpoint,int maxOWsize,List prev_sum_trees,List prev_sum_trees_mat,NumericVector y_scaled,int num_splits,int gridsize,bool zero_split,
                              unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split
){
  // Rcout << "Get to start of get_best_trees_sum_exact. \n";
  
  List eval_model;
  NumericVector lik_list;
  List best_subset;
  int overall_size=1000;
  List overall_trees(overall_size);
  NumericVector overall_lik2;
  IntegerVector overall_parent2;
  List overall_mat(overall_size);
  
  List overall_predvecs(overall_size);
  
  
  int overall_count=0;  
  std::vector<int> overall_parent(overall_size);
  std::vector<double> overall_lik(overall_size);
  NumericVector test_preds;
  
  // Rcout <<"line 7302.\n";
  
  //////   //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  if(zero_split==1){
    for(int q=0; q<parent.size();q++){
      //Rcout << "q= "<< q << ". \n";
      //Rcout << "length of parent = "<< parent.size() << ". \n";
      //Rcout << "length of prev_sum_trees = "<< prev_sum_trees.size() << ". \n";
      
      // if(parent.size()!= prev_sum_trees.size() )throw std::range_error("length of parent !=length of prev_sum_trees ");
      
      SEXP s_temp = prev_sum_trees[parent[q]];
      //Rcout <<"line 1999.\n";
      
      if(is<List>(s_temp)){
        //Rcout << "s is a list. \n";
        List sum_trees2_temp=prev_sum_trees[parent[q]];
        List sum_trees_mat2_temp=prev_sum_trees_mat[parent[q]];
        sum_trees2_temp.push_back(tree_table[0]);
        sum_trees_mat2_temp.push_back(tree_mat[0]);
        //double lik_temp=sumtree_likelihood_function2_exact(y_scaled,sum_trees2_temp,sum_trees_mat2_temp,y_scaled.size(),a,nu,lambda);  
        
        
        List lik_listtemp=sumtree_likelihood_function2_exact(y_scaled,sum_trees2_temp,sum_trees_mat2_temp,y_scaled.size(),a,nu,lambda); 
        double lik_temp=as<double>(lik_listtemp[0]);
        NumericVector temppredoutput=lik_listtemp[1];
        
        double tree_prior_temp=1;
        //int p_other=0;
        //NumericVector other_int_nodes;
        //Rcout << "Get to loop over t. \n";
        
        for(int t=0;t<sum_trees2_temp.size();t++){
          NumericMatrix tree=sum_trees2_temp[t];
          NumericMatrix mat=sum_trees_mat2_temp[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior_temp*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
        //Rcout << "Finish Loop. \n";
        
        //double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other)*log(D1.n_rows);  
        double BIC=-2*(lik_temp+log(tree_prior_temp));  
        
        ////Rcout << "Get to fill in overall_. \n";
        
        overall_trees[overall_count]= tree_table[0];
        overall_mat[overall_count]= tree_mat[0];
        overall_parent[overall_count]=parent[q];
        //Rcout << "parent[q] = " << parent[q] <<  ".\n";
        //Rcout << "When q=" << q << " parent[q]=" << parent[q] << ". overall_count =" << overall_count << ".\n";
        
        //double lik_temp=likelihood_function(resids[0],tree_table[0],tree_mat[0],a,mu,nu,lambda);
        //double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
        //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
        overall_lik[overall_count]= BIC;
        //Rcout << "Get to end of adding no split trees. \n";
        
        overall_predvecs[overall_count]=temppredoutput;
        
        
        overall_count++;
      }else{
        //Rcout << "s is not a list. \n";
        NumericMatrix sum_trees2_temp=prev_sum_trees[parent[q]];
        NumericMatrix sum_trees_mat2=prev_sum_trees_mat[parent[q]];
        List st(2);
        List st_mat(2);
        st[0]=sum_trees2_temp;
        st[1]=tree_table[0];
        st_mat[0]=sum_trees_mat2;
        st_mat[1]=tree_mat[0];
        // return(st);
        //double lik_temp=sumtree_likelihood_function2_exact(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        
        List lik_listtemp=sumtree_likelihood_function2_exact(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        double lik_temp=as<double>(lik_listtemp[0]);
        NumericVector temppredoutput=lik_listtemp[1];
        
        double tree_prior_temp=1;
        //int p_other=0;
        //NumericVector other_int_nodes;
        for(int t=0;t<st.size();t++){
          NumericMatrix tree=st[t];
          NumericMatrix mat=st_mat[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior_temp*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
        
        //double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other)*log(D1.n_rows);  
        double BIC=-2*(lik_temp+log(tree_prior_temp));  
        
        
        overall_trees[overall_count]= tree_table[0];
        overall_mat[overall_count]= tree_mat[0];
        overall_parent[overall_count]=parent[q];
        //Rcout << "When q=" << q << " parent[q]=" << parent[q] << ". overall_count =" << overall_count << ".\n";
        //double lik_temp=likelihood_function(resids[0],tree_table[0],tree_mat[0],a,mu,nu,lambda);
        //double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
        //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
        overall_lik[overall_count]= BIC;
        
        overall_predvecs[overall_count]=temppredoutput;
        
        overall_count++;
        //Rcout << "Get to end of adding no split trees. \n";
        
      }
      //Rcout <<"line 2073.\n";
      //Rcout << "overall_count = " << overall_count << ".\n";
      //Rcout << "overall_size = " << overall_size << ".\n";
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        //Rcout <<"line 2081.\n";
        
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
        overall_predvecs=resize_bigger(overall_predvecs,overall_size);
        
      }
    }
    //int overall_count=0; 
    //Rcout <<"line 2082.\n";
    //Rcout << "overall_count = " << overall_count << ".\n";
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    overall_predvecs=resize(overall_predvecs,overall_count);
    
    //eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
    eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                  overall_predvecs);
    
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
    overall_predvecs=eval_model[4];
    //Rcout <<"line 2094.\n";
    
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=order_(overall_lik2);
      owindices=owindices-1;
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);
      List temp_otrees(maxOWsize);
      List temp_omat(maxOWsize);
      IntegerVector temp_oparent(maxOWsize);
      List temp_opreds(maxOWsize);
      //Rcout <<"line 2106.\n";
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){  
        temp_olik[t]=overall_lik2[owindices[t]];
        temp_otrees[t]=overall_trees[owindices[t]];
        temp_omat[t]= overall_mat[owindices[t]];
        temp_oparent[t]=overall_parent2[owindices[t]];
        temp_opreds[t]=overall_predvecs[owindices[t]];
      }
      
      overall_lik2=temp_olik;
      overall_trees=temp_otrees;
      overall_mat=temp_omat;
      overall_count=overall_trees.size();
      overall_parent2=temp_oparent;
      overall_predvecs=temp_opreds;
    }
    //Rcout <<"line 2122.\n";
    
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }
  }
  
  // Rcout <<"Get past zero split tree block of code.\n";
  // Rcout <<"line 7498.\n";
  
  //Rcout <<"num_splits= " << num_splits << " .\n";
  
  for(int j=0;j<num_splits;j++){
    int lsize=1000;
    List table_subset_curr_round(lsize);
    std::vector<double> lik_subset_curr_round(lsize);
    List mat_subset_curr_round(lsize);
    std::vector<int> parent_curr_round(lsize);
    
    List predvecs_curr_round(lsize);
    
    // Rcout <<"line 7510.\n";
    
    
    int count=0;
    //// Rcout <<"LENGTH OF TREE TABLE LIST = " << tree_table.size() << " !!!!!!!!!.\n";
    //Rcout <<"LENGTH OF cp_mat_list LIST = " << cp_mat_list.size() << " !!!!!!!!!.\n";
    for(int i=0;i<tree_table.size();i++){
      // if(first_round==1){
      //   parent=-1;
      //   //NumericMatrix temp_list=cp_mat_list[0];
      //   best_subset=get_best_split_exact(resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
      //                              lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
      //                              min_num_obs_for_split,min_num_obs_after_split//,first_round
      //                                ); 
      //   
      // }else{
      if(err_list[i]==0){
        //NumericMatrix test_tree=tree_table[i];
        //NumericMatrix test_treemat=tree_mat[i];
        //NumericMatrix test_cpmat= cp_mat_list[parent[i]];
        //need to append current tree_table[i] to its parent sum_of_trees   
        
        // Rcout << " i = " << i << ".\n";
        // Rcout << " j = " << j << ".\n";
        // Rcout <<"line 7535.\n";
        
        if(split_rule_node==1){
          if(j==0){
            best_subset=get_best_split_sum_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                                 parent[i],cp_mat_list[i],
                                                                      alpha,beta,maxOWsize,//first_round,
                                                                      prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                      min_num_obs_for_split,min_num_obs_after_split);   
          }else{
            best_subset=get_best_split_sum_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                                 parent[i],cp_mat_list[i],
                                                                      alpha,beta,maxOWsize,//first_round,
                                                                      prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                      min_num_obs_for_split,min_num_obs_after_split);
          }
        }else{
          // Rcout << " line 7552. no update .\n";
          
          //Rcout << " i = " << i << ".\n";
          // Rcout << " parent[i] = " << parent[i] << ".\n";
          //NumericMatrix temptesttable = tree_table[i];
          //Rcout << " temptesttable.nrow() = " << temptesttable.nrow() << ".\n";
          //Rcout << " temptesttable.ncol() = " << temptesttable.ncol() << ".\n";
          //NumericMatrix temptestmat = tree_mat[i];
          
          //Rcout << " temptestmat.nrow() = " << temptestmat.nrow() << ".\n";
          //Rcout << " temptestmat.ncol() = " << temptestmat.ncol() << ".\n";
          
          //NumericMatrix tempcpmat = cp_mat_list[parent[i]];
          //Rcout << " tempcpmat.nrow() = " << tempcpmat.nrow() << ".\n";
          //Rcout << " tempcpmat.ncol() = " << tempcpmat.ncol() << ".\n";
          
          best_subset=get_best_split_sum_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,
                                               D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                               parent[i],cp_mat_list[parent[i]],
                                                                    alpha,beta,maxOWsize,//first_round,
                                                                    prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                    min_num_obs_for_split,min_num_obs_after_split);
        }
        
        // Rcout << " line 7576. no update .\n";
        
        
        // return(best_subset);
      }else if(err_list[i]==1){
        ////Rcout << "CONTINUE. error list.  \n";
        continue;
      }else{
        //Rcout << " i = " << i << ".\n";
        //Rcout << " j = " << j << ".\n";
        
        throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");
        
        List ret_list(6);
        ret_list[0]=9999;
        ret_list[1]=err_list[i];
        ret_list[2]=i;
        ret_list[3]=j;
        ret_list[4]=tree_table;
        ret_list[5]=err_list;
        return(ret_list);
        throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");
      }
      //}
      // Rcout << "Get past get_best_split_exact. \n";
      
      if(best_subset.size()==1){
        //Rcout << "CONTINUE. \n";
        continue;
      }
      // List temp_trees=best_subset[4];
      // List temp_mat=best_subset[6];
      // lik_list=best_subset[5];
      // IntegerVector temp_parent=best_subset[7];
      // Rcout << " line 7610. no update .\n";
      
      List temp_trees=best_subset[0];
      List temp_mat=best_subset[2];
      lik_list=best_subset[1];
      IntegerVector temp_parent=best_subset[3];
      List temp_predlist=best_subset[4];
      //Rcout << " line 7487. no update .\n";
      
      if(temp_parent.size()!= temp_trees.size()){
        throw std::range_error("there should be a parent for each tree!!!");
      }
      if(lik_list.size()==0){
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      
      if(min(lik_list)<lowest_BIC){
        lowest_BIC=min(lik_list);
      }
      for(int k=0;k<temp_trees.size();k++){
        table_subset_curr_round[count]=temp_trees[k];
        lik_subset_curr_round[count]=lik_list[k];
        mat_subset_curr_round[count]=temp_mat[k];
        parent_curr_round[count]=temp_parent[k];
        
        predvecs_curr_round[count]=temp_predlist[k];
        
        count++;
        
        if(count==(lsize-1)){
          lsize=lsize*2;
          table_subset_curr_round=resize_bigger(table_subset_curr_round,lsize);
          mat_subset_curr_round=resize_bigger(mat_subset_curr_round,lsize);
          lik_subset_curr_round.resize(lsize);
          parent_curr_round.resize(lsize);
          predvecs_curr_round=resize_bigger(predvecs_curr_round,lsize);
          
        }
      }
      //Rcout << " line 7519. no update .\n";
      
    }
    // Rcout << " line 7652. no update .\n";
    
    table_subset_curr_round=resize(table_subset_curr_round,count);
    mat_subset_curr_round=resize(mat_subset_curr_round,count);
    lik_subset_curr_round.resize(count);
    parent_curr_round.resize(count);
    predvecs_curr_round=resize(predvecs_curr_round,count);
    
    
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0
      break;															// break out of the for-loop,
    }
    
    List eval_modeltemp=evaluate_model_occams_window_exact(as<NumericVector>(wrap(lik_subset_curr_round)),
                                                           lowest_BIC,
                                                           log(c),
                                                           table_subset_curr_round,
                                                           mat_subset_curr_round,
                                                           as<IntegerVector>(wrap(parent_curr_round)),
                                                           predvecs_curr_round);
    
    
    lik_subset_curr_round=Rcpp::as<std::vector<double>>(eval_modeltemp[0]);
    table_subset_curr_round=eval_modeltemp[1];
    mat_subset_curr_round=eval_modeltemp[2];
    //overall_count=overall_trees.size();
    parent_curr_round=Rcpp::as<std::vector<int>>(eval_modeltemp[3]);
    
    predvecs_curr_round=eval_modeltemp[4];
    
    
    if(table_subset_curr_round.size()==0){
      //Rcout << "BREAK. \n";
      break;
    }
    //Rcout << "inner loop round WITHOUT CONTINUE OR BREAK. \n";
    
    for(int k=0;k<table_subset_curr_round.size();k++){
      overall_trees[overall_count]=table_subset_curr_round[k];
      overall_lik[overall_count]=lik_subset_curr_round[k];
      overall_mat[overall_count]=mat_subset_curr_round[k];
      overall_parent[overall_count]=parent_curr_round[k];
      
      overall_predvecs[overall_count]=predvecs_curr_round[k];
      
      
      overall_count++;
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
        
        overall_predvecs=resize_bigger(overall_predvecs,overall_size);
        
      }
    }
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    
    overall_predvecs=resize(overall_predvecs,overall_count);
    
    
    
    if(less_greedy==1){
      
    }else{
      //Rcout << "overall_parent[0] BEFORE OW EVALUATION = " << overall_parent[0] << ".\n";
      //eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
      //
      eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                    overall_predvecs);
      
      
      overall_lik2=eval_model[0];
      overall_trees=eval_model[1];
      overall_mat=eval_model[2];
      overall_count=overall_trees.size();
      overall_parent2=eval_model[3];
      overall_predvecs=eval_model[4];
      
      //Rcout << "overall_parent2[0] AFTER OW EVALUATION" << overall_parent2[0] << ".\n";
      //Rcout << "overall_parent2[0] AFTER OW EVALUATION" << overall_parent2[0] << ".\n";
      
      //add in check to see if OW accepted more than the top maxOW models...
      if(overall_lik2.size()>maxOWsize){
        //Rcout << "MORE THAN MAXOWSIZE!!!!!!!!!!!" << overall_parent2[0] << ".\n";
        
        //find the maxOWsize best models and continue with those!
        IntegerVector owindices=orderforOW(overall_lik2);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        List temp_opreds(maxOWsize);
        
        //now only select those elements
        for(int t=0;t<maxOWsize;t++){  
          temp_olik[t]=overall_lik2[owindices[t]];
          temp_otrees[t]=overall_trees[owindices[t]];
          temp_omat[t]= overall_mat[owindices[t]];
          temp_oparent[t]=overall_parent2[owindices[t]];
          temp_opreds[t]=overall_predvecs[owindices[t]];
        }
        
        overall_lik2=temp_olik;
        overall_trees=temp_otrees;
        overall_mat=temp_omat;
        overall_count=overall_trees.size();
        overall_parent2=temp_oparent;
        overall_predvecs=temp_opreds;
        
        //Rcout << "overall_parent2[0] AFTER REMOVING EXTRA MODELS AND REARRAGING = " << overall_parent2[0] << ".\n";
        
      }
      
      
    }
    
    // Rcout << "Line 7754.\n";
    
    tree_table=table_subset_curr_round;
    IntegerVector temp1(table_subset_curr_round.size(),0);
    err_list=temp1;
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }
    tree_mat=mat_subset_curr_round;
    parent=parent_curr_round;
    
    if(split_rule_node==1){
      NumericVector temp_preds;
      List updated_curr_preds;
      NumericVector new_mean;
      lowest_BIC=min(overall_lik2);
      
      
      //NumericMatrix curr_resids(resids.nrow(),resids.ncol());
      
      
      List temp(table_subset_curr_round.size());
      
      cp_mat_list=temp;
      
      
      
      for(int k=0;k<table_subset_curr_round.size();k++){
        //NumericVector terminal_nodes;
        
        if(parent_curr_round[k]==-1){
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a);
        }else{    
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a);
        }  
        
        //terminal_nodes=find_term_nodes(table_subset_curr_round[k]);
        updated_curr_preds=update_predictions(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,D1.n_rows);
        NumericVector test_res;
        
        if(parent_curr_round[k]==-1){
          test_res=resids(_,0); 
        }else{
          test_res=resids(_,parent_curr_round[k]);
        }
        
        NumericVector curr_test_res=updated_curr_preds[1];
        //Rcout << "Line 2548. j = " << j << ". \n";
        
        // if(parent_curr_round[k]==-1){
        //   curr_resids(_,0)=test_res-curr_test_res;
        // }else{
        //   curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;
        // }
        NumericVector temp_curr_resids=test_res-curr_test_res;
        
        List cp_mat_list1;
        if(gridpoint==0){
          cp_mat_list1=make_pelt_cpmat(wrap(D1),temp_curr_resids,pen,num_cp);
        }else{
          cp_mat_list1=make_gridpoint_cpmat(wrap(D1),temp_curr_resids,gridsize,num_cp);
        }
        
        cp_mat_list[k]=cp_mat_list1[0];
        
      }
      
      
      
      //List temp(0);
      
      //cp_mat_list=temp;
      
      // for(int f=0;f<curr_resids.ncol();f++){
      //   List cp_mat_list1;
      //   if(gridpoint==0){
      //     cp_mat_list1=make_pelt_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
      //   }else{
      //     cp_mat_list1=make_gridpoint_cpmat(wrap(D1),curr_resids(_,f),gridsize,num_cp);
      //   }
      //   
      //   cp_mat_list.push_back(cp_mat_list1[0]);      
      // }
      
    }   //end of if-statement split_rule_node==1
  }
  //Rcout << "Get to end of outer loop. \n";
  //Rcout << "overall_trees.size()= " << overall_trees.size() << " \n";
  //Rcout << "overall_mat.size()= " << overall_mat.size() << " \n";
  //Rcout << "overall_lik.size()= " << overall_lik.size() << " \n";
  //Rcout << "overall_parent.size()= " << overall_parent.size() << " \n";
  //Rcout << "overall_parent2.size()= " << overall_parent2.size() << " \n";
  
  overall_trees=resize(overall_trees,overall_count);
  overall_mat=resize(overall_mat,overall_count); 
  overall_lik.resize(overall_count);
  overall_parent.resize(overall_count);
  overall_predvecs=resize(overall_predvecs,overall_count);
  // Rcout <<"line 7865.\n";
  
  
  if(less_greedy==1){
    
    eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                  overall_predvecs);
    
    
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
    overall_predvecs=eval_model[4];
  }
  
  // Rcout << "Line 7882";
  
  NumericVector temp_preds;
  NumericVector temp_true_preds;
  
  List updated_preds;
  NumericVector new_mean;  
  NumericMatrix overall_test_preds(test_data.nrow(),overall_trees.size());  
  NumericMatrix overallpreds(D1.n_rows,overall_trees.size());
  lowest_BIC=min(overall_lik2);
  
  NumericMatrix overallpreds_total(D1.n_rows,overall_trees.size());
  
  
  //Rcout << "Get to start of update mean loop. \n";
  //Rcout << "overall_trees.size()= " << overall_trees.size() << " \n";
  //Rcout << "overall_mat.size()= " << overall_mat.size() << " \n";
  //Rcout << "overall_lik.size()= " << overall_lik.size() << " \n";
  //Rcout << "overall_parent.size()= " << overall_parent.size() << " \n";
  //Rcout << "overall_parent2.size()= " << overall_parent2.size() << " \n";
  
  for(int k=0;k<overall_trees.size();k++){
    //NumericVector terminal_nodes;
    
    if(overall_parent2[k]==-1){
      //Rcout << "within update mean var. overall_parent2[k]= " << overall_parent2[k] << " \n";
      //Rcout << "k= " << k << " \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,0),a);
    }else{
      //Rcout << "within update mean var. overall_parent2[k]= " << overall_parent2[k] << " \n";
      //Rcout << "k= " << k << " \n";
      
      //Rcout << "a= " << a << " \n";
      SEXP tree_checktype = overall_trees[k];
      if(is<NumericMatrix>(tree_checktype)){
        //Rcout << "tree_checktype is a NumericMatrix. \n";
      }else{
        //Rcout << "tree_checktype is NOT a NumericMatrix. \n";
      }
      SEXP mat_checktype = overall_mat[k];
      if(is<NumericMatrix>(mat_checktype)){
        //Rcout << "mat_checktype is a NumericMatrix. \n";
      }else{
        //Rcout << "mat_checktype is NOT a NumericMatrix. \n";
      }
      //  SEXP vec_checktype = resids(_,overall_parent2[k]);
      //  if(is<NumericVector>(vec_checktype)){
      //Rcout << "vec_checktype is a NumericMatrix. \n";
      //  }else{
      //Rcout << "vec_checktype is NOT a NumericMatrix. \n";
      //  }
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a);
    }
    //Rcout << "Get past update mean_var. k= " << k << " \n";
    
    //terminal_nodes=find_term_nodes(overall_trees[k]);
    updated_preds=update_predictions(overall_trees[k],overall_mat[k],new_mean,D1.n_rows);
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs(test_data,overall_trees[k]//,new_mean
    );
    temp_preds=updated_preds[1];
    overallpreds(_,k)=temp_preds;
    if(is_test_data)overall_test_preds(_,k)=test_preds;
    
    
    temp_true_preds=overall_predvecs[k];
    overallpreds_total(_,k)=temp_true_preds;    
    
  }
  //Rcout << "Get to end of update mean loop. \n";
  
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
  //arma::colvec predicted_values=sum(M1,1);
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  //arma::colvec predicted_test_values=sum(M2,1);
  //Rcout << "Get to end of get_est_trees_sum. \n";
  // List ret(8);
  // ret[0]=overall_lik2;
  // ret[1]=overall_trees;
  // ret[2]=overall_mat;
  // ret[3]=predicted_values;
  // ret[4]=overall_parent2;
  // ret[5]=wrap(M1);
  // ret[6]=lowest_BIC;
  // ret[7]=wrap(M2);
  List ret(8);
  ret[0]=overall_lik2;
  ret[1]=overall_trees;
  ret[2]=overall_mat;
  ret[3]=overall_parent2;
  ret[4]=wrap(M1);
  ret[5]=lowest_BIC;
  ret[6]=wrap(M2);
  ret[7]=overallpreds_total;
  return(ret);
}
//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_sum_update_splits_exact(double less_greedy, double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                                            arma::mat& D1,NumericMatrix resids,double a,double mu,double nu,double lambda,
                                            double c,double sigma_mu,List tree_table,List tree_mat,double lowest_BIC,//int first_round,
                                            IntegerVector parent,List cp_mat_list,IntegerVector err_list,NumericMatrix test_data,double alpha,double beta,bool is_test_data,double pen,int num_cp,bool split_rule_node,bool gridpoint,int maxOWsize,List prev_sum_trees,List prev_sum_trees_mat,NumericVector y_scaled,int num_splits,int gridsize,bool zero_split,
                                            unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split
){
  // Rcout << "Get to start of get_best_trees_sum_update_splits_exact. \n";
  
  List eval_model;
  NumericVector lik_list;
  List best_subset;
  int overall_size=1000;
  List overall_trees(overall_size);
  NumericVector overall_lik2;
  IntegerVector overall_parent2;
  List overall_mat(overall_size);
  
  List overall_predvecs(overall_size);
  
  
  int overall_count=0;  
  std::vector<int> overall_parent(overall_size);
  std::vector<double> overall_lik(overall_size);
  NumericVector test_preds;
  
  // Rcout << "Get to line 8006. \n";
  
  //////   //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  //////   //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  if(zero_split==1){
    for(int q=0; q<parent.size();q++){
      //Rcout << "q= "<< q << ". \n";
      //// Rcout << "length of parent = "<< parent.size() << ". \n";
      //// Rcout << "length of prev_sum_trees = "<< prev_sum_trees.size() << ". \n";
      
      // if(parent.size()!= prev_sum_trees.size() )throw std::range_error("length of parent !=length of prev_sum_trees ");
      
      SEXP s_temp = prev_sum_trees[parent[q]];
      
      // Rcout <<"line 8020.\n";
      
      if(is<List>(s_temp)){
        //Rcout << "s is a list. \n";
        List sum_trees2_temp=prev_sum_trees[parent[q]];
        List sum_trees_mat2_temp=prev_sum_trees_mat[parent[q]];
        sum_trees2_temp.push_back(tree_table[0]);
        sum_trees_mat2_temp.push_back(tree_mat[0]);
        //double lik_temp=sumtree_likelihood_function2_exact(y_scaled,sum_trees2_temp,sum_trees_mat2_temp,y_scaled.size(),a,nu,lambda);  
        
        
        List lik_listtemp=sumtree_likelihood_function2_exact(y_scaled,sum_trees2_temp,sum_trees_mat2_temp,y_scaled.size(),a,nu,lambda); 
        double lik_temp=as<double>(lik_listtemp[0]);
        NumericVector temppredoutput=lik_listtemp[1];
        
        double tree_prior_temp=1;
        //int p_other=0;
        //NumericVector other_int_nodes;
        //Rcout << "Get to loop over t. \n";
        
        for(int t=0;t<sum_trees2_temp.size();t++){
          NumericMatrix tree=sum_trees2_temp[t];
          NumericMatrix mat=sum_trees_mat2_temp[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior_temp*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
        //Rcout << "Finish Loop. \n";
        
        //double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other)*log(D1.n_rows);  
        double BIC=-2*(lik_temp+log(tree_prior_temp));  
        
        ////Rcout << "Get to fill in overall_. \n";
        
        overall_trees[overall_count]= tree_table[0];
        overall_mat[overall_count]= tree_mat[0];
        overall_parent[overall_count]=parent[q];
        //Rcout << "parent[q] = " << parent[q] <<  ".\n";
        //Rcout << "When q=" << q << " parent[q]=" << parent[q] << ". overall_count =" << overall_count << ".\n";
        
        //double lik_temp=likelihood_function(resids[0],tree_table[0],tree_mat[0],a,mu,nu,lambda);
        //double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
        //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
        overall_lik[overall_count]= BIC;
        //Rcout << "Get to end of adding no split trees. \n";
        
        overall_predvecs[overall_count]=temppredoutput;
        
        
        overall_count++;
      }else{
        //Rcout << "s is not a list. \n";
        NumericMatrix sum_trees2_temp=prev_sum_trees[parent[q]];
        NumericMatrix sum_trees_mat2=prev_sum_trees_mat[parent[q]];
        List st(2);
        List st_mat(2);
        st[0]=sum_trees2_temp;
        st[1]=tree_table[0];
        st_mat[0]=sum_trees_mat2;
        st_mat[1]=tree_mat[0];
        // return(st);
        //double lik_temp=sumtree_likelihood_function2_exact(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        
        List lik_listtemp=sumtree_likelihood_function2_exact(y_scaled,st,st_mat,y_scaled.size(),a,nu,lambda);
        double lik_temp=as<double>(lik_listtemp[0]);
        NumericVector temppredoutput=lik_listtemp[1];
        
        double tree_prior_temp=1;
        //int p_other=0;
        //NumericVector other_int_nodes;
        for(int t=0;t<st.size();t++){
          NumericMatrix tree=st[t];
          NumericMatrix mat=st_mat[t];
          //other_int_nodes = find_term_nodes(tree);
          //p_other+=other_int_nodes.size();
          tree_prior_temp*=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree,mat,alpha,beta);
        }
        
        //double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other)*log(D1.n_rows);  
        double BIC=-2*(lik_temp+log(tree_prior_temp));  
        
        
        overall_trees[overall_count]= tree_table[0];
        overall_mat[overall_count]= tree_mat[0];
        overall_parent[overall_count]=parent[q];
        //Rcout << "When q=" << q << " parent[q]=" << parent[q] << ". overall_count =" << overall_count << ".\n";
        //double lik_temp=likelihood_function(resids[0],tree_table[0],tree_mat[0],a,mu,nu,lambda);
        //double tree_prior_temp=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,tree_table[0],tree_mat[0],alpha,beta);
        //double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(D1.n_rows);
        overall_lik[overall_count]= BIC;
        
        overall_predvecs[overall_count]=temppredoutput;
        
        overall_count++;
        //Rcout << "Get to end of adding no split trees. \n";
        
      }
      //Rcout <<"line 2073.\n";
      //Rcout << "overall_count = " << overall_count << ".\n";
      //Rcout << "overall_size = " << overall_size << ".\n";
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        //Rcout <<"line 2081.\n";
        
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
        overall_predvecs=resize_bigger(overall_predvecs,overall_size);
        
      }
    }
    //int overall_count=0; 
    //Rcout <<"line 2082.\n";
    //Rcout << "overall_count = " << overall_count << ".\n";
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    overall_predvecs=resize(overall_predvecs,overall_count);
    
    //eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
    eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                  overall_predvecs);
    
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
    overall_predvecs=eval_model[4];
    //Rcout <<"line 2094.\n";
    
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=order_(overall_lik2);
      owindices=owindices-1;
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);
      List temp_otrees(maxOWsize);
      List temp_omat(maxOWsize);
      IntegerVector temp_oparent(maxOWsize);
      List temp_opreds(maxOWsize);
      //Rcout <<"line 2106.\n";
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){  
        temp_olik[t]=overall_lik2[owindices[t]];
        temp_otrees[t]=overall_trees[owindices[t]];
        temp_omat[t]= overall_mat[owindices[t]];
        temp_oparent[t]=overall_parent2[owindices[t]];
        temp_opreds[t]=overall_predvecs[owindices[t]];
      }
      
      overall_lik2=temp_olik;
      overall_trees=temp_otrees;
      overall_mat=temp_omat;
      overall_count=overall_trees.size();
      overall_parent2=temp_oparent;
      overall_predvecs=temp_opreds;
    }
    //Rcout <<"line 2122.\n";
    
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }
  }
  
  // Rcout <<"get sum trees update splits. Get past zero split tree block of code.\n";
  // Rcout <<"line 8204.\n";
  
  
  for(int j=0;j<num_splits;j++){
    int lsize=1000;
    List table_subset_curr_round(lsize);
    std::vector<double> lik_subset_curr_round(lsize);
    List mat_subset_curr_round(lsize);
    std::vector<int> parent_curr_round(lsize);
    
    
    List predvecs_curr_round(lsize);
    
    
    int count=0;
    // Rcout << "begin loop over splits. j == " << j << ".\n";
    
    //Rcout <<"LENGTH OF TREE TABLE LIST = " << tree_table.size() << " !!!!!!!!!.\n";
    for(int i=0;i<tree_table.size();i++){
      // Rcout << "begin loop over models. j == " << j << ". i == "<< i << ".\n";
      // Rcout <<"line 8224.\n";
      
      // if(first_round==1){
      //   parent=-1;
      //   //NumericMatrix temp_list=cp_mat_list[0];
      //   best_subset=get_best_split_exact(resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),
      //                              lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,
      //                              min_num_obs_for_split,min_num_obs_after_split//,first_round
      //                                ); 
      //   
      // }else{
      if(err_list[i]==0){
        //NumericMatrix test_tree=tree_table[i];
        //NumericMatrix test_treemat=tree_mat[i];
        //NumericMatrix test_cpmat= cp_mat_list[parent[i]];
        //need to append current tree_table[i] to its parent sum_of_trees   
        
        // Rcout <<"line 8241.\n";
        
        if(split_rule_node==1){
          if(j==0){
            best_subset=get_best_split_sum_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                                 parent[i],cp_mat_list[parent[i]],
                                                                      alpha,beta,maxOWsize,//first_round,
                                                                      prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                      min_num_obs_for_split,min_num_obs_after_split);   
          }else{
            // Rcout << " line 8251. with update .\n";
            
            // Rcout << " i = " << i << ".\n";
            // Rcout << " parent[i] = " << parent[i] << ".\n";
            best_subset=get_best_split_sum_2_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                                   parent[i],cp_mat_list[i],
                                                                        alpha,beta,maxOWsize,//first_round,
                                                                        prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                        min_num_obs_for_split,min_num_obs_after_split);
          }
        }else{
          throw std::range_error("get_best_trees_update_splits should only apply when split_rule_node==1");
          best_subset=get_best_split_sum_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,
                                               parent[i],cp_mat_list[i],
                                                                    alpha,beta,maxOWsize,//first_round,
                                                                    prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i,
                                                                    min_num_obs_for_split,min_num_obs_after_split);
        }
        
        
        
        // return(best_subset);
      }else if(err_list[i]==1){
        //Rcout << "j == << " << j << "  \n";
        
        throw std::range_error("err_list[i] is 0");
        
        //Rcout << "CONTINUE. error list.  \n";
        continue;
      }else{
        List ret_list(6);
        ret_list[0]=9999;
        ret_list[1]=err_list[i];
        ret_list[2]=i;
        ret_list[3]=j;
        ret_list[4]=tree_table;
        ret_list[5]=err_list;
        //return(ret_list);
        throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");
      }
      //}
      // Rcout << "Get past get_best_split_exact. j == " << j << ".\n";
      // Rcout << " line 8293. with update .\n";
      
      if(best_subset.size()==1){
        //Rcout << "CONTINUE. \n";
        continue;
      }
      // List temp_trees=best_subset[4];
      // List temp_mat=best_subset[6];
      // lik_list=best_subset[5];
      // IntegerVector temp_parent=best_subset[7];
      List temp_trees=best_subset[0];
      List temp_mat=best_subset[2];
      lik_list=best_subset[1];
      IntegerVector temp_parent=best_subset[3];
      List temp_predlist=best_subset[4];
      if(temp_parent.size()!= temp_trees.size()){
        throw std::range_error("there should be a parent for each tree!!!");
      }
      if(lik_list.size()==0){
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      
      if(min(lik_list)<lowest_BIC){
        lowest_BIC=min(lik_list);
      }
      for(int k=0;k<temp_trees.size();k++){
        table_subset_curr_round[count]=temp_trees[k];
        lik_subset_curr_round[count]=lik_list[k];
        mat_subset_curr_round[count]=temp_mat[k];
        parent_curr_round[count]=temp_parent[k];
        
        predvecs_curr_round[count]=temp_predlist[k];
        
        count++;
        
        if(count==(lsize-1)){
          lsize=lsize*2;
          table_subset_curr_round=resize_bigger(table_subset_curr_round,lsize);
          mat_subset_curr_round=resize_bigger(mat_subset_curr_round,lsize);
          lik_subset_curr_round.resize(lsize);
          parent_curr_round.resize(lsize);
          predvecs_curr_round=resize_bigger(predvecs_curr_round,lsize);
          
        }
      }
    }
    table_subset_curr_round=resize(table_subset_curr_round,count);
    mat_subset_curr_round=resize(mat_subset_curr_round,count);
    lik_subset_curr_round.resize(count);
    parent_curr_round.resize(count);
    predvecs_curr_round=resize(predvecs_curr_round,count);
    
    
    if(table_subset_curr_round.size()==0){
      //Rcout << "BREAK. \n";
      break;
    }
    
    List eval_modeltemp=evaluate_model_occams_window_exact(as<NumericVector>(wrap(lik_subset_curr_round)),
                                                           lowest_BIC,
                                                           log(c),
                                                           table_subset_curr_round,
                                                           mat_subset_curr_round,
                                                           as<IntegerVector>(wrap(parent_curr_round)),
                                                           predvecs_curr_round);
    
    lik_subset_curr_round=Rcpp::as<std::vector<double>>(eval_modeltemp[0]);
    table_subset_curr_round=eval_modeltemp[1];
    mat_subset_curr_round=eval_modeltemp[2];
    //overall_count=overall_trees.size();
    parent_curr_round=Rcpp::as<std::vector<int>>(eval_modeltemp[3]);
    predvecs_curr_round=eval_modeltemp[4];
    
    if(table_subset_curr_round.size()==0){
      //Rcout << "BREAK. \n";
      break;
    }
    //Rcout << "inner loop round WITHOUT CONTINUE OR BREAK. \n";
    
    for(int k=0;k<table_subset_curr_round.size();k++){
      overall_trees[overall_count]=table_subset_curr_round[k];
      overall_lik[overall_count]=lik_subset_curr_round[k];
      overall_mat[overall_count]=mat_subset_curr_round[k];
      overall_parent[overall_count]=parent_curr_round[k];
      
      overall_predvecs[overall_count]=predvecs_curr_round[k];
      
      
      overall_count++;
      
      if(overall_count==(overall_size-1)){
        overall_size=overall_size*2;
        overall_trees=resize_bigger(overall_trees,overall_size);
        overall_lik.resize(overall_size);
        overall_mat=resize_bigger(overall_mat,overall_size);
        overall_parent.resize(overall_size);
        
        overall_predvecs=resize_bigger(overall_predvecs,overall_size);
        
      }
    }
    overall_trees=resize(overall_trees,overall_count);
    overall_lik.resize(overall_count);
    overall_mat=resize(overall_mat,overall_count);
    overall_parent.resize(overall_count);
    
    overall_predvecs=resize(overall_predvecs,overall_count);
    
    // Rcout << " line 8381. with update .\n";
    
    if(less_greedy==1){
      
    }else{
      //Rcout << "overall_parent[0] BEFORE OW EVALUATION = " << overall_parent[0] << ".\n";
      //eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));
      //
      eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                    overall_predvecs);
      
      
      overall_lik2=eval_model[0];
      overall_trees=eval_model[1];
      overall_mat=eval_model[2];
      overall_count=overall_trees.size();
      overall_parent2=eval_model[3];
      overall_predvecs=eval_model[4];
      
      //Rcout << "overall_parent2[0] AFTER OW EVALUATION" << overall_parent2[0] << ".\n";
      //Rcout << "overall_parent2[0] AFTER OW EVALUATION" << overall_parent2[0] << ".\n";
      
      //add in check to see if OW accepted more than the top maxOW models...
      if(overall_lik2.size()>maxOWsize){
        //Rcout << "MORE THAN MAXOWSIZE!!!!!!!!!!!" << overall_parent2[0] << ".\n";
        
        //find the maxOWsize best models and continue with those!
        IntegerVector owindices=orderforOW(overall_lik2);
        owindices=owindices-1;
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);
        List temp_otrees(maxOWsize);
        List temp_omat(maxOWsize);
        IntegerVector temp_oparent(maxOWsize);
        List temp_opreds(maxOWsize);
        
        //now only select those elements
        for(int t=0;t<maxOWsize;t++){  
          temp_olik[t]=overall_lik2[owindices[t]];
          temp_otrees[t]=overall_trees[owindices[t]];
          temp_omat[t]= overall_mat[owindices[t]];
          temp_oparent[t]=overall_parent2[owindices[t]];
          temp_opreds[t]=overall_predvecs[owindices[t]];
        }
        
        overall_lik2=temp_olik;
        overall_trees=temp_otrees;
        overall_mat=temp_omat;
        overall_count=overall_trees.size();
        overall_parent2=temp_oparent;
        overall_predvecs=temp_opreds;
        
        //Rcout << "overall_parent2[0] AFTER REMOVING EXTRA MODELS AND REARRAGING = " << overall_parent2[0] << ".\n";
        
      }
      
    }
    
    
    tree_table=table_subset_curr_round;
    IntegerVector temp1(table_subset_curr_round.size(),0);
    err_list=temp1;
    if(overall_trees.size()<overall_size-1){
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size);
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }else{
      overall_size=2*overall_size;
      overall_trees=resize_bigger(overall_trees,overall_size);
      overall_mat=resize_bigger(overall_mat,overall_size); 
      overall_lik.resize(overall_size);
      overall_parent.resize(overall_size);
      overall_predvecs=resize_bigger(overall_predvecs,overall_size);
      
    }
    tree_mat=mat_subset_curr_round;
    parent=parent_curr_round;
    // Rcout << " line 8460. with update .\n";
    
    if(split_rule_node==1){
      NumericVector temp_preds;
      List updated_curr_preds;
      NumericVector new_mean;
      lowest_BIC=min(overall_lik2);
      
      
      //NumericMatrix curr_resids(resids.nrow(),resids.ncol());
      
      
      List temp(table_subset_curr_round.size());
      
      cp_mat_list=temp;
      
      
      
      for(int k=0;k<table_subset_curr_round.size();k++){
        //NumericVector terminal_nodes;
        
        if(parent_curr_round[k]==-1){
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a);
        }else{    
          new_mean=update_mean_var(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a);
        }  
        
        NumericVector terminal_nodes=find_term_nodes(table_subset_curr_round[k]);
        updated_curr_preds=update_predictions(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,D1.n_rows);
        NumericVector test_res;
        
        if(parent_curr_round[k]==-1){
          test_res=resids(_,0); 
        }else{
          test_res=resids(_,parent_curr_round[k]);
        }
        
        //NumericVector curr_test_res=updated_curr_preds[1];
        //Rcout << "Line 2548. j = " << j << ". \n";
        
        // if(parent_curr_round[k]==-1){
        //   curr_resids(_,0)=test_res-curr_test_res;
        // }else{
        //   curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;
        // }
        //NumericVector temp_curr_resids=test_res-curr_test_res;
        
        //Loop over terminal nodes
        List Tempcpmatlist(terminal_nodes.size());
        for(int nodeind=0;nodeind<terminal_nodes.size();nodeind++){
          
          //IntegerVector term_obs=wrap(find_term_obs(mat_subset_curr_round[k],terminal_nodes[nodeind]));
          arma::uvec term_obsarma=find_term_obs(mat_subset_curr_round[k],terminal_nodes[nodeind]);
          //IntegerVector term_obs=wrap(term_obsarma);
          IntegerVector term_obs=wrap(arma::conv_to<arma::ivec>::from(term_obsarma));
          
          
          
          //Rcout << "term_obs = " << term_obs << ".\n";
          //Rcout << "test_res= " << test_res << ".\n";
          
          //Rcout << "test_res[term_obs] = " << test_res[term_obs] << ".\n";
          NumericVector tempsubset = test_res[term_obs];
          //Rcout << "tempsubset = " << tempsubset << ".\n";
          
          
          List cp_mat_list1;
          arma::mat tempdata_subset = D1.rows(term_obsarma);
          if(gridpoint==0){
            //cp_mat_list1=make_pelt_cpmat(wrap(tempdata_subset),temp_curr_resids[term_obs],pen,num_cp);
            cp_mat_list1=make_pelt_cpmat(wrap(tempdata_subset),tempsubset,pen,num_cp);
          }else{
            //cp_mat_list1=make_gridpoint_cpmat(wrap(tempdata_subset),temp_curr_resids[term_obs],gridsize,num_cp);
            cp_mat_list1=make_gridpoint_cpmat(wrap(tempdata_subset),tempsubset,gridsize,num_cp);
          }
          Tempcpmatlist[nodeind]=cp_mat_list1[0];
        }
        cp_mat_list[k]=Tempcpmatlist;
        
      }
      
      
      
      //List temp(0);
      
      //cp_mat_list=temp;
      
      // for(int f=0;f<curr_resids.ncol();f++){
      //   List cp_mat_list1;
      //   if(gridpoint==0){
      //     cp_mat_list1=make_pelt_cpmat(wrap(D1),curr_resids(_,f),pen,num_cp);
      //   }else{
      //     cp_mat_list1=make_gridpoint_cpmat(wrap(D1),curr_resids(_,f),gridsize,num_cp);
      //   }
      //   
      //   cp_mat_list.push_back(cp_mat_list1[0]);      
      // }
      
    }   //end of if-statement split_rule_node==1
  }
  // Rcout << "Get to end of outer loop. \n";
  //Rcout << "overall_trees.size()= " << overall_trees.size() << " \n";
  //Rcout << "overall_mat.size()= " << overall_mat.size() << " \n";
  //Rcout << "overall_lik.size()= " << overall_lik.size() << " \n";
  //Rcout << "overall_parent.size()= " << overall_parent.size() << " \n";
  //Rcout << "overall_parent2.size()= " << overall_parent2.size() << " \n";
  
  overall_trees=resize(overall_trees,overall_count);
  overall_mat=resize(overall_mat,overall_count); 
  overall_lik.resize(overall_count);
  overall_parent.resize(overall_count);
  overall_predvecs=resize(overall_predvecs,overall_count);
  
  // Rcout << " line 8570. with update .\n";
  
  if(less_greedy==1){
    
    eval_model=evaluate_model_occams_window_exact(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)),
                                                  overall_predvecs);
    
    
    overall_lik2=eval_model[0];
    overall_trees=eval_model[1];
    overall_mat=eval_model[2];
    overall_count=overall_trees.size();
    overall_parent2=eval_model[3];
    overall_predvecs=eval_model[4];
  }
  
  
  NumericVector temp_preds;
  NumericVector temp_true_preds;
  
  List updated_preds;
  NumericVector new_mean;  
  NumericMatrix overall_test_preds(test_data.nrow(),overall_trees.size());  
  NumericMatrix overallpreds(D1.n_rows,overall_trees.size());
  lowest_BIC=min(overall_lik2);
  
  NumericMatrix overallpreds_total(D1.n_rows,overall_trees.size());
  
  
  //Rcout << "Get to start of update mean loop. \n";
  //Rcout << "overall_trees.size()= " << overall_trees.size() << " \n";
  //Rcout << "overall_mat.size()= " << overall_mat.size() << " \n";
  //Rcout << "overall_lik.size()= " << overall_lik.size() << " \n";
  //Rcout << "overall_parent.size()= " << overall_parent.size() << " \n";
  //Rcout << "overall_parent2.size()= " << overall_parent2.size() << " \n";
  
  for(int k=0;k<overall_trees.size();k++){
    //NumericVector terminal_nodes;
    
    if(overall_parent2[k]==-1){
      //Rcout << "within update mean var. overall_parent2[k]= " << overall_parent2[k] << " \n";
      //Rcout << "k= " << k << " \n";
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,0),a);
    }else{
      //Rcout << "within update mean var. overall_parent2[k]= " << overall_parent2[k] << " \n";
      //Rcout << "k= " << k << " \n";
      
      //Rcout << "a= " << a << " \n";
      SEXP tree_checktype = overall_trees[k];
      if(is<NumericMatrix>(tree_checktype)){
        //Rcout << "tree_checktype is a NumericMatrix. \n";
      }else{
        //Rcout << "tree_checktype is NOT a NumericMatrix. \n";
      }
      SEXP mat_checktype = overall_mat[k];
      if(is<NumericMatrix>(mat_checktype)){
        //Rcout << "mat_checktype is a NumericMatrix. \n";
      }else{
        //Rcout << "mat_checktype is NOT a NumericMatrix. \n";
      }
      //  SEXP vec_checktype = resids(_,overall_parent2[k]);
      //  if(is<NumericVector>(vec_checktype)){
      //Rcout << "vec_checktype is a NumericMatrix. \n";
      //  }else{
      //Rcout << "vec_checktype is NOT a NumericMatrix. \n";
      //  }
      new_mean=update_mean_var(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a);
    }
    //Rcout << "Get past update mean_var. k= " << k << " \n";
    
    //terminal_nodes=find_term_nodes(overall_trees[k]);
    updated_preds=update_predictions(overall_trees[k],overall_mat[k],new_mean,D1.n_rows);
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs(test_data,overall_trees[k]//,new_mean
    );
    temp_preds=updated_preds[1];
    overallpreds(_,k)=temp_preds;
    if(is_test_data)overall_test_preds(_,k)=test_preds;
    
    
    temp_true_preds=overall_predvecs[k];
    overallpreds_total(_,k)=temp_true_preds;    
    
  }
  // Rcout << "Get to end of update mean loop. \n";
  
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
  //arma::colvec predicted_values=sum(M1,1);
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  //arma::colvec predicted_test_values=sum(M2,1);
  //// Rcout << "Get to end of get_est_trees_sum. \n";
  // List ret(8);
  // ret[0]=overall_lik2;
  // ret[1]=overall_trees;
  // ret[2]=overall_mat;
  // ret[3]=predicted_values;
  // ret[4]=overall_parent2;
  // ret[5]=wrap(M1);
  // ret[6]=lowest_BIC;
  // ret[7]=wrap(M2);
  List ret(8);
  ret[0]=overall_lik2;
  ret[1]=overall_trees;
  ret[2]=overall_mat;
  ret[3]=overall_parent2;
  ret[4]=wrap(M1);
  ret[5]=lowest_BIC;
  ret[6]=wrap(M2);
  ret[7]=overallpreds_total;
  return(ret);
}
//######################################################################################################################//
// [[Rcpp::export]]
NumericVector scale_response(double a,double b,double c,double d,NumericVector y){
  NumericVector y_scaled = -((-b*c+a*d)/(-a+b))+((-c+d)*y/(-a+b));
  
  return(y_scaled);
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericVector get_original(double low,double high,double sp_low,double sp_high,NumericVector sum_preds){
  NumericVector original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);
  
  return(original_y);
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec get_original_arma(double low,double high,double sp_low,double sp_high,arma::vec sum_preds){
  arma::vec original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);
  
  return(original_y);
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec get_original_TE_arma(double low,double high,double sp_low,double sp_high,arma::vec sum_preds){
  arma::vec original_y=sum_preds*((-low+high)/(-sp_low+sp_high));
  
  return(original_y); // reverse scaling of predictions of scaled variable (??)
}
//######################################################################################################################//

// [[Rcpp::export]]
double get_original_TE_double(double low,double high,double sp_low,double sp_high,double sum_preds){
  double original_y=sum_preds*((-low+high)/(-sp_low+sp_high)); 
  
  return(original_y); // reverse scaling of predictions of scaled variable (??)
}
//###########################################################################################################################//

// [[Rcpp::export]]

List get_termobs_test_data(NumericMatrix test_data,NumericMatrix tree_data) {
  // Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.
  
  //test_data is a nxp matrix with the same variable names as the training data the model was built on
  
  //tree_data is the tree table with the tree information i.e. split points and split variables and terminal node mean values
  
  //term_node_means is a vector storing the terminal node mean values
  
  arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
  arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);
  //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);
  NumericVector terminal_nodes=find_term_nodes(tree_data);
  //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
  //NumericVector tree_predictions;
  
  //now for each internal node find the observations that belong to the terminal nodes
  
  //NumericVector predictions(test_data.nrow());
  List term_obs(terminal_nodes.size());
  if(terminal_nodes.size()==1){
    //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
    //predictions=rep(nodemean,test_data.nrow());
    IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
    term_obs[0]= temp_obsvec;
  }
  else{
    for(int i=0;i<terminal_nodes.size();i++){
      //arma::mat subdata=testd;
      int curr_term=terminal_nodes[i];
      int row_index;
      int term_node=terminal_nodes[i];
      
      if(curr_term % 2==0){
        //term node is left daughter
        row_index=terminal_nodes[i];
      }else{
        //term node is right daughter
        row_index=terminal_nodes[i]-1;
      }
      
      //save the left and right node data into arma uvec
      
      arma::vec left_nodes=arma_tree.col(0);
      arma::vec right_nodes=arma_tree.col(1);
      arma::mat node_split_mat;    
      node_split_mat.set_size(0,3);
      
      while(row_index!=1){
        //for each terminal node work backwards and see if the parent node was a left or right node
        //append split info to a matrix 
        int rd=0;
        arma::uvec parent_node=arma::find(left_nodes == term_node);
        
        if(parent_node.size()==0){
          parent_node=arma::find(right_nodes == term_node);
          rd=1;
        }
        
        //want to cout parent node and append to node_split_mat
        
        node_split_mat.insert_rows(0,1);
        node_split_mat(0,0)=tree_data(parent_node[0],2);
        node_split_mat(0,1)=tree_data(parent_node[0],3);
        node_split_mat(0,2)=rd;     
        row_index=parent_node[0]+1;
        term_node=parent_node[0]+1;
      }
      
      //once we have the split info, loop through rows and find the subset indexes for that terminal node!
      //then fill in the predicted value for that tree
      //double prediction = tree_data(term_node,5);
      arma::uvec pred_indices;
      int split= node_split_mat(0,0)-1;
      arma::vec tempvec = testd.col(split);
      double temp_split = node_split_mat(0,1);
      
      if(node_split_mat(0,2)==0){
        pred_indices = arma::find(tempvec <= temp_split);
      }else{
        pred_indices = arma::find(tempvec > temp_split);      
      }
      
      arma::uvec temp_pred_indices;
      arma::vec data_subset = testd.col(split);
      data_subset=data_subset.elem(pred_indices);
      
      //now loop through each row of node_split_mat
      int n=node_split_mat.n_rows;
      
      for(int j=1;j<n;j++){
        int curr_sv=node_split_mat(j,0);
        double split_p = node_split_mat(j,1);
        
        data_subset = testd.col(curr_sv-1);
        data_subset=data_subset.elem(pred_indices);
        
        if(node_split_mat(j,2)==0){
          //split is to the left
          temp_pred_indices=arma::find(data_subset <= split_p);
        }else{
          //split is to the right
          temp_pred_indices=arma::find(data_subset > split_p);
        }
        pred_indices=pred_indices.elem(temp_pred_indices);
        
        if(pred_indices.size()==0){
          continue;
        }
        
      }
      
      //double nodemean=tree_data(terminal_nodes[i]-1,5);
      IntegerVector predind=as<IntegerVector>(wrap(arma::conv_to<arma::ivec>::from(pred_indices)));
      //predictions[predind]= nodemean;
      term_obs[i]=predind;
    } 
  }
  //List ret(1);
  //ret[0] = term_obs;
  
  //ret[0] = terminal_nodes;
  //ret[1] = term_obs;
  //ret[2] = predictions;
  return(term_obs);
}

//###########################################################################################################################//

// [[Rcpp::export]]

arma::field<arma::uvec> get_termobs_test_data_fields(NumericMatrix test_data,NumericMatrix tree_data) {
  // Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.
  
  //test_data is a nxp matrix with the same variable names as the training data the model was built on
  
  //tree_data is the tree table with the tree information i.e. split points and split variables and terminal node mean values
  
  //term_node_means is a vector storing the terminal node mean values
  
  arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
  arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);
  //NumericVector internal_nodes=find_internal_nodes_gs(tree_data);
  NumericVector terminal_nodes=find_term_nodes(tree_data);
  //arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
  //NumericVector tree_predictions;
  
  //now for each internal node find the observations that belong to the terminal nodes
  
  //NumericVector predictions(test_data.nrow());
  //List term_obs(terminal_nodes.size());
  arma::field<arma::uvec> term_obsF(terminal_nodes.size());
  
  if(terminal_nodes.size()==1){
    //double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
    //predictions=rep(nodemean,test_data.nrow());
    IntegerVector temp_obsvec = seq_len(test_data.nrow())-1;
    term_obsF(0)= Rcpp::as<arma::uvec>(temp_obsvec);
  }
  else{
    for(int i=0;i<terminal_nodes.size();i++){
      //arma::mat subdata=testd;
      int curr_term=terminal_nodes[i];
      int row_index;
      int term_node=terminal_nodes[i];
      
      if(curr_term % 2==0){
        //term node is left daughter
        row_index=terminal_nodes[i];
      }else{
        //term node is right daughter
        row_index=terminal_nodes[i]-1;
      }
      
      //save the left and right node data into arma uvec
      
      arma::vec left_nodes=arma_tree.col(0);
      arma::vec right_nodes=arma_tree.col(1);
      arma::mat node_split_mat;    
      node_split_mat.set_size(0,3);
      
      while(row_index!=1){
        //for each terminal node work backwards and see if the parent node was a left or right node
        //append split info to a matrix 
        int rd=0;
        arma::uvec parent_node=arma::find(left_nodes == term_node);
        
        if(parent_node.size()==0){
          parent_node=arma::find(right_nodes == term_node);
          rd=1;
        }
        
        //want to cout parent node and append to node_split_mat
        
        node_split_mat.insert_rows(0,1);
        node_split_mat(0,0)=tree_data(parent_node[0],2);
        node_split_mat(0,1)=tree_data(parent_node[0],3);
        node_split_mat(0,2)=rd;     
        row_index=parent_node[0]+1;
        term_node=parent_node[0]+1;
      }
      
      //once we have the split info, loop through rows and find the subset indexes for that terminal node!
      //then fill in the predicted value for that tree
      //double prediction = tree_data(term_node,5);
      arma::uvec pred_indices;
      int split= node_split_mat(0,0)-1;
      arma::vec tempvec = testd.col(split);
      double temp_split = node_split_mat(0,1);
      
      if(node_split_mat(0,2)==0){
        pred_indices = arma::find(tempvec <= temp_split);
      }else{
        pred_indices = arma::find(tempvec > temp_split);      
      }
      
      arma::uvec temp_pred_indices;
      arma::vec data_subset = testd.col(split);
      data_subset=data_subset.elem(pred_indices);
      
      //now loop through each row of node_split_mat
      int n=node_split_mat.n_rows;
      
      for(int j=1;j<n;j++){
        int curr_sv=node_split_mat(j,0);
        double split_p = node_split_mat(j,1);
        
        data_subset = testd.col(curr_sv-1);
        data_subset=data_subset.elem(pred_indices);
        
        if(node_split_mat(j,2)==0){
          //split is to the left
          temp_pred_indices=arma::find(data_subset <= split_p);
        }else{
          //split is to the right
          temp_pred_indices=arma::find(data_subset > split_p);
        }
        pred_indices=pred_indices.elem(temp_pred_indices);
        
        if(pred_indices.size()==0){
          continue;
        }
        
      }
      
      //double nodemean=tree_data(terminal_nodes[i]-1,5);
      //IntegerVector predind=as<IntegerVector>(wrap(arma::conv_to<arma::ivec>::from((pred_indices)));
      //predictions[predind]= nodemean;
      //term_obsF(i)=Rcpp::as<arma::uvec>(predind);
      term_obsF(i)=pred_indices;
    } 
  }
  //List ret(1);
  //ret[0] = term_obs;
  
  //ret[0] = terminal_nodes;
  //ret[1] = term_obs;
  //ret[2] = predictions;
  return(term_obsF);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_termobs_testdata_overall(List overall_sum_trees,NumericMatrix test_data){
  
  //List overall_term_nodes_trees(overall_sum_trees.size());
  List overall_term_obs_trees(overall_sum_trees.size());
  //List overall_predictions(overall_sum_trees.size());
  
  for(int i=0;i<overall_sum_trees.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = overall_sum_trees[i];
    
    NumericVector test_preds_sum_tree;
    if(is<List>(s)){
      //if current set of trees contains more than one tree...usually does!
      List sum_tree=overall_sum_trees[i];        
      //save all info in list of list format the same as the trees.       
      //List term_nodes_trees(sum_tree.size());
      List term_obs_trees(sum_tree.size());
      //NumericMatrix predictions(num_obs,sum_tree.size());
      
      for(int k=0;k<sum_tree.size();k++){          
        NumericMatrix tree_table=sum_tree[k];
        List tree_info=get_termobs_test_data(test_data, tree_table) ;
        //NumericVector term_nodes=tree_info[0];
        //term_nodes_trees[k]=term_nodes;
        term_obs_trees[k]=tree_info;
        //umericVector term_preds=tree_info[2];
        //predictions(_,k)=term_preds;
      } 
      //overall_term_nodes_trees[i]=term_nodes_trees;
      overall_term_obs_trees[i]= term_obs_trees;
      //overall_predictions[i]=predictions;
    }else{
      NumericMatrix sum_tree=overall_sum_trees[i];
      List tree_info=get_termobs_test_data(test_data, sum_tree) ;
      //overall_term_nodes_trees[i]=tree_info[0];
      List term_obs_trees(1);
      term_obs_trees[0]=tree_info ;
      //NumericVector term_preds=tree_info[2];
      //NumericVector predictions=term_preds;   
      overall_term_obs_trees[i]= term_obs_trees;
      //overall_predictions[i]=predictions;
    }  
  }    
  //List ret(1);
  //ret[0]=overall_term_nodes_trees;
  //ret[0]=overall_term_obs_trees;
  //ret[2]=overall_predictions;
  return(overall_term_obs_trees);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]

arma::field<arma::field<arma::field<arma::uvec>>> get_termobs_testdata_fields_overall(List overall_sum_trees,NumericMatrix test_data){
  
  //List overall_term_nodes_trees(overall_sum_trees.size());
  //  List overall_term_obs_trees(overall_sum_trees.size());
  //List overall_predictions(overall_sum_trees.size());
  
  arma::field<arma::field<arma::field<arma::uvec>>> OutputF(overall_sum_trees.size()) ;
  
  for(int i=0;i<overall_sum_trees.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = overall_sum_trees[i];
    
    NumericVector test_preds_sum_tree;
    if(is<List>(s)){
      //if current set of trees contains more than one tree...usually does!
      List sum_tree=overall_sum_trees[i];        
      //save all info in list of list format the same as the trees.       
      //List term_nodes_trees(sum_tree.size());
      //List term_obs_trees(sum_tree.size());
      arma::field<arma::field<arma::uvec>> term_obs_treesF(sum_tree.size());
      //NumericMatrix predictions(num_obs,sum_tree.size());
      
      for(int k=0;k<sum_tree.size();k++){          
        NumericMatrix tree_table=sum_tree[k];
        arma::field<arma::uvec> tree_infoF=get_termobs_test_data_fields(test_data, tree_table) ;
        //NumericVector term_nodes=tree_info[0];
        //term_nodes_trees[k]=term_nodes;
        term_obs_treesF(k)=tree_infoF;
        //umericVector term_preds=tree_info[2];
        //predictions(_,k)=term_preds;
      } 
      //overall_term_nodes_trees[i]=term_nodes_trees;
      OutputF(i)= term_obs_treesF;
      //overall_predictions[i]=predictions;
    }else{
      NumericMatrix sum_tree=overall_sum_trees[i];
      arma::field<arma::uvec> tree_infoF=get_termobs_test_data_fields(test_data, sum_tree) ;
      //overall_term_nodes_trees[i]=tree_info[0];
      //List term_obs_trees(1);
      arma::field<arma::field<arma::uvec>> term_obs_treesF(1);
      term_obs_treesF(0)=tree_infoF ;
      //NumericVector term_preds=tree_info[2];
      //NumericVector predictions=term_preds;   
      OutputF(i)= term_obs_treesF;
      //overall_predictions[i]=predictions;
    }  
  }    
  //List ret(1);
  //ret[0]=overall_term_nodes_trees;
  //ret[0]=overall_term_obs_trees;
  //ret[2]=overall_predictions;
  return(OutputF);
}
//###########################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat get_J_test(List curr_termobs,NumericVector tree_term_nodes, int n){
  //this function will make a binary nxb matrix where each column assigns observations to terminal nodes
  
  //create J matrix with correct dimensions and fill with zeros
  arma::mat Jmat(n, tree_term_nodes.size());
  Jmat.zeros();
  
  //for each terminal node get the observations associated with it and set column
  for(int i=0;i<tree_term_nodes.size();i++){
    //double tn=tree_term_nodes[i];
    //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
    //assign term_obs to the correct index of J
    IntegerVector term_obs2=curr_termobs[i];
    NumericVector obs_col(n);
    obs_col[term_obs2]=1;
    arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
    Jmat.col(i)= colmat;
    
    // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
    // colmat.elem(term_obs).fill(1);
    // Jmat.col(i)= colmat;
  }
  return(Jmat);
}

//###########################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat get_W_test(List sum_treetable ,
                     List termobs_testdata_onemodel,
                     int n){
  //this will take in a list of obs to node matrices for each tree in the sum make the J matrix assigning observations to terminal nodes
  //J is an nxb_j matrix. It will then iteratively append the J matrix to itself to get the overall W matrix which has dimensions nxsumb_j
  
  //create empty matrix to which we will append individual J matrices
  arma::mat W(n,0);
  int upsilon=0;
  for(int j=0;j<sum_treetable.size();j++){
    
    NumericMatrix curr_tree=sum_treetable[j];
    List curr_termobs=termobs_testdata_onemodel[j];
    NumericVector tree_term_nodes=find_term_nodes(curr_tree);
    int b_j=tree_term_nodes.size();
    //will make J as we go in BART-BMA no need to create it again here....
    arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
    W.insert_cols(upsilon,Jmat);
    upsilon+=b_j;
  }
  
  return(W);
  // ret[1]=mu_vec;
}


//###########################################################################################################################//


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector preds_bbma_lin_alg_outsamp(List overall_sum_trees,
                                         List overall_sum_mat,
                                         NumericVector y,
                                         NumericVector BIC_weights,
                                         int num_iter,int burnin,int num_obs,int num_test_obs,
                                         double a,double sigma,double mu_mu,double nu,
                                         double lambda,//List resids,
                                         NumericMatrix test_data){
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y); 
  
  List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size()); 
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size()); 
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  for(int i=0;i<overall_sum_trees.size();i++){
    
    
    arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    double b=Wmat.n_cols;
    arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
    arma::mat y_arma(num_obs,1);
    y_arma.col(0)=yvec;
    //get exponent
    //double expon=(n+nu)/2;
    //get y^Tpsi^{-1}y
    // arma::mat psi_inv=psi.i();
    
    
    //arma::mat yty=y_arma.t()*y_arma;
    
    //get t(y)inv(psi)J
    //arma::mat ytW=y_arma.t()*Wmat;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=Wmat.t()*Wmat;
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    //get t(J)inv(psi)y
    arma::mat third_term=Wmat.t()*y_arma;
    
    //get m^TV^{-1}m
    //arma::mat mvm= ytW*sec_term_inv*third_term;
    
    arma::vec preds_temp_arma= W_tilde*sec_term_inv*third_term;
    
    
    
    
    
    
    
    
    double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    preds_all_models_arma.col(i)=preds_temp_arma*weight;
    
    
  }
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(preds_all_models_arma,1);
  
  
  NumericVector orig_preds=get_original(min(y),max(y),-0.5,0.5,wrap(predicted_values)) ;
  
  
  
  return(orig_preds);
}

//###########################################################################################################################//


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector preds_bbma_lin_alg_insamp(List overall_sum_trees,List overall_sum_mat,
                                        NumericVector y,NumericVector BIC_weights,
                                        int num_iter,int burnin,int num_obs,
                                        double a,double sigma,double mu_mu,double nu,
                                        double lambda//,List resids
){
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y); 
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size()); 
  arma::mat preds_all_models_arma(num_obs,BIC_weights.size()); 
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  for(int i=0;i<overall_sum_trees.size();i++){
    
    arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    //arma::mat W_tilde=W(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    double b=Wmat.n_cols;
    arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
    arma::mat y_arma(num_obs,1);
    y_arma.col(0)=yvec;
    //get exponent
    //double expon=(n+nu)/2;
    //get y^Tpsi^{-1}y
    // arma::mat psi_inv=psi.i();
    
    
    //arma::mat yty=y_arma.t()*y_arma;
    
    //get t(y)inv(psi)J
    //arma::mat ytW=y_arma.t()*Wmat;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=Wmat.t()*Wmat;
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    //get t(J)inv(psi)y
    arma::mat third_term=Wmat.t()*y_arma;
    
    //get m^TV^{-1}m
    //arma::mat mvm= ytW*sec_term_inv*third_term;
    
    arma::vec preds_temp_arma= Wmat*sec_term_inv*third_term;
    
    
    
    
    
    
    
    
    double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    preds_all_models_arma.col(i)=preds_temp_arma*weight;
    
    
  }
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(preds_all_models_arma,1);
  
  NumericVector orig_preds=get_original(min(y),max(y),-0.5,0.5,wrap(predicted_values)) ;
  
  
  
  //return(wrap(predicted_values));
  return(orig_preds);
}



//###########################################################################################################################//


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List mean_vars_lin_alg_outsamp(List overall_sum_trees,
                               List overall_sum_mat,
                               NumericVector y,
                               NumericVector BIC_weights,
                               int num_iter,int burnin,int num_obs,int num_test_obs,
                               double a,double sigma,double mu_mu,double nu,
                               double lambda,//List resids,
                               NumericMatrix test_data){
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y); 
  
  List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size()); 
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size()); 
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  List covar_matrices(BIC_weights.size());
  NumericVector model_weights(BIC_weights.size());
  
  for(int i=0;i<overall_sum_trees.size();i++){
    //Rcout << "Line 4111";
    
    
    arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    double b=Wmat.n_cols;
    arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
    arma::mat y_arma(num_obs,1);
    y_arma.col(0)=yvec;
    //get exponent
    //double expon=(n+nu)/2;
    //get y^Tpsi^{-1}y
    // arma::mat psi_inv=psi.i();
    
    
    arma::mat yty=y_arma.t()*y_arma;
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*Wmat;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=Wmat.t()*Wmat;
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    //get t(J)inv(psi)y
    arma::mat third_term=Wmat.t()*y_arma;
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    arma::mat w_tilde_M_inv =  W_tilde*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    //arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    arma::mat I_test(num_test_obs,num_test_obs);
    I_test=I_test.eye();
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    
    
    //Rcout << "Line 4162";
    
    
    double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    //Rcout << "Line 4167";
    
    model_weights[i]=weight;
    //Rcout << "Line 4170";
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*weight;
    preds_all_models_arma.col(i)=preds_temp_arma;
    //Rcout << "Line 4173";
    
    covar_matrices[i]= wrap(covar_t);
    //Rcout << "Line 4176";
    
    
    
    
  }
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  
  //NumericVector orig_preds=get_original(min(y),max(y),-0.5,0.5,wrap(predicted_values)) ;
  
  List ret(4);
  ret[0]=wrap(predicted_values);
  ret[1]=model_weights;
  ret[2]=wrap(preds_all_models_arma.t());
  ret[3]=covar_matrices;
  
  return(ret);
}
//###########################################################################################################################//


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List mean_vars_lin_alg_insamp(List overall_sum_trees,
                              List overall_sum_mat,
                              NumericVector y,
                              NumericVector BIC_weights,
                              int num_iter,int burnin,int num_obs,int num_test_obs,
                              double a,double sigma,double mu_mu,double nu,
                              double lambda//,List resids,
){
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y); 
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size()); 
  arma::mat preds_all_models_arma(num_obs,BIC_weights.size()); 
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  List covar_matrices(BIC_weights.size());
  NumericVector model_weights(BIC_weights.size());
  
  for(int i=0;i<overall_sum_trees.size();i++){
    
    
    arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    double b=Wmat.n_cols;
    arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
    arma::mat y_arma(num_obs,1);
    y_arma.col(0)=yvec;
    //get exponent
    //double expon=(n+nu)/2;
    //get y^Tpsi^{-1}y
    // arma::mat psi_inv=psi.i();
    
    
    arma::mat yty=y_arma.t()*y_arma;
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*Wmat;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=Wmat.t()*Wmat;
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    //get t(J)inv(psi)y
    arma::mat third_term=Wmat.t()*y_arma;
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    
    arma::mat W_M_inv = Wmat*sec_term_inv;
    arma::vec preds_temp_arma= W_M_inv*third_term;
    
    //arma::mat I_test(num_test_obs,num_test_obs);
    //I_test=I_test.eye();
    
    arma::mat covar_t=as_scalar((1/(nu+num_obs))*(nu*lambda+yty-mvm))*(W_M_inv*(Wmat.t()));
    
    
    
    
    double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    model_weights[i]=weight;
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*weight;
    preds_all_models_arma.col(i)=preds_temp_arma;
    covar_matrices[i]= covar_t;
    
    
    
    
  }
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  
  //NumericVector orig_preds=get_original(min(y),max(y),-0.5,0.5,wrap(predicted_values)) ;
  
  List ret(4);
  ret[0]=wrap(predicted_values);
  ret[1]=model_weights;
  ret[2]=wrap(preds_all_models_arma.t());
  ret[3]=covar_matrices;
  
  return(ret);
}


//###########################################################################################################################//


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List BART_BMA_sumLikelihood(double less_greedy,double spike_tree, int s_t_hyperprior, double p_s_t, double a_s_t, double b_s_t, double num_obs, double num_vars, double lambda_poisson,
                            NumericMatrix data,NumericVector y,double start_mean,double start_sd,
                            double a,double mu,double nu,double lambda,double c,
                            double sigma_mu,double pen,int num_cp,NumericMatrix test_data,int num_rounds,
                            double alpha,double beta,bool split_rule_node,bool gridpoint,int maxOWsize,int num_splits,int gridsize, bool zero_split,bool only_max_num_trees,
                            unsigned int min_num_obs_for_split, unsigned int min_num_obs_after_split, int exact_residuals){
  
  
  bool is_test_data=0;
  if(test_data.nrow()>0){
    is_test_data=1;
  }
  if(y.size() !=data.nrow()){
    if(y.size()<data.nrow()){
      throw std::range_error("Response length is smaller than the number of observations in the data"); 
    }else{
      throw std::range_error("Response length is greater than the number of observations in the data"); 
    }
  }
  //check test data has the same number of variables as training data
  if(test_data.nrow()>0 && (data.ncol() != test_data.ncol())){
    throw std::range_error("Test data and training data must have the same number of variables. BART BMA assumes variables are in the same order."); 
  }
  
  //check value of c is greater than 1!
  //	if(c<1){
  //		throw std::range_error("Value of Occam's Window has to be greater than 0."); 
  //	}
  if(num_cp<0 || num_cp>100){
    throw std::range_error("Value of num_cp should be a value between 1 and 100."); 
  }
  //NumericMatrix treetable=start_tree(mu,sigma_mu);
  NumericMatrix treetable=start_tree2();
  NumericMatrix treemat=start_matrix(data.nrow());
  //initialize the tree table and matrix
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  //initialize the tree table and matrix
  arma::mat D1(data.begin(), data.nrow(), data.ncol(), false);
  double n=D1.n_rows;
  //	double lik=likelihood_function(y_scaled,treetable,treemat,a,mu,nu,lambda);
  double lik=as<double>(likelihood_function2_exact(y_scaled,treetable,treemat,a,mu,nu,lambda)[0]);
  double tree_prior=get_tree_prior(spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,treetable,treemat,alpha,beta);
  //double lowest_BIC=-2*(lik+log(tree_prior))+1*log(n);
  double lowest_BIC=-2*(lik+log(tree_prior));
  // Rcout << "Initial lowest BIC = " <<  lowest_BIC << ".\n"; 
  List best_subset;  
  List tree_table;
  List tree_mat;
  tree_table.push_back(treetable);
  tree_mat.push_back(treemat);
  List CART_BMA;
  arma::mat r;
  arma::colvec yarma=clone(y_scaled);
  r.insert_cols(0,yarma);
  NumericMatrix resids=wrap(r);
  //int first_round;
  int first_round_break=0;
  
  //List overall_trees(num_rounds);
  List overall_mat;
  List overall_lik;
  NumericMatrix prev_round_preds;  
  NumericVector prev_round_BIC;
  NumericVector prev_round_BIC2;
  arma::mat prev_round_preds2;
  NumericMatrix prev_round_test_preds;
  arma::mat prev_round_test_preds2;
  arma::mat overall_overall_sum_test_preds;
  arma::colvec predicted_test_values;
  List prev_sum_trees;
  List prev_sum_tree_resids;
  List prev_sum_trees_mat;  
  List cp_mat_list;
  //int oo_size=300;
  //List overall_overall_sum_trees(oo_size);
  //List overall_overall_sum_tree_resids(oo_size);
  //List overall_overall_sum_trees_mat(oo_size);
  //List overall_overall_sum_BIC(oo_size);
  
  List overall_overall_sum_trees;
  List overall_overall_sum_tree_resids;
  List overall_overall_sum_trees_mat;
  NumericVector overall_overall_sum_BIC;
  
  //int oo_count=0;
  arma::mat overall_overall_sum_preds;
  IntegerVector prev_par;
  arma::colvec predicted_values;
  
  for(int j=0;j<num_rounds;j++){
    if(j>0){
      if(only_max_num_trees==1){
        lowest_BIC=100000;
      }
    }
    int overall_size=300;
    List overall_sum_trees(overall_size);
    List overall_sum_trees_mat(overall_size);
    List overall_sum_tree_resids(overall_size);
    int overall_count=0;
    IntegerVector parent;   
    NumericVector curr_round_lik;
    List curr_round_trees;
    List curr_round_mat;
    NumericVector curr_BIC;
    IntegerVector curr_round_parent;
    NumericVector overall_sum_BIC;
    arma::mat overall_sum_preds;    
    arma::mat overall_sum_test_preds;
    
    if(j==0){
      parent.push_back(0);
      //first_round=1;
    }//else{
    //  first_round=0;
    //}
    //// Rcout << "resids = " << resids << ".\n";
    
    List resids_cp_mat(resids.ncol());
    int resids_count=0;
    std::vector<int> err_list(resids.ncol());
    //get best splits
    for(int f=0;f<resids.ncol();f++){
      if(gridpoint==0){
        cp_mat_list=make_pelt_cpmat(data,resids(_,f),pen,num_cp);
      }else{
        cp_mat_list=make_gridpoint_cpmat(data,resids(_,f),gridsize,num_cp);
        //NumericVector tempresidvec=resids(_,f);
        // cp_mat_list=make_gridpoint_cpmat_arma(D1,Rcpp::as<arma::vec>(tempresidvec),gridsize,num_cp);
        
      }
      resids_cp_mat[resids_count]=cp_mat_list[0];
      err_list[resids_count]=cp_mat_list[1];
      resids_count++;
    }
    resids_cp_mat=resize(resids_cp_mat,resids_count);
    err_list.resize(resids_count);     
    parent=seq_len(tree_table.size())-1;
    if(is_true(all(as<IntegerVector>(wrap(err_list))==1))){
      if(j==0){
        throw std::range_error("No split points could be found to grow trees");
      }else{
        throw std::range_error("No trees can be grown for the number of iterations desired, as no splits were found.Please try fewer iterations.");
      }
    } 
    
    //Rcout << "Get to before get_best_trees in outer loop number " << j << " . \n";
    
    
    if(exact_residuals==1){
      //get current set of trees.
      if(j==0){
        if(split_rule_node==1){
          CART_BMA=get_best_trees_update_splits_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1, resids, a,mu,nu,lambda,c,sigma_mu,tree_table,tree_mat,lowest_BIC,//first_round,
                                                      parent,resids_cp_mat,//as<IntegerVector>(wrap(err_list)),
                                                      test_data,alpha,beta,is_test_data,pen,num_cp,split_rule_node,gridpoint,maxOWsize,num_splits,gridsize,zero_split,
                                                      min_num_obs_for_split, min_num_obs_after_split);
          
        }else{
          CART_BMA=get_best_trees_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1, resids, a,mu,nu,lambda,c,sigma_mu,tree_table,tree_mat,lowest_BIC,//first_round,
                                        parent,resids_cp_mat,//as<IntegerVector>(wrap(err_list)),
                                        test_data,alpha,beta,is_test_data,pen,num_cp,split_rule_node,gridpoint,maxOWsize,num_splits,gridsize,zero_split,
                                        min_num_obs_for_split, min_num_obs_after_split);
        }
        
      }else{
        
        // Rcout << "Line 9426. Length of tree_table = " << tree_table.size() << " . \n";
        // Rcout << "Line 9427. Length of prev_sum_trees = " << prev_sum_trees.size() << " . \n";
        
        //if j >0 then sum of trees become a list so need to read in list and get likelihood for each split point and terminal node
        if(split_rule_node==1){
          CART_BMA=get_best_trees_sum_update_splits_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1, resids, a,mu,nu,lambda,c,sigma_mu,tree_table,tree_mat,lowest_BIC,//first_round,
                                                          parent,resids_cp_mat,as<IntegerVector>(wrap(err_list)),test_data,alpha,beta,is_test_data,pen,num_cp,split_rule_node,gridpoint,maxOWsize,prev_sum_trees,prev_sum_trees_mat,y_scaled,num_splits,gridsize,zero_split,
                                                          min_num_obs_for_split, min_num_obs_after_split);
          
          
        }else{
          CART_BMA=get_best_trees_sum_exact(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1, resids, a,mu,nu,lambda,c,sigma_mu,tree_table,tree_mat,lowest_BIC,//first_round,
                                            parent,resids_cp_mat,as<IntegerVector>(wrap(err_list)),test_data,alpha,beta,is_test_data,pen,num_cp,split_rule_node,gridpoint,maxOWsize,prev_sum_trees,prev_sum_trees_mat,y_scaled,num_splits,gridsize,zero_split,
                                            min_num_obs_for_split, min_num_obs_after_split);
        }
      }
    }else{
      if(j==0){
        if(split_rule_node==1){
          CART_BMA=get_best_trees_update_splits(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1, resids, a,mu,nu,lambda,c,sigma_mu,tree_table,tree_mat,lowest_BIC,//first_round,
                                                parent,resids_cp_mat,//as<IntegerVector>(wrap(err_list)),
                                                test_data,alpha,beta,is_test_data,pen,num_cp,split_rule_node,gridpoint,maxOWsize,num_splits,gridsize,zero_split,
                                                min_num_obs_for_split, min_num_obs_after_split);
          
        }else{
          CART_BMA=get_best_trees(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1, resids, a,mu,nu,lambda,c,sigma_mu,tree_table,tree_mat,lowest_BIC,//first_round,
                                  parent,resids_cp_mat,//as<IntegerVector>(wrap(err_list)),
                                  test_data,alpha,beta,is_test_data,pen,num_cp,split_rule_node,gridpoint,maxOWsize,num_splits,gridsize,zero_split,
                                  min_num_obs_for_split, min_num_obs_after_split);
        }
        
      }else{
        
        // Rcout << "Line 9459. Length of tree_table = " << tree_table.size() << " . \n";
        // Rcout << "Line 9460. Length of prev_sum_trees = " << prev_sum_trees.size() << " . \n";
        
        //if j >0 then sum of trees become a list so need to read in list and get likelihood for each split point and terminal node
        if(split_rule_node==1){
          CART_BMA=get_best_trees_sum_update_splits(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1, resids, a,mu,nu,lambda,c,sigma_mu,tree_table,tree_mat,lowest_BIC,//first_round,
                                                    parent,resids_cp_mat,as<IntegerVector>(wrap(err_list)),test_data,alpha,beta,is_test_data,pen,num_cp,split_rule_node,gridpoint,maxOWsize,prev_sum_trees,prev_sum_trees_mat,y_scaled,num_splits,gridsize,zero_split,
                                                    min_num_obs_for_split, min_num_obs_after_split);
          
          
        }else{
          CART_BMA=get_best_trees_sum(less_greedy,spike_tree, s_t_hyperprior, p_s_t, a_s_t,b_s_t,num_obs,num_vars,lambda_poisson,D1, resids, a,mu,nu,lambda,c,sigma_mu,tree_table,tree_mat,lowest_BIC,//first_round,
                                      parent,resids_cp_mat,as<IntegerVector>(wrap(err_list)),test_data,alpha,beta,is_test_data,pen,num_cp,split_rule_node,gridpoint,maxOWsize,prev_sum_trees,prev_sum_trees_mat,y_scaled,num_splits,gridsize,zero_split,
                                      min_num_obs_for_split, min_num_obs_after_split);
        }
      }
    }
    
    
    
    //Rcout << "Get past get_best_trees in outer loop number " << j << " . \n";
    
    curr_round_lik=CART_BMA[0];
    
    // Rcout << "Line 9625 outer loop number " << j << " . \n";
    
    curr_round_trees=CART_BMA[1];      
    //// Rcout << "Line 9628 outer loop number " << j << " . \n";
    curr_round_mat=CART_BMA[2];
    //Rcout << "Line 9630 outer loop number " << j << " . \n";
    curr_round_parent=CART_BMA[3];
    //Rcout << "Line 9632 outer loop number " << j << " . \n";
    NumericMatrix curr_round_preds=CART_BMA[4];
    //Rcout << "Line 9634 outer loop number " << j << " . \n";
    NumericMatrix curr_round_test_preds=CART_BMA[6];
    //Rcout << "Line 9636 outer loop number " << j << " . \n";
    curr_BIC=CART_BMA[5];
    // Rcout << "Line 9638 outer loop number " << j << " . \n";
    
    //Predictions from whole sum-of-tree models
    //Matrix dimensions: number of training observations by number of models
    //Each column is a separate model
    NumericMatrix curr_true_preds;
    if(exact_residuals==1){
      NumericMatrix temp_truepredmat = CART_BMA[7];
      curr_true_preds=temp_truepredmat;
    }
    
    // Rcout << "Line 9640 outer loop number " << j << " . \n";
    
    //
    //Rcout << "curr_round_preds = " << curr_round_preds << ".\n";
    
    //Rcout << "curr_true_preds = " << curr_true_preds << ".\n";
    
    if(curr_round_lik.size()==0) {
      //Rcout << "Break because curr_round_lik.size()==0 in outer loop number " << j << " . \n";
      if(j==0){
        first_round_break=1;
      }else{
        overall_overall_sum_trees=prev_sum_trees;
        overall_overall_sum_trees_mat=prev_sum_trees_mat;
        overall_overall_sum_tree_resids=prev_sum_tree_resids;
        overall_overall_sum_BIC=prev_round_BIC2;
        overall_overall_sum_preds=prev_round_preds2;
        if(is_test_data==1) overall_overall_sum_test_preds= prev_round_test_preds2;
      }
      
      break;
    } 
    
    if(curr_BIC[0]<lowest_BIC){
      lowest_BIC=curr_BIC[0];
    }
    tree_table=List();
    tree_mat=List();
    int lsize=curr_round_lik.size();
    tree_table=List(lsize);
    tree_mat=List(lsize);
    NumericMatrix temp_preds(n,curr_round_lik.size());
    NumericMatrix temp_test_preds(test_data.nrow(),curr_round_lik.size());
    NumericMatrix temp_resids(n,curr_round_lik.size());
    NumericVector temp_parent(curr_round_lik.size());
    NumericVector temp_BIC(curr_round_lik.size());
    List temp_sum_trees(lsize);
    List temp_sum_tree_resids(lsize);
    List temp_sum_trees_mat(lsize);
    
    
    // Rcout << "Line 9681 outer loop number " << j << " . \n";
    
    
    
    int count=0; 
    for(int k=0;k<curr_round_lik.size();k++){
      //tree_table[count]=start_tree(mu,sigma_mu);
      tree_table[count]=start_tree2();
      tree_mat[count]=start_matrix(n);
      if(j==0){
        //Rcout << "Line 4864. j=  " << j << " . \n";
        
        //
        if(exact_residuals==1){
          temp_preds(_,k)=curr_true_preds(_,k);
        }else{
          temp_preds(_,k)=curr_round_preds(_,k);
        }
        if(is_test_data==1) temp_test_preds(_,k)=curr_round_test_preds(_,k);
        
        temp_resids(_,k)=y_scaled-temp_preds(_,k);
        
        temp_parent[k]=-1;
        temp_BIC[k]=curr_round_lik[k];
        temp_sum_trees[count]=curr_round_trees[k];
        temp_sum_trees_mat[count]=curr_round_mat[k];
        temp_sum_tree_resids[count]=resids(_,0);
        
        //Rcout << "Line 4875. j=  " << j << " . \n";
        
      }else{
        NumericVector curr_temp_pred=curr_round_preds(_,k) + prev_round_preds(_,curr_round_parent[k]);
        NumericVector curr_temp_test_pred;
        
        if(is_test_data==1) {
          curr_temp_test_pred=curr_round_test_preds(_,k) + prev_round_test_preds(_,curr_round_parent[k]);
          temp_test_preds(_,k) = curr_temp_test_pred;
        }
        temp_BIC[k]=curr_round_lik[k];
        
        if(exact_residuals==1){
          temp_preds(_,k)=curr_true_preds(_,k);
          temp_resids(_,k)=y_scaled-curr_true_preds(_,k);
          
        }else{
          temp_preds(_,k)=curr_round_preds(_,k);
          temp_resids(_,k)=y_scaled-curr_temp_pred;
        }
        
        
        temp_parent[k] = k;
        temp_sum_trees[count]=curr_round_trees[k];
        temp_sum_trees_mat[count]=curr_round_mat[k];
        temp_sum_tree_resids[count]=resids(_,curr_round_parent[k]);        
      }
      count++;     
    }
    // Rcout << "Line 9738. j=  " << j << " . \n";
    
    // Removing unnecessary if statement. Same condition results in a break from the loop above
    //if(curr_round_lik.size()==0){
    //  throw std::range_error("No trees chosen in last round");
    //}
    for(int k=0;k<curr_round_lik.size();k++){
      int size_mat=300;
      List sum_of_trees(size_mat);
      List sum_of_tree_resids(size_mat);
      List sum_of_trees_mat(size_mat);
      int count=0;
      
      if(curr_round_parent[k]==-1){
      }else{
        if(j==1){
          NumericMatrix other_tree=prev_sum_trees[curr_round_parent[k]];
          NumericVector other_resids=prev_sum_tree_resids[curr_round_parent[k]];
          sum_of_trees[count]= other_tree;  
          sum_of_tree_resids[count]=other_resids-curr_round_preds(_,k);
          NumericMatrix other_mat=prev_sum_trees_mat[curr_round_parent[k]];
          sum_of_trees_mat[count]=other_mat;
          count++;
          
          if(count==(size_mat-1)){
            size_mat=size_mat*2;
            sum_of_trees=resize_bigger(sum_of_trees,size_mat);
            sum_of_tree_resids=resize_bigger(sum_of_tree_resids,size_mat);
            sum_of_trees_mat=resize_bigger(sum_of_trees_mat,size_mat);
          }
        }else{
          List other_tree=prev_sum_trees[curr_round_parent[k]];
          List other_tree_resids=prev_sum_tree_resids[curr_round_parent[k]];  
          List other_mat=prev_sum_trees_mat[curr_round_parent[k]];
          for(int f=0;f<other_tree.size();f++){
            if(is<NumericMatrix>(other_tree[f])){
            }else{
              throw std::range_error(" tree is not a numeric matrix!");
            }
            NumericMatrix treetoadd=other_tree[f];
            if(is<NumericVector>(other_tree_resids[f])){
            }else{
              throw std::range_error("other resids not a numeric matrix!");
            }
            NumericVector resids_prevroundtemp=other_tree_resids[f];
            NumericVector residstoadd=resids_prevroundtemp-curr_round_preds(_,k);
            
            if(is<NumericMatrix>(other_mat[f])){
              
            }else{
              throw std::range_error(" other mat not a numeric matrix!");
            }
            NumericMatrix mattoadd=other_mat[f];
            sum_of_trees[count]=treetoadd;
            sum_of_tree_resids[count]=residstoadd;
            sum_of_trees_mat[count]=mattoadd;
            count++;
            
            if(count==(size_mat-1)){
              size_mat=size_mat*2;
              sum_of_trees=resize_bigger(sum_of_trees,size_mat);
              sum_of_tree_resids=resize_bigger(sum_of_tree_resids,size_mat);
              sum_of_trees_mat=resize_bigger(sum_of_trees_mat,size_mat);
            }
          }
        }
        sum_of_trees[count]=temp_sum_trees[k];
        sum_of_tree_resids[count]=temp_sum_tree_resids[k];
        sum_of_trees_mat[count]=temp_sum_trees_mat[k];
        count++;
        
        if(count==(size_mat-1)){
          size_mat=size_mat*2;
          sum_of_trees=resize_bigger(sum_of_trees,size_mat);
          sum_of_trees_mat=resize_bigger(sum_of_trees_mat,size_mat);
          sum_of_tree_resids=resize_bigger(sum_of_tree_resids,size_mat);
        }
      }
      sum_of_trees=resize(sum_of_trees,count);
      sum_of_trees_mat=resize(sum_of_trees_mat,count);
      sum_of_tree_resids=resize(sum_of_tree_resids,count);
      
      if(curr_round_parent[k]!=-1){
        overall_sum_trees[overall_count]=sum_of_trees;
        overall_sum_tree_resids[overall_count]=sum_of_tree_resids;
        overall_sum_trees_mat[overall_count]=sum_of_trees_mat;
        overall_sum_BIC=temp_BIC;
        overall_sum_preds=Rcpp::as<arma::mat>(temp_preds);
        if(is_test_data==1) overall_sum_test_preds=Rcpp::as<arma::mat>(temp_test_preds);
        overall_count++;
        if(overall_count==(overall_size-1)){
          overall_size=overall_size*2;
          overall_sum_trees=resize_bigger(overall_sum_trees,overall_size);
          overall_sum_tree_resids=resize_bigger(overall_sum_tree_resids,overall_size);
          overall_sum_trees_mat=resize_bigger(overall_sum_trees_mat,overall_size);
        }
      }  
    }
    
    
    //// Rcout << "Line 4990 curr_round_lik.size()=" << curr_round_lik.size() << ".\n";
    //// Rcout << "Line 4991 overall_sum_trees.size()=" << overall_sum_trees.size() << ".\n";
    //// Rcout << "Line 4992 overall_count=" << overall_count << ".\n";
    // Rcout << "Line 9841. j=  " << j << " . \n";
    
    
    if(only_max_num_trees==0){
      //Rcout << "Get to checking for models from previous rounds in outer loop number " << j << " . \n";
      //check if there were any trees from the previous round that didn't have daughter trees grown.
      //create vector to count number of possible parents for previous round
      if(j>0){
        IntegerVector prev_par_no_child=match(prev_par,curr_round_parent);
        if(any(is_na(prev_par_no_child))){
          IntegerVector t4=ifelse(is_na(prev_par_no_child),1,0);
          for(int h=0;h<prev_par_no_child.size();h++){
            if(t4[h]==1){
              if(prev_round_BIC2[h]-lowest_BIC<=log(c)){
                SEXP s = prev_sum_trees[h];
                if(is<List>(s)){
                  List tree_no_child=prev_sum_trees[h];
                  List resids_no_child=prev_sum_tree_resids[h];
                  List treemat_no_child=prev_sum_trees_mat[h];
                  overall_sum_trees[overall_count]=tree_no_child;
                  overall_sum_tree_resids[overall_count]=resids_no_child;
                  overall_sum_trees_mat[overall_count]=treemat_no_child;
                  overall_count++;
                  
                  if(overall_count==(overall_size-1)){
                    overall_size=overall_size*2;
                    overall_sum_trees=resize_bigger(overall_sum_trees,overall_size);
                    overall_sum_tree_resids=resize_bigger(overall_sum_tree_resids,overall_size);
                    overall_sum_trees_mat=resize_bigger(overall_sum_trees_mat,overall_size);
                  }
                  double BIC_to_add=prev_round_BIC2[h];
                  overall_sum_BIC.push_back(BIC_to_add);
                  overall_sum_preds.insert_cols(overall_sum_preds.n_cols,prev_round_preds2.col(h));
                  if(is_test_data==1) overall_sum_test_preds.insert_cols(overall_sum_test_preds.n_cols,prev_round_test_preds2.col(h));
                }else{
                  NumericMatrix tree_no_child=prev_sum_trees[h];
                  NumericVector resids_no_child= prev_sum_tree_resids[h];
                  NumericMatrix treemat_no_child=prev_sum_trees_mat[h];
                  overall_sum_trees[overall_count]=tree_no_child;
                  overall_sum_tree_resids[overall_count]=resids_no_child;
                  overall_sum_trees_mat[overall_count]=treemat_no_child;
                  overall_count++;
                  
                  if(overall_count==(overall_size-1)){
                    overall_size=overall_size*2;
                    overall_sum_trees=resize_bigger(overall_sum_trees,overall_size);
                    overall_sum_tree_resids=resize_bigger(overall_sum_tree_resids,overall_size);
                    overall_sum_trees_mat=resize_bigger(overall_sum_trees_mat,overall_size);
                  }
                  
                  double BIC_to_add=prev_round_BIC2[h];
                  overall_sum_BIC.push_back(BIC_to_add);
                  
                  overall_sum_preds.insert_cols(overall_sum_preds.n_cols,prev_round_preds2.col(h));
                  if(is_test_data==1) overall_sum_test_preds.insert_cols(overall_sum_test_preds.n_cols,prev_round_test_preds2.col(h));
                }
              }
            }
          }
        }
      }
      // Rcout << "Get to end of checking for models from previous rounds in outer loop number " << j << " . \n";
      
    }
    
    
    //Rcout << "Line 2828 curr_round_lik.size()=" << curr_round_lik.size() << ".\n";
    //Rcout << "Line 2828 overall_sum_trees.size()=" << overall_sum_trees.size() << ".\n";
    //Rcout << "Line 2828 overall_count=" << overall_count << ".\n";
    
    //Rcout << "temp_preds = " << temp_preds << ".\n";
    // Rcout << "Line 9912. j=  " << j << " . \n";
    
    prev_round_preds=temp_preds;
    if(is_test_data==1) prev_round_test_preds=temp_test_preds;
    prev_round_BIC=temp_BIC;
    prev_round_BIC2=temp_BIC;
    prev_round_preds2=Rcpp::as<arma::mat>(temp_preds);
    if(is_test_data==1) prev_round_test_preds2=Rcpp::as<arma::mat>(temp_test_preds);
    resids=temp_resids;
    parent=temp_parent;    
    overall_sum_trees=resize(overall_sum_trees,overall_count);
    overall_sum_tree_resids=resize(overall_sum_tree_resids,overall_count);
    overall_sum_trees_mat=resize(overall_sum_trees_mat,overall_count);     
    
    // Rcout << "Line 9926. j=  " << j << " . \n";
    
    if(j==0){  
      prev_sum_trees=temp_sum_trees;
      prev_sum_tree_resids=temp_sum_tree_resids;
      //NumericMatrix test=prev_sum_trees[0];
      prev_sum_trees_mat=temp_sum_trees_mat; 
      //overall_sum_trees=resize(temp_sum_trees,temp_sum_trees.size());
      overall_sum_trees=temp_sum_trees;
      //overall_sum_tree_resids=resize(temp_sum_tree_resids,temp_sum_tree_resids.size());
      overall_sum_tree_resids=temp_sum_tree_resids;
      overall_sum_trees_mat=temp_sum_trees_mat; 
      //overall_sum_trees_mat=resize(temp_sum_trees_mat,temp_sum_trees.size());
      overall_sum_BIC=temp_BIC;
      overall_sum_preds= Rcpp::as<arma::mat>(temp_preds);
      if(is_test_data==1) overall_sum_test_preds= Rcpp::as<arma::mat>(temp_test_preds);
    }else{
      
      //Rcout << "Line 2856. Length of tree_table = " << tree_table.size() << " . \n";
      //Rcout << "Line 2857. Length of prev_sum_trees = " << prev_sum_trees.size() << " . \n";
      //Rcout << "Line 2858. Length of overall_sum_trees = " << overall_sum_trees.size() << " . \n";
      
      
      prev_sum_trees=overall_sum_trees;
      prev_sum_tree_resids=overall_sum_tree_resids;
      prev_sum_trees_mat=overall_sum_trees_mat;
      prev_round_BIC2=overall_sum_BIC;
      prev_round_preds2=overall_sum_preds;
      if(is_test_data==1) prev_round_test_preds2=overall_sum_test_preds;
    }
    // overall_overall_sum_trees[oo_count]=overall_sum_trees;
    // overall_overall_sum_tree_resids[oo_count]=overall_sum_tree_resids;
    // overall_overall_sum_trees_mat[oo_count]=overall_sum_trees_mat;
    // overall_overall_sum_BIC[oo_count]=overall_sum_BIC;
    // oo_count ++;
    // if(oo_count==(oo_size-1)){
    //   oo_size=oo_size*2;
    //   overall_overall_sum_trees=resize_bigger(overall_overall_sum_trees,oo_size);
    //   overall_overall_sum_tree_resids=resize_bigger(overall_overall_sum_tree_resids,oo_size);
    //   overall_overall_sum_trees_mat=resize_bigger(overall_overall_sum_trees_mat,oo_size);
    //   overall_overall_sum_BIC=resize_bigger(overall_overall_sum_BIC,oo_size);
    // } 
    if(j==num_rounds-1){
      overall_overall_sum_trees=overall_sum_trees;
      overall_overall_sum_tree_resids=overall_sum_tree_resids;
      overall_overall_sum_trees_mat=overall_sum_trees_mat;
      overall_overall_sum_BIC=overall_sum_BIC;
      overall_overall_sum_preds=overall_sum_preds;
      if(is_test_data==1) overall_overall_sum_test_preds=overall_sum_test_preds;
    }
    
    //overall_trees[j]=curr_round_trees;
    //overall_mat.push_back(curr_round_mat);
    overall_lik.push_back(curr_round_lik);
    prev_par=seq_len(overall_sum_trees.size())-1;
    // Rcout << "Get to end of outer loop number " << j << " . \n";
    
  }
  // Rcout << "Get past outer loop \n";
  
  if(first_round_break==1){
    throw std::range_error("BART-BMA didnt find any suitable model for the data. Maybe limit for Occam's window is too small. Alternatively, try using more observations or change parameter values.");
  }
  
  //overall_overall_sum_trees=resize(overall_overall_sum_trees,oo_count);
  //overall_overall_sum_tree_resids=resize(overall_overall_sum_tree_resids,oo_count);
  //overall_overall_sum_trees_mat=resize(overall_overall_sum_trees_mat,oo_count);
  //overall_overall_sum_BIC=resize(overall_overall_sum_BIC,oo_count);
  
  
  
  // Rcout << "Get to defining end_BIC \n";
  //NumericVector end_BIC=overall_overall_sum_BIC[overall_overall_sum_BIC.size()-1] ;
  NumericVector end_BIC=overall_overall_sum_BIC ;
  // 
  // 
  // //Rcout << "Get to defining past defining end_BIC \n";
  // NumericMatrix overallpreds(n,end_BIC.size());
  // NumericMatrix overall_test_preds(test_data.nrow(),end_BIC.size());
  // NumericVector post_weights(end_BIC.size());
  // //Rcout << "Get to loop that defines temp_pred M1 \n";
  // for(int k=0;k<end_BIC.size();k++){
  //   NumericMatrix oosp=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_preds));
  //   NumericVector temp_preds=oosp(_,k);
  //   NumericVector temp_test_preds;
  //   if(is_test_data==1){
  //     NumericMatrix oostp=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_test_preds));
  //     temp_test_preds=oostp(_,k); 
  //   }
  //   //Rcout << "Get to defining orig_temp_pred \n";
  //   NumericVector orig_temp_preds=get_original(min(y),max(y),-0.5,0.5,temp_preds) ;
  //   NumericVector BICi=-0.5*end_BIC;
  //   double max_BIC=max(BICi);
  //   double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
  //   post_weights[k]=weight;
  //   overallpreds(_,k) = temp_preds*weight;
  //   if(is_test_data==1){
  //     overall_test_preds(_,k) = temp_test_preds*weight;
  //   }
  // }
  // //Rcout << "Get to defining M1 \n";
  // arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);
  // predicted_values=sum(M1,1);
  // arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  // 
  //if(is_test_data==1) predicted_test_values=sum(M2,1);
  if(overall_lik.size()==0){
    throw std::range_error("BART-BMA didnt find any suitable model for the data. Maybe limit for Occam's window is too small.");
  }else{
    
    //NumericVector orig_preds=get_original(min(y),max(y),-0.5,0.5,wrap(predicted_values)) ;
    //NumericVector orig_test_preds;
    //if(is_test_data==1){
    //  orig_test_preds=get_original(min(y),max(y),-0.5,0.5,wrap(predicted_test_values)) ;
    //}
    
    NumericVector orig_preds=preds_bbma_lin_alg_insamp(overall_overall_sum_trees,
                                                       overall_overall_sum_trees_mat,
                                                       y,
                                                       end_BIC,
                                                       100,
                                                       10,
                                                       y.size() ,
                                                       a,
                                                       0,
                                                       0,
                                                       nu,
                                                       lambda//,List resids
    );
    
    NumericVector orig_test_preds;
    if(is_test_data==1){
      orig_test_preds=preds_bbma_lin_alg_outsamp(overall_overall_sum_trees,
                                                 overall_overall_sum_trees_mat,
                                                 y,
                                                 end_BIC,
                                                 100,
                                                 10,
                                                 y.size() ,
                                                 test_data.nrow(),
                                                 a,
                                                 0, 0, nu,
                                                 lambda,//List resids,
                                                 test_data);
    }
    
    //NumericVector minmax(2);
    //minmax[0]=min(y);
    //minmax[1]=max(y);
    // Rcout << "Get to defining ret \n";
    if(is_test_data==1){
      List ret(6);
      ret[0] = orig_preds;
      ret[1] = overall_overall_sum_trees;
      ret[2] =overall_overall_sum_trees_mat;
      ret[3] = end_BIC;
      ret[4] = orig_test_preds;
      ret[5] =overall_overall_sum_tree_resids;
      return(ret);
    }else{
      List ret(5);
      ret[0] = orig_preds;
      ret[1] = overall_overall_sum_trees;
      ret[2] = overall_overall_sum_trees_mat;
      ret[3] = end_BIC;
      ret[4] =overall_overall_sum_tree_resids;
      return(ret);
    }
  }
}


//###########################################################################################################################//

// [[Rcpp::export]]
Rcpp::NumericVector Quantile(Rcpp::NumericVector x, Rcpp::NumericVector probs) {
  // implementation of type 7
  const size_t n=x.size(), np=probs.size();
  if (n==0) return x;
  if (np==0) return probs;
  Rcpp::NumericVector index = (n-1.)*probs, y=x.sort(), x_hi(np), qs(np);
  Rcpp::NumericVector lo = Rcpp::floor(index), hi = Rcpp::ceiling(index);
  Rcpp::NumericVector probs_names=Rcpp::round(probs*100., 6);
  Rcpp::CharacterVector qs_names(np);
  std::string _name, zero=".000000";
  
  for (size_t i=0; i<np; ++i) {
    _name = std::to_string(static_cast<double>(probs_names[i]));
    if(_name=="100.000000") _name=std::string("100");
    else if(_name.substr(1)==zero) _name=_name.substr(0, 1);
    else if(_name.substr(2)==zero) _name=_name.substr(0, 2);
    else if(_name.substr(3)==zero.substr(2)) _name=_name.substr(0, 3);
    else if(_name.substr(4)==zero.substr(3)) _name=_name.substr(0, 4);
    else if(_name.substr(5)==zero.substr(4)) _name=_name.substr(0, 5);
    else if(_name.substr(6)==zero.substr(5)) _name=_name.substr(0, 6);
    else if(_name.substr(7)==zero.substr(6)) _name=_name.substr(0, 7);
    else if(_name.substr(4)==zero.substr(2)) _name=_name.substr(0, 4);
    else if(_name.substr(5)==zero.substr(3)) _name=_name.substr(0, 5);
    else if(_name.substr(6)==zero.substr(4)) _name=_name.substr(0, 6);
    else if(_name.substr(7)==zero.substr(5)) _name=_name.substr(0, 7);
    else if(_name.substr(8)==zero.substr(6)) _name=_name.substr(0, 8);
    /*  // why does this not work?
    else {
    for(size_t j=3; j<8; ++j) {
    if(_name.substr(j)==zero.substr(j-1)) _name=_name.substr(0, j);
    else if(_name.substr(j+1)==zero.substr(j-1)) _name=_name.substr(0, j+1);
    }
    }
    */
    qs_names[i] = _name + std::string("%");
    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];
    if ((index[i]>lo[i]) && (x_hi[i] != qs[i])) {
      double h;
      h = index[i]-lo[i];
      qs[i] = (1.-h)*qs[i] + h*x_hi[i];
    }
  }
  
  qs.names()=qs_names;
  return qs;
  }


//###########################################################################################################################//
#include <algorithm>
#include <cmath>
#include <vector>

template<typename T>
static inline double Lerp(T v0, T v1, T t)
{
  return (1 - t)*v0 + t*v1;
}

//###########################################################################################################################//
template<typename T>
static inline std::vector<T> Quantile2(const std::vector<T>& inData, const std::vector<T>& probs)
{
  // function from https://stackoverflow.com/questions/11964552/finding-quartiles
  if (inData.empty())
  {
    return std::vector<T>();
  }
  
  if (1 == inData.size())
  {
    return std::vector<T>(1, inData[0]);
  }
  
  std::vector<T> data = inData;
  std::sort(data.begin(), data.end());
  std::vector<T> quantiles;
  
  for (size_t i = 0; i < probs.size(); ++i)
  {
    T poi = Lerp<T>(-0.5, data.size() - 0.5, probs[i]);
    
    size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
    size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));
    
    T datLeft = data.at(left);
    T datRight = data.at(right);
    
    T quantile = Lerp<T>(datLeft, datRight, poi - left);
    
    quantiles.push_back(quantile);
  }
  
  return quantiles;
}

//###########################################################################################################################//
#include <boost/math/distributions/students_t.hpp>

// [[Rcpp::depends(BH)]]
// [[Rcpp::export]]
double mixt_eval_cdf(double x_val, double d_o_f, std::vector<double> mean_vec, std::vector<double> var_vec, std::vector<double> weights_vec, double quant_val) {
  
  boost::math::students_t dist(d_o_f);
  
  double ret_val=0;
  for(unsigned int i=0; i < weights_vec.size();i++){
    if(var_vec[i]>0){
      double tempx = (x_val-mean_vec[i])/sqrt(var_vec[i]);
      ret_val += weights_vec[i]*boost::math::cdf(dist,tempx);
    }else{//in some cases (for ITEs) there is zero variance, and can't divide by zero
      //Rcout << " \n \n VARIANCE = " << var_vec[i] << ".\n \n" ;
      if(x_val>=mean_vec[i]){ //if no variance, cdf is zero below mean, and one above mean
        ret_val += weights_vec[i];
      } // no else statement because add zero if below mean (when there is zero variance)
    }
  }
  
  return (ret_val-quant_val) ;  // approximation
}

//###########################################################################################################################//
// [[Rcpp::export]]
double rootmixt(double d_o_f, double a, double b,
                std::vector<double> mean_vec,
                std::vector<double> var_vec,
                std::vector<double> weights_vec, double quant_val, double root_alg_precision){
  
  static const double EPS = root_alg_precision;//1e-15; // 1×10^(-15)
  
  double fa = mixt_eval_cdf(a, d_o_f, mean_vec, var_vec, weights_vec,quant_val), fb = mixt_eval_cdf(b, d_o_f, mean_vec, var_vec, weights_vec,quant_val);
  
  // if either f(a) or f(b) are the root, return that
  // nothing else to do
  if (fa == 0) return a;
  if (fb == 0) return b;
  
  //Rcout << "quant_val = " << quant_val << ".\n";
  
  //Rcout << "fa = " << fa << ".\n";
  //Rcout << "fb = " << fb << ".\n";
  
  // this method only works if the signs of f(a) and f(b)
  
  if(fa * fb >= 0){throw std::range_error("fa * fb >= 0.This method only works if the signs of f(a) and f(b) are different.");}
  do {
    // calculate fun at the midpoint of a,b
    // if that's the root, we're done
    // prefer:
    double midpt = (a + b) / 2;
    double fmid = mixt_eval_cdf(midpt, d_o_f, mean_vec, var_vec, weights_vec,quant_val);
    
    if (fmid == 0) return midpt;
    
    // adjust our bounds to either [a,midpt] or [midpt,b]
    // based on where fmid ends up being. I'm pretty
    // sure the code in the question is wrong, so I fixed it
    if (fa * fmid < 0) { // fmid, not f1
      fb = fmid; 
      b = midpt;
    }
    else {
      fa = fmid; 
      a = midpt; 
    }
  } while (b-a > EPS); // only loop while
  // a and b are sufficiently far
  // apart
  //Rcout << "a = " << a << ".\n";
  //Rcout << "b = " << b << ".\n";
  return (a + b) / 2;  // approximation
}


//###########################################################################################################################//

std::vector<double> mixt_find_boundsQ(double d_o_f, std::vector<double> mean_vec, std::vector<double> var_vec, double quant_val) {
  //boost::math::students_t dist1(d_o_f);
  
  std::vector<double> tempbounds(mean_vec.size());
  
  for(unsigned int i=0; i< mean_vec.size();i++){
    //tempbounds[i]= mean_vec[i]+sqrt(var_vec[i])*boost::math::quantile(dist1,quant_val);
    tempbounds[i]= mean_vec[i]+sqrt(var_vec[i])*quant_val;
  }
  
  std::vector<double> ret(2);
  ret[0]= *std::min_element(tempbounds.begin(), tempbounds.end()); 
  ret[1]= *std::max_element(tempbounds.begin(), tempbounds.end()); 
  
  return(ret) ; 
}
//###########################################################################################################################//

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_exact_outsamp(List overall_sum_trees,
                             List overall_sum_mat,
                             NumericVector y,
                             NumericVector BIC_weights,//double min_possible,double max_possible,
                             int num_obs,int num_test_obs,
                             double a,double sigma,double mu_mu,double nu,
                             double lambda,//List resids,
                             NumericMatrix test_data, double lower_prob, double upper_prob, int num_cores,
                             double root_alg_precision){
  
  
  List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat t_vars_arma(num_test_obs,BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  //int num_models= BIC_weights.size();
  
  
  
  //arma::mat draws_for_preds(0,num_test_obs);
  //arma::mat draws_for_preds(num_iter,num_test_obs);
  
  //  {
  //#pragma omp for schedule(dynamic,1) 
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  arma::mat yty=y_arma.t()*y_arma;
  
  
  //#pragma omp parallel num_threads(1)
  //#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    
    //Rcout << "Line 4413. i= "<< i << ".\n";
    
    arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    double b=Wmat.n_cols;
    
    //get exponent
    //double expon=(n+nu)/2;
    //get y^Tpsi^{-1}y
    // arma::mat psi_inv=psi.i();
    
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*Wmat;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=Wmat.t()*Wmat;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    //arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    //arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    arma::mat third_term=Wmat.t()*y_arma;
    
    
    //arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    
    
    arma::mat w_tilde_M_inv =  W_tilde*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    //arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    //arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    arma::mat I_test(num_test_obs,num_test_obs);
    I_test=I_test.eye();
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    
    
    //arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    //arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights[i];
    preds_all_models_arma.col(i)=preds_temp_arma;
    
    // Rcout << "Line 4468. i= "<< i << ".\n";
    
    
    t_vars_arma.col(i)=covar_t.diag();
    
  }
  
  //}
  //#pragma omp barrier  
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  NumericMatrix output(3, num_test_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  
  //std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  
  
  
  typedef std::vector<double> stdvec;
  std::vector<double> weights_vec= as<stdvec>(post_weights);
  
  boost::math::students_t dist2(nu+num_obs);
  double lq_tstandard= boost::math::quantile(dist2,lower_prob);
  double med_tstandard= boost::math::quantile(dist2,0.5); //This is just 0 ??
  double uq_tstandard= boost::math::quantile(dist2,upper_prob);
  
  if(weights_vec.size()==1){
    
    for(int i=0;i<num_test_obs;i++){
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      boost::math::students_t dist2(nu+num_obs);
      
      
      output(0,i)= tempmeans[0]+sqrt(tempvars[0])*lq_tstandard;
      output(1,i)= tempmeans[0]+sqrt(tempvars[0])*med_tstandard;
      output(2,i)= tempmeans[0]+sqrt(tempvars[0])*uq_tstandard;
      
      
    }
  }else{
    
    for(int i=0;i<num_test_obs;i++){
      //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
      
      //The two lines below can be removed by replacing the standard vectors with armadillo vectors in teh functions below
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      
      std::vector<double> bounds_lQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, lq_tstandard);
      
      output(0,i)=rootmixt(nu+num_obs,  
             bounds_lQ[0]-0.0001,  
             bounds_lQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, lower_prob,root_alg_precision);
      
      std::vector<double> bounds_med = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, med_tstandard);
      
      output(1,i)=rootmixt(nu+num_obs,  
             bounds_med[0]-0.0001,  
             bounds_med[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, 0.5,root_alg_precision);
      
      std::vector<double> bounds_uQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, uq_tstandard);
      
      
      output(2,i)=rootmixt(nu+num_obs,  
             bounds_uQ[0]-0.0001,  
             bounds_uQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, upper_prob,root_alg_precision);
      
      
    }  
  }
  
  
  
  for(int i=0;i<output.ncol();i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output(_,i)=get_original(min(y),max(y),-0.5,0.5, output(_,i));
  }  
  
  List ret(2);
  ret[0]= output;
  ret[1]= wrap(get_original_arma(min(y),max(y),-0.5,0.5,predicted_values));
  
  
  return(ret);
  
}
//###########################################################################################################################//

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_exact_outsamp_par(List overall_sum_trees,
                                 List overall_sum_mat,
                                 NumericVector y,
                                 NumericVector BIC_weights,//double min_possible,double max_possible,
                                 int num_obs,int num_test_obs,
                                 double a,double sigma,double mu_mu,double nu,
                                 double lambda,//List resids,
                                 NumericMatrix test_data, double lower_prob, double upper_prob, int num_cores,
                                 double root_alg_precision){
  
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,test_data);
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat t_vars_arma(num_test_obs,BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  arma::vec post_weights_arma = as<arma::vec>(post_weights);
  
  //int num_models= BIC_weights.size();
  
  
  
  
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  arma::mat yty=y_arma.t()*y_arma;
  
  
  arma::mat I_test(num_test_obs,num_test_obs);
  I_test=I_test.eye();
  
  //create field (armadillo list) of models
  //each model is a field (armadillo list) of trees represented by matrices
  arma::field<arma::field<arma::mat>> modelsF(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    List temp_tree_list = overall_sum_trees[i];
    //Rcout << "Line 5661.\n";
    
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      //Rcout << "Line 5663.\n";
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    //Rcout << "Line 5669.\n";
    
    modelsF(i)=treesF;
  }
  
  
  arma::field<arma::field<arma::mat>> matsF(overall_sum_mat.size());
  for(int i=0;i<overall_sum_mat.size();i++){
    List temp_tree_list = overall_sum_mat[i];
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    matsF(i)=treesF;
  }
  
  
  
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    
    
    
    
    arma::field<arma::mat> tree_list = modelsF(i);
    
    arma::mat W(num_obs,0);
    int upsilon=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat curr_tree=(modelsF(i))(j);
      arma::mat curr_obs_nodes=(matsF(i))(j);
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      int b_j=term_nodes.n_elem;
      //begin J function
      
      //will make J as we go in BART-BMA no need to create it again here....
      // arma::mat Jmat=J(curr_obs_nodes,tree_term_nodes);
      arma::mat tree_matrix_temp = (matsF(i))(j);
      arma::mat Jmat(tree_matrix_temp.n_rows, b_j);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(int q=0;q<b_j;q++){
        //double tn=term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        
        //begin find_term_obs
        
        //arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
        //for reference arma_tree_mat == matsF(i)(j) == tree_matrix_temp
        
        arma::uvec term_obs;
        
        for(unsigned int r=0;r<tree_matrix_temp.n_cols;r++){
          //arma::vec colmat=arma_tree_mat.col(r);
          arma::vec colmat=tree_matrix_temp.col(r);
          term_obs=arma::find(colmat==term_nodes[q]);
          if(term_obs.size()>0){
            break;
          }
        }
        
        //end find_term_obs
        
        
        //assign term_obs to the correct index of J
        //NumericVector term_obs2=as<NumericVector>(wrap(term_obs));
        //NumericVector obs_col(obs_to_nodes_temp.nrow());
        arma::vec obs_col= arma::zeros<arma::vec>(tree_matrix_temp.n_rows);
        //Rcout << "Line 5747.\n";
        obs_col.elem(term_obs)= arma::ones<arma::vec>(term_obs.n_elem);
        //Rcout << "Line 5749.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(i)= colmat;
      }
      
      
      
      
      W.insert_cols(upsilon,Jmat);
      upsilon+=b_j;
    }
    
    
    //Rcout << "Line 5759.\n";
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    ////////////////////////////////////////
    
    arma::mat W_tilde(num_test_obs,0);
    int upsilon2=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_test_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_test_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde.insert_cols(upsilon2,Jmat);
      upsilon2+=b_j;
    }
    
    ////////////////////////////////////////
    //Rcout << "Line 5819.\n";
    
    
    double b=W.n_cols;
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*W;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=W.t()*W;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    //arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    //arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    arma::mat third_term=W.t()*y_arma;
    
    
    //arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    
    
    arma::mat w_tilde_M_inv =  W_tilde*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    //arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    //arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    
    
    //arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    //arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights_arma[i];
    
    preds_all_models_arma.col(i)=preds_temp_arma;
    
    
    
    t_vars_arma.col(i)=covar_t.diag();
    
  }
  
  //}
#pragma omp barrier  
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  arma::mat output(3, num_test_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  
  //std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  
  
  
  typedef std::vector<double> stdvec;
  std::vector<double> weights_vec= as<stdvec>(post_weights);
  
  boost::math::students_t dist2(nu+num_obs);
  double lq_tstandard= boost::math::quantile(dist2,lower_prob);
  double med_tstandard= boost::math::quantile(dist2,0.5); //This is just 0 ??
  double uq_tstandard= boost::math::quantile(dist2,upper_prob);
  
  
  if(weights_vec.size()==1){
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_test_obs;i++){
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      //boost::math::students_t dist2(nu+num_obs);
      
      
      output(0,i)= tempmeans[0]+sqrt(tempvars[0])*lq_tstandard;
      output(1,i)= tempmeans[0]+sqrt(tempvars[0])*med_tstandard;
      output(2,i)= tempmeans[0]+sqrt(tempvars[0])*uq_tstandard;
      
      
    }
#pragma omp barrier  
  }else{
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_test_obs;i++){
      //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      
      std::vector<double> bounds_lQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, lq_tstandard);
      
      output(0,i)=rootmixt(nu+num_obs,  
             bounds_lQ[0]-0.0001,  
             bounds_lQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, lower_prob,root_alg_precision);
      
      
      std::vector<double> bounds_med = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, med_tstandard);
      
      output(1,i)=rootmixt(nu+num_obs,  
             bounds_med[0]-0.0001,  
             bounds_med[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, 0.5,root_alg_precision);
      
      std::vector<double> bounds_uQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, uq_tstandard);
      
      output(2,i)=rootmixt(nu+num_obs,  
             bounds_uQ[0]-0.0001,  
             bounds_uQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, upper_prob,root_alg_precision);
      
      
    }  
#pragma omp barrier  
  }
  
  arma::mat output_rescaled(output.n_rows, output.n_cols);
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(unsigned int i=0;i<output.n_cols;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output_rescaled.col(i)=get_original_arma(min(y),max(y),-0.5,0.5, output.col(i));
    
    
  }  
#pragma omp barrier 
  
  List ret(2);
  ret[0]= wrap(output_rescaled);
  ret[1]= wrap(get_original_arma(min(y),max(y),-0.5,0.5,predicted_values));
  
  
  return(ret);
  
}

//###########################################################################################################################//

#include <RcppArmadilloExtensions/rmultinom.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_lin_alg_outsamp(List overall_sum_trees,
                               List overall_sum_mat,
                               NumericVector y,
                               NumericVector BIC_weights,
                               int num_iter,int burnin,int num_obs,int num_test_obs,
                               double a,double sigma,double mu_mu,double nu,
                               double lambda,//List resids,
                               NumericMatrix test_data, double lower_prob, double upper_prob){
  
  
  List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  //int num_models= BIC_weights.size();
  
  IntegerVector num_its_to_sample = RcppArmadillo::rmultinom(num_iter,post_weights);
  //IntegerVector num_its_sum = cumsum(num_its_to_sample);
  
  
  arma::mat draws_for_preds(0,num_test_obs);
  //arma::mat draws_for_preds(num_iter,num_test_obs);
  
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  //get exponent
  //double expon=(n+nu)/2;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  
  
  arma::mat yty=y_arma.t()*y_arma;
  
  arma::mat I_test(num_test_obs,num_test_obs);
  I_test=I_test.eye();
  
  for(int i=0;i<overall_sum_trees.size();i++){
    
    // Rcout << "Line 4413. i= "<< i << ".\n";
    
    arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    double b=Wmat.n_cols;
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*Wmat;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=Wmat.t()*Wmat;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    //arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    //arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    arma::mat third_term=Wmat.t()*y_arma;
    
    
    //arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    
    
    arma::mat w_tilde_M_inv =  W_tilde*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    //arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    //arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    
    
    //arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    //arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights[i];
    preds_all_models_arma.col(i)=preds_temp_arma;
    
    // Rcout << "Line 4468. i= "<< i << ".\n";
    
    
    //int num_its_to_sample = round(weight*(num_iter));
    
    int num_its_temp = num_its_to_sample[i];
    
    // arma::uword m = num_test_obs;
    // arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    // U = sqrt(nu+num_obs / U);
    // arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    // Rcout << "Line 4483. i= "<< i << ".\n";
    
    
    
    
    //arma::uword m = num_test_obs;
    arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    //arma::vec U = Rcpp::rchisq(num_its_temp,nu+num_obs);
    U = sqrt((nu+num_obs) / U);
    //arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    arma::mat Y = arma::randn(num_its_temp, num_test_obs);
    arma::mat temp_LY = Y*arma::chol(covar_t);
    temp_LY.each_col() %= U;
    arma::mat temp_draws = arma::repmat(preds_temp_arma, 1, num_its_temp).t() +  temp_LY;
    
    
    
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    
    // Rcout << "number of rows of preds_temp_arma = " << preds_temp_arma.n_rows << ".\n";
    // Rcout << "number of rows of temp_LY = " << temp_LY.n_rows << ".\n";
    // 
    // Rcout << "number of rows of draws_for_preds = " << draws_for_preds.n_rows << ".\n";
    // Rcout << "number of rows of temp_draws = " << temp_draws.n_rows << ".\n";
    //
    
    
    draws_for_preds = join_cols(draws_for_preds,temp_draws);
    // if(i==0){
    //   draws_for_preds.rows(0,num_its_sum[i]-1)=temp_draws;
    // }else{
    //   draws_for_preds.rows(num_its_sum[i-1],num_its_sum[i]-1)=temp_draws;
    // }
    
  }
  
  
  //arma::colvec predicted_values;
  // Rcout << "Line 4491";
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  NumericMatrix output(3, num_test_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  
  
  
  // Rcout << "Line 4492";
  //Rcout << "probs_for_quantiles = " << probs_for_quantiles << ".\n";
  typedef std::vector<double> stdvec;
  
  for(int i=0;i<num_test_obs;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    std::vector<double> tempcol= arma::conv_to<stdvec>::from(draws_for_preds.col(i));
    std::vector<double> tempquant= Quantile2(tempcol, probs_for_quantiles);
    NumericVector tempforoutput = wrap(tempquant);
    output(_,i)= tempforoutput;
    
  }  
  
  for(int i=0;i<output.ncol();i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output(_,i)=get_original(min(y),max(y),-0.5,0.5, output(_,i));
  }  
  
  
  List ret(2);
  ret[0]= output;
  ret[1]= wrap(get_original_arma(min(y),max(y),-0.5,0.5,predicted_values));
  
  
  return(ret);
  
}

//###########################################################################################################################//


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_lin_alg_insamp(List overall_sum_trees,
                              List overall_sum_mat,
                              NumericVector y,
                              NumericVector BIC_weights,
                              int num_iter,int burnin,int num_obs,
                              double a,double sigma,double mu_mu,double nu,
                              double lambda,//List resids,
                              double lower_prob, double upper_prob){
  
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_obs,BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  //int num_models= BIC_weights.size();
  
  IntegerVector num_its_to_sample = RcppArmadillo::rmultinom(num_iter,post_weights);
  
  
  
  arma::mat draws_for_preds(0,num_obs);
  
  for(int i=0;i<overall_sum_trees.size();i++){
    
    
    arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    double b=Wmat.n_cols;
    NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
    arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
    arma::mat y_arma(num_obs,1);
    y_arma.col(0)=yvec;
    //get exponent
    //double expon=(n+nu)/2;
    //get y^Tpsi^{-1}y
    // arma::mat psi_inv=psi.i();
    
    
    arma::mat yty=y_arma.t()*y_arma;
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*Wmat;
    
    
    //get t(J)inv(psi)J
    arma::mat WtW=Wmat.t()*Wmat;
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);
    //get t(J)inv(psi)y
    arma::mat third_term=Wmat.t()*y_arma;
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    arma::vec preds_temp_arma= Wmat*sec_term_inv*third_term;
    
    //arma::mat I_test(num_test_obs,num_test_obs);
    //I_test=I_test.eye();
    
    arma::mat covar_t=as_scalar((nu*lambda+yty-mvm)/(nu+num_obs))*(Wmat*sec_term_inv*Wmat.t());
    
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    //weighted_preds_all_models_arma.col(i)=preds_temp_arma*weight;
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights[i];
    
    preds_all_models_arma.col(i)=preds_temp_arma;
    
    
    
    //int num_its_to_sample = round(weight*(num_iter));
    int num_its_temp = num_its_to_sample[i];
    
    arma::uword m = num_obs;
    arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    U = sqrt(nu+num_obs / U);
    arma::mat Y = mvnrnd( arma::vec(m, arma::fill::zeros), covar_t,num_its_temp);
    //arma::mat temp_draws(m, num_its_temp);
    //for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //  temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    //}
    
    
    arma::rowvec Ut=U.t();
    Y.each_row() %= Ut;
    arma::mat temp_draws = arma::repmat(preds_temp_arma, 1, num_its_temp) +  Y;
    
    
    
    //draws_for_preds = join_cols(draws_for_preds,temp_draws.t());
    //if(i==0){
    //  draws_for_preds.rows(0,num_its_sum_arma[i]-1)=temp_draws.t();
    //}else{
    //  draws_for_preds.rows(num_its_sum_arma[i-1],num_its_sum_arma[i]-1)=temp_draws.t();
    //}
    
    draws_for_preds = join_cols(draws_for_preds,temp_draws.t());
    
  }
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  NumericMatrix output(3, num_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  
  
  
  typedef std::vector<double> stdvec;
  
  for(int i=0;i<num_obs;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    std::vector<double> tempcol= arma::conv_to<stdvec>::from(draws_for_preds.col(i));
    std::vector<double> tempquant= Quantile2(tempcol, probs_for_quantiles);
    NumericVector tempforoutput = wrap(tempquant);
    output(_,i)= tempforoutput;
    
  } 
  
  for(int i=0;i<output.ncol();i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output(_,i)=get_original(min(y),max(y),-0.5,0.5, output(_,i));
  }  
  
  
  List ret(2);
  ret[0]= output;
  ret[1]= wrap(get_original_arma(min(y),max(y),-0.5,0.5,predicted_values));
  
  
  return(ret);
}

//###########################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_chol_attempt_outsamp(List overall_sum_trees,
                                    List overall_sum_mat,
                                    NumericVector y,
                                    NumericVector BIC_weights,
                                    int num_iter,int burnin,int num_obs,int num_test_obs,
                                    double a,double sigma,double mu_mu,double nu,
                                    double lambda,//List resids,
                                    NumericMatrix test_data, double lower_prob, double upper_prob){
  
  
  List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  //int num_models= BIC_weights.size();
  
  IntegerVector num_its_to_sample = RcppArmadillo::rmultinom(num_iter,post_weights);
  
  
  arma::mat draws_for_preds(0,num_test_obs);
  
  for(int i=0;i<overall_sum_trees.size();i++){
    
    // Rcout << "Line 4413. i= "<< i << ".\n";
    
    arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    double b=Wmat.n_cols;
    NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
    arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
    arma::mat y_arma(num_obs,1);
    y_arma.col(0)=yvec;
    //get exponent
    //double expon=(n+nu)/2;
    //get y^Tpsi^{-1}y
    // arma::mat psi_inv=psi.i();
    
    
    arma::mat yty=y_arma.t()*y_arma;
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*Wmat;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=Wmat.t()*Wmat;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    //arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    // arma::mat third_term=Wmat.t()*y_arma;
    
    
    arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    //arma::mat mvm= ytW*sec_term_inv*third_term;
    
    
    
    //arma::mat w_tilde_M_inv =  W_tilde*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    //arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    arma::mat I_test(num_test_obs,num_test_obs);
    I_test=I_test.eye();
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    
    
    arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights[i];
    preds_all_models_arma.col(i)=preds_temp_arma;
    
    // Rcout << "Line 4468. i= "<< i << ".\n";
    
    
    //int num_its_to_sample = round(weight*(num_iter));
    
    int num_its_temp = num_its_to_sample[i];
    
    // arma::uword m = num_test_obs;
    // arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    // U = sqrt(nu+num_obs / U);
    // arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    // Rcout << "Line 4483. i= "<< i << ".\n";
    
    
    
    
    //arma::uword m = num_test_obs;
    arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    //arma::vec U = Rcpp::rchisq(num_its_temp,nu+num_obs);
    U = sqrt((nu+num_obs) / U);
    //arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    arma::mat Y = arma::randn(num_its_temp, num_test_obs);
    arma::mat temp_LY = Y*arma::chol(covar_t);
    temp_LY.each_col() %= U;
    arma::mat temp_draws = arma::repmat(preds_temp_arma, 1, num_its_temp).t() +  temp_LY;
    
    
    
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    
    // Rcout << "number of rows of preds_temp_arma = " << preds_temp_arma.n_rows << ".\n";
    // Rcout << "number of rows of temp_LY = " << temp_LY.n_rows << ".\n";
    // 
    // Rcout << "number of rows of draws_for_preds = " << draws_for_preds.n_rows << ".\n";
    // Rcout << "number of rows of temp_draws = " << temp_draws.n_rows << ".\n";
    // 
    draws_for_preds = join_cols(draws_for_preds,temp_draws);
    
  }
  
  //arma::colvec predicted_values;
  // Rcout << "Line 4491";
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  NumericMatrix output(3, num_test_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  
  
  
  typedef std::vector<double> stdvec;
  
  for(int i=0;i<num_test_obs;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    std::vector<double> tempcol= arma::conv_to<stdvec>::from(draws_for_preds.col(i));
    std::vector<double> tempquant= Quantile2(tempcol, probs_for_quantiles);
    NumericVector tempforoutput = wrap(tempquant);
    output(_,i)= tempforoutput;
    
  }  
  
  for(int i=0;i<output.ncol();i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output(_,i)=get_original(min(y),max(y),-0.5,0.5, output(_,i));
  }  
  
  
  
  List ret(2);
  ret[0]= output;
  ret[1]= wrap(get_original_arma(min(y),max(y),-0.5,0.5,predicted_values));
  
  
  return(ret);
}



//###########################################################################################################################//

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_lin_alg_parallel_outsamp(List overall_sum_trees,
                                        List overall_sum_mat,
                                        NumericVector y,
                                        NumericVector BIC_weights,
                                        int num_iter,int burnin,int num_obs,int num_test_obs,
                                        double a,double sigma,double mu_mu,double nu,
                                        double lambda,//List resids,
                                        NumericMatrix test_data, double lower_prob, double upper_prob, int num_cores){
  
  
  List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  //int num_models= BIC_weights.size();
  
  IntegerVector num_its_to_sample = RcppArmadillo::rmultinom(num_iter,post_weights);
  IntegerVector num_its_sum = cumsum(num_its_to_sample);
  
  
  //arma::mat draws_for_preds(0,num_test_obs);
  arma::mat draws_for_preds(num_iter,num_test_obs);
  
  //  {
  //#pragma omp for schedule(dynamic,1) 
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  arma::mat yty=y_arma.t()*y_arma;
  
  
  //#pragma omp parallel num_threads(1)
  //#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    
    //Rcout << "Line 4413. i= "<< i << ".\n";
    
    arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    double b=Wmat.n_cols;
    
    //get exponent
    //double expon=(n+nu)/2;
    //get y^Tpsi^{-1}y
    // arma::mat psi_inv=psi.i();
    
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*Wmat;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=Wmat.t()*Wmat;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    //arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    //arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    arma::mat third_term=Wmat.t()*y_arma;
    
    
    //arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    
    
    arma::mat w_tilde_M_inv =  W_tilde*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    //arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    //arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    arma::mat I_test(num_test_obs,num_test_obs);
    I_test=I_test.eye();
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    
    
    //arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    //arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights[i];
    preds_all_models_arma.col(i)=preds_temp_arma;
    
    // Rcout << "Line 4468. i= "<< i << ".\n";
    
    
    //int num_its_to_sample = round(weight*(num_iter));
    
    int num_its_temp = num_its_to_sample[i];
    
    // arma::uword m = num_test_obs;
    // arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    // U = sqrt(nu+num_obs / U);
    // arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    // Rcout << "Line 4483. i= "<< i << ".\n";
    
    
    
    
    //arma::uword m = num_test_obs;
    arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    //arma::vec U = Rcpp::rchisq(num_its_temp,nu+num_obs);
    U = sqrt((nu+num_obs) / U);
    //arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    arma::mat Y = arma::randn(num_its_temp, num_test_obs);
    arma::mat temp_LY = Y*arma::chol(covar_t);
    temp_LY.each_col() %= U;
    arma::mat temp_draws = arma::repmat(preds_temp_arma, 1, num_its_temp).t() +  temp_LY;
    
    
    
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    
    // Rcout << "number of rows of preds_temp_arma = " << preds_temp_arma.n_rows << ".\n";
    // Rcout << "number of rows of temp_LY = " << temp_LY.n_rows << ".\n";
    // 
    // Rcout << "number of rows of draws_for_preds = " << draws_for_preds.n_rows << ".\n";
    // Rcout << "number of rows of temp_draws = " << temp_draws.n_rows << ".\n";
    //
    
    
    //draws_for_preds = join_cols(draws_for_preds,temp_draws);
    if(i==0){
      draws_for_preds.rows(0,num_its_sum[i]-1)=temp_draws;
    }else{
      draws_for_preds.rows(num_its_sum[i-1],num_its_sum[i]-1)=temp_draws;
    }
    
  }
  
  //}
  //#pragma omp barrier  
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  NumericMatrix output(3, num_test_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  
  
  
  typedef std::vector<double> stdvec;
  
  for(int i=0;i<num_test_obs;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    std::vector<double> tempcol= arma::conv_to<stdvec>::from(draws_for_preds.col(i));
    std::vector<double> tempquant= Quantile2(tempcol, probs_for_quantiles);
    NumericVector tempforoutput = wrap(tempquant);
    output(_,i)= tempforoutput;
    
  }  
  
  for(int i=0;i<output.ncol();i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output(_,i)=get_original(min(y),max(y),-0.5,0.5, output(_,i));
  }  
  
  List ret(2);
  ret[0]= output;
  ret[1]= wrap(get_original_arma(min(y),max(y),-0.5,0.5,predicted_values));
  
  
  return(ret);
  
}


//###########################################################################################################################//

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_lin_alg_fields_outsamp(List overall_sum_trees,
                                      List overall_sum_mat,
                                      NumericVector y,
                                      NumericVector BIC_weights,
                                      int num_iter,int burnin,int num_obs,int num_test_obs,
                                      double a,double sigma,double mu_mu,double nu,
                                      double lambda,//List resids,
                                      NumericMatrix test_data, double lower_prob, double upper_prob, int num_cores){
  
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,test_data);
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  
  arma::vec post_weights_arma = as<arma::vec>(post_weights);
  //int num_models= BIC_weights.size();
  
  IntegerVector num_its_to_sample = RcppArmadillo::rmultinom(num_iter,post_weights);
  IntegerVector num_its_sum = cumsum(num_its_to_sample);
  arma::vec num_its_to_sample_arma = as<arma::vec>(num_its_to_sample);
  arma::vec num_its_sum_arma = as<arma::vec>(num_its_sum);
  
  //arma::mat draws_for_preds(0,num_test_obs);
  //arma::mat draws_for_preds(num_iter,num_test_obs);
  arma::mat draws_for_preds(num_iter,num_test_obs);
  
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  //get exponent
  //double expon=(n+nu)/2;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  
  //Rcout << "Line 5649.\n";
  
  arma::mat yty=y_arma.t()*y_arma;
  
  arma::mat I_test(num_test_obs,num_test_obs);
  I_test=I_test.eye();
  
  //create field (armadillo list) of models
  //each model is a field (armadillo list) of trees represented by matrices
  arma::field<arma::field<arma::mat>> modelsF(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    List temp_tree_list = overall_sum_trees[i];
    //Rcout << "Line 5661.\n";
    
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      //Rcout << "Line 5663.\n";
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    //Rcout << "Line 5669.\n";
    
    modelsF(i)=treesF;
  }
  
  
  arma::field<arma::field<arma::mat>> matsF(overall_sum_mat.size());
  for(int i=0;i<overall_sum_mat.size();i++){
    List temp_tree_list = overall_sum_mat[i];
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    matsF(i)=treesF;
  }
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    
    //Rcout << "Line 5684. i= "<< i << ".\n";
    
    //arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    
    
    arma::field<arma::mat> tree_list = modelsF(i);
    
    arma::mat W(num_obs,0);
    int upsilon=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat curr_tree=(modelsF(i))(j);
      arma::mat curr_obs_nodes=(matsF(i))(j);
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      int b_j=term_nodes.n_elem;
      //begin J function
      
      //will make J as we go in BART-BMA no need to create it again here....
      // arma::mat Jmat=J(curr_obs_nodes,tree_term_nodes);
      arma::mat tree_matrix_temp = (matsF(i))(j);
      arma::mat Jmat(tree_matrix_temp.n_rows, b_j);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(int q=0;q<b_j;q++){
        //double tn=term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        
        //begin find_term_obs
        
        //arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
        //for reference arma_tree_mat == matsF(i)(j) == tree_matrix_temp
        
        arma::uvec term_obs;
        
        for(unsigned int r=0;r<tree_matrix_temp.n_cols;r++){
          //arma::vec colmat=arma_tree_mat.col(r);
          arma::vec colmat=tree_matrix_temp.col(r);
          term_obs=arma::find(colmat==term_nodes[q]);
          if(term_obs.size()>0){
            break;
          }
        }
        
        //end find_term_obs
        
        
        //assign term_obs to the correct index of J
        //NumericVector term_obs2=as<NumericVector>(wrap(term_obs));
        //NumericVector obs_col(obs_to_nodes_temp.nrow());
        arma::vec obs_col= arma::zeros<arma::vec>(tree_matrix_temp.n_rows);
        //Rcout << "Line 5747.\n";
        obs_col.elem(term_obs)= arma::ones<arma::vec>(term_obs.n_elem);
        //Rcout << "Line 5749.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(i)= colmat;
      }
      
      
      
      
      W.insert_cols(upsilon,Jmat);
      upsilon+=b_j;
    }
    
    
    //Rcout << "Line 5759.\n";
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    ////////////////////////////////////////
    
    arma::mat W_tilde(num_test_obs,0);
    int upsilon2=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_test_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_test_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde.insert_cols(upsilon2,Jmat);
      upsilon2+=b_j;
    }
    
    ////////////////////////////////////////
    //Rcout << "Line 5819.\n";
    
    
    double b=W.n_cols;
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*W;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=W.t()*W;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    //arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    //arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    arma::mat third_term=W.t()*y_arma;
    
    
    //arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    
    
    arma::mat w_tilde_M_inv =  W_tilde*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    //arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    //arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    
    
    //arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    //arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights_arma[i];
    
    preds_all_models_arma.col(i)=preds_temp_arma;
    
    // Rcout << "Line 4468. i= "<< i << ".\n";
    
    
    //int num_its_to_sample = round(weight*(num_iter));
    
    int num_its_temp = num_its_to_sample_arma[i];
    
    // arma::uword m = num_test_obs;
    // arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    // U = sqrt(nu+num_obs / U);
    // arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    // Rcout << "Line 4483. i= "<< i << ".\n";
    
    
    
    
    //arma::uword m = num_test_obs;
    arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    //arma::vec U = Rcpp::rchisq(num_its_temp,nu+num_obs);
    U = sqrt((nu+num_obs) / U);
    //arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    arma::mat Y = arma::randn(num_its_temp, num_test_obs);
    arma::mat temp_LY = Y*arma::chol(covar_t);
    temp_LY.each_col() %= U;
    arma::mat temp_draws = arma::repmat(preds_temp_arma, 1, num_its_temp).t() +  temp_LY;
    
    
    
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    
    // Rcout << "number of rows of preds_temp_arma = " << preds_temp_arma.n_rows << ".\n";
    // Rcout << "number of rows of temp_LY = " << temp_LY.n_rows << ".\n";
    // 
    // Rcout << "number of rows of draws_for_preds = " << draws_for_preds.n_rows << ".\n";
    // Rcout << "number of rows of temp_draws = " << temp_draws.n_rows << ".\n";
    //
    
    
    //draws_for_preds = join_cols(draws_for_preds,temp_draws);
    
    //draws_for_preds = join_cols(draws_for_preds,temp_draws);
    if(i==0){
      draws_for_preds.rows(0,num_its_sum_arma[i]-1)=temp_draws;
    }else{
      draws_for_preds.rows(num_its_sum_arma[i-1],num_its_sum_arma[i]-1)=temp_draws;
    }
    
  }
  
#pragma omp barrier  
  
  //arma::colvec predicted_values;
  // Rcout << "Line 4491";
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  NumericMatrix output(3, num_test_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  
  
  
  // Rcout << "Line 4492";
  //Rcout << "probs_for_quantiles = " << probs_for_quantiles << ".\n";
  typedef std::vector<double> stdvec;
  
  for(int i=0;i<num_test_obs;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    std::vector<double> tempcol= arma::conv_to<stdvec>::from(draws_for_preds.col(i));
    std::vector<double> tempquant= Quantile2(tempcol, probs_for_quantiles);
    NumericVector tempforoutput = wrap(tempquant);
    output(_,i)= tempforoutput;
    
  }  
  
  for(int i=0;i<output.ncol();i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output(_,i)=get_original(min(y),max(y),-0.5,0.5, output(_,i));
  }  
  
  
  List ret(2);
  ret[0]= output;
  ret[1]= wrap(get_original_arma(min(y),max(y),-0.5,0.5,predicted_values));
  
  
  return(ret);
  
}
//###########################################################################################################################//
// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_chol_parallel_outsamp(List overall_sum_trees,
                                     List overall_sum_mat,
                                     NumericVector y,
                                     NumericVector BIC_weights,
                                     int num_iter,int burnin,int num_obs,int num_test_obs,
                                     double a,double sigma,double mu_mu,double nu,
                                     double lambda,//List resids,
                                     NumericMatrix test_data, double lower_prob,
                                     double upper_prob, int num_cores){
  
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,test_data);
  
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  arma::vec post_weights_arma = as<arma::vec>(post_weights);
  
  //int num_models= BIC_weights.size();
  
  IntegerVector num_its_to_sample = RcppArmadillo::rmultinom(num_iter,post_weights);
  IntegerVector num_its_sum = cumsum(num_its_to_sample);
  
  arma::vec num_its_to_sample_arma = as<arma::vec>(num_its_to_sample);
  arma::vec num_its_sum_arma = as<arma::vec>(num_its_sum);
  
  //arma::mat draws_for_preds(0,num_test_obs);
  arma::mat draws_for_preds(num_iter,num_test_obs);
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  //get exponent
  //double expon=(n+nu)/2;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  
  //Rcout << "Line 5649.\n";
  
  arma::mat yty=y_arma.t()*y_arma;
  
  arma::mat I_test(num_test_obs,num_test_obs);
  I_test=I_test.eye();
  
  
  
  arma::field<arma::field<arma::mat>> modelsF(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    List temp_tree_list = overall_sum_trees[i];
    //Rcout << "Line 5661.\n";
    
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      //Rcout << "Line 5663.\n";
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    //Rcout << "Line 5669.\n";
    
    modelsF(i)=treesF;
  }
  
  
  arma::field<arma::field<arma::mat>> matsF(overall_sum_mat.size());
  for(int i=0;i<overall_sum_mat.size();i++){
    List temp_tree_list = overall_sum_mat[i];
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    matsF(i)=treesF;
  }
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    
    //Rcout << "Line 5684. i= "<< i << ".\n";
    
    //arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    
    
    arma::field<arma::mat> tree_list = modelsF(i);
    
    arma::mat W(num_obs,0);
    int upsilon=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat curr_tree=(modelsF(i))(j);
      arma::mat curr_obs_nodes=(matsF(i))(j);
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      int b_j=term_nodes.n_elem;
      //begin J function
      
      //will make J as we go in BART-BMA no need to create it again here....
      // arma::mat Jmat=J(curr_obs_nodes,tree_term_nodes);
      arma::mat tree_matrix_temp = (matsF(i))(j);
      arma::mat Jmat(tree_matrix_temp.n_rows, b_j);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(int q=0;q<b_j;q++){
        //double tn=term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        
        //begin find_term_obs
        
        //arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
        //for reference arma_tree_mat == matsF(i)(j) == tree_matrix_temp
        
        arma::uvec term_obs;
        
        for(unsigned int r=0;r<tree_matrix_temp.n_cols;r++){
          //arma::vec colmat=arma_tree_mat.col(r);
          arma::vec colmat=tree_matrix_temp.col(r);
          term_obs=arma::find(colmat==term_nodes[q]);
          if(term_obs.size()>0){
            break;
          }
        }
        
        //end find_term_obs
        
        
        //assign term_obs to the correct index of J
        //NumericVector term_obs2=as<NumericVector>(wrap(term_obs));
        //NumericVector obs_col(obs_to_nodes_temp.nrow());
        arma::vec obs_col= arma::zeros<arma::vec>(tree_matrix_temp.n_rows);
        //Rcout << "Line 5747.\n";
        obs_col.elem(term_obs)= arma::ones<arma::vec>(term_obs.n_elem);
        //Rcout << "Line 5749.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(i)= colmat;
      }
      
      
      
      
      W.insert_cols(upsilon,Jmat);
      upsilon+=b_j;
    }
    
    
    //Rcout << "Line 6224.\n";
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    ////////////////////////////////////////
    
    arma::mat W_tilde(num_test_obs,0);
    int upsilon2=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_test_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_test_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde.insert_cols(upsilon2,Jmat);
      upsilon2+=b_j;
    }
    
    ////////////////////////////////////////
    //Rcout << "Line 5819.\n";
    
    
    double b=W.n_cols;
    //arma::vec yvec=Rcpp::as<arma::vec>(y);
    //arma::mat y_arma(num_obs,1);
    //y_arma.col(0)=yvec;
    //get exponent
    //double expon=(n+nu)/2;
    //get y^Tpsi^{-1}y
    // arma::mat psi_inv=psi.i();
    
    
    //arma::mat yty=y_arma.t()*y_arma;
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*W;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=W.t()*W;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    //arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    // arma::mat third_term=Wmat.t()*y_arma;
    
    
    arma::mat coeff = solve(L_mat.t(), solve(L_mat, W.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    //arma::mat mvm= ytW*sec_term_inv*third_term;
    
    
    
    //arma::mat w_tilde_M_inv =  W_tilde*sec_term_inv;
    
    //Rcout << "Line 6348.\n";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    //arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    arma::mat mvm= y_arma.t()*W*coeff;
    
    
    
    //arma::mat I_test(num_test_obs,num_test_obs);
    //I_test=I_test.eye();
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    
    
    arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    //Rcout << "Line 6394. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights_arma[i];
    preds_all_models_arma.col(i)=preds_temp_arma;
    
    //Rcout << "Line 6403. i= "<< i << ".\n";
    
    
    //int num_its_to_sample = round(weight*(num_iter));
    
    int num_its_temp = num_its_to_sample_arma[i];
    
    // arma::uword m = num_test_obs;
    // arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    // U = sqrt(nu+num_obs / U);
    // arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    // Rcout << "Line 4483. i= "<< i << ".\n";
    
    
    
    
    //arma::uword m = num_test_obs;
    arma::vec U = arma::chi2rnd(nu+num_obs,num_its_temp);
    //arma::vec U = Rcpp::rchisq(num_its_temp,nu+num_obs);
    U = sqrt((nu+num_obs) / U);
    //arma::mat Y = arma::mvnrnd( arma::colvec(m, arma::fill::zeros), covar_t,num_its_temp);
    arma::mat Y = arma::randn(num_its_temp, num_test_obs);
    arma::mat temp_LY = Y*arma::chol(covar_t);
    temp_LY.each_col() %= U;
    arma::mat temp_draws = arma::repmat(preds_temp_arma, 1, num_its_temp).t() +  temp_LY;
    
    //Rcout << "Line 6435. i= "<< i << ".\n";
    
    
    // arma::mat temp_draws(m, num_its_temp);
    // // Rcout << "Line 4478. i= "<< i << ".\n";
    // 
    // for ( arma::uword i = 0; i < temp_draws.n_cols; ++i ) {
    //   temp_draws.col(i) = preds_temp_arma + Y.col(i) * U[i];
    // }
    
    // Rcout << "number of rows of preds_temp_arma = " << preds_temp_arma.n_rows << ".\n";
    // Rcout << "number of rows of temp_LY = " << temp_LY.n_rows << ".\n";
    // 
    // Rcout << "number of rows of draws_for_preds = " << draws_for_preds.n_rows << ".\n";
    // Rcout << "number of rows of temp_draws = " << temp_draws.n_rows << ".\n";
    // 
    //draws_for_preds = join_cols(draws_for_preds,temp_draws);
    
    if(i==0){
      draws_for_preds.rows(0,num_its_sum_arma[i]-1)=temp_draws;
    }else{
      draws_for_preds.rows(num_its_sum_arma[i-1],num_its_sum_arma[i]-1)=temp_draws;
    }
    
    //Rcout << "Line 6459. i= "<< i << ".\n";
    
    
  }
  
  //arma::colvec predicted_values;
  //Rcout << "Line 6463.\n";
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  NumericMatrix output(3, num_test_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  
  
  
  typedef std::vector<double> stdvec;
  
  for(int i=0;i<num_test_obs;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    std::vector<double> tempcol= arma::conv_to<stdvec>::from(draws_for_preds.col(i));
    std::vector<double> tempquant= Quantile2(tempcol, probs_for_quantiles);
    NumericVector tempforoutput = wrap(tempquant);
    output(_,i)= tempforoutput;
    
  }  
  
  for(int i=0;i<output.ncol();i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output(_,i)=get_original(min(y),max(y),-0.5,0.5, output(_,i));
  }  
  
  
  List ret(2);
  ret[0]= output;
  ret[1]= wrap(get_original_arma(min(y),max(y),-0.5,0.5,predicted_values));
  
  
  return(ret);
}



//###########################################################################################################################//

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List mean_vars_lin_alg_parallel_outsamp(List overall_sum_trees,
                                        List overall_sum_mat,
                                        NumericVector y,
                                        NumericVector BIC_weights,
                                        int num_iter,int burnin,int num_obs,int num_test_obs,
                                        double a,double sigma,double mu_mu,double nu,
                                        double lambda,//List resids,
                                        NumericMatrix test_data, int num_cores){
  
  //NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y); 
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,test_data);
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size()); 
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size()); 
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  //create field (armadillo list) of models
  //each model is a field (armadillo list) of trees represented by matrices
  arma::field<arma::field<arma::mat>> modelsF(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    List temp_tree_list = overall_sum_trees[i];
    //Rcout << "Line 5661.\n";
    
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      //Rcout << "Line 5663.\n";
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    //Rcout << "Line 5669.\n";
    
    modelsF(i)=treesF;
  }
  
  
  arma::field<arma::field<arma::mat>> matsF(overall_sum_mat.size());
  for(int i=0;i<overall_sum_mat.size();i++){
    List temp_tree_list = overall_sum_mat[i];
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    matsF(i)=treesF;
  }
  
  
  
  //List covar_matrices(BIC_weights.size());
  //NumericVector model_weights(BIC_weights.size());
  arma::vec model_weightsv(BIC_weights.size());
  arma::field<arma::mat> covar_matricesF(BIC_weights.size());
  
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  //get exponent
  //double expon=(n+nu)/2;
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  
  
  arma::mat yty=y_arma.t()*y_arma;
  
  
  arma::mat I_test(num_test_obs,num_test_obs);
  I_test=I_test.eye();
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    //Rcout << "Line 4111";
    
    
    //Rcout << "Line 5684. i= "<< i << ".\n";
    
    //arma::mat Wmat=W(overall_sum_trees[i],overall_sum_mat[i],num_obs);
    
    
    
    
    arma::field<arma::mat> tree_list = modelsF(i);
    
    arma::mat W(num_obs,0);
    int upsilon=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat curr_tree=(modelsF(i))(j);
      arma::mat curr_obs_nodes=(matsF(i))(j);
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      int b_j=term_nodes.n_elem;
      //begin J function
      
      //will make J as we go in BART-BMA no need to create it again here....
      // arma::mat Jmat=J(curr_obs_nodes,tree_term_nodes);
      arma::mat tree_matrix_temp = (matsF(i))(j);
      arma::mat Jmat(tree_matrix_temp.n_rows, b_j);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(int q=0;q<b_j;q++){
        //double tn=term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        
        //begin find_term_obs
        
        //arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
        //for reference arma_tree_mat == matsF(i)(j) == tree_matrix_temp
        
        arma::uvec term_obs;
        
        for(unsigned int r=0;r<tree_matrix_temp.n_cols;r++){
          //arma::vec colmat=arma_tree_mat.col(r);
          arma::vec colmat=tree_matrix_temp.col(r);
          term_obs=arma::find(colmat==term_nodes[q]);
          if(term_obs.size()>0){
            break;
          }
        }
        
        //end find_term_obs
        
        
        //assign term_obs to the correct index of J
        //NumericVector term_obs2=as<NumericVector>(wrap(term_obs));
        //NumericVector obs_col(obs_to_nodes_temp.nrow());
        arma::vec obs_col= arma::zeros<arma::vec>(tree_matrix_temp.n_rows);
        //Rcout << "Line 5747.\n";
        obs_col.elem(term_obs)= arma::ones<arma::vec>(term_obs.n_elem);
        //Rcout << "Line 5749.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(i)= colmat;
      }
      
      
      
      
      W.insert_cols(upsilon,Jmat);
      upsilon+=b_j;
    }
    
    
    //Rcout << "Line 5759.\n";
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    ////////////////////////////////////////
    
    arma::mat W_tilde(num_test_obs,0);
    int upsilon2=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_test_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_test_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde.insert_cols(upsilon2,Jmat);
      upsilon2+=b_j;
    }
    
    ////////////////////////////////////////
    
    double b=W.n_cols;
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*W;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=W.t()*W;
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    //get t(J)inv(psi)y
    arma::mat third_term=W.t()*y_arma;
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    arma::mat w_tilde_M_inv =  W_tilde*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    //arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    arma::vec preds_temp_arma= W_tilde*sec_term_inv*(W.t())*y_arma;
    
    
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    
    
    //Rcout << "Line 4162";
    
    
    double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    //Rcout << "Line 4167";
    
    model_weightsv[i]=weight;
    //Rcout << "Line 4170";
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*weight;
    preds_all_models_arma.col(i)=preds_temp_arma;
    //Rcout << "Line 4173";
    
    covar_matricesF(i)= covar_t;
    //Rcout << "Line 4176";
    
    
    
    
  }
#pragma omp barrier  
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  
  //NumericVector orig_preds=get_original(min(y),max(y),-0.5,0.5,wrap(predicted_values)) ;
  
  
  List ret(4);
  ret[0]=wrap(predicted_values);
  ret[1]=wrap(model_weightsv);
  ret[2]=wrap(preds_all_models_arma.t());
  ret[3]=wrap(covar_matricesF);
  
  return(ret);
}
//###########################################################################################################################//

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_ITE_outsamp_par(List overall_sum_trees,
                               List overall_sum_mat,
                               NumericVector y,
                               NumericVector BIC_weights,//double min_possible,double max_possible,
                               int num_obs,int num_test_obs,
                               double a,double sigma,double mu_mu,double nu,
                               double lambda,//List resids,
                               NumericMatrix test_data, double lower_prob, double upper_prob, int num_cores,
                               double root_alg_precision,
                               NumericMatrix training_data){
  
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  
  NumericVector onevec1(training_data.nrow(),1.0);
  NumericMatrix Test1data= cbind(onevec1,test_data);
  
  NumericVector zerovec1(training_data.nrow(),0.0);
  NumericMatrix Test0data= cbind(zerovec1,test_data);
  
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata1_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,Test1data);
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata0_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,Test0data);
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat t_vars_arma(num_test_obs,BIC_weights.size());
  
  arma::vec cate_means_arma(BIC_weights.size());
  arma::vec cate_means_weighted_arma(BIC_weights.size());
  
  arma::vec cate_vars_arma(BIC_weights.size());
  
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  arma::vec post_weights_arma = as<arma::vec>(post_weights);
  
  //int num_models= BIC_weights.size();
  
  
  arma::vec averagingvec=(1/double(num_test_obs))*arma::ones<arma::vec>(num_test_obs);
  
  
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  arma::mat yty=y_arma.t()*y_arma;
  
  
  //arma::mat I_test(num_test_obs,num_test_obs);
  //I_test=I_test.eye();
  
  //create field (armadillo list) of models
  //each model is a field (armadillo list) of trees represented by matrices
  arma::field<arma::field<arma::mat>> modelsF(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    List temp_tree_list = overall_sum_trees[i];
    //Rcout << "Line 5661.\n";
    
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      //Rcout << "Line 5663.\n";
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    //Rcout << "Line 5669.\n";
    
    modelsF(i)=treesF;
  }
  
  
  arma::field<arma::field<arma::mat>> matsF(overall_sum_mat.size());
  for(int i=0;i<overall_sum_mat.size();i++){
    List temp_tree_list = overall_sum_mat[i];
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    matsF(i)=treesF;
  }
  
  
  
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    
    
    
    
    arma::field<arma::mat> tree_list = modelsF(i);
    
    arma::mat W(num_obs,0);
    int upsilon=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat curr_tree=(modelsF(i))(j);
      arma::mat curr_obs_nodes=(matsF(i))(j);
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      int b_j=term_nodes.n_elem;
      //begin J function
      
      //will make J as we go in BART-BMA no need to create it again here....
      // arma::mat Jmat=J(curr_obs_nodes,tree_term_nodes);
      arma::mat tree_matrix_temp = (matsF(i))(j);
      arma::mat Jmat(tree_matrix_temp.n_rows, b_j);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(int q=0;q<b_j;q++){
        //double tn=term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        
        //begin find_term_obs
        
        //arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
        //for reference arma_tree_mat == matsF(i)(j) == tree_matrix_temp
        
        arma::uvec term_obs;
        
        for(unsigned int r=0;r<tree_matrix_temp.n_cols;r++){
          //arma::vec colmat=arma_tree_mat.col(r);
          arma::vec colmat=tree_matrix_temp.col(r);
          term_obs=arma::find(colmat==term_nodes[q]);
          if(term_obs.size()>0){
            break;
          }
        }
        
        //end find_term_obs
        
        
        //assign term_obs to the correct index of J
        //NumericVector term_obs2=as<NumericVector>(wrap(term_obs));
        //NumericVector obs_col(obs_to_nodes_temp.nrow());
        arma::vec obs_col= arma::zeros<arma::vec>(tree_matrix_temp.n_rows);
        //Rcout << "Line 5747.\n";
        obs_col.elem(term_obs)= arma::ones<arma::vec>(term_obs.n_elem);
        //Rcout << "Line 5749.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(i)= colmat;
      }
      
      
      
      
      W.insert_cols(upsilon,Jmat);
      upsilon+=b_j;
    }
    
    
    //Rcout << "Line 5759.\n";
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    ////////////////////////////////////////
    
    arma::mat W_tilde1(num_test_obs,0);
    int upsilon2=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata1_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_test_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_test_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde1.insert_cols(upsilon2,Jmat);
      upsilon2+=b_j;
    }
    
    ////////////////////////////////////////
    
    ////////////////////////////////////////
    
    arma::mat W_tilde0(num_test_obs,0);
    int upsilon3=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata0_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_test_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_test_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde0.insert_cols(upsilon3,Jmat);
      upsilon3+=b_j;
    }
    
    ////////////////////////////////////////
    
    //Rcout << "Line 5819.\n";
    
    
    double b=W.n_cols;
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*W;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=W.t()*W;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    //arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    //arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    arma::mat third_term=W.t()*y_arma;
    
    
    //arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    arma::mat Treat_diff = W_tilde1-W_tilde0;
    
    arma::mat w_tilde_M_inv =  Treat_diff*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    //arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    //arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    arma::mat covar_t=temp_scal*(w_tilde_M_inv*(Treat_diff.t()));
    
    arma::mat catevartemp=temp_scal*(averagingvec.t()*w_tilde_M_inv*(Treat_diff.t())*averagingvec);
    
    //arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    //arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights_arma[i];
    
    preds_all_models_arma.col(i)=preds_temp_arma;
    t_vars_arma.col(i)=covar_t.diag();
    
    cate_means_arma(i)=as_scalar(averagingvec.t()*preds_temp_arma);
    cate_means_weighted_arma(i)=cate_means_arma(i)*post_weights_arma(i); 
    
    cate_vars_arma(i)=as_scalar(catevartemp);
    
  }
  
  //}
#pragma omp barrier  
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  double cate_pred=sum(cate_means_weighted_arma);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  arma::mat output(3, num_test_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  
  //std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  arma::mat cate_ints(3, 1);
  
  
  
  typedef std::vector<double> stdvec;
  std::vector<double> weights_vec= as<stdvec>(post_weights);
  
  boost::math::students_t dist2(nu+num_obs);
  double lq_tstandard= boost::math::quantile(dist2,lower_prob);
  double med_tstandard= boost::math::quantile(dist2,0.5); //This is just 0 ??
  double uq_tstandard= boost::math::quantile(dist2,upper_prob);
  
  
  if(weights_vec.size()==1){
    
    cate_ints(0,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*lq_tstandard;
    cate_ints(1,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*med_tstandard;
    cate_ints(2,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*uq_tstandard;   
    
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_test_obs;i++){
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      //boost::math::students_t dist2(nu+num_obs);
      
      //Rcout << "Line 13812 tempvars" << t_vars_arma.row(i) << ".\n";
      //Rcout << "tempmeans" << preds_all_models_arma.row(i) << ".\n";
      
      
      output(0,i)= tempmeans[0]+sqrt(tempvars[0])*lq_tstandard;
      output(1,i)= tempmeans[0]+sqrt(tempvars[0])*med_tstandard;
      output(2,i)= tempmeans[0]+sqrt(tempvars[0])*uq_tstandard;
      
      
    }
#pragma omp barrier  
  }else{
    std::vector<double> tempmeans_cate= arma::conv_to<stdvec>::from(cate_means_arma);
    std::vector<double> tempvars_cate= arma::conv_to<stdvec>::from(cate_vars_arma);
    
    std::vector<double> bounds_lQ_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, lq_tstandard);
    
    //Rcout << "line 13828 cate_vars_arma = " << cate_vars_arma << ".\n";
    //Rcout << "cate_means_arma = " << cate_means_arma << ".\n";
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(0,0)= rootmixt(nu+num_obs,  
              bounds_lQ_cate[0]-0.0001,
              bounds_lQ_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, lower_prob,root_alg_precision);
    
    std::vector<double> bounds_med_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, med_tstandard);
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(1,0)= rootmixt(nu+num_obs,  
              bounds_med_cate[0]-0.0001,  
              bounds_med_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, 0.5, root_alg_precision);
    
    std::vector<double> bounds_uQ_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, uq_tstandard);
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(2,0)= rootmixt(nu+num_obs,  
              bounds_uQ_cate[0]-0.0001,  
              bounds_uQ_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, upper_prob, root_alg_precision);
    
    //Rcout << "line 13871 cate_ints = " << cate_ints << ".\n";
    
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_test_obs;i++){
      //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
      
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      //Rcout << "Line 13859. tempvars" << t_vars_arma.row(i) << ".\n";
      //Rcout << "tempmeans" << preds_all_models_arma.row(i) << ".\n";
      
      std::vector<double> bounds_lQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, lq_tstandard);
      
      output(0,i)=rootmixt(nu+num_obs,  
             bounds_lQ[0]-0.0001,  
             bounds_lQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, lower_prob,root_alg_precision);
      
      
      std::vector<double> bounds_med = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, med_tstandard);
      
      output(1,i)=rootmixt(nu+num_obs,  
             bounds_med[0]-0.0001,  
             bounds_med[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, 0.5,root_alg_precision);
      
      std::vector<double> bounds_uQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, uq_tstandard);
      
      output(2,i)=rootmixt(nu+num_obs, 
             bounds_uQ[0]-0.0001,  
             bounds_uQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, upper_prob,root_alg_precision);
      
      
    }  
#pragma omp barrier  
  }
  
  arma::mat output_rescaled(output.n_rows, output.n_cols);
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(unsigned int i=0;i<output.n_cols;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output_rescaled.col(i)=get_original_TE_arma(min(y),max(y),-0.5,0.5, output.col(i));
    
    
  }  
#pragma omp barrier  
  
  arma::mat cate_ints_rescaled=get_original_TE_arma(min(y),max(y),-0.5,0.5, cate_ints.col(0));
  
  
  List ret(4);
  ret[0]= wrap(output_rescaled);
  ret[1]= wrap(get_original_TE_arma(min(y),max(y),-0.5,0.5,predicted_values));
  ret[2]= get_original_TE_double(min(y),max(y),-0.5,0.5,cate_pred);
  ret[3]= wrap(cate_ints_rescaled);
  
  return(ret);
  
}

//###########################################################################################################################//

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_ITE_insamp_par(List overall_sum_trees,
                              List overall_sum_mat,
                              NumericVector y,
                              NumericVector BIC_weights,//double min_possible,double max_possible,
                              int num_obs,//int num_test_obs,
                              double a,double sigma,
                              double mu_mu,double nu,
                              double lambda,//List resids,NumericMatrix test_data, 
                              double lower_prob, double upper_prob, 
                              int num_cores,
                              double root_alg_precision,
                              NumericMatrix training_data){
  
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  
  NumericVector onevec1(training_data.nrow(),1.0);
  NumericMatrix Training1data= cbind(onevec1,training_data);
  
  NumericVector zerovec1(training_data.nrow(),0.0);
  NumericMatrix Training0data= cbind(zerovec1,training_data);
  
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata1_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,Training1data);
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata0_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,Training0data);
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_obs,BIC_weights.size());
  arma::mat t_vars_arma(num_obs,BIC_weights.size());
  
  arma::vec cate_means_arma(BIC_weights.size());
  arma::vec cate_means_weighted_arma(BIC_weights.size());
  
  arma::vec cate_vars_arma(BIC_weights.size());
  
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  arma::vec post_weights_arma = as<arma::vec>(post_weights);
  
  //int num_models= BIC_weights.size();
  
  
  arma::vec averagingvec=(1/double(num_obs))*arma::ones<arma::vec>(num_obs);
  
  
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  arma::mat yty=y_arma.t()*y_arma;
  
  
  //arma::mat I_test(num_test_obs,num_test_obs);
  //I_test=I_test.eye();
  
  //create field (armadillo list) of models
  //each model is a field (armadillo list) of trees represented by matrices
  arma::field<arma::field<arma::mat>> modelsF(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    List temp_tree_list = overall_sum_trees[i];
    //Rcout << "Line 5661.\n";
    
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      //Rcout << "Line 5663.\n";
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    //Rcout << "Line 5669.\n";
    
    modelsF(i)=treesF;
  }
  
  
  arma::field<arma::field<arma::mat>> matsF(overall_sum_mat.size());
  for(int i=0;i<overall_sum_mat.size();i++){
    List temp_tree_list = overall_sum_mat[i];
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    matsF(i)=treesF;
  }
  
  
  
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    
    
    
    
    arma::field<arma::mat> tree_list = modelsF(i);
    
    arma::mat W(num_obs,0);
    int upsilon=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat curr_tree=(modelsF(i))(j);
      arma::mat curr_obs_nodes=(matsF(i))(j);
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      int b_j=term_nodes.n_elem;
      //begin J function
      
      //will make J as we go in BART-BMA no need to create it again here....
      // arma::mat Jmat=J(curr_obs_nodes,tree_term_nodes);
      arma::mat tree_matrix_temp = (matsF(i))(j);
      arma::mat Jmat(tree_matrix_temp.n_rows, b_j);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(int q=0;q<b_j;q++){
        //double tn=term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        
        //begin find_term_obs
        
        //arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
        //for reference arma_tree_mat == matsF(i)(j) == tree_matrix_temp
        
        arma::uvec term_obs;
        
        for(unsigned int r=0;r<tree_matrix_temp.n_cols;r++){
          //arma::vec colmat=arma_tree_mat.col(r);
          arma::vec colmat=tree_matrix_temp.col(r);
          term_obs=arma::find(colmat==term_nodes[q]);
          if(term_obs.size()>0){
            break;
          }
        }
        
        //end find_term_obs
        
        
        //assign term_obs to the correct index of J
        //NumericVector term_obs2=as<NumericVector>(wrap(term_obs));
        //NumericVector obs_col(obs_to_nodes_temp.nrow());
        arma::vec obs_col= arma::zeros<arma::vec>(tree_matrix_temp.n_rows);
        //Rcout << "Line 5747.\n";
        obs_col.elem(term_obs)= arma::ones<arma::vec>(term_obs.n_elem);
        //Rcout << "Line 5749.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(i)= colmat;
      }
      
      
      
      
      W.insert_cols(upsilon,Jmat);
      upsilon+=b_j;
    }
    
    
    //Rcout << "Line 5759.\n";
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    ////////////////////////////////////////
    
    arma::mat W_tilde1(num_obs,0);
    int upsilon2=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata1_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde1.insert_cols(upsilon2,Jmat);
      upsilon2+=b_j;
    }
    
    ////////////////////////////////////////
    
    ////////////////////////////////////////
    
    arma::mat W_tilde0(num_obs,0);
    int upsilon3=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata0_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde0.insert_cols(upsilon3,Jmat);
      upsilon3+=b_j;
    }
    
    ////////////////////////////////////////
    
    //Rcout << "Line 5819.\n";
    
    
    double b=W.n_cols;
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*W;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=W.t()*W;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    //arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    //arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    arma::mat third_term=W.t()*y_arma;
    
    
    //arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    arma::mat Treat_diff = W_tilde1-W_tilde0;
    
    arma::mat w_tilde_M_inv =  Treat_diff*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    //arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    //arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    arma::mat covar_t=temp_scal*(w_tilde_M_inv*(Treat_diff.t()));
    
    arma::mat catevartemp=temp_scal*(averagingvec.t()*w_tilde_M_inv*(Treat_diff.t())*averagingvec);
    
    //arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    //arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights_arma(i);
    
    preds_all_models_arma.col(i)=preds_temp_arma;
    t_vars_arma.col(i)=covar_t.diag();
    
    cate_means_arma(i)=as_scalar(averagingvec.t()*preds_temp_arma);
    cate_means_weighted_arma(i)=cate_means_arma(i)*post_weights_arma(i); 
    
    cate_vars_arma(i)=as_scalar(catevartemp);
    
  }
  
  //}
#pragma omp barrier  
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  double cate_pred=sum(cate_means_weighted_arma);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  arma::mat output(3, num_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  
  //std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  arma::mat cate_ints(3, 1);
  
  
  
  typedef std::vector<double> stdvec;
  std::vector<double> weights_vec= as<stdvec>(post_weights);
  
  boost::math::students_t dist2(nu+num_obs);
  double lq_tstandard= boost::math::quantile(dist2,lower_prob);
  double med_tstandard= boost::math::quantile(dist2,0.5); //This is just 0 ??
  double uq_tstandard= boost::math::quantile(dist2,upper_prob);
  
  //Rcout << "cate_vars_arma line 13795" << cate_vars_arma << ".\n";
  //Rcout << "cate_means_arma" << cate_means_arma << ".\n";
  
  if(weights_vec.size()==1){
    
    cate_ints(0,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*lq_tstandard;
    cate_ints(1,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*med_tstandard;
    cate_ints(2,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*uq_tstandard;   
    
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_obs;i++){
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      //boost::math::students_t dist2(nu+num_obs);
      
      //Rcout << "Line 13812 tempvars" << t_vars_arma.row(i) << ".\n";
      //Rcout << "tempmeans" << preds_all_models_arma.row(i) << ".\n";
      
      
      output(0,i)= tempmeans[0]+sqrt(tempvars[0])*lq_tstandard;
      output(1,i)= tempmeans[0]+sqrt(tempvars[0])*med_tstandard;
      output(2,i)= tempmeans[0]+sqrt(tempvars[0])*uq_tstandard;
      
      
    }
#pragma omp barrier  
  }else{
    std::vector<double> tempmeans_cate= arma::conv_to<stdvec>::from(cate_means_arma);
    std::vector<double> tempvars_cate= arma::conv_to<stdvec>::from(cate_vars_arma);
    
    std::vector<double> bounds_lQ_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, lq_tstandard);
    
    //Rcout << "line 13828 cate_vars_arma = " << cate_vars_arma << ".\n";
    //Rcout << "cate_means_arma = " << cate_means_arma << ".\n";
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(0,0)= rootmixt(nu+num_obs,  
              bounds_lQ_cate[0]-0.0001,
              bounds_lQ_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, lower_prob,root_alg_precision);
    
    std::vector<double> bounds_med_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, med_tstandard);
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(1,0)= rootmixt(nu+num_obs,  
              bounds_med_cate[0]-0.0001,  
              bounds_med_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, 0.5, root_alg_precision);
    
    std::vector<double> bounds_uQ_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, uq_tstandard);
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(2,0)= rootmixt(nu+num_obs,  
              bounds_uQ_cate[0]-0.0001,  
              bounds_uQ_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, upper_prob, root_alg_precision);
    
    //Rcout << "line 13871 cate_ints = " << cate_ints << ".\n";
    
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_obs;i++){
      //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
      
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      //Rcout << "Line 13859. tempvars" << t_vars_arma.row(i) << ".\n";
      //Rcout << "tempmeans" << preds_all_models_arma.row(i) << ".\n";
      
      std::vector<double> bounds_lQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, lq_tstandard);
      
      output(0,i)=rootmixt(nu+num_obs,  
             bounds_lQ[0]-0.0001,  
             bounds_lQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, lower_prob,root_alg_precision);
      
      
      std::vector<double> bounds_med = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, med_tstandard);
      
      output(1,i)=rootmixt(nu+num_obs,  
             bounds_med[0]-0.0001,  
             bounds_med[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, 0.5,root_alg_precision);
      
      std::vector<double> bounds_uQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, uq_tstandard);
      
      output(2,i)=rootmixt(nu+num_obs, 
             bounds_uQ[0]-0.0001,  
             bounds_uQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, upper_prob,root_alg_precision);
      
      
    }  
#pragma omp barrier  
  }
  arma::mat output_rescaled(output.n_rows, output.n_cols);
  
  double min_y = min(y);
  double max_y = max(y);
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(unsigned int i=0;i<output.n_cols;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output_rescaled.col(i)=get_original_TE_arma(min_y,max_y,-0.5,0.5, output.col(i));
    
    
  }  
#pragma omp barrier  
  
  arma::mat cate_ints_rescaled=get_original_TE_arma(min_y,max_y,-0.5,0.5, cate_ints.col(0));
  
  
  List ret(4);
  ret[0]= wrap(output_rescaled);
  ret[1]= wrap(get_original_TE_arma(min_y,max_y,-0.5,0.5,predicted_values));
  ret[2]= get_original_TE_double(min_y,max_y,-0.5,0.5,cate_pred);
  ret[3]= wrap(cate_ints_rescaled);
  
  return(ret);
  
  
}

//###########################################################################################################################//

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_ITE_CATT_outsamp_par(List overall_sum_trees,
                                    List overall_sum_mat,
                                    NumericVector y,
                                    NumericVector BIC_weights,//double min_possible,double max_possible,
                                    int num_obs,int num_test_obs,
                                    double a,double sigma,double mu_mu,double nu,
                                    double lambda,//List resids,
                                    NumericMatrix test_data, double lower_prob, double upper_prob, int num_cores,
                                    double root_alg_precision,
                                    NumericMatrix training_data,NumericVector ztest){
  
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  
  NumericVector onevec1(training_data.nrow(),1.0);
  NumericMatrix Test1data= cbind(onevec1,test_data);
  
  NumericVector zerovec1(training_data.nrow(),0.0);
  NumericMatrix Test0data= cbind(zerovec1,test_data);
  
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata1_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,Test1data);
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata0_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,Test0data);
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_test_obs,BIC_weights.size());
  arma::mat t_vars_arma(num_test_obs,BIC_weights.size());
  
  arma::vec cate_means_arma(BIC_weights.size());
  arma::vec cate_means_weighted_arma(BIC_weights.size());
  arma::vec cate_vars_arma(BIC_weights.size());
  
  arma::vec catt_means_arma(BIC_weights.size());
  arma::vec catt_means_weighted_arma(BIC_weights.size());
  arma::vec catt_vars_arma(BIC_weights.size());
  
  arma::vec catnt_means_arma(BIC_weights.size());
  arma::vec catnt_means_weighted_arma(BIC_weights.size());
  arma::vec catnt_vars_arma(BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  arma::vec post_weights_arma = as<arma::vec>(post_weights);
  arma::vec ztest_arma = as<arma::vec>(ztest);
  
  
  //int num_models= BIC_weights.size();
  
  
  arma::vec averagingvec=(1/double(num_test_obs))*arma::ones<arma::vec>(num_test_obs);
  
  arma::vec catt_averagingvec=(1/arma::sum(ztest_arma))*ztest_arma;
  arma::vec catnt_averagingvec=(1/(num_test_obs-arma::sum(ztest_arma)))*(1-ztest_arma);
  
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  arma::mat yty=y_arma.t()*y_arma;
  
  
  //arma::mat I_test(num_test_obs,num_test_obs);
  //I_test=I_test.eye();
  
  //create field (armadillo list) of models
  //each model is a field (armadillo list) of trees represented by matrices
  arma::field<arma::field<arma::mat>> modelsF(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    List temp_tree_list = overall_sum_trees[i];
    //Rcout << "Line 5661.\n";
    
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      //Rcout << "Line 5663.\n";
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    //Rcout << "Line 5669.\n";
    
    modelsF(i)=treesF;
  }
  
  
  arma::field<arma::field<arma::mat>> matsF(overall_sum_mat.size());
  for(int i=0;i<overall_sum_mat.size();i++){
    List temp_tree_list = overall_sum_mat[i];
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    matsF(i)=treesF;
  }
  
  
  
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    
    
    
    
    arma::field<arma::mat> tree_list = modelsF(i);
    
    arma::mat W(num_obs,0);
    int upsilon=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat curr_tree=(modelsF(i))(j);
      arma::mat curr_obs_nodes=(matsF(i))(j);
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      int b_j=term_nodes.n_elem;
      //begin J function
      
      //will make J as we go in BART-BMA no need to create it again here....
      // arma::mat Jmat=J(curr_obs_nodes,tree_term_nodes);
      arma::mat tree_matrix_temp = (matsF(i))(j);
      arma::mat Jmat(tree_matrix_temp.n_rows, b_j);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(int q=0;q<b_j;q++){
        //double tn=term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        
        //begin find_term_obs
        
        //arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
        //for reference arma_tree_mat == matsF(i)(j) == tree_matrix_temp
        
        arma::uvec term_obs;
        
        for(unsigned int r=0;r<tree_matrix_temp.n_cols;r++){
          //arma::vec colmat=arma_tree_mat.col(r);
          arma::vec colmat=tree_matrix_temp.col(r);
          term_obs=arma::find(colmat==term_nodes[q]);
          if(term_obs.size()>0){
            break;
          }
        }
        
        //end find_term_obs
        
        
        //assign term_obs to the correct index of J
        //NumericVector term_obs2=as<NumericVector>(wrap(term_obs));
        //NumericVector obs_col(obs_to_nodes_temp.nrow());
        arma::vec obs_col= arma::zeros<arma::vec>(tree_matrix_temp.n_rows);
        //Rcout << "Line 5747.\n";
        obs_col.elem(term_obs)= arma::ones<arma::vec>(term_obs.n_elem);
        //Rcout << "Line 5749.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(i)= colmat;
      }
      
      
      
      
      W.insert_cols(upsilon,Jmat);
      upsilon+=b_j;
    }
    
    
    //Rcout << "Line 5759.\n";
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    ////////////////////////////////////////
    
    arma::mat W_tilde1(num_test_obs,0);
    int upsilon2=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata1_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_test_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_test_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde1.insert_cols(upsilon2,Jmat);
      upsilon2+=b_j;
    }
    
    ////////////////////////////////////////
    
    ////////////////////////////////////////
    
    arma::mat W_tilde0(num_test_obs,0);
    int upsilon3=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata0_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_test_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_test_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde0.insert_cols(upsilon3,Jmat);
      upsilon3+=b_j;
    }
    
    ////////////////////////////////////////
    
    //Rcout << "Line 5819.\n";
    
    
    double b=W.n_cols;
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*W;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=W.t()*W;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    //arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    //arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    arma::mat third_term=W.t()*y_arma;
    
    
    //arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    arma::mat Treat_diff = W_tilde1-W_tilde0;
    
    arma::mat w_tilde_M_inv =  Treat_diff*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    //arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    //arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    arma::mat covar_t=temp_scal*(w_tilde_M_inv*(Treat_diff.t()));
    
    arma::mat catevartemp=temp_scal*(averagingvec.t()*w_tilde_M_inv*(Treat_diff.t())*averagingvec);
    
    arma::mat cattvartemp=temp_scal*(catt_averagingvec.t()*w_tilde_M_inv*(Treat_diff.t())*catt_averagingvec);
    arma::mat catntvartemp=temp_scal*(catnt_averagingvec.t()*w_tilde_M_inv*(Treat_diff.t())*catnt_averagingvec);
    
    //arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    //arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights_arma[i];
    
    preds_all_models_arma.col(i)=preds_temp_arma;
    t_vars_arma.col(i)=covar_t.diag();
    
    cate_means_arma(i)=as_scalar(averagingvec.t()*preds_temp_arma);
    cate_means_weighted_arma(i)=cate_means_arma(i)*post_weights_arma(i); 
    
    cate_vars_arma(i)=as_scalar(catevartemp);
    
    
    catt_means_arma(i)=as_scalar(catt_averagingvec.t()*preds_temp_arma);
    catt_means_weighted_arma(i)=catt_means_arma(i)*post_weights_arma(i); 
    
    catt_vars_arma(i)=as_scalar(cattvartemp);
    
    
    catnt_means_arma(i)=as_scalar(catnt_averagingvec.t()*preds_temp_arma);
    catnt_means_weighted_arma(i)=catnt_means_arma(i)*post_weights_arma(i); 
    
    catnt_vars_arma(i)=as_scalar(catntvartemp);
    
    
  }
  
  //}
#pragma omp barrier  
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  double cate_pred=sum(cate_means_weighted_arma);
  double catt_pred=sum(catt_means_weighted_arma);
  double catnt_pred=sum(catnt_means_weighted_arma);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  arma::mat output(3, num_test_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  
  //std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  arma::mat cate_ints(3, 1);
  arma::mat catt_ints(3, 1);
  arma::mat catnt_ints(3, 1);
  
  
  
  typedef std::vector<double> stdvec;
  std::vector<double> weights_vec= as<stdvec>(post_weights);
  
  boost::math::students_t dist2(nu+num_obs);
  double lq_tstandard= boost::math::quantile(dist2,lower_prob);
  double med_tstandard= boost::math::quantile(dist2,0.5); //This is just 0 ??
  double uq_tstandard= boost::math::quantile(dist2,upper_prob);
  
  
  if(weights_vec.size()==1){
    
    cate_ints(0,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*lq_tstandard;
    cate_ints(1,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*med_tstandard;
    cate_ints(2,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*uq_tstandard;   
    
    catt_ints(0,0)= catt_means_arma(0)+sqrt(catt_vars_arma(0))*lq_tstandard;
    catt_ints(1,0)= catt_means_arma(0)+sqrt(catt_vars_arma(0))*med_tstandard;
    catt_ints(2,0)= catt_means_arma(0)+sqrt(catt_vars_arma(0))*uq_tstandard;   
    
    catnt_ints(0,0)= catnt_means_arma(0)+sqrt(catnt_vars_arma(0))*lq_tstandard;
    catnt_ints(1,0)= catnt_means_arma(0)+sqrt(catnt_vars_arma(0))*med_tstandard;
    catnt_ints(2,0)= catnt_means_arma(0)+sqrt(catnt_vars_arma(0))*uq_tstandard;   
    
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_test_obs;i++){
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      //boost::math::students_t dist2(nu+num_obs);
      
      //Rcout << "Line 13812 tempvars" << t_vars_arma.row(i) << ".\n";
      //Rcout << "tempmeans" << preds_all_models_arma.row(i) << ".\n";
      
      
      output(0,i)= tempmeans[0]+sqrt(tempvars[0])*lq_tstandard;
      output(1,i)= tempmeans[0]+sqrt(tempvars[0])*med_tstandard;
      output(2,i)= tempmeans[0]+sqrt(tempvars[0])*uq_tstandard;
      
      
    }
#pragma omp barrier  
  }else{
    std::vector<double> tempmeans_cate= arma::conv_to<stdvec>::from(cate_means_arma);
    std::vector<double> tempvars_cate= arma::conv_to<stdvec>::from(cate_vars_arma);
    
    std::vector<double> bounds_lQ_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, lq_tstandard);
    
    //Rcout << "line 13828 cate_vars_arma = " << cate_vars_arma << ".\n";
    //Rcout << "cate_means_arma = " << cate_means_arma << ".\n";
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(0,0)= rootmixt(nu+num_obs,  
              bounds_lQ_cate[0]-0.0001,
              bounds_lQ_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, lower_prob,root_alg_precision);
    
    std::vector<double> bounds_med_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, med_tstandard);
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(1,0)= rootmixt(nu+num_obs,  
              bounds_med_cate[0]-0.0001,  
              bounds_med_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, 0.5, root_alg_precision);
    
    std::vector<double> bounds_uQ_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, uq_tstandard);
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(2,0)= rootmixt(nu+num_obs,  
              bounds_uQ_cate[0]-0.0001,  
              bounds_uQ_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, upper_prob, root_alg_precision);
    
    //Rcout << "line 13871 cate_ints = " << cate_ints << ".\n";
    
    
    //
    
    std::vector<double> tempmeans_catt= arma::conv_to<stdvec>::from(catt_means_arma);
    std::vector<double> tempvars_catt= arma::conv_to<stdvec>::from(catt_vars_arma);
    
    std::vector<double> bounds_lQ_catt = mixt_find_boundsQ( nu+num_obs, tempmeans_catt, tempvars_catt, lq_tstandard);
    
    
    catt_ints(0,0)= rootmixt(nu+num_obs,  
              bounds_lQ_catt[0]-0.0001,
              bounds_lQ_catt[1]+0.0001,
              tempmeans_catt,
              tempvars_catt,
              weights_vec, lower_prob,root_alg_precision);
    
    std::vector<double> bounds_med_catt = mixt_find_boundsQ( nu+num_obs, tempmeans_catt, tempvars_catt, med_tstandard);
    
    //Rcout << "bounds_lQ_catt[0] = " << bounds_lQ_catt[0] << ".\n"; 
    //Rcout << "bounds_lQ_catt[1] = " << bounds_lQ_catt[1] << ".\n";
    
    catt_ints(1,0)= rootmixt(nu+num_obs,  
              bounds_med_catt[0]-0.0001,  
              bounds_med_catt[1]+0.0001,
              tempmeans_catt,
              tempvars_catt,
              weights_vec, 0.5, root_alg_precision);
    
    std::vector<double> bounds_uQ_catt = mixt_find_boundsQ( nu+num_obs, tempmeans_catt, tempvars_catt, uq_tstandard);
    
    //Rcout << "bounds_lQ_catt[0] = " << bounds_lQ_catt[0] << ".\n"; 
    //Rcout << "bounds_lQ_catt[1] = " << bounds_lQ_catt[1] << ".\n";
    
    catt_ints(2,0)= rootmixt(nu+num_obs,  
              bounds_uQ_catt[0]-0.0001,  
              bounds_uQ_catt[1]+0.0001,
              tempmeans_catt,
              tempvars_catt,
              weights_vec, upper_prob, root_alg_precision);
    
    
    
    //
    
    
    std::vector<double> tempmeans_catnt= arma::conv_to<stdvec>::from(catnt_means_arma);
    std::vector<double> tempvars_catnt= arma::conv_to<stdvec>::from(catnt_vars_arma);
    
    std::vector<double> bounds_lQ_catnt = mixt_find_boundsQ( nu+num_obs, tempmeans_catnt, tempvars_catnt, lq_tstandard);
    
    
    
    catnt_ints(0,0)= rootmixt(nu+num_obs,  
               bounds_lQ_catnt[0]-0.0001,
               bounds_lQ_catnt[1]+0.0001,
               tempmeans_catnt,
               tempvars_catnt,
               weights_vec, lower_prob,root_alg_precision);
    
    std::vector<double> bounds_med_catnt = mixt_find_boundsQ( nu+num_obs, tempmeans_catnt, tempvars_catnt, med_tstandard);
    
    //Rcout << "bounds_lQ_catnt[0] = " << bounds_lQ_catnt[0] << ".\n"; 
    //Rcout << "bounds_lQ_catnt[1] = " << bounds_lQ_catnt[1] << ".\n";
    
    catnt_ints(1,0)= rootmixt(nu+num_obs,  
               bounds_med_catnt[0]-0.0001,  
               bounds_med_catnt[1]+0.0001,
               tempmeans_catnt,
               tempvars_catnt,
               weights_vec, 0.5, root_alg_precision);
    
    std::vector<double> bounds_uQ_catnt = mixt_find_boundsQ( nu+num_obs, tempmeans_catnt, tempvars_catnt, uq_tstandard);
    
    //Rcout << "bounds_lQ_catnt[0] = " << bounds_lQ_catnt[0] << ".\n"; 
    //Rcout << "bounds_lQ_catnt[1] = " << bounds_lQ_catnt[1] << ".\n";
    
    catnt_ints(2,0)= rootmixt(nu+num_obs,  
               bounds_uQ_catnt[0]-0.0001,  
               bounds_uQ_catnt[1]+0.0001,
               tempmeans_catnt,
               tempvars_catnt,
               weights_vec, upper_prob, root_alg_precision);
    
    
    
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_test_obs;i++){
      //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
      
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      //Rcout << "Line 13859. tempvars" << t_vars_arma.row(i) << ".\n";
      //Rcout << "tempmeans" << preds_all_models_arma.row(i) << ".\n";
      
      std::vector<double> bounds_lQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, lq_tstandard);
      
      output(0,i)=rootmixt(nu+num_obs,  
             bounds_lQ[0]-0.0001,  
             bounds_lQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, lower_prob,root_alg_precision);
      
      
      std::vector<double> bounds_med = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, med_tstandard);
      
      output(1,i)=rootmixt(nu+num_obs,  
             bounds_med[0]-0.0001,  
             bounds_med[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, 0.5,root_alg_precision);
      
      std::vector<double> bounds_uQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, uq_tstandard);
      
      output(2,i)=rootmixt(nu+num_obs, 
             bounds_uQ[0]-0.0001,  
             bounds_uQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, upper_prob,root_alg_precision);
      
      
    }  
#pragma omp barrier  
  }
  
  arma::mat output_rescaled(output.n_rows, output.n_cols);
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(unsigned int i=0;i<output.n_cols;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output_rescaled.col(i)=get_original_TE_arma(min(y),max(y),-0.5,0.5, output.col(i));
    
    
  }  
#pragma omp barrier  
  
  arma::mat cate_ints_rescaled=get_original_TE_arma(min(y),max(y),-0.5,0.5, cate_ints.col(0));
  arma::mat catt_ints_rescaled=get_original_TE_arma(min(y),max(y),-0.5,0.5, catt_ints.col(0));
  arma::mat catnt_ints_rescaled=get_original_TE_arma(min(y),max(y),-0.5,0.5, catnt_ints.col(0));
  
  
  List ret(8);
  ret[0]= wrap(output_rescaled);
  ret[1]= wrap(get_original_TE_arma(min(y),max(y),-0.5,0.5,predicted_values));
  ret[2]= get_original_TE_double(min(y),max(y),-0.5,0.5,cate_pred);
  ret[3]= wrap(cate_ints_rescaled);
  ret[4]= get_original_TE_double(min(y),max(y),-0.5,0.5,catt_pred);
  ret[5]= wrap(catt_ints_rescaled);
  ret[6]= get_original_TE_double(min(y),max(y),-0.5,0.5,catnt_pred);
  ret[7]= wrap(catnt_ints_rescaled);
  
  
  return(ret);
  
}

//###########################################################################################################################//

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List pred_ints_ITE_CATT_insamp_par(List overall_sum_trees,
                                   List overall_sum_mat,
                                   NumericVector y,
                                   NumericVector BIC_weights,//double min_possible,double max_possible,
                                   int num_obs,//int num_test_obs,
                                   double a,double sigma,
                                   double mu_mu,double nu,
                                   double lambda,//List resids,NumericMatrix test_data, 
                                   double lower_prob, double upper_prob, 
                                   int num_cores,
                                   double root_alg_precision,
                                   NumericMatrix training_data,
                                   NumericVector ztrain){
  
  
  //List termobs_testdata_overall= get_termobs_testdata_overall(overall_sum_trees,test_data);
  
  NumericVector onevec1(training_data.nrow(),1.0);
  NumericMatrix Training1data= cbind(onevec1,training_data);
  
  NumericVector zerovec1(training_data.nrow(),0.0);
  NumericMatrix Training0data= cbind(zerovec1,training_data);
  
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata1_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,Training1data);
  arma::field<arma::field<arma::field<arma::uvec>>> termobs_testdata0_overallF=get_termobs_testdata_fields_overall(overall_sum_trees,Training0data);
  
  //NumericMatrix preds_all_models(num_test_obs,BIC_weights.size());
  arma::mat preds_all_models_arma(num_obs,BIC_weights.size());
  arma::mat weighted_preds_all_models_arma(num_obs,BIC_weights.size());
  arma::mat t_vars_arma(num_obs,BIC_weights.size());
  
  arma::vec cate_means_arma(BIC_weights.size());
  arma::vec cate_means_weighted_arma(BIC_weights.size());
  arma::vec cate_vars_arma(BIC_weights.size());
  
  arma::vec catt_means_arma(BIC_weights.size());
  arma::vec catt_means_weighted_arma(BIC_weights.size());
  arma::vec catt_vars_arma(BIC_weights.size());
  
  arma::vec catnt_means_arma(BIC_weights.size());
  arma::vec catnt_means_weighted_arma(BIC_weights.size());
  arma::vec catnt_vars_arma(BIC_weights.size());
  
  // for all sums of trees
  
  NumericVector BICi=-0.5*BIC_weights;
  double max_BIC=max(BICi);
  
  
  NumericVector post_weights(BIC_weights.size());
  
  for(int k=0;k<BIC_weights.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    post_weights[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  arma::vec post_weights_arma = as<arma::vec>(post_weights);
  
  
  arma::vec ztrain_arma = as<arma::vec>(ztrain);
  
  
  
  arma::vec averagingvec=(1/double(num_obs))*arma::ones<arma::vec>(num_obs);
  
  arma::vec catt_averagingvec=(1/arma::sum(ztrain_arma))*ztrain_arma;
  arma::vec catnt_averagingvec=(1/(num_obs-arma::sum(ztrain_arma)))*(1-ztrain_arma);
  
  
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y);
  arma::vec yvec=Rcpp::as<arma::vec>(y_scaled);
  arma::mat y_arma(num_obs,1);
  y_arma.col(0)=yvec;
  arma::mat yty=y_arma.t()*y_arma;
  
  
  //arma::mat I_test(num_test_obs,num_test_obs);
  //I_test=I_test.eye();
  
  //create field (armadillo list) of models
  //each model is a field (armadillo list) of trees represented by matrices
  arma::field<arma::field<arma::mat>> modelsF(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    List temp_tree_list = overall_sum_trees[i];
    //Rcout << "Line 5661.\n";
    
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      //Rcout << "Line 5663.\n";
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    //Rcout << "Line 5669.\n";
    
    modelsF(i)=treesF;
  }
  
  
  arma::field<arma::field<arma::mat>> matsF(overall_sum_mat.size());
  for(int i=0;i<overall_sum_mat.size();i++){
    List temp_tree_list = overall_sum_mat[i];
    arma::field<arma::mat> treesF(temp_tree_list.size());
    for(int q=0;q<temp_tree_list.size();q++){
      arma::mat temp_tree_mat = Rcpp::as<arma::mat>(temp_tree_list[q]);
      treesF(q)=temp_tree_mat;
    }
    matsF(i)=treesF;
  }
  
  
  
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(int i=0;i<overall_sum_trees.size();i++){
    
    
    
    
    arma::field<arma::mat> tree_list = modelsF(i);
    
    arma::mat W(num_obs,0);
    int upsilon=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat curr_tree=(modelsF(i))(j);
      arma::mat curr_obs_nodes=(matsF(i))(j);
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      int b_j=term_nodes.n_elem;
      //begin J function
      
      //will make J as we go in BART-BMA no need to create it again here....
      // arma::mat Jmat=J(curr_obs_nodes,tree_term_nodes);
      arma::mat tree_matrix_temp = (matsF(i))(j);
      arma::mat Jmat(tree_matrix_temp.n_rows, b_j);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(int q=0;q<b_j;q++){
        //double tn=term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        
        //begin find_term_obs
        
        //arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
        //for reference arma_tree_mat == matsF(i)(j) == tree_matrix_temp
        
        arma::uvec term_obs;
        
        for(unsigned int r=0;r<tree_matrix_temp.n_cols;r++){
          //arma::vec colmat=arma_tree_mat.col(r);
          arma::vec colmat=tree_matrix_temp.col(r);
          term_obs=arma::find(colmat==term_nodes[q]);
          if(term_obs.size()>0){
            break;
          }
        }
        
        //end find_term_obs
        
        
        //assign term_obs to the correct index of J
        //NumericVector term_obs2=as<NumericVector>(wrap(term_obs));
        //NumericVector obs_col(obs_to_nodes_temp.nrow());
        arma::vec obs_col= arma::zeros<arma::vec>(tree_matrix_temp.n_rows);
        //Rcout << "Line 5747.\n";
        obs_col.elem(term_obs)= arma::ones<arma::vec>(term_obs.n_elem);
        //Rcout << "Line 5749.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(i)= colmat;
      }
      
      
      
      
      W.insert_cols(upsilon,Jmat);
      upsilon+=b_j;
    }
    
    
    //Rcout << "Line 5759.\n";
    
    
    //arma::mat W_tilde=get_W_test(overall_sum_trees[i],termobs_testdata_overall[i],num_test_obs);
    
    
    ////////////////////////////////////////
    
    arma::mat W_tilde1(num_obs,0);
    int upsilon2=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata1_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde1.insert_cols(upsilon2,Jmat);
      upsilon2+=b_j;
    }
    
    ////////////////////////////////////////
    
    ////////////////////////////////////////
    
    arma::mat W_tilde0(num_obs,0);
    int upsilon3=0;
    for(unsigned int j=0;j<modelsF(i).n_elem;j++){
      
      arma::mat  curr_tree=(modelsF(i))(j);
      arma::field<arma::uvec> curr_termobs=(termobs_testdata0_overallF(i))(j);
      
      //begin find termnodes
      //NumericVector tree_term_nodes=find_term_nodes(curr_tree);
      
      //arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
      arma::vec colmat=curr_tree.col(4);
      arma::uvec term_nodes=arma::find(colmat==-1);
      term_nodes=term_nodes+1;
      
      //return(wrap(term_nodes));
      //Rcout << "Line 5782.\n";
      
      //end find termnodes
      
      int b_j=term_nodes.size();
      //will make J as we go in BART-BMA no need to create it again here....
      //arma::mat Jmat=get_J_test(curr_termobs,tree_term_nodes,n);
      //arma::mat Jmat=get_J_test(curr_termobs,term_nodes,num_test_obs);
      
      //begin J test function
      arma::mat Jmat(num_obs, term_nodes.n_elem);
      Jmat.zeros();
      
      //for each terminal node get the observations associated with it and set column
      for(unsigned int q=0;q<term_nodes.n_elem;q++){
        //double tn=tree_term_nodes[q];
        //arma::uvec term_obs=find_term_obs(obs_to_nodes_temp,tn);
        //assign term_obs to the correct index of J
        arma::uvec term_obs2=curr_termobs(q);
        arma::vec obs_col= arma::zeros<arma::vec>(num_obs);
        //Rcout << "Line 5809.\n";
        obs_col.elem(term_obs2)=arma::ones<arma::vec>(term_obs2.n_elem);
        //Rcout << "Line 5811.\n";
        //arma::vec colmat=Rcpp::as<arma::vec>(obs_col);
        Jmat.col(q)= obs_col;
        
        // arma::vec  colmat=arma::zeros(Jmat.n_rows) ;// colmattest(Jmat.n_rows,0);
        // colmat.elem(term_obs).fill(1);
        // Jmat.col(q)= colmat;
      }
      //return(Jmat);
      
      //end J test function
      
      W_tilde0.insert_cols(upsilon3,Jmat);
      upsilon3+=b_j;
    }
    
    ////////////////////////////////////////
    
    //Rcout << "Line 5819.\n";
    
    
    double b=W.n_cols;
    
    
    //get t(y)inv(psi)J
    arma::mat ytW=y_arma.t()*W;
    
    
    //get t(J)inv(psi)J  
    arma::mat WtW=W.t()*W;
    
    
    
    
    
    //get jpsij +aI
    arma::mat aI(b,b);
    aI=a*aI.eye();
    arma::mat sec_term=WtW+aI;
    
    
    //arma::mat L_mat = arma::chol(sec_term,"lower");
    
    //arma::mat L_inv = arma::inv(L_mat);
    //arma::mat L_inv = arma::inv(arma::chol(sec_term));
    //arma::mat L_inv = arma::inv((L_mat));
    //arma::mat L_inv_t = arma::trans(arma::inv(trimatu(L_mat)));
    //arma::mat L_inv_t = arma::trans(arma::inv((L_mat)));
    
    
    
    //arma::mat sec_term_inv=sec_term.i();
    arma::mat sec_term_inv=inv_sympd(sec_term);  
    
    
    //get t(J)inv(psi)y
    arma::mat third_term=W.t()*y_arma;
    
    
    //arma::mat coeff = solve(L_mat.t(), solve(L_mat, Wmat.t()*y_arma));
    //arma::mat coeff = L_inv.t()*L_inv*Wmat.t()*y_arma;
    
    
    
    //get m^TV^{-1}m
    arma::mat mvm= ytW*sec_term_inv*third_term;
    
    arma::mat Treat_diff = W_tilde1-W_tilde0;
    
    arma::mat w_tilde_M_inv =  Treat_diff*sec_term_inv;
    
    //Rcout << "Line 4151";
    
    
    
    // //Obtain (lower triangular?) matrix t(L) by Cholesky decomposition such that sec_term_inv=L*t(L)
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // //obtain the log of the root of the determinant
    // double rootisum = arma::sum(log(rooti.diag()));
    // 
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    // 
    // arma::mat rel=(b/2)*log(a)+rootisum -expon*log(nu*lambda - arma::sum(LtWtY%LtWtY) +yty);
    // 
    // 
    
    
    // arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sec_term))));
    // arma::mat LtWtY= rooti*(Wmat.t()*y);
    
    
    
    arma::vec preds_temp_arma= w_tilde_M_inv*third_term;
    
    // arma::vec preds_temp_arma= W_tilde*sec_term_inv*(Wmat.t())*y_arma;
    
    //arma::vec preds_temp_arma= W_tilde*coeff;
    //arma::mat mvm= coeff*y_arma;
    //arma::mat mvm= y_arma.t()*Wmat*coeff;
    
    
    
    
    
    arma::mat temp_for_scal = ((nu*lambda+yty-mvm)/(nu+num_obs));
    double temp_scal= as_scalar(temp_for_scal) ;
    //Rcout << "Line 4156";
    //arma::mat covar_t=temp_scal*(I_test+w_tilde_M_inv*(W_tilde.t()));
    arma::mat covar_t=temp_scal*(w_tilde_M_inv*(Treat_diff.t()));
    
    arma::mat catevartemp=temp_scal*(averagingvec.t()*w_tilde_M_inv*(Treat_diff.t())*averagingvec);
    arma::mat cattvartemp=temp_scal*(catt_averagingvec.t()*w_tilde_M_inv*(Treat_diff.t())*catt_averagingvec);
    arma::mat catntvartemp=temp_scal*(catnt_averagingvec.t()*w_tilde_M_inv*(Treat_diff.t())*catnt_averagingvec);
    
    //arma::mat W_tilde_L_inv_t= W_tilde*L_inv_t; 
    //arma::mat covar_t=temp_scal*(I_test+W_tilde_L_inv_t*(W_tilde_L_inv_t.t()));
    
    
    
    // Rcout << "Line 4459. i= "<< i << ".\n";
    
    
    
    //double weight=exp(BICi[i]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    weighted_preds_all_models_arma.col(i)=preds_temp_arma*post_weights_arma[i];
    
    preds_all_models_arma.col(i)=preds_temp_arma;
    t_vars_arma.col(i)=covar_t.diag();
    
    cate_means_arma(i)=as_scalar(averagingvec.t()*preds_temp_arma);
    cate_means_weighted_arma(i)=cate_means_arma(i)*post_weights_arma(i); 
    
    cate_vars_arma(i)=as_scalar(catevartemp);
    
    
    catt_means_arma(i)=as_scalar(catt_averagingvec.t()*preds_temp_arma);
    catt_means_weighted_arma(i)=catt_means_arma(i)*post_weights_arma(i); 
    
    catt_vars_arma(i)=as_scalar(cattvartemp);
    
    
    catnt_means_arma(i)=as_scalar(catnt_averagingvec.t()*preds_temp_arma);
    catnt_means_weighted_arma(i)=catnt_means_arma(i)*post_weights_arma(i); 
    
    catnt_vars_arma(i)=as_scalar(catntvartemp);
    
  }
  
  //}
#pragma omp barrier  
  
  //arma::colvec predicted_values;
  
  //arma::mat M1(preds_all_models.begin(), preds_all_models.nrow(), preds_all_models.ncol(), false);
  arma::colvec predicted_values=sum(weighted_preds_all_models_arma,1);
  
  double cate_pred=sum(cate_means_weighted_arma);
  double catt_pred=sum(catt_means_weighted_arma);
  double catnt_pred=sum(catnt_means_weighted_arma);
  
  //NumericMatrix draws_wrapped= wrap(draws_for_preds);
  arma::mat output(3, num_obs);
  //NumericVector probs_for_quantiles =  NumericVector::create(lower_prob, 0.5, upper_prob);
  
  //std::vector<double> probs_for_quantiles {lower_prob, 0.5, upper_prob};
  arma::mat cate_ints(3, 1);
  arma::mat catt_ints(3, 1);
  arma::mat catnt_ints(3, 1);
  
  
  typedef std::vector<double> stdvec;
  std::vector<double> weights_vec= as<stdvec>(post_weights);
  
  boost::math::students_t dist2(nu+num_obs);
  double lq_tstandard= boost::math::quantile(dist2,lower_prob);
  double med_tstandard= boost::math::quantile(dist2,0.5); //This is just 0 ??
  double uq_tstandard= boost::math::quantile(dist2,upper_prob);
  
  //Rcout << "cate_vars_arma line 13795" << cate_vars_arma << ".\n";
  //Rcout << "cate_means_arma" << cate_means_arma << ".\n";
  
  if(weights_vec.size()==1){
    
    cate_ints(0,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*lq_tstandard;
    cate_ints(1,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*med_tstandard;
    cate_ints(2,0)= cate_means_arma(0)+sqrt(cate_vars_arma(0))*uq_tstandard;
    
    catt_ints(0,0)= catt_means_arma(0)+sqrt(catt_vars_arma(0))*lq_tstandard;
    catt_ints(1,0)= catt_means_arma(0)+sqrt(catt_vars_arma(0))*med_tstandard;
    catt_ints(2,0)= catt_means_arma(0)+sqrt(catt_vars_arma(0))*uq_tstandard;   
    
    catnt_ints(0,0)= catnt_means_arma(0)+sqrt(catnt_vars_arma(0))*lq_tstandard;
    catnt_ints(1,0)= catnt_means_arma(0)+sqrt(catnt_vars_arma(0))*med_tstandard;
    catnt_ints(2,0)= catnt_means_arma(0)+sqrt(catnt_vars_arma(0))*uq_tstandard; 
    
    
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_obs;i++){
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      //boost::math::students_t dist2(nu+num_obs);
      
      //Rcout << "Line 13812 tempvars" << t_vars_arma.row(i) << ".\n";
      //Rcout << "tempmeans" << preds_all_models_arma.row(i) << ".\n";
      
      
      output(0,i)= tempmeans[0]+sqrt(tempvars[0])*lq_tstandard;
      output(1,i)= tempmeans[0]+sqrt(tempvars[0])*med_tstandard;
      output(2,i)= tempmeans[0]+sqrt(tempvars[0])*uq_tstandard;
      
      
    }
#pragma omp barrier  
  }else{
    std::vector<double> tempmeans_cate= arma::conv_to<stdvec>::from(cate_means_arma);
    std::vector<double> tempvars_cate= arma::conv_to<stdvec>::from(cate_vars_arma);
    
    std::vector<double> bounds_lQ_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, lq_tstandard);
    
    //Rcout << "line 13828 cate_vars_arma = " << cate_vars_arma << ".\n";
    //Rcout << "cate_means_arma = " << cate_means_arma << ".\n";
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(0,0)= rootmixt(nu+num_obs,  
              bounds_lQ_cate[0]-0.0001,
              bounds_lQ_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, lower_prob,root_alg_precision);
    
    std::vector<double> bounds_med_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, med_tstandard);
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(1,0)= rootmixt(nu+num_obs,  
              bounds_med_cate[0]-0.0001,  
              bounds_med_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, 0.5, root_alg_precision);
    
    std::vector<double> bounds_uQ_cate = mixt_find_boundsQ( nu+num_obs, tempmeans_cate, tempvars_cate, uq_tstandard);
    
    //Rcout << "bounds_lQ_cate[0] = " << bounds_lQ_cate[0] << ".\n"; 
    //Rcout << "bounds_lQ_cate[1] = " << bounds_lQ_cate[1] << ".\n";
    
    cate_ints(2,0)= rootmixt(nu+num_obs,  
              bounds_uQ_cate[0]-0.0001,  
              bounds_uQ_cate[1]+0.0001,
              tempmeans_cate,
              tempvars_cate,
              weights_vec, upper_prob, root_alg_precision);
    
    //Rcout << "line 13871 cate_ints = " << cate_ints << ".\n";
    
    
    
    
    //
    
    std::vector<double> tempmeans_catt= arma::conv_to<stdvec>::from(catt_means_arma);
    std::vector<double> tempvars_catt= arma::conv_to<stdvec>::from(catt_vars_arma);
    
    std::vector<double> bounds_lQ_catt = mixt_find_boundsQ( nu+num_obs, tempmeans_catt, tempvars_catt, lq_tstandard);
    
    
    catt_ints(0,0)= rootmixt(nu+num_obs,  
              bounds_lQ_catt[0]-0.0001,
              bounds_lQ_catt[1]+0.0001,
              tempmeans_catt,
              tempvars_catt,
              weights_vec, lower_prob,root_alg_precision);
    
    std::vector<double> bounds_med_catt = mixt_find_boundsQ( nu+num_obs, tempmeans_catt, tempvars_catt, med_tstandard);
    
    //Rcout << "bounds_lQ_catt[0] = " << bounds_lQ_catt[0] << ".\n"; 
    //Rcout << "bounds_lQ_catt[1] = " << bounds_lQ_catt[1] << ".\n";
    
    catt_ints(1,0)= rootmixt(nu+num_obs,  
              bounds_med_catt[0]-0.0001,  
              bounds_med_catt[1]+0.0001,
              tempmeans_catt,
              tempvars_catt,
              weights_vec, 0.5, root_alg_precision);
    
    std::vector<double> bounds_uQ_catt = mixt_find_boundsQ( nu+num_obs, tempmeans_catt, tempvars_catt, uq_tstandard);
    
    //Rcout << "bounds_lQ_catt[0] = " << bounds_lQ_catt[0] << ".\n"; 
    //Rcout << "bounds_lQ_catt[1] = " << bounds_lQ_catt[1] << ".\n";
    
    catt_ints(2,0)= rootmixt(nu+num_obs,  
              bounds_uQ_catt[0]-0.0001,  
              bounds_uQ_catt[1]+0.0001,
              tempmeans_catt,
              tempvars_catt,
              weights_vec, upper_prob, root_alg_precision);
    
    
    
    //
    
    
    std::vector<double> tempmeans_catnt= arma::conv_to<stdvec>::from(catnt_means_arma);
    std::vector<double> tempvars_catnt= arma::conv_to<stdvec>::from(catnt_vars_arma);
    
    std::vector<double> bounds_lQ_catnt = mixt_find_boundsQ( nu+num_obs, tempmeans_catnt, tempvars_catnt, lq_tstandard);
    
    
    
    catnt_ints(0,0)= rootmixt(nu+num_obs,  
               bounds_lQ_catnt[0]-0.0001,
               bounds_lQ_catnt[1]+0.0001,
               tempmeans_catnt,
               tempvars_catnt,
               weights_vec, lower_prob,root_alg_precision);
    
    std::vector<double> bounds_med_catnt = mixt_find_boundsQ( nu+num_obs, tempmeans_catnt, tempvars_catnt, med_tstandard);
    
    //Rcout << "bounds_lQ_catnt[0] = " << bounds_lQ_catnt[0] << ".\n"; 
    //Rcout << "bounds_lQ_catnt[1] = " << bounds_lQ_catnt[1] << ".\n";
    
    catnt_ints(1,0)= rootmixt(nu+num_obs,  
               bounds_med_catnt[0]-0.0001,  
               bounds_med_catnt[1]+0.0001,
               tempmeans_catnt,
               tempvars_catnt,
               weights_vec, 0.5, root_alg_precision);
    
    std::vector<double> bounds_uQ_catnt = mixt_find_boundsQ( nu+num_obs, tempmeans_catnt, tempvars_catnt, uq_tstandard);
    
    //Rcout << "bounds_lQ_catnt[0] = " << bounds_lQ_catnt[0] << ".\n"; 
    //Rcout << "bounds_lQ_catnt[1] = " << bounds_lQ_catnt[1] << ".\n";
    
    catnt_ints(2,0)= rootmixt(nu+num_obs,  
               bounds_uQ_catnt[0]-0.0001,  
               bounds_uQ_catnt[1]+0.0001,
               tempmeans_catnt,
               tempvars_catnt,
               weights_vec, upper_prob, root_alg_precision);
    
    
    
#pragma omp parallel num_threads(num_cores)
#pragma omp for
    for(int i=0;i<num_obs;i++){
      //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
      
      std::vector<double> tempmeans= arma::conv_to<stdvec>::from(preds_all_models_arma.row(i));
      std::vector<double> tempvars= arma::conv_to<stdvec>::from(t_vars_arma.row(i));
      
      //Rcout << "Line 13859. tempvars" << t_vars_arma.row(i) << ".\n";
      //Rcout << "tempmeans" << preds_all_models_arma.row(i) << ".\n";
      
      std::vector<double> bounds_lQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, lq_tstandard);
      
      output(0,i)=rootmixt(nu+num_obs,  
             bounds_lQ[0]-0.0001,  
             bounds_lQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, lower_prob,root_alg_precision);
      
      
      std::vector<double> bounds_med = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, med_tstandard);
      
      output(1,i)=rootmixt(nu+num_obs,  
             bounds_med[0]-0.0001,  
             bounds_med[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, 0.5,root_alg_precision);
      
      std::vector<double> bounds_uQ = mixt_find_boundsQ( nu+num_obs, tempmeans, tempvars, uq_tstandard);
      
      output(2,i)=rootmixt(nu+num_obs, 
             bounds_uQ[0]-0.0001,  
             bounds_uQ[1]+0.0001,
             tempmeans,
             tempvars,
             weights_vec, upper_prob,root_alg_precision);
      
      
    }  
#pragma omp barrier  
  }
  
  arma::mat output_rescaled(output.n_rows, output.n_cols);
  
  
#pragma omp parallel num_threads(num_cores)
#pragma omp for
  for(unsigned int i=0;i<output.n_cols;i++){
    //output(_,i)=Quantile(draws_wrapped(_,i), probs_for_quantiles);
    
    output_rescaled.col(i)=get_original_TE_arma(min(y),max(y),-0.5,0.5, output.col(i));
    
    
  }  
#pragma omp barrier  
  
  arma::mat cate_ints_rescaled=get_original_TE_arma(min(y),max(y),-0.5,0.5, cate_ints.col(0));
  arma::mat catt_ints_rescaled=get_original_TE_arma(min(y),max(y),-0.5,0.5, catt_ints.col(0));
  arma::mat catnt_ints_rescaled=get_original_TE_arma(min(y),max(y),-0.5,0.5, catnt_ints.col(0));
  
  
  List ret(8);
  ret[0]= wrap(output_rescaled);
  ret[1]= wrap(get_original_TE_arma(min(y),max(y),-0.5,0.5,predicted_values));
  ret[2]= get_original_TE_double(min(y),max(y),-0.5,0.5,cate_pred);
  ret[3]= wrap(cate_ints_rescaled);
  ret[4]= get_original_TE_double(min(y),max(y),-0.5,0.5,catt_pred);
  ret[5]= wrap(catt_ints_rescaled);
  ret[6]= get_original_TE_double(min(y),max(y),-0.5,0.5,catnt_pred);
  ret[7]= wrap(catnt_ints_rescaled);
  
  
  
  
  return(ret);
  
}
