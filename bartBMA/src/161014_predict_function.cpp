//####################################################################################//
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_internal_nodes_pred(NumericMatrix treetable){
  NumericVector internal_nodes;
  
  for(int l=0;l<treetable.nrow();l++){    
    if(treetable(l,4)==1){
      internal_nodes.push_back(l+1);
    }
  }
  
  NumericVector internal_nodes_sort = clone(internal_nodes);
  std::sort(internal_nodes.begin(), internal_nodes.end());
  
  return(internal_nodes_sort);
}
//####################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_nodes_pred(NumericMatrix tree_table){
  NumericVector terminal_nodes;
  
  for(int i=0;i<tree_table.nrow();i++){
    if(tree_table(i,4)==-1){
      terminal_nodes.push_back(i+1);
    }
  }
  return(terminal_nodes);
}  
//###################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector get_original_pred(double low,double high,double sp_low,double sp_high,NumericVector sum_preds){
  NumericVector original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);
  return(original_y);
}
//###################################################################################//
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector bartBMA_get_testdata_term_obs_pred(NumericMatrix test_data,NumericMatrix tree_data,NumericVector term_node_means) {
  // Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.
  
  //test_data is a nxp matrix with the same variable names as the training data the model was built on...
  //should have an error check for this later on or at least check it's the same dimension
  
  //tree_data is the tree table with the tree information i.e. split points and split variables and terminal node mean values
  
  //term_node_means is a vector storing the terminal node mean values
  
  arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
  arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);
  NumericVector internal_nodes=find_internal_nodes_pred(tree_data);
  NumericVector terminal_nodes=find_term_nodes_pred(tree_data);
  arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
  NumericVector tree_predictions;
  
  //now for each internal node find the observations that belong to the terminal nodes
  
  NumericVector predictions(test_data.nrow());
  if(terminal_nodes.size()==1){
    double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
    predictions=rep(nodemean,test_data.nrow());
  }
  else{
    for(int i=0;i<terminal_nodes.size();i++){
      arma::mat subdata=testd;
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
        }else{
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
      //arma::uvec col_indices=seq_len(testd.n_cols);
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
      double nodemean=tree_data(terminal_nodes[i]-1,5);
      IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
      predictions[predind]= nodemean;
      
    } 
  }
  return(predictions);
}

//##############################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

List get_BART_BMA_test_predictions(NumericMatrix test_data,NumericVector BIC,List sum_trees,NumericVector y_minmax){
  //this will take in a set of sum of trees and loop through each tree in each set.
  //for each tree in each set:
  //call bartBMA_get_testdata_term_obs_pred() and get the predicted values
  //add up values for each set
  //at end weight summed predictions by BIC and put back on original scale
  NumericMatrix sum_tree_preds(test_data.nrow(),sum_trees.size());
  
  for(int i=0;i<sum_trees.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = sum_trees[i];
    
    NumericVector test_preds_sum_tree;
    
    if(is<List>(s)){
      //if set of trees is a list then loop through the list and get predicted values
      List tree_set=sum_trees[i];
      for(int j=0;j<tree_set.size();j++){
        
        NumericMatrix tree_data=tree_set[j];
        NumericVector term_nodes= find_term_nodes_pred(tree_data);
        NumericVector term_node_means;
        for(int k=0;k<term_nodes.size();k++){
          term_nodes[k]=term_nodes[k]-1;
          term_node_means.push_back(tree_data(term_nodes[k],5));
        }
        
        NumericVector test_preds_tree;
        
        if(j==0){          
          test_preds_tree=bartBMA_get_testdata_term_obs_pred(test_data,tree_data,term_node_means);
          test_preds_sum_tree=test_preds_tree;
        }else{
          test_preds_tree=bartBMA_get_testdata_term_obs_pred(test_data,tree_data,term_node_means);
          test_preds_sum_tree=test_preds_sum_tree+test_preds_tree;
        }
        
      }
    }else{
      
      //else there is only one tree in the list element not a sum of trees
      NumericMatrix tree_data=sum_trees[i];
      NumericVector term_nodes= find_term_nodes_pred(tree_data);
      NumericVector term_node_means;
      
      for(int k=0;k<term_nodes.size();k++){
        term_nodes[k]=term_nodes[k]-1;
        term_node_means.push_back(tree_data(term_nodes[k],5));
      }
      
      test_preds_sum_tree=bartBMA_get_testdata_term_obs_pred(test_data,tree_data,term_node_means);
      
    }
    //now have the summed preds for the sum_of_trees add it to column of summed preds
    sum_tree_preds(_,i)= test_preds_sum_tree; 
  }
  
  //now have predictions for each sum_of_trees. Next need to weight each tree prediction by posterior probability and add up.
  NumericMatrix overall_test_preds(sum_tree_preds.nrow(),sum_tree_preds.ncol());  
  for(int k=0;k<BIC.size();k++){  
    NumericVector temp_test_preds=sum_tree_preds(_,k);
    // NumericVector orig_temp_preds=get_original_pred(min(y),max(y),-0.5,0.5,temp_preds) ;
    //double weight=-0.5*BIC[k]/sum(-0.5*BIC);
    
    NumericVector BICi=-0.5*BIC;
    double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    
    overall_test_preds(_,k) = temp_test_preds*weight;
  }  
  
  //sum over all the weighted predictions;
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);
  arma::colvec test_predictions=sum(M2,1);  
  NumericVector orig_test_preds=get_original_pred(y_minmax[0],y_minmax[1],-0.5,0.5,wrap(test_predictions)) ;
  
  //next put predictions on the original scale
  List ret(2);
  ret(0)=orig_test_preds;
  ret(1)= sum_tree_preds;
  return(ret);  
}
