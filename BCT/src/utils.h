#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <random>
#include <map>
#include <set>
#include <iterator>
#include <Rcpp.h>
#include <math.h>



using namespace std;

//Important inputs that are used throught the functions

int m; // alphabet size
int D; // max depth
long double beta; // prior hyper-parameter
long double alpha; // prior hyper-parameter (is function of beta- defined above)
unsigned int k_max; //top-k trees for k-bct algorithm
vector <short> zeros;
vector <short> xn; // input string transformed into vector
map <char, short> encoder; // each input character is mapped to a unique integer
map <short, char> decoder; // stores the decoding rules 

#include "node.h"
#include "Tree_properties.h" //definitions for nodes and trees

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// FUNCTION PROTOTYPES

//1. Functions for initialisation/preprocessing
void init_tree(tree& T);                  // initialises the tree for bct/kbct
void preproc(vector <node *> init);             // does the preprocessing calculations needed for k-bct

void set_param(string &s, int depth, short kmax);  //initialises the parameters
void set_global_parameters(string &s, int depth, short kmax); //used for initialising the global parameters in bct (where kmax is not needed)
void set_global_parameters(string &s, int depth, short kmax, double b); //used for initialising the global parameters in kbct 



//2. Functions used for building Tmax
void update(tree& T, short s, vector <short> ct);    // updates tree for each new symbol
void occur(node * N, short s);                       // node occurrence updates
void insert(tree& T, vector <short> ct2, short ch);  // inserts node with context ct2 to T, links with existing child


//3. Important-main algorithms
long double ctw(tree& T);                      // CTW algorithm
long double bct(tree& T);                      // finds MAP tree
void ctw_bct(tree& T, Tree_properties & tp);                         // runs CTW and BCT together


void kbct_forw(tree &T, vector <node *> init);                      // forward pass of kbct, takes improper and gives improper tree, calculates vectors lm, c for each node
void kbct_back(vector <node *> init, tree T, vector <tree> &trees); // backward pass of kbct
void kbct(tree& T, vector <tree> &trees, vector <node *> init, vector<double> & odds); //runs ctw and kbct together

vector <Tree_properties> build_kbct(); // used in R
double compute_mle(tree &T) ; // computes the likelihood of the MLE tree
long double build_ctw_rcpp(); // performs only the ctw algorithm for the R ctw function



Tree_properties build_bct(); // performs the bct algorithm and transmits the results in a form so that the R environment can interpret
vector <Tree_properties> build_kbct(); // performs the kbct algorithm and transmits the results in a form so that the R environment can interpret
Rcpp::NumericVector compute_log_loss(vector <short> xn, int train_size); // computes the log-loss incurred by prediction 
Rcpp::List online_predict( int train_size); // predicts the most probable character
Rcpp::List map_param(); // finds the MAP model and outputs the probability vectors for each leaf

//4.The remaining secondary functions, which are auxiliary things (but are needed to perform the rest in C++)

int show_leaves(tree T);                       // prints leaves of tree and returns number of leaves
void collect_leaves(tree &T, Tree_properties &tp); // collects the important results after the bct/kbct is performed
void label(tree& T);                           // sets contexts to the nodes of a tree
void copytree(tree  &T, tree &T2);                         // copies only non-deleted nodes to store proper pruned tree
void copy(tree &Tout, tree &Tin);


void comb_initial3(int d, vector <node *> init);                        //used in preprocessing calculations
void comb(int d, int k, tree &T, vector <node *> init);                 // combinations for k-bct forward pass
vector<vector<double> > cartesian_prod(const vector<vector<double>>& v);  //combinations for doubles
vector<vector<short> > cartesian_prod_int(const vector<vector<short>>& v); //combinations for shorts
void makeproper(tree& T); // takes an improper tree and makes it proper


//5. BIC and AIC scores
void counts(tree &T);         // finds counts for proper tree T in sequence xn
double compute_mle(tree &T);                         // finds log-maximum likelihood (input is proper tree of any depth, with count vectors calculated already by function "counts")
void compute_bic_aic_mle(tree &T, double &bic, double &aic, double &mle);     // finds bic, aic, mle of (proper) tree T, sequence xn -> used in bct/kbct R functions



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//FUNCTION DEFINITIONS


void update(tree& T, short s, vector <short> ct) {
  // updates an existing tree 
  
  node * temp = T[0][0]; // start from root
  
  occur(temp, s);
  
  
  for (int j = 0; j < D; j++) {
    
    if (temp->child[ct[j]] != NULL) { // if child node exists in tree
      
      temp = temp->child[ct[j]];    // move to child
      occur(temp, s);               // child occurred
    }
    
    else {                            // create children up to depth D
      
      vector <short> ct2 = ct;        // context of node to be created
      short ch = 0;                   // shows which child exists
      
      for (int k = 0; k < D - j; k++) {
        
        
        
        //insert node with context ct2 to tree
        insert(T, ct2, ch);
        
        occur(T[ct2.size()].back(), s); // inserted nodes occurs
        
        ch = ct2.back();
        ct2.pop_back();
        
      }
      
      j = D + 5;  // break out of previous loop if a child doesn't exist
      temp->child[ch] = T[ct2.size() + 1].back();
    }
    
  }
  
}
// ======================================================================================================
void occur(node * N, short s) {
  
  N->a[s]++; //update count vectors a_s
  int M = 0;
  for (int i = 0; i < m; i++) { M = M + N->a[i]; }
  
  N->le = 1.0* N->le + log2(1.0* N->a[s] - 0.5) - log2(0.5*m + 1.0* M - 1.0); //update log-estimated probability
  
}
// ======================================================================================================
void insert(tree& T, vector <short> ct2, short ch) {
  
  int d = ct2.size();            // node depth
  node * init = new node(m);        // initialise a node pointer
  T[d].push_back(init);         // add a node to tree at corresponding depth
  if (d == D) {
    
    T[d].back()->leaf = 1;   // no children if leaf
    
  }
  
  else {                     // set address of children given by ch
    
    T[d].back()->child[ch] = T[d + 1].back();
  }
  
}
// ======================================================================================================
int show_leaves(tree T) { // uses node class with contexts, returns number of leaves
  
  int n_leaves = 0;
  
  for (int u = 0; u < D + 1; u++) {
    for (unsigned int v = 0; v < T[u].size(); v++)
    {
      
      if (T[u][v]->leaf == 1) {
        n_leaves++;
      }
    }
  }
  
  return n_leaves;
  
}
// ======================================================================================================
void collect_leaves(tree &T, Tree_properties &tp){
  //Takes a tree and stores in tp.context the leaves'contexts of the tree
  
  tp.n_leaves = 0;
  tp.max_depth = 0;
  
  for (int u = 0; u < D + 1; u++) {
    for (unsigned int v = 0; v < T[u].size(); v++)
    {
      
      if (T[u][v]->leaf == 1) {
        
        string su;
        for (unsigned int m = 0; m < T[u][v]->s.size(); m++)
        {
          
          su.push_back(decoder[T[u][v]->s[m]]);
          
        }
        tp.context.push_back(su);
        tp.n_leaves++;
      }
      
      tp.max_depth = u;
    }
  }
}
// ======================================================================================================
void init_tree(tree& T) { //initialise tree to have D+1 rows (i.e. for depth 0 up to D) and add root node to depth0
  
  T.clear();
  vector <node*> row;
  for (int i = 0; i < D + 1; i++) { T.push_back(row); }
  
  T[0].push_back(new node(m)); //only add root node which always exists
  
  
  if (D == 0) {
    T[0][0]->leaf = 1;
  }
  
}
// ======================================================================================================
long double ctw(tree& T) {                   // algorithm takes improper tree and finds prior-predictive likelihood
  
  for (int d = D; d > -1; d--) {           // loop over levels
    
    for (unsigned int k = 0; k < T[d].size(); k++) {  // loop over nodes of each level
      
      if (d == D) {                   // if at max depth, node is a leaf
        T[d][k]->lw = T[d][k]->le;
        
      }
      else {                         // if node is not a leaf
        
        long double sum = 0;
        
        for (int ch = 0; ch < m; ch++) {
          
          if (T[d][k]->child[ch] != NULL) {           // if child #ch exists
            sum = sum + T[d][k]->child[ch]->lw;     // calculate sum of Le(s)
            
          }
          
        }
        
        //calculate weighted log-prob in two cases as explained in notes for numerical precision
        
        long double delta = T[d][k]->le - sum + log2(beta) - log2(1.0 - beta);
        if (delta < 30) {
          
          T[d][k]->lw = log2(1.0 - beta) + sum + log2(1.0 + pow(2.0, delta));
          
        }
        else {
          T[d][k]->lw = log2(beta) + T[d][k]->le + log2(exp(1))*(pow(2.0, -delta) - pow(2.0, -2.0*delta - 1));
          
        }
        
      }
    }
  }
  
  return T[0][0]->lw;              //output value of weighted prob Pw at root
}
// ======================================================================================================
long double bct(tree& T) {                   // algorithm takes improper tree finds maximal probability at root and makes tree proper, then finds bct
  // lw here is used to store the maximal probablity, not the weighted one
  
  // First forward pass (leaves to root) to calculate maximal probabilities Pm at every node
  
  if (D == 0) { // if iid data
    
    return  T[0][0]->le;              //output value of max prob at root
  }
  
  for (int d = D; d > -1; d--) {           // loop over levels
    
    
    
    for (unsigned int k = 0; k < T[d].size(); k++) {  // loop over initially existing nodes of each level
      
      if (d == D) {                   // if at max depth, node is a leaf (if and only if for improper tree)
        T[d][k]->lw = T[d][k]->le;
        
      }
      
      else {                         // if node is not a leaf
        
        long double sum = 0;
        
        for (short ch = 0; ch < m; ch++) {
          
          if (T[d][k]->child[ch] == NULL) {           // if child #ch does not exist, it is equivalent with
            
            if (d < D - 1) {
              sum = sum + log2(beta);
            }
            
          }
          
          else {                                        // if child ch exists
            
            sum = sum + T[d][k]->child[ch]->lw;       // sum of log-probs at children
            
          }
          
        }
        
        // calculate maximal log-prob as explained in notes
        
        
        
        if (log2(1.0 - 1.0*beta) + sum > log2(beta) + T[d][k]->le) { // maximum achieved by children term
          
          T[d][k]->lw = log2(1.0 - 1.0*beta) + sum;                // set max prob of node
          
          
        }
        
        else {                                                        // maximum achived by curent node
          
          T[d][k]->lw = log2(beta) + T[d][k]->le;                  // set max prob of node and mark to be pruned
          T[d][k]->leaf = 1;
          
          for (short ch = 0; ch < m; ch++) {   // for child # ch of each node
            
            
            if (T[d][k]->child[ch] != NULL) {                       // if child exists
              
              T[d][k]->child[ch]->a[0] = -1;                     // mark child to be destructed
              T[d][k]->child[ch] = 0;                            // destruct connection with child
              
            }
            
          }
        }
        
        
      }
    }
  }
  
  // Then backward pass (root to leaves), to prune tree and destroy the required nodes
  // Use a[0] =-1 to mark nodes that need to be pruned
  
  int *length;
  length = new int[D+1];
  //int length[D + 1];
  
  for (int d = 0; d < D + 1; d++) {
    length[d] = T[d].size();
  }
  
  
  
  for (int d = 0; d < D + 1; d++) { // root to leaves now
    
    int check = 0;
    
    for (int k = 0; k < length[d]; k++) {
      
      if (T[d][k]->a[0] == -1) {  // node was marked to be deleted
        
        for (short ch = 0; ch < m; ch++) {   // for child # ch of each node
          
          if (T[d][k]->child[ch] != NULL) {                       // if child exists
            
            T[d][k]->child[ch]->a[0] = -1;                     // mark child to be destructed
            
          }
        }
        
      }
      
      else {
        
        if (T[d][k]->leaf == 0) {               // if child not a leaf and not deleted
          
          for (short ch = 0; ch < m; ch++) {
            
            if (T[d][k]->child[ch] == NULL) {
              
              
              node * init = new node(m);                // insert child to make tree proper
              T[d][k]->child[ch] = init;             // connect child with parent node
              T[d + 1].push_back(init);              // store at appropriate tree depth
              init->leaf = 1;                        // denote it leaf
              
              if (d < D - 1) {
                init->lw = log2(beta);            // set maximal prob for leaf at depth < D, if leaf is at depth D then logP=0;
              }
              
              
            }
          }
        }
      }
      
      if (T[d][k]->a[0] == -1) { //node doesn't exist in reality
        check++;
      }
    }
    if (check == length[d]) { //then no need to look at higher depths
      return T[0][0]->lw;              // output value of weighted prob Pw at root
    }
  }
  delete length;
  return T[0][0]->lw;              // output value of weighted prob Pw at root
  
}
// ======================================================================================================
void label(tree& T) { // takes as input proper tree with no contexts and writes their context in node.s
  
  for (int d = 0; d < D + 1; d++) {
    
    for (unsigned int k = 0; k < T[d].size(); k++) {
      
      if (T[d][k]->leaf == 0) {
        
        for (short ch = 0; ch < m; ch++) {
          
          
          T[d][k]->child[ch]->s = T[d][k]->s;
          T[d][k]->child[ch]->s.push_back(ch);
        }
        
        
      }
    }
    
  }
  
}



void label_inproper(tree& T) { // takes as input proper tree with no contexts and writes their context in node.s
  
  for (int d = 0; d < D + 1; d++) {
    
    for (unsigned int k = 0; k < T[d].size(); k++) {
      
      if (T[d][k]->leaf == 0) {
        
        for (short ch = 0; ch < m; ch++) {
          
          if(T[d][k]->child[ch] != NULL){
            T[d][k]->child[ch]->s = T[d][k]->s;
            T[d][k]->child[ch]->s.push_back(ch);
          }
        }
        
        
      }
    }
    
  }
  
}


// ======================================================================================================
void copytree(tree &T, tree &T2) { // takes tree from BCT which doesn't delete nodes (only marks them so virtually they dont exist)
  // copies only non-marked nodes to another  tree
  // initialise tree of depth D
  for (int d = 0; d < D + 1; d++) {
    
    unsigned int check = 0;
    
    for (unsigned int k = 0; k < T[d].size(); k++) {
      
      if (T[d][k]->a[0] > -1) {
        
        T2[d].push_back(T[d][k]);
      }
      else {
        check++;
      }
    }
    
    if (check == T[d].size()) {
      return;
    }
  }
  
}
// ======================================================================================================
void preproc(vector <node *> init) {   // preprocessing needed for k-bct, gives init[0] for nodes at d=1,...,init[D-1] for d=D
  // no need to include root here as it is always in the data (obviously)
  
  
  init[D - 1]->c.push_back(zeros);       // for d=D , c=0 and lm[0]=0 from construction
  
  for (short d = D - 2; d > -1; d--) {
    
    init[d]->lm[0] = log2(beta);      // for smaller depth first add c=0 with p=logbeta
    init[d]->c.push_back(zeros);
    
    
    comb_initial3(d, init);           // then find all combinations and keep the top k of them
    
  }
  
  
}
// ======================================================================================================
vector<vector<double> > cartesian_prod(const vector<vector<double> >& v) { // cartesian product of vectors
  vector<vector<double> > s = { {} };                                  // returns a matrix with all combinations
  for (auto& u : v) {
    vector<vector<double> > r;
    for (auto& x : s) {
      for (auto y : u) {
        r.push_back(x);
        r.back().push_back(y);
      }
    }
    s.swap(r);
  }
  return s;
}
// ======================================================================================================
vector<vector<short> > cartesian_prod_int(const vector<vector<short> >& v) {   // finds cartesian product of vector of short ints
  vector<vector<short> > s = { {} };                                       // stores in matrix
  for (auto& u : v) {
    vector<vector<short> > r;
    for (auto& x : s) {
      for (auto y : u) {
        r.push_back(x);
        r.back().push_back(y);
      }
    }
    s.swap(r);
  }
  return s;
}
// ======================================================================================================
void comb(int d, int k, tree &T, vector <node *> init) {                   // finds combinations of children of node T[d][k]
  // computes the products for lm and keeps top k
  vector<vector<double> > v;  // stores combinations of probabilities
  vector<vector<short> > c;  // stores respective indices
  
  
  for (short ch = 0; ch < m; ch++) {
    
    
    node * pointer;
    
    if (T[d][k]->child[ch] == NULL) {
      
      pointer = init[d];                         // if child does not exists (not occured) use preprocessing node
      
    }
    
    else {
      
      pointer = T[d][k]->child[ch];                 // else use children of node
      
    }
    
    v.push_back(pointer->lm);
    
    vector <short> temp;
    
    for (unsigned short j = 1; j < pointer->lm.size() + 1; j++) { // use cartesian products of {1,2,3}... to find position vectors
      
      temp.push_back(j);
    }
    
    c.push_back(temp);
  }
  
  vector<vector<double> > v2 = cartesian_prod(v);          // cartesian products : all combinations
  c = cartesian_prod_int(c);
  
  for (unsigned int i = 0; i < v2.size(); i++) {                  //loop combiantions, keep what is necessary
    
    double sum = log2(1.0 - beta);
    
    for (unsigned short j = 0; j < v2[i].size(); j++) {
      sum = sum + v2[i][j];                          // sum of log-lm= prod-lm
    }
    
    vector <vector <short> > tempc;                    // use temporal variables to keep sort appropriately
    
    vector <double> temp;
    
    if (sum > T[d][k]->lm[T[d][k]->lm.size() - 1]) {  // if sum> min of ordered list lm
      
      int j = T[d][k]->lm.size() - 1;
      bool test = 0;
      
      while (j > 0) {
        
        if ((sum > T[d][k]->lm[j]) && (sum <= T[d][k]->lm[j - 1])) {  // if > of prev and < next, inlude there
          // stick element below the first one which is equal to sum
          temp.push_back(sum);
          tempc.push_back(c[i]);
          
          short limit = 0;
          
          if (T[d][k]->lm.size() < k_max) {
            limit = 1;                       // if the size of lm<k then add without replacement
          }
          
          for (unsigned int q = 0; q < T[d][k]->lm.size() - j - 1 + limit; q++) { //set temp
            temp.push_back(T[d][k]->lm[j + q]);
            tempc.push_back(T[d][k]->c[j + q]);
          }
          
          if (T[d][k]->lm.size() < k_max) {   // if size<k then add to list (here initialise space)
            T[d][k]->lm.push_back(0);
            T[d][k]->c.push_back({ 0 });
          }
          
          for (unsigned int q = 0; q < temp.size(); q++) {
            
            T[d][k]->lm[j + q] = temp[q];    // use temp for assignments
            T[d][k]->c[j + q] = tempc[q];
          }
          
          j = 1;            // used to break out of loop when added this combination
          test = 1;         // test used to add element at the top of the list
        }
        
        j--;
        
      }
      
      if (test == 0) {         // if sum> all lm then add it at the top of the list
        // code identical with before at index 0
        
        temp.push_back(sum);
        tempc.push_back(c[i]);
        
        short limit = 0;
        
        if (T[d][k]->lm.size() < k_max) {
          limit = 1;
        }
        
        for (unsigned int q = 0; q < T[d][k]->lm.size() - 1 + limit; q++) {
          temp.push_back(T[d][k]->lm[q]);
          tempc.push_back(T[d][k]->c[q]);
        }
        
        if (T[d][k]->lm.size() < k_max) {
          T[d][k]->lm.push_back(0);
          T[d][k]->c.push_back({ 0 });
        }
        
        for (unsigned int q = 0; q < temp.size(); q++) {
          
          T[d][k]->lm[q] = temp[q];
          T[d][k]->c[q] = tempc[q];
        }
        
      }
    }
    
    else { // if element <= last element of list, only include it at the bottom if size<k
      
      if (T[d][k]->lm.size() < k_max) {
        T[d][k]->lm.push_back(sum);
        T[d][k]->c.push_back(c[i]);
      }
    }
    
    
    
    
    
  }
  
  
}
// ======================================================================================================
void kbct_forw(tree &T, vector <node *> init) { // forward pass of kbct algorithm
  
  for (int d = D; d > -1; d--) {               // loop for leaves to root
    
    for (unsigned int k = 0; k < T[d].size(); k++) {
      
      if (d == D) {
        
        T[d][k]->lm[0] = T[d][k]->le;    // at leaves set lm[0]=le with position vector c=0
        T[d][k]->c.push_back(zeros);
        
      }
      
      
      else {
        
        T[d][k]->lm[0] = log2(beta) + T[d][k]->le; // for d<D first add the c=0 lm=le combination
        T[d][k]->c.push_back(zeros);               // if sth is equal to that it will be stuck below it
        // intuitevely keep c=0 higher to "reward" pruning at ties
        
        comb(d, k, T, init);                       // find combinations, sort and keep top k of them in list
        
      }
    }
  }
  
  
}
// ======================================================================================================
void kbct_back(vector <node *> init, tree T, vector <tree> &trees) {   // backward loop of kbct
  // builds top k trees using the lm and c's
  
  for (unsigned int i = 0; i < k_max; i++) { // loop over k-top trees
    if (T[0][0]->c[i][0] != 0) {  // if the c at the root is not zero, then don't prune there and add nodes at d=1
      
      node * temp = new node(m);  // create new node
      *temp = *T[0][0];        // initialise from T (so have children)
      trees[i][0][0] = temp;   // add to new tree
      
      for (int ch = 0; ch < m; ch++) { // always add m children to have proper tree output
        
        if (trees[i][0][0]->child[ch] != NULL) { // if child exists, initialise new node like the child and
          node * temp2 = new node(m);             // add it to the tree at depth 1
          *temp2 = *trees[i][0][0]->child[ch];
          trees[i][1].push_back(temp2);
          trees[i][0][0]->child[ch] = temp2; // connect it to the root of the new tree
        }
        
        else {  // if child doesn't exist, initialise new node from preprocessing step and add to tree, connect to root
          node * newnode = new node(m);
          *newnode = *init[0];
          trees[i][1].push_back(newnode);
          trees[i][0][0]->child[ch] = newnode;
        }
      }
      
      for (int d = 0; d < D - 1; d++) { // after adding the m nodes at d=1, check about pruning iteratively
        
        for (unsigned int k = 0; k < trees[i][d].size(); k++) {
          
          if (trees[i][d][k]->leaf == 0) {   // if node is not a leaf
            
            for (int j = 0; j < m; j++) {  // for all its children (all of them exist so next line probably not necessary)
              
              if (trees[i][d][k]->child[j]->leaf == 0) { // probably not necessary
                
                int index = 0;
                if (d == 0) { index = i; }
                short t = trees[i][d][k]->c[index][j] - 1; // index t as in notes (-1 because of c++ indexing from 0 and not from 1)
                
                
                // check if appropriate c is not 0, then no pruning and adding all children
                
                if (trees[i][d][k]->child[j]->c[t][0] != 0) {
                  
                  
                  trees[i][d][k]->child[j]->c[0] = trees[i][d][k]->child[j]->c[t]; // after node examined and children added, take t of next step from examined node
                  // (node will not be examined again so this is ok)
                  
                  for (int ch = 0; ch < m; ch++) {
                    
                    if (trees[i][d][k]->child[j]->child[ch] != NULL) {  // if child exists
                      
                      node * temp3 = new node(m);                        // create new child
                      *temp3 = *trees[i][d][k]->child[j]->child[ch];  // initialise from tree T
                      trees[i][d + 2].push_back(temp3);               // add it to tree
                      trees[i][d][k]->child[j]->child[ch] = temp3;    // connect it to parent
                    }
                    
                    else { //if child doesn't exist initialise from preprocessing and then as before
                      
                      node * newnode2 = new node(m);
                      *newnode2 = *init[d + 1];
                      if (d == D - 2) {
                        newnode2->leaf = 1; // the adding tree at depth D, so mark it a leaf
                      }
                      trees[i][d + 2].push_back(newnode2);
                      trees[i][d][k]->child[j]->child[ch] = newnode2;
                    }
                  }
                }
                
                else { // if pruning at node
                  
                  trees[i][d][k]->child[j]->leaf = 1; // mark it to be a leaf and make all children pointer =0 so they don't exist
                  for (int ch = 0; ch < m; ch++) { trees[i][d][k]->child[j]->child[ch] = 0; }
                  
                }
              }
            }
          }
        }
      }
    }
    
    else { // else if T[0][0] c[i]=0 then keep only root node
      
      node * newnode3 = new node(m);
      *newnode3 = *T[0][0];
      newnode3->leaf = 1;
      for (int ch = 0; ch < m; ch++) { newnode3->child[ch] = 0; }
      trees[i][0][0] = newnode3;
    }
  }
}
// ======================================================================================================
void kbct(tree& T, vector <tree> &trees, vector <node *> init, vector <Tree_properties> & tp_vec) { // call all kbct functions together
  
  long double pwl = ctw(T);
  
  kbct_forw(T, init);
  
  kbct_back(init, T, trees);
  
  
  
  for (unsigned int i = 0; i < k_max; i++) {
    
    Tree_properties tp;
    
    label(trees[i]);
    
    collect_leaves(trees[i], tp);
    
    long double prior = log2(pow(alpha, (tp.n_leaves - 1.0))*pow(beta, (tp.n_leaves - trees[i][D].size())));// log-prior
    
    
    tp. prior = pow(2, prior);
    tp.log_prior = prior;
    
    
    long double pml = T[0][0]->lm[i];
    
    long double posterior = pml - pwl; // log-posterior
    
    tp.posterior = pow(2, posterior);
    tp.log_posterior = posterior;
    
    tp.odd_posterior = pow(2, T[0][0]->lm[0] - pml);
    
    
    double bic = 0;
    double aic = 0;
    double mle = 0;
    compute_bic_aic_mle(trees[i], bic, aic, mle);
    tp.bic = bic;
    tp.aic = aic;
    tp.mle = mle;
    
    tp_vec.push_back(tp);
  }
}
// ======================================================================================================
void comb_initial3(int d, vector <node *> init) {    // finds combinations from preprocessing stage
  vector<vector<double> > v;  // stores combinations of probabilities
  vector<vector<short> > c;  // stores respective indices
  
  for (short ch = 0; ch < m; ch++) {
    
    node * pointer;
    pointer = init[d + 1];
    
    v.push_back(pointer->lm);
    
    vector <short> temp;
    
    for (unsigned short j = 1; j < pointer->lm.size() + 1; j++) {
      
      temp.push_back(j);
    }
    
    c.push_back(temp);
  }
  
  vector<vector<double> > v2 = cartesian_prod(v);
  c = cartesian_prod_int(c);
  
  for (unsigned int i = 0; i < v2.size(); i++) {
    
    double sum = log2(1.0 - beta);
    
    for (unsigned short j = 0; j < v2[i].size(); j++) {
      sum = sum + v2[i][j];
    }
    
    vector <vector <short> > tempc;
    
    vector <double> temp;
    
    if (sum > init[d]->lm[init[d]->lm.size() - 1]) {
      
      int j = init[d]->lm.size() - 1;
      bool test = 0;
      
      while (j > 0) {
        
        if ((sum > init[d]->lm[j]) && (sum <= init[d]->lm[j - 1])) {
          
          temp.push_back(sum);
          tempc.push_back(c[i]);
          
          short limit = 0;
          
          if (init[d]->lm.size() < k_max) {
            limit = 1;
          }
          
          for (unsigned int q = 0; q < init[d]->lm.size() - j - 1 + limit; q++) {
            temp.push_back(init[d]->lm[j + q]);
            tempc.push_back(init[d]->c[j + q]);
          }
          
          if (init[d]->lm.size() < k_max) {
            init[d]->lm.push_back(0);
            init[d]->c.push_back({ 0 });
          }
          
          for (unsigned int q = 0; q < temp.size(); q++) {
            
            init[d]->lm[j + q] = temp[q];
            init[d]->c[j + q] = tempc[q];
          }
          
          j = 1;
          test = 1;
        }
        
        j--;
        
      }
      
      if (test == 0) {
        
        temp.push_back(sum);
        tempc.push_back(c[i]);
        
        short limit = 0;
        
        if (init[d]->lm.size() < k_max) {
          limit = 1;
        }
        
        for (unsigned int q = 0; q < init[d]->lm.size() - 1 + limit; q++) {
          temp.push_back(init[d]->lm[q]);
          tempc.push_back(init[d]->c[q]);
        }
        
        if (init[d]->lm.size() < k_max) {
          init[d]->lm.push_back(0);
          init[d]->c.push_back({0});
        }
        
        for (unsigned int q = 0; q < temp.size(); q++) {
          
          init[d]->lm[q] = temp[q];
          init[d]->c[q] = tempc[q];
        }
        
      }
    }
    
    else {
      
      if (init[d]->lm.size() < k_max) {
        init[d]->lm.push_back(sum);
        init[d]->c.push_back(c[i]);
      }
    }
  }
}
// ======================================================================================================
Rcpp::NumericVector compute_log_loss(vector<short> xn, int train_size) {
  
  //This function computes the (cumulative) log-loss at each time-step in the test set.
  //The input xn here is the whole dataset and train_size is the size of the training set.
  
  Rcpp::NumericVector log_loss(xn.size() - train_size);   // Store the the log-loss at each time-step
  
  
  //1. First perform CTW for the training sequence
  
  tree T;
  init_tree(T);
  
  for (int i = D; i < train_size; i++) {
    
    short s = xn[i];          // current symbol
    vector <short> ct(D);     // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1];  // sets context
      
    }
    
    update(T, s, ct);              // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
  }
  
  long double init_ctw = ctw(T);     // Store the prior-predictive likelihood in the training set
  
  
  //2.Then evaluate the log-loss incurred by prediction in the test set by performing CTW sequentially
  
  
  for (unsigned int i = train_size; i < xn.size(); i++) {
    
    short s = xn[i];                           // current symbol
    vector <short> ct(D);                      // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1];     // sets context
      
    }
    
    update(T, s, ct);              // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
    
    // Here to perform CTW we only have to look at the D+1 contexts
    
    vector <node*> nodes_ct; //pointers to these nodes in Tmax (nodes already exist)
    
    nodes_ct.push_back(T[0][0]);
    
    for (int j = 1; j < D + 1; j++) {
      
      nodes_ct.push_back(nodes_ct[j - 1]->child[ct[j - 1]]); // this sets node pointers
      
    }
    for (int d = D; d > -1; d--) {           // loop over levels
      if (d == D) {                   // if at max depth, node is a leaf
        nodes_ct[d]->lw = nodes_ct[d]->le;
      }
      
      else {                         // if node is not a leaf
        
        long double sum = 0;
        
        for (int ch = 0; ch < m; ch++) {
          
          if (nodes_ct[d]->child[ch] != NULL) {           // if child #ch exists
            sum = sum + nodes_ct[d]->child[ch]->lw;     // calculate sum of Le(s)
            
          }
          
        }
        //calculate weighted log-prob in two cases (for numerical precision)
        
        long double delta = nodes_ct[d]->le - sum + log2(beta) - log2(1.0 - beta);
        if (delta < 30) {
          
          nodes_ct[d]->lw = log2(1.0 - beta) + sum + log2(1.0 + pow(2.0, delta));
          
        }
        else {
          nodes_ct[d]->lw = log2(beta) + nodes_ct[d]->le + log2(exp(1))*(pow(2.0, -delta) - pow(2.0, -2.0*delta - 1));
          
        }
        
      }
      
    }
    
    //END of sequential CTW
    
    log_loss[i- train_size] = (init_ctw - T[0][0]->lw); // this evaluates: -log P(x_{n+i}| x_1^n) = P(x_1^n)/P(x_1^{n+i)
    // i.e. the cumulative log-loss up to symbol i of the test set
  }
  for (int i = 0; i < log_loss.size(); i++) { log_loss[i] = log_loss[i] * log(2.0)/(i+1);
  } // convert log2 to ln
  
  return log_loss;
}
// ======================================================================================================
void set_param(string &s, int depth, short kmax){ // with kmax; it is used in the kbct function
  
  D = depth;
  k_max = kmax;
  
  xn.clear();
  zeros.clear();
  encoder.clear();
  decoder.clear();
  for (char c : s) {
    if(encoder.find(c) == encoder.end()){  // e.g. encodes the first enctountered character as 0, the second which is different from the first one to 1, etc...
      encoder.insert({c,encoder.size()});  
    }
    xn.push_back(encoder[c]);
  }
  for (map<char, short>::iterator i = encoder.begin(); i != encoder.end(); ++i)   // contructs the decoder
    decoder[i->second] = i->first;
  
  m = encoder.size(); // alphabet size
  for(int i=0;i<m;i++){
    zeros.push_back(0);
  }
  
}
// ======================================================================================================
void set_param(string &s, int depth){ // just as before, but without the kmax parameter
  
  D = depth;
  
  xn.clear();
  zeros.clear();
  encoder.clear();
  decoder.clear();
  for (char c : s) {
    if(encoder.find(c) == encoder.end()){
      encoder.insert({c,encoder.size()});  
    }
    xn.push_back(encoder[c]);
  }
  for (map<char, short>::iterator i = encoder.begin(); i != encoder.end(); ++i)
    decoder[i->second] = i->first;
  
  m = encoder.size();
  for(int i=0;i<m;i++){
    zeros.push_back(0);
  }
  
}
// ======================================================================================================
void set_global_parameters(string &s, int depth, short kmax){ // after reading and encoding the input dataset, the global parameters are fixed 
  set_param(s, depth, kmax);
  beta = 1 - pow(2, -(m - 1)); // default value for the prior hyper-parameter
  alpha = pow((1.0 - beta), (1.0 / (m - 1.0)));
}
// ======================================================================================================
void set_global_parameters(string &s, int depth, short kmax, double b){
  set_param(s, depth, kmax);
  if(b>0 && b<1)
    beta = b;   // for a custom beta.
  else
    beta = 1 - pow(2, -(m - 1));
  alpha = pow((1.0 - beta), (1.0 / (m - 1.0)));
}
// ======================================================================================================
vector <Tree_properties> build_kbct() {
  
  
  vector <Tree_properties> tp_vec;
  
  //1. Initialise tree to store Tmax
  tree T;
  init_tree(T);
  vector <tree> trees(k_max, T); //initialise top-k trees
  
  
  
  //2. Find pre-processing nodes at each depth
  
  vector <node *> init; // pointers for nodes of pre-processing stage (not needed  for root node)
  
  if (D > 0) {
    
    for (short d = 0; d < D; d++) {
      init.push_back(new node(m)); // initialise them
    }
    preproc(init);
    
  }
  
  //3. Update for each sequence symbol to build Tmax and calculate estimated probabilities
  
  for (unsigned int i = D; i < xn.size(); i++) {
    
    short s = xn[i]; // current symbol
    vector <short> ct(D); // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1]; // sets context
      
    }
    
    update(T, s, ct); // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
  }
  
  //  cout << "Tree was built" << endl;
  
  
  //4. Run CTW and k-BCT
  
  vector<double> odds(k_max, 0.0); // stores posterior odds for top-k trees
  kbct(T, trees, init, tp_vec); // runs CTW followed by k-BCT
  
  // cleaning process: deleting the Tmax from the memory
  for (int d = T.size() - 1; d >= 0; d--) {
    for (unsigned int k = 0; k < T[d].size(); k++) {
      delete [] T[d][k]->child;
      delete [] T[d][k]->a;
      delete T[d][k];
    }
  }
  // cleaning process: deleting the all modelsfrom the memory
  for (unsigned int i=0; i<k_max;i++){
    for (int d = trees[i].size() - 1; d >= 0; d--) {
      for (unsigned int k = 0; k < trees[i][d].size(); k++) {
        delete [] trees[i][d][k]->child;
        delete [] trees[i][d][k]->a;
        delete trees[i][d][k];
      }
    }
  }
  return tp_vec;
  
}
// ======================================================================================================
void compute_bic_aic_mle(tree &T, double &bic, double &aic, double &mle) { 
  // takes a tree and computes the aic, bic scores along with the maximum likelihood of the model
  
  int n_leaves = show_leaves(T);
  
  counts(T);
  mle = compute_mle(T);
  bic = -2 * mle + (n_leaves * (m - 1))* log(xn.size() - D);
  aic = -2 * mle + (n_leaves * (m - 1)) * 2;
  
  
}
// ======================================================================================================
void counts(tree &T) { // finds counts of proper tree T in sequence xn
  
  for (int d = 0; d < D + 1; d++) {  // first reset all counts to zero 
    
    for (unsigned int k = 0; k < T[d].size(); k++) {
      
      for (int ch = 0; ch < m; ch++) {
        
        T[d][k]->a[ch] = 0;
      }
    }
  }
  
  int max;
  for (int i = 0; i < D + 1; i++) {
    if (T[i].size() > 0) {
      max = i;
    }
  }
  max = D;  // use D or max depth for initial value, comment this line out to use max depth
  
  // update for each sequence symbol
  for (unsigned int i = max; i < xn.size(); i++) {
    
    
    short s = xn[i];          // current symbol
    vector <short> ct(D);     // current context
    
    
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1];  // sets context
      
    }
    
    node * temp = T[0][0];
    T[0][0]->a[s]++;
    
    for (int j = 0; j < D; j++) {
      
      if (temp->leaf == 0) {
        
        temp = temp->child[ct[j]];   //proper tree so if not a leaf all children exist
        temp->a[s]++;
      }
      
      else {
        j = D + 5;
      }
    }
  }
}
// ======================================================================================================
double compute_mle(tree &T) {  //finds mle of proper tree T 
  
  double sum = 0;
  
  for (int d = 0; d < D + 1; d++) {
    
    for (unsigned int k = 0; k < T[d].size(); k++) { // loop over all leaves of the tree
      
      if (T[d][k]->leaf == 1) {
        
        int M = 0;
        
        for (int j = 0; j < m; j++) {
          M = M + T[d][k]->a[j];
        }
        
        for (int j = 0; j < m; j++) {
          
          if (T[d][k]->a[j] != 0) {
            
            sum = sum + 1.0 * T[d][k]->a[j] * log(1.0*T[d][k]->a[j] / M);
          }
          
        }
        
      }
    }
    
  }
  return sum;
}
// ======================================================================================================
void ctw_bct(tree& T, Tree_properties & tp) {
  
  // finds the MAP tree and computes exactly the posterior
  long double pwl = ctw(T);
  
  long double pml = bct(T);
  
  tree T2;
  init_tree(T2);
  T2[0].pop_back();
  copytree(T, T2);
  
  label(T2); //add contexts to nodes
  
  collect_leaves(T2, tp);
  long double prior = log2(pow(alpha, (tp.n_leaves - 1.0)) * pow(beta, (tp.n_leaves - T2[D].size()))); // log-prior
  
  tp.prior = pow(2, prior);
  tp.log_prior = prior;
  
  long double posterior = pml - pwl; // log-posterior
  
  tp.posterior = pow(2, posterior);
  tp.log_posterior = posterior;
  
  double bic = 0;
  double aic = 0;
  double mle = 0;
  compute_bic_aic_mle(T2, bic, aic, mle);
  tp.bic = bic;
  tp.aic = aic;
  tp.mle = mle;
  tp.odd_posterior = 0;
}
// ======================================================================================================
Tree_properties build_bct() {
  
  Tree_properties tp;
  
  // Initialise tree to store Tmax
  tree T;
  init_tree(T);
  
  
  // Update for each sequence symbol to build Tmax and calculate estimated probabilities
  
  for (unsigned int i = D; i < xn.size(); i++) {
    
    short s = xn[i]; // current symbol
    vector <short> ct(D); // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1]; // sets context
      
    }
    
    update(T, s, ct); // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
  }
  
  
  
  ctw_bct(T, tp);
  
  
  // Releases the memory occupied by the tree 
  for (int d = T.size() - 1; d >= 0; d--) {
    for (unsigned int k = 0; k < T[d].size(); k++) {
      delete [] T[d][k]->child;
      delete [] T[d][k]->a;
      delete T[d][k];
    }
  }
  
  // Run CTW and BCT
  
  return tp;
  
}
// ======================================================================================================
void copy(tree &Tout, tree &Tin) { // takes as input a tree and gives an identical copy of the tree (so that can change the 2nd one on its own)
  // also orders nodes in correct order in tree
  // works for both proper and improper trees, and trees with contexts or without contexts
  
  
  for (int i = 0; i < D; i++) {
    
    for (unsigned int j = 0; j < Tout[i].size(); j++) {
      
      for (int ch = 0; ch < m; ch++) {
        
        if (Tout[i][j]->child[ch] != NULL) {
          
          node * temp = new node(m);
          Tout[i + 1].push_back(temp);
          *temp = *Tout[i][j]->child[ch];
          Tout[i][j]->child[ch] = temp;
        }
      }
    }
  }
  
}
// ======================================================================================================
Rcpp::List online_predict( int train_size) {
  
  //The input xn here is the whole dataset and train_size is the size of the training set, i.e. test set size is xn.size()-train_size
  //This function performs sequential prediction
  
  Rcpp::CharacterVector prediction;
  //1. First perform CTW for the training sequence
  
  tree T;
  init_tree(T);
  
  // constructs the alphabet
  vector <char> alphabet;
  
  for (map<char, short>::iterator i = encoder.begin(); i != encoder.end(); ++i){
    alphabet.push_back(i->first);
  }
  //sorts the alphabet in alphabetical order
  sort(alphabet.begin(), alphabet.end()); 
  
  
  
  for (int i = D; i < train_size; i++) {
    
    short s = xn[i];          // current symbol
    vector <short> ct(D);     // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1];  // sets context
      
    }
    
    update(T, s, ct);              // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
  }
  
  long double init_ctw = ctw(T);     // Store the prior-predictive likelihood in the training set
  long double temp_ctw = init_ctw;
  
  //2.Then try m different values to find conditional distrution of xn+1/xn
  
  Rcpp::List total_probs;
  for (unsigned int i = train_size; i < xn.size(); i++) {
    
    Rcpp::NumericVector probs(m);
    int max_i = 0; long double max_p = 0.0;
    
    for (int t = 0; t < m; t++) { //t is the tested value of xn+1
      node new_node(m);
      vector <node*> row(1);
      row[0] = &new_node;
      
      tree Temp(D + 1, row);
      init_tree(Temp);
      
      Temp[0][0] = new node(m);
      *Temp[0][0] = *T[0][0];
      
      copy(Temp, T);
      
      short s = t;                                // current symbol
      vector <short> ct(D);                      // current context
      
      for (int j = 0; j < D; j++) {
        
        ct[j] = xn[i - j - 1];     // sets context
        
      }
      
      update(Temp, s, ct);              // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
      
      // Here to perform CTW we only have to look at the D+1 contexts
      
      vector <node*> nodes_ct; //pointers to these nodes in Tmax (nodes already exist)
      
      nodes_ct.push_back(Temp[0][0]);
      
      for (int j = 1; j < D + 1; j++) {
        
        nodes_ct.push_back(nodes_ct[j - 1]->child[ct[j - 1]]); // this sets node pointers
        
      }
      
      
      for (int d = D; d > -1; d--) {           // loop over levels
        
        
        if (d == D) {                   // if at max depth, node is a leaf
          nodes_ct[d]->lw = nodes_ct[d]->le;
        }
        
        else {                         // if node is not a leaf
          
          long double sum = 0;
          
          for (int ch = 0; ch < m; ch++) {
            
            if (nodes_ct[d]->child[ch] != NULL) {           // if child #ch exists
              sum = sum + nodes_ct[d]->child[ch]->lw;     // calculate sum of Le(s)
              
            }
            
          }
          
          //calculate weighted log-prob in two cases (for numerical precision)
          
          long double delta = nodes_ct[d]->le - sum + log2(beta) - log2(1.0 - beta);
          if (delta < 30) {
            
            nodes_ct[d]->lw = log2(1.0 - beta) + sum + log2(1.0 + pow(2.0, delta));
            
          }
          else {
            nodes_ct[d]->lw = log2(beta) + nodes_ct[d]->le + log2(exp(1))*(pow(2.0, -delta) - pow(2.0, -2.0*delta - 1));
          }
        }
      }

      long double prob = pow(2, Temp[0][0]->lw - temp_ctw);
      probs[encoder[alphabet[t]]] = prob;
      
      
      //END of sequential CTW
      // releases the memory
      for (int d = Temp.size() - 1; d >= 0; d--) {
        for (unsigned int k = 0; k < Temp[d].size(); k++) {
          delete [] Temp[d][k]->child;
          delete [] Temp[d][k]->a;
          delete Temp[d][k];
        }
      }
      
      
      
      if (prob > max_p) {
        max_i = t; max_p = prob;
      }
      
    }
    
    short s = xn[i];                           // current symbol
    vector <short> ct(D);                      // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1];     // sets context
      
    }
    
    update(T, s, ct);              // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
    
    // Here to perform CTW we only have to look at the D+1 contexts
    
    vector <node*> nodes_ct; //pointers to these nodes in Tmax (nodes already exist)
    
    nodes_ct.push_back(T[0][0]);
    
    for (int j = 1; j < D + 1; j++) {
      
      nodes_ct.push_back(nodes_ct[j - 1]->child[ct[j - 1]]); // this sets node pointers
      
    }
    
    
    for (int d = D; d > -1; d--) {           // loop over levels
      
      
      if (d == D) {                   // if at max depth, node is a leaf
        nodes_ct[d]->lw = nodes_ct[d]->le;
      }
      
      else {                         // if node is not a leaf
        
        long double sum = 0;
        
        for (int ch = 0; ch < m; ch++) {
          
          if (nodes_ct[d]->child[ch] != NULL) {           // if child #ch exists
            sum = sum + nodes_ct[d]->child[ch]->lw;     // calculate sum of Le(s)
            
          }
          
        }
        
        //calculate weighted log-prob in two cases (for numerical precision)
        
        long double delta = nodes_ct[d]->le - sum + log2(beta) - log2(1.0 - beta);
        if (delta < 30) {
          
          nodes_ct[d]->lw = log2(1.0 - beta) + sum + log2(1.0 + pow(2.0, delta));
          
        }
        else {
          nodes_ct[d]->lw = log2(beta) + nodes_ct[d]->le + log2(exp(1))*(pow(2.0, -delta) - pow(2.0, -2.0*delta - 1));
        }
      }
      
    }
    
    //END of sequential CTW
    temp_ctw = T[0][0]->lw;
    prediction.push_back(decoder[max_i]); // outputs the most likely character
    total_probs.push_back(probs); // append to the total_probs the distribution
    
  }
  // releases memory
  for (int d = T.size() - 1; d >= 0; d--) {
    for (unsigned int k = 0; k < T[d].size(); k++) {
      delete [] T[d][k]->child;
      delete [] T[d][k]->a;
      delete T[d][k];
    }
  }
  total_probs.push_back(prediction, "Prediction");
  //Rcpp::List output = List::create(Named("Predictions") = prediction, _["Probability_vectors"] = total_probs);
  
  return total_probs;
}
// ======================================================================================================
long double mle_tree() { // Computes the likelihood of the Maximum likelihood model
  
  tree T;
  init_tree(T);
  //theta is map that stores parameters for each context
  
  Rcpp::List theta;
  //1. Build improper Tmax which is the MLE tree
  for (unsigned int i = D; i < xn.size(); i++) {
    
    short s = xn[i];          // current symbol
    vector <short> ct(D);     // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1];  // sets context
      
    }
    
    update(T, s, ct);              // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
  }
  
  
  
  //2.Associate parameter vectors to leaves of Tmax based on counts and make tree proper
  int length_D = T[D].size();// initial nodes of depth D that occurred in xn
  
  Rcpp::NumericVector theta_s; 
  for (int i = 0; i < m; i++) { theta_s.push_back(0); }
  
  
  
  
  for (int d = 0; d < D; d++) { //for nodes at depth <D
    
    for (unsigned int k = 0; k < T[d].size(); k++) {
      
      if (T[d][k]->leaf == 0) { // i.e. for all internal nodes
        
        for (short ch = 0; ch < m; ch++) {
          
          if (T[d][k]->child[ch] != NULL) { // if a child does exist label it
            
            T[d][k]->child[ch]->s = T[d][k]->s;
            T[d][k]->child[ch]->s.push_back(ch);
            
          }
          
          else { //if it doesn't exist, add it and label it
            
            node * init = new node(m);
            T[d][k]->child[ch] = init;             // connect child with parent node
            T[d + 1].push_back(init);              // store at appropriate tree depth
            init->leaf = 1;                        // denote it leaf (then it's not considered next)
            init->s = T[d][k]->s;                  // adds its context
            init->s.push_back(ch);
            
            //and take distirbution from its parent
            
            float sum = 0; // for M_s
            for (int j = 0; j < m; j++) {
              
              theta_s[j] = 1.0 * T[d][k]->a[j];
              sum = sum + 1.0 * T[d][k]->a[j]; //this is M_s at the end of loop
              
            }
            
            for (int j = 0; j < m; j++) { //normalise
              theta_s[j] = 1.0 * theta_s[j] / sum; // these are estimated parameters
            }
            
            //theta[init->s] = theta_s;  // set to theta dictionary for children
            
            Rcpp::String chv;
            
            for (unsigned int j = 0; j<(init->s).size(); j++)
              chv.push_back(decoder[init->s[j]]);
            
            theta.push_back(theta_s, chv);
            
          }
          
        }
      }
    }
  }
  
  
  
  //for nodes at depth D (only need for those that have occurred in xn, others have 0 count vector)
  
  long double mle = 0;
  
  for (int k = 0; k < length_D; k++) {
    
    
    float sum = 0;
    for (int j = 0; j < m; j++) {
      
      // theta_s[j] = 1.0 * T[D][k]->a[j];
      sum = sum + 1.0 * T[D][k]->a[j];
      
    }
    
    for (int j = 0; j < m; j++) { //find mle
      if (T[D][k]->a[j] > 0) {
        mle = mle + 1.0 * T[D][k]->a[j] * (log(1.0 * T[D][k]->a[j]) - log(sum));
      }
      
    }
    Rcpp::String chv;
    
  }
  
  Rcpp::Rcout << "log-ML is: "<<endl;
  for (int d = T.size() - 1; d >= 0; d--) {
    for (unsigned int k = 0; k < T[d].size(); k++) {
      delete [] T[d][k]->child;
      delete [] T[d][k]->a;
      delete T[d][k];
    }
  }
  return mle;
  
}
// ======================================================================================================
string vec2str(vector <short> &s){
  
  string chv(s.size(), '0');
  for (unsigned int i = 0; i<s.size(); i++){
    chv[i] = decoder[s[i]];
  }
  return chv;
}
// ======================================================================================================
void makeproper(tree& T) { // takes improper tree and makes it proper
  
  
  for (int d = 0; d < D + 1; d++) {
    for (unsigned int k = 0; k < T[d].size(); k++) {
      
      if (T[d][k]->leaf == 0) {
        
        for (int ch = 0; ch < m; ch++) {
          
          if (T[d][k]->child[ch] == NULL) {
            
            node * init = new node(m);
            T[d][k]->child[ch] = init;             // connect child with parent node
            T[d + 1].push_back(init);              // store at appropriate tree depth
            init->leaf = 1;                        // denote it leaf
            init->s = T[d][k]->s;                  // if contexts used
            init->s.push_back(ch);
          }
        }
      }
    }
  }
  
}
// ======================================================================================================
long double build_ctw_rcpp() { 
  
  // Initialise tree to store Tmax
  tree T;
  init_tree(T);
  
  
  // Update for each sequence symbol to build Tmax and calculate estimated probabilities
  
  for (unsigned int i = D; i < xn.size(); i++) {
    
    short s = xn[i]; // current symbol
    vector <short> ct(D); // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1]; // sets context
      
    }
    
    update(T, s, ct); // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
  }
  
  
  long double prob = ctw(T);
  
  for (int d = T.size() - 1; d >= 0; d--) {
    for (unsigned int k = 0; k < T[d].size(); k++) {
      delete [] T[d][k]->child;
      delete [] T[d][k]->a;
      delete T[d][k];
    }
  }
  // Run CTW and BCT
  
  return prob;
  
}
// ======================================================================================================
map<string, vector<int>> dictionary_counts() { // construct a dictionary to look for existence of nodes in Tmax
  
  vector <char> alphabet;
  
  // construct the alphabet and order it; 
  for (map<char, short>::iterator i = encoder.begin(); i != encoder.end(); ++i){
    alphabet.push_back(i->first);
  }
  
  sort(alphabet.begin(), alphabet.end()); 
  
  tree Tmax;
  init_tree(Tmax);
  
  for (unsigned int i = D; i < xn.size(); i++) {
    
    short s = xn[i]; // current symbol
    vector <short> ct(D); // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1]; // sets context
      
    }
    
    update(Tmax, s, ct); // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
  }
  
 // makeproper(Tmax);
  label_inproper(Tmax);  
  
  // using only one operation, by looking at [context] element
  map <string, vector<int> > dict;
  vector<int> vec(m);
  for (int i = 0; i < m; i++) {
    vec[encoder[alphabet[i]]] = (Tmax[0][0]->a[i]);// no need to encode root node as it always exists, but can have empty context
  }
  
  dict["Root"] = vec;
  
  for (int d = 1; d < D + 1; d++) {
    
    for (unsigned int k = 0; k < Tmax[d].size(); k++) {
      vector<int> vec2(m);
      for (int i = 0; i < m; i++) {
        vec2[encoder[alphabet[i]]] = Tmax[d][k]->a[i];// no need to encode root node as it always exists, but can have empty context
      }
      string chv = vec2str(Tmax[d][k]->s);
      dict[chv] = vec2;
    }
  }
  // releases memory
  for (int d = Tmax.size() - 1; d >= 0; d--) {
    for (unsigned int k = 0; k < Tmax[d].size(); k++) {
      delete [] Tmax[d][k]->child;
      delete [] Tmax[d][k]->a;
      delete Tmax[d][k];
    }
  }
  
  return(dict);
}
// ======================================================================================================
Rcpp::List map_param() {
  // calculates the parameters of each leaf within the MAP model
  
  // Initialise tree to store Tmax
  tree T;
  init_tree(T);
  
  // Update for each sequence symbol to build Tmax and calculate estimated probabilities
  
  for (unsigned int i = D; i < xn.size(); i++) {
    
    short s = xn[i]; // current symbol
    vector <short> ct(D); // current context
    
    for (int j = 0; j < D; j++) {
      
      ct[j] = xn[i - j - 1]; // sets context
      
    }
    
    update(T, s, ct); // Updates sequentially Tmax when symbol "s" follows context "ct" in xn
  }
  
  
  vector <char> alphabet;
  
  for (map<char, short>::iterator i = encoder.begin(); i != encoder.end(); ++i){
    alphabet.push_back(i->first);
  }
  
  sort(alphabet.begin(), alphabet.end()); 
  
  
  Rcpp::List theta;
  
  long double pml = bct(T);
  pml = pml+0;
  
  tree T2;
  init_tree(T2);
  T2[0].pop_back();
  copytree(T, T2); // in T2 is stored the MAP tree
  label(T2); //add contexts to nodes
  
  for (int d = 0; d < D + 1; d++) {
    for (unsigned int k = 0; k < T2[d].size(); k++) {
      if (T2[d][k]->leaf == 1) {
        
        double sum = 0;
        Rcpp::NumericVector param(m);
        
        for (int j = 0; j < m; j++) {
          param[encoder[alphabet[j]]] = (1.0 * T2[d][k]->a[j] + 0.5);
          sum = sum + 1.0 * T2[d][k]->a[j] + 0.5;
        }
        for (int j = 0; j < m; j++) { //normalize
          
          param[j] = param[j] / sum;
          
        }
        
        string chv = vec2str(T2[d][k]->s);
        theta.push_back(param, chv);
      }
    }
  }
  //releases memory
  for (int d = T.size() - 1; d >= 0; d--) {
    for (unsigned int k = 0; k < T[d].size(); k++) {
      delete [] T[d][k]->child;
      delete [] T[d][k]->a;
      delete T[d][k];
    }
  }
  
  return(theta);
}
// ======================================================================================================

#endif // UTILS_H
