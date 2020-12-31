#ifndef NODE_H
#define NODE_H
#include <vector>
#include <algorithm>
using namespace std;

typedef vector<vector<short> > matrix;

class node {
public:
  
  int m; // Alphabet size
  vector<short>  s; // Node context
  int * a; // Occurrences for each j in alphabet
  
  double le; // Logarithm of estimated probability Pe(as)
  double lw; // Logarithm of weighted probability probability
  vector <double> lm; // List  of log-maximal probabilities
  matrix c; // Position vectors for k-bct
  
  bool leaf; // Indicates if node is a leaf
  node ** child; // Pointers for children nodes
  node(int alphabetsize);
  node& operator=(const node& other);
  
};

typedef node* nodeptr;



node::node(int alphabetsize) {
  m = alphabetsize;
  le = lw = 0;
  lm.push_back(0.0);
  leaf = 0;
  
  a = new int[alphabetsize];
  memset(a, 0, sizeof(int)*alphabetsize);
  child = new nodeptr[alphabetsize];
  memset(child, 0, sizeof(nodeptr)*alphabetsize);
  
  
}

node& node::operator=(const node& other) // copy assignment
{
  if (this != &other) { 
    this->s=other.s;
    this->le=other.le;
    this->lw=other.lw;
    this->lm=other.lm;
    this->c = other.c;
    this->leaf = other.leaf;
    std::copy(other.a, other.a + other.m, this->a);
    std::copy(other.child, other.child + other.m, this->child);
  }
  return *this;
}
// Introduce a structure for a tree: group of node pointers
//      description: tree[m] stores pointers of nodes at depth m
typedef vector <vector <node *> > tree;
#endif // NODE_H
