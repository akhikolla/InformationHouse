#include <Rcpp.h>

#define TRUE 1
#define FALSE 0

#define MAX(a,b) (((a)<(b))?(b):(a))

using namespace Rcpp;

//' Find the shortest paths to other vertices
//' 
//' @param dist a matrix of distances between nodes
//' @param src an integer vertex ID 
//' @param node_costs a list of costs, in order, of all nodes represented in the 
//'     sociomatrix, all are assumed 0 if unspecified
//' @return a numeric vector, entry \emph{i} of which is the vertex immediately preceeding 
//'     vertex \emph{i} in the shortest path leading to that vertex. Full paths must be constructed
//'     recursively.
// [[Rcpp::export]]
IntegerVector dijkstra_nodes(NumericMatrix dist, int src, NumericVector node_costs){
  /* Important Note: Indices will all be shifted down by one
  * from their usual R values, so dist(0,2) (this is a method
  * call to the NumericMatrix class accessing index [0][2])
  * will return the connection from vertex 1 to vertex 3.
  */
  src = src - 1; // Shift down R-style input
  
  int nv = dist.nrow();
  IntegerVector prev(nv,NA_INTEGER);
  int on;
  int some_more_unvisited;
  int up_next;
  std::vector<int> visited(nv,FALSE);
  double src_to_on;
  double cost_of_on;
  double on_to_i;
  double dist_to_beat;
  std::vector<double> short_src_to_(nv, R_PosInf);
  
  short_src_to_[src] = 0;
  on = src; // Move to src vertex
  
  do{ //make sure this is robust to disconnected graphs
    Rcpp::checkUserInterrupt();
    visited[on] = TRUE; // Mark current vertex as visited
    
    src_to_on = short_src_to_[on];
    
    some_more_unvisited = FALSE; // Wishful Thinking 
    
    dist_to_beat = R_PosInf; // Haven't Checked distances to other vertices yet
    up_next = NA_INTEGER; // Not Yet Known 
    
    for(int i = 0; i < nv; i++){ // Iterate through Unvisited Vertices and consider paths through the one we're On
      if(!visited[i]){
        some_more_unvisited = TRUE; // Dang... The loop is evaluated at least once - There will be more vertices to visit after the one we're On
        on_to_i = dist(on,i);
        if(on == src){cost_of_on = 0;}else{cost_of_on = node_costs(on);}
        if(src_to_on + cost_of_on + on_to_i < short_src_to_[i]){ // Is path through "On" to i any shorter?
          short_src_to_[i] = src_to_on + cost_of_on + on_to_i;
          prev[i] = on;
        } 
        if(short_src_to_[i] < dist_to_beat){ // Keep Track of Closest Vertex to get Next Step.
          dist_to_beat = short_src_to_[i];
          up_next = i;
        }
      }
    }
    on = up_next;
  }while(some_more_unvisited &&  (up_next >= 0) ); // NA_INTEGER is stored as smallest int
  
  // Shift back to R-style indices
  for(int i = 0; i < nv; i++){
    if(i != src) prev[i]++;
  }
  
  return prev;
}

//' Find the shortest L-Inf norm paths to other vertices
//' 
//' @param dist a matrix of distances between nodes
//' @param src an integer vertex ID 
//' @return a numeric vector, entry \emph{i} of which is the vertex immediately preceeding 
//'   vertex \emph{i} in the shortest path leading to that vertex. Full paths must be constructed
//'   recursively.
// [[Rcpp::export]]
IntegerVector dijkstra_inf(NumericMatrix dist, int src){
  /* Important Note: Indices will all be shifted down by one
  * from their usual R values, so dist(0,2) (this is a method
  * call to the NumericMatrix class accessing index [0][2])
  * will return the connection from vertex 1 to vertex 3.
  */
  src = src - 1; // shift down R-style input
  
  int nv = dist.nrow();
  IntegerVector prev(nv,NA_INTEGER);
  int on;
  int some_more_unvisited;
  int up_next;
  std::vector<int> visited(nv,FALSE);
  double src_to_on;
  double on_to_i;
  double dist_to_beat;
  std::vector<double> short_src_to_(nv, R_PosInf);
  
  short_src_to_[src] = 0;
  on = src; // Move to src vertex
  
  do{ //make sure this is robust to disconnected graphs
    Rcpp::checkUserInterrupt();
    visited[on] = TRUE; // Mark current vertex as visited
    
    src_to_on = short_src_to_[on];
    
    some_more_unvisited = FALSE; // Wishful Thinking 
    
    dist_to_beat = R_PosInf; // Haven't Checked distances to other vertices yet
    up_next = NA_INTEGER; // Not Yet Known 
    
    for(int i = 0; i < nv; i++){ // Iterate through Unvisited Vertices and consider paths through the one we're On
      if(!visited[i]){
        some_more_unvisited = TRUE; // Dang... The loop is evaluated at least once - There will be more vertices to visit after the one we're On
        on_to_i = dist(on,i);
        if(MAX(src_to_on, on_to_i) < short_src_to_[i]){ // Is path through "On" to i any shorter? && Make sure it exists (NA is -1)
          short_src_to_[i] = MAX(src_to_on, on_to_i);
          prev[i] = on;
        } 
        if(short_src_to_[i] < dist_to_beat){ // Keep Track of Closest Vertex to get Next Step.
          dist_to_beat = short_src_to_[i];
          up_next = i;
        }
      }
    }
    on = up_next;
  }while(some_more_unvisited &&  (up_next >= 0) ); // NA_INTEGER is stored as smallest int
  
  // Shift back to R-style indices
  for(int i = 0; i < nv; i++){
    if(i != src) prev[i]++;
  }
  return prev;
}
