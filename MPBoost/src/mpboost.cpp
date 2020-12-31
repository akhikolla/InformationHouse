#include <Rcpp.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdlib>

// Adding information to vertices ( 3 new properties):
// level
struct vertex_info1_t {
  typedef boost::vertex_property_tag kind;
};

// n2
struct vertex_info2_t {
  typedef boost::vertex_property_tag kind;

};

// V(.,.)
struct vertex_info3_t {
  typedef boost::vertex_property_tag kind;
};

typedef boost::property<vertex_info1_t, int, boost::property<vertex_info2_t, int, boost::property<vertex_info3_t, boost::multiprecision::cpp_int>>> VertexProperties;

typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, VertexProperties, EdgeWeightProperty> DirectedGraph;

typedef boost::graph_traits<DirectedGraph>::edge_iterator edge_iter;

typedef boost::graph_traits<DirectedGraph>::out_edge_iterator out_edge_iterator;

typedef boost::graph_traits<DirectedGraph>::vertex_iterator vertex_iter;

// [[Rcpp::export()]]
Rcpp::IntegerVector mpboost(int N1, int N2, int MTI) {
  int N = N1 + N2;
  boost::multiprecision::cpp_int temp;
  DirectedGraph g((N + 1)*(N + 2)/2);
  boost::property_map<DirectedGraph, vertex_info1_t>::type VertexInfo1Map = get(vertex_info1_t(), g);
  boost::property_map<DirectedGraph, vertex_info2_t>::type VertexInfo2Map = get(vertex_info2_t(), g);
  boost::property_map<DirectedGraph, vertex_info3_t>::type VertexInfo3Map = get(vertex_info3_t(), g);
  boost::property_map<DirectedGraph, boost::edge_weight_t>::type EdgeWeightMap = get(boost::edge_weight_t(), g);
  out_edge_iterator oei, oei_end;

  for (int i = 0; i != N; ++i) {
    for (int j = i*(i + 1)/2; j != i*(i + 1)/2 + i + 1; ++j) {
      if (fabs(i - static_cast<double>(N)/N2*(j - i*(i + 1)/2)) <= MTI) {
	add_edge(j, j + i + 1, 0, g); 
	add_edge(j, j + i + 2, 0, g);
	put(VertexInfo1Map, j, i); // Info 1: level
	put(VertexInfo2Map, j, j - i*(i + 1)/2); // Info 2: number of n2 at each level
      }
    }
  }
  // Upper level (i=N)
  int vertn = (N + 1)*(N + 2)/2 - (N + 1) + N2;
  put(VertexInfo1Map, vertn, N);
  put(VertexInfo2Map, vertn, N2);
  // initialize V(N,N2)=1
  put(VertexInfo3Map, vertn, 1);

  // Compute V(.,.) in lower levels
  for (int i = N - 1; i != -1; --i) {
    for (int j = i*(i + 1)/2; j != i*(i + 1)/2 + i + 1; ++j) {
      temp = 0;
      for (tie(oei, oei_end) = out_edges(j, g); oei != oei_end; ++oei)
	temp += VertexInfo3Map[target(*oei, g)];
      put(VertexInfo3Map, j, temp);
    }
  }
  
  // Compute weights
  for (std::pair<edge_iter, edge_iter> edgePair = edges(g); edgePair.first != edgePair.second; ++edgePair.first) {
    if (VertexInfo3Map[source(*edgePair.first, g)] > 0) {  // This condition is to avoid nan's coming from 0/0...
      boost::multiprecision::cpp_rational tempr(VertexInfo3Map[target(*edgePair.first, g)], VertexInfo3Map[source(*edgePair.first, g)]);
      put(EdgeWeightMap, *edgePair.first, tempr.convert_to<double>());
    }
  }
  
  // Vertex information
  std::pair<vertex_iter, vertex_iter> vertexPair;
  // Generate sequence of x_i taking the weigths into account and starting from vertex (0,0)
  Rcpp::IntegerVector x;
  double bern;
  unsigned int vertnext = 0;
  for (vertexPair = vertices(g); vertexPair.first != vertexPair.second; ++vertexPair.first) {
    if ((VertexInfo1Map[vertnext] != N) & (vertnext == *vertexPair.first)) {
      bern = R::rbinom(1, EdgeWeightMap[*out_edges(*vertexPair.first, g).first]);
      if (bern == 1) {
	x.push_back(1);
	vertnext = target(*out_edges(*vertexPair.first, g).first, g);
      }
      else {
	x.push_back(2);
	vertnext = target(*(++out_edges(*vertexPair.first, g).first), g);
      }
    }
  }
  return x;
}
