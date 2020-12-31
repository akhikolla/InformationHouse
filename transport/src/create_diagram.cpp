#ifdef WITH_CGAL
  #include <stdio.h>
  #include <iostream>
  //#include <fstream>
  #include <cassert>
  //#include <time.h>
  #include <sstream>
#endif

#include <Rcpp.h>

using namespace Rcpp;

#ifdef WITH_CGAL

// define the kernel
typedef double numberType;
#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<numberType>  Kernel;

// includes for defining the Voronoi diagram adaptor
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Apollonius_graph_adaptation_traits_2.h>
#include <CGAL/Apollonius_graph_adaptation_policies_2.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

// define the Voronoi diagram adaptor
typedef CGAL::Apollonius_graph_traits_2<Kernel> Traits;
typedef CGAL::Apollonius_graph_2<Traits> Apo_graph;
typedef CGAL::Apollonius_graph_adaptation_traits_2<Apo_graph> AT;
typedef CGAL::Apollonius_graph_degeneracy_removal_policy_2<Apo_graph> AP;
typedef CGAL::Voronoi_diagram_2<Apo_graph,AT,AP> Diagram;

typedef CGAL::Apollonius_site_2<Kernel> Site;
typedef Site::Point_2 Point;
typedef Site::Weight Weight;
typedef Diagram::Face Face;
typedef Diagram::Halfedge Halfedge;
typedef Diagram::Ccb_halfedge_circulator Ccb_halfedge_circulator;
typedef Diagram::Face_iterator Face_iterator;



void insertSite(Apo_graph* ag, FILE* f, double x, double y, double weight) {
  Point p(x, y);
  Weight w(weight);
  Site s(p, w);
  ag->insert(s);
  if (f != NULL) {
    fprintf(f, "%.30f %.30f %.30f\n", p.x(), p.y(), w);
  }
}

// inserts sites outside of the visible area so that each hyperbola inside has
// finite endpoints
void insertBoundingSites(Apo_graph* ag, double maxX, double maxY, double maxWeight) {
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      // I don't use x and y in the header of the for loop directly because of
      // to avoid the comparison of unprecise double values.
      double x = -maxWeight - maxX + i * (2*maxWeight + 3*maxX);
      double y = -maxWeight - maxY + j * (2*maxWeight + 3*maxY);
      Point p(x, y);
      Weight w(0);
      Site s(p, w);
      ag->insert(s);
    }
  }
}

// generates the Apollonius graph of n random points with random weights in the
// rectangle width x height
// the generated data is saved to f
/* Apo_graph* generateApo_graphRandom(int n, double width, double height, std::string sitesPath) {
  Apo_graph* ag = new Apo_graph();
  FILE* sitesFile = fopen(sitesPath.c_str(), "w");

  double maxWeight = 0;

  for (int i = 0; i < n; i++) {
    double x = width * rand()/RAND_MAX;
    double y = height * rand()/RAND_MAX;
    double weight = std::min(width, height)/6 * rand()/RAND_MAX;
    maxWeight = std::max(maxWeight, weight);
    insertSite(ag, sitesFile, x, y, weight);
  }

  insertBoundingSites(ag, width, height, maxWeight);

  fclose(sitesFile);

	return ag;
}
*/
 
// takes a sites file as input and thus doesn't need to create one
Apo_graph* generateApo_graphCustom(NumericMatrix& sites) {
  // note: speed gain of pass-by-ref is minor here
  // (somehow since sites is Rcpp object wrapping a SEXP, which is a pointer)
  if (sites.ncol() != 3) {
    stop("Matrix 'sites' should have three columns x, y, weights");
  }
    
  double maxWeight = 0;
  double maxX = 0;
  double maxY = 0;

  maxX = Rcpp::max(sites(_,0));
  maxY = Rcpp::max(sites(_,1));
  maxWeight = Rcpp::max(sites(_,2));

  Apo_graph* ag = new Apo_graph();

  insertBoundingSites(ag, maxX, maxY, maxWeight);

  int numSites = sites.nrow();

  for (int i = 0; i < numSites; ++i) {
    insertSite(ag, NULL, sites(i,0), sites(i,1), sites(i,2));
  }

  return ag;
}

/* Apo_graph* generateApo_graphGrid(std::string weightsPath, std::string sitesPath) {
  FILE* sitesFile = fopen(sitesPath.c_str(), "w");

  std::ifstream weightsStream(weightsPath.c_str());
  int rows = -1, cols = 0;
  std::vector<double> weights;
  double maxWeight = 0;

  while (weightsStream) {
    rows++;
    std::string currLine;
    getline(weightsStream, currLine);

    std::stringstream currLineStream;
    currLineStream << currLine;
    double currWeight;

    while (currLineStream >> currWeight) {
      weights.push_back(currWeight);
      maxWeight = std::max(maxWeight, currWeight);
      if (rows == 0) {
        cols++;
      }
    }
  }

  Apo_graph* ag = new Apo_graph();

  insertBoundingSites(ag, 1, 1, maxWeight);

  // the coordinates of the points are divided by normalizingFactor to make them
  // be contained in the square [0,1]x[0,1]
  int normalizingFactor = std::max(rows, cols);

  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      double x = (col + 0.5) / normalizingFactor;
      double y = (rows - row - 0.5) / normalizingFactor;
      double weight = weights[row*cols + col];
      insertSite(ag, sitesFile, x, y, weight);
    }
  }

  fclose(sitesFile);

  return ag;
}
 */

NumericMatrix formatIntersections(Apo_graph* ag) {

  Diagram diagram(*ag);
  NumericVector intersection_vec(0);
  // DS: For now I append to a vector and transform
  // it then into a i times 10 matrix
  // Probably we know the (maybe only maximal) number of rows
  // of the matrix beforehand: then it would be more efficient to set up
  // such a matrix here and fill it in below

  // To prevent edges from being considered from both sides, we collect the
  // already visited ones in a set.
  std::set<Halfedge> edges;
  Halfedge currHalfedge, currOppositeHalfedge;
  // the two sites adjacent to a halfedge
  Site origin, neighbor;
  // the two endpoints of a halfedge
  Point source, target;

  int numNonTrivEdges = 0, numEdges = 0, numHalfedges = 0;

  for (Face_iterator face_it = diagram.faces_begin();
    face_it != diagram.faces_end(); face_it++) {
    Face currFace = *face_it;
    origin = currFace.dual()->site();

    Ccb_halfedge_circulator edge_ci, edge_ci_end;
    edge_ci = edge_ci_end = currFace.ccb();
    do {
      numHalfedges++;
      currHalfedge = *edge_ci;
      currOppositeHalfedge = *currHalfedge.opposite();
      // the current edge hasn't been visited yet
      if (edges.find(currOppositeHalfedge) == edges.end()) {
        numEdges++;

        neighbor = currOppositeHalfedge.face()->dual()->site();

        if (!currHalfedge.has_source() || !currHalfedge.has_target()) {
          continue;
        }
        
        numNonTrivEdges++;
        intersection_vec.push_back(origin.point().x());
        intersection_vec.push_back(origin.point().y());
        intersection_vec.push_back(origin.weight());
        intersection_vec.push_back(neighbor.point().x());
        intersection_vec.push_back(neighbor.point().y());
        intersection_vec.push_back(neighbor.weight());
 
        // the first endpoint of the hyperbola
        if (currHalfedge.has_source()) {
          source = (*currHalfedge.source()).point();
          intersection_vec.push_back(source.x());
          intersection_vec.push_back(source.y());
        } else {
          intersection_vec.push_back(R_NaN);
          intersection_vec.push_back(R_NaN);
        }

        // the second endpoint of the hyperbola
        if (currHalfedge.has_target()) {
          target = (*currHalfedge.target()).point();
          intersection_vec.push_back(target.x());
          intersection_vec.push_back(target.y());
        } else {
          intersection_vec.push_back(R_NaN);
          intersection_vec.push_back(R_NaN);
        }
      }
      edges.insert(currHalfedge);
    } while(++edge_ci != edge_ci_end);
  }
  // Rcpp is column major, therefore:
  NumericMatrix intersections(10,numNonTrivEdges,intersection_vec.begin());
  return(transpose(intersections));
}

#endif

// [[Rcpp::export]]
NumericMatrix create_diagram(NumericMatrix sites) {

#ifndef WITH_CGAL
  Rcpp::stop("C++ function create_diagram called, but CGAL not available.");
  return(NumericMatrix(1,1));
}
#endif
#ifdef WITH_CGAL
  Apo_graph* ag;

  ag = generateApo_graphCustom(sites);
  return(formatIntersections(ag));
}

#endif
