#include <Rcpp.h>
#include "vector.h"
#include "quadTree.h"

using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix nbody(NumericMatrix pos, NumericVector strength, double theta, double min_dist, double max_dist, double alpha) {
  size_t i;
  NumericMatrix res(pos.nrow(), pos.ncol());
  QuadTree<2> tree(theta, min_dist, max_dist, alpha);
  std::deque<Body<2> *> bodies;
  for (i = 0; i < pos.nrow(); ++i) {
    Body<2> * body = new Body<2>();
    body->pos.coord[0] = pos(i, 0);
    body->pos.coord[1] = pos(i, 1);
    body->strength = strength[i];
    bodies.push_back(body);
  }

  tree.insertBodies(bodies);

  for (i = 0; i < pos.nrow(); ++i) {
    tree.updateBodyForce(bodies[i]);
  }
  for (i = 0; i < pos.nrow(); ++i) {
    res(i, 0) = bodies[i]->force.coord[0];
    res(i, 1) = bodies[i]->force.coord[1];
  }

  return res;
}
