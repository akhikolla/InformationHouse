/*
 * Author: Valentin Hartmann
 */

#include "Target.h"

#include <math.h>
#include <stdlib.h> // rand
#include <set>
#include <limits>
#include <Rcpp.h>

Target::Target(
  const std::vector<std::pair<double, double> >& points,
  const std::vector<double>& masses)
  : points(points), masses(masses), grid(false) {
    // the masses are normalized such that the accumulated mass is 1
    double accMass = 0;
    for (std::vector<double>::iterator it = this->masses.begin(); it != this->masses.end(); ++it) {
      accMass += *it;
    }
    for (std::vector<double>::iterator it = this->masses.begin(); it != this->masses.end(); ++it) {
      *it = *it / accMass;
    }
}

Target::Target(
  const std::vector<std::pair<double, double> >& points,
  const std::vector<double>& masses,
  const double accMass)
  : points(points), masses(masses), grid(false) {
    // the masses are normalized such that the accumulated mass is 1
    for (std::vector<double>::iterator it = this->masses.begin(); it != this->masses.end(); ++it) {
      *it = *it / accMass;
    }
}

Target::Target(
  const int rows,
  const int cols,
  const std::vector<double>& masses,
  const double accMass)
  : masses(masses), grid(true) {
    // the masses are normalized such that the accumulated mass is 1
    for (std::vector<double>::iterator it = this->masses.begin(); it != this->masses.end(); ++it) {
      *it = *it / accMass;
    }

    points.resize(masses.size());
    // the target measure's support is normalized to be contained in [0,1]x[0,1]
    double normalizingFactor = std::max(rows, cols);
    for (int i = rows - 1; i >= 0; --i) {
      for (int j = 0; j < cols; j++) {
        points[(rows - 1 - i)*cols + j] = std::make_pair((j + 0.5) / normalizingFactor, (i + 0.5) / normalizingFactor);
      }
    }
}

std::pair<Target*, std::vector<std::vector<int> > > Target::coarsen(int reduction) {
  // NOTE: new... denotes the attributes of the newly created measure (with
  // smaller support) to prevent possible confusion that could occur when using
  // e.g. points and this->points

  int newCount = points.size() / reduction;
  std::vector<std::pair<double, double> > newPoints(newCount);
  std::vector<double> newMasses(newCount);

  typedef std::vector<std::vector<int> > Assignments;
  // the vector to be returned
  Assignments assignments(newCount);
  for (Assignments::iterator it = assignments.begin(); it != assignments.end(); ++it) {
    *it = std::vector<int>();
  }

  // ***** initialization of Lloyd's algorithm *****

  // initialize Lloyd's algorithm with random points from the original target
  // using a std::set to enforce uniqueness
  std::set<std::pair<double, double> > pointSet;
  while (pointSet.size() < newCount) {
    int index = Rcpp::runif(1, 0.0, 1.0)[0] * ((double) points.size());
    pointSet.insert(points[index]);
  }

  std::vector<std::pair<double, double> >::iterator newPointsIt = newPoints.begin();
  for (std::set<std::pair<double, double> >::iterator pointSetIt = pointSet.begin();
    pointSetIt != pointSet.end(); ++pointSetIt, ++newPointsIt) {
    *newPointsIt = *pointSetIt;
  }

  // ***********************************************


  // ***** Lloyd's algorithm *****

  // the points from the last iteration of the algorithm to check whether it has
  // converged yet
  std::vector<std::pair<double, double> > pointsLastIteration;

  // // the number of points that changed during the last iteration
  // int diff;

  do {
    pointsLastIteration = newPoints;

    // clear assignments
    for (Assignments::iterator it = assignments.begin(); it != assignments.end(); ++it) {
      std::vector<int>().swap(*it);
    }

    // * maximization step *
    for (int i = 0; i < points.size(); ++i) {
      // the currently nearest point of the new measure
      int nearest;
      // the distance to the currently nearest point of the new measure
      double nearestDist = std::numeric_limits<double>::infinity();
      for (int j = 0; j < newPoints.size(); ++j) {
        // using the squared Euclidean distance instead of the usual Euclidean
        // distance speeds up computation and doesn't change the result
        double currDist = pow(newPoints[j].first - points[i].first, 2) +
                        pow(newPoints[j].second - points[i].second, 2);
        if (currDist < nearestDist) {
          nearestDist = currDist;
          nearest = j;
        }
      }
      assignments[nearest].push_back(i);
    }

    // * expectation step *
    // the points of the original target the current point is computed from
    std::vector<int> currAssignedPoints;

    for (int i = 0; i < newCount; i++) {
      newPoints[i] = std::make_pair(0, 0);
      newMasses[i] = 0;
      currAssignedPoints = assignments[i];
      // to prevent dividing by 0 below
      if (currAssignedPoints.size()) {
        for (std::vector<int>::iterator it = currAssignedPoints.begin();
          it != currAssignedPoints.end(); ++it) {
          // multiply with the mass since we solve a *weighted* k-means problem
          newPoints[i].first += points[*it].first * masses[*it];
          newPoints[i].second += points[*it].second * masses[*it];
          newMasses[i] += masses[*it];
        }
        // if it was an unweighted problem:
        // newPoints[i].first /= (currAssignedPoints.size() * newMasses[i]);
        newPoints[i].first /= newMasses[i];
        newPoints[i].second /= newMasses[i];
      }
    }

  //   diff = 0;
  //   std::vector<std::pair<double, double> >::iterator newPointsIt = newPoints->begin();
  //   for (std::vector<std::pair<double, double> >::iterator pointsLastIterationIt = pointsLastIteration->begin(); pointsLastIterationIt != pointsLastIteration->end(); ++pointsLastIterationIt, ++newPointsIt) {
  //     if (*newPointsIt != *pointsLastIterationIt) {
  //       diff++;
  //     }
  //   }
  // } while (diff > 100000);
  } while (newPoints != pointsLastIteration);

  // *****************************

  return std::make_pair(new Target(newPoints, newMasses), assignments);
}
