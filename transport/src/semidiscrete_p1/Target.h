/*
 * Author: Valentin Hartmann
 */

#ifndef TARGET_H
#define TARGET_H

#include <vector>

class Target {
  protected:
    /*
     * NOTE: We don't use a std::map<std::pair<double, double>, double>* which
     * contains the points and the associated masses since we need an order for
     * the weight vector in the Voronoi diagram.
     */

    // the support of the target measure
    std::vector<std::pair<double, double> > points;

    // the target's masses
    // the accumulated mass is always 1
    std::vector<double> masses;

    // Is the target created from a grid?
    const bool grid;

  public:
    Target(
      const std::vector<std::pair<double, double> >& points,
      const std::vector<double>& masses);

    Target(
      const std::vector<std::pair<double, double> >& points,
      const std::vector<double>& masses,
      const double accMass);

    /*
     * Constructs a Target from a given (pixel) grid.
     * The support will be contained in [0,1]x[0,1].
     * @param rows:     the number of rows of the grid
     * @param cols:     the number of columns of the grid
     * @param accMass:  the sum of all masses
     */
    Target(
      const int rows,
      const int cols,
      const std::vector<double>& masses,
      const double accMass);

    const std::vector<std::pair<double, double> >& getPoints() const {return points;}
    const std::vector<double>& getMasses() const {return masses;}

    /*
     * Creates a measure with the cardinal number of its support equal to
     * count / reduction (rounded down).
     * @return Target*: the created measure
     * @return std::vector<std::vector<int> >:  a vector of count/reduction
                        vectors where the vector vector[i] contains the points
                        of the original measure that are closer to the i-th
                        point of the created measure than to all other points in
                        the support of the latter (where distance means
                        Euclidean distance)
     */
    std::pair<Target*, std::vector<std::vector<int> > > coarsen(int reduction);
};

#endif
