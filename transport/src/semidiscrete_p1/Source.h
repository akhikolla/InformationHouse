/*
 * Author: Valentin Hartmann
 */

#ifndef SOURCE_H
#define SOURCE_H

#include "semidiscrete_p1/Target.h"
#include "semidiscrete_p1/lbfgs.h"
#include <math.h>
#include <vector>
#include <map>
#include <Rcpp.h>

// the number type
typedef double NT;
#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<NT>  Kernel;

#include <CGAL/Apollonius_graph_hierarchy_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>


// ***** configuration of CGAL and CGAL typedefs *****
typedef CGAL::Apollonius_graph_traits_2<Kernel> Traits;
typedef CGAL::Apollonius_graph_2<Traits> ApoGraph;
typedef CGAL::Apollonius_site_2<Kernel> Site;
typedef Site::Point_2 Point;
typedef Site::Weight Weight;
typedef ApoGraph::Vertex_handle VertexHandle;

class Source {
  protected:
    // the pixel's masses, i.e. colors, in the order
    // [0][0], [0][1], [0][2], ..., [1][0], [1][1], ...
    std::vector<double> pixelMasses;

    // the pixel's masses, normalized to compensate the refinement factor
    std::vector<double> refinedPixelMasses;

    // the number of rows, columns of the picture
    const int rows;
    const int cols;

    // the pixels will be split into 2^refinement subpixels
    int refinement;

    // the target measure
    const Target* target;

    // the parameters for the L-BFGS algorithm
    lbfgs_parameter_t* params;

    // the weight vector which will be optimized by the L-BFGS algorithm
    lbfgsfloatval_t* weights;

    // the Wasserstein distance between the source and the target measure for
    // the weight vector of the current optimization step
    double wasserstein;

    // the squared Wasserstein distance from the last function evaluation;
    // wasserstein will only be updated with this value if the weight vector of
    // this function evaluation becomes the next step's weight vector
    // as an explanation: lbfgs usually performs several function evaluations
    // with different weight vectors for each optimization step
    double wassersteinUpdate;

    // The following two variables are important if the optimization process
    // stops before the stopping criterion is reached because in this case the
    // final weight vector is not necessarily the one at which minTransported is
    // reached.

    // the minimum L^1 norm of the gradient reached so far in the optimization
    // process
    double minMistransported;
    // the weight vector at which minMistransported is reached
    lbfgsfloatval_t* minWeights = NULL;

  public:
    // the optimization terminates when the L^1 norm of the gradient is < EPS
    constexpr static double EPS = 5e-3;

    // the linesearch algorithm used for the optimization process
    // const static int LINESEARCH = LBFGS_LINESEARCH_DEFAULT;
    const static int LINESEARCH = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
    // const static int LINESEARCH = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    // const static int LINESEARCH = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;

    // the maximum number of trials for the line search
    const static int MAX_LINESEARCH = 200;

    // the number of previous iteration results used to approximate the inverse
    // hessian matrix
    const static int NUM_CORRECTIONS = 6;

    // the factor with which the L2-summand for stabilizing the optimization is
    // multiplied
    // higher value -> more stable but less precise
    const double REGULARIZATION_STRENGTH;
    
    /*
     * @param pixelMasses:    the pixel's masses, i.e., colors
     * @param rows, cols:     the number of rows, columns of the picture
     * @param accPixelMasses: the accumulated mass of the pixels, i.e., \mu(R^2)
                              (if known in advance, will be computed otherwise)
     */
    Source(
      const std::vector<double>& pixelMasses,
      const int rows,
      const int cols,
      const double accPixelMass,
      const double REGULARIZATION_STRENGTH);
    Source(
      const std::vector<double>& pixelMasses,
      const int rows,
      const int cols,
      const double REGULARIZATION_STRENGTH);

    ~Source();

    lbfgsfloatval_t* getWeights() const {return weights;}

    double getMinMistransported() const {return minMistransported;}

    lbfgsfloatval_t* getMinWeights() const {return minWeights;}

    double getWasserstein() const {return wasserstein;}

    /*
     * Starts the optimization.
    * @param target:          the target measure;
                              its support must be contained in the rectangle
                              [0,rows]x[0,cols]
    * @param startWeights:    the initial weight vector from which the
                              optimization will be started
    * @param refinementRatio: the pixels of the source measure will be split
                              into subpixels such that
                              |supp(source)| / |supp(target)| >= refinementRatio
                              0 means no refinement
    * @return:                the status code of lbfgs
    */
    int optimize(
      const Target* const target,
      lbfgsfloatval_t* const startWeights,
      const int refinementRatio = 0);

    /* Prints the transport plan for the given data to the file "filename".
     * @param target:   the target measure
     * @param weights:  the weight vector
     * @param filename: the name of the file the transport plan is printed to
     */
    void createTransportPlan(
      const Target* const target,
      const lbfgsfloatval_t* const weights,
      Rcpp::NumericMatrix transportplan);

  protected:
    // Initializes the parameters for the L-BFGS algorithm.
    void initLbfgs();

    // Initializes the floating point values since that can't be done in the
    // .h-file.
    void initFloats();

    /*
     * Creates and fills the Apollonius graph with the current target measure
     * and weight vector.
     * @param x:              the weight vector
     * @param targetHandles:  same as target->getPoints() but with the
                              VertexHandles belonging to the points;
                              used for faster access to the indeces
     */
    ApoGraph* createApoGraph(
      const lbfgsfloatval_t* const x,
      std::map<Point/*center of the Voronoi cell*/,
      int/*index*/>* const targetHandles);


    // Required by lbfgs. Should evaluate the function to be optimized and its
    // gradient.
    static lbfgsfloatval_t _evaluate(
      void* instance,
      const lbfgsfloatval_t* x,
      lbfgsfloatval_t* g,
      const int n,
      const lbfgsfloatval_t step)
    {return reinterpret_cast<Source*>(instance)->evaluate(x, g, n, step, Rcpp::NumericMatrix(1,1));}

    /*
     * @param transportplan: the transport plan is printed to NumericMatrix
     *                       transportplan (3 columns: from-index, to-index, mass) 
     *                       unless it is NumericMatrix(1,1), in which case we need the g
     */
    lbfgsfloatval_t evaluate(
      const lbfgsfloatval_t* x,
      lbfgsfloatval_t* g,
      const int n,
      const lbfgsfloatval_t step,
      Rcpp::NumericMatrix transportplan);

    // Callback of lbfgs. Called whenever an optimization step is finished.
    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls)
    {return reinterpret_cast<Source*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);}

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls);

    /*
     * Computes the integral of ||x-p|| over the given square pixel.
     * @param x1, x2: first and second coordinate of the pixel's top left vertex
     * @param size:   width/height of the pixel
     */
    static double distIntegral(
      double x1, double x2,
      double p1, double p2,
      double size = 1);

    /*
     * Approximates distIntegral. Note the different parameters.
     * @param x1, x2:   first and second coordinate of the pixel's bottom left
     *                  vertex
     * @param sizeHalf: HALF the size of the pixel
     */
    static double distIntegralApprox(
      double x1, double x2,
      double p1, double p2,
      double sizeHalf = 0.5);

  private:
    // ***** helper functions for distIntegral *****

    static double normPlus(double a1, double a2) {return sqrt(a1*a1 + a2*a2);}
    static double euclideanDist(double a1, double a2, double b1, double b2) {
      return normPlus(a1 - b1, a2 - b2);
    }

    /*
     * @param x: the weight vector
     * @param g: the gradient of Phi
     * @param n: the length of the weight and gradient vector
     * @return: the amount of mistransported mass for the current transport
     * partition
     */
    double computeMistransported(
      const lbfgsfloatval_t *x,
      const lbfgsfloatval_t *g,
      int n);

    /*
     * @param pos:          the position at which the primitive is computed
     * @param posl2, posu2: corresponds to l1l2, l1u2 resp. or u1l2, u1u2 resp.
                            where l1u2 = normPlus(l1, u2) and l and u correspond
                            to the lower and upper integration limits
     */
    static double primitive(
      double pos,
      double posl2, double posu2,
      double l2, double u2);

    // *********************************************

    /*
     * common code of
     * Source(const std::vector<double>& pixelMasses, int rows, int cols, double accPixelMass)
     * and
     * Source(const std::vector<double>& pixelMasses, int rows, int cols)
     */
    void constructFromVector(
      const std::vector<double>& pixelMasses,
      const double accPixelWeight);
};

#endif
