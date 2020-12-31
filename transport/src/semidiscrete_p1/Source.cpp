/*
 * Author: Valentin Hartmann
 */

#include "Source.h"
#include <stdlib.h> // NULL
#include <map>
#include <set>
#include <vector>
#include <Rcpp.h>

/*
 * debug output
 */
// #define PRINT_PROGRESS

/*
 * stores minMistransported and minWeights
 * increases computation time of evaluate()
 */
// #define STORE_INTERMEDIATE_RESULTS

Source::Source(
  const std::vector<double>& pixelMasses,
  const int rows,
  const int cols,
  const double accPixelMass,
  const double REGULARIZATION_STRENGTH)
  : rows(rows), cols(cols), REGULARIZATION_STRENGTH(REGULARIZATION_STRENGTH)
  {
    initLbfgs();
    constructFromVector(pixelMasses, accPixelMass);
}

Source::Source(
  const std::vector<double>& pixelMasses,
  const int rows,
  const int cols,
  const double REGULARIZATION_STRENGTH)
  : rows(rows), cols(cols), REGULARIZATION_STRENGTH(REGULARIZATION_STRENGTH)
  {
    // the accumulated mass of the pixels, i.e. \mu(R^2)
    double accPixelMass = 0;
    for (std::vector<double>::const_iterator it = pixelMasses.begin(); it != pixelMasses.end(); ++it) {
      accPixelMass += *it;
    }

    initLbfgs();
    constructFromVector(pixelMasses, accPixelMass);
}

Source::~Source() {
  if (minWeights) {
    lbfgs_free(minWeights);
  }
}

void Source::initFloats() {
  wasserstein = 0;
}

void Source::initLbfgs() {
  params = new lbfgs_parameter_t();
  lbfgs_parameter_init(params);
  // we have our own convergence criterion
  params->epsilon = 0;
  params->linesearch = LINESEARCH;
  params->max_linesearch = MAX_LINESEARCH;
  params->m = NUM_CORRECTIONS;
}

void Source::constructFromVector(
  const std::vector<double>& pixelMasses,
  const double accPixelMass) {
  this->pixelMasses.resize(pixelMasses.size());
  for (int i = 0; i < pixelMasses.size(); i++) {
    this->pixelMasses[i] = pixelMasses[i] / accPixelMass;
  }

  refinedPixelMasses.resize(pixelMasses.size());
}


void Source::createTransportPlan(
  const Target* const target,
  const lbfgsfloatval_t* const weights,
  Rcpp::NumericMatrix transportplan)
  {
    evaluate(weights, NULL, target->getPoints().size(), 0, transportplan);
}


int Source::optimize(
  const Target* const target,
  lbfgsfloatval_t* const startWeights,
  const int refinementRatio)
  {
    minMistransported = 2;
    if (minWeights) {
      lbfgs_free(minWeights);
      minWeights = NULL;
    }

    this->target = target;
    int numTargetPoints = target->getPoints().size();

    if (refinementRatio > 0) {
      double refinementFactor = refinementRatio * numTargetPoints / (double) pixelMasses.size();
      refinement = ceil(sqrt(refinementFactor));

      // normalize the pixelMasses
      // (due to the refinement their number and therefore sum is multiplied by
      // refinement^2)
      refinementFactor = pow(refinement, 2);
      if (refinement == 1) {
        refinementFactor = 1;
      }
      for (int i = 0; i < pixelMasses.size(); i++) {
        refinedPixelMasses[i] = pixelMasses[i] / refinementFactor;
      }
    }

    weights = startWeights;

    int lbfgsReturn = 0;
    if (numTargetPoints < 10) {
      // For small instances Wolfe works better.
      params->linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
      lbfgsReturn = lbfgs(numTargetPoints, weights, NULL, _evaluate, _progress, this, params);
      params->linesearch = LINESEARCH;
    } else {
      lbfgsReturn = lbfgs(numTargetPoints, weights, NULL, _evaluate, _progress, this, params);
    }
    
    return lbfgsReturn;
}

int Source::progress(
  const lbfgsfloatval_t *x,
  const lbfgsfloatval_t *g,
  const lbfgsfloatval_t fx,
  const lbfgsfloatval_t xnorm,
  const lbfgsfloatval_t gnorm,
  const lbfgsfloatval_t step,
  int n,
  int k,
  int ls)
  {
    wasserstein = wassersteinUpdate;

    double mistransported = computeMistransported(x, g, n);
  
    #ifdef PRINT_PROGRESS
    printf("Iteration %d:\n", k);
    printf("Euclidean norms: weights: %f, gradient: %f\n", xnorm, gnorm);
    printf("Phi(x): %f\n", fx);
    printf("Step size: %f\n", step);
    printf("mistransported: %f\n", mistransported);
    printf("\n");
    #endif

    // approximation is sufficiently precise; cancel the optimization process
    // if (gnorm < EPS) {
    if (mistransported < EPS) {
      return 1;
    }
    return 0;
}

lbfgsfloatval_t Source::evaluate(
  const lbfgsfloatval_t* x,
  lbfgsfloatval_t* g,
  const int n,
  const lbfgsfloatval_t step,
  Rcpp::NumericMatrix transportplan)
  {
    // std::cout << "evaluation" << std::endl;

    double phi = 0;

    if (transportplan.nrow() == 1 && transportplan.ncol() == 1) {
      // can only happen if evaluate is called from _evaluate
      for (int i = 0; i < n; i++) {
        double currTargetMass = target->getMasses()[i];
        g[i] = -currTargetMass;
        phi -= currTargetMass * x[i];
      }
    }

    std::map<std::pair<int/*source*/, int/*target*/>, double/*mass*/> transportPlan;

    // same as target->getPoints() but with the VertexHandles belonging to the points;
    // used for faster access to the indices
    std::map<Point/*center of the Voronoi cell*/, int/*index*/> targetHandles;

    ApoGraph* graph = createApoGraph(x, &targetHandles);

    // the cell in which the left neighbor of the current square lies
    VertexHandle leftCell;
    // the cell in which the top neighbor of the current square lies
    // (for the case that we reached started at the beginning of a new row)
    VertexHandle topCell = NULL;
    // the cell in which the current square lies
    VertexHandle currCell;

    wassersteinUpdate = 0;
    // following four variables: increase since the values don't need to be
    // recomputed in every iteration
    int refinedRows = refinement * rows;
    int refinedCols = refinement * cols;
    // the source measure's support is normalized to be contained in [0,1]x[0,1]
    double normalizingFactor = std::max(refinedRows, refinedCols);
    double pixelSize = 1 / normalizingFactor;
    double pixelSizeHalf = pixelSize / 2;
    for (int i = refinedRows - 1; i >= 0; --i) {
      for (int j = 0; j < refinedCols; j++) {
        // the center of the current square
        Point currSquare((j + 0.5) / normalizingFactor, (i + 0.5) / normalizingFactor);
        if (j == 0) {
          currCell = graph->nearest_neighbor(currSquare, topCell);
          topCell = currCell;
          leftCell = currCell;
        } else {
          currCell = graph->nearest_neighbor(currSquare, leftCell);
          leftCell = currCell;
        }

        if (transportplan.nrow() == 1 && transportplan.ncol() == 1) {
          // the center of currCell
          Point currCenter = currCell->site().point();

          // the index of currCell
          int currIndex = targetHandles[currCenter];

          // The expression in the index can't be simplified since it relies on
          // rounding down.
          double massCurrPixel = refinedPixelMasses
            [/*keep this as it is*/((refinedRows - 1 - i)/refinement)*cols + j/refinement];
          g[currIndex] += massCurrPixel;
          phi += x[currIndex] * massCurrPixel;
          // double topOfPixel = (i + 1) / normalizingFactor;
          double bottomOfPixel = i / normalizingFactor;
          double leftOfPixel = j / normalizingFactor;
          // wassersteinUpdate += massCurrPixel *
          //   distIntegral(leftOfPixel, topOfPixel,
          //     currCenter.x(), currCenter.y(), pixelSize);
          wassersteinUpdate += massCurrPixel *
            distIntegralApprox(leftOfPixel, bottomOfPixel,
              currCenter.x(), currCenter.y(), pixelSizeHalf);
        } else {
          int source = ((refinedRows - 1 - i)/refinement)*cols + j/refinement;
          int target = targetHandles[currCell->site().point()];
          double mass = refinedPixelMasses[source];
          transportPlan[std::make_pair(source, target)] += mass;
        }
      }
    }
    if (transportplan.nrow() == 1 && transportplan.ncol() == 1) {
      phi -= wassersteinUpdate;
    } else {
      int i=0;
      for (std::map<std::pair<int, int>, double>::iterator it = transportPlan.begin();
        it != transportPlan.end(); ++it){
        transportplan(i,1) = it->first.first;
        transportplan(i,2) = it->first.second;
        transportplan(i,3) = it->second;
        i++;
        if (i == transportplan.nrow()) break;
        // fix later: we currently preassign the number of rows that transportplan must have;
        // of course with number of points in the refined source times number of target points
        // we should be on the safe side
      }
    }

    delete graph;

    #ifdef STORE_INTERMEDIATE_RESULTS
    double mistransported = computeMistransported(x, g, n);
    if (mistransported < minMistransported) {
      minMistransported = mistransported;
      if (minWeights) {
        lbfgs_free(minWeights);
      }
      minWeights = lbfgs_malloc(n);
      for (int i = 0; i < n; ++i) {
        minWeights[i] = weights[i];
      }
    }
    #endif

    for (int i = 0; i < n; i++) {
      g[i] += REGULARIZATION_STRENGTH * 2 * x[i];
      phi += REGULARIZATION_STRENGTH * pow(x[i], 2);
    }
    
    // double factor = 1;
    // phi += x[0]>0 ? factor * x[0] : -factor * x[0];
    // g[0] += x[0] > 0 ? factor : -factor;

    // phi += mistransported;

    // printf("\nevaluate\n");
    // printf("Phi(x): %f\n", phi);
    // printf("step: %f\n", step);
    // printf("g[0]: %f, g[1]: %f, g[2]: %f, g[3]: %f\n", g[0], g[1], g[2], g[3]);
    // printf("x[0]: %f, x[1]: %f, x[2]: %f, x[3]: %f\n", x[0], x[1], x[2], x[3]);

    return phi;
}

ApoGraph* Source::createApoGraph(
  const lbfgsfloatval_t* const x,
  std::map<Point/*center of the Voronoi cell*/,
  int/*index*/>* const targetHandles)
  {
  ApoGraph* graph = new ApoGraph();
  int numTargets = target->getPoints().size();
  for (int i = 0; i < numTargets; i++) {
    Point p(target->getPoints()[i].first, target->getPoints()[i].second);
    Site s(p, x[i]);
    VertexHandle vertexCreated = graph->insert(s);
    // otherwise the vertex is not visible and therefore won't be the nearest
    // neighbor of any point
    if (vertexCreated != VertexHandle(NULL)) {
      (*targetHandles)[vertexCreated->site().point()] = i;
    }
  }

  return graph;
}

double Source::distIntegral(
  double x1, double x2,
  double p1, double p2,
  double size)
  {
    // the lower integration limits
    double l1 = x1 - p1, l2 = x2 - p2;
    // the upper integration limits
    double u1 = l1 + size, u2 = l2 + size;
    double l1l2 = normPlus(l1, l2), l1u2 = normPlus(l1, u2),
    u1l2 = normPlus(u1, l2), u1u2 = normPlus(u1, u2);
    return primitive(u1, u1l2, u1u2, l2, u2) - primitive(l1, l1l2, l1u2, l2, u2);
}

double Source::computeMistransported(
  const lbfgsfloatval_t *x,
  const lbfgsfloatval_t *g,
  int n)
  {
    double mistransported = 0;
    for (int i = 0; i < n; ++i) {
      double without_regularizer = g[i] - REGULARIZATION_STRENGTH * 2 * x[i];
      if (without_regularizer < 0) {
        mistransported -= without_regularizer;
      } else {
        mistransported += without_regularizer;
      }
    }
    // Mistransported mass gets counted twice: Once were it is missing and once
    // were it is in surplus.
    mistransported /= 2;
    return mistransported;
}

double Source::primitive(
  double pos,
  double posl2, double posu2,
  double l2, double u2)
  {
    double log_posu2_pos, log_posl2_pos, log_posu2_u2, log_posl2_l2;
    log_posu2_pos = log(posu2 + pos);
    log_posl2_pos = log(posl2 + pos);
    log_posu2_u2 = log(posu2 + u2);
    log_posl2_l2 = log(posl2 + l2);
    if (log_posu2_pos == log(0))
      log_posu2_pos = 0;
    if (log_posl2_pos == log(0))
      log_posl2_pos = 0;
    if (log_posu2_u2 == log(0))
      log_posu2_u2 = 0;
    if (log_posl2_l2 == log(0))
      log_posl2_l2 = 0;
    return ((u2*(pos*posu2 + u2*u2*log_posu2_pos) - l2*(pos*posl2 + l2*l2*log_posl2_pos))/2
    + (pos*(u2*posu2 - l2*posl2 + 2*pos*pos*(log_posu2_u2 - log_posl2_l2)) - u2*u2*u2*log_posu2_pos + l2*l2*l2*log_posl2_pos)/6)/2;
}

/*
 * Approximates distIntegral by assuming ||x-p|| is ||c-p|| everywhere on
 * the pixel where c is the center of the pixel.
 */
double Source::distIntegralApprox(
  double x1, double x2,
  double p1, double p2,
  double sizeHalf) {
    return euclideanDist(x1 + sizeHalf, x2 + sizeHalf, p1, p2);
}
