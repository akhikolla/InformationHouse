/*
 * Author: Valentin Hartmann
 * Adapted for R-packages by Dominic Schuhmacher
 * 
 * (The code if WITH_CGAL is defined is identical to src/opt_transport.cpp 
 * in the unpublished R-package apollonius v0.0-3)
 */

#ifdef WITH_CGAL
  #include "semidiscrete_p1/Source.h"
  #include "semidiscrete_p1/Target.h"

  #include "semidiscrete_p1/lbfgs.h"

  #include <vector>
  #include <time.h>
  // #include <fstream>
  #include <string>
  #include <sstream>
#endif
  
#include <Rcpp.h>

using namespace Rcpp;  
  
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The remaining code only WITH_CGAL  
// except that we have to provide an empty semidiscrete_p1 function also 
// if WITH_CGAL is undefined
#ifdef WITH_CGAL  

/*
 * prints the computed Wasserstein distance
 */
// #define PRINT_WASSERSTEIN_DISTANCE

// ***** debug output *****

/*
 * prints the current level of the coarsening of \nu
 */
// #define PRINT_LEVELS

/*
 * prints the return code of L-BFGS for each coarsening of \nu including the
 * original \nu
 * 1 = minimization was successful
 */
// #define PRINT_LBFGS_RETURN

/*
 * prints the return code of L-BFGS only for the original \nu
 * 1 = minimization was successful
 */
// #define PRINT_FINAL_LBFGS_RETURN

/*
 * prints the error return codes of lbfgs
 */
// #define LBFGS_ERRORS

// ************************

struct MeasureData {
  std::vector<double> masses;
  double accMass;
};

struct GridMeasureData : MeasureData {
  int rows;
  int cols;
};

struct GeneralMeasureData : MeasureData {
  std::vector<std::pair<double, double> > points;
};

/*
 * @param refinement:         the refinement ratio for the source measure
 * @param targetReduction:    |supp(\nu_l)| = |supp(\nu_{l-1})| / targetReduction
 * @param sourcePath:         the path of the file containing the source masses
 * @param targetPath:         the path of the file containing the target masses
 * @param normalizingFactor:  in the case that targetPath contains a list of
 *                            points and corresponding masses this is the factor
 *                            by which all coordinates are divided to make them
 *                            being contained in [0,1]x[0,1]
 *                            The caller must make sure that this is actually
 *                            the case.
 * @param weightsPath:        the file the weights vector should be printed to
 * @param transportPlanPath:  the file the transport plan should be printed to;
 *                            if omitted, the transport plan isn't printed
 * @param regularizationStrength:  the factor with which the L2-summand for
 *                            stabilizing the optimization is multiplied
 *                            higher value -> more stable but less precise
 */
int computeWeightVector(
  int refinement,
  int targetReduction,
  NumericMatrix source_mat,
  NumericMatrix target_mat,
  NumericVector weights_vec,
  bool verbose,
  double* wasserstein_dist,
  NumericMatrix transportplan,
  double normalizingFactor,
  double regularizationStrength);
/*
 * @param filepath: a file containing the masses of a measure on a grid
 * @return: a corresponding GridMeasureData object
 */
GridMeasureData* NumericMatrix2GridMeasureData(NumericMatrix measure_mat);

/*
 * The points are assumed to have non-negative coordinates.
 * Note that this function normalizes the points to be contained in [0,1]x[0,1]
 * in such a way that at least one point lies on the line (x,1) or (1,y).
 * @param filepath:           a file containing the points and corresponding
 *                            masses of a measure in a list
 * @param normalizingFactor:  the factor by which the coordinates should be
 *                            divided
 *                            After division they need to be contained in
 *                            [0,1]x[0,1].
 * @return: a corresponding GeneralMeasureData object
 */
GeneralMeasureData* NumericMatrix2GenMeasureData(NumericMatrix measure_mat, double normalizingFactor);

/*
 * Normalizes a vector such that its smallest entry is 0 by adding a constant to
 * all entries.PRINT_LBFGS_RETURN
 * @param vector:     the vector to be normalized
 * @param numEntries: the number of entries of vector
 */
void normalizeVector(lbfgsfloatval_t* vector, int numEntries);

/*
 * Prints a hint on how to use the command line version of the program.
 */
 // void printUsageHint();
 
#endif

// [[Rcpp::export]]
List semidiscrete_p1(NumericMatrix source_mat, NumericMatrix target_mat, 
                     bool verbose = false, bool target_in_genpos = true,
                     double regularization_strength = 0.0,
                     NumericMatrix transportplan = NumericMatrix(1,1)) {
  //                          const Nullable<NumericMatrix>& transportplan_storage = R_NilValue) {
  //                          leads to too many complications.
  //                          Now: If a single 0 is passed -> no transportplan
  
  // regularization_strength is regularizationStrength in other functions
  
    // We assume that also if target is a list of points in general position
  // these positions are in [0,1] x [0,1]
  // so no scaling factor (a.k.a. normalizingFactor). This is R's job.
  //srand(time(NULL));

  // ***************************************************************************

#ifndef WITH_CGAL
  Rcpp::stop("C++ function semidiscrete_p1 called, but CGAL not available.");
  List res = List::create(0, 0, 0);
  return(res);
}
#endif
#ifdef WITH_CGAL  
  
  int refinement = 1000;
  int targetReduction = 5;
  
  double normalizingFactor = -1.0;
  if (target_in_genpos) {
    normalizingFactor = 1.0;
  }

  NumericVector weights_vec(target_mat.nrow());
  double wasserstein_dist = 0.0;
  double* pwasserstein_dist = &wasserstein_dist;
  int lbfgs_return;
  
  lbfgs_return = computeWeightVector(refinement, targetReduction, source_mat, target_mat,
                      weights_vec, verbose, pwasserstein_dist, transportplan, 1.0, regularization_strength);
    // 1.0 is for the normalizingFactor. 
  List res = List::create(weights_vec, wasserstein_dist, lbfgs_return);
  return(res);
}

// void printUsageHint() {
//   printf("Usage: opt_transport <-l/-g> <source> <target> [<normalizing factor for list target>] <weights> [<transport plan>]\n");
// }

int computeWeightVector(
  int refinement,
  int targetReduction,
  NumericMatrix source_mat,
  NumericMatrix target_mat,
  NumericVector weights_vec,
  bool verbose,
  double* pwasserstein_dist,
  NumericMatrix transportplan,
  double normalizingFactor,
  double regularizationStrength) {

  bool targetAsList;
  if (normalizingFactor < 0) {
    targetAsList = false;
  } else {
    targetAsList = true;
  }

  GridMeasureData* sourceData = NumericMatrix2GridMeasureData(source_mat);
  Source source(sourceData->masses, sourceData->rows, sourceData->cols,
    sourceData->accMass, regularizationStrength);

//  Rcout << "Source: " << std::endl;
//  Rcout << "First and last mass: " << sourceData->masses[0] <<
//    "  " << sourceData->masses[1022] << "  " << sourceData->masses[1023] << std::endl;
//  Rcout << "No. of rows: " << sourceData->rows << std::endl;
//  Rcout << "No. of cols: " << sourceData->cols << std::endl;
//  Rcout << "Total mass: " << sourceData->accMass << std::endl;
  
  MeasureData* targetData;
  Target* target;

  if (targetAsList) {
    targetData = NumericMatrix2GenMeasureData(target_mat, normalizingFactor);
    GeneralMeasureData* generalTargetData = static_cast<GeneralMeasureData*>(
      targetData);
    target = new Target(generalTargetData->points, generalTargetData->masses,
      generalTargetData->accMass);
 
//  Rcout << "Target: " << std::endl;   
//  Rcout << "First point: " << generalTargetData->points[0].first << 
//    "  " << generalTargetData->points[0].second << std::endl;
//  Rcout << "Last point: " << generalTargetData->points[9].first << 
//    "  " << generalTargetData->points[9].second << std::endl;
//  Rcout << "First and last mass: " << generalTargetData->masses[0] <<
//   "  " << generalTargetData->masses[8] << "  " << generalTargetData->masses[9] << std::endl;
//  Rcout << "Total mass: " << generalTargetData->accMass << std::endl;
 
  } else {
    targetData = NumericMatrix2GridMeasureData(target_mat);
    GridMeasureData* gridTargetData = static_cast<GridMeasureData*>(targetData);
    target = new Target(gridTargetData->rows, gridTargetData->cols,
      gridTargetData->masses, gridTargetData->accMass);
  }


  // ***** the coarsening of the target measure *****

  int numberOfMeasures = log(targetData->masses.size()) / log(targetReduction);
  // if we increased numberOfMeasures also for <= 2,
  // we would have |supp(\nu_L)| == 1
  if (targetData->masses.size() / pow(5, numberOfMeasures) >= 2) {
    ++numberOfMeasures;
  }
  typedef std::vector<std::pair<Target*, std::vector<std::vector<int> > > > Coarsening;
  Coarsening coarsening(numberOfMeasures);
  coarsening[0] = std::make_pair(target, std::vector<std::vector<int> >());
  for (int l = 1; l < coarsening.size(); ++l) {
    coarsening[l] = coarsening[l - 1].first->coarsen(targetReduction);
  }

  // the start weights of \nu_L
  int supportSizeNuL = coarsening.back().first->getPoints().size();
  lbfgsfloatval_t* startWeights = lbfgs_malloc(supportSizeNuL);
  for (int i = 0; i < supportSizeNuL; ++i) {
    startWeights[i] = 0;
  }

  // ***************************************************


  // ***** the minimization of \Phi *****

  int lbfgs_return;

  for (int l = coarsening.size() - 1; l > 0; --l) {
    if (verbose)
      Rprintf("***********\nLevel of Coarsening: %d\n***********\n\n", l);

    lbfgs_return = source.optimize(coarsening[l].first, startWeights, refinement);

    if (verbose)
      Rprintf("L-BFGS return code: %d\n\n", lbfgs_return);

    int lengthOldStartWeights = coarsening[l].first->getPoints().size();
    lbfgsfloatval_t* oldStartWeights = lbfgs_malloc(lengthOldStartWeights);
    for (int i = 0; i < lengthOldStartWeights; ++i) {
      oldStartWeights[i] = startWeights[i];
    }

    lbfgs_free(startWeights);
    startWeights = lbfgs_malloc(coarsening[l - 1].first->getPoints().size());

    for (int i = 0; i < coarsening[l].second.size(); ++i) {
      // the points of \nu_{l-1} that belong to the i-th point of \nu_l
      std::vector<int> currAssignedPoints = coarsening[l].second[i];
      // the weight of the i-th point of \nu_l
      double currWeight = oldStartWeights[i];
      for (int j = 0; j < currAssignedPoints.size(); ++j) {
        startWeights[currAssignedPoints[j]] = currWeight;
      }
    }
    lbfgs_free(oldStartWeights);
  }

  if (verbose)
    Rprintf("***********\nLevel of Coarsening: 0\n***********\n\n");

  lbfgs_return =  source.optimize(coarsening[0].first, startWeights, refinement);

  if (verbose)
    Rprintf("L-BFGS return code: %d\n", lbfgs_return);

  *pwasserstein_dist = source.getWasserstein();
  
  if (verbose)
    Rprintf("Wasserstein distance: %.20f\n", *pwasserstein_dist);
  
  checkUserInterrupt();

  // ************************************


  // ***** fill inte the NumericVector weights_vec *****

  // int precision = 10;

  normalizeVector(startWeights, target->getPoints().size());

  // weights_vec will contain a weight vector (as Matrix row major!)
  // no matter whether target is genpos or grid
  int len = target->getPoints().size();
  if (len != weights_vec.length()) {
    Rcout << "Provided length: " << weights_vec.length() << ";   Computed length: "
          << len << std::endl;
    stop("The two lengths for NumericVector weights_vec do not match");
  }
  for (int i = 0; i < len; i++) {
    weights_vec[i] = startWeights[i];
  }

  if (!(transportplan.nrow()==1 && transportplan.ncol()==1)) {
    source.createTransportPlan(coarsening[0].first, startWeights, transportplan);
  }

  lbfgs_free(startWeights);

  for (Coarsening::iterator it = coarsening.begin(); it != coarsening.end(); ++it) {
    delete (*it).first;
  }
  delete targetData;
  return lbfgs_return;
}


GridMeasureData* NumericMatrix2GridMeasureData(NumericMatrix measure_mat) {
 
  GridMeasureData* m = new GridMeasureData;
  m->accMass = 0;
  m->rows = measure_mat.nrow(), m->cols = measure_mat.ncol();
  double currMass;
  
  for (int i = 0; i < measure_mat.nrow(); ++i) {
  for (int j = 0; j < measure_mat.ncol(); ++j) {
    currMass = measure_mat(i,j);
    m->masses.push_back(currMass);
    m->accMass += currMass;
  }
  }
  
  return m;
}

GeneralMeasureData* NumericMatrix2GenMeasureData(NumericMatrix measure_mat, 
                                                 double normalizingFactor) {
  GeneralMeasureData* m = new GeneralMeasureData;
  m->accMass = 0;

  for (int i = 0; i < measure_mat.nrow(); ++i) {
    m->points.push_back(std::make_pair(measure_mat(i,0), measure_mat(i,1)));
    m->masses.push_back(measure_mat(i,2));
    m->accMass += measure_mat(i,2);
  }

  if (normalizingFactor != 1) {
  // normalize the points to be contained in [0,1]x[0,1]
  Rcout << "nu is internally normalized" << std::endl;
  for (std::vector<std::pair<double, double> >::iterator it = m->points.begin();
    it != m->points.end(); ++it) {
    *it = std::make_pair(it->first / normalizingFactor,
      it->second / normalizingFactor);
  }
  }
   
  return m;
}

void normalizeVector(lbfgsfloatval_t* vector, int numEntries) {
  lbfgsfloatval_t min = vector[0];
  for (int i = 1; i < numEntries; ++i) {
    min = std::min(min, vector[i]);
  }
  for (int i = 0; i < numEntries; ++i) {
    vector[i] -= min;
  }
}

#endif
