#include <cstdlib>
#include <iostream>
#include <Rcpp.h>

#include"ShortCutSolver/MultiScaleSolver.h"
#ifdef WITH_CPLEX
  #include "ShortCutSolver_CPLEX/Interfaces-CPLEX.h"
#endif
#ifndef WITH_CPLEX
  #include "ShortCutSolver_SparseSimplex/Interfaces-SparseSimplex.h"
#endif



using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int* reverseArray(int *data, int len) {
	int *result=(int*) malloc(sizeof(int)*len);
	for(int i=0;i<len;i++) {
		result[len-1-i]=data[i];
	}
	return result;
}

// [[Rcpp::export]]
SEXP SolveHierarchicalTransport(SEXP x, SEXP y, SEXP xydepth, SEXP xydimensions, SEXP compdepth,
				SEXP measureScaleVecPre, SEXP keepBasisVecPre, SEXP refineBasisVecPre,
				SEXP layerCoarsestVecPre, SEXP verboseVecPre,
		                SEXP assignment, SEXP udummy, SEXP vdummy) {
	
  // much of what follows is still in the old C / .C style
  // just emulated with Rcpp, but we try to treat the new retbasis IntegerVector
  // Rcpp-style, i.e. we do not pass it from R
        Rcpp::NumericVector dataX(x);
	Rcpp::NumericVector dataY(y);
	Rcpp::IntegerVector datadepth(xydepth);
     	Rcpp::IntegerVector datadim(xydimensions);
        Rcpp::IntegerVector cdepthvec(compdepth);
        Rcpp::NumericVector measureScaleVec(measureScaleVecPre);
        Rcpp::IntegerVector keepBasisVec(keepBasisVecPre);
        Rcpp::IntegerVector refineBasisVec(refineBasisVecPre);
        Rcpp::IntegerVector coupling(assignment);
	Rcpp::NumericVector u(udummy);
	Rcpp::NumericVector v(vdummy);
	Rcpp::IntegerVector layerCoarsestVec(layerCoarsestVecPre);
	Rcpp::LogicalVector verboseVec(verboseVecPre);
	
	TDoubleMatrix muxx = {dataX.begin(), *datadepth.begin(), datadim.begin()};
	TDoubleMatrix muyy = {dataY.begin(), *datadepth.begin(), datadim.begin()};
	TDoubleMatrix muxxRev={muxx.data, muxx.depth, reverseArray(muxx.dimensions,muxx.depth)};
	TDoubleMatrix muyyRev={muyy.data, muyy.depth, reverseArray(muyy.dimensions,muyy.depth)};
	TDoubleMatrix *muX = &muxxRev, *muY = &muyyRev;
	

	
	int cdepth;
	int keepBasis;
	int refineBasis;
	double measureScale;
	int layerCoarsest;
	bool verbose;
	
	cdepth = *cdepthvec.begin();
	keepBasis=*keepBasisVec.begin();
	refineBasis=*refineBasisVec.begin();
	measureScale=*measureScaleVec.begin();
        layerCoarsest=*layerCoarsestVec.begin();
	verbose=*verboseVec.begin();
	
	int msg;

	// new things
	// int layerCoarsest=1; // at which refinement level does multiScale solving begin?
	verbose_mode=verbose; // regulating verbosity
	

	eprintf("entered c++ code\n");
	eprintf("constructing setup object\n");
	#ifdef WITH_CPLEX
		eprintf("\tusing cplex\n");
		TMultiScaleSetupCPLEX<TMultiScaleSetupW2Grid> MultiScaleSetup(muX,muY,cdepth);
	#else
		eprintf("\tusing custom simplex\n");
		TMultiScaleSetupSparseSimplex<TMultiScaleSetupW2Grid> MultiScaleSetup({muX,muY,cdepth},measureScale,keepBasis,refineBasis);
	#endif
	
	msg=MultiScaleSetup.Setup();
	if(msg!=0) {
		eprintf("error during setup: %d\n",msg);
	        return Rcpp::wrap(Rcpp::List::create(msg));
	}

	eprintf("initializing multiscale solver\n");
	TMultiScaleSolver MultiScaleSolver(
			MultiScaleSetup.FactoryCostFunctionProvider,
			MultiScaleSetup.FactoryCouplingHandlerExt,
			MultiScaleSetup.FactorySolverInterface,
			MultiScaleSetup.FactoryShieldGenerator,
			MultiScaleSetup.HPX, MultiScaleSetup.HPY,
			MultiScaleSetup.muXH, MultiScaleSetup.muYH,
			layerCoarsest,
			//TShortCutSolver::VCHECK_DUAL
			TShortCutSolver::VCHECK_PRIMAL
			);
	
	// delay destruction of shortcut solver objects until destruction of MultiScaleSolver
	// required, to extract basis after solving is complete
	MultiScaleSolver.autoDeletePointers=false;
	
	eprintf("solving\n");
	msg=MultiScaleSolver.solve();
	if(msg!=0) {
		eprintf("error during solving: %d\n",msg);
	        return Rcpp::wrap(Rcpp::List::create(msg));
	}
	eprintf("solved\n");
	
	eprintf("exporting result\n");

	int xres=MultiScaleSetup.xres;
	int yres=MultiScaleSetup.yres;
	
	int xit,yit;
	for(xit=0;xit<xres;xit++) {
		u[xit]=MultiScaleSolver.alpha[xit];
	}
	for(yit=0;yit<yres;yit++) {
		v[yit]=MultiScaleSolver.beta[yit];
	}
	TVarListHandler *xVarsFinal=MultiScaleSolver.xVarsFinal;
	double *muFinal=MultiScaleSolver.muFinal;
	
	// copy optimal coupling to "coupling" array, to return to R
	TVarListHandler::TIterator it=xVarsFinal->iterationInitialize();
	while(xVarsFinal->iterate(&it)) {
		// TODO: rounding of double CPLEX couplings to integer precision may be problematic
		coupling[it.y*xres+it.x]=(int) round(muFinal[it.offset]);
	}
	
	///////////////////////////////////////////////////////////////////////////////////
	// if desired: basis extraction
	eprintf("basis extraction\n");
	
	TVarListHandler *basis;
	#ifdef WITH_CPLEX
		// CPLEX
		TSolverInterfaceCPLEX *solverInterface=(TSolverInterfaceCPLEX*) MultiScaleSolver.solverInterface;
		solverInterface->extractBasisVarList(&basis);
	#else
		// Custom simplex
		TSolverInterfaceSparseSimplex *solverInterface=(TSolverInterfaceSparseSimplex*) MultiScaleSolver.solverInterface;
		solverInterface->extractBasisVarList(&basis);
	#endif

	Rcpp::NumericMatrix retbasis(basis->total,2);
	it=basis->iterationInitialize();
	int i=0;
	while(basis->iterate(&it)) {
	  retbasis(i,0) = it.x+1;  // adapt to R-indices
	  retbasis(i,1) = it.y+1;  // adapt to R-indices
	  // eprintf("%d %d\n",it.x,it.y);
	  i++;
	}
	eprintf("%d basis vectors extracted\n",basis->total);

	// print basis
	// eprintf("printing basis\n");
	// eprintf("%d\n",basis->total);
	// TVarListHandler::TIterator it=basis->iterationInitialize();
	// it=basis->iterationInitialize();
	// while(basis->iterate(&it)) {
	//	eprintf("%d %d\n",it.x,it.y);
	// }

	delete basis;

	// end basis extraction
	////////////////////////////////////////////////////////////////////


	

	// clean up reversed dimension arrays
	free(muxxRev.dimensions);
	free(muyyRev.dimensions);
 
        return Rcpp::wrap(Rcpp::List::create(0,x,y,coupling,retbasis,u,v));
	
}

