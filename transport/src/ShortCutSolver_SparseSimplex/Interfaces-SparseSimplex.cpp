#include"Interfaces-SparseSimplex.h"


/////////////////////////////////////////////////////////////////////////////////////////////
// TSolverInterfaceSparseSimplex Methods
/////////////////////////////////////////////////////////////////////////////////////////////


TSolverInterfaceSparseSimplex::TSolverInterfaceSparseSimplex(TCouplingHandlerExtBase *_couplingHandler,
		TSparseSimplexSolverBase *_solver, double *_alpha, double *_beta,
		bool _copyBasis,
		bool _deleteSolverOnDestroy) {
		
	alpha=_alpha;
	beta=_beta;
	couplingHandler=_couplingHandler;
	solver = _solver;

	copyBasis=_copyBasis;
	deleteSolverOnDestroy=_deleteSolverOnDestroy;
	
}


TSolverInterfaceSparseSimplex::~TSolverInterfaceSparseSimplex() {
	if(deleteSolverOnDestroy && (solver!=NULL)) {
		delete solver;
	}
}

int TSolverInterfaceSparseSimplex::solve() {
	int msg;
	msg=solver->solve();
	return msg;
}

int TSolverInterfaceSparseSimplex::extractBasisVarList(TVarListHandler **basis) {
	// pointer to last varList for convenience
	TVarListHandler *xVars=couplingHandler->getXVars();
	// xres variable for convenience
	int xres=solver->xres;
	
	// create new varListHandler to store basis
	TVarListHandler *result=new TVarListHandler();
	result->setupEmpty(xres);


	// iterate over varList
	TVarListHandler::TIterator it=xVars->iterationInitialize();
	while(xVars->iterate(&it)) {
		// test if basis entry is 1
		if(solver->basis[it.y*xres+it.x]==1) {
			// if yes: add to basis list (can ignore duplicate check)
			result->addToLine(it.x,it.y,false);
		}
	}
	
	// write result to return variable
	*basis=result;
	
	return 0;
}


int TSolverInterfaceSparseSimplex::prepareUpdate(TVarListHandler *newXVars) {
	// go through var list and create new var list that captures basis position
	if (copyBasis) {
		int xres=couplingHandler->getXres();
		//int yres=couplingHandler->getYres();
		TVarListHandler *xVars=couplingHandler->getXVars();
		eprintf("\t\tpreparing next update: add old basis to newXVars\n");
		
		// iterate over varList
		TVarListHandler::TIterator it=xVars->iterationInitialize();
		while(xVars->iterate(&it)) {
			// test if basis entry is 1
			if(solver->basis[it.y*xres+it.x]==1) {
				// if yes: add old basis entries to newXVars
				newXVars->addToLine(it.x,it.y,true);
			}
		}

	}
	return 0;
}

int TSolverInterfaceSparseSimplex::update() {
	if (copyBasis) {
		// communicate to solver object that basis array contains valid basis.
		solver->basisstartgiven=1;
	}
	return 0;
}

double TSolverInterfaceSparseSimplex::getObjective() {
	return solver->getObjective();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Factory Class for MultiScaleSolver
////////////////////////////////////////////////////////////////////////////////////////////////////////


TFactorySolverInterfaceSparseSimplex::TFactorySolverInterfaceSparseSimplex(
		TFactoryShieldGeneratorBase *_FactoryShieldGenerator,
		bool _keepBasis, bool _refineBasis,
		THierarchicalPartition *_HPX, THierarchicalPartition *_HPY,
		double **_muXH, double **_muYH) {
		
	FactoryShieldGenerator=_FactoryShieldGenerator;
	keepBasis=_keepBasis;
	refineBasis=_refineBasis;

	HPX=_HPX;
	HPY=_HPY;
	muXH=_muXH;
	muYH=_muYH;
	
	storedOldBasis=false;
	xVarsC=NULL;
	basisC=NULL;
	piC=NULL;
	xVarsF=NULL;
	basisF=NULL;
	piF=NULL;
}


int TFactorySolverInterfaceSparseSimplex::generate(
		__attribute__((unused)) int layerId,
		TCouplingHandlerSparse *couplingHandler,
		TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface,
		double *muX, double *muY, double *alpha, double *beta,
		TSolverInterface **result) {
		
	int msg;
//	TVarListHandler *xVarsF=NULL;
//	bool *basisF=NULL;
//	double *piF=NULL;

//	/////////////////////////////////////////////////////////////////////////////////////////////////
//	if(!storedOldBasis) {
//		// shield initially given set of variables
//		eprintf("\t\tfirst shielding.\n\t\t\ttotal variables: %d\n",couplingHandlerInterface->getXVars()->total);
//		TVarListHandler *newXVars;
//		newXVars=new TVarListHandler(couplingHandlerInterface->getXVars());
//	
//		TShieldGeneratorBase *shieldGenerator;
//		msg=FactoryShieldGenerator->generate(layerId,&shieldGenerator);
//		if(msg!=0) {
//			return msg;
//		}
//		shieldGenerator->generateShield(newXVars,newXVars);
//		delete shieldGenerator;
//		newXVars->sort();
//		couplingHandler->updateXVars(newXVars,false);
//		eprintf("\t\t\tnew total variables: %d\n",couplingHandlerInterface->getXVars()->total);
//	/////////////////////////////////////////////////////////////////////////////////////////////////
//	} else {
//	/////////////////////////////////////////////////////////////////////////////////////////////////
//		eprintf("\t\trefining coarse basis.\n");
//		
//		int layerC=layerId-1; // coarse level id
//		// refine var list		
//		xVarsF=refineVarList(
//				HPX, HPY,
//				xVarsC, layerC, true);
//		// compute refined basis
//		msg=MultiScaleRefineBasis(
//				HPX, HPY,
//				xVarsC,
//				basisC, piC,
//				muXH[layerId], muYH[layerId], xVarsF,
//				layerC,
//				&basisF,&piF);
//		if(msg!=0) {
//			return msg;
//		};
//		// merging basis: copy current var list, then merge refined var list depending on basisF
//		eprintf("\t\t\tbefore merging. total variables: %d\n",couplingHandlerInterface->getXVars()->total);
//		TVarListHandler *newXVars;
//		newXVars=new TVarListHandler(couplingHandlerInterface->getXVars());
//		newXVars->mergeSelected(xVarsF,basisF);		
//		newXVars->sort();
//		couplingHandler->updateXVars(newXVars,false);
//		eprintf("\t\t\tafter merging: %d\n",couplingHandlerInterface->getXVars()->total);		
//	}
//	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// create integer marginals
	int *muXint=(int*) malloc(couplingHandler->xres*sizeof(int));
	int *muYint=(int*) malloc(couplingHandler->yres*sizeof(int));
	// sums of marginals to make consistency check
	int muXintSum=0, muYintSum=0;

	// assign truncated marginals
	const double muScale=1;
	for(int x=0;x<couplingHandler->xres;x++) {
		muXint[x]=(int)(round(muX[x]/muScale));
		muXintSum+=muXint[x];
	}
	for(int y=0;y<couplingHandler->yres;y++) {
		muYint[y]=(int)(round(muY[y]/muScale));
		muYintSum+=muYint[y];
	}
	// check if marginals have same mass
	if(muXintSum!=muYintSum) {
		eprintf("ERROR: marginals have different mass after truncation: %d vs %d\n",muXintSum,muYintSum);
		return ERR_SPARSESIMPLEX_TRUNCATION;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////////////////////
	// create actual subsolver and interface objects
	TSparseSimplexSolver<TCouplingHandlerSparse> *subSolver;
	subSolver=new TSparseSimplexSolver<TCouplingHandlerSparse>(couplingHandler,
				muXint,muYint, alpha, beta, true);

	TSolverInterfaceSparseSimplex *solverInterface;
	solverInterface= new TSolverInterfaceSparseSimplex(
			couplingHandlerInterface, subSolver,
			alpha, beta,
			keepBasis,true
			);
		

	// initialize sub solver
	msg=subSolver->setup();
	if(msg!=0) {
		return msg;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////////////////////
	if(storedOldBasis) {
		// copy refined basis into basis array of subSolver
		eprintf("\t\tcopying basis\n");
		
		int xres=xVarsF->res;
		TVarListHandler::TIterator it=xVarsF->iterationInitialize();
		while(xVarsF->iterate(&it)) {
			if(basisF[it.offset]) {
				// if this entry is part of refined basis

				// write to basis array
				subSolver->basis[it.y*xres+it.x]=1;
				// write to assignment
				subSolver->assignment[it.y*xres+it.x]=(int)(round(piF[it.offset]/muScale));
				
			}		
		}
		
		// communicate to solver object that basis array contains valid basis.
		subSolver->basisstartgiven=1;

		
		// clean up


		storedOldBasis=false;
		delete xVarsC;
		free(basisC);
		free(piC);

		delete xVarsF;
		free(basisF);
		free(piF);
	}
	
	*result=solverInterface;
	return 0;

}


int TFactorySolverInterfaceSparseSimplex::prepareRefinement(
		__attribute__((unused)) int layerId,
		__attribute__((unused)) TSolverInterface *solverInterface
		) {
	
	if(!refineBasis) {
		// if refineBasis is not activated, do nothing
		return 0;
	}
	
	
	eprintf("\t\textract coarse basis and coupling\n");
	TSolverInterfaceSparseSimplex *solverInterfaceSpSi=(TSolverInterfaceSparseSimplex*) solverInterface;
	TSparseSimplexSolver<TCouplingHandlerSparse> *solver=(TSparseSimplexSolver<TCouplingHandlerSparse>*) solverInterfaceSpSi->solver;
	
	xVarsC=new TVarListHandler(solverInterfaceSpSi->couplingHandler->getXVars());
	piC=(double*) malloc(sizeof(double)*xVarsC->total);
	basisC=(bool*) malloc(sizeof(bool)*xVarsC->total);
	
	int xres=xVarsC->res;
	TVarListHandler::TIterator it=xVarsC->iterationInitialize();
	while(xVarsC->iterate(&it)) {
		basisC[it.offset]=(bool) solver->basis[it.y*xres+it.x]; // cast: int -> bool
		piC[it.offset]=(double) solver->assignment[it.y*xres+it.x]; // cast: int -> double
	}
	// indicate that basis was extracted
	storedOldBasis=true;
	return 0;
}


int TFactorySolverInterfaceSparseSimplex::customizeRefinement(int layerId, TVarListHandler *xVars) {
	int msg;
	/////////////////////////////////////////////////////////////////////////////////////////////////
	if(!storedOldBasis) {
		// shield initially given set of variables
		eprintf("\t\tfirst shielding.\n\t\t\ttotal variables: %d\n",xVars->total);

		TShieldGeneratorBase *shieldGenerator;
		msg=FactoryShieldGenerator->generate(layerId,&shieldGenerator);
		if(msg!=0) {
			return msg;
		}
		shieldGenerator->generateShield(xVars,xVars);
		delete shieldGenerator;
		xVars->sort();
		eprintf("\t\t\tnew total variables: %d\n",xVars->total);
	/////////////////////////////////////////////////////////////////////////////////////////////////
	} else {
	/////////////////////////////////////////////////////////////////////////////////////////////////
		eprintf("\t\trefining coarse basis.\n");
		
		int layerC=layerId-1; // coarse level id
		// refine var list		
		xVarsF=refineVarList(
				HPX, HPY,
				xVarsC, layerC, true);
		// compute refined basis
		msg=MultiScaleRefineBasis(
				HPX, HPY,
				xVarsC,
				basisC, piC,
				muXH[layerId], muYH[layerId], xVarsF,
				layerC,
				&basisF,&piF);
		if(msg!=0) {
			return msg;
		};
		// merging basis: copy current var list, then merge refined var list depending on basisF
		eprintf("\t\t\tbefore merging. total variables: %d\n",xVars->total);
		xVars->mergeSelected(xVarsF,basisF);		
		xVars->sort();
		eprintf("\t\t\tafter merging: %d\n",xVars->total);		
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Derivation of MultiScaleSetup Class
////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TReferenceSetup>
TMultiScaleSetupSparseSimplex<TReferenceSetup>::~TMultiScaleSetupSparseSimplex() {
	if(FactorySolverInterface!=NULL) {
		delete FactorySolverInterface;
	}
}


template<class TReferenceSetup>
int TMultiScaleSetupSparseSimplex<TReferenceSetup>::SetupSolverSpecificPrepropcessing() {
	// truncate measures to measureScale
	return MeasureToolsTruncateMeasures(muX,muY,xres,yres,measureScale);	
	return 0;
}


template<class TReferenceSetup>
int TMultiScaleSetupSparseSimplex<TReferenceSetup>::SetupFactorySolverInterface() {
	FactorySolverInterface = new TFactorySolverInterfaceSparseSimplex(FactoryShieldGenerator,keepBasis,refineBasis,HPX,HPY,muXH,muYH);
	return 0;
}


// instantantiate solver class templates

template class TSparseSimplexSolver<TCouplingHandlerSemiDense>;
template class TSparseSimplexSolver<TCouplingHandlerSparse>;

template class TMultiScaleSetupSparseSimplex<TMultiScaleSetupW2Grid>;

