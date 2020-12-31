#include"Interfaces-CPLEX.h"


/////////////////////////////////////////////////////////////////////////////////////////////
// TSolverInterfaceCPLEX Methods
/////////////////////////////////////////////////////////////////////////////////////////////


TSolverInterfaceCPLEX::TSolverInterfaceCPLEX(TCouplingHandlerExtBase *_couplingHandler,
		TCPLEXNetSolverBase *_solver,
		bool _initializeBases, bool _deleteSolverOnDestroy) {
	couplingHandler=_couplingHandler;
	solver = _solver;

	initializeBases=_initializeBases;
	deleteSolverOnDestroy=_deleteSolverOnDestroy;
	
	oldBasis=NULL;
}

TSolverInterfaceCPLEX::TSolverInterfaceCPLEX(TCouplingHandlerExtBase *_couplingHandler,
		TCPLEXNetSolverBase *_solver,
		 double *_alpha, double *_beta,
			bool _initializeBases, bool _deleteSolverOnDestroy) {
	alpha=_alpha;
	beta=_beta;
	couplingHandler=_couplingHandler;
	solver = _solver;

	initializeBases=_initializeBases;
	deleteSolverOnDestroy=_deleteSolverOnDestroy;
}

TSolverInterfaceCPLEX::~TSolverInterfaceCPLEX() {
	if(deleteSolverOnDestroy && (solver!=NULL)) {
		delete solver;
	}
	
	if(oldBasis!=NULL) {
		delete oldBasis;
	}
}

int TSolverInterfaceCPLEX::solve() {
	int msg;
	msg=solver->solve();
	return msg;
}

int TSolverInterfaceCPLEX::extractBasisVarList(TVarListHandler **basis) {
	// pointer to last varList for convenience
	TVarListHandler *xVars=couplingHandler->getXVars();
	
	// create new varListHandler to store basis
	TVarListHandler *result=new TVarListHandler();
	result->setupEmpty(solver->xres);
	// varlist signal on last varlist to extract basis info
	TVarListSignal<int> *arcStatus=new TVarListSignal<int>(xVars,0);
	// extracting
	int msg;
	msg=solver->extractBase(arcStatus->signal,NULL);
	if(msg!=0) {
		return msg;
	}
	// iterate over basis info
	TVarListHandler::TIterator it=xVars->iterationInitialize();
	while(xVars->iterate(&it)) {
		// if variable was part of basis
		if(arcStatus->signal[it.offset]==CPX_BASIC) {
			// add to basis, ignore duplicate check
			result->addToLine(it.x,it.y,false);
		}
	}
	// cleanup
	delete arcStatus;
	// write result to return variable
	*basis=result;
	
	return 0;
}

int TSolverInterfaceCPLEX::prepareUpdate(__attribute__((unused)) TVarListHandler *newXVars) {
	if(initializeBases) {
	
		// code that can be replaced by extractBasisVarList
			
//		// pointer to last varList for convenience
//		TVarListHandler *xVars=couplingHandler->getXVars();
//		// var list that stores basis variables
//		oldBasis=new TVarListHandler();
//		oldBasis->setupEmpty(solver->xres);
//		// varlist signal on last varlist to extract basis info
//		TVarListSignal<int> *arcStatus=new TVarListSignal<int>(xVars,0);
//		// extracting
//		int msg;
//		msg=solver->extractBase(arcStatus->signal,NULL);
//		if(msg!=0) {
//			return msg;
//		}
//		// iterate over basis info
//		TVarListHandler::TIterator it=xVars->iterationInitialize();
//		while(xVars->iterate(&it)) {
//			// if variable was part of basis
//			if(arcStatus->signal[it.offset]==CPX_BASIC) {
//				// add to basis, ignore duplicate check
//				oldBasis->addToLine(it.x,it.y,false);
//				// add to newXVars, with duplicate check
//				newXVars->addToLine(it.x,it.y,true);
//			}
//		}
//		// cleanup
//		delete arcStatus;

		// extract basis as varList and write to oldBasis
		extractBasisVarList(&oldBasis);
		// merge oldBasis into newXVars
		newXVars->merge(oldBasis);


	}
	return 0;
}

int TSolverInterfaceCPLEX::update() {
	int msg;
	msg=solver->deleteArcs();
	if(msg!=0) {
		return msg;
	}
	msg=solver->setupArcs();
	if(msg!=0) {
		return msg;
	}

	if(initializeBases) {
		// initialize pseudo node status
		// observed from cplex: always first pseudo arc is chosen to be part of basis
		int *nodeStatus=(int*) calloc((solver->xres+solver->yres),sizeof(int));
		nodeStatus[0]=1;
		
		TVarListSignal<int> *newArcStatus=new TVarListSignal<int>(couplingHandler->getXVars(),0);
		// need this to be able to write to newArcStatus by absolute indices
		newArcStatus->computeOffsets();
		// iterate over old basis
		TVarListHandler::TIterator it=oldBasis->iterationInitialize();
		while(oldBasis->iterate(&it)) {
			// set all corresponding entries of newArcStatus to CPX_BASIC
			newArcStatus->write(it.x,it.y,CPX_BASIC);
		}
		msg=solver->setupBase(newArcStatus->signal,nodeStatus);
		free(nodeStatus);
		delete newArcStatus;
		delete oldBasis;
		oldBasis=NULL;
		if(msg!=0) {
				return msg;
		}
	}

	return 0;
}


double TSolverInterfaceCPLEX::getObjective() {
	return solver->getObjective();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Factory Class for MultiScaleSolver
////////////////////////////////////////////////////////////////////////////////////////////////////////

TFactorySolverInterfaceCPLEX::TFactorySolverInterfaceCPLEX() {
}


int TFactorySolverInterfaceCPLEX::generate(
		__attribute__((unused)) int layerId,
		TCouplingHandlerSparse *couplingHandler,
		TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface,
		double *muX, double *muY, double *alpha, double *beta,
		TSolverInterface **result) {
		
	int msg;
	
	TCPLEXNetSolver<TCouplingHandlerSparse> *subSolver;
	subSolver=new TCPLEXNetSolver<TCouplingHandlerSparse>(couplingHandler,
			muX,muY,alpha,beta);
			
	
	TSolverInterfaceCPLEX *solverInterface;
	solverInterface= new TSolverInterfaceCPLEX(
			couplingHandlerInterface, subSolver,
			alpha, beta,
			true, true
			);
	
	// initialize sub solver
	msg=subSolver->setup();
	if(msg!=0) {
		return msg;
	}
	

	*result=solverInterface;
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Derivation of MultiScaleSetup Class
////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class TReferenceSetup>
TMultiScaleSetupCPLEX<TReferenceSetup>::~TMultiScaleSetupCPLEX() {
	if(FactorySolverInterface!=NULL) {
		delete FactorySolverInterface;
	}
}



template<class TReferenceSetup>
int TMultiScaleSetupCPLEX<TReferenceSetup>::SetupSolverSpecificPrepropcessing() {
	// flip sign of muY
	// skip this for now. will be done inside solver
	//for(int y=0;y<yres;y++) {
	//	muY[y]=-muY[y];
	//}

	return 0;
}


template<class TReferenceSetup>
int TMultiScaleSetupCPLEX<TReferenceSetup>::SetupFactorySolverInterface() {
	FactorySolverInterface = new TFactorySolverInterfaceCPLEX;
	return 0;
}



// instantantiate solver class templates

template class TCPLEXNetSolver<TCouplingHandlerSemiDense>;
template class TCPLEXNetSolver<TCouplingHandlerSparse>;


template class TMultiScaleSetupCPLEX<TMultiScaleSetupW2Grid>;

