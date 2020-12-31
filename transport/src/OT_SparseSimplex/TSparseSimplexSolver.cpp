#include"TSparseSimplexSolver.h"



TSparseSimplexSolverBase::TSparseSimplexSolverBase(int *_muX, int *_muY, double *_alpha, double *_beta, bool _deleteMarginals) {
	xres=0;
	yres=0;

	muX=_muX;
	muY=_muY;
	
	alpha=_alpha;
	beta=_beta;

	objective=0.;
	
	assignment=NULL;
	basis=NULL;

	basisstartgiven=0;

	setupStatus=false;
	solutionStatus=false;
	
	deleteMarginals=_deleteMarginals;
}


int TSparseSimplexSolverBase::solve() {
	return MSG_NOT_IMPLEMENTED;
}

int TSparseSimplexSolverBase::setup() {
	if(setupStatus) {
		return 0;
	}
	
	// create empty assignment and basis matrices
	assignment=(int*) malloc(xres*yres*sizeof(int));
	basis=(int*) malloc(xres*yres*sizeof(int));
	
	// initialize with 0
	setupBaseZero();
	// indicate that basis does not contain feasible basis
	basisstartgiven=0;
	
	setupStatus=true;
	return 0;
}

int TSparseSimplexSolverBase::setupBaseZero() {
	// initialize with 0
	for(int y=0;y<yres;y++) {
		for(int x=0;x<xres;x++) {
			assignment[y*xres+x]=0;
			basis[y*xres+x]=0;
		}
	}
	return 0;
}

int TSparseSimplexSolverBase::cleanup() {
	if(!setupStatus) {
		return 0;
	}
	
	if(assignment!=NULL) {
		free(assignment);
	}
	if(basis!=NULL) {
		free(basis);
	}
	setupStatus=false;
	return 0;
}
	


TSparseSimplexSolverBase::~TSparseSimplexSolverBase() {
	if(deleteMarginals) {
		free(muX);
		free(muY);
	}
	cleanup();
}


double TSparseSimplexSolverBase::getObjective() {
	return objective;
}
