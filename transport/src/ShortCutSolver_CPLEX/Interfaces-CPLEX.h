#ifndef Interfaces_CPLEX_H_
#define Interfaces_CPLEX_H_

#include<cstdlib>

#include<OT_CPLEX/TCPLEXNetSolver.h>
#include<ShortCutSolver/Interfaces.h>
#include<ShortCutSolver/MultiScaleSolver.h>

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Short Cut Solver Interface
////////////////////////////////////////////////////////////////////////////////////////////////////////


class TSolverInterfaceCPLEX: public TSolverInterface {
public:
	TCPLEXNetSolverBase *solver;
	TCouplingHandlerExtBase *couplingHandler;
	bool deleteSolverOnDestroy;
	// whether or not to destroy the solver on destroy
	bool initializeBases;
	// initializeBases=true: optimal bases are extracted from CPLEX NET solver and
	//    re-sent to solver (on the variables that are kept) before re-solving.
	TVarListHandler *oldBasis;
	TSolverInterfaceCPLEX(TCouplingHandlerExtBase *_couplingHandler,
			TCPLEXNetSolverBase *_solver,
			bool _initializeBases, bool _deleteSolverOnDestroy=false);
	TSolverInterfaceCPLEX(TCouplingHandlerExtBase *_couplingHandler,
			TCPLEXNetSolverBase *_solver, double *_alpha, double *_beta,
			bool _initializeBases, bool _deleteSolverOnDestroy=false);
	~TSolverInterfaceCPLEX();

	int solve();
	int prepareUpdate(__attribute__((unused)) TVarListHandler *newXVars);
	int update();

	double getObjective();

	int extractBasisVarList(TVarListHandler **basis);

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Factory Class for MultiScaleSolver
////////////////////////////////////////////////////////////////////////////////////////////////////////


class TFactorySolverInterfaceCPLEX : public TFactorySolverInterfaceBase {
public:
	TFactorySolverInterfaceCPLEX();

	virtual int generate(
			__attribute__((unused)) int layerId,
			TCouplingHandlerSparse *couplingHandler,
			TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface,
			double *muX, double *muY, double *alpha, double *beta,
			TSolverInterface **result);
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Derivation of MultiScaleSetup Class
////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class TReferenceSetup>
class TMultiScaleSetupCPLEX : public TReferenceSetup {
public:
	using TReferenceSetup::muY;
	using TReferenceSetup::yres;

	using TReferenceSetup::FactorySolverInterface;
	using TReferenceSetup::TReferenceSetup;
	
	~TMultiScaleSetupCPLEX();
	
	//TMultiScaleSetupCPLEX(typename TReferenceSetup::ConstructorArgument arg) : TReferenceSetup(arg) {};
	virtual int SetupSolverSpecificPrepropcessing();
	virtual int SetupFactorySolverInterface();
};


#endif
