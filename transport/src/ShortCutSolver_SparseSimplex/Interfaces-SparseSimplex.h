#ifndef Interfaces_SparseSimplex_H_
#define Interfaces_SparseSimplex_H_

#include<cstdlib>

#include<OT_SparseSimplex/TSparseSimplexSolver.h>
#include<ShortCutSolver/Interfaces.h>
#include<ShortCutSolver/MultiScaleSolver.h>

using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Short Cut Solver Interface
////////////////////////////////////////////////////////////////////////////////////////////////////////


class TSolverInterfaceSparseSimplex: public TSolverInterface {
public:
	TSparseSimplexSolverBase *solver;
	TCouplingHandlerExtBase *couplingHandler;
	bool copyBasis;
	bool deleteSolverOnDestroy;
	// whether or not to destroy the solver on destroy
	TSolverInterfaceSparseSimplex(TCouplingHandlerExtBase *_couplingHandler,
			TSparseSimplexSolverBase *_solver, double *_alpha, double *_beta,
			bool _copyBasis,
			bool _deleteSolverOnDestroy);
	~TSolverInterfaceSparseSimplex();

	int solve();
	int prepareUpdate(TVarListHandler *newXVars);
	int update();

	double getObjective();

	int extractBasisVarList(TVarListHandler **basis);

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Factory Class for MultiScaleSolver
////////////////////////////////////////////////////////////////////////////////////////////////////////


class TFactorySolverInterfaceSparseSimplex : public TFactorySolverInterfaceBase {
public:
	bool keepBasis; // whether to recycle basis during multiple solves on same hierarchy scale
	bool refineBasis; // whether to construct refined basis from old coarse basis during hierarchy scale refinement

	// reference to a shieldGenerator factory class, in case the "initial shielding" is required at each hierarchy level
	// (if refineBasis=false)
	TFactoryShieldGeneratorBase *FactoryShieldGenerator;

	
	// details about hierarchical problem structure
	// is required for basis refinement
	THierarchicalPartition *HPX;
	THierarchicalPartition *HPY;
	double **muXH;
	double **muYH;

	
	bool storedOldBasis; // whether a coarse old basis is currently stored
	// variables that store old basis
	TVarListHandler *xVarsC;
	double *piC;
	bool *basisC;
	// variables that store refined basis
	TVarListHandler *xVarsF;
	bool *basisF;
	double *piF;

	
	TFactorySolverInterfaceSparseSimplex(TFactoryShieldGeneratorBase *_FactoryShieldGenerator,
			bool _keepBasis, bool _refineBasis,
			THierarchicalPartition *_HPX, THierarchicalPartition *_HPY,
			double **_muXH, double **_muYH);
	virtual int generate(
			__attribute__((unused)) int layerId,
			TCouplingHandlerSparse *couplingHandler,
			TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface,
			double *muX, double *muY, double *alpha, double *beta,
			TSolverInterface **result);
	
	virtual int prepareRefinement(
			__attribute__((unused)) int layerId,
			__attribute__((unused)) TSolverInterface *solverInterface
			);
	virtual int customizeRefinement(
			__attribute__((unused)) int layerId,
			__attribute__((unused)) TVarListHandler *xVars
			);
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Derivation of MultiScaleSetup Class
////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class TReferenceSetup>
class TMultiScaleSetupSparseSimplex : public TReferenceSetup {
public:
	using TReferenceSetup::muY;
	using TReferenceSetup::muX;
	using TReferenceSetup::xres;
	using TReferenceSetup::yres;
	using TReferenceSetup::HPX;
	using TReferenceSetup::HPY;
	using TReferenceSetup::muXH;
	using TReferenceSetup::muYH;
	using TReferenceSetup::FactoryShieldGenerator;
	using TReferenceSetup::FactorySolverInterface;
	using TReferenceSetup::TReferenceSetup;
	double measureScale;
	int keepBasis;
	int refineBasis;
	
	TMultiScaleSetupSparseSimplex(typename TReferenceSetup::ConstructorArgument arg,
			double _measureScale, int _keepBasis, int _refineBasis) : TReferenceSetup(arg) {
		measureScale=_measureScale;
		keepBasis=_keepBasis;
		refineBasis=_refineBasis;
	};
	
	~TMultiScaleSetupSparseSimplex();
	
	virtual int SetupSolverSpecificPrepropcessing();
	virtual int SetupFactorySolverInterface();
};



#endif
