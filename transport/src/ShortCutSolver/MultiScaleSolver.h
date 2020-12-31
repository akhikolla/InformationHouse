#ifndef MultiScaleSolver_H_
#define MultiScaleSolver_H_

#include<Common/PythonTypes.h>
#include<Common/ErrorCodes.h>
#include<Common/tools.h>
#include<Common/Verbose.h>

#include<Common/GridTools.h>
#include<Common/THierarchyBuilder.h>
#include<Common/THierarchicalPartition.h>
#include<Common/TCostFunctionProvider.h>
#include<Common/TCostFunctionProvider-Dynamic.h>
#include<ShortCutSolver/Interfaces.h>
#include<ShortCutSolver/TShieldGenerator.h>
#include<ShortCutSolver/TShortCutSolver.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Factory Classes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// First define all sorts of factory classes for cost function providers, coupling handlers and so forth

/////////////////////////////
// CostFunctionProvider
/////////////////////////////

// Base
class TFactoryCostFunctionProviderBase {
public:
	virtual ~TFactoryCostFunctionProviderBase() {};
	virtual int generate(__attribute__((unused)) int layerId, TCostFunctionProviderBase **result) {
		*result=NULL;
		return ERR_BASE_NOTIMPLEMENTED;
	}
};

// Squared Euclidean Distance Dynamic
class TFactoryCostFunctionProvider_Dynamic : public TFactoryCostFunctionProviderBase {
public:
	int *xresH, *yresH;
	double **xposH, **yposH;
	int dim;
	TFactoryCostFunctionProvider_Dynamic(int *_xresH, int *_yresH, double **_xposH, double **_yposH, int _dim) {
		xresH=_xresH;
		yresH=_yresH;
		xposH=_xposH;
		yposH=_yposH;
		dim=_dim;
	}
	virtual int generate(int layerId, TCostFunctionProviderBase **result) {
		*result=new TCostFunctionProvider_Dynamic(xresH[layerId], yresH[layerId], xposH[layerId], yposH[layerId], dim);
		return 0;
	}
};


/////////////////////////////
// ShieldGenerator
/////////////////////////////

// Base
class TFactoryShieldGeneratorBase {
public:
	virtual ~TFactoryShieldGeneratorBase() {};
	virtual int generate(__attribute__((unused)) int layerId, TShieldGeneratorBase **result) {
		*result=NULL;
		return ERR_BASE_NOTIMPLEMENTED;
	}
};

// Squared Euclidean Distance Dynamic
class TFactoryShieldGeneratorGrid_SqrEuclidean : public TFactoryShieldGeneratorBase {
public:
	int *xDimH, *yDimH;
	int dim;
	TFactoryShieldGeneratorGrid_SqrEuclidean(int *_xDimH, int *_yDimH, int _dim) {
		xDimH=_xDimH;
		yDimH=_yDimH;
		dim=_dim;
	}
	virtual int generate(int layerId, TShieldGeneratorBase **result) {
		*result= new TShieldGeneratorGrid_SqrEuclidean(dim,xDimH+dim*layerId,yDimH+dim*layerId);
		return 0;
	}
};


/////////////////////////////
// CouplingHandlerExt
/////////////////////////////
// For now this is fixed to the type TCouplingHandlerExt<TCouplingHandlerSparse>
// which, as of now, is the only relevant one for all practical purposes
// Base
class TFactoryCouplingHandlerExt {
public:
	int generate(__attribute__((unused)) int layerId,
			int xres, int yres,
			TCostFunctionProviderBase *costFunctionProvider, TVarListHandler *xVars,
			TCouplingHandlerExt<TCouplingHandlerSparse> ** result
			) {
		TCouplingHandlerSparse *couplingHandler;
		couplingHandler = new TCouplingHandlerSparse(xres,yres, costFunctionProvider, xVars);
		TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface;
		couplingHandlerInterface = new TCouplingHandlerExt<TCouplingHandlerSparse>(couplingHandler,true);
		
		*result=couplingHandlerInterface;
		return 0;
	}
};

/////////////////////////////
// SubSolver
/////////////////////////////
//
//


// Base
class TFactorySolverInterfaceBase {
public:
	virtual ~TFactorySolverInterfaceBase() {};
	virtual int generate(
			__attribute__((unused)) int layerId,
			__attribute__((unused)) TCouplingHandlerSparse *couplingHandler,
			__attribute__((unused)) TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface,
			__attribute__((unused)) double *muX, __attribute__((unused)) double *muY,
			__attribute__((unused)) double *alpha, __attribute__((unused)) double *beta,
			TSolverInterface **result) {
		*result=NULL;
		eprintf("not implemented in subsolver factory\n");
		return ERR_BASE_NOTIMPLEMENTED;
	}
	virtual int prepareRefinement(
			__attribute__((unused)) int layerId,
			__attribute__((unused)) TSolverInterface *solverInterface
			) {
		return 0;
	}
	virtual int customizeRefinement(
			__attribute__((unused)) int layerId,
			__attribute__((unused)) TVarListHandler *xVars
			) {
		return 0;
	}
};




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MultiScale Solver
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TMultiScaleSolver {
public:
	// Factories
	TFactoryCostFunctionProviderBase *FactoryCostFunctionProvider;
	TFactoryCouplingHandlerExt *FactoryCouplingHandlerExt;
	TFactorySolverInterfaceBase *FactorySolverInterface;
	TFactoryShieldGeneratorBase *FactoryShieldGenerator;
	// Hierarchical Partitions
	THierarchicalPartition *HPX, *HPY;
	// Hierarchical Mases
	double **muXH,**muYH;
	int VCHECK;
	
	// pointers to ShortCutSolver components
	TCostFunctionProviderBase *costFunctionProvider;
	TCouplingHandlerExt<TCouplingHandlerSparse> *couplingHandlerInterface;
	TSolverInterface *solverInterface;
	TShieldGeneratorBase *shieldGenerator;
	TShortCutSolver *ShortCutSolver;
	// whether to delete pointers after solving
	bool autoDeletePointers;
	
	// after solving, these variables store optimal primal and dual variables
	TVarListHandler *xVarsFinal;
	double *muFinal;
	double *alpha, *beta;
	double objective; // value of optimal objective


	int layerCoarsest; // index of first layer to consider for hierarchical solver

	
	TMultiScaleSolver(
			TFactoryCostFunctionProviderBase *_FactoryCostFunctionProvider,
			TFactoryCouplingHandlerExt *_FactoryCouplingHandlerExt,
			TFactorySolverInterfaceBase *_FactorySolverInterface,
			TFactoryShieldGeneratorBase *_FactoryShieldGenerator,
			// Hierarchical Partitions
			THierarchicalPartition *_HPX, THierarchicalPartition *_HPY,
			// Hierarchical Mases
			double **_muXH, double **_muYH,
			// coarsest level to start with solving
			int _layerCoarsest,
			// Violation check method
			int _VCHECK
			);
			
	~TMultiScaleSolver();
	void cleanupShortCutComponents();
	int solve();

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hierarchical Solver Initialization
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class TMultiScaleSetupBase {
public:
	TDoubleMatrix *posX, *posY; // point clouds for marginal positions
	double *muX, *muY; // pointers to marginal masses
	int xres, yres; // integers for total cardinality of marginal points
	int dim; // dimensionality of marginal positions
	
	// hierarchy
	int depth; // parameter for controlling nr of layers in hierarchical partition
	THierarchyBuilder *HBX,*HBY; // hierarchy builder classes
	THierarchicalPartition *HPX,*HPY; // hierarchical partition classes
	double **posXH, **posYH; // hierarchical positions
	double **muXH,**muYH; // hierarchical masses
	int *xresH, *yresH; // hierarchical marginal cardinality
		
	// factories
	TFactoryCostFunctionProviderBase *FactoryCostFunctionProvider;
	TFactoryShieldGeneratorBase *FactoryShieldGenerator;
	TFactoryCouplingHandlerExt *FactoryCouplingHandlerExt;
	TFactorySolverInterfaceBase *FactorySolverInterface;
	
	//////////////////////////////////////////////////////////////////////////////
	
	TMultiScaleSetupBase(TDoubleMatrix *_posX, TDoubleMatrix *_posY, double *_muX, double *_muY, int _depth);
	virtual ~TMultiScaleSetupBase();

	int BasicSetup();
	int BasicMeasureChecks();
	virtual int SetupSolverSpecificPrepropcessing() { return ERR_BASE_NOTIMPLEMENTED; };
	int SetupHierarchicalPartition(double *mu, double *pos, int res, int dim, int depth,
			THierarchyBuilder **_HB, THierarchicalPartition **_HP, double ***_posH, double ***_muH, int **_resH);
	int SetupHierarchicalPartitions();
	
	
	// methods for initializing factories
	virtual int SetupFactoryCostFunctionProvider() { return ERR_BASE_NOTIMPLEMENTED; };
	virtual int SetupFactoryShieldGenerator() { return ERR_BASE_NOTIMPLEMENTED; };
	virtual int SetupFactoryCouplingHandlerExt();
	virtual int SetupFactorySolverInterface() { return ERR_BASE_NOTIMPLEMENTED; };


	virtual int Setup();
	virtual int SetupNoSolver();
	
	
};


class TMultiScaleSetupW2Grid : public TMultiScaleSetupBase {
public:
	// constructor argument datatype
	// needed to simplify inheritance for subsolver templates
	struct ConstructorArgument {
		TDoubleMatrix *muXGrid, *muYGrid;
		int depth;
	};
	// multidimensional arrays of marginal measures
	TDoubleMatrix *muXGrid, *muYGrid;
	// grid dimensions of each hierarchy level. required for shield generators
	int *xDimH,*yDimH;
	
	TMultiScaleSetupW2Grid(TDoubleMatrix *_muXGrid, TDoubleMatrix *_muYGrid, int _depth);
	// delegate constructor
	TMultiScaleSetupW2Grid(ConstructorArgument arg) :
			TMultiScaleSetupW2Grid(arg.muXGrid,arg.muYGrid,arg.depth) {};
	// destructor
	virtual ~TMultiScaleSetupW2Grid();
	
	virtual int SetupFactoryCostFunctionProvider();
	virtual int SetupFactoryShieldGenerator();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr double MultiScaleRefineBasis_MassTolerance=1E-12;


int MultiScaleRefineBasis(THierarchicalPartition *HPX, THierarchicalPartition *HPY,
		TVarListHandler *xVarsC,
		bool *basisC, double *piC,
		double *muXF, double *muYF, TVarListHandler *xVarsF,
		int layerC,
		bool **basisFRes, double **piFRes
		);
int MultiScaleRefineBasis_NWCinCell(
		int *xList, int *yList, double *muXF, double *muYF, double *muXFSpent, double *muYFSpent,
		int xresLoc, int yresLoc,
		int *xActive, int *yActive,
		TVarListSignal<bool> *basisFSignal, TVarListSignal<double> *piFSignal,
		double m
		);

#endif
