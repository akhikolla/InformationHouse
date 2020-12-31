#include"MultiScaleSolver.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MultiScale Solver
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TMultiScaleSolver::TMultiScaleSolver(
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
		) {
	FactoryCostFunctionProvider=_FactoryCostFunctionProvider;
	FactoryCouplingHandlerExt=_FactoryCouplingHandlerExt;
	FactorySolverInterface=_FactorySolverInterface;
	FactoryShieldGenerator=_FactoryShieldGenerator;
	HPX=_HPX;
	HPY=_HPY;
	muXH=_muXH;
	muYH=_muYH;
	
	layerCoarsest=_layerCoarsest;
	
	VCHECK=_VCHECK;
	
	costFunctionProvider=NULL;
	couplingHandlerInterface=NULL;
	solverInterface=NULL;
	shieldGenerator=NULL;
	ShortCutSolver=NULL;

	
	
	xVarsFinal=NULL;
	muFinal=NULL;
	alpha=NULL;
	beta=NULL;
	
	objective=0.;
	
	// whether to auto delete pointers to shortcut components upon completing solving
	autoDeletePointers=true;
	// for now simply initialize this to true.
	// can be overwritten manually after construction, if desired

}

TMultiScaleSolver::~TMultiScaleSolver() {

	cleanupShortCutComponents();

	if(xVarsFinal!=NULL) {
		delete xVarsFinal;
	}
	if(muFinal!=NULL) {
		free(muFinal);
	}
	if(alpha!=NULL) {
		free(alpha);
	}
	if(beta!=NULL) {
		free(beta);
	}
}

void TMultiScaleSolver::cleanupShortCutComponents() {
	if(costFunctionProvider!=NULL) {
		delete costFunctionProvider;
		costFunctionProvider=NULL;
	}
	if(couplingHandlerInterface!=NULL) {
		delete couplingHandlerInterface;
		couplingHandlerInterface=NULL;;
	}
	if(solverInterface!=NULL) {
		delete solverInterface;
		solverInterface=NULL;;
	}
	if(shieldGenerator!=NULL) {
		delete shieldGenerator;
		shieldGenerator=NULL;;
	}
	if(ShortCutSolver!=NULL) {
		delete ShortCutSolver;
		ShortCutSolver=NULL;
	}
}

int TMultiScaleSolver::solve() {
	int msg;
	
	// pointer for varLists
	// store which variables are used in current sparse problem (i.e. the current neighbourhood in the paper)
	TVarListHandler *xVars=NULL,*xVarsCoarse=NULL;

	// indexing which layers will be solved.
	// layer 0 is coarsest layer with only one node per hierarchical partition
	// the last layer (with highest possibl index) corresponds to original problem
	int layerFinest=HPX->nLayers-1; // index of finest layer for hierarchical solver
	// (solving layerCoarsest=0 is pointless since it is just one node per hierarchical dimension)

	// index of layer that is currently solved
	int layerId=layerCoarsest;
		
	// while not at finest layer
	while(layerId<=layerFinest) {
		// start solving
		eprintf("current layer: %d\n",layerId);
	
		// cardinality of marginals at current layer
		int xres=HPX->layers[layerId]->nCells;
		int yres=HPY->layers[layerId]->nCells;

		//////////////////////////////////////////////////////////////////////////////////////////////
		// choose initial varList for ShortCutSolver
		if(layerId==layerCoarsest) {
			eprintf("\tfull var list\n");
			// initialize full var list when on coarsest level, i.e. solve dense problem
			xVars=GetFullVarList(xres,yres);
		} else {
			eprintf("\trefining var list\n");
			// refine var list
			// if coarser layer has been solved, generate new initial var list by refining
			xVars=refineVarList(HPX, HPY, xVarsCoarse, layerId-1, true);
			delete xVarsCoarse;
			eprintf("\ttotal new variables: %d\n",xVars->total);
			// allow sub solver to customize refinement
			msg=FactorySolverInterface->customizeRefinement(layerId,xVars);
			if(msg!=0) {
				return msg;
			}
		}

		//////////////////////////////////////////////////////////////////////////////////////////////
		// cost function provider
		// this class models the cost function that is used in the problem.
		msg=FactoryCostFunctionProvider->generate(layerId,&costFunctionProvider);
		if(msg!=0) {
			return msg;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////
		// coupling handler & interface class
		// this is the data structure that stores the current coupling
		// (and the cost function of the currently selected variables)
		// gets a reference to the costFunctionProvider to obtain required cost function values
		// couplingHandler is the actual data structure handler
		// couplingHandlerInterface is a slightly extended interface class with some more methods
		// that will be used by the ShortCutSolver to communicate with the coupling handler
		eprintf("\tcoupling handler interface\n");
		msg=FactoryCouplingHandlerExt->generate(layerId,xres,yres, costFunctionProvider, xVars,&couplingHandlerInterface);
		if(msg!=0) {
			return msg;
		}
		TCouplingHandlerSparse *couplingHandler;
		couplingHandler=couplingHandlerInterface->couplingHandler;


		//////////////////////////////////////////////////////////////////////////////////////////////
		// sub solver
		eprintf("\tsubsolver\n");
		// allocate dual variables
		alpha=(double*) malloc(sizeof(double)*xres);
		beta=(double*) malloc(sizeof(double)*yres);
		// generate instance subsolver & interface class
		msg=FactorySolverInterface->generate(
				layerId,
				couplingHandler,couplingHandlerInterface,muXH[layerId],muYH[layerId],alpha,beta,
				&solverInterface);
		if(msg!=0) {
			return msg;
		}
		
		
		
		//////////////////////////////////////////////////////////////////////////////////////////////
		// shielding
		// class for constructing shielding neighbourhoods
		eprintf("\tshielding generator\n");	
		msg=FactoryShieldGenerator->generate(layerId,&shieldGenerator);
		if(msg!=0) {
			return msg;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////
		// ShortCutSolver
		// gets references to couplingHanderInterface, solverInterface and shieldGenerator
		// implements algorithm 4.1 of paper
		eprintf("\tShortCut solver\n");
		ShortCutSolver= new TShortCutSolver(couplingHandlerInterface,solverInterface,shieldGenerator,
				VCHECK
				);

		msg=ShortCutSolver->initialize();
		if (msg!=0) {
			return msg;
		}


		eprintf("\tsolving\n");
		// maximal number of iterations in algorithm 4.1.
		// 100 should be way more than enough. in practice it takes about 3-4 per level
		msg=ShortCutSolver->step(100);
		if(msg!=0) {
			return msg;
		}

		// quick status update in terminal
		eprintf("\tstatus:\n");
		eprintf("\t\tsolved: %d\n",ShortCutSolver->report.solved);
		eprintf("\t\tsteps: %d\n",ShortCutSolver->report.steps);
		eprintf("\t\tobjective: %f\n",ShortCutSolver->report.objective);

		// something went wrong. problem at current level is not solved
		if(ShortCutSolver->report.solved!=1) {
			return ERR_MULTISCALE_SHORTCUTUNSOLVED;
		}

		
		if(layerId==layerFinest) {
			// if at finest layer, prepare export of results
			xVarsFinal=couplingHandler->xVars;
			muFinal=(double*) malloc(sizeof(double)*xVarsFinal->total);
			doubleArrayCopy(couplingHandler->mu,muFinal,xVarsFinal->total);
			objective=ShortCutSolver->report.objective;
		} else {
			// otherwise prepare refinement and clean up
			
			// prepare refinement for sub solver object
			FactorySolverInterface->prepareRefinement(layerId,solverInterface);
			
			// extract support of optimal solution, will be refined at next level
			xVarsCoarse=couplingHandlerInterface->getSupport();

			// then do some cleaning
			free(alpha);
			free(beta);		
			delete couplingHandler->xVars;
		}
	
		if((layerId<layerFinest)||autoDeletePointers) {
			// clean up stuff at this level
			cleanupShortCutComponents();
		}
			
		layerId++;
	}

	return 0;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hierarchical Solver Initialization
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TMultiScaleSetupBase::TMultiScaleSetupBase(TDoubleMatrix *_posX, TDoubleMatrix *_posY, double *_muX, double *_muY,
		int _depth) {
	posX=_posX;
	posY=_posY;
	muX=_muX;
	muY=_muY;
	depth=_depth;
	
	FactoryCostFunctionProvider=NULL;
	FactoryShieldGenerator=NULL;
	FactoryCouplingHandlerExt=NULL;
	FactorySolverInterface=NULL;
	
	HBX=NULL;
	HBY=NULL;

}

TMultiScaleSetupBase::~TMultiScaleSetupBase() {
	if(HBX!=NULL) {

		free(xresH);		
		HBX->freeSignal(muXH,HBX->layers.size());
		HBX->freeSignal(posXH,HBX->layers.size());
		delete HPX;
		delete HBX;
	}
	if(HBY!=NULL) {

		free(yresH);		
		HBY->freeSignal(muYH,HBY->layers.size());
		HBY->freeSignal(posYH,HBY->layers.size());
		delete HPY;
		delete HBY;
	
		
	}
	
	if(FactoryCouplingHandlerExt!=NULL) {
		delete FactoryCouplingHandlerExt;
	}
	
	
}

int TMultiScaleSetupBase::BasicSetup() {
	if((posX->depth!=2) || (posY->depth!=2)) {
		eprintf("ERROR: marginal point clouds must be 2d arrays.\n");
		return ERR_PREP_INIT_POINTCLOUD2D;
	}
	if(posX->dimensions[1]!=posY->dimensions[1]) {
		return ERR_PREP_INIT_DIMMISMATCH;
	}
	xres=posX->dimensions[0];
	yres=posY->dimensions[0];
	dim=posX->dimensions[1];
	return 0;
}

int TMultiScaleSetupBase::BasicMeasureChecks() {

	// sanity check: muX and muY must be strictly positive
	if(doubleArrayMin(muX,xres)<=0.) {
		eprintf("ERROR: minimum of muX is not strictly positive.\n");
		return ERR_PREP_INIT_MUXNEG;
	}
	if(doubleArrayMin(muY,yres)<=0.) {
		eprintf("ERROR: minimum of muY is not strictly positive.\n");
		return ERR_PREP_INIT_MUYNEG;
	}

	// TODO: more sanity checks?

	return 0;
}


int TMultiScaleSetupBase::SetupHierarchicalPartition(double *mu, double *pos, int res, int dim, int depth,
		THierarchyBuilder **_HB, THierarchicalPartition **_HP, double ***_posH, double ***_muH, int **_resH) {
	
	// create hierarchical partition
	THierarchyBuilder *HB = new THierarchyBuilder(pos,res,dim, THierarchyBuilder::CM_Grid, depth);
		
	// convert to format used by hierarchical solver
	THierarchicalPartition *HP=HB->convert();
	
	// hierarchical positions
	// position of nodes in coarser hierarchy levels
	double **posH=HB->allocateDoubleSignal(dim);
	HB->getSignalPos(posH);
	
	// hierarchical masses
	// combined masses of nodes in coarser hierarchy levels
	double **muH=HB->allocateDoubleSignal(1);
	HP->computeHierarchicalMasses(mu,muH);

	// get hierarchical cardinalities of each marginal on each level
	int *resH=HB->getResH();

	*_HB=HB;
	*_HP=HP;
	*_posH=posH;
	*_muH=muH;
	*_resH=resH;
	
	return 0;
	
}

int TMultiScaleSetupBase::SetupHierarchicalPartitions() {
	int msg;
	msg=SetupHierarchicalPartition(muX, posX->data, xres, dim, depth,
			&HBX, &HPX, &posXH, &muXH, &xresH);
	if(msg!=0) {
		return msg;
	}
	
	msg=SetupHierarchicalPartition(muY, posY->data, yres, dim, depth,
			&HBY, &HPY, &posYH, &muYH, &yresH);
	if(msg!=0) {
		return msg;
	}
		
	return 0;
}

int TMultiScaleSetupBase::SetupFactoryCouplingHandlerExt() {
	FactoryCouplingHandlerExt=new TFactoryCouplingHandlerExt;
	return 0;
}


int TMultiScaleSetupBase::Setup() {
	int msg;

	msg=BasicSetup();
	if(msg!=0) { return msg; }
	msg=BasicMeasureChecks();
	if(msg!=0) { return msg; }
	msg=SetupSolverSpecificPrepropcessing();
	if(msg!=0) { return msg; }
	msg=SetupHierarchicalPartitions();
	if(msg!=0) { return msg; }
	
	msg=SetupFactoryCostFunctionProvider();
	if(msg!=0) { return msg; }
	msg=SetupFactoryShieldGenerator();
	if(msg!=0) { return msg; }
	msg=SetupFactoryCouplingHandlerExt();
	if(msg!=0) { return msg; }
	msg=SetupFactorySolverInterface();
	if(msg!=0) { return msg; }

	return 0;

}

int TMultiScaleSetupBase::SetupNoSolver() {
	int msg;

	msg=BasicSetup();
	if(msg!=0) { return msg; }
	msg=BasicMeasureChecks();
	if(msg!=0) { return msg; }
	msg=SetupHierarchicalPartitions();
	if(msg!=0) { return msg; }
	
	msg=SetupFactoryCostFunctionProvider();
	if(msg!=0) { return msg; }
	msg=SetupFactoryCouplingHandlerExt();
	if(msg!=0) { return msg; }

	return 0;

}


//////////////////////////////////////////

// W2 Grid

TMultiScaleSetupW2Grid::TMultiScaleSetupW2Grid(TDoubleMatrix *_muXGrid, TDoubleMatrix *_muYGrid, int _depth) :
	TMultiScaleSetupBase(NULL, NULL, _muXGrid->data, _muYGrid->data, _depth) {
	muXGrid=_muXGrid;
	muYGrid=_muYGrid;
	posX=GridToolsGetGridMatrix(muXGrid->depth,muXGrid->dimensions);
	posY=GridToolsGetGridMatrix(muYGrid->depth,muYGrid->dimensions);
	
}

TMultiScaleSetupW2Grid::~TMultiScaleSetupW2Grid() {
	freeTDoubleMatrix(posX);
	freeTDoubleMatrix(posY);
	if(FactoryCostFunctionProvider!=NULL) {
		delete FactoryCostFunctionProvider;
	}
	if(FactoryShieldGenerator!=NULL) {
		free(xDimH);
		free(yDimH);
		delete FactoryShieldGenerator;
	}
}

int TMultiScaleSetupW2Grid::SetupFactoryCostFunctionProvider() {
	FactoryCostFunctionProvider=new TFactoryCostFunctionProvider_Dynamic(xresH,yresH,posXH,posYH,dim);
	return 0;
}

int TMultiScaleSetupW2Grid::SetupFactoryShieldGenerator() {
	xDimH=HBX->getDimH(muXGrid->dimensions);
	yDimH=HBY->getDimH(muYGrid->dimensions);
	FactoryShieldGenerator = new TFactoryShieldGeneratorGrid_SqrEuclidean(xDimH,yDimH,dim);
	return 0;
}

//////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Basis Refinement
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// from coarse optimal coupling and simplex basis construct initial basis for finer problem
// HPX and HPY are respective hierarchical partitions
// xVarsC is coarse active varList, piC and basisC are respective entries for optimal coupling and basis
// muXF and muYF are marginals on fine level, xVarsF is refined varlist of xVarsC
// layerC indicates coarse level id
// basisFRes and piFRes are targets for refined basis and coupling
int MultiScaleRefineBasis(THierarchicalPartition *HPX, THierarchicalPartition *HPY,
		TVarListHandler *xVarsC, bool *basisC, double *piC,
		double *muXF, double *muYF, TVarListHandler *xVarsF,
		int layerC,
		bool **basisFRes, double **piFRes
		) {
		
	int msg;
	
	int xresF, yresF;
	int xresC, yresC;
	
	xresF=HPX->layers[layerC+1]->nCells;
	yresF=HPY->layers[layerC+1]->nCells;
	xresC=HPX->layers[layerC]->nCells;
	yresC=HPY->layers[layerC]->nCells;
	
//	// test if coarse input basis has correct number of non-zero entries
//	int basisCSum=0;
//	for(int i=0;i<xVarsC->total;i++) {
//		basisCSum+=basisC[i];
//	}
//	if(basisCSum!=xresC+yresC-1) {
//		eprintf("ERROR: coarse input basis has wrong number of entries: sum(basisC)=%d, xresC+yresC-1=%d\n",basisCSum,xresC+yresC-1);
//		return ERR_MULTISCALE_BASISREFINEMENT_COARSEBASISENTRIES;
//	}	

	
	
	int x,y,yIndex,rowLen;

	// variables that keep track how much weight we have already assigned to each vertex	
	double *muXFSpent=(double*) malloc(sizeof(double)*xresF);
	for(x=0;x<xresF;x++) {
		muXFSpent[x]=0.;
	}
	double *muYFSpent=(double*) malloc(sizeof(double)*yresF);
	for(y=0;y<yresF;y++) {
		muYFSpent[y]=0.;
	}
	
	// variable for refined basis and coupling
	bool *basisF=(bool*) malloc(sizeof(bool)*xVarsF->total);
	double *piF=(double*) malloc(sizeof(double)*xVarsF->total);
	for(x=0;x<xVarsF->total;x++) {
		basisF[x]=false;
		piF[x]=0.;
	}
	// set up TVarListSignal<> to allow write access to basisF via absolute indices
	// same for piF
	TVarListSignal<bool> *basisFSignal=new TVarListSignal<bool>(xVarsF,basisF);
	basisFSignal->computeOffsets();
	TVarListSignal<double> *piFSignal=new TVarListSignal<double>(xVarsF,piF);
	piFSignal->computeOffsets();
	
	// lists for active nodes
	int *xActive=(int*) malloc(sizeof(int)*xresC);
	for(x=0;x<xresC;x++) {
		xActive[x]=0;
	}
	int *yActive=(int*) malloc(sizeof(int)*yresC);
	for(y=0;y<yresC;y++) {
		yActive[y]=0;
	}
	
	//
	int inputPointerOffset=0; // offset to current element in basisC and piC
	for(x=0;x<xresC;x++) {
		rowLen=xVarsC->lenList->at(x);
		int xresLoc=HPX->layers[layerC]->nChildren[x];
		for(yIndex=0;yIndex<rowLen;yIndex++) {
			y=xVarsC->varList[x]->at(yIndex);
			int yresLoc=HPY->layers[layerC]->nChildren[y];
			if(basisC[inputPointerOffset]==1) {
				msg=MultiScaleRefineBasis_NWCinCell(
						HPX->layers[layerC]->children[x],
						HPY->layers[layerC]->children[y],
						muXF, muYF, muXFSpent, muYFSpent,
						xresLoc, yresLoc,
						xActive+x, yActive+y,
						basisFSignal, piFSignal,
						piC[inputPointerOffset]
						);
				if(msg!=0) {
					return msg;
				}

			}
			inputPointerOffset++;
		}
		
	}
	
	//////////////////////////////////////////////////////////////////
	// some sanity checks
	
//	// test if all mass was spent
//	for(x=0;x<xresF;x++) {
//		if(abs(muXF[x]-muXFSpent[x])>=MultiScaleRefineBasis_MassTolerance) {
//			eprintf("ERROR: X marginal of constructed basis solution differs from muX");
//			return ERR_MULTISCALE_BASISREFINEMENT_MARGX;
//		}
//	}
//	for(y=0;y<yresF;y++) {
//		if(abs(muYF[y]-muYFSpent[y])>=MultiScaleRefineBasis_MassTolerance) {
//			eprintf("ERROR: Y marginal of constructed basis solution differs from muY");
//			return ERR_MULTISCALE_BASISREFINEMENT_MARGY;
//		}
//	}

//	// test if refined coupling actually has full mass
//	double piFSum=0;
//	double muXSum=0;
//	for(x=0;x<xresF;x++) {
//		muXSum+=muXF[x];
//	}
//	for(int i=0;i<xVarsF->total;i++) {
//		piFSum+=piF[i];
//	}
//	if(muXSum!=piFSum) {
//		eprintf("ERROR: refined coupling has wrong total mass: sum(piF)=%f, sum(muXF)=%f\n",piFSum,muXSum);
//		return ERR_MULTISCALE_BASISREFINEMENT_BASISENTRIES;
//	}	


//	// test if basis has correct number of non-zero entries
//	int basisFSum=0;
//	for(int i=0;i<xVarsF->total;i++) {
//		if(basisF[i]) {
//			basisFSum+=1;
//		}
//	}
//	if(basisFSum!=xresF+yresF-1) {
//		eprintf("ERROR: refined basis has wrong number of entries: sum(basisF)=%d, xresF+yresF-1=%d\n",basisFSum,xresF+yresF-1);
//		
//		// print full basis
//		inputPointerOffset=0; // offset to current element in basisC and piC
//		for(x=0;x<xresF;x++) {
//			rowLen=xVarsF->lenList->at(x);
//			for(yIndex=0;yIndex<rowLen;yIndex++) {
//				y=xVarsF->varList[x]->at(yIndex);
//				if(basisF[inputPointerOffset]) {
//					eprintf("{%d,%d},",x,y);
//				}
//				inputPointerOffset++;
//			}
//		}
//		eprintf("\n");
//		
//		
//		return ERR_MULTISCALE_BASISREFINEMENT_BASISENTRIES;
//	}	
		
	//////////////////////////////////////////////////////////////////

	
	// free all variables
	free(muXFSpent);
	free(muYFSpent);
	free(xActive);
	free(yActive);
	delete basisFSignal;
	delete piFSignal;
	
	*basisFRes=basisF;
	*piFRes=piF;
	
	return 0;
	
}

int MultiScaleRefineBasis_NWCinCell(
		int *xList, int *yList, double *muXF, double *muYF, double *muXFSpent, double *muYFSpent,
		int xresLoc, int yresLoc,
		int *xActive, int *yActive,
		TVarListSignal<bool> *basisFSignal, TVarListSignal<double> *piFSignal,
		double m
		) {
	if(m<MultiScaleRefineBasis_MassTolerance) {
		// if no mass is to be distributed in this cell, only a basis link is to be set
		// this is arbitrarily set at (0,0) in the xList,yList coordinates
		basisFSignal->write(xList[0],yList[0],true);
		return 0;
	} else {
		// else start local north west corner rule
		double q=m;
		int x=*xActive;
		int y=*yActive;
		double dX;
		double dY;
		double delta;
		while((x<xresLoc) && (y<yresLoc) && (q>MultiScaleRefineBasis_MassTolerance)) {
			// determine mass that is left on each active vertex
			dX=muXF[xList[x]]-muXFSpent[xList[x]];
			dY=muYF[yList[y]]-muYFSpent[yList[y]];
			
			// check that this is non-zero
			if(dX<=MultiScaleRefineBasis_MassTolerance) {
				eprintf("ERROR: active x was depleted in basis refinement.\n");
				return ERR_MULTISCALE_BASISREFINEMENT_DEPLETED;
			}
			if(dY<=MultiScaleRefineBasis_MassTolerance) {
				eprintf("ERROR: active y was depleted in basis refinement.\n");
				return ERR_MULTISCALE_BASISREFINEMENT_DEPLETED;
			}
			
			// some mass will be put onto (x,y). therefore, this will be a basis field
			basisFSignal->write(xList[x],yList[y],true);
			// determine mass for this field and update refined coupling
			delta=min(q,min(dX,dY));
			piFSignal->write(xList[x],yList[y],delta);			
			
			
			// check if auxiliary edges are required			
			if(q>delta+MultiScaleRefineBasis_MassTolerance) {
				// delta is induced by dX or dY, there will be another step on this tile
				// only need one aux edge (does not matter which), if both are increased and tile is not completed
				if(abs(dX-dY)<MultiScaleRefineBasis_MassTolerance) {
					if((x<xresLoc-1)&&(y<yresLoc-1)) {
						basisFSignal->write(xList[x],yList[y+1],true);
					}
					// if exactly one of the two is true, the while loop will terminate with q>0 and signal an error
				}
			} else {
				// delta is induced by q. this is the final step on this tile
				// any side that is increased and allows an aux edge, gets one
				if((dX<=delta+MultiScaleRefineBasis_MassTolerance) && (x<xresLoc-1)) {
					basisFSignal->write(xList[x+1],yList[y],true);
				}
				if((dY<=delta+MultiScaleRefineBasis_MassTolerance) && (y<yresLoc-1)) {
					basisFSignal->write(xList[x],yList[y+1],true);
				}
			}


			// update spent masses and q
			muXFSpent[xList[x]]+=delta;
			muYFSpent[yList[y]]+=delta;
			q-=delta;
			// check if active vertices must be increased (if vertex mass is depleted)
			if(dX<=delta+MultiScaleRefineBasis_MassTolerance) {
				x++;
			}
			if(dY<=delta+MultiScaleRefineBasis_MassTolerance) {
				y++;
			}


				
		}
		// check if all budget was spent
		if(q>=MultiScaleRefineBasis_MassTolerance) {
			eprintf("ERROR: not all mass was spent in local north west corner rule during basis refinement.\n");
			return ERR_MULTISCALE_BASISREFINEMENT_EXCESS;
		}
		*xActive=x;
		*yActive=y;
		return 0;
	}	
}
