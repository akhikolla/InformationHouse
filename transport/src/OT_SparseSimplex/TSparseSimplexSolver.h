#ifndef TSparseSimplexSolver_H_
#define TSparseSimplexSolver_H_

#include<cstdlib>
#include<iostream>
#include<time.h>
extern "C" {
#include<OT_SparseSimplex/sparsesimplex.h>
}
#include<Common/TCouplingHandler.h>
#include<Common/tools.h>


using namespace std;

class TSparseSimplexSolverBase {
public:
	static const int MSG_NOT_IMPLEMENTED = -207;

	bool setupStatus; // true if setup was run
	bool solutionStatus; // true if problem has been run successfully.

	int *muX, *muY;
	double objective;
	int xres, yres;
	int basisstartgiven;
	double *alpha, *beta; // dual variables

	// arrays to store assignment and basis. for now dense. to be replaced by sparse data structure.
	int *assignment;
	int *basis;
	
	// variable to keep tack whether class instance is responsible for freeimg memory of muX and muY
	bool deleteMarginals;

	TSparseSimplexSolverBase(int *_muX, int *_muY, double *_alpha, double *_beta, bool _deleteMarginals);

	virtual ~TSparseSimplexSolverBase();

	virtual int solve();
	
	virtual int setupBaseZero();
	virtual int setup();
	virtual int cleanup();

	virtual double getObjective();


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class TCouplingHandlerType>
class TSparseSimplexSolver : public TSparseSimplexSolverBase {
public:

	TCouplingHandlerType *CouplingHandler;

	TSparseSimplexSolver(TCouplingHandlerType *_CouplingHandler, int *_muX, int *_muY, double *_alpha, double *_beta, bool _deleteMarginals);

	int solve();


};






template <class TCouplingHandlerType>
TSparseSimplexSolver<TCouplingHandlerType>::TSparseSimplexSolver(TCouplingHandlerType *_CouplingHandler,
		int *_muX, int *_muY, double *_alpha, double *_beta,  bool _deleteMarginals) :
		TSparseSimplexSolverBase(_muX, _muY, _alpha, _beta, _deleteMarginals) {
	CouplingHandler=_CouplingHandler;
	xres=CouplingHandler->xres;
	yres=CouplingHandler->yres;
}


template <class TCouplingHandlerType>
int TSparseSimplexSolver<TCouplingHandlerType>::solve() {

	
	int x,y,rowLen,yIndex;
	// allocate dense cost matrix and fill with values from coupling handler.
	// COLUMN MAJOR!
	double *costm=(double*) malloc(xres*yres*sizeof(double));
		
	
	// the standard interface for CouplingHandlers is: access elements row wise.
	// iterate over all rows
	for(x=0;x<xres;x++) {
		// extract number of entries in this row
		rowLen=CouplingHandler->getRowLength(x);
		// iterate over number of entries
		for(yIndex=0;yIndex<rowLen;yIndex++) {
			// extract position of entry
			y=CouplingHandler->getColNr(x,yIndex);
			// write cost function value to cost matrix in colum major format
			costm[y*xres+x]=CouplingHandler->getCRow(x,yIndex);
		}
	}
	
	// create compressed sparse row-like representation of neighbourhood
	int *channels_byrow_over; // list of entries per row
	int **channels_byrow;  // list of lists of entries in each row

	channels_byrow_over=(int*) malloc(xres*sizeof(int));
	channels_byrow=(int**) malloc(xres*sizeof(int*));
	for(x=0;x<xres;x++) {
		rowLen=CouplingHandler->getRowLength(x);
		channels_byrow_over[x]=rowLen;
		channels_byrow[x]=(int*) malloc(rowLen*sizeof(int));
		for(yIndex=0;yIndex<rowLen;yIndex++) {
			channels_byrow[x][yIndex]=CouplingHandler->getColNr(x,yIndex);
		}
	}
	
	
	// START CHANGE DOMINIC
	int dense=0;    // D: If we set dense=1 the dense starting solution is computed
	                //    If this is a feasible solution everything works as before
	int time1,time2;
	time1=clock();
	eprintf("\t\tcalling sparsesimplex. startgiven: %d, nr of vars: %d\n",basisstartgiven,CouplingHandler->total);
	sparsesimplex(xres, yres, muX, muY, costm, channels_byrow_over, channels_byrow, assignment, basis, alpha, beta, basisstartgiven, dense);
	time2=clock();
	eprintf("\t\ttime: %d\n",(time2-time1));
	
	
//	// print full basis
//	cout << "basis:" << endl;
//	for(y=0;y<yres;y++) {
//		cout << "{ " << basis[y*xres];
//		for(x=1;x<xres;x++) {
//			cout << ", " << basis[y*xres+x];
//		}
//		cout << "}," << endl;
//	}


//	// print full assignment
//	cout << "assignment:" << endl;
//	for(y=0;y<yres;y++) {
//		cout << "{ " << assignment[y*xres];
//		for(x=1;x<xres;x++) {
//			cout << ", " << assignment[y*xres+x];
//		}
//		cout << "}," << endl;
//	}

//	// print marginals
//	cout << "muX:" << endl;
//	cout << "{ " << muX[0];
//	for(x=1;x<xres;x++) {
//		cout << ", " << muX[x];
//	}
//	cout << "}" << endl;
//	cout << "muY:" << endl;
//	cout << "{ " << muY[0];
//	for(y=1;y<yres;y++) {
//		cout << ", " << muY[y];
//	}
//	cout << "}" << endl;


	// just pretend everything goes well for now.
	solutionStatus=true;

	// copy optimal assignment back to coupling handler
	// compute (primal) objective in same loop
	objective=0;
	for(x=0;x<xres;x++) {
		rowLen=CouplingHandler->getRowLength(x);
		for(yIndex=0;yIndex<rowLen;yIndex++) {
			y=CouplingHandler->getColNr(x,yIndex);
			CouplingHandler->setMuRow(x,yIndex,assignment[y*xres+x]);
			objective+=CouplingHandler->getCRow(x,yIndex)*CouplingHandler->getMuRow(x,yIndex);
		}
	}
	

	free(costm);
	for(x=0;x<xres;x++) {
		free(channels_byrow[x]);
	}
	free(channels_byrow);
	free(channels_byrow_over);
		
	return 0;

}


#endif
