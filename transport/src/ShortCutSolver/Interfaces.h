#ifndef Interfaces_H_
#define Interfaces_H_

#include<cstdlib>
#include<Common/TCouplingHandler.h>
#include<Common/TVarListHandler.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Defines prototypes for interface classes for TCouplingHandler and the subsolver for the sparse
// transport problems, which will be used by the TShortCutSolver class.
////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coupling Handler Interface
////////////////////////////////////////////////////////////////////////////////////////////////////////


class TCouplingHandlerExtBase {
	public:
	static constexpr double MASS_TOLERANCE=1E-12;
	static constexpr double SLACK_TOLERANCE=1E-5;
	static constexpr double SLACK_INFINITY=1E12;

	TCouplingHandlerExtBase();
	virtual ~TCouplingHandlerExtBase();

	virtual int getXres();
	virtual int getYres();
	virtual TVarListHandler* getXVars();
	virtual TVarListHandler* getSupport();
	virtual TVarListSignal<double>* getSupportSignal();
	virtual void updateXVars(TVarListHandler *_newXVars, bool keepXVars);
	virtual bool dualViolationCheck(bool doProjection, double *alpha, double *beta);
};



template<class TCouplingHandlerType>
class TCouplingHandlerExt : public TCouplingHandlerExtBase {
	public:
	TCouplingHandlerType *couplingHandler;
	bool deleteHandlerOnDestroy;

	TCouplingHandlerExt();
	TCouplingHandlerExt(TCouplingHandlerType *_couplingHandler, bool _deleteHandlerOnDestroy=false);
	~TCouplingHandlerExt();

	int getXres();
	int getYres();
	TVarListHandler* getXVars();

	TVarListHandler* getSupport();
	TVarListSignal<double>* getSupportSignal();
	void updateXVars(TVarListHandler *_newXVars, bool keepXVars);
	bool dualViolationCheck(bool doProjection, double *alpha, double *beta);
};




////////////////////////////////////////////////////////////////////////////////////////////////////////
// SubSolver Interface
////////////////////////////////////////////////////////////////////////////////////////////////////////

class TSolverInterface {
public:
	static const int MSG_NOT_IMPLEMENTED=-10001;
	double *alpha,*beta;
	TSolverInterface() {
		alpha=NULL;
		beta=NULL;
	}
	virtual ~TSolverInterface() {}

	virtual int solve();
	virtual int prepareUpdate(__attribute__((unused)) TVarListHandler *newXVars);
	virtual int update();
	virtual double getObjective();
};


#endif
