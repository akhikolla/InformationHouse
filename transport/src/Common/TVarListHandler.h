#ifndef TVARLISTHANDLER_H_
#define TVARLISTHANDLER_H_

#include<cstdlib>
#include<vector>

#include<Common/ErrorCodes.h>
#include<Common/Verbose.h>
#include<Common/MergeSort.h>

using namespace std;

class TVarListHandler {
public:
	int res,total;
	vector<int> *lenList;
	vector<int> **varList;
	TVarListHandler();
	TVarListHandler(TVarListHandler *base);
	TVarListHandler(int *indices, int *indptr, int _res, int _total, bool sorted=false);
	virtual ~TVarListHandler();
	void clear();
	void setupEmpty(int _res);
	void fillViaTranspose(TVarListHandler *transpose, int yres);
	void fillFromCSRIndexList(int *indices, int *indptr, int _res, int _total);
	void writeToCSRIndexList(int *indices, int *indptr);
	int merge(TVarListHandler *addition);
	int mergeSelected(TVarListHandler *addition, bool *selection);
	void addToLine(int x, vector<int>* yCandidates);
	void addToLine(int x, int yCandidate);
	void addToLine(int x, int yCandidate, bool testDuplicate);
	void sort();
	void sort(int x);
	
	static bool LowerEq(int a, int b);

	// iteration tools
	struct TIterator {
		int x,yIndex,y,offset;
		bool iterationInitialized;
	};
	TIterator iterationInitialize();
	bool iterate(TIterator *it);

};


template <class T>
class TVarListSignal {
public:
	TVarListHandler *varList;
	T *signal;
	bool internalSignal;
	
	int *offsets;
	bool computedOffsets;
	// pointer offsets for each row of varList in internalSignal, to simplify read/write access
	// analogous to TCouplingHandler
	
	TVarListSignal(TVarListHandler *_varList, T *_signal);
	TVarListSignal(TVarListHandler *_varList, T init);
	~TVarListSignal();
	
	void transcribeSorted(TVarListSignal<T> *src, T defaultValue);
	
	void computeOffsets(); // compute offsets for pointers
	void write(int x, int y, T value); // write to signal by using actual indices
	void writeFromTranspose(TVarListSignal<T> *transpose); // write to signal by copying transposed signal
};


TVarListHandler* GetFullVarList(int xres, int yres);
bool VarListTools_HasEmptyRows(TVarListHandler *vars);

#endif
