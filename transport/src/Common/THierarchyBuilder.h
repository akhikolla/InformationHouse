#ifndef THierarchyBuilder_H_
#define THierarchyBuilder_H_

#include<cstdlib>
#include<cmath>
#include<vector>
#include<Common/THierarchicalPartition.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class THierarchyBuilderNode {
public:
	vector<int> children;
	vector<int> leaves;
	vector<int> posCode;
	int parent;
};

class THierarchyBuilderLayer {
public:
	vector<THierarchyBuilderNode> nodes;
};	


class THierarchyBuilder {
public:
	static constexpr int CM_Tree=0;
	static constexpr int CM_Grid=1;
	static constexpr double boxTolerance=1.E-10;
	
	double *points;
	int nPoints; // number of points at lowest level
	int dim;
	vector<double> boxLo, boxHi; // lower and higher bounds for box around point cloud;
	vector<THierarchyBuilderLayer> layers;
	int childMode;
	
	THierarchyBuilder(double *_points, int _nPoints, int _dim,
			int _childMode, int partitionDepth);
	
	void setBox();
	void reset();
	void refine();
	vector<vector<int> > getChildrenPosCodes(int layerId, int nodeId);
	vector<vector<int> > getChildrenLeaves(int layerId, int nodeId, vector<vector<int> > childrenPosCodes);
	void getRelPosCodeFromIndex(int index, int dim, int *posCode);
	void getOffsetPosCode(int *relPosCode, int *parentPosCode, int dim);
	bool isInBox(double *coord, int *posCode, int dim, int layerId);
	void addAtomicLayer();
	
	static double max(double *x, int n, int step, int offset);
	static double min(double *x, int n, int step, int offset);
	
	THierarchicalPartition* convert();
	double** allocateDoubleSignal(int sigdim, int lBottom=0);
	void freeSignal(double **signal, int lBottom);
	void getSignalPos(double **signal);
	int* getDimH(int *finestDims);
	int* getResH();
};


#endif
