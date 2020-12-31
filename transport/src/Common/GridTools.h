#ifndef GridTools_H_
#define GridTools_H_

#include<cmath>

#include<Common/TVarListHandler.h>
#include<Common/ErrorCodes.h>
#include<Common/PythonTypes.h>
#include<Common/tools.h>

using namespace std;

void GridToolsGetStrides(int dim, int *dims, int *strides);
int GridToolsGetIdFromPos(int dim, int *pos, int *strides);
void GridToolsGetPosFromId(int dim, int id, int *pos, int *strides);

void GridToolsGetNeighbours(int dim, int *dims, TVarListHandler *neighbours);
void GridToolsGetNeighbours_Torus(int dim, int *dims, int torusDim, TVarListHandler *neighbours);
void GridToolsGetNeighbours_Torus_iterateXVariables(TVarListHandler *neighbours, int *pos, int *dims, int *strides,
		int dim, int torusDim, int d);


int GridToolsGetTotalPoints(int depth, int *dimensions);
double* GridToolsGetGrid(int depth, int *dimensions);
TDoubleMatrix* GridToolsGetGridMatrix(int depth, int *dimensions);


double MeasureToolsTruncateMeasure(double *muX, int xres, double measureScale);
int MeasureToolsTruncateMeasures(double *muX, double *muY, int xres, int yres, double measureScale);


#endif

