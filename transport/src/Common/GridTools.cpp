#include"GridTools.h"


void GridToolsGetStrides(int dim, int *dims, int *strides) {
	int d;
	strides[dim-1]=1;
	for(d=dim-2;d>=0;d--) {
		strides[d]=strides[d+1]*dims[d+1];
	}
}

int GridToolsGetIdFromPos(int dim, int *pos, int *strides) {
	int d, result;
	result=0;
	for(d=0;d<dim;d++) {
		result+=pos[d]*strides[d];
	}
	return result;
}

void GridToolsGetPosFromId(int dim, int id, int *pos, int *strides) {
	int d;
	pos[0]=id/strides[0];
	for(d=1;d<dim;d++) {
		pos[d]=(id%strides[d-1])/strides[d];
	}
}

//////////////////////////////

void GridToolsGetNeighbours(int dim, int *dims, TVarListHandler *neighbours) {
		GridToolsGetNeighbours_Torus(dim,dims,0,neighbours);
}

void GridToolsGetNeighbours_Torus(int dim, int *dims, int torusDim, TVarListHandler *neighbours) {
	/* dim: dimension of grid, dims: points per axis, torusDim: dimensions 0 to torusDim-1 are considered cyclical. */
	int *strides;
	strides=(int*) malloc(sizeof(int)*dim);
	int *pos;
	pos=(int*) malloc(sizeof(int)*dim);

	GridToolsGetStrides(dim,dims,strides);

	GridToolsGetNeighbours_Torus_iterateXVariables(neighbours, pos, dims, strides, dim, torusDim, 0);

	free(strides);
	free(pos);
}

void GridToolsGetNeighbours_Torus_iterateXVariables(TVarListHandler *neighbours, int *pos, int *dims, int *strides,
		int dim, int torusDim, int d) {
	int x;
	if(d<dim) {
		// as long as we have not reached the last axis, launch iteration along next axis
		for(x=0;x<dims[d];x++) {
			pos[d]=x;
			GridToolsGetNeighbours_Torus_iterateXVariables(neighbours,pos,dims,strides,dim,torusDim,d+1);
		}
	} else {
		// once a specific grid point is reached, add points
		//addVariables_Shields(xVars,xMap,xPos);
		int xId,d;
		xId=GridToolsGetIdFromPos(dim,pos,strides);
		for(d=0;d<dim;d++) {
			if(pos[d]>0) {
				neighbours->addToLine(xId,xId-strides[d]);
			} else {
				if(d<torusDim) {
					neighbours->addToLine(xId,xId+strides[d]*(dims[d]-1));
				}
			}
			if(pos[d]<dims[d]-1) {
				neighbours->addToLine(xId,xId+strides[d]);
			} else {
				if(d<torusDim) {
					neighbours->addToLine(xId,xId-strides[d]*(dims[d]-1));
				}
			}
		}

	}
}




int GridToolsGetTotalPoints(int depth, int *dimensions) {
	int totalSize=1;
	for(int d=0;d<depth;d++) {
		totalSize*=dimensions[d];
	}
	return totalSize;
}

double* GridToolsGetGrid(int depth, int *dimensions) {
	double *result;
	// total number of grid points
	int totalSize=GridToolsGetTotalPoints(depth,dimensions);
	
	// allocate grid points
	result=(double*) malloc(sizeof(double)*totalSize*depth);
	
	// set grid points
	// iterate along each dimension. set coordinate of that dimension for all matching points
	for(int d=0;d<depth;d++) {
		// set combined dimension cardinality variables
		// for all coarser dimensions
		int cardCoarse=GridToolsGetTotalPoints(d,dimensions);

		// for finer dimensions
		int cardFine=GridToolsGetTotalPoints(depth-d-1,dimensions+d+1);
		
		for(int i=0;i<cardCoarse;i++) {
			for(int j=0;j<dimensions[d];j++) {
				for(int k=0;k<cardFine;k++) {
					result[i*dimensions[d]*cardFine*depth+j*cardFine*depth+k*depth+d]=j;
				}
			}
		}
	}
	return result;
}

TDoubleMatrix* GridToolsGetGridMatrix(int depth, int *dimensions) {
	TDoubleMatrix *result=(TDoubleMatrix*) malloc(sizeof(TDoubleMatrix));
	result->data=GridToolsGetGrid(depth,dimensions);
	result->depth=2;
	result->dimensions=(int*) malloc(sizeof(int)*2);
	result->dimensions[0]=GridToolsGetTotalPoints(depth,dimensions);
	result->dimensions[1]=depth;
	return result;
	
}

////////////////////////////////////////



double MeasureToolsTruncateMeasure(double *muX, int xres, double measureScale) {
	// truncate a single measure to integer values, , by dividing through measureScale and rounding.
	// return total truncated mass
	double result=0;
	for(int i=0;i<xres;i++) {
		muX[i]=round(muX[i]/measureScale);
		result+=muX[i];
	}
	return result;
}

int MeasureToolsTruncateMeasures(double *muX, double *muY, int xres, int yres, double measureScale) {
	// truncate two measures to integer values, by dividing through measureScale and rounding.
	// small differences in total mass are then compensated by adjusting the last entry
	double muXsum,muYsum;
	muXsum=MeasureToolsTruncateMeasure(muX,xres,measureScale);
	muYsum=MeasureToolsTruncateMeasure(muY,yres,measureScale);
	
//	// try to compensate mass difference by adjusting last entries
//	if(muXsum<muYsum) {
//		muX[xres-1]+=muYsum-muXsum;
//	} else {
//		muY[yres-1]+=muXsum-muYsum;
//	}
	// this will only work if error due to truncation is small compared to mass in single bins
	// (in particular the last bin) if this difference is too large
	// the compensation will lead to a negative mass value which results in an error.

	// hopefully better:
	// compensate mass difference by iterating over measuer with fewer mass and increasing entries until balance is restored
	double *muLess;
	int resLess;
	int increments;
	if(muXsum<muYsum) {
		muLess=muX;
		resLess=xres;
		increments=round(muYsum-muXsum);
	} else {
		muLess=muY;
		resLess=yres;
		increments=round(muXsum-muYsum);
	}
	int offset=0;
	while(increments>0) {
		muLess[offset]+=1.;
		increments--;
		offset++;
		if(offset>=resLess) {
			offset=0;
		}
	}

	// if somee entries have been rounded to zero, throw error
	if(doubleArrayMin(muX,xres)<=0.) {
		return ERR_PREP_TRUNC_MUXNEG;
	}
	if(doubleArrayMin(muY,yres)<=0.) {
		return ERR_PREP_TRUNC_MUYNEG;
	}

	// return 0 to indicate successful truncation
	return 0;	
}
