#include"tools.h"

void doubleArrayCopy(double *a, double *b, int n) {
	int i;
	for(i=0;i<n;i++) {
		b[i]=a[i];
	}
}

void doubleArrayScale(double *a, double b, int n) {
	int i;
	for(i=0;i<n;i++) {
		a[i]=a[i]*b;
	 }
}


double doubleArrayMin(double *data, int res) {
	// compute minimum of array
	double result=data[0];
	for(int i=1;i<res;i++) {
		if(result>data[i]) {
			result=data[i];
		}
	}
	return result;
}

double doubleArrayMax(double *data, int res) {
	// compute maximum of array
	double result=data[0];
	for(int i=1;i<res;i++) {
		if(result<data[i]) {
			result=data[i];
		}
	}
	return result;
}


void freeTDoubleMatrix(TDoubleMatrix *a) {
	free(a->dimensions);
	free(a->data);
	free(a);
}

