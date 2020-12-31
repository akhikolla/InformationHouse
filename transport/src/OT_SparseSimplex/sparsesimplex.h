#ifndef SPARSESIMPLEX_H_
#define SPARSESIMPLEX_H_

// #include"OT_SparseSimplex/sparsebasicfeasible.h"
// apprently #includes belong into header file only if header needs them directly (otherwise this increases compilation)

int sparsesimplex(int mm, int nn, int *a, int *b, double *costm, int *channels_byrow_over, int **channels_byrow,
		  int *assignment, int *basis, double *u, double *v, int startgiven, int dense);
#endif
