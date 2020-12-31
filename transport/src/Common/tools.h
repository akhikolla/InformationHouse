#ifndef __tools_H
#define __tools_H

#include<Common/PythonTypes.h>
#include<cstdlib>

void doubleArrayCopy(double *a, double *b, int n);
void doubleArrayScale(double *a, double b, int n);

double doubleArrayMin(double *data, int res);
double doubleArrayMax(double *data, int res);

void freeTDoubleMatrix(TDoubleMatrix *a);


#endif
