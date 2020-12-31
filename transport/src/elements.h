/*

   elements.h
   (mainly defines basic operations on arrays)

   $Revision: 0.1 $   $Date: 2013/09/08 14:48:00 $

   Code by Dominic Schuhmacher
   
*/

#ifndef ELEMENTS_H
#define ELEMENTS_H

/* in auction */
double arraysum(double *a, int n);
int arrayargmax(double *a, int n);
double arraymax(double *a, int n);
double arraysec(double *a, int n, int arg);
double arraymin(double *a, int n);

/* in primaldual */
int arrayisum(int *a, int n);
int arrayimin(int *a, int n);

#endif

