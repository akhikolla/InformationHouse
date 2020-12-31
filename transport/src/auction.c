/*

   auction.c
   (was Wasserstein/Code/auction2.c)

   $Revision: 0.4 $   $Date: 2013/09/08 14:42:00 $

   Code by Dominic Schuhmacher
   
*/


/* n >= 2 is assumed throughout !!!!!!!!! */

/* #include <stdio.h> */
/* #include <stdlib.h> */
#include <R.h>
#include <math.h>
#include <R_ext/Utils.h>

#include "elements.h"

typedef struct State {
  int n; 
  double epsbid;   /* in this new version: the current eps */
  int nofassigned;  /* number of assigned persons */
  int *pers_to_obj;  /* -1 means unassigned */
  int *obj_to_pers;  /* -1 means unassigned */           
  double *price;   
  int *desiremat;        /* matrix of desires */ 
  double *persvalue; /* desire minus price of current person */
                     /* last two only used in bid, but maybe better
                        to reserve memory once and for all */
  /* int *flowmatrix;        matrix of flows */
} State;

#define DESIRE(I,J,STATE,NVALUE) ((STATE)->desiremat)[(NVALUE) * (J) + (I)]
#define MIN(A,B) ((A)<(B) ? (A) : (B))

void bid(State *state, int i);
/* double arraysum(double *a, int n);
int arrayargmax(double *a, int n);
double arraymax(double *a, int n);
double arraysec(double *a, int n, int arg);
double arraymin(double *a, int n); */
/* void printit(State *state); */



/* ------------ The main function ----------------------------- */

void auction(int *desirem, int *nn, int *pers_to_obj, double *price, int *kk, double *eps)
{
  int i,j,r; /* indices */
   int k,n;
   State state;

   /* inputs */
   state.n = n = *nn;
   k = *kk;    /* length of eps, only needed in outside loop */
   state.pers_to_obj = pers_to_obj;      /* n vector: person i gets which object */
   state.price = price;  /* n vector: price of object j */
   state.desiremat = desirem;  /* n x n vector: desire of person i for object j */  

   /* scratch space */ 
   state.obj_to_pers = (int *) R_alloc((long) n, sizeof(int));
   state.persvalue = (double *) R_alloc((long) n, sizeof(double));
   /* With Calloc we can get extra memory in addition to Rs possibilities
      (shouldn't be needed here), but then we need to "Free" it  */
   /*   state.flowmatrix = (int *) R_alloc((long) (n1 * n2), sizeof(int));
   state.arcmatrix = (int *) R_alloc((long) (n1 * n2), sizeof(int));
   state.collectvals = (int *) R_alloc((long) (n1 * n2), sizeof(int)); */

   for (r = 0; r < k; r++) {
     state.epsbid = eps[r];
     state.nofassigned = 0;
     for (j = 0; j < n; j++) {
       state.pers_to_obj[j] = -1;
       state.obj_to_pers[j] = -1;
     }
     /* At start everything is unassigned */

     while (state.nofassigned < n) {
       /*  printit(&state); */
       R_CheckUserInterrupt();
       for (i = 0; i < n; i++) {
         if (state.pers_to_obj[i] == -1) {
           bid(&state, i);  /* bid does assigning and unassigning and changes nofassigned */
         }
       }
     }
   }

}


/* ------------ Functions called by auction ------------------------- */

void bid(State *state, int person) {
  int j;
  int n;
  int bidfor, oldpers;
  double bidamount;

  n = state->n;
  for (j = 0; j < n; j++) {
    state->persvalue[j] = DESIRE(person,j,state,n) - state->price[j];
  }
  bidfor = arrayargmax(state->persvalue, n);
  bidamount = state->persvalue[bidfor] - arraysec(state->persvalue,n,bidfor) + state->epsbid;
  /* here we get a float result, the rest are int results */
  oldpers = state->obj_to_pers[bidfor];
  if (oldpers == -1) {
    state->nofassigned++; 
  }
  else {
    state->pers_to_obj[oldpers] = -1; 
  }
  state->pers_to_obj[person] = bidfor;
  state->obj_to_pers[bidfor] = person;
  state->price[bidfor] = state->price[bidfor] + bidamount;
}


/* Sum of double array */
double arraysum(double *a, int n) {
   int i;
   double asum = 0.0;
   for (i = 0; i < n; i++)
      asum += a[i];
   return(asum);
}

/* Gives first index that maximizes array */
int arrayargmax(double *a, int n) {
  int i, arg;
  double amax;
  arg = 0;
  amax = a[0];
  for (i = 1; i < n; i++)
    if (a[i] > amax) {
    arg = i;
    amax = a[i];
  }
  return(arg);
}

/* Maximal element of a non-negative double array */
double arraymax(double *a, int n) {
  int i;
  double amax;
  amax = a[0];
  for (i = 1; i < n; i++)
    if (a[i] > amax) amax = a[i];
  return(amax);
}

/* Second largest element of a non-negative integer array
   knowing the largest is at index arg */
double arraysec(double *a, int n, int arg) {
  int i;
  double amax;
  if (arg > 0) amax = a[0];
  else amax = a[1]; 
  for (i = 0; i < arg; i++)
    if (a[i] > amax) amax = a[i];
  for (i = arg+1; i < n; i++)
    if (a[i] > amax) amax = a[i]; 
  return(amax);
}

/* Minimal element of a non-negative double array */
double arraymin(double *a, int n) {
  int i;
  double amin;
  amin = a[0];
  for (i = 0; i < n; i++)
    if (a[i] < amin) amin = a[i];
  return(amin);
}


/* void printit(State *state)
{
  int i=0,n=0;

  n = state->n;

  Rprintf("Current state: \n");

  Rprintf("nofassigned:  %d \n", state->nofassigned);

  Rprintf("pers_to_obj:  ");
  for (i = 0; i < n; i++) {
    Rprintf("%d ", state->pers_to_obj[i]);
  }
  Rprintf("\n");

  Rprintf("obj_to_pers:  ");
  for (i = 0; i < n; i++) {
    Rprintf("%d ", state->obj_to_pers[i]);
  }
  Rprintf("\n");

  Rprintf("price:  ");
  for (i = 0; i < n; i++) {
    Rprintf("%2.9lf ", state->price[i]);
  }
  Rprintf("\n");

  Rprintf("persvalue:  ");
  for (i = 0; i < n; i++) {
    Rprintf("%2.9lf ", state->persvalue[i]);
  }
  Rprintf("\n\n\n");
}
*/
