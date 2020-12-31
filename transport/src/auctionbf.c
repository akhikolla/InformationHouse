/*

   auctionbf.c
   (was /Users/schumi/Documents/Research/Papers_present/Wasserstein/Code/auction4.c)

   $Revision: 0.4 $   $Date: 2013/09/08 14:43:00 $

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
  int backwards; /* 0 if we should do forward auction, 1 if we should do backward auction */ 
  int nofassigned;  /* number of assigned persons */
  int *pers_to_obj;  /* -1 means unassigned */
  int *obj_to_pers;  /* -1 means unassigned */           
  double *price;   
  double *profit;
  int *desiremat;        /* matrix of desires */ 
  double *persvalue; /* desire minus price of current person in forward phase */
  double *objvalue; /* desire minus profit of current object in reverse phase */
                     /* last three only used in bid, but maybe better
                        to reserve memory once and for all */
} State;

#define DESIRE(I,J,STATE,NVALUE) ((STATE)->desiremat)[(NVALUE) * (J) + (I)]
#define DESIREMAIN(I,J,STATE,NVALUE) ((STATE).desiremat)[(NVALUE) * (J) + (I)]
#define MIN(A,B) ((A)<(B) ? (A) : (B))

void bidbf(State *state, int person);
void lurebf(State *state, int obj);
/* double arraysum(double *a, int n);
int arrayargmax(double *a, int n);
double arraymax(double *a, int n);
double arraysec(double *a, int n, int arg);
double arraymin(double *a, int n); */
/* void printit(State *state); */



/* ------------ The main function ----------------------------- */

void auctionbf(int *desirem, int *nn, int *pers_to_obj, double *price, double *profit, int *kk, double *eps)
{
   int i,j,r; /* indices */
   int k,n;
   State state;

   /* inputs */
   state.n = n = *nn;
   k = *kk;    /* length of eps, only needed in outside loop */
   state.pers_to_obj = pers_to_obj;      /* n vector: person i gets which object */
   state.price = price;    /* n vector: price of object j */
   state.profit = profit;  /* n vector: profit of person i */
   state.desiremat = desirem;  /* n x n vector: desire of person i for object j */  

   /* scratch space */ 
   state.obj_to_pers = (int *) R_alloc((long) n, sizeof(int));
   state.persvalue = (double *) R_alloc((long) n, sizeof(double));
   state.objvalue = (double *) R_alloc((long) n, sizeof(double));
   /* With Calloc we can get extra memory in addition to Rs possibilities
      (shouldn't be needed here), but then we need to "Free" it  */
   /*   state.flowmatrix = (int *) R_alloc((long) (n1 * n2), sizeof(int));
   state.arcmatrix = (int *) R_alloc((long) (n1 * n2), sizeof(int));
   state.collectvals = (int *) R_alloc((long) (n1 * n2), sizeof(int)); */

   /* Prices start at what the R-function supplied (usually 0) */
   /* Profits are set to the rowwise max that satisfies eps-CS */
   for (i = 0; i < n; i++) {
     for (j = 0; j < n; j++) {
       state.persvalue[j] = DESIREMAIN(i,j,state,n);
     }
     state.profit[i] = arrayargmax(state.persvalue, n);
   }

   for (r = 0; r < k; r++) {
     state.backwards = 0;
     state.epsbid = eps[r];
     /* At start everything is unassigned */
     state.nofassigned = 0;
     for (j = 0; j < n; j++) {
       state.pers_to_obj[j] = -1;
       state.obj_to_pers[j] = -1;
     }
     
     while (state.nofassigned < n) {
       /*  printit(&state); */
       R_CheckUserInterrupt();
       if (state.backwards == 0) {
         /* printit(&state); */
         for (i = 0; i < n; i++) {
           if (state.pers_to_obj[i] == -1) {
             /* Rprintf("Bid \n"); */
             bidbf(&state, i);  /* bid does assigning and unassigning and changes nofassigned */
           }
         }
       } else {
         /* printit(&state); */
         for (j = 0; j < n; j++) {
           if (state.obj_to_pers[j] == -1) {
             /* Rprintf("Lure \n"); */        
             lurebf(&state, j);  /* lure does assigning and unassigning and changes nofassigned */
           }
         }
       }

     }   /* eof while */
   }     /* eof eps-scaling for-loop */
}


/* ------------ Functions called by auction ------------------------- */

void bidbf(State *state, int person) {
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
    state->backwards = 1;
  }
  else {
    state->pers_to_obj[oldpers] = -1; 
  }
  state->pers_to_obj[person] = bidfor;
  state->obj_to_pers[bidfor] = person;
  state->price[bidfor] = state->price[bidfor] + bidamount;
  /* new forward/reverse auction algo */
  state->profit[person] = DESIRE(person,bidfor,state,n) - state->price[bidfor];
}


/* like bidbf, but for reverse auction */
void lurebf(State *state, int obj) {
  int i;
  int n;
  int lurepno, oldobj;
  double lureamount;

  n = state->n;
  for (i = 0; i < n; i++) {
    state->objvalue[i] = DESIRE(i,obj,state,n) - state->profit[i];
  }
  lurepno = arrayargmax(state->objvalue, n);
  lureamount = state->objvalue[lurepno] - arraysec(state->objvalue,n,lurepno) + state->epsbid;
  /* here we get a float result, the rest are int results */
  oldobj = state->pers_to_obj[lurepno];
  if (oldobj == -1) {
    state->nofassigned++;
    state->backwards = 0; 
  }
  else {
    state->obj_to_pers[oldobj] = -1; 
  }
  state->obj_to_pers[obj] = lurepno;
  state->pers_to_obj[lurepno] = obj;
  state->profit[lurepno] = state->profit[lurepno] + lureamount;
  /* new forward/reverse auction algo */
  state->price[obj] = DESIRE(lurepno,obj,state,n) - state->profit[lurepno];
}


/* void printit(State *state)
{
  int i=0,n=0;

  n = state->n;

  Rprintf("Current state: \n");

  Rprintf("backwards: %d \n", state->backwards);

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

  Rprintf("profit:  ");
  for (i = 0; i < n; i++) {
    Rprintf("%2.9lf ", state->profit[i]);
  }
  Rprintf("\n");

  Rprintf("persvalue:  ");
  for (i = 0; i < n; i++) {
    Rprintf("%2.9lf ", state->persvalue[i]);
  }
  Rprintf("\n");

  Rprintf("objvalue:  ");
  for (i = 0; i < n; i++) {
    Rprintf("%2.9lf ", state->objvalue[i]);
  }
  Rprintf("\n");

  Rprintf("\n\n\n");
}
*/ 
