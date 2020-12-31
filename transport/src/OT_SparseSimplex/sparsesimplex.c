/*
   WARNINGS: It seems the degenerate case is not catered for!?
   In view of Luenberger, p.50 I think we should
*/


/*

   sparsesimplex.c

   $Revision: 0.0.1 $   $Date: 2017/06/24 23:43:00 $

   By Dominic Schuhmacher, based on revsimplex.c 0.3.1; original version ported
   from java code for class Shortlist2014 by Carsten Gottschlich
   (the present file is without shortlist components)
   (CG: in the comments refers to correspondences with original code,
   mainly used to give the original variable names where I have changed them) 

*/


/*  Call by

   .C("revsimplex", as.integer(m), as.integer(n), as.integer(apos), as.integer(bpos),
	          as.double(dd), assignment = as.integer(initassig), basis = as.integer(initbasis), startgiven = as.integer(startgiven),
	          DUP=TRUE, PACKAGE="transport")

*/

/* About using long instead of int below (from R-ext, p.17):
Do be very careful with passing arguments between R, C and FORTRAN code. In particular, long in C will be 32-bit on most R platforms (including those mostly used by the CRAN maintainers), but 64-bit on many modern Unix and Linux platforms. It is rather unlikely that the use of long in C code has been thought through: if you need a longer type than int you should use a configure test for a C99 type such as int_fast64_t (and failing that, long long) and typedef your own type to be long or long long, or use another suitable type (such as size_t). Note that integer in FORTRAN corresponds to int in C on all R platforms.

On my system int is 32 bit (I guess on every system), long is 64 bit (and longlong also).
 */
 
#include"sparsesimplex.h"
#include"OT_SparseSimplex/sparsebasicfeasible.h"

#include <math.h>
#include <R.h>
#include <R_ext/Utils.h>

/* R_ext/Utils.h is included to provide the function   void R_CheckUserInterrupt(void) 

  Maybe useful at some point:
  #include <Rmath.h>
  #include <R_ext/Lapack.h>

  Usually not needed:
  #include <stdio.h>
  #include <stdlib.h>
*/

#define THRESHOLD (-1.0e-6)
/* seems a bit far from zero to me */

/* Note: the State structure includes many "temporary variables", because I suspect
   that the many calls to Calloc inside the loops would slow down the programm
   AND above all because this simplifies debugging (just print state).
   On the other hand a more dynamic allocation of the byrow / bycol stuff
   could lead to a substantial reduction of memory requirements */
typedef struct State {
  int m, n;   /* lengths of mass vectors (CG: p, c) */
  int *a, *b;   /* mass vectors (CG: producer, CG: consumer) */ 
  double *costm;    /* m x n matrix of costs */
  int **channels_byrow;  /* where the mass can flow: new and sparser data structure
                            adapt for basis_byrow as well */
  int *channels_byrow_over; /* first index that contains garbage */
  int *assignment;   /* m x n matrix giving the mass that is assigned from i to j */

  int *basis;    /* m x n (CG: isBasis, my philosophy is that it should be is_basis_vector) */
                 /* I don't use isInitialBasis */

  int *basis_byrow;         /* m x n (CG: producerBasis) */
  int *basis_byrow_over;   /* first index that contains garbage (CG: producerBasisCount) */
  int *basis_bycol;         /* n x m (CG: consumerBasis) */
  int *basis_bycol_over;   /* first index that contains garbage (CG: consumerBasisCount) */ 
  /* basis_byrow gives the col-indices of basis vectors by row
     i-th row contains at positions 0 to basis_byrow_over-1 the col-indices of the basis vectors 
     in i-th row of basis. basis_bycol likewise
  */

  int indi;      /* i and j index of basic variable to add or remove from basis*/
  int indj;      /* (CG: indexU, indexV) */
  int maxdim;    /* max(m,n) (CG: maxLen) */
  int iter;   /* number of the current iteration */ 

  /* For computing u/v and determinung new basic variable */
  int next_row;    /* next row for which relative costs are considered in rowmostnegative
                      pivot finding, is not temporary but persists over iterations
                      (CG currentProducer) */
  double *u, *v;   /* dual variables ("prices") */
  int *is_computed_u, *is_computed_v;   /* which entries of u and v have been newly
                                        computed (0/1) */
  int *list;    /* indices (without distinction if row or column index) to which prices
                     u, v have been assigned */   
  int *is_row;    /* marks the row indices in the above list (1)              (CG: isProd) */
  int over;       /* the first pos in list that contains garbage
                     (sometimes last pos. that still makes sense) (CG: count) */
  /*  int curr;          up to which position has list been processed             (CG: done) */

  /* For finding circle and redistributing mass */
  int *circlea;   /* length m+n (-1 would be enough) */
  int *circleb;   /* gives the index pairs in the circle (or path while under construction) */
  int circ_over;  /* the last pos in the circle list that contains something useful */
  int *candlist; /* list of promissing row or col indices for extending the path */
  int *rem_curr;
  int *rem_next_branch;
  int *rem_do_rowscan;
} State;

/* We use col major convention (non-standard in C)
   Note that: col major is used in primaldual.c
              col major is standard in R              */
/* #define COSTM(I,J,STATE,NROW) ((STATE)->costm)[(NROW) * (J) + (I)]
#define ASSIG(I,J,STATE,NROW) ((STATE)->assignment)[(NROW) * (J) + (I)]
#define BASIS(I,J,STATE,NROW) ((STATE)->basis)[(NROW) * (J) + (I)]
#define IBASIS(I,J,STATE,NROW) ((STATE)->initialbasis)[(NROW) * (J) + (I)]
#define BYROW(I,J,STATE,NROW) ((STATE)->basis_byrow)[(NROW) * (J) + (I)]
#define BYCOL(I,J,STATE,NROW) ((STATE)->basis_bycol)[(NROW) * (J) + (I)]
*/
#define MAT(NAME, I, J) (state->NAME)[state->m * (J) + (I)]  // col-major as in R !
// there is probably room for improvement here: several steps in the algo work row-wise
// that C is row-major has no influence because we make C go col-major here
#define TMAT(NAME, I, J) (state->NAME)[state->n * (J) + (I)]
/* letzteres für basis_bycol, die n*m ist */
#define MAX(A,B) ((A)>(B) ? (A) : (B))
#define MIN(A,B) ((A)<(B) ? (A) : (B))

int spa_update_transport_rowmostneg(State *state);
/*void init_channels_byrow(int *channelm, State *state);  to transform matrix channel representation
                                                           into ragged matrix channel representation */
void spa_init_sparse_solution(State *state);
void spa_init_solution(State *state);
void spa_init_helpers(State *state);          /*  (CG: initGraph, and more) */
int spa_new_basic_variable_rowmostneg(State *state);
void spa_find_circle(State *state);
void spa_move_mass(State *state);
void spa_add_to_basis(State *state);
void spa_remove_from_basis(State *state);

/* double compute_total_costs();
int get_iterations();
int count_flows();
int count_basis_entries(); */

/* for debugging */
void spa_printit(State *state);
void spa_pricediag(State *state);
void spa_circlediag(State *state);
void spa_printfvec(int n, double *a);
void spa_printvec(int n, int *a);

/* for inner call where sparsesimplex is used for getting starting solution */
int sparsesimplex2(int mm, int nn, int *a, int *b, double *costm, int *channels_byrow_over, int **channels_byrow,
		   int *assignment, int *basis, double *u, double *v, int startgiven, int dense);
int spa_update_transport_firstavail(State *state);
int spa_new_basic_variable_firstavail(State *state);

/* ------------ The main function ----------------------------- */

int sparsesimplex(int mm, int nn, int *a, int *b, double *costm, int *channels_byrow_over, int **channels_byrow,
		  int *assignment, int *basis, double *u, double *v, int startgiven, int dense)
/* inputs */
/* int mm, nn;*/
/* int *a, *b;*/
/* double *costm;*/
/* int *channels_byrow_over;   // „Vektor“ der Länge *mm der Nachbarschaftsgrößen der Quellen (Anzahl „Channels“) */
/* int **channels_byrow;       // Pointer von *mm Pointern auf (ersten Eintrag)*/
/*                             // von Vektoren der Länge channels_byrow_over[0],…,channels_byrow_over[mm-1]*/
/* int *startgiven;*/
/* outputs */
/* int *assignment;*/
/* int *basis;*/
{
  
  int i,j,m,n,maxdim;
  State state;
  int is_optimal;

  /* inputs/outputs */
  state.m = m = mm;
  state.n = n = nn;
  state.a = a;
  state.b = b;
  state.costm = costm;
  /*  state.channels_byrow_over = nnchannel; */

  state.assignment = assignment;
  state.basis = basis;

  state.maxdim = maxdim = MAX(m,n);
  state.iter = 0;
  state.next_row = 0;

  state.channels_byrow_over = channels_byrow_over;
  state.channels_byrow = channels_byrow;

  /*  for (int i=0; i < state.m; i++) {
    spa_printvec(state.channels_byrow_over[i], state.channels_byrow[i]);
  }
  Rprintf("\n"); */
  
  /* scratch space */
  /* I doubt that casting to long does something (positive) */
  state.basis_byrow = (int *) R_alloc((long) (m * n), sizeof(int));
  state.basis_byrow_over = (int *) R_alloc((long) m, sizeof(int));
  state.basis_bycol = (int *) R_alloc((long) (n * m), sizeof(int));
  state.basis_bycol_over = (int *) R_alloc((long) n, sizeof(int));
  state.u = (double *) R_alloc((long) m, sizeof(double));
  state.v = (double *) R_alloc((long) n, sizeof(double));
  state.is_computed_u = (int *) R_alloc((long) m, sizeof(int));
  state.is_computed_v = (int *) R_alloc((long) n, sizeof(int));
  state.list = (int *) R_alloc((long) (m + n), sizeof(int));
  state.is_row = (int *) R_alloc((long) (m + n), sizeof(int));
  state.circlea = (int *) R_alloc((long) (m + n), sizeof(int));
  state.circleb = (int *) R_alloc((long) (m + n), sizeof(int));
  state.candlist = (int *) R_alloc((long) maxdim, sizeof(int));
  state.rem_curr = (int *) R_alloc((long) maxdim, sizeof(int));
  state.rem_next_branch = (int *) R_alloc((long) maxdim, sizeof(int));
  state.rem_do_rowscan = (int *) R_alloc((long) maxdim, sizeof(int));
  
  if (startgiven == 0) {
    if (dense == 0) {
      spa_init_sparse_solution(&state); 
    } else {
      spa_init_solution(&state); // by modified row minimum rule (much faster)
    }
  } 
  spa_init_helpers(&state);

  is_optimal = 0;

  /* Maybe add condition for count_update */
  while (is_optimal == 0) {
    R_CheckUserInterrupt();
    state.iter++;
    is_optimal = spa_update_transport_rowmostneg(&state);
  }

  for (i = 0; i < m; i++) {
    u[i] = state.u[i];
  }
  for (j = 0; j < n; j++) {
    v[j] = state.v[j];
  }

  return 0;
}


/*    */
int spa_update_transport_rowmostneg(State *state)
{
  int found;

  found = spa_new_basic_variable_rowmostneg(state);
  /* implicitly returns state->indi and state->indj, i.e.
     i and j index of new basic variable */

  if (found == 0) {
    return(1);
  }

  spa_add_to_basis(state);

  spa_find_circle(state);

  spa_move_mass(state);

  spa_remove_from_basis(state);

  /* printit(state); */
  /*  error("test interrupt");  */

  return(0);  
}


/* Bernhard's program gives this info directly */
/* void init_channels_byrow(int *channelm, State *state)
{
  int i, j, k;

  state->channels_byrow = (int **) R_alloc(state->m, sizeof(int *));    
  for (i = 0; i < state->m; i++) {
      state->channels_byrow[i] = (int *) R_alloc(state->channels_byrow_over[i], sizeof(int));
      k = 0;
  for (j = 0; j < state->n; j++) {
    if (channelm[state->m * j + i] == 1) {
      state->channels_byrow[i][k] = j;
      k++;
    }
  }
  }
  } */



 // Currently very inefficient because we need the space for everything twice
 // In the long run: reserve one extra dimension for everything before calling sparsesimplex
void spa_init_sparse_solution(State *state)
{
  int i, j, k, allok;
  int totsum, rowsum, colsum;
  int m2, n2;
  int *a2, *b2;
  double *costm2;
  int *channels_byrow_over2;   // „Vektor“ der Länge *mm der Nachbarschaftsgrößen der Quellen (Anzahl „Channels“) 
  int **channels_byrow2;       // Pointer von *mm Pointern auf (ersten Eintrag)
                              // von Vektoren der Länge channels_byrow_over[0],…,channels_byrow_over[mm-1]
  /* outputs */
  int *assignment2;
  int *basis2;
  double *u2, *v2;

  // Realloc could be a good solution here
  // Freeing is not 100% fool proof: user interrupt (and error) does not Free the memory
  m2 = state->m + 1;
  n2 = state->n + 1;
  a2 = (int *) Calloc((long) m2, int);
  b2 = (int *) Calloc((long) n2, int);
  costm2 = (double *) Calloc((long) (m2 * n2), double);
  channels_byrow_over2 = (int *) Calloc((long) m2, int);
  channels_byrow2 = (int **) Calloc((long) m2, int *);
  assignment2 = (int *) Calloc((long) (m2 * n2), int);
  basis2 = (int *) Calloc((long) (m2 * n2), int);
  u2 = (double *) Calloc((long) m2, double);
  v2 = (double *) Calloc((long) n2, double);
  for (i = 0; i < state->m; i++) {
    channels_byrow2[i] = (int *) Calloc(state->channels_byrow_over[i]+1, int);      
  }
  channels_byrow2[state->m] = (int *) Calloc(n2, int);
  
  // Rprintf("After memory allocation in init_sparse_solution \n");
  // copy what we can from actual problem
  totsum = 0;
  for (i = 0; i < state->m; i++) {
    totsum += state->a[i];  // dangerous, we actually do not need that sum(a) is a representable number otherwise
                            // introduce a check that it is representable
    a2[i] = state->a[i];
    channels_byrow_over2[i] = state->channels_byrow_over[i]+1;
    u2[i] = 0; 
  } 
  for (j = 0; j < state->n; j++) {
    b2[j] = state->b[j];
    v2[j] = 0;
  }
  for (i = 0; i < state->m; i++) {    
  for (j = 0; j < state->n; j++) {
    costm2[m2 * j + i] = 0;   //  change m * j +i to m * i + j globally
    assignment2[m2 * j + i] = 0;
    basis2[m2 * j + i] = 0;
  }
  for (k = 0; k < state->channels_byrow_over[i]; k++) {
    channels_byrow2[i][k] = state->channels_byrow[i][k];
  }
  }
  u2[state->m] = 0;
  v2[state->n] = 0;

  // Rprintf("totsum: %d \n", totsum);

  
  // extra values for extended problem
  a2[state->m] = totsum;
  b2[state->n] = totsum;
  channels_byrow_over2[state->m] = n2;
  
  for (i = 0; i < state->m; i++) {  
    costm2[m2 * state->n + i] = 1.0;   //  change m * j +i to m * i + j globally
    assignment2[m2 * state->n + i] = state->a[i];
    basis2[m2 * state->n + i] = 1;
    channels_byrow2[i][state->channels_byrow_over[i]] = state->n;
  }  
  for (j = 0; j < state->n; j++) {  
    costm2[m2 * j + state->m] = 1.0;   //  change m * j +i to m * i + j globally
    assignment2[m2 * j + state->m] = state->b[j];
    basis2[m2 * j + state->m] = 1;
    channels_byrow2[state->m][j] = j;
  }
  
  costm2[m2 * state->n + state->m] = 0.0;
  assignment2[m2 * state->n + state->m] = 0;
  basis2[m2 * state->n + state->m] = 1;
  channels_byrow2[state->m][state->n] = state->n;

  //  for (i = 0; i < m2; i++) {
  //  printvec(channels_byrow_over2[i], channels_byrow2[i]);
  //  }

  i = sparsesimplex2(m2, n2, a2, b2, costm2, channels_byrow_over2, channels_byrow2,
		assignment2, basis2, u2, v2, 1, 0);
  // the only difference is that new basic variable is determined by firstavail
  // rather than modrowmin. There is much to be gained by replacing this
  // by a procedure that finds a starting solution in a more targeted way 

  if (assignment2[m2 * state->n + state->m] != totsum) {
    error("sparse problem has no feasible solution");
  }


  // cut off extra row and column and count remaining basis vectors
  totsum = 0;
  for (i = 0; i < state->m; i++) {    
  for (j = 0; j < state->n; j++) {
    MAT(assignment,i,j) = assignment2[m2 * j + i];   
    MAT(basis,i,j) = basis2[m2 * j + i];
    if (basis2[m2 * j + i] == 1) {
      totsum++;
    }
  }
  }
  
  // check basis
    for (i = 0; i < state->m; i++) {
    for (j = 0; j < state->n; j++) {
      allok = 0;
      if (MAT(basis,i,j) == 1) {
	for (k = 0; k < state->channels_byrow_over[i]; k++) {
	  if (state->channels_byrow[i][k] == j) {
            allok = 1;
	  }
	}
	if (allok == 0) {
          error("starting solution for sparse simplex has illegal basis entry");
	}
      }
    }
    }

    for (i = 0; i < state->m; i++) {
      rowsum = 0;
      for (j = 0; j < state->n; j++) {
        rowsum += MAT(assignment,i,j);
      }
      if (rowsum != state->a[i]) {
        Rprintf("assignment:  \n");
        for (i = 0; i < state->m; i++) {
	Rprintf("%d:  ", state->a[i]);  
        for (j = 0; j < state->n; j++) {
        Rprintf("%d ", MAT(assignment,i,j));
        }
        Rprintf("\n");
        }
        Rprintf("\n");
        error("starting assignment for sparse simplex has wrong rowsum");
      }
    }

    for (j = 0; j < state->n; j++) {
    colsum = 0;
    for (i = 0; i < state->m; i++) {
      colsum += MAT(assignment,i,j);
    }
    if (colsum != state->b[j]) {
      error("starting assignment for sparse simplex has wrong colsum");
    }
    }

  if (totsum != state->m + state->n -1) {
    /* Rprintf("basis:  \n");
    for (i = 0; i < state->m; i++) {
    for (j = 0; j < state->n; j++) {
      Rprintf("%d ", MAT(basis,i,j));
    }
    Rprintf("\n");
    }             */
    // Rprintf("\n");
    Rprintf("%d basis entries, should be %d \n",totsum, state->m + state->n - 1); 
    Rprintf("Trying to fix this... ");
    dedewithchannels(state->m, state->n, totsum, state->basis,
		     state->channels_byrow_over, state->channels_byrow);
    Rprintf("Success! \n");
  }
  
  Free(a2);
  Free(b2);
  Free(costm2);
  Free(channels_byrow_over2);
  for (i = 0; i < m2; i++) {
    Free(channels_byrow2[i]);
  }
  Free(channels_byrow2);
  Free(assignment2);
  Free(basis2);
  Free(u2);
  Free(v2); 
}



/* initializes assignment and basis by using the "modified row minimum rule" 
   (and taking care of potential degeneracies:
   create an additional basis vector between i where degeneracy occurs and any j
   that has still mass left, but has no more mass left when next degeneracy occurs) */
void spa_init_solution(State *state)
{
  int i, j, jstar;
  int mass;
  int degennum, degeni, degenj;
  int m,n,numleft;
  double mini;
  int *aleft, *bleft;   
  int *aisleft, *bisleft;    
  int *degenisj;

  m = state->m;
  n = state->n;

  aleft = (int *) Calloc((long) m, int);
  bleft = (int *) Calloc((long) n, int);
  aisleft = (int *) Calloc((long) m, int);
  bisleft = (int *) Calloc((long) n, int);
  degenisj = (int *) Calloc((long) n, int);

  for (i = 0; i < m; i++) {
  for (j = 0; j < n; j++) {
    MAT(assignment,i,j) = 0;
    MAT(basis,i,j) = 0;
  }
  }

  /* in all such instances it might be more efficient
     to copy first one array and then the other (needs only increments of 1) */
  for (i = 0; i < m; i++) {
    aleft[i] = state->a[i];
    aisleft[i] = 1;
  }

  for (j = 0; j < n; j++) {
    bleft[j] = state->b[j];
    bisleft[j] = 1;
    //    degenisj[j] = 0;
  } 

  degennum = 0;
  numleft = m;
  while (numleft > 0) {
    R_CheckUserInterrupt();
    for (i = 0; i < m; i++) {      
      if (aisleft[i] == 1) {
	mini = R_PosInf;
	for (j = 0; j < n; j++) {
	  if (bisleft[j] == 1) {
            if (MAT(costm,i,j) < mini) {
              mini = MAT(costm,i,j);
	      jstar = j;
	    }
	  }
	}
        mass = MIN(aleft[i],bleft[jstar]);
	aleft[i] -= mass;
        bleft[jstar] -= mass;
        MAT(assignment,i,jstar) = mass;
        MAT(basis,i,jstar) = 1;
        if (aleft[i] == 0) {
          aisleft[i] = 0;
	  numleft--;
        }
        if (bleft[jstar] == 0) {
          bisleft[jstar] = 0;
        }
        if (aleft[i] == 0 && bleft[jstar] == 0) {   /* degenerate case */
          if (degennum > 0) {
	    for (j = 0; j < n; j++) {
              if (degenisj[j] == 1 && bisleft[j] == 0) {
                degenj = j;
	        MAT(basis, degeni, degenj) = 1;
		break;
	      }
	    }
	  }
	  degeni = i;
	  for (j = 0; j < n; j++) {
            degenisj[j] = bisleft[j];
          } 
          degennum++;
	}         
      }
    }
  }
  /* if (degennum > 0) {
    warning("Starting solution is degenerate. Nothing to worry about!");
  } */
  Free(aleft);
  Free(bleft);
  Free(aisleft);
  Free(bisleft);
  Free(degenisj);
}


/* Initializes the various helper variables in state */
void spa_init_helpers(State *state)
{
  int i, j;

  /* initializations to zero */
  /* Some of this might be superfluous, but
     R-int says setting to zero by R_alloc cannot be trusted */
  for (i = 0; i < state->m; i++) {
    state->basis_byrow_over[i] = 0;
  }
  for (j = 0; j < state->n; j++) {
    state->basis_bycol_over[j] = 0;
  }
  /*   for (k = 0; k < m+n; k++) {
    state->circlea[k] = 0;
    state->circleb[k] = 0;
  }
  for (k = 0; k < maxdim; k++) {
    state->rem_curr[k] = 0;
    state->rem_next_branch[k] = 0;
    state->rem_rowscan[k] = 0;
  }
  */  

  /* more interesting initializations */
  for (i = 0; i < state->m; i++) {
  for (j = 0; j < state->n; j++) {
    if (MAT(basis,i,j) == 1) {
      MAT(basis_byrow,i,state->basis_byrow_over[i]) = j;
      state->basis_byrow_over[i]++;
      TMAT(basis_bycol,j,state->basis_bycol_over[j]) = i;
      state->basis_bycol_over[j]++;
    }
  }
  }

}


int spa_new_basic_variable_rowmostneg(State *state)
{
  int i,j,l,newind;
  int count;
  double bestred;
  double redcost; 
  /* "local copies" of the essential objects in state */
  int m,n;  
  int *list;    /* indices (without distinction if row or column index) to which prices
                     u, v have been assigned */   
  int *is_row;    /* marks the row indices in the above list (1)              (CG: isProd) */
  int over=0;       /* the first pos in list that contains garbage
                     (sometimes last pos. that still makes sense) (CG: count) */
  int curr=0;       /* up to which position has list been processed             (CG: done) */
  double *u, *v;   /* dual variables ("prices") */
  int *is_computed_u, *is_computed_v;   /* which entries of u and v have been newly                                                            computed (0/1) */
  m = state->m;
  n = state->n;
  list = state->list;
  is_row = state->is_row;
  u = state->u;
  v = state->v;
  is_computed_u = state->is_computed_u;
  is_computed_v = state->is_computed_v;

  for (i = 0; i < m; i++) {
    is_computed_u[i] = 0;
  }
  for (j = 0; j < n; j++) {
    is_computed_v[j] = 0;
  }

  /* compute u, v based on 
     costm[i,j] = u[i]+v[j] 
     whenever i,j indices of a basis-vector */
  
  u[0] = 0.0;
  is_computed_u[0] = 1;
  list[over] = 0;
  is_row[over] = 1;
  over++;

  while (curr < over) {
    if (is_row[curr] == 1) {
      for (j = 0; j < state->basis_byrow_over[list[curr]]; j++) {
        newind = MAT(basis_byrow, list[curr], j);
        if (is_computed_v[newind] == 0) {
          v[newind] = MAT(costm, list[curr], newind) - u[list[curr]];
          is_computed_v[newind] = 1;
          list[over] = newind;
          is_row[over] = 0;
          over++;
        }
      }
    }
    else {
      for (i = 0; i < state->basis_bycol_over[list[curr]]; i++) {
        newind = TMAT(basis_bycol, list[curr], i);
        if (is_computed_u[newind] == 0) {
          u[newind] = MAT(costm, newind, list[curr]) - v[list[curr]];
          is_computed_u[newind] = 1;
          list[over] = newind;
          is_row[over] = 1;
          over++;
        }
      }      
    }
    curr++;
  }

  
  bestred = 0.0;
  /* we break the count-loop and return to caller
     if the THRESHOLD is understepped in any row;
     state->next_row contains the index of the next row not considered */
  for (count = 0; count < m; count++) {
    i = state->next_row;
    for (l = 0; l < state->channels_byrow_over[i]; l++) {
      j = state->channels_byrow[i][l];
      /* we don't check that (i,j) is not 
         already in basis; theoretically the redcost is 0 in that case anyway;
         practically it should be zero also because the dual costs
         come exactly from taking diffs at all basis vectors
         remains the question: is computing n+m-1 times extra the redcost
         and comparing to bestred really more efficient than checking n*m times
         if (i,j) is 1 in state->basis */
      /* In the sparse version I changed this as we are checking only
         much fewer than n*m times if (i,j) is 1 in state->basis */
      if (MAT(basis,i,j) == 0) {	
        redcost = MAT(costm,i,j) - u[i] - v[j];
        if (redcost < bestred) {
          bestred = redcost;
          state->indi = i;
          state->indj = j;
        }
      }
    }
    state->next_row++;
    if (state->next_row == m) {
      state->next_row = 0;
    }
    if (bestred < THRESHOLD) {
      state->over = over;  /* just for diagnostics ??! */
      return(1);
    }
  }
  
  state->over = over;
  return(0);
}



void spa_find_circle(State *state)
{
  int i, j;
  int lasti, lastj;
  int candi, candj;
  int ncand;
  int curr, curr_fork, next_branch;
  int do_rowscan, finished;
  /* "local copies" of the essential objects in state */
  int indi, indj;
  int *circlea, *circleb;
  int *candlist;
  int *rem_curr;
  int *rem_next_branch;
  int *rem_do_rowscan;

  indi = state->indi;
  indj = state->indj; 

  circlea = state->circlea;
  circleb = state->circleb;
  candlist = state->candlist;
  rem_curr = state->rem_curr;
  rem_next_branch = state->rem_next_branch;
  rem_do_rowscan = state->rem_do_rowscan;   

  circlea[0] = indi;
  circleb[0] = indj;

  curr = 1;           /*  CG: countCircle  */         
  curr_fork = -1;      /*  CG: countFork    */
                       /*  slight adaptation: curr_fork is always CG's countFork - 1
                           i.e. index of current fork, -1 if no fork */ 
                      /*  index for circle; circlea[curr_fork], circleb[curr_fork]
                          gibt uns den tatsächlichen "Knoten" (Indexpaar in Basis), wo
                          die Fork stattfindet */ 
  finished = 0;

  do_rowscan = 1;     /*  CG: searchAlongTheSameProducer */
  next_branch = 0;    /*  CG: forkToTryNext    */
                      /*  this is often actually the "current branch" */ 

  while (finished == 0) {

    /* ROWSCAN */
    if (do_rowscan == 1) {
      lasti = circlea[curr-1];
      lastj = circleb[curr-1];
      ncand = 0;
      
      /* count candidates and compile candlist */
      for (j = 0; j < state->basis_byrow_over[lasti]; j++) {
        candj = MAT(basis_byrow, lasti, j);
        if (state->basis_bycol_over[candj] > 1 && candj != lastj) {
          candlist[ncand] = candj;
          /* statt nur den ersten zu nehmen
             nicht klar, ob genügend forks auftreten, dass das was bringt,
             aber kosten tut's kaum was */ 
          ncand++;
        }
      }

      if (ncand == 0) {
        /* dead end (leaf), try new branch from last fork */
        curr = rem_curr[curr_fork];
        do_rowscan = rem_do_rowscan[curr_fork];
        next_branch = rem_next_branch[curr_fork];
      }
      else if (ncand == 1) {
        /* easy, add only choice to circle */
        circlea[curr] = lasti;
        circleb[curr] = candlist[0];
        curr++;
        do_rowscan = 0;
        next_branch = 0;
      }
      else {
	/* if (ncand >= 2)  */
        if (next_branch == 0) {
          circlea[curr] = lasti;
          circleb[curr] = candlist[0];  
          
          curr_fork++;
          rem_curr[curr_fork] = curr;
     	  rem_do_rowscan[curr_fork] = 1;
          rem_next_branch[curr_fork] = 1;
          curr++;
          
          do_rowscan = 0;
          next_branch = 0;
        }
        else if (next_branch >= ncand) {
          /* dead end (all branches starting at current fork)
             try next branch */
          curr_fork--;
          curr = rem_curr[curr_fork];
          do_rowscan = rem_do_rowscan[curr_fork];
          next_branch = rem_next_branch[curr_fork];
        }
        else {
          rem_next_branch[curr_fork]++;
          circlea[curr] = lasti;
          circleb[curr] = candlist[next_branch];
          curr++;
          do_rowscan = 0;
          next_branch = 0;
        }
      }
    }
    else {

    /* COLSCAN */
      lasti = circlea[curr-1];
      lastj = circleb[curr-1];
      ncand = 0;

      /* count candidates and compile candlist */
      for (i = 0; i < state->basis_bycol_over[lastj]; i++) {
        candi = TMAT(basis_bycol, lastj, i);
        if (candi == indi && curr > 2) {
	  /* changed > 3 to > 2 on 15/10/07
             as it seems a direct loop with four station has curr = 3 here */
          /* This works because by basis triangularity there are no circles in the
             original basis (i.e. without adding the (indi,indj)-vector */
          finished = 1;
          break;  /* the for loop */
        }
        if (state->basis_byrow_over[candi] > 1 && candi != lasti) {
          candlist[ncand] = candi;
          /* statt nur den ersten zu nehmen
             nicht klar, ob genügend forks auftreten, dass das was bringt,
             aber kosten tut's kaum was */ 
          ncand++;
        }
      }
      if (finished == 1) {
        break;  /* the outer while loop, i.e. the condition for this loop
                   might as well be empty */
      }

      if (ncand == 0) {
        /* dead end (leaf), try new branch from last fork */
        curr = rem_curr[curr_fork];
        do_rowscan = rem_do_rowscan[curr_fork];
        next_branch = rem_next_branch[curr_fork];
      }
      else if (ncand == 1) {
        /* easy, add only choice to circle */
        circlea[curr] = candlist[0];
        circleb[curr] = lastj;
        curr++;
        do_rowscan = 1;
        next_branch = 0;
      }
      else {
	/* if (ncand >= 2)  */
        if (next_branch == 0) {
          circlea[curr] = candlist[0];
          circleb[curr] = lastj;  
          
          curr_fork++;
          rem_curr[curr_fork] = curr;
     	  rem_do_rowscan[curr_fork] = 0;
          rem_next_branch[curr_fork] = 1;
          curr++;
          
          do_rowscan = 1;
          next_branch = 0;
        }
        else if (next_branch >= ncand) {
          /* dead end (all branches starting at current fork)
             try next branch */
          curr_fork--;
          curr = rem_curr[curr_fork];
          do_rowscan = rem_do_rowscan[curr_fork];
          next_branch = rem_next_branch[curr_fork];
        }
        else {
          rem_next_branch[curr_fork]++;
          circlea[curr] = candlist[next_branch];
          circleb[curr] = lastj;
          curr++;
          do_rowscan = 1;
          next_branch = 0;
        }
      }
    }  /* end of COLSCAN */ 

    /*      circlediag(state);
      Rprintf("ncand: %d \n", ncand);
      Rprintf("curr_fork/next_branch/curr: %d / %d / %d \n", curr_fork, next_branch, curr); 
      Rprintf("do_rowscan: %d \n", do_rowscan);   */
      /* Rprintf("%d, ", k); */

  }  /* end of while */
  
  /*  Rprintf("\n FINISHED! \n");
  circlediag(state);
  Rprintf("ncand: %d \n", ncand);
  Rprintf("curr_fork/next_branch/curr: %d / %d / %d \n", curr_fork, next_branch, curr); 
  Rprintf("do_rowscan: %d \n \n", do_rowscan);  */

  state->circ_over = curr;  /* is always <= m+n (the "all basis vectors in the circle" case) */
}


void spa_move_mass(State *state)
{
  int k;
  int mass, newmass, argmin;
  /* "local copies" of the essential objects in state */  
  int *circlea, *circleb;
  int circ_over;

  circ_over = state->circ_over; 

  circlea = state->circlea;
  circleb = state->circleb;

  /* minimize mass over "all minus signs" */
  mass = MAT(assignment, circlea[1], circleb[1]);
  argmin = 1;    
  for (k = 3; k < circ_over; k += 2) {
    newmass = MAT(assignment, circlea[k], circleb[k]);
    if (newmass < mass) {
      mass = newmass;
      argmin = k;
    }
  }

  if (mass > 0) {
    for (k = 0; k < circ_over; k += 2) {
      MAT(assignment, circlea[k], circleb[k]) += mass;
      MAT(assignment, circlea[k+1], circleb[k+1]) -= mass;
      /* should work, because circ_over should always be even */
    }
  }
  /*  else {
    Rprintf("%d ", state->iter);
  }                               deleted 0.6 only  */

  /* index of basis vector which shall be removed */
  state->indi = circlea[argmin];
  state->indj = circleb[argmin];
}


void spa_add_to_basis(State *state)
{
  int indi, indj;

  indi = state->indi;
  indj = state->indj;  

  MAT(basis, indi, indj) = 1;

  MAT(basis_byrow, indi, state->basis_byrow_over[indi]) = indj;
  state->basis_byrow_over[indi]++;

  TMAT(basis_bycol, indj, state->basis_bycol_over[indj]) = indi;
  state->basis_bycol_over[indj]++;
}


void spa_remove_from_basis(State *state)
{
  int i,j;
  int indi, indj;

  indi = state->indi;
  indj = state->indj;  

  MAT(basis, indi, indj) = 0;

  if (state->basis_byrow_over[indi] == 1) {
    state->basis_byrow_over[indi] = 0;
  }
  else {
    for (j = 0; j < state->basis_byrow_over[indi]; j++) {
      if (MAT(basis_byrow, indi, j) == indj) {
        MAT(basis_byrow, indi, j) = MAT(basis_byrow, indi, state->basis_byrow_over[indi]-1);
        state->basis_byrow_over[indi]--;
        break;
      }
    }
  }

  if (state->basis_bycol_over[indj] == 1) {
    state->basis_bycol_over[indj] = 0;
  }
  else {
    for (i = 0; i < state->basis_bycol_over[indj]; i++) {
      if (TMAT(basis_bycol, indj, i) == indi) {
        TMAT(basis_bycol, indj, i) = TMAT(basis_bycol, indj, state->basis_bycol_over[indj]-1);
        state->basis_bycol_over[indj]--;
        break;
      }
    }
  }
  
}




void spa_printit(State *state)
{
  int i,j;
  int m,n;

  m = state->m;
  n = state->n;

  Rprintf("Current state: \n");

  Rprintf("dim:  %d  %d \n", m, n);
  Rprintf("maxdim:  %d \n", state->maxdim);

  Rprintf("a:  ");
  for (i = 0; i < m; i++) {
    Rprintf("%d ", state->a[i]);
  }
  Rprintf("\n");

  Rprintf("b:  ");
  for (j = 0; j < n; j++) {
    Rprintf("%d ", state->b[j]);
  }
  Rprintf("\n");

  Rprintf("costm:  \n");
  for (i = 0; i < m; i++) {
  for (j = 0; j < n; j++) {
    Rprintf("%2.6lf ", MAT(costm,i,j));
  }
  Rprintf("\n");
  }
  Rprintf("\n");

  Rprintf("assignment:  \n");
  for (i = 0; i < m; i++) {
  for (j = 0; j < n; j++) {
    Rprintf("%d ", MAT(assignment,i,j));
  }
  Rprintf("\n");
  }
  Rprintf("\n");

  Rprintf("basis:  \n");
  for (i = 0; i < m; i++) {
  for (j = 0; j < n; j++) {
    Rprintf("%d ", MAT(basis,i,j));
  }
  Rprintf("\n");
  }
  Rprintf("\n");

  Rprintf("basis_byrow:  \n");
  for (i = 0; i < m; i++) {
  for (j = 0; j < state->basis_byrow_over[i]; j++) {
    Rprintf("%d ", MAT(basis_byrow,i,j));
  }
  Rprintf("\n");
  }
  Rprintf("\n");

  Rprintf("basis_bycol:  \n");
  for (j = 0; j < n; j++) {
  for (i = 0; i < state->basis_bycol_over[j]; i++) {
    Rprintf("%d ", TMAT(basis_bycol,j,i));
  }
  Rprintf("\n");
  }
  Rprintf("\n");

  Rprintf("next entry in/out:  %d  %d \n", state->indi, state->indj);
  Rprintf("\n\n\n");
}


void spa_pricediag(State *state)
{
  int i;
  int m,n;

  m = state->m;
  n = state->n;

  Rprintf("u:  ");
  for (i = 0; i < m; i++) {
    Rprintf("%2.9lf ", state->u[i]);
  }
  Rprintf("\n");
  Rprintf("is_computed_u:  ");
  for (i = 0; i < m; i++) {
    Rprintf("%d ", state->is_computed_u[i]);
  }
  Rprintf("\n");

  Rprintf("v:  ");
  for (i = 0; i < n; i++) {
    Rprintf("%2.9lf ", state->v[i]);
  }
  Rprintf("\n");
  Rprintf("is_computed_v:  ");
  for (i = 0; i < n; i++) {
    Rprintf("%d ", state->is_computed_v[i]);
  }
  Rprintf("\n");

  Rprintf("list:  ");
  for (i = 0; i < m+n; i++) {
    Rprintf("%d ", state->list[i]);
  }
  Rprintf("\n");
  Rprintf("is_row:  ");
  for (i = 0; i < m+n; i++) {
    Rprintf("%d ", state->is_row[i]);
  }
  Rprintf("\n");  
}


void spa_circlediag(State *state)
{
  int i;
  int m,n;

  m = state->m;
  n = state->n;

  Rprintf("circlea: ");
  for (i = 0; i < m+n; i++) {
    Rprintf("%d ", state->circlea[i]);
  }
  Rprintf("\n");
  Rprintf("circleb: ");
  for (i = 0; i < m+n; i++) {
    Rprintf("%d ", state->circleb[i]);
  }
  Rprintf("\n");

  Rprintf("candlist:  ");
  for (i = 0; i < state->maxdim; i++) {
    Rprintf("%d ", state->candlist[i]);
  }
  Rprintf("\n");

  Rprintf("rem_curr:  ");
  for (i = 0; i < state->maxdim; i++) {
    Rprintf("%d ", state->rem_curr[i]);
  }
  Rprintf("\n");

  Rprintf("rem_next_branch:  ");
  for (i = 0; i < state->maxdim; i++) {
    Rprintf("%d ", state->rem_next_branch[i]);
  }
  Rprintf("\n");
  Rprintf("rem_do_rowscan:  ");
  for (i = 0; i < state->maxdim; i++) {
    Rprintf("%d ", state->rem_do_rowscan[i]);
  }
  Rprintf("\n");  
}


void spa_printfvec(int n, double *a)
{
  int i;

  Rprintf("\n");
  for (i = 0; i < n; i++) {
    Rprintf("%2.9lf ", a[i]);
  }
  Rprintf("\n");
}


void spa_printvec(int n, int *a)
{
  int i;

  Rprintf("\n");
  for (i = 0; i < n; i++) {
    Rprintf("%d ", a[i]);
  }
  Rprintf("\n");
  }   




/* ------------ Alternative functions for determining starting solution ----------------------------- */
/* ------------ Used for second call to sparsesimplex ----------------------------- */

int sparsesimplex2(int mm, int nn, int *a, int *b, double *costm, int *channels_byrow_over, int **channels_byrow,
		  int *assignment, int *basis, double *u, double *v, int startgiven, int dense)
/* inputs */
/* int mm, nn;*/
/* int *a, *b;*/
/* double *costm;*/
/* int *channels_byrow_over;   // „Vektor“ der Länge *mm der Nachbarschaftsgrößen der Quellen (Anzahl „Channels“) */
/* int **channels_byrow;       // Pointer von *mm Pointern auf (ersten Eintrag)*/
/*                             // von Vektoren der Länge channels_byrow_over[0],…,channels_byrow_over[mm-1]*/
/* int *startgiven;*/
/* outputs */
/* int *assignment;*/
/* int *basis;*/
{
  
  int i,j,m,n,maxdim;
  State state;
  int is_optimal;

  /* inputs/outputs */
  state.m = m = mm;
  state.n = n = nn;
  state.a = a;
  state.b = b;
  state.costm = costm;
  /*  state.channels_byrow_over = nnchannel; */

  state.assignment = assignment;
  state.basis = basis;

  state.maxdim = maxdim = MAX(m,n);
  state.iter = 0;
  state.next_row = 0;

  state.channels_byrow_over = channels_byrow_over;
  state.channels_byrow = channels_byrow;

  /*  for (int i=0; i < state.m; i++) {
    printvec(state.channels_byrow_over[i], state.channels_byrow[i]);
  }
  Rprintf("\n"); */
  
  /* scratch space */
  /* I doubt that casting to long does something (positive) */
  state.basis_byrow = (int *) R_alloc((long) (m * n), sizeof(int));
  state.basis_byrow_over = (int *) R_alloc((long) m, sizeof(int));
  state.basis_bycol = (int *) R_alloc((long) (n * m), sizeof(int));
  state.basis_bycol_over = (int *) R_alloc((long) n, sizeof(int));
  state.u = (double *) R_alloc((long) m, sizeof(double));
  state.v = (double *) R_alloc((long) n, sizeof(double));
  state.is_computed_u = (int *) R_alloc((long) m, sizeof(int));
  state.is_computed_v = (int *) R_alloc((long) n, sizeof(int));
  state.list = (int *) R_alloc((long) (m + n), sizeof(int));
  state.is_row = (int *) R_alloc((long) (m + n), sizeof(int));
  state.circlea = (int *) R_alloc((long) (m + n), sizeof(int));
  state.circleb = (int *) R_alloc((long) (m + n), sizeof(int));
  state.candlist = (int *) R_alloc((long) maxdim, sizeof(int));
  state.rem_curr = (int *) R_alloc((long) maxdim, sizeof(int));
  state.rem_next_branch = (int *) R_alloc((long) maxdim, sizeof(int));
  state.rem_do_rowscan = (int *) R_alloc((long) maxdim, sizeof(int));
  
  if (startgiven == 0) {
    if (dense == 0) {
      spa_init_sparse_solution(&state); 
    } else {
      spa_init_solution(&state); // by modified row minimum rule (much faster)
    }
  } 
  spa_init_helpers(&state);

  is_optimal = 0;

  /* Maybe add condition for count_update */
  while (is_optimal == 0) {
    R_CheckUserInterrupt();
    state.iter++;
    is_optimal = spa_update_transport_firstavail(&state);
  }

  for (i = 0; i < m; i++) {
    u[i] = state.u[i];
  }
  for (j = 0; j < n; j++) {
    v[j] = state.v[j];
  }

  return 0;
}


/*    */
int spa_update_transport_firstavail(State *state)
{
  int found;

  found = spa_new_basic_variable_firstavail(state);
  /* implicitly returns state->indi and state->indj, i.e.
     i and j index of new basic variable */

  if (found == 0) {
    return(1);
  }

  spa_add_to_basis(state);

  spa_find_circle(state);

  spa_move_mass(state);

  spa_remove_from_basis(state);

  /* printit(state); */
  /*  error("test interrupt");  */

  return(0);  
}


int spa_new_basic_variable_firstavail(State *state)
{
  int i,j,l,newind;
  int count;
  double bestred;
  double redcost; 
  /* "local copies" of the essential objects in state */
  int m,n;  
  int *list;    /* indices (without distinction if row or column index) to which prices
                     u, v have been assigned */   
  int *is_row;    /* marks the row indices in the above list (1)              (CG: isProd) */
  int over=0;       /* the first pos in list that contains garbage
                     (sometimes last pos. that still makes sense) (CG: count) */
  int curr=0;       /* up to which position has list been processed             (CG: done) */
  double *u, *v;   /* dual variables ("prices") */
  int *is_computed_u, *is_computed_v;   /* which entries of u and v have been newly                                                            computed (0/1) */
  m = state->m;
  n = state->n;
  list = state->list;
  is_row = state->is_row;
  u = state->u;
  v = state->v;
  is_computed_u = state->is_computed_u;
  is_computed_v = state->is_computed_v;

  for (i = 0; i < m; i++) {
    is_computed_u[i] = 0;
  }
  for (j = 0; j < n; j++) {
    is_computed_v[j] = 0;
  }

  /* compute u, v based on 
     costm[i,j] = u[i]+v[j] 
     whenever i,j indices of a basis-vector */
  
  u[0] = 0.0;
  is_computed_u[0] = 1;
  list[over] = 0;
  is_row[over] = 1;
  over++;

  while (curr < over) {
    if (is_row[curr] == 1) {
      for (j = 0; j < state->basis_byrow_over[list[curr]]; j++) {
        newind = MAT(basis_byrow, list[curr], j);
        if (is_computed_v[newind] == 0) {
          v[newind] = MAT(costm, list[curr], newind) - u[list[curr]];
          is_computed_v[newind] = 1;
          list[over] = newind;
          is_row[over] = 0;
          over++;
        }
      }
    }
    else {
      for (i = 0; i < state->basis_bycol_over[list[curr]]; i++) {
        newind = TMAT(basis_bycol, list[curr], i);
        if (is_computed_u[newind] == 0) {
          u[newind] = MAT(costm, newind, list[curr]) - v[list[curr]];
          is_computed_u[newind] = 1;
          list[over] = newind;
          is_row[over] = 1;
          over++;
        }
      }      
    }
    curr++;
  }

  
  bestred = 0.0;
  /* we break the count-loop and return to caller
     if the THRESHOLD is understepped in any row;
     state->next_row contains the index of the next row not considered */
  for (count = 0; count < m; count++) {
    i = state->next_row;
    for (l = 0; l < state->channels_byrow_over[i]; l++) {
      j = state->channels_byrow[i][l];
      /* we don't check that (i,j) is not 
         already in basis; theoretically the redcost is 0 in that case anyway;
         practically it should be zero also because the dual costs
         come exactly from taking diffs at all basis vectors
         remains the question: is computing n+m-1 times extra the redcost
         and comparing to bestred really more efficient than checking n*m times
         if (i,j) is 1 in state->basis */
      /* In the sparse version I changed this as we are checking only
         much fewer than n*m times if (i,j) is 1 in state->basis */
      if (MAT(basis,i,j) == 0) {	
        redcost = MAT(costm,i,j) - u[i] - v[j];
        if (redcost < bestred) {
          bestred = redcost;
          state->indi = i;
          state->indj = j;
        }
        if (bestred < THRESHOLD) {
          state->over = over;  /* just for diagnostics ??! */
          return(1);
        }
      }
    }
    state->next_row++;
    if (state->next_row == m) {
      state->next_row = 0;
    }
  }
    
  state->over = over;
  return(0);
}

