/*
   WARNINGS: It seems the degenerate case is not catered for!?
   In view of Luenberger, p.50 I think we should

   A small advantage is to be expected from taking distances as ints,
   by Luenberger, u, v are ints as well
*/


/*

   shortsimplex.c

   $Revision: 0.4.1 $   $Date: 2015/10/21 20:56:00 $

   Ported by Dominic Schuhmacher
   based on java code for class Shortlist2014 by Carsten Gottschlich
   (CG: in the comments refers to correspondences with original code,
   mainly used to give the original variable names where I have changed them) 

*/


/*  Call by

   .C("shortsimplex", mm = as.integer(m), nn = as.integer(n), a = as.integer(a), b = as.integer(b),
       costm = as.double(costm), assignment = as.integer(assignment), basis = as.integer(basis),
       DUP=TRUE) 
   # DUP kann man später auf falsch setzen um Zeit zu sparen SOLLTE MAN ABER NUR IN NOTFAELLEN!!! 
   # in package mit zusaetzlicher Option PACKAGE = "blabla"   

   a and b must NOT contain zeroes! (was ok in revsimplex I think)

*/

/* About using long instead of int below (from R-ext, p.17):
Do be very careful with passing arguments between R, C and FORTRAN code. In particular, long in C will be 32-bit on most R platforms (including those mostly used by the CRAN maintainers), but 64-bit on many modern Unix and Linux platforms. It is rather unlikely that the use of long in C code has been thought through: if you need a longer type than int you should use a configure test for a C99 type such as int_fast64_t (and failing that, long long) and typedef your own type to be long or long long, or use another suitable type (such as size_t). Note that integer in FORTRAN corresponds to int in C on all R platforms.

On my system int is 32 bit (I guess on every system), long is 64 bit (and longlong also).
 */

#include <math.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <stdlib.h>
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
/* FOR SHORTLIST STUFF: in vectors needed for initialization phases 1 and 2 
   memory is reserved by Calloc and Freed right afterwards */
typedef struct State {
  int shl_s, shl_s1, shl_k;  /* the shortlist parameters s and k, and s1 = s-1 */ 
  int shl_nabs_p; /* maximal number of searched lists, absolute value 
                      instead of relative p. CG: maxSearchedList */
  int *shl_byrow; /* the shortlist */

  int m, n;   /* lengths of mass vectors (CG: p, c) */
  int *a, *b;   /* mass vectors (CG: producer, CG: consumer) */ 
  double *costm;    /* m x n matrix of costs */
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
/* Note: the other dimension of the matrix does not matter
   ---> apply to shl_byrow */
#define MAT(NAME, I, J) (state->NAME)[state->m * (J) + (I)]
#define TMAT(NAME, I, J) (state->NAME)[state->n * (J) + (I)]
/* letzteres für basis_bycol, die n*m ist */
#define MAX(A,B) ((A)>(B) ? (A) : (B))
#define MIN(A,B) ((A)<(B) ? (A) : (B))

void init_shortlist(State *state);        /*  (CG: initShortListWithQuicksort) */
void partial_qsort(double *a, int *ind, int firsti, int lasti, int s1);
void init_assignment(State *state);       /*  (CG:initAssignment...)           */
void init_basis(State *state);            /*  (CG:findBasisForAssignment...)   */
void find_first_unconnected(State *state, int *firsti, int *firstj);
void label_connected(State *state, int firsti, int firstj);
void shl_init_helpers(State *state);          /*  (CG: initGraph, and more)        */
int shl_update_transport_rowmostneg(State *state);
int update_transport_shortlist(State *state);
int shl_new_basic_variable_rowmostneg(State *state);
int new_basic_variable_shortlist(State *state);
void shl_find_circle(State *state);
void shl_move_mass(State *state);
void shl_add_to_basis(State *state);
void shl_remove_from_basis(State *state);

/* double compute_total_costs();
int get_iterations();
int count_flows();
int count_basis_entries(); */

/* for debugging */
void shl_printit(State *state);
/* void pricediag(State *state); */
/* void circlediag(State *state); */
void shl_printfvec(int n, double *a);
void shl_printvec(int n, int *a);
void shl_printmat(int m, int n, int *a);


/* ------------ The main function ----------------------------- */

/* rename into revsimplex */
void shortsimplex(ss, kk, pp, mm, nn, a, b, costm, assignment, basis)
  /* inputs */
  int *ss, *kk;
  double *pp;
  int *mm, *nn;
  int *a, *b;
  double *costm;
  /* outputs */
  int *assignment;
  int *basis;
{
  int m,n,maxdim;
  State state;
  int is_optimal;

  /* inputs/outputs */
  state.shl_s = *ss;
  state.shl_s1 = *ss - 1;
  state.shl_k = *kk;   /* careful: never use k alone (is index in various places) */
  state.m = m = *mm;
  state.n = n = *nn;
  state.shl_nabs_p = trunc(*pp * m);
  state.shl_nabs_p = MAX(1, state.shl_nabs_p);
  state.a = a;
  state.b = b;
  state.costm = costm;
  
  state.assignment = assignment;
  state.basis = basis;

  state.maxdim = maxdim = MAX(m,n);
  state.iter = 0;
  state.next_row = 0;

  /* scratch space */
  /* I doubt that casting to long does something (positive) */
  state.shl_byrow = (int *) R_alloc((long) (m * state.shl_s), sizeof(int));
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

/* --- Phase 1 ----------------------- */
  init_shortlist(&state);
  /* shl_printmat(state.m, state.shl_s, state.shl_byrow);   Leave for diagnostics */
/* --- Phase 2 ----------------------- */
  init_assignment(&state);
  /* shl_printmat(state.m, state.n, state.assignment);    Leave for diagnostics */
  init_basis(&state);
  /* shl_printmat(state.m, state.n, state.basis);     Leave for diagnostics */
/* --- Phase 3 ----------------------- */
  shl_init_helpers(&state);
  is_optimal = 0;
  /* Maybe add a counter */
  while (is_optimal == 0) {
    R_CheckUserInterrupt();
    state.iter++;
    is_optimal = update_transport_shortlist(&state); 
                        /* of course also by row. most neg. within shortlist entries */
  }

  /* shl_printit(&state);     Leave for diagnostics */
  /* error("test interrupt"); */
/* --- Phase 4 ----------------------- */
  is_optimal = 0;
  /* Maybe add a counter */
  while (is_optimal == 0) {
    R_CheckUserInterrupt();
    state.iter++;
    is_optimal = shl_update_transport_rowmostneg(&state);
  }
/* ----------------------------------- */

  /* Rprintf("\n");   deleted 0.6 only  */
  /* line break for output in move_mass */
}


/*    */
int shl_update_transport_rowmostneg(State *state)
{
  int found;

  found = shl_new_basic_variable_rowmostneg(state);
  /* implicitly returns state->indi and state->indj, i.e.
     i and j index of new basic variable */

  if (found == 0) {
    return(1);
  }

  shl_add_to_basis(state);

  shl_find_circle(state);

  shl_move_mass(state);

  shl_remove_from_basis(state);

  /* shl_printit(state); */
  /*  error("test interrupt");  */

  return(0);  
}


int update_transport_shortlist(State *state)
{
  int found;

  found = new_basic_variable_shortlist(state);
  /* implicitly returns state->indi and state->indj, i.e.
     i and j index of new basic variable */

  if (found == 0) {
    return(1);
  }

  shl_add_to_basis(state);

  shl_find_circle(state);

  shl_move_mass(state);

  shl_remove_from_basis(state);

  /* shl_printit(state); */
  /*  error("test interrupt");  */

  return(0);  
}


/* Do partial quicksort for each row */
void init_shortlist(State *state)
{
  int i, j, m, n;
  double *costfromrow;    /* +/- CG: cost */
  int *indfromrow;    /* +/- CG: index */
                      /* i-th row contains the first shl_s entries of what we obtain in R by
                         order(costm[i,]) - 1   (minus 1 since indices in R start at 0) */
  m = state->m;
  n = state->n;

  costfromrow = (double *) Calloc((long) n, double);
  indfromrow = (int *) Calloc((long) n, int);
  
  for (i = 0; i < m; i++) {

    for (j = 0; j < n; j++) {
      indfromrow[j] = j;
      costfromrow[j] = MAT(costm,i,j);
    }
    /* copying two length n vectors here costs more
       than copying the two length shl_s vectors in Carsten's implementation */

    partial_qsort(costfromrow, indfromrow, 0, n-1, state->shl_s1);

    for (j = 0; j < state->shl_s; j++) {
      MAT(shl_byrow,i,j) = indfromrow[j];
    }

  }

  Free(costfromrow);
  Free(indfromrow);
}


void partial_qsort(double *a, int *ind, int firsti, int lasti, int s1)
{
  int piv, i, j, midi, tindex;  /* t as in "temporary" */
  double pivcost, temp;
  
  if (firsti < lasti) {  /* two or more elements to sort */
    i = firsti;
    j = lasti-1;

    /* choose pivot as "median of three" */
    midi = firsti + trunc((lasti - firsti)/2);
    if (a[midi] < a[firsti] && a[firsti] < a[lasti]) {
      piv = firsti;
    }
    else if (a[midi] < a[lasti] && a[lasti] < a[firsti]) {
      piv = lasti;
    } 
    else {
      piv = midi;
    }     

    /* Rprintf("firsti, piv, lasti: %d %d %d \n", firsti, midi, lasti); */
    pivcost = a[piv]; 

    /* move <=elements to the left of pivot and >elements to the right of pivot
       and update ind accordingly 
      -------------------------------------------*/
    /* swap pivot element with last element for getting it out of the way
       not entirely sure if this is the most efficient form of in-place sorting: */
    a[piv] = a[lasti];
    a[lasti] = pivcost;
    tindex = ind[piv];
    ind[piv] = ind[lasti];
    ind[lasti] = tindex; 

    for ( ; ; ) {  
      while (a[i] <= pivcost && i < lasti) {
        i++;
      }
      while (j >= 0 && a[j] > pivcost) {
        /* j can get to -1 and this is by design, but then a[j] will 
           not be evaluated due to the short circuiting of && */
        j--;
      }      
      if (i < j) {
        /* swap a[i] with a[j] and update indices accordingly */
        temp = a[i];
        a[i] = a[j];
        a[j] = temp;
        tindex = ind[i];
        ind[i] = ind[j];
        ind[j] = tindex;
      } else break;
    }

    /* swap pivot element into right place */
    a[lasti] = a[j+1];
    a[j+1] = pivcost;
    tindex = ind[lasti];
    ind[lasti] = ind[j+1];
    ind[j+1] = tindex;
    /*-------------------------------------------*/
 
    /* recursion; the if makes it partial... */
    partial_qsort(a, ind, firsti, j, s1);
    if (j+1 < s1) {
      partial_qsort(a, ind, j+2, lasti, s1);
    } 
  } 
}


/* return position (1,2 or 3)
   int medianofthree(double a, double b, double c) */


/* Initializes assignment and basis with shortlist
   anything passed from R to C is overwritten;
   CG:initAssignment... */
void init_assignment(State *state)
{
  int i, j, k;
  int assigmass, totmass, mass;
  int indexJ;  /* the selected colindex where we transport mass to */
  double mini;   /* 2015/10/21: in fact the large homogeneous examples get quite a bit slower
                    by this, but this leads to proper shortlist init assignment, whereas before
                    the init assignment was really random due to interpreting temporary minima
                    in rows of costm as ints. */

  int m,n;  
  int *aleft, *bleft;   /* row/col masses left (CG: prod, cons) */
  int *adone, *bdone;      /* 0-1, which rows/cols are exhausted
                              (CG: producerDone, consumerDone) */

  m = state->m;
  n = state->n;

  aleft = (int *) Calloc((long) m, int);
  bleft = (int *) Calloc((long) n, int);
  adone = (int *) Calloc((long) m, int);
  bdone = (int *) Calloc((long) n, int);

  for (i = 0; i < m; i++) {
  for (j = 0; j < n; j++) {
    MAT(assignment,i,j) = 0;
  }
  }

  assigmass = 0;
  totmass = 0;
  for (i = 0; i < state->m; i++) {
    totmass += state->a[i];
  }

  for (i = 0; i < m; i++) {
    aleft[i] = state->a[i];
    adone[i] = 0;
  }
  for (j = 0; j < n; j++) {
    bleft[j] = state->b[j];
    bdone[j] = 0;
  } 

  while (assigmass < totmass) {
     
    for (i = 0; i < m; i++) {
      if (adone[i] == 0) {

        indexJ = -1;
        for (k = 0; indexJ == -1 && k < state->shl_s; k++) {
          if (bdone[MAT(shl_byrow,i,k)] == 0) {
            indexJ = MAT(shl_byrow,i,k);
          }
        }

        if (indexJ == -1) {  /* we didn't get rid of all the rowmass in row i 
                                by assigning to cols from the shortlist */
          mini = R_PosInf;  
	  
          for (j = 0; j < n; j++) {
            if (bdone[j] == 0) {
              if (MAT(costm,i,j) < mini) {
                mini = MAT(costm,i,j);
                indexJ = j;
              }
            }
          }
        }

        mass = MIN(aleft[i],bleft[indexJ]); 
        MAT(assignment,i,indexJ) += mass;
        assigmass += mass;
        aleft[i] -= mass;
        bleft[indexJ] -= mass;
        if (aleft[i] == 0) {
          adone[i] = 1;
        }
        if (bleft[indexJ] == 0) {
          bdone[indexJ] = 1;
        }
        
      }
    }

  }
  
  Free(aleft);
  Free(bleft);
  Free(adone);
  Free(bdone);
  /*  MAT(basis,i,j) = */
  /* FIXBASIS (DE-DEGENERATE) */
}


/* Initializes basis based on assignment and fixes degeneracy
   anything passed from R to C is overwritten;
   CG:findBasisForAssignment... */
/* This functions is for computing non-degenerate basis based
   on assignment NO MATTER where the assignment comes from 
   assuming it contains <= m+n-1 positive entries */
void init_basis(State *state)
{
  int i, j;
  int ncurr;   /* current number of basis entries */
  int ntarget; /* target number of basis entries = m+n-1 */
  int m,n;  
  int firsti, firstj, secondi, secondj;
  /* in the case of degeneracy: basis entries that currently need a connection */ 

  m = state->m;
  n = state->n;
  ntarget = m + n - 1;
  ncurr = 0;
  for (i = 0; i < m; i++) {
  for (j = 0; j < n; j++) {
    if (MAT(assignment,i,j) > 0) {
      MAT(basis,i,j) = 1;
      ncurr++;
    }
    else {
      MAT(basis,i,j) = 0;
    }
  }
  }

  if (ncurr > ntarget) {
    error("the computed initial 'basis' has too many entries");
  }

  if (ncurr == ntarget) {
    return;
  }

  /* instead of CG:isConnected is set connected basis entries to 2. 
     saves some memory and might even be a little bit more efficient, although
     we have to undo it in the end */
  Rprintf("Initial solution based on shortlist is degenerate. Adding %d basis vector(s)... ", ntarget - ncurr);

  /* Find any basis entry */
  find_first_unconnected(state, &firsti, &firstj);
  label_connected(state, firsti, firstj);

  while (ncurr < ntarget) {
    find_first_unconnected(state, &secondi, &secondj);
    /* connect the second tree to the first and mark all of its entries as connected */
    MAT(basis,firsti,secondj) = 2;
    ncurr++;
    label_connected(state, secondi, secondj);
  }

  /* remove the marks */
  for (i = 0; i < m; i++) {
  for (j = 0; j < n; j++) {
    if (MAT(basis,i,j) == 2) {
      MAT(basis,i,j) = 1;
    }
  }
  }

  Rprintf("done.\n");  
}


/* this function is mainly outsourced so that we can (double) break gracefully
   when the entry is found */
void find_first_unconnected(State *state, int *firsti, int *firstj)
{
  int i,j;
  for (i = 0; i < state->m; i++) {
  for (j = 0; j < state->n; j++) {
    if (MAT(basis,i,j) == 1) {
      *firsti = i;
      *firstj = j;
      return;
    }
  }
  }
  error("no unconnected basis entry found in 'find_first_unconnected'");
}


/* mark all basis entry that connect (via same row/same column relation) to
   the entry (firsti, firstj) as connected (=2).
   Instead of CG:list we recycle state->circlea, state-> circleb here
   (the name is misleading, of course there are no circles here) */
void label_connected(State *state, int firsti, int firstj)
{
  int i,j,i0,j0;
  int todo, curr;  /* CG: countTodo, countDone */
  int *circlea, *circleb;   

  todo = 1;
  curr = 0;
  circlea = state->circlea;
  circleb = state->circleb;
  circlea[0] = firsti;
  circleb[0] = firstj;  
  MAT(basis, circlea[0], circleb[0]) = 2;

  while (curr < todo) {
    
    i0 = circlea[curr];
    j0 = circleb[curr];
    curr++;

    for (i = 0; i < state->m; i++) {
      if (MAT(basis,i,j0) == 1) {
        circlea[todo] = i;
        circleb[todo] = j0;
        todo++;
        MAT(basis,i,j0) = 2;
      }
    }

    for (j = 0; j < state->n; j++) {
      if (MAT(basis,i0,j) == 1) {
        circlea[todo] = i0;
        circleb[todo] = j;
        todo++;
        MAT(basis,i0,j) = 2;
      }
    }

  }
}


/* Initializes the various helper variables in state */
void shl_init_helpers(State *state)
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



int shl_new_basic_variable_rowmostneg(State *state)
{
  int i,j,newind;
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
  int curr=0;       /* is NOT in state; up to which position has list been processed (CG: done) */
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

  bestred = 0;
  /* we break the count-loop and return to caller
     if the THRESHOLD is understep in any row;
     state->next_row contains the index of the next row not considered */
  for (count = 0; count < m; count++) {
    i = state->next_row;
    for (j = 0; j < n; j++) {
      redcost = MAT(costm,i,j) - u[i] - v[j];
      if (redcost < bestred) {
        bestred = redcost;
        state->indi = i;
        state->indj = j;
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



int new_basic_variable_shortlist(State *state)
{
  int i,j,k,newind;
  int countlists, countcands;
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

  bestred = 0;
  countcands = 0;

  for (countlists = 0; countlists < state->shl_nabs_p; countlists++) {
    i = state->next_row;
    for (k = 0; k < state->shl_s; k++) {
      j = MAT(shl_byrow,i,k);
      if (MAT(basis,i,j) == 0) {
        redcost = MAT(costm,i,j) - u[i] - v[j];
        if (redcost < 0) {
          countcands++;
          if (redcost < bestred) {
            bestred = redcost;
            state->indi = i;
            state->indj = j;
          }
        }
      }
    }
    state->next_row++;
    if (state->next_row == m) {
      state->next_row = 0;
    }

    if (countcands >= state->shl_k) {
      break;  /* the outer for-loop */
    }
    
  }
   
  state->over = over;  /* just for diagnostics ??! */
  if (bestred < THRESHOLD) {
    return(1);
  }
  else {
    return(0);  /* current transport plan is optimal */
  } 
}



void shl_find_circle(State *state)
{
  int i, j, k;
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

  k = 0;
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
        /* KOMMEN WIR HIER UEBERHAUPT VORBEI? WIR SCHAUEN JA IMMER
           EINS VORAUS, UM DIESEN FALL FRÜHZEITIG ABZUFANGEN! */
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
        if (candi == indi && curr > 3) {
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
      k++;
      /* Rprintf("%d, ", k); */

  }  /* end of while */
  
  /*  Rprintf("\n FINISHED! \n");
  circlediag(state);
  Rprintf("ncand: %d \n", ncand);
  Rprintf("curr_fork/next_branch/curr: %d / %d / %d \n", curr_fork, next_branch, curr); 
  Rprintf("do_rowscan: %d \n \n", do_rowscan);  */

  state->circ_over = curr;  /* is always <= m+n (the "all basis vectors in the circle" case) */
}


void shl_move_mass(State *state)
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


void shl_add_to_basis(State *state)
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


void shl_remove_from_basis(State *state)
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




void shl_printit(State *state)
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
  Rprintf("\n\n");

  Rprintf("SHORTLIST STUFF\n\n");
  Rprintf("shortlist parameters (s,k,nabs_p): %d %d %d \n\n", state->shl_s,state->shl_k,state->shl_nabs_p);
  Rprintf("shortlist:  \n");
  for (i = 0; i < m; i++) {
  for (j = 0; j < state->shl_s; j++) {
    Rprintf("%d ", MAT(shl_byrow,i,j));
  }
  Rprintf("\n");
  }
  Rprintf("\n");

  Rprintf("\n\n\n");
}

/*
void pricediag(State *state)
{
  int i,j;
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
*/

/*
void circlediag(State *state)
{
  int i,j;
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
*/

void shl_printfvec(int n, double *a)
{
  int i;

  Rprintf("\n");
  for (i = 0; i < n; i++) {
    Rprintf("%2.9lf ", a[i]);
  }
  Rprintf("\n");
}

void shl_printvec(int n, int *a)
{
  int i;

  Rprintf("\n");
  for (i = 0; i < n; i++) {
    Rprintf("%d ", a[i]);
  }
  Rprintf("\n");
}   

void shl_printmat(int m, int n, int *a)
{
  int i,j;

  for (i = 0; i < m; i++) {
  for (j = 0; j < n; j++) {
    Rprintf("%d ", a[m * j + i]);
  }
  Rprintf("\n");
  }
  Rprintf("\n");
}

