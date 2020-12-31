#include"sparsebasicfeasible.h"

#include <math.h>
#include <R.h>
#include <R_ext/Utils.h>

// the next two are no forward-declared in header
// void dedewithchannels(int m, int n, int no_of_ones, int *basis, int *channels_byrow_over, int **channels_byrow);
// void findblocks(int m, int n, int *basis, int *nblock, int *rowblock, int *colblock);   for
// void printvec(int n, int *a);

#define MAX(A,B) ((A)>(B) ? (A) : (B))
#define MIN(A,B) ((A)<(B) ? (A) : (B))
#define CMAT(NAME, I, J) (NAME)[m * (J) + (I)]
#define ALTMAT(NAME, I, J) (NAME)[nbl * (J) + (I)]

// int arrayimin(int *a, int n)
  

// vermutllich koennen wir nicht davon ausgehen, dass channels_byrow sortiert ist, oder ?
//
// Output: basis
void dedewithchannels(int m, int n, int no_of_ones, int *basis, int *channels_byrow_over, int **channels_byrow)
{
  int i,j,k,l,g,h;
  int alarm;
  int rowsum, colsum;
  //  int *oldbasis;

  int *nblock;
  int *rowblock;
  int *colblock;

  int nbl;
  int *gchannels; 
  int *gwhichrows, *gwhichcols;
  int *chainlist;
  int chainlist_over;
  int rowg, colg, newg;
 

  rowblock = (int *) Calloc((long) m, int);
  colblock = (int *) Calloc((long) n, int);
  //  oldbasis = (int *) Calloc((long) (m * n), int);
  
  //  for (i = 0; i < m; i++) {    
  //  for (j = 0; j < n; j++) {
  //    CMAT(oldbasis,i,j) =  CMAT(basis,i,j);
  //  }
  //  }
  
  // check if given basis is legal
  for (i = 0; i < m; i++) {
  for (j = 0; j < n; j++) {
    if (CMAT(basis,i,j) == 1) {
      alarm = 1;
      for (k = 0; k < channels_byrow_over[i]; k++) {
        if (channels_byrow[i][k] == j) {
	    alarm = 0;
	}
      }
      if (alarm == 1) {
        error("dedegeneration failed. Degenerated basis has entry at %d %d, channels allow no transport there", i,j); 
      }
    }
  }
  }

  // basis rows and cols of zeros get a 1 in first entry 
  for (i = 0; i < m; i++) {
    rowsum = 0;
    for (j = 0; j < n; j++) {
      rowsum += CMAT(basis,i,j);
    }
    if (rowsum == 0) {
      CMAT(basis,i,1) = 1; 
    }
  }

  for (j = 0; j < n; j++) {
    colsum = 0;
    for (i = 0; i < m; i++) {
      colsum += CMAT(basis,i,j);
    }
    if (colsum == 0) {
      CMAT(basis,1,j) = 1; 
    }
  }

  nbl = 0;
  nblock = &nbl;
  findblocks(m, n, basis, nblock, rowblock, colblock);
  /*  Rprintf("%d blocks found: \n", nbl);
      printvec(m, rowblock);
      printvec(n, colblock); */
  // if (nbl == 0) error("pointer failure");  
  if (no_of_ones + nbl != m + n) {
    error("dedegeneration failed. Basis appears to have cycles. (This analysis is based on no_of_ones passed and occurred after fixing zero rows and cols.)");
  } 
  if (nbl == 1) {
    return;
  }

  // abstracting from blocks to entries:
  gchannels = (int *) Calloc(nbl * nbl, int);   // between which blocks is transport possible
  gwhichrows = (int *) Calloc(nbl * nbl, int);  // if transport is possible between pair of blocks, from which source to which target
  gwhichcols = (int *) Calloc(nbl * nbl, int);  // from which source (...rows) to which target (...cols)
  chainlist = (int *) Calloc(nbl, int);         // in which order do we chain blocks to 0-block (first entry 0)

  for (g = 0; g < nbl; g++) {
    chainlist[g] = 0;
  for (h = 0; h < nbl; h++) {
    ALTMAT(gchannels,g,h) = 0;
    ALTMAT(gwhichrows,g,h) = 0;
    ALTMAT(gwhichcols,g,h) = 0;
  }
  }
  
  for (i = 0; i < m; i++) {
  for (k = 0; k < channels_byrow_over[i]; k++) {
    j = channels_byrow[i][k];
    g = rowblock[i];
    h = colblock[j];
    if (g != h) {
      ALTMAT(gchannels,g,h) = 1;
      ALTMAT(gwhichrows,g,h) = i;  // contrary to the R-function we take last (i,j)
      ALTMAT(gwhichcols,g,h) = j;  // (for R-function first (i,j)); never mind
    }
  }
  }

  // chain blocks to 0-block via ones in rows or cols of gchannel
  chainlist[0] = 0;
  chainlist_over = 1;

  /*  Rprintf("gchannels:\n");
  for (g = 0; g < nbl; g++) {
  for (h = 0; h < nbl; h++) {
    Rprintf("%d ", ALTMAT(gchannels,g,h));
  }
  Rprintf("\n");
  }
  Rprintf("\n"); */

  /* Rprintf("gwhichrows:\n");
  for (g = 0; g < nbl; g++) {
  for (h = 0; h < nbl; h++) {
    Rprintf("%d ", ALTMAT(gwhichrows,g,h));
  }
  Rprintf("\n");
  }
  Rprintf("\n"); */

  /* Rprintf("gwhichcols:\n");
  for (g = 0; g < nbl; g++) {
  for (h = 0; h < nbl; h++) {
    Rprintf("%d ", ALTMAT(gwhichcols,g,h));
  }
  Rprintf("\n");
  }
  Rprintf("\n"); */
  
  for (l = 0; l < nbl-1; l++) {
    newg = -1;
    for (k = 0; k < chainlist_over; k++) {
      for (h = 0; h < nbl; h++) {
        if (ALTMAT(gchannels,chainlist[k],h) == 1) {
          rowg = chainlist[k];
	  colg = h;
          newg = h;
          break;
        }
      }
      if (newg > -1) break;
    }
    if (newg == -1) {
      for (k = 0; k < chainlist_over; k++) {
        for (g = 0; g < nbl; g++) {
          if (ALTMAT(gchannels,g,chainlist[k]) == 1) {
            rowg = g;
            colg = chainlist[k];
            newg = g;
            break;
          }
        }
        if (newg > -1) break;
      }
    }
    if (newg == -1) {
      error("dedegeneration failed. Channels don't allow for an extended basis that is not degenerated");
    }

    // set the connecting entry in basis that is responsible for the 1 in at entry (rowg,colg) in gchannel to 1
    CMAT(basis,ALTMAT(gwhichrows,rowg,colg),ALTMAT(gwhichcols,rowg,colg)) = 1;
    chainlist[chainlist_over] = newg;
    chainlist_over++;
    for (k = 0; k < chainlist_over; k++) {
      ALTMAT(gchannels,chainlist[k],newg) = 0;
      ALTMAT(gchannels,newg,chainlist[k]) = 0;  // entry (newg,newg) set to zero twice ;-)
    }    
  }
  
  Free(gchannels);
  Free(gwhichrows);
  Free(gwhichcols);
  Free(chainlist);
  
  Free(rowblock);
  Free(colblock);
}




// we use a simpler data structure for the blocks compared to the R function
// lose the order of the spanning tree (note that we didn't have the spanning tree
// itself in the R function either)
//
// Checks just what is connected to what if we interpret basis as
// "bipartite adjacency matrix" (i.e. row and cols stand for different nodes)
void findblocks(int m, int n, int *basis, int *nblock, int *rowblock, int *colblock)
{
  int i,j,k;
  int searchforrows;
  int nr2go, nc2go, nrlist_over, nclist_over;
  int *rows2go, *cols2go;
  int *rowlist, *collist;
  
  rows2go = (int *) Calloc(m, int);
  rowlist = (int *) Calloc(m, int);
  cols2go = (int *) Calloc(n, int);
  collist = (int *) Calloc(n, int);
  
  for (i = 0; i < m; i++) {
    rowblock[i] = -1;   // different from R: which block the i-th row is assigned to
                        // -1 means not assigned
    rows2go[i] = 1;  // different from R: 1 if this row index is to go,
                     // 0 if we have dealt with it already
    rowlist[i] = -1;   // i-th element current list of row indices (-1 means NA)
  }

  for (j = 0; j < n; j++) {
    colblock[j] = -1;
    cols2go[j] = 1;
    collist[j] = -1;
  }
  
  *nblock = -1;  // no of current block
  nr2go = m;    // no of 1s in rows2go
  nc2go = n;
  while (nr2go > 0 || nc2go > 0) {
    searchforrows = 1;
    nrlist_over = 0;  // index one larger than the largest that matters
    nclist_over = 0;
    for (j = 0; j < n; j++) {
      if (cols2go[j] == 1) {
	collist[nclist_over] = j;  
	nclist_over++;
	cols2go[j] = 0;
	nc2go--;
	break;
      }      
    }
    
    // the following can only happen if there are rows without basis vectors
    // (and hence zero mass rows if the basis comes from a transport)
    // but was previously uncaught (especially in the R-functions)
    if (nclist_over == 0) {
      for (i = 0; i < m; i++) {
	if (rows2go[i] == 1) {
	  (*nblock)++;
	  rowblock[i] = *nblock; 
          rows2go[i] = 0;
	  nr2go--;
	}
      }
      break;  // outer while
    }

    (*nblock)++;
    colblock[collist[nclist_over-1]] = *nblock;

    while (nrlist_over > 0 || nclist_over > 0) {
      if (searchforrows == 1) {
        for (k = 0; k < nclist_over; k++) {
          for (i = 0; i < m; i++) {
            if (CMAT(basis,i,collist[k]) == 1 && rows2go[i] == 1) {
	      rowlist[nrlist_over] = i;
	      nrlist_over++;
              rowblock[i] = *nblock;
	      rows2go[i] = 0;
	      nr2go--;
	    }
	  }   
	}
        nclist_over = 0;
        searchforrows = 0;
      } else {
        for (k = 0; k < nrlist_over; k++) {
          for (j = 0; j < n; j++) {
            if (CMAT(basis,rowlist[k],j) == 1 && cols2go[j] == 1) {
	      collist[nclist_over] = j;
	      nclist_over++;
	      colblock[j] = *nblock;
	      cols2go[j] = 0;
	      nc2go--;
	    }
	  }
	}
	nrlist_over = 0;
	searchforrows = 1;
      }
    }   // end inner while
    
  }     // end outer while

  (*nblock)++;  // nblock was pointing on the last block no written, so now +1 
  Free(rows2go);
  Free(rowlist);
  Free(cols2go);
  Free(collist);
  //  return 0;
}



/* void printvec(int n, int *a)
{
  int i;

  Rprintf("\n");
  for (i = 0; i < n; i++) {
    Rprintf("%d ", a[i]);
  }
  Rprintf("\n");
  }   
*/
