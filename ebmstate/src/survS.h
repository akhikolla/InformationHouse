/* [From the package 'survival' (v. 2.44-1.1) by Terry Therneau]
 */

/* 
**  This file started out with support for Splus, then morphed to allow
**  either R or Splus (based on ifdef lines), and now is R only.
*/
#include "R.h"
#include "Rinternals.h"
#include <R_ext/Utils.h>  

/* typedef int Sint; */  /* no longer needed */

/*
** Memory defined with ALLOC is removed automatically by S.
**  That with "Calloc" I have to remove myself.  Use the
**  latter for objects that need to to persist between calls.
*/
#define ALLOC(a,b)  R_alloc(a,b)

