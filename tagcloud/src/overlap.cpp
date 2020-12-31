#include "Rcpp.h"
#include <cstdlib>
#include <time.h>

using namespace Rcpp;

/* check whether row no. row_in overlaps with any other row in the boxes_in
 * matrix */
// [[Rcpp::export(name="row.overlap")]]
LogicalVector row_overlap( NumericVector row_in, NumericMatrix boxes ){

	int row = as<int>(row_in);
  double x1, y1, w1, h1 ;
  double x2, y2, w2, h2 ;

  row-- ;
	if( boxes.nrow() < 2 || boxes.nrow() < row ) return wrap(false);

  x1 = boxes( row, 0 ) ; y1 = boxes( row, 1 ) ; w1 = boxes( row, 2 ) ; h1 = boxes( row, 3 ) ;

  for( int i = 0 ; i < boxes.nrow() ; i++ ) {
    if( i == row ) continue ;

    x2 = boxes( i, 0 ) ; y2 = boxes( i, 1 ) ; w2 = boxes( i, 2 ) ; h2 = boxes( i, 3 ) ;

    if( ! ( ( x1 + w1 < x2 ) | ( x2 + w2 < x1 ) | 
            ( y1 + h1 < y2 ) | ( y2 + h2 < y1 ) ) ) {
      return( wrap( true ) ) ;
    }

  }

	return wrap(false);
}

/* report all overlaps within a matrix. Returns a square matrix with as
 * many rows and columns as there are rows in boxes_in */
RcppExport SEXP all_overlaps( SEXP boxes_in ){
	NumericMatrix boxes( boxes_in ) ;
	NumericMatrix ret( boxes.nrow(), boxes.nrow() ) ;
  double x1, y1, w1, h1 ;
  double x2, y2, w2, h2 ;
	bool overlap= true;


  for( int i = 0 ; i < boxes.nrow() ; i++ ) {
    x1 = boxes( i, 0 ) ; y1 = boxes( i, 1 ) ; w1 = boxes( i, 2 ) ; h1 = boxes( i, 3 ) ;

    for( int j = i  ; j < boxes.nrow() ; j++ ) {
      if( i == j ) {
        ret( i, j ) = 1 ;
        continue ;
      }
      x2 = boxes( j, 0 ) ; y2 = boxes( j, 1 ) ; w2 = boxes( j, 2 ) ; h2 = boxes( j, 3 ) ;

      overlap= false ;
      if( ! ( ( x1 + w1 < x2 ) | ( x2 + w2 < x1 ) | 
              ( y1 + h1 < y2 ) | ( y2 + h2 < y1 ) ) ) {
        overlap= true ;
      }

      if(overlap) {
        ret( i, j ) = 1 ;
        ret( j, i ) = 1 ;
      } else {
        ret( i, j ) = 0 ;
        ret( j, i ) = 0 ;
      }

    }
  }

	return wrap(ret);
}


/* reports whether any two rows in the matrix boxes_in overlap. Returns
 * true or false */
// [[Rcpp::export(name="any.overlap")]]
LogicalVector any_overlap( NumericMatrix boxes ){
  double x1, y1, w1, h1 ;
  double x2, y2, w2, h2 ;

  if( boxes.nrow() < 2 )
	return wrap(false);

  for( int i = 0 ; i < boxes.nrow() - 1 ; i++ ) {
    x1 = boxes( i, 0 ) ; y1 = boxes( i, 1 ) ; w1 = boxes( i, 2 ) ; h1 = boxes( i, 3 ) ;

    for( int j = i + 1 ; j < boxes.nrow() ; j++ ) {
      x2 = boxes( j, 0 ) ; y2 = boxes( j, 1 ) ; w2 = boxes( j, 2 ) ; h2 = boxes( j, 3 ) ;

      if( ! ( ( x1 + w1 < x2 ) | ( x2 + w2 < x1 ) | 
              ( y1 + h1 < y2 ) | ( y2 + h2 < y1 ) ) ) {
        return( wrap( true ) ) ;
      }
    }
  }

	return wrap(false);
}


/* checks whether the vector box_in overlaps with any of the rows in
 * boxes_in */
// [[Rcpp::export(name="is.overlap")]]
LogicalVector is_overlap( NumericVector box, NumericMatrix boxes ) {

  double x1=box[0], y1=box[1], w1=box[2], h1=box[3] ;
	NumericVector bnds;
	double x2, y2, w2, h2;

	for (int i=0;i < boxes.nrow();i++) {
    x2 = boxes(i,0) ;
    y2 = boxes(i,1) ;
    w2 = boxes(i,2) ;
    h2 = boxes(i,3) ;

    if( ! ( ( x1 + w1 < x2 ) | ( x2 + w2 < x1 ) | 
            ( y1 + h1 < y2 ) | ( y2 + h2 < y1 ) ) ) {
      return( wrap( true ) ) ;
    }

	}

	return wrap(false);
}

/* move in a spiral */
// [[Rcpp::export(name="run.spiral")]]
NumericVector spiral( List params,
                      NumericMatrix boxes ) {

  NumericVector ret( 3 ) ;

	bool overlap= true;
  double x2, y2, w2, h2 ;
  double x =  as<double>(params["x"]),
         y =  as<double>(params["y"]),
         w =  as<double>(params["w"]),
         h =  as<double>(params["h"]),
         r =  as<double>(params["r"]),
         angle  = as<double>(params["angle"]),
         astep  = as<double>(params["astep"]),
         rstep  = as<double>(params["rstep"]),
         asp    = as<double>(params["aspect"]),
         maxr   = as<double>(params["maxr"]) ;
  int    dir  = as<int>(params["dir"]),
         max_iter  = as<int>(params["max.iter"]) ;

  while( max_iter > 0 ) {

    /* test the current overlap. None? we're done */
    overlap = false ;
    for (int i=0;i < boxes.nrow();i++) {
      x2 = boxes(i,0) ;
      y2 = boxes(i,1) ;
      w2 = boxes(i,2) ;
      h2 = boxes(i,3) ;

      if( ! ( ( x + w < x2 ) | ( x2 + w2 < x ) | 
              ( y + h < y2 ) | ( y2 + h2 < y ) ) ) {
        /* there is an overlap */
        overlap = true ;
        break ;
      } 

    }

    if( ! overlap ) {
      ret[0] = x ;
      ret[1] = y ;
      ret[2] = 1 ;
      return( ret ) ;
    }

    angle += dir * astep ;
    r += rstep ;
    x = r * cos( angle ) / 2 ;
    y = r * sin( angle ) / 2 * asp ;

    if( r > 2 * maxr ) {
      return( ret ) ;
    }

  }
  
  return( ret ) ;
}


/* move in an ulam spiral */
// [[Rcpp::export(name="run.ulam")]]
NumericVector ulam( List params,
                      NumericMatrix boxes ) {

  NumericVector ret( 3 ) ;

	bool overlap= true;
  double x2, y2, w2, h2 ;
  double x =  as<double>(params["x"]),
         y =  as<double>(params["y"]),
         w =  as<double>(params["w"]),
         h =  as<double>(params["h"]),
         rstep = as<double>(params["rstep"]),
         asp   = as<double>(params["aspect"]),
         maxr  = as<double>(params["maxr"]) ;
  double dr = 0, r = rstep ;
  int    dir1  = as<int>(params["dir1"]),
         dir2  = as<int>(params["dir2"]),
         max_iter  = as<int>(params["max.iter"]),
         tmp ;

  NumericVector foo( 1 ) ;
  //srand( ( unsigned ) time( NULL ) ); 

  while( max_iter > 0 ) {

    /* test the current overlap. None? we're done */
    overlap = false ;
    for (int i=0;i < boxes.nrow();i++) {
      x2 = boxes(i,0) ;
      y2 = boxes(i,1) ;
      w2 = boxes(i,2) ;
      h2 = boxes(i,3) ;

      if( ! ( ( x + w < x2 ) | ( x2 + w2 < x ) | 
              ( y + h < y2 ) | ( y2 + h2 < y ) ) ) {
        overlap = true ;
        break ;
      } 

    }

    if( ! overlap ) {
      ret[0] = x ;
      ret[1] = y ;
      ret[2] = 1 ;
      return( ret ) ;
    }

    dr += rstep ;
    x += dir1 * rstep ;
    y += dir2 * rstep * asp ;

    if( dr > r ) {
      foo = runif( 1, 0, 1 ) ;
      dr = rstep * 0.5 * (*foo.begin()) ;
      r += rstep ;
      tmp = dir1 ;
      dir1 = -dir2 ;
      dir2 = tmp ;
    }

    if( r > 2 * maxr ) {
      return( ret ) ;
    }

  }
  

  return( ret ) ;
}
