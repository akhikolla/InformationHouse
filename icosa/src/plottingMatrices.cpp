#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix edgeMatTri_(NumericMatrix v, NumericMatrix e) {
	int oldRows = e.nrow();
	int newRows = oldRows*2;

	NumericMatrix endObj(newRows, 3);
	int firstPoint;
	int secondPoint;

	for(int i=0; i<oldRows; i++){
		firstPoint = e(i,0);
		secondPoint = e(i,1);

		endObj(i*2,_) = v(firstPoint,_);
		endObj(i*2+1,_) = v(secondPoint,_);
	}

	return endObj;
}


// [[Rcpp::export]]
NumericMatrix triMatTri_(NumericMatrix v, NumericMatrix f) {
	int oldRows = f.nrow();
	int newRows = oldRows*3;

	NumericMatrix endObj(newRows, 3);
	int firstPoint;
	int secondPoint;
	int thirdPoint;

	for(int i=0; i<oldRows; i++){
		firstPoint = f(i,0);
		secondPoint = f(i,1);
		thirdPoint = f(i,2);

		endObj(i*3,_) = v(firstPoint,_);
		endObj(i*3+1,_) = v(secondPoint,_);
		endObj(i*3+2,_) = v(thirdPoint,_);
	}

	return endObj;
}

// [[Rcpp::export]]
NumericVector pointLayerColorOrder_(NumericMatrix f){
	int rows= f.nrow();
	
	NumericVector endObj(rows*3);
	
	for(int i=0; i<rows;i++){
		endObj(i*3) = f(i,0);
		endObj(i*3+1) = f(i,1);
		endObj(i*3+2) = f(i,2);
		
	}
	
	return endObj;
}
