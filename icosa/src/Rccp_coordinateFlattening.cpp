#include <Rcpp.h>
using namespace Rcpp;

// this is a small utility function to rearrange the vertices
// so that they form tetrahedra for every  face,
// in an xxxxyyyyzzzz format. In the final next function only the pointer moves.
// assumes the origin is (0,0,0)

// [[Rcpp::export]]
NumericVector xxxxyyyyzzzz_(NumericMatrix v, NumericMatrix f) {
	//the number of faces in a
	int nFaces = f.nrow();

	// the resulting vector with thhe coordinates
	NumericVector endObj(12*nFaces);

	//track the position
	int counter =0;

	NumericVector indexVec(3);

	NumericVector tempVec(3);

	//loop for the faces
	for(int i=0; i<nFaces;i++){

		indexVec =f(i,_);

		//loop for the coordinates
		for(int j=0; j<3;j++){

			//do a small loop to get the vector
			for(int d=0; d<3;d++){
			  int inVal= indexVec(d);
			  tempVec(d)= v(inVal,j);
			}

			//loop for the coordinates
			for(int l=0;l<4;l++){
				if(l<3){
					endObj(counter) = tempVec(l);
				}else{
					endObj(counter) = 0;
				}

				counter++;
			}

		}
	}

  return endObj;
}

// [[Rcpp::export]]
NumericVector xyz1xyz1xyz1xyz1_(NumericMatrix v, NumericMatrix f) {
	//the number of faces in a
	int nFaces = f.nrow();

	// the resulting vector with thhe coordinates
	NumericVector endObj(16*nFaces);

	//track the position
	int counter =0;

	NumericVector indexVec(3);

	int inVal=0;

	//loop for the faces - i
	for(int i=0; i<nFaces;i++){

		indexVec =f(i,_);

		//loop for the indices and origin
		for(int j=0; j<4;j++){
			// the indices of vectors
			if(j<3){
				inVal= indexVec(j);

				for(int k=0;k<4; k++){
					// for xyz
					if(k<3){
						endObj(counter) = v(inVal, k);
					// for 1
					}else{
						endObj(counter) = 1;
					}
					counter++;
				}
			// for the origin
			}else{
				for(int d=0;d<4;d++){
					if(d<3){
						endObj(counter) =0;
					}else{
						endObj(counter) =1;
					}
					counter++;
				}


			}
		}
	}

  return endObj;
}

// [[Rcpp::export]]
NumericVector xyz1 (NumericMatrix q){
	int rows = q.nrow();

	int newVecSize = rows*4;

	NumericVector newVec(newVecSize);

	for(int i=0; i<(rows);i++){
		newVec(i*4) = q(i,0);
		newVec(i*4+1) = q(i,1);
		newVec(i*4+2) = q(i,2);
		newVec(i*4+3) = 1;
	}

	return newVec;

}
