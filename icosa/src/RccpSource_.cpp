#include <Rcpp.h>
#include <math.h>
#include <cstdlib>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Collapse_(NumericVector vect){
	int orig = vect.size();
	NumericVector tempVect(orig);

	//initialize
	tempVect(0) = vect(0);
	int counter=1;
	for(int i=1; i<orig;i++){
		if(vect(i)!=vect(i-1)){
			tempVect(counter) = vect(i);
			counter++;
		}
	}


	// copy
	NumericVector newVect(counter);
	for(int i=0; i<counter; i++){
		newVect(i) = tempVect(i);
	}


	return newVect;

}



double dist(NumericVector point1, NumericVector point2){
	double endObj;
	endObj = sqrt(pow(point1(0)-point2(0),2.0)+pow(point1(1)-point2(1),2.0)+pow(point1(2)-point2(2),2.0));
	return endObj;
}

// [[Rcpp::export]]
double ArcDist_(NumericVector coord1, NumericVector coord2, NumericVector origin, bool method) {
	double endObj = 0;

	double bLength = dist(coord1, origin);
	double cLength = dist(coord2, origin);
	double aLength = dist(coord2, coord1);

	double alpha = 0;

	double inside = (pow(bLength,2.0)+pow(cLength,2.0)-pow(aLength,2.0))/(2*bLength*cLength);
	alpha = acos(inside);

	// disnance
	if(method){
		endObj = alpha*bLength;
	}else{
		endObj = alpha;
	};

	return endObj;
}


// [[Rcpp::export]]
NumericMatrix SymmetricArcDistMat_(NumericMatrix xyzMat, NumericVector origin, bool method) {

	int entries = xyzMat.nrow();

	NumericMatrix endObj(entries,entries);

	for(int i=0; i<entries; i++){
		for(int j=(i+1); j<entries;j++){
		  double value = ArcDist_(xyzMat(i,_),xyzMat(j,_), origin, method);
			endObj(i,j) = value;
			endObj(j,i) = endObj(i,j);
		}

	}

	return endObj;
}


// [[Rcpp::export]]
NumericMatrix ArcDistMat_(NumericMatrix gridPoints,NumericMatrix queries, NumericVector origin, bool method) {

	int rows = gridPoints.nrow();
	int cols = queries.nrow();

	NumericMatrix endObj(rows,cols);

	for(int i=0; i<rows; i++){
		for(int j=0; j<cols;j++){
			double value = ArcDist_(gridPoints(i,_),queries(j,_), origin, method);
			endObj(i,j) = value;
		}

	}
	return endObj;
}

// [[Rcpp::export]]
NumericVector ArcDistMany_(NumericMatrix p0, NumericMatrix p1, NumericVector origin, bool method){
	int rows=p0.nrow();

	NumericVector endObj(rows);

	for(int i=0;i<rows;i++){
		endObj(i) = ArcDist_(p0(i,_),p1(i,_), origin, method);
	}

	return endObj;

}

// [[Rcpp::export]]
int whichMinVector_(NumericVector one) {
	int length = one.size();

	double minVal = one(0);
	int endIndex =0;
	double tempVal;

	for(int i=1; i<length; i++){
		tempVal = one(i);
		if(tempVal<=minVal){
			endIndex = i;
			minVal = tempVal;
		}

	}
	return endIndex;
}

// [[Rcpp::export]]
NumericVector Aggregate_(NumericMatrix gridPoints,NumericMatrix queries, NumericVector origin){

	int rows = gridPoints.nrow();
	int cols = queries.nrow();

	NumericVector endObj(cols);

	bool method = 1;

	double minVal=0;

	// for every query
	for(int j=0; j<cols; j++){

		// a proximity check is necessary
		for(int i=0; i<rows;i++){
			double value = ArcDist_(gridPoints(i,_),queries(j,_), origin, method);
			if(i==0){
				minVal = value;
				endObj(j) = i;
			}
			if(value<minVal){

				endObj(j) = i;
				minVal = value;
			}

		}


	}

	return endObj;

}

// [[Rcpp::export]]
NumericVector edges_(NumericMatrix v, NumericMatrix e, NumericVector origin, bool method){
	// the number of edges
	int nEdge = e.nrow();
	//vector containing the lengths
	NumericVector endObj(nEdge);
	// temporary variable storing the indices
	int *tempIndex = new int[2];

	for(int i=0; i<nEdge; i++){
		*tempIndex = e(i,0);
		*(tempIndex+1) = e(i,1);

		endObj(i) = ArcDist_(v(*tempIndex,_),v(*(tempIndex+1),_), origin, method);
	}

	delete[] tempIndex;
	return endObj;

}

// calculate the dihedral angle from dihedral angle opposite of 'a'
double DihedralAngle(double a, double b, double c){
	double temp = acos((cos(a)-(cos(b)*cos(c)))/(sin(b)*sin(c)));
	return temp;

}

 // the surface of one spherical triangle
// [[Rcpp::export]]
 double SphericalTriangleSurface_(NumericVector coord1, NumericVector coord2, NumericVector coord3, NumericVector origin, double pi){

	double radius = dist(origin, coord1);

 	double ang12 = ArcDist_(coord1, coord2, origin,0);
 	double ang13 = ArcDist_(coord1, coord3, origin,0);
 	double ang23 = ArcDist_(coord2, coord3, origin,0);

 	// dihedral angle in opposite of a
 	double d12 = DihedralAngle(ang12,ang13,ang23);
 	double d13 = DihedralAngle(ang13,ang12,ang23);
 	double d23 = DihedralAngle(ang23,ang12,ang13);

 	double surface = pow(radius,2.0)*((d12+d13+d23)-pi);

 	return surface;


 }

 // [[Rcpp::export]]
 NumericVector spherTriSurfs(NumericMatrix v, NumericMatrix f, NumericVector origin, double pi){
 	// the number of faces in the grid(highest res.)
	int nFaces = f.nrow();

	//  the final container vecctor
 	NumericVector endObj(nFaces);

	// repeat surface calculation for every triangle
	for(int i=0;i<nFaces;i++){
		endObj(i) = SphericalTriangleSurface_(
			v(f(i,0),_),
			v(f(i,1),_),
			v(f(i,2),_),
			origin,
			pi
		);

	}

	return endObj;
 }


// [[Rcpp::export]]
double surfConvHullTri(NumericMatrix v, NumericVector cent, NumericVector origin, double pi){
	// warning! use only if the great circle archs do not intersect, otherwise you get area inflation!

	// the number of triangles 
	int nTri = v.nrow();

	//  the total area covered by the shape
	double endSum=0;

	// repeat surface calculation for every triangle
	for(int i=1;i<nTri;i++){
		endSum += SphericalTriangleSurface_(
			v(i-1,_),
			v(i,_),
			cent,
			origin,
			pi
		);
	}

	// and then add up the last triangle
	endSum += SphericalTriangleSurface_(
		v(nTri-1,_),
		v(0,_),
		cent,
		origin,
		pi
	);

	return endSum;
}


// [[Rcpp::export]]
NumericMatrix SplitArc_(NumericVector coord1, NumericVector coord2, NumericVector center, int breaks, bool onlyNew){
  NumericMatrix tempObj(breaks,3);

	//define the origin
	NumericVector origin(3);
	origin(0) = 0;
	origin(1) = 0;
	origin(2) = 0;

	// transform the coordinates to the origin
	NumericVector tCoord1(3);
	NumericVector tCoord2(3);
	for(int i=0;i<3;i++){
		tCoord1(i)=coord1(i)-center(i);
		tCoord2(i)=coord2(i)-center(i);
	}

	// given that there are no breaks
	if(breaks!=0){
		//a temporary object

		double bLength = dist(tCoord1, origin);
		double cLength = dist(tCoord2, origin);
		double aLength = dist(tCoord1, tCoord2);

		double alpha = acos((pow(bLength,2.0)+pow(cLength,2.0)-pow(aLength,2.0))/(2*bLength*cLength));

		// to be rewritten in each loop
		double theta;

		for(int i=0; i<breaks; i++){
			theta = alpha/(breaks+1)*(i+1);
			tempObj(i,0) = center(0)+(sin(alpha-theta)*tCoord1(0)+sin(theta)*tCoord2(0))/sin(alpha);
			tempObj(i,1) = center(1)+(sin(alpha-theta)*tCoord1(1)+sin(theta)*tCoord2(1))/sin(alpha);
			tempObj(i,2) = center(2)+(sin(alpha-theta)*tCoord1(2)+sin(theta)*tCoord2(2))/sin(alpha);

		}
	}

	// do you want to return the first and the last values or only the new coords?
	if(onlyNew){

		return tempObj;

	}else{
		//declare
		NumericMatrix endObj(breaks+2, 3);

		if(breaks==0){
			endObj(0,_) = coord1;
			endObj(1,_) = coord2;

		}else{

			//the first point
			endObj(0,0) = coord1(0);
			endObj(0,1) = coord1(1);
			endObj(0,2) = coord1(2);

			//copy in the New matrix
			for(int i=0; i<breaks;i++){
				endObj(i+1,0) = tempObj(i,0);
				endObj(i+1,1) = tempObj(i,1);
				endObj(i+1,2) = tempObj(i,2);
			}

			int rows = endObj.nrow();
			//the last point
			endObj(rows-1,0) = coord2(0);
			endObj(rows-1,1) = coord2(1);
			endObj(rows-1,2) = coord2(2);
		}
		return endObj;
	}


}

// [[Rcpp::export]]
NumericMatrix GreatCircle_(NumericVector coord1, NumericVector coord2, NumericVector origin, int breaks, double pi){
  NumericMatrix tempObj(breaks,3);

	// given that there are no breaks
	if(breaks!=0){
		//a temporary object

		double bLength = dist(coord1, origin);
		double cLength = dist(coord2, origin);
		double aLength = dist(coord1, coord2);

		double alpha = acos((pow(bLength,2.0)+pow(cLength,2.0)-pow(aLength,2.0))/(2*bLength*cLength));

		// to be rewritten in each loop
		double theta;

		for(int i=0; i<breaks; i++){
			theta = alpha/(breaks+1)*(i+1);
			theta = 2*pi/(breaks+1)*(i+1);
			tempObj(i,0) = origin(0)+(sin(alpha-theta)*coord1(0)+sin(theta)*coord2(0))/sin(alpha);
			tempObj(i,1) = origin(1)+(sin(alpha-theta)*coord1(1)+sin(theta)*coord2(1))/sin(alpha);
			tempObj(i,2) = origin(2)+(sin(alpha-theta)*coord1(2)+sin(theta)*coord2(2))/sin(alpha);
		}
	}

	// do you want to return the first and the last values or only the new coords?

	//declare
	NumericMatrix endObj(breaks+1, 3);

	if(breaks==0){
		endObj(0,_) = coord1;


	}else{

		//the first point
		endObj(0,0) = coord1(0);
		endObj(0,1) = coord1(1);
		endObj(0,2) = coord1(2);

		//copy in the New matrix
		for(int i=0; i<breaks;i++){
			endObj(i+1,0) = tempObj(i,0);
			endObj(i+1,1) = tempObj(i,1);
			endObj(i+1,2) = tempObj(i,2);
		}


	}

	return endObj;


}

// [[Rcpp::export]]
NumericVector NeighbourOfOneFace_(NumericMatrix faces, int faceNo){

	int nF = faces.rows();
	NumericVector actFace = faces(faceNo,_);

	LogicalMatrix tempRes(nF,3);

	int rowRes;

	//No face can have more than 6 neighbours
	NumericVector semiEndObj(6);
	int counter=0;

	// do consequtive lookups
	//all three coordinates
	//every row in a column
	for(int i=0;i<nF;i++){
		//all vertex columns
		for(int j=0;j<3;j++){
			//all three coordinates
			for(int k=0;k<3;k++){
				if(actFace(k)==faces(i,j)){
					tempRes(i,j) = 1;

				}


			}

		}
		rowRes= tempRes(i,0)+tempRes(i,1)+tempRes(i,2);
		if(rowRes>1){
			semiEndObj(counter) = i;
			counter++;
		}

	}

	// clean up:
	NumericVector endObj(counter);
	for(int i=0; i<counter;i++){
		endObj(i) = semiEndObj(i);
	}

	return endObj;

}

// [[Rcpp::export]]
NumericMatrix DirectNeighboursTri_(NumericMatrix faces){
	int rows = faces.rows();

	NumericMatrix endObj (rows,4);
	NumericVector currentFace;

	int counter;

	for(int i=0; i<rows;i++){
		currentFace = NeighbourOfOneFace_(faces, i);
		counter=0;
		for(int j=0;j<currentFace.size();j++){
			endObj(i,counter)=currentFace(j);
			counter++;

		}

	}

	return endObj;

}


// function to extract all the neighbouring faces in hierarchies
// will accept the faces matrix in a 3 column format

// [[Rcpp::export]]
NumericMatrix AllNeighboursTri_(NumericMatrix allFaces, NumericVector div){
	// the number of divisions
	int nDiv = div.size();

	// the total number of end-res faces
	NumericVector faceInRes(nDiv);
	int nF=1;
	for(int i=0;i<nDiv;i++){
		nF=nF*div(i);
		faceInRes(i) = nF;
	}

	// offsets between the indices of different resolutions in the table
	NumericVector offset (nDiv+1);
	offset(0) = 0;
	int faceNo=1;

	for(int d=0;d<nDiv;d++){
		faceNo = div(d)*faceNo;
		offset(d+1)= offset(d)+faceNo;
	}

	//start with the icosahedron
	NumericMatrix prevResMat(nF,4); // set matrix size so that it no longer needs to be adjusted

	// get the starting neighbours (icosahedron)
	NumericMatrix icosaMat(20,3);
	for(int i=0;i<20;i++){
		icosaMat(i,_) = allFaces(i,_);
	}

	// the actual object containing the neighbouring faces
	NumericMatrix startNeighbours;
	startNeighbours = DirectNeighboursTri_(icosaMat);

	// add an exception for the icosahedron trigrid
	if(nDiv==1){
		return startNeighbours;
	}

	// put these results to the object that will be used by the loops
	for(int i=0;i<20;i++){
		prevResMat(i,_) = startNeighbours(i,_);
	}


	// vector container of the potential neighbour face indices
	NumericVector potentialNeighbours(nF);

	//the matrix containing the results of the div loop (including the final)
	NumericMatrix newResMat(nF,4);

	// changing variables
	int oldIndex;
	int divOffSet;
	int counter;
	int newInd=-99;
	int faceCounter;

	NumericVector facesToBeSearched(4);
//	NumericMatrix faceMat;
	NumericVector tempRes(4);

	// iterate for all divisions
	for(int i=1;i<nDiv;i++){
		faceCounter=0;
//		Rcout << "\n";
//		Rcout << "starting div loop: " << i << "\n";


		// count the faces to know how many faces we know the neighbours of

		divOffSet = offset(i);
//		Rcout << "divOffSet: " << divOffSet << "\n";

		//iterate for all faces in that division
		for(int j=0; j<faceInRes(i); j++){
//			Rcout << "\n";
//			Rcout << "starting face loop: " << j << "\n";

			//the active face belongs to a face in the previous divRound - integer division
			oldIndex = j/div(i);

			// the Neighbours of that face will be on the line oldIndex
			// all the coarse resolution faces that needs to be searched are these
			facesToBeSearched= prevResMat(oldIndex,_);
//			Rcout << "oldIndex: " << oldIndex << "\n";

//			Rcout << "facesToBeSearched: " << facesToBeSearched << "\n";

			// for all the 4 coarse faces
			counter=0;
			for(int k=0;k<4;k++){
				// every coarse face is divided to div(i) subfaces
				for (int l=0; l<div(i);l++){
					potentialNeighbours(counter) = facesToBeSearched(k)*div(i)+l;
					counter++;
				}
			}

//			Rcout << "potentialNeighbours: " << potentialNeighbours << "\n";

			// now the indices of the faces to be searched for are in potentialNeighbours
			// we must store these in a new object
			NumericMatrix faceMat(counter,3);
			// saturate this with the vertex indices
			for(int k=0;k<counter;k++){
				faceMat(k,_) = allFaces(potentialNeighbours(k)+divOffSet,_);

				// the index of the face we are interested in the current round has changed in faceMat
				if(j==potentialNeighbours(k)){
					newInd = k;
				}
			}

//			Rcout << "faceMat: " << faceMat << "\n";
//			Rcout << "newInd: " << newInd << "\n";
			// now search for the active face and store it
			tempRes = NeighbourOfOneFace_(faceMat,newInd);

//			Rcout << "tempRes: " << tempRes << "\n";

			// this vector contains indices that point to the faceMat table
			// reconstruct the original indices
			for(int k=0;k<4;k++){
				newResMat(j,k) =  potentialNeighbours(tempRes(k));
			}
			NumericVector tempRow = newResMat(j,_);

//			Rcout << "original indices: " << tempRow << "\n";
			faceCounter++;
		}

		// produce a neighbourmatrix that will be used later
		// never use NumeriMatrix1 = NumericMatrix2 - pointer conflict!!
		for(int j=0;j<faceCounter;j++){
			prevResMat(j,_) = newResMat(j,_);
		}

//		Rcout << "newResMat: " << newResMat << "\n";
	}

	return newResMat;
}




// [[Rcpp::export]]
NumericVector GetPatch_(NumericMatrix faceNeighbours, LogicalVector activeFaces, int startFace, int maxRound){
	int nF =faceNeighbours.rows();

	// the procedure will start with this
	NumericVector start;

	start=faceNeighbours(startFace,_);

	// declare container for the results of one call
	NumericVector currentRound;

	// the faces in the patch
	NumericVector checked(nF); //max set by the number of faces
	int checkedNo = start.size();
	int nLook = start.size();

//	Rcout << "start: " << start << "\n";

	//initialize - the checked thing
	for(int i=0; i<checkedNo;i++){
		checked(i) = start(i);
	}

	//initialize the lookup vector
	NumericVector Lookup(nF); //max set by the number of faces
	for(int i=0; i < checkedNo;i++) {
		Lookup(i) = start(i);
	}


	int totalNewFaces=1;

	NumericVector NewLookup(nF);

	// will be used to set the maximum amount of trials too
	int roundCount =0;

	int newFaces=1;

	int present=0;

	int checkStart;
	int tempStart=0;

	NumericVector indices(nF);
	//while there are more rounds
	while(totalNewFaces>0 && roundCount <maxRound){
//		Rcout << "\n";
//		Rcout << "started round " << roundCount << "\n";

		// the number of faces found in this round
		totalNewFaces =0;

		// for every face in the current round
		for(int i=0;i<nLook;i++){

			// new faces found for this face
			newFaces=0;

			// look up

			currentRound = faceNeighbours(Lookup(i),_);
//			Rcout << "\n";

//			Rcout << "lookin up face: " << Lookup(i) << "\n";
//			Rcout << "Faces found in current round: "<< currentRound << "\n";

			// calculate the number of faces that needs to be checked
			tempStart = indices(roundCount);
			if(tempStart<0){
				checkStart=0;
			}else {
				checkStart = tempStart;
			}
//			Rcout << "checking stored faces between " << checkStart << " and " << checkedNo << "\n";

			//for every potential new face
			for(int j=0;j<currentRound.size();j++){
				// check for the presence of currentRound(j) in the checked vector
				present =0;



				for(int k=checkStart;k<checkedNo;k++){
					if(currentRound(j) == checked(k)){
						present++;
					}
				}


				//if it is not present and it is a face that is considered active
				if(present==0 && activeFaces(currentRound(j))){
					//add to the lookup of the next round
					NewLookup(totalNewFaces) = currentRound(j);
					// put to the checked
					checked(checkedNo+newFaces) = currentRound(j);
					newFaces++;

					// the total number of faces found in this round
					totalNewFaces++;
				}

			}
			// the number of faces checked for additional faces
			checkedNo= checkedNo+newFaces;
//			Rcout << "checkedNo: " << checkedNo << "\n";

//			Rcout << "newFaces: " << newFaces << "\n";
		}

//		Rcout << "\n";
//		Rcout << "NewLookup: " << NewLookup << "\n";
//		Rcout << "totalNewFaces: " << totalNewFaces << "\n";

		for(int i=0;i<totalNewFaces;i++){
			Lookup(i) = NewLookup(i);
		}

		indices(roundCount+2)= checkedNo-totalNewFaces;
//		Rcout << "indices value stored in the round: " << checkedNo-totalNewFaces << "\n";

		nLook = totalNewFaces;
		roundCount++;
	}

	NumericVector endObj(checkedNo);
	for(int i=0;i<checkedNo;i++){
		endObj(i) = checked(i);
	}

	return endObj;

}

// [[Rcpp::export]]
NumericMatrix Partition_(NumericMatrix faceNeighbours, LogicalVector activeFaces, int maxRound){

	// the number of faces
	int nF=faceNeighbours.rows();

	//final object will be a matrix: 1st col: face row, 2nd col: partition number
	NumericMatrix tempObj(nF,2);

	// numeric vector containing a single partition
	NumericVector currentPart;

	// the number of partitions
	int counter=0;

	// the number of faces that are assigned to a partition
	int nCheckedFaces=0;

	// the number of faces in the current partition
	int nNewPart;

	// logical vector indicating whether a face is assigned or not
	bool*checked=new bool[nF];

	// start looking up for partitions with this
	int startFace = 0;

	//for the index check
	bool checkSwitch;

	while(0<=startFace){
//		Rcout << "start of part loop: " << counter << "\n";


		//get the next partition
		currentPart = GetPatch_(faceNeighbours, activeFaces, startFace, maxRound);
		nNewPart = currentPart.size();

		for(int i=0;i<nNewPart;i++){
			tempObj(nCheckedFaces+i,0) = currentPart(i);
			tempObj(nCheckedFaces+i,1) = counter;

			//check the faces
			for(int j=0;j<nF;j++){
				if(currentPart(i)==j){
					checked[j] = 1;
				}
			}
		}

//		Rcout << "found partition: " << nNewPart << "faces" << "\n";

		// count the checked faces
		nCheckedFaces = nCheckedFaces+nNewPart;

		//remaining faces to be checked
		startFace = -99;

		checkSwitch=1;

		for(int index=0; index<nF; index++){
			if(checkSwitch){
				if(checked[index]==0 && activeFaces(index)){
					//store once
					startFace=index;
					checkSwitch=0;
				}

			}

		}
		//if the index is still within the bounds of the matrix

//		Rcout << "will start new search from : " << startFace << "\n";

//		Rcout << "end of part loop: " << counter << "\n";

		//finished partitions
		counter++;
	}

	delete[] checked;

	// get rid of the other values
	NumericMatrix endObj(nCheckedFaces,2);
	for(int i=0;i<nCheckedFaces;i++){
		endObj(i,_) = tempObj(i,_);
	}

	return endObj;

}





//// [[Rcpp::export]]
//NumericMatrix NormalGreatCircle_(NumericVector coord1, NumericVector coord2, NumericVector origin, int breaks, double pi){
// NumericMatrix tempObj;
//
//	// given that there are no breaks
//	if(breaks!=0){
//		//a temporary object
//			double radius = dist(coord1, origin);
//			NumericMatrix oldGC = GreatCircle_(coord1, coord2, origin, 3, pi);
//
//			NumericVector u(3);
//			NumericVector v(3);
//
//
//			// for all the coordinates
//			for(int i=0;i<3;i++){
//				u(i) = oldGC(1,i)-oldGC(0,i);
//				v(i) = oldGC(2,i)-oldGC(0,i);
//			}
//
//			NumericVector cp(3);
//
//			cp(0) = u(1)*v(2)-u(2)*v(1);
//			cp(1) = u(2)*v(0)-u(0)*v(2);
//			cp(2) = u(0)*v(1)-u(1)*v(0);
//
//			NumericVector pointOnPerp(3);
//
//			double rat= sqrt(pow(cp(0),2.0)+pow(cp(1),2.0)+pow(cp(2),2.0))*radius;
//
//			pointOnPerp(0) = cp(0)/rat;
//			pointOnPerp(1) = cp(1)/rat;
//			pointOnPerp(2) = cp(2)/rat;
//
//		tempObj = GreatCircle_(oldGC(0,_),pointOnPerp,origin, breaks, pi);
//	}
//
//	// do you want to return the first and the last values or only the new coords?
//
//	//declare
//	NumericMatrix endObj(breaks+1, 3);
//
//	if(breaks==0){
//		endObj(0,_) = coord1;
//
//
//	}else{
//
//		endObj = tempObj;
//
//	}
//
//	return endObj;
//
//}

// [[Rcpp::export]]
NumericMatrix Refine2d_(NumericMatrix From, int breaks){
	int origRows = From.rows();

	NumericMatrix newMat ((origRows-1)*(breaks+1)*2,2);
	NumericVector vect(2);
	NumericVector p0(2);
	NumericVector p1(2);


	int counter =0;

	for(int i=0; i<(origRows-1); i++){
		p0 =  From(i,_);
		p1 =  From(i+1,_);
		vect(0) = p1(0) - p0(0);
		vect(1) = p1(1) - p0(1);

		for(int j=0; j<(breaks+1); j++){
			newMat(counter,0) = p0(0)+vect(0)/(breaks+1)*j;
			newMat(counter,1) = p0(1)+vect(1)/(breaks+1)*j;
			counter++;
		}

	}

	NumericMatrix endMat(counter, 2);
	for(int i=0; i<counter; i++){
		endMat(i,0) = newMat(i,0);
		endMat(i,1) = newMat(i,1);
	}



	return endMat;
}


// [[Rcpp::export]]
NumericVector SphericalTriangleCenter_(NumericVector v0, NumericVector v1, NumericVector v2, NumericVector origin){
	NumericVector endObj(3);
	NumericVector planeCenter(3);

	for(int i=0; i<3; i++){
		planeCenter(i) = (v0(i)+v1(i)+v2(i))/3;
	}

	double distOri = dist(planeCenter, origin);
	double radius = dist(v0, origin);

	for(int i=0; i<3; i++){
		endObj(i) = (planeCenter(i)-origin(i))/distOri*radius;
	}

	return endObj;
}

// [[Rcpp::export]]
NumericMatrix EdgesFromPoints_(NumericMatrix verts, NumericVector howMany, NumericVector origin){
	// the number of points
	int nPoints = verts.nrow();

	// we do not know how many edge-pairs will be (do not want to calc.)
	double*tempObjCol1 = new double[nPoints*15];
	double*tempObjCol2 = new double[nPoints*15];

	//the number of distances from a point
	double *dists = new double[nPoints];
	double *sorted = new double[nPoints];

	//storing
	int counter =0;


	//for every point
	for(int i=0; i<nPoints;i++){
		// calculate the distance to every other point
		for(int j=0; j<nPoints;j++){
			dists[j] = ArcDist_(verts(i,_), verts(j,_), origin, 1);
			sorted[j]=dists[j];
		}
		//sort the distances
	  std::sort(sorted,sorted+nPoints);
		//the closest howMany(i)+1 elements are forming edges+self
		for(int k=1; k<=howMany(i); k++){
			//look up the same values
			for(int l=0; l<nPoints;l++){
				//when found, store the pair
				if(sorted[k]==dists[l]){
					tempObjCol1[counter] = i;
					tempObjCol2[counter] = l;
					counter++;
				}
			}
		}

	}

	delete[] dists;
	delete[] sorted;


	//clean up the results
	double*pair = new double[2];
	for(int i =0; i<counter; i++){
		*pair = tempObjCol1[i];
		*(pair+1) = tempObjCol2[i];
		//organize them
		std::sort(pair, pair+2);
		tempObjCol1[i]= pair[0];
		tempObjCol2[i]= pair[1];
	}

	delete[] pair;

//	return newTemp;


	//paste in the final object - max determined by prev loop

	int max= counter;

	//get rid of the duplicates
	int*col1 = new int[max*2];
	int*col2 = new int[max*2];

	//initialize
	counter=1;
	col1[0] = tempObjCol1[0];
	col2[0] = tempObjCol2[0];


	int present;
	for(int i=0; i<max;i++){
		//check whether the pair is already in the stack
		present=0;
		for(int k=0; k<counter;k++){
			if(col1[k]==tempObjCol1[i] && col2[k]==tempObjCol2[i]){
				present++;
			}
		}
		if(present==0){
			col1[counter] = tempObjCol1[i];
			col2[counter] = tempObjCol2[i];
			counter++;
		}

	}

	int size = counter;
	NumericMatrix endObj(size,2);

	for(int i=0;i<size; i++){
			endObj(i,0)= col1[i];
			endObj(i,1)= col2[i];
	}


	delete[] col1;
	delete[] col2;



	// delete the tempobj columns
	delete[] tempObjCol1;
	delete[] tempObjCol2;

	return endObj;

}


// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
   NumericVector y = clone(x);
   std::sort(y.begin(), y.end());
   return y;
}


// [[Rcpp::export]]
NumericMatrix EdgesToFaces_(NumericMatrix edges){

	int nPoints = edges.nrow();

	//2 member array with point indices of the current row
	int*actEdge = new int[2];

	//pairs indices where either one is present
	bool*b1 = new bool[nPoints];
	bool*b2 = new bool[nPoints];

	//the number of TRUE values in b1 and b2
	int sumb1;
	int sumb2;

	//sizes of the u1 and u2 arrays
	int u1Size;
	int u2Size;
	int doubleSize;
	int chNewSize;

	//recurring temporary values
	int present;
	int counter;

	//temporary return value
	NumericMatrix endObj(nPoints*2, 3);

	//count the number of faces
	int globalCounter =0;

	for(int i=0;i<nPoints;i++){
		actEdge[0] = edges(i,0);
		actEdge[1] = edges(i,1);


	//the first point
		//1. b1 logical vector
		sumb1 = 0;
		for(int j=0;j<nPoints;j++){
			b1[j]= actEdge[0]==edges(j,0) || actEdge[0]==edges(j,1);
			if(b1[j]==1) sumb1++;
		}

		//2. u1 array
		int*u1Multi = new int[sumb1*2];
		counter=0;

		// concatenate the two parts where one of the rows is tru
		for(int j=0;j<nPoints;j++){
			if(b1[j]){
				u1Multi[counter] = edges(j,0);
				counter++;
				u1Multi[counter] = edges(j,1);
				counter++;
			}
		}

		// the unique function in R
		//get rid of the double values
		int*u1 = new int[sumb1*2];

		//initialize
		u1[0] = u1Multi[0];

		counter =1;
		for(int j=1; j<sumb1*2;j++){
			//look up the copy vector whether the new element is present
			present =0;
			for(int k=0; k<counter;k++){
				if(u1Multi[j]==u1[k]) present++;
			}

			//if the value was not found
			if(present==0){
				u1[counter] = u1Multi[j];
				counter++;
			}

		}

		// clean up the things that lead to u1
		delete[] u1Multi;
		u1Size = counter;

	/////////////////////////////////////
	//the second point
		//1. b2 logical vector
		sumb2 = 0;
		for(int j=0;j<nPoints;j++){
			b2[j]= actEdge[1]==edges(j,0) || actEdge[1]==edges(j,1);
			if(b2[j]==1) sumb2++;
		}

		//2. u1 array
		int*u2Multi = new int[sumb2*2];
		counter=0;

		// concatenate the two parts where one of the rows is tru
		for(int j=0;j<nPoints;j++){
			if(b2[j]){
				u2Multi[counter] = edges(j,0);
				counter++;
				u2Multi[counter] = edges(j,1);
				counter++;
			}
		}

		// the unique function in R
		//get rid of the double values
		int*u2 = new int[sumb2*2];

		//initialize
		u2[0] = u2Multi[0];

		counter =1;
		for(int j=1; j<sumb2*2;j++){
			//look up the copy vector whether the new element is present
			present =0;
			for(int k=0; k<counter;k++){
				if(u2Multi[j]==u2[k]) present++;
			}

			//if the value was not found
			if(present==0){
				u2[counter] = u2Multi[j];
				counter++;
			}

		}

		// clean up the things that lead to u2
		delete[] u2Multi;
		u2Size = counter;

	////Double faces
		// present
		bool*logDoub = new bool[u1Size];

		// the values in the chDouble array
		counter=0;

		int *chDouble = new int[u1Size];
		for(int j=0; j<u1Size;j++){
			//initialize
			logDoub[j] =0;
			for(int k=0;k<u2Size;k++){
				//if the u1[j] element is present in u2
				if(u1[j]==u2[k]){
					logDoub[j]=1;
				}
			}
			//if the current element is present u1[j]
			if(logDoub[j]){
				//copy it over
				chDouble[counter] = u1[j];
				counter++;
			}
		}
		doubleSize= counter;
		delete[] logDoub;


	// get the new points in the faces
		int*chNew = new int[doubleSize];
		counter =0;
		for(int j=0;j<doubleSize;j++){
			if(!(chDouble[j]== actEdge[0] || chDouble[j]== actEdge[1])){
				chNew[counter] = chDouble[j];
				counter++;
			}
		}
		chNewSize = counter;

	//store the results (new faces)
		endObj(globalCounter,0) = actEdge[0];
		endObj(globalCounter,1) = actEdge[1];
		endObj(globalCounter,2) = chNew[0];
		endObj(globalCounter,_) = stl_sort(endObj(globalCounter,_));
		globalCounter++;

		// in case the opposing face is still present
		if(chNewSize>1){
			endObj(globalCounter,0) = actEdge[0];
			endObj(globalCounter,1) = actEdge[1];
			endObj(globalCounter,2) = chNew[1];
			endObj(globalCounter,_) = stl_sort(endObj(globalCounter,_));
			globalCounter++;
		}



		delete[] u1;
		delete[] u2;
		delete[] chDouble;
		delete[] chNew;

	}
	//get rid of the unneeded values
	NumericMatrix endObj2(globalCounter, 3);
	for(int i=0;i<globalCounter;i++){
		endObj2(i,_) = endObj(i,_);
	}

	delete[] actEdge;
	delete[] b1;
	delete[] b2;

	int max=globalCounter;

	//get rid of the duplicates
	int*col1 = new int[max*2];
	int*col2 = new int[max*2];
	int*col3 = new int[max*2];

	//initialize
	counter=1;
	col1[0] = endObj2(0,0);
	col2[0] = endObj2(0,1);
	col3[0] = endObj2(0,2);



	for(int i=0; i<max;i++){
		//check whether the pair is already in the stack
		present=0;
		for(int k=0; k<counter;k++){
			if(col1[k]==endObj2(i,0) && col2[k]==endObj2(i,1) && col3[k]==endObj2(i,2)){
				present++;
			}
		}
		if(present==0){
			col1[counter] = endObj2(i,0);
			col2[counter] = endObj2(i,1);
			col3[counter] = endObj2(i,2);
			counter++;
		}

	}

	NumericMatrix endObj3(counter,3);
	for(int i=0;i<counter;i++){
		endObj3(i,0) = col1[i];
		endObj3(i,1) = col2[i];
		endObj3(i,2) = col3[i];

	}

	delete[] col1;
	delete[] col3;
	delete[] col2;

	return endObj3;
}

// [[Rcpp::export]]
List TriangleTesselation_(NumericVector v0, NumericVector v1, NumericVector v2, NumericVector origin, int lineBreak){

	// array of points
		int currentRowNo = 3;

	//A. primary splitting
	//1.the first edge [1] - [2]
		NumericMatrix split12;
		split12 = SplitArc_(v0,v1, origin, lineBreak, 0);

	// 2. split the second edge [1] - [3]
		NumericMatrix split13;
		split13 = SplitArc_(v0,v2, origin, lineBreak, 0);

	// 3. split the third edge [2] - [3]
		NumericMatrix split23;
		split23 = SplitArc_(v1,v2, origin, lineBreak, 0);

		int splitSize = split12.rows();
		int newRows = splitSize*3;

		// create a new container
			int newSize = currentRowNo+newRows;
			double*newPointsCol0 = new double[newSize];
			double*newPointsCol1 = new double[newSize];
			double*newPointsCol2 = new double[newSize];

		// fill the container
			// with the original values at the beginning
			// x coordinates
			newPointsCol0[0] = v0(0);
			newPointsCol0[1] = v1(0);
			newPointsCol0[2] = v2(0);

			// y coordinates
			newPointsCol1[0] = v0(1);
			newPointsCol1[1] = v1(1);
			newPointsCol1[2] = v2(1);

			// z coordinates
			newPointsCol2[0] = v0(2);
			newPointsCol2[1] = v1(2);
			newPointsCol2[2] = v2(2);

		// with the new values
		// values of the split12
			int offset =3;
			for(int i=0;i<splitSize;i++){
				newPointsCol0[i+offset] = split12(i,0);
				newPointsCol1[i+offset] = split12(i,1);
				newPointsCol2[i+offset] = split12(i,2);
			}

		// values of the split13
		offset = 3+splitSize;
			for(int i=0;i<splitSize;i++){
				newPointsCol0[i+offset] = split13(i,0);
				newPointsCol1[i+offset] = split13(i,1);
				newPointsCol2[i+offset] = split13(i,2);
			}

		// values of the split23
			offset = 3+splitSize*2;
			for(int i=0;i<splitSize;i++){
				newPointsCol0[i+offset] = split23(i,0);
				newPointsCol1[i+offset] = split23(i,1);
				newPointsCol2[i+offset] = split23(i,2);
			}

		//get rid of the duplicates
			double*col1 = new double[newSize*2];
			double*col2 = new double[newSize*2];
			double*col3 = new double[newSize*2];

			//initialize
			int counter=1;
			col1[0] = newPointsCol0[0];
			col2[0] = newPointsCol1[0];
			col3[0] = newPointsCol2[0];

			int present;

			for(int i=0; i<newSize;i++){
				//check whether the pair is already in the stack
				present=0;
				for(int k=0; k<counter;k++){
					if(col1[k]==newPointsCol0[i] && col2[k]==newPointsCol1[i] && col3[k]==newPointsCol2[i]){
						present++;
					}
				}
				if(present==0){
					col1[counter] = newPointsCol0[i];
					col2[counter] = newPointsCol1[i];
					col3[counter] = newPointsCol2[i];
					counter++;
				}

			}
			int uniqueSize = counter;

			delete[] newPointsCol0;
			delete[] newPointsCol1;
			delete[] newPointsCol2;

		// the number of points bound by two and four edges
			int bound4 = uniqueSize;


	//B. secondary splitting
	//only if lineBreak >1
		if(lineBreak>1){

		//1. containers - oversized on purpose!!
			int sumPoints=0;
			for(int i=0;i<(lineBreak+1);i++){
				sumPoints=sumPoints+i;
			}

			double*npCol0 = new double[sumPoints*3];
			double*npCol1 = new double[sumPoints*3];
			double*npCol2 = new double[sumPoints*3];

		//2. connect the two split edges [1-2] - [1-3]
			NumericMatrix newLine;

			NumericVector point1(3);
			NumericVector point2(3);

			offset = 0;
			for(int i=0;i<lineBreak;i++){
				// split the points
				point1 = split12(i+1,_);
				point2 = split13(i+1,_);
				newLine = SplitArc_(point1, point2, origin, i, 1);
				newRows=newLine.rows();
				if(newRows>0){
					//copy the new rows
					for(int j=0;j<newRows;j++){
						npCol0[j+offset] = newLine(j,0);
						npCol1[j+offset] = newLine(j,1);
						npCol2[j+offset] = newLine(j,2);
					}
					offset= offset+newRows;
				}
			}

		//3. connect split lines [2-1] - [2-3]
			//reverse the line 1-2
			NumericMatrix split21(splitSize,3);
			for(int i=0;i<splitSize;i++){
				split21(splitSize-i-1,_) = split12(i,_);
			}

			for(int i=0;i<lineBreak;i++){
				// split the points
				point1 = split21(i+1,_);
				point2 = split23(i+1,_);
				newLine = SplitArc_(point1, point2, origin, i, 1);

				newRows=newLine.rows();
				if(newRows>0){
					//copy the new rows
					for(int j=0;j<newRows;j++){
						npCol0[j+offset] = newLine(j,0);
						npCol1[j+offset] = newLine(j,1);
						npCol2[j+offset] = newLine(j,2);
					}
					offset= offset+newRows;
				}
			}

		//4. connect split lines [3-1] - [3-2]
			//reverse the line 1-3
			NumericMatrix split31(splitSize,3);
			for(int i=0;i<splitSize;i++){
				split31(splitSize-i-1,_) = split13(i,_);
			}

			//reverse the line 2-3
			NumericMatrix split32(splitSize,3);
			for(int i=0;i<splitSize;i++){
				split32(splitSize-i-1,_) = split23(i,_);
			}


			for(int i=0;i<lineBreak;i++){
				// split the points
				point1 = split31(i+1,_);
				point2 = split32(i+1,_);
				newLine = SplitArc_(point1, point2, origin, i, 1);

				newRows=newLine.rows();
				if(newRows>0){
					//copy the new rows
					for(int j=0;j<newRows;j++){
						npCol0[j+offset] = newLine(j,0);
						npCol1[j+offset] = newLine(j,1);
						npCol2[j+offset] = newLine(j,2);
					}
					offset= offset+newRows;
				}
			}

			int pointNo = offset/3;
			int matSize = lineBreak-1;

		//5. averaging the results
		//5.A. replicating the pattern with indices- to know what points to average
			//a. matrix 1
		//	int m1 [matSize][matSize];
			int *m1 = new int[matSize * matSize];
			int count = 0;
			for(int i=0; i<matSize;i++){
				for (int j=0;j<(i+1);j++){
				//	m1[i][j] = count;
					m1[j*matSize+i] = count;
					count++;
				}
			}

			//b. matrix 2
			// -/1 creating
		//	int m2[matSize][matSize];
			int *m2 = new int[matSize * matSize];
			for(int j=0;j<matSize;j++){
				for(int i=matSize-1-j;i<matSize;i++){
				//	m2[i][j] = count;
					m2[j*matSize+i] = count;
					count++;
				}

			}

			//-/2 transforming
		//	int tm2[matSize][matSize];
			int *tm2 = new int[matSize * matSize];
			for(int i=0; i<matSize;i++){
				for(int j=0;j<(i+1);j++){
				//	tm2[i][j] = m2[i][matSize-i+j-1];
					tm2[j*matSize+i] = m2[(matSize-i+j-1)*matSize+i];
				}
			}

			delete[] m2;
			
			//c. matrix 3
		//	int m3[matSize][matSize];
			int *m3 = new int[matSize * matSize];
			for(int j=matSize-1;j> -1; j--){
				for(int i=j; i<matSize; i++){
				//	m3[i][j]=count;
					m3[j*matSize+i]=count;
					count++;
				}
			}

		// 5.B. "average" the corresponding coordinates
			//temporary containers of indices
			int index1;
			int index2;
			int index3;

			// containers for the new columns (average npCol)
			// size should be exactly pointNo, but who knows?
			double *aNpCol0 = new double[pointNo+15];
			double *aNpCol1 = new double[pointNo+15];
			double *aNpCol2 = new double[pointNo+15];
			int storeCounter =0;

			// temporary containers for the new vectors
			NumericVector tempV0(3);
			NumericVector tempV1(3);
			NumericVector tempV2(3);
			NumericVector tempCenter(3);

			// replicate the matrix form and go through the coordinates
			for(int i=0; i<matSize;i++){
				for (int j=0;j<(i+1);j++){
					//the first temporary vertex
				//	index1= m1[i][j];
					index1= m1[j*matSize+i];
					tempV0(0) = npCol0[index1];
					tempV0(1) = npCol1[index1];
					tempV0(2) = npCol2[index1];

					//the second temporary vertex
				//	index2= tm2[i][j];
					index2= tm2[j*matSize+i];
					tempV1(0) = npCol0[index2];
					tempV1(1) = npCol1[index2];
					tempV1(2) = npCol2[index2];

					//the third temporary vertex
				//	index3= m3[i][j];
					index3= m3[j*matSize+i];
					tempV2(0) = npCol0[index3];
					tempV2(1) = npCol1[index3];
					tempV2(2) = npCol2[index3];

					//calculate the center
					tempCenter = SphericalTriangleCenter_(tempV0, tempV1, tempV2, origin);

					// store its coordinates
					aNpCol0[storeCounter] = tempCenter(0);
					aNpCol1[storeCounter] = tempCenter(1);
					aNpCol2[storeCounter] = tempCenter(2);
					storeCounter++;

				}
			}
			
			delete[] m1;
			delete[] tm2;
			delete[] m3;
			

			delete[] npCol0;
			delete[] npCol1;
			delete[] npCol2;

			// add up the two parts
			int totalPoints = uniqueSize+storeCounter;
			NumericMatrix vertices(totalPoints,3);
			//first part
			for(int i=0;i<uniqueSize;i++){
				vertices(i,0) = col1[i];
				vertices(i,1) = col2[i];
				vertices(i,2) = col3[i];
			}

			delete[] col1;
			delete[] col2;
			delete[] col3;

			//second part
			for(int i=0;i<storeCounter;i++){
				vertices(uniqueSize+i,0) = aNpCol0[i];
				vertices(uniqueSize+i,1) = aNpCol1[i];
				vertices(uniqueSize+i,2) = aNpCol2[i];

			}

			delete[] aNpCol0;
			delete[] aNpCol1;
			delete[] aNpCol2;

	//		return vertices;


		//C. register the number of edges that will belong to the points
			int nVert= vertices.rows();

			NumericVector allEdge(nVert);

			allEdge (0) = 2;
			allEdge (1) = 2;
			allEdge (2) = 2;

			for(int i=3;i<bound4;i++){
				allEdge(i) = 4;
			}
			if(bound4!=nVert){
				for(int i=bound4;i<nVert;i++){
					allEdge(i) = 6;
				}
			}
//			return allEdge;
		//D. contstruction of the edges
			NumericMatrix edges;
			edges = EdgesFromPoints_(vertices, allEdge, origin);


		//E. contstruction of the faces
			NumericMatrix faces;
			faces = EdgesToFaces_(edges);

			return Rcpp::List::create(Rcpp::Named("vertices") = vertices,
                          Rcpp::Named("faces") = faces,
						  Rcpp::Named("edgeVert") = bound4);

		}else{

			NumericMatrix vertices(uniqueSize,3);
			for(int i=0;i<uniqueSize;i++){
				vertices(i,0)= col1[i];
				vertices(i,1)= col2[i];
				vertices(i,2)= col3[i];
			}

			delete[] col1;
			delete[] col2;
			delete[] col3;

	//		return vertices;

		//C. register the number of edges that will belong to the points
			int nVert= vertices.rows();

			NumericVector allEdge(nVert);

			allEdge (0) = 2;
			allEdge (1) = 2;
			allEdge (2) = 2;

			for(int i=3;i<bound4;i++){
				allEdge(i) = 4;
			}
			if(bound4!=nVert){
				for(int i=bound4;i<nVert;i++){
					allEdge(i) = 6;
				}
			}

//			return allEdge;
		//D. contstruction of the edges
			NumericMatrix edges;
			edges = EdgesFromPoints_(vertices, allEdge, origin);


		//E. contstruction of the faces
			NumericMatrix faces;
			faces = EdgesToFaces_(edges);

			return Rcpp::List::create(Rcpp::Named("vertices") = vertices,
                          Rcpp::Named("faces") = faces,
						  Rcpp::Named("edgeVert") = bound4);
		}

}

// [[Rcpp::export]]
NumericVector SizeEstimate_(NumericVector tesselation){
	int nDegree=tesselation.size();

	NumericVector endObj(2);
	// the number of vertices
	endObj(0) = 12;

	// the number of faces
	endObj(1) = 20;

	int t;
	int oneTri;
	int prod;

	NumericVector prev(2);
	for(int i=0; i<nDegree;i++){

		//results of the previous loop
		prev(0) = endObj(0);
		prev(1) = endObj(1);

		// the Number of vertices - face edges will be duplicates!!
		t=tesselation[i];
		oneTri=0;
		for(int j=0;j<(t+2);j++){
			oneTri=oneTri+j;
		}

		endObj(0)= endObj(0)+prev(1)*oneTri;

		// the number of faces
		prod=20;
		for(int j=0; j<i+1;j++){
			prod= prod*pow(tesselation[j],2.0);

		}
		endObj(1)= prev(1)+prod;
	}

	return endObj;
	///will overestimate!!! vertices by a lot,

}



// [[Rcpp::export]]
List IcosahedronTesselation_(NumericMatrix oldV, NumericMatrix oldF, NumericVector tesselation, NumericVector origin){

	//the number of rounds (degrees)
	int nDegree = tesselation.size();

	// tesselation in the active round
	int tessel;

	// estimate (over) the necessary memory size of the skeleton object
	NumericVector TempEst = SizeEstimate_(tesselation);

	// first: estimated number of vertices (overestimation!)
	int estimV = TempEst(0)+150;

	// second: estimated number of faces
	int estimF = TempEst(1);

	// the end variables
	NumericMatrix v(estimV,3);
	NumericMatrix f(estimF,5);

	//paste in the starting values
	// vertices
	for(int i=0;i<12;i++){
		v(i,_) = oldV(i,_);
	}
	// faces
	for(int i=0;i<20;i++){
		f(i,_) = oldF(i,_);
	}

	//indexing the faces that are to be tesselated
	int nStartF=0;
	int nEndF=20;

	// number of faces to be tesselated
	int nFaceDiff;

	// variables to hold the indices of the vertices that will be tesselated
	int p0ind;
	int p1ind;
	int p2ind;

	// coordinates of tesselated points
	NumericVector v0;
	NumericVector v1;
	NumericVector v2;

	// number of new vertices and faces
	int nNewVert;
	int nNewFace;

	//number of rows in the v and f matrices
	int nF= 20;
	int nV= 12;

	//auxilliary declarations
	int present;
	int which;
	int found;
	int nReplace;
	int index;
	int counter;

	//loop for degrees of tesslelation
	for(int degree=0; degree<nDegree; degree++){
//	int degree=0;

		//number of breaking points
		tessel = tesselation(degree)-1;

		//the number of faces to be tesselated
		nFaceDiff= nEndF-nStartF;
//	nFaceDiff=18;
		// loop for every face to be tesselated
		for(int fa=0;fa<nFaceDiff;fa++){
	//	int fa=0;

			// the vertex indices
			p0ind = f(nStartF+fa,0);
			p1ind = f(nStartF+fa,1);
			p2ind = f(nStartF+fa,2);

			// the vertex coordinates
			v0 = v(p0ind,_);
			v1 = v(p1ind,_);
			v2 = v(p2ind,_);

		//--// tesselation of the face
			List OneFace = TriangleTesselation_(v0,v1,v2, origin, tessel);

			// get the components
			NumericMatrix vertices = OneFace(0);
			NumericMatrix faces = OneFace(1);
	//		NumericMatrix origFaces= OneFace(1);

			// the parameters
			nNewVert = vertices.rows();
			nNewFace = faces.rows();

		//--// search for the already existing vertices  (+10 to avoid out of bounds)
			int*replThis = new int[nNewVert+10];
			int*withThis = new int[nNewVert+10];

			counter=0;
			for(int i=0;i<nNewVert;i++){
				// check whether the new vertex in question is present
				present =0;
				which= -99;
				for(int j=0; j<nV;j++){
					if(vertices(i,0)==v(j,0) &&
						vertices(i,1)==v(j,1) &&
						vertices(i,2)==v(j,2)){
							present=1;
							which = j;
					}
				}

				//if it is present
				if(present>0){
					withThis[counter] = which;
					replThis[counter] = i;
					counter++;
				}
			}

			nReplace = counter;

		//--// create a lookup table for the vertex indices
			found =0;

			int*oldVertex = new int[nNewVert];
			int*newVertex = new int[nNewVert];

			counter=0;
			// loop through all the vertices
			for(int i=0;i<nNewVert;i++){
				oldVertex[i] = i;
				present =0;
				for(int j=0;j<nReplace;j++){
					if(i==replThis[j]){
						present++;
						index = j;
					}
				}
				//if present
				if(present>0){
					newVertex[i] = withThis[index];
					found++;

				// if not present before
				}else{
					newVertex[i] = i - found+nV;
					// copy over the new vertex coordinates to the rest
					v(nV+counter,_) = vertices(i,_);
					counter++;

				}

			}
			// the new number of vertices
			nV = counter+nV;

			// clean up the mess
			delete[] replThis;
			delete[] withThis;

		//--// replace vertex indices in the facetables

			// for the rows
			for(int j=0;j<nNewFace;j++){
				//for the columns
				for(int i=0;i<3;i++){
					int toBeChanged=1;
					for(int t=0; t<nNewVert; t++){
						if(faces(j,i)==oldVertex[t] && toBeChanged){
							faces(j,i)=newVertex[t];
							toBeChanged=0;
						}
					}
				}

				// before passing to the next, copy over the face
				f(nF+j,0) = faces(j,0);
				f(nF+j,1) = faces(j,1);
				f(nF+j,2) = faces(j,2);
				f(nF+j,3) = degree;
				f(nF+j,4) = fa;


			}
			// the number of faces in the f table
			nF=nF+nNewFace;

			delete[] oldVertex;
			delete[] newVertex;


		}


		//the next tesselation order - do not run
		nStartF=nEndF;
		nEndF = nF;

	}

	// clean up the variables
	NumericMatrix endV(nV,3);
	NumericMatrix endF(nF,5);

	// first the vertices
	for(int i=0;i<nV;i++){
		endV(i,_) = v(i,_);
	}

	// and second the faces
	for(int i=0;i<nF;i++){
		endF(i,_) = f(i,_);
	}



	return Rcpp::List::create(Rcpp::Named("v") = endV,
                          Rcpp::Named("f") = endF,
						 Rcpp::Named("nStartF") = nStartF);

}

// function after the tesselation:
// expanding the faces matrix to an edges matrix

// [[Rcpp::export]]
NumericMatrix expandFacesToEdges_(NumericMatrix faces){
	int rows= faces.rows();

	NumericMatrix fullEdges(rows*3,2);

	int counter = 0;
	for(int i=0;i<rows;i++){
		// 1 - 2
		fullEdges(counter,0) = faces(i,0);
		fullEdges(counter,1) = faces(i,1);
		counter++;

		// 1 - 3
		fullEdges(counter,0) = faces(i,0);
		fullEdges(counter,1) = faces(i,2);
		counter++;

		// 2 - 3
		fullEdges(counter,0) = faces(i,1);
		fullEdges(counter,1) = faces(i,2);
		counter++;

	}

	// can get rid of the duplicates in R

	return fullEdges;
}

// function to order the grid cells in rings

// [[Rcpp::export]]
List orderTriGrid_(NumericMatrix faces, NumericMatrix neigh, NumericVector startFaces, NumericVector startVert, int nBelts, int nV){

	int nF = faces.rows();

	NumericVector faceOrder(nF);
	NumericVector vertexOrder(nV);

	// the looping variables
	NumericVector prevVerts(nV); // oversized!
	NumericVector newVerts(nV);
	int nPrevVert=5;
	int nvCounter;

	// face results
	NumericVector newFaces(nF);
	NumericVector prevFaces(nF);
	NumericVector prev2Faces(nF);

	int nPrev2Face=0;
	int nPrevFace=5;
	int nfCounter;

	// starting faces
	for(int j=0;j<nPrevFace;j++){
		faceOrder(j) = startFaces(j);
		prevFaces(j) = startFaces(j);
	}

	//starting vertices -all
	for(int k=0;k<(nPrevVert+1); k++){
		vertexOrder(k) = startVert(k);
	}

	// the previous round
	for(int k=0;k<nPrevVert; k++){
		prevVerts(k) = startVert(k+1);
	}

	//total counters
	int fCounter=5;
	int vCounter=6;

	int focusFace=-9;
	int present;
	int present1;
	int present2;
	int firstFace=-9;
	int prevFace=-9;
	int befPrevFace=-99;
	int vLoopCounter1;
	int vLoopCounter2;
	int newVertInd;

	NumericVector potVert(3);
	NumericVector allNeighbours(4);

	bool loopSwitch;
	bool outWhile;

	//start of the belts
	NumericVector belts(nBelts);
	belts(0) = 0;

	// start looping from the second belt
	for(int i=1;i<nBelts; i++)
	{
		belts(i) = fCounter;
//		Rcout << "\nbelt round: " << i << "\n";

//		Rcout << "the vertices found in th previous round are: " << prevVerts << "\n\n";
		//reset the counters
		nvCounter=0;
		nfCounter=0;

		// everything here is for the first face
		// neighbours of the last face in the previous loop
		allNeighbours= neigh(faceOrder(fCounter-1),_);

		outWhile=1;

		while(outWhile){
			//check which of these faces is actually new - first face of the new row
			for(int k=0;k<4;k++){
				// in the previous round
				present1=0;
				for(int j=0; j<nPrevFace; j++){
					if(allNeighbours(k) == prevFaces(j)){
						present1++;
					}
				}
				// in the two rounds before
				present2=0;
				if(i>2){
					for(int j=0; j<nPrev2Face; j++){
						if(allNeighbours(k) == prev2Faces(j)){
							present2++;
						}
					}
				}

				if(present1==0 && present2==0){
					firstFace = allNeighbours(k);
					outWhile=0;
				}else{
					// this is for the weird face!
					if(present2>0){
						allNeighbours = neigh(faceOrder(fCounter-2),_); // then the checks fun again
//						Rcout << "\nstepping back\n";
					}

				}
			}
		}

		// save the first face of the belt!
			//rowdata
			newFaces(nfCounter) = firstFace;
			nfCounter++;

			//total data
			faceOrder(fCounter)=firstFace;
			fCounter++;

			for(int j=0;j<3;j++){
				present=0;
				for(int k=0;k<nPrevVert;k++){
					if(faces(firstFace,j)==prevVerts(k)){
						present++;
					}
				}
				// the new vertex
				if(present==0){
					newVertInd = faces(firstFace,j);
					vertexOrder(vCounter) = newVertInd;
					vCounter++;

					// pass to the next loop
					newVerts(nvCounter) = newVertInd;
					nvCounter++;
				}

			}

//		Rcout << "first face in the belt is " << firstFace << "\n";

		prevFace = firstFace;
		loopSwitch=1;

		vLoopCounter1=0;
		vLoopCounter2=1;

		//order the faces within the loop
		// while there are new faces
		while(loopSwitch){
			if(vLoopCounter1==nPrevVert) vLoopCounter2=0;

			// get the neighbours of the first face
			allNeighbours = neigh(prevFace,_);

//			Rcout << "The neighbours of " << prevFace << " are: " << allNeighbours << "\n";
//			Rcout << "These will be tested for vertices:" << prevVerts[vLoopCounter1] << " and " << prevVerts[vLoopCounter2] << "\n";

			// which of these is new?
			for(int j=0;j<4;j++){

				// the new face cannot be in the previous round
				present=0;

				for(int k=0;k<nPrevFace;k++){
					if(prevFaces(k)==allNeighbours(j)){
						present++;
					}
				}
				// and itself or the one before
				if(allNeighbours(j)==prevFace || allNeighbours(j)== befPrevFace){
					present++;
				}


				// and the face before prevFace

				if(present==0){
					// and should include either prevVerts[vLoopCounter1] or prevVerts[vLoopCounter2]

					present1=0;
					present2=0;
					for(int k=0;k<3;k++){
						if(faces(allNeighbours(j),k) == prevVerts[vLoopCounter1]){
							present1++;

						}
						if(faces(allNeighbours(j),k) == prevVerts[vLoopCounter2]){
							present2++;

						}
					}
					if(present1>0){
						focusFace = allNeighbours(j);
					}
					if(present2>0){
						focusFace = allNeighbours(j);
						// shift the vertex lookup values
						vLoopCounter1++;
						vLoopCounter2++;
					}
				}
			}

//			Rcout << "next found face is:  " << focusFace << "\n";

			potVert= faces(focusFace,_);
			// look up all three vertices - which is new?
			for(int k=0;k<3;k++){
				//new only if it was not present in the previous round
				present=0;
				for(int j=0;j<nPrevVert;j++){
					if(potVert(k)==prevVerts(j)){
						present++;
					}

				}

				//or if it is not present in the current round
				for(int j=0;j<nvCounter;j++){
					if(potVert(k)==newVerts(j)){
						present++;
					}

				}

				//there is a new vertex
				if(present==0){
					newVertInd = potVert(k);

//					Rcout << "next found vertex is:  " << newVertInd << "\n";

					//final output
					vertexOrder(vCounter) = newVertInd;
					vCounter++;

					// pass to the next loop
					newVerts(nvCounter) = newVertInd;
					nvCounter++;
				}
			}

			//loop breaking condition - the first face is the neighbour of the current focus face
			allNeighbours = neigh(focusFace,_);

			// after a couple of faces has been found in the belt
			if(nfCounter>3){
				present=0;
				for(int k=0;k<4; k++){
					if(allNeighbours(k) == firstFace){
						loopSwitch=0;
					}
				}

			}


			// pass to the next loop
				//faces
				befPrevFace = prevFace;
				prevFace=focusFace;
				newFaces(nfCounter) = focusFace;
				nfCounter++;


			//save to the final
				faceOrder(fCounter)=focusFace;
				fCounter++;

		}

		// copy over the values
		for(int j=0;j<nvCounter;j++){
			prevVerts(j) = newVerts(j);
		}

		for(int j=0;j<nPrevFace; j++){
			prev2Faces(j) = prevFaces(j);
		}

		for(int j=0;j<nfCounter;j++){
			prevFaces(j) = newFaces(j);
		}

		// shift the counters
		nPrev2Face = nPrevFace;
		nPrevFace = nfCounter;
		nPrevVert = nvCounter;
	}

	return Rcpp::List::create(Rcpp::Named("faceOrder") = faceOrder,
                          Rcpp::Named("vertexOrder") = vertexOrder,
						  Rcpp::Named("belts") = belts);

}

// this is a function to look up the spherical centers!

// [[Rcpp::export]]
NumericMatrix allTriangleCenters_(NumericMatrix v, NumericMatrix f, NumericVector origin){

	// the rows in faces
	int nFaces = f.rows();

	// the final object of coordinates
	NumericMatrix endObj(nFaces,3);

	// temporary variables for the indices of points
	int pointInd0;
	int pointInd1;
	int pointInd2;

	NumericVector center(3);

	// get the centers
	for(int i=0; i<nFaces; i++){
		pointInd0 = f(i,0);
		pointInd1 = f(i,1);
		pointInd2 = f(i,2);

		center = SphericalTriangleCenter_(v(pointInd0,_),v(pointInd1,_),v(pointInd2,_), origin);

		endObj(i,_) = center;

	}

	return endObj;


}

// function to create the lookup table for the hexagrid

// [[Rcpp::export]]
NumericMatrix CreateHexaSubfaces_(NumericMatrix n, NumericMatrix f, int nV){

	// number of triangles
	int nF = n.rows();
	//number of subfaces
	NumericMatrix endObj(nF*6,5);

	//container to store the neighbours of focal face
	NumericVector actNeigh(4);

	//container to store the vertices of the focal neighbour
	NumericVector neighVerts(3);

	NumericVector tempObj(3);

	int present;
	int counter =0;
	// basic loop for all the faces of the triangular grid
	for(int i=0;i<nF;i++){

		//neighbours of the focal cell
		actNeigh = n(i,_);

		// iterate for all the  neighbours
		for(int j=0;j<4;j++){
			// not for the current face - 3 remains
			if(actNeigh(j)!=i){
				// vertices in the current neighbour
				neighVerts = f(actNeigh(j),_);

				// for all the vertices that are shared with the focal triangle
				for(int k=0;k<3;k++){
					// which is shared? - there will be two
					present=0;
					for(int l=0;l<3;l++){
						if(f(i,l)==neighVerts(k)){
							present++;

						}
					}
					//if it is shared
					if(present>0){

						// the first column is the vertex that is shared has to be < nV
						tempObj(0) = neighVerts(k);

						// facecenter of the neighbour triangle (new vertex of hexagrid)
						// shifted to the second part of the new $v
						tempObj(1) = actNeigh(j)+nV;

						// facecenter of the focal triangle (new vertex of hexagrid)
						// shifted to the second part of the new $v
						tempObj(2) = i+nV;

						// sort the indices - for later uniquing
						tempObj = stl_sort(tempObj);
						for(int u=0;u<3;u++){
							endObj(counter, u)=tempObj(u);
						}

						//indicate the "tesselation level"
						endObj(counter, 3) = -6;

						//which original face does the subface belong to
						endObj(counter, 4) = i;


						counter++;
					}


				}

			}
		}


	}

	return endObj;


}

//function to calculate the faces table of a hexagrid from the subfaces

// [[Rcpp::export]]
NumericMatrix HexaFaces_(NumericMatrix fOrd){
	// the number of rows in the
	int rows=  fOrd.rows();

	//the final internal faces table//12 pentagons!
	NumericMatrix endObj(((rows+12)/6), 6);

	// face information
	NumericMatrix tempObj(6,3);

	//starting
	int prevFace=0;
	int counter=0;
	int prevCounter;

	bool presBool1;
	bool presBool2;

	for(int i=0;i<rows;i++){
		// a new face is found
		if(fOrd(i,0)!=prevFace){
			prevCounter=counter;
			counter=0;

			// do something here with tempObj (the previous faces)
			//look up every position
			//initialize the row
			endObj(prevFace,0) = tempObj(0,1);
			endObj(prevFace,1) = tempObj(0,2);

			LogicalVector temp(prevCounter);
			temp(0)= 1;
			//look up the rest
			for(int j=2;j<prevCounter;j++){
				for(int k=1;k<prevCounter;k++){
					// do not consider the 0
					presBool1 = tempObj(k,1)==endObj(prevFace,j-1);
					presBool2 = tempObj(k,2)==endObj(prevFace, j-1);

					if(presBool1 && !(temp(k))){
						// select the other
						endObj(prevFace, j) = tempObj(k, 2);
						temp(k) = 1;
					}
					if(presBool2 && !(temp(k))){
						// select the other
						endObj(prevFace, j) = tempObj(k, 1);
						temp(k) = 1;
					}
				}

			}


			//store in tempObj the info
			tempObj(counter,_)=fOrd(i,_);
			counter++;

		}else{
			//do not do anything special
			tempObj(counter,_)=fOrd(i,_);
			counter++;
		}

		prevFace=fOrd(i,0);
	}
	prevCounter=counter;

	//run the last round
	//look up every position
	//initialize the row
	endObj(prevFace,0) = tempObj(0,1);
	endObj(prevFace,1) = tempObj(0,2);

	LogicalVector temp(prevCounter);
	temp(0)= 1;
	//look up the rest
	for(int j=2;j<prevCounter;j++){
		for(int k=1;k<prevCounter;k++){
			// do not consider the 0
			presBool1 = tempObj(k,1)==endObj(prevFace,j-1);
			presBool2 = tempObj(k,2)==endObj(prevFace, j-1);

			if(presBool1 && !(temp(k))){
				// select the other
				endObj(prevFace, j) = tempObj(k, 2);
				temp(k) = 1;
			}
			if(presBool2 && !(temp(k))){
				// select the other
				endObj(prevFace, j) = tempObj(k, 1);
				temp(k) = 1;
			}
		}

	}


	return endObj;
}

// [[Rcpp::export]]
NumericMatrix RetrieveIndexMat_(NumericVector indices){
	int all = indices.size();
	NumericVector entries(all);

	NumericMatrix endObj(all, 12);

	int maxInd=0;
	int actInd;
	for(int i=0; i<all; i++){
		actInd= indices(i)-1;
		if(actInd>maxInd){
			maxInd = actInd;
		}

		endObj(actInd,entries(actInd)) = i+1;

		entries(actInd)++;
	}

	NumericMatrix endObj2((maxInd+1), 12);
	for(int i=0;i<(maxInd+1);i++){
		endObj2(i,_) = endObj(i,_);

	}

	return endObj2;
}



double InvWeight_(double* val,double* prob, int n){
	double endVal=0;
	double divVal=0;

	for(int i=0;i<n;i++){
		endVal=endVal+val[i]*(1/prob[i]);
		divVal=divVal+(1/prob[i]);
	}

	endVal=endVal*(1/divVal);

	return endVal;

}

double Mean_(NumericVector val){
	double endVal=0;
	double counter=0;
	int n=val.size();

	for(int i=0;i<n;i++){
		endVal=endVal+val(i);
		counter++;
	}

	endVal= endVal/counter;

	return endVal;
}

// [[Rcpp::export]]
NumericVector InverseWeightByFaceCenter_(NumericMatrix fcNew, NumericVector loc, NumericMatrix n, NumericMatrix fcOld, NumericVector values, NumericVector origin, int deg){
	// number of faces in the new grid
	int newN = fcNew.rows();

	// index of the current neighbour
	int actNeigh;

	// final  container
	NumericVector newVals(newN);

	// temporary vectors
	NumericVector dists(newN);
	NumericVector vals(newN);
	NumericVector allNeighbours(newN);
	NumericVector neighCounts(newN);

	NumericVector startNeighbours;
	int counter=0;
	int newNeighbour;
	int present;

	// for every new face(cetner)
	for(int i=0;i<newN;i++){
//		Rcout << "\nnew point no " << i << "\n";

		counter=0;

		startNeighbours =n(loc(i),_);
		for(int j=0;j<4;j++){
			allNeighbours(j)=startNeighbours(j);
			counter++;
		}

		neighCounts(0)=4;

//		Rcout << "startneighbours: " << startNeighbours << "\n";

		// for every degree of neighbourhood
		for(int k=0; k<deg;k++){
//			Rcout << "starting degree loop: " << k << "\n";

			// index of the neighbour
			for(int j=0;j<neighCounts(k);j++){
//				Rcout << "neighbour search of  " << allNeighbours(j) << "\n";
				// for all the new neighbours
				for(int l=0;l<4;l++){

					newNeighbour = n(allNeighbours(j),l);

//					Rcout << "potential new neighbour:  " << newNeighbour << "\n";
					// check whether this was not present previously
					present=0;
					for(int m=0;m<counter;m++){
						if(newNeighbour == allNeighbours(m)){
							present++;
						}

					}

					if(present == 0){
						allNeighbours(counter)=newNeighbour;
						counter++;

//						Rcout << "not present before:  " << newNeighbour << "\n";
					}
				}


			}
			neighCounts(k+1) = counter;

		}

//		Rcout << "starting  the distance and value extraction\n";
//		Rcout << "counter: "  << counter << "\n";
//		if(i==0){
//			Rcout << "allNeighbours: "  << allNeighbours << "\n";
//		}
		double* dists = new double[counter];
		double* vals = new double[counter];

		for(int j=0;j<counter;j++){
			actNeigh= allNeighbours(j);
			// distance between the new face center and the center of the older
			dists[j] = ArcDist_(fcNew(i,_),fcOld(actNeigh,_), origin, 0);
			//values that will make up the new one
			vals[j] = values(actNeigh);
		}


//		Rcout << "weighting\n";
		newVals(i) = InvWeight_(vals, dists, counter);
		delete[] dists;
		delete[] vals;
	}

	return newVals;
}

// [[Rcpp::export]]
NumericVector OccupiedCellUpSampling_(NumericVector values, NumericVector loc){
	int n= loc.size();
	NumericVector endVal(n);
	for(int i=0;i<n;i++){

		endVal(i)=values(loc(i));

	}
	return endVal;

}



// [[Rcpp::export]]
NumericMatrix ExpandBoundariesToCols_(NumericMatrix f, NumericMatrix v, NumericVector res, NumericVector origin, int pent){
	int nF = f.rows();
	int nColF = f.cols();

	// supposed final final object
	int overPoints = max(res);

	// what the edge looks like, when resolution is the same for all faces
	NumericMatrix templateMatrix(overPoints,3);

	// and fill with NAs
	//	std::fill( templateMatrix.begin(), templateMatrix.end(), NumericVector::get_na() ) ;
	std::fill( templateMatrix.begin(), templateMatrix.end(), -80000 ) ;

	// this will be copied over to
	NumericMatrix nowUsed;

	NumericMatrix semiEndObj(nColF*(overPoints+4)*10,nF);

	NumericMatrix tempMat;

	int counter=0;
	int first=0;
	int second=0;

	bool onlyNew=1;

	int nRun;

	NumericVector ID(nF);

	for(int i=0; i<nF; i++){
		counter= 0;

		// the exception for the first 'pent' pentagons in a hexagrid
		if(nColF==6 && i<pent){
			nRun = nColF-1;
		}else{
			nRun = nColF;
		}

		//for every edge = no. of points
		for(int j=0; j<nRun;j++){
			//first vertex in the edge
			first=f(i,j);

			// second vertex in the edge
			if(j<(nRun-1)){
				second= f(i,j+1);

			}else{
				second=f(i,0);
			}

			//divide the edge
			tempMat=SplitArc_(v(first,_), v(second,_), origin, res(i), onlyNew);

			// now this changes with every run, but all we have to do is to fill it with NAs
			nowUsed = templateMatrix;
			std::fill( nowUsed.begin(), nowUsed.end(), -80000 ) ;

			//fill in the NA-matrix with actual values
			for(int u=0;u<tempMat.rows();u++){
				nowUsed(u,_) = tempMat(u,_);
			}

			// refer to the copied matrix, that should constrain the structure
			for(int k=0;k<3;k++){
				semiEndObj(counter,i)=v(first,k);
				counter++;
				for(int m=0;m<nowUsed.rows();m++){
					semiEndObj(counter,i) = nowUsed(m,k);
					counter++;
				}


			}


		}
		semiEndObj(counter,i) = v(second,0);
		counter++;
		semiEndObj(counter,i) = v(second,1);
		counter++;
		semiEndObj(counter,i) = v(second,2);
		counter++;

		ID(i) = i;
	}

	NumericMatrix endObj(counter+1,nF);

	for(int i=0;i<counter;i++){
		endObj(i,_) = semiEndObj(i,_);
	}
	endObj(counter,_) = ID;

	return endObj;
}



// cell shape change function (C++)
// will operate in 2d
// [[Rcpp::export]]
double ShapeTri_(NumericVector p0, NumericVector p1, NumericVector p2){
	NumericVector allDist(3);

	allDist(0) = dist(p0, p1);
	allDist(1)= dist(p1, p2);
	allDist(2) = dist(p2, p0);

	// sort the distances of the cells
	 NumericVector allDist2 = stl_sort(allDist);

//	 Rcout << allDist2 << "\n";

	double Shape;

	Shape = (allDist2(0)*allDist2(1))/pow(allDist2(2), 2.0);

	return Shape;
}

// [[Rcpp::export]]
NumericVector AllShapeTri_(NumericMatrix v, NumericMatrix f){
	int nF = f.rows();
	NumericVector endObj(nF);

	for(int i=0;i<nF; i++){
		endObj(i) = ShapeTri_(
			v(f(i,0),_),
			v(f(i,1),_),
			v(f(i,2),_)
		);

	}

	return endObj;
}

// [[Rcpp::export]]
NumericMatrix expandEdges_(NumericMatrix eExp, NumericVector center, double res){
	// the number of edges
	int nEdges = eExp.rows()/2;

	// allocate some memory to the final matrix
	NumericMatrix mTemp(nEdges*(res+2)*3+10,3);

	// count the number of entries in the temporary matrix
	int counter =0;
	NumericVector p0(3);
	NumericVector p1(3);

	NumericMatrix split;

	for(int i=0; i<nEdges;i++){
		// the start and end points of the new edge
		p0 = eExp(2*i,_);
		p1 = eExp(2*i+1,_);

		// Split the edges
		split = SplitArc_(p0,p1, center, res, 0);

		// for every startpoint of an edge
		for(int j=0;j<(split.rows()-1);j++){
			mTemp(counter,_) = split(j,_);
			counter++;
			mTemp(counter,_) = split(j+1,_);
			counter++;
		}

	}

	//copy over to the final data matrix
	NumericMatrix endObj(counter, 3);
	for(int i=0; i<counter;i++){
		endObj(i,_) = mTemp(i,_);
	}

	return endObj;

}


 // [[Rcpp::export]]
 NumericMatrix ExpandEdgesByFacesTri_(NumericMatrix origV, NumericMatrix origSubF, NumericVector center, int breaks){
 	int nFaces = origSubF.rows();

	// the temporary data holder of the total split
 	NumericMatrix tempObj(nFaces*(breaks+6)*3,4);

	// first point in the split
	NumericVector p0(3);

	// second point in the split
 	NumericVector p1(3);

	// the newly split edge
 	NumericMatrix tempEdge;

	// count the total number of points in the new table
 	int counter = 0;

 	// do this procedure for every single face
 	for(int i=0; i<nFaces;i++){
 		// for all 3 points
 		for(int j=0; j<3;j++){
 			p0 = origV(origSubF(i,j),_);

 			if(j==2){
 				p1 = origV(origSubF(i,0),_);
 			} else {
 				p1 = origV(origSubF(i,j+1),_);
 			}

 			// do the splitArc
 			tempEdge = SplitArc_(p0, p1, center, breaks, 1);

 			// save the results
 			// the first point
 			for(int l=0; l<3;l++){
				tempObj(counter,l) = p0(l);
			}
 			tempObj(counter,3) = i;
			counter++;

			// the new points
			for(int k=0;k<breaks;k++){
				for(int l=0;l<3;l++){
					tempObj(counter,l) = tempEdge(k,l);
				}
				tempObj(counter,3) = i;
				counter++;
			}

		}

 	}

 	// return only the meaningful results
 	NumericMatrix endObj(counter,4);

 	for(int i=0;i<counter; i++){
 	  for(int l=0;l<4;l++){
 		  endObj(i,l) = tempObj(i,l);
 	  }
 	}

 	return endObj;

 }


 // [[Rcpp::export]]
 NumericMatrix edgeListFromNeighbours_(NumericMatrix outN){
	// define the parameters
	int nCol=outN.cols();
	int nRow=outN.rows();

	// the edgelist will be:
	int nColNew = nRow*(nCol-1);
	NumericMatrix edgeList(nColNew,2);

	int counter =0;
	int actVal =-9;
	int fromVal=-9;

	// fill in this edgelist
	for(int i=0; i<nRow; i++){
		for(int j=0; j<4; j++){
			actVal= outN(i,j);
			fromVal=i+1;
			if(fromVal!=actVal){
				if(fromVal<actVal){
					edgeList(counter, 0) = fromVal;
					edgeList(counter, 1) = actVal;
				}else{
					edgeList(counter, 1) = fromVal;
					edgeList(counter, 0) = actVal;
				}

				counter++;

			}

		}

	}

	return edgeList;


 }

 // [[Rcpp::export]]
NumericVector EvenInterpolation_(NumericMatrix xyz, NumericVector origin, double critValue){
	// 1. calculate the distances between the points
	int nRow = xyz.rows();
	NumericVector distances(nRow);
	int p0;
	int p1;
	for(int i=0; i<(nRow);i++){
		p0 = i;
		p1= i+1;
		if(p0==(nRow-1)){
			p1 = 0;
		}
		distances(i) = ArcDist_(xyz(p0,_), xyz(p1,_), origin, 0);
	}

	// 2. estimate the number of points necessary based on this
	int nEstimate =0;
	double q=-1.0;
	int q1=-1;
	for(int i=0; i<distances.size();i++){
		//if the distance is larger then the critical value
		if(distances(i)>critValue){
			//divide to a number of parts
				q = distances(i)/critValue;
				q1 = (int)round(q);
				nEstimate = nEstimate+q1;
		} else {
		//store
			nEstimate++;
		}

	}

	// 3. get the new points
	NumericMatrix tempObj(nEstimate*2,3);
	tempObj(0,_) = xyz(0,_);

	NumericMatrix tempLine;

	int counter =1;
	for(int i=0;i<nRow;i++){
		p0 = i;
		p1= i+1;
		if(p0==(nRow-1)){
			p1 = 0;
		}
		if(distances(i)>critValue){
			q = distances(i)/critValue;
			q1 = (int)round(q)-1;

			// the new points
			tempLine = SplitArc_(xyz(p0,_), xyz(p1,_), origin, q1, 1);

			// just a regular copy of the new points
			for(int j=0;j<q1;j++){
				tempObj(counter,_) = tempLine(j,_);
				counter++;
			}

			// and the final point
			tempObj(counter,_) = xyz(p1,_);
			counter++;

		} else {
		//store
			tempObj(counter,_) = xyz(p1,_);
			counter++;
		}

	}

	// 4. clean up
	NumericMatrix endObj(counter, 3);
	for(int i=0;i<counter;i++){
		endObj(i,_) = tempObj(i,_);

	}

	return endObj;

}


// [[Rcpp::export]]
NumericMatrix centroidPoints_(NumericMatrix coords, NumericMatrix combin, NumericVector center, int breaks){
	// the number of combinations!
	int nComb = combin.rows();

	// estimate the size of the resulting object
	NumericMatrix tempObj(nComb*(breaks+3),3);

	int counter=0;

	NumericMatrix tempLine;
	NumericVector p0(3);
	NumericVector p1(3);

	for(int i=0;i<nComb;i++){
		p0= coords(combin(i,0),_);
		p1= coords(combin(i,1),_);

		tempLine = SplitArc_(p0, p1, center, breaks, 0);
		for(int j=0;j<(breaks+2);j++){
			tempObj(counter,_) = tempLine(j,_);
			counter++;
		}


	}

	// copy over
	NumericMatrix endObj(counter, 3);
	for(int i=0;i<counter;i++){
		endObj(i,_) = tempObj(i,_);
	}

	return endObj;
}

// [[Rcpp::export]]
NumericMatrix projectCloseToPoint_(NumericMatrix coords, NumericVector toPoint, NumericVector center, int breaks){
	int nPoints = coords.rows();
	NumericMatrix newPoints(nPoints, 3);

	for(int i=0;i<nPoints;i++){
		newPoints(i,_)= SplitArc_(toPoint, coords(i,_) , center, breaks, 1)(0,_);

	}

	return newPoints;

}
