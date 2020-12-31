#include <math.h>

// function toi calculate the determinant of a 4*4 matrix

extern "C" void Determinant4x4( double *v0,double *v1,double *v2,double *v3 , double *det){
 *det =
	v0[3]*v1[2]*v2[1]*v3[0] - v0[2]*v1[3]*v2[1]*v3[0] -
    v0[3]*v1[1]*v2[2]*v3[0] + v0[1]*v1[3]*v2[2]*v3[0] +

    v0[2]*v1[1]*v2[3]*v3[0] - v0[1]*v1[2]*v2[3]*v3[0] -
    v0[3]*v1[2]*v2[0]*v3[1] + v0[2]*v1[3]*v2[0]*v3[1] +

    v0[3]*v1[0]*v2[2]*v3[1] - v0[0]*v1[3]*v2[2]*v3[1] -
    v0[2]*v1[0]*v2[3]*v3[1] + v0[0]*v1[2]*v2[3]*v3[1] +

    v0[3]*v1[1]*v2[0]*v3[2] - v0[1]*v1[3]*v2[0]*v3[2] -
    v0[3]*v1[0]*v2[1]*v3[2] + v0[0]*v1[3]*v2[1]*v3[2] +

    v0[1]*v1[0]*v2[3]*v3[2] - v0[0]*v1[1]*v2[3]*v3[2] -
    v0[2]*v1[1]*v2[0]*v3[3] + v0[1]*v1[2]*v2[0]*v3[3] +

    v0[2]*v1[0]*v2[1]*v3[3] - v0[0]*v1[2]*v2[1]*v3[3] -
    v0[1]*v1[0]*v2[2]*v3[3] + v0[0]*v1[1]*v2[2]*v3[3];

}

// function to determine whether a point is inside the tetrahedron
// vertices and query objects shoudl be xyz1 coordinates.
// dets array of 5 double values-rewritten
// based on:

extern "C" void PointInTetrahedron_(double *vertices,double *query,double *dets,int *result){


	double *v0 = vertices;
	double *v1 = vertices+4;
	double *v2 = vertices+8;
	double *v3 = vertices+12;


    Determinant4x4(v0, v1, v2, v3, dets);
    Determinant4x4(query, v1, v2, v3, dets+1);
    Determinant4x4(v0, query, v2, v3, dets+2);
    Determinant4x4(v0, v1, query, v3, dets+3);
    Determinant4x4(v0, v1, v2, query, dets+4);
    /**
    If by chance the Determinant det0 is 0, then your tetrahedron is degenerate (the points are coplanar).
    If any other Di=0, then P lies on boundary i (boundary i being that boundary formed by the three points other than Vi).
    If the sign of any Di differs from that of D0 then P is outside boundary i.
    If the sign of any Di equals that of D0 then P is inside boundary i.
    If P is inside all 4 boundaries, then it is inside the tetrahedron.
    As a check, it must be that D0 = D1+D2+D3+D4.
    The pattern here should be clear; the computations can be extended to simplicies of any dimension. (The 2D and 3D case are the triangle and the tetrahedron).
    If it is meaningful to you, the quantities bi = Di/D0 are the usual barycentric coordinates.
    Comparing signs of Di and D0 is only a check that P and Vi are on the same side of boundary i.
    */
	if(dets[0]==0){
		// no tetrahedron - bad!
		*result = -9;
	}else{
		//default: the point is not in the tetrahedron
		*result = -1;

		//the first determinant is negative
		if (dets[0] < 0)
        {
            if ((dets[1] < 0) && (dets[2] < 0) && (dets[3] < 0) && (dets[4] < 0))
            {
              *result = 1;
            }
        }
		//the first determinant is positive
        if (dets[0] > 0)
        {
            if ((dets[1] > 0) && (dets[2] > 0) && (dets[3] > 0) && (dets[4] > 0))
            {
                *result = 1;
            }
        }

		//more: the point lies on the boundary planes of the tetrahedron
		// unfortunately this is not enough, because other face-tetrahedra
		// are bounded by these planes - additional tests are needed
		// for the central triangles

	}

}

// hierarchical searching
extern "C" void _locateTriangle_(
	double *allVertices,
	int *divs,
	int *nDivs,
	double *queries,
	int *nQrs,
	int *queryIndex,
	int *faceIndex,
	int *offset,
	int *faceContainer,
	int *foundMiddle,
	int *tempF
	){
	// function starts from here

	// index for the last
	int resIndex=0;

	// calculate the offset of degrees in the faces table
//	int *offset= new int[*nDivs+1];

	// it starts with one
	offset[0] = 0;
	int faceNo=1;

	for(int d=0;d<*nDivs;d++){
		faceNo = divs[d]*faceNo;
		offset[d+1]= offset[d]+faceNo;
	}



	// declaration of arrays of indices
	// start of the active degree
	double*divStart;
	// start of the subfacegroup belonging to the
	double*divFaceStart;
	// start of the face-tetrahedron
	double*verts;
	// declaration of data point arrays
	double *actQuery;


	//the temporary determinants array - for the tetrahedron check - will be rewritten in every loop
	double *dets = new double[5];

	// single points can end up being on multiple faces -  this needs to be taken care of
	int nFoundPoints;
	int counterMid;
	int nFoundCounter;

	// 6 - theoretical maximum number of faces a point can sit on
	int*foundPoints = new int[12];
//	int*foundMiddle = new int[12];
	int*foundPointsInner = new int[12];


	// iterate lookup for every queried (q) point
	for(int q=0;q<*nQrs; q++){
//-	int q=0;
		foundPoints[0] = 0;
		nFoundPoints=1;

		// step over the queried points in every loop
		actQuery = queries+q*4;

		// iterate for every degree (d)
		for(int d=0; d<*nDivs;d++){
//-		int d=0;
			// shift the start to the start of the degree
			divStart = allVertices+(offset[d])*16;


			//empty container to store the results of the lookups
			//overwritten in every search

//				int*faceContainer = new int[divs[d]];

			counterMid=0;

			// iterate for every results (where 0 and 1 is present)
			for(int r=0;r<nFoundPoints;r++){
//-				int r=0;
				// shift the start to the first subface of the face,
				// where the point was found in the previous degree
				divFaceStart= divStart+(foundPoints[r]*divs[d])*16;

	//-			divFaceStart=allVertices;

				// in which face is the point found
				nFoundCounter=0;

				// iterate for every relevant face (f) - determined by d
				for(int f=0; f<divs[d];f++){
		//int f=0;
					verts= divFaceStart+f*16;

					PointInTetrahedron_(verts, actQuery, dets, (faceContainer+f));

					// check the container immediately for results
					if(*(faceContainer+f)==1){
						foundPointsInner[nFoundCounter] = f+foundPoints[r]*divs[d];
						nFoundCounter++;

					}
					
					//tempF
					if(d==0 && q==8){
						tempF[f] = faceContainer[f];
						
					}

				}


				// store the faces in a common container for this
				// degree
				for(int o=0;o<nFoundCounter;o++){
					foundMiddle[counterMid] = foundPointsInner[o];
					counterMid++;
				}


			}

			//end of looping - if in the last division
			if(d==((*nDivs)-1)){
				for(int r2=0;r2<counterMid; r2++){

					queryIndex[resIndex] = q;
					faceIndex[resIndex] = foundMiddle[r2];
					resIndex++;
				}

			}


			// prepare for the next loop
			// do checking for all the faces found previously
			for(int p=0;p<counterMid;p++){
				foundPoints[p] = foundMiddle[p];
			}

			nFoundPoints = counterMid;

//			delete faceContainer;

		}

	}
	// clean up

	delete[] foundPoints;
//	delete[] foundMiddle;
	delete[] foundPointsInner;
	delete[] dets;

}

