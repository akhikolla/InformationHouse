# S4 class definitions
# the core class! -  no need to export
obj3d <- setClass(
	#class
	"obj3d",
	
	slots=c(
		vertices = "matrix",
		faces = "matrix",
		edges = "matrix",
		center = "numeric",
		length = "numeric"
	)
)

#the icosahedron
icosahedron <- setClass(
	"icosahedron",
	
	slots = c(
		edgeLength = "numeric",
		skeleton = "list",
		r = "numeric",
		sp="ANY"
	),
	
	contain="obj3d"
)

#constructor of icosahedron
setMethod(
	"initialize",
	signature = "icosahedron",
	definition = function (.Object, r=FALSE,a=FALSE){
		# calculate the missing thing
		if(missing(a))
		{
			.Object@r<-r
			.Object@edgeLength<-r/(sin(2*pi/5))
		}
		if(missing(r))
		{
			.Object@edgeLength<-a
			.Object@r<-a*(sin(2*pi/5))
		}
		
		phi<-0.5*(1+sqrt(5))
		
		P1<-c(0,1,phi)/2
		P2<-c(0,1,-phi)/2
		P3<-c(0,-1,phi)/2
		P4<-c(0,-1,-phi)/2
		
		P5<-c(1,phi,0)/2
		P6<-c(1,-phi,0)/2
		P7<-c(-1,phi,0)/2
		P8<-c(-1,-phi,0)/2
		
		P9<-c(phi,0,1)/2
		P10<-c(phi,0,-1)/2
		P11<-c(-phi,0,1)/2
		P12<-c(-phi,0,-1)/2
		
		vertices<-rbind(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12)*.Object@edgeLength
		
		#maximum precision achieved by rotation
		rotVal<-0.55357435889704526002
		vertices<-t(apply(vertices, 1, rotateOnePoint, angles=c(rotVal,0,0), origin=c(0,0,0)))

		colnames(vertices)<-c("x","y","z")
		
		#the edges
		edges<-NULL
		for(i in 1:nrow(vertices))
		{
			for (j in 1:nrow(vertices))
			{
				A<-unlist(vertices[i,])
				B<-unlist(vertices[j,])
				mTwo<-matrix(c(A,B),ncol=3, byrow=TRUE)
				d<-stats::dist(mTwo)
				if(round(d-.Object@edgeLength,10)==0)
				{
					E<-c(rownames(vertices)[i], rownames(vertices)[j])
					E<-sort(E)
					edges<-rbind(edges, E)
				}
			
			
			
			}
			
		}
		edges<-unique(edges)
		rownames(edges)<-paste("E", 1:nrow(edges), sep=(""))
		
		faces<-NULL
		#the faces
		for(i in 1:nrow(edges))
		{
			actEdge<-unlist(edges[i,])
			
			#P1
			b1<-edges[,1]%in%actEdge[1] | edges[,2]%in%actEdge[1]
			
			u1<-unique(c(edges[b1,1], edges[b1,2]))
			
			#P3
			b2<-edges[,1]%in%actEdge[2] | edges[,2]%in%actEdge[2]
			
			u2<-unique(c(edges[b2,1], edges[b2,2]))
			#double faces
				chDoub<-u1[u1%in%u2]
			#new points
				chNew<-chDoub[!chDoub%in%actEdge]
			
			#new lines in the faces
			#the first new point
			F<-sort(c(actEdge, chNew[1]))
			faces<-rbind(faces,F)
			
			#the second new point
			F<-sort(c(actEdge, chNew[2]))
			faces<-rbind(faces,F)
			
				
		}
		faces<-unique(faces)
		rownames(faces)<-paste("F",1:nrow(faces), sep="")
		
		.Object@faces<-faces
		.Object@edges<-edges
		.Object@vertices<-vertices
		.Object@center <- c(0,0,0)
		.Object@length<- c(
			"vertices"=nrow(vertices),
			"edges"=nrow(edges),
			"faces"=nrow(faces))
			
		# the skeleton:
			#vertices
			v<-vertices
			rownames(v)<-NULL
			
			#faces
			f<-matrix(NA, nrow=20,ncol=5)
			fTemp<-as.numeric(unlist(lapply(strsplit(as.character(faces),"P"), function(x){x[2]})))
			f[,1:3]<-matrix(fTemp, ncol=3)-1
			f[,4] <- rep(-1,nrow(f))
			f[,5] <- rep(-1,nrow(f))
		
	
		.Object@skeleton<-list(v=v, f=f)
		
		return(.Object)
	}
)


#' A triangular icosahedral grid
#' 
#' \code{trigrid()} creates a triangular grid based on the
#'    tessellation of an icosahedron.
#'
#' The grid structure functions as a frame for data graining, plotting and spatial
#'	calculations. Data can be stored in layers that are linked to the grid object. In the current version only the 
#'	\code{facelayer} class is implemented, which allows the user to render data to the cells
#'	of the grid, which are usually referred to as faces. 
#'	The grid 'user interface' is made up of four primary tables: the \code{@vertices} table for the coordinates of the vertices,
#' 	the \code{faceCenters} for the coordinates of the centers of faces,
#'	the \code{faces} and the \code{edges} tables that contain which vertices form which faces and edges respectively.
#'	In these tables, the faces and vertices are sorted to form spirals that go from the north pole in a counter-clockwise
#'	direction. In case grid subsetting is performed these tables get truncated.
#'	
#'	At finer resolutions, the large number of spatial elements render all calculations resource demanding and slow, 
#'	therefore the hierarchical structure created during the tessellation procedure is retained for efficient implementation.
#'	These data are stored in a list in the slot \code{@skeleton} and are 0-indexed integer tables for Rccp-based functions. \code{$v} 
#'	stores vertex, \code{$f} the edge, and \code{$e} contains the edge data for plotting and calculations. In these tables
#'	the original hierarchy based orderings of the units are retained, during subsetting, additional vectors are used to indicate
#'	deactivation of these units. Any sort of meddling with the \code{@skeleton} object will lead to unexpected behavior.
#'	
#' @slot vertices Matrix of the vertex XYZ coordinates.
#'
#' @slot faces Matrix of the verticies forming the faces.
#'
#' @slot edges Matrix of the vertices forming the edges.
#'
#'	@slot tessellation Contains the tessellation vector.
#'
#'	@slot orientation Contains the grid orientation in xyz 3d space, values in radian relative to the (0,1,0) direction.
#'
#'	@slot center is the xyz coordinates of the grids origin/center.
#'
#'	@slot div vector contains the number of faces that a single face of the previous tessellation level is decomposed to.
#'
#'	@slot faceCenters contains the xyz coordinates of the centers of the faces on the surface of the sphere.	
#'	@slot belts Vector of integers indicating the belt the face belongs to.
#' @slot edgeLength the length of an average edge in km and degrees.
#' @slot graph an 'igraph' class graph object.
#' @slot length integer vector of length=3. The number of vertices, edges and faces in this order.
#' @slot proj4string a CRS class object indicating the model in the PROJ.4 system
#' @slot r the radius of the grid
#' @slot sp The SpatialPolygons representation of the grid. If missing, it can be created with newsp().
#' @slot skeleton data tables with sequential indexing for the C functions.
#'
#'
#' @param tessellation (\code{numeric}) An integer vector with the tessellation values. Each number
#'    describes the number of new edges replacing one original edge. Multiple series of tessellations
#' 	  are possible this way. The total tessellation is the product of the tessellation vector. 
#'	  Higher values result in more uniform cell sizes, but the larger number of tessellation series
#'	  increases the speed of lookup functions.
#'
#' @param sp (\code{logical}) Flag indicating whether the \code{\link[sp]{SpatialPolygons}} class representation of the grid
#'	should be added to the object when the grid is calculated. If set to \code{TRUE} the \code{SpPolygons()} function will be run with with the resolution parameter set to 25. The 
#'  resulting object will be stored in slot \code{@sp}. As the calculation of this object can substantially increase the grid creation time,
#'	 by default this argument has a value of \code{FALSE}. The \code{\link[sp]{SpatialPolygons}} class representation can be added on demand by running the function \code{newsp}.
#'
#' @param graph (\code{logical}) Flag indicating whether the \code{'igraph'} class representation of the grid
#'	should be added to the object when the grid is calculated. This argument defaults to \code{TRUE} because this option has only minor performance load on the grid 
#' constructor function. For familiarization with the
#' object structure, however, setting this parameter to \code{FALSE} might help, as invoking \code{\link[utils]{str}} on the \code{'igraph'} class slot of the class might flood the console.
#'
#' @param radius (\code{numeric}) The radius of the grid. Defaults to the authalic radius of Earth.
#' @param center (\code{numeric}) The origin of the grid in the reference Cartesian coordinate system. Defaults to (0,0,0).
#'
#' @return A triangular grid object, with class \code{trigrid}.
#' @examples
#' # single tessellation value
#' g <- trigrid(c(8))
#' g
#' # series of tessellations
#' g1 <- trigrid(c(2,3,4))
#' g1
#' @exportClass trigrid
trigrid<-setClass(
	"trigrid",
	contain="icosahedron",
	slots=c(
		tessellation="numeric",
		div="numeric",
		faceCenters="matrix",
		orientation="numeric",
		belts="numeric",
		proj4string="CRS",
		graph="ANY"
	)
)


#' @export trigrid
setMethod(
	"initialize",
	signature="trigrid",
	definition=function(.Object, tessellation=1, sp=FALSE, graph=TRUE, radius=authRadius, center=origin){
		
		# trial variables
	#	tessellation<-c(2,2,2,2)
	#	authRadius<-6371
	#	degree<-1
	#	fa<-4
		##check the tessellation vector
		if(!sum(tessellation%%1)==0 | !is.numeric(tessellation))
		{
			stop("Invalid tessellation vector. Please enter a vector of integers only.")
		}
		if(prod(tessellation)==0){
			stop("Invalid tessellation vector. Please enter a vector of integers without 0s.")
		}
		

		# create a basic icosahedron
			icosa<-icosahedron(r=radius)
		
			# the final trigrid object
			.Object@tessellation <-tessellation
			.Object@r <- radius
			.Object@center <- icosa@center
			
			# add the CRS for spatial transformations
			#supress scientific notation
			options(scipen=999)
			.Object@proj4string <- sp::CRS(paste("+proj=longlat +a=", round(radius*1000), " +b=", round(radius*1000), sep=""))
			
			
		# extract the skeleton	
			f=icosa@skeleton$f
			v=icosa@skeleton$v
		
		# in case a tessellation is happening
		if(prod(tessellation)>1){
			# invoke the c++ function to create a grid
			newGrid<-.Call(Cpp_icosa_IcosahedronTesselation_, 
				v,
				f, 
				tessellation, 
				icosa@center)
		}else{
			newGrid <-list(f=f, v=v)
		}
			
		
		# additional data
			# the divisions
			div <- c(20, tessellation^2)
			
			# in case of the icosahedron
			if(prod(tessellation)==1){
				div<-20
			}
			
			.Object@div <- div
	
		# calculate the edges
			# boolean variables for the subsetting
			if(prod(tessellation)>1){
				aF <- newGrid$f[,4]==(length(tessellation)-1)
			# for special case of the icosahedron
			}else{
				aF <- newGrid$f[,4]==-1
			}
			offSet<-min(which(aF))-1
			
			aV <- rep(T, nrow(newGrid$v))
		
			faces<-subset(newGrid$f, aF)[,1:3]
			# edges
			e<-unique(.Call(Cpp_icosa_expandFacesToEdges_, faces))
			aE <- rep(T, nrow(e))
			# the neighbours of the faces
			n<-.Call(Cpp_icosa_AllNeighboursTri_, newGrid$f[,1:3], div)
			
		
		# the R grid UI
		# vertices
			vertices<- newGrid$v
			colnames(vertices) <- c("x","y","z")
			
		
		# the centers of the triangles
			faceCenters<-.Call(Cpp_icosa_allTriangleCenters_, newGrid$v, faces, icosa@center)
			colnames(faceCenters)<-c("x","y","z")
			
	
		#ordering the face and vertex data
			#vertex at the very north
				topVert<- which(max(vertices[,3])==vertices[,3])-1
			
			#top five faces
				firstFiveFace<-order(faceCenters[,3], decreasing=TRUE)[1:5]
				#ordered by longitude
				longOrd<-order(CarToPol(faceCenters[firstFiveFace,], norad=TRUE, origin=icosa@center)[,"long"])
				startFaces<-firstFiveFace[longOrd]
			
			#starting vertices
				startF<-faces[startFaces,]
				
				startVert<-numeric(6)
				startVert[1]<-topVert
				startVert[2]<-startF[1,c(F,startF[1,2:3]%in%startF[5,2:3])]
				
				for(i in 1:4){
					startVert[i+2]<-startF[i,!startF[i,]%in%c(startVert)]
				}
				
			# arguments for the ordering function	
				nBelts<-prod(tessellation)*3
				# in case the  thing is an icosahedron
				
				nV<-nrow(vertices)
				startFaces<-startFaces-1
			
			#call
				ordering<-.Call(Cpp_icosa_orderTriGrid_, faces, n, startFaces, startVert, nBelts, nV)
			
			# the belts
				# where do the belts start
				beltStartsAndEnd<-c(ordering$belts+1, nrow(faces))
				
				#empty container
				belts<-rep(NA, nrow(faces))
				
				#for every belt
				for(i in 1:(length(beltStartsAndEnd)-1)){
					
					if(i<(length(beltStartsAndEnd)-1)){
						actInd<-beltStartsAndEnd[i]:(beltStartsAndEnd[i+1]-1)
					}else{
						actInd<-beltStartsAndEnd[i]:(beltStartsAndEnd[i+1])
					}
					#store
					belts[actInd]<-i
				}

				.Object@belts<-belts
			#output formatting (R indices)
				#UI-skeleton translation
				faceOrder<-ordering$faceOrder+1 # good
				names(faceOrder)<-paste("F", 1:nrow(faces), sep="") #good
				
				vertexOrder<-ordering$vertexOrder+1
				names(vertexOrder)<-paste("P", 1:nrow(vertices), sep="") #good
				
				#skeleton-UI translation (R indicies)
				faceInvertOrder<-rep(0, length(faceOrder))
				faceInvertOrder[faceOrder]<-1:length(faceOrder)
				
				vertexInvertOrder<-rep(0, length(vertexOrder))
				vertexInvertOrder[vertexOrder]<-1:length(vertexOrder)
				
				
				aF[aF]<-faceInvertOrder #not good
				aV[aV]<-vertexInvertOrder#not good
				
			
			# the skeleton
				.Object@skeleton <- list(v=newGrid$v,aV=aV, uiV=vertexOrder, e=e, aE=aE,  f=newGrid$f, aF=aF, uiF=faceOrder, n=n, offsetF=offSet)
#				.Object@skeleton <- list(v=newGrid$v,aV=aV, e=e, aE=aE,  f=newGrid$f, aF=aF)
			
			#using the output to order
			vertices<-vertices[vertexOrder,]
			rownames(vertices) <- paste("P",1:nrow(vertices) ,sep="")
		
			# the edges (outer ordering)
			tempE<-aV[as.numeric(e+1)]
			
			edges<-matrix(paste("P",tempE, sep=""), ncol=2)
			rownames(edges)<-paste("E", 1:nrow(edges), sep="")
			
			# the faces
			faces<-faces[faceOrder,]
			#translate the vertex information!
			facesNum<-as.numeric(faces)+1
			faces2<-aV[facesNum]
			
			faces<-matrix(paste("P",faces2, sep=""), ncol=3)
			rownames(faces)<-paste("F", 1:nrow(faces), sep="")
		
			#the faceCenters
			faceCenters<-faceCenters[faceOrder,]
			rownames(faceCenters)<-paste("F", 1:nrow(faceCenters), sep="")
			
			options(scipen=0)
			
			#add to the object
			.Object@faces <- faces
			.Object@vertices <- vertices
			.Object@edges <- edges
			.Object@faceCenters <- faceCenters
		
		#length attribute
		.Object@length<- c(
			"vertices"=nrow(.Object@vertices),
			"edges"=nrow(.Object@edges),
			"faces"=nrow(.Object@faces))
			
		.Object@edgeLength <- c(
			mean(.Call(Cpp_icosa_edges_, newGrid$v, e, icosa@center, 1)),
			mean(.Call(Cpp_icosa_edges_, newGrid$v, e, icosa@center, 0))/pi*180)
			
		names(.Object@edgeLength) <- c("km", "deg")
		
		.Object@orientation<-c(0,0,0)
		
		# add the igraph of the grid
		dummy<-NA
		.Object@graph<-dummy
		
		if(graph==TRUE){
				.Object@graph<-gridgraph(.Object)
		}
		
		#2d grid!
		if(sp==TRUE){
				.Object@sp<-SpPolygons(.Object, res=25)
		}else{
			dummy<-NA
			
			.Object@sp<-dummy
		}
		# correct the coordinates with the center data
		for(vc in 1:3){
			.Object@vertices[,vc]<-.Object@vertices[,vc]+center[vc]
			.Object@faceCenters[,vc]<-.Object@faceCenters[,vc]+center[vc]
			.Object@skeleton$v[,vc]<-.Object@skeleton$v[,vc]+center[vc]
		}
		.Object@center <- center
		
		return(.Object)
	
	}
)

	
#' Construct a penta-hexagonal icosahedral grid
#' 
#' The \code{hexagrid} function constrcucts a hexa-pentagonal grid based on the inversion of a 
#'    tessellated icosahedron.
#'
#' Inherits from the \code{trigrid} class.
#'
#' The grid structure functions as a frame for data graining, plotting and
#'	calculations. Data can be stored in layers that are linked to the grid object. In the current version only the 
#'	\code{\link{facelayer}} class is implemented which allows the user to render data to the cells
#'	of the grid which are called faces. 
#' 	The grid 'user interface' is made up of four primary tables: the \code{@vertices} table for the coordinates of the vertices,
#' 	the \code{faceCenters} for the coordinates of the centers of faces,
#'	the \code{faces} and the \code{edges} tables that contain which vertices form which faces and edges respectively.
#'	In these tables, the faces and vertices are sorted to form spirals that go from the north pole in a counter-clockwise
#'	direction. In case grid subsetting is performed these tables get truncated.
#'	
#'	At finer resolutions, the large number of spatial elements render all calculations very resource demanding and slow, 
#'	therefore the hierarchical structure created during the tessellation procedure is retained for efficient implementations.
#'	These data are stored in a list in the slot \code{@skeleton} and are 0-indexed integer tables for Rccp-based functions. \code{$v} 
#'	stores vertex, \code{$f} the edge, and \code{$e} contains the edge data for plotting and calculations. In these tables
#'	the original hierarchy based orderings of the units are retained, during subsetting, additional vectors are used to indicate
#'	deactivation of these units. Any sort of meddling with the @skeleton object will lead to unexpected behavior.
#'	
#' @slot vertices Matrix of the vertex coordinates.
#'
#' @slot faces Matrix of the verticies forming the faces
#'
#' @slot edges Matrix of the vertices forming the edges.
#'
#'	@slot tessellation Contains the tessellation vector.
#'
#'	@slot orientation Contains the grid orientation in xyz 3d space, values in radian.
#'
#'	@slot center The xyz coordinates of the grid's origin/center.
#'
#'	@slot div Contains the number of faces that a single face of the previous tessellation level is decomposed to.
#'
#'	@slot faceCenters Contains the xyz coordinates of the centers of the faces on the surface of the sphere.	
#'
#'
#' @param tessellation (\code{numeric}) An integer vector with the tessellation values. Each number
#'    describes the number of new edges replacing one original edge. Multiple series of tessellations
#' 	  are possible this way. The total tessellation is the product of the tessellation vector. 
#'	  Higher values result in more uniform cell sizes, but the larger number of tessellation series,
#'	  increases the speed of lookup functions.
#'
#' @param sp (\code{logical}) Flag indicating whether the \code{\link[sp]{SpatialPolygons}} class representation of the grid
#'	should be added to the object when the grid is calculated. If set to true the \code{\link{SpPolygons}} function will be run with with the resolution parameter set to \code{25}. The 
#'  resulting object will be stored in slot \code{@sp}. As the calculation of this object can increase the grid creation time substantially
#'	 by default this argument has a value \code{FALSE}. This can be added on demand by running the function \code{\link{newsp}}.
#'
#' @param graph (\code{logical}) Flag indicating whether the \code{\link[igraph:aaa-igraph-package]{igraph}} class representation of the grid
#'	should be added to the object when the grid is calculated. This argument defaults to \code{TRUE} because this option has only minor performance load on the grid 
#' constructor function. For familiarization with the
#' object structure, however, setting this parameter to \code{FALSE} might help, as invoking \code{\link[utils]{str}} on the 'igraph' class slot of the class might flood the console.
#'
#' @param radius (\code{numeric}) The radius of the grid. Defaults to the authalic radius of Earth.
#' @param center (\code{numeric}) The origin of the grid in the reference Cartesian coordinate system. Defaults to \code{c(0,0,0)}.
#'
#'
#' @return A hexagonal grid object, with class \code{hexagrid}.
#' @examples
#' g <- hexagrid(c(8), sp=TRUE)
#' g1 <- hexagrid(c(2,3,4))
#' @exportClass hexagrid
hexagrid<-setClass(
	"hexagrid",
	contain="trigrid",
)



#' @export hexagrid
setMethod(
	"initialize",
	signature="hexagrid",
	definition=function(.Object, tessellation=1, sp=FALSE, graph=TRUE, center=origin, radius=authRadius){
			
		tGrid<-trigrid(tessellation, radius=radius)
		# v part of the skeleton and the active vertices
			# new $v: first part original vertices of the trigrid (faceCenters)
			# of the hexagrid
			# second par: original face centers of the trigrid (in the order of n and f)
			# vertices of the new hexagrid
			v<-rbind(tGrid@skeleton$v,tGrid@faceCenters[tGrid@skeleton$aF[as.logical(tGrid@skeleton$aF)],])
			rownames(v)<-NULL
			aV<-rep(FALSE, nrow(v))
			
			aV[(nrow(tGrid@skeleton$v)+1):length(aV)]<- tGrid@skeleton$aF[as.logical(tGrid@skeleton$aF)]
			
		# the number of original vertices
			nOrigV<-nrow(tGrid@skeleton$v)
			
			#vertices member (ordered in the trigrid)
			vertices<-tGrid@faceCenters
			options(scipen=999)
			rownames(vertices)<-paste("P", 1:nrow(vertices), sep="")
			
			#uiV - the R indexes of the vertices in $v
			uiV<-rep(0, nrow(vertices))
			uiV[aV[as.logical(aV)]]<-1:length(uiV)+nOrigV
			names(uiV)<-rownames(vertices)
			
		#the faces
			#triangular faces in the last tessellation
			fTemp<-tGrid@skeleton$f[(tGrid@skeleton$offsetF+1):nrow(tGrid@skeleton$f),1:3]
			
			#table of the subfaces:
			#first column: original vertices of trigrid, the internal order of hexagonal faces
			sF<-.Call(Cpp_icosa_CreateHexaSubfaces_, tGrid@skeleton$n, fTemp, nOrigV)
			
			#total faces table required for lookups
			f<-rbind(tGrid@skeleton$f, sF)
			
			# links the subfaces to their UI orders (0-no subface, 1: F1 in the ui)
			aSF<-rep(0, nrow(f))
			aSF[(nrow(tGrid@skeleton$f)+1):nrow(f)]<-tGrid@skeleton$aV[sF[,1]+1]
			
		#edges 
			#internal representation of edges - plotting
			e<-unique(sF[,2:3])
			
			#outer
			# the vertices with outer indices
			edges<-matrix(paste("P", aV[as.numeric(e)+1], sep=""), ncol=2)
			rownames(edges)<-paste("E", 1:nrow(edges), sep="")
			
			#activation - for subsetting
			aE<-rep(T, nrow(e))
			
		#faces matrix
			#the first column implies the face the subfaces belong to
			fOrdered<-unique(sF[order(sF[,1]),1:3])
			
			#which subfaces (R indices of rows in $f table) belong to which face
			#ordered like the internal representation of trigrids vertices!
			facesInternal<-.Call(Cpp_icosa_HexaFaces_, fOrdered)
			facesInternal[1:12,6]<-NA
			vF<-facesInternal
			
			#replace the vertex indices to the outer indices (names)
			facesExpandOut<-aV[as.numeric(facesInternal)+1]
			facesExpandOut[!is.na(facesExpandOut)]<-paste("P", facesExpandOut[!is.na(facesExpandOut)], sep="")
			
			#create a matrix
			facesExpandOut<-matrix(facesExpandOut, ncol=6)
			
			#reorder the faces from north to south
			faces<-facesExpandOut
			faces[tGrid@skeleton$aV,]<-facesExpandOut
			
			# suppress scientific notation
			rownames(faces)<-paste("F", 1:nrow(faces), sep="")
			
			#uiF - will be used for subsetting
			indices<-tGrid@skeleton$aV[sF[,1]+1]
			
			uiF<-.Call(Cpp_icosa_RetrieveIndexMat_, indices)
			uiF[uiF==0]<-NA
			
			#add the offset (R indices)
			uiF<-uiF+nrow(tGrid@skeleton$f)
			rownames(uiF)<-paste("F", 1:nrow(uiF), sep="")
			
			# when doing subsetting:
			# 1. select the faces by names
			# 2. flatten with as.numeric()
			# 3. subset $f for rows
			# 4. unique the output -> than you can use the indices for $v
			
		#face centers
			faceCenters<-tGrid@vertices
			rownames(faceCenters) <- paste("F", 1:nrow(faceCenters), sep="")
			
			# the centers of the faces for plotting ($plotV)
			#ordered as tGrird@skeleton$v
			transFace<-t(faces)
			flatF<-as.character(transFace)
			nas<-is.na(flatF)
			flatF[nas]<-"P1"
			
			vertsF<-vertices[flatF,]
			vertsF[nas,]<-NA
			rownames(vertsF)<-NULL
			indF<-rep(1:nrow(faces),each=6)
			x<-tapply(vertsF[,1], indF, mean, na.rm=TRUE)
			y<-tapply(vertsF[,2], indF, mean, na.rm=TRUE)
			z<-tapply(vertsF[,3], indF, mean, na.rm=TRUE)
			
			#still the outer order
			fcPlot<-cbind(x,y,z)
			
			#the inner order (which is necessary for plotting)
			plotV<-v
			plotV[tGrid@skeleton$uiV,]<-fcPlot[,]
			
			
			#the inner order
			
			# total skeleton
			aF<-tGrid@skeleton$aV
			
			skeleton<-list(f=f, vF=vF, aF=aF,aSF=aSF, uiF=uiF, v=v, aV=aV, uiV=uiV,plotV=plotV, e=e, aE=aE, edgeTri=tGrid@edges)
		
	
			.Object@tessellation <-tessellation
			.Object@div <-tGrid@div
			
		# the face belts of the hexagrid
			# the total number of subdivisons
			prodTess<-prod(tessellation)
			
			# the ascending part of the face indices
			beltIndicesPart<-rep(NA, prodTess+1)
			beltIndicesPart[1]<-1
			
			# increasing by 5 in every belt
			for(i in 2:length(beltIndicesPart)){
				beltIndicesPart[i]<-(i-1)*5
			}
			beltIndicesMid<-rep(beltIndicesPart[i], prodTess-1)
			
			beltIndices<-c(beltIndicesPart, beltIndicesMid, rev(beltIndicesPart))
			beltIndices<-cumsum(beltIndices)
			
			# transform the indices to belt numbers
			belts<-rep(length(beltIndices), max(beltIndices))
			for(i in (length(beltIndices)-1):1){
				belts[beltIndices[i]:1]<-i
			}
			.Object@belts<-belts
		
		
		#user interface
			.Object@faces <- faces
			.Object@vertices <- vertices
			.Object@edges <- edges
			.Object@skeleton <- skeleton
			.Object@faceCenters <- faceCenters
			
			#edge Length calcualtion
			.Object@edgeLength <- c(
				mean(.Call(Cpp_icosa_edges_, skeleton$v, skeleton$e, tGrid@center, 1)),
				mean(.Call(Cpp_icosa_edges_, skeleton$v, skeleton$e, tGrid@center, 0))/pi*180)
				names(.Object@edgeLength) <- c("km", "deg")
			
			.Object@length<- c(
				"vertices"=nrow(.Object@vertices),
				"edges"=nrow(.Object@edges),
				"faces"=nrow(.Object@faces))
				
		
			#temporary solutions
			.Object@r <- tGrid@r
			.Object@proj4string<- tGrid@proj4string
			.Object@orientation<-c(0,0,0)
			.Object@center <- center
		
		
		# calculate the graph representation
		dummy<-NA
		.Object@graph<-dummy		
		if(graph==TRUE){
			.Object@graph<-gridgraph(.Object)
		}		
		
		if(sp==TRUE){
				.Object@sp<-SpPolygons(.Object, res=25)
		}else{
			dummy<-NA
			.Object@sp<-dummy
		}
		
		# translate the 3d information with the center
		for(vc in 1:3){
			.Object@vertices[,vc] <- .Object@vertices[,vc]+center[vc]
			.Object@faceCenters[,vc] <- .Object@faceCenters[,vc]+center[vc]
			.Object@skeleton$v[,vc] <- .Object@skeleton$v[,vc]+center[vc]
			.Object@skeleton$plotV[,vc] <- .Object@skeleton$plotV[,vc]+center[vc]
		}
		options(scipen=0)	
		return(.Object)
		
	}
)
	
		





