#' Tessellation guide to \code{\link{hexagrid}} objects
#' 
#' The table includes basic properties of \code{\link{hexagrid}}s described with specific tessellation parameters
#' 
#' @format A \code{data.frame} with 129 observations and 18 variables:
#' \describe{
#' 	\item{\code{total}}{The total tessellation of the grid, the number of points inserted between icosahedron vertices along an edge.}
#' 	\item{\code{level1}}{Level 1 tessellation.}
#' 	\item{\code{level2}}{Level 2 tessellation - second value of the tessellation vector. }
#' 	\item{\code{level3}}{Level 3 tessellation - third value of the tessellation vector. }
#' 	\item{\code{level4}}{Level 4 tessellation - four value of the tessellation vector. }
#' 	\item{\code{faces}}{The number of faces in the grid.}
#' 	\item{\code{vertices}}{The number of vertices in the grid.}
#' 	\item{\code{meanEdgeLength_deg}}{Mean edge length in degrees.}
#' 	\item{\code{sdEdgeLength_deg}}{Standard deviation of edge length in degrees.}
#' 	\item{\code{meanEdgeLength_km}}{Mean edge length in kilometers.}
#' 	\item{\code{sdEdgeLength_km}}{Standard devation of edge length in kilometers.}
#' 	\item{\code{meanArea_km2}}{Mean face area in square-kilometers.}
#' 	\item{\code{sdArea_km2}}{Standard deviation of face area in square-kilometers.}
#' 	\item{\code{time}}{Time to compute grid with an Intel Xeon E-1650 prcessor.}
#' 	\item{\code{time_sp}}{Time to compute grid with an Intel Xeon E-1650 prcessor, with the 'sp' member.}
#' 	\item{\code{size}}{The size of the grid in bytes.}
#' 	\item{\code{size_sp}}{The size of the grid object in bytes, with the 'sp' member.}
#' 	\item{\code{timeLocate_5000}}{Time to locate 5000 points with an Intel Xeon E-1650 processor in seconds.}
#' }
"tessguide"