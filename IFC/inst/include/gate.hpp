/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2020 Yohann Demont                                              
                                                                                
  It is part of IFC package, please cite:                                       
  -IFC: An R Package for Imaging Flow Cytometry                                 
  -YEAR: 2020                                                                   
  -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,              
                      Jean-Pierre Marolleau, Loïc Garçon,                       
                      INSERM, UPD, CHU Amiens                                   
                                                                                
                                                                                
  DISCLAIMER:                                                                   
  -You are using this package on your own risk!                                 
  -We do not guarantee privacy nor confidentiality.                             
  -This program is distributed in the hope that it will be useful, but WITHOUT  
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or         
  FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or  
  contributors be liable for any direct, indirect, incidental, special,         
  exemplary, or consequential damages (including, but not limited to,           
  procurement of substitute goods or services; loss of use, data, or profits;   
  or business interruption) however caused and on any theory of liability,      
  whether in contract, strict liability, or tort (including negligence or       
  otherwise) arising in any way out of the use of this software, even if        
  advised of the possibility of such damage.                                    
                                                                                
  You should have received a copy of the GNU General Public License             
  along with IFC. If not, see <http://www.gnu.org/licenses/>.                   
*/

#ifndef IFC_GATE_HPP
#define IFC_GATE_HPP

#include <Rcpp.h>
using namespace Rcpp;

//' @title Polygon Closing
//' @name close_polygon
//' @description
//' Pushes 1st row at the end of a matrix.
//' @param M NumericMatrix.
//' @keywords internal
Rcpp::NumericMatrix close_polygon (const Rcpp::NumericMatrix M) {
  R_len_t i, n = M.nrow();
  Rcpp::NumericMatrix MM = Rcpp::no_init_matrix(n + 1, M.ncol());
  for(i = 0; i < n; i++) MM(i, _) = M(i, _);
  MM(i, _) = M(0, _);
  return MM;
}

//' @title Point in Polygon
//' @name trigo_pnt_in_poly
//' @description
//' This function works out if 2D point lies within the boundaries of a
//' defined polygon. \cr \cr \bold{Note:} Points that lie on the boundaries of
//' the polygon or vertices are assumed to be within the polygon.
//' 
//' The algorithm implements a sum of the angles made between the test point and
//' each pair of points making up the polygon. The point is interior if the sum
//' is 2pi, otherwise, the point is exterior if the sum is 0. This works for
//' simple and complex polygons (with holes) given that the hole is defined with
//' a path made up of edges into and out of the hole. \cr \cr This sum of angles
//' is not able to consistently assign points that fall on vertices or on the
//' boundary of the polygon. The algorithm defined here assumes that points
//' falling on a boundary or polygon vertex are part of the polygon.
//' @param pnt NumericVector, x and y coordinates of the points.
//' @param poly NumericMatrix, a 2-column matrix defining the locations (x and y) of vertices of the polygon of interest.
//' @param epsilon double, threshold value.  Default is 0.000000000001
//' @author Jeremy VanDerWal, Lorena Falconi, Stephanie Januchowski, Luke Shoo and Collin Storlie \email{jjvanderwal@@gmail.com}.
//' @source \url{https://github.com/jjvanderwal/SDMTools}
//' @keywords internal
bool trigo_pnt_in_poly (const Rcpp::NumericVector pnt,
                        const Rcpp::NumericMatrix poly,
                        const double epsilon = 0.000000000001 ) {
  //define some other variables
  R_len_t i, n = poly.nrow();
  double x, x1, x2, y, y1, y2, dy, dx, dd, angle = 0.0;
  //cycle through the polygon vertices and sum the angles
  for (i = 0; i < n - 1; i++) {
    //define the points
    x1 = poly(i, 0); x2 = poly((i + 1), 0); x = pnt[0];
    y1 = poly(i, 1); y2 = poly((i + 1), 1); y = pnt[1];
    //check if point is vertix
    if (x == x1 && y == y1) { angle = PI + 1; break; }
    //check if point is on border line between 2 points
    if (x == x1 && x == x2) { if ((y1 <= y && y <= y2) || (y1 >= y && y >= y2)) { angle = PI + 1; break; } } // check point between two horizontal points
    if (y == y1 && y == y2) { if ((x1 <= x && x <= x2) || (x1 >= x && x >= x2)) { angle = PI + 1; break; } } // check point between two verticle points
    dy = (y1==y2) ? -9999:(y1-y)/(y1-y2); //check if the relative change in x == relative change in y
    dx = (x1==x2) ? -9999:(x1-x)/(x1-x2); //check if the relative change in x == relative change in y
    dd = dy-dx; dd = (dd<0) ? -dd:dd;
    if (dd < epsilon && dy>0 && dy<1) { angle = PI + 1; break; } // if dx == dy and dy is between 0 & 1 ... point is on the border line
    // && dy > 0 && dy < 1
    //if not a vertex or on border lines... sum the angles
    double dtheta = std::atan2(y2 - y, x2 - x) - std::atan2(y1 - y, x1 - x);
    while (dtheta > PI) dtheta -= 2 * PI;
    while (dtheta < -PI) dtheta += 2 * PI;
    angle += dtheta;
  }
  return (std::fabs(angle) >= PI);
}

//' @title Point in Ellipse
//' @name pnt_in_ell
//' @description
//' Checks if 2D point lies within a defined ellipse.
//' @param pnt NumericVector, x and y coordinates of the point.
//' @param ell NumericVector, coordinates of the ellipse.
//' @keywords internal
bool pnt_in_ell (const Rcpp::NumericVector pnt,
                 const Rcpp::NumericVector ell) {
  return (( ((pnt[0] - ell[2]) * (pnt[0] - ell[2]) / ell[0] / ell[0]) + ((pnt[1] - ell[3]) * (pnt[1] - ell[3]) / ell[1] / ell[1])) <= 1);
}

//' @title Ellipse Boundaries to Coordinates
//' @name cpp_ell_coord
//' @description
//' Transforms ellipse boundaries to usefull coordinates.
//' @param bound_x NumericVector, x-boundaries of the ellipse.
//' @param bound_y NumericVector, y-boundaries of the ellipse.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector hpp_ell_coord (const Rcpp::NumericVector bound_x,
                                   const Rcpp::NumericVector bound_y) {
  double xmin, ymin, xmax, ymax;
  Rcpp::NumericVector out(4); 
  xmin = min(bound_x);
  ymin = min(bound_y);
  xmax = max(bound_x);
  ymax = max(bound_y);
  out[0] = (xmax - xmin) / 2; // rx
  out[1] = (ymax - ymin) / 2; // ry
  out[2] = xmin + out[0];     // cx
  out[3] = ymin + out[1];     // cy
  return out;
}

//' @title Point in Gate
//' @name cpp_pnt_in_gate
//' @description
//' This function checks if points lie in a polygon or ellipse.
//' @param pnts NumericMatrix, a 2-columns matrix with (x and y) coordinates of the points of interest.
//' @param gate NumericMatrix, a 2-columns matrix defining polygon vertices or ellipse boundaries.
//' @param algorithm int, used for computation. Default is 1.\cr
//' 1: Trigonometry.\cr
//' 2: Special case = axes-aligned rectangle.\cr
//' 3: Special case = axes-aligned ellipse.
//' @param epsilon double, epsilon threshold value. Default is 0.000000000001
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalVector hpp_pnt_in_gate (const Rcpp::NumericMatrix pnts,
                                     const Rcpp::NumericMatrix gate,
                                     const int algorithm = 1,
                                     const double epsilon = 0.000000000001 ) {
  R_len_t k, n = gate.nrow(), L = pnts.nrow();
  Rcpp::NumericMatrix P;
  if( (gate(0,0) == gate(n - 1,0)) && (gate(0,1) == gate(n - 1,1)) ) {
    P = gate;
  } else {
    P = close_polygon(gate);
  }
  if((L < 1) || (P.nrow() < 2) || (gate.ncol() != 2) || (pnts.ncol() != 2)) Rcpp::stop("hpp_pnt_in_gate: Bad dimension in pnt_in_poly inputs");
  if((algorithm < 0) || (algorithm > 3)) Rcpp::stop("hpp_pnt_in_gate: 'algorithm' should be 1(Trigonometry), 2(Special case = axes-aligned rectangle) or 3(Special case = axes-aligned ellipse)");
  Rcpp::NumericVector xran = range(P(_,0));
  Rcpp::NumericVector yran = range(P(_,1));
  Rcpp::LogicalVector C(L, false);
  
  switch(algorithm) {
  case 1:
    for(k = 0; k < L; k++) {
      if((pnts(k,0) >= xran[0]) && (pnts(k,0) <= xran[1]) && (pnts(k,1) >= yran[0]) && (pnts(k,1) <= yran[1]))
        C[k] = trigo_pnt_in_poly(pnts(k, _), P, epsilon);
    }
    break;
  case 2:
    for(k = 0; k < L; k++) {
      if((pnts(k,0) >= xran[0]) && (pnts(k,0) <= xran[1]) && (pnts(k,1) >= yran[0]) && (pnts(k,1) <= yran[1]))
        C[k] = true;
    }
    break;
  case 3:
    Rcpp::NumericVector ell = hpp_ell_coord(xran, yran);
    for(k = 0; k < L; k++) {
      if((pnts(k,0) >= xran[0]) && (pnts(k,0) <= xran[1]) && (pnts(k,1) >= yran[0]) && (pnts(k,1) <= yran[1]))
        C[k] = pnt_in_ell(pnts(k, _), ell);
    }
    break;
  }
  // case 4:
  //   for(k = 0; k < L; k++) {
  //     if((pnts(k,0) >= xran[0]) && (pnts(k,0) <= xran[1]) && (pnts(k,1) >= yran[0]) && (pnts(k,1) <= yran[1]))
  //       C[k] = ray_pnt_in_poly(pnts(k, _), P);
  //   }
  //   break;
  // case 5:
  //   for(k = 0; k < L; k++) {
  //     if((pnts(k,0) >= xran[0]) && (pnts(k,0) <= xran[1]) && (pnts(k,1) >= yran[0]) && (pnts(k,1) <= yran[1]))
  //       C[k] = wn_pnt_in_poly(pnts(k, _), P);
  //   }
  //   break;
  return C;
}

#endif
