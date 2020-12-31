/*
 * Implementation of an algorithm which computes least-squares assignments
 * between two (discrete) measures in R^2, where the source measure
 * is supposed to be given on the nxm-grid { (i+0.5,j+0.5) | 0<=i<m, 0<=j<n } and 
 * the target measure can be arbitrary, as long it is discrete with a finite number of support 
 * points within the nxm grid. The source measure is 
 * interpreted as a pixel density on the squares { [i,j]x[i+1,j+1] | 0<=i,j<n }. 
 * This solves optimal transport between those distributions with respect 
 * to the squared euclidean distance as the cost function.
 *
 * See
 *
 * [1] Aurenhammern, Hoffmann, Aronov: Minkowsky-Type
 * Theorems and Least-Squares Clustering, Algorithmica 20, 1 (1998), 61-76.
 *
 * [2] Merigot: A multiscale approach to optimal transport, Computer
 * Graphics Forum 30, 5 (2011), 1584-1592.
 *
 * Code by Bjoern Baehre
 */

#include <R.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include "power_diagram.h"

#define AHA_TRANSPORT_MEMORY 32000

// some global variables which make it
// easier to use our c-functions in R
Triangulation aha_rt;
double *aha_x;
double *aha_y;
double *aha_transport_from;
double *aha_transport_to;
double *aha_transport_mass;
double aha_rect[4];
double aha_pert;
// image size = aha_n*aha_m
int aha_n;
int aha_m;
double *aha_dphi_val;

// global variables used to raster a cell
int *aha_edge_pixel; // n*m values, marking pixels intersected by edges
double *aha_area; // n*m values, area of pixel intersected with power cell
int *aha_ixmin;
int *aha_ixmax;
int aha_iymin;
int aha_iymax;

/* determines n weighted points (x,y,m) which are close to the n0 points (x0,y0,m0)
 * 'in the Wasserstein sense' by using the (weighted) standard k-means algorithm,
 * which runs as long as at least on point of (x,y) is moved a distance greater than eps;
 * p maps the points of (x0,y0) to their corresponding cluster */
void decompose_c(int *n, double *x, double *y, double *m, int *n0, double *x0, double *y0,
               double *m0, int *p, double *eps)
{
    int i,j;
    double *cluster_x = Calloc(*n0*sizeof(double),double);
    double *cluster_y = Calloc(*n0*sizeof(double),double);
    double *cluster_m = Calloc(*n0*sizeof(double),double);
    double a,min;
    int imin;
    double e;

    // weighted k-means algorithm
    do {
        // reset size and weight of clusters
        for (i=0;i<*n;i++) {
            cluster_x[i] = 0;
            cluster_y[i] = 0;
            cluster_m[i] = 0;
        }

        // assign each point of (x0,y0,m0) to a cluster in such a way
        // that the squared distances are minimized
        for (i=0;i<*n0;i++) {
            imin = 0;
            min = ((x0[i]-x[0])*(x0[i]-x[0])+(y0[i]-y[0])*(y0[i]-y[0]));
            for (j=1;j<*n;j++) {
                a = ((x0[i]-x[j])*(x0[i]-x[j])+(y0[i]-y[j])*(y0[i]-y[j]));
                if (a<min) {
                    min = a;
                    imin = j;
                }
            }
            p[i] = imin;
            cluster_x[imin] += m0[i]*x0[i];
            cluster_y[imin] += m0[i]*y0[i];
            cluster_m[imin] += m0[i];
        }

        e = 0.0;
        // compute new cluster coordinates, using the weighted mean
        for (i=0;i<*n;i++) {
            if (cluster_m[i]>0) {
                cluster_x[i] = cluster_x[i]/(cluster_m[i]);
                cluster_y[i] = cluster_y[i]/(cluster_m[i]);
                a = (cluster_x[i]-x[i])*(cluster_x[i]-x[i])+
                    (cluster_y[i]-y[i])*(cluster_y[i]-y[i]);
                e = (a>e ? a : e);
                x[i] = cluster_x[i];
                y[i] = cluster_y[i];
            }

            m[i] = cluster_m[i];
        }

    } while(e>(*eps)*(*eps));

    Free(cluster_x);
    Free(cluster_y);
    Free(cluster_m);
}

/* area of [0,1]x[0,1] intersected with halfspace to the left of (x1,y1)->(x2,y2),
 * assuming that the edge is neither horizontal nor vertical */
double pixel_edge_area(double x1, double y1, double x2, double y2)
{
    double a, b, rx1, rx2, ry1, ry2, area;
    int onx1, onx2, ony1, ony2;
    a = (y2-y1)/(x2-x1);
    b = y1-a*x1;

    rx1 = -b/a;
    rx2 = (1-b)/a;
    ry1 = b;
    ry2 = a+b;

    onx1 = (0 <= rx1 && rx1 <= 1);
    onx2 = (0 <= rx2 && rx2 <= 1);
    ony1 = (0 <= ry1 && ry1 <= 1);
    ony2 = (0 <= ry2 && ry2 <= 1);

    if (ony1 && ony2) {
        area = (ry1+ry2)/2.0;
        if (x2 > x1) area = 1.0-area;
    } else if (onx1 && onx2) {
        area = 1.0 - (rx1+rx2)/2.0;
        if (y2 > y1) area = 1.0-area;
    } else if (onx1 && ony2) {
        area = (1.0-rx1)*ry2/2.0;
        if (x1 < x2) area = 1.0-area;
    } else if (ony1 && onx2) {
        area = (1.0-ry1)*rx2/2.0;
        if (x2 < x1) area = 1.0-area;
    } else if (ony1 && onx1) {
        area = ry1*rx1/2.0;
        if (x1 < x2) area = 1.0-area;
    } else { // if (onx2 && ony2) {
        area = (1.0-ry2)*(1.0-rx2)/2.0;
        if (x2 < x1) area = 1.0-area;
    }

    return area;
}

/* corresponds to the integral of ||x||^2 over the triangle given by (v1,v2,v3) */
double triangle_integral(double v1x, double v1y, double v2x, double v2y, double v3x, double v3y)
{
    double a = v2x-v1x;
    double b = v3x-v1x;
    double c = v2y-v1y;
    double d = v3y-v1y;
    return fabs((a*d-b*c))*((a*a+b*b+c*c+d*d+a*b+c*d)/12.0 + 
            (v1x*(a+b)+v1y*(c+d))/3.0+(v1x*v1x+v1y*v1y)/2.0);
} 


/* given an edge (x1,y1)->(x2,y2), determines integers 0<=ix1<=ix2<max s.t.
 * ix1<=x1,x2<=ix2 but ix1+1>(x1 or x2) and ix2-1<(x1 or x2);
 * to use it for y-coordinates, call pixel_range(&iy1, &iy2, y1, y2, x2, x1, max') */
void pixel_range(int *ix1, int *ix2, double x1, double x2, double y1, double y2, int max)
{
    int tx;

    // special case: vertical edge with integer x-coordinate;
    // depending on the y-coordinates, we want to consider 
    // either the pixels to the left or to the right
    if (x1==x2 && x1==round(x1)) {
        *ix1 = (int) round(x1);
        if (y1<y2) *ix1 -= 1;
        *ix2 = *ix1;
    } else {
        *ix1 = (int) floor(x1);
        *ix2 = (int) floor(x2);
    }

    // make sure that ix1<=ix2
    if (*ix2<*ix1) {
        tx = *ix1;
        *ix1 = *ix2;
        *ix2 = tx;
    }

    // make sure that 0<=ix1<=ix2<max
    if (*ix1 < 0) *ix1 = 0;
    if (*ix2 < 0) *ix2 = 0;
    if (*ix1 >= max) *ix1 = max-1;
    if (*ix2 >= max) *ix2 = max-1;
}

/* compute the pixels which intersect the power cell of v;
 * the result will be saved in global variables
 * int *aha_ixmin;
 * int *aha_ixmax;
 * int aha_iymin;
 * int aha_iymax;
 * int *aha_edge_pixel;
 * double *aha_area
 */
void raster_cell(Site *v, int n, double *x, double *y)
{
    if (n<3) return;
    double a, b, x1, x2, y1, y2, ty1, ty2, ymin=R_PosInf, ymax=R_NegInf;
    int ix1, ix2, iy1, iy2, ity1, ity2, i, j, k;

    // determine min/max of cell y-coordinates
    for (i=0;i<n;i++) {
        if (y[i]<ymin) ymin = y[i];
        if (y[i]>ymax) ymax = y[i];
    }

    // reset ixmin/ixmax
    pixel_range(&aha_iymin,&aha_iymax,ymin,ymax,0,0,aha_n);
    for (i=aha_iymin;i<=aha_iymax;i++) {
        aha_ixmin[i] = aha_m-1;
        aha_ixmax[i] = 0;
    }

    // determine x-coordinates of contributing pixels (with y-coordinates
    // between iymin and iymax), iterate over edges
    for (k=0;k<n;k++) {
        // get edge coordinates
        x1 = x[k];
        y1 = y[k];
        x2 = x[(k<n-1 ? k+1 : 0)];
        y2 = y[(k<n-1 ? k+1 : 0)];
        pixel_range(&ix1,&ix2,x1,x2,y1,y2,aha_m);
        pixel_range(&iy1,&iy2,y1,y2,x2,x1,aha_n);
        // vertical edge
        if (x1==x2) {
            for (j=iy1;j<=iy2;j++) {
                ++(aha_edge_pixel[j*aha_m+ix1]);
                aha_area[j*aha_m+ix1] = (y1 > y2 ? ix1+1.0-x1 : x1-ix1);
                if (ix1<aha_ixmin[j]) aha_ixmin[j] = ix1;
                if (ix1>aha_ixmax[j]) aha_ixmax[j] = ix1;
            }
        // horizontal edge
        } else if (y1==y2) {
            for (i=ix1;i<=ix2;i++) {
                ++(aha_edge_pixel[iy1*aha_m+i]);
                aha_area[iy1*aha_m+i] = (x1 < x2 ? iy1+1.0-y1 : y1-iy1);
            }
            if (ix1<aha_ixmin[iy1]) aha_ixmin[iy1] = ix1;
            if (ix2>aha_ixmax[iy1]) aha_ixmax[iy1] = ix2;
        } else { 
            // y=a*x+b
            a = (y2-y1)/(x2-x1);
            b = y1-a*x1;

            for (i=ix1;i<=ix2;i++) {
            	ty1 = a*i+b;
            	ty2 = a*(i+1)+b;

            	if (i==ix1)
                	ty1 = (x1<x2 ? y1 : y2);
            	if (i==ix2)
                	ty2 = (x1<x2 ? y2 : y1);

            	pixel_range(&ity1,&ity2,ty1,ty2,0,0,aha_n);
                for (j=ity1;j<=ity2;j++) {
                    ++(aha_edge_pixel[j*aha_m+i]);
                    aha_area[j*aha_m+i] = pixel_edge_area(x1-i,y1-j,x2-i,y2-j);
                    if (i<aha_ixmin[j]) aha_ixmin[j] = i;
                    if (i>aha_ixmax[j]) aha_ixmax[j] = i;
                }
            }
        }
    }
}

void cell_integral(Site *s, double *source_measure, double *integral, int cost, int weight, int exact, int reset)
{
    int ix,iy,j,k;
    double c=0,area,x,y;
    *integral = 0;

    for (iy=aha_iymin;iy<=aha_iymax;iy++) {
        for (ix=aha_ixmin[iy];ix<=aha_ixmax[iy];ix++) {
            if (!exact) {
                if (cost)
                    c = (ix+0.5-s->x)*(ix+0.5-s->x)+(iy+0.5-s->y)*(iy+0.5-s->y);
                if (aha_edge_pixel[iy*aha_m+ix] > 1) {
                    if (reset) aha_edge_pixel[iy*aha_m+ix] = 0;
                    power_cell(&k,aha_x,aha_y,&aha_rt,s,ix,iy,ix+1.0,iy+1.0);
                    area = polygon_area(k,aha_x,aha_y);
                } else if (aha_edge_pixel[iy*aha_m+ix] == 1) {
                    if (reset) aha_edge_pixel[iy*aha_m+ix] = 0;
                    area = aha_area[iy*aha_m+ix];
                } else {
                    area = 1;
                }
                if (weight && cost) *integral += (c-s->w)*area*source_measure[ix*aha_n+aha_n-1-iy];
                else if (cost) *integral += c*area*source_measure[ix*aha_n+aha_n-1-iy];
                else *integral += area*source_measure[ix*aha_n+aha_n-1-iy];
            } else {
                if (aha_edge_pixel[iy*aha_m+ix] > 0) {
                    if (reset) aha_edge_pixel[iy*aha_m+ix] = 0;
                    power_cell(&k,aha_x,aha_y,&aha_rt,s,ix,iy,ix+1.0,iy+1.0);
                    if (k<3) continue;
                    if (cost) {
                        c = 0;
                        for (j=1;j<k-1;j++) {
                            c += triangle_integral(aha_x[0]-s->x,aha_y[0]-s->y,aha_x[j]-s->x,
                                    aha_y[j]-s->y,aha_x[j+1]-s->x,aha_y[j+1]-s->y);
                        }
                    }
                    if (weight && cost) *integral += (c-s->w*polygon_area(k,aha_x,aha_y))*source_measure[ix*aha_n+aha_n-1-iy];
                    else if (cost) *integral += c*source_measure[ix*aha_n+aha_n-1-iy];
                    else *integral += polygon_area(k,aha_x,aha_y)*source_measure[ix*aha_n+aha_n-1-iy];
                } else {
                    x=ix-s->x;
                    y=iy-s->y;
                    if (cost && weight) *integral += source_measure[ix*aha_n+aha_n-1-iy]*(x*(1+x)+(y*(1+y))+2.0/3.0-s->w);
                    else if (cost) *integral += source_measure[ix*aha_n+aha_n-1-iy]*(x*(1+x)+(y*(1+y))+2.0/3.0);
                    else *integral += source_measure[ix*aha_n+aha_n-1-iy];
                }
            }
        }
    }
}

/* objective function */
void aha_phi(int *n, double *x, double *y, double *w, double *source_measure, double *target_measure, int *exact, double *res)
{
    int i,k;
    double integral;
    *res = 0;

    triangulate(&aha_rt, *n, x, y, w, aha_rect[0], aha_rect[1], aha_rect[2], aha_rect[3], aha_pert);

    for (i=0;i<aha_rt.size;i++) {
        power_cell(&k,aha_x,aha_y,&aha_rt,&aha_rt.sites[i],aha_rect[0],aha_rect[1],aha_rect[2],aha_rect[3]);
        if (k>2) {
            raster_cell(&aha_rt.sites[i],k,aha_x,aha_y);
            cell_integral(&aha_rt.sites[i],source_measure,&integral,1,1,*exact,0);
            *res += w[i]*target_measure[i] + integral;

            cell_integral(&aha_rt.sites[i],source_measure,&integral,0,0,*exact,1);
            aha_dphi_val[i] = target_measure[i] - integral;
        } else {
            *res += w[i]*target_measure[i];
            aha_dphi_val[i] = target_measure[i];
        }
    }
}

/* gradient of objective function */
void aha_dphi(int *n, double *x, double *y, double *w, double *source_measure, double *target_measure, int *exact, double *res)
{
    int i;

    for (i=0;i<aha_rt.size;i++) {
        res[i] = aha_dphi_val[i];
    }
}


/* Computes the L^2-Wasserstein distance */
void aha_wasserstein(int *n, double *x, double *y, double *w, double *source_measure, double *res)
{
    int i,k;
    double integral;

    *res = 0.0;
    triangulate(&aha_rt, *n, x, y, w, aha_rect[0], aha_rect[1], aha_rect[2], aha_rect[3], aha_pert);

    for (i=0;i<aha_rt.size;i++) {
        power_cell(&k,aha_x,aha_y,&aha_rt,&aha_rt.sites[i],aha_rect[0],aha_rect[1],aha_rect[2],aha_rect[3]);
        if (k>2) {
            raster_cell(&aha_rt.sites[i],k,aha_x,aha_y);
            cell_integral(&aha_rt.sites[i],source_measure,&integral,1,0,1,1);
            *res += integral;
        }
    }

    *res = sqrt(*res);
}

/* computes the transport induced by (x,y,w) w.r.t. the source measure 
 * and write the result into global variables
 * double *aha_transport_from;
 * double *aha_transport_to;
 * double *aha_transport_mass;
 */
void aha_compute_transport(int *n, double *x, double *y, double *w, double *source_measure, int *res)
{
    int from, to, k, j, ix, iy, memory = 1;
    double integral;

    *res = 0; // number of transportations, i.e. number of (from,to) value pairs
    triangulate(&aha_rt, *n, x, y, w, aha_rect[0], aha_rect[1], aha_rect[2], aha_rect[3], aha_pert);
    // to each vertex transport the mass which is covered by its power cell
    for (to=0;to<aha_rt.size;to++) {
        power_cell(&k,aha_x,aha_y,&aha_rt,&aha_rt.sites[to],aha_rect[0],aha_rect[1],aha_rect[2],aha_rect[3]);
        if (k>2) {
            raster_cell(&aha_rt.sites[to],k,aha_x,aha_y);
            // intersect each pixel with the power cell of the current vertex
            // and multiply the resulting area with the pixel's weight
            // the get the amount of mass which is transported from the pixel
            // to the current vertex
            for (iy=aha_iymin;iy<=aha_iymax;iy++) {
                for (ix=aha_ixmin[iy];ix<=aha_ixmax[iy];ix++) {
                    from = ix*aha_n+aha_n-1-iy;
                    if (aha_edge_pixel[iy*aha_m+ix] > 0) {
                        aha_edge_pixel[iy*aha_m+ix] = 0;
                        power_cell(&j,aha_x,aha_y,&aha_rt,&aha_rt.sites[to],ix,iy,ix+1.0,iy+1.0);
                        integral = polygon_area(j,aha_x,aha_y)*source_measure[from];
                    } else {
                        integral = source_measure[from];
                    }
                    if (integral==0) continue;
                    aha_transport_from[*res] = from; 
                    aha_transport_to[*res] = to; 
                    aha_transport_mass[*res] = integral;
                    ++(*res);
                    if ((*res)>=memory*AHA_TRANSPORT_MEMORY) {
                        ++memory;
                        aha_transport_from = Realloc(aha_transport_from, memory*AHA_TRANSPORT_MEMORY*sizeof(double),double);
                        aha_transport_to = Realloc(aha_transport_to, memory*AHA_TRANSPORT_MEMORY*sizeof(double),double);
                        aha_transport_mass = Realloc(aha_transport_mass, memory*AHA_TRANSPORT_MEMORY*sizeof(double),double);
                    }
                }
            }
        }
    }
}

/* copy the results from aha_compute_transport for usage in R */
void aha_get_transport(int *size, double *from, double *to, double *mass)
{
    int i;
    for (i=0;i<*size;i++) {
        from[i] = aha_transport_from[i];
        to[i] = aha_transport_to[i];
        mass[i] = aha_transport_mass[i];
    }
}

void aha_init(int *n, int *m, double *rect)
{
    int k;
    aha_n = *n;
    aha_m = *m;
    aha_rect[0] = rect[0];
    aha_rect[1] = rect[1];
    aha_rect[2] = rect[2];
    aha_rect[3] = rect[3];
    aha_x = Calloc(((aha_n+1)*(aha_m+1)+8)*sizeof(double),double);
    aha_y = Calloc(((aha_n+1)*(aha_m+1)+8)*sizeof(double),double);
    aha_ixmin = Calloc(aha_n*sizeof(int),int);
    aha_ixmax = Calloc(aha_n*sizeof(int),int);
    aha_edge_pixel = Calloc(aha_n*aha_m*sizeof(int),int);
    aha_area = Calloc(aha_n*aha_m*sizeof(double),double);
    aha_dphi_val = Calloc(aha_n*aha_m*sizeof(double),double);
    aha_transport_from = Calloc(AHA_TRANSPORT_MEMORY*sizeof(double),double);
    aha_transport_to = Calloc(AHA_TRANSPORT_MEMORY*sizeof(double),double);
    aha_transport_mass = Calloc(AHA_TRANSPORT_MEMORY*sizeof(double),double);

    aha_pert = 0;

    for(k=0;k<aha_n*aha_m;k++) {
        aha_edge_pixel[k] = 0;
        aha_area[k] = 0.0;
    }
    init_triangulation(&aha_rt);
}

void aha_free()
{
    Free(aha_x);
    Free(aha_y);
    Free(aha_ixmin);
    Free(aha_ixmax);
    Free(aha_edge_pixel);
    Free(aha_area);
    Free(aha_dphi_val);
    Free(aha_transport_from);
    Free(aha_transport_to);
    Free(aha_transport_mass);
    free_triangulation(&aha_rt);
}

/*
 * Code for power diagrams.
 */

#define PD_DEFAULT_MEMSIZE 32000

double *pd_x;
double *pd_y;

void compute_power_diagram(int *cell_size, int *n, double *x, double *y, double *w, double *rect)
{
    Triangulation rt;
    static int memory = 1;
    int i, j, k=0, total=0;
    double xmin, ymin, xmax, ymax, pert, xt, yt, a;
    double *cx = Calloc((*n+4)*sizeof(double),double);
    double *cy = Calloc((*n+4)*sizeof(double),double);

    pd_x = Calloc(memory*PD_DEFAULT_MEMSIZE*sizeof(double),double);
    pd_y = Calloc(memory*PD_DEFAULT_MEMSIZE*sizeof(double),double);

    xmin = x[0];
    xmax = x[0];
    ymin = y[0];
    ymax = y[0];
    for (i=0;i<*n;i++) {
        xmin = (x[i]<xmin ? x[i] : xmin);
        xmax = (x[i]>xmax ? x[i] : xmax);
        ymin = (y[i]<xmin ? y[i] : ymin);
        ymax = (y[i]>ymax ? y[i] : ymax);
    }

    pert = R_PosInf;
    for (i=0;i<*n;i++) {
        for (j=i+1;j<*n;j++) { 
            xt = fabs(x[i]-x[j]);
            yt = fabs(y[i]-y[j]);
            a = (xt > yt ? xt : yt);
            if (a < pert) pert = a;
        }
    }

    init_triangulation(&rt);
    triangulate(&rt, *n, x, y, w, xmin-pert, ymin-pert, xmax+pert, ymax+pert, pert*1e-5);

    for (i=0;i<rt.size;i++) {
        power_cell(&cell_size[i],cx,cy,&rt,&rt.sites[i],rect[0],rect[1],rect[2],rect[3]);
        total += cell_size[i];
        if (total > memory*PD_DEFAULT_MEMSIZE) {
            pd_x = Realloc(pd_x, (++memory)*PD_DEFAULT_MEMSIZE*sizeof(double),double);
            pd_y = Realloc(pd_y, (++memory)*PD_DEFAULT_MEMSIZE*sizeof(double),double);
        }
        for(j=0;j<cell_size[i];j++) {
            pd_x[k] = cx[j];
            pd_y[k++] = cy[j];
        }
    }

    free_triangulation(&rt);
    Free(cx);
    Free(cy);
}

void get_power_diagram(int *size, double *x, double *y)
{
    int i;
    for(i=0;i<(*size);i++) {
        x[i] = pd_x[i];
        y[i] = pd_y[i];
    }
    Free(pd_x);
    Free(pd_y);
}
