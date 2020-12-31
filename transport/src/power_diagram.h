/*
 * Functions to compute the power-diagram of a set of
 * weighted points in R^2.
 *
 * See 
 *
 * [1] H. Edelsbrunner, N. R. Shah: Incremental Topological 
 * Flipping Works for Regular Triangulations.  
 * Algorithmica (1996) 15: 223-241.
 * 
 * [2] O. Devillers: The Delaunay Hierarchy.
 * Internat. J. Found. Comput. Sci., 13: 163-180, 2002
 *
 * With ideas from J. R. Shewchuk's 'Triangle'
 * and CGAL's 'Regular_triangulation_2' package
 *
 * Code by Bjoern Baehre
 */

#ifndef POWER_DIAGRAM_H
#define POWER_DIAGRAM_H

#include <R.h>
#include <Rmath.h>

#define MALLOC(size,t) Calloc(size * sizeof(t),t)
#define REALLOC(p,size,t) Realloc(p, size * sizeof(t),t)
#define FREE(p) Free(p)

#include <stdlib.h>
#include <stdio.h>
#include <time.h> 
#include <math.h>
#include <float.h> 
//#include <mpfr.h>

// memory for triangles and link facets is claimed in blocks
#define MEMORY_BLOCKSIZE 64000
#define LINK_FACETS_BLOCKSIZE 1024
// for the meaning of 'level', see [2];
#define MAX_LEVEL 6

typedef struct Site 
{
    double x;
    double y;
    double w;
    int level; // see [2]

    // edges pointing to v
    struct Edge *level_edge;
} 
Site;

// the data structures for triangles and edge handlers
// are in part adapted from J. R. Shewchuk's 'Triangle'
typedef struct Triangle 
{
    // vertices are always listed in counterclockwise order
    Site *vertex[3];
    struct Triangle *neighbour[3];
    // the orientation of each neighbour from
    // the perspective of this triangle;
    // only used by edge handlers
    int neighbour_edge[3];
    int level; // see [2]
} 
Triangle;

// edge handler, used to navigate and operate
// on a triangulation via macros;
// e.g. an edge handler pointing to a triangle abc
// with shift=1 represents the edge bc;
// incrementing the value of shift (modulo 3) means
// iterating over the edges of the triangle in 
// counterclockwise oder 
typedef struct Edge
{
    Triangle *triangle;
    int shift;
} 
Edge;

typedef struct Triangulation 
{
    // it is convienient to have an artificial
    // triangle of 'infinite size'
    Site s_inf1, s_inf2, s_inf3;
    Triangle t_inf;

    // linked-list of memory blocks containing triangles
    struct MBlock {
        Triangle *triangles;
        struct MBlock *next, *prev;
        int i;
    } *root, **current_block;

    Site *sites;
    int size; // number of vertices
    // site of highest level where point location begins
    Site *high_site;
    int highest_level;

    Edge *level_edges; // see [2]
    Edge *link_facets; // see [1]
    int stack_i; // link facets are treated as a stack
} 
Triangulation;

void power_cell(int *, double *, double *, Triangulation *, 
        Site *, double, double, double, double);

int __macro_i = 0;
Edge __macro_e = { NULL, 0 };

// if the edge handler e points to a triangle
// abc and represents the edge ab, then
// ORG(e) = a, DEST(e) = b, APEX(e) = c
#define ORG(e) \
    ((e).triangle->vertex[((e).shift+1)%3])

#define DEST(e) \
    ((e).triangle->vertex[((e).shift+2)%3])

#define APEX(e) \
    ((e).triangle->vertex[(e).shift])

#define SET_ORG(e,v) \
    ((e).triangle)->vertex[((e).shift+1)%3] = (v)

#define SET_DEST(e,v) \
    ((e).triangle)->vertex[((e).shift+2)%3] = (v)

#define SET_APEX(e,v) \
    ((e).triangle)->vertex[(e).shift] = (v)

#define SET_ALL(e,o,d,a) \
    SET_ORG(e,o); \
    SET_DEST(e,d); \
    SET_APEX(e,a)

// if e points to abc and represents the edge ab,
// then SWITCH(e) will point to the neighbour of abc at ab,
// representing the edge ba
#define SWITCH(e) \
    __macro_i = (e).shift; \
    (e).shift = ((e).triangle)->neighbour_edge[(e).shift]; \
    (e).triangle = ((e).triangle)->neighbour[__macro_i]

// next edge in couterclockwise oder
#define NEXT(e) \
    (e).shift = ((e).shift+1)%3

// next edge in clockwise oder
#define PREV(e) \
    (e).shift = ((e).shift+2)%3

// f becomes a copy of e
#define COPY(e,f) \
    (f).triangle = (e).triangle; \
    (f).shift = (e).shift

// equivalent to COPY(e,f); SWITCH(f);
#define TO_SWITCHED(e,f) \
    (f).triangle = (e).triangle->neighbour[(e).shift]; \
    (f).shift = (e).triangle->neighbour_edge[(e).shift]

// equivalent to COPY(e,f); NEXT(f);
#define TO_NEXT(e,f) \
    (f).triangle = (e).triangle; \
    (f).shift = ((e).shift+1)%3

// equivalent to COPY(e,f); PREV(f);
#define TO_PREV(e,f) \
    (f).triangle = (e).triangle; \
    (f).shift = ((e).shift+2)%3

// the triangles e and f point to will become
// neighbours so that SWITCH(e)==f 
// and SWITCH(f)==e holds
#define GLUE(e,f) \
    (e).triangle->neighbour[(e).shift] = (f).triangle; \
    (f).triangle->neighbour[(f).shift] = (e).triangle; \
    (e).triangle->neighbour_edge[(e).shift] = (f).shift; \
    (f).triangle->neighbour_edge[(f).shift] = (e).shift

// equivalent to SWITCH(e); GLUE(e,f); SWITCH(e);
#define SWITCH_GLUE(e,f) \
    TO_SWITCHED(e,__macro_e); \
    GLUE(__macro_e, f)

// vector macros, just for convenience
#define ORIENTATION(a,b,c) \
    orientation((a)->x,(a)->y,(b)->x,(b)->y,(c)->x,(c)->y)

// checks if v is one of the three infinite vertices of rt
#define INF(rt,v) \
    ((v)==&((rt)->s_inf1) || (v)==&((rt)->s_inf2) || (v)==&((rt)->s_inf3))

/* returns 1 if abc are arranged in counterclockwise order,
 * -1 if they are arranged in clockwise order,
 *  0 if they are collinear */
int orientation(double ax, double ay, double bx, double by, double cx, double cy)
{
    double ori = (ax-cx)*(by-cy)-(bx-cx)*(ay-cy);
    if (ori>0) return 1;
    if (ori<0) return -1;
    return 0;

    /*
    int res;
    mpfr_t a,b,c,d;
    mpfr_inits2(512,a,b,c,d,(mpfr_ptr)0); 
    mpfr_set_d(a,ax,MPFR_RNDN);
    mpfr_set_d(b,by,MPFR_RNDN);
    mpfr_set_d(c,bx,MPFR_RNDN);
    mpfr_set_d(d,ay,MPFR_RNDN);

    mpfr_sub_d(a,a,cx,MPFR_RNDN);
    mpfr_sub_d(b,b,cy,MPFR_RNDN);
    mpfr_sub_d(c,c,cx,MPFR_RNDN);
    mpfr_sub_d(d,d,cy,MPFR_RNDN);

    mpfr_mul(a,a,b,MPFR_RNDN);
    mpfr_mul(c,c,d,MPFR_RNDN);
    mpfr_sub(a,a,c,MPFR_RNDN);

    res = mpfr_sgn(a);
    mpfr_clears(a,b,c,d,(mpfr_ptr)0); 
    return res;
    */
}

void orthogonal_center(double *x, double *y, Site *a, Site *b, Site *c)
{
    double D = 2*(c->y*(b->x-a->x)+c->x*(a->y-b->y)-a->y*b->x+a->x*b->y);
    double A = b->x*b->x+b->y*b->y-a->x*a->x-a->y*a->y-b->w+a->w;
    double B = c->x*c->x+c->y*c->y-a->x*a->x-a->y*a->y-c->w+a->w;
    *x = ((c->y-a->y)*A+(a->y-b->y)*B)/D;
    *y = ((a->x-c->x)*A+(b->x-a->x)*B)/D;
}

double pdist(double x, double y, Site *s)
{
    return (x-s->x)*(x-s->x)+(y-s->y)*(y-s->y)-s->w;
}

/* computes the area of a convex polygon (x0,y0)...(xn-1,yn-1) */
double polygon_area(int n, double *x, double *y)
{
    if (n<3) return 0.0;
    int i;
    double a = x[n-1]*y[0] - x[0]*y[n-1];

    for (i=0;i<n-1;i++)
        a += x[i]*y[i+1] - x[i+1]*y[i];
    return 0.5*a;
}

/* create a new triangle for rt so that e.triangle points to it */
void create_triangle(Triangulation *rt, Edge *e, int level)
{
    struct MBlock *b = *rt->current_block;

    // if current memory block is full, create a new one
    if (b->i == MEMORY_BLOCKSIZE) {
        if (b->next == NULL) {
            b->next = MALLOC(1, struct MBlock);
            if (b->next == NULL) {
                error("ERROR power_diagram.h: malloc failed\n");
            }
            b->next->triangles = MALLOC(MEMORY_BLOCKSIZE, Triangle);
            if ((b->next)->triangles == NULL) {
                error("ERROR power_diagram.h: malloc failed\n");
            }
            b->next->next = NULL;
            b->next->prev = b;
        }
        b->next->i = 0;
        rt->current_block = &b->next;
    }

    Triangle *t = &(*rt->current_block)->triangles[((*rt->current_block)->i)++];
    (*t) = (Triangle) { { NULL, NULL, NULL }, { &(rt->t_inf), &(rt->t_inf), &(rt->t_inf) }, { 0, 0, 0 }, level };
    (*e) = (Edge) { t, 0 };
}

void push_link(Triangulation *rt, Edge e)
{
    static int k = 1;
    if (rt->stack_i == k*LINK_FACETS_BLOCKSIZE) {
        ++k;
        rt->link_facets = REALLOC(rt->link_facets, k*LINK_FACETS_BLOCKSIZE, Edge);
    }
    Edge *f = &rt->link_facets[(rt->stack_i)++];
    COPY(e,*f);
}

Edge *pop_link(Triangulation *rt)
{
    if (rt->stack_i == 0)
        return NULL;
    return &rt->link_facets[--(rt->stack_i)];
}

/* relative to the origin of e */
int order(Edge e)
{
    Site *s = DEST(e);
    int i=1;
    PREV(e);
    SWITCH(e);
    // optimization: we are only interested in orders lower than 5
    while (i<5 && s != DEST(e)) {
        ++i;
        PREV(e);
        SWITCH(e);
    } 
    return i;
}

// The following flipping operations are explained in [1].
// In general, they remove or create triangles in the current
// triangulation, using edge handlers to glue them together

/* inserts v into the triangle associated with e */
void flip1_3(Triangulation *rt, Site *v, Edge e, int level)
{
    Edge f,g,h;
    int i;

    e.triangle->level = -1;
    for (i=0;i<3;i++) {
        create_triangle(rt, &f, level);
        SET_ALL(f, ORG(e), DEST(e), v);
        SWITCH_GLUE(e,f);
        COPY(f,DEST(f)->level_edge[level]);
        push_link(rt, f);
        if (i==0) {
            COPY(f,g);
            PREV(g);
            NEXT(f);
        } else {
            PREV(f);
            GLUE(f,h);
            PREV(f);
        }
        push_link(rt, f);
        COPY(f,h);
        NEXT(e);
    }
    GLUE(f,g);
    COPY(f,v->level_edge[level]);

    if (rt->highest_level == level)
        rt->high_site = v;
}

/* inserts v into e */
void flip2_4(Triangulation *rt, Site *v, Edge e, int level)
{
    Edge f,g,h;
    int i;

    g.triangle = NULL;
    h.triangle = NULL;

    e.triangle->level = -1;
    NEXT(e);
    for (i=0;i<4;i++) {
        create_triangle(rt, &f, level);
        SET_ALL(f,ORG(e),DEST(e),v);
        SWITCH_GLUE(e,f);
        COPY(f,DEST(f)->level_edge[level]);
        push_link(rt, f);
        if (i==0) {
            COPY(f,g);
            PREV(g);
            NEXT(f);
        } else {
            PREV(f);
            GLUE(f,h);
            PREV(f);
        }
        push_link(rt, f);
        COPY(f,h);
        NEXT(e);
        if (i==1) {
            SWITCH(e);
            e.triangle->level = -1;
            NEXT(e);
        }
    }
    GLUE(f,g);
    COPY(f,v->level_edge[level]);

    if (rt->highest_level == level)
        rt->high_site = v;
}

void flip2_2(Triangulation *rt, Edge e, int level)
{
    Edge f,g,h;
    int i;

    TO_SWITCHED(e,f);
    create_triangle(rt, &g, level);
    create_triangle(rt, &h, level);
    SET_ALL(g, APEX(e), APEX(f), DEST(e));
    SET_ALL(h, APEX(f), APEX(e), DEST(f));
    GLUE(g,h);

    for (i=0;i<2;i++) {
        NEXT(g);
        NEXT(h);
        PREV(e);
        PREV(f);
        if (i==0) { SWITCH_GLUE(e,h); } else { SWITCH_GLUE(f,h); }
        if (i==0) { SWITCH_GLUE(f,g); } else { SWITCH_GLUE(e,g); }
        COPY(g,DEST(g)->level_edge[level]);
        COPY(h,DEST(h)->level_edge[level]);
        push_link(rt, g);
        push_link(rt, h);
    }
    e.triangle->level = -1;
    f.triangle->level = -1;
}

/* removes the origin of e */
void flip3_1(Triangulation *rt, Edge e, int level)
{
    Edge f;
    int i;
    Site *v = ORG(e);
    v->level = -1;

    create_triangle(rt, &f, level);
    for (i=0;i<3;i++) {
        NEXT(e);
        SET_ORG(f,ORG(e));
        SWITCH_GLUE(e,f);
        e.triangle->level = -1;
        push_link(rt,f);
        NEXT(e);
        SWITCH(e);
        NEXT(f);
    }
    COPY(f,DEST(f)->level_edge[level]);
    NEXT(f);
    COPY(f,DEST(f)->level_edge[level]);
    NEXT(f);
    COPY(f,DEST(f)->level_edge[level]);
    COPY(f,v->level_edge[level]);

    if (level==rt->highest_level && v==rt->high_site) 
        rt->high_site = DEST(f);
}

/* removes the origin of e */
void flip4_2(Triangulation *rt, Edge e, int level)
{
    Edge f,u,l;
    Site *v = ORG(e);
    int i;

    v->level = -1;

    TO_PREV(e,f);
    SWITCH(f);
    NEXT(f);
    SWITCH(e);
    PREV(e);
    create_triangle(rt, &u, level);
    create_triangle(rt, &l, level);
    SET_ORG(u,ORG(e));
    SET_ORG(l,ORG(f));
    for (i=0;i<2;i++) {
        SET_DEST(u,DEST(e));
        SET_DEST(l,DEST(f));
        SWITCH_GLUE(e,u);
        SWITCH_GLUE(f,l);
        push_link(rt,u);
        push_link(rt,l);
        e.triangle->level = -1;
        f.triangle->level = -1;
        NEXT(e);
        SWITCH(e);
        NEXT(e);
        NEXT(f);
        SWITCH(f);
        NEXT(f);
        NEXT(u);
        NEXT(l);
    }
    GLUE(u,l);
    push_link(rt,u);
    COPY(u,DEST(u)->level_edge[level]);
    COPY(l,DEST(l)->level_edge[level]);
    NEXT(u);
    NEXT(l);
    COPY(u,DEST(u)->level_edge[level]);
    COPY(l,DEST(l)->level_edge[level]);
    COPY(l,v->level_edge[level]);

    if (level==rt->highest_level && v==rt->high_site) 
        rt->high_site = DEST(l);
}

/* checks if e is locally regular */
int locally_regular(Triangulation *rt, Edge e) 
{
    Edge f;
    Site *a, *b, *c, *d;

    TO_SWITCHED(e,f);

    a = APEX(e);
    b = ORG(e);
    c = DEST(e);
    d = APEX(f);

    if (INF(rt,b) && INF(rt,c))
        return 1;
    if (INF(rt,a) && INF(rt,d))
        return 1;
    if (INF(rt,c) || INF(rt,b))
        return 0;
    if (INF(rt,a) || INF(rt,d))
        return 1;

    double a1, a2, b1, b2, c1, c2;
    a1 = b->x-a->x;
    a2 = c->x-b->x;
    b1 = b->y-a->y;
    b2 = c->y-b->y;
    c1 = a->x*a->x+a->y*a->y-b->x*b->x-b->y*b->y-a->w+b->w;
    c2 = b->x*b->x+b->y*b->y-c->x*c->x-c->y*c->y-b->w+c->w;

    return (d->x-a->x)*(c2*b1-c1*b2)+(d->y-a->y)*(c1*a2-c2*a1)
           < (a1*b2-a2*b1)*(d->x*d->x+d->y*d->y-a->x*a->x-a->y*a->y+a->w-d->w);

    /*
    mpfr_t a1,a2,b1,b2,c1,c2,c3,m1,m2,acc;
    int s;

    mpfr_inits2(512,a1,a2,b1,b2,c1,c2,c3,m1,m2,acc,(mpfr_ptr)0); 
    mpfr_set_d(a1,b->x,MPFR_RNDN);
    mpfr_sub_d(a1,a1,a->x,MPFR_RNDN);
    mpfr_set_d(a2,c->x,MPFR_RNDN);
    mpfr_sub_d(a2,a2,b->x,MPFR_RNDN);
    mpfr_set_d(b1,b->y,MPFR_RNDN);
    mpfr_sub_d(b1,b1,a->y,MPFR_RNDN);
    mpfr_set_d(b2,c->y,MPFR_RNDN);
    mpfr_sub_d(b2,b2,b->y,MPFR_RNDN);

    // c1 = a->x*a->x+a->y*a->y-b->x*b->x-b->y*b->y-a->w+b->w
    mpfr_set_d(c1,a->x,MPFR_RNDN);
    mpfr_mul(c1,c1,c1,MPFR_RNDN);
    mpfr_set_d(m1,a->y,MPFR_RNDN);
    mpfr_mul(m1,m1,m1,MPFR_RNDN);
    mpfr_add(c1,c1,m1,MPFR_RNDN);
    mpfr_set_d(m1,b->x,MPFR_RNDN);
    mpfr_mul(m1,m1,m1,MPFR_RNDN);
    mpfr_sub(c1,c1,m1,MPFR_RNDN);
    mpfr_set_d(m1,b->y,MPFR_RNDN);
    mpfr_mul(m1,m1,m1,MPFR_RNDN);
    mpfr_sub(c1,c1,m1,MPFR_RNDN);
    mpfr_sub_d(c1,c1,a->w,MPFR_RNDN);
    mpfr_add_d(c1,c1,b->w,MPFR_RNDN);

    // c2 = b->x*b->x+b->y*b->y-c->x*c->x-c->y*c->y-b->w+c->w
    mpfr_set_d(c2,b->x,MPFR_RNDN);
    mpfr_mul(c2,c2,c2,MPFR_RNDN);
    mpfr_set_d(m1,b->y,MPFR_RNDN);
    mpfr_mul(m1,m1,m1,MPFR_RNDN);
    mpfr_add(c2,c2,m1,MPFR_RNDN);
    mpfr_set_d(m1,c->x,MPFR_RNDN);
    mpfr_mul(m1,m1,m1,MPFR_RNDN);
    mpfr_sub(c2,c2,m1,MPFR_RNDN);
    mpfr_set_d(m1,c->y,MPFR_RNDN);
    mpfr_mul(m1,m1,m1,MPFR_RNDN);
    mpfr_sub(c2,c2,m1,MPFR_RNDN);
    mpfr_sub_d(c2,c2,b->w,MPFR_RNDN);
    mpfr_add_d(c2,c2,c->w,MPFR_RNDN);

    // c3 = d->x*d->x+d->y*d->y-a->x*a->x-a->y*a->y+a->w-d->w
    mpfr_set_d(c3,d->x,MPFR_RNDN);
    mpfr_mul(c3,c3,c3,MPFR_RNDN);
    mpfr_set_d(m1,d->y,MPFR_RNDN);
    mpfr_mul(m1,m1,m1,MPFR_RNDN);
    mpfr_add(c3,c3,m1,MPFR_RNDN);
    mpfr_set_d(m1,a->x,MPFR_RNDN);
    mpfr_mul(m1,m1,m1,MPFR_RNDN);
    mpfr_sub(c3,c3,m1,MPFR_RNDN);
    mpfr_set_d(m1,a->y,MPFR_RNDN);
    mpfr_mul(m1,m1,m1,MPFR_RNDN);
    mpfr_sub(c3,c3,m1,MPFR_RNDN);
    mpfr_sub_d(c3,c3,d->w,MPFR_RNDN);
    mpfr_add_d(c3,c3,a->w,MPFR_RNDN);

    // return (a1*b2-a2*b1)*c3
    //        - (d->x-a->x)*(c2*b1-c1*b2)-(d->y-a->y)*(c1*a2-c2*a1) > 0
    mpfr_mul(m1,a1,b2,MPFR_RNDN);
    mpfr_mul(m2,a2,b1,MPFR_RNDN);
    mpfr_sub(m1,m1,m2,MPFR_RNDN);
    mpfr_mul(acc,m1,c3,MPFR_RNDN);

    mpfr_mul(m1,c2,b1,MPFR_RNDN);
    mpfr_mul(m2,c1,b2,MPFR_RNDN);
    mpfr_sub(m1,m1,m2,MPFR_RNDN);
    mpfr_set_d(m2,d->x,MPFR_RNDN);
    mpfr_sub_d(m2,m2,a->x,MPFR_RNDN);
    mpfr_mul(m1,m1,m2,MPFR_RNDN);
    mpfr_sub(acc,acc,m1,MPFR_RNDN);

    mpfr_mul(m1,c1,a2,MPFR_RNDN);
    mpfr_mul(m2,c2,a1,MPFR_RNDN);
    mpfr_sub(m1,m1,m2,MPFR_RNDN);
    mpfr_set_d(m2,d->y,MPFR_RNDN);
    mpfr_sub_d(m2,m2,a->y,MPFR_RNDN);
    mpfr_mul(m1,m1,m2,MPFR_RNDN);
    mpfr_sub(acc,acc,m1,MPFR_RNDN);

    s = mpfr_sgn(acc);
    mpfr_clears(a1,a2,b1,b2,c1,c2,c3,m1,m2,acc,(mpfr_ptr)0); 
    return s>0;
    */
}

/* flips e if possible */
void flip(Triangulation *rt, Edge e, int level) {

    Edge f;
    TO_SWITCHED(e,f);
    int ori_base_e = ORIENTATION(APEX(e),ORG(e),APEX(f));
    int ori_base_f = ORIENTATION(APEX(f),ORG(f),APEX(e));

    // 2-2-flippable?
    if (ori_base_e > 0 && ori_base_f > 0) {
        flip2_2(rt, e, level);
        return;
    }

    if (INF(rt,ORG(e)) || INF(rt,DEST(e)))
        return;

    // 3-1-flippable?
    int order_e = order(e);
    if (order_e == 3 && ori_base_e < 0) {
        flip3_1(rt, e, level);
        return;
    }

    int order_f = order(f);
    if (order_f == 3 && ori_base_f < 0) {
        flip3_1(rt, f, level);
        return;
    }

    // 4-2-flippable?
    if (order_e == 4 && ori_base_e == 0) {
        flip4_2(rt, e, level);
        return;
    }

    if (order_f == 4 && ori_base_f == 0) {
        flip4_2(rt, f, level);
        return;
    }
}

/* Initializes a new level of triangulation, see [2] */
void create_new_level(Triangulation *rt, Site *s)
{
    Edge e, f;
    int level;

    level = ++(rt->highest_level);
    ++(rt->s_inf1.level);
    ++(rt->s_inf2.level);
    ++(rt->s_inf3.level);

    rt->t_inf = 
        (Triangle) { { &rt->s_inf1, &rt->s_inf2, &rt->s_inf3}, 
                     { &rt->t_inf, &rt->t_inf, &rt->t_inf }, 
                     { 0, 0, 0 }, 
                     level };

    f.triangle = &rt->t_inf;
    f.shift = 0;

    // create bounding triangle
    create_triangle(rt, &e, level);
    SET_ALL(e, &rt->s_inf1, &rt->s_inf2, &rt->s_inf3);
    COPY(e, rt->s_inf2.level_edge[level]);
    GLUE(e,f);
    NEXT(e);
    PREV(f);
    COPY(e, rt->s_inf3.level_edge[level]);
    GLUE(e,f);
    NEXT(e);
    PREV(f);
    COPY(e, rt->s_inf1.level_edge[level]);
    GLUE(e,f);

    // add s as first point
    ++(s->level);
    flip1_3(rt, s, e, level);
    rt->stack_i = 0;
}

/* Inserts vertex v into triangle t, assuming that v lies inside t */
void insert_site(Triangulation *rt, Site *s, Triangle *t, int level, int is_on_border)
{
    Edge e;
    e.triangle = t;

    // insert v
    if (is_on_border >= 0) {
        e.shift = is_on_border;
        flip2_4(rt, s, e, level);
    } else {
        e.shift = 0;
        flip1_3(rt, s, e, level);
    }
    // flip irregular edges
    Edge *f;
    //int debug0 = 0;
    while ((f = pop_link(rt)) != NULL) {
        /*
        ++debug0;
        if (debug0>10000) {
            error("debug0");
            break;
        }
        */
        // only flip if f is still an edge of the current triangulation
        if (f->triangle->level==level && !locally_regular(rt, *f))
            flip(rt, *f, level);
    }
}

/* Adds site v to the triangulation. If the site turs out to be redundant (or hidden, i.e.
 * it's power cell is empty) it will immediately be removed during the flipping process.
 * To find the triangle to which v must be added, line search is used, starting at site 'start'. */
void add_site(Triangulation *rt, Site *s, Site *start, int level)
{
    Triangle *t;
    Edge e;
    int is_on_border, ori1, ori2, ori3;

    COPY(start->level_edge[level], e);
    //int debug1=0,debug2=0,debug3=0; 

    // 'start' might be redundant in this level but still active
    // at the level above, i.e. e.triangle could be removed;
    // the level edges are updated during deletion so we can
    // update e here until e.triangle is part of the current 
    // triangulation
    while (e.triangle->level==-1) {
        /*
        ++debug2;
        if (debug2>10000) {
            error("debug2");
            return;
        }
        */
        COPY(DEST(e)->level_edge[level], e);
    }
    start = DEST(e);

    // find an edge of the current triangulation which intersects the edge (start->v)
    while (1) {
        /*
        ++debug1;
        if (debug1>10000) {
            error("debug1");
            return;
        }
        */

        ori1 = ORIENTATION(DEST(e),ORG(e),s);
        if (ori1 > 0) {
            SWITCH(e);
            PREV(e);
            continue;
        }

        ori3 = ORIENTATION(APEX(e),DEST(e),s);
        if (ori3 > 0) {
            NEXT(e);
            SWITCH(e);
            continue;
        }

        ori2 = ORIENTATION(ORG(e),APEX(e),s);
        if (ori2 > 0) {
            PREV(e);
            SWITCH(e);
            break;
        }

        // now e.triangle must contain s
        if (ori1 == 0) 
            is_on_border = e.shift;
        else if (ori2 == 0) 
            is_on_border = (e.shift+2)%3;
        else if (ori3 == 0) 
            is_on_border = (e.shift+1)%3;
        else
            is_on_border = -1;

        t = e.triangle;
        goto found;
    }

    // walk along start->s to find the triangle containing s
    while (1) {
        /*
        ++debug3;
        if (debug3>10000) {
            error("debug3");
            return;
        }
        */

        ori1 = ORIENTATION(DEST(e),APEX(e),s);
        ori2 = ORIENTATION(APEX(e),ORG(e),s);
        if ((ori1 >= 0) &&
            (ori2 >= 0)) {
            if (ori1 == 0)
                is_on_border = (e.shift+2)%3;
            else if (ori2 == 0)
                is_on_border = (e.shift+1)%3;
            else
                is_on_border = -1;

            t = e.triangle;
            break;
        } else if (ORIENTATION(start,APEX(e),s) > 0) {
            PREV(e);
            SWITCH(e);
        } else {
            NEXT(e);
            SWITCH(e);
        }
    }

    found:
    if (level > 0) {
        //TODO: take near site?
        add_site(rt, s, APEX(e), level-1);
    }

    if (level == s->level) {

        insert_site(rt, s, t, level, is_on_border);

        // increment the level of s with a probability of 1/40 (see [2])
        if (s->level > -1 && runif(0,1)<1.0/40.0 && s->level < MAX_LEVEL-1) {
            if (level==rt->highest_level) {
                do {
                    create_new_level(rt, s);
                } while (runif(0,1)<1.0/40.0 && s->level < MAX_LEVEL-1);
            } else {
                ++(s->level);
            }
        }
    }
}

/* Computest the regular triangulation of points (x,y,w). It is necessary to call
 * init_triangulation(rt) once; this function can then be used repeatedly; call
 * free_memory(rt) when finished.
 * Input data should not contain dublicates and must lie in the bounding box [rx1,rx2]x[ry1,ry2], 
 * rx1<rx2, ry1<ry2. */
void triangulate(Triangulation *rt, int sites, double *x, double *y, double *w, 
                 double rx1, double ry1, double rx2, double ry2,
                 double pert)
{
    Site *s;
    int i;

    // reset triangulation
    rt->current_block = &rt->root;
    rt->root->i = 0; 
    rt->stack_i = 0;

    rt->highest_level = -1;
    rt->s_inf1.level = -1;
    rt->s_inf2.level = -1;
    rt->s_inf3.level = -1;

    if (sites<1)
        return;

    if (sites > rt->size) {
        rt->sites = REALLOC(rt->sites, sites, Site);
        rt->level_edges = REALLOC(rt->level_edges, sites*MAX_LEVEL, Edge);
    }
    rt->size = sites;

    // set up the bounding triangle at 'infinity'
    double xi = 3e10*(rx2-rx1 > ry2-ry1 ? rx2-rx1 : ry2-ry1);
    rt->s_inf1.x = rx1-xi;
    rt->s_inf1.y = ry1-xi;
    rt->s_inf1.w = 0;
    rt->s_inf2.x = rx1+xi;
    rt->s_inf2.y = ry1;
    rt->s_inf2.w = 0;
    rt->s_inf3.x = rx1;
    rt->s_inf3.y = ry1+xi;
    rt->s_inf3.w = 0;

    GetRNGstate();
    // add sites one by one
    for (i=0; i<sites; i++) {
        s = &rt->sites[i];
        s->x = x[i] + runif(-pert,pert);
        s->y = y[i] + runif(-pert,pert);
        s->w = w[i];
        s->level = 0;
        s->level_edge = &rt->level_edges[i*MAX_LEVEL];
        if (i==0) {
            s->level = -1;
            create_new_level(rt, s);
        } else {
            add_site(rt, s, rt->high_site, rt->highest_level);
        }
    }
    PutRNGstate();

    // DEBUG
    /*
    int n,j;
    double a=0.0, cx[10000], cy[10000], ox, oy;
    Edge e;

    // check sum of triangle areas
    for (i=0; i<sites; i++) { 
        power_cell(&n,cx,cy,rt,&rt->sites[i],rx1,ry1,rx2,ry2);
        a += polygon_area(n,cx,cy);
    }

    if (a!=(rx2-rx1)*(ry2-ry1)) {
        printf("pd area diff = %lf\n",a-(rx2-rx1)*(ry2-ry1));
    }

    // check regularity
    for (i=0; i<sites; i++) { 
        e = rt->sites[i].level_edge[0];

        if (rt->sites[i].level == -1)
            continue;
        if (INF(rt,DEST(e)) || INF(rt,ORG(e)) || INF(rt,APEX(e)))
            continue;

        orthogonal_center(&ox,&oy,ORG(e),DEST(e),APEX(e));

        for (j=0; j<sites; j++) { 
            if (APEX(e)==&rt->sites[j] ||
                DEST(e)==&rt->sites[j] ||
                ORG(e)==&rt->sites[j])
                continue;

            if (pdist(ox,oy,&rt->sites[j]) <= pdist(ox,oy,ORG(e)) ||
                pdist(ox,oy,&rt->sites[j]) <= pdist(ox,oy,DEST(e)) ||
                pdist(ox,oy,&rt->sites[j]) <= pdist(ox,oy,APEX(e))) {
                printf("irregular triangle found\n");
                goto end;
            }
        }
    }
    end:;
    */
    // END DEBUG
}

/* Initialization */
void init_triangulation(Triangulation *rt)
{
    rt->s_inf1 = (Site) { 0, 0, 0, -1, NULL }; // (-inf,-inf)
    rt->s_inf2 = (Site) { 0, 0, 0, -1, NULL }; // (+inf,0)
    rt->s_inf3 = (Site) { 0, 0, 0, -1, NULL }; // (0,+inf)
    rt->s_inf1.level_edge = MALLOC(MAX_LEVEL, Edge);
    rt->s_inf2.level_edge = MALLOC(MAX_LEVEL, Edge);
    rt->s_inf3.level_edge = MALLOC(MAX_LEVEL, Edge);
    rt->sites = MALLOC(1, Site);
    rt->root = MALLOC(1, struct MBlock);
    rt->root->triangles = MALLOC(MEMORY_BLOCKSIZE, Triangle);
    rt->link_facets = MALLOC(LINK_FACETS_BLOCKSIZE, Edge);
    rt->level_edges = MALLOC(MAX_LEVEL, Edge);

    rt->size = 0;

    rt->root->i = 0;
    rt->root->prev = NULL;
    rt->root->next = NULL;
    rt->current_block = &rt->root;

    rt->stack_i = 0;
}

/* Free memory */
void free_triangulation(Triangulation *rt)
{
    struct MBlock *b = *rt->current_block;

    while (b->next != NULL) {
        b = b->next;
    }

    while (b != rt->root) {
        FREE(b->triangles);
        b = b->prev;
        FREE(b->next);
    }

    FREE(rt->root->triangles);
    FREE(rt->root);
    FREE(rt->sites);
    FREE(rt->link_facets);
    FREE(rt->level_edges);
    FREE(rt->s_inf1.level_edge);
    FREE(rt->s_inf2.level_edge);
    FREE(rt->s_inf3.level_edge);
}

/* Computes the intersection of (x1,y1)->(x2,y2) and the line given by ax+by+c = 0,
 * assuming it exists */
void isect_edge_line(double *x, double *y, double x1, double y1, double x2, double y2,
                     double a, double b, double c)
{
    double A = -(c+a*x1+b*y1)/(a*(x2-x1)+b*(y2-y1));
    *x = A*(x2-x1)+x1;
    *y = A*(y2-y1)+y1;
}

/* Intersects the polygon (n,x,y) with the halfspace given by ax+by+c <= 0 */
void isect_polygon_halfspace(int *n, double *x, double *y, double a, double b, double c)
{
    int i, j=0;

    x[*n] = x[0];
    y[*n] = y[0];
    x[*n+1] = x[1];
    y[*n+1] = y[1];
    x[*n+2] = x[2];
    y[*n+2] = y[2];
    // iterate over edges
    for (i=2;i<*n+2;i++) {
        // remove origin vertex of the edge if it is not part of the halfspace
        if (a*x[i]+b*y[i]+c > 0) {
            // if the destination vertex is part of the halfspace,
            // it might be necessary to add an intersection vertex
            if (a*x[i+1]+b*y[i+1]+c < 0) {
                isect_edge_line(&x[j],&y[j],x[i],y[i],x[i+1],y[i+1],a,b,c);
                ++j;
            }
        } else {
            x[j] = x[i];
            y[j] = y[i];
            ++j;
            // remove destination vertex if it is not part of the halfspace
            if (a*x[i+1]+b*y[i+1]+c > 0) {
                // since in this situation the origin vertex is part of
                // the halfspace and the destination vertex is not,
                // it might be necessary to add an intersection vertex
                if (a*x[i]+b*y[i]+c != 0) {
                    isect_edge_line(&x[j],&y[j],x[i],y[i],x[i+1],y[i+1],a,b,c);
                    ++j;
                }
            }
        }
    }
    *n = j;
}

/* Computes the power cell of s intersected with the rectangle [rx1,rx2]x[ry1,ry2] */
void power_cell(int *n, double *x, double *y, Triangulation *rt, Site *s, 
                double rx1, double ry1, double rx2, double ry2)
{
    Edge e = s->level_edge[0];
    Site *t, *start = ORG(e);
    double a,b,c;

    if (s->level == -1) {
        *n = 0;
        return;
    }

    // initialize rectange polygon
    *n = 4;
    x[0] = rx1;
    y[0] = ry1;
    x[1] = rx2;
    y[1] = ry1;
    x[2] = rx2;
    y[2] = ry2;
    x[3] = rx1;
    y[3] = ry2;

    do {
        if (!INF(rt,DEST(e))) {
            t = ORG(e);
            a = 2*(t->x-s->x);
            b = 2*(t->y-s->y);
            c = s->x*s->x + s->y*s->y - t->x*t->x - t->y*t->y + t->w - s->w;
            isect_polygon_halfspace(n,x,y,a,b,c);
        }
        SWITCH(e);
        PREV(e);
    } while(ORG(e) != start);
}
#endif
