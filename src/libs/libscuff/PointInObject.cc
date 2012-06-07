/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* Fast tree-based point-in-object testing, by Steven G. Johnson (3/2012).

   We test whether a point (x,y,z) is in an object by counting the number
   of triangular panels intersected by the ray from (x,y,z) to (x,y,-infinty):
   this is EVEN if the point IS NOT in the object, and ODD if it IS.
   (This works even for non-convex objects and objects with holes.)

   Testing ray intersections boils down to (quickly) finding which triangles
   contain (x,y) when projected onto the z=0 plane.  We do this by a kd-tree
   of triangles: a tree structure which recursively bisects the xy plane
   by x=# or y=# lines (for some #). 

   The preprocessing time is O(N log N), where N is the number of panels,
   and the storage is O(N).  Subsequent inclusion tests should be O(log N).
*/

#include <stdlib.h>
#include <string.h>

#include "libscuff.h"

namespace scuff {

/*************************************************************************/

/* Test whether (x,y) is in the triangle with vertices (xt[i],yt[i]),
   sorted in ascending order of yt[i].

   Based on the PNPOLY algorithm by W. Randolph Franklin
   (http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html),
   simplified to the case of a triangle, and modified to exploit y-sorted
   vertices to avoid divisions.   The inequalities are carefully ordered
   to handle points on edges in a consistent fashion for adjacent triangles.

   Note that vertices are stored in single precision (Kahan-style "fiducial"
   representation), but computations are done in double precision) */

#define PAST_EDGE(x,y,x0,y0,x1,y1) (((y0>y) != (y1>y)) &&		\
				    (x-x0)*(y1-y0) < (x1-x0)*(y-y0))
static int point_in_tri(double x, double y,
			const float xt[3], const float yt[3])
{
  /* true if (x,y) is to the right of an odd number of edges */
  return (PAST_EDGE(x,y,xt[0],yt[0],xt[1],yt[1]) +
	  PAST_EDGE(x,y,xt[1],yt[1],xt[2],yt[2]) +
	  PAST_EDGE(x,y,xt[0],yt[0],xt[2],yt[2])) % 2;
}

/* Test whether (x,y,z) is above the given triangle: whether the
   ray from (x,y,z) to (x,y,-infinity) intersects the triangle. */
static int point_above_tri(double x, double y, double z,
		    const float xt[3], const float yt[3], const float zt[3])
{
  if (point_in_tri(x,y,xt,yt)) {
    double d1x = xt[1] - xt[0];
    double d1y = yt[1] - yt[0];
    double d1z = zt[1] - zt[0];
    double d2x = xt[2] - xt[0];
    double d2y = yt[2] - yt[0];
    double d2z = zt[2] - zt[0];
    double cx = d1y * d2z - d1z * d2y; /* (d1 x d2).x */
    double cy = d1z * d2x - d1x * d2z; /* (d1 x d2).y */
    double cz = d1x * d2y - d1y * d2x; /* (d1 x d2).z */
    return ((x - xt[0]) * cx + (y - yt[0]) * cy + (z - zt[0]) * cz) * cz > 0;
  }
  return 0;
}

/*************************************************************************/
/* Fast (linear-time) array partitioning */

/* swap size bytes between a_ and b_ */
static void swap(void *a_, void *b_, size_t size)
{
     if (a_ == b_) return;
     {
          size_t i, nlong = size / sizeof(long);
          long *a = (long *) a_, *b = (long *) b_;
          for (i = 0; i < nlong; ++i) {
               long c = a[i];
               a[i] = b[i];
               b[i] = c;
          }
	  a_ = (void*) (a + nlong);
	  b_ = (void*) (b + nlong);
     }
     {
          size_t i;
          char *a = (char *) a_, *b = (char *) b_;
          size = size % sizeof(long);
          for (i = 0; i < size; ++i) {
               char c = a[i];
               a[i] = b[i];
               b[i] = c;
          }
     }
}

/* like qsort_r, but reorders the array so that the first WHICH (> 0, < NMEMB)
   elements are <= (according to the given comparison function) than
   the last NMEMB-WHICH elements, but the elements are not otherwise
   necessarily in order.  This requires O(nmemb * size) time.  */
static void partition(void *base_, size_t nmemb, size_t size, size_t which,
		      void *thunk,
		      int (*compar)(void *, const void *, const void *))
{
     char *base = (char *) base_;
     if (nmemb < 10) { /* just use O(nmemb^2) sort for small enough nmemb */
          size_t idx[10] = {0,1,2,3,4,5,6,7,8,9};
	  size_t i, j;
	  for (i = 0; i+1 < nmemb; ++i)
	       for (j = i+1; j < nmemb; ++j)
		    if (compar(thunk, base+idx[i]*size, base+idx[j]*size) > 0){
			 int idx_ = idx[i];
			 idx[i] = idx[j];
			 idx[j] = idx_;
		    }
	  for (i = 0; i < which; ++i)
	       if (idx[i] >= which) {
		    for (j = which; j < size && idx[j] >= which; ++j) ;
		    swap(base+i*size, base+j*size, size);
		    idx[j] = idx[i];
	       }
     }
     else {
	  size_t i, pivot, npart;
	  /* pick median of first/middle/last elements as pivot */
	  {
	       const char *a = base, *b = base + (nmemb/2)*size, 
		    *c = base + (nmemb-1)*size;
	       pivot = compar(thunk,a,b) < 0
		    ? (compar(thunk,b,c) < 0 ? nmemb/2 :
		       (compar(thunk,a,c) < 0 ? nmemb-1 : 0))
		    : (compar(thunk,a,c) < 0 ? 0 :
		       (compar(thunk,b,c) < 0 ? nmemb-1 : nmemb/2));
	  }
	  /* partition array */
	  swap(base + pivot*size, base + (nmemb-1) * size, size);
	  pivot = (nmemb - 1) * size;
	  for (i = npart = 0; i < nmemb-1; ++i)
	       if (compar(thunk, base+i*size, base+pivot) <= 0)
		    swap(base+i*size, base+(npart++)*size, size);
	  swap(base+npart*size, base+pivot, size);
	  /* recursive partition one of the two partitions */
	  if (which < npart)
	       partition(base, npart, size, which, thunk, compar);
	  else if (which > ++npart) {
	       partition(base+npart*size, nmemb-npart, size, 
			 which-npart, thunk, compar);
	  }
     }
}

/***********************************************************************/

typedef struct {
  float bmin[2], bmax[2]; /* corners of bounding box */
  float x[3], y[3], z[3]; /* corners of triangle, sorted by y[i] */
} boxtri;

struct kdtri_s {
  int dim; /* which dimension to split on (0 or 1) */
  double div; /* divide into boxtri's <= div and >= div */
  struct kdtri_s *le, *gt; /* subtrees */
  size_t n;
  boxtri *B; /* array of n boxtries in this leaf if le == gt == 0 */
  float bmin[3], bmax[3]; /* bounding box of B */
};

/***********************************************************************/

/* compare two boxtries A and B along the dim direction.  A </> B
   if the DIM components of all the corners of A are </> the DIM components
   of all corners of B, otherwise A == B for purposes of this ordering.
   Note that this is only a partial ordering, since "==" is not transitive. */
static int boxtri_compar(void *dim_, const void *a_, const void *b_)
{
     int dim = *((int *) dim_);
     const boxtri *a = (const boxtri *) a_;
     const boxtri *b = (const boxtri *) b_;
     return (a->bmax[dim] < b->bmin[dim] ? -1
	     : (b->bmax[dim] < a->bmin[dim] ? +1 : 0));
}

static void get_partition(kdtri t, size_t *nle, size_t *ngt) {
     size_t i, n = t->n;
     double div;
     boxtri *B = t->B;
     int dim = t->dim;

     partition(B, n, sizeof(boxtri), n/2, &dim, boxtri_compar);
     div = t->div = 0.5 * (B[n/2].bmin[dim] + B[n/2].bmax[dim]);

     *nle = *ngt = 0;
     for (i = 0; i < n; ++i) {
	  if (B[i].bmin[dim] <= div)
	       *nle += 1;
	  if (B[i].bmax[dim] > div)
	       *ngt += 1;
     }
}

static size_t limax2(size_t a, size_t b) { return a > b ? a : b; }

static void best_partition(kdtri t, size_t *nle, size_t *ngt)
{
     int best_dim;
     size_t best_red, red;

     t->dim = 0;
     get_partition(t, nle, ngt);
     best_dim = t->dim;
     best_red = limax2(*nle, *ngt);

     t->dim = 1;
     get_partition(t, nle, ngt);
     if ((red = limax2(*nle, *ngt)) < best_red) {
	  best_dim = t->dim;
	  best_red = red;
     }
     
     if (best_dim != t->dim) {
	  t->dim = best_dim;
	  get_partition(t, nle, ngt);
     }
}

#define MIN2(a,b) (a < b ? a : b)
#define MAX2(a,b) (a > b ? a : b)
#define MIN3(a,b,c) (a < b ? MIN2(a,c) : MIN2(b,c))
#define MAX3(a,b,c) (a > b ? MAX2(a,c) : MAX2(b,c))

/* compute bounding box */
static void kdtri_compute_bounds(kdtri t)
{
  if (!t) return;
  if (t->n > 0) {
    size_t i, n = t->n;

    t->bmin[0] = t->B[0].bmin[0];
    t->bmin[1] = t->B[0].bmin[1];
    t->bmin[2] = MIN3(t->B[0].z[0],t->B[0].z[1],t->B[0].z[2]);

    t->bmax[0] = t->B[0].bmax[0];
    t->bmax[1] = t->B[0].bmax[1];
    t->bmax[2] = MAX3(t->B[0].z[0],t->B[0].z[1],t->B[0].z[2]);

    for (i = 1; i < n; ++i) { /* union of bounding boxes */
      boxtri *b = t->B + i;
      float z;

      t->bmin[0] = MIN2(t->bmin[0], b->bmin[0]);
      t->bmin[1] = MIN2(t->bmin[1], b->bmin[1]);
      z = MIN3(b->z[0],b->z[1],b->z[2]);
      t->bmin[2] = MIN2(t->bmin[2], z);

      t->bmax[0] = MAX2(t->bmax[0], b->bmax[0]);
      t->bmax[1] = MAX2(t->bmax[1], b->bmax[1]);
      z = MAX3(b->z[0],b->z[1],b->z[2]);
      t->bmax[2] = MAX2(t->bmax[2], z);
    }
  }
  else
    t->bmax[0] = t->bmax[1] = t->bmax[2] =
      t->bmax[0] = t->bmax[1] = t->bmax[2] = 0;
}

/* create a new t tree with a boxtri array pointed to by B */
static kdtri kdtri_create(boxtri *B, size_t n)
{
     kdtri t = (kdtri) malloc(sizeof(struct kdtri_s));
     if (t) {
       t->B = B;
       t->n = n;
       t->le = t->gt = NULL;
       t->dim = -1;
       t->div = 0;
       kdtri_compute_bounds(t);
     }
     return t;
}

/* deallocate tree, including the boxtri arrays */
void kdtri_destroy(kdtri t)
{
     if (t) {
	  kdtri_destroy(t->le);
	  kdtri_destroy(t->gt);
	  free(t->B);
	  free(t);
     }
}

/* partition the tree */
static int kdtri_partition(kdtri t)
{
     size_t n;
     if (t && (n = t->n) > 10) { /* only partition large leaves */
	  boxtri *B0;
	  int dim;
	  double div;
	  size_t nle, ngt;

	  best_partition(t, &nle, &ngt);
	  dim = t->dim;
	  div = t->div;
	  B0 = t->B;

	  if (limax2(nle, ngt) > (3*n)/4)
	       return EXIT_SUCCESS; /* subdivision doesn't gain us enough */

	  if (nle > 0) {
	       size_t i, ile;
	       boxtri *B = (boxtri *) malloc(sizeof(boxtri) * nle);
	       if (!B) return EXIT_FAILURE;
	       for (i = ile = 0; i < n; ++i)
		    if (B0[i].bmin[dim] <= div)
			 B[ile++] = B0[i];
	       t->le = kdtri_create(B, nle);
	       if (!(t->le) || kdtri_partition(t->le) != EXIT_SUCCESS)
		 return EXIT_FAILURE;
	  }
	  if (ngt > 0) {
	       size_t i, igt;
	       boxtri *B = (boxtri *) malloc(sizeof(boxtri) * ngt);
	       if (!B) return EXIT_FAILURE;
	       for (i = igt = 0; i < n; ++i)
		    if (B0[i].bmax[dim] > div)
			 B[igt++] = B0[i];
	       t->gt = kdtri_create(B, ngt);
	       if (!(t->gt) || kdtri_partition(t->gt) != EXIT_SUCCESS)
		 return EXIT_FAILURE;
	  }

	  free(t->B);
	  t->B = NULL;
     }
     return EXIT_SUCCESS;
}

/***********************************************************************/

/* some tree statistics: */
static unsigned imax2(unsigned a, unsigned b) { return (a > b ? a : b); }
unsigned kdtri_maxdepth(kdtri t)
{
     if (!t) return 0;
     return imax2(kdtri_maxdepth(t->le), kdtri_maxdepth(t->gt)) + 1;
}
size_t kdtri_maxleaf(kdtri t)
{
     if (!t) return 0;
     if (!t->le) return t->n;
     return imax2(kdtri_maxleaf(t->le), kdtri_maxleaf(t->gt));
}
double kdtri_meandepth(kdtri t)
{
     if (!t) return 0.0;
     return 1.0 + 0.5 * (kdtri_meandepth(t->le) + kdtri_meandepth(t->gt));
}
double kdtri_meanleaf(kdtri t)
{
     if (!t) return 0.0;
     if (!t->le) return 1.0 * t->n;
     return 0.5 * (kdtri_meanleaf(t->le) + kdtri_meanleaf(t->gt));
}

/***********************************************************************/

/* return true if bounding box of b contains p */
static int boxtri_contains(const boxtri *b, const double p[2])
{
     return (p[0] >= b->bmin[0] && p[0] <= b->bmax[0] &&
	     p[1] >= b->bmin[1] && p[1] <= b->bmax[1]);
}

#if 0
/* return pointer to an element in t that contains p, or
   NULL if none, where contains tells whether a point is in an element */
static boxtri *kdtri_find(kdtri t, const double p[2])
{
     if (!t) return NULL;
     if (t->le) { /* search subtrees */
	  if (p[t->dim] <= t->div)
	       return kdtri_find(t->le, p, contains);
	  else
	       return kdtri_find(t->gt, p, contains);
     }
     else { /* we are at a leaf node: search directly */
	  size_t i, n = t->n;
	  boxtri *B = t->B;
	  for (i = 0; i < n; ++i)
	       if (boxtri_contains(B+i, p)
		   && point_in_tri(p[0],p[1], B[i].x,B[i].y))
		    return B + i;
	  return NULL;
     }
}
#endif

/* return count of the number of triangles in t lying below p, i.e.
   intersecting the ray from p to (p[0],p[1],-infinity). */
static int kdtri_count_below(kdtri t, const double p[3])
{
  if (!t) return 0;
  if (t->le) { /* search subtrees */
    if (p[t->dim] <= t->div)
      return kdtri_count_below(t->le, p);
    else
      return kdtri_count_below(t->gt, p);
  }
  else { /* we are at a leaf node: search directly */
    size_t i, n = t->n;
    boxtri *B = t->B;
    int count = 0;
    for (i = 0; i < n; ++i)
      if (boxtri_contains(B+i, p)
	  && point_above_tri(p[0],p[1],p[2], B[i].x,B[i].y,B[i].z))
	++count;
    return count;
  }
}

/* Return true iff the object that was used to create kdtri
   contains the point p. */
static int kdtri_object_contains(kdtri t, const double p[3])
{
  if (!t ||
      p[0] < t->bmin[0] || p[0] > t->bmax[0] ||
      p[1] < t->bmin[1] || p[1] > t->bmax[1] ||
      p[2] < t->bmin[2] || p[2] > t->bmax[2])
    return 0;

  /* contains p iff an odd number of triangles lie below p */
  return kdtri_count_below(t, p) % 2;
}

/***********************************************************************/

/* Create a (tree-partitioned) kdtri object for the panels of O, using
   O(NumPanels) storage and O(NumPanels*log(NumPanels)) time.   Returns
   NULL if we ran out of memory. */
static kdtri kdtri_create_from_object(const RWGObject *O)
{
  int i;
  boxtri *B;
  kdtri t = NULL;

  if (!O || !O->NumPanels) goto done;
  if (!(B = (boxtri *) malloc(sizeof(boxtri) * O->NumPanels))) goto done;

  for (i = 0; i < O->NumPanels; ++i) {
    int j;
    for (j = 0; j < 3; ++j) {
      int k, vj = O->Panels[i]->VI[j];
      B[i].x[j] = O->Vertices[3*vj];
      B[i].y[j] = O->Vertices[3*vj+1];
      B[i].z[j] = O->Vertices[3*vj+2];
      for (k = 0; k < j; ++k) /* sort in y order */
	if (B[i].y[k] > B[i].y[j]) {
	  float x = B[i].x[k], y = B[i].y[k], z = B[i].z[k];
	  B[i].x[k] = B[i].x[j];
	  B[i].y[k] = B[i].y[j];
	  B[i].z[k] = B[i].z[j];
	  B[i].x[j] = x; B[i].y[j] = y; B[i].z[j] = z;
	}
    }

    /* compute bounding box */
    B[i].bmin[0] = MIN3(B[i].x[0], B[i].x[1], B[i].x[2]);
    B[i].bmin[1] = MIN3(B[i].y[0], B[i].y[1], B[i].y[2]);
    B[i].bmax[0] = MAX3(B[i].x[0], B[i].x[1], B[i].x[2]);
    B[i].bmax[1] = MAX3(B[i].y[0], B[i].y[1], B[i].y[2]);
  }

  if (!(t = kdtri_create(B, (size_t) O->NumPanels))) {
    free(B);
    goto done;
  }

  if (kdtri_partition(t) != EXIT_SUCCESS) {
    kdtri_destroy(t); t = NULL;
    goto done;
  }

 done:
  return t;
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void RWGObject::InitkdPanels(bool reinit, int LogLevel)
{
  if (kdPanels) {
    if (!reinit) return;
    kdtri_destroy(kdPanels); // destroy pre-existing tree
    kdPanels = NULL;
  }
  if (!NumPanels) return;
  kdPanels = kdtri_create_from_object(this);
  if (!kdPanels) ErrExit("out of memory when creating panel kd-tree for %s",
			 Label);

  if (LogLevel >= SCUFF_VERBOSELOGGING) {
    Log("kdPanels(%s): depth mean %g/max %u, leaf mean %g/max %lu", Label,
	kdtri_meandepth(kdPanels), kdtri_maxdepth(kdPanels),
	kdtri_meanleaf(kdPanels), (unsigned long) kdtri_maxleaf(kdPanels));
  }
}

/***********************************************************************/

/* return whether object contains X */
bool RWGObject::Contains(const double X[3])
{
  InitkdPanels(false);
  return kdtri_object_contains(kdPanels, X);
}

/***********************************************************************/

/* return whether object contains O; assumes that surfaces
   are non-intersecting (so that it either contains every
   vertex or no vertices of O).

   Note: the Vertices list contains interior reference points in addition
   to actual surface vertices, so we must use a vertex in the Panels list
   rather that just the first vertex in the Vertices list. */
bool RWGObject::Contains(const RWGObject *O)
{
  if (!O || O->NumPanels <= 0) return false;

  /* return false for all PEC objects.  This is because PEC objects
     are actually treated as infinitesimally thin PEC surfaces (which 
     may not even be closed surfaces), so they have no interior. */
  if (O->MP->IsPEC()) return false;

  int vi = O->Panels[0]->VI[0]; // the first vertex of the first panel
#if 0 // debugging: exhaustively check all vertices
  bool cont = Contains(O->Vertices + 3*vi);
  for (int np = 0; np < O->NumPanels; ++np)
    for (int j = 0; j < 3; ++j) {
      int vj = O->Panels[np]->VI[j];
      if (Contains(O->Vertices + 3*vj) != cont) {
	ErrExit("object %s and %s have intersecting surfaces? (%g,%g,%g) is %s %s but (%g,%g,%g) %s", Label, O->Label, O->Vertices[3*vi],O->Vertices[3*vi+1],O->Vertices[3*vi+2], cont ? "is in" : "is not in", Label, O->Vertices[3*vj],O->Vertices[3*vj+1],O->Vertices[3*vj+2], cont ? "is not" : "is");
      }
    }
#endif
  return Contains(O->Vertices + 3*vi);
}

/***********************************************************************/

int RWGGeometry::GetObjectIndex(const double X[3]) {
  // find the innermost object containing X
  for (int i = NumObjects - 1; i >= 0; --i) // innermost to outermost order
    if (Objects[i]->Contains(X))
      return i;
  return -1; // not in any object
}

RWGObject *RWGGeometry::GetObject(const double X[3]) {
  int i = GetObjectIndex(X);
  return i < 0 ? NULL : Objects[i];
}


} // namespace scuff
