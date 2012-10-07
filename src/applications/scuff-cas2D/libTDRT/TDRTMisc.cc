/*
 * TDRTMisc.cc -- miscellaneous routines for libcas2D 
 */

#include <stdlib.h>
#include <math.h>
#include "libTDRT.h"

double VecDot(double *v1, double *v2)
{ return v1[0]*v2[0] + v1[1]*v2[1]; }

double VecNorm(double *v)
{ return sqrt(v[0]*v[0] + v[1]*v[1]); }

double VecNormalize(double *v)
{ double vv=VecNorm(v);
  v[0]/=vv;
  v[1]/=vv;
  return vv;
}


/* v3 = v1 - v2 */
void VecSub(double *v1, double *v2, double *v3)
{ v3[0]=v1[0] - v2[0];
  v3[1]=v1[1] - v2[1];
}

void VecScale(double *v, double alpha)
{ v[0]*=alpha;
  v[1]*=alpha;
}

void VecScaleAdd(double *v1, double alpha, double *v2, double *v3)
{ v3[0]=v1[0] + alpha*v2[0];
  v3[1]=v1[1] + alpha*v2[1];
}

void VecPlusEquals(double *v1, double alpha, double *v2)
{ v1[0]+= alpha*v2[0];
  v1[1]+= alpha*v2[1];
}

double VecDistance(double *v1, double *v2)
{ double v3[2];

  v3[0]=v1[0]-v2[0];
  v3[1]=v1[1]-v2[1];
  return sqrt(v3[0]*v3[0] + v3[1]*v3[1]);
}


double VecD2(double *v1, double *v2)
{ double v3[2];

  v3[0]=v1[0]-v2[0];
  v3[1]=v1[1]-v2[1];

  return (v3[0]*v3[0] + v3[1]*v3[1]);
}
