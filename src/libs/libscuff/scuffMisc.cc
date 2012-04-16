/*
 * scuffMisc.cc  -- miscellaneous stuff for libscuff
 *
 * homer reid  -- 10/2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "libscuff.h"

namespace scuff {

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* vector routines                                              */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/* v <= 0 */
void VecZero(double v[3])
{ memset(v,0,3*sizeof(double));
}

/* v2 = v1 */
double *VecCopy(const double v1[3], double v2[3]) {
  v2[0] = v1[0]; v2[1] = v1[1]; v2[2] = v1[2];
  return v2;
}

/* v2 = alpha * v1 */
double *VecScale(const double v1[3], double alpha, double v2[3])
{
  v2[0] = alpha * v1[0];
  v2[1] = alpha * v1[1];
  v2[2] = alpha * v1[2];
  return v2;
}

/* v *= alpha */
double *VecScale(double v[3], double alpha)
{ v[0]*=alpha;
  v[1]*=alpha;
  v[2]*=alpha;
  return v;
}

/* v3 = v1 + alpha*v2 */
double *VecScaleAdd(const double v1[3], double alpha, const double v2[3], double v3[3])
{ v3[0]=v1[0] + alpha*v2[0];
  v3[1]=v1[1] + alpha*v2[1];
  v3[2]=v1[2] + alpha*v2[2];
  return v3;
}

/* v3 = alpha*v1 + beta*v2 */
double *VecLinComb(double alpha, const double v1[3], double beta, const double v2[3], double v3[3])
{ v3[0]=alpha*v1[0] + beta*v2[0];
  v3[1]=alpha*v1[1] + beta*v2[1];
  v3[2]=alpha*v1[2] + beta*v2[2];
  return v3;
}

/* v3 = v1 + v2 */
double *VecAdd(const double v1[3], const double v2[3], double v3[3])
{ v3[0]=v1[0] + v2[0];
  v3[1]=v1[1] + v2[1];
  v3[2]=v1[2] + v2[2];
  return v3;
}

/* v3 = v1 - v2 */
double *VecSub(const double v1[3], const double v2[3], double v3[3])
{ v3[0]=v1[0] - v2[0];
  v3[1]=v1[1] - v2[1];
  v3[2]=v1[2] - v2[2];
  return v3;
}

/* v1 += alpha*v2 */
double *VecPlusEquals(double v1[3], double alpha, const double v2[3])
{ v1[0]+=alpha*v2[0];
  v1[1]+=alpha*v2[1];
  v1[2]+=alpha*v2[2];
  return v1;
}

/* v3 = v1 \times v2 */
double *VecCross(const double v1[3], const double v2[3], double v3[3])
{ v3[0]=v1[1]*v2[2] - v1[2]*v2[1];
  v3[1]=v1[2]*v2[0] - v1[0]*v2[2];
  v3[2]=v1[0]*v2[1] - v1[1]*v2[0];
  return v3;
}

/* return v1 \dot v2 */
double VecDot(const double v1[3], const double v2[3])
{ return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }

double VecNorm2(const double v[3])
{ return VecDot(v,v); }

double VecNorm(const double v[3])
{ return sqrt(VecDot(v,v)); }

double VecDistance(const double v1[3], const double v2[3])
{ double v3[3];
  VecSub(v1,v2,v3);
  return VecNorm(v3);
}

double VecDistance2(const double v1[3], const double v2[3])
{ double v3[3];
  VecSub(v1,v2,v3);
  return VecDot(v3,v3);
}

double VecNormalize(double v[3])
{ double d=VecNorm(v);  
  v[0]/=d;
  v[1]/=d;
  v[2]/=d;
  return d;
}

void SixVecPlus(const cdouble V1[6], const cdouble Alpha,
                const cdouble V2[6], cdouble V3[6])
{ 
  for(int n=0; n<6; n++)
   V3[n] = V1[n] + Alpha*V2[n];
}

void SixVecPlusEquals(const cdouble V1[6], const cdouble Alpha,
                      const cdouble V2[6])
{ 
  for(int n=0; n<6; n++)
   V1[n] += Alpha*V2[n];
}

void SixVecPlusEquals(const cdouble V1[6], const cdouble V2[6])
 { SixVecPlusEquals(V1, 1.0, V2); }

} // namespace scuff
