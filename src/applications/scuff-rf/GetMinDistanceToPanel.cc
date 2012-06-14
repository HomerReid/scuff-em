#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libRWG.h>

int GMDTPVerbose=0;

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetMinDistanceToPanel(double *V1, double *V2, double *V3,
                             double *X, double *XMin, int *OnBoundary)
{ 
  double XmV[3], A[3], B[3], uMin, vMin;

  VecSub(X, V1, XmV);
  VecSub(V2, V1, A);
  VecSub(V3, V2, B);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double A2=VecNorm2(A);
  double B2=VecNorm2(B);
  double AdB=VecDot(A, B);
  double XmVdA=VecDot(XmV, A);
  double XmVdB=VecDot(XmV, B);
  double Denominator = A2*B2 - AdB*AdB;

  uMin = (B2*XmVdA - AdB*XmVdB) / Denominator;
  vMin = (A2*XmVdB - AdB*XmVdA) / Denominator;

if (GMDTPVerbose)
 { printf("%e %e %e = \n", VecDistance(X,V1), VecDistance(X,V2), VecDistance(X,V3));
   printf("uMin, vMin=%e,%e\n",uMin,vMin);
 };

  /***************************************************************/
  /* if the point found that way lies inside the triangle, then  */
  /* we accept it as the optimal point                           */
  /***************************************************************/
  if (uMin>=0.0 && uMin<=1.0 && vMin>=0.0 && vMin<=uMin)
   { 
     *OnBoundary=0;
     XMin[0] = V1[0] + uMin*A[0] + vMin*B[0];
     XMin[1] = V1[1] + uMin*A[1] + vMin*B[1];
     XMin[2] = V1[2] + uMin*A[2] + vMin*B[2];
     return VecDistance(X, XMin);
   };

  /***************************************************************/
  /* otherwise the minimum is attained on one of the three edges */
  /***************************************************************/
  *OnBoundary=1;
  double Distance, MinDistance=1e9;
  double XP[3];

  /* look on the line V1--V2 */
  VecSub(X, V1, XmV );
  VecSub(V2, V1, A);
  uMin = VecDot(XmV , A) / VecNorm2(A);
  if (uMin<0.0) uMin=0.0;
  if (uMin>1.0) uMin=1.0;
  XP[0] = V1[0] + uMin*A[0];
  XP[1] = V1[1] + uMin*A[1];
  XP[2] = V1[2] + uMin*A[2];
  Distance=VecDistance(X, XP);
  if (Distance<MinDistance)
   { MinDistance=Distance;
     memcpy(XMin, XP, 3*sizeof(double));
   };

  /* look on the line V2--V3 */
  VecSub(X, V2, XmV );
  VecSub(V3, V2, A);
  uMin = VecDot(XmV , A) / VecNorm2(A);
  if (uMin<0.0) uMin=0.0;
  if (uMin>1.0) uMin=1.0;
  XP[0] = V2[0] + uMin*A[0];
  XP[1] = V2[1] + uMin*A[1];
  XP[2] = V2[2] + uMin*A[2];
  Distance=VecDistance(X, XP);
  if (Distance<MinDistance)
   { MinDistance=Distance;
     memcpy(XMin, XP, 3*sizeof(double));
   };

  /* look on the line V3--V1 */
  VecSub(X, V3, XmV );
  VecSub(V1, V3, A);
  uMin = VecDot(XmV , A) / VecNorm2(A);
  if (uMin<0.0) uMin=0.0;
  if (uMin>1.0) uMin=1.0;
  XP[0] = V3[0] + uMin*A[0];
  XP[1] = V3[1] + uMin*A[1];
  XP[2] = V3[2] + uMin*A[2];
  Distance=VecDistance(X, XP);
  if (Distance<MinDistance)
   { MinDistance=Distance;
     memcpy(XMin, XP, 3*sizeof(double));
   };

  return MinDistance;

}
