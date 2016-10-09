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

/*
 * Vector.cc   -- simple vector manipulation routines (formerly
 *                defined in scuffMisc.cc)
 * homer reid  -- 10/2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "libhrutil.h"

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* vector routines                                              */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/* v <= 0 */
void VecZero(double *v, int N )
{ memset(v,0,N*sizeof(double));
}

/* v2 = v1 */
double *VecCopy(const double *v1, double *v2, int N) 
{
  memcpy(v2, v1, N*sizeof(double));
  return v2;
}

/* v3 = Alpha*v1 + Beta*v2 */
double *VecLinComb(double Alpha, const double *v1, double Beta, const double *v2, double *v3, int N)
{ 
  for(int n=0; n<N; n++)
   v3[n]=Alpha*v1[n] + Beta*v2[n];
  return v3;
}

/* v2 = Alpha * v1 */
double *VecScale(const double *v1, double Alpha, double *v2, int N)
{ return VecLinComb(Alpha, v1, 0.0, v1, v2, N); }

/* v *= Alpha */
double *VecScale(double *v, double Alpha, int N)
{ return VecLinComb(Alpha, v, 0.0, v, v, N); }

/* v3 = v1 + Alpha*v2 */
double *VecScaleAdd(const double *v1, double Alpha, const double *v2, double *v3, int N)
{ return VecLinComb(1.0, v1, Alpha, v2, v3, N); }

/* v3 = v1 + v2 */
double *VecAdd(const double *v1, const double *v2, double *v3, int N)
{ return VecLinComb(1.0, v1, 1.0, v2, v3, N); }

/* v3 = v1 - v2 */
double *VecSub(const double *v1, const double *v2, double *v3, int N)
{ return VecLinComb(1.0, v1, -1.0, v2, v3, N); } 

/* v1 += Alpha*v2 */
double *VecPlusEquals(double *v1, double Alpha, const double *v2, int N)
{ return VecLinComb(1.0, v1, Alpha, v2, v1, N); }

/* return v1 \dot v2 */
double VecDot(const double *v1, const double *v2, int N)
{ double Dot=0.0;
  for(int n=0; n<N; n++) Dot+=v1[n]*v2[n];
  return Dot;
}

double VecNorm2(const double *v, int N)
{ return VecDot(v,v,N); }

double VecNorm(const double *v, int N)
{ return sqrt(VecDot(v,v,N)); }

double VecNormalize(double *v, int N)
{ double d=VecNorm(v, N); 
  VecScale(v, 1.0/d, N);
  return d;
}

double VecDistance2(const double *v1, const double *v2, int N)
{ double d=0.0;
  for(int n=0; n<N; n++) 
   d += (v1[n]-v2[n])*(v1[n]-v2[n]);
  return d;
}

double VecDistance(const double *v1, const double *v2, int N)
{ return sqrt(VecDistance2(v1,v2,N));
}

// note cross product is for 3D vectors only

/* v3 = v1 \times v2 */
double *VecCross(const double *v1, const double *v2, double *v3)
{ v3[0]=v1[1]*v2[2] - v1[2]*v2[1];
  v3[1]=v1[2]*v2[0] - v1[0]*v2[2];
  v3[2]=v1[0]*v2[1] - v1[1]*v2[0];
  return v3;
}

/***************************************************************/
/* complex vector arithmetic ***********************************/
/***************************************************************/
/* v3 = Alpha*v1 + Beta*v2 */
cdouble *VecLinComb(cdouble Alpha, const cdouble *v1, cdouble Beta, const cdouble *v2, cdouble *v3, int N)
{ 
  for(int n=0; n<N; n++)
   v3[n]=Alpha*v1[n] + Beta*v2[n];
  return v3;
}

/* v2 = Alpha * v1 */
cdouble *VecScale(const cdouble *v1, cdouble Alpha, cdouble *v2, int N)
{ return VecLinComb(Alpha, v1, 0.0, v1, v2, N); }

/* v *= Alpha */
cdouble *VecScale(cdouble *v, cdouble Alpha, int N)
{ return VecLinComb(Alpha, v, 0.0, v, v, N); }

/* v3 = v1 + Alpha*v2 */
cdouble *VecScaleAdd(const cdouble *v1, cdouble Alpha, const cdouble *v2, cdouble *v3, int N)
{ return VecLinComb(1.0, v1, Alpha, v2, v3, N); }

/* v3 = v1 + v2 */
cdouble *VecAdd(const cdouble *v1, const cdouble *v2, cdouble *v3, int N)
{ return VecLinComb(1.0, v1, 1.0, v2, v3, N); }

/* v3 = v1 - v2 */
cdouble *VecSub(const cdouble *v1, const cdouble *v2, cdouble *v3, int N)
{ return VecLinComb(1.0, v1, -1.0, v2, v3, N); } 

/* v1 += Alpha*v2 */
cdouble *VecPlusEquals(cdouble *v1, cdouble Alpha, const cdouble *v2, int N)
{ return VecLinComb(1.0, v1, Alpha, v2, v1, N); }

/* v1 += Alpha*v2 */
cdouble *VecPlusEquals(cdouble *v1, cdouble Alpha, const double *v2, int N)
{ for(int n=0; n<N; n++)
   v1[n] += Alpha*v2[n];
  return v1;
}

/***************************************************************/
/* finite-difference derivatives of vector-valued functions    */
/***************************************************************/
cdouble GetDivCurl(VVFunction VVFun, void *UserData, int Order,
                   double X[3], double Delta, cdouble CurlF[3])
{
  bool SecondOrder = (Order==2);

  // dV[Mu][Nu] = d_\mu F_\nu
  cdouble V0[3], dV[3][3];
  
  if (!SecondOrder)
   VVFun(UserData, X, V0);

  for(int Mu=0; Mu<3; Mu++)
   { 
     double XP[3];
     XP[0]=X[0];
     XP[1]=X[1];
     XP[2]=X[2];

     double DX = (X[Mu]==0.0) ? Delta : Delta*fabs(X[Mu]);

     cdouble dVP[3];
     XP[Mu] += DX;
     VVFun(UserData, XP, dVP);

     cdouble dVMBuffer[3], *dVM, Denom;
     if (SecondOrder)
      { XP[Mu] -= 2.0*DX;
        VVFun(UserData, XP, dVMBuffer);
        dVM=dVMBuffer;
        Denom = 2.0*DX;
      }
     else
      { dVM=V0;
        Denom=1.0*DX;
      };
     
     dV[Mu][0] = (dVP[0] - dVM[0]) / Denom;
     dV[Mu][1] = (dVP[1] - dVM[1]) / Denom;
     dV[Mu][2] = (dVP[2] - dVM[2]) / Denom;

   };

  CurlF[0] = dV[1][2] - dV[2][1];
  CurlF[1] = dV[2][0] - dV[0][2];
  CurlF[2] = dV[0][1] - dV[1][0];

  return dV[0][0] + dV[1][1] + dV[2][2];
}

/***************************************************************/
/* matrix arithmetic *******************************************/
/***************************************************************/
bool Matrix2x2_Inverse(double *a[2], double ainv[2][2])
{
  double det = a[0][0]*a[1][1] - a[0][1]*a[1][0];
  double detinv = 1 / det;
  // store in local vars to work even if a == ainv
  double ainv00 = detinv * a[1][1];
  double ainv11 = detinv * a[0][0];
  double ainv01 = -detinv * a[0][1];
  double ainv10 = -detinv * a[1][0];
  ainv[0][0] = ainv00;
  ainv[1][1] = ainv11;
  ainv[0][1] = ainv01;
  ainv[1][0] = ainv10;
  return det != 0.0;
}

/***************************************************************/
/* single-precision comparisons of double-precision numbers    */
/***************************************************************/
bool EqualFloat(const double a, const double b) 
{ return ( float(a) == float(b) ); }

bool EqualFloat(const cdouble a, const cdouble b) 
{ return      ( float(real(a)) == float(real(b)) )
          &&  ( float(imag(a)) == float(imag(b)) );
}

bool VecEqualFloat(const double *a, const double *b) 
{
   return (     float(a[0]) == float(b[0])
            &&  float(a[1]) == float(b[1])
            &&  float(a[2]) == float(b[2])
          );
}

bool VecClose(const double *a, const double *b, double abstol)
{
  return fabs(a[0]-b[0]) + fabs(a[1]-b[1]) + fabs(a[2]-b[2]) <= 3*abstol;
}   

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Compare(cdouble *V1, cdouble *V2, int N,
             const char *str1, const char *str2)
{ 
  printf(" n | %-25s | %-25s | RD      | Ratio\n",str1,str2);
  for(int n=0; n<N; n++)
   printf("%2i | (%+.4e,%+.4e) | (%+.4e,%+.4e) | %.1e | %.3e\n",n,
    real(V1[n]),imag(V1[n]), real(V2[n]),imag(V2[n]),
    RD(V1[n],V2[n]), abs(V1[n]/V2[n]));
  printf("\n");
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Compare(double *V1, double *V2, int N,
             const char *str1, const char *str2)
{ 
  printf(" n | %-12s | %-12s | RD      | Ratio\n",str1,str2);
  for(int n=0; n<N; n++)
   printf("%2i | %+12.4e | %+12.4e | %.1e | %.3e\n",n,
    V1[n],V2[n],RD(V1[n],V2[n]),fabs(V1[n]/V2[n]));
  printf("\n");
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Compare(double *V1, double *V2, double *V3, int N,
             const char *str1, const char *str2, const char *str3)
{ 
  printf(" n | %-12s | %-12s | RD      | %-12s | RD\n",str1,str2,str3);
  for(int n=0; n<N; n++)
   printf("%2i | %+12.4e | %+12.4e | %.1e | %+12.4e | %.1e\n",n,
    V1[n], V2[n], RD(V1[n],V2[n]), V3[n], RD(V1[n],V3[n]));
  printf("\n");
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
void Compare(cdouble *V1, cdouble *V2, cdouble *V3, int N,
             const char *str1, const char *str2, const char *str3)
{ 
  printf(" n | %-25s | %-25s | RD      | %-25s | RD\n",str1,str2,str3);
  for(int n=0; n<N; n++)
   printf("%2i | (%+.4e,%+.4e) | (%+.4e,%+.4e) | %.1e | (%+.4e, %+.4e) | %.1e\n",n,
    real(V1[n]),imag(V1[n]),
    real(V2[n]),imag(V2[n]),
    RD(V1[n],V2[n]), 
    real(V3[n]),imag(V3[n]),
    RD(V1[n],V3[n]));
  printf("\n");
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DofprintVec(FILE *f, void *v, bool Complex, int Length,
                 const char *fmt, bool CR)
{
  char format[100];
  snprintf(format,100,"%s ",fmt);

  if (Complex)
   { cdouble *cv=(cdouble *)v;
     for(int n=0; n<Length; n++)
      fprintf(f,format,real(cv[n]),imag(cv[n]));
   }
  else
   { double *rv=(double *)v;
     for(int n=0; n<Length; n++)
      fprintf(f,format,rv[n]);
   };

  if (CR) 
   fprintf(f,"\n");
}

void fprintVec(FILE *f, double *v, int Length, const char *format)
{  DofprintVec(f, (void *)v, false, Length, format, false); }

void fprintVecCR(FILE *f, double *v, int Length, const char *format)
{  DofprintVec(f, (void *)v, false, Length, format, true); }

void fprintVec(FILE *f, cdouble *v, int Length, const char *format)
{  DofprintVec(f, (void *)v, true, Length, format, false); }

void fprintVecCR(FILE *f, cdouble *v, int Length, const char *format)
{  DofprintVec(f, (void *)v, true, Length, format, true); }
