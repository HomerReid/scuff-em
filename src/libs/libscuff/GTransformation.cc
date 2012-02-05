/*
 * GTransformation.cc -- a very simple mechanism for handling geometric 
 *                    -- transformations (displacements and rotations)
 *                        
 * homer reid         -- 11/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "GTransformation.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void ConstructRotationMatrix(double *ZHat, double Theta, double M[3][3]);

/***************************************************************/
/* a transformation that displaces through vector DX           */
/***************************************************************/
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, double *DX)
{
  GTransformation *NGT;

  if (GT==0)
   { 
     NGT=(GTransformation *)malloc(sizeof(*NGT));

     memcpy(NGT->DX, DX, 3*sizeof(double));

     memset(NGT->M[0], 0, 3*sizeof(double));
     memset(NGT->M[1], 0, 3*sizeof(double));
     memset(NGT->M[2], 0, 3*sizeof(double));
     NGT->M[0][0]=NGT->M[1][1]=NGT->M[2][2]=1.0;

     NGT->Type=GTRANSFORMATION_DISPLACEMENT;

     return NGT;
   }
  else
   { GT->DX[0] += DX[0];
     GT->DX[1] += DX[1];
     GT->DX[2] += DX[2];
     GT->Type |= GTRANSFORMATION_DISPLACEMENT;
     return GT;
   };

}

GTransformation *CreateGTransformation(double *DX)
 { CreateOrAugmentGTransformation(0, DX); }

/***************************************************************/
/* a transformation that rotates through Theta degrees         */
/* (DEGREES, NOT RADIANS) about an axis that passes through    */
/* the origin and the point with coordinates ZHat[0..2]        */
/***************************************************************/
GTransformation *CreateOrAugmentGTransformation(GTransformation *GT, 
                                                double *ZHat, double Theta)
{
  GTransformation *NGT;

  if (GT==0)
   { 
     NGT=(GTransformation *)malloc(sizeof(*NGT));

     memset(NGT->DX, 0, 3*sizeof(double));

     ConstructRotationMatrix(ZHat, Theta, NGT->M);

     NGT->Type=GTRANSFORMATION_ROTATION;

     return NGT;
   }
  else
   { double MP[3][3], NewDX[3], NewM[3][3];
     int i, j, k;

     ConstructRotationMatrix(ZHat, Theta, MP);

     for(i=0; i<3; i++)
      for(NewDX[i]=0.0, j=0; j<3; j++)
       NewDX[i] += MP[i][j]*(GT->DX[j]);
     memcpy(GT->DX, NewDX, 3*sizeof(double));

     for(i=0; i<3; i++)
      for(j=0; j<3; j++)
       for(NewM[i][j]=0.0, k=0; k<3; k++)
        NewM[i][j] += MP[i][k]*GT->M[k][j];

     for(i=0; i<3; i++)
      memcpy(GT->M[i],NewM[i],3*sizeof(double));

     GT->Type |= GTRANSFORMATION_ROTATION;
     return GT;
   };
}

GTransformation *CreateGTransformation(double *ZHat, double Theta)
 { CreateOrAugmentGTransformation(0, ZHat, Theta); }

/***************************************************************/
/* GT -> DeltaGT*GT ********************************************/
/***************************************************************/
void AugmentGTransformation(GTransformation *DeltaGT, GTransformation *GT)
{ 
  int i, j, k;
  double NewM[3][3], NewDX[3];
 
  if ( DeltaGT->Type & GTRANSFORMATION_ROTATION )
   { 
     for(i=0; i<3; i++)
      for(j=0; j<3; j++)
       for(NewM[i][j]=0.0, k=0; k<3; k++)
        NewM[i][j] += DeltaGT->M[i][k] * GT->M[k][j];

     for(i=0; i<3; i++)
      memcpy(GT->M[i], NewM[i], 3*sizeof(double));

     for(i=0; i<3; i++)
      for(NewDX[i]=0.0, j=0; j<3; j++)
       NewDX[i] += DeltaGT->M[i][j]*GT->DX[j];

     memcpy(GT->DX, NewDX, 3*sizeof(double));

   };

  GT->DX[0] += DeltaGT->DX[0];
  GT->DX[1] += DeltaGT->DX[1];
  GT->DX[2] += DeltaGT->DX[2];
   
} 

/***************************************************************/
/* apply the transformation (in-place) to a list of NX points  */
/* X[3*nx+0, 3*nx+1, 3*nx+2] = cartesian coords of point #nx   */
/***************************************************************/
void ApplyGTransformation(GTransformation *GT, double *X, int NX)
{ 
  if (GT==0) 
   return;

  int nx; 
  int i, j;
  double XP[3];

  if ( GT->Type & GTRANSFORMATION_ROTATION )
   { 
     for(nx=0; nx<NX; nx++)
      { memcpy(XP, GT->DX, 3*sizeof(double));
        for(i=0; i<3; i++) 
         for(j=0; j<3; j++)      
          XP[i] += GT->M[i][j] * X[3*nx+j];
        memcpy(X+3*nx, XP, 3*sizeof(double));
      };
   }
  else
   { for(nx=0; nx<NX; nx++)
      { X[3*nx+0] += GT->DX[0];
        X[3*nx+1] += GT->DX[1];
        X[3*nx+2] += GT->DX[2];
      };
   };
}

void ApplyGTransformation(GTransformation *GT, double *X)
{ ApplyGTransformation(GT, X, 1); }

/***************************************************************/
/* like the above, but operate out-of-place.                   */
/***************************************************************/
void ApplyGTransformation(GTransformation *GT, double *X, double *XP, int NX)
{ 
  if (GT==0) 
   { memcpy(XP, X, 3*NX*sizeof(double));
     return;
   };

  int nx; 
  int i, j;
  if ( GT->Type & GTRANSFORMATION_ROTATION )
   { 
     for(nx=0; nx<NX; nx++)
      { memcpy(XP+3*nx, GT->DX, 3*sizeof(double));
        for(i=0; i<3; i++) 
         for(j=0; j<3; j++)      
          XP[3*nx+i] += GT->M[i][j] * X[3*nx+j];
      };
   }
  else
   { for(nx=0; nx<NX; nx++)
      { XP[3*nx+0] = X[3*nx+0] + GT->DX[0];
        XP[3*nx+1] = X[3*nx+1] + GT->DX[1];
        XP[3*nx+2] = X[3*nx+2] + GT->DX[2];
      };
   };
}

void ApplyGTransformation(GTransformation *GT, double *X, double *XP)
 { ApplyGTransformation(GT, X, XP, 1); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void UnApplyGTransformation(GTransformation *GT, double *X, int NX)
{
  int i, j, nx;
  double XP[3];

  for(nx=0; nx<NX; nx++)
   { X[3*nx+0] -= GT->DX[0];
     X[3*nx+1] -= GT->DX[1];
     X[3*nx+2] -= GT->DX[2];
   };

  if ( GT->Type & GTRANSFORMATION_ROTATION )
   { for(nx=0; nx<NX; nx++)
      { 
        for(i=0; i<3; i++)
         for(XP[i]=0.0, j=0; j<3; j++)
          XP[i] += GT->M[j][i]*X[3*nx+j];

        memcpy(X+3*nx, XP, 3*sizeof(double));

      };
   };
}

void ResetGTransformation(GTransformation *GT)
{ 
  memset(GT->DX, 0, 3*sizeof(double));
  memset(GT->M[0], 0, 3*sizeof(double));
  memset(GT->M[1], 0, 3*sizeof(double));
  memset(GT->M[2], 0, 3*sizeof(double));

  GT->M[0][0]=GT->M[1][1]=GT->M[2][2]=1.0;

  GT->Type = GTRANSFORMATION_DISPLACEMENT;

}
   

/***************************************************************/
/* Construct the 3x3 matrix that represents a rotation of      */
/* Theta degrees (note we interpret Theta in DEGREES, NOT      */
/* RADIANS!) about the axis specified by cartesian coordinates */
/* ZHat.                                                       */
/* Algorithm:                                                  */
/*  1. Construct matrix M1 that rotates Z axis into alignment  */
/*     with ZHat.                                              */
/*  2. Construct matrix M2 that rotates through Theta about    */ 
/*     Z axis.                                                 */
/*  3. Construct matrix M=M1^{-1}*M2*M1=M1^T*M2*M1.            */  
/***************************************************************/
static void ConstructRotationMatrix(double *ZHat, double Theta, double M[3][3])
{ 
  int Mu, Nu, Rho;
  double ct, st, cp, sp, CT, ST;
  double M2M1[3][3], M2[3][3], M1[3][3];

  // first normalize ZHat
  double nZHat=sqrt(ZHat[0]*ZHat[0] + ZHat[1]*ZHat[1] + ZHat[2]*ZHat[2]);
  double NZHat[3];
  NZHat[0]=ZHat[0] / nZHat;
  NZHat[1]=ZHat[1] / nZHat;
  NZHat[2]=ZHat[2] / nZHat;
 
  /* construct M1 */
  ct=NZHat[2];
  st=sqrt(1.0-ct*ct);
  cp= ( st < 1.0e-8 ) ? 1.0 : NZHat[0] / st;
  sp= ( st < 1.0e-8 ) ? 0.0 : NZHat[1] / st;
  M1[0][0]=ct*cp;  M1[0][1]=ct*sp;   M1[0][2]=-st;
  M1[1][0]=-sp;    M1[1][1]=cp;      M1[1][2]=0.0;
  M1[2][0]=st*cp;  M1[2][1]=st*sp;   M1[2][2]=ct;

  /* construct M2 */
  CT=cos(Theta*M_PI/180.0);
  ST=sin(Theta*M_PI/180.0);
  M2[0][0]=CT;      M2[0][1]=-ST;     M2[0][2]=0.0;
  M2[1][0]=ST;      M2[1][1]=CT;      M2[1][2]=0.0;
  M2[2][0]=0.0;     M2[2][1]=0.0;     M2[2][2]=1.0;

  /* M2M1 <- M2*M1 */
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    for(M2M1[Mu][Nu]=0.0, Rho=0; Rho<3; Rho++)
     M2M1[Mu][Nu] += M2[Mu][Rho] * M1[Rho][Nu];

  /* M <- M1^T * M2M1 */
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    for(M[Mu][Nu]=0.0, Rho=0; Rho<3; Rho++)
     M[Mu][Nu] += M1[Rho][Mu]*M2M1[Rho][Nu];
}
