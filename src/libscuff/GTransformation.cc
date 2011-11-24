/*
 * GTransformation.cc -- a very simple mechanism for 
 *                    -- handling geometric transformations
 *                    -- (displacements and rotations)
 *                        
 * homer reid         -- 11/2011
 */

#include "GTransformation.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ConstructRotationMatrix(double *ZHat, double Theta, double M[3][3]);

/***************************************************************/
/* constructor 1 ***********************************************/
/***************************************************************/
GTransformation *CreateGTransformation(double *DX)
{ GTransformation *GT=(GTransformation *)malloc(sizeof *GT);
  GT->Type=GTRANSFORMATION_DISPLACEMENT;
  memcpy(GT->DX, DX, 3*sizeof(double));
}

/***************************************************************/
/* constructor 2 ***********************************************/
/***************************************************************/
GTransformation *CreateGTransformation(double *ZHat, double Theta)
{ GTransformation *GT=(GTransformation *)malloc(sizeof *GT);
  GT->Type=GTRANSFORMATION_ROTATION;
  ConstructRotationMatrix(ZHat, Theta, GT->M);
}

/***************************************************************/
/* apply the transformation (in-place) to a list of NX points  */
/* X[3*nx+0, 3*nx+1, 3*nx+2] = cartesian coords of point #nx   */
/***************************************************************/
void ApplyGTransformation(GTransformation *GT, double *X, int NX)
{ 
  int nx; 
  if (GT->Type==GTRANSFORMATION_DISPLACEMENT)
   { for(nx=0; nx<NX; nx++)
      { X[3*nx+0] += GT->DX[0];
        X[3*nx+1] += GT->DX[1];
        X[3*nx+2] += GT->DX[2];
      };
   }
  else
   { int Mu, Nu;
     double XP[3];
     for(nx=0; nx<NX; nx++)
      { memset(XP, 0, 3*sizeof(double));
        for(Mu=0; Mu<3; Mu++) 
         for(Nu=0; Nu<3; Nu++)      
          XP[Mu] += GT->M[Mu][Nu] * X[3*nx+Nu];
        memcpy(X+3*nx, XP, 3*sizeof(double));
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
  int nx; 
  if (GT->Type==GTRANSFORMATION_DISPLACEMENT)
   { for(nx=0; nx<NX; nx++)
      { XP[3*nx+0] = X[3*nx+0] + GT->DX[0];
        XP[3*nx+1] = X[3*nx+1] + GT->DX[1];
        XP[3*nx+2] = X[3*nx+2] + GT->DX[2];
      };
   }
  else
   { int Mu, Nu;
     for(nx=0; nx<NX; nx++)
      { memset(XP+3*nx, 0, 3*sizeof(double));
        for(Mu=0; Mu<3; Mu++) 
         for(Nu=0; Nu<3; Nu++)      
          XP[3*nx + Mu] += GT->M[Mu][Nu] * X[3*nx+Nu];
      };
   };

}

void ApplyGTransformation(GTransformation *GT, double *X, double *XP)
{ ApplyGTransformation(GT, X, XP, 1); }

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

  VecNormalize(ZHat);
 
  /* construct M1 */
  ct=ZHat[2];
  st=sqrt(1.0-ct*ct);
  cp= ( st < 1.0e-8 ) ? 1.0 : ZHat[0] / st;
  sp= ( st < 1.0e-8 ) ? 0.0 : ZHat[1] / st;
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
