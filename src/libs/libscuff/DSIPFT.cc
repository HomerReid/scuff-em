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

/***************************************************************/
/* DSIPFT.cc  -- SCUFF-EM code for computing power, force, and */
/*            -- torque (PFT) using the 'displaced surface     */
/*            -- integral' (DSI) method                        */
/*                                                             */
/* Homer Reid -- 1/2014                                        */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libTriInt.h>
#include <config.h>
 
namespace scuff{

void GetReducedFarFields(RWGGeometry *G, const int ns, const int ne,
                         const double X0[3], const cdouble k,
                         cdouble e[3], cdouble h[3]);

#define II cdouble(0.0,1.0) 
#define TENTHIRDS 3.33333333333333333333333

#define SIPOWER   0
#define SIXFORCE  1
#define SIYFORCE  2
#define SIZFORCE  3
#define SIXTORQUE 4
#define SIYTORQUE 5
#define SIZTORQUE 6
#define NUMPFT    7

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetReducedFields(RWGGeometry *G, int ns, int ne,
                      double X0[3],  cdouble k,
                      cdouble e[3], cdouble h[3])
{
  cdouble a[3], curla[3], Gradp[3];

  G->Surfaces[ns]->GetReducedPotentials(ne, X0, k, 0, a, curla, Gradp);

  cdouble k2=k*k;

  e[0] = a[0] + Gradp[0]/k2;
  e[1] = a[1] + Gradp[1]/k2;
  e[2] = a[2] + Gradp[2]/k2;
  h[0] = curla[0];
  h[1] = curla[1];
  h[2] = curla[2];

}


/***************************************************************/
/* Compute the G and C dyadic Green's functions, retaining only*/
/* far-field contributions.                                    */
/***************************************************************/
void CalcGCFarField(double R[3], cdouble k,
                    cdouble GMuNu[3][3], cdouble CMuNu[3][3])
{ 
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];

  if (r2==0.0)
   return;

  double r=sqrt(r2);

  cdouble ExpFac=exp(II*k*r)/(4.0*M_PI*r);

  GMuNu[0][0] =               ExpFac * ( 1.0 - R[0]*R[0]/r2 );
  GMuNu[0][1] = GMuNu[1][0] = ExpFac * (     - R[0]*R[1]/r2 );
  GMuNu[0][2] = GMuNu[2][0] = ExpFac * (     - R[0]*R[2]/r2 );
  GMuNu[1][1] =               ExpFac * ( 1.0 - R[1]*R[1]/r2 );
  GMuNu[1][2] = GMuNu[2][1] = ExpFac * (     - R[1]*R[2]/r2 );
  GMuNu[2][2] =               ExpFac * ( 1.0 - R[2]*R[2]/r2 );

  CMuNu[0][0]=CMuNu[1][1]=CMuNu[2][2]=0.0;

  CMuNu[0][1] = ExpFac * R[2]/r;
  CMuNu[1][2] = ExpFac * R[0]/r;
  CMuNu[2][0] = ExpFac * R[1]/r;

  CMuNu[1][0] = -CMuNu[0][1];
  CMuNu[2][1] = -CMuNu[1][2];
  CMuNu[0][2] = -CMuNu[2][0];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFarFieldGCIntegrals(RWGSurface *S, int ne,
                            double *X, cdouble k,
                            cdouble GInt[3], cdouble CInt[3])
{
  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  RWGEdge *E    = S->Edges[ne];
  double *QP    = S->Vertices + 3*E->iQP;
  double *V1    = S->Vertices + 3*E->iV1;
  double *V2    = S->Vertices + 3*E->iV2;
  double *QM    = S->Vertices + 3*E->iQM;
  double Length = E->Length;

  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/
  double AP[3], BP[3], AM[3], BM[3];
  for(int Mu=0; Mu<3; Mu++)
   { AP[Mu] = V1[Mu] - QP[Mu];
     AM[Mu] = V1[Mu] - QM[Mu];
     BP[Mu] = V2[Mu] - QP[Mu];
     BM[Mu] = V2[Mu] - QM[Mu];
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumPts;
  double *TCR = GetTCR(9, &NumPts);
  GInt[0]=GInt[1]=GInt[2]=CInt[0]=CInt[1]=CInt[2]=0.0;
  for(int np=0, ncp=0; np<NumPts; np++)
   { 
     double u=TCR[ncp++];
     double v=TCR[ncp++]; 
     double w=TCR[ncp++] * Length;

     /***************************************************************/
     /* set X and F=X-Q *********************************************/
     /***************************************************************/
     double XmQP[3], XP[3], XmXP[3], XmQM[3], XM[3], XmXM[3];
     for(int Mu=0; Mu<3; Mu++)
      { XmQP[Mu] = u*AP[Mu] + v*BP[Mu];
          XP[Mu] = XmQP[Mu] + QP[Mu];
        XmXP[Mu] = X[Mu] - XP[Mu];
        XmQM[Mu] = u*AM[Mu] + v*BM[Mu];
          XM[Mu] = XmQM[Mu] + QM[Mu];
        XmXM[Mu] = X[Mu] - XM[Mu];
      };

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     cdouble GP[3][3], CP[3][3], GM[3][3], CM[3][3];
     CalcGCFarField(XmXP, k, GP, CP);
     CalcGCFarField(XmXM, k, GM, CM);
     GInt[0]=GInt[1]=GInt[2]=CInt[0]=CInt[1]=CInt[2]=0.0;
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { GInt[Mu] += w*(   GP[Mu][Nu]*XmQP[Mu]
                         - GM[Mu][Nu]*XmQM[Mu]
                       );

         CInt[Mu] += w*(   CP[Mu][Nu]*XmQP[Mu]
                         - CM[Mu][Nu]*XmQM[Mu]
                       );
       };

   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetFarFields(RWGGeometry *G, IncField *IF, HVector *KN,
                      cdouble Omega, HMatrix *XMatrix)
{ 
  (void)G;
  (void)IF;
  (void)KN;
  (void)Omega;
  (void)XMatrix;
#if 0
  Log("Computing far fields at %i points...",XMatrix->NR);
//#ifdef USE_OPENMP
  int NumThreads=GetNumThreads();
//#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
//#endif
  for(int nenx=0; nenx<NE*NX; nenx++)
   { 
     int ne=nenx / NX;
     int nx=nenx % NX;

     GetFarFields(S, ne, X, k, GInt, CInt);

     HMatrix *GetFSVMatrix(RWGSurface *S, HMatrix *CRMatrix,
                      cdouble Omega, cdouble Eps, cdouble Mu,
                      bool FarField)

#endif
  return 0;
}

/***************************************************************/
/* CR matrix stands for 'Cubature rule matrix.'                */
/*                                                             */
/* A CR matrix is an NCx7 matrix, where the nth row has entries*/
/*  w_n x_n y_n z_n nx_n ny_n nz_n                             */
/* where w is the cubature weight, (x,y,z) is the cubature pt, */
/* and nx_n, ny_n, nz_n is the outward-pointing surface normal.*/
/*                                                             */
/* If BS is non-null, then it describes a closed triangulated  */
/* surface, and the cubature rule evaluates an integral over   */
/* this surface by doing a one-point cubature over each panel. */
/*                                                             */
/* Otherwise, the cubature rule describes a cubature rule with */
/* NumPoints cubature points over a sphere of radius R.        */
/*                                                             */
/* If Lebedev is true, this is a Lebedev cubature rule. In     */
/* this case, NumPoints must be one of the numbers of cubature */
/* points supported by the GetLebedevRule() routine in         */
/* libTriInt.                                                  */
/*                                                             */
/* Otherwise (Lebedev==false) the cubature rule is a product   */
/* rule with a Clenshaw-Curtis grid in the Theta direction and */
/* an evenly-spaced grid in the Phi direction, and NumPoints   */
/* should be an odd integer between 9 and 99 inclusive.        */
/*                                                             */
/* If OTGT and/or GT are non-null, then each cubature point and*/
/* normal vector is transformed by OTGT (first) and/or GT.     */
/***************************************************************/
HMatrix *GetCRMatrix(char *BSMesh, double R, int NumPoints, bool Lebedev,
                     GTransformation *OTGT, GTransformation *GT)
{
  HMatrix *CRMatrix;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (BSMesh)
   { 
     RWGSurface *BS=new RWGSurface(BSMesh);
     CRMatrix = new HMatrix(BS->NumPanels, 7);
     for(int np=0; np<BS->NumPanels; np++)
      { CRMatrix->SetEntry(np, 0, BS->Panels[np]->Area);
        CRMatrix->SetEntry(np, 1, BS->Panels[np]->Centroid[0]);
        CRMatrix->SetEntry(np, 2, BS->Panels[np]->Centroid[1]);
        CRMatrix->SetEntry(np, 3, BS->Panels[np]->Centroid[2]);
        CRMatrix->SetEntry(np, 4, BS->Panels[np]->ZHat[0]);
        CRMatrix->SetEntry(np, 5, BS->Panels[np]->ZHat[1]);
        CRMatrix->SetEntry(np, 6, BS->Panels[np]->ZHat[2]);
      };
     delete BS;
   }
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  else if (Lebedev)
   { double *LRule = GetLebedevRule(NumPoints);
     if (LRule==0) ErrExit("no Lebedev rule with %i points",NumPoints);
     CRMatrix = new HMatrix(NumPoints, 7);
     for(int np=0; np<NumPoints; np++)
      { CRMatrix->SetEntry(np,0, R*R*LRule[4*np + 3]);
        CRMatrix->SetEntry(np,1, R*LRule[4*np + 0]);
        CRMatrix->SetEntry(np,2, R*LRule[4*np + 1]);
        CRMatrix->SetEntry(np,3, R*LRule[4*np + 2]);
        CRMatrix->SetEntry(np,4, LRule[4*np + 0]);
        CRMatrix->SetEntry(np,5, LRule[4*np + 1]);
        CRMatrix->SetEntry(np,6, LRule[4*np + 2]);
      };
   }
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  else
   {
     int NumThetaPoints=NumPoints;
     int NumPhiPoints=2*(NumPoints+1);
     int TotalPoints = NumPhiPoints*(NumThetaPoints-2) + 2;
     CRMatrix = new HMatrix(TotalPoints, 7);
     double dPhi= 2.0*M_PI / NumPhiPoints;

     double *CCRule = GetCCRule(NumThetaPoints);
     if (CCRule==0) ErrExit("invalid value of NumPoints in GetCRMatrix");
     for(int npTheta=0; npTheta<NumThetaPoints; npTheta++)
      { 
        double CosTheta = CCRule[2*npTheta+0];
        double w = CCRule[2*npTheta+1];
        double SinTheta = sqrt(1.0-CosTheta*CosTheta);

        // north, south pole
        if ( npTheta==0 || npTheta==(NumThetaPoints-1) )
         { int Index  = (npTheta==0) ? 0 : TotalPoints-1;
           CRMatrix->SetEntry(Index, 0, 2.0*M_PI*R*R*w);
           CRMatrix->SetEntry(Index, 1, 0.0);
           CRMatrix->SetEntry(Index, 2, 0.0);
           CRMatrix->SetEntry(Index, 3, R*CosTheta);
           CRMatrix->SetEntry(Index, 4, 0.0);
           CRMatrix->SetEntry(Index, 5, 0.0);
           CRMatrix->SetEntry(Index, 6, CosTheta);
         }
        else
         { for(int npPhi=0; npPhi<NumPhiPoints; npPhi++)
            { int Index  = 1 + NumPhiPoints*(npTheta-1) + npPhi;
              double Phi = npPhi * dPhi;
              CRMatrix->SetEntry(Index, 0, R*R*w*dPhi);
              CRMatrix->SetEntry(Index, 1, R*SinTheta*cos(Phi));
              CRMatrix->SetEntry(Index, 2, R*SinTheta*sin(Phi));
              CRMatrix->SetEntry(Index, 3, R*CosTheta);
              CRMatrix->SetEntry(Index, 4, SinTheta*cos(Phi));
              CRMatrix->SetEntry(Index, 5, SinTheta*sin(Phi));
              CRMatrix->SetEntry(Index, 6, CosTheta);
            };
         };
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (OTGT || GT)
   { 
     for(int nr=0; nr<CRMatrix->NR; nr++)
      { 
        double X[3], nHat[3], XP[3];

        CRMatrix->GetEntriesD(nr,"1:3",X);
        CRMatrix->GetEntriesD(nr,"4:6",nHat);
        VecScaleAdd(X, 1.0, nHat, XP);

        if (OTGT) 
         { OTGT->Apply(X); OTGT->Apply(XP); }
        if (GT) 
         { GT->Apply(X); GT->Apply(XP); }

        VecSub(XP, X, nHat);

        CRMatrix->SetEntriesD(nr,"1:3",X);
        CRMatrix->SetEntriesD(nr,"4:6",nHat);
      };
   };

  return CRMatrix;

}

/***************************************************************/
/* Compute a matrix of field six-vectors.                      */
/*                                                             */
/* Each column of an FSVMatrix is a 6-component vector giving  */
/* the E and H fields, at a single point in space, due to a    */
/* single RWG basis function.                                  */
/*                                                             */
/* The rows of FSV matrix are indexed by pairs                 */
/*                                                             */
/*  (evaluation point index, basis function index)             */
/*                                                             */
/* More specifically, if nx is the index of an evaluation point*/
/* and ne is the index of an internal edge in an RWGSurface,   */
/* then the fields at point #nx due to a unit-strength electric*/
/* current in basis function #ne are in row                    */
/*                                                             */
/*  nx*(2*NE) + 2*ne + 0                                       */
/*                                                             */
/* while the fields at pnt #nx due to a unit-strength magnetic */ 
/* current in basis function #ne are in row                    */ 
/*                                                             */
/*  nx*(2*NE) + 2*ne + 0.                                      */
/*                                                             */
/* Here NE is the total number of edges.                       */
/*                                                             */
/* If SurfaceIndex is >=0, then only the contributions of      */
/* edges on that surface S are taken into account. Otherwise,  */
/* all edges on all surfaces in the geometry are taken into    */
/* account.                                                    */
/*                                                             */
/* Note that this routine currently assumes all evaluation     */
/* points lie in the exterior region of the geometry. Extending*/
/* to consider points lying inside bodies is possible but not  */
/* implemented at present.                                     */
/***************************************************************/
HMatrix *GetFSVMatrix(RWGGeometry *G, int SurfaceIndex,
                      HMatrix *XMatrix, const char *ColDesc,
                      cdouble Omega, bool FarField=false,
                      HMatrix *FSVMatrix=0)
{ 
  int NBF;
  if (SurfaceIndex>=0)
   NBF = G->Surfaces[SurfaceIndex]->NumBFs;
  else
   NBF = G->TotalBFs;

  int NX    = XMatrix->NR;
  int NBFNX = NBF*NX;

  /***************************************************************/
  /* reallocate FSVMatrix as necessary ***************************/
  /***************************************************************/
  if ( FSVMatrix && ( (FSVMatrix->NR!=6) || FSVMatrix->NC!=NBFNX) )
   { Warn("wrong-size FSVMatrix passed to GetFSVMatrix (reallocating)");
     delete FSVMatrix;
     FSVMatrix=0;
   };
  if (FSVMatrix==0)
   FSVMatrix = new HMatrix(6, NBFNX, LHM_COMPLEX);

  /***************************************************************/
  /* get material parameters of exterior region ******************/
  /***************************************************************/
  cdouble EpsRel, MuRel;
  G->RegionMPs[0]->GetEpsMu(Omega, &EpsRel, &MuRel);
  cdouble k     = Omega*sqrt(EpsRel*MuRel);
  cdouble ZRel  = sqrt(MuRel/EpsRel);
  cdouble IKZ   = II*k*ZVAC*ZRel;
  cdouble IKOZ  = II*k/(ZVAC*ZRel);
  
  /***************************************************************/
  /* loop over all evaluation points and all edges on all        */
  /* relevant surfaces                                           */
  /***************************************************************/
  Log("Precomputing FSV matrices (%i columns)",NBFNX);
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     if (SurfaceIndex>=0 && SurfaceIndex!=ns)
      continue;

     RWGSurface *S = G->Surfaces[ns];
     int Offset    = (SurfaceIndex>=0) ? 0 : G->BFIndexOffset[ns];
     int NE        = S->NumEdges;
     int NENX      = NE*NX;

#ifdef USE_OPENMP
     int NumThreads=GetNumThreads();
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
     for(int nenx=0; nenx<NENX; nenx++)
      { 
        int ne=nenx/NX;
        int nx=nenx%NX;

        if (nx==0) LogPercent(ne,NE,5);

        // coordinates of nxth evaluation point
        double X[3];
        XMatrix->GetEntriesD(nx,ColDesc,X);

        // reduced fields due to this edge at this eval pnt
        cdouble e[3], h[3];
        if (FarField)
         GetReducedFarFields(G, ns, ne, X, k, e, h);
        else
         GetReducedFields(G, ns, ne, X, k, e, h);


        // full fields due to this edge, populated with 
        // unit strength as an electric or magnetic current
        for(int Mu=0; Mu<3; Mu++)
         { 
           int nbf = Offset + ( (S->IsPEC) ? ne : 2*ne );
           int Column = nx*NBF + nbf;
           FSVMatrix->SetEntry( 0+Mu, Column, IKZ*e[Mu]);
           FSVMatrix->SetEntry( 3+Mu, Column,     h[Mu]);
           if (S->IsPEC) continue;
           FSVMatrix->SetEntry( 0+Mu, Column+1,     -h[Mu]);
           FSVMatrix->SetEntry( 3+Mu, Column+1, IKOZ*e[Mu]);
         };

      }; // for (int nenx=0...)

   }; // for(int ns=0; ns<G->NumSurfaces; ns++)

  return FSVMatrix;
}

/***************************************************************/
/* return the 3x3 matrices that are sandwiched between E and H */
/* three-vectors to yield the Poynting vector and the Maxwell  */
/* stress tensor.                                              */
/***************************************************************/
static double LeviCivita[3][3][3]=
{ { { 0.0,  0.0,  0.0 }, {  0.0, 0.0, +1.0 }, {  0.0, -1.0, +0.0 }  },
  { { 0.0,  0.0, -1.0 }, {  0.0, 0.0,  0.0 }, { +1.0, +0.0, +0.0 }  },
  { { 0.0, +1.0,  0.0 }, { -1.0, 0.0,  0.0 }, {  0.0,  0.0,  0.0 }  }
};
void GetNMatrices(double nHat[3], double X[3], double XTorque[3],
                  double NMatrix[NUMPFT][3][3])
{
  /***************************************************************/
  /* first matrix is the power matrix ****************************/
  /***************************************************************/
  NMatrix[SIPOWER][0][0] = NMatrix[SIPOWER][1][1] = NMatrix[SIPOWER][2][2] = 0.0;
  NMatrix[SIPOWER][1][2] = nHat[0]; NMatrix[SIPOWER][2][1] = -nHat[0];
  NMatrix[SIPOWER][2][0] = nHat[1]; NMatrix[SIPOWER][0][2] = -nHat[1];
  NMatrix[SIPOWER][0][1] = nHat[2]; NMatrix[SIPOWER][1][0] = -nHat[2];

  /***************************************************************/
  /* next three matrices are force matrices                      */
  /***************************************************************/
  NMatrix[SIXFORCE][0][0]                           =  nHat[0]; 
  NMatrix[SIXFORCE][1][1]                           = -nHat[0];
  NMatrix[SIXFORCE][2][2]                           = -nHat[0];
  NMatrix[SIXFORCE][0][1] = NMatrix[SIXFORCE][1][0] =  nHat[1];
  NMatrix[SIXFORCE][0][2] = NMatrix[SIXFORCE][2][0] =  nHat[2];
  NMatrix[SIXFORCE][1][2] = NMatrix[SIXFORCE][2][1] =  0.0;

  NMatrix[SIYFORCE][0][0]                           = -nHat[1]; 
  NMatrix[SIYFORCE][1][1]                           =  nHat[1];
  NMatrix[SIYFORCE][2][2]                           = -nHat[1];
  NMatrix[SIYFORCE][0][1] = NMatrix[SIYFORCE][1][0] =  nHat[0];
  NMatrix[SIYFORCE][1][2] = NMatrix[SIYFORCE][2][1] =  nHat[2];
  NMatrix[SIYFORCE][0][2] = NMatrix[SIYFORCE][2][0] =  0.0;

  NMatrix[SIZFORCE][0][0]                           = -nHat[2]; 
  NMatrix[SIZFORCE][1][1]                           = -nHat[2];
  NMatrix[SIZFORCE][2][2]                           =  nHat[2];
  NMatrix[SIZFORCE][0][2] = NMatrix[SIZFORCE][2][0] =  nHat[0];
  NMatrix[SIZFORCE][1][2] = NMatrix[SIZFORCE][2][1] =  nHat[1];
  NMatrix[SIZFORCE][0][1] = NMatrix[SIZFORCE][1][0] =  0.0;

  /***************************************************************/
  /* last three matrices are the torque matrices                 */
  /***************************************************************/
  double D[3], DxN[3];
  D[0] = X[0]-XTorque[0];
  D[1] = X[1]-XTorque[1];
  D[2] = X[2]-XTorque[2];
  DxN[0] = D[1]*nHat[2] - D[2]*nHat[1]; 
  DxN[1] = D[2]*nHat[0] - D[0]*nHat[2];
  DxN[2] = D[0]*nHat[1] - D[1]*nHat[0];

  for(int a=0; a<3; a++)
   for(int b=0; b<3; b++)
    for(int i=0; i<3; i++)
     { NMatrix[SIXTORQUE+i][a][b] = (a==b) ? -DxN[i] : 0.0;
       for(int j=0; j<3; j++)
        NMatrix[SIXTORQUE+i][a][b] 
         +=   LeviCivita[i][j][a]*D[j]*nHat[b]
            + LeviCivita[i][j][b]*D[j]*nHat[a];
     };

}

/***************************************************************/
/* 'hermitian vector/matrix/vector product'                    */
/***************************************************************/
double HVMVP(cdouble V1[3], double M[3][3], cdouble V2[3])
{ 
  double Sum=0.0;

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    Sum += real( conj(V1[Mu]) * M[Mu][Nu] * V2[Nu]);

  return Sum;
}

/***************************************************************/
/* Get power, force, and torque by the displaced               */
/* surface-integral method.                                    */
/***************************************************************/
void RWGGeometry::GetDSIPFT(HVector *KN, IncField *IF, cdouble Omega,
                            double PFT[NUMPFT],
                            char *BSMesh, double R, int NumPoints,
                            bool Lebedev, bool FarField,
                            GTransformation *OTGT,
                            GTransformation *GT)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (BSMesh)
   Log("Computing SIPFT over bounding surface %s...",BSMesh);
  else
   Log("Computing SIPFT: (R,NPts,Lebedev)=(%e,%i,%s)",
        R, NumPoints, Lebedev ? "true" : "false");

  /***************************************************************/
  /* get cubature-rule matrix ************************************/
  /***************************************************************/
  HMatrix *CRMatrix = GetCRMatrix(BSMesh, R, NumPoints, Lebedev,
                                  OTGT, GT);

  /***************************************************************/
  /* we assume that all cubature points lie in the same region   */
  /* of the scuff geometry, so we use the first point in the rule*/
  /* to determine which region that is and look up its eps/mu    */
  /***************************************************************/
  double X0[3];
  CRMatrix->GetEntriesD(0,"1:3",X0);
  int RegionIndex=GetRegionIndex(X0); 
  cdouble EpsRel, MuRel;
  RegionMPs[ RegionIndex ] -> GetEpsMu(Omega, &EpsRel, &MuRel);
  double EpsAbs = TENTHIRDS * real(EpsRel) / ZVAC;
  double  MuAbs = TENTHIRDS * real(MuRel) * ZVAC;

  double XTorque[3] = {0.0, 0.0, 0.0};
  if (OTGT) OTGT->Apply(XTorque);
  if (GT) GT->Apply(XTorque);

  /***************************************************************/
  /* get the total fields at the cubature points                 */
  /***************************************************************/
  HMatrix *XMatrix = CRMatrix->ExtractEntries("[:,1:3]");
  HMatrix *FMatrix;
  if (FarField)
   FMatrix = GetFarFields(this, IF, KN, Omega, XMatrix);
  else
   FMatrix = GetFields(IF, KN, Omega, XMatrix);

  /***************************************************************/
  /* loop over points in the cubature rule                       */
  /***************************************************************/
  memset(PFT, 0, NUMPFT*sizeof(double));
  for(int nr=0; nr<CRMatrix->NR; nr++)
   { 
     double w, X[3], nHat[3];
     w = CRMatrix->GetEntryD(nr, 0);
     CRMatrix->GetEntriesD(nr, "1:3", X);
     CRMatrix->GetEntriesD(nr, "4:6", nHat);

     double NMatrix[NUMPFT][3][3];
     GetNMatrices(nHat, X, XTorque, NMatrix);

     cdouble E[3], H[3];
     FMatrix->GetEntries(nr, "0:2", E);
     FMatrix->GetEntries(nr, "3:5", H);

     PFT[SIPOWER] += 0.25 * w * (  HVMVP(E, NMatrix[SIPOWER], H)
                                  -HVMVP(H, NMatrix[SIPOWER], E)
                                );

     for(int n=SIXFORCE; n<=SIZTORQUE; n++)
      PFT[n] += 0.25 * w * ( EpsAbs*HVMVP(E, NMatrix[n], E)
                             +MuAbs*HVMVP(H, NMatrix[n], H)
                           );
   };

  delete FMatrix; 
  delete XMatrix; 
  delete CRMatrix;
  
}

/***************************************************************/
/* This routine gets the contributions of a single pair of     */
/* RWG basis functions (weighted with unit strength) to the    */
/* DSIPFT.                                                     */
/*                                                             */
/* Entries is an array with enough room for 4xNUMPFT cdoubles. */
/*                                                             */
/* On return,                                                  */
/*  Entries[4*nq + 0..3] = (EE, EM, ME, MM) entries of nq-th   */
/*                         SIPFT matrix, where nq=0..7 ranges  */
/*                         over all PFT quantities.            */
/***************************************************************/
void GetDSIPFTMatrixEntries(RWGSurface *S, int neA, int neB,
                            HMatrix *CRMatrix, HMatrix *FSVMatrix,
                            cdouble EpsRel, cdouble MuRel,
                            double *XTorque, cdouble *SEntries)
{
  /***************************************************************/
  /* loop over cubature points                                   */
  /***************************************************************/
  memset(SEntries, 0, 4*NUMPFT*sizeof(cdouble));
  int NE = S->NumEdges;
  int NC = CRMatrix->NR;
  cdouble EpsAbs = TENTHIRDS * EpsRel / ZVAC;
  cdouble MuAbs  = TENTHIRDS * MuRel * ZVAC;
  for (int nc=0; nc<NC; nc++)
   {  
     double w, X[3], nHat[3];
     w = CRMatrix->GetEntryD(nc, 0);
     CRMatrix->GetEntriesD(nc, "1:3", X);
     CRMatrix->GetEntriesD(nc, "4:6", nHat);

     /***************************************************************/
     /* get 3x3 N matrices at this cubature point *******************/
     /***************************************************************/
     double NMatrix[NUMPFT][3][3];
     GetNMatrices(nHat, X, XTorque, NMatrix);

     /***************************************************************/
     /* fetch 'F' six-vectors for this cubature point           *****/
     /***************************************************************/
     int ColumnKA = nc*(2*NE) + (2*neA) + 0;
     int ColumnNA = nc*(2*NE) + (2*neA) + 1;
     int ColumnKB = nc*(2*NE) + (2*neB) + 0;
     int ColumnNB = nc*(2*NE) + (2*neB) + 1;
     cdouble FSVKA[6], FSVNA[6], FSVKB[6], FSVNB[6];
     for(int n=0; n<6; n++)
      { FSVKA[n] = FSVMatrix->GetEntry(n, ColumnKA);
        FSVNA[n] = FSVMatrix->GetEntry(n, ColumnNA);
        FSVKB[n] = FSVMatrix->GetEntry(n, ColumnKB);
        FSVNB[n] = FSVMatrix->GetEntry(n, ColumnNB);
      };
 
     /***************************************************************/
     /* entries of power matrix               ***********************/
     /***************************************************************/
     for(int m=0; m<3; m++)
      for(int n=0; n<3; n++)
       { SEntries[ 0 ] += 0.5*w*conj(FSVKA[m]) * NMatrix[SIPOWER][m][n] * FSVKB[3+n];
         SEntries[ 1 ] += 0.5*w*conj(FSVKA[m]) * NMatrix[SIPOWER][m][n] * FSVNB[3+n];
         SEntries[ 2 ] += 0.5*w*conj(FSVNA[m]) * NMatrix[SIPOWER][m][n] * FSVKB[3+n];
         SEntries[ 3 ] += 0.5*w*conj(FSVNA[m]) * NMatrix[SIPOWER][m][n] * FSVNB[3+n];
       }; 

     /***************************************************************/
     /* entries of force and torque matrices ************************/
     /***************************************************************/
     for(int nSIFT=0; nSIFT<NUMPFT-1; nSIFT++)
      { int nSIPFT=nSIFT+1;
        cdouble UpperVMVP[4]={0.0, 0.0, 0.0, 0.0};
        cdouble LowerVMVP[4]={0.0, 0.0, 0.0, 0.0};
        for(int m=0; m<3; m++)
         for(int n=0; n<3; n++)
          { 
            UpperVMVP[0] += conj(FSVKA[m]) * NMatrix[nSIPFT][m][n] * FSVKB[n];
            UpperVMVP[1] += conj(FSVKA[m]) * NMatrix[nSIPFT][m][n] * FSVNB[n];
            UpperVMVP[2] += conj(FSVNA[m]) * NMatrix[nSIPFT][m][n] * FSVKB[n];
            UpperVMVP[3] += conj(FSVNA[m]) * NMatrix[nSIPFT][m][n] * FSVNB[n];

            LowerVMVP[0] += conj(FSVKA[3+m]) * NMatrix[nSIPFT][m][n] * FSVKB[3+n];
            LowerVMVP[1] += conj(FSVKA[3+m]) * NMatrix[nSIPFT][m][n] * FSVNB[3+n];
            LowerVMVP[2] += conj(FSVNA[3+m]) * NMatrix[nSIPFT][m][n] * FSVKB[3+n];
            LowerVMVP[3] += conj(FSVNA[3+m]) * NMatrix[nSIPFT][m][n] * FSVNB[3+n];
          };
        SEntries[ 4*nSIPFT + 0 ] += 0.25*w*( EpsAbs*UpperVMVP[0] + MuAbs*LowerVMVP[0] );
        SEntries[ 4*nSIPFT + 1 ] += 0.25*w*( EpsAbs*UpperVMVP[1] + MuAbs*LowerVMVP[1] );
        SEntries[ 4*nSIPFT + 2 ] += 0.25*w*( EpsAbs*UpperVMVP[2] + MuAbs*LowerVMVP[2] );
        SEntries[ 4*nSIPFT + 3 ] += 0.25*w*( EpsAbs*UpperVMVP[3] + MuAbs*LowerVMVP[3] );

      }; //for(int nSIFT=0; nSIFT<NUMPFT-1; nSIFT++)

   }; //for (int nc=0; nc<NC; nc++)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::GetDSIPFTTrace(int SurfaceIndex, cdouble Omega,
                                 HVector *KNVector, HMatrix *SigmaMatrix,
                                 double PFT[7], double **ByEdge,
                                 char *BSMesh, double R, int NumPoints,
                                 bool Lebedev, bool FarField)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGSurface *S=Surfaces[SurfaceIndex];
  int Offset = BFIndexOffset[SurfaceIndex];
  if (BSMesh)
   Log("Computing SIPFT for surface %s over bounding surface %s...",
        S->MeshFileName, BSMesh);
  else
   Log("Computing SIPFT for surface %s: (R,NPts,Lebedev)=(%e,%i,%s)",
        S->MeshFileName, R, NumPoints, Lebedev ? "true" : "false");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (S->IsPEC) 
   { 
     Warn("DSIPFT() not implemented for PEC bodies");
     memset(PFT, 0, 7*sizeof(double)); 
     return;
   };

  cdouble Eps, Mu;
  RegionMPs[ S->RegionIndices[0] ] -> GetEpsMu(Omega, &Eps, &Mu);

  double XTorque[3] = {0.0, 0.0, 0.0};
  if (S->OTGT) S->OTGT->Transform(XTorque);
  if (S->GT) S->GT->Transform(XTorque);

  /*--------------------------------------------------------------*/
  /*- initialize edge-by-edge contributions to zero --------------*/
  /*--------------------------------------------------------------*/
  int NE = S->NumEdges;
  if (ByEdge)
   { for(int nq=0; nq<NUMPFT; nq++)
      if (ByEdge[nq])
       memset(ByEdge[nq],0,NE*sizeof(double));
   };

  /*--------------------------------------------------------------*/
  /*- fetch cubature rule and precompute field six-vectors       -*/
  /*--------------------------------------------------------------*/
  HMatrix *CRMatrix  = GetCRMatrix(BSMesh, R, NumPoints, Lebedev, S->OTGT, S->GT);
  HMatrix *FSVMatrix = GetFSVMatrix(this, SurfaceIndex, CRMatrix, "1:3",
                                    Omega, FarField);

  /*--------------------------------------------------------------*/
  /*- loop over all pairs of edges -------------------------------*/
  /*--------------------------------------------------------------*/
  double PAbs=0.0, Fx=0.0, Fy=0.0, Fz=0.0, Taux=0.0, Tauy=0.0, Tauz=0.0;
#ifndef USE_OPENMP
  if (LogLevel>=SCUFF_VERBOSE2) Log(" no multithreading...");
#else
  int NumThreads=GetNumThreads();
  if (LogLevel>=SCUFF_VERBOSE2) Log(" using %i OpenMP threads",NumThreads);
#pragma omp parallel for schedule(dynamic,1), 		\
                         num_threads(NumThreads)	\
                         reduction(+:PAbs, Fx, Fy, Fz, Taux, Tauy, Tauz)
#endif
  for(int neaneb=0; neaneb<NE*NE; neaneb++)
    { 
      int nea = neaneb / NE;
      int neb = neaneb % NE;
      if (neb<nea) continue;

      if (neb==nea) LogPercent(nea,NE,10);

      /*--------------------------------------------------------------*/
      /*- get SIPFT contributions from this pair of basis functions---*/
      /*--------------------------------------------------------------*/
      cdouble Entries[NUMPFT*4];
      GetDSIPFTMatrixEntries(S, nea, neb, CRMatrix, FSVMatrix,
                             Eps, Mu, XTorque, Entries);

      /*--------------------------------------------------------------*/
      /*- extract the surface-current coefficient either from the KN -*/
      /*- vector or the Sigma matrix                                 -*/
      /*--------------------------------------------------------------*/
      cdouble KK, KN, NK, NN;
      if (KNVector)
       { 
         cdouble kAlpha =       KNVector->GetEntry(Offset+2*nea+0);
         cdouble nAlpha = -ZVAC*KNVector->GetEntry(Offset+2*nea+1);
         cdouble kBeta  =       KNVector->GetEntry(Offset+2*neb+0);
         cdouble nBeta  = -ZVAC*KNVector->GetEntry(Offset+2*neb+1);

         KK = conj(kAlpha) * kBeta;
         KN = conj(kAlpha) * nBeta;
         NK = conj(nAlpha) * kBeta;
         NN = conj(nAlpha) * nBeta;
       }
      else
       { KK = SigmaMatrix->GetEntry(Offset+2*neb+0, Offset+2*nea+0);
         KN = SigmaMatrix->GetEntry(Offset+2*neb+1, Offset+2*nea+0);
         NK = SigmaMatrix->GetEntry(Offset+2*neb+0, Offset+2*nea+1);
         NN = SigmaMatrix->GetEntry(Offset+2*neb+1, Offset+2*nea+1);
       };

      /*--------------------------------------------------------------*/
      /*- get the contributions of this edge pair to all quantities   */
      /*--------------------------------------------------------------*/
      double DeltaPFT[NUMPFT];
      double Weight = (nea==neb) ? 1.0 : 2.0;
      for(int nq=0; nq<NUMPFT; nq++)
       DeltaPFT[nq] = Weight * real(  KK*Entries[4*nq+0] + KN*Entries[4*nq+1]
                                     +NK*Entries[4*nq+2] + NN*Entries[4*nq+3]
                                   );
      
      /*--------------------------------------------------------------*/
      /*- accumulate contributions to full sums ----------------------*/
      /*--------------------------------------------------------------*/
      PAbs +=  DeltaPFT[0];
      Fx   +=  DeltaPFT[1];
      Fy   +=  DeltaPFT[2];
      Fz   +=  DeltaPFT[3];
      Taux +=  DeltaPFT[4];
      Tauy +=  DeltaPFT[5];
      Tauz +=  DeltaPFT[6];

      /*--------------------------------------------------------------*/
      /*- accumulate contributions to by-edge sums                    */
      /*--------------------------------------------------------------*/
      if (ByEdge)
       {  
         #pragma omp critical (ByEdge)
          { if (ByEdge[0]) ByEdge[0][nea] += DeltaPFT[0];
            if (ByEdge[1]) ByEdge[1][nea] += DeltaPFT[1];
            if (ByEdge[2]) ByEdge[2][nea] += DeltaPFT[2];
            if (ByEdge[3]) ByEdge[3][nea] += DeltaPFT[3];
            if (ByEdge[4]) ByEdge[4][nea] += DeltaPFT[4];
            if (ByEdge[5]) ByEdge[5][nea] += DeltaPFT[5];
            if (ByEdge[6]) ByEdge[6][nea] += DeltaPFT[6];
          };
       };

    }; // end of multithreaded loop

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  delete FSVMatrix;
  delete CRMatrix;

  PFT[0]=PAbs;
  PFT[1]=Fx;
  PFT[2]=Fy;
  PFT[3]=Fz;
  PFT[4]=Taux;
  PFT[5]=Tauy;
  PFT[6]=Tauz;
    
}

/***************************************************************/
/* Evaluate trace formulas for the spatially-resolved fluxes   */
/* at individual points in space.                              */
/*                                                             */
/* XMatrix is an NXx3 matrix storing the cartesian coordinates */
/* of the evaluation points.                                   */
/*                                                             */
/* On return, FMatrix is an NXx21 matrix whose columns are the */
/* components of the average Poynting vector (PV) and Maxwell  */
/* stress tensor (MST) at each evaluation point.               */
/*                                                             */
/* FMatrix[nx, 0..2]  = PV_{x,y,z};                            */
/* FMatrix[nx, 3..11] = MST_{xx}, MST_{xy}, ..., MST_{zz}      */
/* FMatrix[nx,12..20] = rxMST_{xx}, rxMST_{xy}, ..., rxMST_{zz}*/
/***************************************************************/
HMatrix *RWGGeometry::GetSRFluxTrace(HMatrix *XMatrix, cdouble Omega,
                                     HVector *KNVector, HMatrix *SigmaMatrix,
                                     HMatrix *FMatrix, bool FarField)
{
  /***************************************************************/
  /* (re)allocate FMatrix as necessary ***************************/
  /***************************************************************/
  int NX = XMatrix->NR;
  if ( FMatrix && ( (FMatrix->NR != NX) || (FMatrix->NC != 3*NUMPFT) ) )
   { Warn("Wrong-size FMatrix in GetSRFluxTrace (reallocating...)");
     delete FMatrix;
     FMatrix=0;
   };
  if (FMatrix==0)
   FMatrix = new HMatrix(NX, 3*NUMPFT, LHM_REAL);

  Log("Computing spatially-resolved fluxes at %i evaluation points...",NX);

  // assume all points lie in the exterior region
  cdouble EpsRel, MuRel;
  RegionMPs[ 0 ] -> GetEpsMu(Omega, &EpsRel, &MuRel);
  double EpsAbs = TENTHIRDS * real(EpsRel) / ZVAC;
  double  MuAbs = TENTHIRDS * real(MuRel) * ZVAC;

  // assume all surfaces are non-PEC
  int NE = TotalBFs / 2;

  // compute torque about the origin of coordinates
  double XTorque[3]={0.0, 0.0, 0.0};

  /***************************************************************/
  /* precompute field six-vectors for all evaluation points      */
  /***************************************************************/
  HMatrix *FSVMatrix = GetFSVMatrix(this, -1, XMatrix, "0:2", Omega, FarField);

  /***************************************************************/
  /*- loop over all basis functions and all pairs of eval points */
  /***************************************************************/
//#ifndef USE_OPENMP
  if (LogLevel>=SCUFF_VERBOSE2) Log(" no multithreading...");
//#else
  int NumThreads=GetNumThreads();
  if (LogLevel>=SCUFF_VERBOSE2) Log(" using %i OpenMP threads",NumThreads);
//#pragma omp parallel for collapse(3), schedule(dynamic,1), num_threads(NumThreads)
//#endif
  for(int nx=0; nx<NX; nx++)
   for(int nea=0; nea<NE; nea++)
    for(int neb=0; neb<NE; neb++)
     { 
       if (nx==0) LogPercent(nea,NE,10);

       /*--------------------------------------------------------------*/
       /* extract the surface-current coefficients either from the KN  */
       /* vector or the Sigma matrix                                   */
       /*--------------------------------------------------------------*/
       cdouble KK, KN, NK, NN;
       if (KNVector)
        { 
          cdouble kAlpha =       KNVector->GetEntry(2*nea+0);
          cdouble nAlpha = -ZVAC*KNVector->GetEntry(2*nea+1);
          cdouble kBeta  =       KNVector->GetEntry(2*neb+0);
          cdouble nBeta  = -ZVAC*KNVector->GetEntry(2*neb+1);
  
          KK = conj(kAlpha) * kBeta;
          KN = conj(kAlpha) * nBeta;
          NK = conj(nAlpha) * kBeta;
          NN = conj(nAlpha) * nBeta;
        }
       else
        { KK = SigmaMatrix->GetEntry(2*neb+0, 2*nea+0);
          KN = SigmaMatrix->GetEntry(2*neb+1, 2*nea+0);
          NK = SigmaMatrix->GetEntry(2*neb+0, 2*nea+1);
          NN = SigmaMatrix->GetEntry(2*neb+1, 2*nea+1);
        };

       /*--------------------------------------------------------------*/
       /* The E-field due to basis function Alpha is                   */
       /*  kAlpha*EKAlpha + nAlpha*ENAlpha.                            */
       /* The H-field due to basis function Alpha is                   */
       /*  kAlpha*HKAlpha + nAlpha*HNAlpha.                            */
       /*--------------------------------------------------------------*/
       int ColumnKA = nx*(2*NE) + (2*nea) + 0;
       int ColumnNA = nx*(2*NE) + (2*nea) + 1;
       int ColumnKB = nx*(2*NE) + (2*neb) + 0;
       int ColumnNB = nx*(2*NE) + (2*neb) + 1;
       cdouble EKAlpha[3], HKAlpha[3], ENAlpha[3], HNAlpha[3];
       cdouble EKBeta[3], HKBeta[3], ENBeta[3], HNBeta[3];
       FSVMatrix->GetEntries( "0:2", ColumnKA, EKAlpha);
       FSVMatrix->GetEntries( "3:5", ColumnKA, HKAlpha);
       FSVMatrix->GetEntries( "0:2", ColumnNA, ENAlpha);
       FSVMatrix->GetEntries( "3:5", ColumnNA, HNAlpha);
       FSVMatrix->GetEntries( "0:2", ColumnKB, EKBeta);
       FSVMatrix->GetEntries( "3:5", ColumnKB, HKBeta);
       FSVMatrix->GetEntries( "0:2", ColumnNB, ENBeta);
       FSVMatrix->GetEntries( "3:5", ColumnNB, HNBeta);

       /*--------------------------------------------------------------*/
       /*- get the contributions of this edge pair to all quantities   */
       /*--------------------------------------------------------------*/
       double X[3];
       XMatrix->GetEntriesD(nx, "0:2", X);
       for(int Mu=0; Mu<3; Mu++)
        { 
          // N-matrices for the Muth cartesian direction
          double NMatrix[NUMPFT][3][3];
          double nHat[3]={0.0, 0.0, 0.0};
          nHat[Mu]=1.0;
          GetNMatrices(nHat, X, XTorque, NMatrix);
          
          double DeltaPFT[7];

          // power
          DeltaPFT[0] 
           =0.25*real( KK*(  HVMVP(EKAlpha, NMatrix[SIPOWER], HKBeta)
                            -HVMVP(HKAlpha, NMatrix[SIPOWER], EKBeta) )
                      +KN*(  HVMVP(EKAlpha, NMatrix[SIPOWER], HNBeta)
                            -HVMVP(HKAlpha, NMatrix[SIPOWER], ENBeta) )
                      +NK*(  HVMVP(ENAlpha, NMatrix[SIPOWER], HKBeta)
                            -HVMVP(HNAlpha, NMatrix[SIPOWER], EKBeta) )
                      +NN*(  HVMVP(ENAlpha, NMatrix[SIPOWER], HNBeta)
                            -HVMVP(HNAlpha, NMatrix[SIPOWER], ENBeta) )
                     );
 
          // force and torque 
          for(int nq=1; nq<NUMPFT; nq++)
           { DeltaPFT[nq]
              = 0.25*real( KK*( EpsAbs*HVMVP(EKAlpha, NMatrix[nq], EKBeta)
                                +MuAbs*HVMVP(HKAlpha, NMatrix[nq], HKBeta) )
                          +KN*( EpsAbs*HVMVP(EKAlpha, NMatrix[nq], ENBeta)
                                +MuAbs*HVMVP(HKAlpha, NMatrix[nq], HNBeta) )
                          +NK*( EpsAbs*HVMVP(ENAlpha, NMatrix[nq], EKBeta)
                                +MuAbs*HVMVP(HNAlpha, NMatrix[nq], HKBeta) )
                          +NN*( EpsAbs*HVMVP(ENAlpha, NMatrix[nq], ENBeta)
                                +MuAbs*HVMVP(HNAlpha, NMatrix[nq], HNBeta) )
                         );
           };
 
         #pragma omp critical (Accumulate)
          {
            FMatrix->AddEntry(nx, 0  + Mu, DeltaPFT[0]); // PV[Mu]
            FMatrix->AddEntry(nx, 3  + 3*Mu + 0, DeltaPFT[1]); // MST[Mu][0]
            FMatrix->AddEntry(nx, 3  + 3*Mu + 1, DeltaPFT[2]); // MST[Mu][1]
            FMatrix->AddEntry(nx, 3  + 3*Mu + 2, DeltaPFT[3]); // MST[Mu][2]
            FMatrix->AddEntry(nx, 12 + 3*Mu + 0, DeltaPFT[4]); // rxMST[Mu][0]
            FMatrix->AddEntry(nx, 12 + 3*Mu + 1, DeltaPFT[5]); // rxMST[Mu][1]
            FMatrix->AddEntry(nx, 12 + 3*Mu + 2, DeltaPFT[6]); // rxMST[Mu][2]
          };
 
        }; // for(int Mu=0; Mu<3; Mu++)

    }; // for(nxne==...)

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  delete FSVMatrix;
  return FMatrix;

} // routine GetSRFluxTrace

}
