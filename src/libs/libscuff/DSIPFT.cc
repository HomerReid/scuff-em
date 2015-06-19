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

#ifdef USE_OPENMP
 #include <omp.h>
#endif
 
namespace scuff{

void GetReducedFarFields(RWGSurface *S, const int ne,
                         const double X0[3], const cdouble k,
                         cdouble e[3], cdouble h[3]);

#define II cdouble(0.0,1.0) 
#define TENTHIRDS 3.33333333333333333333333

#define SIPOWER   0
#define SIPSCAT   1
#define SIXFORCE  2
#define SIYFORCE  3
#define SIZFORCE  4
#define SIXTORQUE 5
#define SIYTORQUE 6
#define SIZTORQUE 7

#define NUMSRFLUX 12

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetReducedFields(RWGSurface *S, int ne,
                      double X0[3],  cdouble k,
                      cdouble e[3], cdouble h[3])
{
  cdouble a[3], curla[3], Gradp[3];

  S->GetReducedPotentials(ne, X0, k, 0, a, curla, Gradp);

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
       { GInt[Mu] += w*( GP[Mu][Nu]*XmQP[Nu] - GM[Mu][Nu]*XmQM[Nu] );
         CInt[Mu] += w*( CP[Mu][Nu]*XmQP[Nu] - CM[Mu][Nu]*XmQM[Nu] );
       };

   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetFarFields(RWGGeometry *G, IncField *IF, HVector *KN,
                      cdouble Omega, HMatrix *XMatrix)
{ 
  int NX = XMatrix->NR;
  Log("Computing far fields at %i points...",NX);
  HMatrix *FMatrix = new HMatrix(NX, 6, LHM_COMPLEX);

  cdouble EpsOut, MuOut;
  G->RegionMPs[0]->GetEpsMu(Omega, &EpsOut, &MuOut);
  cdouble k    = sqrt(EpsOut*MuOut)*Omega;
  cdouble ZRel = sqrt(MuOut/EpsOut);
  cdouble IKZ  = II*k*ZVAC*ZRel;
  cdouble IKOZ = II*k/(ZVAC*ZRel);

#ifdef USE_OPENMP
  int NumThreads=GetNumThreads();
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nx=0; nx<NX; nx++)
   { 
     double X[3];
     XMatrix->GetEntriesD(nx, "0:2", X);

     if ( !(G->PointInRegion(0,X)) )
      ErrExit("%s:%i: points must lie in exterior region",__FILE__,__LINE__);

     cdouble EH[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
     for(int ns=0; ns<G->NumSurfaces; ns++)
      { 
        RWGSurface *S=G->Surfaces[ns];
        int Offset = G->BFIndexOffset[ns];

        for(int ne=0; ne<S->NumEdges; ne++)
         { 
           cdouble kAlpha, nAlpha;
           if (S->IsPEC)
            { kAlpha = KN->GetEntry(Offset + ne);
              nAlpha = 0.0;
            }
           else
            { kAlpha =       KN->GetEntry(Offset + 2*ne + 0);
              nAlpha = -ZVAC*KN->GetEntry(Offset + 2*ne + 1);
            };

           cdouble e[3], h[3];
           GetReducedFarFields(S, ne, X, k, e, h);

           EH[0] += IKZ*kAlpha*e[0] - nAlpha*h[0];
           EH[1] += IKZ*kAlpha*e[1] - nAlpha*h[1];
           EH[2] += IKZ*kAlpha*e[2] - nAlpha*h[2];
           EH[3] +=     kAlpha*h[0] + IKOZ*nAlpha*e[0];
           EH[4] +=     kAlpha*h[1] + IKOZ*nAlpha*e[1];
           EH[5] +=     kAlpha*h[2] + IKOZ*nAlpha*e[2];

         }; // for(int ne=0; ne<S->NumEdges; ne++)
      }; // for(int ns=0; ns<G->NumSurfaces; ns++)

     if (IF) 
      { cdouble EHI[6];
        IF->GetFields(X, EHI);
        for(int Mu=0; Mu<6; Mu++) 
         EH[Mu] += EHI[Mu];
      };

     FMatrix->SetEntries(nx, "0:5", EH);

   }; // for(int nx=0; nx<XMatrix->NR; nx++)

  return FMatrix;
}

/***************************************************************/
/* SCR matrix stands for 'surface cubature rule matrix.'       */
/*                                                             */
/* An SCR matrix is an NCx7 matrix. The nth row has entries    */
/*                                                             */
/*    x_n y_n z_n nx_n ny_n nz_n w_n                           */
/*                                                             */
/* where (x,y,z) is the nth cubature point, (nx, ny, nz) are   */
/* the components of the outward-pointing normal at that point,*/
/* and w_n is the cubature weight.                             */
/*                                                             */
/* If BS is non-null, then it describes a closed triangulated  */
/* surface, and the cubature rule evaluates an integral over   */
/* this surface by doing a one-point cubature over each panel. */
/*                                                             */
/* Otherwise, the cubature rule describes a cubature rule with */
/* NumPoints cubature points over a sphere of radius R.        */
/*                                                             */
/* If UseCCQ is false, this is a Lebedev cubature rule. In     */
/* this case, NumPoints must be one of the numbers of cubature */
/* points supported by the GetLebedevRule() routine in         */
/* libTriInt.                                                  */
/*                                                             */
/* Otherwise (UseCCQ==true) the cubature rule is a product     */
/* rule with a Clenshaw-Curtis grid in the Theta direction and */
/* an evenly-spaced grid in the Phi direction, and NumPoints   */
/* should be an odd integer between 9 and 99 inclusive.        */
/*                                                             */
/* If GT is non-null, then each cubature point and             */
/* normal vector is transformed by GT.                         */
/***************************************************************/
HMatrix *GetSCRMatrix(char *BSMesh, double R, int NumPoints, 
                      bool UseCCQ, GTransformation *GT)
{
  HMatrix *SCRMatrix;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (BSMesh)
   { 
     RWGSurface *BS=new RWGSurface(BSMesh);
     SCRMatrix = new HMatrix(BS->NumPanels, 7);
     for(int np=0; np<BS->NumPanels; np++)
      { 
        double *X0   = BS->Panels[np]->Centroid;
        double *ZHat = BS->Panels[np]->ZHat;
        
        // define Sign= \pm 1 such that Sign*ZHat is the
        // outward-pointing surface normal
        double Sign=1.0;
        double XP[3];
        double LengthScale= BS->Panels[np]->Radius;
        VecScaleAdd(X0,0.1*LengthScale,ZHat,XP);
        if ( BS->Contains(XP) )
         Sign=-1.0;

        SCRMatrix->SetEntry(np, 0, X0[0]);
        SCRMatrix->SetEntry(np, 1, X0[1]);
        SCRMatrix->SetEntry(np, 2, X0[2]);
        SCRMatrix->SetEntry(np, 3, Sign*ZHat[0]);
        SCRMatrix->SetEntry(np, 4, Sign*ZHat[1]);
        SCRMatrix->SetEntry(np, 5, Sign*ZHat[2]);
        SCRMatrix->SetEntry(np, 6, BS->Panels[np]->Area);
      };
     delete BS;
   }
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  else if (!UseCCQ)
   { double *LRule = GetLebedevRule(NumPoints);
     if (LRule==0) ErrExit("no Lebedev rule with %i points",NumPoints);
     SCRMatrix = new HMatrix(NumPoints, 7);
     for(int np=0; np<NumPoints; np++)
      { SCRMatrix->SetEntry(np,0, R*LRule[4*np + 0]);
        SCRMatrix->SetEntry(np,1, R*LRule[4*np + 1]);
        SCRMatrix->SetEntry(np,2, R*LRule[4*np + 2]);
        SCRMatrix->SetEntry(np,3, LRule[4*np + 0]);
        SCRMatrix->SetEntry(np,4, LRule[4*np + 1]);
        SCRMatrix->SetEntry(np,5, LRule[4*np + 2]);
        SCRMatrix->SetEntry(np,6, R*R*LRule[4*np + 3]);
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
     SCRMatrix = new HMatrix(TotalPoints, 7);
     double dPhi= 2.0*M_PI / NumPhiPoints;

     double *CCRule = GetCCRule(NumThetaPoints);
     if (CCRule==0) ErrExit("invalid value of NumPoints in GetSCRMatrix");
     for(int npTheta=0; npTheta<NumThetaPoints; npTheta++)
      { 
        double CosTheta = CCRule[2*npTheta+0];
        double w = CCRule[2*npTheta+1];
        double SinTheta = sqrt(1.0-CosTheta*CosTheta);

        // north, south pole
        if ( npTheta==0 || npTheta==(NumThetaPoints-1) )
         { int Index  = (npTheta==0) ? 0 : TotalPoints-1;
           SCRMatrix->SetEntry(Index, 0, 0.0);
           SCRMatrix->SetEntry(Index, 1, 0.0);
           SCRMatrix->SetEntry(Index, 2, R*CosTheta);
           SCRMatrix->SetEntry(Index, 3, 0.0);
           SCRMatrix->SetEntry(Index, 4, 0.0);
           SCRMatrix->SetEntry(Index, 5, CosTheta);
           SCRMatrix->SetEntry(Index, 6, 2.0*M_PI*R*R*w);
         }
        else
         { for(int npPhi=0; npPhi<NumPhiPoints; npPhi++)
            { int Index  = 1 + NumPhiPoints*(npTheta-1) + npPhi;
              double Phi = npPhi * dPhi;
              SCRMatrix->SetEntry(Index, 0, R*SinTheta*cos(Phi));
              SCRMatrix->SetEntry(Index, 1, R*SinTheta*sin(Phi));
              SCRMatrix->SetEntry(Index, 2, R*CosTheta);
              SCRMatrix->SetEntry(Index, 3, SinTheta*cos(Phi));
              SCRMatrix->SetEntry(Index, 4, SinTheta*sin(Phi));
              SCRMatrix->SetEntry(Index, 5, CosTheta);
              SCRMatrix->SetEntry(Index, 6, R*R*w*dPhi);
            };
         };
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (GT)
   { 
     for(int nr=0; nr<SCRMatrix->NR; nr++)
      { 
        double X[3], nHat[3], XP[3];

        SCRMatrix->GetEntriesD(nr,"0:2",X);
        SCRMatrix->GetEntriesD(nr,"3:5",nHat);
        VecScaleAdd(X, 1.0, nHat, XP);

        GT->Apply(X); GT->Apply(XP);

        VecSub(XP, X, nHat);

        SCRMatrix->SetEntriesD(nr,"0:2",X);
        SCRMatrix->SetEntriesD(nr,"3:5",nHat);
      };
   };

  return SCRMatrix;

}

/***************************************************************/
/* Compute a matrix of field six-vectors.                      */
/*                                                             */
/* Each column of an FSVMatrix is a 6-component vector giving  */
/* the E and H fields, at a single point in space, due to a    */
/* single RWG basis function.                                  */
/*                                                             */
/* If SurfaceIndex is >=0, then only the contributions of      */
/* RWG functions on surface #SurfaceIndex are taken into       */
/* account. Otherwise, all RWG functions on all surfaces in    */
/* the geometry are taken into account.                        */
/*                                                             */
/* The cartesian coordinates of the evaluation points are the  */
/* first three columns of XMatrix.                             */
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
/*  nx*(NBF) + 2*ne + 0                                        */
/*                                                             */
/* while the fields at pnt #nx due to a unit-strength magnetic */ 
/* current in basis function #ne are in row                    */ 
/*                                                             */
/*  nx*(NBF) + 2*ne + 1.                                       */
/*                                                             */
/* Here NBF is the total number of basis functions (in the     */
/* surface, if SurfaceIndex>=0, or in the entire geometry      */
/* (if SurfaceIndex==-1).                                      */
/***************************************************************/
HMatrix *GetFSVMatrix(RWGGeometry *G, int SurfaceIndex,
                      HMatrix *XMatrix, cdouble Omega,
                      bool FarField=false, HMatrix *FSVMatrix=0)
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
        XMatrix->GetEntriesD(nx,"0:2",X);

        // figure out what region this point is in and 
        // whether or not the surface in question bounds that
        // region
        int RegionIndex=G->GetRegionIndex(X);
        double Sign;
        if (S->RegionIndices[0]==RegionIndex)
         Sign=1.0;
        else if (S->RegionIndices[1]==RegionIndex)
         Sign=-1.0;
        else  
         continue; // surface does not contribute to fields at X

        cdouble EpsRel, MuRel;
        G->RegionMPs[RegionIndex]->GetEpsMu(Omega, &EpsRel, &MuRel);
        cdouble k     = Omega*sqrt(EpsRel*MuRel);
        cdouble ZRel  = sqrt(MuRel/EpsRel);
        cdouble IKZ   = II*k*ZVAC*ZRel;
        cdouble IKOZ  = II*k/(ZVAC*ZRel);

        // reduced fields due to this edge at this eval pnt
        cdouble e[3], h[3];
        if (FarField)
         GetReducedFarFields(S, ne, X, k, e, h);
        else
         GetReducedFields(S, ne, X, k, e, h);

        // full fields due to this edge, populated with
        // unit strength as an electric or magnetic current
        for(int Mu=0; Mu<3; Mu++)
         { 
           int nbf = Offset + ( (S->IsPEC) ? ne : 2*ne );
           int Column = nx*NBF + nbf;
           FSVMatrix->SetEntry( 0+Mu, Column,    Sign*IKZ*e[Mu]);
           FSVMatrix->SetEntry( 3+Mu, Column,        Sign*h[Mu]);
           if (S->IsPEC) continue;
           FSVMatrix->SetEntry( 0+Mu, Column+1,     -Sign*h[Mu]);
           FSVMatrix->SetEntry( 3+Mu, Column+1, Sign*IKOZ*e[Mu]);
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
                  double NMatrix[NUMPFT][3][3], 
                  bool *NeedQuantity=0)
{
  if ( !NeedQuantity || NeedQuantity[SIPOWER] )
   { NMatrix[SIPOWER][0][0] = NMatrix[SIPOWER][1][1] = NMatrix[SIPOWER][2][2] = 0.0;
     NMatrix[SIPOWER][1][2] = nHat[0]; NMatrix[SIPOWER][2][1] = -nHat[0];
     NMatrix[SIPOWER][2][0] = nHat[1]; NMatrix[SIPOWER][0][2] = -nHat[1];
     NMatrix[SIPOWER][0][1] = nHat[2]; NMatrix[SIPOWER][1][0] = -nHat[2];
   };

  if ( !NeedQuantity || NeedQuantity[SIXFORCE] )
   { NMatrix[SIXFORCE][0][0]                           =  nHat[0]; 
     NMatrix[SIXFORCE][1][1]                           = -nHat[0];
     NMatrix[SIXFORCE][2][2]                           = -nHat[0];
     NMatrix[SIXFORCE][0][1] = NMatrix[SIXFORCE][1][0] =  nHat[1];
     NMatrix[SIXFORCE][0][2] = NMatrix[SIXFORCE][2][0] =  nHat[2];
     NMatrix[SIXFORCE][1][2] = NMatrix[SIXFORCE][2][1] =  0.0;
   };

  if ( !NeedQuantity || NeedQuantity[SIYFORCE] )
   { NMatrix[SIYFORCE][0][0]                           = -nHat[1]; 
     NMatrix[SIYFORCE][1][1]                           =  nHat[1];
     NMatrix[SIYFORCE][2][2]                           = -nHat[1];
     NMatrix[SIYFORCE][0][1] = NMatrix[SIYFORCE][1][0] =  nHat[0];
     NMatrix[SIYFORCE][1][2] = NMatrix[SIYFORCE][2][1] =  nHat[2];
     NMatrix[SIYFORCE][0][2] = NMatrix[SIYFORCE][2][0] =  0.0;
   };

  if ( !NeedQuantity || NeedQuantity[SIZFORCE] )
   { NMatrix[SIZFORCE][0][0]                           = -nHat[2];
     NMatrix[SIZFORCE][1][1]                           = -nHat[2];
     NMatrix[SIZFORCE][2][2]                           =  nHat[2];
     NMatrix[SIZFORCE][0][2] = NMatrix[SIZFORCE][2][0] =  nHat[0];
     NMatrix[SIZFORCE][1][2] = NMatrix[SIZFORCE][2][1] =  nHat[1];
     NMatrix[SIZFORCE][0][1] = NMatrix[SIZFORCE][1][0] =  0.0;
   };

  if (      NeedQuantity
       && ( NeedQuantity[SIXTORQUE]==false )
       && ( NeedQuantity[SIYTORQUE]==false )
       && ( NeedQuantity[SIZTORQUE]==false )
     ) return;

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

  for(int i=0; i<3; i++)
   { 
     if ( NeedQuantity && (NeedQuantity[SIXTORQUE+i]==false) ) 
      continue;

     for(int a=0; a<3; a++)
      for(int b=0; b<3; b++)
       { NMatrix[SIXTORQUE+i][a][b] = (a==b) ? -DxN[i] : 0.0;
         for(int j=0; j<3; j++)
          NMatrix[SIXTORQUE+i][a][b] 
           +=   LeviCivita[i][j][a]*D[j]*nHat[b]
              + LeviCivita[i][j][b]*D[j]*nHat[a];
       };
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
void GetDSIPFT(RWGGeometry *G, cdouble Omega, double *kBloch,
               HVector *KN, IncField *IF,
               double PFT[NUMPFT],
               char *BSMesh, double R, int NumPoints,
               bool UseCCQ, bool FarField,
               char *PlotFileName, GTransformation *GT)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (BSMesh)
   Log("Computing DSIPFT over bounding surface %s...",BSMesh);
  else
   Log("Computing DSIPFT: (R,NPts,Lebedev)=(%e,%i,%s)",
        R, NumPoints, UseCCQ ? "false" : "true");

  /***************************************************************/
  /* get cubature-rule matrix ************************************/
  /***************************************************************/
  HMatrix *SCRMatrix = GetSCRMatrix(BSMesh, R, NumPoints, UseCCQ, GT);

  /***************************************************************/
  /* we assume that all cubature points lie in the same region   */
  /* of the scuff geometry, so we use the first point in the rule*/
  /* to determine which region that is and look up its eps/mu    */
  /***************************************************************/
  double X0[3];
  SCRMatrix->GetEntriesD(0,"0:2",X0);
  int RegionIndex=G->GetRegionIndex(X0); 
  cdouble EpsRel, MuRel;
  G->RegionMPs[ RegionIndex ] -> GetEpsMu(Omega, &EpsRel, &MuRel);
  double EpsAbs = TENTHIRDS * real(EpsRel) / ZVAC;
  double  MuAbs = TENTHIRDS * real(MuRel) * ZVAC;

  double XTorque[3] = {0.0, 0.0, 0.0};
  if (GT) GT->Apply(XTorque);

  /***************************************************************/
  /* get the scattered and total fields at the cubature points   */
  /***************************************************************/
  HMatrix *FMatrixScat;
  if (FarField && kBloch==0)
   FMatrixScat = GetFarFields(G, 0, KN, Omega, SCRMatrix);
  else
   FMatrixScat = G->GetFields(0, KN, Omega, kBloch, SCRMatrix);

  HMatrix *FMatrix;
  if (IF==0)
   { 
     FMatrix = FMatrixScat;
   }
  else
   { FMatrix = G->GetFields(IF, 0, Omega, kBloch, SCRMatrix);
     FMatrix->AddBlock(FMatrixScat, 0, 0);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double **ByPanel=0;
  RWGSurface *BS=0;
  if (BSMesh && PlotFileName)
   { BS=new RWGSurface(BSMesh);
     if (GT) BS->Transform(GT);
     ByPanel = (double **)mallocEC(NUMPFT*sizeof(double *));
     ByPanel[0] = (double *)mallocEC(NUMPFT*(BS->NumPanels)*sizeof(double));
     for(int nq=1; nq<NUMPFT; nq++)
      ByPanel[nq] = ByPanel[nq-1] + BS->NumPanels;
   };

  /***************************************************************/
  /* loop over points in the cubature rule                       */
  /***************************************************************/
  memset(PFT, 0, NUMPFT*sizeof(double));
  for(int nr=0; nr<SCRMatrix->NR; nr++)
   { 
     double w, X[3], nHat[3];
     SCRMatrix->GetEntriesD(nr, "0:2", X);
     SCRMatrix->GetEntriesD(nr, "3:5", nHat);
     w = SCRMatrix->GetEntryD(nr, 6);

     double NMatrix[NUMPFT][3][3];
     GetNMatrices(nHat, X, XTorque, NMatrix);

     cdouble ES[3], HS[3], E[3], H[3];
     FMatrixScat->GetEntries(nr, "0:2", ES);
     FMatrixScat->GetEntries(nr, "3:5", HS);
     FMatrix->GetEntries(nr, "0:2", E);
     FMatrix->GetEntries(nr, "3:5", H);

     // absorbed power 
     double dP = -0.25 * w * (  HVMVP(E, NMatrix[SIPOWER], H)
                               -HVMVP(H, NMatrix[SIPOWER], E)
                             );
     PFT[SCUFF_PABS] += dP;
     if (ByPanel) ByPanel[SCUFF_PABS][ nr ] = dP;

     // scattered power
     dP = 0.25 * w * (  HVMVP(ES, NMatrix[SIPOWER], HS)
                       -HVMVP(HS, NMatrix[SIPOWER], ES)
                     );
     PFT[SCUFF_PSCAT] += dP;
     if (ByPanel) ByPanel[SCUFF_PSCAT][ nr ] = dP;

     // force and torque
     double dFT[NUMPFT];
     for(int nq=2; nq<NUMPFT; nq++)
      { dFT[nq] = 0.25 * w * ( EpsAbs*HVMVP(E, NMatrix[nq], E)
                               +MuAbs*HVMVP(H, NMatrix[nq], H)
                             );
        PFT[nq] += dFT[nq];
        if (ByPanel) ByPanel[nq][ nr ] = dFT[nq];
      };
   };

  if (FMatrix!=FMatrixScat) delete FMatrix;
  delete FMatrixScat;
  delete SCRMatrix;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (ByPanel)
   { 
     static const char *PFTNames[8]
      ={"PAbs","PScat","FX","FY","FZ","TX","TY","TZ"};

     FILE *f=fopen(PlotFileName,"a");
     for(int nq=0; nq<NUMPFT; nq++)
      BS->PlotScalarDensity(ByPanel[nq], false, PlotFileName,
                            "%s(%s)",PFTNames[nq],z2s(Omega));
     fclose(f);
     
     free(ByPanel[0]);
     free(ByPanel);
     delete BS;
   };

}

/***************************************************************/
/* This routine gets the contributions of a single pair of     */
/* RWG basis functions to the DSIPFT.                          */
/*                                                             */
/* DeltaPFT is an array with enough room for NUMPFT doubles.   */
/*                                                             */
/* If NeedQuantity[nq] == false then the calculation of that   */
/* quantity is skipped.                                        */
/***************************************************************/
void GetEdgeEdgeDSIPFT(RWGGeometry *G,
                       int nsa, int nea, int nsb, int neb,
                       cdouble KK, cdouble KN, cdouble NK, cdouble NN,
                       HMatrix *SCRMatrix, HMatrix *FSVMatrix,
                       cdouble EpsRel, cdouble MuRel,
                       double *XTorque,
                       double DeltaPFT[NUMPFT],
                       bool NeedQuantity[NUMPFT],
                       double **ByPanel)
{
  memset(DeltaPFT, 0, NUMPFT*sizeof(double));
 
  // MicroByPanel is an array that stores the contributions
  // of each panel on a meshed bounding surface to each
  // PFT quantity
  double *MicroByPanel=0;
  if (ByPanel)
   MicroByPanel = new double[ NUMPFT*(SCRMatrix->NR) ];

  double Weight = (nsa==nsb && neb==nea) ? 1.0 : 2.0;

  /***************************************************************/
  /* loop over cubature points                                   */
  /***************************************************************/
  int OffsetA    = G->BFIndexOffset[nsa];
  int OffsetB    = G->BFIndexOffset[nsb];
  int NETot      = G->TotalEdges;

  int NC         = SCRMatrix->NR;
  double EpsAbs  = TENTHIRDS * real(EpsRel) / ZVAC;
  double MuAbs   = TENTHIRDS * real(MuRel) * ZVAC;
  for (int nc=0; nc<NC; nc++)
   {  
     double X[3], nHat[3];
     SCRMatrix->GetEntriesD(nc, "0:2", X);
     SCRMatrix->GetEntriesD(nc, "3:5", nHat);
     double w = Weight * SCRMatrix->GetEntryD(nc, 6);

     /***************************************************************/
     /* get 3x3 N matrices at this cubature point *******************/
     /***************************************************************/
     double NMatrix[NUMPFT][3][3];
     GetNMatrices(nHat, X, XTorque, NMatrix, NeedQuantity);

     /***************************************************************/
     /* fetch 'F' six-vectors for this cubature point           *****/
     /***************************************************************/
     int ColumnKA = nc*(2*NETot) + OffsetA + (2*nea) + 0;
     int ColumnNA = nc*(2*NETot) + OffsetA + (2*nea) + 1;
     int ColumnKB = nc*(2*NETot) + OffsetB + (2*neb) + 0;
     int ColumnNB = nc*(2*NETot) + OffsetB + (2*neb) + 1;
     cdouble FSVKA[6], FSVNA[6], FSVKB[6], FSVNB[6];
     for(int n=0; n<6; n++)
      { FSVKA[n] = FSVMatrix->GetEntry(n, ColumnKA);
        FSVNA[n] = FSVMatrix->GetEntry(n, ColumnNA);
        FSVKB[n] = FSVMatrix->GetEntry(n, ColumnKB);
        FSVNB[n] = FSVMatrix->GetEntry(n, ColumnNB);
      };

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     double EpsEE[3][3], MuHH[3][3], EH[3][3];
     for(int m=0; m<3; m++)
      for(int n=0; n<3; n++)
       { EpsEE[m][n] = EpsAbs * real( KK*conj(FSVKA[m])*FSVKB[n]
                                     +KN*conj(FSVKA[m])*FSVNB[n]
                                     +NK*conj(FSVNA[m])*FSVKB[n]
                                     +NN*conj(FSVNA[m])*FSVNB[n]
                                    );
     
         MuHH[m][n] = MuAbs * real(  KK*conj(FSVKA[3+m])*FSVKB[3+n]
                                    +KN*conj(FSVKA[3+m])*FSVNB[3+n]
                                    +NK*conj(FSVNA[3+m])*FSVKB[3+n] 
                                    +NN*conj(FSVNA[3+m])*FSVNB[3+n]
                                  );
 
         EH[m][n] = real( KK*conj(FSVKA[m])*FSVKB[3+n]
                         +KN*conj(FSVKA[m])*FSVNB[3+n]
                         +NK*conj(FSVNA[m])*FSVKB[3+n]
                         +NN*conj(FSVNA[m])*FSVNB[3+n]
                        );
       };
 
     /***************************************************************/
     /* absorbed power **********************************************/
     /***************************************************************/
     if ( NeedQuantity[SIPOWER] )
      { 
        // this could be accelerated by exploiting the
        // structure of the NMatrix
        double MicroDelta=0.0;
        for(int m=0; m<3; m++)
         for(int n=0; n<3; n++)
          MicroDelta += 0.5*w*NMatrix[SIPOWER][m][n]*EH[m][n];

        DeltaPFT[SCUFF_PABS] += MicroDelta;
        if (ByPanel) MicroByPanel[SCUFF_PABS*NC + nc]=MicroDelta;
      };

     DeltaPFT[SCUFF_PSCAT]=0.0; // no scattered power in the DSI formalism

     /***************************************************************/
     /* entries of force and torque matrices ************************/
     /***************************************************************/
     for(int nq=SIXFORCE; nq<=SIZTORQUE; nq++)
      { 
        if ( NeedQuantity[nq]==false ) continue;

        double MicroDelta=0.0;
        for(int m=0; m<3; m++)
         for(int n=0; n<3; n++)
          MicroDelta+=0.25*w*NMatrix[nq][m][n]*(EpsEE[m][n] + MuHH[m][n]);

        DeltaPFT[nq] += MicroDelta;
        if (ByPanel) MicroByPanel[nq*NC + nc]=MicroDelta;

      }; //for(int nSIFT=0; nSIFT<NUMPFT-1; nSIFT++)

   }; //for (int nc=0; nc<NC; nc++)

  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  if (ByPanel)
   { 
#ifdef USE_OPENMP
#pragma omp critical(Accumulate)
#endif
      {
        for(int nq=0; nq<NUMPFT; nq++)
         { if (NeedQuantity[nq]==false)
            continue;
           for(int nc=0; nc<NC; nc++)
            ByPanel[nq][nc] += MicroByPanel[nq*NC + nc];
         };
      };
     delete[] MicroByPanel;
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetDSIPFTTrace(RWGGeometry *G, cdouble Omega,
                    HMatrix *RytovMatrix,
                    double PFT[NUMPFT], bool NeedQuantity[NUMPFT],
                    char *BSMesh, double R, int NumPoints,
                    bool UseCCQ, bool FarField,
                    char *PlotFileName, GTransformation *GT)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Log("Computing DSIPFT: ");
  if (BSMesh)
   LogC(" (BS mesh %s)",BSMesh);
  else
   LogC(" (R=%e, NC=%i, Lebedev=%i, FarField=%i)",R,NumPoints,!UseCCQ,FarField);

  // assume that the integration surface lies in the exterior medium
  cdouble Eps, Mu;
  G->RegionMPs[ 0 ] -> GetEpsMu(Omega, &Eps, &Mu);

  double XTorque[3] = {0.0, 0.0, 0.0};
  if (GT) GT->Transform(XTorque);

  /*--------------------------------------------------------------*/
  /*- fetch cubature rule and precompute field six-vectors       -*/
  /*--------------------------------------------------------------*/
  HMatrix *SCRMatrix = GetSCRMatrix(BSMesh, R, NumPoints, UseCCQ, GT);
  HMatrix *FSVMatrix = GetFSVMatrix(G, -1, SCRMatrix, Omega, FarField);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double **ByPanel=0;
  RWGSurface *BS=0;
  if (BSMesh && PlotFileName)
   { BS=new RWGSurface(BSMesh);
     if (GT) BS->Transform(GT);
     ByPanel = (double **)mallocEC(NUMPFT*sizeof(double *));
     ByPanel[0] = (double *)mallocEC(NUMPFT*(BS->NumPanels)*sizeof(double));
     for(int nq=1; nq<NUMPFT; nq++)
      ByPanel[nq] = ByPanel[nq-1] + BS->NumPanels;
   };

  /*--------------------------------------------------------------*/
  /*- loop over all pairs of edges on all surfaces for which the  */
  /*- outer region is the exterior region of the vacuum           */
  /*--------------------------------------------------------------*/
  int NS=G->NumSurfaces;
  double PAbs=0.0, Fx=0.0, Fy=0.0, Fz=0.0, Taux=0.0, Tauy=0.0, Tauz=0.0;
  for(int nsa=0; nsa<NS; nsa++)
   for(int nsb=nsa; nsb<NS; nsb++)
    { 
      if (G->Surfaces[nsa]->RegionIndices[0]!=0) continue;
      if (G->Surfaces[nsb]->RegionIndices[0]!=0) continue;
      int NEA=G->Surfaces[nsa]->NumEdges;
      int NEB=G->Surfaces[nsb]->NumEdges;
      int OffsetA=G->BFIndexOffset[nsa];
      int OffsetB=G->BFIndexOffset[nsb];

#ifndef USE_OPENMP
  if (G->LogLevel>=SCUFF_VERBOSE2) Log(" no multithreading...");
#else
  int NumThreads=GetNumThreads();
  if (G->LogLevel>=SCUFF_VERBOSE2) Log(" using %i OpenMP threads",NumThreads);
#pragma omp parallel for schedule(dynamic,1), 		\
                         num_threads(NumThreads),	\
                         collapse(2),			\
                         reduction(+:PAbs, Fx, Fy, Fz, Taux, Tauy, Tauz)
#endif
    for(int nea=0; nea<NEA; nea++)
     for(int neb=0; neb<NEB; neb++)
      { 
        if (nsb==nsa && neb<nea) continue;

        if (neb==nea) LogPercent(OffsetA+2*nea,G->TotalBFs,10);

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        cdouble KK = RytovMatrix->GetEntry(OffsetB+2*neb+0, OffsetA+2*nea+0);
        cdouble KN = RytovMatrix->GetEntry(OffsetB+2*neb+1, OffsetA+2*nea+0);
        cdouble NK = RytovMatrix->GetEntry(OffsetB+2*neb+0, OffsetA+2*nea+1);
        cdouble NN = RytovMatrix->GetEntry(OffsetB+2*neb+1, OffsetA+2*nea+1);

        /*--------------------------------------------------------------*/
        /*- get DSIPFT contributions from this pair of basis functions -*/
        /*--------------------------------------------------------------*/
        double DeltaPFT[NUMPFT];
        GetEdgeEdgeDSIPFT(G, nsa, nea, nsb, neb, KK, KN, NK, NN,
                          SCRMatrix, FSVMatrix, Eps, Mu, XTorque,
                          DeltaPFT, NeedQuantity, ByPanel);

        /*--------------------------------------------------------------*/
        /*- accumulate contributions to full sums ----------------------*/
        /*--------------------------------------------------------------*/
        PAbs +=  DeltaPFT[SCUFF_PABS];
        Fx   +=  DeltaPFT[SCUFF_XFORCE];
        Fy   +=  DeltaPFT[SCUFF_YFORCE];
        Fz   +=  DeltaPFT[SCUFF_ZFORCE];
        Taux +=  DeltaPFT[SCUFF_XTORQUE];
        Tauy +=  DeltaPFT[SCUFF_YTORQUE];
        Tauz +=  DeltaPFT[SCUFF_ZTORQUE];

      }; // for(int nea=... for (int neb=...

    }; // for (int nsa... for (int nsb=...

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  delete FSVMatrix;
  delete SCRMatrix;

  PFT[SCUFF_PABS]=PAbs;
  PFT[SCUFF_PSCAT]=0.0;
  PFT[SCUFF_XFORCE]=Fx;
  PFT[SCUFF_YFORCE]=Fy;
  PFT[SCUFF_ZFORCE]=Fz;
  PFT[SCUFF_XTORQUE]=Taux;
  PFT[SCUFF_YTORQUE]=Tauy;
  PFT[SCUFF_ZTORQUE]=Tauz;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (ByPanel)
   { 
     static const char *PFTNames[8]
      ={"PAbs","PScat", "FX","FY","FZ","TX","TY","TZ"};

     FILE *f=fopen(PlotFileName,"a");
     for(int nq=0; nq<NUMPFT; nq++)
      BS->PlotScalarDensity(ByPanel[nq], false, PlotFileName,
                            "%s(%s)",PFTNames[nq],z2s(Omega));
     fclose(f);
     
     free(ByPanel[0]);
     free(ByPanel);
     delete BS;
   };
    
}

/***************************************************************/
/* Compute a matrix of reduced-field six-vectors.              */
/*                                                             */
/* Each column of ehMatrix is a 6-component vector giving      */
/* the reduced fields e and h, at a single point in space, due */
/* single RWG basis function.                                  */
/*                                                             */
/* If SurfaceIndex is >=0, then only the contributions of      */
/* RWG functions on surface #SurfaceIndex are taken into       */
/* account. Otherwise, all RWG functions on all surfaces in    */
/* the geometry are taken into account.                        */
/*                                                             */
/* The cartesian coordinates of the evaluation points are the  */
/* first three columns of XMatrix.                             */
/*                                                             */
/* The columns of ehMatrix are indexed by pairs                */
/*                                                             */
/*  (evaluation point index, basis function index)             */
/*                                                             */
/* More specifically, if nx is the index of an evaluation pt   */
/* and ne is the index of an internal edge in an RWGSurface,   */
/* then the components of the reduced fields at point #nx due  */
/* to basis function #ne populated with unit strength are in   */
/*                                                             */
/*  e_[0,1,2] = ehMatrix[ 0:2, nx*NE + ne ]                    */
/*  h_[0,1,2] = ehMatrix[ 3:5, nx*NE + ne ]                    */
/*                                                             */
/* Here NE is the total number of internal edges (in the       */
/* surface, if SurfaceIndex>=0, or in the entire geometry      */
/* if SurfaceIndex==-1).                                       */
/***************************************************************/
HMatrix *Get_ehMatrix(RWGGeometry *G, int SurfaceIndex,
                      HMatrix *XMatrix, cdouble Omega,
                      bool FarField=false, HMatrix *ehMatrix=0)
{ 
  int NX = XMatrix->NR;
  int NET; // total number of edges
  int NumSurfaces;
  if (SurfaceIndex>=0)
   { NET = G->Surfaces[SurfaceIndex]->NumEdges;
     NumSurfaces=1;
   }
  else
   { 
     NET = G->TotalEdges;
     NumSurfaces=G->NumSurfaces;
   };

  /***************************************************************/
  /* reallocate the ehMatrix as necessary ************************/
  /***************************************************************/
  int NETNX = NET*NX;
  if ( ehMatrix && ( (ehMatrix->NR!=6) || ehMatrix->NC!=NETNX) )
   { Warn("wrong-size ehMatrix passed to Get_ehMatrix (reallocating)");
     delete ehMatrix;
     ehMatrix=0;
   };
  if (ehMatrix==0)
   ehMatrix = new HMatrix(6, NETNX, LHM_COMPLEX);
  
  /***************************************************************/
  /* loop over all evaluation points and all edges on all        */
  /* relevant surfaces                                           */
  /***************************************************************/
  Log("Precomputing ehMatrix (%i columns)",NETNX);
#ifdef USE_OPENMP
  int NumThreads=GetNumThreads();
  if (G->LogLevel>=SCUFF_VERBOSE2)
   Log("OpenMP multithreading (%i threads)",NumThreads);
#pragma omp parallel for collapse(3), schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nx=0; nx<NX; nx++)
   for(int ns=0; ns<NumSurfaces; ns++)
    for(int ne=0; ne<NET; ne++)
     { 
       RWGSurface *S;
       int Offset; 
       if (SurfaceIndex==-1) 
        { S=G->Surfaces[ns];
          Offset=G->EdgeIndexOffset[ns];
        }
       else
        { S=G->Surfaces[SurfaceIndex];
          Offset=0;
        };
       int NE = S->NumEdges;
       if (ne>=NE) continue;

       if (ns==0 && ne==0) LogPercent(nx,NX,10);

       // coordinates of nxth evaluation point
       double X[3];
       XMatrix->GetEntriesD(nx,"0:2",X);

       // figure out what region this point is in and 
       // whether or not the surface in question bounds that
       // region
       int RegionIndex=G->GetRegionIndex(X);
       double Sign;
       if (S->RegionIndices[0]==RegionIndex)
        Sign=1.0;
       else if (S->RegionIndices[1]==RegionIndex)
        Sign=-1.0;
       else  
        continue; // surface does not contribute to fields at X

       // get reduced fields due to this edge at this eval pnt
       cdouble EpsRel, MuRel;
       G->RegionMPs[RegionIndex]->GetEpsMu(Omega, &EpsRel, &MuRel);
       cdouble k = Omega*sqrt(EpsRel*MuRel);
       cdouble e[3], h[3];
       if (FarField)
        GetReducedFarFields(S, ne, X, k, e, h);
       else
        GetReducedFields(S, ne, X, k, e, h);

       int nc = nx*NET + Offset + ne;
       for(int Mu=0; Mu<3; Mu++)
        { ehMatrix->SetEntry( 0+Mu, nc, Sign*e[Mu]);
          ehMatrix->SetEntry( 3+Mu, nc, Sign*h[Mu]);
        };

   }; //multithreaded loop

  return ehMatrix;
}

/***************************************************************/
/* This routine gets the contributions of a single pair of     */
/* RWG basis functions to the average Poynting vector and      */
/* Maxwell stress tensor at a given point x.                   */
/***************************************************************/
void GetEdgeEdgeSRFlux(RWGGeometry *G, int nsa, int nea, int nsb, int neb,
                       HMatrix *ehMatrix, int nx,
                       cdouble KK, cdouble KN, cdouble NK, cdouble NN,
                       cdouble Omega, cdouble EpsRel, cdouble MuRel,
                       double SRFlux[NUMSRFLUX])
{
  /***************************************************************/
  /* fetch reduced fields for this cubature point                */
  /***************************************************************/
  int Offseta     = G->EdgeIndexOffset[nsa];
  int Offsetb     = G->EdgeIndexOffset[nsb];
  int NET         = G->TotalEdges;
  int ncAlpha     = nx*NET+Offseta+nea;
  cdouble *eAlpha = ehMatrix->ZM + 6*ncAlpha;
  cdouble *hAlpha = eAlpha+3;
  int ncBeta      = nx*NET+Offsetb+neb;
  cdouble *eBeta  = ehMatrix->ZM + 6*ncBeta;
  cdouble *hBeta  = eBeta+3;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double EpsEE[3][3], MuHH[3][3], EH[3][3], HE[3][3];
  cdouble k    = sqrt(EpsRel*MuRel)*Omega;
  cdouble ZRel = sqrt(MuRel/EpsRel);
  cdouble  IKZ = II*k*ZVAC*ZRel;
  cdouble IKOZ = II*k/(ZVAC*ZRel);
  double EpsAbs = TENTHIRDS * real(EpsRel) / ZVAC;
  double  MuAbs = TENTHIRDS * real(MuRel) * ZVAC;
  for(int m=0; m<3; m++)
   for(int n=0; n<3; n++)
    { 
      EpsEE[m][n] = EpsAbs * real( KK*conj(IKZ)*IKZ*conj(eAlpha[m])*eBeta[n]
                                  -KN*conj(IKZ)    *conj(eAlpha[m])*hBeta[n]
                                  -NK          *IKZ*conj(hAlpha[m])*eBeta[n]
                                  +NN              *conj(hAlpha[m])*hBeta[n]
                                 );
     
      MuHH[m][n]  = MuAbs * real(  KK                *conj(hAlpha[m])*hBeta[n]
                                  +KN           *IKOZ*conj(hAlpha[m])*eBeta[n]
                                  +NK*conj(IKOZ)     *conj(eAlpha[m])*hBeta[n]
                                  +NN*conj(IKOZ)*IKOZ*conj(eAlpha[m])*eBeta[n]
                                );
 
      EH[m][n] = real( KK*conj(IKZ)     *conj(eAlpha[m])*hBeta[n]
                      +KN*conj(IKZ)*IKOZ*conj(eAlpha[m])*eBeta[n]
                      -NK               *conj(hAlpha[m])*hBeta[n]
                      -NN          *IKOZ*conj(hAlpha[m])*eBeta[n]
                     );

      HE[m][n] = real( KK           *IKZ*conj(hAlpha[m])*eBeta[n]
                      -KN               *conj(hAlpha[m])*hBeta[n]
                      +NK*conj(IKOZ)*IKZ*conj(eAlpha[m])*eBeta[n]
                      -NN*conj(IKOZ)    *conj(eAlpha[m])*hBeta[n]
                     );
    };
  double Trace  =   EpsEE[0][0] + EpsEE[1][1] + EpsEE[2][2]
                  +  MuHH[0][0] +  MuHH[1][1] +  MuHH[2][2];
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  // P (poynting flux)
  SRFlux[0] = 0.25*(EH[1][2] - EH[2][1] - HE[1][2] + HE[2][1]);
  SRFlux[1] = 0.25*(EH[2][0] - EH[0][2] - HE[2][0] + HE[0][2]);
  SRFlux[2] = 0.25*(EH[0][1] - EH[1][0] - HE[0][1] + HE[1][0]);

  // T (maxwell stress tensor)
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { 
      double TMuNu = 0.25*(  EpsEE[Mu][Nu] + MuHH[Mu][Nu]
                           + EpsEE[Nu][Mu] + MuHH[Nu][Mu]
                          );
      if (Mu==Nu) TMuNu -= 0.25*Trace;
      SRFlux[ 3 + 3*Mu + Nu ] = TMuNu;
    };
}

/***************************************************************/
/* Evaluate trace formulas for the spatially-resolved fluxes   */
/* at individual points in space.                              */
/*                                                             */
/* XMatrix is an NXx3 matrix storing the cartesian coordinates */
/* of the evaluation points.                                   */
/*                                                             */
/* On return, FMatrix is an NXx12 matrix whose columns are the */
/* components of the average Poynting vector (PV) and Maxwell  */
/* stress tensor (MST) at each evaluation point.               */
/*                                                             */
/* FMatrix[nx, 0..2]  = PV_{x,y,z};                            */
/* FMatrix[nx, 3..11] = MST_{xx}, MST_{xy}, ..., MST_{zz}      */
/***************************************************************/

// compute a unique index for the contribution of thread #nt to the
// SRFlux quantity #nq at spatial point #nx
int GetSRFluxIndex(int nt, int nx, int nq, int NX)
{ return nt*NX*NUMSRFLUX + nx*NUMSRFLUX + nq; }

HMatrix *GetSRFlux(RWGGeometry *G, HMatrix *XMatrix, cdouble Omega,
                   HVector *KNVector, HMatrix *RytovMatrix,
                   HMatrix *FMatrix, bool FarField)
{ 
  /***************************************************************/
  /* (re)allocate FMatrix as necessary ***************************/
  /***************************************************************/
  int NX = XMatrix->NR;
  if ( FMatrix && ( (FMatrix->NR != NX) || (FMatrix->NC != NUMSRFLUX) ) )
   { Warn("Wrong-size FMatrix in GetSRFluxTrace (reallocating...)");
     delete FMatrix;
     FMatrix=0;
   };
  if (FMatrix==0)
   FMatrix = new HMatrix(NX, NUMSRFLUX, LHM_REAL);

  Log("Computing spatially-resolved fluxes at %i evaluation points...",NX);

  /***************************************************************/
  /* precompute 'reduced fields' e,h of all basis functions at   */
  /* all evaluation points                                       */
  /***************************************************************/
  HMatrix *ehMatrix = Get_ehMatrix(G, -1, XMatrix, Omega, FarField);

  /***************************************************************/
  /* prefetch material properties for all regions at this freq   */
  /***************************************************************/
  int NR=G->NumRegions;
  cdouble *RegionEpsMu = new cdouble[ 2*NR ];
  for(int nr=0; nr<G->NumRegions; nr++)
   { cdouble EpsRel, MuRel;
     G->RegionMPs[nr]->GetEpsMu(Omega, &EpsRel, &MuRel);
     RegionEpsMu[2*nr + 0] = EpsRel;
     RegionEpsMu[2*nr + 1] =  MuRel;
   };

  /***************************************************************/
  /* allocate per-thread storage to avoid costly synchronization */
  /* primitives in the multithreaded loop                        */
  /* [note the array is automatically zeroed by mallocEC()]      */
  /***************************************************************/
  int NumThreads=1;
#ifdef USE_OPENMP
  NumThreads=GetNumThreads();
#endif
  double *DeltaSRFlux=(double *)mallocEC(NumThreads*NX*NUMSRFLUX*sizeof(double));
   
  bool UseSymmetry=true;

  /***************************************************************/
  /*- loop over all basis functions and all pairs of eval points */
  /***************************************************************/
  int NET = G->TotalEdges;
  Log("Assembling SR flux...");
#ifdef USE_OPENMP  
  if (G->LogLevel>=SCUFF_VERBOSE2) 
   Log(" using %i OpenMP threads",NumThreads);
#pragma omp parallel for collapse(3), schedule(dynamic,1), num_threads(NumThreads)
#endif  
  for(int nx=0; nx<NX; nx++)
   for(int neta=0; neta<NET; neta++)
    for(int netb=0; netb<NET; netb++)
     { 
       if (UseSymmetry && netb<neta)
        continue;

       if (neta==0 && netb==0)
        LogPercent(nx,NX,100);

       // break up total edge indices into (surfaceIndex, edgeIndex) pairs
       int nsa=0, nea=neta;
       while(nea > G->Surfaces[nsa]->NumEdges)
        nea -= G->Surfaces[nsa++]->NumEdges;
       int nsb=0, neb=netb;
       while(neb > G->Surfaces[nsb]->NumEdges)
        neb -= G->Surfaces[nsb++]->NumEdges;

       RWGSurface *Sa = G->Surfaces[nsa];
       RWGSurface *Sb = G->Surfaces[nsb];
  
       // identify the region in which the eval point lies
       // and skip if the given surfaces do not contribute
       // to the fields in that region
       double X[3];
       XMatrix->GetEntriesD(nx, "0:2", X);
       int RegionIndex = G->GetRegionIndex(X);
       if (    (Sa->RegionIndices[0] != RegionIndex )
            && (Sa->RegionIndices[1] != RegionIndex )
          ) continue;
       if (    (Sb->RegionIndices[0] != RegionIndex )
            && (Sb->RegionIndices[1] != RegionIndex )
          ) continue;

       cdouble EpsRel=RegionEpsMu[2*RegionIndex+0];
       cdouble  MuRel=RegionEpsMu[2*RegionIndex+1];

       /*--------------------------------------------------------------*/
       /* extract the surface-current coefficients either from the KN  */
       /* vector or the Sigma matrix                                   */
       /*--------------------------------------------------------------*/
       bool IsPECa = Sa->IsPEC;
       bool IsPECb = Sb->IsPEC;
       int nbfa = G->BFIndexOffset[nsa] + ( (IsPECa) ? nea : 2*nea );
       int nbfb = G->BFIndexOffset[nsb] + ( (IsPECb) ? neb : 2*neb );
       cdouble KK, KN, NK, NN;
       if (KNVector)
        { 
          cdouble kAlpha =       KNVector->GetEntry(nbfa);
          cdouble kBeta  =       KNVector->GetEntry(nbfb);
          cdouble nAlpha, nBeta = 0.0;
          if (!IsPECa)
           nAlpha = -ZVAC*KNVector->GetEntry(nbfa+1);
          if (!IsPECb)
           nBeta  = -ZVAC*KNVector->GetEntry(nbfb+1);
  
          KK = conj(kAlpha) * kBeta;
          KN = conj(kAlpha) * nBeta;
          NK = conj(nAlpha) * kBeta;
          NN = conj(nAlpha) * nBeta;
        }
       else
        { KK = RytovMatrix->GetEntry(nbfb+0, nbfa+0);
          KN = RytovMatrix->GetEntry(nbfb+1, nbfa+0);
          NK = RytovMatrix->GetEntry(nbfb+0, nbfa+1);
          NN = RytovMatrix->GetEntry(nbfb+1, nbfa+1);
        };

       /*--------------------------------------------------------------*/
       /*- get the contributions of this edge pair --------------------*/
       /*--------------------------------------------------------------*/
       double SRFlux[NUMSRFLUX];
       GetEdgeEdgeSRFlux(G, nsa, nea, nsb, neb, ehMatrix, nx,
                         KK, KN, NK, NN, Omega, EpsRel, MuRel, SRFlux);

       /*--------------------------------------------------------------*/
       /*- accumulate the contributions of this edge pair              */
       /*--------------------------------------------------------------*/
       int nt=0;
#ifdef USE_OPENMP
       nt=omp_get_thread_num();
#endif
       double Weight = (neta==netb) ? 1.0 : 2.0;
       if (!UseSymmetry) Weight=1;

       for(int nq=0; nq<NUMSRFLUX; nq++)
        DeltaSRFlux[ GetSRFluxIndex(nt, nx, nq, NX) ] += Weight*SRFlux[nq];

    }; // for(nx==...)

  /*--------------------------------------------------------------*/
  /*- sum contributions of all threads ---------------------------*/
  /*--------------------------------------------------------------*/
  FMatrix->Zero();
  for(int nx=0; nx<NX; nx++)
   for(int nq=0; nq<NUMSRFLUX; nq++)
    for(int nt=0; nt<NumThreads; nt++)
     FMatrix->AddEntry(nx, nq, DeltaSRFlux[ GetSRFluxIndex(nt, nx, nq, NX)] );

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  free(DeltaSRFlux);
  delete ehMatrix;
  delete[] RegionEpsMu;

  return FMatrix;

} // routine GetSRFlux

} // namespace scuff
