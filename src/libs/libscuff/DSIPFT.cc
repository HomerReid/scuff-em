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
#include "PFTOptions.h"

#ifdef USE_OPENMP
 #include <omp.h>
#endif
 
namespace scuff{

#define II cdouble(0.0,1.0) 

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

// unique index for the contribution of thread #nt to the
// SRFlux quantity #nq at spatial point #nx
inline int GetSRFluxIndex(int NX, int nt, int nx, int nq)
{ return nt*NX*NUMSRFLUX + nx*NUMSRFLUX + nq; }

HMatrix *GetSRFluxTrace(RWGGeometry *G, HMatrix *XMatrix, cdouble Omega,
                        HMatrix *DRMatrix, HMatrix *FMatrix)
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

  /***************************************************************/
  /* FIXME ? *****************************************************/
  /***************************************************************/
  for(int ns=0; ns<G->NumSurfaces; ns++)
   if (G->Surfaces[ns]->IsPEC)
    ErrExit("GetSRFluxTrace not implemented for PEC bodies");
  bool IsPECA=false;
  bool IsPECB=false;

  Log("Computing spatially-resolved fluxes at %i evaluation points...",NX);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NBF = G->TotalBFs;
  static int NBFSave = 0, NXSave=0;
  static HMatrix *RFMatrix=0;
  if (NBFSave!=NBF || NXSave<NX)
   { NBFSave = NBF;
     NXSave  = NX;
     if (RFMatrix) delete RFMatrix;
     RFMatrix = new HMatrix(NBF, 6*NX, LHM_COMPLEX);
   };
  G->GetRFMatrix(Omega, 0, XMatrix, RFMatrix);

  /***************************************************************/
  /* allocate per-thread storage to avoid costly synchronization */
  /* primitives in the multithreaded loop                        */
  /* [note the array is automatically zeroed by mallocEC()]      */
  /***************************************************************/
  int NumThreads=1;
#ifdef USE_OPENMP
  NumThreads=GetNumThreads();
#endif

  size_t DeltaSRFluxSize = NumThreads*NX*NUMSRFLUX*sizeof(cdouble);
  static size_t DeltaSRFluxSizeSave=0;
  static cdouble *DeltaSRFlux = 0;
  if (DeltaSRFluxSizeSave<DeltaSRFluxSize)
   { 
     DeltaSRFluxSizeSave=DeltaSRFluxSize;
     DeltaSRFlux=(cdouble *)reallocEC(DeltaSRFlux,DeltaSRFluxSize);
   };
  memset(DeltaSRFlux, 0, DeltaSRFluxSize);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  G->UpdateCachedEpsMuValues(Omega);
#ifdef USE_OPENMP
  LogC("(%i threads)",NumThreads);
#pragma omp parallel for schedule(dynamic,1),		\
                         num_threads(NumThreads)
#endif
  for(int nx=0; nx<NX; nx++)
   {
     double X[3];
     XMatrix->GetEntriesD(nx,"0:2",X);
     int nr=G->GetRegionIndex(X);
     double  MuAbs = TENTHIRDS*real(G->MuTF[nr] )*ZVAC;
     double EpsAbs = TENTHIRDS*real(G->EpsTF[nr])/ZVAC;
    
     cdouble *EKN[3], *HKN[3];
     for(int Mu=0; Mu<3; Mu++)
      { EKN[Mu] = RFMatrix->ZM + NBF*(6*nx + 3*0 + Mu);
        HKN[Mu] = RFMatrix->ZM + NBF*(6*nx + 3*1 + Mu);
      };

     for(int neaTot=0; neaTot<G->TotalEdges; neaTot++)
      for(int nebTot=0; nebTot<G->TotalEdges; nebTot++)
       { 
         int nsa, nea, nsb, neb, KNIndexA, KNIndexB;
         G->ResolveEdge(neaTot, &nsa, &nea, &KNIndexA);
         G->ResolveEdge(nebTot, &nsb, &neb, &KNIndexB);

         cdouble Bilinears[4];
         GetKNBilinears(0, DRMatrix,
                        IsPECA, KNIndexA, IsPECB, KNIndexB, Bilinears);
         cdouble KK=Bilinears[0];
         cdouble KN=Bilinears[1]/(-1.0*ZVAC);
         cdouble NK=Bilinears[2]/(-1.0*ZVAC);
         cdouble NN=Bilinears[3]/(ZVAC*ZVAC);
 
         cdouble EE[3][3], EH[3][3], HH[3][3];
         for(int Mu=0; Mu<3; Mu++)
          for(int Nu=0; Nu<3; Nu++)
           { EE[Mu][Nu] =  KK*conj(EKN[Mu][KNIndexA+0])*EKN[Nu][KNIndexB+0]
                          +KN*conj(EKN[Mu][KNIndexA+0])*EKN[Nu][KNIndexB+1]
                          +NK*conj(EKN[Mu][KNIndexA+1])*EKN[Nu][KNIndexB+0]
                          +NN*conj(EKN[Mu][KNIndexA+1])*EKN[Nu][KNIndexB+1];

             EH[Mu][Nu] =  KK*conj(EKN[Mu][KNIndexA+0])*HKN[Nu][KNIndexB+0]
                          +KN*conj(EKN[Mu][KNIndexA+0])*HKN[Nu][KNIndexB+1]
                          +NK*conj(EKN[Mu][KNIndexA+1])*HKN[Nu][KNIndexB+0]
                          +NN*conj(EKN[Mu][KNIndexA+1])*HKN[Nu][KNIndexB+1];

             HH[Mu][Nu] =  KK*conj(HKN[Mu][KNIndexA+0])*HKN[Nu][KNIndexB+0]
                          +KN*conj(HKN[Mu][KNIndexA+0])*HKN[Nu][KNIndexB+1]
                          +NK*conj(HKN[Mu][KNIndexA+1])*HKN[Nu][KNIndexB+0]
                          +NN*conj(HKN[Mu][KNIndexA+1])*HKN[Nu][KNIndexB+1];
           };

         cdouble Trace, PV[3], MST[3][3];
         Trace = EpsAbs*(EE[0][0] + EE[1][1] + EE[2][2])
                 +MuAbs*(HH[0][0] + HH[1][1] + HH[2][2]);

         PV[0] = 0.5*( EH[1][2] - EH[2][1] );
         PV[1] = 0.5*( EH[2][0] - EH[0][2] );
         PV[2] = 0.5*( EH[0][1] - EH[1][0] );

         for(int Mu=0; Mu<3; Mu++)
          for(int Nu=0; Nu<3; Nu++)
           MST[Mu][Nu] = 0.5*(EpsAbs*EE[Mu][Nu] + MuAbs*HH[Mu][Nu]);
         MST[0][0] -= 0.25*Trace;
         MST[1][1] -= 0.25*Trace;
         MST[2][2] -= 0.25*Trace;

         int nt=0;
#ifdef USE_OPENMP
         nt = omp_get_thread_num();
#endif
         int Index=GetSRFluxIndex(NX, nt, nx, 0);
         for(int Mu=0; Mu<3; Mu++)
          DeltaSRFlux[Index++] += PV[Mu];
         for(int Mu=0; Mu<3; Mu++)
          for(int Nu=0; Nu<3; Nu++)
           DeltaSRFlux[Index++] += MST[Mu][Nu];

       }; // for(int nbfA=0 ... for(int nbfB=0...

   }; //for(int nx=0; nx<NX; nx++)
         
  /*--------------------------------------------------------------*/
  /*- sum contributions of all threads ---------------------------*/
  /*--------------------------------------------------------------*/
  FMatrix->Zero();
  for(int nx=0; nx<NX; nx++)
   for(int nq=0; nq<NUMSRFLUX; nq++)
    for(int nt=0; nt<NumThreads; nt++)
     FMatrix->AddEntry(nx, nq, real(DeltaSRFlux[ GetSRFluxIndex(NX, nt, nx, nq)] ));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  return FMatrix;

} // routine GetSRFlux

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
/* Otherwise, the cubature rule describes a Lebedev cubature   */
/* rule with DSIPoints cubature points on a sphere of radius R.*/
/*                                                             */
/* If GT1 and/or GT2 are non-null, then each cubature point    */
/* normal vector is transformed by GT1 (first) and then GT2.   */
/***************************************************************/
HMatrix *GetSCRMatrix(char *DSIMesh, double DSIRadius, int DSIPoints,
                      GTransformation *GT1=0, GTransformation *GT2=0)
{
  HMatrix *SCRMatrix;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (DSIMesh)
   { 
     RWGSurface *BS=new RWGSurface(DSIMesh);
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
  else
   { double *LRule = GetLebedevRule(DSIPoints);
     if (LRule==0) ErrExit("no Lebedev rule with %i points",DSIPoints);
     SCRMatrix = new HMatrix(DSIPoints, 7);
     double R=DSIRadius, R2=R*R;
     for(int np=0; np<DSIPoints; np++)
      { SCRMatrix->SetEntry(np,0, R*LRule[4*np + 0]);
        SCRMatrix->SetEntry(np,1, R*LRule[4*np + 1]);
        SCRMatrix->SetEntry(np,2, R*LRule[4*np + 2]);
        SCRMatrix->SetEntry(np,3, LRule[4*np + 0]);
        SCRMatrix->SetEntry(np,4, LRule[4*np + 1]);
        SCRMatrix->SetEntry(np,5, LRule[4*np + 2]);
        SCRMatrix->SetEntry(np,6, R2*LRule[4*np + 3]);
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  GTransformation *GTList[2];
  GTList[0]=GT1;
  GTList[1]=GT2;
  for(int nGT=0; nGT<2; nGT++)
   { 
     GTransformation *GT = GTList[nGT];
     if (GT==0) continue;
     
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
  if ( !NeedQuantity || NeedQuantity[PFT_PABS] )
   { NMatrix[PFT_PABS][0][0] = NMatrix[PFT_PABS][1][1] = NMatrix[PFT_PABS][2][2] = 0.0;
     NMatrix[PFT_PABS][1][2] = nHat[0]; NMatrix[PFT_PABS][2][1] = -nHat[0];
     NMatrix[PFT_PABS][2][0] = nHat[1]; NMatrix[PFT_PABS][0][2] = -nHat[1];
     NMatrix[PFT_PABS][0][1] = nHat[2]; NMatrix[PFT_PABS][1][0] = -nHat[2];
   };

  if ( !NeedQuantity || NeedQuantity[PFT_XFORCE] )
   { NMatrix[PFT_XFORCE][0][0]                           =  nHat[0]; 
     NMatrix[PFT_XFORCE][1][1]                           = -nHat[0];
     NMatrix[PFT_XFORCE][2][2]                           = -nHat[0];
     NMatrix[PFT_XFORCE][0][1] = NMatrix[PFT_XFORCE][1][0] =  nHat[1];
     NMatrix[PFT_XFORCE][0][2] = NMatrix[PFT_XFORCE][2][0] =  nHat[2];
     NMatrix[PFT_XFORCE][1][2] = NMatrix[PFT_XFORCE][2][1] =  0.0;
   };

  if ( !NeedQuantity || NeedQuantity[PFT_YFORCE] )
   { NMatrix[PFT_YFORCE][0][0]                           = -nHat[1]; 
     NMatrix[PFT_YFORCE][1][1]                           =  nHat[1];
     NMatrix[PFT_YFORCE][2][2]                           = -nHat[1];
     NMatrix[PFT_YFORCE][0][1] = NMatrix[PFT_YFORCE][1][0] =  nHat[0];
     NMatrix[PFT_YFORCE][1][2] = NMatrix[PFT_YFORCE][2][1] =  nHat[2];
     NMatrix[PFT_YFORCE][0][2] = NMatrix[PFT_YFORCE][2][0] =  0.0;
   };

  if ( !NeedQuantity || NeedQuantity[PFT_ZFORCE] )
   { NMatrix[PFT_ZFORCE][0][0]                           = -nHat[2];
     NMatrix[PFT_ZFORCE][1][1]                           = -nHat[2];
     NMatrix[PFT_ZFORCE][2][2]                           =  nHat[2];
     NMatrix[PFT_ZFORCE][0][2] = NMatrix[PFT_ZFORCE][2][0] =  nHat[0];
     NMatrix[PFT_ZFORCE][1][2] = NMatrix[PFT_ZFORCE][2][1] =  nHat[1];
     NMatrix[PFT_ZFORCE][0][1] = NMatrix[PFT_ZFORCE][1][0] =  0.0;
   };

  if (      NeedQuantity
       && ( NeedQuantity[PFT_XTORQUE]==false )
       && ( NeedQuantity[PFT_YTORQUE]==false )
       && ( NeedQuantity[PFT_ZTORQUE]==false )
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
     if ( NeedQuantity && (NeedQuantity[PFT_XTORQUE+i]==false) ) 
      continue;

     for(int a=0; a<3; a++)
      for(int b=0; b<3; b++)
       { NMatrix[PFT_XTORQUE+i][a][b] = (a==b) ? -DxN[i] : 0.0;
         for(int j=0; j<3; j++)
          NMatrix[PFT_XTORQUE+i][a][b] 
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
               HVector *KN, IncField *IF, double PFT[NUMPFT],
               char *DSIMesh, double DSIRadius, int DSIPoints,
               bool FarField, char *PlotFileName, 
               GTransformation *GT1, GTransformation *GT2)
{
  (void) FarField;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (DSIMesh)
   Log("Computing DSIPFT over bounding surface %s...",DSIMesh);
  else
   Log("Computing DSIPFT: (Radius,# Pts)=(%e,%i)",DSIRadius,DSIPoints);

  /***************************************************************/
  /* get cubature-rule matrix ************************************/
  /***************************************************************/
  HMatrix *SCRMatrix = GetSCRMatrix(DSIMesh, DSIRadius, DSIPoints, GT1, GT2);

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
  if (GT1) GT1->Apply(XTorque);
  if (GT2) GT2->Apply(XTorque);

  /***************************************************************/
  /* get the scattered and total fields at the cubature points   */
  /***************************************************************/
  HMatrix *FMatrixScat;
/*
  if (FarField && kBloch==0)
   FMatrixScat = GetFarFields(G, 0, KN, Omega, SCRMatrix);
  else
   FMatrixScat = G->GetFields(0, KN, Omega, kBloch, SCRMatrix);
*/
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
  if (DSIMesh && PlotFileName)
   { BS=new RWGSurface(DSIMesh);
     if (GT1) BS->Transform(GT1);
     if (GT2) BS->Transform(GT2);
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
     double dP = -0.25 * w * (  HVMVP(E, NMatrix[PFT_PABS], H)
                               -HVMVP(H, NMatrix[PFT_PABS], E)
                             );
     PFT[PFT_PABS] += dP;
     if (ByPanel) ByPanel[PFT_PABS][ nr ] = dP;

     // scattered power
     dP = 0.25 * w * (  HVMVP(ES, NMatrix[PFT_PABS], HS)
                       -HVMVP(HS, NMatrix[PFT_PABS], ES)
                     );
     PFT[PFT_PSCAT] += dP;
     if (ByPanel) ByPanel[PFT_PSCAT][ nr ] = dP;

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
/***************************************************************/
/***************************************************************/
void GetDSIPFTTrace(RWGGeometry *G, cdouble Omega, HMatrix *DRMatrix,
                    double PFT[NUMPFT], bool NeedQuantity[NUMPFT],
                    char *DSIMesh, double DSIRadius, int DSIPoints,
                    bool FarField, char *PlotFileName,
                    GTransformation *GT1, GTransformation *GT2)
{
  (void) FarField;
  (void) PlotFileName;
  (void) NeedQuantity;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Log("Computing DSIPFT: ");
  if (DSIMesh)
   LogC(" (BS mesh %s)",DSIMesh);
  else
   LogC(" (R=%e, NC=%i)",DSIRadius,DSIPoints);

  double XTorque[3] = {0.0, 0.0, 0.0};
  if (GT1) GT1->Transform(XTorque);
  if (GT2) GT2->Transform(XTorque);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *SCRMatrix = GetSCRMatrix(DSIMesh, DSIRadius, DSIPoints, GT1, GT2);
  HMatrix *SRMatrix  = GetSRFluxTrace(G, SCRMatrix, Omega, DRMatrix);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  memset(PFT, 0, NUMPFT*sizeof(double));
  double MeanPV[3]={0.0, 0.0, 0.0};
  double MeanMST[9]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double Area=0.0;
  for(int nx=0; nx<SCRMatrix->NR; nx++)
   {
     double X[3], nHat[3], Weight;
     SCRMatrix->GetEntriesD(nx, "0:2", X);
     SCRMatrix->GetEntriesD(nx, "3:5", nHat);
     SCRMatrix->GetEntriesD(nx, "6",   &Weight);

     double XmXT[3];
     VecSub(X, XTorque, XmXT);
     
     double PV[3], MST[9];
     SRMatrix->GetEntriesD(nx, "0:2",  PV);
     SRMatrix->GetEntriesD(nx, "3:11", MST);

     Area += Weight;
     VecPlusEquals(MeanPV,Weight,PV,3);
     VecPlusEquals(MeanMST,Weight,MST,9);
   
     double F[3]={0.0, 0.0, 0.0};
     for(int Mu=0; Mu<3; Mu++)
      { PFT[PFT_PABS] += Weight*PV[Mu]*nHat[Mu];
        for(int Nu=0; Nu<3; Nu++)
         F[Mu] += Weight*MST[3*Mu+Nu]*nHat[Nu];
      }; 

     for(int Mu=0; Mu<3; Mu++)
      { int MP1 = (Mu+1)%3, MP2=(Mu+2)%3;
        PFT[PFT_XFORCE  + Mu] += F[Mu];
        PFT[PFT_XTORQUE + Mu] += XmXT[MP1]*F[MP2] - XmXT[MP2]*F[MP1];
      };
   };
  VecScale(MeanPV, 1.0/Area, 3);
  VecScale(MeanMST, 1.0/Area, 9);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (DSIMesh && PlotFileName)
   { 
     static const char *PFTNames[NUMPFT]
      ={"PAbs","PScat","FX","FY","FZ","TX","TY","TZ"};

     RWGSurface *BS=new RWGSurface(DSIMesh);
     FILE *f=fopen(PlotFileName,"a");
    for(int SubtractMean=0; SubtractMean<2; SubtractMean++)
     for(int nq=0; nq<NUMPFT; nq++)
      { 
        if (nq==PFT_PSCAT) continue;

        const char *Suffix = (SubtractMean==0 ? "" : "Bar");
         fprintf(f,"View \"%s%s(Omega=%s)\" {\n",PFTNames[nq],Suffix,z2s(Omega));
        for(int np=0; np<BS->NumPanels; np++)
         { 
           RWGPanel *P=BS->Panels[np];
           double *X0   = BS->Panels[np]->Centroid;
           double *nHat = BS->Panels[np]->ZHat;
           double *TV[3]; // triangle vertices
           TV[0]=BS->Vertices + 3*P->VI[0];
           TV[1]=BS->Vertices + 3*P->VI[1];
           TV[2]=BS->Vertices + 3*P->VI[2];

           double XmXT[3];
           VecSub(X0, XTorque, XmXT);
           double PV[3], MST[9];
           SRMatrix->GetEntriesD(np, "0:2",  PV);
           SRMatrix->GetEntriesD(np, "3:11", MST);

           if (SubtractMean)
            { VecPlusEquals(PV,  -1.0,  MeanPV, 3);
              VecPlusEquals(MST, -1.0, MeanMST, 9);
            };

           double PFlux, FFlux[3]; 
           PFlux    = VecDot(PV, nHat);
           FFlux[0] = VecDot(MST + 0*3, nHat);
           FFlux[1] = VecDot(MST + 1*3, nHat);
           FFlux[2] = VecDot(MST + 2*3, nHat);

           double Val;
           if (nq==PFT_PABS)
            Val = PFlux;
           else if (PFT_XFORCE<=nq && nq<=PFT_ZFORCE)
            Val = FFlux[nq-PFT_XFORCE];
           else if (PFT_XTORQUE<=nq && nq<=PFT_ZTORQUE)
            { int Mu  = nq - PFT_XTORQUE;
              int MP1 = (Mu+1)%3; 
              int MP2 = (Mu+2)%3;
              Val = XmXT[MP1]*FFlux[MP2] - XmXT[MP2]*FFlux[MP1];
            };

           fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                      TV[0][0], TV[0][1], TV[0][2],
                      TV[1][0], TV[1][1], TV[1][2],
                      TV[2][0], TV[2][1], TV[2][2],
                      Val,Val,Val);
         }; // for(int np=0; np<BS->NumPanels; np++)
       fprintf(f,"};\n");
       fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=1;\n");
   
      }; // for(int SubtractMean=0; SubtractMean<2; SubtractMean++), for(int nq=0; nq<NUMPFT; nq++)

     delete BS;
     fclose(f);

   }; // if (DSIMesh && PlotFileName)

  delete SCRMatrix;
  delete SRMatrix;

}

} // namespace scuff
