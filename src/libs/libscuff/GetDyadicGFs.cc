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
 * GetDyadicGF.cc  -- compute the scattering parts of the electric and
 *                    magnetic dyadic Green's functions
 *
 * homer reid      -- 5/2012
 *                 -- 20151212 overhaul and major acceleration
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "PanelCubature.h"

#define II cdouble(0.0,1.0)

namespace scuff { 

cdouble GetG(double R[3], cdouble k, cdouble *dG, cdouble *ddG);

RWGSurface *ResolveNE(RWGGeometry *G, int neFull,
                      int *pns, int *pne, int *pKNIndex);

void GetReducedFields_Nearby(RWGSurface *S, const int ne,
                             const double X0[3], const cdouble k,
                             cdouble e[3], cdouble h[3]);

/***************************************************************/
/* G[0..2] = + <b_ai(x) | G_ij(x,x0)> for j=x,y,z              */
/* C[3..5] = + <b_ai(x) | C_ij(x,x0)> for j=x,y,z              */
/***************************************************************/
typedef struct GCX0Data 
 { cdouble k;
   double *X0;
   GBarAccelerator *GBA;
 } GCX0Data;

void GCX0Integrand(double X[3], double b[3],
                   void *UserData, double W, double *Integral)
{
  GCX0Data *Data       = (GCX0Data *)UserData;
  cdouble k            = Data->k;
  double *X0           = Data->X0;
  GBarAccelerator *GBA = Data->GBA;

  double XmX0[3];
  XmX0[0] = X[0] - X0[0];
  XmX0[1] = X[1] - X0[1];
  XmX0[2] = X[2] - X0[2];

  // Gij = \delta_{ij}*G0 + ddGBar[3*i + j]/k2;
  // Cxy = dGBar[z]/(ik) + cyclic permutations
  cdouble G0, dG[3], ddG[9];
  if (GBA)
   G0=GetGBar(XmX0, GBA, dG, ddG);
  else
   G0=GetG(XmX0, k, dG, ddG);

  cdouble k2=k*k, ik=II*k;
  cdouble *GC = (cdouble *)Integral;
  for(int i=0; i<3; i++)
   { GC[i] += W*(G0 + ddG[3*i+i]/k2)*b[i];
     for(int j=i+1; j<3; j++)
      { GC[i] += W*(ddG[3*i+j]/k2)*b[j]; 
        GC[j] += W*(ddG[3*j+i]/k2)*b[i];
      };
   };

  GC[3+0] += W*(b[1]*dG[2] - b[2]*dG[1]) / ik;
  GC[3+1] += W*(b[2]*dG[0] - b[0]*dG[2]) / ik;
  GC[3+2] += W*(b[0]*dG[1] - b[1]*dG[0]) / ik;

}

/***************************************************************/
/* All evaluation points are assumed to lie in the same region.*/
/* VMatrix is a matrix of dimensions NBF x 6*NX.               */
/***************************************************************/
cdouble GetVMatrix(RWGGeometry *G, cdouble Omega, double *kBloch,
                   HMatrix *XMatrix, HMatrix *VMatrix,
                   int ColumnOffset=0)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double rRelOuterThreshold=4.0;
  double rRelInnerThreshold=0.0;
  int LowOrder=7;
  int HighOrder=20;
  char *s1=getenv("SCUFF_RREL_OUTER_THRESHOLD");
  char *s2=getenv("SCUFF_RREL_INNER_THRESHOLD");
  char *s3=getenv("SCUFF_LOWORDER");
  char *s4=getenv("SCUFF_HIGHORDER");
  if (s1) sscanf(s1,"%le",&rRelOuterThreshold);
  if (s2) sscanf(s2,"%le",&rRelInnerThreshold);
  if (s3) sscanf(s3,"%i",&LowOrder);
  if (s4) sscanf(s4,"%i",&HighOrder);
  if (s1||s2||s3||s4)
   Log("({O,I}rRelThreshold | LowOrder | HighOrder)=(%e,%e,%i,%i)",
       rRelOuterThreshold,rRelInnerThreshold,LowOrder,HighOrder);

  /***************************************************************/
  /* get region index from first point                           */
  /***************************************************************/
  double XFirst[3];
  XFirst[0]=XMatrix->GetEntryD(0,ColumnOffset+0);
  XFirst[1]=XMatrix->GetEntryD(0,ColumnOffset+1);
  XFirst[2]=XMatrix->GetEntryD(0,ColumnOffset+2);
  int RegionIndex = G->GetRegionIndex(XFirst);
  if (RegionIndex==-1) // PEC
   { VMatrix->Zero();
     return 0.0;  
   };
  cdouble nn = G->RegionMPs[RegionIndex]->GetRefractiveIndex(Omega);
  cdouble k  = Omega * nn;

  GBarAccelerator *GBA=0;
  if (G->LDim>0)
   GBA=G->CreateRegionGBA(RegionIndex, Omega, kBloch, XMatrix);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NE=G->TotalEdges;
  int NX=XMatrix->NR;
  int NENX=NE*NX;
#ifndef USE_OPENMP
  Log("Computing VMatrix entries at %i points",NX);
#else
  int NumThreads=GetNumThreads();
  Log("Computing VMatrix entries (%i threads) at %i points",NumThreads,NX);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nenx=0; nenx<NENX; nenx++)
   { 
     int nx     = nenx / NE;
     int neFull = nenx % NE;

     int ns, ne, Offset;
     ResolveNE(G, neFull, &ns, &ne, &Offset);
     RWGSurface *S = G->Surfaces[ns];
     RWGEdge *E    = S->Edges[ne];
   
     double Sign=0.0;
     if      (G->Surfaces[ns]->RegionIndices[0]==RegionIndex) 
      Sign=+1.0;
     else if (G->Surfaces[ns]->RegionIndices[1]==RegionIndex)
      Sign=-1.0;
     else 
      continue;

     double X0[3];
     X0[0]=XMatrix->GetEntryD(nx,ColumnOffset+0);
     X0[1]=XMatrix->GetEntryD(nx,ColumnOffset+1);
     X0[2]=XMatrix->GetEntryD(nx,ColumnOffset+2);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     cdouble GC[6];
     GCX0Data MyData, *Data=&MyData;
     Data->k   = k;
     Data->X0  = X0;
     Data->GBA = GBA;
// note: this is not the right metric to be using...
#define IDIM 12
     double rRel = VecDistance(X0, E->Centroid) / E->Radius;

     if (rRel >= rRelOuterThreshold)
      { 
        GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
                       IDIM, LowOrder, (double *)GC);
      }
     else if (rRel>=rRelInnerThreshold)
      { 
        GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
                       IDIM, HighOrder, (double *)GC);
      }
     else
      { 
        GetReducedFields_Nearby(S, ne, X0, k, GC+0, GC+3);
        GC[3] /= (II*k);
        GC[4] /= (II*k);
        GC[5] /= (II*k);

        if (GBA)
         { cdouble GC1[6], GC2[6];
           int Order=4;
           GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
                          IDIM, Order, (double *)GC1);
           Data->GBA = 0;
           GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
                          IDIM, Order, (double *)GC2);
           for(int Mu=0; Mu<6; Mu++) GC[Mu] += (GC1[Mu] - GC2[Mu]);
         };
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     VMatrix->SetEntry(Offset, 6*nx + 0, Sign*GC[0]);
     VMatrix->SetEntry(Offset, 6*nx + 1, Sign*GC[1]);
     VMatrix->SetEntry(Offset, 6*nx + 2, Sign*GC[2]);
     VMatrix->SetEntry(Offset, 6*nx + 3, Sign*GC[3]);
     VMatrix->SetEntry(Offset, 6*nx + 4, Sign*GC[4]);
     VMatrix->SetEntry(Offset, 6*nx + 5, Sign*GC[5]);
     if ( !(G->Surfaces[ns]->IsPEC) )
      { VMatrix->SetEntry(Offset+1, 6*nx + 0, -Sign*GC[3]);
        VMatrix->SetEntry(Offset+1, 6*nx + 1, -Sign*GC[4]);
        VMatrix->SetEntry(Offset+1, 6*nx + 2, -Sign*GC[5]);
        VMatrix->SetEntry(Offset+1, 6*nx + 3,  Sign*GC[0]);
        VMatrix->SetEntry(Offset+1, 6*nx + 4,  Sign*GC[1]);
        VMatrix->SetEntry(Offset+1, 6*nx + 5,  Sign*GC[2]);
      };
   };

  if (GBA) 
   DestroyGBarAccelerator(GBA);

  return k;
}

/***************************************************************/
/* 20151212 new routine for computing dyadic GFs that uses a   */
/* much faster strategy and computes DGFs at many points       */
/* simultaneously.                                             */
/* There are two options:                                      */
/*  (a) XMatrix has dimension NX x 3. In this case the source  */
/*      and destination points are the same, with cartesian    */
/*      components X[nx,0:2].                                  */
/*  (b) XMatrix has dimension NX x 6. In this case the source  */
/*      and destination points can be different, and we have   */
/*       X[nx,0:2] = destination point                         */
/*       X[nx,3:5] = source points                             */
/***************************************************************/
HMatrix *RWGGeometry::GetDyadicGFs(cdouble Omega, double *kBloch,
                                   HMatrix *XMatrix, HMatrix *M,
                                   HMatrix *GMatrix)
{ 
  int NBF = TotalBFs;
  int NX  = XMatrix->NR;
  Log("Getting DGFs at %i eval points...",NX);

  /*--------------------------------------------------------------*/
  /*- allocate storage for VSource, VDest matrices. I keep these -*/
  /*  on hand as statically-allocated buffers on the assumption  -*/
  /*  that the routine will be called many times with the same   -*/
  /*  number of evaluation points, for example in Brillouin-zone -*/
  /*  integrations                                               -*/
  /*--------------------------------------------------------------*/
  static HMatrix *VSource=0, *VDest=0;
  if ( VSource==0 || VSource->NR!=NBF || VSource->NC!=(6*NX) )
   { 
     if (VSource) delete VSource;
     if (VDest)   delete VDest;
     VSource=new HMatrix(NBF, 6*NX, LHM_COMPLEX);
     VDest=new HMatrix(NBF, 6*NX, LHM_COMPLEX);
   };

  /*--------------------------------------------------------------*/
  /*- allocate an output matrix of the right size if necessary   -*/
  /*--------------------------------------------------------------*/
  if (GMatrix==0 || GMatrix->NR!=NX || GMatrix->NC!=18)
   { 
     if (GMatrix) 
      { Warn("wrong-size GMatrix passed to GetDyadicGFs (reallocating)");
        delete GMatrix;
      };
     GMatrix=new HMatrix(NBF, 18, LHM_COMPLEX);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log("Fetching VMatrix...");
  cdouble k=GetVMatrix(this, Omega, kBloch, XMatrix, VDest);
  if (XMatrix->NC>=6)
   GetVMatrix(this, Omega, kBloch, XMatrix, VSource, 3);
  else
   VSource->Copy(VDest);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log(" LUSolving...");
  M->LUSolve(VSource);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log(" Computing VMVPs...");
  for(int nx=0; nx<NX; nx++)
   for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
     { cdouble GE=0.0, GM=0.0;
       for(int nbf=0; nbf<NBF; nbf++)
        { GE+=VDest->GetEntry(nbf, 6*nx+0+i) * VSource->GetEntry(nbf, 6*nx+0+j);
          GM+=VDest->GetEntry(nbf, 6*nx+3+i) * VSource->GetEntry(nbf, 6*nx+3+j);
        };
       GMatrix->SetEntry(nx, 0 + 3*i + j, -II*k*GE);
       GMatrix->SetEntry(nx, 9 + 3*i + j, +II*k*GM);
     };

  return GMatrix;

}

/***************************************************************/
/* wrapper to support legacy interface to GetDyadicGFs         */
/***************************************************************/
void RWGGeometry::GetDyadicGFs(double XEval[3], double XSource[3],
                               cdouble Omega, double *kBloch,
                               HMatrix *M, HVector *KN,
                               cdouble GEScat[3][3], cdouble GMScat[3][3],
                               cdouble GETot[3][3], cdouble GMTot[3][3])
{ 
  (void )KN; // not used, retained for backward compatibility

  double XBuffer[6];
  memcpy(XBuffer+0, XEval,   3*sizeof(double));
  memcpy(XBuffer+3, XSource, 3*sizeof(double));
  HMatrix XMatrix(1,6,LHM_COMPLEX,LHM_NORMAL,XBuffer);

  if ( (LBasis && !kBloch) || (!LBasis && kBloch) )
   ErrExit("%s:%i: incorrect kBloch specification",__FILE__,__LINE__);

  cdouble GBuffer[18];
  HMatrix GMatrix(1,18,LHM_COMPLEX,LHM_NORMAL,GBuffer);
  
  GetDyadicGFs(Omega, kBloch, &XMatrix, M, &GMatrix);

  for(int i=0; i<3; i++)
   for(int j=0; j<3; j++)
    { GEScat[i][j] = GMatrix.GetEntry(0, 0 + 3*i + j );
      GMScat[i][j] = GMatrix.GetEntry(0, 9 + 3*i + j );
    };
  
  /***************************************************************/
  /* add direction contributions of point source                 */
  /***************************************************************/
  cdouble EpsRel, MuRel;
  int nr=GetRegionIndex(XSource);
  RegionMPs[nr]->GetEpsMu(Omega, &EpsRel, &MuRel);
  cdouble k=Omega*sqrt(EpsRel*MuRel);
  cdouble EFactor = k*k/EpsRel;
  cdouble MFactor = k*k/MuRel;

  cdouble P[3]={1.0, 0.0, 0.0};
  PointSource PS(XSource, P);
  if (LDim>0)
   { PS.SetLattice(LBasis);
     PS.SetkBloch(kBloch);
   };

  for(int j=0; j<3; j++)
   { 
     cdouble EH[6];

     // set point source to point in the jth direction
     memset(P, 0, 3*sizeof(cdouble));
     P[j]=1.0;
     PS.SetP(P);

     PS.SetType(LIF_ELECTRIC_DIPOLE);
     PS.GetFields(XEval, EH);
     GETot[0][j] = GEScat[0][j] + EH[0] / EFactor;
     GETot[1][j] = GEScat[1][j] + EH[1] / EFactor;
     GETot[2][j] = GEScat[2][j] + EH[2] / EFactor;

     PS.SetType(LIF_MAGNETIC_DIPOLE);
     PS.GetFields(XEval, EH);
     GMTot[0][j] = GMScat[0][j] + EH[3] / MFactor;
     GMTot[1][j] = GMScat[1][j] + EH[4] / MFactor;
     GMTot[2][j] = GMScat[2][j] + EH[5] / MFactor;
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::GetDyadicGFs(double X[3], cdouble Omega, double *kBloch,
                               HMatrix *M, HVector *KN,
                               cdouble GEScat[3][3], cdouble GMScat[3][3])
{
  cdouble Dummy[3][3];
  GetDyadicGFs(X, X, Omega, kBloch, M, KN, GEScat, GMScat, Dummy, Dummy);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::GetDyadicGFs2(cdouble Omega, double *kBloch,
                                    HMatrix *XMatrix, HMatrix *M,
                                    HMatrix *GMatrix)
{ 
  int NBF = TotalBFs;
  int NX  = XMatrix->NR;
  Log("Getting DGFs at %i eval points...",NX);

  /*--------------------------------------------------------------*/
  /*- allocate storage for VSource, VDest matrices. I keep these -*/
  /*  on hand as statically-allocated buffers on the assumption  -*/
  /*  that the routine will be called many times with the same   -*/
  /*  number of evaluation points, for example in Brillouin-zone -*/
  /*  integrations.                                              -*/
  /*--------------------------------------------------------------*/
  static HMatrix *VSource=0, *VDest=0;
  if ( VSource==0 || VSource->NR!=NBF || VSource->NC!=(6*NX) )
   { 
     if (VSource) delete VSource;
     if (VDest)   delete VDest;
     VSource=new HMatrix(NBF, 6*NX, LHM_COMPLEX);
     VDest=new HMatrix(NBF, 6*NX, LHM_COMPLEX);
   };

  /*--------------------------------------------------------------*/
  /*- allocate an output matrix of the right size if necessary   -*/
  /*--------------------------------------------------------------*/
  if (GMatrix==0 || GMatrix->NR!=NX || GMatrix->NC!=18)
   { 
     if (GMatrix) 
      { Warn("wrong-size GMatrix passed to GetDyadicGFs (reallocating)");
        delete GMatrix;
      };
     GMatrix=new HMatrix(NBF, 18, LHM_COMPLEX);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log("Fetching VMatrix...");
  cdouble k=GetVMatrix(this, Omega, kBloch, XMatrix, VDest);
  if (XMatrix->NC>=6)
   GetVMatrix(this, Omega, kBloch, XMatrix, VSource, 3);
  else
   VSource->Copy(VDest);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log(" LUSolving...");
  M->LUSolve(VSource);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Log(" Computing VMVPs...");
  for(int nx=0; nx<NX; nx++)
   for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
     { cdouble GE=0.0, GM=0.0;
       for(int nbf=0; nbf<NBF; nbf++)
        { GE+=VDest->GetEntry(nbf, 6*nx+0+i) * VSource->GetEntry(nbf, 6*nx+0+j);
          GM+=VDest->GetEntry(nbf, 6*nx+3+i) * VSource->GetEntry(nbf, 6*nx+3+j);
        };
       GMatrix->SetEntry(nx, 0 + 3*i + j, -II*k*GE);
       GMatrix->SetEntry(nx, 9 + 3*i + j, +II*k*GM);
     };

  return GMatrix;

}

} // namespace scuff
