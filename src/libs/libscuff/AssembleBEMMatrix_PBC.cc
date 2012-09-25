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
 * AssembleBEMMatrix_PBC.cc -- the periodic-boundary-condition version of 
 *                          -- the BEM matrix assembly routine 
 *                          --
 * homer reid               -- 8/2011 -- 8/2012  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libhrutil.h>
#include <libTriInt.h>
#include <libMDInterp.h>
#include <libscuff.h>
#include <libscuffInternals.h>

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#define II cdouble(0.0,1.0)

namespace scuff{

/***************************************************************/
/* this is a panel-panel cubature routine similar to           */
/* GetPPIs_Cubature, but with the difference that the usual    */
/* helmholtz GF is replaced by its periodic version (with the  */
/* contributions of the innermost 9 grid cells excluded) as    */
/* computed by the Interpolator object.                        */
/***************************************************************/
void GetAB9PanelPanelInteraction(double **Va, double *Qa,
                                 double **Vb, double *Qb,
                                 cdouble k, Interp3D *Interpolator,
                                 cdouble GC[2])
{ 
  /***************************************************************/
  /* preliminary setup for numerical cubature, similar to code in*/
  /* PanelPanelInteractions.cc.                                  */
  /***************************************************************/
  double *V0, A[3], B[3], *Q;
  V0=Va[0];
  VecSub(Va[1], Va[0], A);
  VecSub(Va[2], Va[0], B);
  Q=Qa;

  double *V0P, AP[3], BP[3], *QP;
  V0P=Vb[0];
  VecSub(Vb[1], Vb[0], AP);
  VecSub(Vb[2], Vb[0], BP);
  QP=Qb;

  /***************************************************************/
  /* TCR ('triangle cubature rule') points to a vector of 3N     */
  /* doubles (for an N-point cubature rule).                     */
  /* TCR[3*n,3*n+1,3*n+2]=(u,v,w), where (u,v)                   */
  /* are the (x,y) coordinates of the nth quadrature point and w */
  /* is its weight.                                              */
  /* note we use the same quadrature rule for both the source    */
  /* and destination triangles.                                  */
  /***************************************************************/
  double *TCR;
  int NumPts;
  TCR=GetTCR(RWGGeometry::PBCCubatureOrder, &NumPts);

  /***************************************************************/
  /* outer loop **************************************************/
  /***************************************************************/
  double hDot, hNabla=4.0; // note hNabla is constant throughout
  cdouble hPlus, ik = II*k, ik2=ik*ik;
  int np, ncp, npp, ncpp, i;
  double u, v, w, up, vp, wp;
  double X[3], F[3], XP[3], FP[3], R[3], FxFP[3];
  int ZFlipped;
  double PhiVD[16];
  cdouble GBar, GradGBar[3];
  cdouble GCInner[2];
  memset(GC,0,2*sizeof(cdouble));
  for(np=ncp=0; np<NumPts; np++) 
   { 
     u=TCR[ncp++]; v=TCR[ncp++]; w=TCR[ncp++];

     /***************************************************************/
     /* set X and F=X-Q *********************************************/
     /***************************************************************/
     for(i=0; i<3; i++)
      { X[i] = V0[i] + u*A[i] + v*B[i];
        F[i] = X[i] - Q[i];
      };

     /***************************************************************/
     /* inner loop to calculate value of inner integrand ************/
     /***************************************************************/
     memset(GCInner,0,2*sizeof(cdouble));
     for(npp=ncpp=0; npp<NumPts; npp++)
      { 
        up=TCR[ncpp++]; vp=TCR[ncpp++]; wp=TCR[ncpp++];

        /***************************************************************/ 
        /* set XP and FP=XP-QP *****************************************/
        /***************************************************************/
        for(i=0; i<3; i++)
         { XP[i] = V0P[i] + up*AP[i] + vp*BP[i];
           FP[i] = XP[i] - QP[i];
           R[i] = X[i] - XP[i];
         };

        /***************************************************************/
        /* use interpolators to get the value of the periodic GF at R  */
        /***************************************************************/
        if ( R[2] < 0.0 )
         { R[2] *= -1.0;
           ZFlipped=1;
         }
        else
         ZFlipped=0;

        Interpolator->EvaluatePlus(R[0], R[1], R[2], PhiVD);

        GBar = cdouble(PhiVD[0],PhiVD[8+0]);
        GradGBar[0] = cdouble(PhiVD[1],PhiVD[8+1]);
        GradGBar[1] = cdouble(PhiVD[2],PhiVD[8+2]);
        GradGBar[2] = cdouble(PhiVD[3],PhiVD[8+3]);

        if (ZFlipped)
         GradGBar[2]*=-1.0; // flip the sign of dG/dz 
      
        /***************************************************************/
        /* inner integrand  ********************************************/
        /***************************************************************/
        /* compute h factors */
        hDot=VecDot(F, FP);
        hPlus=hDot + hNabla/ik2;
        VecCross(F, FP, FxFP);
   
        GCInner[0] += wp * hPlus * GBar;   
        GCInner[1] += wp * (  FxFP[0]*GradGBar[0] 
                             +FxFP[1]*GradGBar[1] 
                             +FxFP[2]*GradGBar[2] );

      }; /* for(npp=ncpp=0; npp<NumPts; npp++) */

     /*--------------------------------------------------------------*/
     /*- accumulate contributions to outer integral                  */
     /*--------------------------------------------------------------*/
     GC[0]+=w*GCInner[0];
     GC[1]+=w*GCInner[1];

   }; // for(np=ncp=0; np<nPts; np++)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetAB9EdgeEdgeInteractions(RWGSurface *Sa, int nea, RWGSurface *Sb, int neb,
                                cdouble k, Interp3D *Interpolator, cdouble *GC)
{

  RWGEdge *Ea = Sa->Edges[nea];
  RWGEdge *Eb = Sb->Edges[neb];

  double *Va[3], *Qa, *Vb[3], *Qb;
  cdouble GCPP[2], GCPM[2], GCMP[2], GCMM[2];

  Va[1] = Sa->Vertices + 3*(Ea->iV1);
  Va[2] = Sa->Vertices + 3*(Ea->iV2);

  Vb[1] = Sb->Vertices + 3*(Eb->iV1);
  Vb[2] = Sb->Vertices + 3*(Eb->iV2);

  /*- PP ---------------------------------------------------------*/ 
  Va[0] = Qa = Sa->Vertices + 3*(Ea->iQP);
  Vb[0] = Qb = Sb->Vertices + 3*(Eb->iQP);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
RWGPanel *Pa = Sa->Panels[Ea->iPPanel];
RWGPanel *Pb = Sb->Panels[Eb->iPPanel];
Va[0] = Sa->Vertices + 3*(Pa->VI[0]);
Va[1] = Sa->Vertices + 3*(Pa->VI[1]);
Va[2] = Sa->Vertices + 3*(Pa->VI[2]);
Vb[0] = Sb->Vertices + 3*(Pb->VI[0]);
Vb[1] = Sb->Vertices + 3*(Pb->VI[1]);
Vb[2] = Sb->Vertices + 3*(Pb->VI[2]);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  GetAB9PanelPanelInteraction(Va, Qa, Vb, Qb, k, Interpolator, GCPP);

  /*- PM ---------------------------------------------------------*/ 
  Va[0] = Qa = Sa->Vertices + 3*(Ea->iQP);
  Vb[0] = Qb = Sb->Vertices + 3*(Eb->iQM);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
Pa = Sa->Panels[Ea->iPPanel];
Pb = Sb->Panels[Eb->iMPanel];
Va[0] = Sa->Vertices + 3*(Pa->VI[0]);
Va[1] = Sa->Vertices + 3*(Pa->VI[1]);
Va[2] = Sa->Vertices + 3*(Pa->VI[2]);
Vb[0] = Sb->Vertices + 3*(Pb->VI[0]);
Vb[1] = Sb->Vertices + 3*(Pb->VI[1]);
Vb[2] = Sb->Vertices + 3*(Pb->VI[2]);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  GetAB9PanelPanelInteraction(Va, Qa, Vb, Qb, k, Interpolator, GCPM);

  /*- MP ---------------------------------------------------------*/ 
  Va[0] = Qa = Sa->Vertices + 3*(Ea->iQM);
  Vb[0] = Qb = Sb->Vertices + 3*(Eb->iQP);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
Pa = Sa->Panels[Ea->iMPanel];
Pb = Sb->Panels[Eb->iPPanel];
Va[0] = Sa->Vertices + 3*(Pa->VI[0]);
Va[1] = Sa->Vertices + 3*(Pa->VI[1]);
Va[2] = Sa->Vertices + 3*(Pa->VI[2]);
Vb[0] = Sb->Vertices + 3*(Pb->VI[0]);
Vb[1] = Sb->Vertices + 3*(Pb->VI[1]);
Vb[2] = Sb->Vertices + 3*(Pb->VI[2]);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  GetAB9PanelPanelInteraction(Va, Qa, Vb, Qb, k, Interpolator, GCMP);

  /*- MM ---------------------------------------------------------*/ 
  Va[0] = Qa = Sa->Vertices + 3*(Ea->iQM);
  Vb[0] = Qb = Sb->Vertices + 3*(Eb->iQM);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
Pa = Sa->Panels[Ea->iMPanel];
Pb = Sb->Panels[Eb->iMPanel];
Va[0] = Sa->Vertices + 3*(Pa->VI[0]);
Va[1] = Sa->Vertices + 3*(Pa->VI[1]);
Va[2] = Sa->Vertices + 3*(Pa->VI[2]);
Vb[0] = Sb->Vertices + 3*(Pb->VI[0]);
Vb[1] = Sb->Vertices + 3*(Pb->VI[1]);
Vb[2] = Sb->Vertices + 3*(Pb->VI[2]);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  GetAB9PanelPanelInteraction(Va, Qa, Vb, Qb, k, Interpolator, GCMM);

  double GPreFac = Ea->Length * Eb->Length;
  cdouble CPreFac = GPreFac / (II*k);
  GC[0] = GPreFac * (GCPP[0] - GCPM[0] - GCMP[0] + GCMM[0]);
  GC[1] = CPreFac * (GCPP[1] - GCPM[1] - GCMP[1] + GCMM[1]);
  
}

/***************************************************************/
/* add the matrix elements corresponding to a given pair of    */
/* RWG basis functions to the BEM matrix                       */
/***************************************************************/
void AddMEs(HMatrix *M, RWGSurface *Sa, int nea, RWGSurface *Sb, int neb, 
            int RowOffset, int ColOffset, cdouble PreFac[3], cdouble GC[2])
{
  int X, Y;
  int Symmetric=0;

  if ( Sa->IsPEC && Sb->IsPEC )
   { 
     X=RowOffset + nea;
     Y=ColOffset + neb;  

     M->AddEntry( X, Y, PreFac[0]*GC[0] );
   }
  else if ( Sa->IsPEC && !Sb->IsPEC )
   { 
     X=RowOffset + nea;
     Y=ColOffset + 2*neb;  

     M->AddEntry( X, Y,   PreFac[0]*GC[0] );
     M->AddEntry( X, Y+1, PreFac[1]*GC[1] );

   }
  else if ( !Sa->IsPEC && Sb->IsPEC )
   {
     X=RowOffset + 2*nea;
     Y=ColOffset + neb;  

     M->AddEntry( X,   Y, PreFac[0]*GC[0] );
     M->AddEntry( X+1, Y, PreFac[1]*GC[1] );
   }
  else // ( !Sa->IsPEC && !Sb->IsPEC )
   { 
     X=RowOffset + 2*nea;
     Y=ColOffset + 2*neb;  

     M->AddEntry( X, Y,   PreFac[0]*GC[0]);
     M->AddEntry( X, Y+1, PreFac[1]*GC[1]);
     if ( !Symmetric || (nea!=neb) )
      M->AddEntry( X+1, Y, PreFac[1]*GC[1]);
     M->AddEntry( X+1, Y+1, PreFac[2]*GC[0]);

   };


}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   RWGGeometry *G;
   HMatrix *M;
   int nt, nTask;
 } ThreadData;

/***************************************************************/
/* thread routine for AddOuterCellContributions ****************/
/***************************************************************/
void *AOCC_Thread(void *data)
{
  /***************************************************************/
  /* extract local copies of fields in argument structure        */
  /***************************************************************/
  ThreadData *TD=(ThreadData *)data;
  RWGGeometry *G       = TD->G;
  HMatrix *M           = TD->M;

  cdouble *EpsTF       = G->EpsTF;
  cdouble *MuTF        = G->MuTF;
  cdouble Omega        = G->StoredOmega;

  /***************************************************************/
  /* other local variables ***************************************/
  /***************************************************************/
  int nt=0;
  RWGSurface *Sa, *Sb;
  int RowOffset, ColOffset;

  double Signs[2];
  int CommonRegions[2];
  int NumCommonRegions;

  int nr1, nr2;
  cdouble k1, k2;
  Interp3D *Interpolator1, *Interpolator2;
  cdouble PreFac1[3], PreFac2[3];
  cdouble GC[2];

  for(int nsa=0; nsa<G->NumSurfaces; nsa++)
   for(int nsb=0; nsb<G->NumSurfaces; nsb++)
    { 
      Sa=G->Surfaces[nsa];
      Sb=G->Surfaces[nsb];

      NumCommonRegions=CountCommonRegions(Sa, Sb, CommonRegions, Signs);
      if (NumCommonRegions==0) 
       continue;

      RowOffset = G->BFIndexOffset[nsa];
      ColOffset = G->BFIndexOffset[nsb];

      nr1=CommonRegions[0];
      Interpolator1=G->GBarAB9Interpolators[nr1];
      k1=csqrt2(EpsTF[nr1]*MuTF[nr1])*Omega;
      PreFac1[0] =  Signs[0]*II*MuTF[nr1]*Omega;
      PreFac1[1] = -Signs[0]*II*k1;
      PreFac1[2] = -Signs[0]*II*EpsTF[nr1]*Omega;

      if ( NumCommonRegions==2 )
       { nr2=CommonRegions[1];
         Interpolator2=G->GBarAB9Interpolators[nr2];
         k2=csqrt2(EpsTF[nr2]*MuTF[nr2])*Omega;
         PreFac2[0] =  Signs[1]*II*MuTF[nr2]*Omega;
         PreFac2[1] = -Signs[1]*II*k2;
         PreFac2[2] = -Signs[1]*II*EpsTF[nr2]*Omega;
       }
      else
       Interpolator2=NULL;

      if (Interpolator1==NULL && Interpolator2==NULL)
       continue;       // this can happen if neither region is extended

      for(int nea=0; nea<Sa->NumEdges; nea++)  
       for(int neb=0; neb<Sb->NumEdges; neb++)
        { 
          nt++;
          if (nt==TD->nTask) nt=0;
          if (nt!=TD->nt) continue;

          if (neb==0) LogPercent(nea, Sa->NumEdges);

          // contribution from first medium if it is extended
          if (Interpolator1)
           { GetAB9EdgeEdgeInteractions(Sa, nea, Sb, neb, k1, Interpolator1, GC);
             AddMEs(M, Sa, nea, Sb, neb, RowOffset, ColOffset, PreFac1, GC);
           };

          // contribution from second medium if it is extended
          if (Interpolator2)
           { GetAB9EdgeEdgeInteractions(Sa, nea, Sb, neb, k2, Interpolator2, GC);
             AddMEs(M, Sa, nea, Sb, neb, RowOffset, ColOffset, PreFac2, GC);
           };
          
        }; // for(nea=0 ... for(neb=0 ... 

    }; // for (nsa=0 ... for(nsb=0 ... 

  return 0;

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::AddOuterCellContributions(double kBloch[MAXLATTICE], HMatrix *M)
{ 

  Log(" Adding contributions of outer cells...");

  /*--------------------------------------------------------------*/
  /*- initialize GBarAB9 interpolation tables for all regions at -*/
  /* the present frequency and bloch wavevector                  -*/
  /*--------------------------------------------------------------*/
  GBarData MyGBarData, *GBD=&MyGBarData;
  GBD->kBloch = kBloch;
  GBD->ExcludeInner9=true;
  GBD->E=-1.0;
  GBD->LBV[0]=LatticeBasisVectors[0];
  GBD->LBV[1]=LatticeBasisVectors[1];
  for(int nr=0; nr<NumRegions; nr++)
   {
     if (GBarAB9Interpolators[nr]==0) 
      continue;
     Log("  Initializing interpolator for region %i (%s)...",nr,RegionLabels[nr]);
     GBD->k = csqrt2(EpsTF[nr]*MuTF[nr])*StoredOmega;
     GBarAB9Interpolators[nr]->ReInitialize(GBarVDPhi3D, (void *)GBD);
   };
  
  Log("  Computing GBarAB9 panel--panel integrals...");

  /*--------------------------------------------------------------*/
  /*- fire off threads -------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nt;
  int nThread=GetNumThreads();
#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[nThread], *TD;
  pthread_t *Threads = new pthread_t[nThread];
  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDs[nt]);
     TD->G = this;
     TD->M = M;
     TD->nt=nt;
     TD->nTask=nThread;
     if (nt+1 == nThread)
       AOCC_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, AOCC_Thread, (void *)TD);
   }
  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);
  delete[] Threads;
  delete[] TDs;

#else 
  int nTask;
#ifndef USE_OPENMP
  nThread=nTask=1;
#else
  nTask=nThread*100;
#pragma omp parallel for schedule(dynamic,1), num_threads(nThread)
#endif
  for(nt=0; nt<nTask; nt++)
   { 
     ThreadData TD1;
     TD1.G = this;
     TD1.M = M;
     TD1.nt=nt;
     TD1.nTask=nTask;
     AOCC_Thread((void *)&TD1);
   };
#endif
  
} 

/***************************************************************/
/* This routine handles the contributions of the innermost     */
/* lattice cells to the BEM matrix.                            */
/*                                                             */
/* If UsePBCAcceleration is true, then the contributions are   */
/* simply stored in the MPP, MPM, MPZ, MZP, and MZZ blocks     */
/* inside the RWGGeometry class and are not stamped into the   */
/* actual BEM matrix. In this case the M input parameter       */
/* should be set to NULL.                                      */
/*                                                             */
/* If UsePBCAcceleration is false, then the contributions are  */
/* stamped into the BEM matrix (stored in the M parameter)     */ 
/* with the proper phase factors depending on the bloch        */
/* wavevector, as they are assembled.                          */
/*                                                             */
/* Here Mij = the interaction of the unit cell geometry with   */
/*            the periodic copy of that geometry at i*L1 + j*L2*/
/*            where L1, L2 are the two lattice vectors and     */
/*            where i,j = P, M, Z for +1, -1, 0                */
/*                                                             */
/* NOTE: There is one tricky thing here. The way we compute    */
/*       the portion of the Mij block that accounts for the    */
/*       interactions between two surfaces Sa and Sb is simply */
/*       to displace the entire surface Sb through a           */
/*       translation vector i*L1 + j*L2. Everything else about */
/*       the surfaces is unchanged, including the information  */
/*       about which regions they bound. The difficulty is     */
/*       that if Sa and Sb interact through a medium in the    */
/*       unit cell that is compact, i.e. not extended (imagine */
/*       Sa and Sb are the upper and lower hemispherical       */
/*       surfaces of a sphere that is fully contained within   */
/*       the unit cell), then the contribution to MPP from     */
/*       these two surfaces comes only from the interaction    */
/*       through the exterior medium. However, the             */
/*       AssembleBEMMatrixBlock() routine doesn't know this,   */
/*       and will erroneously compute a contribution coming    */
/*       from the region interior to the sphere as well.       */
/*       To avoid this spurious contribution, we detect such   */
/*       cases and zero out the material properties of the     */
/*       compact region (in the case described here it would   */
/*       be the region inside the sphere), then restore the    */
/*       material properties after the call to                 */
/*       AssembleBEMMatrixBlock. Kinda klugey, but it works.   */
/***************************************************************/
void RWGGeometry::AssembleInnerCellBlocks(double *kBloch, HMatrix *M)
{
  /*- note: PFPP = 'phase factor, plus-plus'                      */
  /*-       PFPZ = 'phase factor, plus-zero'                      */
  /*- etc.                                                        */
  cdouble PFPP, PFPM, PFPZ, PFZP;
  if (M)
   { 
     double *LBV[2];
     LBV[0]=LatticeBasisVectors[0];
     LBV[1]=LatticeBasisVectors[1];
     PFPP = exp( II* (kBloch[0]*(LBV[0][0] + LBV[1][0]) + kBloch[1]*(LBV[0][1] + LBV[1][1]) ) );
     PFPM = exp( II* (kBloch[0]*(LBV[0][0] - LBV[1][0]) + kBloch[1]*(LBV[0][1] - LBV[1][1]) ) );
     PFPZ = exp( II* (kBloch[0]*(LBV[0][0]            ) + kBloch[1]*(LBV[0][1]            ) ) );
     PFZP = exp( II* (kBloch[0]*(            LBV[1][0]) + kBloch[1]*(            LBV[1][1]) ) );
   };

  Log(" Assembling inner matrix blocks at Omega=%s",z2s(StoredOmega));

  Log(" MZZ block...");
  AssembleBEMMatrix(StoredOmega, M ? M : MZZ );

  // if we are not using PBC acceleration, then we use MZZ as a 
  // scratch matrix for assembling all inner cell blocks.
  HMatrix *MScratch = MZZ;

  int nr1, nr2=0;
  int NumCommonRegions, CommonRegionIndices[2];
  double Signs[2];
  bool Region1Contributes, Region2Contributes;
  int SaveZeroed1, SaveZeroed2=0;

  int ns, nsp;
  ABMBArgStruct MyABMBArgStruct, *Args=&MyABMBArgStruct;
  InitABMBArgs(Args);
  Args->G = this;
  Args->Omega = StoredOmega;
  Args->Symmetric=0;

  double Displacement[3];
  Displacement[2]=0.0;
  Args->Displacement=Displacement;

  Log(" MPP block ...");
  Displacement[0]=LatticeBasisVectors[0][0] + LatticeBasisVectors[1][0];
  Displacement[1]=LatticeBasisVectors[0][1] + LatticeBasisVectors[1][1];
  Args->B = M ? MScratch : MPP;
  for(ns=0; ns<NumSurfaces; ns++)
   for(nsp=0; nsp<NumSurfaces; nsp++)
    { 
      Args->Sa=Surfaces[ns];
      Args->Sb=Surfaces[nsp];
      Args->Symmetric=0;
      Args->RowOffset=BFIndexOffset[ns];
      Args->ColOffset=BFIndexOffset[nsp];

      NumCommonRegions=CountCommonRegions(Args->Sa, Args->Sb, CommonRegionIndices, Signs);
      if (NumCommonRegions==0) 
       continue;

      // see NOTE above
      nr1=CommonRegionIndices[0]; 
      Region1Contributes = (RegionIsExtended[MAXLATTICE*nr1+0] && RegionIsExtended[MAXLATTICE*nr1+1]);
      if (NumCommonRegions==2)
       { nr2=CommonRegionIndices[1];
         Region2Contributes = (RegionIsExtended[MAXLATTICE*nr2+0] && RegionIsExtended[MAXLATTICE*nr2+1]);
       }
      else
       Region2Contributes=false;

      if ( Region1Contributes==false  && Region2Contributes==false )
       continue;

      if ( Region1Contributes==false )
       { SaveZeroed1=RegionMPs[nr1]->Zeroed; 
         RegionMPs[nr1]->Zero();
       };
      if ( NumCommonRegions==2 && Region2Contributes==false )
       { SaveZeroed2=RegionMPs[nr2]->Zeroed; 
         RegionMPs[nr2]->Zero();
       };

      AssembleBEMMatrixBlock(Args);

      if ( Region1Contributes==false )
       RegionMPs[nr1]->Zeroed=SaveZeroed1;
      if ( NumCommonRegions==2 && Region2Contributes==false )
       RegionMPs[nr2]->Zeroed=SaveZeroed2;

    };
  if (M)
   { for(int nr=0; nr<M->NR; nr++)
      for(int nc=0; nc<M->NR; nc++)
       M->AddEntry(nr, nc,       PFPP* MScratch->GetEntry(nr, nc) +
                            conj(PFPP)*MScratch->GetEntry(nc, nr) );
   };

  Log(" MPM block ...");
  Displacement[0]=LatticeBasisVectors[0][0] - LatticeBasisVectors[1][0];
  Displacement[1]=LatticeBasisVectors[0][1] - LatticeBasisVectors[1][1];
  Args->B = M ? MScratch : MPM;
  for(ns=0; ns<NumSurfaces; ns++)
   for(nsp=0; nsp<NumSurfaces; nsp++)
    { Args->Sa=Surfaces[ns];
      Args->RowOffset=BFIndexOffset[ns];
      Args->Sb=Surfaces[nsp];
      Args->ColOffset=BFIndexOffset[nsp];

      NumCommonRegions=CountCommonRegions(Args->Sa, Args->Sb, CommonRegionIndices, Signs);
      if (NumCommonRegions==0) 
       continue;

      // see NOTE above
      nr1=CommonRegionIndices[0]; 
      Region1Contributes = (RegionIsExtended[MAXLATTICE*nr1+0] && RegionIsExtended[MAXLATTICE*nr1+1]);
      if (NumCommonRegions==2)
       { nr2=CommonRegionIndices[1];
         Region2Contributes = (RegionIsExtended[MAXLATTICE*nr2+0] && RegionIsExtended[MAXLATTICE*nr2+1]);
       }
      else
       Region2Contributes=false;

      if ( Region1Contributes==false  && Region2Contributes==false )
       continue;

      if ( Region1Contributes==false )
       { SaveZeroed1=RegionMPs[nr1]->Zeroed; 
         RegionMPs[nr1]->Zero();
       };
      if ( NumCommonRegions==2 && Region2Contributes==false )
       { SaveZeroed2=RegionMPs[nr2]->Zeroed; 
         RegionMPs[nr2]->Zero();
       };

      AssembleBEMMatrixBlock(Args);

      if ( Region1Contributes==false )
       RegionMPs[nr1]->Zeroed=SaveZeroed1;
      if ( NumCommonRegions==2 && Region2Contributes==false )
       RegionMPs[nr2]->Zeroed=SaveZeroed2;

    };
  if (M)
   { for(int nr=0; nr<M->NR; nr++)
      for(int nc=0; nc<M->NR; nc++)
       M->AddEntry(nr, nc,       PFPM *MScratch->GetEntry(nr, nc) +
                            conj(PFPM)*MScratch->GetEntry(nc, nr) );
   };


  Log(" MPZ block ...");
  Displacement[0]=LatticeBasisVectors[0][0]; 
  Displacement[1]=LatticeBasisVectors[0][1];
  Args->B = M ? MScratch : MPZ;
  for(ns=0; ns<NumSurfaces; ns++)
   for(nsp=0; nsp<NumSurfaces; nsp++)
    { Args->Sa=Surfaces[ns];
      Args->RowOffset=BFIndexOffset[ns];
      Args->Sb=Surfaces[nsp];
      Args->ColOffset=BFIndexOffset[nsp];

      NumCommonRegions=CountCommonRegions(Args->Sa, Args->Sb, CommonRegionIndices, Signs);
      if (NumCommonRegions==0) 
       continue;

      // see NOTE above
      nr1=CommonRegionIndices[0]; 
      Region1Contributes = RegionIsExtended[MAXLATTICE*nr1+0];
      if (NumCommonRegions==2)
       { nr2=CommonRegionIndices[1];
         Region2Contributes = RegionIsExtended[MAXLATTICE*nr2+0];
       }
      else
       Region2Contributes=false;

      if ( Region1Contributes==false && Region2Contributes==false )
       continue;

      if ( Region1Contributes==false )
       { SaveZeroed1=RegionMPs[nr1]->Zeroed; 
         RegionMPs[nr1]->Zero();
       };
      if ( NumCommonRegions==2 && Region2Contributes==false )
       { SaveZeroed2=RegionMPs[nr2]->Zeroed; 
         RegionMPs[nr2]->Zero();
       };

      AssembleBEMMatrixBlock(Args);

      if ( Region1Contributes==false )
       RegionMPs[nr1]->Zeroed=SaveZeroed1;
      if ( NumCommonRegions==2 && Region2Contributes==false )
       RegionMPs[nr2]->Zeroed=SaveZeroed2;

    };
  if (M)
   { for(int nr=0; nr<M->NR; nr++)
      for(int nc=0; nc<M->NR; nc++)
       M->AddEntry(nr, nc,       PFPZ *MScratch->GetEntry(nr, nc) +
                            conj(PFPZ)*MScratch->GetEntry(nc, nr) );
   };




  Log(" MZP block ...");
  Displacement[0]=LatticeBasisVectors[1][0]; 
  Displacement[1]=LatticeBasisVectors[1][1];
  Args->B = M ? MScratch : MZP;
  for(ns=0; ns<NumSurfaces; ns++)
   for(nsp=0; nsp<NumSurfaces; nsp++)
    { 
      Args->Sa=Surfaces[ns];
      Args->RowOffset=BFIndexOffset[ns];
      Args->Sb=Surfaces[nsp];
      Args->ColOffset=BFIndexOffset[nsp];

      NumCommonRegions=CountCommonRegions(Args->Sa, Args->Sb, CommonRegionIndices, Signs);
      if (NumCommonRegions==0) 
       continue;

      // see NOTE above
      nr1=CommonRegionIndices[0]; 
      Region1Contributes = RegionIsExtended[MAXLATTICE*nr1+1];
      if (NumCommonRegions==2)
       { nr2=CommonRegionIndices[1];
         Region2Contributes = RegionIsExtended[MAXLATTICE*nr2+1];
       }
      else
       Region2Contributes=false;

      if ( Region1Contributes==false && Region2Contributes==false )
       continue;

      if ( Region1Contributes==false )
       { SaveZeroed1=RegionMPs[nr1]->Zeroed; 
         RegionMPs[nr1]->Zero();
       };
      if ( NumCommonRegions==2 && Region2Contributes==false )
       { SaveZeroed2=RegionMPs[nr2]->Zeroed; 
         RegionMPs[nr2]->Zero();
       };

      AssembleBEMMatrixBlock(Args);

      if ( Region1Contributes==false )
       RegionMPs[nr1]->Zeroed=SaveZeroed1;
      if ( NumCommonRegions==2 && Region2Contributes==false )
       RegionMPs[nr2]->Zeroed=SaveZeroed2;

    };
  if (M)
   { for(int nr=0; nr<M->NR; nr++)
      for(int nc=0; nc<M->NR; nc++)
       M->AddEntry(nr, nc,       PFZP* MScratch->GetEntry(nr, nc) +
                            conj(PFZP)*MScratch->GetEntry(nc, nr) );
   };

 // 20120924: the zeroing/unzeroing of MPs may have wreaked havoc 
 //           on the internally-stored Eps/Mu arrays, so i will 
 //           quickly restore them here. 
 //           it must be said that this zeroing/unzeroing          
 //           paradigm is fairly kludgy, and it would be nice to 
 //           redo it more elegantly at some point.
 UpdateCachedEpsMuValues(StoredOmega);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::AssembleBEMMatrix(cdouble Omega, double kBloch[MAXLATTICE], HMatrix *M)
{
  if (NumLatticeBasisVectors!=2)
   ErrExit("%s:%i: not yet implemented",__FILE__,__LINE__);
  Log("Assembling PBC BEM matrix at (Omega,Kx,Ky)=(%s,%g,%g) (Mem=%lu)",z2s(Omega),kBloch[0],kBloch[1],GetMemoryUsage()/1048576);

  /*--------------------------------------------------------------*/
  /*- (re)allocate the matrix as necessary -----------------------*/
  /*--------------------------------------------------------------*/
  if ( M!=0 && (M->NR!=TotalBFs || M->NC!=TotalBFs ) )
   { M=0;
     Warn("wrong-sized HMatrix passed to AssembleBEMMatrix() (reallocating...)");
   };

  if(M==0)
   M=new HMatrix(TotalBFs, TotalBFs, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /*- make sure cached values of epsilon and mu are correct for   */
  /*- this frequency; also, if we are using PBCAcceleration and   */
  /*- the frequency has changed since we last computed the inner  */
  /*- cell blocks, we need to recompute those.                    */
  /*--------------------------------------------------------------*/
  if( StoredOmega != Omega )
   { UpdateCachedEpsMuValues(Omega);
     if (UsePBCAcceleration)
      AssembleInnerCellBlocks(0, 0);
   };

  /*--------------------------------------------------------------*/
  /*- stamp in the contributions of the innermost 9 lattice cells,*/
  /*- by looking them up in the MPP, MPM, ... internal storage    */
  /*- blocks if those are present, and otherwise by computing them*/
  /*- on the fly.                                                 */
  /*- note: PFPP = 'phase factor, plus-plus'                      */
  /*-       PFPZ = 'phase factor, plus-zero'                      */
  /*- etc.                                                        */
  /*--------------------------------------------------------------*/
  if ( UsePBCAcceleration )
   { 
     double *LBV[2];
     LBV[0]=LatticeBasisVectors[0];
     LBV[1]=LatticeBasisVectors[1];
     cdouble PFPP = exp( II* ( kBloch[0]*(LBV[0][0] + LBV[1][0]) + kBloch[1]*(LBV[0][1] + LBV[1][1]) ) );
     cdouble PFPM = exp( II* ( kBloch[0]*(LBV[0][0] - LBV[1][0]) + kBloch[1]*(LBV[0][1] - LBV[1][1]) ) );
     cdouble PFPZ = exp( II* ( kBloch[0]*(LBV[0][0]            ) + kBloch[1]*(LBV[0][1]            ) ) );
     cdouble PFZP = exp( II* ( kBloch[0]*(            LBV[1][0]) + kBloch[1]*(            LBV[1][1]) ) );
     for(int nr=0; nr<M->NR; nr++)
      for(int nc=0; nc<M->NR; nc++)
       M->SetEntry(nr, nc,  PFPP*MPP->GetEntry(nr, nc) + conj(PFPP)*MPP->GetEntry(nc,nr)
                          + PFPM*MPM->GetEntry(nr, nc) + conj(PFPM)*MPM->GetEntry(nc,nr)
                          + PFPZ*MPZ->GetEntry(nr, nc) + conj(PFPZ)*MPZ->GetEntry(nc,nr)
                          + PFZP*MZP->GetEntry(nr, nc) + conj(PFZP)*MZP->GetEntry(nc,nr)
                          + MZZ->GetEntry(nr,nc) 
                  );
   }
  else
   AssembleInnerCellBlocks(kBloch, M);

  /*--------------------------------------------------------------*/
  /*- add the contribution of lattice cells beyond the innermost 9*/
  /*--------------------------------------------------------------*/
  AddOuterCellContributions(kBloch, M);

  return M;

}

} //namespace scuff{
