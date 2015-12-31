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
 * You should have received a copy of the GNU General Public License * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * GetFields.cc  -- libscuff class methods for computing scattered
 *               -- electric and magnetic fields 
 *
 * homer reid    -- 3/2007  -- 11/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libTriInt.h>
#include <libMDInterp.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "FieldGrid.h"
#include "PanelCubature.h"

#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#define NUMFIELDS 6 // Ex, Ey, Ez, Hx, Hy, Hz
#define MAXFUNC 40

namespace scuff {

RWGSurface *ResolveNE(RWGGeometry *G, int neFull,
                      int *pns, int *pne, int *pKNIndex);

void GetReducedPotentials_Nearby(RWGSurface *S, const int ne,
                                 const double X0[3], const cdouble k,
                                 cdouble *p, cdouble a[3],
                                 cdouble dp[3], cdouble da[3][3],
                                 cdouble ddp[3][3], cdouble dcurla[3][3]);

void GetReducedFields_Nearby(RWGSurface *S, const int ne,
                             const double X0[3], const cdouble k,
                             cdouble e[3], cdouble h[3]);

typedef struct GCX0Data
 { cdouble k;
   double *X0;
   GBarAccelerator *GBA;
 } GCX0Data;

void GCX0Integrand(double X[3], double b[3],
                   void *UserData, double W, double *Integral);

#define II cdouble(0,1)

/*******************************************************************/
/* 'GetReducedPotentialsIntegrand' is the integrand routine passed */
/* to TriInt to compute the reduced potentials of a single panel.  */
/*                                                                 */
/* ('reduced' in the context of this file means 'without certain   */
/*   constant prefactors.')                                        */
/*                                                                 */
/* The GBarInterp field may be used to specify an optional         */
/* Interp3D interpolator for the periodic Green's function. If     */
/* this field is NULL, then the computation is done using the usual*/
/* non-periodic Green's function, (e^{iKR}/4\piR).                 */
/*                                                                 */
/* integrand values:                                               */
/*  f[0],  f[1]  == real , imag  (X-Q)_x            Phi(|X-X0|)    */
/*  f[2],  f[3]  == real , imag  (X-Q)_y            Phi(|X-X0|)    */
/*  f[4],  f[5]  == real , imag  (X-Q)_z            Phi(|X-X0|)    */
/*  f[6],  f[7]  == real , imag  [(X-Q) x (X-X0)]_x Psi(|X-X0|)    */
/*  f[8],  f[9]  == real , imag  [(X-Q) x (X-X0)]_y Psi(|X-X0|)    */
/*  f[10], f[11] == real , imag  [(X-Q) x (X-X0)]_z Psi(|X-X0|)    */
/*  f[12], f[13] == real , imag  -2(X-X0)_x         Psi(|X-X0|)    */
/*  f[14], f[15] == real , imag  -2(X-X0)_y         Psi(|X-X0|)    */
/*  f[16], f[17] == real , imag  -2(X-X0)_z         Psi(|X-X0|)    */
/*******************************************************************/
typedef struct GRPIntegrandData
 { 
   double *Q;            // RWG basis function source/sink vertex  
   double PreFac;        // RWG basis function prefactor 
   const double *X0;     // field evaluation point 
   cdouble K;            // \sqrt{Eps*Mu} * frequency
   GBarAccelerator *GBA; // optional accelerator for PBC geometries
   bool ForceFullEwald;
 } GRPIData;

static void GRPIntegrand(double *X, void *parms, double *f)
{ 
  GRPIntegrandData *GRPID=(GRPIntegrandData *)parms;
  double *Q            = GRPID->Q;
  double PreFac        = GRPID->PreFac;
  const double *X0     = GRPID->X0;
  cdouble K            = GRPID->K;
  GBarAccelerator *GBA = GRPID->GBA;

  /* get the value of the RWG basis function at X */
  double fRWG[3];
  VecSub(X,Q,fRWG);
  VecScale(fRWG,PreFac);

  double R[3];
  VecSub(X,X0,R);

  cdouble *zf=(cdouble *)f;
  
  if (GBA==0)
   {
     /* compute the scalar functions Phi and Psi  */
     /* (Phi = e^{ikr}/4*pi*r, and Psi is defined */
     /* such that \nabla Phi = r*Psi              */
     double fxR[3];
     VecCross(fRWG,R,fxR);
     double r=VecNorm(R);
     cdouble Phi = exp(II*K*r) / (4.0*M_PI*r);
     cdouble Psi = (II*K - 1.0/r) * Phi / r;
  
     /* assemble integrand components */
     zf[0]= fRWG[0] * Phi;
     zf[1]= fRWG[1] * Phi;
     zf[2]= fRWG[2] * Phi;
     zf[3]= fxR[0] * Psi;
     zf[4]= fxR[1] * Psi;
     zf[5]= fxR[2] * Psi;
     zf[6]= -2.0 * PreFac * R[0] * Psi;
     zf[7]= -2.0 * PreFac * R[1] * Psi;
     zf[8]= -2.0 * PreFac * R[2] * Psi;
   }
  else
   { 
     cdouble G, dG[3];
     if (GRPID->ForceFullEwald)
      { cdouble GBarVD[8];
        GBarVDEwald(R, GBA->k, GBA->kBloch, GBA->LBV, GBA->LDim,
                    -1.0, true, GBarVD);
        G=GBarVD[0];
        dG[0]=GBarVD[1];
        dG[1]=GBarVD[2];
        dG[2]=GBarVD[3];
      } 
     else
      G=GetGBar(R, GBA, dG);

     /* assemble integrand components */
     zf[0]= fRWG[0] * G;
     zf[1]= fRWG[1] * G;
     zf[2]= fRWG[2] * G;
     zf[3]= (fRWG[1] * dG[2] - fRWG[2] * dG[1]);
     zf[4]= (fRWG[2] * dG[0] - fRWG[0] * dG[2]);
     zf[5]= (fRWG[0] * dG[1] - fRWG[1] * dG[0]);
     zf[6]= -2.0 * PreFac * dG[0];
     zf[7]= -2.0 * PreFac * dG[1];
     zf[8]= -2.0 * PreFac * dG[2];
   };

} 

/***************************************************************/
/* Compute the 'reduced potentials' of a single RWG basis      */
/* function.                                                   */
/*                                                             */
/* The 'reduced potentials' are dimensionless analogues of the */
/* scalar and vector potentials of the source distributions    */
/* described by the RWG function populated with unit strength: */
/*                                                             */
/*    p(x) = \int G(x,y) \nabla \cdot f(y) dy                  */
/*  a_i(x) = \int G(x,y) f_i(y) dy                             */
/*                                                             */
/* where G(x,y) is the scalar green's function and f(y) is the */
/* vector-valued RWG current at y.                             */
/***************************************************************/
void RWGSurface::GetReducedPotentials(int ne, const double *X, cdouble K,
                                      GBarAccelerator *GBA, 
                                      cdouble *a, cdouble *Curla, cdouble *Gradp,
                                      bool ForceFullEwald)
{
  double *QP, *V1, *V2, *QM;
  double PArea, MArea;
  RWGEdge *E;
  GRPIntegrandData MyGRPIData, *GRPID=&MyGRPIData;
  cdouble IP[9], IM[9];

  /* get edge vertices */
  E=Edges[ne];
  QP=Vertices + 3*(E->iQP);
  V1=Vertices + 3*(E->iV1);
  V2=Vertices + 3*(E->iV2);
  PArea=Panels[E->iPPanel]->Area;

  if (E->iQM == -1)
   QM=0;
  else
   { QM=Vertices + 3*(E->iQM);
     MArea=Panels[E->iMPanel]->Area;
   };

  /* set up data structure passed to GRPIntegrand */
  GRPID->X0=X;
  GRPID->K=K;
  GRPID->GBA = GBA;
  GRPID->ForceFullEwald=ForceFullEwald;

  /* contribution of positive panel */
  GRPID->Q=QP;
  GRPID->PreFac = E->Length / (2.0*PArea);
  TriIntFixed(GRPIntegrand, 18, (void *)GRPID, QP, V1, V2, 25, (double *)IP);

  /* contribution of negative panel if present */
  if (QM)
   { GRPID->Q=QM;
     GRPID->PreFac = E->Length / (2.0*MArea);
     TriIntFixed(GRPIntegrand, 18, (void *)GRPID, V1, V2, QM, 25, (double *)IM);
   }
  else
   memset(IM, 0, 9*sizeof(cdouble));

  for(int Mu=0; Mu<3; Mu++) 
   { a[Mu]     = IP[Mu]   - IM[Mu];
     Curla[Mu] = IP[Mu+3] - IM[Mu+3];
     Gradp[Mu] = IP[Mu+6] - IM[Mu+6];
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetScatteredFields(RWGGeometry *G, const double *X, const int RegionIndex,
                        HVector *KN, const cdouble Omega, GBarAccelerator *GBA,
                        cdouble EHS[6])
{ 
  memset(EHS, 0, 6*sizeof(cdouble));

  cdouble Eps=G->EpsTF[RegionIndex];
  cdouble Mu=G->MuTF[RegionIndex];
  cdouble iwe=II*Omega*Eps;
  cdouble iwu=II*Omega*Mu;
  cdouble K=csqrt2(Eps*Mu)*Omega;

  RWGSurface *S;
  int i, ne, Offset;
  cdouble KAlpha, NAlpha, a[3], Curla[3], Gradp[3];
  double Sign;
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     S=G->Surfaces[ns];
     Offset=G->BFIndexOffset[ns];

     /*****************************************************************/
     /* figure out the sign of the contribution of currents on this   */
     /* surface to the field at the evaluation point.                 */
     /*****************************************************************/
     if ( S->RegionIndices[0] == RegionIndex )
      Sign=+1.0;
     else if ( S->RegionIndices[1] == RegionIndex )
      Sign=-1.0;
     else
      continue; // in this case S does not contribute to field at eval pt

     /***************************************************************/ 
     /* now loop over all basis functions on the surface to         */ 
     /* get contributions to the field at the evaluation point.     */ 
     /***************************************************************/
     for(ne=0; ne<S->NumEdges; ne++)
      { 
        if ( S->IsPEC )
         { 
           KAlpha = Sign*KN->GetEntry( Offset + ne );
           NAlpha = 0.0;
         }
        else
         { KAlpha = Sign*KN->GetEntry( Offset + 2*ne + 0 );
           NAlpha = Sign*KN->GetEntry( Offset + 2*ne + 1 );
         };
      
        if (RWGGeometry::UseGetFieldsV2P0==false)
         S->GetReducedPotentials(ne, X, K, GBA, a, Curla, Gradp);
        else
         { 
           RWGEdge *E = S->Edges[ne];
           double rRel = VecDistance(X, E->Centroid ) / E->Radius;
           if (rRel > 3.0) 
            S->GetReducedPotentials(ne, X, K, GBA, a, Curla, Gradp);
           else if (GBA==0)
            { cdouble p, da[3][3], ddp[3][3], dcurla[3][3];
              GetReducedPotentials_Nearby(S, ne, X, K,
                                          &p, a, Gradp, da, ddp, dcurla);
              Curla[0] = da[1][2]-da[2][1];
              Curla[1] = da[2][0]-da[0][2];
              Curla[2] = da[0][1]-da[1][0];
            }
           else if (GBA->LDim==1)
            {
              S->GetReducedPotentials(ne, X, K, GBA, a, Curla, Gradp, true);

              double XmL[3];
              XmL[1]=X[1];
              XmL[2]=X[2];
              for(int nx=-1; nx<=1; nx++)
               { double L = nx*GBA->LBV[0][0];
                 XmL[0]=X[0] - L;
                 cdouble aDelta[3], GradpDelta[3];
                 cdouble p, da[3][3], ddp[3][3], dcurla[3][3];
                 GetReducedPotentials_Nearby(S, ne, XmL, K,
                                             &p, aDelta, GradpDelta, da, ddp, dcurla);
                 cdouble PhaseFactor = exp(II*L*(GBA->kBloch[0]));
                 for(int n=0; n<3; n++)
                  { int np1=(n+1)%3, np2=(n+2)%3;
                    a[n]     += PhaseFactor*aDelta[n];
                    Gradp[n] += PhaseFactor*GradpDelta[n];
                    Curla[n] += PhaseFactor*(da[np1][np2]-da[np2][np1]);
                  };
               };
            };
         };

        for(i=0; i<3; i++)
         { EHS[i]   += ZVAC*( KAlpha*(iwu*a[i] - Gradp[i]/iwe) + NAlpha*Curla[i] );
           EHS[i+3] += -1.0*NAlpha*(iwe*a[i] - Gradp[i]/iwu) + KAlpha*Curla[i];
         };

      }; // for (ne=0 ... 

    }; // for(ns=0 ... 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   int nt, NumTasks;

   RWGGeometry *G;
   HMatrix *XMatrix;
   HMatrix *FMatrix;
   HVector *KN;
   IncField *IF;
   cdouble Omega;
   GBarAccelerator **RegionGBAs;
   ParsedFieldFunc **PFFuncs;
   int NumFuncs;

 } ThreadData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void *GetFields_Thread(void *data)
{
  ThreadData *TD=(ThreadData *)data;

#ifdef USE_PTHREAD
  SetCPUAffinity(TD->nt);
#endif

  /***************************************************************/
  /* fields unpacked from thread data structure ******************/
  /***************************************************************/
  RWGGeometry *G                 = TD->G;
  HMatrix *XMatrix               = TD->XMatrix;
  HMatrix *FMatrix               = TD->FMatrix;
  HVector *KN                    = TD->KN;
  IncField *IFList               = TD->IF;
  cdouble Omega                  = TD->Omega;
  GBarAccelerator **RegionGBAs   = TD->RegionGBAs;
  ParsedFieldFunc **PFFuncs      = TD->PFFuncs;
  int NumFuncs                   = TD->NumFuncs;

  /***************************************************************/
  /* loop over all eval points (all rows of the XMatrix)         */
  /***************************************************************/
  int nt=0;
  for(int nr=0; nr<XMatrix->NR; nr++)
   { 
     nt++;
     if (nt==TD->NumTasks) nt=0;
     if (nt!=TD->nt) continue;

     double X[3];
     X[0]=XMatrix->GetEntryD(nr, 0);
     X[1]=XMatrix->GetEntryD(nr, 1);
     X[2]=XMatrix->GetEntryD(nr, 2);
   
     cdouble EH[6];
     memset(EH, 0, 6*sizeof(cdouble));

     int RegionIndex = G->GetRegionIndex(X);
     if (RegionIndex<0 || G->RegionMPs[RegionIndex]->IsPEC())
      continue;

     cdouble Eps = G->EpsTF[RegionIndex];
     cdouble Mu  = G->MuTF[RegionIndex];
     GBarAccelerator *GBA = RegionGBAs ? RegionGBAs[RegionIndex] : 0;
    
     /*--------------------------------------------------------------*/
     /*- get scattered fields at X                                   */
     /*--------------------------------------------------------------*/
     if (KN)
      GetScatteredFields(G, X, RegionIndex, KN, Omega, GBA, EH);

     /*--------------------------------------------------------------*/
     /*- add incident fields by summing contributions of all        -*/
     /*- IncFields whose sources lie in the same region as X        -*/
     /*--------------------------------------------------------------*/
     if (IFList)
      for(IncField *IF=IFList; IF; IF=IF->Next)
       if ( IF->RegionIndex == RegionIndex )
        { cdouble dEH[6];
          IF->GetFields(X, dEH);
          VecPlusEquals(EH, 1.0, dEH, 6);
        };

     /*--------------------------------------------------------------*/
     /*- compute field functions ------------------------------------*/
     /*--------------------------------------------------------------*/
     double dA[3]={1.0, 0.0, 0.0};
     for(int nf=0; nf<NumFuncs; nf++)
      FMatrix->SetEntry(nr, nf, PFFuncs[nf]->Eval(X, dA, EH, Eps, Mu));

   }; // for (nr=0; nr<XMatrix->NR; nr++)

  return 0;

} 

/***************************************************************/
/* set kBloch=NULL for non-PBC geometries **********************/
/***************************************************************/
HMatrix *RWGGeometry::GetFields(IncField *IF, HVector *KN, 
                                cdouble Omega, double *kBloch, 
                                HMatrix *XMatrix, HMatrix *FMatrix, char *FuncString)
{ 
  int NumThreads = GetNumThreads();
 
  /***************************************************************/
  /* preprocess the Functions string to count the number of      */
  /* comma-separated function strings and verify that each string*/
  /* is a valid function                                         */
  /***************************************************************/
  char *FCopy;
  char *Funcs[MAXFUNC];
  int nf, NumFuncs;

  if (FuncString==NULL)
   FCopy=strdupEC("Ex,Ey,Ez,Hx,Hy,Hz"); // default is cartesian field components
  else
   FCopy=strdupEC(FuncString);

  NumFuncs=Tokenize(FCopy, Funcs, MAXFUNC, ",");

  ParsedFieldFunc **PFFuncs = new ParsedFieldFunc *[NumFuncs];
  for(nf=0; nf<NumFuncs; nf++)
   PFFuncs[nf] = new ParsedFieldFunc(Funcs[nf]);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( XMatrix==0 || XMatrix->NC<3 || XMatrix->NR==0 )
   ErrExit("wrong-size XMatrix (%ix%i) passed to GetFields",XMatrix->NR,XMatrix->NC);
  if (FMatrix==0) 
   FMatrix=new HMatrix(XMatrix->NR, NumFuncs, LHM_COMPLEX);
  else if ( (FMatrix->NR != XMatrix->NR) || (FMatrix->NC!=NumFuncs) ) 
   { Warn(" ** warning: wrong-size FMatrix passed to GetFields(); allocating new matrix");
     FMatrix=new HMatrix(XMatrix->NR, NumFuncs, LHM_COMPLEX);
   };

  if (LogLevel >= SCUFF_VERBOSELOGGING)
   Log("Computing fields at %i evaluation points...",XMatrix->NR);

  /***************************************************************/
  /* the incident fields will most likely have been updated at   */
  /* the current frequency already by an earlier call to         */
  /* AssembleRHSVector(), but someone might call GetFields()     */
  /* to get information on just the incident fields before       */
  /* before setting up and solving the BEM problem, so we should */
  /* do this just to make sure.                                  */
  /***************************************************************/
  UpdateIncFields(IF, Omega, kBloch);

  /***************************************************************/
  /* For the periodic-boundary-condition case, we need to        */
  /* initialize accelerator objects to accelerate computation of */
  /* the periodic Green's function in each extended region of    */
  /* the geometry.                                               */
  /***************************************************************/
  GBarAccelerator **RegionGBAs=0;
  if (KN && LBasis)
   { RegionGBAs=(GBarAccelerator **)mallocEC(NumRegions*sizeof(GBarAccelerator *));
     for(int nr=0; nr<NumRegions; nr++)
      if ( ! ( RegionMPs[nr]->IsPEC() ) )
       RegionGBAs[nr]=CreateRegionGBA(nr, Omega, kBloch, XMatrix);
   };

  /***************************************************************/
  /* fire off threads                                            */
  /***************************************************************/
  int nt;

  // set up an instance of ThreadData containing all fields
  // that are common to all threads, which we can subsequently
  // copy wholesale to initialize new ThreadData structures.
  // Note in this version of the routine there are no periodic 
  // boundary conditions and hence no interpolator object for 
  // the periodic Green's function (so GBarInterp=0).
  ThreadData ReferenceTD; 
  ReferenceTD.G=this;
  ReferenceTD.XMatrix = XMatrix;
  ReferenceTD.FMatrix = FMatrix;
  ReferenceTD.KN=KN;
  ReferenceTD.IF=IF;
  ReferenceTD.Omega=Omega;
  ReferenceTD.RegionGBAs=RegionGBAs;
  ReferenceTD.PFFuncs=PFFuncs;
  ReferenceTD.NumFuncs=NumFuncs;

#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[NumThreads], *TD;
  pthread_t *Threads = new pthread_t[NumThreads];
  ReferenceTD.NumTasks=NumThreads;
  for(nt=0; nt<NumThreads; nt++)
   { 
     TD=&(TDs[nt]);
     memcpy(TD, &ReferenceTD, sizeof(ThreadData));
     TD->nt=nt;

     if (nt+1 == NumThreads )
       GetFields_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, GetFields_Thread, (void *)TD);
   }
  for(nt=0; nt<NumThreads-1; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDs;

#else 
#ifndef USE_OPENMP
  NumThreads=ReferenceTD.NumTasks=1;
#else
  ReferenceTD.NumTasks=NumThreads*100;
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(nt=0; nt<ReferenceTD.NumTasks; nt++)
   { 
     ThreadData TD1;
     memcpy(&TD1, &ReferenceTD, sizeof(ThreadData));
     TD1.nt=nt;
     GetFields_Thread((void *)&TD1);
   };
#endif

  /***************************************************************/
  /* deallocate temporary storage ********************************/
  /***************************************************************/
  free(FCopy);
  for(nf=0; nf<NumFuncs; nf++)
   delete PFFuncs[nf];
  delete[] PFFuncs;

  if (KN && LBasis)
   { for(int nr=0; nr<NumRegions; nr++)
      if (RegionGBAs[nr])
       DestroyGBarAccelerator(RegionGBAs[nr]);
     free(RegionGBAs);
   };

  return FMatrix;

}
  
/***************************************************************/
/* simpler interface to GetFields with only a single eval point*/
/***************************************************************/
void RWGGeometry::GetFields(IncField *IF, HVector *KN, cdouble Omega, 
                            double *kBloch, double *X, cdouble *EH)
{
  HMatrix XMatrix(1, 3, LHM_REAL, LHM_NORMAL, (void *)X);
  HMatrix FMatrix(1, 6, LHM_COMPLEX, LHM_NORMAL, (void *)EH);
  GetFields(IF, KN, Omega, kBloch, &XMatrix, &FMatrix, 0);
} 

/***************************************************************/
/* non-PBC entry points to GetFields                           */
/***************************************************************/
HMatrix *RWGGeometry::GetFields(IncField *IF, HVector *KN, cdouble Omega, 
                                HMatrix *XMatrix, HMatrix *FMatrix, char *FuncString)
{ return GetFields(IF, KN, Omega, 0, XMatrix, FMatrix, FuncString); }

void RWGGeometry::GetFields(IncField *IF, HVector *KN, cdouble Omega, 
                            double *X, cdouble *EH)
{ GetFields(IF, KN, Omega, 0, X, EH); }

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void GetDistanceToPanel(double X0[3], RWGSurface *S, int np)
{
  RWGPanel *P=S->Panels[np];
  double XdN=VecDot(X0, P->ZHat);

  double XPerp[3], XPar[3];
  l
}

void GetDistanceToBF(double X0[3], RWGSurface *S, int ne)
{
  RWGEdge *E=S->Edges[ne];
  double Distance = GetDistanceToPanel(X0, S, E->iPPanel);
  if (E->iMPanel!=-1)
   Distance = fmin(Distance, GetDistanceToPanel(X0, S, E->iMPanel);
   
}
#endif

/***************************************************************/
/* ehMatrix is a matrix whose columns may be dot-producted with*/
/* the KN vector (BEM system solution vector) to yield         */
/* components of the scattered E and H fields.                 */
/* More specifically, for Mu=0...5, the (6*nx + Mu)th column   */
/* of ehMatrix is dotted into KN to yield the Muth component   */
/* of the field six-vector F=\{ E \choose H \}.                */
/***************************************************************/
HMatrix *GetehMatrix(RWGGeometry *G, cdouble Omega, double *kBloch,
                     HMatrix *XMatrix, HMatrix *ehMatrix=0,
                     int ColumnOffset=0)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NE  = G->TotalEdges;
  int NBF = G->TotalBFs;
  int NX  = XMatrix->NR;
  if (     ehMatrix==0 
       || (ehMatrix->NR != NBF) 
       || (ehMatrix->NC != 6*NX) 
     )
   { if (ehMatrix) 
      { Warn("wrong-size ehMatrix passed to GetehMatrix; reallocating");
        delete ehMatrix;
      };
     ehMatrix = new HMatrix(NBF, 6*NX, LHM_COMPLEX);
   };
  ehMatrix->Zero();

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
  /***************************************************************/
  /***************************************************************/
  cdouble *ZRels   = new cdouble[G->NumRegions];
  cdouble *ks      = new cdouble[G->NumRegions];
  for(int nr=0; nr<G->NumRegions; nr++)
   { cdouble EpsRel, MuRel;
     G->RegionMPs[nr]->GetEpsMu(Omega, &EpsRel, &MuRel);
     ZRels[nr] = sqrt(MuRel/EpsRel);
     ks[nr]    = sqrt(MuRel*EpsRel) * Omega;
   };

  /***************************************************************/
  /* For the periodic-boundary-condition case, we need to        */
  /* initialize accelerator objects to accelerate computation of */
  /* the periodic Green's function in each extended region of    */
  /* the geometry.                                               */
  /***************************************************************/
  GBarAccelerator **RegionGBAs=0;
  int NumRegions=G->NumRegions;
  if (G->LBasis)
   { RegionGBAs=
      (GBarAccelerator **)mallocEC(NumRegions*sizeof(RegionGBAs[0]));
     for(int nr=0; nr<NumRegions; nr++)
      if ( ! ( G->RegionMPs[nr]->IsPEC() ) )
       RegionGBAs[nr]=G->CreateRegionGBA(nr, Omega, kBloch, XMatrix);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NENX=NE*NX;
#ifndef USE_OPENMP
  Log("Computing ehMatrix entries at %i points",NX);
#else
  int NumThreads=GetNumThreads();
  Log("Computing ehMatrix entries (%i threads) at %i points",NumThreads,NX);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nenx=0; nenx<NENX; nenx++)
   { 
     int nx     = nenx / NE;
     int neFull = nenx % NE;

     int ns, ne, nbf;
     ResolveNE(G, neFull, &ns, &ne, &nbf);
     RWGSurface *S = G->Surfaces[ns];
     RWGEdge *E    = S->Edges[ne];

     double X[3];
     X[0]=XMatrix->GetEntryD(nx,ColumnOffset+0);
     X[1]=XMatrix->GetEntryD(nx,ColumnOffset+1);
     X[2]=XMatrix->GetEntryD(nx,ColumnOffset+2);
     int RegionIndex = G->GetRegionIndex(X);
     if (RegionIndex==-1) continue; // inside a closed PEC surface
   
     double Sign=0.0;
     if      (G->Surfaces[ns]->RegionIndices[0]==RegionIndex) 
      Sign=+1.0;
     else if (G->Surfaces[ns]->RegionIndices[1]==RegionIndex)
      Sign=-1.0;
     else 
      continue;
   
     cdouble k    = ks[RegionIndex];
     cdouble ZRel = ZRels[RegionIndex];

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     cdouble GC[6];
     GCX0Data MyData, *Data=&MyData;
     Data->X0  = X;
     Data->k   = ks[RegionIndex];
     Data->GBA = RegionGBAs ? RegionGBAs[RegionIndex] : 0;

     double rRel = VecDistance(X, E->Centroid) / E->Radius;
     const int IDim=12;
     if (rRel >= rRelOuterThreshold)
      { 
        GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
                       IDim, LowOrder, (double *)GC);
      }
     else if (rRel>=rRelInnerThreshold)
      { 
        GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
                       IDim, HighOrder, (double *)GC);
      }
     else
      { 
        GetReducedFields_Nearby(S, ne, X, k, GC+0, GC+3);
        GC[3] /= (-II*k);
        GC[4] /= (-II*k);
        GC[5] /= (-II*k);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
cdouble OldGCLow[6], OldGCHigh[6];
GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
               IDim, LowOrder, (double *)OldGCLow);
GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
               IDim, HighOrder, (double *)OldGCHigh);

#pragma omp critical
{
  double MaxRDLow=0.0, MaxRDHigh=0.0;
  for(int n=0; n<6; n++)
   { MaxRDLow = fmax(MaxRDLow, RD(OldGCLow[n], GC[n]));
     MaxRDHigh = fmax(MaxRDHigh, RD(OldGCHigh[n], GC[n]));
   };
  
 printf("%+g %.5e 0 %i %e ",X[2],rRel,RegionIndex,MaxRDLow);
 for(int n=0; n<6; n++)
  printf("%+12.4e %+12.4e ",real(OldGCLow[n]),imag(OldGCLow[n]));
 printf("\n");
 printf("%+g %.5e 1 %i %e ",X[2],rRel,RegionIndex,MaxRDHigh);
 for(int n=0; n<6; n++)
  printf("%+12.4e %+12.4e ",real(OldGCHigh[n]),imag(OldGCHigh[n]));
 printf("\n");
 printf("%+g %.5e 2 %i %e ",X[2],rRel,RegionIndex,0.0);
 for(int n=0; n<6; n++)
  printf("%+12.4e %+12.4e ",real(GC[n]),imag(GC[n]));
 printf("\n");
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

        if (RegionGBAs)
         { cdouble GC1[6], GC2[6];
           int Order=4;
           GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
                          IDim, Order, (double *)GC1);
           Data->GBA = 0;
           GetBFCubature2(G, ns, ne, GCX0Integrand, (void *)Data,
                          IDim, Order, (double *)GC2);
           for(int Mu=0; Mu<6; Mu++) 
            GC[Mu] += (GC1[Mu] - GC2[Mu]);
         };
      };

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#if 0
E = ik*Z0 * Zr * k*g + ik*n*c
  = ik*Z0 * Zr * k*g - ik*Z0*nScuff*c
H =        -ik * k*c + (ik/(Z0*Zr)) * n*c
H =        -ik * k*c - (ik/Zr) *nScuff*c
#endif
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     cdouble *GG=GC+0, *CC=GC+3;
     cdouble EKFactor =      Sign*II*k*ZRel*ZVAC;
     cdouble HKFactor = -1.0*Sign*II*k;
     cdouble ENFactor = -1.0*Sign*II*k*ZVAC;
     cdouble HNFactor = -1.0*Sign*II*k/ZRel;

     ehMatrix->SetEntry(nbf, 6*nx + 0, EKFactor * GG[0] );
     ehMatrix->SetEntry(nbf, 6*nx + 1, EKFactor * GG[1] );
     ehMatrix->SetEntry(nbf, 6*nx + 2, EKFactor * GG[2] );
     ehMatrix->SetEntry(nbf, 6*nx + 3, HKFactor * CC[0] );
     ehMatrix->SetEntry(nbf, 6*nx + 4, HKFactor * CC[1] );
     ehMatrix->SetEntry(nbf, 6*nx + 5, HKFactor * CC[2] );

     if ( !(S->IsPEC) )
      { ehMatrix->SetEntry(nbf+1, 6*nx + 0, ENFactor * CC[0] );
        ehMatrix->SetEntry(nbf+1, 6*nx + 1, ENFactor * CC[1] );
        ehMatrix->SetEntry(nbf+1, 6*nx + 2, ENFactor * CC[2] );
        ehMatrix->SetEntry(nbf+1, 6*nx + 3, HNFactor * GG[0] );
        ehMatrix->SetEntry(nbf+1, 6*nx + 4, HNFactor * GG[1] );
        ehMatrix->SetEntry(nbf+1, 6*nx + 5, HNFactor * GG[2] );
      };

   }; // for(int nenx=0; nenx<NENX; nenx++)

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (RegionGBAs)
   { for(int nr=0; nr<NumRegions; nr++)
      if (RegionGBAs[nr])
       DestroyGBarAccelerator(RegionGBAs[nr]);
     free(RegionGBAs);
   };

  delete[] ZRels;
  delete[] ks;

  return ehMatrix;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RWGGeometry::GetFields2(IncField *IFList, HVector *KN,
                                 cdouble Omega, double *kBloch,
                                 HMatrix *XMatrix, HMatrix *FMatrix)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( XMatrix==0 || XMatrix->NC<3 || XMatrix->NR==0 )
   ErrExit("wrong-size XMatrix (%ix%i) passed to GetFields",
            XMatrix->NR,XMatrix->NC);

  int NX=XMatrix->NR;
  if (LogLevel >= SCUFF_VERBOSELOGGING)
   Log("Computing fields at %i evaluation points...",NX);

  /***************************************************************/
  /* (re)allocate output matrix as necessary *********************/
  /***************************************************************/
  if (FMatrix==0 || FMatrix->NR!=NX || FMatrix->NC!=NUMFIELDS)
   { if (FMatrix)
      { Warn(" ** warning: wrong-size FMatrix passed to GetFields(); reallocating");
        delete FMatrix;
      };
     FMatrix=new HMatrix(NX, NUMFIELDS, LHM_COMPLEX);
   };
  FMatrix->Zero();

  /***************************************************************/
  /* the incident fields will most likely have been updated at   */
  /* the current frequency already by an earlier call to         */
  /* AssembleRHSVector(), but someone might call GetFields()     */
  /* to get information on just the incident fields before       */
  /* before setting up and solving the BEM problem, so we should */
  /* do this just to make sure.                                  */
  /***************************************************************/
  if (IFList)
   UpdateIncFields(IFList, Omega, kBloch);

  /***************************************************************/
  /* get contributions of surface currents if present ************/
  /***************************************************************/
  if (KN)
   {
     HMatrix *ehMatrix = GetehMatrix(this, Omega, kBloch, XMatrix);
     HMatrix KNMatrix(1, TotalBFs, LHM_COMPLEX, LHM_NORMAL, (void *)KN->ZV);
     HMatrix *FMatrixT = new HMatrix(1, 6*NX, LHM_COMPLEX);
     KNMatrix.Multiply(ehMatrix, FMatrixT);
     for(int nx=0; nx<NX; nx++)
      for(int Mu=0; Mu<6; Mu++)
       FMatrix->SetEntry(nx, Mu, FMatrixT->GetEntry(0, 6*nx + Mu));

     delete ehMatrix;
     delete FMatrixT;
   };

  /***************************************************************/
  /* add contributions of incident fields if present *************/
  /***************************************************************/
  if (IFList)
   for(int nx=0; nx<NX; nx++)
    { 
      double X[3];
      XMatrix->GetEntriesD(nx,":",X);
      int RegionIndex = GetRegionIndex(X);
      if (RegionIndex==-1) continue; // inside a closed PEC surface

      for(IncField *IF=IFList; IF; IF=IF->Next)
       if ( IF->RegionIndex == RegionIndex )
        { cdouble EH[6];
          IF->GetFields(X, EH);
          for(int Mu=0; Mu<6; Mu++)
           FMatrix->AddEntry(nx, Mu, EH[Mu]);
        };
    };

  return FMatrix;
         
}

} // namespace scuff
