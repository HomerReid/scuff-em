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

#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#define MAXFUNC 50

namespace scuff {

void GetReducedPotentials_Nearby(RWGSurface *S, const int ne,
                                 const double X0[3], const cdouble k,
                                 cdouble *p, cdouble a[3],
                                 cdouble dp[3], cdouble da[3][3],
                                 cdouble ddp[3][3], cdouble dcurla[3][3]);

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
/*  f[12], f[13] == real , imag  -2(X-X0)_x          Psi(|X-X0|)   */
/*  f[14], f[15] == real , imag  -2(X-X0)_y          Psi(|X-X0|)   */
/*  f[16], f[17] == real , imag  -2(X-X0)_z          Psi(|X-X0|)   */
/*******************************************************************/
typedef struct GRPIntegrandData
 { 
   double *Q;            // RWG basis function source/sink vertex  
   double PreFac;        // RWG basis function prefactor 
   const double *X0;     // field evaluation point 
   cdouble K;            // \sqrt{Eps*Mu} * frequency
   GBarAccelerator *GBA; // optional accelerator for PBC geometries
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
                                      cdouble *a, cdouble *Curla, cdouble *Gradp)
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
         { RWGEdge *E = S->Edges[ne];
           double rRel = VecDistance(X, E->Centroid ) / E->Radius;
           if (rRel > 5.0) 
            S->GetReducedPotentials(ne, X, K, GBA, a, Curla, Gradp);
           else
            { cdouble p, da[3][3], ddp[3][3], dcurla[3][3];
              GetReducedPotentials_Nearby(S, ne, X, K,
                                          &p, a, Gradp, da, ddp, dcurla);
              Curla[0] = da[1][2]-da[2][1];
              Curla[1] = da[2][0]-da[0][2];
              Curla[2] = da[0][1]-da[1][0];
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
     if (G->RegionMPs[RegionIndex]->IsPEC())
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
      { for(IncField *IF=IFList; IF; IF=IF->Next)
         if ( IF->RegionIndex == RegionIndex )
          { cdouble dEH[6];
            IF->GetFields(X, dEH);
            SixVecPlusEquals(EH, 1.0, dEH);
          };
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
  if ( XMatrix==0 || XMatrix->NC!=3 || XMatrix->NR==0 )
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
  if (KN && LDim>0)
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

  if (KN && LDim>0)
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

} // namespace scuff
