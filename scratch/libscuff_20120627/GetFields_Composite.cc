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

#include "libscuff.h"
#include "RWGComposite.h"
#include "FieldGrid.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#define USE_OPENMP
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#define MAXFUNC 50

namespace scuff {

#define II cdouble(0,1)

/*******************************************************************/
/* 'GetReducedPotentialsIntegrand' is the integrand routine passed */
/* to TriInt to compute the reduced potentials of a single panel.  */
/*                                                                 */
/* integrand values:                                               */
/*  f[0],  f[1]  == real , imag  (X-Q)_x            Phi(|X-X0|)    */
/*  f[2],  f[3]  == real , imag  (X-Q)_y            Phi(|X-X0|)    */
/*  f[4],  f[5]  == real , imag  (X-Q)_z            Phi(|X-X0|)    */
/*  f[6],  f[7]  == real , imag  [(X-Q) x (X-X0)]_x Psi(|X-X0|)    */
/*  f[8],  f[9]  == real , imag  [(X-Q) x (X-X0)]_y Psi(|X-X0|)    */
/*  f[10], f[11] == real , imag  [(X-Q) x (X-X0)]_z Psi(|X-X0|)    */
/*  f[12], f[13] == real , imag  2(X-X0)_x          Psi(|X-X0|)    */
/*  f[14], f[15] == real , imag  2(X-X0)_y          Psi(|X-X0|)    */
/*  f[16], f[17] == real , imag  2(X-X0)_z          Psi(|X-X0|)    */
/*******************************************************************/
typedef struct GRPIntegrandData
 { 
   double *Q;         // RWG basis function source/sink vertex  
   double PreFac;     // RWG basis function prefactor 
   const double *X0;  // field evaluation point 
   cdouble K;         // \sqrt{Eps*Mu} * frequency
 } GRPIData;

static void GRPIntegrand(double *X, void *parms, double *f)
{ 
  GRPIntegrandData *GRPID=(GRPIntegrandData *)parms;
  double *Q         = GRPID->Q;
  double PreFac     = GRPID->PreFac;
  const double *X0  = GRPID->X0;
  cdouble K         = GRPID->K;

  /* get the value of the RWG basis function at X */
  double fRWG[3];
  VecSub(X,Q,fRWG);
  VecScale(fRWG,PreFac);
  
  /* compute the scalar helmholtz GFs */
  double XmX0[3], fxR[3];
  VecSub(X,X0,XmX0);
  VecCross(fRWG,XmX0,fxR);
  double r=VecNorm(XmX0);
  cdouble Phi = exp(II*K*r) / (4.0*M_PI*r);
  cdouble Psi = (II*K - 1.0/r) * Phi / r;
  
  /* assemble integrand components */
  cdouble *zf=(cdouble *)f;
  zf[0]= fRWG[0] * Phi;
  zf[1]= fRWG[1] * Phi;
  zf[2]= fRWG[2] * Phi;
  zf[3]= fxR[0] * Psi;
  zf[4]= fxR[1] * Psi;
  zf[5]= fxR[2] * Psi;
  zf[6]= -2.0 * PreFac * XmX0[0] * Psi;
  zf[7]= -2.0 * PreFac * XmX0[1] * Psi;
  zf[8]= -2.0 * PreFac * XmX0[2] * Psi;

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
double MyQR10[]=
 { 1.304673574141418e-02, 3.333567215434143e-02, 
   6.746831665550773e-02, 7.472567457529027e-02, 
   1.602952158504878e-01, 1.095431812579910e-01,
   2.833023029353764e-01, 1.346333596549959e-01, 
   4.255628305091844e-01, 1.477621123573765e-01, 
   5.744371694908156e-01, 7.166976970646236e-01, 
   8.397047841495122e-01, 9.325316833444923e-01,
   9.869532642585859e-01, 1.477621123573765e-01,
   1.346333596549959e-01, 1.095431812579910e-01, 
   7.472567457529027e-02, 3.333567215434143e-02
 };

void GetEdgeContributionToGradp(const double *X0, double *V1, double *V2, 
                                cdouble K, cdouble *Gradp)
{
  int np;
  double u, w, X0mV1[3], V2mV1[3], X0mX[3], r;
  cdouble Psi;

// R = X0 - X = X0 - ( V1 + u*(V2-V1) )

  int NumPts=10;
  double *QR=MyQR10;

  VecSub(V2, V1, V2mV1);
  VecSub(X0, V1, X0mV1);
  Gradp[0]=Gradp[1]=Gradp[2]=0.0;

  for(np=0; np<NumPts; np++)
   { 
     u=QR[2*np+0]; w=QR[2*np+1];

     X0mX[0] = X0mV1[0] + u*V2mV1[0];
     X0mX[1] = X0mV1[1] + u*V2mV1[1];
     X0mX[2] = X0mV1[2] + u*V2mV1[2];

     r=VecNorm(X0mX);

     Psi = (II*K - 1.0/r) * exp(II*K*r) / (4.0*M_PI*r*r);

     Gradp[0] += w*X0mX[0]*Psi;
     Gradp[1] += w*X0mX[1]*Psi;
     Gradp[2] += w*X0mX[2]*Psi;
     
   };

  cdouble PreFac = -1.0*VecNorm(V2mV1);
  Gradp[0] *= PreFac;
  Gradp[1] *= PreFac;
  Gradp[2] *= PreFac;
 
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
/***************************************************************/
void GetReducedPotentials(RWGEdge *E, RWGPanel **Panels, 
                          double *Vertices,
                          const double *X, cdouble K,
                          cdouble *a, cdouble *Curla,
                          cdouble *Gradp)
{
  double *QP, *V1, *V2, *QM;
  double PArea, MArea;
  int mu;
  GRPIntegrandData MyGRPIData, *GRPID=&MyGRPIData;
  cdouble IP[9], IM[9];
  cdouble Gradp_Edge[3];

  /* get edge vertices */
  QP=Vertices + 3*(E->iQP);
  V1=Vertices + 3*(E->iV1);
  V2=Vertices + 3*(E->iV2);
  PArea=Panels[E->iPPanel]->Area;

  /* set up data structure passed to GRPIntegrand */
  GRPID->X0=X;
  GRPID->K=K;

  /* contribution of positive panel */
  GRPID->Q=QP;
  GRPID->PreFac = E->Length / (2.0*PArea);
  TriIntFixed(GRPIntegrand, 18, (void *)GRPID, QP, V1, V2, 25, (double *)IP);

  /* contribution of negative panel if present */
  if ( E->iQM!=-1 )
   { 
     MArea=Panels[E->iMPanel]->Area;
     QM=Vertices + 3*(E->iQM);
     GRPID->Q=QM;
     GRPID->PreFac = E->Length / (2.0*MArea);
     TriIntFixed(GRPIntegrand, 18, (void *)GRPID, V1, V2, QM, 25, (double *)IM);
     memset(Gradp_Edge, 0, 3*sizeof(cdouble));
   }
  else // if there is no negative panel then there is a line-charge edge 
   { memset(IM, 0, 9*sizeof(cdouble));
     GetEdgeContributionToGradp(X, V1, V2, K, Gradp_Edge); 
   };

  for(mu=0; mu<3; mu++) 
   { a[mu]     = IP[mu]   - IM[mu];
     Curla[mu] = IP[mu+3] - IM[mu+3];
     Gradp[mu] = IP[mu+6] - IM[mu+6] + Gradp_Edge[mu];
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetScatteredFields(RWGComposite *C, 
                        int SubRegion, 
                        const double *X, 
                        HVector *KN,
                        const cdouble Omega, 
                        const cdouble Eps,
                        const cdouble Mu,
                        cdouble EHS[6])
{ 
  memset(EHS, 0, 6*sizeof(cdouble));

  cdouble iwe=II*Omega*Eps;
  cdouble iwu=II*Omega*Mu;
  cdouble K=csqrt2(Eps*Mu)*Omega;

  int nps, nte, Offset, i;
  double Sign;
  cdouble KAlpha, NAlpha, a[3], Curla[3], Gradp[3];
  PartialSurface *PS;
  RWGEdge *E;
  for(nps=0; nps<C->NumPartialSurfaces; nps++)
   { 
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     if ( C->PSSubRegions[2*nps + 0] == SubRegion )
      Sign=1.0;
     else if ( C->PSSubRegions[2*nps + 1] == SubRegion )
      Sign=-1.0;
     else
      continue;

     /***************************************************************/
     /* now loop over panels on object's surface to get             */
     /* contributions to field at evaluation point.                 */
     /***************************************************************/
     PS=C->PartialSurfaces[nps];
     Offset = C->BFIndexOffset[nps];
     for(nte=0; nte<PS->NumTotalEdges; nte++)
      { 
        KAlpha = Sign*KN->GetEntry( Offset + 2*nte + 0 );
        NAlpha = Sign*KN->GetEntry( Offset + 2*nte + 1 );

        if ( nte<PS->NumEdges )
         E=PS->Edges[nte];
        else
         E=PS->HEdges[ nte - PS->NumEdges ];

        GetReducedPotentials(E, PS->Panels, C->Vertices, 
                             X, K, a, Curla, Gradp);

        for(i=0; i<3; i++)
         { EHS[i]   += ZVAC*( KAlpha*(iwu*a[i] - Gradp[i]/iwe) + NAlpha*Curla[i] );
           EHS[i+3] += -1.0*NAlpha*(iwe*a[i] - Gradp[i]/iwu) + KAlpha*Curla[i];
         };

      }; // for (ne=0 ... 

    }; // for(no=0 ... 
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   int nt, nTask;

   RWGGeometry *G;
   RWGComposite *C;
   int SubRegion;
   HMatrix *XMatrix;
   HMatrix *FMatrix;
   HVector *KN;
   IncField *IF;
   cdouble Omega;
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
 // RWGGeometry *G              = TD->G;
  RWGComposite *C             = TD->C;
  int SubRegion               = TD->SubRegion;
  HMatrix *XMatrix            = TD->XMatrix;
  HMatrix *FMatrix            = TD->FMatrix;
  HVector *KN                 = TD->KN;
  IncField *IFList            = TD->IF;   
  cdouble Omega               = TD->Omega;
  ParsedFieldFunc **PFFuncs   = TD->PFFuncs;
  int NumFuncs                = TD->NumFuncs;

  /***************************************************************/
  /* other local variables ***************************************/
  /***************************************************************/
  double X[3];
  cdouble EH[6], dEH[6];
  cdouble Eps, Mu;
  double dA[3]={0.0, 0.0, 0.0};
  IncField *IF;

  /***************************************************************/
  /* loop over all eval points (all rows of the XMatrix)         */
  /***************************************************************/
  int nt=0;
  for(int nr=0; nr<XMatrix->NR; nr++)
   { 
     nt++;
     if (nt==TD->nTask) nt=0;
     if (nt!=TD->nt) continue;

     X[0]=XMatrix->GetEntryD(nr, 0);
     X[1]=XMatrix->GetEntryD(nr, 1);
     X[2]=XMatrix->GetEntryD(nr, 2);

     //ObjectIndex = G->GetObjectIndex(X);
     //if (ObjectIndex==-1)
     // G->ExteriorMP->GetEpsMu(Omega, &Eps, &Mu);
     //else
     // G->Objects[ObjectIndex]->MP->GetEpsMu(Omega, &Eps, &Mu);
     Eps = C->EpsTF[SubRegion];
     Mu  = C->MuTF[SubRegion];
    
     /*--------------------------------------------------------------*/
     /*- get scattered fields at X                                   */
     /*--------------------------------------------------------------*/
     if (KN)
      GetScatteredFields(C, SubRegion, X, KN, Omega, Eps, Mu, EH);
     else
      memset(EH, 0, 6*sizeof(cdouble));

     /*--------------------------------------------------------------*/
     /*- add incident fields by summing contributions of all        -*/
     /*- IncFields whose sources lie in the same region as X        -*/
     /*--------------------------------------------------------------*/
if ( IFList && SubRegion==0 )
 { IFList->GetFields(X, dEH);
   SixVecPlusEquals(EH, 1.0, dEH);
 };

#if 0
     if (IFList)
      { for(IF=IFList; IF; IF=IF->Next)
         if ( IF->ObjectIndex == SubRegion /*ObjectIndex*/ )
          { IF->GetFields(X, dEH);
            SixVecPlusEquals(EH, 1.0, dEH);
          };
      };
#endif

     /*--------------------------------------------------------------*/
     /*- compute field functions ------------------------------------*/
     /*--------------------------------------------------------------*/
     for(int nf=0; nf<NumFuncs; nf++)
      FMatrix->SetEntry(nr, nf, PFFuncs[nf]->Eval(X, dA, EH, Eps, Mu));

   }; // for (nr=0; nr<XMatrix->NR; nr++)

  return 0;

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetFields(RWGComposite *C, int SubRegion,
                   IncField *IF, HVector *KN,
                   cdouble Omega, HMatrix *XMatrix,
                   HMatrix *FMatrix, char *FuncString,
                   int nThread)
{ 
  if (nThread <= 0) nThread = GetNumThreads();
 
  /***************************************************************/
  /* preprocess the Functions string to count the number of      */
  /* comma-separated function strings and verify that each string*/
  /* is a valid function                                         */
  /***************************************************************/
  char *FCopy;
  char *Funcs[MAXFUNC];
  int nf, NumFuncs;

  if (FuncString==NULL)
   FCopy=strdup("Ex,Ey,Ez,Hx,Hy,Hz"); // default is cartesian field components
  else
   FCopy=strdup(FuncString);

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

  /***************************************************************/
  /* the incident fields will most likely have been updated at   */
  /* the current frequency already by an earlier call to         */
  /* AssembleRHSVector(), but somewhat might call GetFields()    */
  /* to get information on just the incident fields before       */
  /* before setting up and solving the BEM problem, so we should */
  /* do this just to make sure.                                  */
  /***************************************************************/
  //UpdateIncFields(IF, Omega);

  /***************************************************************/
  /* fire off threads                                            */
  /***************************************************************/
  int nt;

  // set up an instance of ThreadData containing all fields
  // that are common to all threads, which we can subsequently
  // copy wholesale to initialize new ThreadData structures
  ThreadData ReferenceTD; 
  // ReferenceTD.G=this;
  ReferenceTD.C=C;
  ReferenceTD.SubRegion=SubRegion;
  ReferenceTD.XMatrix = XMatrix;
  ReferenceTD.FMatrix = FMatrix;
  ReferenceTD.KN=KN;
  ReferenceTD.IF=IF;
  ReferenceTD.Omega=Omega;
  ReferenceTD.PFFuncs=PFFuncs;
  ReferenceTD.NumFuncs=NumFuncs;

#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[nThread], *TD;
  pthread_t *Threads = new pthread_t[nThread];
  ReferenceTD.nTask=nThread;
  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDs[nt]);
     memcpy(TD, &ReferenceTD, sizeof(ThreadData));
     TD->nt=nt;

     if (nt+1 == nThread)
       GetFields_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, GetFields_Thread, (void *)TD);
   }
  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDs;

#else 
#ifndef USE_OPENMP
  nThread=ReferenceTD.nTask=1;
#else
  ReferenceTD.nTask=nThread*100;
#pragma omp parallel for schedule(dynamic,1), num_threads(nThread)
#endif
  for(nt=0; nt<ReferenceTD.nTask; nt++)
   { 
     ThreadData TD1;
     memcpy(&TD1, &ReferenceTD, sizeof(ThreadData));
     TD1.nt=nt;
     GetFields_Thread((void *)&TD1);
   };
#endif

  free(FCopy);
  for(nf=0; nf<NumFuncs; nf++)
   delete PFFuncs[nf];
  delete[] PFFuncs;

  return FMatrix;

}

/***************************************************************/
/* simple interface to GetFields *******************************/
/***************************************************************/
#if 0
void GetFields(RWGComposite *C, int SubRegion,
               IncField *IF, HVector *KN, cdouble Omega, double *X,
               cdouble *EH, int nThread)
{
  HMatrix XMatrix(1, 3, LHM_REAL, LHM_NORMAL, (void *)X);

  HMatrix FMatrix(1, 6, LHM_COMPLEX, LHM_NORMAL, (void *)EH);

  GetFields(C, SubRegion, IF, KN, Omega, &XMatrix, &FMatrix, 0, nThread);
} 
#endif


} // namespace scuff
