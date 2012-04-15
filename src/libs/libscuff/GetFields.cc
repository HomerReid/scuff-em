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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

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
   double *Q;           // RWG basis function source/sink vertex 
   double PreFac;       // RWG basis function prefactor 
   const double *X0;    // field evaluation point 
   cdouble K;           // \sqrt{Eps*Mu} * frequency
 } GRPIData;

static void GRPIntegrand(double *X, void *parms, double *f)
{ 
  GRPIntegrandData *GRPID=(GRPIntegrandData *)parms;
  double fRWG[3], XmX0[3], fxR[3];
  cdouble Phi, Psi, *L=(cdouble *)f;
  double r;
  cdouble K;
  double PreFac=GRPID->PreFac;

  /* get the value of the RWG basis function at XP */
  VecSub(X,GRPID->Q,fRWG);
  VecScale(fRWG,PreFac);
  
  VecSub(X,GRPID->X0,XmX0);
  VecCross(fRWG,XmX0,fxR);
  r=VecNorm(XmX0);

  K=GRPID->K;
  
  Phi= exp(II*K*r) / (4.0*M_PI*r);
  Psi= (II*K - 1.0/r) * Phi / r;
  
  L[0]= fRWG[0] * Phi;
  L[1]= fRWG[1] * Phi;
  L[2]= fRWG[2] * Phi;
  L[3]= fxR[0] * Psi;
  L[4]= fxR[1] * Psi;
  L[5]= fxR[2] * Psi;
  L[6]= -2.0 * PreFac * XmX0[0] * Psi;
  L[7]= -2.0 * PreFac * XmX0[1] * Psi;
  L[8]= -2.0 * PreFac * XmX0[2] * Psi;

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
void RWGObject::GetReducedPotentials(int ne, const double *X, cdouble K,
                                     cdouble *a, cdouble *Curla,
                                     cdouble *Gradp)
{
  double *QP, *V1, *V2, *QM;
  double PArea, MArea;
  RWGEdge *E;
  int mu;
  GRPIntegrandData MyGRPIData, *GRPID=&MyGRPIData;
  cdouble IP[9], IM[9];

  /* get edge vertices */
  E=Edges[ne];
  QP=Vertices + 3*(E->iQP);
  V1=Vertices + 3*(E->iV1);
  V2=Vertices + 3*(E->iV2);
  QM=Vertices + 3*(E->iQM);
  PArea=Panels[E->iPPanel]->Area;
  MArea=Panels[E->iMPanel]->Area;

  /* set up data structure passed to GRPIntegrand */
  GRPID->X0=X;
  GRPID->K=K;

  /* contribution of positive panel */
  GRPID->Q=QP;
  GRPID->PreFac = E->Length / (2.0*PArea);
  TriIntFixed(GRPIntegrand, 18, (void *)GRPID, QP, V1, V2, 25, (double *)IP);

  /* contribution of negative panel */
  GRPID->Q=QM;
  GRPID->PreFac = E->Length / (2.0*MArea);
  TriIntFixed(GRPIntegrand, 18, (void *)GRPID, V1, V2, QM, 25, (double *)IM);

  for(mu=0; mu<3; mu++) 
   { a[mu]     = IP[mu]   - IM[mu];
     Curla[mu] = IP[mu+3] - IM[mu+3];
     Gradp[mu] = IP[mu+6] - IM[mu+6];
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   int nt, nThread;

   RWGGeometry *G;   
   const double *X;            /* eval point */
   cdouble Omega;
   cdouble Eps, Mu;
   RWGObject *ObjectInQuestion;
   HVector *KN;
   cdouble EH[6];

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
  RWGGeometry *G              = TD->G;
  const double *X             = TD->X;
  cdouble Omega               = TD->Omega;
  cdouble Eps                 = TD->Eps;
  cdouble Mu                  = TD->Mu;
  RWGObject *ObjectInQuestion = TD->ObjectInQuestion;
  HVector *KN                 = TD->KN;
  cdouble *EH                 = TD->EH;

  /***************************************************************/
  /* other local fields ******************************************/
  /***************************************************************/
  double *DKN;
  cdouble *ZKN;
  cdouble iwe=II*Omega*Eps;
  cdouble iwu=II*Omega*Mu;
  cdouble K=csqrt2(Eps*Mu)*Omega;
  if ( KN->RealComplex==LHM_REAL )
   { DKN=KN->DV;
     ZKN=NULL;
   }
  else
   { 
     DKN=NULL;
     ZKN=KN->ZV;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGObject *O;
  int i, ne, no, Type, Offset;
  cdouble KAlpha, NAlpha, a[3], Curla[3], Gradp[3];
  double Sign;
  int nt=0;
  memset(EH, 0, 6*sizeof(cdouble));
  for(no=0, O=G->Objects[0]; no<G->NumObjects; O=G->Objects[++no])
   { 
     /******************************************************************/
     /* figure out the sign of the contribution of currents on this    */
     /* object's surface to the field at the evaluation point.         */
     /* note that this code does the correct thing if ObjectInQuestion */
     /* is NULL (i.e. evaluation point in external medium)             */
     /*****************************************************************/
     if (O==ObjectInQuestion)
      Sign=-1.0;
     else if (O->ContainingObject==ObjectInQuestion)
      Sign=+1.0;
     else
      continue; // in this case O does not contribute to field at eval pt

     Type=O->MP->Type;
     Offset=G->BFIndexOffset[O->Index];

     /***************************************************************/
     /* now loop over panels on object's surface to get             */
     /* contributions to field at evaluation point.                 */
     /***************************************************************/
     for(ne=0; ne<O->NumEdges; ne++)
      { 
        nt++;
        if (nt==TD->nThread) nt=0;
        if (nt!=TD->nt) continue;

        if ( Type==MP_PEC && ZKN!=NULL )
         { KAlpha=ZKN[ Offset + ne ];
           NAlpha=0.0;
         }
        else if ( Type==MP_PEC && DKN!=NULL )
         { KAlpha=DKN[ Offset + ne ];
           NAlpha=0.0; 
         }
        else if ( Type!=MP_PEC && ZKN!=NULL )
         { KAlpha=ZKN[ Offset + 2*ne ];
           NAlpha=ZKN[ Offset + 2*ne + 1 ];
         }
        else if ( Type!=MP_PEC && DKN!=NULL )
         { KAlpha=DKN[ Offset + 2*ne ];
           NAlpha=DKN[ Offset + 2*ne + 1 ];
         };

        KAlpha*=Sign;
        NAlpha*=Sign;
      
        O->GetReducedPotentials(ne, X, K, a, Curla, Gradp);

        for(i=0; i<3; i++)
         { EH[i]   += ZVAC*( KAlpha*(iwu*a[i] - Gradp[i]/iwe) + NAlpha*Curla[i] );
           EH[i+3] += -1.0*NAlpha*(iwe*a[i] - Gradp[i]/iwu) + KAlpha*Curla[i];
         };

      }; // for (ne=0 ... 

    }; // for(no=0 ... 

  return 0;
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void *GetFields_Thread(ThreadData *TD)
{
  ThreadData *TD=(ThreadData *)data;

#ifdef USE_PTHREAD
  SetCPUAffinity(TD->nt);
#endif

  /***************************************************************/
  /* fields unpacked from thread data structure ******************/
  /***************************************************************/
  RWGGeometry *G              = TD->G;
  const double *X             = TD->X;
  cdouble Omega               = TD->Omega;
  cdouble Eps                 = TD->Eps;
  cdouble Mu                  = TD->Mu;
  RWGObject *ObjectInQuestion = TD->ObjectInQuestion;
  HVector *KN                 = TD->KN;
  cdouble *EH                 = TD->EH;

  /***************************************************************/
  /* other local fields ******************************************/
  /***************************************************************/
  double *DKN;
  cdouble *ZKN;
  cdouble iwe=II*Omega*Eps;
  cdouble iwu=II*Omega*Mu;
  cdouble K=csqrt2(Eps*Mu)*Omega;
  if ( KN->RealComplex==LHM_REAL )
   { DKN=KN->DV;
     ZKN=NULL;
   }
  else
   { 
     DKN=NULL;
     ZKN=KN->ZV;
   };

  /***************************************************************/
  /* loop over all eval points (all rows of the XMatrix)         */
  /***************************************************************/
  for(nr=0; nr<XMatrix->NR; nr++) 
   { 
     X[0]=
    
     /*--------------------------------------------------------------*/
     /*- get scattered fields at eval point -------------------------*/
     /*--------------------------------------------------------------*/
     for(no=0;  

     /*--------------------------------------------------------------*/
     /*- add incident fields   --------------------------------------*/
     /*--------------------------------------------------------------*/

     /*--------------------------------------------------------------*/
     /*- compute field functions ------------------------------------*/
     /*--------------------------------------------------------------*/

   }; // for (nr=0; nr<XMatrix->NR; nr++)

} 


/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::GetFields(const IncField *IF, const HVector *KN, 
                            const cdouble Omega, const HMatrix *XMatrix, 
                            const char *Functions, HMatrix *FMatrix, 
                            int nThread)
{ 
  if (nThread <= 0) nThread = GetNumThreads();
  
  /***************************************************************/
  /* fire off threads                                            */
  /***************************************************************/
// ''modern'' attempt at multithreading begins here  

  int nt, nc;

  // set up an instance of ThreadData containing all fields
  // that are common to all threads, which we can subsequently
  // copy wholesale to initialize new ThreadData structures
  ThreadData ReferenceTD; 
  ReferenceTD.G=this;
  ReferenceTD.X=X;
  ReferenceTD.Omega=Omega;
  ReferenceTD.Eps=Eps;
  ReferenceTD.Mu=Mu;
  ReferenceTD.ObjectInQuestion=ObjectInQuestion;
  ReferenceTD.KN=KN;

  memset(EH, 0, 6*sizeof(cdouble));

#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[nThread], *TD;
  pthread_t *Threads = new pthread_t[nThread];

  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDs[nt]);
     memcpy(TD, &ReferenceTD, sizeof(ThreadData));
     TD->nt=nt;
     TD->nThread=nThread;

     if (nt+1 == nThread)
       GetFields_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, GetFields_Thread, (void *)TD);
   }
  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);

  /***************************************************************/
  /* sum contributions from all threads                          */
  /***************************************************************/
  for(nt=0; nt<nThread; nt++)
   for(nc=0; nc<6; nc++)
    EH[nc]+=TDs[nt].EH[nc];

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
     memcpy(&TD1, &ReferenceTD, sizeof(ThreadData));
     TD1.nt=nt;
     TD1.nThread=nTask;
     GetFields_Thread((void *)&TD1);
     for(nc=0; nc<6; nc++)
      EH[nc]+=TD1.EH[nc];
   };
#endif

// end of ''modern'' attempt at multithreading

}

/***************************************************************/
/* simple (old) interface to GetFields *************************/
/***************************************************************/
void GetFields(const IncField *IF, const HVector *KN, 
               const cdouble Omega, const double *X, 
               cdouble *EH, int nThread)
{
  HMatrix XMatrix(1, 3, LHM_REAL, LHM_NORMAL, (void *)X);

  HMatrix FMatrix(1, 6, LHM_COMPLEX, LHM_NORMAL, (void *)EH);

  GetFields(IF, KN, Omega, XMatrix, FMatrix, nThread);
} 


} // namespace scuff
