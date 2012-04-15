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
/* integrand passed to TriInt to compute contributions to electric */
/* and magnetic fields from a single panel                         */
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
typedef struct EHFieldIntegrandData
 { 
   double *Q;           // RWG basis function source/sink vertex 
   double PreFac;       // RWG basis function prefactor 
   const double *X0;    // field evaluation point 
   cdouble K;           // \sqrt{Eps*Mu} * frequency
 } EHFIData;

static void EHFieldIntegrand(double *X, void *parms, double *f)
{ 
  EHFieldIntegrandData *EHFID=(EHFieldIntegrandData *)parms;
  double fRWG[3], XmX0[3], fxR[3];
  cdouble Phi, Psi, *L=(cdouble *)f;
  double r;
  cdouble K;
  double PreFac=EHFID->PreFac;

  /* get the value of the RWG basis function at XP */
  VecSub(X,EHFID->Q,fRWG);
  VecScale(fRWG,PreFac);
  
  VecSub(X,EHFID->X0,XmX0);
  VecCross(fRWG,XmX0,fxR);
  r=VecNorm(XmX0);

  K=EHFID->K;
  
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
/***************************************************************/
/***************************************************************/
void RWGObject::GetReducedPotentials(int ne, const double *X, cdouble K,
                                     cdouble *a, cdouble *Curla,
                                     cdouble *Gradp)
{
  double *QP, *V1, *V2, *QM;
  double PArea, MArea;
  RWGEdge *E;
  int mu;
  EHFieldIntegrandData MyEHFIData, *EHFID=&MyEHFIData;
  cdouble IP[9], IM[9];

  /* get edge vertices */
  E=Edges[ne];
  QP=Vertices + 3*(E->iQP);
  V1=Vertices + 3*(E->iV1);
  V2=Vertices + 3*(E->iV2);
  QM=Vertices + 3*(E->iQM);
  PArea=Panels[E->iPPanel]->Area;
  MArea=Panels[E->iMPanel]->Area;

  /* set up data structure passed to EFieldIntegrand */
  EHFID->X0=X;
  EHFID->K=K;

  /* contribution of positive panel */
  EHFID->Q=QP;
  EHFID->PreFac = E->Length / (2.0*PArea);
  TriIntFixed(EHFieldIntegrand, 18, (void *)EHFID, QP, V1, V2, 25, (double *)IP);

  /* contribution of negative panel */
  EHFID->Q=QM;
  EHFID->PreFac = E->Length / (2.0*MArea);
  TriIntFixed(EHFieldIntegrand, 18, (void *)EHFID, V1, V2, QM, 25, (double *)IM);

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
/* Get scattered fields at point X.                            */
/* If ObjectIndex = -1, X is assumed to lie in the external    */
/* region. Otherwise, X is assumed to lie in the interior of   */
/* object #ObjectIndex.                                        */
/***************************************************************/
void RWGGeometry::GetFields(const double X[3], int ObjectIndex, cdouble Omega,
                            HVector *KN, cdouble EH[6], int nThread)
{ 
  if (nThread <= 0) nThread = GetNumThreads();
  
  /***************************************************************/
  /* switch off to determine whether we are in the external      */
  /* medium, inside an object, or otherwise                      */
  /***************************************************************/
  RWGObject *ObjectInQuestion;
  cdouble Eps, Mu; 

  if ( ObjectIndex==-1 )                                /* in external medium */
   { ObjectInQuestion=0;
     ExteriorMP->GetEpsMu(Omega, &Eps, &Mu);
   }
  else if ( ObjectIndex<-1 || ObjectIndex>=NumObjects ) /* invalid object */
   { fprintf(stderr,"\n*\n* WARNING: invalid object selected in GetFields\n*\n");
     memset(EH,0,6*sizeof(cdouble));
     return;
   }
  else if ( Objects[ObjectIndex]->MP->IsPEC() )         /* inside a PEC object */
   { memset(EH,0,6*sizeof(cdouble));      /* fields vanish in a PEC body */
     return;
   }
  else                                                  /* inside a non-PEC object*/
   { ObjectInQuestion=Objects[ObjectIndex];
     ObjectInQuestion->MP->GetEpsMu(Omega, &Eps, &Mu);
   };

  /***************************************************************/
  /* fire off threads                                            */
  /***************************************************************/

#if 0 

  this is the code as of 3/25/2012, commented out 
  but retained in the file the for the time being 
  while i experiment with new ways to do the multithreading 

#ifdef USE_PTHREAD
  pthread_t *Threads = new pthread_t[nThread];
  ThreadData *TDS = new ThreadData[nThread], *TD;
  cdouble *PartialEH = new cdouble[6*nThread]; 
#else
  ThreadData TD1;
  double EH0r=0,EH1r=0,EH2r=0,EH3r=0,EH4r=0,EH5r=0;
  double EH0i=0,EH1i=0,EH2i=0,EH3i=0,EH4i=0,EH5i=0;
#endif
  int nt;

#ifdef USE_OPENMP
#pragma omp parallel for private(TD1), schedule(static,1), reduction(+:EH0r,EH1r,EH2r,EH3r,EH4r,EH5r,EH0i,EH1i,EH2i,EH3i,EH4i,EH5i), num_threads(nThread)
#endif
  for(nt=0; nt<nThread; nt++)
   { 
#ifdef USE_PTHREAD
     TD=&(TDS[nt]);
#else
     ThreadData *TD=&TD1;
#endif
     TD->nt=nt;
     TD->nThread=nThread;

     TD->G=this;
     TD->X=X;
     TD->Omega=Omega;
     TD->Eps=Eps;
     TD->Mu=Mu;
     TD->ObjectInQuestion=ObjectInQuestion;
     TD->KN=KN;

#ifdef USE_PTHREAD
     TD->EH=PartialEH + 6*nt;
     if (nt+1 == nThread)
       GetFields_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, GetFields_Thread, (void *)TD);
#else
     cdouble PartialEH[6];
     TD->EH=PartialEH;
     GetFields_Thread((void *)TD);
     // annoyance: openmp doesn't support reductions on complex types
     EH0r += real(PartialEH[0]); EH0i += imag(PartialEH[0]);
     EH1r += real(PartialEH[1]); EH1i += imag(PartialEH[1]);
     EH2r += real(PartialEH[2]); EH2i += imag(PartialEH[2]);
     EH3r += real(PartialEH[3]); EH3i += imag(PartialEH[3]);
     EH4r += real(PartialEH[4]); EH4i += imag(PartialEH[4]);
     EH5r += real(PartialEH[5]); EH5i += imag(PartialEH[5]);
#endif
   }

#ifdef USE_PTHREAD
  /***************************************************************/
  /* wait for threads to complete                                */
  /***************************************************************/
  for(nt=0; nt<nThread-1; nt++)
      pthread_join(Threads[nt],0);
#endif

#ifdef USE_PTHREAD
  /***************************************************************/
  /* sum contributions from all threads *                        */
  /***************************************************************/
  memset(EH,0,6*sizeof(cdouble));
  for(nt=0; nt<nThread; nt++)
   { EH[0]+=PartialEH[6*nt + 0]; 
     EH[1]+=PartialEH[6*nt + 1]; 
     EH[2]+=PartialEH[6*nt + 2]; 
     EH[3]+=PartialEH[6*nt + 3]; 
     EH[4]+=PartialEH[6*nt + 4]; 
     EH[5]+=PartialEH[6*nt + 5]; 
   };

  delete[] Threads;
  delete[] TDS;
  delete[] PartialEH;
#else
  EH[0] = cdouble(EH0r,EH0i);
  EH[1] = cdouble(EH1r,EH1i);
  EH[2] = cdouble(EH2r,EH2i);
  EH[3] = cdouble(EH3r,EH3i);
  EH[4] = cdouble(EH4r,EH4i);
  EH[5] = cdouble(EH5r,EH5i);
#endif

#endif // end of block commented out on 3/25/2012

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
/* entry point to GetFields() in which the caller identifies   */
/* the object inside which the evaluation point lies by its    */
/* label (as assigned using the LABEL keyword in the .scuffgeo */
/* file) or using the keywords "EXTERIOR" or "MEDIUM" for the  */
/* exterior medium.                                            */
/***************************************************************/
void RWGGeometry::GetFields(const double X[3], 
                            const char *ObjectLabel,
                            cdouble Omega,
                            HVector *KN, cdouble EH[6], int nThread)
{
  if (!ObjectLabel)
   ErrExit("%s:%i:internal error",__FILE__,__LINE__);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int ObjectIndex;
  if ( !strcasecmp(ObjectLabel,"EXTERIOR") || !strcasecmp(ObjectLabel,"MEDIUM") )
   ObjectIndex=-1;
  else
   { for (ObjectIndex=0; ObjectIndex<NumObjects; ObjectIndex++)
      if ( !strcasecmp(ObjectLabel,Objects[ObjectIndex]->Label) )
       break;
   }

  if (ObjectIndex==NumObjects)
   ErrExit("unknown object label %s in GetFields()",ObjectLabel);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GetFields(X, ObjectIndex, Omega, KN, EH, nThread);

}

/***************************************************************/
/* entry point to GetFields() with autodetection of where the  */
/* evaluation point lies.                                      */
/***************************************************************/
void RWGGeometry::GetFields(const double X[3], cdouble Omega,
                            HVector *KN, cdouble EH[6], int nThread)
{
  GetFields(X, GetObjectIndex(X), Omega, KN, EH, nThread);
}

/***************************************************************/
/* Autodetection of object where a given point lies; assumes   */
/* objects have been topologically sorted (so that if object A */
/* contains object B, then B comes after A).                  */
/***************************************************************/

int RWGGeometry::GetObjectIndex(const double X[3]) {
  // find the innermost object containing X
  for (int i = NumObjects - 1; i >= 0; --i) // innermost to outermost order
    if (Objects[i]->Contains(X))
      return i;
  return -1; // not in any object
}

RWGObject *RWGGeometry::GetObject(const double X[3]) {
  int i = GetObjectIndex(X);
  return i < 0 ? NULL : Objects[i];
}

} // namespace scuff
