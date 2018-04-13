/*
 * Assemble.cc  -- libTDRT routine for assembling the T and U blocks 
 *              -- of the BEM matrix
 * 
 * homer reid   -- 10/2006 -- 10/2010
 * 
 * note: the dimensions of the vector space in which
 *       the BEM matrix operates are organized as follows:
 *
 * EFIE: 
 *  V_{2m}   =  electric currents in f_m^\parallel
 *  V_{2m+1} =  electric currents in f_m^z
 *
 * PMCHW:
 *  V_{4m}   =  electric currents in f_m^\parallel
 *  V_{4m+1} =  magnetic currents in f_m^\parallel
 *  V_{4m+2} =  electric currents in f_m^z
 *  V_{4m+3} =  magnetic currents in f_m^z
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhmat.h>
#include <libMatProp.h>
#include <libhrutil.h>

#include "libTDRT.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#define II cdouble(0,1)

#define FORMULATION_EFIE  0
#define FORMULATION_PMCHW 1

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   /* these fields only used by AssembleU */
   TDRTObject *Oa, *Ob;
   HMatrix *U, *dUdX, *dUdY;

   /* these fields only used by AssembleT */
   TDRTObject *O;
   HMatrix *T;

   /* these fields used by both */
   int nt, NumTasks;
   double Xi, q;
   double EpsOut, MuOut, EpsIn, MuIn;
   StaticSSIDataTable *SSSIDT;

 } ThreadData;

/***************************************************************/
/* AssembleU_Thread          *********************************/
/***************************************************************/
static void *AssembleU_Thread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  /* local copies of fields in ThreadData structure */
  TDRTObject *Oa       = TD->Oa;
  TDRTObject *Ob       = TD->Ob;
  HMatrix *U           = TD->U;
  HMatrix *dUdX        = TD->dUdX;
  HMatrix *dUdY        = TD->dUdY;
  double Xi            = TD->Xi;
  double q             = TD->q;
  double EpsOut        = TD->EpsOut;
  double MuOut         = TD->MuOut;
  double EpsIn         = TD->EpsIn;

  StaticSSIDataTable *SSSIDT = TD->SSSIDT;

  /* other local variables */
  int niva, NIVa, nivb, NIVb, nt, IndexA, IndexB;
  double Kappa, Kappa2, Z, KZ, KoZ;
  int Formulation;
  LFBuffer LBuf, *L=&LBuf;
  LFBuffer dLdXBuf, *dLdX=&dLdXBuf;
  LFBuffer dLdYBuf, *dLdY=&dLdYBuf;

  NIVa=Oa->NumIVs;
  NIVb=Ob->NumIVs;
  Kappa2=EpsOut*MuOut*Xi*Xi;
  Kappa=sqrt(Kappa2);
  Z=sqrt(MuOut/EpsOut);
  KZ=Kappa*Z;
  KoZ=Kappa/Z;

  /* figure out whether we are doing EFIE or PMCHW.      */
  /* note that EpsIn==-1.0 is code for EFIE, but that    */
  /* EpsIn, MuIn are otherwise not used                  */
  if ( EpsIn==-1.0 )
   Formulation=FORMULATION_EFIE;
  else  
   Formulation=FORMULATION_PMCHW;

  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/
#ifdef USE_PTHREAD
#ifdef _GNU_SOURCE
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(TD->nt,&cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
#endif
  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/

  /****************************************************************/ 
  /*- loop over all control points (interior vertices) on both   -*/
  /*- objects; for each pair of control points, compute the      -*/
  /*- L-functions for all basis functions associated with the    -*/
  /*- pair of control points and stamp them into their proper    -*/
  /*- slots in the U matrix                                      -*/
  /****************************************************************/  
  nt=0;
  for(niva=0; niva<NIVa; niva++)
   for(nivb=0; nivb<NIVb; nivb++)
    { 
      nt++;
      if (nt==TD->NumTasks) nt=0;
      if (nt!=TD->nt) continue;

      /* get the L functions for all basis functions associated with */
      /* this pair of control points                                 */
      ComputeLFunctions(Oa, niva, Ob, nivb, Kappa, q, SSSIDT, L, dLdX, dLdY);

      /* stamp the L-functions into their appropriate places in the */
      /* U matrix                                                   */
      if ( Formulation == FORMULATION_EFIE )
       { 
          IndexA=2*niva;  IndexB=2*nivb;

          U->SetEntry( IndexA,   IndexB,   -KZ*(L->LPPBullet + L->LPPNabla/Kappa2) );
          U->SetEntry( IndexA,   IndexB+1, -KZ*(L->LPZBullet + L->LPZNabla/Kappa2) );
          U->SetEntry( IndexA+1, IndexB,   -KZ*(L->LZPBullet + L->LZPNabla/Kappa2) );
          U->SetEntry( IndexA+1, IndexB+1, -KZ*(L->LZZBullet + L->LZZNabla/Kappa2) );

          if (dUdX)
           { dUdX->SetEntry( IndexA,   IndexB,   -KZ*(dLdX->LPPBullet + dLdX->LPPNabla/Kappa2) );
             dUdX->SetEntry( IndexA,   IndexB+1, -KZ*(dLdX->LPZBullet + dLdX->LPZNabla/Kappa2) );
             dUdX->SetEntry( IndexA+1, IndexB,   -KZ*(dLdX->LZPBullet + dLdX->LZPNabla/Kappa2) );
             dUdX->SetEntry( IndexA+1, IndexB+1, -KZ*(dLdX->LZZBullet + dLdX->LZZNabla/Kappa2) );
           };

          if (dUdY)
           { dUdY->SetEntry( IndexA,   IndexB,   -KZ*(dLdY->LPPBullet + dLdY->LPPNabla/Kappa2) );
             dUdY->SetEntry( IndexA,   IndexB+1, -KZ*(dLdY->LPZBullet + dLdY->LPZNabla/Kappa2) );
             dUdY->SetEntry( IndexA+1, IndexB,   -KZ*(dLdY->LZPBullet + dLdY->LZPNabla/Kappa2) );
             dUdY->SetEntry( IndexA+1, IndexB+1, -KZ*(dLdY->LZZBullet + dLdY->LZZNabla/Kappa2) );
           };

       }
      else /* PMCHW */
       { 
          IndexA=4*niva;  IndexB=4*nivb;
          U->SetEntry( IndexA,   IndexB,    -KZ*(L->LPPBullet + L->LPPNabla/Kappa2) );
          U->SetEntry( IndexA,   IndexB+1,  -II*L->LPPTimes );
          U->SetEntry( IndexA+1, IndexB,    -II*L->LPPTimes );
          U->SetEntry( IndexA+1, IndexB+1,  KoZ*(L->LPPBullet + L->LPPNabla/Kappa2) );

          IndexA=4*niva;  IndexB=4*nivb+2;
          U->SetEntry( IndexA,   IndexB,    -KZ*(L->LPZBullet + L->LPZNabla/Kappa2) );
          U->SetEntry( IndexA,   IndexB+1,  -II*L->LPZTimes );
          U->SetEntry( IndexA+1, IndexB,    -II*L->LPZTimes );
          U->SetEntry( IndexA+1, IndexB+1,  KoZ*(L->LPZBullet + L->LPZNabla/Kappa2) );

          IndexA=4*niva+2; IndexB=4*nivb;
          U->SetEntry( IndexA,   IndexB,    -KZ*(L->LZPBullet + L->LZPNabla/Kappa2) );
          U->SetEntry( IndexA,   IndexB+1,  -II*L->LZPTimes );
          U->SetEntry( IndexA+1, IndexB,    -II*L->LZPTimes );
          U->SetEntry( IndexA+1, IndexB+1,  KoZ*(L->LZPBullet + L->LZPNabla/Kappa2) );

          IndexA=4*niva+2; IndexB=4*nivb+2;
          U->SetEntry( IndexA,   IndexB,    -KZ*(L->LZZBullet + L->LZZNabla/Kappa2) );
          U->SetEntry( IndexA,   IndexB+1,  -L->LZZTimes );
          U->SetEntry( IndexA+1, IndexB,    -L->LZZTimes );
          U->SetEntry( IndexA+1, IndexB+1,  KoZ*(L->LZZBullet + L->LZZNabla/Kappa2) );

          if (dUdX)
           { 
             IndexA=4*niva;  IndexB=4*nivb;
             dUdX->SetEntry( IndexA,   IndexB,    -KZ*(dLdX->LPPBullet + dLdX->LPPNabla/Kappa2) );
             dUdX->SetEntry( IndexA,   IndexB+1,  -II*dLdX->LPPTimes );
             dUdX->SetEntry( IndexA+1, IndexB,    -II*dLdX->LPPTimes );
             dUdX->SetEntry( IndexA+1, IndexB+1,  KoZ*(dLdX->LPPBullet + dLdX->LPPNabla/Kappa2) );
   
             IndexA=4*niva;  IndexB=4*nivb+2;
             dUdX->SetEntry( IndexA,   IndexB,    -KZ*(dLdX->LPZBullet + dLdX->LPZNabla/Kappa2) );
             dUdX->SetEntry( IndexA,   IndexB+1,  -II*dLdX->LPZTimes );
             dUdX->SetEntry( IndexA+1, IndexB,    -II*dLdX->LPZTimes );
             dUdX->SetEntry( IndexA+1, IndexB+1,  KoZ*(dLdX->LPZBullet + dLdX->LPZNabla/Kappa2) );
   
             IndexA=4*niva+2; IndexB=4*nivb;
             dUdX->SetEntry( IndexA,   IndexB,    -KZ*(dLdX->LZPBullet + dLdX->LZPNabla/Kappa2) );
             dUdX->SetEntry( IndexA,   IndexB+1,  -II*dLdX->LZPTimes );
             dUdX->SetEntry( IndexA+1, IndexB,    -II*dLdX->LZPTimes );
             dUdX->SetEntry( IndexA+1, IndexB+1,  KoZ*(dLdX->LZPBullet + dLdX->LZPNabla/Kappa2) );
   
             IndexA=4*niva+2; IndexB=4*nivb+2;
             dUdX->SetEntry( IndexA,   IndexB,    -KZ*(dLdX->LZZBullet + dLdX->LZZNabla/Kappa2) );
             dUdX->SetEntry( IndexA,   IndexB+1,  -dLdX->LZZTimes );
             dUdX->SetEntry( IndexA+1, IndexB,    -dLdX->LZZTimes );
             dUdX->SetEntry( IndexA+1, IndexB+1,  KoZ*(dLdX->LZZBullet + dLdX->LZZNabla/Kappa2) );
           }; // if (dUdX)

          if (dUdY)
           { 
             IndexA=4*niva;  IndexB=4*nivb;
             dUdY->SetEntry( IndexA,   IndexB,    -KZ*(dLdY->LPPBullet + dLdY->LPPNabla/Kappa2) );
             dUdY->SetEntry( IndexA,   IndexB+1,  -II*dLdY->LPPTimes );
             dUdY->SetEntry( IndexA+1, IndexB,    -II*dLdY->LPPTimes );
             dUdY->SetEntry( IndexA+1, IndexB+1,  KoZ*(dLdY->LPPBullet + dLdY->LPPNabla/Kappa2) );
   
             IndexA=4*niva;  IndexB=4*nivb+2;
             dUdY->SetEntry( IndexA,   IndexB,    -KZ*(dLdY->LPZBullet + dLdY->LPZNabla/Kappa2) );
             dUdY->SetEntry( IndexA,   IndexB+1,  -II*dLdY->LPZTimes );
             dUdY->SetEntry( IndexA+1, IndexB,    -II*dLdY->LPZTimes );
             dUdY->SetEntry( IndexA+1, IndexB+1,  KoZ*(dLdY->LPZBullet + dLdY->LPZNabla/Kappa2) );
    
             IndexA=4*niva+2; IndexB=4*nivb;
             dUdY->SetEntry( IndexA,   IndexB,    -KZ*(dLdY->LZPBullet + dLdY->LZPNabla/Kappa2) );
             dUdY->SetEntry( IndexA,   IndexB+1,  -II*dLdY->LZPTimes );
             dUdY->SetEntry( IndexA+1, IndexB,    -II*dLdY->LZPTimes );
             dUdY->SetEntry( IndexA+1, IndexB+1,  KoZ*(dLdY->LZPBullet + dLdY->LZPNabla/Kappa2) );
   
             IndexA=4*niva+2; IndexB=4*nivb+2;
             dUdY->SetEntry( IndexA,   IndexB,    -KZ*(dLdY->LZZBullet + dLdY->LZZNabla/Kappa2) );
             dUdY->SetEntry( IndexA,   IndexB+1,  -dLdY->LZZTimes );
             dUdY->SetEntry( IndexA+1, IndexB,    -dLdY->LZZTimes );
             dUdY->SetEntry( IndexA+1, IndexB+1,  KoZ*(dLdY->LZZBullet + dLdY->LZZNabla/Kappa2) );
           }; // if (dUdY)

       }; // if EFIE ... else PMCHW ... 

    }; // for(niva=0; ... for (nivb=0; ... 

  return 0;

}
 
/***************************************************************/  
/***************************************************************/  
/***************************************************************/
void AssembleU(TDRTObject *Oa, TDRTObject *Ob, double Xi, double q,
               double EpsOut, double MuOut, int NumThreads,
               StaticSSIDataTable *SSSIDT,
               HMatrix *U, HMatrix *dUdX, HMatrix *dUdY)
{ 
  /* for the U matrices we don't need the values of epsilon and mu */
  /* inside the objects, but we do need to know if they are PEC or */
  /* not; we communicate this information to the _Thread routine   */
  /* by setting EpsIn=-1.0 for the PEC case and =0.0 otherwise     */
  double EpsIn;
  if ( Oa->MP->IsPEC() && Ob->MP->IsPEC() )
   EpsIn=-1.0;
  else if ( !(Oa->MP->IsPEC()) && !(Ob->MP->IsPEC()) )
   EpsIn=0.0;  
  else
   ErrExit("mixed PEC / dielectric geometries not supported");

  Log("Assembling U(%s,%s) at (Xi,q)=(%e,%e)",Oa->Label,Ob->Label,Xi,q);
#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[NumThreads], *TD;
  pthread_t *Threads = new pthread_t[NumThreads];
  Log(" POSIX multithreading (%i threads)",NumThreads);
  for(int nt=0; nt<NumThreads; nt++)
   { 
     TD=&(TDs[nt]);
     TD->nt=nt;
     TD->NumTasks=NumThreads;

     TD->SSSIDT=SSSIDT;
     TD->Oa=Oa;
     TD->Ob=Ob;
     TD->Xi=Xi;
     TD->q=q;
     TD->EpsOut=EpsOut;
     TD->MuOut=MuOut;
     TD->EpsIn=EpsIn;
     TD->U=U;
     TD->dUdX=dUdX;
     TD->dUdY=dUdY;

     if (nt+1 == NumThreads)
       AssembleU_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, AssembleU_Thread, (void *)TD);
   }

  for(int nt=0; nt<NumThreads-1; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDs;

#else 
#ifndef USE_OPENMP
  int NumTasks=NumThreads=1;
#else
  int NumTasks=NumThreads*100;
  Log(" OpenMP multithreading (%i/%i threads/tasks)",NumThreads/NumTasks);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nt=0; nt<NumTasks; nt++)
   { 
     ThreadData TD1;
     TD1.nt=nt;
     TD1.NumTasks=NumTasks;

     TD1.SSSIDT=SSSIDT;
     TD1.Oa=Oa;
     TD1.Ob=Ob;
     TD1.Xi=Xi;
     TD1.q=q;
     TD1.EpsOut=EpsOut;
     TD1.MuOut=MuOut;
     TD1.EpsIn=EpsIn;
     TD1.U=U;
     TD1.dUdX=dUdX;
     TD1.dUdY=dUdY;

     AssembleU_Thread((void *)&TD1);
   };
#endif

}



/***************************************************************/
/* NOTE: because the T matrix is symmetric, we only fill in    */
/* its upper triangle, and it is stored in packed format       */ 
/* inside the HMatrix class. also, in this routine, we compute */
/* the matrix entries in 2x2 blocks at a time (corresponding   */
/* to a single pair of basis functions.) when we go to stamp   */
/* a given 2x2 block of entries into the matrix, there is a    */
/* subtlety that affects blocks on the diagonal.               */
/* consider for example what happens when we go to stamp in    */
/* the 2x2 block corresponding to the interactions of basis    */
/* function #7 with itself. this block corresponds to elements */
/*  (14,14) (14,15)                                            */
/*  (15,14) (15,15)                                            */
/* of the larger matrix. but notice that (15,14) lies in the   */
/* lower triangle of the matrix. hence, if we try to set this  */
/* entry to a value, that value will override the value that   */
/* is already set for (14,15) (because of the way HMatrix is   */
/* implemented.) hence, whenever we go to stamp in a 2x2 block */
/* corresponding to ne==nep in the loop below, we need NOT to  */
/* stamp in the (2*ne+1, 2*nep) element.                       */
/***************************************************************/
/***************************************************************/
/* AssembleT_Thread                                            */
/***************************************************************/
static void *AssembleT_Thread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  /* local copies of fields in ThreadData structure */
  TDRTObject *O       = TD->O;
  double Xi           = TD->Xi;
  double q            = TD->q;
  double EpsOut       = TD->EpsOut;
  double MuOut        = TD->MuOut;
  double EpsIn        = TD->EpsIn;
  double MuIn         = TD->MuIn;
  HMatrix *T          = TD->T;

  StaticSSIDataTable *SSSIDT = TD->SSSIDT;

  /* other local variables */
  int niva, nivb, NIV, nt, IndexA, IndexB;
  int Formulation;

  double KappaOut2, KappaOut, ZOut, KZOut, KoZOut;
  double KappaIn2=0.0, KappaIn=0.0, ZIn=0.0, KZIn=0.0, KoZIn=0.0;
  LFBuffer LBuf, *L=&LBuf;

  NIV=O->NumIVs;

  KappaOut2=EpsOut*MuOut*Xi*Xi;
  KappaOut=sqrt(KappaOut2);
  ZOut=sqrt(MuOut/EpsOut);
  KZOut=KappaOut*ZOut;
  KoZOut=KappaOut/ZOut;

  /* figure out what formulation we are using.      */ 
  /* (EpsIn=-1.0 is code for EFIE; otherwise PMCHW) */
  if ( EpsIn==-1.0 )
   Formulation=FORMULATION_EFIE;
  else
   { Formulation=FORMULATION_PMCHW; 
     KappaIn2=EpsIn*MuIn*Xi*Xi;
     KappaIn=sqrt(KappaIn2);
     ZIn=sqrt(MuIn/EpsIn);
     KZIn=KappaIn*ZIn;
     KoZIn=KappaIn/ZIn;
   };

  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/
#ifdef USE_PTHREAD
#ifdef _GNU_SOURCE
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(TD->nt,&cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
#endif
  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/

  nt=0;
  for(niva=0; niva<NIV; niva++)
   for(nivb=niva; nivb<NIV; nivb++)
    { 
      nt++;
      if (nt==TD->NumTasks) nt=0;
      if (nt!=TD->nt) continue;

  
      if ( Formulation == FORMULATION_EFIE )
       { 
          ComputeLFunctions(O, niva, O, nivb, KappaOut, q, SSSIDT, L, 0, 0);

          IndexA=2*niva;
          IndexB=2*nivb;

          T->SetEntry( IndexA,   IndexB,   -KZOut*(L->LPPBullet + L->LPPNabla/KappaOut2) );
          T->SetEntry( IndexA,   IndexB+1, -KZOut*(L->LPZBullet + L->LPZNabla/KappaOut2) );
          if (nivb!=niva)
           T->SetEntry( IndexA+1, IndexB,  -KZOut*(L->LZPBullet + L->LZPNabla/KappaOut2) );
          T->SetEntry( IndexA+1, IndexB+1, -KZOut*(L->LZZBullet + L->LZZNabla/KappaOut2) );

       }
      else
       { 
          ComputeLFunctions(O, niva, O, nivb, KappaOut, q, SSSIDT, L, 0, 0);

          IndexA=4*niva; IndexB=4*nivb;
          T->SetEntry( IndexA,   IndexB,   -KZOut*(L->LPPBullet + L->LPPNabla/KappaOut2) );
          T->SetEntry( IndexA,   IndexB+1, -II*L->LPPTimes);
          if (nivb!=niva)
           T->SetEntry( IndexA+1, IndexB,  -II*L->LPPTimes);
          T->SetEntry( IndexA+1, IndexB+1, KoZOut*(L->LPPBullet + L->LPPNabla/KappaOut2) );

          IndexA=4*niva; IndexB=4*nivb+2;
          T->SetEntry( IndexA,   IndexB,   -KZOut*(L->LPZBullet + L->LPZNabla/KappaOut2) );
          T->SetEntry( IndexA,   IndexB+1, -II*L->LPZTimes);
          T->SetEntry( IndexA+1, IndexB,   -II*L->LPZTimes);
          T->SetEntry( IndexA+1, IndexB+1, KoZOut*(L->LPZBullet + L->LPZNabla/KappaOut2) );

          if (nivb!=niva)
           { IndexA=4*niva+2; IndexB=4*nivb;
             T->SetEntry( IndexA,   IndexB,   -KZOut*(L->LZPBullet + L->LZPNabla/KappaOut2) );
             T->SetEntry( IndexA,   IndexB+1, -II*L->LZPTimes);
             T->SetEntry( IndexA+1, IndexB,   -II*L->LZPTimes);
             T->SetEntry( IndexA+1, IndexB+1, KoZOut*(L->LZPBullet + L->LZPNabla/KappaOut2) );
           };

          IndexA=4*niva+2; IndexB=4*nivb+2;
          T->SetEntry( IndexA,   IndexB,   -KZOut*(L->LZZBullet + L->LZZNabla/KappaOut2) );
          T->SetEntry( IndexA,   IndexB+1, -L->LZZTimes);
          if (nivb!=niva)
           T->SetEntry( IndexA+1, IndexB,  -L->LZZTimes);
          T->SetEntry( IndexA+1, IndexB+1, KoZOut*(L->LZZBullet + L->LZZNabla/KappaOut2) );

          ComputeLFunctions(O, niva, O, nivb, KappaIn, q, SSSIDT, L, 0, 0);

          IndexA=4*niva; IndexB=4*nivb;
          T->AddEntry( IndexA,   IndexB,   -KZIn*(L->LPPBullet + L->LPPNabla/KappaIn2) );
          T->AddEntry( IndexA,   IndexB+1, -II*L->LPPTimes);
          if (nivb!=niva)
           T->AddEntry( IndexA+1, IndexB,  -II*L->LPPTimes);
          T->AddEntry( IndexA+1, IndexB+1, KoZIn*(L->LPPBullet + L->LPPNabla/KappaIn2) );

          IndexA=4*niva; IndexB=4*nivb+2;
          T->AddEntry( IndexA,   IndexB,   -KZIn*(L->LPZBullet + L->LPZNabla/KappaIn2) );
          T->AddEntry( IndexA,   IndexB+1, -II*L->LPZTimes);
          T->AddEntry( IndexA+1, IndexB,   -II*L->LPZTimes);
          T->AddEntry( IndexA+1, IndexB+1, KoZIn*(L->LPZBullet + L->LPZNabla/KappaIn2) );

          if (nivb!=niva)
           { IndexA=4*niva+2; IndexB=4*nivb;
             T->AddEntry( IndexA,   IndexB,   -KZIn*(L->LZPBullet + L->LZPNabla/KappaIn2) );
             T->AddEntry( IndexA,   IndexB+1, -II*L->LZPTimes);
             T->AddEntry( IndexA+1, IndexB,   -II*L->LZPTimes);
             T->AddEntry( IndexA+1, IndexB+1, KoZIn*(L->LZPBullet + L->LZPNabla/KappaIn2) );
           };

          IndexA=4*niva+2; IndexB=4*nivb+2;
          T->AddEntry( IndexA,   IndexB,   -KZIn*(L->LZZBullet + L->LZZNabla/KappaIn2) );
          T->AddEntry( IndexA,   IndexB+1, -L->LZZTimes);
          if (nivb!=niva)
           T->AddEntry( IndexA+1, IndexB,  -L->LZZTimes);
          T->AddEntry( IndexA+1, IndexB+1, KoZIn*(L->LZZBullet + L->LZZNabla/KappaIn2) );

       }; //  if ( Formulation == FORMULATION_EFIE ) ... else ... 

    }; // for (niva = ... ) for (nivb = ... )

  return 0;

}
 
/***************************************************************/  
/***************************************************************/  
/***************************************************************/
void AssembleT(TDRTObject *O, double Xi, double q,
               double EpsOut, double MuOut, int NumThreads, 
               StaticSSIDataTable *SSSIDT, HMatrix *T)
{ 
  double EpsIn, MuIn;
  if ( O->MP->IsPEC() )
   EpsIn=-1.0;
  else
   { cdouble zEps, zMu;
     O->MP->GetEpsMu( cdouble(0,Xi), &zEps, &zMu);
     EpsIn = real(zEps);
     MuIn  = real(zMu);
   };
  Log("Assembling T(%s) at (Xi,q)=(%e,%e)",O->Label,Xi,q);

#ifdef USE_PTHREAD
  Log(" POSIX multithreading (%i threads)",NumThreads);
  ThreadData *TDs = new ThreadData[NumThreads], *TD;
  pthread_t *Threads = new pthread_t[NumThreads];
  for(int nt=0; nt<NumThreads; nt++)
   { 
     TD=&(TDs[nt]);
     TD->nt=nt;
     TD->NumTasks=NumThreads;

     TD->SSSIDT=SSSIDT;
     TD->O=O;
     TD->Xi=Xi;
     TD->q=q;
     TD->EpsOut=EpsOut;
     TD->MuOut=MuOut;
     TD->EpsIn=EpsIn;
     TD->MuIn=MuIn;
     TD->T=T;

     if (nt+1 == NumThreads)
       AssembleT_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, AssembleT_Thread, (void *)TD);
   }

  for(int nt=0; nt<NumThreads-1; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDs;

#else 
#ifndef USE_OPENMP
  int NumTasks=NumThreads=1;
#else
  int NumTasks=NumThreads*100;
  Log(" OpenMP multithreading (%i/%i threads/tasks)",NumThreads/NumTasks);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nt=0; nt<NumTasks; nt++)
   { 
     ThreadData TD1;
     TD1.nt=nt;
     TD1.NumTasks=NumTasks;

     TD1.SSSIDT=SSSIDT;
     TD1.O=O;
     TD1.Xi=Xi;
     TD1.q=q;
     TD1.EpsOut=EpsOut;
     TD1.MuOut=MuOut;
     TD1.EpsIn=EpsIn;
     TD1.MuIn=MuIn;
     TD1.T=T;

     AssembleT_Thread((void *)&TD1);
   };
#endif

}
