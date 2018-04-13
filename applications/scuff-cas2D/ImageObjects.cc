/*
 * ImageObjects.cc -- Casimir2D code module that handles image objects
 *                 -- for implicit treatment of a metallic ground plane
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhmat.h>
#include <libMatProp.h>
#include <libhrutil.h>

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#include "libTDRT.h"
#include "scuff-cas2D.h"

#define FORMULATION_EFIE  0
#define FORMULATION_PMCHW 1

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   /* these fields used only by AddImageContributionsToU */
   TDRTObject *Oa, *ObImage;
   HMatrix *U, *dUdX, *dUdY;

   /* these fields used only by AssembleTI */
   TDRTObject *O, *OImage;
   HMatrix *TI, *dTIdY;

   /* these fields used by both */
   int nt, NumTasks;
   double Xi, q;
   double EpsOut, MuOut, EpsIn, MuIn;
   StaticSSIDataTable *SSSIDT;

 } ThreadData;

/***************************************************************/
/* AddImageContributionToU_Thread ******************************/
/***************************************************************/
static void *AIC2U_Thread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  /* local copies of fields in ThreadData structure */
  TDRTObject *Oa       = TD->Oa;
  TDRTObject *ObImage  = TD->ObImage;
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
  NIVb=ObImage->NumIVs;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
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
  cdouble MSIGN = -1.0;
  for(niva=0; niva<NIVa; niva++)
   for(nivb=0; nivb<NIVb; nivb++)
    { 
      nt++;
      if (nt==TD->NumTasks) nt=0;
      if (nt!=TD->nt) continue;

      /* get the L functions for all basis functions associated with */
      /* this pair of control points                                 */
      ComputeLFunctions(Oa, niva, ObImage, nivb, Kappa, q, SSSIDT, L, dLdX, dLdY);

      /* add the contributions of the image object to the matrix elements.  */
      /* how it works:                                                      */
      /*  (a) contributions from electric currents on the image panel enter */
      /*      with a minus sign                                             */
      /*  (b) contributions from magnetic currents on the image panel enter */
      /*      with a plus sign .                                            */
      if ( Formulation == FORMULATION_EFIE )
       { 
          IndexA=2*niva;  IndexB=2*nivb;

          U->AddEntry( IndexA,   IndexB,   MSIGN*(-KZ*(L->LPPBullet + L->LPPNabla/Kappa2)) );
          U->AddEntry( IndexA,   IndexB+1, MSIGN*(-KZ*(L->LPZBullet + L->LPZNabla/Kappa2)) );
          U->AddEntry( IndexA+1, IndexB,   MSIGN*(-KZ*(L->LZPBullet + L->LZPNabla/Kappa2)) );
          U->AddEntry( IndexA+1, IndexB+1, MSIGN*(-KZ*(L->LZZBullet + L->LZZNabla/Kappa2)) );

          if (dUdX)
           { dUdX->AddEntry( IndexA,   IndexB,   MSIGN*(-KZ*(dLdX->LPPBullet + dLdX->LPPNabla/Kappa2)) );
             dUdX->AddEntry( IndexA,   IndexB+1, MSIGN*(-KZ*(dLdX->LPZBullet + dLdX->LPZNabla/Kappa2)) );
             dUdX->AddEntry( IndexA+1, IndexB,   MSIGN*(-KZ*(dLdX->LZPBullet + dLdX->LZPNabla/Kappa2)) );
             dUdX->AddEntry( IndexA+1, IndexB+1, MSIGN*(-KZ*(dLdX->LZZBullet + dLdX->LZZNabla/Kappa2)) );
           };

          if (dUdY)
           { dUdY->AddEntry( IndexA,   IndexB,   MSIGN*(-KZ*(dLdY->LPPBullet + dLdY->LPPNabla/Kappa2)) );
             dUdY->AddEntry( IndexA,   IndexB+1, MSIGN*(-KZ*(dLdY->LPZBullet + dLdY->LPZNabla/Kappa2)) );
             dUdY->AddEntry( IndexA+1, IndexB,   MSIGN*(-KZ*(dLdY->LZPBullet + dLdY->LZPNabla/Kappa2)) );
             dUdY->AddEntry( IndexA+1, IndexB+1, MSIGN*(-KZ*(dLdY->LZZBullet + dLdY->LZZNabla/Kappa2)) );
           };

       }
      else /* PMCHW */
       { 
          IndexA=4*niva;  IndexB=4*nivb;
          U->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(L->LPPBullet + L->LPPNabla/Kappa2)) );
          U->AddEntry( IndexA,   IndexB+1,        (-II*L->LPPTimes)                  );
          U->AddEntry( IndexA+1, IndexB,    MSIGN*(-II*L->LPPTimes)                  );
          U->AddEntry( IndexA+1, IndexB+1,        (KoZ*(L->LPPBullet + L->LPPNabla/Kappa2)) );

          IndexA=4*niva;  IndexB=4*nivb+2;
          U->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(L->LPZBullet + L->LPZNabla/Kappa2)) );
          U->AddEntry( IndexA,   IndexB+1,        (-II*L->LPZTimes)                  );
          U->AddEntry( IndexA+1, IndexB,    MSIGN*(-II*L->LPZTimes)                  );
          U->AddEntry( IndexA+1, IndexB+1,        (KoZ*(L->LPZBullet + L->LPZNabla/Kappa2)) );

          IndexA=4*niva+2; IndexB=4*nivb;
          U->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(L->LZPBullet + L->LZPNabla/Kappa2)) );
          U->AddEntry( IndexA,   IndexB+1,        (-II*L->LZPTimes)                  );
          U->AddEntry( IndexA+1, IndexB,    MSIGN*(-II*L->LZPTimes)                  );
          U->AddEntry( IndexA+1, IndexB+1,          KoZ*(L->LZPBullet + L->LZPNabla/Kappa2) );

          IndexA=4*niva+2; IndexB=4*nivb+2;
          U->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(L->LZZBullet + L->LZZNabla/Kappa2)) );
          U->AddEntry( IndexA,   IndexB+1,        (-L->LZZTimes)                            );
          U->AddEntry( IndexA+1, IndexB,    MSIGN*(-L->LZZTimes)                            );
          U->AddEntry( IndexA+1, IndexB+1,        (KoZ*(L->LZZBullet + L->LZZNabla/Kappa2)) );

          if (dUdX)
           { 
             IndexA=4*niva;  IndexB=4*nivb;
             dUdX->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(dLdX->LPPBullet + dLdX->LPPNabla/Kappa2) ));
             dUdX->AddEntry( IndexA,   IndexB+1,        (-II*dLdX->LPPTimes )                    );
             dUdX->AddEntry( IndexA+1, IndexB,    MSIGN*(-II*dLdX->LPPTimes )                    );
             dUdX->AddEntry( IndexA+1, IndexB+1,        (KoZ*(dLdX->LPPBullet + dLdX->LPPNabla/Kappa2)) );
   
             IndexA=4*niva;  IndexB=4*nivb+2;
             dUdX->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(dLdX->LPZBullet + dLdX->LPZNabla/Kappa2)) );
             dUdX->AddEntry( IndexA,   IndexB+1,        (-II*dLdX->LPZTimes )                    );
             dUdX->AddEntry( IndexA+1, IndexB,    MSIGN*(-II*dLdX->LPZTimes )                    );
             dUdX->AddEntry( IndexA+1, IndexB+1,        (KoZ*(dLdX->LPZBullet + dLdX->LPZNabla/Kappa2)) );
   
             IndexA=4*niva+2; IndexB=4*nivb;
             dUdX->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(dLdX->LZPBullet + dLdX->LZPNabla/Kappa2)) );
             dUdX->AddEntry( IndexA,   IndexB+1,        (-II*dLdX->LZPTimes )                    );
             dUdX->AddEntry( IndexA+1, IndexB,    MSIGN*(-II*dLdX->LZPTimes )                    );
             dUdX->AddEntry( IndexA+1, IndexB+1,        (KoZ*(dLdX->LZPBullet + dLdX->LZPNabla/Kappa2)) );
   
             IndexA=4*niva+2; IndexB=4*nivb+2;
             dUdX->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(dLdX->LZZBullet + dLdX->LZZNabla/Kappa2)) );
             dUdX->AddEntry( IndexA,   IndexB+1,        (-dLdX->LZZTimes )                              );
             dUdX->AddEntry( IndexA+1, IndexB,    MSIGN*(-dLdX->LZZTimes )                              );
             dUdX->AddEntry( IndexA+1, IndexB+1,        (KoZ*(dLdX->LZZBullet + dLdX->LZZNabla/Kappa2)) );
           }; // if (dUdX)

          if (dUdY)
           { 
             IndexA=4*niva;  IndexB=4*nivb;
             dUdY->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(dLdY->LPPBullet + dLdY->LPPNabla/Kappa2) ));
             dUdY->AddEntry( IndexA,   IndexB+1,        (-II*dLdY->LPPTimes )                    );
             dUdY->AddEntry( IndexA+1, IndexB,    MSIGN*(-II*dLdY->LPPTimes )                    );
             dUdY->AddEntry( IndexA+1, IndexB+1,        (KoZ*(dLdY->LPPBullet + dLdY->LPPNabla/Kappa2)) );
   
             IndexA=4*niva;  IndexB=4*nivb+2;
             dUdY->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(dLdY->LPZBullet + dLdY->LPZNabla/Kappa2)) );
             dUdY->AddEntry( IndexA,   IndexB+1,        (-II*dLdY->LPZTimes )                    );
             dUdY->AddEntry( IndexA+1, IndexB,    MSIGN*(-II*dLdY->LPZTimes )                    );
             dUdY->AddEntry( IndexA+1, IndexB+1,        (KoZ*(dLdY->LPZBullet + dLdY->LPZNabla/Kappa2)) );
   
             IndexA=4*niva+2; IndexB=4*nivb;
             dUdY->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(dLdY->LZPBullet + dLdY->LZPNabla/Kappa2)) );
             dUdY->AddEntry( IndexA,   IndexB+1,        (-II*dLdY->LZPTimes )                           );
             dUdY->AddEntry( IndexA+1, IndexB,    MSIGN*(-II*dLdY->LZPTimes )                           );
             dUdY->AddEntry( IndexA+1, IndexB+1,        (KoZ*(dLdY->LZPBullet + dLdY->LZPNabla/Kappa2)) );
   
             IndexA=4*niva+2; IndexB=4*nivb+2;
             dUdY->AddEntry( IndexA,   IndexB,    MSIGN*(-KZ*(dLdY->LZZBullet + dLdY->LZZNabla/Kappa2)) );
             dUdY->AddEntry( IndexA,   IndexB+1,        (-dLdY->LZZTimes )                              );
             dUdY->AddEntry( IndexA+1, IndexB,    MSIGN*(-dLdY->LZZTimes )                              );

           }; // if (dUdY)

       }; // if EFIE ... else PMCHW ... 

    }; // for(niva=0; ... for (nivb=0; ... 

  return 0;

}

/***************************************************************/  
/* AddImageContributionToU -- modeled after AssembleU routine  */  
/* in libTDRT                                                  */  
/***************************************************************/
void AddImageContributionToU(TDRTObject *Oa, TDRTObject *ObImage, 
                             double Xi, double q,
                             double EpsOut, double MuOut, int NumThreads, 
                             StaticSSIDataTable *SSSIDT,
                             HMatrix *U, HMatrix *dUdX, HMatrix *dUdY)
{ 
  /* for the U matrices we don't need the values of epsilon and mu */
  /* inside the objects, but we do need to know if they are PEC or */
  /* not; we communicate this information to the _Thread routine   */
  /* by setting EpsIn=-1.0 for the PEC case and =0.0 otherwise     */
  double EpsIn;
  if ( Oa->MP->IsPEC() && ObImage->MP->IsPEC() )
   EpsIn=-1.0;
  else if ( !(Oa->MP->IsPEC()) && !(ObImage->MP->IsPEC()) )
   EpsIn=0.0;  
  else
   ErrExit("mixed PEC / dielectric geometries not supported");

#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[NumThreads], *TD;
  pthread_t *Threads = new pthread_t[NumThreads];
  for(int nt=0; nt<NumThreads; nt++)
   { 
     TD=&(TDS[nt]);
     TD->nt=nt;
     TD->NumTasks=NumThreads;

     TD->SSSIDT=SSSIDT;

     TD->Oa=Oa;
     TD->ObImage=ObImage;
     TD->Xi=Xi;
     TD->q=q;
     TD->EpsOut=EpsOut;
     TD->MuOut=MuOut;
     TD->EpsIn=EpsIn;
     TD->U=U;
     TD->dUdX=dUdX;
     TD->dUdY=dUdY;

     if (nt+1 == NumThreads)
      AIC2U_Thread( (void *)TD );
     else
      pthread_create( &(Threads[nt]), 0, AIC2U_Thread, (void *)TD);
   };

  for(int nt=0; nt<NumThreads; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDs;
#else 
#ifndef USE_OPENMP
  int NumTasks=NumThreads=1;
#else
  int NumTasks=NumThreads*100;
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nt=0; nt<NumTasks; nt++)
   { 
     ThreadData TD1;
     TD1.nt=nt;
     TD1.NumTasks=NumTasks;

     TD1.SSSIDT=SSSIDT;

     TD1.Oa=Oa;
     TD1.ObImage=ObImage;
     TD1.Xi=Xi;
     TD1.q=q;
     TD1.EpsOut=EpsOut;
     TD1.MuOut=MuOut;
     TD1.EpsIn=EpsIn;
     TD1.U=U;
     TD1.dUdX=dUdX;
     TD1.dUdY=dUdY;

     AIC2U_Thread((void *)&TD1);
   };
#endif

}


/***************************************************************/
/* AssembleTI_Thread          **********************************/
/***************************************************************/
static void *AssembleTI_Thread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  /* local copies of fields in ThreadData structure */
  TDRTObject *O       = TD->O;
  TDRTObject *OImage  = TD->OImage;
  HMatrix *TI         = TD->TI;
  HMatrix *dTIdY      = TD->dTIdY;
  double Xi           = TD->Xi;
  double q            = TD->q;
  double EpsOut       = TD->EpsOut;
  double MuOut        = TD->MuOut;
  double EpsIn        = TD->EpsIn;

  StaticSSIDataTable *SSSIDT = TD->SSSIDT;

  /* other local variables */
  int nt;
  LFBuffer LBuf, *L=&LBuf;
  LFBuffer dLdYBuf, *dLdY=(dTIdY==0 ? 0 : &dLdYBuf);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Kappa2, Kappa, ZOut, KZ, KoZ;
  Kappa2=EpsOut*MuOut*Xi*Xi;
  Kappa=sqrt(Kappa2);
  ZOut=sqrt(MuOut/EpsOut);
  KZ=Kappa*ZOut;
  KoZ=Kappa/ZOut;

  /*--------------------------------------------------------------*/
  /* figure out what formulation we are using.                    */ 
  /*--------------------------------------------------------------*/
  /* (EpsIn=-1.0 is code for EFIE; otherwise PMCHW) */
  int Formulation;
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

  nt=0;
  cdouble MSIGN = -1.0;
  int niva, IndexA, nivb, IndexB, NIV=O->NumIVs;
  for(niva=0; niva<NIV; niva++)
   for(nivb=niva; nivb<NIV; nivb++)
    { 
      nt++;
      if (nt==TD->NumTasks) nt=0;
      if (nt!=TD->nt) continue;

      ComputeLFunctions(O, niva, OImage, nivb, Kappa, q, SSSIDT, L, 0, dLdY);
  
      if ( Formulation == FORMULATION_EFIE )
       { 

          IndexA=2*niva;
          IndexB=2*nivb;

          TI->SetEntry( IndexA,   IndexB,   MSIGN*(-KZ*(L->LPPBullet + L->LPPNabla/Kappa2) ) );
          TI->SetEntry( IndexA,   IndexB+1, MSIGN*(-KZ*(L->LPZBullet + L->LPZNabla/Kappa2) ) );
          if (nivb!=niva)
           TI->SetEntry( IndexA+1, IndexB,  MSIGN*(-KZ*(L->LZPBullet + L->LZPNabla/Kappa2) ) );
          TI->SetEntry( IndexA+1, IndexB+1, MSIGN*(-KZ*(L->LZZBullet + L->LZZNabla/Kappa2) ) );

          if (dTIdY)
           { dTIdY->SetEntry( IndexA,   IndexB,   MSIGN*(KZ*(dLdY->LPPBullet + dLdY->LPPNabla/Kappa2)) );
             dTIdY->SetEntry( IndexA,   IndexB+1, MSIGN*(KZ*(dLdY->LPZBullet + dLdY->LPZNabla/Kappa2)) );
             if (nivb!=niva)
              dTIdY->SetEntry( IndexA+1, IndexB,  MSIGN*(KZ*(dLdY->LZPBullet + dLdY->LZPNabla/Kappa2)) );
             dTIdY->SetEntry( IndexA+1, IndexB+1, MSIGN*(KZ*(dLdY->LZZBullet + dLdY->LZZNabla/Kappa2)) );
           };

       }
      else
       { 
          IndexA=4*niva; IndexB=4*nivb;
          TI->SetEntry( IndexA,   IndexB,   MSIGN*(-KZ*(L->LPPBullet + L->LPPNabla/Kappa2)) );
          TI->SetEntry( IndexA,   IndexB+1,       (-II*L->LPPTimes)                        ) ;
          if (nivb!=niva)
           TI->SetEntry( IndexA+1, IndexB,  MSIGN*(-II*L->LPPTimes)                        );
          TI->SetEntry( IndexA+1, IndexB+1,       (KoZ*(L->LPPBullet + L->LPPNabla/Kappa2)) );

          IndexA=4*niva; IndexB=4*nivb+2;
          TI->SetEntry( IndexA,   IndexB,   MSIGN*(-KZ*(L->LPZBullet + L->LPZNabla/Kappa2)) );
          TI->SetEntry( IndexA,   IndexB+1,       (-II*L->LPZTimes)                        ); 
          TI->SetEntry( IndexA+1, IndexB,   MSIGN*(-II*L->LPZTimes)                        );
          TI->SetEntry( IndexA+1, IndexB+1,       (KoZ*(L->LPZBullet + L->LPZNabla/Kappa2)) );

          if (nivb!=niva)
           { IndexA=4*niva+2; IndexB=4*nivb;
             TI->SetEntry( IndexA,   IndexB,   MSIGN*(-KZ*(L->LZPBullet + L->LZPNabla/Kappa2)) );
             TI->SetEntry( IndexA,   IndexB+1,       (-II*L->LZPTimes)                        );
             TI->SetEntry( IndexA+1, IndexB,   MSIGN*(-II*L->LZPTimes)                        );
             TI->SetEntry( IndexA+1, IndexB+1,       (KoZ*(L->LZPBullet + L->LZPNabla/Kappa2)) );
           };

          IndexA=4*niva+2; IndexB=4*nivb+2;
          TI->SetEntry( IndexA,   IndexB,   MSIGN*(-KZ*(L->LZZBullet + L->LZZNabla/Kappa2)) );
          TI->SetEntry( IndexA,   IndexB+1,       (-L->LZZTimes)                                  );
          if (nivb!=niva)
           TI->SetEntry( IndexA+1, IndexB,  MSIGN*(-L->LZZTimes)                                  );
          TI->SetEntry( IndexA+1, IndexB+1,       (KoZ*(L->LZZBullet + L->LZZNabla/Kappa2)) );

          if (dTIdY)
           {
             IndexA=4*niva; IndexB=4*nivb;
             dTIdY->SetEntry( IndexA,   IndexB,   MSIGN*(-KZ*(dLdY->LPPBullet + dLdY->LPPNabla/Kappa2)) );
             dTIdY->SetEntry( IndexA,   IndexB+1,       (-II*dLdY->LPPTimes)                            ) ;
             if (nivb!=niva)
              dTIdY->SetEntry( IndexA+1, IndexB,  MSIGN*(-II*dLdY->LPPTimes)                            );
             dTIdY->SetEntry( IndexA+1, IndexB+1,       (KoZ*(dLdY->LPPBullet + dLdY->LPPNabla/Kappa2)) );
   
             IndexA=4*niva; IndexB=4*nivb+2;
             dTIdY->SetEntry( IndexA,   IndexB,   MSIGN*(-KZ*(dLdY->LPZBullet + dLdY->LPZNabla/Kappa2)) );
             dTIdY->SetEntry( IndexA,   IndexB+1,       (-II*dLdY->LPZTimes)                            ); 
             dTIdY->SetEntry( IndexA+1, IndexB,   MSIGN*(-II*dLdY->LPZTimes)                            );
             dTIdY->SetEntry( IndexA+1, IndexB+1,       (KoZ*(dLdY->LPZBullet + dLdY->LPZNabla/Kappa2)) );
   
             if (nivb!=niva)
              { IndexA=4*niva+2; IndexB=4*nivb;
                dTIdY->SetEntry( IndexA,   IndexB,   MSIGN*(-KZ*(dLdY->LZPBullet + dLdY->LZPNabla/Kappa2)) );
                dTIdY->SetEntry( IndexA,   IndexB+1,       (-II*dLdY->LZPTimes)                            );
                dTIdY->SetEntry( IndexA+1, IndexB,   MSIGN*(-II*dLdY->LZPTimes)                            );
                dTIdY->SetEntry( IndexA+1, IndexB+1,       (KoZ*(dLdY->LZPBullet + dLdY->LZPNabla/Kappa2)) );
              };
   
             IndexA=4*niva+2; IndexB=4*nivb+2;
             dTIdY->SetEntry( IndexA,   IndexB,   MSIGN*(-KZ*(dLdY->LZZBullet + dLdY->LZZNabla/Kappa2)) );
             dTIdY->SetEntry( IndexA,   IndexB+1,       (-dLdY->LZZTimes)                                     );
             if (nivb!=niva)
              dTIdY->SetEntry( IndexA+1, IndexB,  MSIGN*(-dLdY->LZZTimes)                                     );
             dTIdY->SetEntry( IndexA+1, IndexB+1,       (KoZ*(L->LZZBullet + L->LZZNabla/Kappa2))       );

           };

       }; //  if ( Formulation == FORMULATION_EFIE ) ... else ... 

    }; // for (niva = ... ) for (nivb = ... )

  return 0;

}
 
/******************************************************************/
/* AssembleTI: Assemble the contribution to the T block of the    */
/* BEM matrix coming from the image object alone. This is         */
/* really doing the same thing as AssembleUab but with the added  */
/* bonus that the matrix is symmetric.                            */
/******************************************************************/
void AssembleTI(TDRTObject *O, TDRTObject *OImage, double Xi, double q,
                double EpsOut, double MuOut, int NumThreads,
                StaticSSIDataTable *SSSIDT, HMatrix *TI, HMatrix *dTIdY)
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double EpsIn, MuIn;
  if ( O->MP->IsPEC() )
   EpsIn=-1.0;
  else
   { cdouble zEps, zMu;
     O->MP->GetEpsMu( cdouble(0,Xi), &zEps, &zMu );
     EpsIn=real(zEps);
     MuIn=real(zMu);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[NumThreads], *TD;
  pthread_t *Threads = new pthread_t[NumThreads];
  for(int nt=0; nt<NumThreads; nt++)
   { 
     TD=&(TDS[nt]);
     TD->nt=nt;
     TD->NumTasks=NumThreads;

     TD->SSSIDT=SSSIDT;
     TD->O=O;
     TD->OImage=OImage;
     TD->Xi=Xi;
     TD->q=q;
     TD->EpsOut=EpsOut;
     TD->MuOut=MuOut;
     TD->EpsIn=EpsIn;
     TD->MuIn=MuIn;
     TD->TI=TI;
     TD->dTIdY=dTIdY;

     if ( (nt+1) == NumThreads )
      AssembleTI_Thread( (void *)&TD );
     else
      pthread_create( &(Threads[nt]), 0, AssembleTI_Thread, (void *)TD);
   };

  for(int nt=0; nt<NumThreads; nt++)
   pthread_join(Threads[nt],0);

  delete[] Threads;
  delete[] TDs;
#else
#ifndef USE_OPENMP
  int NumTasks=NumThreads=1;
#else
  int NumTasks=NumThreads*100;
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nt=0; nt<NumTasks; nt++)
   { 
     ThreadData TD1;
     TD1.nt=nt;
     TD1.NumTasks=NumTasks;

     TD1.SSSIDT=SSSIDT;
     TD1.O=O;
     TD1.OImage=OImage;
     TD1.Xi=Xi;
     TD1.q=q;
     TD1.EpsOut=EpsOut;
     TD1.MuOut=MuOut;
     TD1.EpsIn=EpsIn;
     TD1.MuIn=MuIn;
     TD1.TI=TI;
     TD1.dTIdY=dTIdY;

     AssembleTI_Thread((void *)&TD1);
   };
#endif

}


/***************************************************************/
/***************************************************************/
/***************************************************************/
TDRTObject *CreateImageObject(TDRTObject *O)
{ 
  int nv;
  char buffer[1000];

  TDRTObject *IO=(TDRTObject *)malloc(sizeof(TDRTObject));

  /*--------------------------------------------------------------*/
  /*- copy vertices and invert y components ----------------------*/
  /*--------------------------------------------------------------*/
  IO->NumVertices=O->NumVertices;
  IO->Vertices=(double *)memdup(O->Vertices, 2*O->NumVertices*sizeof(double));  
  for (nv=0; nv<IO->NumVertices; nv++)
   IO->Vertices[2*nv + 1] *= -1.0;

  /*--------------------------------------------------------------*/
  /*- copy Segments, IVs, Neighbors arrays. these are just copied */
  /*- verbatim since they only store indices into the Vertices    */
  /*- table.                                                      */
  /*--------------------------------------------------------------*/
  IO->NumSegments=O->NumSegments;
  IO->NumIVs=O->NumIVs;
  IO->NumBFs=O->NumBFs;

  IO->Segments=(int *)memdup(O->Segments,2*O->NumSegments * sizeof(int)); 
  IO->IVs=(int *)memdup(O->IVs,O->NumIVs * sizeof(int)); 
  IO->Neighbors=(int *)memdup(O->Neighbors,2*O->NumIVs * sizeof(int)); 
 
  memset(IO->Displacement,0,2*sizeof(double));
  IO->Rotation=0;

  snprintf(buffer,1000,"%s.image",O->MeshFileName);
  IO->MeshFileName=strdup(buffer);
  snprintf(buffer,1000,"%s_image",O->Label);
  IO->Label=strdup(buffer);
  IO->UnTransLine[0]=0;

  IO->MP = O->MP;

  return IO;

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteImageObjectPPMeshes(C2DWorkspace *W, const char *FileName, const char *Tag)
{ 
  FILE *f;
  TDRTObject *O;
  double *VA, *VB, Val;
  int i, no, ns;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  f=fopen(FileName,"a");
  if (!f) return;
  fprintf(f,"View.ShowElement=1;\n");
  fprintf(f,"View.ShowScale=0;\n");
  fprintf(f,"View \"%s_Images\" {\n",Tag);
  
  /***************************************************************/
  /* plot all line segments on all image objects. ****************/
  /***************************************************************/
  for(no=0, O=W->ImageObjects[no]; no<W->G->NumObjects; O=W->ImageObjects[++no])
   for(ns=0; ns<O->NumSegments; ns++)
     { 
       VA=O->Vertices + 2*O->Segments[2*ns];
       VB=O->Vertices + 2*O->Segments[2*ns+1];

       Val=(double)(no+1);

       fprintf(f,"SL(%e,%e,%e,%e,%e,%e) {%e,%e};\n",
                  VA[0], VA[1], 0.0, VB[0], VB[1], 0.0, Val, Val);
     };

   
  fprintf(f,"};\n");
  fclose(f);
}
