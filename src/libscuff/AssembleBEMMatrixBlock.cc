/*
 * AssembleBEMMatrixBlock.cc -- libscuff routine for assembling a single 
 *                              block of the BEM matrix (i.e. the       
 *                              interactions of two objects in the geometry)
 * 
 * homer reid                -- 10/2006 -- 10/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <libhmat.h>
#include <libhrutil.h>

#include "libscuff.h"

#define II cdouble(0,1)

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   ABMBArgStruct *Args;
   int nt;

 } ThreadData;

/***************************************************************/
/* 'AssembleBMatrixBlockThread'                                */
/***************************************************************/
static void *ABMBThread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  /* extract local copies of fields in argument structure */
  ABMBArgStruct *Args  = TD->Args;
  RWGGeometry *G       = Args->G;
  RWGObject *Oa        = Args->Oa;
  RWGObject *Ob        = Args->Ob;
  cdouble Frequency    = Args->Frequency;
  int nThread          = Args->nThread;
  int NumTorqueAxes    = Args->NumTorqueAxes;
  double *GammaMatrix  = Args->GammaMatrix;
  int RowOffset        = Args->RowOffset;
  int ColOffset        = Args->ColOffset;
  int ForceSymmetric   = Args->ForceSymmetric;
  HMatrix *B           = Args->B;
  HMatrix **GradB      = Args->GradB;
  HMatrix **dBdTheta   = Args->dBdTheta;
  double Sign          = Args->Sign;
  cdouble EpsA         = Args->EpsA;
  cdouble EpsB         = Args->EpsB;
  double MuA           = Args->MuA;
  double MuB           = Args->MuB;
  int OaIsPEC          = Args->OaIsPEC;
  int ObIsPEC          = Args->ObIsPEC;

  /* other local variables */
  int nt, nea, NEa, BFMultA, neb, NEb, BFMultB;

  NEa=Oa->NumEdges;
  BFMultA = Oa->MP->IsPEC() ? 1 : 2;

  NEb=Ob->NumEdges;
  BFMultB = Ob->MP->IsPEC() ? 1 : 2;


  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/
#ifdef _GNU_SOURCE
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(TD->nt,&cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  EEPreFac = Sign*II*

  /***************************************************************/
  /* loop over all internal edges on both objects.               */
  /*                                                             */
  /* note: if the block matrix B is symmetric, or if the caller  */
  /* specified ForceSymmetric, then the second loop starts at    */
  /* the current value of the first loop counter; otherwise the  */
  /* second loop starts at 0                                     */
  /***************************************************************/
  int nt=0;
  if ( ForceSymmetric || (B->StorageType!=LHM_NORMAL) )
   Symmetric=1;
  else
   Symmetric=0; 
  for(nea=0; nea<NEa; nea++)
   for(neb=Symmetric*nea; neb<NEb; neb++)
    { 
      nt++;
      if (nt==TD->nThread) nt=0;
      if (nt!=TD->nt) continue;

      /*--------------------------------------------------------------*/
      /*- contributions of first medium (EpsA, MuA)  -----------------*/
      /*--------------------------------------------------------------*/
      EEIArgs->nea  = nea;
      EEIArgs->neb  = neb;
      EEIArgs->Eps  = EpsA;
      EEIArgs->Mu   = MuA;
      GetEEIs(&EEIArgs);

      if ( OaIsPEC && ObIsPEC )
       { B->SetEntry( RowOffset + nea, ColOffset + neb, EEPreFacA*GInt )
         for(Mu=0; Mu<NumGradientComponents; Mu++)
          GradB[Mu]->SetEntry( RowOffset + nea, ColOffset + neb, EEPrefacA*GradGInt[Mu]);
         for(Mu=0; Mu<NumTorqueAxes; Mu++)
          dBdTheta[Mu]->SetEntry( RowOffset + nea, ColOffset + neb, EEPrefacA*dGIntdTheta[Mu]);
       }
      else if ( OaIsPEC && !ObIsPEC )
       { B->SetEntry( RowOffset + nea, ColOffset + 2*neb,   EEPreFacA*GInt )
         B->SetEntry( RowOffset + nea, ColOffset + 2*neb+1, EMPreFacA*CInt )
         for(Mu=0; Mu<NumGradientComponents; Mu++)
          { GradB[Mu]->SetEntry( RowOffset + nea, ColOffset + 2*neb,   EEPrefacA*GradGInt[Mu]);
            GradB[Mu]->SetEntry( RowOffset + nea, ColOffset + 2*neb+1, EMPrefacA*GradCInt[Mu]);
          };
         for(Mu=0; Mu<NumTorqueAxes; Mu++)
          { dBdTheta[Mu]->SetEntry( RowOffset + nea, ColOffset + 2*neb,   EEPrefacA*dGIntdTheta[Mu]);
            dBdTheta[Mu]->SetEntry( RowOffset + nea, ColOffset + 2*neb+1, EMPrefacA*dCIntdTheta[Mu]);
          { dBdTheta[Mu]->SetEntry( RowOffset + nea, ColOffset + 2*neb,   EEPrefacA*dGIntdTheta[Mu]);
       }

    };

  return 0;

}

/***************************************************************/  
/***************************************************************/  
/***************************************************************/
void RWGGeometry::AssembleU(int noa, int nob, double Frequency, int RealFreq,
                            int NumTorqueAxes, double *GammaMatrix, int nThread, 
                            HMatrix *U, 
                            HMatrix *dUdX, HMatrix *dUdY, HMatrix *dUdZ, 
                            HMatrix *dUdTheta1, HMatrix *dUdTheta2, HMatrix *dUdTheta3)
{ 
  int i, j, nt;

  ThreadData TDS[nThread], *TD;
  pthread_t Threads[nThread];
  static void **pLFWs=0;
  static int nThreadSave=0;

  /***************************************************************/
  /* update the cached values of Epsilon and Mu inside the two   */
  /* objects and inside (this)                                   */
  /***************************************************************/
  MP->GetEpsMu(Frequency, RealFreq, &EpsThisFreq, &MuThisFreq);
  Objects[noa]->MP->GetEpsMu(Frequency, RealFreq, 
                             &(Objects[noa]->EpsThisFreq), 
                             &(Objects[noa]->MuThisFreq));
  Objects[nob]->MP->GetEpsMu(Frequency, RealFreq, 
                             &(Objects[nob]->EpsThisFreq), 
                             &(Objects[nob]->MuThisFreq));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (nThreadSave!=nThread)
   { nThreadSave=nThread;
     pLFWs=(void **)RWGMalloc(nThread*sizeof(void *));
     for(nt=0; nt<nThread; nt++)
      pLFWs[nt]=CreateLFWorkspace();
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDS[nt]);
     TD->nt=nt;
     TD->nThread=nThread;

     TD->G=this;
     TD->Oa=Objects[noa];
     TD->Ob=Objects[nob];
     TD->Frequency=Frequency;
     TD->RealFreq=RealFreq;
     TD->NumTorqueAxes=NumTorqueAxes;
     TD->GammaMatrix=GammaMatrix;
     TD->U=U;
     TD->dUdX=dUdX;
     TD->dUdY=dUdY;
     TD->dUdZ=dUdZ;
     TD->dUdTheta1=dUdTheta1;
     TD->dUdTheta2=dUdTheta2;
     TD->dUdTheta3=dUdTheta3;
     TD->pLFW=pLFWs[nt];

     pthread_create( &(Threads[nt]), 0, AssembleU_Thread, (void *)TD);
   };

  for(nt=0; nt<nThread; nt++)
   pthread_join(Threads[nt],0);

}


/***************************************************************/
/* AssembleT_Thread                                            */
/***************************************************************/
static void *AssembleT_Thread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  /* local copies of fields in ThreadData structure */
  int nThread         = TD->nThread;
  RWGGeometry *G      = TD->G;
  RWGObject *O        = TD->O;
  void *pLFW          = TD->pLFW;
  double Frequency    = TD->Frequency;
  int RealFreq        = TD->RealFreq;
  HMatrix *T          = TD->T;

  /* other local variables */
  int nt, ne, nep, NE;
  int BFMult=O->MP->IsPEC() ? 1 : 2;
  int PerCent=0;

  NE=O->NumEdges;

  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/
#ifdef _GNU_SOURCE
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(TD->nt,&cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#endif
  /*--------------------------------------------------------------*/
  /*- EXPERIMENTAL -----------------------------------------------*/
  /*--------------------------------------------------------------*/

  nt=0;
  for(ne=0; ne<NE; ne++)
   for(nep=ne; nep<NE; nep++)
    { 
      nt++;
      if (nt==TD->nThread) nt=0;
      if (nt!=TD->nt) continue;

      for(PerCent=0; PerCent<9; PerCent++)
       if ( nep==ne &&  (ne == (PerCent*NE)/10) )
        Log("%i0 %% (%i/%i)...",PerCent,ne,NE);

      G->StampMatrixElements(pLFW, O, ne, BFMult*ne, O, nep, BFMult*nep,
                             Frequency, RealFreq, 0, 0,
                             T, 0, 0, 0, 0, 0, 0);

    };
  
  return 0;

}
 
/***************************************************************/  
/***************************************************************/  
/***************************************************************/
void AssembleBEMMatrixBlock(ABMBArgStruct *Args)
{ 
  /***************************************************************/
  /* figure out the container relationship between the objects   */
  /* to get values for Sign, EpsA, MuA, EpsB, MuB.               */
  /* note: EpsA, MuA are the material properties of the medium   */
  /*       through which the two objects interact                */
  /*       if the two objects are identical, then EpsB, MuB are  */
  /*       the material properties of the medium interior to the */
  /*       object; otherwise, EpsB is set to 0.0.                */
  /***************************************************************/
  RWGGeometry *G=Args->Ga;
  cdouble Frequency=Args->Frequency;
  RWGObject *Oa=Args->Oa;
  RWGObject *Ob=Args->Ob;
  if ( Ob->ContainingObject == Oa )
   { Args->Sign=-1.0;
     Oa->MP->GetEpsMu(Frequency, &(Args->EpsA), &(Args->MuA) );
   }
  else if ( Oa->ContainingObject == Ob )
   { Args->Sign=-1.0;
     Ob->MP->GetEpsMu(Frequency, &(Args->EpsA), &(Args->MuA) );
   } 
  else if ( Oa->ContainingObject == Ob->ContainingObject )
   { Args->Sign=1.0;
     if (Oa->ContainingObject==0)
      G->MP->GetEpsMu(Frequency, &(Args->EpsA), &(Args->MuA) );
     else 
      Oa->ContainingObject->MP->GetEpsMu(Frequency, &(Args->EpsA), &(Args->MuA) );
   };

  Args->OaIsPEC = Oa->MP->IsPEC();
  Args->ObIsPEC = Ob->MP->IsPEC();
  if ( !(Args->OaIsPEC) && !(Args->ObIsPEC) && Oa==Ob ) 
   Oa->MP->GetEpsMu(Frequency, &(Args->EpsB), &(Args->MuB) );
  else
   Args->EpsB=0.0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int nt, nThread=Args->nThread;
  ThreadData TDs[nThread];
  pthread_t Threads[nThread];

  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDS[nt]);
     TD->nt=nt;
     TD->nThread=nThread;

     TD->G=this;
     TD->O=Objects[no];
     TD->Frequency=Frequency;
     TD->RealFreq=RealFreq;
     TD->T=T;
     TD->pLFW=pLFWs[nt]; 

     pthread_create( &(Threads[nt]), 0, AssembleT_Thread, (void *)TD);
   };

  for(nt=0; nt<nThread; nt++)
   pthread_join(Threads[nt],0);

}
