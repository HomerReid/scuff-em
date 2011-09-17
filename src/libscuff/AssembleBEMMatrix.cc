/*
 * AssembleM.cc -- libscuff routine for assembling BEM matrices
 * 
 * homer reid   -- 10/2006 -- 6/2011
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
   /* these fields only used by AssembleUb */
   RWGObject *Oa, *Ob;
   HMatrix *U, *dUdX, *dUdY, *dUdZ; 
   HMatrix *dUdTheta1, *dUdTheta2, *dUdTheta3;
   int NumTorqueAxes;
   double *GammaMatrix;

   /* these fields only used by AssembleT */
   RWGObject *O;
   HMatrix *T;

   /* these fields used by both */
   RWGGeometry *G;
   int nt, nThread;
   double Frequency;
   int RealFreq;    /* =1 for real frequency, 0 for imag */
   void *pLFW;

 } ThreadData;

/***************************************************************/
/* AssembleU_Thread          *********************************/
/***************************************************************/
static void *AssembleU_Thread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  /* local copies of fields in ThreadData structure */
  RWGGeometry *G       = TD->G;
  RWGObject *Oa        = TD->Oa;
  RWGObject *Ob        = TD->Ob;
  HMatrix *U           = TD->U;
  HMatrix *dUdX        = TD->dUdX;
  HMatrix *dUdY        = TD->dUdY;
  HMatrix *dUdZ        = TD->dUdZ;
  HMatrix *dUdTheta1   = TD->dUdTheta1;
  HMatrix *dUdTheta2   = TD->dUdTheta2;
  HMatrix *dUdTheta3   = TD->dUdTheta3;
  double NumTorqueAxes = TD->NumTorqueAxes;
  double *GammaMatrix  = TD->GammaMatrix;
  int nThread          = TD->nThread;
  void *pLFW           = TD->pLFW;
  double Frequency     = TD->Frequency;
  int RealFreq         = TD->RealFreq;

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

  nt=0;
  for(nea=0; nea<NEa; nea++)
   for(neb=0; neb<NEb; neb++)
    { 
      nt++;
      if (nt==TD->nThread) nt=0;
      if (nt!=TD->nt) continue;

      G->StampMatrixElements(pLFW, 
                             Oa, nea, BFMultA*nea, 
                             Ob, neb, BFMultB*neb, 
                             Frequency, RealFreq, NumTorqueAxes, GammaMatrix,
                             U, dUdX, dUdY, dUdZ, 
                             dUdTheta1, dUdTheta2, dUdTheta3);
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
void RWGGeometry::AssembleT(int no, double Frequency, int RealFreq,
                            int nThread, HMatrix *T)
{ 
  int nt;

  ThreadData TDS[nThread], *TD;
  pthread_t Threads[nThread];
  static void **pLFWs=0;
  static int nThreadSave=0;

  /***************************************************************/
  /* update the values of Epsilon and Mu cached within the object*/
  /* and within (this)                                           */
  /***************************************************************/
  MP->GetEpsMu(Frequency, RealFreq, &EpsThisFreq, &MuThisFreq);
  Objects[no]->MP->GetEpsMu(Frequency, RealFreq, 
                            &(Objects[no]->EpsThisFreq), 
                            &(Objects[no]->MuThisFreq));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (nThreadSave!=nThread)
   { nThreadSave=nThread;
     pLFWs=(void **)RWGMalloc(nThread*sizeof(void *));
     for(nt=0; nt<nThread; nt++)
      pLFWs[nt]=CreateLFWorkspace();
   };

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

/***************************************************************/
/* Assemble the full BEM matrix. *******************************/
/***************************************************************/
void RWGGeometry::AssembleBEMMatrix(double Frequency, int RealFreq,
                                    int nThread, HMatrix *M)
{ 
  int no, nop, N, NP, RowOffset, ColOffset;
  cdouble Eps;
  double Mu;
  HMatrix *T, *U;

  if (nThread<=0)
   ErrExit("AssembleBEMMatrix called with nThread=%i",nThread);

  MP->GetEpsMu(Frequency, RealFreq, &Eps, &Mu);
  Log("Assembling BEM matrix at frequency %g (medium eps,mu=%g+%gi,%g)",
       Frequency,real(Eps),imag(Eps),Mu);


  M->Zero();

  /***************************************************************/
  /* T blocks ****************************************************/
  /***************************************************************/
  for(no=0; no<NumObjects; no++)
   { 
     if ( Mate[no]!=-1 ) continue;

     Objects[no]->MP->GetEpsMu(Frequency, RealFreq, &Eps, &Mu);
     Log(" T[%i] (%s) (eps,mu=%g+%gi,%g)",no,Objects[no]->MeshFileName,
           real(Eps),imag(Eps),Mu);

     N=Objects[no]->NumBFs;
     if (RealFreq==LHM_REAL)
      T=new HMatrix(N,N,LHM_REAL,LHM_SYMMETRIC);
     else
      T=new HMatrix(N,N,LHM_COMPLEX,LHM_SYMMETRIC);

     AssembleT(no, Frequency, RealFreq, nThread, T);
  
     M->InsertBlock(T,BFIndexOffset[no], BFIndexOffset[no]);

     for(nop=no+1; nop<NumObjects; nop++)
      if ( Mate[nop]==no )
       M->InsertBlock(T,BFIndexOffset[nop], BFIndexOffset[nop]);

     delete T;

   };

  /***************************************************************/
  /* U blocks.                                                   */
  /***************************************************************/
  for(no=0; no<NumObjects; no++)
   for(nop=no+1; nop<NumObjects; nop++)
    { 
      N=Objects[no]->NumBFs;
      NP=Objects[nop]->NumBFs;
      U=new HMatrix(N,NP,(RealFreq ? LHM_COMPLEX : LHM_REAL),LHM_NORMAL);

      Log(" U[%i,%i]",no,nop);
      AssembleU(no, nop, Frequency, RealFreq,
                0, 0, nThread, U, 0, 0, 0, 0, 0, 0);

      RowOffset=BFIndexOffset[no];
      ColOffset=BFIndexOffset[nop];

      M->InsertBlock(U,RowOffset,ColOffset);

      delete U;
    };

} 

/***************************************************************/
/* Allocate a new HMatrix to store the BEM matrix.    **********/
/***************************************************************/
HMatrix *RWGGeometry::AllocateBEMMatrix(int RealFreq)
{ 
  HMatrix *M;

  if (RealFreq)
   M=new HMatrix(TotalBFs,TotalBFs,LHM_COMPLEX,LHM_SYMMETRIC);
  else
   M=new HMatrix(TotalBFs,TotalBFs,LHM_REAL,LHM_SYMMETRIC);

  return M;

} 
