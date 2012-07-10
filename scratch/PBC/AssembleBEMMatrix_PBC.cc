/*
 *
 */

#include <libscuff.h>
#include <libscuffInternals.h>
#include <PBCGeometry.h>

namespace scuff{

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddMEs(HMatrix *M, 
            RWGObject *Oa, int nea, RWGObject *Ob, int neb, 
            int RowOffset, int ColOffset, 
            cdouble PreFac[3], cdouble GC[2])
{
  int OaIsPEC = Oa->MP->IsPEC();
  int ObIsPEC = Ob->MP->IsPEC();
  int X, Y;
  int Symmetric=0;

  if ( OaIsPEC && ObIsPEC )
   { 
     X=RowOffset + nea;
     Y=ColOffset + neb;  

     M->AddEntry( X, Y, PreFac[0]*GC[0] );
   }
  else if ( OaIsPEC && !ObIsPEC )
   { 
     X=RowOffset + nea;
     Y=ColOffset + 2*neb;  

     B->AddEntry( X, Y,   PreFac[0]*GC[0] );
     B->AddEntry( X, Y+1, PreFac[1]*GC[1] );

   }
  else if ( !OaIsPEC && ObIsPEC )
   {
     X=RowOffset + 2*nea;
     Y=ColOffset + neb;  

     B->AddEntry( X,   Y, PreFac[0]*GC[0] );
     B->AddEntry( X+1, Y, PreFac[1]*GC[1] );
   }
  else // ( !OaIsPEC && !ObIsPEC )
   { 
     X=RowOffset + 2*nea;
     Y=ColOffset + 2*neb;  

     B->AddEntry( X, Y,   PreFac[0]*GC[0]);
     B->AddEntry( X, Y+1, PreFac[1]*GC[1]);
     if ( !Symmetric || (nea!=neb) )
      B->AddEntry( X+1, Y, PreFac[1]*GC[1]);
     B->AddEntry( X+1, Y+1, PreFac[2]*GC[0]);

   };


}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   PBCGeometry *PG;
   int nt, nTask;
 };

/***************************************************************/
/* thread routine for AddOuterContributions ********************/
/***************************************************************/
void AOC_Thread()
{
  /***************************************************************/
  /* extract local copies of fields in argument structure */
  /***************************************************************/
  ThreadData *TD=(ThreadData *)data;
  PBCGeometry *PG      = Args->PG;
  RWGGeometry *G       = PG->G;
  cdouble Omega        = PG->CurrenOmega;

  /***************************************************************/
  /* other local variables ***************************************/
  /***************************************************************/
  int nt=0;
  RWGObject *Oa, *Ob;
  int RowOffset, ColOffset;

  cdouble EpsExt=PG->EpsTF[0];
  cdouble MuExt=PG->MuTF[0];
  cdouble kExt=csqrt2(EpsExt*MuExt)*Omega;
  cdouble ExteriorPreFac[3];
  ExteriorPreFac[0] = II*MuExt*Omega;
  ExteriorPreFac[1] = -II*kExt;
  ExteriorPreFac[2] = -1.0*II*EpsExt*Omega;

  cdouble EpsInt, MuInt, kInt, InteriorPreFac[3];

  Interp3D *ExteriorInterpolator = PG->GBarAB9_Exterior;
  Interp3D *InteriorInterpolator;

  for(int noa=0; noa<G->NumObjects; noa++)
   for(int nob=0; nob<G->NumObjects; nob++)
    { 
      Oa=G->Objects[noa];
      RowOffset = G->BFIndexOffset[noa];
      Ob=G->Objects[nob];
      ColOffset = G->BFIndexOffset[nob];

      if ( noa==nob && PG->GBarAB9_Interior[noa] )
       { InteriorInterpolator = GBarGB9_Interior[noa];
         EpsInt=PG->EpsTF[noa+1];
         MuInt=PG->MuTF[noa+1];
         kInt=csqrt2(EpsInt*MuInt)*Omega;
         InteriorPreFac[0] = II*MuInt*Omega;
         InteriorPreFac[1] = -II*kInt;
         InteriorPreFac[2] = -1.0*II*EpsInt*Omega;
       }
      else
       InteriorInterpolator=0;

      for(int nea=0; nea<Oa->NumEdges; nea++)  
       for(int neb=0; neb<Ob->NumEdges; neb++)
        { 
          nt++;
          if (nt==TD->nTask) nt=0;
          if (nt!=TD->nt) continue;

          // contribution from exterior medium 
          GetAB9MatrixElements(Oa, nea, Ob, neb, k, ExteriorInterpolator, GC);
          AddMEs(M, Oa, nea, Ob, neb, RowOffset, ColOffset, ExteriorPreFac, GC);
          
          // contribution from interior medium if present 
          if (InteriorInterpolator)
           { 
             GetAB9MatrixElements(Oa, nea, Ob, neb, k, InteriorInterpolator, GC);
             AddMEs(M, Oa, nea, Ob, neb, RowOffset, ColOffset, InteriorPreFac, GC);
           };

        }; // for(nea=0 ... for(neb=0 ... 

    }; // for (noa=0 ... for(nob=0 ... 

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PBCGeometry::AddOuterCellContributions(double *P)
{ 
  /*--------------------------------------------------------------*/
  /*- initialize the interpolation tables at the present frequency*/
  /*- and bloch wavevector                                        */
  /*--------------------------------------------------------------*/
  GBarData MyGBarData, *GBD=&MyGBarData;

  cdouble Eps, Mu;
  GBarAB9Exterior->ReInitialize(nThread, GBarVDPhi3D,

  /*--------------------------------------------------------------*/
  /*- fire off threads -------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nThread=GetNumThreads();
#ifdef USE_PTHREAD
  ThreadData *TDs = new ThreadData[nThread], *TD;
  pthread_t *Threads = new pthread_t[nThread];
  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDs[nt]);
     TD->PG = this;
     TD->nt=nt;
     TD->nTask=nThread;
     if (nt+1 == nThread)
       AOC_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, AOC_Thread, (void *)TD);
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
     TD1.PG=this;
     TD1.nt=nt;
     TD1.nTask=nTask;
     AOC_Thread((void *)&TD1);
   };
#endif
  
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PBCGeometry::AssembleInnerCellBlocks()
{
  Log("Assembling inner matrix blocks at Omega=%s\n",z2s(CurrentOmega));

  Log(" MZZ block...");
  G->AssembleBEMMatrix(Omega, MZZ);

  int no, nop;

  ABMBArgStruct MyABMBArgStruct, *Args=&MyABMBArgStruct;
  Args->G = G; 
  Args->Omega = CurrentOmega;
  Args->nThread = GetNumThreads();
  Args->Symmetric=0;

  double Displacement[3];
  Displacement[2]=0.0;
  Args->Displacement=Displacement;

  Log("Assembling MPP block ...");
  Displacement[0]=LBV[0][0] + LBV[1][0];
  Displacement[1]=LBV[0][1] + LBV[1][1];
  Args->B=MPP;
  for(no=0; no<G->NumObjects; no++)
   for(nop=0; nop<G->NumObjects; nop++)
    { 
      Args->Oa=G->Objects[no];
      Args->RowOffset=G->BFIndexOffset[no];
      Args->Ob=G->Objects[nop];
      Args->ColOffset=G->BFIndexOffset[nop];
      Args->Symmetric=0;
      Args->B=MPP;

      // explain me
      if ( no==nop && !(Oa->MP->IsPEC()) && NumStraddlers[2*no+0]==0 || NumStraddlers[2*no+1]==0 )
       Oa->MP->Zero()

      AssembleBEMMatrixBlock(Args);

      if ( no==nop && !(Oa->MP->IsPEC()) && NumStraddlers[2*no+0]==0 || NumStraddlers[2*no+1]==0 )
       Oa->MP->UnZero()

    };

  Log("Assembling MPM block ...");
  Displacement[0]=LBV[0][0] - LBV[1][0];
  Displacement[1]=LBV[0][1] - LBV[1][1];
  Args->B=MPM;
  for(no=0; no<G->NumObjects; no++)
   for(nop=0; nop<G->NumObjects; nop++)
    { Args->Oa=G->Objects[no];
      Args->RowOffset=G->BFIndexOffset[no];
      Args->Ob=G->Objects[nop];
      Args->ColOffset=G->BFIndexOffset[nop];

      // explain me
      if ( no==nop && !(Oa->MP->IsPEC()) && NumStraddlers[2*no+0]==0 || NumStraddlers[2*no+1]==0 )
       Oa->MP->Zero()

      AssembleBEMMatrixBlock(Args);

      if ( no==nop && !(Oa->MP->IsPEC()) && NumStraddlers[2*no+0]==0 || NumStraddlers[2*no+1]==0 )
       Oa->MP->UnZero()

    };

  Log("Assembling MPZ block ...");
  Displacement[0]=LBV[0][0]; 
  Displacement[1]=LBV[0][1];
  Args->B=MPZ;
  for(no=0; no<G->NumObjects; no++)
   for(nop=0; nop<G->NumObjects; nop++)
    { Args->Oa=G->Objects[no];
      Args->RowOffset=G->BFIndexOffset[no];
      Args->Ob=G->Objects[nop];
      Args->ColOffset=G->BFIndexOffset[nop];

      // explain me
      if ( no==nop && !(Oa->MP->IsPEC()) && NumStraddlers[2*no+0]==0 )
       Oa->MP->Zero()

      AssembleBEMMatrixBlock(Args);

      if ( no==nop && !(Oa->MP->IsPEC()) && NumStraddlers[2*no+0]==0 )
       Oa->MP->UnZero()

    };

  Log("Assembling MZP block ...");
  Displacement[0]=LBV[1][0]; 
  Displacement[1]=LBV[1][1];
  Args->B=MZP;
  for(no=0; no<G->NumObjects; no++)
   for(nop=0; nop<G->NumObjects; nop++)
    { Args->Oa=G->Objects[no];
      Args->RowOffset=G->BFIndexOffset[no];
      Args->Ob=G->Objects[nop];
      Args->ColOffset=G->BFIndexOffset[nop];

      // explain me
      if ( no==nop && !(Oa->MP->IsPEC()) && NumStraddlers[2*no+1]==0 )
       Oa->MP->Zero()

      AssembleBEMMatrixBlock(Args);

      if ( no==nop && !(Oa->MP->IsPEC()) && NumStraddlers[2*no+1]==0 )
       Oa->MP->UnZero()

    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *PBCGeometry::AssembleBEMMatrix(cdouble Omega, double P[2], HMatrix *M)
{
  /*--------------------------------------------------------------*/
  /*- (re)allocate the matrix as necessary -----------------------*/
  /*--------------------------------------------------------------*/
  if ( M!=0 && M->NR!=TotalBFs || M->NC!=TotalBFs )
   { M=0;
     Warn("wrong-sized HMatrix passed to AssembleBEMMatrix() (reallocating...)");
   };

  if(M==0)

   M=new HMatrix(TotalBFs, TotalBFs, LHM_COMPLEX);

  /*--------------------------------------------------------------*/
  /*- assemble contributions of innermost lattice cells if the   -*/
  /*- frequency has changed since the last call                  -*/
  /*--------------------------------------------------------------*/
  if( CurrentOmega != Omega )
   { CurrentOmega=Omega;
     G->ExteriorMP->GetEpsMu(Omega, PG->EpsTF+0, PG->MuTF+0);
     for(no=0; no<G->NumObjects; no++)
      G->Objects[no]->MP->GetEpsMu(Omega, PG->EpsTF+no+1, PG->MuTF+no+1);
     AssembleInnerCellBlocks();
   };

  /*--------------------------------------------------------------*/
  /*- stamp in the contributions of the innermost lattice cells   */
  /*- note: PFPP = 'phase factor, plus-plus'                      */
  /*-       PFPZ = 'phase factor, plus-zero'                      */
  /*- etc.                                                        */
  /*--------------------------------------------------------------*/
  cdouble PFPP = exp( II* ( P[0]*(LBV[0][0] + LBV[1][0]) + P[1]*(LBV[0][1] + LBV[1][1]) ) );
  cdouble PFPM = exp( II* ( P[0]*(LBV[0][0] - LBV[1][0]) + P[1]*(LBV[0][1] - LBV[1][1]) ) );
  cdouble PFPZ = exp( II* ( P[0]*(LBV[0][0]            ) + P[1]*(LBV[0][1]            ) ) );
  cdouble PFZP = exp( II* ( P[0]*(            LBV[1][0]) + P[1]*(            LBV[1][1]) ) );
  for(int nr=0; nr<M->NR; nr++)
   for(int nc=0; nc<M->NR; nc++)
    M->SetEntry(nr, nc,  PFPP*MPP->GetEntry(nr, nc) + conj(PFPP)*MPP->GetEntry(nc,nr)
                       + PFPM*MPM->GetEntry(nr, nc) + conj(PFPM)*MPM->GetEntry(nc,nr)
                       + PFPZ*MPZ->GetEntry(nr, nc) + conj(PFPZ)*MPZ->GetEntry(nc,nr)
                       + PFZP*MZP->GetEntry(nr, nc) + conj(PFZP)*MZP->GetEntry(nc,nr)
                       + MZZ->GetEntry(nr,nc) 
               );
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  AddOuterCellContributions(P);

}

} //namespace scuff{
