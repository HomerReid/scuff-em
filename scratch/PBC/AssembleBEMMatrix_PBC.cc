/*
 *
 */

#include <libscuff.h>
#include <libscuffInternals.h>
#include <PBCGeometry.h>

namespace scuff{

#define ABSURD 1.234567e89

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddOuterContributions_Thread()
{
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PBCGeometry::GetOuterContributions(double *P)
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nThread=GetNumThreads();
  GBarData MyGBarData, *GBD=&MyGBarData;

  cdouble Eps, Mu;
  GBarAB9Exterior->ReInitialize(nThread, GBarVDPhi3D,

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(no=0; no<
  
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PBCGeometry::GetInnerContributions()
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
      AssembleBEMMatrixBlock(Args);
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
      AssembleBEMMatrixBlock(Args);
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
      AssembleBEMMatrixBlock(Args);
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
      AssembleBEMMatrixBlock(Args);
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
     GetInnerContributions();
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
  AddOuterContributions(P);

}

} // namespace scuff
