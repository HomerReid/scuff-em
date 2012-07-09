/*
 *
 */

namespace scuff{

#define ABSURD 1.234567e89

bool EqualFloat(const double *a, const double *b) 
{
   return ( float(a) == float(b) );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
/***************************************************************/
void SetOmega(RWGGeometry *G, PBCAccelerator *PBCA, cdouble Omega)
{

  // the change in Omega value invalidates
  PBCA->CurrentP[0]=ABSURD;
  PBCA->CurrentP[1]=ABSURD;
  
}

void SetP(RWGGeometry *G, PBCAccelerator *PBCA, cdouble P)
{
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleBEMMatrix_PBC(RWGGeometry *G, cdouble Omega, 
                           double **LBV, double *P, PBCAccelerator *PBCA, 
                           HMatrix *M)
{
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if( PBCA->CurrentOmega != Omega )
   SetOmega(G, PBCA, Omega);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if( !EqualFloat(PBCA->CurrentP[0], P[0]) || !EqualFloat(PBCA->CurrentP[1], P[1]) )
   SetP(G, PBCA, P);

  /*--------------------------------------------------------------*/
  /*- stamp in the contributions of the innermost 9 lattice cells */
  /*- note: PFPP = 'phase factor, plus-plus'                      */
  /*-       PFPZ = 'phase factor, plus-zero'                      */
  /*- etc.                                                        */
  /*--------------------------------------------------------------*/
  cdouble PFPP = exp( II* ( P[0]*(LBV[0][0] + LBV[1][0]) + P[1]*(LBV[0][1] + LBV[1][1]) ) );
  cdouble PFPM = exp( II* ( P[0]*(LBV[0][0] - LBV[1][0]) + P[1]*(LBV[0][1] - LBV[1][1]) ) );
  cdouble PFPZ = exp( II* ( P[0]*(LBV[0][0]            ) + P[1]*(LBV[0][1]            ) ) );
  cdouble PFZP = exp( II* ( P[0]*(            LBV[1][0]) + P[1]*(            LBV[1][1]) ) );
  HMatrix *MPP=PBCA->MPP;
  HMatrix *MPM=PBCA->MPM;
  HMatrix *MPZ=PBCA->MPZ;
  HMatrix *MZP=PBCA->MZP;
  HMatrix *MZZ=PBCA->MZZ;

  for(nr=0; nr<M->NR; nr++)
   for(nc=0; nc<M->NR; nc++)
    M->SetEntry(nr, nc,  PFPP*MPP->GetEntry(nr, nc) + conj(PFPP)*MPP->GetEntry(nc,nr)
                       + PFPM*MPM->GetEntry(nr, nc) + conj(PFPM)*MPM->GetEntry(nc,nr)
                       + PFPZ*MPZ->GetEntry(nr, nc) + conj(PFPZ)*MPZ->GetEntry(nc,nr)
                       + PFZP*MZP->GetEntry(nr, nc) + conj(PFZP)*MZP->GetEntry(nc,nr)
                       + MZZ->GetEntry(nr,nc) 
               );
  
  /*--------------------------------------------------------------*/
  /*- compute the contributions of the outer 9 grid cells         */
  /*--------------------------------------------------------------*/
  
  
  

}

} // namespace scuff
