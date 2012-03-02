/*
 * FrequencyIntegrand.cc -- evaluate the spectral density of 
 *                       -- power transfer at a single frequency
 *
 * homer reid            -- 2/2012
 *
 */

#include "scuff-heat.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PlotFlux(ScuffHeatData *SHD, cdouble Omega)
{ 
  RWGGeometry *G = SHD->G;

  /***************************************************************/
  /* sanity check ************************************************/
  /***************************************************************/
  int no;
  for(no=0; no<G->NumObjects; no++)
   if ( G->Objects[no]->MP->IsPEC() )
    ErrExit("flux plot not available for geometries containing PEC objects");

  /***************************************************************/
  /* allocate a vector with enough slots to store one double     */
  /* value per panel in the mesh                                 */
  /***************************************************************/
  double *PFV=(double *)malloc( (G->TotalPanels)*sizeof(double) ); // panel flux vector
  memset(PFV, 0, (G->TotalPanels)*sizeof(double));

  /***************************************************************/
  /* fill in the panel flux vector with the flux on each panel   */
  /***************************************************************/
  RWGObject *O;
  RWGEdge *E;
  int ne, BFIndex, PanelIndex;
  double Value;
  for(no=0; no<G->NumObjects; no++)
   for(O=G->Objects[no], ne=0; ne<O->NumEdges; ne++)
    { 
      E=O->Edges[ne];
      
      BFIndex = G->BFIndexOffset[no] + 2*ne;
      Value = 0.5 * ( SHD->DV->GetEntryD(BFIndex + 0) + SHD->DV->GetEntryD(BFIndex + 1) ); 

      PanelIndex = G->PanelIndexOffset[no] + E->iPPanel;
      PFV[ PanelIndex ] += 0.5*Value / (O->Panels[E->iPPanel]->Area);

      PanelIndex = G->PanelIndexOffset[no] + E->iMPanel;
      PFV[ PanelIndex ] += 0.5*Value / (O->Panels[E->iMPanel]->Area);

    };
  
  /***************************************************************/
  /* create a GMSH postprocessing file containing a single 'view'*/
  /* that plots the flux on each panel                           */
  /***************************************************************/
  FILE *f=vfopen("%s.%g.flux","w",GetFileBase(G->GeoFileName),real(Omega));
  fprintf(f,"View \"Flux{\n");
  int np;
  RWGPanel *P;
  double *PV[3];
  for(PanelIndex=no=0; no<G->NumObjects; no++)
   for(O=G->Objects[no], np=0, P=O->Panels[0]; np<O->NumPanels; np++, PanelIndex++)
   {
      P=O->Panels[np];
      PV[0]=O->Vertices + 3*P->VI[0];
      PV[1]=O->Vertices + 3*P->VI[1];
      PV[2]=O->Vertices + 3*P->VI[2];

      fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                 PV[0][0], PV[0][1], PV[0][2],
                 PV[1][0], PV[1][1], PV[1][2],
                 PV[2][0], PV[2][1], PV[2][2],
                 PFV[PanelIndex],PFV[PanelIndex],PFV[PanelIndex]);
    };
  fprintf(f,"};\n\n");
  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFrequencyIntegrand(ScuffHeatData *SHD, cdouble Omega, double *FI)
{
  RWGGeometry *G = SHD->G;
  int nThread    = SHD->nThread;
  HMatrix *M0    = SHD->M0;
  HMatrix *M1    = SHD->M1;
  HMatrix *M2    = SHD->M2;

  /***************************************************************/
  /* assemble the three separate blocks that contribute to the   */
  /* BEM matrix:                                                 */
  /*  0: contributions of environment only                       */
  /*  1: contributions of object 1 only                          */
  /*  2: contributions of objects 2...N only                     */
  /***************************************************************/

  /*--------------------------------------------------------------*/
  /*- 0: contributions of environment only -----------------------*/
  /*--------------------------------------------------------------*/
  int no;
  for(no=0; no<G->NumObjects; no++)
   G->Objects[no]->MP->Zero();

  Log("Assembling M0 matrix...");
  G->AssembleBEMMatrix(Omega, nThread, M0);

  // multiply by the 'S' matrix 
  int nr, nc;
  for(nr=1; nr<M0->NR; nr+=2)
   for(nc=0; nc<M0->NC; nc++)
    M0->SetEntry(nr, nc, -1.0*M0->GetEntry(nr,nc));

  /*--------------------------------------------------------------*/
  /*- 1: contributions of object 1 only    -----------------------*/
  /*--------------------------------------------------------------*/
  G->ExteriorMP->Zero();
  G->Objects[0]->MP->UnZero();

  Log("Assembling M1 matrix...");
  G->AssembleBEMMatrix(Omega, nThread, M1);

  for(nr=1; nr<M1->NR; nr+=2)
   for(nc=0; nc<M1->NC; nc++)
    M1->SetEntry(nr, nc, -1.0*M1->GetEntry(nr,nc));

  /*--------------------------------------------------------------*/
  /*- 2: contributions of objects 2--N only.                      */
  /*-                                                             */
  /*-  NOTE: if we only have a single object then we set M2 = M0; */
  /*-        the calculation them amounts to computing the heat   */
  /*-        transfer the single body to the environment.         */
  /*--------------------------------------------------------------*/
  if (G->NumObjects==1)
   M2->Copy(M0);
  else // (G->NumObjects>1)
   { 
     G->Objects[0]->MP->Zero();
     for(no=1; no<G->NumObjects; no++)
      G->Objects[no]->MP->UnZero();

     Log("Assembling M2 matrix...");
     G->AssembleBEMMatrix(Omega, nThread, M2);

     for(nr=1; nr<M2->NR; nr+=2)
      for(nc=0; nc<M2->NC; nc++)
       M2->SetEntry(nr, nc, -1.0*M2->GetEntry(nr,nc));
   };
  
  // undo the zeroing out of the environment and object 1
  G->ExteriorMP->UnZero();
  G->Objects[0]->MP->UnZero();

  /***************************************************************/
  /* assemble and LU-factorize the full BEM matrix (for now we   */
  /* do this in-place using the M0 matrix)                       */
  /***************************************************************/
  if (G->NumObjects==1)
   { for(nr=0; nr<M0->NR; nr++)
      for(nc=0; nc<M0->NC; nc++)
       M0->AddEntry(nr, nc, M1->GetEntry(nr,nc) );
   }
  else
   { for(nr=0; nr<M0->NR; nr++)
      for(nc=0; nc<M0->NC; nc++)
       M0->AddEntry(nr, nc, M1->GetEntry(nr,nc) + M2->GetEntry(nr,nc) );
   };

  Log("LU-factorizing MTotal...");
  M0->LUFactorize();

  /***************************************************************/
  /* set M1 = sym(M1) and M2 = sym(M2)                           */
  /*  (note that the loop runs over only the upper triangle of   */
  /*  the matrices)                                              */
  /***************************************************************/
  cdouble Sym, SymT;
  for(nr=0; nr<M1->NR; nr++)
   for(nc=nr; nc<M1->NC; nc++)
    { 
      Sym = 0.5*(M1->GetEntry(nr, nc) + conj(M1->GetEntry(nc, nr)));
      M1->SetEntry(nr, nc, Sym );
      if(nc>nr) M1->SetEntry(nc, nr, conj(Sym) );
    }; 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (SHD->PlotFlux)
   {
     for(nr=0; nr<M2->NR; nr++)
      for(nc=nr; nc<M2->NC; nc++)
       { Sym  = conj(M2->GetEntry(nr, nc));
         SymT = conj(M2->GetEntry(nc, nr));
         M2->SetEntry(nc, nr, Sym );
         if (nc>nr) M2->SetEntry(nr, nc, Sym );
       };

     Log("LU-solving M1...");
     M0->LUSolve(M1,'N');
     Log("LU-solving M2...");
     M0->LUSolve(M2,'C');
     Log("Multipliying...");
     M1->Multiply(M2, M0);
     for(nr=0; nr<M1->NR; nr++)
      SHD->DV->SetEntry(nr, M1->GetEntry(nr,nr));
     Log("Plotting flux vector...");
     PlotFlux(SHD, Omega);

     *FI=0.0;
   }
  else
   {
     for(nr=0; nr<M2->NR; nr++)
      for(nc=nr; nc<M2->NC; nc++)
       { Sym = 0.5*(M2->GetEntry(nr, nc) + conj(M2->GetEntry(nc, nr)));
         M2->SetEntry(nr, nc, Sym );
         if(nc>nr) M2->SetEntry(nc, nr, conj(Sym) );
       };

     /***************************************************************/
     /* set M1 <= M^{-1'} * M1 **************************************/
     /* set M2 <= M^{-1}  * M2 **************************************/
     /***************************************************************/
     Log("LU-solving M1...");
     M0->LUSolve(M1,'C');
     Log("LU-solving M2...");
     M0->LUSolve(M2,'N');
   
     /***************************************************************/
     /* set M0 = M1*M2                                              */
     /***************************************************************/
     Log("Multiplying M1*M2...");
     M1->Multiply(M2, M0);
   
     /***************************************************************/
     /* the value of the frequency integrand is now the trace of M0 */
     /***************************************************************/
     *FI = real( M0->GetTrace() ) / 8.0 ;
     Log("...done!");

     /***************************************************************/
     /* write the result to the frequency-resolved output file ******/
     /***************************************************************/
     FILE *f=fopen(SHD->ByOmegaFile, "a");
     fprintf(f,"%s %e\n",z2s(Omega),*FI);
     fclose(f);
   };
   
}
