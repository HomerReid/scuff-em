/***************************************************************/
/***************************************************************/
/***************************************************************/
void CreateFluxPlot(SNEQData *SNEQD, cdouble Omega, char *Tag)
{ 
  RWGGeometry *G = SNEQD->G;

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
  double *PFV=(double *)mallocEC( (G->TotalPanels)*sizeof(double) ); // panel flux vector
  memset(PFV, 0, (G->TotalPanels)*sizeof(double));

  /***************************************************************/
  /* fill in the panel flux vector with the flux on each panel   */
  /***************************************************************/
  RWGObject *O;
  RWGEdge *E;
  int ne, BFIndex, PanelIndex;
  double Value;
  int Offset = (G->NumObjects == 1) ? 0 : G->Objects[0]->NumBFs;
  for(no=(G->NumObjects==1) ? 0 : 1; no<G->NumObjects; no++)
   for(O=G->Objects[no], ne=0; ne<O->NumEdges; ne++)
    { 
      E=O->Edges[ne];
      
      BFIndex = G->BFIndexOffset[no] - Offset + 2*ne;
      Value = 0.5 * ( SNEQD->DV->GetEntryD(BFIndex + 0) + SNEQD->DV->GetEntryD(BFIndex + 1) ); 

      PanelIndex = G->PanelIndexOffset[no] + E->iPPanel;
      PFV[ PanelIndex ] += 0.5*Value / (O->Panels[E->iPPanel]->Area);

      PanelIndex = G->PanelIndexOffset[no] + E->iMPanel;
      PFV[ PanelIndex ] += 0.5*Value / (O->Panels[E->iMPanel]->Area);

    };
  
  /***************************************************************/
  /* create a GMSH postprocessing file containing a single 'view'*/
  /* that plots the flux on each panel                           */
  /***************************************************************/
  FILE *f=vfopen("%s.%g.%s.flux.pp","w",GetFileBase(G->GeoFileName),real(Omega),Tag);
  fprintf(f,"View \"Flux\"{\n");
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

      if ( PFV[PanelIndex]!=0.0 )
       fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                  PV[0][0], PV[0][1], PV[0][2],
                  PV[1][0], PV[1][1], PV[1][2],
                  PV[2][0], PV[2][1], PV[2][2],
                  PFV[PanelIndex],PFV[PanelIndex],PFV[PanelIndex]);
    };
  fprintf(f,"};\n\n");
  fclose(f);


}
