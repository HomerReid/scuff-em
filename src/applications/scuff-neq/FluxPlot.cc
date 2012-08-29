/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/***************************************************************/
/***************************************************************/
/***************************************************************/
void CreateFluxPlot(SNEQData *SNEQD, cdouble Omega, char *Tag)
{ 
  RWGGeometry *G = SNEQD->G;

  /***************************************************************/
  /* sanity check ************************************************/
  /***************************************************************/
  int ns;
  for(ns=0; ns<G->NumSurfaces; ns++)
   if ( G->Surfaces[ns]->MP->IsPEC() )
    ErrExit("flux plot not available for geometries containing PEC surfaces");

  /***************************************************************/
  /* allocate a vector with enough slots to store one double     */
  /* value per panel in the mesh                                 */
  /***************************************************************/
  double *PFV=(double *)mallocEC( (G->TotalPanels)*sizeof(double) ); // panel flux vector
  memset(PFV, 0, (G->TotalPanels)*sizeof(double));

  /***************************************************************/
  /* fill in the panel flux vector with the flux on each panel   */
  /***************************************************************/
  RWGSurface *S;
  RWGEdge *E;
  int ne, BFIndex, PanelIndex;
  double Value;
  int Offset = (G->NumSurfaces == 1) ? 0 : G->Surfaces[0]->NumBFs;
  for(ns=(G->NumSurfaces==1) ? 0 : 1; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], ne=0; ne<S->NumEdges; ne++)
    { 
      E=S->Edges[ne];
      
      BFIndex = G->BFIndexOffset[ns] - Offset + 2*ne;
      Value = 0.5 * ( SNEQD->DV->GetEntryD(BFIndex + 0) + SNEQD->DV->GetEntryD(BFIndex + 1) ); 

      PanelIndex = G->PanelIndexOffset[ns] + E->iPPanel;
      PFV[ PanelIndex ] += 0.5*Value / (S->Panels[E->iPPanel]->Area);

      PanelIndex = G->PanelIndexOffset[ns] + E->iMPanel;
      PFV[ PanelIndex ] += 0.5*Value / (S->Panels[E->iMPanel]->Area);

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
  for(PanelIndex=ns=0; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], np=0, P=S->Panels[0]; np<S->NumPanels; np++, PanelIndex++)
   {
      P=S->Panels[np];
      PV[0]=S->Vertices + 3*P->VI[0];
      PV[1]=S->Vertices + 3*P->VI[1];
      PV[2]=S->Vertices + 3*P->VI[2];

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
