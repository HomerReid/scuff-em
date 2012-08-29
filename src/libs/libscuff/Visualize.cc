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

/*
 * Visualize.cc
 * 
 * homer reid -- 11/2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "libscuff.h"

namespace scuff {

#define II cdouble(0,1)

#define MAXSTR 1000

/************************************************************/
/* subroutines for emitting GMSH postprocessing code        */
/************************************************************/
/* vector point (otherwise known as 'arrow') */
void WriteVP(double *X, double *V, FILE *f)
{
  fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",X[0],X[1],X[2],V[0],V[1],V[2]);
}

/* scalar triangle */
void WriteST(double **VV, double Val, FILE *f)
{
  fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
             VV[0][0], VV[0][1], VV[0][2],
             VV[1][0], VV[1][1], VV[1][2],
             VV[2][0], VV[2][1], VV[2][2], 
             Val,Val,Val);
}

/***************************************************************/
/* WritePPMesh routine. Writes geometry to a .pp file suitable */
/* for opening as a GMSH postprocessing file.                  */
/* Note calling convention and file handling are different     */
/* from those of WriteGPMesh.                                  */
/* If PlotNormals is nonzero, we also plot the normal to each  */
/* panel.                                                      */
/***************************************************************/
void RWGGeometry::WritePPMesh(const char *FileName, const char *Tag, int PlotNormals)
{ 
  FILE *f;
  RWGSurface *S;
  RWGPanel *P;
  char buffer[MAXSTR], *p;
  double *PV[3], Val;
  int ns, np;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  strncpy(buffer,FileName,996);
  p=strrchr(buffer,'.');
  if ( !p || strcmp(p,".pp") )
   strcat(buffer,".pp");

  f=fopen(buffer,"a");
  if (!f) 
   { fprintf(stderr,"warning: could not open file %s \n",FileName);
     return;
   };

  /***************************************************************/
  /* plot all panels on all objects  *****************************/
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n",Tag);
  for(ns=0, S=Surfaces[0]; ns<NumSurfaces; S=Surfaces[++ns])
   for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
    { 
      PV[0]=S->Vertices + 3*P->VI[0];
      PV[1]=S->Vertices + 3*P->VI[1];
      PV[2]=S->Vertices + 3*P->VI[2];
      Val=(double)(ns+1);
      fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                 PV[0][0], PV[0][1], PV[0][2],
                 PV[1][0], PV[1][1], PV[1][2],
                 PV[2][0], PV[2][1], PV[2][2],
                 Val,Val,Val);
    };
  fprintf(f,"};\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=1;\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowScale=0;\n");

  /***************************************************************/
  /* additionally plot panel normals if that was also requested **/
  /***************************************************************/
  if (PlotNormals)
   { 
     fprintf(f,"View \"%s.Normals\" {\n",Tag);
     for(ns=0, S=Surfaces[0]; ns<NumSurfaces; S=Surfaces[++ns])
      for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
       { 
          Val=fmax(VecNorm(P->Centroid), 2.0*P->Area);
          fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                     P->Centroid[0], P->Centroid[1], P->Centroid[2], 
                     Val*P->ZHat[0], Val*P->ZHat[1], Val*P->ZHat[2]);
       };
     fprintf(f,"};\n");
   }; /* if (PlotNormals) */

  fclose(f);
}

/***************************************************************/
/* WriteGPMesh routine. This writes out the panel set as a     */
/* collection of lines. Useful for visualizing in gnuplot.     */
/*                                                             */
/* The mesh (with normals) may be plotted with the gnuplot cmd */
/*  splot 'file.gpmsh'  w l                                    */
/* while                                                       */
/*  splot 'file.gpmsh' i 0:10000:2 w l                         */
/* will show the mesh sans normals.                            */
/*                                                             */
/***************************************************************/
void RWGGeometry::WriteGPMesh(const char *format, ...)
{ 
  va_list ap;
  char FileName[MAXSTR], *p;
  FILE *f;
  RWGSurface *S;
  RWGPanel *P;
  double *PV[3];
  int i, ns, np;

  va_start(ap,format);
  vsnprintfEC(FileName,997,format,ap);
  va_end(ap);

  p=strrchr(FileName,'.');
  if ( !p || strcmp(p,".gp") )
   strcat(FileName,".gp");

  f=fopen(FileName,"w");
  if (!f) 
   { fprintf(stderr,"warning: could not open file %s \n",FileName);
     return;
   };
  
  /***************************************************************/
  /* plot mesh ***************************************************/
  /***************************************************************/
  for(ns=0, S=Surfaces[0]; ns<NumSurfaces; S=Surfaces[++ns])
   for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
     { 
       PV[0]=S->Vertices + 3*P->VI[0];
       PV[1]=S->Vertices + 3*P->VI[1];
       PV[2]=S->Vertices + 3*P->VI[2];

       /* plot edges */
       for(i=0; i<3; i++)
        fprintf(f,"%e %e %e\n",PV[i][0],PV[i][1], PV[i][2]);
       fprintf(f,"%e %e %e\n",PV[0][0],PV[0][1],PV[0][2]);
       fprintf(f,"\n\n");
     };
  fclose(f);
}


/***************************************************************/
/* This is the same as WriteGPMesh except that it writes more  */
/* information to the plot file.                               */
/***************************************************************/
void RWGGeometry::WriteGPMeshPlus(const char *format, ...)
{ 
  va_list ap;
  char FileName[MAXSTR];
  FILE *f;
  RWGPanel *P;
  RWGEdge *E;
  RWGSurface *S;
  double *V, *PV[3];
  double rArea;
  char LabelFileName[200];
  int i, ns, np, ne, nv, LabelIndex;

  va_start(ap,format);
  vsnprintfEC(FileName,MAXSTR,format,ap);
  va_end(ap);
  
  /***************************************************************/
  /* plot mesh ***************************************************/
  /***************************************************************/
  f=fopen(FileName,"w");
  for(ns=0, S=Surfaces[0]; ns<NumSurfaces; S=Surfaces[++ns])
   for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
     { 
       PV[0]=S->Vertices + 3*P->VI[0];
       PV[1]=S->Vertices + 3*P->VI[1];
       PV[2]=S->Vertices + 3*P->VI[2];

       /* plot edges */
       for(i=0; i<3; i++)
        fprintf(f,"%e %e %e\n",PV[i][0],PV[i][1], PV[i][2]);
       fprintf(f,"%e %e %e\n",PV[0][0],PV[0][1], PV[0][2]);
       fprintf(f,"\n\n");
     };
  fclose(f);

  /***************************************************************/
  /* plot normals ************************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".normals");
  f=fopen(LabelFileName,"w");
  for(ns=0, S=Surfaces[0]; ns<NumSurfaces; S=Surfaces[++ns])
   for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
     { 
       fprintf(f,"%e %e %e\n",P->Centroid[0], P->Centroid[1], P->Centroid[2]);
       rArea=sqrt(P->Area);
       VecPlusEquals(P->Centroid,rArea,P->ZHat);
       fprintf(f,"%e %e %e\n",P->Centroid[0], P->Centroid[1], P->Centroid[2]);
       VecPlusEquals(P->Centroid,-rArea,P->ZHat);
       fprintf(f,"\n\n");
     };
  fclose(f);

  /***************************************************************/
  /* plot edge labels ********************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".edgelabels");
  f=fopen(LabelFileName,"w");
  LabelIndex=1;
  for(ns=0, S=Surfaces[0]; ns<NumSurfaces; S=Surfaces[++ns])
   for(ne=0, E=S->Edges[0]; ne<S->NumEdges; E=S->Edges[++ne])
    fprintf(f,"set label %i \"%i\" at %e,%e,%e\n",LabelIndex++,
               ne,E->Centroid[0], E->Centroid[1], E->Centroid[2]);
  fclose(f);

  /***************************************************************/
  /* plot panel labels *******************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".panellabels");
  f=fopen(LabelFileName,"w");
  for(ns=0, S=Surfaces[0]; ns<NumSurfaces; S=Surfaces[++ns])
   for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
    fprintf(f,"set label %i \"%i\" at %e,%e,%e tc lt 2\n",LabelIndex++,
               np, P->Centroid[0], P->Centroid[1], P->Centroid[2]);
  fclose(f);

  /***************************************************************/
  /* plot vertex labels ******************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".vertexlabels");
  f=fopen(LabelFileName,"w");
  for(ns=0, S=Surfaces[0]; ns<NumSurfaces; S=Surfaces[++ns])
   for(nv=0, V=S->Vertices; nv<S->NumVertices; nv++, V+=3)
    fprintf(f,"set label %i \"%i\" at %e,%e,%e tc lt 3\n",LabelIndex++,
               nv,V[0],V[1],V[2]);
  fclose(f);

}


/***************************************************************/
/* WritePPMesh routine for RWGSurfaces.                         */
/***************************************************************/
void RWGSurface::WritePPMesh(const char *FileName, const char *Tag, int PlotNormals)
{ 
  FILE *f;
  RWGPanel *P;
  char buffer[MAXSTR], *p;
  double *PV[3], Val;
  int np;

  /***************************************************************/
  /* attempt to open .pp file ************************************/
  /***************************************************************/
  strncpy(buffer,FileName,996);
  p=strrchr(buffer,'.');
  if ( !p || strcmp(p,".pp") )
   strcat(buffer,".pp");

  f=fopen(buffer,"a");
  if (!f) 
   { fprintf(stderr,"warning: could not open file %s \n",FileName);
     return;
   };

  fprintf(f,"View \"%s\" {\n",Tag);
  
  /***************************************************************/
  /* plot all panels on the object   *****************************/
  /***************************************************************/
  for(np=0, P=Panels[0]; np<NumPanels; P=Panels[++np])
   { 
     PV[0]=Vertices + 3*P->VI[0];
     PV[1]=Vertices + 3*P->VI[1];
     PV[2]=Vertices + 3*P->VI[2];

     Val=(double)(Index+1);

     fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                PV[0][0], PV[0][1], PV[0][2],
                PV[1][0], PV[1][1], PV[1][2],
                PV[2][0], PV[2][1], PV[2][2],
                Val,Val,Val);
     };

   
  fprintf(f,"};\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=1;\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowScale=0;\n");

  /***************************************************************/
  /* additionally plot panel normals if that was also requested **/
  /***************************************************************/
  if (PlotNormals)
   { 
     fprintf(f,"View \"%s.Normals\" {\n",Tag);
     for(np=0, P=Panels[0]; np<NumPanels; P=Panels[++np])
      { 
        Val=fmax(VecNorm(P->Centroid), 2.0*P->Area);
        fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                   P->Centroid[0], P->Centroid[1], P->Centroid[2], 
                   Val*P->ZHat[0], Val*P->ZHat[1], Val*P->ZHat[2]);
      };
     fprintf(f,"};\n");
   }; /* if (PlotNormals) */

  fclose(f);
}

/***************************************************************/
/* WritePPMeshLabels: create a GMSH post-processing view       */
/* containing text strings for an object's panel indices.      */
/***************************************************************/
#define LS_PANELINDICES        1
#define LS_INTERIOREDGEINDICES 2
#define LS_EXTERIOREDGEINDICES 4
#define LS_VERTEXINDICES       8
void RWGSurface::WritePPMeshLabels(const char *FileName,
                                  const char *Tag, 
                                  int WhichLabels)
{ 
  FILE *f;
  RWGPanel *P;
  RWGEdge *E;
  char buffer[MAXSTR], *p;
  int np, ne, nv;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  strncpy(buffer,FileName,996);
  p=strrchr(buffer,'.');
  if ( !p || strcmp(p,".pp") )
   strcat(buffer,".pp");

  f=fopen(buffer,"a");
  if (!f) 
   { fprintf(stderr,"warning: could not open file %s \n",FileName);
     return;
   };

  /***************************************************************/
  /* panel indices ***********************************************/
  /***************************************************************/
  if (WhichLabels & LS_PANELINDICES)
   {
     if (Tag)
      fprintf(f,"View \"%s_Panels\" {\n",Tag);
     else
      fprintf(f,"View \"Panels\" {\n");
     for(np=0, P=Panels[0]; np<NumPanels; P=Panels[++np])
      fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",
                 P->Centroid[0],P->Centroid[1],P->Centroid[2],P->Index);
     fprintf(f,"};\n");
     fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=0;\n");
     fprintf(f,"View[PostProcessing.NbViews-1].ShowScale=0;\n");
   };

  /***************************************************************/
  /* interior edge indices ***************************************/
  /***************************************************************/
  if (WhichLabels & LS_INTERIOREDGEINDICES)
   {
     if (Tag)
      fprintf(f,"View \"%s.InteriorEdges\" {\n",Tag);
     else
      fprintf(f,"View \"InteriorEdges\" {\n");
     for(ne=0, E=Edges[0]; ne<NumEdges; E=Edges[++ne])
      fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",
                 E->Centroid[0],E->Centroid[1],E->Centroid[2],E->Index);
     fprintf(f,"};\n");
     fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=0;\n");
     fprintf(f,"View[PostProcessing.NbViews-1].ShowScale=0;\n");
   };

  /***************************************************************/
  /* exterior edge indices ***************************************/
  /***************************************************************/
  if (WhichLabels & LS_EXTERIOREDGEINDICES)
   {
     if (Tag)
      fprintf(f,"View \"%s.ExteriorEdges\" {\n",Tag);
     else
      fprintf(f,"View \"ExteriorEdges\" {\n");
     for(ne=0, E=ExteriorEdges[0]; ne<NumExteriorEdges; E=ExteriorEdges[++ne])
      fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",
                 E->Centroid[0],E->Centroid[1],E->Centroid[2],-(E->Index+1));
     fprintf(f,"};\n");
     fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=0;\n");
     fprintf(f,"View[PostProcessing.NbViews-1].ShowScale=0;\n");
   };

  /***************************************************************/
  /* vertex indices          *************************************/
  /***************************************************************/
  if (WhichLabels & LS_VERTEXINDICES)
   {
     if (Tag)
      fprintf(f,"View \"%s.Vertices\" {\n",Tag);
     else
      fprintf(f,"View \"Vertices\" {\n");
     for(nv=0; nv<NumVertices; nv++)
      fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",
                 Vertices[3*nv+0], Vertices[3*nv+1], Vertices[3*nv+2],nv);
     fprintf(f,"};\n");
     fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=0;\n");
     fprintf(f,"View[PostProcessing.NbViews-1].ShowScale=0;\n");
   };


  fclose(f);
}

void RWGSurface::WritePPMeshLabels(const char *FileName, const char *Tag)
{ WritePPMeshLabels(FileName, Tag,   LS_PANELINDICES 
                                   | LS_INTERIOREDGEINDICES
                                   | LS_EXTERIOREDGEINDICES
                                   | LS_VERTEXINDICES );
}

/***************************************************************/
/* WriteGPMesh routine for RWGSurfaces.                        */    
/*                                                             */
/* The mesh (with normals) may be plotted with the gnuplot cmd */
/*  splot 'file.gpmsh'  w l                                    */
/* while                                                       */
/*  splot 'file.gpmsh' i 0:10000:2 w l                         */
/* will show the mesh sans normals.                            */
/***************************************************************/
void RWGSurface::WriteGPMesh(const char *format, ...)
{ 
  va_list ap;
  char FileName[MAXSTR];
  FILE *f;
  RWGPanel *P;
  RWGEdge *E;
  double *V, *PV[3];
  double rArea;
  char LabelFileName[200];
  int i, np, ne, nv, LabelIndex;

  va_start(ap,format);
  vsnprintfEC(FileName,MAXSTR,format,ap);
  va_end(ap);
  
  /***************************************************************/
  /* plot mesh ***************************************************/
  /***************************************************************/
  f=fopen(FileName,"w");
  for(np=0, P=Panels[0]; np<NumPanels; P=Panels[++np])
   { 
     PV[0]=Vertices + 3*P->VI[0];
     PV[1]=Vertices + 3*P->VI[1];
     PV[2]=Vertices + 3*P->VI[2];

     /* plot edges */
     for(i=0; i<3; i++)
      fprintf(f,"%e %e %e\n",PV[i][0],PV[i][1],PV[i][2]);
     fprintf(f,"%e %e %e\n",PV[0][0],PV[0][1],PV[0][2]);
     fprintf(f,"\n\n");
   };
  fclose(f);

  /***************************************************************/
  /* plot normals ************************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".normals");
  f=fopen(LabelFileName,"w");
  for(np=0, P=Panels[0]; np<NumPanels; P=Panels[++np])
   { 
     fprintf(f,"%e %e %e\n",P->Centroid[0], P->Centroid[1], P->Centroid[2]);
     rArea=sqrt(P->Area);
     VecPlusEquals(P->Centroid,rArea,P->ZHat);
     fprintf(f,"%e %e %e\n",P->Centroid[0], P->Centroid[1], P->Centroid[2]);
     VecPlusEquals(P->Centroid,-rArea,P->ZHat);
     fprintf(f,"\n\n");
   };
  fclose(f);

  /***************************************************************/
  /* plot edge labels ********************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".edgelabels");
  f=fopen(LabelFileName,"w");
  LabelIndex=1;
  for(ne=0, E=Edges[0]; ne<NumEdges; E=Edges[++ne])
   fprintf(f,"set label %i \"%i\" at %e,%e,%e\n",LabelIndex++,
              ne,E->Centroid[0], E->Centroid[1], E->Centroid[2]);
  fclose(f);

  /***************************************************************/
  /* plot panel labels *******************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".panellabels");
  f=fopen(LabelFileName,"w");
  for(np=0, P=Panels[0]; np<NumPanels; P=Panels[++np])
   fprintf(f,"set label %i \"%i\" at %e,%e,%e tc lt 2\n",LabelIndex++,
              np, P->Centroid[0], P->Centroid[1], P->Centroid[2]);
  fclose(f);

  /***************************************************************/
  /* plot vertex labels ******************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".vertexlabels");
  f=fopen(LabelFileName,"w");
  for(nv=0, V=Vertices; nv<NumVertices; nv++)
   fprintf(f,"set label %i \"%i\" at %e,%e,%e tc lt 3\n",LabelIndex++,
              nv,V[3*nv],V[3*nv+1],V[3*nv+2]);
  fclose(f);

}

/***************************************************************/
/* Emit GMSH postprocessing code for visualizing the current   */
/* distribution described by a single vector of surface-current*/
/* expansion coefficients.                                     */
/***************************************************************/
void RWGGeometry::PlotSurfaceCurrents(const char *SurfaceLabel,
                                      HVector *KN, cdouble Omega, 
                                      const char *format, ...)
{ 
  int ns, np, ne;
  FILE *f;
  int MM, Offset, NeedMagnetic;
  RWGSurface *S, *WhichSurface;
  RWGPanel *P;
  RWGEdge *E;
  double *PV[3];
  double XmQ[3];
  cdouble Weight, Rho, J[3];

  char FileName[MAXSTR];
  va_list ap;

  cdouble iw = II*Omega;

  if (SurfaceLabel)
   WhichSurface=GetSurfaceByLabel(SurfaceLabel);
  else
   WhichSurface=NULL;
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  va_start(ap,format);
  vsnprintfEC(FileName,MAXSTR,format,ap);
  va_end(ap);
  f=fopen(FileName,"w");
  if (!f) return;

  /***************************************************************/
  /* plot electric charge densities of all basis functions       */
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","Electric Charge Density");
  for(ns=0; ns<NumSurfaces; ns++)
   { 
     S=Surfaces[ns];
     if (WhichSurface && S!=WhichSurface) 
      continue;

     for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
      { 
        MM=S->IsPEC ? 1 : 2;
        Offset=BFIndexOffset[ns];

        /* */
        Rho=0.0;
        for(ne=0; ne<S->NumEdges; ne++)
         if ( S->Edges[ne]->iPPanel == np )
          Rho += S->Edges[ne]->Length * KN->GetEntry(Offset + MM*ne);
         else if ( S->Edges[ne]->iMPanel == np )
          Rho -= S->Edges[ne]->Length * KN->GetEntry(Offset + MM*ne);
        Rho /= (P->Area * iw);

        PV[0]=S->Vertices + 3*P->VI[0];
        PV[1]=S->Vertices + 3*P->VI[1];
        PV[2]=S->Vertices + 3*P->VI[2];

        fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                   PV[0][0], PV[0][1], PV[0][2],
                   PV[1][0], PV[1][1], PV[1][2],
                   PV[2][0], PV[2][1], PV[2][2],
                   real(Rho),real(Rho),real(Rho));
      };
   }; 
  fprintf(f,"};\n");

  /***************************************************************/
  /* plot electric current densities at centroids of all panels  */
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","Electric Current");
  for(ns=0; ns<NumSurfaces; ns++)
   { 
     S=Surfaces[ns];
     if (WhichSurface && S!=WhichSurface) 
      continue;

     for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
      { 
        MM=S->IsPEC ? 1 : 2;
        Offset=BFIndexOffset[ns];
 
        memset(J,0,3*sizeof(cdouble));
        for(ne=0, E=S->Edges[0]; ne<S->NumEdges; E=S->Edges[++ne])
         { if ( E->iPPanel == np )
            { VecSub( P->Centroid, S->Vertices + 3*(E->iQP), XmQ);
              Weight = E->Length * KN->GetEntry(Offset + MM*ne);
              J[0] += Weight * XmQ[0];
              J[1] += Weight * XmQ[1];
              J[2] += Weight * XmQ[2];
            }
           else if ( E->iMPanel == np )
            { VecSub( P->Centroid, S->Vertices + 3*(E->iQM), XmQ);
              Weight = E->Length * KN->GetEntry(Offset + MM*ne);
              J[0] -= Weight * XmQ[0];
              J[1] -= Weight * XmQ[1];
              J[2] -= Weight * XmQ[2];
            };
         };
        J[0] /= (2.0*P->Area);
        J[1] /= (2.0*P->Area);
        J[2] /= (2.0*P->Area);
 
        fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                   P->Centroid[0],P->Centroid[1],P->Centroid[2],
                   abs(J[0]),abs(J[1]),abs(J[2]));
 
      };
   };
  fprintf(f,"};\n");

  /***************************************************************/
  /* we only need magnetic charge/current displays iff any       */
  /* object is non-PEC                                           */
  /***************************************************************/
  NeedMagnetic=0;
  for(ns=0; NeedMagnetic==0 && ns<NumSurfaces; ns++)
   { 
     S=Surfaces[ns];
     if (WhichSurface && S!=WhichSurface)
      continue;
     if ( S->IsPEC )
      continue;
     NeedMagnetic=1;
   };

  if (NeedMagnetic==0)
   { fclose(f); 
     return;
   };

  /***************************************************************/
  /* plot magnetic charge densities of all basis functions       */
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","Magnetic Charge Density");
  for(ns=0; ns<NumSurfaces; ns++)
   { 
     S=Surfaces[ns];
     if (WhichSurface && S!=WhichSurface) 
      continue;
     if ( S->IsPEC ) 
      continue;
     MM=2;
     Offset=BFIndexOffset[ns];

     for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
       { 
         Rho=0.0;
         for(ne=0; ne<S->NumEdges; ne++)
          if ( S->Edges[ne]->iPPanel == np )
           Rho += S->Edges[ne]->Length * KN->GetEntry(Offset + MM*ne + 1);
          else if ( S->Edges[ne]->iMPanel == np )
           Rho -= S->Edges[ne]->Length * KN->GetEntry(Offset + MM*ne + 1);
         Rho /= (P->Area * iw);
         Rho*=(-1.0*ZVAC);

         PV[0]=S->Vertices + 3*P->VI[0];
         PV[1]=S->Vertices + 3*P->VI[1];
         PV[2]=S->Vertices + 3*P->VI[2];

         fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                    PV[0][0], PV[0][1], PV[0][2],
                    PV[1][0], PV[1][1], PV[1][2],
                    PV[2][0], PV[2][1], PV[2][2],
                    real(Rho),real(Rho),real(Rho));
        };
   };
  fprintf(f,"};\n");

  /***************************************************************/
  /* plot electric current densities at centroids of all panels  */
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","Magnetic Current");
  for(ns=0; ns<NumSurfaces; ns++)
   { 
     S=Surfaces[ns];
     if (WhichSurface && S!=WhichSurface) 
      continue;
     if ( S->IsPEC )
      continue;
     MM=2;
     Offset=BFIndexOffset[ns];

     for(np=0, P=S->Panels[0]; np<S->NumPanels; P=S->Panels[++np])
      {
        /* */
        memset(J,0,3*sizeof(cdouble));
        for(ne=0, E=S->Edges[0]; ne<S->NumEdges; E=S->Edges[++ne])
         { if ( E->iPPanel == np )
            { VecSub( P->Centroid, S->Vertices + 3*(E->iQP), XmQ);
              Weight = E->Length * KN->GetEntry(Offset + MM*ne + 1);
              Weight *= -1.0*ZVAC;
              J[0] += Weight * XmQ[0];
              J[1] += Weight * XmQ[1];
              J[2] += Weight * XmQ[2];
            }
           else if ( E->iMPanel == np )
            { VecSub( P->Centroid, S->Vertices + 3*(E->iQM), XmQ);
              Weight = E->Length * KN->GetEntry(Offset + MM*ne + 1);
              Weight *= -1.0*ZVAC;
              J[0] -= Weight * XmQ[0];
              J[1] -= Weight * XmQ[1];
              J[2] -= Weight * XmQ[2];
            };
         };
        J[0] /= (2.0*P->Area);
        J[1] /= (2.0*P->Area);
        J[2] /= (2.0*P->Area);
 
        fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                   P->Centroid[0],P->Centroid[1],P->Centroid[2],
                   abs(J[0]),abs(J[1]),abs(J[2]));

      };
   };
  fprintf(f,"};\n");
  fclose(f);

}

/***************************************************************/
/* an alternative entry point to PlotSurfaceCurrents in which  */
/* the called doesn't specify a surface label; in this case    */
/* currents on all surfaces are plotted.                       */
/***************************************************************/
void RWGGeometry::PlotSurfaceCurrents(HVector *KN, cdouble Omega, 
                                      const char *format, ...)
{ 
  va_list ap;
  char FileName[MAXSTR];
  va_start(ap,format);
  vsnprintfEC(FileName,MAXSTR,format,ap);
  va_end(ap);
  
  PlotSurfaceCurrents(0, KN, Omega, FileName);
}

} // namespace scuff
