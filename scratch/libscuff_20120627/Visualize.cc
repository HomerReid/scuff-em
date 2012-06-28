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
  RWGObject *O;
  RWGPanel *P;
  char buffer[1000], *p;
  double *PV[3], Val;
  int no, np;

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
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
    { 
      PV[0]=O->Vertices + 3*P->VI[0];
      PV[1]=O->Vertices + 3*P->VI[1];
      PV[2]=O->Vertices + 3*P->VI[2];
      Val=(double)(no+1);
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
     for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
      for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
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

void RWGGeometry::WritePPMesh(const char *FileName, const char *Tag)
 { WritePPMesh(FileName, Tag, 0); }

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
  char FileName[1000], *p;
  FILE *f;
  RWGObject *O;
  RWGPanel *P;
  double *PV[3];
  int i, no, np;

  va_start(ap,format);
  vsnprintf(FileName,997,format,ap);
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
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
     { 
       PV[0]=O->Vertices + 3*P->VI[0];
       PV[1]=O->Vertices + 3*P->VI[1];
       PV[2]=O->Vertices + 3*P->VI[2];

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
  char FileName[1000];
  FILE *f;
  RWGPanel *P;
  RWGEdge *E;
  RWGObject *O;
  double *V, *PV[3];
  double rArea;
  char LabelFileName[200];
  int i, no, np, ne, nv, LabelIndex;

  va_start(ap,format);
  vsnprintf(FileName,1000,format,ap);
  va_end(ap);
  
  /***************************************************************/
  /* plot mesh ***************************************************/
  /***************************************************************/
  f=fopen(FileName,"w");
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
     { 
       PV[0]=O->Vertices + 3*P->VI[0];
       PV[1]=O->Vertices + 3*P->VI[1];
       PV[2]=O->Vertices + 3*P->VI[2];

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
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
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
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(ne=0, E=O->Edges[0]; ne<O->NumEdges; E=O->Edges[++ne])
    fprintf(f,"set label %i \"%i\" at %e,%e,%e\n",LabelIndex++,
               ne,E->Centroid[0], E->Centroid[1], E->Centroid[2]);
  fclose(f);

  /***************************************************************/
  /* plot panel labels *******************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".panellabels");
  f=fopen(LabelFileName,"w");
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
    fprintf(f,"set label %i \"%i\" at %e,%e,%e tc lt 2\n",LabelIndex++,
               np, P->Centroid[0], P->Centroid[1], P->Centroid[2]);
  fclose(f);

  /***************************************************************/
  /* plot vertex labels ******************************************/
  /***************************************************************/
  strncpy(LabelFileName,FileName,180);
  strcat(LabelFileName,".vertexlabels");
  f=fopen(LabelFileName,"w");
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(nv=0, V=O->Vertices; nv<O->NumVertices; nv++, V+=3)
    fprintf(f,"set label %i \"%i\" at %e,%e,%e tc lt 3\n",LabelIndex++,
               nv,V[0],V[1],V[2]);
  fclose(f);

}


/***************************************************************/
/* WritePPMesh routine for RWGObjects.                         */
/***************************************************************/
void RWGObject::WritePPMesh(const char *FileName, const char *Tag, int PlotNormals)
{ 
  FILE *f;
  RWGPanel *P;
  char buffer[1000], *p;
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

void RWGObject::WritePPMesh(const char *FileName, const char *Tag)
 { WritePPMesh(FileName, Tag, 0); }

/***************************************************************/
/* WritePPMeshLabels: create a GMSH post-processing view       */
/* containing text strings for an object's panel indices.      */
/***************************************************************/
#define LS_PANELINDICES        1
#define LS_INTERIOREDGEINDICES 2
#define LS_EXTERIOREDGEINDICES 4
#define LS_VERTEXINDICES       8
void RWGObject::WritePPMeshLabels(const char *FileName,
                                  const char *Tag, 
                                  int WhichLabels)
{ 
  FILE *f;
  RWGPanel *P;
  RWGEdge *E;
  char buffer[1000], *p;
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

void RWGObject::WritePPMeshLabels(const char *FileName, const char *Tag)
{ WritePPMeshLabels(FileName, Tag,   LS_PANELINDICES 
                                   | LS_INTERIOREDGEINDICES
                                   | LS_EXTERIOREDGEINDICES
                                   | LS_VERTEXINDICES );
}

/***************************************************************/
/* WriteGPMesh routine for RWGObjects.                        */    
/*                                                             */
/* The mesh (with normals) may be plotted with the gnuplot cmd */
/*  splot 'file.gpmsh'  w l                                    */
/* while                                                       */
/*  splot 'file.gpmsh' i 0:10000:2 w l                         */
/* will show the mesh sans normals.                            */
/***************************************************************/
void RWGObject::WriteGPMesh(const char *format, ...)
{ 
  va_list ap;
  char FileName[1000];
  FILE *f;
  RWGPanel *P;
  RWGEdge *E;
  double *V, *PV[3];
  double rArea;
  char LabelFileName[200];
  int i, np, ne, nv, LabelIndex;

  va_start(ap,format);
  vsnprintf(FileName,1000,format,ap);
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
/***************************************************************/
/* Experimental (10/19/09) code for visualizing the current    */
/* distribution described by a single vector.                  */
/***************************************************************/
/***************************************************************/
#if 0
void RWGGeometry::PlotVector(double *KVec, const char *format, ...)
{ 
  int i, no, np, nv, ne, ei, na, NumArrows;
  int Offset;
  RWGObject *O;
  RWGEdge *E;
  RWGPanel *P;
  double J[3];
  char FileBase[1000], *p;
  va_list ap;
  void *pCC;
  double *x, *y, *z, *u, *v, *w, *X, *Y, *Z, *Tri, *C;
 
  double *x2, *y2, *z2, *u2, *v2, *w2;
  int NumArrows2, VertexIndexOffset[10];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  va_start(ap,format);
  vsnprintf(FileBase,1000,format,ap);
  va_end(ap);
  pCC=C2MLOpen(FileBase);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  NumArrows=TotalPanels;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  NumArrows2=Objects[0]->NumVertices;
  VertexIndexOffset[0]=0;
  for(no=1; no<NumObjects; no++) 
   { NumArrows2 += Objects[no]->NumVertices;
     VertexIndexOffset[no]=VertexIndexOffset[no-1]+Objects[no-1]->NumVertices;
   };
 
  /***************************************************************/
  /* allocate space for vectors exported to matlab              **/
  /***************************************************************/
  x=(double *)mallocEC(NumArrows*sizeof(double));
  y=(double *)mallocEC(NumArrows*sizeof(double));
  z=(double *)mallocEC(NumArrows*sizeof(double));
  u=(double *)mallocEC(NumArrows*sizeof(double));
  v=(double *)mallocEC(NumArrows*sizeof(double));
  w=(double *)mallocEC(NumArrows*sizeof(double));

  x2=(double *)mallocEC(NumArrows2*sizeof(double));
  y2=(double *)mallocEC(NumArrows2*sizeof(double));
  z2=(double *)mallocEC(NumArrows2*sizeof(double));
  u2=(double *)mallocEC(NumArrows2*sizeof(double));
  v2=(double *)mallocEC(NumArrows2*sizeof(double));
  w2=(double *)mallocEC(NumArrows2*sizeof(double));

  X=(double *)mallocEC(3*NumArrows*sizeof(double));
  Y=(double *)mallocEC(3*NumArrows*sizeof(double));
  Z=(double *)mallocEC(3*NumArrows*sizeof(double));
  Tri=(double *)mallocEC(3*NumArrows*sizeof(double));
  C=(double *)mallocEC(NumArrows*sizeof(double));

  /***************************************************************/
  /* first pass to fill in x,y,z, X, Y, Z, Tri vectors ***********/
  /***************************************************************/
  for(no=0; no<NumObjects; no++)
   {
     O=Objects[no];

     for(np=0; np<O->NumPanels; np++)
      { 
        na=PanelIndexOffset[no] + np;

        x[na]=O->Panels[np]->Centroid[0];
        y[na]=O->Panels[np]->Centroid[1];
        z[na]=O->Panels[np]->Centroid[2];

        for(i=0; i<3; i++)
         { X[ 3*na + i ] = O->Vertices[ 3*(O->Panels[np]->VI[i]) + 0 ];
           Y[ 3*na + i ] = O->Vertices[ 3*(O->Panels[np]->VI[i]) + 1 ];
           Z[ 3*na + i ] = O->Vertices[ 3*(O->Panels[np]->VI[i]) + 2 ];
           Tri[ na + i*NumArrows ]= (double ) 3*na + i + 1;
         };

      };

     for(nv=0; nv<O->NumVertices; nv++)
      { 
        na=VertexIndexOffset[no] + nv;
        x2[na]=O->Vertices[3*nv];
        y2[na]=O->Vertices[3*nv + 1];
        z2[na]=O->Vertices[3*nv + 2];
      };
   };
  
  /***************************************************************/
  /* second pass to fill in u, v, w, c vectors  ******************/
  /***************************************************************/
  memset(u,0,NumArrows*sizeof(double));
  memset(v,0,NumArrows*sizeof(double));
  memset(w,0,NumArrows*sizeof(double));
  memset(C,0,NumArrows*sizeof(double));

  for(ei=0; ei<TotalEdges; ei++)
   { 
     ne=GetObjectAndEdgeIndex(ei, &O);
     no=O->Index;
     E=O->Edges[ne];

     /*--------------------------------------------------------------*/
     /*- compute contribution of this edge to the arrows at the     -*/
     /*- centroids of the two panels to which it belongs            -*/
     /*--------------------------------------------------------------*/
     P=O->Panels[E->iPPanel];
     na=PanelIndexOffset[no] + P->Index;
     VecSub(P->Centroid, O->Vertices + 3*(E->iQP), J);
     VecScale(J, KVec[ei] * E->Length / (2.0*P->Area) );
     u[na] += J[0];
     v[na] += J[1];
     w[na] += J[2];
     C[na] += KVec[ei] * E->Length / ( P->Area );

     P=O->Panels[E->iMPanel];
     na=PanelIndexOffset[no] + P->Index;
     VecSub(P->Centroid, O->Vertices + 3*(E->iQM), J);
     VecScale(J, KVec[ei] * E->Length / (2.0*P->Area) );
     u[na] -= J[0];
     v[na] -= J[1];
     w[na] -= J[2];
     C[na] -= KVec[ei] * E->Length / ( P->Area );

     /*--------------------------------------------------------------*/
     /*- compute contribution of this edge to the arrows at its     -*/
     /*- two endpoints                                              -*/
     /*--------------------------------------------------------------*/
     P=O->Panels[E->iPPanel];
     na=VertexIndexOffset[no] + E->iV1;
     VecSub(O->Vertices+3*(E->iV1), O->Vertices + 3*(E->iQP), J);
     VecScale(J, KVec[ei] * E->Length / (2.0*P->Area) );
     u2[na] += J[0];
     v2[na] += J[1];
     w2[na] += J[2];

     na=VertexIndexOffset[no] + E->iV2;
     VecSub(O->Vertices+3*(E->iV2), O->Vertices + 3*(E->iQP), J);
     VecScale(J, KVec[ei] * E->Length / (2.0*P->Area) );
     u2[na] += J[0];
     v2[na] += J[1];
     w2[na] += J[2];
     
   };

  /***************************************************************/
  /* write quiver-plot vectors to hd5 file ***********************/
  /***************************************************************/
  C2MLVector(pCC,x,NumArrows,"x");
  C2MLVector(pCC,y,NumArrows,"y");
  C2MLVector(pCC,z,NumArrows,"z");
  C2MLVector(pCC,u,NumArrows,"u");
  C2MLVector(pCC,v,NumArrows,"v");
  C2MLVector(pCC,w,NumArrows,"w");

  C2MLVector(pCC,X,3*NumArrows,"X");
  C2MLVector(pCC,Y,3*NumArrows,"Y");
  C2MLVector(pCC,Z,3*NumArrows,"Z");
  C2MLMatrix_RM(pCC,Tri,NumArrows,3,"Tri",no);
  C2MLVector(pCC,C,NumArrows,"C");

  C2MLVector(pCC,x2,NumArrows2,"x2");
  C2MLVector(pCC,y2,NumArrows2,"y2");
  C2MLVector(pCC,z2,NumArrows2,"z2");
  C2MLVector(pCC,u2,NumArrows2,"u2");
  C2MLVector(pCC,v2,NumArrows2,"v2");
  C2MLVector(pCC,w2,NumArrows2,"w2");


  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  free(x);
  free(y);
  free(z);
  free(u);
  free(v);
  free(w);

  free(x2);
  free(y2);
  free(z2);
  free(u2);
  free(v2);
  free(w2);

  free(X);
  free(Y);
  free(Z);
  free(Tri);
  free(C);

  C2MLClose(pCC); 

}
#endif

/***************************************************************/
/***************************************************************/
/* Experimental (12/16/09) code for visualizing the current    */
/* distribution described by a single vector.                  */
/***************************************************************/
/***************************************************************/
void RWGGeometry::PlotSurfaceCurrents(HVector *KN, double Frequency, int RealFreq,
                                      const char *format, ...)
{ 
  int no, np, ne;
  FILE *f;
  int MM, Offset, NeedMagnetic;
  RWGObject *O;
  RWGPanel *P;
  RWGEdge *E;
  double *PV[3];
  double XmQ[3];
  cdouble Weight, Rho, J[3];

  char FileName[1000];
  va_list ap;

  cdouble iw = RealFreq ? II*Frequency : -1.0*Frequency;
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  va_start(ap,format);
  vsnprintf(FileName,1000,format,ap);
  va_end(ap);
  f=fopen(FileName,"w");
  if (!f) return;

  /***************************************************************/
  /* plot electric charge densities of all basis functions       */
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","Electric Charge Density");
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
     { 
       MM=O->MP->IsPEC() ? 1 : 2;
       Offset=BFIndexOffset[no];

       /* */
       Rho=0.0;
       for(ne=0; ne<O->NumEdges; ne++)
        if ( O->Edges[ne]->iPPanel == np )
         Rho += O->Edges[ne]->Length * KN->GetEntry(Offset + MM*ne);
        else if ( O->Edges[ne]->iMPanel == np )
         Rho -= O->Edges[ne]->Length * KN->GetEntry(Offset + MM*ne);
       Rho /= (P->Area * iw);

       PV[0]=O->Vertices + 3*P->VI[0];
       PV[1]=O->Vertices + 3*P->VI[1];
       PV[2]=O->Vertices + 3*P->VI[2];

       fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                  PV[0][0], PV[0][1], PV[0][2],
                  PV[1][0], PV[1][1], PV[1][2],
                  PV[2][0], PV[2][1], PV[2][2],
                  real(Rho),real(Rho),real(Rho));
     };
  fprintf(f,"};\n");

  /***************************************************************/
  /* plot electric current densities at centroids of all panels  */
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","Electric Current");
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
     { 
       MM=O->MP->IsPEC() ? 1 : 2;
       Offset=BFIndexOffset[no];

       /* */
       memset(J,0,3*sizeof(cdouble));
       for(ne=0, E=O->Edges[0]; ne<O->NumEdges; E=O->Edges[++ne])
        { if ( E->iPPanel == np )
           { VecSub( P->Centroid, O->Vertices + 3*(E->iQP), XmQ);
             Weight = E->Length * KN->GetEntry(Offset + MM*ne);
             J[0] += Weight * XmQ[0];
             J[1] += Weight * XmQ[1];
             J[2] += Weight * XmQ[2];
           }
          else if ( E->iMPanel == np )
           { VecSub( P->Centroid, O->Vertices + 3*(E->iQM), XmQ);
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
                  real(J[0]),real(J[1]),real(J[2]));

     };
  fprintf(f,"};\n");

  /***************************************************************/
  /* we only need magnetic charge/current displays iff any       */
  /* object is non-PEC                                           */
  /***************************************************************/
  NeedMagnetic=0;
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   if ( !(O->MP->IsPEC()) )
    NeedMagnetic=1;
  if (NeedMagnetic==0)
   { fclose(f); 
     return;
   };

  /***************************************************************/
  /* plot magnetic charge densities of all basis functions       */
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","Magnetic Charge Density");
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   { 
     if ( O->MP->IsPEC() ) 
      continue;
     MM=2;
     Offset=BFIndexOffset[no];

     for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
       { 
         Rho=0.0;
         for(ne=0; ne<O->NumEdges; ne++)
          if ( O->Edges[ne]->iPPanel == np )
           Rho += O->Edges[ne]->Length * KN->GetEntry(Offset + MM*ne + 1);
          else if ( O->Edges[ne]->iMPanel == np )
           Rho -= O->Edges[ne]->Length * KN->GetEntry(Offset + MM*ne + 1);
         Rho /= (P->Area * iw);
         Rho*=(-1.0*ZVAC);

         PV[0]=O->Vertices + 3*P->VI[0];
         PV[1]=O->Vertices + 3*P->VI[1];
         PV[2]=O->Vertices + 3*P->VI[2];

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
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   { 
     if ( O->MP->IsPEC() ) 
      continue;
     MM=2;
     Offset=BFIndexOffset[no];

     for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
      {
        /* */
        memset(J,0,3*sizeof(cdouble));
        for(ne=0, E=O->Edges[0]; ne<O->NumEdges; E=O->Edges[++ne])
         { if ( E->iPPanel == np )
            { VecSub( P->Centroid, O->Vertices + 3*(E->iQP), XmQ);
              Weight = E->Length * KN->GetEntry(Offset + MM*ne + 1);
              Weight *= -1.0*ZVAC;
              J[0] += Weight * XmQ[0];
              J[1] += Weight * XmQ[1];
              J[2] += Weight * XmQ[2];
            }
           else if ( E->iMPanel == np )
            { VecSub( P->Centroid, O->Vertices + 3*(E->iQM), XmQ);
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
                   real(J[0]),real(J[1]),real(J[2]));

      };
   };
  fprintf(f,"};\n");
  fclose(f);

}

} // namespace scuff
