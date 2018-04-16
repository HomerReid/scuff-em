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

#define COMPLETE_VARARGS(Format, Buffer)  \
  va_list ap;                             \
  char Buffer[MAXSTR];                    \
  va_start(ap,Format);                    \
  vsnprintfEC(Buffer,MAXSTR,Format,ap);   \
  va_end(ap);                             \

/************************************************************/
/* subroutines for emitting GMSH postprocessing code        */
/************************************************************/

/* vector point (otherwise known as 'arrow') */
void WriteVP(double *X, double Vx, double Vy, double Vz, FILE *f, double Scale=1.0)
{
  fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",X[0],X[1],X[2],Scale*Vx,Scale*Vy,Scale*Vz);
}

void WriteVP(double *X, double *V, FILE *f, double Scale=1.0)
{ WriteVP(X, V[0], V[1], V[2], f, Scale);
}

void WriteVP(RWGSurface *S, int np, double Vx, double Vy, double Vz, FILE *f, double Scale=1.0)
{ WriteVP(S->Panels[np]->Centroid, Vx, Vy, Vz, f, Scale); }

void WriteVP(RWGSurface *S, int np, double *V, FILE *f, double Scale=1.0)
{ WriteVP(S->Panels[np]->Centroid, V[0], V[1], V[2], f, Scale); }

/* scalar triangle */
void WriteST(double **VV, double Vals[3], FILE *f)
{
  fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
             VV[0][0], VV[0][1], VV[0][2],
             VV[1][0], VV[1][1], VV[1][2],
             VV[2][0], VV[2][1], VV[2][2], 
             Vals[0],Vals[1],Vals[2]);
}

void WriteST(double **VV, double Val, FILE *f)
{ double Vals[3];
  Vals[0]=Vals[1]=Vals[2]=Val;
  WriteST(VV, Vals, f);
}

void WriteST(RWGSurface *S, int np, double Vals[3], FILE *f)
{ double *VV[3];
  VV[0] = S->Vertices + 3*S->Panels[np]->VI[0];
  VV[1] = S->Vertices + 3*S->Panels[np]->VI[1];
  VV[2] = S->Vertices + 3*S->Panels[np]->VI[2];
  WriteST(VV, Vals, f);
}

void WriteST(RWGSurface *S, int np, double Val, FILE *f)
{ double Vals[3];
  Vals[0]=Vals[1]=Vals[2]=Val;
  WriteST(S, np, Vals, f);
}
  
/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::WritePPMesh(const char *FileName, const char *Tag, bool PlotNormals)
{  
  for(int ns=0; ns<NumSurfaces; ns++)
   Surfaces[ns]->WritePPMesh(FileName, Tag, PlotNormals);
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
/* WritePPMesh routine. Writes geometry to a .pp file suitable */
/* for opening as a GMSH postprocessing file.                  */
/* Note calling convention and file handling are different     */
/* from those of WriteGPMesh.                                  */
/* If PlotNormals is nonzero, we also plot the normal to each  */
/* panel.                                                      */
/***************************************************************/
void RWGSurface::WritePPMesh(const char *FileName, const char *Tag, bool PlotNormals)
{ 
  /***************************************************************/
  /* attempt to open .pp file ************************************/
  /***************************************************************/
  char buffer[MAXSTR];
  strncpy(buffer,FileName,MAXSTR-4);
  char *p=strrchr(buffer,'.');
  if ( !p || strcmp(p,".pp") )
   strcat(buffer,".pp");

  FILE *f=fopen(buffer,"a");
  if (!f) 
   { fprintf(stderr,"warning: could not open file %s \n",FileName);
     return;
   };

  char ViewName[MAXSTR];
  if (Tag)
   snprintf(ViewName,MAXSTR,"%s(%s)",Label,Tag);
  else 
   snprintf(ViewName,MAXSTR,"%s",Label);

  /***************************************************************/
  /* plot all panels on the object   *****************************/
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n",ViewName);
  for(int np=0; np<NumPanels; np++)
   WriteST(this, np, (double)(Index+1), f);
  fprintf(f,"};\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=1;\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowScale=0;\n");

  /***************************************************************/
  /* additionally plot panel normals if requested                */
  /***************************************************************/
  if (PlotNormals)
   { fprintf(f,"View \"%s.Normals\" {\n",ViewName);
     for(int np=0; np<NumPanels; np++)
      { RWGPanel *P=Panels[np];
        double Scale=fmax(VecNorm(P->Centroid), 2.0*P->Area);
        WriteVP(this, np, P->ZHat, f, Scale);
      }
     fprintf(f,"};\n");
   }

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
#define LS_BFDIRECTIONS        16
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
     if (NumExteriorEdges>0) 
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

  /***************************************************************/
  /* basis-function directions ***********************************/
  /***************************************************************/
  if (WhichLabels & LS_BFDIRECTIONS)
   {
     if (Tag)
      fprintf(f,"View \"%s.BFDirections\" {\n",Tag);
     else
      fprintf(f,"View \"BFDirections\" {\n");
     for(ne=0, E=Edges[0]; ne<NumEdges; E=Edges[++ne])
      { 
        double *PCentroid, *MCentroid;
        PCentroid = Panels[E->iPPanel]->Centroid;
        if (E->iMPanel==-1)
         MCentroid = E->Centroid;
        else
         MCentroid = Panels[E->iMPanel]->Centroid;

        double Dir[3]; 
        Dir[0] = MCentroid[0] - PCentroid[0];
        Dir[1] = MCentroid[1] - PCentroid[1];
        Dir[2] = MCentroid[2] - PCentroid[2];
        fprintf(f,"VP(%e,%e,%e) {%e, %e, %e};\n",
                 E->Centroid[0],E->Centroid[1],E->Centroid[2],
                 Dir[0],Dir[1],Dir[2]);
      };
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
                                   | LS_VERTEXINDICES
                                   | LS_BFDIRECTIONS
                   );
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
/***************************************************************/
/***************************************************************/
double GetTriangleArea(double *V1, double *V2, double *V3)
{ 
  double A[3], B[3], AxB[3];
  return 0.5*VecNorm(VecCross( VecSub(V2, V1, A), VecSub(V3, V1, B), AxB ));
}

/***************************************************************/
/* Plot a function defined on RWG edges or panels.             */
/*                                                             */
/* If ByEdge==true, Values[ne] is the value associated with    */
/* edge #ne. In this case, the first NumEdges elements of the  */
/* Values array are accessed.                                  */
/*                                                             */
/* If ByEdge==false, Values[np] is the value associated with   */
/* panel #np. In this case, the first NP elements of Values    */
/* are accessed.                                               */
/*                                                             */
/* How it works: Suppose we have a quantity Q defined as a     */
/* sum of contributions from individual edges:                 */
/*                                                             */
/*  Q = \sum_{edges e} Q_e                                     */
/*                                                             */
/* where Q_e is the contribution of edge #e. The RWG basis     */
/* function associated with edge #e is supported on one or two */
/* triangles; letting A_e be the total area of those triangles,*/
/* and imagining that the contribution Q_e is evenly           */
/* distributed over that area, we obtain an effective *surface */
/* density* \rho_e = Q_e/A_e of the quantity Q associated with */
/* the edge.                                                   */
/***************************************************************/
double GetScalarValue(void *Values, int PlotType, int nv)
{ double *dValues=(double *)Values;
  cdouble *zValues=(cdouble *)Values;
  switch(PlotType)
   { case 0: return dValues[nv];
     case 1: return abs(zValues[nv]);
     case 2: return real(zValues[nv]);
     case 3: 
     default: return imag(zValues[nv]);
   }
  return 0.0;
}

/*******************************************************************/
/* Plot a scalar quantity defined on panels or edges.              */
/*                                                                 */
/* PlotType==0 --> Values is an array of doubles                   */
/* PlotType==1 --> Values is an array of cdoubles, plot magnitude  */
/* PlotType==2 --> Values is an array of cdoubles, plot real part  */
/* PlotType==3 --> Values is an array of cdoubles, plot imag part  */
/*                                                                 */
/* AveragingMethod = 0: Values is an array of length S->NumPanels; */
/*                      for each panel simply plot the correspondng*/
/*                      entry of Values (i.e. no averaging done.)  */
/* AveragingMethod = 1: Values is an array of length S->NumPanels; */
/*                      for each panel vertex plot the average     */
/*                      value of all panels sharing that vertex.   */
/* AveragingMethod = 2: Values is an array of length S->NumEdges;  */
/*                      for each panel vertex plot the average     */
/*                      value of all edges sharing that vertex.    */
/*                                                                 */
/* The return value is the surface integral of the quantity over   */
/* the surface.                                                    */
/*******************************************************************/
double RWGSurface::PlotScalarDensity(int PlotType, void *Values,
                                     int AveragingMethod,
                                     const char *UserFileName,
                                     const char *ViewNameFormat, ...)
{ 
  char FileName[MAXSTR];
  strncpy(FileName,UserFileName,MAXSTR-4);
  char *p=strrchr(FileName,'.');
  if ( !p || strcmp(p,".pp") )
   strcat(FileName,".pp");

  FILE *f=fopen(FileName,"a");
  if (!f) 
   { fprintf(stderr,"warning: could not open file %s \n",FileName);
     return 0.0;
   };
  
  COMPLETE_VARARGS(ViewNameFormat, ViewName);

  /***************************************************************/
  /* For each panel vertex, average the contributions of all     */
  /* edges (if ByEdge==true) or of all panels (if ByPanel==false)*/
  /* that share that vertex.                                     */
  /***************************************************************/
  double *ValuePerVertex = 0, *AreaPerVertex = 0;
     int *NumPerVertex   = 0;

  if (AveragingMethod!=0)
   { ValuePerVertex = new double [NumVertices];
     AreaPerVertex  = new double [NumVertices];
     NumPerVertex   = new int    [NumVertices];
     memset(ValuePerVertex, 0, NumVertices*sizeof(double));
     memset(AreaPerVertex,  0, NumVertices*sizeof(double));
     memset(NumPerVertex,   0, NumVertices*sizeof(int));
   }
  if (AveragingMethod==1)
   for(int np=0; np<NumPanels; np++)
    { RWGPanel *P = Panels[np];
      double Area = P->Area;
      double Val  = GetScalarValue(Values, PlotType, np);
      for(int i=0; i<3; i++)
       { ValuePerVertex[P->VI[i]] += Area*Val;
         AreaPerVertex[P->VI[i]]  += Area;
         NumPerVertex[P->VI[i]]    = 1;
       }
    }
  else if (AveragingMethod==2)
   for(int ne=0; ne<NumEdges; ne++)
    { RWGEdge *E = Edges[ne];
      double Area = Panels[E->iPPanel]->Area;
      if (E->iMPanel!=-1) Area+=Panels[E->iMPanel]->Area;
      double Val = GetScalarValue(Values, PlotType, ne);
      ValuePerVertex[E->iV1] += Area*Val;
      ValuePerVertex[E->iV2] += Area*Val;
      AreaPerVertex[E->iV1]  += Area;
      AreaPerVertex[E->iV2]  += Area;
      NumPerVertex[E->iV1]   += 1;
      NumPerVertex[E->iV2]   += 1;
    }
   
  double Integral=0.0;
  fprintf(f,"View \"%s\" {\n",ViewName);
  for(int np=0; np<NumPanels; np++)
   if (AveragingMethod==0)
    { double Value=GetScalarValue(Values, PlotType, np);
      WriteST(this, np, Value, f);
      Integral += Value * Panels[np]->Area;
    }
   else
    { int *VI = Panels[np]->VI;
      double VVals[3];
      for(int i=0; i<3; i++)
       VVals[i] = ValuePerVertex[VI[i]] / (NumPerVertex[VI[i]]*AreaPerVertex[VI[i]]);
      WriteST(this, np, VVals, f);
      Integral += Panels[np]->Area * (VVals[0]+VVals[1]+VVals[2])/3.0;
    }
  fprintf(f,"};\n");
  fclose(f);

  if (ValuePerVertex) delete[] ValuePerVertex;
  if (AreaPerVertex) delete[] AreaPerVertex;
  if (NumPerVertex)   delete[] NumPerVertex;
  return Integral;
}

// alternate entry points to PlotScalarDensity
double RWGSurface::PlotScalarDensity(int PlotType, HMatrix *M, int WhichColumn,
                                     int AveragingMethod, const char *FileName,
                                     const char *ViewNameFormat, ...)
{ COMPLETE_VARARGS(ViewNameFormat, ViewName);
  void *Values = M->GetColumnPointer(WhichColumn);
  return PlotScalarDensity(PlotType, Values, AveragingMethod, FileName, ViewName);
}

double RWGSurface::PlotScalarDensity(double *Values, int AveragingMethod,
                                     const char *FileName, const char *ViewNameFormat, ...)
{ COMPLETE_VARARGS(ViewNameFormat, ViewName);
  return PlotScalarDensity(0, (void *)Values, AveragingMethod, FileName, ViewName);
}

double RWGSurface::PlotScalarDensity(int PlotType, cdouble *Values, int AveragingMethod,
                                     const char *FileName, const char *ViewNameFormat, ...)
{ COMPLETE_VARARGS(ViewNameFormat, ViewName);
  return PlotScalarDensity(PlotType, (void *)Values, AveragingMethod, FileName, ViewName);
}

void GetVectorValue(void *Values[3], int PlotType, int nv, double V[3])
{ for(int Mu=0; Mu<3; Mu++)
   V[Mu] = GetScalarValue(Values[Mu], PlotType, nv);
}

double RWGSurface::PlotVectorDensity(int PlotType, void *Values[3], const char *FileName,
                                     const char *ViewNameFormat, ...)
{ 
  COMPLETE_VARARGS(ViewNameFormat, ViewName);
  
  FILE *f=fopen(FileName,"a");
  if (!f) { Warn("could not open file %s",FileName); return 0.0; }

  double Integral=0.0;
  fprintf(f,"View \"%s\" {\n",ViewName);
  for(int np=0; np<NumPanels; np++)
   { double V[3];
     GetVectorValue(Values, PlotType, np, V);
     WriteVP(this, np, V, f);
     Integral += Panels[np]->Area * VecDot(V, Panels[np]->ZHat);
   }
  fprintf(f,"};\n");
  fclose(f);
  return Integral;
}

double RWGSurface::PlotVectorDensity(double *Values[3], const char *FileName,
                                     const char *ViewNameFormat, ...)
{
  COMPLETE_VARARGS(ViewNameFormat, ViewName);

  void *vValues[3];
  vValues[0]=(void *)Values[0];
  vValues[1]=(void *)Values[1];
  vValues[2]=(void *)Values[2];
  return PlotVectorDensity(0, vValues, FileName, ViewName);
}

double RWGSurface::PlotVectorDensity(int PlotType, HMatrix *M, int WhichColumn,
                                     const char *FileName, const char *ViewNameFormat, ...)
{ 
  if (M->NR != NumPanels)
   { Warn("PlotVectorDensity: wrong-length data (%i!=%i)",M->NR,NumPanels);
     return 0.0;
   }
  if (WhichColumn+2 >= M->NC)
   { Warn("PlotVectorDensity: too few columns (%i,%i)",WhichColumn,M->NC);
     return 0.0;
   }
  void *Values[3];
  Values[0] = M->GetColumnPointer(WhichColumn+0);
  Values[1] = M->GetColumnPointer(WhichColumn+1);
  Values[2] = M->GetColumnPointer(WhichColumn+2);
  COMPLETE_VARARGS(ViewNameFormat, ViewName);
  return PlotVectorDensity(PlotType, Values, FileName, ViewName);
}

/***************************************************************/
/* return true if panel #np on surface #ns is a ``straddler''  */
/* panel. detection of straddler panels is simpleminded: we    */
/* assume that the unit-cell geometry lives in the upper-right */
/* quadrant of the xy plane, so any panel whose centroid does  */
/* not live in that region is a straddler.                     */
/***************************************************************/
bool IsStraddlerPanel(RWGGeometry *G, int ns, int np)
{
  if (G->LBasis==0) return false;
  double *Centroid = G->Surfaces[ns]->Panels[np]->Centroid;
  if ( Centroid[0] < 0.0 ) return true;
  if ( (G->LBasis->NC>=2) && Centroid[1] < 0.0 ) return true;
  return false;
}

/***************************************************************/
/* Emit GMSH postprocessing code for visualizing induced       */
/* charges and currents                                        */
/***************************************************************/
void RWGGeometry::PlotSurfaceCurrents(HMatrix *PSDMatrix,
                                      cdouble Omega, double *kBloch,
                                      const char *FileNameFormat, ...)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  COMPLETE_VARARGS(FileNameFormat, FileName);
  FILE *f=fopen(FileName,"a");
  if (!f) { Warn("could not open file %s",FileName); return; }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char Tag[MAXSTR];
  snprintf(Tag,MAXSTR,"{Omega=%s",z2s(Omega));
  if (LDim > 0 && kBloch==0 )
   ErrExit("%s:%i: missing kBloch for PBC geometry",__FILE__,__LINE__);
  if (LDim>=1) vstrncat(Tag,MAXSTR,",kx=_%g",kBloch[0]);
  if (LDim>=2) vstrncat(Tag,MAXSTR,",ky=_%g",kBloch[1]);
  vstrncat(Tag,MAXSTR,"}");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  const char *QNames[5]={"Electric charge", "Electric current", "Magnetic charge", "Magnetic current", "Poynting flux"};
  const char *ReImStr[2]={"re","im"};
  for(int nv=0; nv<5; nv++)
   for(int ReIm=0; ReIm<2; ReIm++)
    for(int ns=0; ns<NumSurfaces; ns++)
     { 
       RWGSurface *S = Surfaces[ns];
       int Offset    = PanelIndexOffset[ns];
       if (S->IsPEC && nv>=2) continue;

       if (ns==0) fprintf(f,"View \"%s(%s) %s\" {\n",ReImStr[ReIm],QNames[nv],Tag);

       for(int np=0; np<S->NumPanels; np++)
        { if (kBloch && IsStraddlerPanel(this,ns,np)) continue;

          if (nv==0 || nv==2 || nv==4) // electric or magnetic charge
           { cdouble Val=PSDMatrix->GetEntry(Offset+np,(nv==0) ? 4 : (nv==2) ? 8 : 12);
             WriteST(S, np, ReIm ? imag(Val) : real(Val), f);
           }
          else // // electric or magnetic current
           { cdouble Vx=PSDMatrix->GetEntry(Offset+np,(nv==1) ? 5 : 9);
             cdouble Vy=PSDMatrix->GetEntry(Offset+np,(nv==1) ? 6 : 10);
             cdouble Vz=PSDMatrix->GetEntry(Offset+np,(nv==1) ? 7 : 11);
             if (ReIm==0)
              WriteVP(S, np, real(Vx), real(Vy), real(Vz), f);
             else
              WriteVP(S, np, imag(Vx), imag(Vy), imag(Vz), f);
           }
        } // for(int np=0; np<S->NumPanels; np++)

       if (ns==(NumSurfaces-1)) fprintf(f,"};\n");
     } // 
  fclose(f);
}

void RWGGeometry::PlotSurfaceCurrents(HMatrix *PSDMatrix,
                                      cdouble Omega,
                                      const char *FileNameFormat,
                                      ...)
{ 
  COMPLETE_VARARGS(FileNameFormat, FileName);
  PlotSurfaceCurrents(PSDMatrix, Omega, 0, FileName);
}

void RWGGeometry::PlotSurfaceCurrents(HVector *KN, cdouble Omega,
                                      double *kBloch,
                                      const char *FileNameFormat, ...)
{ 
  COMPLETE_VARARGS(FileNameFormat, FileName);
  HMatrix *PSD=GetPanelSourceDensities(Omega, kBloch, KN, 0);
  PlotSurfaceCurrents(PSD, Omega, kBloch, FileName);
  delete PSD;
}

void RWGGeometry::PlotSurfaceCurrents(HVector *KN, cdouble Omega,
                                      const char *FileNameFormat, ...)
{ 
  COMPLETE_VARARGS(FileNameFormat, FileName);
  PlotSurfaceCurrents(KN, Omega, FileName);
}

/***************************************************************/
/* MDFunc is a caller-supplied routine that inputs a list of   */
/* coordinates for NX points (as a 3 x NX HMatrix) and returns */
/* an NX x ND HMatrix whose columns are the values of ND data  */
/* quantities at those points. (MDFunc may optionally return   */
/* information on the names and nature of the data columns.)   */
/*                                                             */
/* PPOptions is an optional string containing GMSH post-       */
/* processing commands to be written to the output file.       */
/*                                                             */
/* If UseCentroids==true, the function is evaluated at the     */
/* panel centroids; otherwise at the vertices.                 */
/*                                                             */
/* The return value is a newly allocated HVector of length ND  */
/* whose entries are the integrals of the data quantities      */
/* over the area of the mesh.                                  */
/*                                                             */
/* If the value of a data quantity at any mesh vertex exceeds  */
/* ABSURD in absolute value, that point is excluded and the    */
/* average over the triangle is determined by the values at    */
/* the remaining vertices.                                     */
/***************************************************************/
#define ABSURD 1.0e100
HVector *RWGSurface::MakeMeshPlot(MeshDataFunc MDFunc, void *UserData,
                                  const char *OutFileBase, const char *PPOptions,
                                  bool UseCentroids)
{
  /***************************************************************/
  /* try to open output file and write user's GMSH PP options    */
  /***************************************************************/
  char PPFileName[100];
  snprintf(PPFileName, 100, "%s.pp",(OutFileBase ? OutFileBase : GetFileBase(MeshFileName)));
  FILE *f=fopen(PPFileName,"w");
  if (!f)
   { Warn("could not create output file %s",PPFileName);
     return 0;
   }
  if (PPOptions) fprintf(f,"%s",PPOptions);
  fclose(f);
     
  /***************************************************************/
  /* read in mesh and put vertex coordinates into XMatrix        */
  /***************************************************************/
  int NX = UseCentroids ? NumPanels : NumVertices;
  HMatrix XMatrix(3, NX, LHM_REAL);
  for(int nx=0; nx<NX; nx++)
   XMatrix.SetEntriesD(":", nx, UseCentroids ? Panels[nx]->Centroid : Vertices + 3*nx);

  /***************************************************************/
  /* call user's function to populate data matrix ****************/
  /***************************************************************/
  const char **DataNames=0;
  HMatrix *DataMatrix=MDFunc(UserData, &XMatrix, &DataNames);
  if ( DataMatrix==0 || DataMatrix->NR!=NX )
   { Warn("Mesh data function returned wrong-size DataMatrix");
     delete DataMatrix;
     return 0;
   }
  int NumData = DataMatrix->NC;
  bool RealData = (DataMatrix->RealComplex==LHM_REAL);
  HVector *Integrals = new HVector(NumData, LHM_REAL);

  /***************************************************************/
  /* By default we assume all data quantities are real-valued    */
  /* scalars.                                                    */
  /* The user's data function can change this by prepending a    */
  /* 4-character prefix to the DataName for a quantity:          */
  /* #RS_ real part of scalar quantity (default)                 */
  /* #IS_ imag part of scalar quantity                           */
  /* #MS_ magnitude of scalar quantity                           */
  /* #RV_ real part of vector quantity                           */
  /* #IV_ imag part of vector quantity                           */
  /* #MV_ magnitude of vector quantity                           */
  /***************************************************************/
  for(int nd=0; nd<NumData; nd++)
   {
     int PlotType = RealData ? 0 : 2;
     bool VectorPlot = false;
     char ndBuffer[10];
     const char *DataName = (DataNames ? DataNames[nd] : 0);
     if (DataName==0)
      { DataName=ndBuffer;
        snprintf(ndBuffer,10,"%i",nd);
      }
     else if (strlen(DataName)>4 && DataName[0]=='#' && DataName[3]=='_')
      { if (toupper(DataName[1])=='M') PlotType=1;
        if (toupper(DataName[1])=='I') PlotType=3;
        if (toupper(DataName[2])=='V') VectorPlot=true;
        DataName+=4;
        if (RealData && PlotType!=0)
         { Warn("Data name %s suggests complex data, but DataMatrix is real-valued",DataNames[nd]);
           PlotType=0;
         }
        if(VectorPlot && (nd+3 > NumData))
         { Warn("Not enough columns for vector data %s (skipping)",DataNames[nd]);
           continue;
         }
        if(VectorPlot && UseCentroids==false)
         { Warn("Vector-valued mesh plots require UseCentroids==true");
           continue;
         }
      }

     if (VectorPlot)
      { double I = PlotVectorDensity(PlotType, DataMatrix, nd, PPFileName, DataName);
        Integrals->SetEntry(nd,I);
        nd+=2;
      }
     else
      { int AveragingMethod = (UseCentroids ? 1 : 2);
        double I = PlotScalarDensity(PlotType, DataMatrix, nd, AveragingMethod, PPFileName, DataName);
        Integrals->SetEntry(nd, I);
      }
   }
  delete DataMatrix;
  return Integrals;
}

/***************************************************************/
/* (nd,nt) entry of returned matrix                            */
/*  = integral of quantity #nd for transform #nt               */
/***************************************************************/
HMatrix *MakeMeshPlot(MeshDataFunc MDFunc, void *UserData,
                      const char *MeshFile, const char *TransFile,
                      const char *CallerFileBase, const char *PPOptions, bool UseCentroids)
{
  char MeshFileCopy[100];
  snprintf(MeshFileCopy,100,"%s",MeshFile);
  
  char PPFileBase[100];
  if (CallerFileBase)
   snprintf(PPFileBase,100,"%s.%s",CallerFileBase,GetFileBase(MeshFileCopy));  
  else
   snprintf(PPFileBase,100,"%s",GetFileBase(MeshFileCopy));
  char FFBuffer[100], *FullFileBase = (TransFile ? FFBuffer : PPFileBase);

  RWGSurface *S=new RWGSurface(MeshFile);
  GTCList GTCs = ReadTransFile(TransFile);
  int NT = GTCs.size();
  HMatrix *Integrals=0;
  for(int nt=0; nt<NT; nt++)
   { if ( GTCs[nt]->GTs.size() > 0 )
      S->Transform(GTCs[nt]->GTs[0]);
     if (TransFile)
      snprintf(FullFileBase,100,"%s.%s",PPFileBase,GTCs[nt]->Tag);
     HVector *I=S->MakeMeshPlot(MDFunc, UserData, FullFileBase, PPOptions, UseCentroids);
     int NumData=I->N;
     if (!Integrals)
      Integrals = new HMatrix(NumData, NT, LHM_REAL);
     for(int nd=0; nd<NumData; nd++)
      Integrals->SetEntry(nd,nt,I->GetEntry(nd));
     delete I;
     S->UnTransform();
   }
  delete S;
  return Integrals;
}

} // namespace scuff
