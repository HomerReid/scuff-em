/*
 * Visualize.cc -- code to produce gmsh and matlab files for visualizing
 *              -- libCas2D geometries under translations and rotations
 *
 * homer reid   -- 11/2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include <libhrutil.h>

#include "libTDRT.h"

/***************************************************************/
/* WriteGPMesh routine. This writes out the panel set as a     */
/* collection of lines. Useful for visualizing in gnuplot.     */
/***************************************************************/
void TDRTGeometry::WriteGPMesh(const char *format, ...)
{ 
  va_list ap;
  char FileName[1000];
  FILE *f;
  TDRTObject *O;
  double *V;
  int i, l, no, ns;

  va_start(ap,format);
  vsnprintf(FileName,1000,format,ap);
  va_end(ap);
  
  l=strlen(FileName);
  if ( l>3 && !strcmp(FileName+l-3,".gp") )
   FileName[l-3]=0;
  
  /***************************************************************/
  /* plot mesh ***************************************************/
  /***************************************************************/
  f=vfopen("%s.gp","w",FileName);
  if (!f) { fprintf(stderr,"warning: could not open file %s.gp\n",RemoveExtension(FileName)); return; }
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(ns=0; ns<O->NumSegments; ns++)
     { 
       V=O->Vertices + 2*O->Segments[2*ns];
       fprintf(f,"%e %e\n",V[0],V[1]);
       V=O->Vertices + 2*O->Segments[2*ns+1];
       fprintf(f,"%e %e\n",V[0],V[1]);
       fprintf(f,"\n\n");
     };
  fclose(f);
}

/***************************************************************/
/* WriteGPMeshPlus routine: Like WriteGPMesh, but also writes  */
/* text labels for plotting in gnuplot.                        */
/***************************************************************/
void TDRTGeometry::WriteGPMeshPlus(const char *format, ...)
{ 
  va_list ap;
  char FileName[1000];
  char buffer[100];
  FILE *f;
  TDRTObject *O;
  double *V, *V1, *V2; 
  int l, no, ns, niv, LabelIndex;
  double Location[2], MidPoint[2], SegmentNormal[2];
  double Length;

  va_start(ap,format);
  vsnprintf(FileName,1000,format,ap);
  va_end(ap);
  
  l=strlen(FileName);
  if ( l>3 && !strcmp(FileName+l-3,".gp") )
   FileName[l-3]=0;

  /***************************************************************/
  /* plot mesh ***************************************************/
  /***************************************************************/
  f=vfopen("%s.gp","w",FileName);
  if (!f) return;
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(ns=0; ns<O->NumSegments; ns++)
     { 
       V=O->Vertices  + 2*O->Segments[2*ns];
       fprintf(f,"%e %e\n",V[0],V[1]);
       V=O->Vertices  + 2*O->Segments[2*ns+1];
       fprintf(f,"%e %e\n",V[0],V[1]);
       fprintf(f,"\n\n");
     };
  fclose(f);

  /***************************************************************/
  /* plot element labels *****************************************/
  /***************************************************************/
  f=vfopen("%s.labels","w",FileName);
  if (!f) return;
  LabelIndex=1;
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(niv=0; niv<O->NumIVs; niv++)
    { 
      /* draw a label 'N+' just above the first line segment   */
      /* within basis function N                               */
      V1=O->Vertices + 2*O->Neighbors[2*niv];
      V2=O->Vertices + 2*O->IVs[niv];
      MidPoint[0]=0.5*(V2[0]+V1[0]);
      MidPoint[1]=0.5*(V2[1]+V1[1]);
      SegmentNormal[0]=(V2[1]-V1[1]);
      SegmentNormal[1]=-(V2[0]-V1[0]);
      VecNormalize(SegmentNormal);
      Length=VecDistance(V2,V1);
      Location[0]=MidPoint[0]+0.25*Length*SegmentNormal[0];
      Location[1]=MidPoint[1]+0.25*Length*SegmentNormal[1];
      fprintf(f,"set label %i \"%i+\" at %e,%e center tc lt 3\n",LabelIndex++,
               niv,Location[0],Location[1]);

      /* draw a label 'N-' just below the second line segment   */
      /* within basis function N                               */
      V1=O->Vertices + 2*O->IVs[niv];
      V2=O->Vertices + 2*O->Neighbors[2*niv+1];
      MidPoint[0]=0.5*(V2[0]+V1[0]);
      MidPoint[1]=0.5*(V2[1]+V1[1]);
      SegmentNormal[0]=(V2[1]-V1[1]);
      SegmentNormal[1]=-(V2[0]-V1[0]);
      VecNormalize(SegmentNormal);
      Length=VecDistance(V2,V1);
      Location[0]=MidPoint[0]-0.25*Length*SegmentNormal[0];
      Location[1]=MidPoint[1]-0.25*Length*SegmentNormal[1];
      fprintf(f,"set label %i \"%i-\" at %e,%e center tc lt 4\n",LabelIndex++,
               niv,Location[0],Location[1]);

    };
  fclose(f);

}

/***************************************************************/
/* WritePPMesh routine. Writes a .pp file suitable for opening */
/* as a GMSH postprocessing file.                              */
/* Note calling convention and file handling are different     */
/* from WriteGPMesh.                                           */
/***************************************************************/
void TDRTGeometry::WritePPMesh(char *FileName, char *Tag)
{ 
  FILE *f;
  TDRTObject *O;
  double *VA, *VB, Val;
  int i, no, ns;

/*
  va_list ap;
  va_start(ap,format);
  vsnprintf(FileName,997,format,ap);
  va_end(ap);
*/

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  f=fopen(FileName,"a");
  if (!f) return;
  fprintf(f,"View.ShowElement=1;\n");
  fprintf(f,"View.ShowScale=0;\n");
  fprintf(f,"View \"%s\" {\n",Tag);
  
  /***************************************************************/
  /* plot all line segments on all objects. **********************/
  /***************************************************************/
  for(no=0, O=Objects[0]; no<NumObjects; O=Objects[++no])
   for(ns=0; ns<O->NumSegments; ns++)
     { 
       VA=O->Vertices + 2*O->Segments[2*ns];
       VB=O->Vertices + 2*O->Segments[2*ns+1];

       Val=(double)(no+1);

       fprintf(f,"SL(%e,%e,%e,%e,%e,%e) {%e,%e};\n",
                  VA[0], VA[1], 0.0, VB[0], VB[1], 0.0, Val, Val);
     };

   
  fprintf(f,"};\n");
  fclose(f);
}
