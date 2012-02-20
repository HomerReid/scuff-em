/*
 * scuff-analyze.cc -- a simple standalone program within the scuff-EM
 *                  -- suite for printing statistics on a geometry or  
 *                  -- on a single object, and optionally for generating
 *                  -- visualization files that identify the internal 
 *                  -- numbering scheme that libscuff uses for vertices, 
 *                  -- panels, and edges
 *
 * homer reid       -- 6/2009 -- 2/2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libhrutil.h>

#include "libscuff.h"

#define MAXSTR 1000

/***************************************************************/
/* subroutine to analyze a single RWG object *******************/
/***************************************************************/
void AnalyzeObject(RWGObject *O, int WriteGPFiles, int WritePPFiles)
{
  double TotalArea, AvgArea;
  int np;
  for(TotalArea=AvgArea=0.0, np=0; np<O->NumPanels; np++)
   { TotalArea+=O->Panels[np]->Area;
     AvgArea+=O->Panels[np]->Area;
   };
  AvgArea/=((double)O->NumPanels);
  
  printf(" Meshfile: %s \n",O->MeshFileName);
  printf(" %i panels\n",O->NumPanels);
  printf(" %i total edges\n",O->NumTotalEdges);
  printf(" %i total basis functions\n",O->NumBFs);
  printf(" %i interior edges\n",O->NumEdges);
  printf(" %i total vertices ",O->NumVertices);
  printf("(after eliminating %i redundant vertices)\n",
            O->NumRedundantVertices);
  printf(" %i interior vertices\n",O->NumInteriorVertices);
  printf(" %i boundary contours\n",O->NumBCs);

  printf("\n");
  printf(" interior vertices - interior edges + panels = euler characteristic\n");
  printf(" %17i - %14i + %6i = %i\n",
            O->NumInteriorVertices, O->NumEdges, O->NumPanels, 
            O->NumInteriorVertices-O->NumEdges+O->NumPanels);
  printf("\n");

  printf(" Total area: %9.7e \n",TotalArea);
  printf(" Avg area: %9.7e // sqrt(Avg Area)=%9.7e\n",AvgArea,sqrt(AvgArea));
  printf(" \n");

  if (WriteGPFiles)
   { O->WriteGPMesh("%s.gp",GetFileBase(O->MeshFileName));
     printf("Mesh visualization data written to GNUPLOT file %s.gp.\n\n",O->MeshFileName);
   };

  if (WritePPFiles)
   { char buffer[MAXSTR];
     snprintf(buffer,MAXSTR,"%s.pp",GetFileBase(O->MeshFileName));
     unlink(buffer);
     O->WritePPMesh(buffer,GetFileBase(O->MeshFileName),1);
     O->WritePPMeshLabels(buffer,GetFileBase(O->MeshFileName));
     printf("Mesh visualization data and panel normals written to GMSH file %s.\n\n",buffer);
   };

}

/***************************************************************/
/* subroutine to analyze an RWGGeometry ************************/
/***************************************************************/
void AnalyzeGeometry(char *GeoFile, int WriteGPFiles, int WritePPFiles)
{ 
  RWGGeometry *G=new RWGGeometry(GeoFile);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  unsigned RAM = G->TotalBFs * G->TotalBFs * sizeof(cdouble);

  printf("***********************************************\n");
  printf("*  GEOMETRY %s \n",GeoFile);
  printf("***********************************************\n");
  printf(" %6i objects\n",G->NumObjects);
  printf(" %6i total basis functions\n",G->TotalBFs);
  printf(" Size of BEM matrix: %.2f GB\n",((double)RAM)/1.0e9);
  printf("\n");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int no;
  RWGObject *O;
  for(no=0, O=G->Objects[0]; no<G->NumObjects; O=G->Objects[++no])
   {
     printf("***********************************************\n");
     printf("*  OBJECT %i: Label = %s\n",no,O->Label);
     printf("***********************************************\n");
     if (G->Mate[no]!=-1)
      printf(" (duplicate of object %i)\n\n",G->Mate[no]+1);
     else
      AnalyzeObject(G->Objects[no],0,0);
   };

  if (WriteGPFiles)
   { G->WriteGPMesh("%s.gp",GetFileBase(GeoFile));
     printf("Geometry visualization data written to GNUPLOT file %s.gp.\n\n",
             GetFileBase(GeoFile));
   };

  if (WritePPFiles)
   { char buffer[MAXSTR];
     snprintf(buffer,MAXSTR,"%s.pp",GetFileBase(GeoFile));
     unlink(buffer);
     G->WritePPMesh(buffer,buffer);
     printf("Geometry visualization data written to GMSH file %s.\n\n",buffer);
   };

}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  char *MeshFile=0;
  int WriteGPFiles=0;
  int WritePPFiles=0;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING, 1, 1, (void *)&GeoFile,      0, "geometry file"},
     {"mesh",               PA_STRING, 1, 1, (void *)&MeshFile,     0, "mesh file"},
     {"WriteGnuplotFiles",  PA_BOOL,   0, 1, (void *)&WriteGPFiles, 0, "write gnuplot visualization files"},
     {"WriteGMSHFiles",     PA_BOOL,   0, 1, (void *)&WritePPFiles, 0, "write GMSH visualization files "},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (GeoFile==0 && MeshFile==0)
   OSUsage(argv[0],OSArray,"either --geometry or --mesh option must be specified");
  if (GeoFile!=0 && MeshFile!=0)
   OSUsage(argv[0],OSArray,"--geometry and --mesh options are mutually exclusive");
   
  /***************************************************************/
  /**************************************************************/
  /***************************************************************/
  if (MeshFile)
   { RWGObject *O=new RWGObject(MeshFile);
     AnalyzeObject(O, WriteGPFiles, WritePPFiles);
   }
  else
   AnalyzeGeometry(GeoFile, WriteGPFiles, WritePPFiles);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
