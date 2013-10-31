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
 * scuff-analyze.cc -- a simple standalone program within the scuff-EM
 *                  -- suite for printing statistics on a geometry or  
 *                  -- on a single mesh, and optionally for generating
 *                  -- visualization files that identify the internal 
 *                  -- numbering scheme that libscuff uses for vertices, 
 *                  -- panels, and edges
 *
 * homer reid       -- 6/2009 -- 2/2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>

#include "libscuff.h"

using namespace scuff;

#define MAXSTR 1000

/***************************************************************/
/* subroutine to analyze a single RWG surface ******************/
/***************************************************************/
void AnalyzeSurface(RWGSurface *S)
{
  double TotalArea=0.0, AvgArea=0.0, AvgArea2=0.0;
  for(int np=0; np<S->NumPanels; np++)
   { TotalArea+=S->Panels[np]->Area;
     AvgArea+=S->Panels[np]->Area;
     AvgArea2+=(S->Panels[np]->Area)*(S->Panels[np]->Area);
   };
  AvgArea/=((double)S->NumPanels);
  AvgArea2/=((double)S->NumPanels);
  double StdDevArea = AvgArea2 - AvgArea*AvgArea;
  
  printf(" Meshfile: %s \n",S->MeshFileName);
  printf(" %i panels\n",S->NumPanels);
  printf(" %i total edges\n",S->NumTotalEdges);
  printf(" %i total basis functions\n",S->NumBFs);
  printf(" %i interior edges\n",S->NumEdges);
  printf(" %i total vertices ",S->NumVertices);
  printf("(after eliminating %i redundant vertices)\n",
            S->NumRedundantVertices);
  printf(" %i interior vertices\n",S->NumInteriorVertices);
  printf(" %i boundary contours\n",S->NumBCs);

  printf("\n");
  printf(" interior vertices - interior edges + panels = euler characteristic\n");
  printf(" %17i - %14i + %6i = %i\n",
            S->NumInteriorVertices, S->NumEdges, S->NumPanels, 
            S->NumInteriorVertices-S->NumEdges+S->NumPanels);
  printf("\n");

  printf(" Total area: %9.7e \n",TotalArea);
  printf(" Avg area: %9.7e // sqrt(Avg Area)=%9.7e\n",AvgArea,sqrt(AvgArea));
  printf(" Standard deviation of panel areas: %9.7e\n",StdDevArea);
  printf(" \n");

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteMeshGPFiles(RWGSurface *S)
{ S->WriteGPMesh("%s.gp",GetFileBase(S->MeshFileName));
  printf("Mesh visualization data written to GNUPLOT file %s.gp.\n\n",S->MeshFileName);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteMeshPPFiles(RWGSurface *S)
{
  char PPFileName[MAXSTR];
  snprintf(PPFileName,MAXSTR,"%s.pp",GetFileBase(S->MeshFileName));
  unlink(PPFileName);
  S->WritePPMesh(PPFileName,GetFileBase(S->MeshFileName),1);
  S->WritePPMeshLabels(PPFileName,GetFileBase(S->MeshFileName));
  printf("Mesh visualization data written to GMSH file %s.\n\n",PPFileName);
}

/***************************************************************/
/* subroutine to analyze an RWGGeometry ************************/
/***************************************************************/
void AnalyzeGeometry(RWGGeometry *G)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double RAM =   ((double)G->TotalBFs) 
               * ((double)G->TotalBFs) 
               * ((double)sizeof(cdouble));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("***********************************************\n");
  printf("*  GEOMETRY %s \n",G->GeoFileName);
  printf("***********************************************\n");
  printf(" %6i surfaces\n",G->NumSurfaces);
  printf(" %6i total basis functions\n",G->TotalBFs);
  printf(" Size of BEM matrix: %.2f GB\n",RAM/1.0e9);
  printf("\n");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("***********************************************\n");
  printf("*  REGIONS: %i\n",G->NumRegions);
  printf("***********************************************\n");
  printf("\n");
  int RLL, MaxRLL = strlen(G->RegionLabels[0]); // 'RLL=region label length'
  for(int nr=1; nr<G->NumRegions; nr++)
   { RLL=strlen(G->RegionLabels[nr]);
     if (RLL > MaxRLL) MaxRLL=RLL;
   };
  if (MaxRLL>150) 
   ErrExit("%s:%i:internal error",__FILE__,__LINE__);
  char fs[200]; // 'format string'

  sprintf(fs,"index | %%-%is | Material properties\n",MaxRLL);
  printf(fs,"Label");

  int nh;
  for(nh=0; nh<(MaxRLL+30); nh++) fs[nh]='-';
  fs[nh]=0;
  printf("%s\n",fs);

  sprintf(fs,"%%5i | %%-%is | %%s\n",MaxRLL);
  for(int nr=0; nr<G->NumRegions; nr++)
   printf(fs,nr,G->RegionLabels[nr],G->RegionMPs[nr]->Name);
  printf("\n");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int ns=0; ns<G->NumSurfaces; ns++)
   {
     RWGSurface *S=G->Surfaces[ns];
     printf("***********************************************\n");
     printf("*  SURFACE %i: Label = %s\n",ns,S->Label);
     printf("***********************************************\n");
     if (G->Mate[ns]!=-1)
      printf(" (duplicate of surface %s)\n\n",G->Surfaces[G->Mate[ns]]->Label);
     else
      { 
        if (S->IsPEC)
         printf("\n PEC surface in region %i (%s) \n\n",
                 S->RegionIndices[0], G->RegionLabels[S->RegionIndices[0]]);
        else
         printf("\n Interface between regions %i (%s) and %i (%s)\n\n",
              S->RegionIndices[0], G->RegionLabels[S->RegionIndices[0]],
              S->RegionIndices[1], G->RegionLabels[S->RegionIndices[1]]);
        AnalyzeSurface(G->Surfaces[ns]);
      };
   };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteGeometryGPFiles(RWGGeometry *G)
{
  G->WriteGPMesh("%s.gp",GetFileBase(G->GeoFileName));
  printf("Geometry visualization data written to GNUPLOT file %s.gp.\n\n",
          GetFileBase(G->GeoFileName));
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteGeometryPPFiles(RWGGeometry *G, int Neighbors)
{
  char PPFileName[MAXSTR];
  snprintf(PPFileName,MAXSTR,"%s.pp",GetFileBase(G->GeoFileName));
  unlink(PPFileName);

  G->WritePPMesh(PPFileName,G->GeoFileName);

  int LDim = G->NumLatticeBasisVectors;
  int n2Mult = (LDim==1) ? 0 : 1;

  if (Neighbors)
   { double DL[3];
     char Tag[20];
     for(int n1=-Neighbors; n1<=Neighbors; n1++)
      for(int n2=-Neighbors*n2Mult; n2<=Neighbors*n2Mult; n2++)
       { 
         if (n1==0 && n2==0) 
          continue;
 
         if (LDim==1)
          { DL[0] = n1*G->LatticeBasisVectors[0][0];
            DL[1] = n1*G->LatticeBasisVectors[0][1];
            DL[2] = n1*G->LatticeBasisVectors[0][2];
            snprintf(Tag,20,"(%i)",n1);
          }
         else
          { DL[0] = n1*G->LatticeBasisVectors[0][0] + n2*G->LatticeBasisVectors[1][0];
            DL[1] = n1*G->LatticeBasisVectors[0][1] + n2*G->LatticeBasisVectors[1][1];
            DL[2] = n1*G->LatticeBasisVectors[0][2] + n2*G->LatticeBasisVectors[1][2];
            snprintf(Tag,20,"(%i,%i)",n1,n2);
          }

         for(int ns=0; ns<G->NumSurfaces; ns++)
          G->Surfaces[ns]->Transform("DISPLACED %e %e %e \n",DL[0],DL[1],DL[2]);

         G->WritePPMesh(PPFileName,Tag);

         for(int ns=0; ns<G->NumSurfaces; ns++)
          G->Surfaces[ns]->UnTransform();
      };
            
   };

  printf("Geometry visualization data written to GMSH file %s.\n\n",PPFileName);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteGeometryLabels(RWGGeometry *G)
{
  int ns, np, ne, npOverall, nbfOverall;
  double Val, *X, *Z;
  RWGSurface *S;

  FILE *f=vfopen("%s.pp","a",GetFileBase(G->GeoFileName));

  /***************************************************************/
  /* panel normals ***********************************************/
  /***************************************************************/
  fprintf(f,"View \"Panel normals\" {\n");
  for(ns=0; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], np=0; np<S->NumPanels; np++)
    { 
      //Val=fmax(VecNorm(P->Centroid), 2.0*P->Area);
      Val=S->Panels[np]->Radius;
      X=S->Panels[np]->Centroid;
      Z=S->Panels[np]->ZHat;
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",X[0],X[1],X[2],Val*Z[0],Val*Z[1],Val*Z[2]);
    };
  fprintf(f,"};\n");

  /***************************************************************/
  /* panel indices ***********************************************/
  /***************************************************************/
  fprintf(f,"View \"Panel indices\" {\n");
  for(npOverall=ns=0; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], np=0; np<S->NumPanels; np++, npOverall++)
    { X=S->Panels[np]->Centroid;
      fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",X[0],X[1],X[2],npOverall);
    };
  fprintf(f,"};\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=0;\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowScale=0;\n");

  /***************************************************************/
  /* BF indices **************************************************/
  /***************************************************************/
  fprintf(f,"View \"Basis function indices\" {\n");
  for(nbfOverall=ns=0; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], ne=0; ne<S->NumEdges; ne++)
    { X=S->Edges[ne]->Centroid;
      if (S->IsPEC)
       { fprintf(f,"T3(%e,%e,%e,0.0) {\"%i\"};\n",X[0],X[1],X[2],nbfOverall);
         nbfOverall+=1;
       }
      else
       { fprintf(f,"T3(%e,%e,%e,0.0) {\"%i,%i\"};\n",X[0],X[1],X[2],nbfOverall,nbfOverall+1);
         nbfOverall+=2;
       };
    };
  fprintf(f,"};\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowElement=0;\n");
  fprintf(f,"View[PostProcessing.NbViews-1].ShowScale=0;\n");
 
  fclose(f);
  
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
  int PhysicalRegion=-1;
  char *TransFile=0;
  int WriteGPFiles=0;
  int WritePPFiles=0;
  int WriteLabels=0;
  int Neighbors=0;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING, 1, 1, (void *)&GeoFile,        0, "geometry file"},
     {"mesh",               PA_STRING, 1, 1, (void *)&MeshFile,       0, "mesh file"},
     {"meshfile",           PA_STRING, 1, 1, (void *)&MeshFile,       0, "mesh file"},
     {"PhysicalRegion",     PA_INT,    1, 1, (void *)&PhysicalRegion, 0, "index of surface within mesh file"},
     {"transfile",          PA_STRING, 1, 1, (void *)&TransFile,      0, "list of transformations"},
     {"WriteGnuplotFiles",  PA_BOOL,   0, 1, (void *)&WriteGPFiles,   0, "write gnuplot visualization files"},
     {"WriteGMSHFiles",     PA_BOOL,   0, 1, (void *)&WritePPFiles,   0, "write GMSH visualization files "},
     {"WriteGMSHLabels",    PA_BOOL,   0, 1, (void *)&WriteLabels,    0, "write GMSH labels"},
     {"Neighbors",          PA_INT,    1, 1, (void *)&Neighbors,      0, "number of neighboring cells to plot"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (GeoFile==0 && MeshFile==0)
   OSUsage(argv[0],OSArray,"either --geometry or --meshfile option must be specified");
  if (GeoFile!=0 && MeshFile!=0)
   ErrExit("--geometry and --meshfile options are mutually exclusive");
  if (PhysicalRegion!=-1 && MeshFile==0)
   ErrExit("--PhysicalRegion option may only be used with --meshfile");
  if (TransFile && GeoFile==0)
   ErrExit("--transfile option may only be used with --geometry");
   
  /***************************************************************/
  /**************************************************************/
  /***************************************************************/
  RWGSurface *S=0;
  RWGGeometry *G=0;
  SetLogFileName("scuff-analyze.log");
  if (MeshFile)
   {
     RWGGeometry::AssignBasisFunctionsToExteriorEdges=false;
     S=new RWGSurface(MeshFile, PhysicalRegion);
     AnalyzeSurface(S);

     if (WriteGPFiles)
      WriteMeshGPFiles(S);

     if (WritePPFiles)
      WriteMeshPPFiles(S);

   }
  else
   { G=new RWGGeometry(GeoFile);

     AnalyzeGeometry(G);

     if (Neighbors!=0 && G->NumLatticeBasisVectors==0 )
      { Warn("--Neighbors option only makes sense for periodic geometries");
        Neighbors=0;
      }

     if (WriteGPFiles)
      WriteGeometryGPFiles(G);

     if (WritePPFiles)
      WriteGeometryPPFiles(G, Neighbors);

     if (WriteLabels)
      WriteGeometryLabels(G);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (TransFile)
   { 
     if (!G) 
      ErrExit("--transfile option may only be used with --geometry");

     int ngtc, NGTC;
     GTComplex **GTCList=ReadTransFile(TransFile, &NGTC);
     char *ErrMsg=G->CheckGTCList(GTCList, NGTC);
     if (ErrMsg)
      ErrExit("file %s: %s",TransFile,ErrMsg);

     char PPFileName[MAXSTR];
     snprintf(PPFileName,MAXSTR,"%s.transformed.pp",GetFileBase(G->GeoFileName));
     unlink(PPFileName);

     for(ngtc=0; ngtc<NGTC; ngtc++) 
      {
        G->Transform(GTCList[ngtc]);
        G->WritePPMesh(PPFileName, GTCList[ngtc]->Tag);
        G->UnTransform();
      };

     printf("Visualizations for %i transforms written to %s.\n",NGTC,PPFileName);
 
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
