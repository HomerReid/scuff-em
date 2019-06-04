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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "libscuff.h"
#include "EquivalentEdgePairs.h"

using namespace scuff;

#define MAXSTR 1000

/***************************************************************/
/* subroutine to analyze a single RWG surface ******************/
/***************************************************************/
void AnalyzeSurface(RWGSurface *S)
{
  double TotalArea=0.0, AvgArea=0.0, AvgArea2=0.0;
  double TotalVolume=0.0;
  for(int np=0; np<S->NumPanels; np++)
   { 
     RWGPanel *P=S->Panels[np];
     double Area = P->Area;

     TotalArea+=Area;
     AvgArea+=Area;
     AvgArea2+=Area*Area;

     double XmX0[3];
     XmX0[0] = P->Centroid[0];
     XmX0[1] = P->Centroid[1];
     XmX0[2] = P->Centroid[2];
     if (S->OTGT) S->OTGT->UnApply(XmX0);
     TotalVolume += VecDot(XmX0, P->ZHat) * Area / 3.0;
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
  if (S->IsClosed)
   printf(" Total volume: %9.7e \n",TotalVolume);
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
  
  if (G->Substrate)
   { printf("***********************************************\n");
     printf("*  SUBSTRATE:\n");
     printf("***********************************************\n");
     G->Substrate->Describe();
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
double GetRegionVolume(RWGGeometry *G, int nr)
{
  double V=0.0;
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     RWGSurface *S=G->Surfaces[ns];

     // open surfaces embedded in a region do
     // not contribute to the volume of that region
//     if (!S->IsClosed) 
//     continue;

     // surfaces that do not bound a region do not 
     // contribute to its volume
     double Sign=0.0;
     if (S->RegionIndices[0]==nr)
      Sign=-1.0;
     else if (S->RegionIndices[1]==nr)
      Sign=+1.0;
     else 
      continue;

     // evaluate surface integral of \vec{x} \dot \vec{nHat}/3
     // over this surface to get its contribution to the 
     // volume of the region it partially bounds
     for(int np=0; np<S->NumPanels; np++)
      { 
        RWGPanel *P=S->Panels[np];
        V += Sign*VecDot(P->Centroid, P->ZHat) * P->Area / 3.0;
      };
   };

  return V;
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

  int LDim = G->LBasis ? G->LBasis->NC : 0;
  int n2Mult = (LDim==1) ? 0 : 1;

  if (Neighbors)
   { double DL[3];
     char Tag[20];
     for(int n1=-Neighbors; n1<=Neighbors; n1++)
      for(int n2=-Neighbors*n2Mult; n2<=Neighbors*n2Mult; n2++)
       { 
         if (n1==0 && n2==0)
          continue;
 
         DL[0] = n1*G->LBasis->GetEntryD(0,0);
         DL[1] = n1*G->LBasis->GetEntryD(1,0);
         DL[2] = n1*G->LBasis->GetEntryD(2,0);

         if (LDim>=2)
          { DL[0] += n2*G->LBasis->GetEntryD(0,1);
            DL[1] += n2*G->LBasis->GetEntryD(1,1);
            DL[2] += n2*G->LBasis->GetEntryD(2,1);
          };

         if (LDim==1)
          snprintf(Tag,20,"(%i)",n1);
         else if (LDim==2)
          snprintf(Tag,20,"(%i,%i)",n1,n2);

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

  /***************************************************************/
  /* BF directions ***********************************************/
  /***************************************************************/
  fprintf(f,"View \"BFDirections\" {\n");
  for(ns=0; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], ne=0; ne<S->NumEdges; ne++)
    { 
      RWGEdge *E = S->Edges[ne];
      double *PCentroid, *MCentroid;
      PCentroid = S->Panels[E->iPPanel]->Centroid;
      if (E->iMPanel==-1)
        MCentroid = E->Centroid;
      else
        MCentroid = S->Panels[E->iMPanel]->Centroid;

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
 
  fclose(f);
  
}

/***************************************************************/
/* main function   *********************************************/
/***************************************************************/  
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* convenience shortcuts that allow the code to be invoked as  */
  /*  % scuff-analyze File.msh                                   */
  /* or                                                          */
  /*  % scuff-analyze File.scuffgeo                              */
  /***************************************************************/
  char *GeoFile=0;
  char *MeshFile=0;
  if (argc>=2)
   { char *Ext = GetFileExtension(argv[1]);
     if (Ext && !strcasecmp(Ext,"scuffgeo"))
      { GeoFile=strdup(argv[1]); 
        argv[1]=0; 
      }
     else if (Ext && !strcasecmp(Ext,"msh"))
      { MeshFile=strdup(argv[1]); 
        argv[1]=0; 
      }
   }

  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  int PhysicalRegion=-1;
  char *TransFile=0;
  bool WriteGPFiles=0;
  bool WritePPFiles=0;
  bool WriteLabels=0;
  int Neighbors=0;
  int EEPs[2]={-1,-1};
  bool RegionVolumes=false;
  double WhichRegion[3]={HUGE_VAL, HUGE_VAL, HUGE_VAL};
  char *EPFile=0;
  /* name, type, # args, max # instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"geometry",           PA_STRING, 1, 1, (void *)&GeoFile,        0, "geometry file"},
     {"mesh",               PA_STRING, 1, 1, (void *)&MeshFile,       0, "mesh file"},
     {"meshfile",           PA_STRING, 1, 1, (void *)&MeshFile,       0, "mesh file"},
     {"PhysicalRegion",     PA_INT,    1, 1, (void *)&PhysicalRegion, 0, "index of surface within mesh file"},
     {"TransFile",          PA_STRING, 1, 1, (void *)&TransFile,      0, "list of transformations"},
     {"EPFile",             PA_STRING, 1, 1, (void *)&EPFile,         0, "list of points"},
     {"WriteGnuplotFiles",  PA_BOOL,   0, 1, (void *)&WriteGPFiles,   0, "write gnuplot visualization files"},
     {"WriteGMSHFiles",     PA_BOOL,   0, 1, (void *)&WritePPFiles,   0, "write GMSH visualization files "},
     {"WriteGMSHLabels",    PA_BOOL,   0, 1, (void *)&WriteLabels,    0, "write GMSH labels"},
     {"Neighbors",          PA_INT,    1, 1, (void *)&Neighbors,      0, "number of neighboring cells to plot"},
     {"EEPs",               PA_INT,    2, 1, (void *)EEPs,            0, "analyze equivalent edge pairs for surfaces xx, xx"},
     {"RegionVolumes",      PA_BOOL,   0, 1, (void *)&RegionVolumes,  0, "compute volumes of closed regions"},
     {"WhichRegion",        PA_DOUBLE, 3, 1, (void *)&WhichRegion,    0, "identify region in which the given point lives"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (GeoFile==0 && MeshFile==0)
   OSUsage(argv[0], VERSION, OSArray,"either --geometry or --meshfile option must be specified");
  if (GeoFile!=0 && MeshFile!=0)
   ErrExit("--geometry and --meshfile options are mutually exclusive");
  if (PhysicalRegion!=-1 && MeshFile==0)
   ErrExit("--PhysicalRegion option may only be used with --meshfile");
  if (TransFile && GeoFile==0)
   ErrExit("--transfile option may only be used with --geometry");
  if (EPFile)
   WritePPFiles=true;
  if (WriteLabels)
   WritePPFiles=true;
   
  /***************************************************************/
  /**************************************************************/
  /***************************************************************/
  RWGSurface *S=0;
  RWGGeometry *G=0;
  SetLogFileName("scuff-analyze.log");
  if (MeshFile)
   {
     RWGGeometry::UseHRWGFunctions=false;
     S=new RWGSurface(MeshFile, PhysicalRegion);
     if (S->ErrMsg)
      ErrExit(S->ErrMsg);
     AnalyzeSurface(S);

     if (WriteGPFiles)
      WriteMeshGPFiles(S);

     if (WritePPFiles)
      WriteMeshPPFiles(S);

   }
  else
   { G=new RWGGeometry(GeoFile);

     AnalyzeGeometry(G);

     if (Neighbors!=0 && G->LBasis==0)
      { Warn("--Neighbors option only makes sense for periodic geometries");
        Neighbors=0;
      }

     if (WriteGPFiles)
      WriteGeometryGPFiles(G);

     if (WritePPFiles)
      WriteGeometryPPFiles(G, Neighbors);

     if (WriteLabels)
      WriteGeometryLabels(G);
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (TransFile)
   { 
     if (!G) 
      ErrExit("--transfile option may only be used with --geometry");

     GTCList GTCs=ReadTransFile(TransFile);
     int NGTC=GTCs.size();
     char *ErrMsg=G->CheckGTCList(GTCs);
     if (ErrMsg)
      ErrExit("file %s: %s",TransFile,ErrMsg);

     char PPFileName[MAXSTR];
     snprintf(PPFileName,MAXSTR,"%s.transformed.pp",GetFileBase(G->GeoFileName));
     unlink(PPFileName);

     for(int ngtc=0; ngtc<NGTC; ngtc++) 
      {
        G->Transform(GTCs[ngtc]);
        G->WritePPMesh(PPFileName, GTCs[ngtc]->Tag);
        G->UnTransform();
      }

     printf("Visualizations for %i transforms written to %s.\n",NGTC,PPFileName);
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (EPFile)
   { HMatrix *XMatrix=new HMatrix(EPFile);
     char PPFileName[MAXSTR];
     if (MeshFile)
      snprintf(PPFileName,MAXSTR,"%s.pp",GetFileBase(S->MeshFileName));
     else
      snprintf(PPFileName,MAXSTR,"%s.pp",GetFileBase(G->GeoFileName));

     FILE *f=fopen(PPFileName,"a");
     fprintf(f,"View \"%s\" {\n",EPFile);
     for(int n=0; n<XMatrix->NR; n++)
      fprintf(f,"SP(%e,%e,%e) {0.0};\n",
                 XMatrix->GetEntryD(n,0),
                 XMatrix->GetEntryD(n,1),
                 XMatrix->GetEntryD(n,2));
     fprintf(f,"};\n\n");
     fclose(f);
  }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (RegionVolumes)
   { printf("Region volumes: \n");
     for(int nr=1; nr<G->NumRegions; nr++)
      { 
        printf("%i (%-20s) ",nr,G->RegionLabels[nr]);
        double RV=GetRegionVolume(G,nr);
        if (RV==0.0)
         printf("unbounded\n");
        else
         printf("%e \n",RV);
      };
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (WhichRegion[0]!=HUGE_VAL)
   { int nr=G->GetRegionIndex(WhichRegion);
     printf("Point {%g,%g,%g} lives in region %s.\n",
             WhichRegion[0],WhichRegion[1],WhichRegion[2],G->RegionLabels[nr]);
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (EEPs[0]!=-1 && EEPs[0]!=-1)
   { char EEPFileName[100];
     if (MeshFile)
      { char GeoFileName[100];
        snprintf(GeoFileName,100,"MESH__%s",GetFileBase(MeshFile));
        G=new RWGGeometry(GeoFileName);
        if (!G) ErrExit("Could not open %s",GeoFileName);
        snprintf(EEPFileName,100,"%s.%i.%i.EEPs",GetFileBase(MeshFile),EEPs[0],EEPs[1]);
      }
     else
      snprintf(EEPFileName,100,"%s.%i.%i.EEPs",GetFileBase(GeoFile),EEPs[0],EEPs[1]);

     RWGSurface *Sa = G->Surfaces[EEPs[0]], *Sb = G->Surfaces[EEPs[1]]; 
     printf("Analyzing equivalent edge pairs between surfaces %i (%s) and %i (%s)...\n",
             EEPs[0],Sa->Label,EEPs[1],Sb->Label);
     EquivalentEdgePairTable EEPTable(G,EEPs[0],EEPs[1]);
     int NEA = Sa->NumEdges, NEB=Sb->NumEdges;
     int NEPairs = NEA * (Sa==Sb ? (NEA+1)/2 : NEB);
     int NumChildPairs  = EEPTable.NumChildren();
     int NumParentPairs = EEPTable.NumParents();
     printf(" Of %u total edge-edge pairs:\n ",NEPairs);
     printf("    %u are children (savings of %.1f %%)\n",NumChildPairs,100.0*((double)NumChildPairs)/((double)NEPairs));
     printf("    %u are parents (%.1f %%)\n",NumParentPairs, 100.0*((double)NumParentPairs) / ((double)NEPairs));
     printf("    %u are unicorns\n",NEPairs - NumParentPairs - NumChildPairs);

     EEPTable.Export(EEPFileName);
     printf("Wrote equivalent-edge pair table to %s.\n",EEPFileName);
   }

  if (MeshFile)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Thank you for your support.\n");

}
