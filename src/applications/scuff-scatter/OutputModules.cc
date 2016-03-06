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
 * OutputModules.cc -- various types of 'output modules' for EMScatter
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libSGJC.h>

#include "scuff-scatter.h"

#define MAXSTR   1000 

#define ABSTOL   0.0
#define RELTOL   5.0e-2
#define MAXEVALS 20000

namespace scuff {
void AddIFContributionsToEMTPFT(RWGGeometry *G, HVector *KN,
                                IncField *IF, cdouble Omega,
                                HMatrix *PFTMatrix, bool Interior);
                }

#define II cdouble(0.0,1.0)

/***************************************************************/
/* compute scattered and total fields at a user-specified list */
/* of evaluation points                                        */
/***************************************************************/
void ProcessEPFile(SSData *SSD, char *EPFileName)
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G  = SSD->G;
  IncField *IF    = SSD->IF;
  HVector  *KN    = SSD->KN;
  cdouble  Omega  = SSD->Omega;
  double *kBloch  = SSD->kBloch;

  /*--------------------------------------------------------------*/
  /*- try to read eval points from file --------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *XMatrix=new HMatrix(EPFileName,LHM_TEXT,"-ncol 3");
  if (XMatrix->ErrMsg)
   { fprintf(stderr,"Error processing EP file: %s\n",XMatrix->ErrMsg);
     delete XMatrix;
     return;
   };

  /*--------------------------------------------------------------*/
  /*- get components of scattered and incident fields            -*/
  /*--------------------------------------------------------------*/
  Log("Evaluating fields at points in file %s...",EPFileName);

  HMatrix *SFMatrix = G->GetFields( 0, KN, Omega, kBloch, XMatrix); // scattered
  HMatrix *IFMatrix = G->GetFields(IF,  0, Omega, kBloch, XMatrix); // incident

  /*--------------------------------------------------------------*/
  /*- create .scattered and .total output files and write fields -*/
  /*--------------------------------------------------------------*/
  SetDefaultCD2SFormat("%+.8e %+.8e ");
  const char *Ext[2]={"scattered","total"};
  for(int ST=0; ST<2; ST++)
   { char OutFileName[MAXSTR];
     snprintf(OutFileName,MAXSTR,"%s.%s",GetFileBase(EPFileName),Ext[ST]);
     FILE *f=CreateUniqueFile(OutFileName,1);
     fprintf(f,"# scuff-scatter run on %s (%s)\n",GetHostName(),GetTimeString());
     fprintf(f,"# columns: \n");
     fprintf(f,"# 1,2,3   x,y,z (evaluation point coordinates)\n");
     fprintf(f,"# 4,5     real, imag Ex\n");
     fprintf(f,"# 6,7     real, imag Ey\n");
     fprintf(f,"# 8,9     real, imag Ez\n");
     fprintf(f,"# 10,11   real, imag Hx\n");
     fprintf(f,"# 12,13   real, imag Hy\n");
     fprintf(f,"# 14,15   real, imag Hz\n");
     for(int nr=0; nr<SFMatrix->NR; nr++)
      { double X[3];
        cdouble EH[6];
        XMatrix->GetEntriesD(nr,":",X);
        SFMatrix->GetEntries(nr,":",EH);
        if (ST==1) 
         for(int nc=0; nc<6; nc++) 
          EH[nc]+=IFMatrix->GetEntry(nr,nc);
        fprintf(f,"%+.8e %+.8e %+.8e ",X[0],X[1],X[2]);
        fprintf(f,"%s %s %s   ",CD2S(EH[0]),CD2S(EH[1]),CD2S(EH[2]));
        fprintf(f,"%s %s %s\n", CD2S(EH[3]),CD2S(EH[4]),CD2S(EH[5]));
      };
     fclose(f);
   };

  delete XMatrix;
  delete SFMatrix;
  delete IFMatrix;

}

/***************************************************************/
/* VisualizeFields() produces a color plot of the E and H      */
/* fields on a user-specified surface mesh for visualization   */
/* in GMSH.                                                    */
/***************************************************************/
static const char *FieldTitles[]=
 {"|Ex|", "|Ey|", "|Ez|", "|E|",
  "|Hx|", "|Hy|", "|Hz|", "|H|",
 };

#define NUMFIELDFUNCS 8

void VisualizeFields(SSData *SSD, char *MeshFileName)
{ 
  /*--------------------------------------------------------------*/
  /*- try to open output file ------------------------------------*/
  /*--------------------------------------------------------------*/
  char GeoFileBase[100], PPFileName[100];
  strncpy(GeoFileBase,GetFileBase(SSD->G->GeoFileName),100);
  snprintf(PPFileName,100,"%s.%s.pp",GeoFileBase,GetFileBase(MeshFileName));
  FILE *f=fopen(PPFileName,"a");
  if (!f) 
   ErrExit("could not open field visualization file %s",PPFileName);
  
  /*--------------------------------------------------------------*/
  /*- try to open user's mesh file -------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=new RWGSurface(MeshFileName);

  Log("Creating flux plot for surface %s...",MeshFileName);
  printf("Creating flux plot for surface %s...\n",MeshFileName);

  /*--------------------------------------------------------------*/
  /*- create an Nx3 HMatrix whose columns are the coordinates of  */
  /*- the flux mesh panel vertices                                */
  /*--------------------------------------------------------------*/
  HMatrix *XMatrix=new HMatrix(S->NumVertices, 3);
  for(int nv=0; nv<S->NumVertices; nv++)
   XMatrix->SetEntriesD(nv, ":", S->Vertices + 3*nv);

  /*--------------------------------------------------------------*/
  /* 20150404 explain me -----------------------------------------*/
  /*--------------------------------------------------------------*/
  int nvRef=S->Panels[0]->VI[0];
  for(int nv=0; nv<S->NumVertices; nv++)
   { 
     bool VertexUsed=false;
     for(int np=0; np<S->NumPanels && !VertexUsed; np++)
      if (     nv==S->Panels[np]->VI[0]
           ||  nv==S->Panels[np]->VI[1]
           ||  nv==S->Panels[np]->VI[2]
         ) VertexUsed=true;

     if (!VertexUsed)
      { 
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
        printf("Replacing %i:{%.2e,%.2e,%.2e} with %i: {%.2e,%.2e,%.2e}\n",
                nv,XMatrix->GetEntryD(nv,0), 
                   XMatrix->GetEntryD(nv,1),
                   XMatrix->GetEntryD(nv,2),
                nvRef,S->Vertices[3*nvRef+0],
                      S->Vertices[3*nvRef+1],
                      S->Vertices[3*nvRef+2]);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

        XMatrix->SetEntriesD(nv, ":", S->Vertices + 3*nvRef);
      };
   };

  /*--------------------------------------------------------------*/
  /*- get the total fields at the panel vertices                 -*/
  /*--------------------------------------------------------------*/
  HMatrix *FMatrix=SSD->G->GetFields(SSD->IF, SSD->KN, SSD->Omega, SSD->kBloch, 
                                     XMatrix);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nff=0; nff<NUMFIELDFUNCS; nff++)
   { 
     fprintf(f,"View \"%s(%s)\" {\n",FieldTitles[nff],z2s(SSD->Omega));

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     for(int np=0; np<S->NumPanels; np++)
      { 
        RWGPanel *P=S->Panels[np];

        double *V[3]; // vertices
        double Q[3];  // quantities
        for(int nv=0; nv<3; nv++)
         { 
           int VI = P->VI[nv];
           V[nv]  = S->Vertices + 3*VI;

           cdouble EH[6];
           FMatrix->GetEntries(VI, ":", EH);
           cdouble *F = (nff>=4) ? EH+0 : EH+3;
           if (nff==3 || nff==7)
            Q[nv] = sqrt(norm(F[0]) + norm(F[1]) + norm(F[2]));
           else
            Q[nv] = abs( F[nff%4] );
         };
           
        fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                   V[0][0], V[0][1], V[0][2],
                   V[1][0], V[1][1], V[1][2],
                   V[2][0], V[2][1], V[2][2],
                   Q[0], Q[1], Q[2]);
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     fprintf(f,"};\n\n");
   };
  fclose(f);

  delete FMatrix;
  delete XMatrix;

  delete S;

}

/***************************************************************/
/* evaluate the induced dipole moments                         */
/***************************************************************/
void GetMoments(SSData *SSD, char *MomentFile)
{
  /***************************************************************/
  /* write file preamble on initial file creation ****************/
  /***************************************************************/
  FILE *f=fopen(MomentFile,"r");
  if (!f)
   { f=fopen(MomentFile,"w");
     if ( !f ) ErrExit("could not open file %s",MomentFile);
     fprintf(f,"# data file columns: \n");
     fprintf(f,"# 01    angular frequency (3e14 rad/sec)\n");
     fprintf(f,"#\n");
     fprintf(f,"# 02    surface label 1 \n");
     fprintf(f,"# 03,04 real,imag px_1 (electric dipole moment, surface 1)\n");
     fprintf(f,"# 05,06 real,imag py_1\n");
     fprintf(f,"# 07,08 real,imag pz_1\n");
     fprintf(f,"# 09,10 real,imag mx_1 (magnetic dipole moment, surface 1)\n");
     fprintf(f,"# 11,12 real,imag my_1\n");
     fprintf(f,"# 13,14 real,imag mz_1\n");
     fprintf(f,"#\n");
     fprintf(f,"# 15    surface label 2 \n");
     fprintf(f,"# 16-17 real,imag px_2 (electric dipole moment, surface 2)\n");
     fprintf(f,"# ...   \n");
     fprintf(f,"# 26,27 real,imag mz_2 (magnetic dipole moment, surface 2\n");
     fprintf(f,"# 28    surface label 3 \n");
     fprintf(f,"# ...   and so on\n");
   };
  fclose(f);

  /***************************************************************/
  /* open the file and write the frequency at the top of the line*/
  /***************************************************************/
  f=fopen(MomentFile,"a");
  
  /***************************************************************/
  /* get dipole moments ******************************************/
  /***************************************************************/
  RWGGeometry *G = SSD->G;
  HVector *KN    = SSD->KN;
  cdouble Omega  = SSD->Omega;

  HVector *PM  = G->GetDipoleMoments(Omega, KN);

  /***************************************************************/
  /* print to file ***********************************************/
  /***************************************************************/
  fprintf(f,"%s ",z2s(Omega));
  for (int ns=0; ns<G->NumSurfaces; ns++)
   { fprintf(f,"%s ",G->Surfaces[ns]->Label);
     for(int Mu=0; Mu<6; Mu++)
      fprintf(f,"%s ",CD2S(PM->GetEntry(6*ns + Mu),"%.8e %.8e "));
   };
  fprintf(f,"\n");
  fflush(f);

  delete PM;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WritePFTFilePreamble(FILE *f)
{
  fprintf(f,"# scuff-scatter run on %s (%s)\n",GetHostName(),GetTimeString());
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1 omega           (rad/sec) \n");
  fprintf(f,"# 2 surface label \n");
  fprintf(f,"# 3 absorbed power  (watts)\n");
  fprintf(f,"# 4 scattered power (watts)\n");
  fprintf(f,"# 5 x-force         (nanonewtons)\n");
  fprintf(f,"# 6 y-force         (nanonewtons)\n");
  fprintf(f,"# 7 z-force         (nanonewtons)\n");
  fprintf(f,"# 8 x-torque        (nanonewtons * microns)\n");
  fprintf(f,"# 9 y-torque        (nanonewtons * microns)\n");
  fprintf(f,"#10 z-torque        (nanonewtons * microns)\n");
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WritePFTFile(SSData *SSD, PFTOptions *PFTOpts, int Method,
                  bool PlotFlux, char *FileName)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGGeometry *G     = SSD->G;
  HVector *KN        = SSD->KN;
  cdouble Omega      = SSD->Omega;
  PFTOpts->RHSVector = SSD->RHS;
  PFTOpts->IF        = SSD->IF;
  PFTOpts->kBloch    = SSD->kBloch;
  PFTOpts->PFTMethod = Method;

  char FileNameBuffer[200];
  if (PlotFlux)
   { snprintf(FileNameBuffer, 200, "%s.pp", FileName);
     PFTOpts->FluxFileName=FileNameBuffer;
   }
  else
   PFTOpts->FluxFileName=0;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *PFTMatrix=G->GetPFTMatrix(KN, Omega, PFTOpts);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double *RegionPFTs=0;
  if (PFTOpts->GetRegionPFTs)
   RegionPFTs=(double *)mallocEC( G->NumRegions*NUMPFT*sizeof(double) );

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  FILE *f=fopen(FileName,"a");
  if (!f) return;
  WritePFTFilePreamble(f);
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     double PFT[NUMPFT];
     PFTMatrix->GetEntriesD(ns, ":", PFT);

     fprintf(f,"%s %s ",z2s(Omega),G->Surfaces[ns]->Label);
     for(int nq=0; nq<NUMPFT; nq++)
      fprintf(f,"%e ",PFT[nq]);
     fprintf(f,"\n");

     if (RegionPFTs)
      { RWGSurface *S=G->Surfaces[ns];
        int nr1=S->RegionIndices[0], nr2=S->RegionIndices[1];
        if (nr1>=0) 
         VecPlusEquals(RegionPFTs + NUMPFT*nr1, +1.0, PFT, NUMPFT);
        if (nr2>=0)
         VecPlusEquals(RegionPFTs + NUMPFT*nr2, -1.0, PFT, NUMPFT);
      };
   };

  if (RegionPFTs)
   { for(int nr=0; nr<G->NumRegions; nr++)
      { fprintf(f,"%s %s ",z2s(Omega),G->RegionLabels[nr]);
        for(int nq=0; nq<NUMPFT; nq++)
         fprintf(f,"%e ",RegionPFTs[nr*NUMPFT+ nq]);
        fprintf(f,"\n");
      };
     free(RegionPFTs);
   };

  fclose(f);
  delete PFTMatrix;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteEMTPFTFiles(SSData *SSD, char *FileBase)
{ 
  RWGGeometry *G = SSD->G;
  HVector *KN    = SSD->KN;
  IncField *IF   = SSD->IF;
  cdouble Omega  = SSD->Omega;

  PFTOptions MyOptions, *Options=&MyOptions;
  InitPFTOptions(Options);
  Options->PFTMethod = SCUFF_PFT_EMT;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NS = SSD->G->NumSurfaces;
  HMatrix *ExtinctionPFT, *ScatteredPFT[SCUFF_EMTPFT_NUMMETHODS];
  ExtinctionPFT=new HMatrix(NS, NUMPFT);
  for(int nm=0; nm<SCUFF_EMTPFT_NUMMETHODS; nm++)
   ScatteredPFT[nm]=new HMatrix(NS, NUMPFT);
  for(int ExtInt=0; ExtInt<2; ExtInt++)
   {
      bool Interior = (ExtInt==1);
      const char *IEString = (Interior ? "Interior" : "Exterior");

      ExtinctionPFT->Zero();
      AddIFContributionsToEMTPFT(G, KN, IF, Omega, ExtinctionPFT, Interior);

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      for(int Method=0; Method<SCUFF_EMTPFT_NUMMETHODS; Method++)
       { 
         Log("Computing scattered EMTPFT (%s%i) at Omega=%s...",IEString,Method,z2s(SSD->Omega));
         Options->Interior=Interior;
         Options->EMTPFTMethod=Method;
         G->GetPFTMatrix(KN, Omega, Options, ScatteredPFT[Method]);
       };

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      FILE *f=vfopen("%s.%cEMTPFT","r",FileBase,IEString[0]);
      bool WritePreamble = (f==0);
      if (f) fclose(f);
      f=vfopen("%s.%cEMTPFT","a",FileBase,IEString[0]);
      if (WritePreamble) WritePFTFilePreamble(f);
      for(int ns=0; ns<NS; ns++)
       { fprintf(f,"%s %s ",z2s(Omega),G->Surfaces[ns]->Label);
         for(int nq=0; nq<NUMPFT; nq++)
          fprintf(f,"%e ",ExtinctionPFT->GetEntryD(ns,nq));
         for(int nm=0; nm<SCUFF_EMTPFT_NUMMETHODS; nm++)
          for(int nq=0; nq<NUMPFT; nq++)
           fprintf(f,"%e ",ScatteredPFT[nm]->GetEntryD(ns,nq));
         fprintf(f,"\n");
       };
      fclose(f);

   }; // for(int Interior=0; Interior<2; Interior++)

  delete ExtinctionPFT;
  for(int nm=0; nm<SCUFF_EMTPFT_NUMMETHODS; nm++)
   delete ScatteredPFT[nm];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WritePSDFile(SSData *SSD, char *PSDFile)
{
  /* write output file preamble if necessary */
  FILE *f=fopen(PSDFile,"r");
  fclose(f);
  if (f==0)
   { f=fopen(PSDFile,"w");
     fprintf(f,"# Data file columns: \n");
     fprintf(f,"# 1:      angular frequency\n");
     fprintf(f,"# 2 3 4:  x, y, z coordinates of panel centroid\n");
     fprintf(f,"# 5:      panel area\n");
     fprintf(f,"# 6,  7   real, imag Sigma_E (electric surface charge density)\n");
     fprintf(f,"# 8,  9   real, imag K_x (electric surface current density)\n");
     fprintf(f,"# 10, 11  real, imag K_y (electric surface current density)\n");
     fprintf(f,"# 12, 13  real, imag K_z (electric surface current density)\n");
     fprintf(f,"# 14, 15  real, imag Sigma_M (magnetic surface charge density)\n");
     fprintf(f,"# 16, 17  real, imag N_x (magnetic surface current density)\n");
     fprintf(f,"# 18, 19  real, imag N_y (magnetic surface current density)\n");
     fprintf(f,"# 20, 21  real, imag N_z (magnetic surface current density)\n");
     fprintf(f,"# 22      inward-directed normal Poynting flux\n");
     fclose(f);
   };

  f=fopen(PSDFile,"a");
  if (!f)
   { Warn("could not open file %s for writing (skipping PSD output)",PSDFile);
     return;
   };

  Log("Computing panel source densities at Omega=%s...",z2s(SSD->Omega));

  static HMatrix *PSDMatrix=0;
  if (PSDMatrix==0)
   PSDMatrix=SSD->G->GetPanelSourceDensities(SSD->Omega, SSD->KN, 0); 
  else
   SSD->G->GetPanelSourceDensities(SSD->Omega, SSD->KN, PSDMatrix); 

  for(int nr=0; nr<PSDMatrix->NR; nr++)
   { fprintf(f,"%s ",z2s(SSD->Omega));
     for(int nc=0; nc<4; nc++)
      fprintf(f,"%e ",real(PSDMatrix->GetEntry(nr,nc)));
     for(int nc=4; nc<=11; nc++)
      fprintf(f,"%e %e ",real(PSDMatrix->GetEntry(nr,nc)),
                         imag(PSDMatrix->GetEntry(nr,nc)));
     for(int nc=12; nc<=12; nc++)
      fprintf(f,"%e ",real(PSDMatrix->GetEntry(nr,nc)));
     fprintf(f,"\n");
   };
  fclose(f);

}
