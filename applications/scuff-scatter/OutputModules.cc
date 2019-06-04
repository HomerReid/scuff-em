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
  char *FileBase  = SSD->FileBase;

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
  char OmegaStr[100];
  snprintf(OmegaStr,100,"%s",z2s(Omega));
  char *TransformLabel=SSD->TransformLabel;
  char *IFLabel=SSD->IFLabel;
  const char *Ext[2]={"scattered","total"};
  for(int ST=0; ST<2; ST++)
   { char OutFileName[MAXSTR];
     snprintf(OutFileName,MAXSTR,"%s.%s.%s",FileBase,GetFileBase(EPFileName),Ext[ST]);
     FILE *f=fopen(OutFileName,"a");
     fprintf(f,"# scuff-scatter " VERSIONSTRING " run on %s (%s)\n",GetHostName(),GetTimeString());
     fprintf(f,"# columns: \n");
     fprintf(f,"# 1,2,3   x,y,z (evaluation point coordinates)\n");
     fprintf(f,"# 4       omega (angular frequency)\n");
     int nc=5;
     if (TransformLabel)
      fprintf(f,"# %i       geometrical transform\n",nc++);
     if (IFLabel)
      fprintf(f,"# %i       incident field\n",nc++);
     fprintf(f,"# %02i,%02i   real, imag Ex\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Ey\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Ez\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Hx\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Hy\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i   real, imag Hz\n",nc,nc+1); nc+=2;
     for(int nr=0; nr<SFMatrix->NR; nr++)
      { double X[3];
        cdouble EH[6];
        XMatrix->GetEntriesD(nr,":",X);
        SFMatrix->GetEntries(nr,":",EH);
        if (ST==1) 
         for(int n=0; n<6; n++) 
          EH[n]+=IFMatrix->GetEntry(nr,n);
        fprintf(f,"%+.8e %+.8e %+.8e ",X[0],X[1],X[2]);
        fprintf(f,"%s ",OmegaStr);
        if (TransformLabel) fprintf(f,"%s ",TransformLabel);
        if (IFLabel) fprintf(f,"%s ",IFLabel);
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
/***************************************************************/
/***************************************************************/
static const char *FieldTitles[]=
 {"Ex", "Ey", "Ez", "|E|",
  "Hx", "Hy", "Hz", "|H|",
 };

#define NUMFIELDFUNCS 8

void WriteFVMesh(SSData *SSD, RWGSurface *S, FILE *f, char *FuncString)
{
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
      { printf("Replacing %i:{%.2e,%.2e,%.2e} with %i: {%.2e,%.2e,%.2e}\n",
                nv,XMatrix->GetEntryD(nv,0), XMatrix->GetEntryD(nv,1), XMatrix->GetEntryD(nv,2),
                nvRef,S->Vertices[3*nvRef+0], S->Vertices[3*nvRef+1], S->Vertices[3*nvRef+2]);
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
  char *TransformLabel=SSD->TransformLabel, *IFLabel=SSD->IFLabel;
  for(int nff=0; nff<NUMFIELDFUNCS; nff++)
   { 
     if (FuncString && !strcasestr(FuncString,FieldTitles[nff]))
      continue;
     fprintf(f,"View \"%s(%s)",FieldTitles[nff],z2s(SSD->Omega));
     if (TransformLabel)
      fprintf(f,"(%s)",TransformLabel);
     if (IFLabel)
      fprintf(f,"(%s)",IFLabel);
     fprintf(f,"\" {\n");

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
      }; // for(int np=0; np<S->NumPanels; np++)

     fprintf(f,"};\n\n");

   };  // for(int nff=0; nff<NUMFIELDFUNCS; nff++)

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  delete FMatrix;
  delete XMatrix;
}

/***************************************************************/
/* VisualizeFields() produces a color plot of the E and H      */
/* fields on a user-specified surface mesh for visualization   */
/* in GMSH.                                                    */
/***************************************************************/
void VisualizeFields(SSData *SSD, char *FVMesh, char *TransFile, char *FuncList)
{ 
  GTCList GTCs=ReadTransFile(TransFile);
  int NT = GTCs.size();
  
  /*--------------------------------------------------------------*/
  /*- try to open user's mesh file -------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=new RWGSurface(FVMesh);

  char *GeoFileBase=SSD->FileBase; 
  char *FVMFileBase=GetFileBase(FVMesh);
  for(int nt=0; nt<NT; nt++)
   {
     char *Tag = GTCs[nt]->Tag;
     char PPFileName[100];
     if (NT>1)
      { 
        snprintf(PPFileName,100,"%s.%s.%s.pp",GeoFileBase,FVMFileBase,Tag);
        Log("Creating flux plot for surface %s, transform %s...",FVMesh,Tag);
      }
     else
      {  snprintf(PPFileName,100,"%s.%s.pp",GeoFileBase,FVMFileBase);
         Log("Creating flux plot for surface %s...",FVMesh);
      }
     FILE *f=fopen(PPFileName,"a");
     if (!f) 
      { Warn("could not open field visualization file %s",PPFileName);
       continue;
      }

     if (GTCs[nt]->GTs.size()>0) S->Transform(GTCs[nt]->GTs[0]);
     WriteFVMesh(SSD, S, f, FuncList);
     S->UnTransform();

     fclose(f);
   }

  delete S;
  DestroyGTCList(GTCs);

}

/***************************************************************/
/* evaluate the induced dipole moments                         */
/***************************************************************/
void GetMoments(SSData *SSD, char *MomentFile)
{
  char *TransformLabel = SSD->TransformLabel;
  char *IFLabel        = SSD->IFLabel;

  /***************************************************************/
  /* write file preamble on initial file creation ****************/
  /***************************************************************/
  FILE *f=fopen(MomentFile,"r");
  if (!f)
   { f=fopen(MomentFile,"w");
     if ( !f ) ErrExit("could not open file %s",MomentFile);
     fprintf(f,"# data file columns: \n");
     fprintf(f,"# 1     angular frequency (3e14 rad/sec)\n");
     int nc=2;
     if (TransformLabel)
      fprintf(f,"# %i     geometrical transform\n",nc++);
     if (IFLabel)
      fprintf(f,"# %i     incident field\n",nc++);
     fprintf(f,"# %i     surface label\n",nc++);
     fprintf(f,"# %02i,%02i real,imag px (electric dipole moment)\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag py \n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag pz \n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag mx (magnetic dipole moment)\n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag my \n",nc,nc+1); nc+=2;
     fprintf(f,"# %02i,%02i real,imag mz \n",nc,nc+1); nc+=2;
     fprintf(f,"#\n");
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
  HMatrix *PM    = G->GetDipoleMoments(Omega, KN);

  /***************************************************************/
  /* print to file ***********************************************/
  /***************************************************************/
  for (int ns=0; ns<G->NumSurfaces; ns++)
   { fprintf(f,"%s ",z2s(Omega));
     if (TransformLabel)
      fprintf(f,"%s ",TransformLabel);
     if (IFLabel)
      fprintf(f,"%s ",IFLabel);
     fprintf(f,"%s ",G->Surfaces[ns]->Label);
     for(int Mu=0; Mu<6; Mu++)
      fprintf(f,"%s ",CD2S(PM->GetEntry(ns,Mu),"%.8e %.8e "));
     fprintf(f,"\n");
   };
  fclose(f);

  delete PM;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WritePFTFilePreamble(FILE *f, char *TransformLabel, char *IFLabel)
{
  fprintf(f,"# scuff-scatter " VERSIONSTRING " run on %s (%s)\n",GetHostName(),GetTimeString());
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1   omega           (rad/sec) \n");
  int nc=2;
  if (TransformLabel)
   fprintf(f,"# %i   geometrical transform\n",nc++);
  if (IFLabel)
   fprintf(f,"# %i   incident field\n",nc++);
  fprintf(f,"#%2i   surface label \n",nc++);             
  fprintf(f,"#%2i   absorbed power  (watts)\n",nc++);
  fprintf(f,"#%2i   scattered power (watts)\n",nc++);
  fprintf(f,"#%2i   x-force         (nanonewtons)\n",nc++);
  fprintf(f,"#%2i   y-force         (nanonewtons)\n",nc++);
  fprintf(f,"#%2i   z-force         (nanonewtons)\n",nc++);
  fprintf(f,"#%2i   x-torque        (nanonewtons * microns)\n",nc++);
  fprintf(f,"#%2i   y-torque        (nanonewtons * microns)\n",nc++);
  fprintf(f,"#%2i   z-torque        (nanonewtons * microns)\n",nc++);
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
  RWGGeometry *G       = SSD->G;
  HVector *KN          = SSD->KN;
  PFTOpts->RHSVector   = SSD->RHS;
  cdouble Omega        = SSD->Omega;
  PFTOpts->kBloch      = SSD->kBloch;
  char *TransformLabel = SSD->TransformLabel;
  char *IFLabel        = SSD->IFLabel;
  PFTOpts->IF          = SSD->IF;
  PFTOpts->PFTMethod   = Method;

  char FileNameBuffer[200];
  if (PlotFlux)
   { snprintf(FileNameBuffer, 200, "%s.pp", FileName);
     PFTOpts->FluxFileName=FileNameBuffer;
   }
  else
   PFTOpts->FluxFileName=0;

  static bool WrotePreamble[SCUFF_PFT_NUMMETHODS];
  static bool Initialized=false;
  if (!Initialized)
   { Initialized=true; 
     memset(WrotePreamble, 0, SCUFF_PFT_NUMMETHODS*sizeof(bool));
   };

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
  if (!WrotePreamble[Method]) 
   { WritePFTFilePreamble(f,TransformLabel,IFLabel);
     WrotePreamble[Method]=true;
   };
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     double PFT[NUMPFT];
     PFTMatrix->GetEntriesD(ns, ":", PFT);

     fprintf(f,"%s ",z2s(Omega));
     if (TransformLabel) 
      fprintf(f,"%s ",TransformLabel);
     if (IFLabel) 
      fprintf(f,"%s ",IFLabel);
     fprintf(f,"%s ",G->Surfaces[ns]->Label);
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
void WritePSDFile(SSData *SSD, char *PSDFile)
{
  /* write output file preamble if necessary */
  char *TransformLabel = SSD->TransformLabel;
  char *IFLabel        = SSD->IFLabel;
  FILE *f=fopen(PSDFile,"r");
  fclose(f);
  if (f==0)
   { f=fopen(PSDFile,"w");
     fprintf(f,"# Data file columns: \n");
     fprintf(f,"# 1:      angular frequency\n");
     int nc=2;
     if (TransformLabel)
      fprintf(f,"# %i:     geometrical transform\n",nc++);
     if (IFLabel)
      fprintf(f,"# %i:     incident field\n",nc++);
     fprintf(f,"# %i %i %i: x, y, z coordinates of panel centroid\n",nc, nc+1, nc+2); nc+=3;
     fprintf(f,"# %i      panel area\n",nc++);
     fprintf(f,"# %i, %i  real, imag Sigma_E (electric surface charge density)\n",nc+1,nc+2); nc+=2;
     fprintf(f,"# %i, %i  real, imag K_x (electric surface current density)\n",nc+1,nc+2); nc+=2;
     fprintf(f,"# %i, %i  real, imag K_y\n",nc+1,nc+2); nc+=2;
     fprintf(f,"# %i, %i  real, imag K_z\n",nc+1,nc+2); nc+=2;
     fprintf(f,"# %i, %i  real, imag Sigma_M (magnetic surface charge density)\n",nc+1,nc+2); nc+=2;
     fprintf(f,"# %i, %i  real, imag N_x (magnetic surface current density)\n",nc+1,nc+2); nc+=2;
     fprintf(f,"# %i, %i  real, imag N_y\n",nc+1,nc+2); nc+=2;
     fprintf(f,"# %i, %i  real, imag N_z\n",nc+1,nc+2); nc+=2;
     fprintf(f,"# %i      inward-directed normal Poynting flux\n",nc++);
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
     if (TransformLabel) 
      fprintf(f,"%s ",TransformLabel);
     if (IFLabel) 
      fprintf(f,"%s ",IFLabel);
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

