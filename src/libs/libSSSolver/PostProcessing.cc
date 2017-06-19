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
 * OutputModules.cc -- implementation of various electrostatic
 *                  -- calculations for the scuff-em electrostatics
 *                  -- module
 *
 * homer reid       -- 5/2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libTriInt.h>

#include "libscuff.h"
#include "SSSolver.h"
#include "StaticSubstrate.h"

namespace scuff {

#define MAXSTR 1000

/***************************************************************/
/* This is a helper routine for VisualizeFields() below.       */
/***************************************************************/
void WriteFVMesh(SSSolver *SSS, RWGSurface *S, HVector *Sigma,
                 StaticExcitation *SE, FILE *f)
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
  HMatrix *PhiE = SSS->GetFields(SE, Sigma, XMatrix);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
const char *FieldTitles[]={"Phi","Ex","Ey","Ez"};
#define NUMFIELDFUNCS 4
  for(int nff=0; nff<NUMFIELDFUNCS; nff++)
   { 
     fprintf(f,"View \"%s",FieldTitles[nff]);
     if (!SSS->SeparateOutputFiles && SSS->TransformLabel)
      fprintf(f,"{%s}",SSS->TransformLabel);
     if (!SSS->SeparateOutputFiles && SSS->ExcitationLabel)
      fprintf(f,"{%s}",SSS->ExcitationLabel);
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
           Q[nv] = PhiE->GetEntryD(VI,nff);
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
  delete PhiE;
  delete XMatrix;
}

/***************************************************************/
/* VisualizeFields() produces a color plot of the potential    */
/* and E field on a user-specified surface mesh for GMSH       */
/* visualization.                                              */
/***************************************************************/
void SSSolver::VisualizeFields(HVector *Sigma,
                               char *FVMeshFile,
                               char *OutFileBase,
                               StaticExcitation *SE,
                               char *TransFile)
{ 
  /*--------------------------------------------------------------*/
  /*- try to open user's mesh file -------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=new RWGSurface(FVMeshFile);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumFVMeshTransforms;
  GTComplex **FVMeshGTCList=ReadTransFile(TransFile, &NumFVMeshTransforms);

  if (!OutFileBase)
   OutFileBase=strdup(GetFileBase(G->GeoFileName));

  char PPFileBase[MAXSTR];
  char *FVMFileBase=GetFileBase(FVMeshFile);
  if (SeparateOutputFiles && (TransformLabel && ExcitationLabel) )
   snprintf(PPFileBase,MAXSTR,"%s.%s.%s.%s",OutFileBase,TransformLabel,ExcitationLabel,FVMFileBase);
  else if (SeparateOutputFiles && TransformLabel)
   snprintf(PPFileBase,MAXSTR,"%s.%s.%s",OutFileBase,TransformLabel,FVMFileBase);
  else if (SeparateOutputFiles && ExcitationLabel)
   snprintf(PPFileBase,MAXSTR,"%s.%s.%s",OutFileBase,ExcitationLabel,FVMFileBase);
  else
   snprintf(PPFileBase,MAXSTR,"%s.%s",OutFileBase,FVMFileBase);
  for(int nt=0; nt<NumFVMeshTransforms; nt++)
   {
     GTComplex *GTC=FVMeshGTCList[nt];
     GTransformation *GT=GTC->GT;
     char *Tag = GTC->Tag;
     char PPFileName[MAXSTR];
     if (NumFVMeshTransforms>1)
      { 
        snprintf(PPFileName,MAXSTR,"%s.%s.pp",PPFileBase,Tag);
        Log("Creating flux plot for surface %s, transform %s...",FVMeshFile,Tag);
      }
     else
      {  snprintf(PPFileName,100,"%s.pp",PPFileBase);
         Log("Creating flux plot for surface %s...",FVMeshFile);
      };
     FILE *f=fopen(PPFileName,"a");
     if (!f) 
      { Warn("could not open field visualization file %s",PPFileName);
       continue;
      };

     if (GT) S->Transform(GT);
     WriteFVMesh(this, S, Sigma, SE, f);
     if (GT) S->UnTransform();

     fclose(f);
   };

  delete S;
  DestroyGTCList(FVMeshGTCList,NumFVMeshTransforms);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *SSSolver::GetCapacitanceMatrix(HMatrix *M,
                                        HVector *Sigma,
                                        HMatrix *CMatrix)
{
  int NS = G->NumSurfaces;
  Log("Computing capacitance matrix for %s",G->GeoFileName);

  /*--------------------------------------------------------------*/
  /*- (re)allocate capacitance matrix as necessary ---------------*/
  /*--------------------------------------------------------------*/
  int NCS=0; // number of conducting surfaces 
  for(int ns=0; ns<NS; ns++)
   if (G->Surfaces[ns]->IsPEC)
    NCS++;
  if (NCS==0)
   { Warn("No conducting surfaces! Aborting capacitance calculation.");
     return 0;
   };
  if (CMatrix==0 || CMatrix->NR!=NCS || CMatrix->NC!=NCS )
   { if (CMatrix) delete CMatrix;
     CMatrix=0;
   };
  if (CMatrix==0)
   CMatrix = new HMatrix(NCS, NCS);

  /*--------------------------------------------------------------*/
  /*- allocate BEM matrix and RHS vector if necessary             */
  /*--------------------------------------------------------------*/
  bool OwnsM = (M==0);
  if (OwnsM)
   { M=AssembleBEMMatrix();
     M->LUFactorize();
   };

  bool OwnsSigma = (Sigma==0);
  if (OwnsSigma)
   Sigma=AllocateRHSVector();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double *Potentials = new double[NS];
  HMatrix *QP = new HMatrix(NS, 4);
  for(int ns=0, ncs=-1; ns<NS; ns++)
   { 
     if ( !(G->Surfaces[ns]->IsPEC) )
      continue;
     ncs++;

     memset(Potentials, 0, NS*sizeof(double));
     Potentials[ns]=1.0;
     AssembleRHSVector(Potentials, Sigma);
     M->LUSolve(Sigma);
     GetCartesianMoments(Sigma, QP);

     for(int nsp=0, ncsp=-1; nsp<NS; nsp++)
      { if ( !(G->Surfaces[nsp]->IsPEC) )
         continue;
        ncsp++;
        CMatrix->SetEntry(ncsp, ncs, QP->GetEntry(nsp,0));
      };
   };
  delete[] Potentials;
  delete QP;

  if (OwnsM) delete M;
  if (OwnsSigma) delete Sigma;

  return CMatrix;
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
void SSSolver::PlotChargeDensity(HVector *Sigma, char *BaseFileName)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char FileName[MAXSTR];
  strncpy(FileName, BaseFileName, MAXSTR-1);
  if (SeparateOutputFiles && (TransformLabel || ExcitationLabel) )
   { char Extension[MAXSTR];
     strncpy(Extension, GetFileExtension(FileName), MAXSTR-1);
     RemoveExtension(FileName);
     if (TransformLabel && ExcitationLabel)
      snprintf(FileName,MAXSTR,"%s.%s.%s.%s",FileName,TransformLabel,ExcitationLabel,Extension);
     else if (TransformLabel)
      snprintf(FileName,MAXSTR,"%s.%s.%s",FileName,TransformLabel,Extension);
     else if (ExcitationLabel)
      snprintf(FileName,MAXSTR,"%s.%s.%s",FileName,ExcitationLabel,Extension);
   };

  FILE *f=fopen(FileName,"a");
  if (!f) return;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S;
  int np;
  fprintf(f,"View \"Surface Charge Density");
  if (!SeparateOutputFiles && TransformLabel)
   fprintf(f,"{%s}",TransformLabel);
  if (!SeparateOutputFiles && ExcitationLabel)
   fprintf(f,"{%s}",ExcitationLabel);
  fprintf(f,"\" {\n");
  for(int ns=0, nbf=0; ns<G->NumSurfaces; ns++)
   for(S=G->Surfaces[ns], np=0; np<S->NumPanels; np++, nbf++)
    { 
      RWGPanel *P=S->Panels[np];
      double *PV[3];
      PV[0]=S->Vertices + 3*P->VI[0];
      PV[1]=S->Vertices + 3*P->VI[1];
      PV[2]=S->Vertices + 3*P->VI[2];
      double Val = Sigma->GetEntryD( nbf );
      fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                   PV[0][0], PV[0][1], PV[0][2],
                   PV[1][0], PV[1][1], PV[1][2],
                   PV[2][0], PV[2][1], PV[2][2],
                   Val, Val, Val);
   }; 
  fprintf(f,"};\n");
  fclose(f);
 
}

} // namespace scuff
