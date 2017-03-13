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
 * OutputModules.cc  -- output modules for scuff-static 
 *
 * homer reid        -- 10/2006 -- 11/2013
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <libscuff.h>
#include <libSpherical.h>
#include <SSSolver.h>
#include <cmatheval.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct USFData 
 {
   void *PhiEvaluator;
   void *EEvaluator[3];
 } USFData;

// 2.0e-20
//#define DELTA 9.5367431640625e-7
#define DELTA 1.0e-6
void UserStaticField(double *x, void *UserData, double PhiE[4])
{
  USFData *SFD = (USFData *)UserData;

  static const char *VariableNames[3] = { "x", "y", "z" };
  cdouble VariableValues[3];
  VariableValues[0] = cdouble(x[0], 0.0 );
  VariableValues[1] = cdouble(x[1], 0.0 );
  VariableValues[2] = cdouble(x[2], 0.0 );

  memset(PhiE, 0, 4*sizeof(double));
  if ( SFD==0 || SFD->PhiEvaluator==0 ) return;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  PhiE[0] = real( cevaluator_evaluate(SFD->PhiEvaluator, 3, 
                                      const_cast<char **>(VariableNames), 
                                      VariableValues) );

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int Mu=0; Mu<3; Mu++)
   { if (SFD->EEvaluator[Mu])
      { PhiE[1+Mu] = real( cevaluator_evaluate(SFD->EEvaluator[Mu], 3, 
                                               const_cast<char **>(VariableNames), 
                                               VariableValues) 
                         );
      }
     else
      { VariableValues[Mu] *= (1.0 + DELTA);
        double Temp = real( cevaluator_evaluate(SFD->PhiEvaluator, 3,
                                                const_cast<char **>(VariableNames), 
                                                VariableValues) 
                          );
        VariableValues[Mu] /= (1.0 + DELTA);
        PhiE[1+Mu] = -(Temp - PhiE[0]) / DELTA; 
      };
   };
   
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PhiEConstant(double *x, void *UserData, double PhiE[4])
{ 
  memset(PhiE, 0, 4.0*sizeof(double));
  int Direction = *((int *)UserData);
  PhiE[0]             = -x[Direction];
  PhiE[1 + Direction] = 1.0;
} 

/***************************************************************/
/* (real-valued) spherical harmonic incident field *************/
/***************************************************************/
typedef struct PhiESphericalData
 { 
   int l, m;

 } PhiESphericalData;

void PhiESpherical(double *x, void *UserData, double PhiE[4])
{
  PhiESphericalData *PESD = (PhiESphericalData *)UserData;
  int l = PESD->l;
  int m = PESD->m;

  double r, Theta, Phi;
  CoordinateC2S(x, &r, &Theta, &Phi);

  PhiE[0] = pow(r,l)*GetRealYlm(l,m,Theta,Phi);

  // FIXME
  PhiE[1] = PhiE[2] = PhiE[3] = 0.0;
}

/***************************************************************/
/* Parse a "potential" file to read a list of conductor        */
/* potentials.                                                 */
/* The file should look something like                         */
/*                                                             */
/*  BottomSurface 3.2                                          */
/*  TopSurface    1.4                                          */
/*                                                             */
/* where the first string is the name of a PEC object/surface  */
/* specified in the .scuffgeo file, and the second string is   */
/* the potential in volts at which that object will be held    */
/* in the scuff-static calculation.                            */
/***************************************************************/
void ParsePotentialFile(RWGGeometry *G, char *PotFile, double *Potentials)
{ 
  FILE *f=fopen(PotFile,"r");
  if (!f)
   ErrExit("could not open file %s",PotFile);

  int LineNum=0;
  char Line[100];
  while (fgets(Line,100,f))
   { 
     LineNum++;

     char *Tokens[3];
     int NumTokens=Tokenize(Line, Tokens, 3);

     if (NumTokens==0 || Tokens[0][0]=='#') 
      continue; // skip blank lines and comments

     if (NumTokens!=2) 
      ErrExit("%s:%i: syntax error",PotFile,LineNum);

     int ns;
     if (G->GetSurfaceByLabel(Tokens[0],&ns)==0)
      ErrExit("%s:%i: unknown surface",PotFile,LineNum,Tokens[0]);
     if ( ! (G->Surfaces[ns]->IsPEC) )
      ErrExit("%s:%i: attempt to assign potential to non-PEC surface %s",Tokens[0]);

     double V;
     if (1!=sscanf(Tokens[1],"%le",&V))
      ErrExit("%s:%i: invalid potential specification",PotFile,LineNum);

     Potentials[ns]=V;
     Log("Setting potential of surface %s to %e volts.\n",Tokens[0],V);
   }; 

  fclose(f); 
}

/***************************************************************/
/* solve the BEM electrostatics problem to fill in Sigma.      */
/* on entry, M is the LU-factorized BEM matrix.                */
/***************************************************************/
void Solve(SSSolver *SSS, HMatrix *M, HVector *Sigma,
           char *PotFile, char *PhiExt, int ConstFieldDirection)
{
  /***************************************************************/
  /* process user's conductor potential file if present          */
  /***************************************************************/
  RWGGeometry *G=SSS->G;
  double *Potentials = new double[G->NumSurfaces];
  memset(Potentials, 0.0, G->NumSurfaces*sizeof(double));
  if (PotFile)
   ParsePotentialFile(G, PotFile, Potentials);

  /***************************************************************/
  /* process user's external-field specification if present      */
  /***************************************************************/
  if (PhiExt)
   { 
     ErrExit("--phiexternal option not yet supported");
     /*
     USFData MyData, *D = &MyData;
     D->PhiEvaluator = cevaluator_create(PhiExt);
     D->EEvaluator[0] = D->EEvaluator[1] = D->EEvaluator[2] = 0;
     SSS->AssembleRHSVector(Potentials, UserStaticField, (void *)D, Sigma);
     */
   }
  else if ( ConstFieldDirection!=-1 )
   SSS->AssembleRHSVector(Potentials, PhiEConstant, &ConstFieldDirection, Sigma);
  else
   SSS->AssembleRHSVector(Potentials, 0, 0, Sigma);

  /***************************************************************/
  /* solve the problem *******************************************/
  /***************************************************************/
  M->LUSolve(Sigma);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WritePolarizabilities(SSSolver *SSS, HMatrix *M, HVector *Sigma, char *FileName)
{
  RWGGeometry *G = SSS->G;
  int NS = G->NumSurfaces;

  HMatrix *PolMatrix = new HMatrix(NS, 9);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *QP        = new HMatrix(NS, 4);
  for(int Mu=0; Mu<3; Mu++)
   { 
     SSS->AssembleRHSVector(0, PhiEConstant, (void *)(&Mu), Sigma);
     M->LUSolve(Sigma);
     SSS->GetCartesianMoments(Sigma, QP);
     for(int ns=0; ns<NS; ns++)
      { PolMatrix->SetEntry(ns, 0*3+Mu, QP->GetEntryD(ns,1));
        PolMatrix->SetEntry(ns, 1*3+Mu, QP->GetEntryD(ns,2));
        PolMatrix->SetEntry(ns, 2*3+Mu, QP->GetEntryD(ns,3));
      };
   };
  delete QP;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  static bool WroteHeader=false;
  FILE *f=fopen(FileName,"a");
  if (!WroteHeader)
   { WroteHeader=true;
     fprintf(f,"# scuff-static run on %s (%s)",GetHostName(),GetTimeString());
     fprintf(f,"# data file columns: \n");
     int nc=1;
     if (SSS->TransformLabel)
      fprintf(f,"# %02i: Transformation \n",nc++);
     fprintf(f,"# %02i: object label \n",nc++);
     fprintf(f,"# %02i: alpha_{xx} \n",nc++);
     fprintf(f,"# %02i: alpha_{xy} \n",nc++);
     fprintf(f,"# %02i: alpha_{xz} \n",nc++);
     fprintf(f,"# %02i: alpha_{yx} \n",nc++);
     fprintf(f,"# %02i: alpha_{yy} \n",nc++);
     fprintf(f,"# %02i: alpha_{yz} \n",nc++);
     fprintf(f,"# %02i: alpha_{zx} \n",nc++);
     fprintf(f,"# %02i: alpha_{zy} \n",nc++);
     fprintf(f,"# %02i: alpha_{zz} \n",nc++);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SSS->TransformLabel) fprintf(f,"%s ",SSS->TransformLabel);
  for(int ns=0; ns<NS; ns++)
   fprintf(f,"%s %e %e %e %e %e %e %e %e %e \n",
              G->Surfaces[ns]->Label,
              PolMatrix->GetEntryD(ns,0),
              PolMatrix->GetEntryD(ns,1),
              PolMatrix->GetEntryD(ns,2),
              PolMatrix->GetEntryD(ns,3),
              PolMatrix->GetEntryD(ns,4),
              PolMatrix->GetEntryD(ns,5),
              PolMatrix->GetEntryD(ns,6),
              PolMatrix->GetEntryD(ns,7),
              PolMatrix->GetEntryD(ns,8));
;
  delete PolMatrix;
 	
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetCapacitanceMatrix(SSSolver *SSS, HMatrix *M,
                              HVector *Sigma, HMatrix *CMatrix)
{
  RWGGeometry *G = SSS->G;
  int NS = G->NumSurfaces;

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
     SSS->AssembleRHSVector(Potentials, 0, 0, Sigma);
     M->LUSolve(Sigma);
     SSS->GetCartesianMoments(Sigma, QP);

     for(int nsp=0, ncsp=-1; nsp<NS; nsp++)
      { if ( !(G->Surfaces[nsp]->IsPEC) )
         continue;
        ncsp++;
        CMatrix->SetEntry(ncsp, ncs, QP->GetEntry(nsp,0));
      };
   };
  delete[] Potentials;
  delete QP;

  return CMatrix;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteCapacitanceMatrix(SSSolver *SSS, HMatrix *M,
                            HVector *Sigma, char *CapFile)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *CapMatrix=GetCapacitanceMatrix(SSS, M, Sigma, 0);
  if (!CapMatrix) return;

  /*--------------------------------------------------------------*/
  /*- write file header the first time ---------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=fopen(CapFile,"a");
  static bool WroteHeader=false;
  RWGGeometry *G=SSS->G;
  if (WroteHeader==false)
   { WroteHeader=true;
     fprintf(f,"# scuff-static run on %s (%s)",GetHostName(),GetTimeString());
     fprintf(f,"# indices of conducting surfaces: ");
     fprintf(f,"# data file columns: \n");
     int NCS=0;
     for(int ns=0; ns<G->NumSurfaces; ns++)
      if (G->Surfaces[ns]->IsPEC)
       fprintf(f,"# %i %s\n",NCS++,G->Surfaces[ns]->Label);
     int nc=1;
     if (SSS->TransformLabel)
      fprintf(f,"# %02i transformation label\n",nc++);
     for(int p=0; p<NCS; p++)
      for(int q=p; q<NCS; q++)
       fprintf(f,"# %02i: C_{%i,%i} \n",nc++,p,q);
   };

  /*--------------------------------------------------------------*/
  /*- write data -------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SSS->TransformLabel)
   fprintf(f,"%s ",SSS->TransformLabel);
  for(int nr=0; nr<CapMatrix->NR; nr++)
   for(int nc=nr; nc<CapMatrix->NC; nc++)
    fprintf(f,"%e ",CapMatrix->GetEntryD(nr,nc));
  fprintf(f,"\n");
  fclose(f);

  delete CapMatrix;
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteCMatrix(SSSolver *SSS, HMatrix *M, HVector *Sigma, 
                  int lMax, char *TextFileName, char *HDF5FileName)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NAlpha = (lMax+1)*(lMax+1);
  HVector *Moments = new HVector(NAlpha);
  HMatrix *CMatrix=new HMatrix(NAlpha, NAlpha);
  for(int l=0, Alpha=0; l<=lMax; l++)
   for(int m=-l; m<=l; m++, Alpha++)
    { 
      /*--------------------------------------------------------------*/
      /* setup and solve the electrostatics problem with an (l,m)     */
      /* spherical-harmonic incident field                            */
      /*--------------------------------------------------------------*/
      PhiESphericalData MyPhiESphericalData, *PESD=&MyPhiESphericalData;
      PESD->l = l;
      PESD->m = m;
      SSS->AssembleRHSVector(0, PhiESpherical, (void *)PESD, Sigma);
      M->LUSolve(Sigma);

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      SSS->GetSphericalMoments(Sigma, lMax, Moments);
      for(int lp=0, AlphaP=0; lp<=lMax; lp++)
       for(int mp=-lp; mp<=lp; mp++, AlphaP++)
        CMatrix->SetEntry(AlphaP,Alpha,Moments->GetEntry(AlphaP));

    };
  delete Moments;

  /***************************************************************/
  /* 20131219 there is unquestionably a very much more efficient */
  /* way to do this...                                           */
  /***************************************************************/
  HMatrix *Gamma=new HMatrix(NAlpha, NAlpha, LHM_COMPLEX);
  Gamma->Zero();
  cdouble OORT2(0.7071067811865475, 0.0);
  cdouble OORT2I(0.0,-0.7071067811865475);
  for(int l=0, Alpha=0; l<=lMax; l++)
   for(int m=-l; m<=l; m++, Alpha++)
    { 
      double Sign = (m%2) ? -1.0 : 1.0;
      if (m>0) 
       { Gamma->SetEntry( Alpha, LM2ALPHA(l,  m), OORT2);
         Gamma->SetEntry( Alpha, LM2ALPHA(l, -m), Sign*OORT2);
       }
      else if (m<0) 
       { Gamma->SetEntry( Alpha, LM2ALPHA(l,  m),  OORT2I);
         Gamma->SetEntry( Alpha, LM2ALPHA(l, -m), -Sign*OORT2I);
       }
      else
       Gamma->SetEntry( Alpha, Alpha, 1.0 );
    };

  HMatrix *CBarMatrix=new HMatrix(NAlpha, NAlpha, LHM_COMPLEX);
  for(int Alpha=0; Alpha<NAlpha; Alpha++)
   for(int AlphaP=0; AlphaP<NAlpha; AlphaP++)
    { cdouble Sum=0.0;
      for(int Beta=0; Beta<NAlpha; Beta++)
       for(int BetaP=0; BetaP<NAlpha; BetaP++)
        Sum +=  Gamma->GetEntry(Alpha,Beta)
               *CMatrix->GetEntry(Beta,BetaP)
               *conj(Gamma->GetEntry(AlphaP,BetaP));
      CBarMatrix->SetEntry(Alpha, AlphaP, Sum);
    };
    
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (TextFileName)
   { FILE *f=fopen(TextFileName,"w");
     if (SSS->TransformLabel) fprintf(f,"%s ",SSS->TransformLabel);
     for(int l=0, Alpha=0; l<=lMax; l++)
      for(int m=-l; m<=l; m++, Alpha++)
       for(int lp=0, AlphaP=0; lp<=lMax; lp++)
        for(int mp=-lp; mp<=lp; mp++, AlphaP++)
         fprintf(f,"%e %e ", real(CBarMatrix->GetEntry(Alpha,AlphaP)),
                             imag(CBarMatrix->GetEntry(Alpha,AlphaP)));
     fprintf(f,"\n");
     fclose(f);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (HDF5FileName)
   { 
     void *pHC = HMatrix::OpenHDF5Context(HDF5FileName);

     for(int Alpha=0; Alpha<NAlpha; Alpha++)
      for(int AlphaP=0; AlphaP<NAlpha; AlphaP++)
       CMatrix->SetEntry(Alpha, AlphaP, 
                          real(CBarMatrix->GetEntry(Alpha, AlphaP)));
     CMatrix->ExportToHDF5(pHC,"CReal");

     for(int Alpha=0; Alpha<NAlpha; Alpha++)
      for(int AlphaP=0; AlphaP<NAlpha; AlphaP++)
       CMatrix->SetEntry(Alpha, AlphaP, 
                          imag(CBarMatrix->GetEntry(Alpha, AlphaP)));
     CMatrix->ExportToHDF5(pHC,"CImag");

     HMatrix::CloseHDF5Context(pHC);
   };

  delete CMatrix;
  delete CBarMatrix;
  delete Gamma;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFields(SSSolver *SSS, HVector *Sigma,
                 char *PhiExt, int ConstFieldDirection,
                 char **EPFiles, int nEPFiles)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *TransformLabel = SSS->TransformLabel;
  char *FileBase       = SSS->FileBase;
  for(int nepf=0; nepf<nEPFiles; nepf++)
   {
     HMatrix *X = new HMatrix(EPFiles[nepf]);
     if (X->ErrMsg)
      ErrExit(X->ErrMsg);

     HMatrix *PhiE;
     if (PhiExt)
      ErrExit("%s:%i: internal error ",__FILE__,__LINE__);
     if ( ConstFieldDirection!=-1 )
      PhiE = SSS->GetFields(PhiEConstant, &ConstFieldDirection, Sigma, X, 0);
     else
      PhiE = SSS->GetFields(0, 0, Sigma, X, 0);

     FILE *f=vfopen("%s.%s.out","w",FileBase,GetFileBase(EPFiles[nepf]));
     fprintf(f,"# scuff-static run on %s (%s)",GetHostName(),GetTimeString());
     fprintf(f,"# data file columns: \n");
     int nc=1;
     if (TransformLabel) 
      fprintf(f,"# %i: transform label\n",nc++);
     fprintf(f,"# %i, %i, %i: x, y, z (evaluation point coordinates)\n",nc,nc+1,nc+2); nc+=3;
     fprintf(f,"# %i:       Phi      (electrostatic potential)\n",nc++);
     fprintf(f,"# %i, %i, %i: Ex,Ey,Ez (electrostatic field components)\n",nc,nc+1,nc+2); nc+=3;
     for(int nr=0; nr<X->NR; nr++)
      { if (TransformLabel) fprintf(f,"%s ",TransformLabel);
        fprintf(f,"%e %e %e %e %e %e %e\n",
                 X->GetEntryD(nr,0), X->GetEntryD(nr,1), X->GetEntryD(nr,2),
                 PhiE->GetEntryD(nr,0), PhiE->GetEntryD(nr,1), 
                 PhiE->GetEntryD(nr,2), PhiE->GetEntryD(nr,3));
      };
     fclose(f);
   }; // for(int nepf=0; nepf<nEPFiles; nepf++)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFVMesh(SSSolver *SSS, RWGSurface *S, HVector *Sigma,
                 char *PhiExt, int ConstFieldDirection,
                 char *TransformLabel, FILE *f)
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
  HMatrix *PhiE=0;
  if (PhiExt)
   ErrExit("%s:%i: internal error ",__FILE__,__LINE__);
  else if ( ConstFieldDirection!=-1 )
   PhiE = SSS->GetFields(PhiEConstant, &ConstFieldDirection, Sigma, XMatrix, 0);
  else
   PhiE = SSS->GetFields(0, 0, Sigma, XMatrix, 0);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
const char *FieldTitles[]={"Phi","Ex","Ey","Ez"};
#define NUMFIELDFUNCS 4
  for(int nff=0; nff<NUMFIELDFUNCS; nff++)
   { 
     fprintf(f,"View \"%s",FieldTitles[nff]);
     if (TransformLabel)
      fprintf(f,"(%s)",TransformLabel);
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
void VisualizeFields(SSSolver *SSS, HVector *Sigma,
                     char *PhiExt, int ConstFieldDirection,
                     char *FVMesh, char *TransFile)
{ 
  /*--------------------------------------------------------------*/
  /*- try to open user's mesh file -------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=new RWGSurface(FVMesh);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumFVMeshTransforms;
  GTComplex **FVMeshGTCList=ReadTransFile(TransFile, &NumFVMeshTransforms);

  char *GeoFileBase=SSS->FileBase;
  char *FVMFileBase=GetFileBase(FVMesh);
  for(int nt=0; nt<NumFVMeshTransforms; nt++)
   {
     GTComplex *GTC=FVMeshGTCList[nt];
     GTransformation *GT=GTC->GT;
     char *Tag = GTC->Tag;
     char PPFileName[100];
     if (NumFVMeshTransforms>1)
      { 
        snprintf(PPFileName,100,"%s.%s.%s.pp",GeoFileBase,FVMFileBase,Tag);
        Log("Creating flux plot for surface %s, transform %s...",FVMesh,Tag);
      }
     else
      {  snprintf(PPFileName,100,"%s.%s.pp",GeoFileBase,FVMFileBase);
         Log("Creating flux plot for surface %s...",FVMesh);
      };
     FILE *f=fopen(PPFileName,"a");
     if (!f) 
      { Warn("could not open field visualization file %s",PPFileName);
       continue;
      };

     if (GT) S->Transform(GT);
     WriteFVMesh(SSS, S, Sigma, PhiExt, ConstFieldDirection, Tag, f);
     if (GT) S->UnTransform();

     fclose(f);
   };

  delete S;
  DestroyGTCList(FVMeshGTCList,NumFVMeshTransforms);

}
