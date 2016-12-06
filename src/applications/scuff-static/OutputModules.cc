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
  if (CMatrix==0 || CMatrix->NR!=NS || CMatrix->NC!=NS)
   { if (CMatrix) delete CMatrix;
     CMatrix=0;
   };
  if (CMatrix==0)
   CMatrix = new HMatrix(NS, NS);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double *Potentials = new double[NS];
  HMatrix *QP = new HMatrix(NS, 4);
  for(int ns=0; ns<NS; ns++)
   { 
     memset(Potentials, 0, NS*sizeof(double));
     Potentials[ns]=1.0;
     SSS->AssembleRHSVector(Potentials, 0, 0, Sigma);
     M->LUSolve(Sigma);
     SSS->GetCartesianMoments(Sigma, QP);
     for(int nsp=0; nsp<NS; nsp++)
      CMatrix->SetEntry(nsp, ns, QP->GetEntry(nsp,0));
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

  /*--------------------------------------------------------------*/
  /*- write file header the first time ---------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=fopen(CapFile,"a");
  static bool WroteHeader=false;
  int NS = SSS->G->NumSurfaces;
  if (WroteHeader==false)
   { WroteHeader=true;
     fprintf(f,"# scuff-static run on %s (%s)",GetHostName(),GetTimeString());
     fprintf(f,"# data file columns: \n");
     int nc=1;
     if (SSS->TransformLabel)
      fprintf(f,"# %02i transformation label\n",nc++);
     for(int p=0; p<NS; p++)
      for(int q=p; q<NS; q++)
       fprintf(f,"# %02i: C_{%i,%i} \n",nc++,p,q);
   };

  /*--------------------------------------------------------------*/
  /*- write data -------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SSS->TransformLabel)
   fprintf(f,"%s ",SSS->TransformLabel);
  for(int nr=0; nr<NS; nr++)
   for(int nc=nr; nc<NS; nc++)
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
/***************************************************************/
/***************************************************************/
void WriteFields(SSSolver *SSS, HMatrix *M, HVector *Sigma,
                 char *PotFile, char *PhiExt, int ConstFieldDirection,
                 char *PlotFile, char **EPFiles, int nEPFiles, char *FileBase)
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

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/ 
  if (PlotFile)
   SSS->PlotChargeDensity(Sigma, PlotFile, 0);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *TransformLabel = SSS->TransformLabel;
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
