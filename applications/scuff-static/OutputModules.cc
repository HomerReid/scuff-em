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
#include <StaticSolver.h>
#include <cmatheval.h>

using namespace scuff;

#define MAXSTR 1000

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFields(StaticSolver *SS, HVector *Sigma, StaticExcitation *SE,
                 char **EPFiles, int nEPFiles, char *OutFileBase)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *TransformLabel  = SS->TransformLabel;
  char *ExcitationLabel = SS->ExcitationLabel;
  bool Separate         = SS->SeparateOutputFiles;
  for(int nepf=0; nepf<nEPFiles; nepf++)
   {
     HMatrix *X = new HMatrix(EPFiles[nepf]);
     if (X->ErrMsg)
      ErrExit(X->ErrMsg);

     HMatrix *PhiE = SS->GetFields(SE, Sigma, X);

     char FileName[MAXSTR];
     if (Separate && TransformLabel && ExcitationLabel)
      snprintf(FileName, MAXSTR, "%s.%s.%s.%s.out",OutFileBase,TransformLabel,ExcitationLabel,GetFileBase(EPFiles[nepf]));
     else if (Separate && TransformLabel)
      snprintf(FileName, MAXSTR, "%s.%s.%s.out",OutFileBase,TransformLabel,GetFileBase(EPFiles[nepf]));
     else if (Separate && ExcitationLabel)
      snprintf(FileName, MAXSTR, "%s.%s.%s.out",OutFileBase,ExcitationLabel,GetFileBase(EPFiles[nepf]));
     else 
      snprintf(FileName, MAXSTR, "%s.%s.out",OutFileBase,GetFileBase(EPFiles[nepf]));
     FILE *f=fopen(FileName,Separate ? "w" : "a");
     fprintf(f,"# scuff-static " VERSIONSTRING " run on %s (%s)",GetHostName(),GetTimeString());
     fprintf(f,"# data file columns: \n");
     int nc=1;
     if (!Separate && TransformLabel) 
      fprintf(f,"# %i: transform label\n",nc++);
     if (!Separate && ExcitationLabel)
      fprintf(f,"# %i: excitation label\n",nc++);
     fprintf(f,"# %i, %i, %i: x, y, z (evaluation point coordinates)\n",nc,nc+1,nc+2); nc+=3;
     fprintf(f,"# %i:       Phi      (electrostatic potential)\n",nc++);
     fprintf(f,"# %i, %i, %i: Ex,Ey,Ez (electrostatic field components)\n",nc,nc+1,nc+2); nc+=3;
     for(int nr=0; nr<X->NR; nr++)
      { if (!Separate && TransformLabel)  fprintf(f,"%s ",TransformLabel);
        if (!Separate && ExcitationLabel) fprintf(f,"%s ",ExcitationLabel);
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
void WritePolarizabilities(StaticSolver *SS, HMatrix *M, HVector *Sigma, char *PolFile)
{
  RWGGeometry *G = SS->G;
  int NS = G->NumSurfaces;

  HMatrix *PolMatrix = new HMatrix(NS, 9);

  ConstantSFData Data;
  StaticField SFs[1] = {ConstantStaticField};
  void *SFData [1]   = {(void *)&Data};

  StaticExcitation SE;
  SE.Label=0;
  SE.Potentials=0;
  SE.SFs    = SFs;
  SE.SFData = SFData;
  SE.NumSFs=1;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *QP        = new HMatrix(NS, 4);
  for(int Mu=0; Mu<3; Mu++)
   { 
     Data.E0[0]=Data.E0[1]=Data.E0[2]=0.0;
     Data.E0[Mu]=1.0;
  
     SS->AssembleRHSVector(&SE, Sigma);
     M->LUSolve(Sigma);
     SS->GetCartesianMoments(Sigma, QP);
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
  FILE *f=fopen(PolFile,"a");
  if (!WroteHeader)
   { WroteHeader=true;
     fprintf(f,"# scuff-static " VERSIONSTRING " run on %s (%s)",GetHostName(),GetTimeString());
     fprintf(f,"# data file columns: \n");
     int nc=1;
     if (SS->TransformLabel)
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
  if (SS->TransformLabel) fprintf(f,"%s ",SS->TransformLabel);
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
void WriteCapacitanceMatrix(StaticSolver *SS, HMatrix *M,
                            HVector *Sigma, char *CapFile)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *CapMatrix=SS->GetCapacitanceMatrix(M, Sigma, 0);
  if (!CapMatrix) return;

  /*--------------------------------------------------------------*/
  /*- write file header the first time ---------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=fopen(CapFile,"a");
  static bool WroteHeader=false;
  RWGGeometry *G=SS->G;
  if (WroteHeader==false)
   { WroteHeader=true;
     fprintf(f,"# scuff-static " VERSIONSTRING " run on %s (%s)",GetHostName(),GetTimeString());
     fprintf(f,"# indices of conducting surfaces: ");
     fprintf(f,"# data file columns: \n");
     int NCS=0;
     for(int ns=0; ns<G->NumSurfaces; ns++)
      if (G->Surfaces[ns]->IsPEC)
       fprintf(f,"# %i %s\n",NCS++,G->Surfaces[ns]->Label);
     int nc=1;
     if (SS->TransformLabel)
      fprintf(f,"# %02i transformation label\n",nc++);
     for(int p=0; p<NCS; p++)
      for(int q=p; q<NCS; q++)
       fprintf(f,"# %02i: C_{%i,%i} \n",nc++,p,q);
   };

  /*--------------------------------------------------------------*/
  /*- write data -------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (SS->TransformLabel)
   fprintf(f,"%s ",SS->TransformLabel);
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
void WriteCMatrix(StaticSolver *SS, HMatrix *M, HVector *Sigma, 
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
      SphericalSFData MyData, *Data=&MyData;
      Data->l = l;
      Data->m = m;
      StaticField SFs[1]={SphericalStaticField};
      StaticExcitation SE={0,0,SFs,(void **)(&Data),1};
      SS->AssembleRHSVector(&SE, Sigma);
      M->LUSolve(Sigma);

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      SS->GetSphericalMoments(Sigma, lMax, Moments);
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
     if (SS->TransformLabel) fprintf(f,"%s ",SS->TransformLabel);
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
