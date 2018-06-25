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
 * OutputModules.cc -- calculations of various types of output quantities
 *                  -- for the scuff-EM scuffSolver module
 * 
 * homer reid       -- 9/2011 - 4/2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libIncField.h>

#include "scuffSolver.h"

#define II cdouble(0.0,1.0)

using namespace scuff;

namespace scuff {
  
/***************************************************************/
/***************************************************************/
/***************************************************************/
void scuffSolver::ProcessEPFile(char *EPFile, char *OutFileName)
{
  /***************************************************************/
  /* read in the list of evaluation points                       */
  /***************************************************************/
  HMatrix *XMatrix=new HMatrix(EPFile);
  if (XMatrix->NC!=3)
   { Warn("ProcessEPFile: file %s should have exactly 3 coordinates per line (skipping)",EPFile);
     return;
   }
  XMatrix->Transpose();

  /***************************************************************/
  /* Get matrix of field components at evaluation points         */
  /***************************************************************/
  HMatrix *PFMatrix = GetFields(XMatrix);

  /***************************************************************/
  /* write field components to output file ***********************/
  /***************************************************************/
  char FileNameBuffer[100];
  if (!OutFileName) 
   { OutFileName = FileNameBuffer;
     snprintf(OutFileName,100,"%s.%s.fields",FileBase,GetFileBase(EPFile));
   }
  FILE *f=fopen(OutFileName,"w");
  if (!f) 
   { Warn("ProcessEPFile: could not open file %s (skipping)",OutFileName); 
     return;
   }

  char TimeString[200];
  time_t MyTime=time(0);
  struct tm *MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);
  fprintf(f,"# scuffSolver::ProcessEPFile ran on %s (%s)\n",getenv("HOST"),TimeString);
  fprintf(f,"# columns:\n");
  fprintf(f,"#   1-3: x, y, z\n");
  fprintf(f,"#   4,5: re, im Ex \n");
  fprintf(f,"#   6,7: re, im Ey \n");
  fprintf(f,"#   8,9: re, im Ez \n");
  fprintf(f,"# 10,11: re, im iwAx\n");
  fprintf(f,"# 12,13: re, im iwAy\n");
  fprintf(f,"# 14,15: re, im iwAz\n");
  fprintf(f,"# 16,17: re, im Phi\n");
  fprintf(f,"# 18,19: re, im -dxPhi \n");
  fprintf(f,"# 20,21: re, im -dyPhi \n");
  fprintf(f,"# 22,23: re, im -dzPhi \n");
  cdouble Omega=G->StoredOmega;
  for(int nx=0; nx<XMatrix->NC; nx++)
   { 
     double *X   = (double *)XMatrix->GetColumnPointer(nx);
     cdouble *PF = (cdouble *)PFMatrix->GetColumnPointer(nx);
     cdouble E[3], iwA[3], mdPhi[3], Phi=PF[_PF_PHI];
     for(int i=0; i<3; i++)
      { iwA[i]   = II*Omega*PF[_PF_AX    + i];
        mdPhi[i] =     -1.0*PF[_PF_DXPHI + i];
        E[i]     = iwA[i] + mdPhi[i];
      }
     fprintVec(f,X);
     fprintVec(f,E);
     fprintVec(f,iwA);
     fprintVec(f,&Phi,1);
     fprintVecCR(f,mdPhi);
   }
  fclose(f);

  delete PFMatrix;
  delete XMatrix;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static const char *DataNames[]=
   { "#RS_re Phi",       "#IS_im Phi",   "#MS_|Phi|",
     "#RV_re E",0,0,     "#IV_im E",0,0, "#MS_|E|",
     "#RS_re Sigma",     "#IS_im Sigma",
     "#RV_re K",0,0,     "#IV_im K",0,0,
     "#RS_re K dot iwA", "#IS_im K dot iwA",
     "#RS_re SigmaPhi",  "#IS_im SigmaPhi",
     "#RS_re Kdot E",    "#RS_im Kdot E"
   };
#define NUMDATA_WITHPSDS    (sizeof(DataNames)/sizeof(DataNames[0]))
#define NUMDATA_WITHOUTPSDS 10
  
HMatrix *RFFieldsMDF(void *UserData, HMatrix *XMatrix, const char ***pDataNames)
{
  scuffSolver *Solver      = (scuffSolver *)UserData;

  RWGGeometry *G        = Solver->G;
  cdouble Omega         = G->StoredOmega;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double zMax=XMatrix->GetEntryD(2,0), zMin=zMax;
  int NX = XMatrix->NC;
  for(int nx=1; nx<NX; nx++)
   { zMax = fmax(zMax, XMatrix->GetEntryD(2,nx) );
     zMin = fmin(zMin, XMatrix->GetEntryD(2,nx) );
   }
  double zSubstrate = 0.0;
  if (G->Substrate && (G->Substrate->NumInterfaces==1))
   zSubstrate = G->Substrate->zInterface[0];
  bool GetPSDs = ( EqualFloat(zMax, zMin) && EqualFloat(zMax, zSubstrate) );
  int NumData = GetPSDs ? NUMDATA_WITHPSDS : NUMDATA_WITHOUTPSDS;
  *pDataNames = DataNames;

  Log("Computing FVMesh data: %i quantities at %i points",NumData,NX);
  HMatrix *PFMatrix = Solver->GetFields(XMatrix);

  /***************************************************************/
  /* Q[0,1,2][nd] = {total, BF contribution, Port contribution}  */
  /*                to data quantity #nd                         */
  /***************************************************************/
  cdouble iw=II*Omega;
  HMatrix *DataMatrix=new HMatrix(NX, NumData, LHM_COMPLEX);
  for(int nx=0; nx<NX; nx++)
   { 
     cdouble iwA[3], Phi, mdPhi[3];
     iwA[0]   = iw*PFMatrix->GetEntry(_PF_AX,nx);
     iwA[1]   = iw*PFMatrix->GetEntry(_PF_AY,nx);
     iwA[2]   = iw*PFMatrix->GetEntry(_PF_AZ,nx);
     Phi      =    PFMatrix->GetEntry(_PF_PHI,nx);
     mdPhi[0] = -1.0*PFMatrix->GetEntry(_PF_DXPHI,nx);
     mdPhi[1] = -1.0*PFMatrix->GetEntry(_PF_DYPHI,nx);
     mdPhi[2] = -1.0*PFMatrix->GetEntry(_PF_DZPHI,nx);

     cdouble E[3];
     VecAdd(iwA, mdPhi, E);
     double EMagnitude = sqrt( norm(E[0]) + norm(E[1]) + norm(E[2]) );

     int nc=0;
     DataMatrix->SetEntry(nx,nc++,Phi);
     DataMatrix->SetEntry(nx,nc++,Phi);
     DataMatrix->SetEntry(nx,nc++,Phi);
     DataMatrix->SetEntry(nx,nc++,E[0]);
     DataMatrix->SetEntry(nx,nc++,E[1]);
     DataMatrix->SetEntry(nx,nc++,E[2]);
     DataMatrix->SetEntry(nx,nc++,E[0]);
     DataMatrix->SetEntry(nx,nc++,E[1]);
     DataMatrix->SetEntry(nx,nc++,E[2]);
     DataMatrix->SetEntry(nx,nc++,EMagnitude);

     if (GetPSDs)
      { 
        double *X = (double *)XMatrix->GetColumnPointer(nx);
        cdouble iwSigmaK[4];
        Solver->EvalSourceDistribution(X, iwSigmaK);

        cdouble K[3], Sigma, KdotiwA, SigmaPhi, KdotE;
        Sigma    = iwSigmaK[0] / iw;
        K[0]     = iwSigmaK[1];
        K[1]     = iwSigmaK[2];
        K[2]     = iwSigmaK[3];
        KdotiwA  = K[0]*iwA[0] + K[1]*iwA[1] + K[2]*iwA[2];
        SigmaPhi = Sigma*Phi;
        KdotE    = K[0]*E[0]   + K[1]*E[1]   + K[2]*E[2];

        DataMatrix->SetEntry(nx,nc++,Sigma);
        DataMatrix->SetEntry(nx,nc++,Sigma);
        DataMatrix->SetEntry(nx,nc++,K[0]);
        DataMatrix->SetEntry(nx,nc++,K[1]);
        DataMatrix->SetEntry(nx,nc++,K[2]);
        DataMatrix->SetEntry(nx,nc++,K[0]);
        DataMatrix->SetEntry(nx,nc++,K[1]);
        DataMatrix->SetEntry(nx,nc++,K[2]);
        DataMatrix->SetEntry(nx,nc++,KdotiwA);
        DataMatrix->SetEntry(nx,nc++,KdotiwA);
        DataMatrix->SetEntry(nx,nc++,SigmaPhi);
        DataMatrix->SetEntry(nx,nc++,SigmaPhi);
        DataMatrix->SetEntry(nx,nc++,KdotE);
        DataMatrix->SetEntry(nx,nc++,KdotE);
      }
   }

  delete PFMatrix;
  return DataMatrix;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
const char *PPOptions=
 //"General.ScaleX   = 3;\n"
 //"General.ScaleY   = 3;\n"
 //"General.ScaleZ   = 3;\n"
 //"General.Axes     = 0;\n"
 "View.Light   = 0;\n"
 "View.Visible = 0;\n";

HMatrix *scuffSolver::ProcessFVMesh(char *FVMesh, char *TransFile, char *OutFileBase)
{
  char FileBaseBuffer[100];
  if (!OutFileBase) 
   { OutFileBase = FileBaseBuffer;
     snprintf(OutFileBase,100,"%s.%s",FileBase,GetFileBase(FVMesh));
    }
  HMatrix *Integrals = MakeMeshPlot(RFFieldsMDF, (void *)this, FVMesh, TransFile, OutFileBase, PPOptions);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#if 0
  char *s = getenv("SCUFF_WRITE_MESH_INTEGRALS");
  if (s && s[0]=='1')
   { 
     bool WithPSDs = (Integrals->NR == NUMDATA_WITHPSDS);
     const char **DataNames = WithPSDs ? DataNamesWithPSDs : DataNamesWithoutPSDs;
     FILE *f=vfopen("/tmp/%s.%s.MeshIntegrals","a",FileBase,FVMesh);
     for(int nt=0; nt<Integrals->NC; nt++)
      { fprintf(f,"\n Data quantities at f=%g, transform %i:\n",real(OMEGA2FREQ*Omega),nt);
        for(int nd=0; nd<Integrals->NR; nd++)
         fprintf(f,"%+8e # %s\n",Integrals->GetEntryD(nd,nt),DataNames[nd]);
      }
     fclose(f);
   }
#endif
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  return Integrals;
}

/***************************************************************/
/* on entry: column #np of KMatrix =                           */
/*  vector of basis-function weights obtained by solving       */
/*  BEM system for unit-strength current input into port #np   */
/*                                                             */
/* Z_{pq} -= I_q * DeltaV_{p}                                  */
/*                                                             */
/* where I_q == 1                                              */
/*       DeltaV_p = V_p^+ -V_p^-                               */
/*       V_p^\pm = average value of Phi at \pm terminal of     */
/*                 port p due to unit current input to port q  */
/***************************************************************/
void scuffSolver::AddMinusIdVTermsToZMatrix(HMatrix *KMatrix, HMatrix *ZMatrix)
{
  int NBF           = G->TotalBFs;
  int NumPortEdges  = PortList->PortEdges.size();
  HMatrix *XMatrix  = new HMatrix(3, NumPortEdges, LHM_REAL);
  for(int npe=0; npe<NumPortEdges; npe++)
   {
     RWGPortEdge *PE = PortList->PortEdges[npe];
     RWGSurface *S   = G->Surfaces[PE->ns];
     RWGEdge *E      = S->GetEdgeByIndex(PE->ne);
     XMatrix->SetEntriesD(":",npe,E->Centroid);
   }
  HMatrix *PFMatrix = new HMatrix(NPFC, NumPortEdges, LHM_COMPLEX);
  
  for(int SourcePort=0; SourcePort<NumPorts; SourcePort++)
   { 
     // get fields at centroids of all port edges due to unit-strength current into SourcePort
     HVector KNSource(NBF, LHM_COMPLEX, KMatrix->GetColumnPointer(SourcePort));
     if (!KN) KN=G->AllocateRHSVector();
     KN->Copy(&KNSource);
     if (PortCurrents==0) PortCurrents=new cdouble[NumPorts];
     memset(PortCurrents, 0, NumPorts*sizeof(cdouble));
     PortCurrents[SourcePort]=1.0;
     GetFields(XMatrix, PFMatrix);
     
     // sum contributions to mean voltage gaps across all destination ports
     for(int npe=0; npe<NumPortEdges; npe++)
      { RWGPortEdge *DestPE = PortList->PortEdges[npe];
        int DestPort        = DestPE->nPort;
        int DestPol         = DestPE->Pol;
        double DestSign     = (DestPol==_PLUS ? 1.0 : -1.0);
	double NDest        = (double)(PortList->Ports[DestPort]->PortEdges[DestPol].size());
        cdouble PhiDest     = PFMatrix->GetEntry(_PF_PHI, npe);
        ZMatrix->AddEntry(DestPort, SourcePort, -1.0*DestSign*PhiDest/NDest);
      }
    }

  delete XMatrix;
  delete PFMatrix;
}

/***************************************************************/
/* compute the NPxNP impedance matrix (NP=number of ports).    */
/***************************************************************/
HMatrix *scuffSolver::GetZMatrix(HMatrix *ZMatrix, HMatrix **pZTerms)
{
  HMatrix *ZTerms[3];
  if (pZTerms)
   { ZTerms[0] = pZTerms[0] = CheckHMatrix(pZTerms[0], NumPorts, NumPorts, LHM_COMPLEX);
     ZTerms[1] = pZTerms[1] = CheckHMatrix(pZTerms[1], NumPorts, NumPorts, LHM_COMPLEX);
     ZTerms[2] = pZTerms[2] = CheckHMatrix(pZTerms[2], NumPorts, NumPorts, LHM_COMPLEX);
   }
  else
   { ZTerms[0] = new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);
     ZTerms[1] = new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);
     ZTerms[2] = new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);
   }

  AssemblePortBFInteractionMatrix(G->StoredOmega);

  // Term1 = -(R^T * W * R) where R = port/BF interaction matrix
  HMatrix *RMatrix=PBFIMatrix, *WRMatrix = new HMatrix(RMatrix);
  M->LUSolve(WRMatrix);
  RMatrix->Multiply(WRMatrix, ZTerms[0], "--transA T");
  ZTerms[0]->Scale(-1.0*ZVAC);

  // Term2 = port--port interaction matrix
  ZTerms[1]->Copy(PPIMatrix);
  ZTerms[1]->Scale(1.0*ZVAC);

  // Term3 = matrix of port voltage gaps
  ZTerms[2]->Zero();
  AddMinusIdVTermsToZMatrix(WRMatrix, ZTerms[2]);

  delete WRMatrix;

  ZMatrix=CheckHMatrix(ZMatrix, NumPorts, NumPorts, LHM_COMPLEX, "GetZMatrix");
  ZMatrix->Zero();
  ZMatrix->Add(ZTerms[0]);
  ZMatrix->Add(ZTerms[1]);
  ZMatrix->Add(ZTerms[2]);

  char *TermFile=0;
  if (CheckEnv("SCUFF_ZMATRIX_TERMFILE",&TermFile))
   { 
     HMatrix *ZSTerms[8];
     ZSTerms[0]=ZMatrix;
     ZSTerms[1]=ZTerms[0];
     ZSTerms[2]=ZTerms[1];
     ZSTerms[3]=ZTerms[2];
     ZSTerms[4]=Z2S(ZMatrix);
     ZSTerms[5]=Z2S(ZTerms[0]);
     ZSTerms[6]=Z2S(ZTerms[1]);
     ZSTerms[7]=Z2S(ZTerms[1]);

     static char Mode[2]="w";
     FILE *f=vfopen(TermFile,Mode,FileBase);
     Mode[0]='a';
     fprintf(f,"%e ",real(G->StoredOmega*OMEGA2FREQ));
     for(int nt=0; nt<8; nt++)
      for(int npp=0; npp<NumPorts*NumPorts; npp++)
       fprintf(f,"%e %e ",abs(ZSTerms[nt]->ZM[npp]), arg(ZSTerms[nt]->ZM[npp]));
     fprintf(f,"\n");
     fclose(f);
     delete ZSTerms[4];
     delete ZSTerms[5];
     delete ZSTerms[6];
     delete ZSTerms[7];
   }

  if (!pZTerms)
   { delete ZTerms[0];
     delete ZTerms[1];
     delete ZTerms[2];
   }

  return ZMatrix;
}

/***************************************************************/
/* routines for transforming impedance parameters into         */
/* scattering parameters                                       */
/*                                                             */
/* formulas: S = (Z-z0*I)*(Z+z0*I)^{-1}                        */
/*           Z = z0* (I-S)^{-1} * (I+S)                        */
/*                                                             */
/* where z0 = characteristic impedance, I=identity matrix      */
/*                                                             */
/* If the second argument is NULL, a new matrix is allocated.  */
/* If the second argument is equal to the first argument, the  */
/* routines operate in-place.                                  */
/***************************************************************/
HMatrix *scuffSolver::Z2S(HMatrix *Z, HMatrix *S, double ZCharacteristic)
{ 
  HMatrix *ZpZC = new HMatrix(Z);
  if (S==0) 
   S = new HMatrix(Z);
  else if (S!=Z)
   S->InsertBlock(Z,0,0);
  for(int nr=0; nr<ZpZC->NR; nr++)
   { ZpZC->AddEntry(nr, nr,  ZCharacteristic);
        S->AddEntry(nr, nr, -ZCharacteristic);
   }
  ZpZC->LUFactorize();
  ZpZC->LUSolve(S);
  delete ZpZC;
  return S;
}

HMatrix *scuffSolver::S2Z(HMatrix *S, HMatrix *Z, double ZCharacteristic)
{
  // 'OmS = 'one minus S'
  HMatrix *OmS=new HMatrix(S);
  HMatrix *OpS=new HMatrix(S);
  OmS->Scale(-1.0);
  for(int nr=0; nr<S->NR; nr++)
   { 
     OpS->AddEntry(nr, nr, 1.0);
     OmS->AddEntry(nr, nr, 1.0);
   }
  OmS->LUFactorize();
  OmS->LUInvert();
  if (Z==0) Z=new HMatrix(S);
  OmS->Multiply(OpS, Z);
  Z->Scale(ZCharacteristic);

  delete OpS;
  delete OmS;
  return Z;
}

} // namespace scuff
