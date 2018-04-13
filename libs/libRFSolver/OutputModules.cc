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
 *                  -- for the scuff-EM RFSolver module
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

#include "RFSolver.h"

#define II cdouble(0.0,1.0)

using namespace scuff;

namespace scuff {
  
/***************************************************************/
/***************************************************************/
/***************************************************************/
void ProcessEPFile(RWGGeometry *G, RWGPortList *PortList, cdouble Omega,
                   char *EPFile, HVector *KN, cdouble *PortCurrents,
                   char *FileBase)
{
  /***************************************************************/
  /* read in the list of evaluation points  **********************/
  /***************************************************************/ 
  HMatrix *XMatrix=new HMatrix(EPFile);
  XMatrix->Transpose();
  int NX = XMatrix->NC;

  HMatrix *PFContributions[2]={0,0};
  char *s=getenv("SCUFF_MOIFIELDS2");
  if (s && s[0]=='1')
   { PFContributions[0] = 0;
     PFContributions[1] = 0;
   }
  else  
   { PFContributions[0] = new HMatrix(NPFC, NX, LHM_COMPLEX);
     PFContributions[1] = new HMatrix(NPFC, NX, LHM_COMPLEX);
   }

  HMatrix *PFMatrix = GetMOIFields2(G, PortList, Omega, XMatrix, KN, PortCurrents, PFContributions);

  /***************************************************************/
  /* write field components to output file ***********************/
  /***************************************************************/
  FILE *f = vfopen("%s.%s.out","w",FileBase,GetFileBase(EPFile));
  char TimeString[200];
  time_t MyTime=time(0);
  struct tm *MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);
  fprintf(f,"# scuff-microstrip ran on %s (%s)\n",getenv("HOST"),TimeString);
  fprintf(f,"# columns:\n");
  fprintf(f,"#   1-3: x, y, z\n");
  fprintf(f,"#   4,5: Ax    (BF)\n");
  fprintf(f,"#   6,7: Ay    (BF)\n");
  fprintf(f,"#   8,9: Az    (BF)\n");
  fprintf(f,"# 10,11: Phi   (BF)\n");
  fprintf(f,"# 12,13: dxPhi (BF)\n");
  fprintf(f,"# 14,15: dyPhi (BF)\n");
  fprintf(f,"# 16,17: dzPhi (BF)\n");
  fprintf(f,"# 18,19: Ax    (Port)\n");
  fprintf(f,"# 20,21: Ay    (Port)\n");
  fprintf(f,"# 22,23: Az    (Port)\n");
  fprintf(f,"# 24,25: Phi   (Port)\n");
  fprintf(f,"# 26,27: dxPhi (Port)\n");
  fprintf(f,"# 28,29: dyPhi (Port)\n");
  fprintf(f,"# 30,31: dzPhi (Port)\n");
  for(int nx=0; nx<NX; nx++)
   { 
     fprintVec(f,(double *)XMatrix->GetColumnPointer(nx));

     if (PFContributions[0])
      { fprintVec(f,(cdouble *)PFContributions[0]->GetColumnPointer(nx),NPFC);
        fprintVecCR(f,(cdouble *)PFContributions[1]->GetColumnPointer(nx),NPFC);
      }
     else
      fprintVecCR(f,(cdouble *)PFMatrix->GetColumnPointer(nx),NPFC);

   }
  fclose(f);
 
  if (PFContributions[0])
   { delete PFContributions[0];
     delete PFContributions[1];
   }
  delete PFMatrix;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct RFMeshData
 { RWGGeometry *G;
   HVector *KN;
   cdouble Omega;
   RWGPortList *PortList;
   cdouble *PortCurrents;
 } RFMeshData;

HMatrix *RFFieldsMDF(void *UserData, HMatrix *XTMatrix, const char ***pDataNames)
{
  static const char *DataNames[]=
   //{ "Ex", "Ey", "Ez", "|E|", "|H|", "Radial Poynting flux" };
   { "re AxFull", "re Ay_Full", "re Az_Full", "re Phi_Full", "re dxPhi_Full", "re dyPhi_Full", "re dzPhi_Full",
     "re AxBF", "re Ay_BF", "re Az_BF", "re Phi_BF", "re dxPhi_BF", "re dyPhi_BF", "re dzPhi_BF",
     "re AxPort", "re Ay_Port", "re Az_Port", "re Phi_Port", "re dxPhi_Port", "re dyPhi_Port", "re dzPhi_Port",
     "abs AxFull", "abs Ay_Full", "abs Az_Full", "abs Phi_Full", "abs dxPhi_Full", "abs dyPhi_Full", "abs dzPhi_Full",
     "abs AxBF", "abs Ay_BF", "abs Az_BF", "abs Phi_BF", "abs dxPhi_BF", "abs dyPhi_BF", "abs dzPhi_BF",
     "abs AxPort", "abs Ay_Port", "abs Az_Port", "abs Phi_Port", "abs dxPhi_Port", "abs dyPhi_Port", "abs dzPhi_Port"
   };
  int NumData = (sizeof(DataNames) / sizeof(DataNames[0]));
  *pDataNames = DataNames;

  RFMeshData *Data      = (RFMeshData *)UserData;
  RWGGeometry *G        = Data->G;
  HVector *KN           = Data->KN;
  cdouble Omega         = Data->Omega;
  RWGPortList *PortList = Data->PortList;
  cdouble *PortCurrents = Data->PortCurrents;

  HMatrix *XMatrix = XTMatrix->Copy();
  XMatrix->Transpose();

  int NX = XMatrix->NC;
  Log("Generating FVMesh data for %i points...",NX);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *PFContributions[2];
  PFContributions[0] = new HMatrix(NPFC, NX, LHM_COMPLEX);
  PFContributions[1] = new HMatrix(NPFC, NX, LHM_COMPLEX);
  HMatrix *PFMatrix = GetMOIFields2(G, PortList, Omega, XMatrix, KN, PortCurrents, PFContributions);

  /***************************************************************/
  /* now go through the list of evaluation points and add the    */
  /* contributions of driven ports to the fields at each point,  */
  /* then fill in the appropriate colums of the DataMatrix.      */
  /***************************************************************/
  HMatrix *DataMatrix=new HMatrix(NX, NumData, LHM_REAL);
  DataMatrix->Zero();
  for(int nx=0; nx<NX; nx++)
   { 
     double X[3];
     XMatrix->GetEntriesD("0:2", nx, X);

     for(int npfc=0; npfc<NPFC; npfc++)
      { DataMatrix->SetEntry(nx,0  + npfc, PFMatrix->GetEntryD(npfc, nx));
        DataMatrix->SetEntry(nx,7  + npfc, PFContributions[0]->GetEntryD(npfc, nx));
        DataMatrix->SetEntry(nx,14 + npfc, PFContributions[1]->GetEntryD(npfc, nx));
        DataMatrix->SetEntry(nx,21 + npfc, abs(PFMatrix->GetEntry(npfc, nx)));
        DataMatrix->SetEntry(nx,28 + npfc, abs(PFContributions[0]->GetEntry(npfc, nx)));
        DataMatrix->SetEntry(nx,35 + npfc, abs(PFContributions[1]->GetEntry(npfc, nx)));
      }
   }

  delete PFMatrix;
  delete PFContributions[0];
  delete PFContributions[1];
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

void ProcessFVMesh(RWGGeometry *G, RWGPortList *PortList, cdouble Omega,
                   char *FVMesh,  HVector *KN, cdouble *PortCurrents,
                   char *FileBase)
{
  RFMeshData MyData, *Data = &MyData;
  Data->G            = G;
  Data->KN           = KN;
  Data->Omega        = Omega;
  Data->PortList     = PortList;
  Data->PortCurrents = PortCurrents;

  char OutFileName[1000];
  snprintf(OutFileName,1000,"%s.%s",FileBase,GetFileBase(FVMesh));
/* return the NPxNP impedance matrix of the geometry (NP=number*/
  MakeMeshPlot(RFFieldsMDF, (void *)Data, FVMesh, PPOptions, OutFileName);
}

/***************************************************************/
/* Z_{pq} -= I_q * DeltaV_{p}                                  */
/*                                                             */
/* where I_q == 1                                              */
/*       DeltaV_p = V_p^+ -V_p^-                               */
/*       V_p^\pm = average value of Phi at \pm terminal of     */
/*                 port p due to unit current input to port q  */
/***************************************************************/
void AddMinusIdVTermsToZMatrix(RWGGeometry *G, RWGPortList *PortList, cdouble Omega,
                               HMatrix *KMatrix, HMatrix *ZMatrix)
{
  int NBF           = G->TotalBFs;
  int NumPorts      = PortList->Ports.size();
  int NumPortEdges  = PortList->PortEdges.size();
  HMatrix *XMatrix  = new HMatrix(3, NumPortEdges, LHM_REAL);
  for(int npe=0; npe<NumPortEdges; npe++)
   {
     RWGPortEdge *PE = PortList->PortEdges[npe];
     RWGSurface *S   = G->Surfaces[PE->ns];
     RWGEdge *E      = S->GetEdgeByIndex(PE->ne);
     XMatrix->SetEntriesD(":",npe,E->Centroid);
   };
  HMatrix *PFMatrix = new HMatrix(NPFC, NumPortEdges, LHM_COMPLEX);
  cdouble *SourcePortCurrents = new cdouble[NumPorts];
  
  for(int SourcePort=0; SourcePort<NumPorts; SourcePort++)
   { 
     // get fields at centroids of all port edges due to unit-strength current into SourcePort
     HVector SourceKVector(NBF, LHM_COMPLEX, KMatrix->GetColumnPointer(SourcePort));
     memset(SourcePortCurrents, 0, NumPorts*sizeof(cdouble));
     SourcePortCurrents[SourcePort]=1.0;
     GetMOIFields(G, PortList, Omega, XMatrix, &SourceKVector, SourcePortCurrents, PFMatrix);
     
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
  delete[] SourcePortCurrents;
}

/***************************************************************/
/* return the NPxNP impedance matrix (NP=number of ports).     */
/*                                                             */
/* M = LU-factorized SIE system matrix                         */
/*                                                             */
/* If ZMatrix is non-null, it should be a caller-allocated     */
/* complex-valued HMatrix of dimension NumPorts x NumPorts.    */
/***************************************************************/
HMatrix *GetZMatrix(RWGGeometry *G, RWGPortList *PortList, cdouble Omega,
                    HMatrix *M, HMatrix *ZMatrix, HMatrix **pZTerms)
{
  int NBF      = G->TotalBFs;
  int NumPorts = PortList->Ports.size();
  
  // sanity check sizes of input arguments
  ZMatrix=CheckHMatrix(ZMatrix, NumPorts, NumPorts, LHM_COMPLEX, "GetZMatrix");

  // int WorkSize = (2*NBF + NumPorts)*NumPorts;
  // HMatrix *Workspace  = CheckHMatrix(pWorkspace ? *pWorkspace : 0, WorkSize, 1, LHM_COMPLEX, "GetZMatrix2");
  // if (pWorkspace) *pWorkspace=Workspace;
  // cdouble *WorkBuffer = Workspace->ZM;
  // HMatrix   RMatrix(NBF,       NumPorts, WorkBuffer + 0*NBF*NumPorts );
  // HMatrix   KMatrix(NBF,       NumPorts, WorkBuffer + 1*NBF*NumPorts );
  // HMatrix PPIMatrix(NumPorts,  NumPorts, WorkBuffer + 2*NBF*NumPorts );

  HMatrix   RMatrix(NBF,       NumPorts, LHM_COMPLEX);
  HMatrix   KMatrix(NBF,       NumPorts, LHM_COMPLEX);
  HMatrix PPIMatrix(NumPorts,  NumPorts, LHM_COMPLEX);

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

  Log("Computing port/BF interaction matrix...");
  GetPortBFInteractionMatrix(G, PortList, Omega, &RMatrix, &PPIMatrix);
  KMatrix.Copy(&RMatrix);
  M->LUSolve(&KMatrix);
  KMatrix.Multiply(&RMatrix, ZTerms[0], "--transA T");
  ZTerms[0]->Scale(-1.0*ZVAC);

  ZTerms[1]->Copy(&PPIMatrix);
  ZTerms[1]->Scale(1.0*ZVAC);

  ZTerms[2]->Zero();
  AddMinusIdVTermsToZMatrix(G, PortList, Omega, &KMatrix, ZTerms[2]);

  ZMatrix->Zero();
  ZMatrix->Add(ZTerms[0]);
  ZMatrix->Add(ZTerms[1]);
  ZMatrix->Add(ZTerms[2]);

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
/* Operate in-place if the second argument is NULL.            */
/***************************************************************/
void Z2S(HMatrix *Z, HMatrix *S, double ZCharacteristic)
{
  HMatrix *ZpZC=new HMatrix(Z);
  HMatrix *ZmZC=new HMatrix(Z);
  for(int nr=0; nr<ZpZC->NR; nr++)
   { ZpZC->AddEntry(nr, nr, ZCharacteristic);
     ZmZC->AddEntry(nr, nr, -ZCharacteristic);
   }
  ZpZC->LUFactorize();
  ZpZC->LUInvert();
  ZmZC->Multiply(ZpZC, S ? S : Z);

  delete ZpZC;
  delete ZmZC;
}

void S2Z(HMatrix *S, HMatrix *Z, double ZCharacteristic)
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
  if (Z==0) Z=S;
  OmS->Multiply(OpS, Z ? Z : S);
  Z->Scale(ZCharacteristic);

  delete OpS;
  delete OmS;

}

} // namespace scuff
