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

static const char *DataNamesWithPSDs[]=
   { "#RS_KdotE",
     "#IS_KdotE",
     "#RS_KdotiwA",
     "#IS_KdotiwA",
     "#RS_SigmaPhi",
     "#IS_SigmaPhi",
     "#RV_re K      (total)", 0, 0,
     "#IV_im K      (total)", 0, 0,
     "#MV_|K|       (total)", 0, 0,
     "#RS_re Sigma  (total)",
     "#IS_im Sigma  (total)",
     "#IS_|Sigma|   (total)",
     "#RV_re E      (total)", 0, 0,
     "#IV_im E      (total)", 0, 0,
     "#MV_|E|       (total)", 0, 0,
     "#RV_re iwA    (total)", 0, 0,
     "#IV_im iwA    (total)", 0, 0,
     "#MV_|iwA|     (total)", 0, 0,
     "#RS_re Phi    (total)",
     "#IS_im Phi    (total)",
     "#IS_|Phi|     (total)",
     "#RV_re -dPhi  (total)", 0, 0,
     "#IV_im -dPhi  (total)", 0, 0,
     "#MV_|dPhi|    (total)", 0, 0,
//
     "#RV_re K      (BF   )", 0, 0,
     "#IV_im K      (BF   )", 0, 0,
     "#MV_|K|       (BF   )", 0, 0,
     "#RS_re Sigma  (BF   )",
     "#IS_im Sigma  (BF   )",
     "#IS_|Sigma|   (BF   )",
     "#RV_re E      (BF   )", 0, 0,
     "#IV_im E      (BF   )", 0, 0,
     "#MV_|E|       (BF   )", 0, 0,
     "#RV_re iwA    (BF   )", 0, 0,
     "#IV_im iwA    (BF   )", 0, 0,
     "#MV_|iwA|     (BF   )", 0, 0,
     "#RS_re Phi    (BF   )",
     "#IS_im Phi    (BF   )",
     "#IS_|Phi|     (BF   )",
     "#RV_re -dPhi  (BF   )", 0, 0,
     "#IV_im -dPhi  (BF   )", 0, 0,
     "#MV_|dPhi|    (BF   )", 0, 0,
//
     "#RV_re K      (Port )", 0, 0,
     "#IV_im K      (Port )", 0, 0,
     "#MV_|K|       (Port )", 0, 0,
     "#RS_re Sigma  (Port )",
     "#IS_im Sigma  (Port )",
     "#IS_|Sigma|   (Port )",
     "#RV_re E      (Port )", 0, 0,
     "#IV_im E      (Port )", 0, 0,
     "#MV_|E|       (Port )", 0, 0,
     "#RV_re iwA    (Port )", 0, 0,
     "#IV_im iwA    (Port )", 0, 0,
     "#MV_|iwA|     (Port )", 0, 0,
     "#RS_re Phi    (Port )",
     "#IS_im Phi    (Port )",
     "#IS_|Phi|     (Port )",
     "#RV_re -dPhi  (Port )", 0, 0,
     "#IV_im -dPhi  (Port )", 0, 0,
     "#MV_|dPhi|    (Port )", 0, 0
   };
#define NUMDATA_WITHPSDS (sizeof(DataNamesWithPSDs)/sizeof(DataNamesWithPSDs[0]))

static const char *DataNamesWithoutPSDs[]=
   { "#RV_re E      (total)", 0, 0,
     "#IV_im E      (total)", 0, 0,
     "#MV_|E|       (total)", 0, 0,
     "#RV_re iwA    (total)", 0, 0,
     "#IV_im iwA    (total)", 0, 0,
     "#MV_|iwA|     (total)", 0, 0,
     "#RS_re Phi    (total)",
     "#IS_im Phi    (total)",
     "#IS_|Phi|     (total)",
     "#RV_re -dPhi  (total)", 0, 0,
     "#IV_im -dPhi  (total)", 0, 0,
     "#MV_|dPhi|    (total)", 0, 0,
//
     "#RV_re E      (BF   )", 0, 0,
     "#IV_im E      (BF   )", 0, 0,
     "#MV_|E|       (BF   )", 0, 0,
     "#RV_re iwA    (BF   )", 0, 0,
     "#IV_im iwA    (BF   )", 0, 0,
     "#MV_|iwA|     (BF   )", 0, 0,
     "#RS_re Phi    (BF   )",
     "#IS_im Phi    (BF   )",
     "#IS_|Phi|     (BF   )",
     "#RV_re -dPhi  (BF   )", 0, 0,
     "#IV_im -dPhi  (BF   )", 0, 0,
     "#MV_|dPhi|    (BF   )", 0, 0,
//
     "#RV_re E      (Port )", 0, 0,
     "#IV_im E      (Port )", 0, 0,
     "#MV_|E|       (Port )", 0, 0,
     "#RV_re iwA    (Port )", 0, 0,
     "#IV_im iwA    (Port )", 0, 0,
     "#MV_|iwA|     (Port )", 0, 0,
     "#RS_re Phi    (Port )",
     "#IS_im Phi    (Port )",
     "#IS_|Phi|     (Port )",
     "#RV_re -dPhi  (Port )", 0, 0,
     "#IV_im -dPhi  (Port )", 0, 0,
     "#MV_|dPhi|    (Port )", 0, 0
   };
#define NUMDATA_WITHOUTPSDS (sizeof(DataNamesWithoutPSDs)/sizeof(DataNamesWithoutPSDs[0]))
  
HMatrix *RFFieldsMDF(void *UserData, HMatrix *XMatrix, const char ***pDataNames)
{
  RFMeshData *Data      = (RFMeshData *)UserData;
  RWGGeometry *G        = Data->G;
  HVector *KNVector     = Data->KN;
  cdouble Omega         = Data->Omega;
  RWGPortList *PortList = Data->PortList;
  cdouble *PortCurrents = Data->PortCurrents;

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
  *pDataNames = GetPSDs ? DataNamesWithPSDs : DataNamesWithoutPSDs;

  Log("Computing FVMesh data: %i quantities at %i points",NumData,NX);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *PFContributions[2];
  PFContributions[0] = new HMatrix(NPFC, NX, LHM_COMPLEX);
  PFContributions[1] = new HMatrix(NPFC, NX, LHM_COMPLEX);
  HMatrix *PFMatrix = GetMOIFields2(G, PortList, Omega, XMatrix, KNVector, PortCurrents, PFContributions);

  /***************************************************************/
  /* Q[0,1,2][nd] = {total, BF contribution, Port contribution}  */
  /*                to data quantity #nd                         */
  /***************************************************************/
  HMatrix *DataMatrix=new HMatrix(NX, NumData, LHM_COMPLEX);
  DataMatrix->Zero();
  cdouble iw=II*Omega;
#define NUMTERMS 3
  for(int nx=0; nx<NX; nx++)
   { 
     cdouble K[NUMTERMS][3], Sigma[3], E[NUMTERMS][3], iwA[NUMTERMS][3], Phi[3], mdPhi[NUMTERMS][3];
     cdouble KdotE, KdotiwA, SigmaPhi;

     for(int Term=1; Term<=2; Term++)
      { iwA[Term][0]   =   iw*PFContributions[Term-1]->GetEntry(_PF_AX,nx);
        iwA[Term][1]   =   iw*PFContributions[Term-1]->GetEntry(_PF_AY,nx);
        iwA[Term][2]   =   iw*PFContributions[Term-1]->GetEntry(_PF_AZ,nx);
        Phi[Term]      =      PFContributions[Term-1]->GetEntry(_PF_PHI,nx);
        mdPhi[Term][0] = -1.0*PFContributions[Term-1]->GetEntry(_PF_DXPHI,nx);
        mdPhi[Term][1] = -1.0*PFContributions[Term-1]->GetEntry(_PF_DYPHI,nx);
        mdPhi[Term][2] = -1.0*PFContributions[Term-1]->GetEntry(_PF_DZPHI,nx);
        E[Term][0]     = iwA[Term][0] + mdPhi[Term][0];
        E[Term][1]     = iwA[Term][1] + mdPhi[Term][1];
        E[Term][2]     = iwA[Term][2] + mdPhi[Term][2];
      }
     VecAdd(E[1], E[2], E[0]);
     VecAdd(iwA[1], iwA[2], iwA[0]);
     Phi[0] = Phi[1] + Phi[2];
     VecAdd(mdPhi[1], mdPhi[2], mdPhi[0]);

     if (GetPSDs)
      { 
        double *X = (double *)XMatrix->GetColumnPointer(nx);
        cdouble KN[6], iwSigmaTau[2];
        EvalSourceDistribution(G, PortList, X, KNVector,     0, KN, iwSigmaTau);
        K[1][0]  = KN[0];
        K[1][1]  = KN[1];
        K[1][2]  = KN[2];
        Sigma[1] = iwSigmaTau[0] / iw;
        EvalSourceDistribution(G, PortList, X, 0, PortCurrents, KN, iwSigmaTau);
        K[2][0]  = KN[0];
        K[2][1]  = KN[1];
        K[2][2]  = KN[2];
        Sigma[2] = iwSigmaTau[0] / iw;

        VecAdd(K[1], K[2], K[0]);
        Sigma[0] = Sigma[1] + Sigma[2];
        
        KdotE    = (K[0][0]*E[0][0] + K[0][1]*E[0][1] + K[0][2]*E[0][2]);
        KdotiwA  = (K[0][0]*iwA[0][0] + K[0][1]*iwA[0][1] + K[0][2]*iwA[0][2]);
        SigmaPhi =  Sigma[0] * Phi[0];
      }
     
     for(int Term=0, nc=0; Term<NUMTERMS; Term++)
       { if (GetPSDs)
          { if (Term==0)
{
            DataMatrix->SetEntry(nx,nc++,KdotE);
            DataMatrix->SetEntry(nx,nc++,KdotE);
            DataMatrix->SetEntry(nx,nc++,KdotiwA);
            DataMatrix->SetEntry(nx,nc++,KdotiwA);
            DataMatrix->SetEntry(nx,nc++,SigmaPhi);
            DataMatrix->SetEntry(nx,nc++,SigmaPhi);
}

            DataMatrix->SetEntry(nx,nc++,K[Term][0]);
            DataMatrix->SetEntry(nx,nc++,K[Term][1]);
            DataMatrix->SetEntry(nx,nc++,K[Term][2]);
            DataMatrix->SetEntry(nx,nc++,K[Term][0]);
            DataMatrix->SetEntry(nx,nc++,K[Term][1]);
            DataMatrix->SetEntry(nx,nc++,K[Term][2]);
            DataMatrix->SetEntry(nx,nc++,K[Term][0]);
            DataMatrix->SetEntry(nx,nc++,K[Term][1]);
            DataMatrix->SetEntry(nx,nc++,K[Term][2]);
            DataMatrix->SetEntry(nx,nc++,Sigma[Term]);
            DataMatrix->SetEntry(nx,nc++,Sigma[Term]);
            DataMatrix->SetEntry(nx,nc++,Sigma[Term]);
          }
         DataMatrix->SetEntry(nx,nc++,E[Term][0]);
         DataMatrix->SetEntry(nx,nc++,E[Term][1]);
         DataMatrix->SetEntry(nx,nc++,E[Term][2]);
         DataMatrix->SetEntry(nx,nc++,E[Term][0]);
         DataMatrix->SetEntry(nx,nc++,E[Term][1]);
         DataMatrix->SetEntry(nx,nc++,E[Term][2]);
         DataMatrix->SetEntry(nx,nc++,E[Term][0]);
         DataMatrix->SetEntry(nx,nc++,E[Term][1]);
         DataMatrix->SetEntry(nx,nc++,E[Term][2]);

         DataMatrix->SetEntry(nx,nc++,iwA[Term][0]);
         DataMatrix->SetEntry(nx,nc++,iwA[Term][1]);
         DataMatrix->SetEntry(nx,nc++,iwA[Term][2]);
         DataMatrix->SetEntry(nx,nc++,iwA[Term][0]);
         DataMatrix->SetEntry(nx,nc++,iwA[Term][1]);
         DataMatrix->SetEntry(nx,nc++,iwA[Term][2]);
         DataMatrix->SetEntry(nx,nc++,iwA[Term][0]);
         DataMatrix->SetEntry(nx,nc++,iwA[Term][1]);
         DataMatrix->SetEntry(nx,nc++,iwA[Term][2]);

         DataMatrix->SetEntry(nx,nc++,Phi[Term]);
         DataMatrix->SetEntry(nx,nc++,Phi[Term]);
         DataMatrix->SetEntry(nx,nc++,Phi[Term]);

         DataMatrix->SetEntry(nx,nc++,mdPhi[Term][0]);
         DataMatrix->SetEntry(nx,nc++,mdPhi[Term][1]);
         DataMatrix->SetEntry(nx,nc++,mdPhi[Term][2]);
         DataMatrix->SetEntry(nx,nc++,mdPhi[Term][0]);
         DataMatrix->SetEntry(nx,nc++,mdPhi[Term][1]);
         DataMatrix->SetEntry(nx,nc++,mdPhi[Term][2]);
         DataMatrix->SetEntry(nx,nc++,mdPhi[Term][0]);
         DataMatrix->SetEntry(nx,nc++,mdPhi[Term][1]);
         DataMatrix->SetEntry(nx,nc++,mdPhi[Term][2]);
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

HMatrix *ProcessFVMesh(RWGGeometry *G, RWGPortList *PortList, cdouble Omega,
                       const char *FVMesh, const char *TransFile,
                       HVector *KN, cdouble *PortCurrents, char *FileBase)
{
  RFMeshData MyData, *Data = &MyData;
  Data->G            = G;
  Data->KN           = KN;
  Data->Omega        = Omega;
  Data->PortList     = PortList;
  Data->PortCurrents = PortCurrents;

  HMatrix *Integrals = MakeMeshPlot(RFFieldsMDF,(void *)Data, FVMesh, TransFile, FileBase, PPOptions);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
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
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  return Integrals;
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
