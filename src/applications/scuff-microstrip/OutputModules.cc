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
 * ProcessEPFile.cc -- scuff-rf code file for evaluating radiated fields 
 *                  -- at user-specified evaluation points 
 * 
 * homer reid       -- 9/2011
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

#include "RWGPorts.h"

#define FREQ2OMEGA (2.0*M_PI/300.0)
#define OMEGA2FREQ (1/(FREQ2OMEGA))
#define II cdouble(0.0,1.0)

using namespace scuff;

// FIXME put me somewhere else
#define NPFC 7
#define _PF_PHI 3

HMatrix *GetMOIFields(RWGGeometry *G, LayeredSubstrate *S, cdouble Omega, HMatrix *XMatrix,
                      HVector *KN, RWGPortList *PortList, cdouble *PortCurrents,
                      HMatrix **PFContributions=0);

void GetMOIRPFMatrices(RWGGeometry *G, LayeredSubstrate *Substrate, RWGPortList *PortList,
                       cdouble Omega, HMatrix *XMatrix,
                       HMatrix **pBFRPFMatrix, HMatrix **pPortRPFMatrix);

void ZToS(HMatrix *ZMatrix, HMatrix *SMatrix);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitZSParmFile(char *FileBase, int NumPorts, char ZS)
{
  FILE *f=vfopen("%s.%cparms","w",FileBase,ZS);
  setbuf(f,0);

  char TimeString[200];
  time_t MyTime=time(0);
  struct tm *MyTm=localtime(&MyTime);
  strftime(TimeString,30,"%D::%T",MyTm);
  fprintf(f,"# scuff-microstrip ran on %s (%s)\n",getenv("HOST"),TimeString);
  fprintf(f,"# columns:\n");
  fprintf(f,"# 1 frequency (GHz)\n");
  for(int ndp=0, nc=2; ndp<NumPorts; ndp++)
   for(int nsp=0; nsp<NumPorts; nsp++, nc+=2)
    fprintf(f,"#%i,%i real,imag %c_{%i%i}\n",nc,nc+1,ZS,ndp,nsp);
  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeSZParms(RWGGeometry *G, RWGPortList *PortList, cdouble Omega,
                    HMatrix *M, HVector *KN, char *FileBase, bool SParms)
{
  int NumPorts = PortList->Ports.size();
  static bool Init=true;
  static HMatrix *XMatrix, *BFRPFMatrix=0, *PortRPFMatrix=0;
  static HVector *BFPFVector, *PortPFVector;
  if (Init)
   { Init=false;

     InitZSParmFile(FileBase, NumPorts, 'Z');
     if(SParms) InitZSParmFile(FileBase, NumPorts, 'S');

     XMatrix = new HMatrix(3, 2*NumPorts, LHM_REAL);
     for(int nPort=0; nPort<PortList->Ports.size(); nPort++)
      { RWGPort *Port = PortList->Ports[nPort];
        XMatrix->SetEntriesD(":",2*nPort + _PLUS,  Port->RefPoint[_PLUS]);
        XMatrix->SetEntriesD(":",2*nPort + _MINUS, Port->RefPoint[_MINUS]);
      }

     int NX = XMatrix->NC;
     BFPFVector   = new HVector(NPFC*NX, LHM_COMPLEX);
     PortPFVector = new HVector(NPFC*NX, LHM_COMPLEX);
   }
  cdouble *PortCurrents = new cdouble[NumPorts];
  HMatrix *ZMatrix      = new HMatrix(NumPorts, NumPorts, LHM_COMPLEX);
  HMatrix *SMatrix      = SParms ? new HMatrix(NumPorts, NumPorts, LHM_COMPLEX) : 0;

  Log("Precomputing reduced potential/field matrices for all port reference points");
  GetMOIRPFMatrices(G, G->Substrate, PortList, Omega, XMatrix, &BFRPFMatrix, &PortRPFMatrix);

  for(int nsPort=0; nsPort<NumPorts; nsPort++)
   { 
     Log(" Computing row %i of the impedance matrix:",nsPort+1);
     memset(PortCurrents, 0, NumPorts*sizeof(cdouble));
     PortCurrents[nsPort]=1.0;
     KN->Zero();
     AddPortContributionsToRHS(G, PortList, PortCurrents, Omega, KN);

     Log("  solving the BEM system");
     M->LUSolve(KN);

     /* note: the entry in the Z-matrix is the complex conjugate */
     /* of the measured port voltage, because the Z-matrix is    */
     /* defined using the usual circuit theory convention in     */
     /* which all quantities have time dependence exp(+iwt),     */
     /* whereas scuff-EM uses the opposite sign convention.      */
     Log("  computing port voltages");
     BFRPFMatrix->Apply(KN, BFPFVector);
     HVector PCVector(NumPorts, PortCurrents);
     PortRPFMatrix->Apply(&PCVector, PortPFVector);
     for(int ndPort=0; ndPort<NumPorts; ndPort++)
      { int nxP = 2*ndPort+0, nxM = 2*ndPort+1;
        cdouble VP = BFPFVector->GetEntry(NPFC*nxP + _PF_PHI) + PortPFVector->GetEntry(NPFC*nxP + _PF_PHI);
        cdouble VM = BFPFVector->GetEntry(NPFC*nxM + _PF_PHI) + PortPFVector->GetEntry(NPFC*nxM + _PF_PHI);
        ZMatrix->SetEntry(ndPort, nsPort, conj(VP-VM) );
      }
   } // for(int nsp=0; nsp<NumPorts; nsp++)

  if (SParms) ZToS(ZMatrix, SMatrix);

  /*--------------------------------------------------------------*/
  /*- write output files                                          */
  /*--------------------------------------------------------------*/
  for(int nZS=0; nZS<(SParms ? 2 : 1); nZS++)
   { char ZS = (nZS==0 ? 'Z' : 'S');
     HMatrix *ZSMatrix = (nZS==0 ? ZMatrix : SMatrix);
     FILE *f=vfopen("%s.%cparms","a",FileBase,ZS);
     fprintf(f,"%e ",real(Omega)*OMEGA2FREQ);
     for(int ndPort=0; ndPort<NumPorts; ndPort++)
      for(int nsPort=0; nsPort<NumPorts; nsPort++)
       fprintf(f,"%e %e ",real(ZSMatrix->GetEntry(ndPort,nsPort)), 
                          imag(ZSMatrix->GetEntry(ndPort,nsPort)));
     fprintf(f,"\n");
     fclose(f);
   }

  delete[] PortCurrents;
  delete ZMatrix;
  if (SMatrix) delete SMatrix;
}
  
/***************************************************************/
/***************************************************************/
/***************************************************************/
void ProcessEPFile(RWGGeometry *G, HVector *KN, cdouble Omega,
                   RWGPortList *PortList, cdouble *PortCurrents,
                   char *EPFile, char *FileBase)
{
  /***************************************************************/
  /* read in the list of evaluation points  **********************/
  /***************************************************************/ 
  HMatrix *XMatrix=new HMatrix(EPFile);
  XMatrix->Transpose();
  int NX = XMatrix->NC;

  HMatrix *PFContributions[2];
  PFContributions[0] = new HMatrix(NPFC, NX);
  PFContributions[1] = new HMatrix(NPFC, NX);
  HMatrix *PFMatrix = GetMOIFields(G, G->Substrate, Omega, XMatrix, KN, PortList, PortCurrents, PFContributions);

  /***************************************************************/
  /* write field components to output file ***********************/
  /***************************************************************/
  FILE *f = vfopen("%s.%s.PF",FileBase,EPFile);
  for(int nx=0; nx<NX; nx++)
   { fprintVec(f,(double *)XMatrix->GetColumnPointer(nx));
     fprintVec(f,(cdouble *)PFContributions[0]->GetColumnPointer(nx),NPFC);
     fprintVecCR(f,(cdouble *)PFContributions[1]->GetColumnPointer(nx),NPFC);
   };
  fclose(f);

  delete PFContributions[0];
  delete PFContributions[1];
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
  PFContributions[0] = new HMatrix(NPFC, NX);
  PFContributions[1] = new HMatrix(NPFC, NX);
  HMatrix *PFMatrix = GetMOIFields(G, G->Substrate, Omega, XMatrix, KN, PortList, PortCurrents, PFContributions);

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

     for(int npfc=0, nd=0; npfc<NPFC; npfc++)
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

void ProcessFVMesh(RWGGeometry *G, HVector *KN, cdouble Omega,
                   RWGPortList *PortList, cdouble *PortCurrents,
                   char *FVMesh, char *FileBase)
{
  RFMeshData MyData, *Data = &MyData;
  Data->G            = G;
  Data->KN           = KN;
  Data->Omega        = Omega;
  Data->PortList     = PortList;
  Data->PortCurrents = PortCurrents;

  char OutFileName[1000];
  snprintf(OutFileName,1000,"%s.%s",FileBase,GetFileBase(FVMesh));
  MakeMeshPlot(RFFieldsMDF, (void *)Data, FVMesh, PPOptions, OutFileName);
}
