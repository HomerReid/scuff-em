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

#include <stdio.h>
#include <stdlib.h>

#include "libscuff.h"
#include "libhmat.h"
#include "CylindricalWave.h"

using namespace scuff;
#define II cdouble (0.0,1.0)

void GetCylMN(cdouble k0, int Nu, double kz, int WaveType,
              double RPZ[3], cdouble MVec[3], cdouble NVec[3]);

static char *FieldFuncs=const_cast<char *>(
 "|Ex|,|Ey|,|Ez|,"
 "sqrt(|Ex|^2+|Ey|^2+|Ez|^2),"
 "|Hx|,|Hy|,|Hz|,"
 "sqrt(|Hx|^2+|Hy|^2+|Hz|^2)");

static const char *FieldTitles[]=
 {"|Ex|", "|Ey|", "|Ez|", "|E|",
  "|Hx|", "|Hy|", "|Hz|", "|H|",
 };

#define NUMFIELDFUNCS 8

/***************************************************************/
/***************************************************************/
/***************************************************************/
void VisualizeFields(RWGGeometry *G, cdouble Omega, double *kBloch,
                     HVector *KN, IncField *IF, char *MeshFileName)
{ 
  /*--------------------------------------------------------------*/
  /*- try to open output file ------------------------------------*/
  /*--------------------------------------------------------------*/
  char GeoFileBase[100], PPFileName[100];
  strncpy(GeoFileBase,GetFileBase(G->GeoFileName),100);
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
  HMatrix *FMatrix=G->GetFields(IF, KN, Omega, kBloch, 
                                XMatrix, 0, FieldFuncs);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nff=0; nff<NUMFIELDFUNCS; nff++)
   { 
     fprintf(f,"View \"%s(%s)\" {\n",FieldTitles[nff],z2s(Omega));

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     for(int np=0; np<S->NumPanels; np++)
      {
        RWGPanel *P=S->Panels[np];
        int iV1 = P->VI[0];  double *V1 = S->Vertices + 3*iV1;
        int iV2 = P->VI[1];  double *V2 = S->Vertices + 3*iV2;
        int iV3 = P->VI[2];  double *V3 = S->Vertices + 3*iV3;

        fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                   V1[0], V1[1], V1[2],
                   V2[0], V2[1], V2[2],
                   V3[0], V3[1], V3[2],
                   FMatrix->GetEntryD(iV1,nff),
                   FMatrix->GetEntryD(iV2,nff),
                   FMatrix->GetEntryD(iV3,nff));
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
/***************************************************************/
/***************************************************************/
void GetScatteringCoefficients(cdouble k0, int Nu, double kz, int Pol,
                               double R, cdouble Epsilon, cdouble ABCD[4])
{
  double RPZ[3]={0.0, 0.0, 0.0};
  RPZ[0]=R;

  cdouble nn=sqrt(Epsilon);

  cdouble MReg[3], NReg[3];
  cdouble MOut[3], NOut[3];
  cdouble MInc[3], NInc[3];
  GetCylMN(nn*k0, Nu, kz, CW_REGULAR,  RPZ, MReg, NReg);
  GetCylMN(   k0, Nu, kz, CW_OUTGOING, RPZ, MOut, NOut);
  GetCylMN(   k0, Nu, kz, CW_INCOMING, RPZ, MInc, NInc);
 
  HMatrix MMatrix(4,4,LHM_COMPLEX);

  MMatrix.SetEntry(0,0,  MReg[1] );
  MMatrix.SetEntry(0,1,  NReg[1] );
  MMatrix.SetEntry(0,2, -MOut[1] );
  MMatrix.SetEntry(0,3, -NOut[1] );

  MMatrix.SetEntry(1,0,  nn*NReg[1] );
  MMatrix.SetEntry(1,1, -nn*MReg[1] );
  MMatrix.SetEntry(1,2,    -NOut[1] );
  MMatrix.SetEntry(1,3,    +MOut[1] );

  MMatrix.SetEntry(2,0,  MReg[2] );
  MMatrix.SetEntry(2,1,  NReg[2] );
  MMatrix.SetEntry(2,2, -MOut[2] );
  MMatrix.SetEntry(2,3, -NOut[2] );

  MMatrix.SetEntry(3,0,  nn*NReg[2] );
  MMatrix.SetEntry(3,1, -nn*MReg[2] );
  MMatrix.SetEntry(3,2,    -NOut[2] );
  MMatrix.SetEntry(3,3,    +MOut[2] );

  double P=0.0, Q=0.0;
  if (Pol==CW_TE2Z)
   P=1.0;
  else
   Q=1.0;

  HVector RHS(4,LHM_COMPLEX);
  RHS.SetEntry(0, P*MInc[1] + Q*NInc[1] );
  RHS.SetEntry(1, P*NInc[1] - Q*MInc[1] );
  RHS.SetEntry(2, P*MInc[2] + Q*NInc[2] );
  RHS.SetEntry(3, P*NInc[2] - Q*MInc[2] );

  MMatrix.LUFactorize();
  MMatrix.LUSolve(&RHS);
  ABCD[0]=RHS.ZV[0];
  ABCD[1]=RHS.ZV[1];
  ABCD[2]=RHS.ZV[2];
  ABCD[3]=RHS.ZV[3];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetExactFields(cdouble Omega, int Nu, double kz, int Pol,
                    double R, cdouble EpsRel, cdouble ABCD[4],
                    HMatrix *XMatrix, HMatrix *FScat, HMatrix *FTot)
{ 
  double P = (Pol==CW_TE2Z ? 1.0 : 0.0);
  double Q = (Pol==CW_TE2Z ? 0.0 : 1.0);
  cdouble A=ABCD[0];
  cdouble B=ABCD[1];
  cdouble C=ABCD[2];
  cdouble D=ABCD[3];

  CylindricalWave CW_MOut(Nu, kz, CW_TE2Z, CW_OUTGOING);
  CylindricalWave CW_NOut(Nu, kz, CW_TM2Z, CW_OUTGOING);
  CylindricalWave CW_MInc(Nu, kz, CW_TE2Z, CW_INCOMING);
  CylindricalWave CW_NInc(Nu, kz, CW_TM2Z, CW_INCOMING);
  CylindricalWave CW_MReg(Nu, kz, CW_TE2Z, CW_REGULAR);
  CylindricalWave CW_NReg(Nu, kz, CW_TM2Z, CW_REGULAR);

  CW_MOut.SetFrequency(Omega);
  CW_NOut.SetFrequency(Omega);
  CW_MInc.SetFrequency(Omega);
  CW_NInc.SetFrequency(Omega);
  CW_MReg.SetFrequencyAndEpsMu(Omega,EpsRel,1.0);
  CW_NReg.SetFrequencyAndEpsMu(Omega,EpsRel,1.0);

  double R2=R*R;
  for (int nr=0; nr<XMatrix->NR; nr++)
   { 
     double XYZ[3];
     XMatrix->GetEntriesD(nr,":",XYZ);

     cdouble EHTot[6], EHScat[6];
     if ( (XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2]) >= R2 )
      {  
        cdouble EHMOut[6], EHNOut[6], EHMInc[6], EHNInc[6];
        CW_MOut.GetFields(XYZ, EHMOut);
        CW_NOut.GetFields(XYZ, EHNOut);
        CW_MInc.GetFields(XYZ, EHMInc);
        CW_NInc.GetFields(XYZ, EHNInc);
        for(int Mu=0; Mu<6; Mu++)
         { EHScat[Mu] = C*EHMOut[Mu] + D*EHNOut[Mu];
           EHTot[Mu]  = EHScat[Mu] + P*EHMInc[Mu] + Q*EHNInc[Mu];
         };
      }
     else
      {  
        cdouble EHMReg[6], EHNReg[6];
        CW_MReg.GetFields(XYZ, EHMReg);
        CW_NReg.GetFields(XYZ, EHNReg);
        for(int Mu=0; Mu<6; Mu++)
         EHScat[Mu] = EHTot[Mu] = A*EHMReg[Mu] + B*EHNReg[Mu];
      };

     if (FScat) FScat->SetEntries(nr,":",EHScat);
     if (FTot)  FTot->SetEntries(nr,":",EHTot);
   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetExactSurfaceCurrents(double XYZ[3],
                             cdouble Omega, int Nu, double kz, int Pol,
                             cdouble EpsRel, cdouble ABCD[4],
                             cdouble KExact[3], cdouble NExact[3])
{
  cdouble A=ABCD[0];
  cdouble B=ABCD[1];

  CylindricalWave CW_MReg(Nu, kz, CW_TE2Z, CW_REGULAR);
  CylindricalWave CW_NReg(Nu, kz, CW_TM2Z, CW_REGULAR);

  CW_MReg.SetFrequencyAndEpsMu(Omega,EpsRel,1.0);
  CW_NReg.SetFrequencyAndEpsMu(Omega,EpsRel,1.0);

  cdouble EHMReg[6], EHNReg[6], EH[6];
  CW_MReg.GetFields(XYZ, EHMReg);
  CW_NReg.GetFields(XYZ, EHNReg);
  for(int Mu=0; Mu<6; Mu++)
   EH[Mu] = A*EHMReg[Mu] + B*EHNReg[Mu];

  double RhoHat[3];
  RhoHat[0] = 0.0;
  RhoHat[1] = XYZ[1];
  RhoHat[2] = XYZ[2];
  VecNormalize(RhoHat);

  for(int Mu=0; Mu<3; Mu++)
   { int Nu=(Mu+1)%3, Rho=(Mu+2)%3;
     KExact[Mu] = -(RhoHat[Nu]*EH[3+Rho] - RhoHat[Rho]*EH[3+Nu]);
     NExact[Mu] =  (RhoHat[Nu]*EH[Rho]   - RhoHat[Rho]*EH[Nu]);
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PlotExactSurfaceCurrents(RWGGeometry *G,
                              cdouble Omega, int Nu,
                              double kz, int Pol,
                              cdouble EpsRel, cdouble ABCD[4])
{
  FILE *f=fopen("ExactSurfaceCurrents.pp","w");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","Exact Electric Current");
  for(int ns=0; ns<G->NumSurfaces; ns++)
   for(int np=0; np<G->Surfaces[ns]->NumPanels; np++)
    { 
      double *XYZ=G->Surfaces[ns]->Panels[np]->Centroid;
      cdouble KExact[3], NExact[3];
      GetExactSurfaceCurrents(XYZ, Omega, Nu, kz, Pol,
                              EpsRel, ABCD, KExact, NExact);
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 XYZ[0],XYZ[1],XYZ[2],
                 real(KExact[0]), real(KExact[1]), real(KExact[2]));
    };
  fprintf(f,"};\n");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","Exact Magnetic Current");
  for(int ns=0; ns<G->NumSurfaces; ns++)
   for(int np=0; np<G->Surfaces[ns]->NumPanels; np++)
    { 
      double *XYZ=G->Surfaces[ns]->Panels[np]->Centroid;
      cdouble KExact[3], NExact[3];
      GetExactSurfaceCurrents(XYZ, Omega, Nu, kz, Pol,
                              EpsRel, ABCD, KExact, NExact);
      fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                 XYZ[0],XYZ[1],XYZ[2],
                 real(NExact[0]), real(NExact[1]), real(NExact[2]));
    };
  fprintf(f,"};\n");
  fclose(f);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PlotSCUFFSurfaceCurrents(HMatrix *PSDMatrix)
{
  FILE *f=fopen("SCUFFSurfaceCurrents.pp","w");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","SCUFF Electric Current");
  for(int np=0; np<PSDMatrix->NR; np++)
   { 
     double XYZ[3];
     cdouble KExact[3], NExact[3];
     PSDMatrix->GetEntriesD(np,"0:2",XYZ);
     PSDMatrix->GetEntries(np,"5:7",KExact);
     PSDMatrix->GetEntries(np,"9:11",NExact);
     fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                XYZ[0],XYZ[1],XYZ[2],
                real(KExact[0]), real(KExact[1]), real(KExact[2]));
   };
  fprintf(f,"};\n");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fprintf(f,"View \"%s\" {\n","SCUFF Magnetic Current");
  for(int np=0; np<PSDMatrix->NR; np++)
   { 
     double XYZ[3];
     cdouble KExact[3], NExact[3];
     PSDMatrix->GetEntriesD(np,"0:2",XYZ);
     PSDMatrix->GetEntries(np,"5:7",KExact);
     PSDMatrix->GetEntries(np,"9:11",NExact);
     fprintf(f,"VP(%e,%e,%e) {%e,%e,%e};\n",
                XYZ[0],XYZ[1],XYZ[2],
                real(NExact[0]), real(NExact[1]), real(NExact[2]));
   };
  fprintf(f,"};\n");

  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void CompareSurfaceCurrents(RWGGeometry *G, HVector *KN,
                            cdouble Omega, int Nu, double kz, int Pol, 
                            cdouble EpsRel, cdouble ABCD[4])
{
  static HMatrix *PSDMatrix=0;
  double kBloch[2]={0.0,0.0};
  kBloch[0]=kz;
  PSDMatrix=G->GetPanelSourceDensities(Omega, kBloch, KN, PSDMatrix);
  PlotSCUFFSurfaceCurrents(PSDMatrix);

  FILE *f=vfopen("%s.Currents","w",GetFileBase(G->GeoFileName));
  double KExactAvg=0.0, KHRAvg=0.0;
  double NExactAvg=0.0, NHRAvg=0.0;
  double TotalArea=0.0;
  for(int np=0; np<PSDMatrix->NR; np++)
   { 
     double XYZ[3], Area;
     cdouble KHR[3], NHR[3];
     PSDMatrix->GetEntriesD(np,"0:2",XYZ);
     Area=PSDMatrix->GetEntryD(np,3);
     PSDMatrix->GetEntries(np,"5:7",KHR);
     PSDMatrix->GetEntries(np,"9:11",NHR);
  
     cdouble KExact[3], NExact[3];
     GetExactSurfaceCurrents(XYZ, Omega, Nu, kz, Pol, EpsRel, ABCD, KExact, NExact);
     
     fprintf(f,"%e %e %e ",XYZ[0],XYZ[1],XYZ[2]);
     fprintf(f,"%s %s %s ",CD2S(KHR[0]),CD2S(KHR[1]),CD2S(KHR[2]));
     fprintf(f,"%s %s %s ",CD2S(NHR[0]),CD2S(NHR[1]),CD2S(NHR[2]));
     fprintf(f,"%s %s %s ",CD2S(KExact[0]),CD2S(KExact[1]),CD2S(KExact[2]));
     fprintf(f,"%s %s %s ",CD2S(NExact[0]),CD2S(NExact[1]),CD2S(NExact[2]));
     fprintf(f,"\n");

     TotalArea += Area;
     KExactAvg += Area*sqrt( norm(KExact[0]) + norm(KExact[1]) + norm(KExact[2]) );
     NExactAvg += Area*sqrt( norm(NExact[0]) + norm(NExact[1]) + norm(NExact[2]) );
     KHRAvg    += Area*sqrt( norm(KHR[0])    + norm(KHR[1])    + norm(KHR[2])    );
     NHRAvg    += Area*sqrt( norm(NHR[0])    + norm(NHR[1])    + norm(NHR[2])    );

   };
  fclose(f);
  KExactAvg/=TotalArea;
  NExactAvg/=TotalArea;
  KHRAvg/=TotalArea;
  NHRAvg/=TotalArea;

  f=fopen("tCS.stats","a");
  fprintf(f,"%e %e %e %i %i ",real(Omega),imag(Omega),kz,Nu,Pol);
  fprintf(f,"%e %e %e ",KExactAvg,KHRAvg,RD(KExactAvg,KHRAvg));
  fprintf(f,"%e %e %e ",NExactAvg,NHRAvg,RD(NExactAvg,NHRAvg));
  fprintf(f,"\n");
  fclose(f);
 
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();

  /*--------------------------------------------------------------*/
  /*- process command-line options -------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;
  double Radius=1.0;
  cdouble Omega=1.0;
  char *OmegaFile=0;
  char *EPFile=0;
  int Nu=0;
  double kz=0.0;
  bool TMPolarization=false;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",  PA_STRING,  1, 1, (void *)&GeoFileName,  0,  ".scuffgeo file"},
     {"Radius",    PA_DOUBLE,  1, 1, (void *)&Radius,       0,  "cylinder radius"},
     {"Omega",     PA_CDOUBLE, 1, 1, (void *)&Omega,        0,  "angular frequency"},
     {"OmegaFile", PA_STRING,  1, 1, (void *)&OmegaFile,    0,  "omega file"},
     {"EPFile",    PA_STRING,  1, 1, (void *)&EPFile,       0,  "list of evaluation points"},
     {"Nu",        PA_INT,     1, 1, (void *)&Nu,           0,  "nu"},
     {"kz",        PA_DOUBLE,  1, 1, (void *)&kz,           0,  "kz"},
     {"TM",        PA_BOOL,    0, 1, (void *)&TMPolarization, 0,  "use TM instead of default TE polarization"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /*--------------------------------------------------------------*/
  /*- process frequency options ----------------------------------*/
  /*--------------------------------------------------------------*/
  HVector *OmegaVector=0;
  if (OmegaFile)
   OmegaVector=new HVector(OmegaFile);
  else if ( Omega!=0.0 )
   { OmegaVector=new HVector(1, LHM_COMPLEX);
     OmegaVector->SetEntry(0,Omega);
   }
  else
   OSUsage(argv[0],OSArray,"either --omega or --omegafile must be specified");

  double kBloch[2]={0.0, 0.0};
  kBloch[0] = kz;

  /*--------------------------------------------------------------*/
  /* create the RWGGeometry from the .scuffgeo file               */
  /*--------------------------------------------------------------*/
  SetLogFileName("%s.log",GetFileBase(GeoFileName));
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);

  /*--------------------------------------------------------------*/
  /* preallocate BEM matrix and RHS vector                        */
  /*--------------------------------------------------------------*/
  HMatrix *M  = G->AllocateBEMMatrix();
  HVector *KN = G->AllocateRHSVector();

  int Pol = TMPolarization ? CW_TM2Z : CW_TE2Z;
  CylindricalWave CW(Nu, kz, Pol, CW_INCOMING);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *XMatrix1=new HMatrix(2,3);
  XMatrix1->SetEntry(0,0,0.0);
  XMatrix1->SetEntry(0,1,0.0);
  XMatrix1->SetEntry(0,2,0.0);
  XMatrix1->SetEntry(1,0,0.0);
  XMatrix1->SetEntry(1,1,2.0);
  XMatrix1->SetEntry(1,2,2.0);
  HMatrix *FMatrix1=new HMatrix(2,6,LHM_COMPLEX);

  HMatrix *XMatrix2=0, *FScat=0, *FTot=0;
  if (EPFile)
   { XMatrix2=new HMatrix(EPFile);
     if (XMatrix2->ErrMsg)
      ErrExit(XMatrix2->ErrMsg);
     FScat=new HMatrix(XMatrix2->NR,6,LHM_COMPLEX);
     FTot=new HMatrix(XMatrix2->NR,6,LHM_COMPLEX);
   };

  SetDefaultCD2SFormat("%.8e %.8e");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)
   {
     Omega=OmegaVector->GetEntry(nOmega);
     Log("Computing at frequency %s...",z2s(Omega));

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     cdouble ABCD[4];
     cdouble EpsRel, MuRel;
     G->RegionMPs[1]->GetEpsMu(Omega, &EpsRel, &MuRel);
     GetScatteringCoefficients(Omega, Nu, kz, Pol, Radius, EpsRel, ABCD);
     FILE *f0=fopen("ABCD.out","w");
     fprintf(f0,"%s %i %e %i %s ",z2s(Omega),Nu,kz,Pol,z2s(EpsRel));
     fprintf(f0,"%s %s %s %s ", CD2S(ABCD[0]),CD2S(ABCD[1]),
                                CD2S(ABCD[2]),CD2S(ABCD[3]));
     fprintf(f0,"\n");
     fclose(f0);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     if (XMatrix2)
      { 
        GetExactFields(Omega, Nu, kz, Pol, Radius, EpsRel, ABCD,
                       XMatrix2, FScat, FTot);

        FILE *fScat=vfopen("%s.exact.scattered","w",GetFileBase(GeoFileName));
        FILE *fTot=vfopen("%s.exact.total","w",GetFileBase(GeoFileName));
        for(int nr=0; nr<XMatrix2->NR; nr++)
         { 
           fprintf(fScat,"%s %e %i %i ",z2s(Omega),kz,Nu,Pol);
           fprintf(fScat,"%e %e %e ",XMatrix2->GetEntryD(nr,0),
                                     XMatrix2->GetEntryD(nr,1),
                                     XMatrix2->GetEntryD(nr,2));
           fprintf(fScat,"%s %s %s ",CD2S(FScat->GetEntry(nr,0)),
                                     CD2S(FScat->GetEntry(nr,1)),
                                     CD2S(FScat->GetEntry(nr,2)));
           fprintf(fScat,"%s %s %s ",CD2S(FScat->GetEntry(nr,3)),
                                     CD2S(FScat->GetEntry(nr,4)),
                                     CD2S(FScat->GetEntry(nr,5)));
           fprintf(fScat,"\n");

           fprintf(fTot,"%s %e %i %i ",z2s(Omega),kz,Nu,Pol);
           fprintf(fTot,"%e %e %e ",XMatrix2->GetEntryD(nr,0),
                                    XMatrix2->GetEntryD(nr,1),
                                    XMatrix2->GetEntryD(nr,2));
           fprintf(fTot,"%s %s %s ",CD2S(FTot->GetEntry(nr,0)),
                                    CD2S(FTot->GetEntry(nr,1)),
                                    CD2S(FTot->GetEntry(nr,2)));
           fprintf(fTot,"%s %s %s ",CD2S(FTot->GetEntry(nr,3)),
                                    CD2S(FTot->GetEntry(nr,4)),
                                    CD2S(FTot->GetEntry(nr,5)));
           fprintf(fTot,"\n");

         };
        fclose(fScat);
        fclose(fTot);

      }; // if (XMatrix2)

     /*--------------------------------------------------------------*/
     /* assemble and factorize the BEM matrix at this frequency      */
     /*--------------------------------------------------------------*/
     G->AssembleBEMMatrix(Omega, kBloch, M);
     M->LUFactorize();

     G->AssembleRHSVector(Omega, kBloch, &CW, KN);
     M->LUSolve(KN);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     CompareSurfaceCurrents(G, KN, Omega, Nu, kz, Pol, EpsRel, ABCD);

     G->PlotSurfaceCurrents(KN, Omega, kBloch, "%s.pp",GetFileBase(GeoFileName));
     PlotExactSurfaceCurrents(G, Omega, Nu, kz, Pol, EpsRel, ABCD);
/*
     VisualizeFields(G, Omega, kBloch, 0, &CW, "Square_1160.msh");
     VisualizeFields(G, Omega, kBloch, KN, 0, "Square_1160.msh");
*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/ 
     G->GetFields(0, KN, Omega, kBloch, XMatrix1, FMatrix1);
     FILE *f1=vfopen("%s.out","a",GetFileBase(GeoFileName));
     fprintf(f1,"%s %e %i %i ",z2s(Omega),kz,Nu,Pol);
     fprintf(f1,"%s %s %s %s %s %s ",CD2S(FMatrix1->GetEntry(0,0)),
                                     CD2S(FMatrix1->GetEntry(0,1)),
                                     CD2S(FMatrix1->GetEntry(0,2)),
                                     CD2S(FMatrix1->GetEntry(0,3)),
                                     CD2S(FMatrix1->GetEntry(0,4)),
                                     CD2S(FMatrix1->GetEntry(0,5)));
     fprintf(f1,"%s %s %s %s %s %s ",CD2S(FMatrix1->GetEntry(1,0)),
                                     CD2S(FMatrix1->GetEntry(1,1)),
                                     CD2S(FMatrix1->GetEntry(1,2)),
                                     CD2S(FMatrix1->GetEntry(1,3)),
                                     CD2S(FMatrix1->GetEntry(1,4)),
                                     CD2S(FMatrix1->GetEntry(1,5)));
     fprintf(f1,"\n");
     fclose(f1);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     GetExactFields(Omega, Nu, kz, Pol, Radius, EpsRel, ABCD,
                    XMatrix1, FMatrix1, 0);
     f1=vfopen("%s.exact","a",GetFileBase(GeoFileName));
     fprintf(f1,"%s %e %i %i ",z2s(Omega),kz,Nu,Pol);
     fprintf(f1,"%s %s %s %s %s %s ",CD2S(FMatrix1->GetEntry(0,0)),
                                     CD2S(FMatrix1->GetEntry(0,1)),
                                     CD2S(FMatrix1->GetEntry(0,2)),
                                     CD2S(FMatrix1->GetEntry(0,3)),
                                     CD2S(FMatrix1->GetEntry(0,4)),
                                     CD2S(FMatrix1->GetEntry(0,5)));
     fprintf(f1,"%s %s %s %s %s %s ",CD2S(FMatrix1->GetEntry(1,0)),
                                     CD2S(FMatrix1->GetEntry(1,1)),
                                     CD2S(FMatrix1->GetEntry(1,2)),
                                     CD2S(FMatrix1->GetEntry(1,3)),
                                     CD2S(FMatrix1->GetEntry(1,4)),
                                     CD2S(FMatrix1->GetEntry(1,5)));
     fprintf(f1,"\n");
     fclose(f1);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     if (XMatrix2)
      { 
        G->GetFields(0, KN, Omega, kBloch, XMatrix2, FScat);
        FILE *f2=vfopen("%s.scuff.scattered","w",GetFileBase(GeoFileName));
        for(int nr=0; nr<XMatrix2->NR; nr++)
         { 
           fprintf(f2,"%s %e %i %i ",z2s(Omega),kz,Nu,Pol);
           fprintf(f2,"%e %e %e ",XMatrix2->GetEntryD(nr,0),
                                  XMatrix2->GetEntryD(nr,1),
                                  XMatrix2->GetEntryD(nr,2));
           fprintf(f2,"%s %s %s ",CD2S(FScat->GetEntry(nr,0)),
                                  CD2S(FScat->GetEntry(nr,1)),
                                  CD2S(FScat->GetEntry(nr,2)));
           fprintf(f2,"%s %s %s ",CD2S(FScat->GetEntry(nr,3)),
                                  CD2S(FScat->GetEntry(nr,4)),
                                  CD2S(FScat->GetEntry(nr,5)));
           fprintf(f2,"\n");
         };
        fclose(f2);

      }; // if (XMatrix2)

   }; //for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)

}
