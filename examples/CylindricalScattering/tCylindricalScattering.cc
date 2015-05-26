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
  XMatrix1->SetEntry(1,1,0.0);
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
