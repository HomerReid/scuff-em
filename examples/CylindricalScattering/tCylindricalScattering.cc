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
  cdouble Omega=1.0;
  char *OmegaFile=0;
  char *EPFile=0;
  int Nu=0;
  double kz=0.0;
  bool TMPolarization=false;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",  PA_STRING,  1, 1, (void *)&GeoFileName,  0,  ".scuffgeo file"},
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
  if (Cache)
   PreloadCache(Cache);

  /*--------------------------------------------------------------*/
  /* preallocate BEM matrix and RHS vector                        */
  /*--------------------------------------------------------------*/
  HMatrix *M  = G->AllocateBEMMatrix();
  HVector *KN = G->AllocateRHSVector();

  int Pol = TM_POLARIZATION ? CW_TM2Z : CW_TE2Z;
  CylindricalWave CW(Nu, kz, Pol, CW_INCOMING);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *XMatrix1=new HMatrix(2,3);
  XMatrix1->SetEntry(0,0,0.0);
  XMatrix1->SetEntry(0,1,0.0);
  XMatrix1->SetEntry(0,2,0.0);
  XMatrix1->SetEntry(1,0,2.0);
  XMatrix1->SetEntry(1,1,0.0);
  XMatrix1->SetEntry(1,2,0.0);

  HMatrix *FMatrix1=new HMatrix(2,6,LHM_COMPLEX);

  HMatrix *XMatrix2=0, *FMatrix2=0;
  if (EPFile)
   { XMatrix2=new HMatrix(EPFile);
     if (XMatrix2 && XMatrix2->ErrMsg)
      ErrExit(XMatrix2->ErrMsg);
     FMatrix2=new HMatrix(XMatrix2->NR,6,LHM_COMPLEX);
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)
   {
     Log("Computing at frequency %s...",z2s(Omega));

     /*--------------------------------------------------------------*/
     /* assemble and factorize the BEM matrix at this frequency      */
     /*--------------------------------------------------------------*/
     Omega=OmegaVector->GetEntry(nOmega);
     G->AssembleBEMMatrix(Omega, kBloch, M);
     M->LUFactorize();

     G->AssembleRHSVector(Omega, kBloch, &CW, KN);
     M->LUSolve(KN);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/ 
     G->GetFields(0, KN, Omega, kBloch, XMatrix1, FMatrix1);
     FILE *f1=vfopen("%s.out1","a",GetFileBase(GeoFileName));
     fprintf(f1,"%s %e %i %i ",z2s(Omega),kz,Nu,Pol);
     fprintf(f1,"%s %s %s ",CD2S(FMatrix1->GetEntry(0,0)),
                            CD2S(FMatrix1->GetEntry(0,1)),
                            CD2S(FMatrix1->GetEntry(0,2)));
     fprintf(f1,"%s %s %s ",CD2S(FMatrix1->GetEntry(1,0)),
                            CD2S(FMatrix1->GetEntry(1,1)),
                            CD2S(FMatrix1->GetEntry(1,2)));
     fprintf(f1,"\n");

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     if (XMatrix2)
      { 
        G->GetFields(0, KN, Omega, kBloch, XMatrix2, FMatrix2);

        FILE *f2=vfopen("%s.out2","w",GetFileBase(GeoFileName));
        for(int nr=0; nr<XMatrix2->NR; nr++)
         { 
           fprintf(f2,"%s %e %i %i ",z2s(Omega),kz,Nu,Pol);
           fprintf(f2,"%e %e %e ",XMatrix2->GetEntry(nr,0),
                                  XMatrix2->GetEntry(nr,1),
                                  XMatrix2->GetEntry(nr,2));
           fprintf(f2,"%s %s %s ",CD2S(FMatrix2->GetEntry(nr,0)),
                                  CD2S(FMatrix2->GetEntry(nr,1)),
                                  CD2S(FMatrix2->GetEntry(nr,2)));
           fprintf(f2,"%s %s %s ",CD2S(FMatrix2->GetEntry(nr,3))
                                  CD2S(FMatrix2->GetEntry(nr,4))
                                  CD2S(FMatrix2->GetEntry(nr,5)));
           fprintf(f2,"\n");
         };
        fclose(f2);
      }; // if (XMatrix2)

   }; //for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)

}
