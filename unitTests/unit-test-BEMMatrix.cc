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
 * unit-test-BEMMatrix.cc -- SCUFF-EM unit test for the full BEM matrix
 * 
 * homer reid             -- 11/2005 -- 10/2011
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"
#include "libscuffInternals.h"

using namespace scuff;

#define REF_FILENAME "unit-test-BEMMatrix.reference.hdf5"
#define FAIL_FILENAME "unit-test-BEMMatrix.failed.hdf5"

#define TESTNAME1  "PEC sphere, low real frequency"
#define TESTNAME2  "PEC sphere, medium real frequency"
#define TESTNAME3  "PEC sphere, imaginary frequency"
#define TESTNAME4  "Dielectric sphere, low real frequency"
#define TESTNAME5  "Dielectric sphere, medium real frequency"
#define TESTNAME6  "Dielectric sphere, imaginary frequency"
#define TESTNAME7  "Two dielectric spheres"
#define TESTNAME8  "Extended PEC plate, kBloch=0"
#define TESTNAME9  "Extended PEC plate, kBloch!=0"
#define TESTNAME10 "Extended dielectric slab, kBloch=0"
#define TESTNAME11 "Extended dielectric slab, kBloch!=0"
#define TESTNAME12 "Sphere slab array"
#define NUMTESTS   12 

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void CompareMatrices(HMatrix *M, HMatrix *MRef,
                     double *AvgRelError, double *MismatchRate)
{ 
  int Mismatches=0;
  double TotalRelError=0.0;
  for(int nr=0; nr<M->NR; nr++)
   for(int nc=nr; nc<M->NC; nc++)
    {  
      cdouble m      = M->GetEntry(nr,nc);
      cdouble mRef   = MRef->GetEntry(nr,nc);
      double Scale   = abs(mRef);
      if (Scale==0.0) Scale=1.0;
      TotalRelError  +=  abs(m-mRef) / Scale;
      if ( !EqualFloat(m, mRef) ) Mismatches++;
    };
  int NumEntries = M->NR * M->NC;
  *AvgRelError = TotalRelError / ((double)NumEntries);
  *MismatchRate = ((double)Mismatches) / ((double)NumEntries);

}
 
/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  SetLogFileName("scuff-test-BEMMatrix.log");
  Log("SCUFF-EM BEM matrix unit test running on %s",GetHostName());
 
  /*--------------------------------------------------------------*/
  /*- Use the --WriteFiles option to create the .hdf5 files that  */
  /*- are used as the comparison for subsequent runs.             */
  /*--------------------------------------------------------------*/
  bool WriteFiles=false;
  bool Test1=false, Test2=false,  Test3=false, Test4=false;
  bool Test5=false, Test6=false,  Test7=false, Test8=false;
  bool Test9=false, Test10=false, Test11=false, Test12=false;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { 
     {"Test1",     PA_BOOL, 0, 1, (void *)&Test1,     0, TESTNAME1},
     {"Test2",     PA_BOOL, 0, 1, (void *)&Test2,     0, TESTNAME2},
     {"Test3",     PA_BOOL, 0, 1, (void *)&Test3,     0, TESTNAME3},
     {"Test4",     PA_BOOL, 0, 1, (void *)&Test4,     0, TESTNAME4},
     {"Test5",     PA_BOOL, 0, 1, (void *)&Test5,     0, TESTNAME5},
     {"Test6",     PA_BOOL, 0, 1, (void *)&Test6,     0, TESTNAME6},
     {"Test7",     PA_BOOL, 0, 1, (void *)&Test7,     0, TESTNAME7},
     {"Test8",     PA_BOOL, 0, 1, (void *)&Test8,     0, TESTNAME8},
     {"Test9",     PA_BOOL, 0, 1, (void *)&Test9,     0, TESTNAME9},
     {"Test10",    PA_BOOL, 0, 1, (void *)&Test10,    0, TESTNAME10},
     {"Test11",    PA_BOOL, 0, 1, (void *)&Test11,    0, TESTNAME11},
     {"Test12",    PA_BOOL, 0, 1, (void *)&Test12,    0, TESTNAME12},
     {"Reference", PA_BOOL, 0, 1, (void *)&WriteFiles, 0, "write reference .hdf5 file"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool AllTests = (argc==1);

  const char *GeoFileNames[NUMTESTS];
  cdouble Omega[NUMTESTS];
  double kBloch[NUMTESTS][2];
  const char *MNames[NUMTESTS];
  const char *TestNames[NUMTESTS];
  int NumTests=0;

  if ( Test1 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "PECSphere_255.scuffgeo";
     Omega[NumTests]        = 0.01;
     MNames[NumTests]       = "M_PECSphere_w0P01";
     TestNames[NumTests]    = TESTNAME1;
     NumTests++;
   };
  if ( Test2 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "PECSphere_255.scuffgeo";
     Omega[NumTests]        = 1.0;
     MNames[NumTests]       = "M_PECSphere_w1P0";
     TestNames[NumTests]    = TESTNAME2;
     NumTests++;
   };
  if ( Test3 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "PECSphere_255.scuffgeo";
     Omega[NumTests]        = 0.1*II;
     MNames[NumTests]       = "M_PECSphere_w0P1I";
     TestNames[NumTests]    = TESTNAME3;
     NumTests++;
   };
  if ( Test4 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "SiSphere_255.scuffgeo";
     Omega[NumTests]        = 0.01;
     MNames[NumTests]       = "M_SiSphere_w0P01";
     TestNames[NumTests]    = TESTNAME4;
     NumTests++;
   };
  if ( Test5 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "SiSphere_255.scuffgeo";
     Omega[NumTests]        = 1.0;
     MNames[NumTests]       = "M_SiSphere_w1P0";
     TestNames[NumTests]    = TESTNAME5;
     NumTests++;
   };
  if ( Test6 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "SiSphere_255.scuffgeo";
     Omega[NumTests]        = 0.1*II;
     MNames[NumTests]       = "M_SiSphere_w0P1I";
     TestNames[NumTests]    = TESTNAME6;
     NumTests++;
   };
  if ( Test7 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "SiSpheres_255.scuffgeo";
     Omega[NumTests]        = 0.1;
     MNames[NumTests]       = "M_SiSpheres_w0P1";
     TestNames[NumTests]    = TESTNAME7;
     NumTests++;
   };
  if ( Test8 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "PECPlate_40.scuffgeo";
     Omega[NumTests]        = 0.1;
     kBloch[NumTests][0]    = 0.0;
     kBloch[NumTests][1]    = 0.0;
     MNames[NumTests]       = "M_PECPlate_w0P1_KZ";
     TestNames[NumTests]    = TESTNAME8;
     NumTests++;
   };
  if ( Test9 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "PECPlate_40.scuffgeo";
     Omega[NumTests]        = 1.1;
     kBloch[NumTests][0]    = 0.7;
     kBloch[NumTests][1]    = 0.9;
     MNames[NumTests]       = "M_PECPlate_w1P1_KNZ";
     TestNames[NumTests]    = TESTNAME9;
     NumTests++;
   };
  if ( Test10 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "SiSlab_40.scuffgeo";
     Omega[NumTests]        = 0.1;
     kBloch[NumTests][0]    = 0.0;
     kBloch[NumTests][1]    = 0.0;
     MNames[NumTests]       = "M_SiSlab_w0P1_KZ";
     TestNames[NumTests]    = TESTNAME10;
     NumTests++;
   };
  if ( Test11 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "SiSlab_40.scuffgeo";
     Omega[NumTests]        = 1.1;
     kBloch[NumTests][0]    = 0.7;
     kBloch[NumTests][1]    = 0.9;
     MNames[NumTests]       = "M_SiSlab_w1P1_KZ";
     TestNames[NumTests]    = TESTNAME11;
     NumTests++;
   };
  if ( Test12 || AllTests || WriteFiles)
   { GeoFileNames[NumTests] = "SphereSlabArray.scuffgeo";
     Omega[NumTests]        = 1.1;
     kBloch[NumTests][0]    = 0.7;
     kBloch[NumTests][1]    = 0.9;
     MNames[NumTests]       = "M_SphereSlabArray";
     TestNames[NumTests]    = TESTNAME12;
     NumTests++;
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  void *pHC = 0;
  if (WriteFiles)
   { pHC = HMatrix::OpenHDF5Context(REF_FILENAME);
     if (pHC==0)
      ErrExit("could not open reference file %s",REF_FILENAME);
   }
  else
   { FILE *f=fopen(REF_FILENAME,"r");
     if (f) 
      fclose(f);
     else
      ErrExit("could not open reference file %s",REF_FILENAME);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool Success=true;
  for(int nt=0; nt<NumTests; nt++)
   { 
     RWGGeometry *G = new RWGGeometry(GeoFileNames[nt]);
     HMatrix *M = G->AllocateBEMMatrix();
     printf("Test %i (%s): \n",nt,TestNames[nt]);

     if (G->LDim==0)
      G->AssembleBEMMatrix(Omega[nt], M);
     else
      G->AssembleBEMMatrix(Omega[nt], kBloch[nt], M);

     if (WriteFiles)
      M->ExportToHDF5(pHC, MNames[nt]);
     else
      { HMatrix *MRef = G->AllocateBEMMatrix();
        MRef->ImportFromHDF5(REF_FILENAME, MNames[nt]);
        if (MRef->ErrMsg)
         Warn("could not find matrix %s in %s (skipping test)",
               MNames[nt],REF_FILENAME);
        double AvgRelError, MismatchRate;
        CompareMatrices(M, MRef, &AvgRelError, &MismatchRate);
        if ( AvgRelError>1.0e-6 || MismatchRate>0.1 )
         { Success=false;       
           printf(" FAILED ");
           if (pHC==0) pHC = HMatrix::OpenHDF5Context(FAIL_FILENAME);
           M->ExportToHDF5(pHC, MNames[nt]);
         }
        else
         printf(" PASSED ");
        printf(" (AvgRelErr = %.1e, Mismatch rate = %.g %%)\n",
                 AvgRelError, 100.0*MismatchRate);
        delete MRef;
      };
     delete M;
   };

  if (pHC)
   { HMatrix::CloseHDF5Context(pHC);
     if (!WriteFiles)
      printf(" Failed matrices written to file %s.\n",FAIL_FILENAME);
   };

  if (Success) 
   exit(0);
  else
   exit(1);

}
