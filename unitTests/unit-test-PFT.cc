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
 * unit-test-PFT.cc -- SCUFF-EM unit tests for PFT calculations
 * 
 * homer reid       -- 11/2005 -- 10/2011
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"
#include "PFTOptions.h"
#include "libscuffInternals.h"
#include "libIncField.h"

using namespace scuff;

#define II cdouble (0.0,1.0)
#define TENTHIRDS 3.33333333333333333334
#define RT1_2     0.70710678118654752440

// values for the WhichCase parameter below 
#define CASE_PEC  0
#define CASE_SIO2 1
#define CASE_GOLD 2

// values for the ComparisonType field below
#define COMP_REL 0
#define COMP_ABS 1

/***************************************************************/
/* absolute and relative error tolerances. Quantities computed */
/* by SCUFF must have absolute or relative errors less than    */
/* this for a test to pass.                                    */
/***************************************************************/
#define RELTOL 0.25
#define ABSTOL 1.0e-10

/***************************************************************/
/* tables of power, force, and torque data for spheres of      */
/* radius R=1um and various material compositions.             */ 
/*                                                             */
/* in each case, powers and forces are for the case of         */
/* illumination by a unit-magnitude z-traveling plane wave     */
/* with left-circular polarization, ie the incident E-field is */
/*  E(x,y,z) = (1/sqrt{2}) (\vec{x} + i \vec{y}) e^{ikz}       */
/*                                                             */
/* data table entries:                                         */
/*  (1) omega (angular frequency, units 3e14 rad/sec)          */
/*  (2) real Eps(omega) (relative permittivity, dim'nless)     */
/*  (3) imag Eps(omega)                                        */
/*  (4) absorbed  power (watts)                                */
/*  (5) scattered power (watts)                                */
/*  (6) Z-force (nanoNewtons)                                  */
/*  (7) Z-torque (nanoNewtons*microns)                         */
/***************************************************************/
#define NUMOMEGAS 3
#define NUMQUANTITIES 7
typedef double PFTData[NUMOMEGAS][NUMQUANTITIES];
double OmegaList[] = { 0.01, 0.1, 1.0 };

PFTData PECPFTData =
 { { 1.0e-2,
     0.0, 0.0,
     0.0, 1.389884e-10, 7.412216e-10, -4.632946e-08
   },
   { 1.0e-1,
     0.0, 0.0,
     0.0, 1.393182e-06, 7.38044e-06, -4.643938e-05
   }, 
   { 1.0,
     0.0, 0.0,
     0.0, 8.488640e-03, 2.132454e-2,  -2.829547e-02
   }
 };

PFTData SIO2PFTData =
 { { 1.0e-2,
     3.634170e+00, 3.794466e-03,
     5.981414e-08, 2.430533e-11, 1.994615e-07, 1.993805e-05
   },
   { 1.0e-1,
     3.778579e+00, 4.786666e-02,
     7.245380e-06, 2.580820e-07, 2.500941e-05, 2.415127e-04
   },
   { 1.0,
     1.429487e+00, 3.647089e-02,
     1.774032e-04, 1.351696e-04, 9.622336e-04, 5.913442e-04
   }
 };

PFTData GoldPFTData =
 { {  1.0e-2,
     -6.610472e+04, 1.172275e+06,
      1.373257e-05, 1.297492e-10, 4.577582e-05, 4.577524e-03
   },
   {  1.0e-1,
     -5.031480e+04, 8.922668e+04,
      3.741193e-05, 1.352849e-06, 1.309154e-04, 1.247064e-03
   },
   {  1.0,
     -2.020863e+03, 3.585437e+02, 
      8.270806e-05, 8.860794e-03, 3.487239e-02, 2.756935e-04 
   }
 };

/***************************************************************/
/* SCUFF data for each frequency and each geometry are written */
/* to these buffers                                            */
/***************************************************************/
double OPFT[8], DSIPFT[8], EPPFTWithout[8], EPPFTWith[8];

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct Test
 { int QIndex;
   double *QScuff;
   const char *Name;
 } Test;

Test Tests[]=
 { { 3, &(OPFT[0]),          "absorbed power (overlap)"              }, /*0*/
   { 4, &(OPFT[1]),          "scattered power (overlap)"             }, /*1*/
   { 5, &(OPFT[4]),          "force (overlap)"                       }, /*2*/
   { 6, &(OPFT[7]),          "torque (overlap)"                      }, /*3*/
   { 3, &(DSIPFT[0]),        "absorbed power (DSI)"                  }, /*4*/
   { 4, &(DSIPFT[1]),        "scattered power (DSI)"                 }, /*5*/
   { 5, &(DSIPFT[4]),        "force (DSI)"                           }, /*6*/
   { 6, &(DSIPFT[7]),        "torque (DSI)"                          }, /*7*/
   { 3, &(EPPFTWithout[0]),  "absorbed power (EP, not precomputed)"  }, /*8*/
   { 4, &(EPPFTWithout[1]),  "scattered power (EP, not precomputed)" }, /*9*/
   { 3, &(EPPFTWith[0]),     "absorbed power (EP, precomputed)"      }, /*10*/
   { 4, &(EPPFTWith[1]),     "scattered power (EP, precomputed)"     }  /*11*/
 };

#define TEST_OPABS           0
#define TEST_OPSCAT          1
#define TEST_OFORCE          2
#define TEST_OTORQUE         3
#define TEST_DSIPABS         4
#define TEST_DSIPSCAT        5
#define TEST_DSIFORCE        6
#define TEST_DSITORQUE       7
#define TEST_EPPABS_WITHOUT  8
#define TEST_EPPSCAT_WITHOUT 9
#define TEST_EPPABS_WITH     10
#define TEST_EPPSCAT_WITH    11

#define NUMTESTS ( sizeof(Tests) / sizeof(Test) )

#define ALLOMEGA -1

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct TestToSkip
 { unsigned int WhichTest;
   int WhichCase;
   int WhichOmega; // == -1 means skip at all frequencies
 } TestToSkip;

TestToSkip TestsToSkip[]=
{ {TEST_OFORCE,    CASE_PEC,   ALLOMEGA},
  {TEST_OTORQUE,   CASE_PEC,   ALLOMEGA},
  {TEST_DSITORQUE, CASE_PEC,   ALLOMEGA},
  {TEST_DSIFORCE,  CASE_PEC,   2,      },
  {TEST_OPSCAT,    CASE_GOLD,  0,      },
  {TEST_OPSCAT,    CASE_GOLD,  1,      },
  {TEST_OFORCE,    CASE_GOLD,  0,      },
  {TEST_OFORCE,    CASE_GOLD,  1,      },
  {TEST_OTORQUE,   CASE_GOLD,  ALLOMEGA},
  {TEST_DSITORQUE, CASE_GOLD,  ALLOMEGA},
  {TEST_OPSCAT,    CASE_SIO2,  0       },
  {TEST_OFORCE,    CASE_SIO2,  0       }
};

#define NUMTESTSTOSKIP (sizeof(TestsToSkip) / sizeof(TestToSkip))

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  InstallHRSignalHandler();
  SetLogFileName("scuff-unit-tests.log");
  Log("SCUFF-EM PFT unit tests running on %s",GetHostName());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool WriteData=false;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"WriteData",  PA_BOOL, 0, 1, (void *)&WriteData, 0, "write data output file"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /***************************************************************/
  /* set up incident field ***************************************/
  /***************************************************************/
  const cdouble E0[3]  = { RT1_2, II*RT1_2, 0.0 };
  const double nHat[3] = { 0.0, 0.0, 1.0 };
  PlaneWave *PW = new PlaneWave(E0, nHat);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTOptions *Options=InitPFTOptions();

  /***************************************************************/
  /* loop over all material geometries ***************************/
  /***************************************************************/
  #define NUMCASES 3
  const char *GeoFileNames[NUMCASES] = { "PECSphere_501.scuffgeo",
                                         "SiO2Sphere_501.scuffgeo",
                                         "GoldSphere_501.scuffgeo"
                                       };
  PFTData *ExactData[NUMCASES]
   = { &PECPFTData, &SIO2PFTData, &GoldPFTData };

  int PassedTests=0, TotalTests=0;
  FILE *DataFile=WriteData ? fopen("unit-test-PFT.data","w") : 0;
  for(int nCase=0; nCase<NUMCASES; nCase++)
   { 
     Log("Testing geometry %s",GeoFileNames[nCase]);
     RWGGeometry *G = new RWGGeometry(GeoFileNames[nCase]);
     G->SetLogLevel(SCUFF_VERBOSELOGGING);
     RWGSurface  *S = G->Surfaces[0];
     HMatrix *TIn   = S->IsPEC ? 0 : G->AllocateBEMMatrix();
     HMatrix *TOut  = G->AllocateBEMMatrix();
     HMatrix *M     = G->AllocateBEMMatrix();
     HVector *KN    = G->AllocateRHSVector();
     HVector *RHS   = G->AllocateRHSVector() ;

     /***************************************************************/
     /* loop over frequencies                                       */
     /***************************************************************/
     for (int nOmega=0; nOmega<NUMOMEGAS; nOmega++)
      { 
        /*--------------------------------------------------------------*/
        /*- set up the BEM scattering problem and solve for the         */
        /*- surface currents                                            */
        /*--------------------------------------------------------------*/
        double Omega=OmegaList[nOmega];
        Log("Solving BEM system at omega=%g",Omega);
        if (S->IsPEC)
         {
           G->AssembleBEMMatrix(Omega, TOut);
           M->Copy(TOut);
         }
        else
         { for(int nr=0; nr<G->NumRegions; nr++)
            G->RegionMPs[nr]->Zero();

           G->RegionMPs[1]->UnZero();
           G->AssembleBEMMatrix(Omega, TIn);
           G->RegionMPs[1]->Zero();

           G->RegionMPs[0]->UnZero();
           G->AssembleBEMMatrix(Omega, TOut);
           G->RegionMPs[0]->Zero();
 
           for(int nr=0; nr<G->NumRegions; nr++)
            G->RegionMPs[nr]->UnZero();

           M->Copy(TOut);
           M->AddBlock(TIn,0,0);
         };
        M->LUFactorize();
        G->AssembleRHSVector(Omega, PW, RHS);
        KN->Copy(RHS);
        M->LUSolve(KN);

        /*--------------------------------------------------------------*/
        /*- compute overlap PFT                                         */
        /*--------------------------------------------------------------*/
        Log("Computing overlap PFT...");
        InitPFTOptions(Options);
        Options->PFTMethod=SCUFF_PFT_OVERLAP;
        Options->RHSVector=RHS;
        G->GetPFT(0, PW, KN, Omega, OPFT, Options);
  
        /*--------------------------------------------------------------*/
        /*- compute DSIPFT                                              */
        /*--------------------------------------------------------------*/
        Log("Computing DSIPFT...");
        InitPFTOptions(Options);
        Options->PFTMethod=SCUFF_PFT_DSI;
        G->GetPFT(0, PW, KN, Omega, DSIPFT, Options);

        /*--------------------------------------------------------------*/
        /*- compute EP power without precomputed matrices               */
        /*--------------------------------------------------------------*/
        Log("Computing EPPFT without...");
        InitPFTOptions(Options);
        Options->PFTMethod=SCUFF_PFT_EP;
        G->GetPFT(0, PW, KN, Omega, EPPFTWithout, Options);

        /*--------------------------------------------------------------*/
        /*- compute EP power with precomputed matrices                  */
        /*--------------------------------------------------------------*/
        Log("Computing EPPFT with...");
        InitPFTOptions(Options);
        Options->PFTMethod=SCUFF_PFT_EP;
        Options->TInterior=TIn;
        Options->TExterior=TOut;
        G->GetPFT(0, PW, KN, Omega, EPPFTWith, Options);

        /*--------------------------------------------------------------*/
        /*- write data to log file if requested-------------------------*/
        /*--------------------------------------------------------------*/
        if (WriteData)
         { 
           if (nOmega==0) fprintf(DataFile,"# geometry %s: \n",G->GeoFileName);
           fprintf(DataFile,"%e %e %e %e %e %e %e %e %e %e %e %e %e \n",Omega,
                             OPFT[0],OPFT[1],OPFT[4],OPFT[7],
                             DSIPFT[0],DSIPFT[1],DSIPFT[4],DSIPFT[7],
                             EPPFTWith[0],EPPFTWith[1],
                             EPPFTWithout[0],EPPFTWithout[1]);
           fflush(DataFile);
         };

        /*--------------------------------------------------------------*/
        /*- apply all tests --------------------------------------------*/
        /*--------------------------------------------------------------*/
        for(unsigned int nt=0; nt<NUMTESTS; nt++)
         { 
           bool SkipThisTest=false;
           for(unsigned int ntts=0; (!SkipThisTest) && ntts<NUMTESTSTOSKIP; ntts++)
            if (    TestsToSkip[ntts].WhichTest==nt
                 && TestsToSkip[ntts].WhichCase==nCase
                 && (    TestsToSkip[ntts].WhichOmega==ALLOMEGA
                      || TestsToSkip[ntts].WhichOmega==nOmega
                    )
               ) SkipThisTest=true;
           if (SkipThisTest) 
            continue;

           TotalTests++;
           const char *TestName = Tests[nt].Name;
           double QExact  = (*ExactData[nCase])[nOmega][Tests[nt].QIndex];
           double QScuff  = *(Tests[nt].QScuff);
           double MyRD = RD(QExact, QScuff);
           Log("Test %s: {exact,scuff,RD}={%+.4e,%+.4e,%+.1e}...",
                TestName,QExact,QScuff,MyRD);

           bool Passed=false;
           if ( QExact==0.0 )
            Passed = abs(QScuff)<ABSTOL;
           else
            Passed = RD( QExact, QScuff ) < RELTOL;

           if (Passed) 
            { PassedTests++;
              LogC("PASSED");
            }
           else
            LogC("FAILED");

         }; // for(unsigned int nt=0; nt<NUMTESTS; nt++)

      }; //for (int nOmega=0; nOmega<NUMOMEGAS; nOmega++)

     delete G; 
     if (TIn) delete TIn;
     delete TOut;
     delete M; 
     delete KN;
     delete RHS;

   }; // for(int nCase=0; nCase<NUMCASES; nCase++)

  if(DataFile)
   fclose(DataFile);
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  Log("%i/%i tests successfully passed.",PassedTests,TotalTests);
  printf("%i/%i tests successfully passed.\n",PassedTests,TotalTests);

  int FailedTests=TotalTests - PassedTests;
  if (FailedTests>0)
   abort();

  return 0;

}
