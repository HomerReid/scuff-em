/*}} Copyright (C) 2005-2011 M. T. Homer Reid
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
 * tGTwiddle -- libSubstrate unit test 
 */

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>
#include <unistd.h>

#include "libhrutil.h"
#include "libSubstrate.h"
#include "libscuff.h"

#define NUMTESTS 3
const char *TestSubstrates[NUMTESTS]=
 { "0.0 CONST_EPS_11.7\n",
   "0.0 CONST_EPS_11.7\n -1.0 VACUUM\n",
   "0.0 CONST_EPS_11.7\n -1.0 GROUNDPLANE\n"
 };
const char *TestNames[NUMTESTS]=
 { "Si half space",
   "Si slab",
   "Si slab with ground plane",
 };

using namespace scuff;

#define II cdouble(0.0,1.0)

/***************************************************************/
/* test satisfaction of boundary conditions at material        */
/* interface #ni                                               */
/***************************************************************/
bool TestInterface(LayeredSubstrate *S, int ni, cdouble Omega, cdouble q2D[2],
                   double Eta=1.0e-6, double Threshold=1.0e-6)
{ 
  
  int NI       = S->NumInterfaces;
  double zGP   = S->zGP;
  double zDest = S->zInterface[ni];

  // set length scale to be minimum width of any substrate layer
  double zScale=1.0;
  for(int n=1; n<NI; n++)
   zScale = fmin(zScale, S->zInterface[n-1] - S->zInterface[n]);

  // material properties of layers above and below interface #ni
  S->UpdateCachedEpsMu(Omega);
  cdouble EpsAbove = S->EpsLayer[ni];
  cdouble EpsBelow = S->EpsLayer[ni+1];
  cdouble MuAbove  = S->MuLayer[ni];
  cdouble MuBelow  = S->MuLayer[ni+1];

  // test with source points in all layers
  bool AllPassed=true;
  for(int nl=0; nl<S->NumLayers; nl++)
   {
     double zSource;
     if (nl==0)
      zSource=S->zInterface[0] + 1.0;
     else if (nl<NI)
      zSource=0.5*(S->zInterface[nl-1] + S->zInterface[nl]);
     else if (!isinf(zGP))
      zSource = 0.5*(S->zInterface[NI-1] + zGP);
     else
      zSource = S->zInterface[nl-1] - 1.0;

     HMatrix GAbove(6,6,LHM_COMPLEX), GBelow(6,6,LHM_COMPLEX);
     S->GetScriptGTwiddle(Omega, q2D, zDest + Eta*zScale, zSource, &GAbove, 0, false, false, true);
     S->GetScriptGTwiddle(Omega, q2D, zDest - Eta*zScale, zSource, &GBelow, 0, false, false, true);

     /***************************************************************/
     /* scale fields appropriately: *********************************/
     /*  (1) scale E/H fields by ZVAC so magnitudes are comparable  */
     /*  (2) scale z components of E/H by Epsilon/Mu to yield D/B   */
     /***************************************************************/
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { GAbove.ScaleEntry(0+Mu, 0+Nu, 1.0/ZVAC);
         GBelow.ScaleEntry(0+Mu, 0+Nu, 1.0/ZVAC);
         GAbove.ScaleEntry(3+Mu, 3+Nu, ZVAC);
         GBelow.ScaleEntry(3+Mu, 3+Nu, ZVAC);
       }

     for(int Nu=0; Nu<6; Nu++)
      { GAbove.ScaleEntry(2, Nu, EpsAbove);
        GBelow.ScaleEntry(2, Nu, EpsBelow);
        GAbove.ScaleEntry(5, Nu, MuAbove);
        GBelow.ScaleEntry(5, Nu, MuBelow);
      }

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     double RefVal=0.0;
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       RefVal = fmax(RefVal, 0.5*( abs(GAbove.GetEntry(Mu,Nu)) + abs(GBelow.GetEntry(Mu,Nu)) ) );

     printf("Interface %i, source in layer %i (zSource=%g): ref val %e \n",ni,nl,zSource,RefVal);
     printf("         | %23s | %23s | RD\n","Above","Below");
     SetDefaultCD2SFormat("{%+.3e,%+.3e}");
     bool ThisPassed=true;
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       { cdouble X=GAbove.GetEntry(Mu,Nu);
         cdouble Y=GBelow.GetEntry(Mu,Nu);
         double RelDiff = abs(X-Y) / RefVal;
         bool MuNuPassed = RelDiff<=Threshold;
         ThisPassed &= MuNuPassed;
         if (!MuNuPassed) printf("\033[1;31m");
         const char EH[3]="EH", xyz[4]="xyz";
         printf("%c%c%c%c | %s | %s | %e\n",EH[Mu/3],EH[Nu/3],xyz[Mu%3],xyz[Nu%3],CD2S(X),CD2S(Y),RelDiff);
         if (!MuNuPassed) printf("\033[0m");
      };
     if (!ThisPassed)
      Warn(" test failed for interface %i with source in layer %i",ni,nl);
     AllPassed&=ThisPassed;

   } // //for(int nl=0; nl<S->NumLayers; nl++)

  return AllPassed;
}

/***************************************************************/
/* test satisfaction of boundary conditions at ground plane    */
/***************************************************************/
bool TestGroundPlane(LayeredSubstrate *S, cdouble Omega, cdouble q2D[2],
                     double Eta=1.0e-6, double Threshold=1.0e-6)
{ 
  if (isinf(S->zGP))
   return true;

  int NI=S->NumInterfaces;
  int NL=S->NumLayers;
  double zGP=S->zGP;

  // set length scale to be minimum width of any substrate layer
  double zScale=1.0;
  for(int n=1; n<NI; n++)
   zScale = fmin(zScale, S->zInterface[n-1] - S->zInterface[n]);

  double zDest = S->zGP + Eta*zScale;
  double zRef = 1.0;
  if (NI>0)
   zRef = 0.5*(S->zGP + S->zInterface[NI-1]);

  bool AllPassed=true;
  for(int nl=0; nl<S->NumLayers; nl++)
   {
     double zSource;
     if (nl==0)
      zSource = (NI>0 ? S->zInterface[0] : zGP) + 1.0;
     else if ( nl < NL )
      zSource=0.5*(S->zInterface[nl-1] + S->zInterface[nl]);
     else
      zSource = 0.5*(zGP + S->zInterface[NI-1]);

     HMatrix GGP(6,6,LHM_COMPLEX), GRef(6,6,LHM_COMPLEX);
     S->GetScriptGTwiddle(Omega, q2D, zDest, zSource, &GGP, 0, false, false, true);
     S->GetScriptGTwiddle(Omega, q2D, zRef,  zSource, &GRef, 0, false, false, true);

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { GGP.ScaleEntry(0+Mu, 0+Nu, 1.0/ZVAC);
         GRef.ScaleEntry(0+Mu, 0+Nu, 1.0/ZVAC);
         GGP.ScaleEntry(3+Mu, 3+Nu, ZVAC);
         GRef.ScaleEntry(3+Mu, 3+Nu, ZVAC);
       }

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     double RefVal=0.0;
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       RefVal = fmax(RefVal, abs(GRef.GetEntry(Mu,Nu)));

     printf("GP, source in layer %i (zSource=%g): ref val %e \n",nl,zSource,RefVal);
     printf("    | %23s | abs/ref \n","G(zGP)");
     for(int Mu=0; Mu<6; Mu++)
      for(int Nu=0; Nu<6; Nu++)
       { if (Mu==2 || Mu==5) continue;
         cdouble X=GGP.GetEntry(Mu,Nu);
         double RelDiff = abs(X) / RefVal;
         if (RelDiff>Threshold) printf("\033[1;31m");
         const char EH[3]="EH", xyz[4]="xyz";
         printf("%c%c%c%c | %s | %e\n",EH[Mu/3],EH[Nu/3],xyz[Mu%3],xyz[Nu%3],CD2S(X),RelDiff);
         if (RelDiff>Threshold) printf("\033[0m");
         AllPassed &= (RelDiff<=Threshold);
      };

   } // //for(int nl=0; nl<S->NumLayers; nl++)

  return AllPassed;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool RunUnitTest(LayeredSubstrate *S, cdouble Omega, cdouble q2D[2], 
                  double Eta=1.0e-6, double Threshold=1.0e-6)
{

  bool AllPassed=true;
  for(int ni=0; ni<S->NumInterfaces; ni++)
   AllPassed &= TestInterface(S, ni, Omega, q2D, Eta, Threshold);
  if (!isinf(S->zGP))
   AllPassed &= TestGroundPlane(S, Omega, q2D, Eta, Threshold);
  return AllPassed;
}

bool RunUnitTest(int nTest, cdouble Omega, cdouble q2D[2], 
                double Eta=1.0e-6, double Threshold=1.0e-6)
{ 
  LayeredSubstrate *S = CreateLayeredSubstrate(TestSubstrates[nTest]);
  bool Status = RunUnitTest(S, Omega, q2D, Eta, Threshold);
  delete S;
  return Status;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  int WhichTest=-1;
  char *SubstrateFile=0;
  char *XDSFile=0;
  cdouble Omega=0.7;
  cdouble q2D[2]={0.1, 0.2};
  double Eta=1.0e-6;
  double Threshold=1.0e-2;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Test",          PA_INT,     1, 1, (void *)&WhichTest,     0, "index of unit test"},
     {"SubstrateFile", PA_STRING,  1, 1, (void *)&SubstrateFile, 0, ".substrate file"},
     {"XDSFile",       PA_STRING,  1, 1, (void *)&XDSFile,       0, "XDS file"},
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,         0, "angular frequency"},
     {"q2D",           PA_CDOUBLE, 2, 1, (void *)q2D,            0, "qx qy"},
     {"Eta",           PA_DOUBLE,  1, 1, (void *)&Eta,           0, ""},
     {"Threshold",     PA_DOUBLE,  1, 1, (void *)&Threshold,     0, ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*--------------------------------------------------------------*/
  // if no specific test or substrate file was specified, run all unit tests
  /*--------------------------------------------------------------*/
  bool AllPassed=true;
  if (SubstrateFile==0 && WhichTest==-1)
   { for(int nTest=0; nTest<NUMTESTS; nTest++)
       AllPassed &= RunUnitTest(nTest, Omega, q2D, Eta, Threshold);
     return AllPassed;
   };

  LayeredSubstrate *S;
  if (WhichTest!=-1)
   S=CreateLayeredSubstrate(TestSubstrates[WhichTest]);
  else
   S=new LayeredSubstrate(SubstrateFile);

  AllPassed=RunUnitTest(S, Omega, q2D, Eta, Threshold);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if(XDSFile)
   { HMatrix *XDSMatrix=new HMatrix(XDSFile);
     FILE *f=fopen("/tmp/tGTwiddle.out","w");
     for(int nx=0; nx<XDSMatrix->NR; nx++)
      { double zDest = XDSMatrix->GetEntryD(nx,2), zSource = XDSMatrix->GetEntryD(nx,5);
        HMatrix GTwiddle(6,6,LHM_COMPLEX);
        S->GetScriptGTwiddle(Omega, q2D, zDest, zSource, &GTwiddle);
        cdouble EpsRel = S->EpsLayer[S->GetLayerIndex(zDest)];
        cdouble MuRel  = S->MuLayer[S->GetLayerIndex(zDest)];
        fprintf(f,"%e %e ",zDest,zSource);
        SetDefaultCD2SFormat("%e %e");
        for(int Mu=0; Mu<6; Mu++)
         for(int Nu=0; Nu<6; Nu++)
          { cdouble Factor = ( (Mu==2) ? EpsRel : (Mu==5) ? MuRel : 1.0);
            fprintf(f,"%s ",CD2S(Factor*GTwiddle.GetEntry(Mu,Nu)));
          }
        cdouble Gamma0Twiddle[6][6];
        memset(Gamma0Twiddle,0,36*sizeof(cdouble));
        if ( S->GetLayerIndex(zDest) == S->GetLayerIndex(zSource) )
         S->GetGamma0Twiddle(Omega, q2D, zDest, zSource, Gamma0Twiddle);
        for(int Mu=0; Mu<6; Mu++)
         for(int Nu=0; Nu<6; Nu++)
          { cdouble Factor = ( (Mu==2) ? EpsRel : (Mu==5) ? MuRel : 1.0);
            fprintf(f,"%s ",CD2S(Factor*Gamma0Twiddle[Mu][Nu]));
          }
        fprintf(f,"\n");
      }
     fclose(f);
   }
 
 return AllPassed;
}
