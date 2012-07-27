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
 * scuff-transmission.cc -- a command-line tool to compute transmission
 *                          through thin films and metamaterial arrays
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libIncField.h>
#include "../PBCGeometry.h"

using namespace scuff;

extern int GetFieldObject;

#define II cdouble (0.0, 1.0)

#define MAXFREQ 10
#define MAXCACHE 10    // max number of cache files for preload

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  EnableAllCPUs();
  InstallHRSignalHandler();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;
  double LBV1[2]={1.0, 0.0};
  double LBV2[2]={0.0, 1.0};
  cdouble OmegaVals[MAXFREQ];	int nOmegaVals;
  char *OmegaFile=0;
  double ThetaMin=0.0;
  double ThetaMax=M_PI/2.0;
  double ThetaPoints=25;
  double XAbove[3]={0.0, 0.0, +2.0};
  double XBelow[3]={0.0, 0.0, -2.0};
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",    PA_STRING,  1, 1,       (void *)&GeoFileName,  0,       ".scuffgeo file"},
     {"LBV1",        PA_DOUBLE,  2, 1,       (void *)LBV1,          0,       "lattice basis vector 1"},
     {"LBV2",        PA_DOUBLE,  2, 1,       (void *)LBV2,          0,       "lattice basis vector 2"},
/**/
     {"Omega",       PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,     &nOmegaVals, "(angular) frequency"},
     {"OmegaFile",   PA_STRING,  1, 1,       (void *)&OmegaFile,    0,       "list of (angular) frequencies"},
/**/
     {"ThetaMin",    PA_DOUBLE,  1, 1,       (void *)&ThetaMin,     0,       "minimum incident angle"},
     {"ThetaMax",    PA_DOUBLE,  1, 1,       (void *)&ThetaMax,     0,       "maximum incident angle"},
     {"ThetaPoints", PA_INT,     1, 1,       (void *)&ThetaPoints,  0,       "number of incident angles"},
/**/
     {"XAbove",      PA_DOUBLE,  3, 1,       (void *)XAbove,        0,       "upper evaluation point"},
     {"XBelow",      PA_DOUBLE,  3, 1,       (void *)XBelow,        0,       "lower evaluation point"},
/**/
     {"Cache",       PA_STRING,  1, 1,       (void *)&Cache,        0,             "read/write cache"},
     {"ReadCache",   PA_STRING,  1, MAXCACHE,(void *)ReadCache,     &nReadCache,   "read cache"},
     {"WriteCache",  PA_STRING,  1, 1,       (void *)&WriteCache,   0,             "write cache"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run calculations                        */
  /*******************************************************************/
  HVector *OmegaVector=0, *OmegaVector0;
  int nFreq, nOV, NumFreqs=0;
  if (OmegaFile)
   { 
     OmegaVector=new HVector(OmegaFile,LHM_TEXT);
     if (OmegaVector->ErrMsg)
      ErrExit(OmegaVector->ErrMsg);
     NumFreqs=OmegaVector->N;
   };

  // now add any individually specified --Omega options
  if (nOmegaVals>0)
   { 
     NumFreqs += nOmegaVals;
     OmegaVector0=OmegaVector;
     OmegaVector=new HVector(NumFreqs, LHM_COMPLEX);
     nFreq=0;
     if (OmegaVector0)
      { for(nFreq=0; nFreq<OmegaVector0->N; nFreq++)
         OmegaVector->SetEntry(nFreq, OmegaVector0->GetEntry(nFreq));
        delete OmegaVector0;
      };
     for(nOV=0; nOV<nOmegaVals; nOV++)
      OmegaVector->SetEntry(nFreq+nOV, OmegaVals[nOV]);
   };

  if ( !OmegaVector || OmegaVector->N==0)
   OSUsage(argv[0], OSArray, "you must specify at least one frequency");

  /*******************************************************************/
  /* process incident-angle-related options to construct a list of   */
  /* incident angles at which to run calculations                    */
  /*******************************************************************/
  if ( ThetaMax<ThetaMin )
   OSUsage(argv[0], OSArray, "--ThetaMin must not be greater than --ThetaMax");
  if ( ThetaMax==ThetaMin )  
   ThetaPoints=1;
  HVector *ThetaVector=LinSpace(ThetaMin, ThetaMax, ThetaPoints);
  
  /*******************************************************************/
  /*- create the RWGGeometry and PBCGeometry structures             -*/
  /*******************************************************************/
  SetLogFileName("%s.log",GetFileBase(GeoFileName));
  RWGGeometry *G=new RWGGeometry(GeoFileName);

  double *LBV[2]={LBV1, LBV2};
  PBCGeometry *PG=new PBCGeometry(G, LBV);

  HMatrix *M  = PG->AllocateBEMMatrix();
  HVector *KN = PG->AllocateRHSVector();

  /*******************************************************************/
  /* preload the scuff cache with any cache preload files the user   */
  /* may have specified                                              */
  /*******************************************************************/
  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");
  if (Cache) 
   WriteCache=Cache;
  for (int nrc=0; nrc<nReadCache; nrc++)
   PreloadCache( ReadCache[nrc] );
  if (Cache)
   PreloadCache( Cache );

  /*******************************************************************/
  /*- create the incident field                                      */
  /*******************************************************************/
  cdouble E0[3]={1.0, 0.0, 0.0};
  double nHat[3]={0.0, 0.0, 1.0};
  PlaneWave PW(E0, nHat);

  cdouble EpsExterior, MuExterior, kExterior;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char OutFileName[1000];
  snprintf(OutFileName,1000,"%s.transmission",GetFileBase(GeoFileName));
  FILE *f=CreateUniqueFile(OutFileName,1);
  setlinebuf(f);
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1: omega \n");
  fprintf(f,"# 2: theta (incident angle) (theta=0 --> normal incidence)\n");
  fprintf(f,"# 3: mag (rPerp) \n");
  fprintf(f,"# 4: phase (rPerp) \n");
  fprintf(f,"# 5: mag (rPar) \n");
  fprintf(f,"# 6: phase (rPar) \n");
  fprintf(f,"# 7: mag (tPerp) \n");
  fprintf(f,"# 8: phase (tPerp) \n");
  fprintf(f,"# 9: mag (tPar) \n");
  fprintf(f,"# 10: phase (tPar) \n");

  /*--------------------------------------------------------------*/
  /*- loop over frequencies and incident angles ------------------*/
  /*--------------------------------------------------------------*/
  double BlochP[2];
  double Theta, SinTheta, CosTheta, SinThetaP, CosThetaP;
  cdouble Omega;
  cdouble EHAbovePerp[6], EHBelowPerp[6], EHAbovePar[6], EHBelowPar[6];
  cdouble ExpFacAbove, ExpFacBelow;
  cdouble rPerp, rPar, tPerp, tPar;
  for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)
   for(int nTheta=0; nTheta<ThetaVector->N; nTheta++)
    { 
      Omega = OmegaVector->GetEntry(nOmega);

      G->ExteriorMP->GetEpsMu(Omega, &EpsExterior, &MuExterior);
      kExterior = csqrt2(EpsExterior*MuExterior)*Omega;

      Theta = ThetaVector->GetEntryD(nTheta);
      SinTheta=sin(Theta);
      CosTheta=cos(Theta);

      // set bloch wavevector and assemble BEM matrix 
      BlochP[0] = real(kExterior)*SinTheta;
      BlochP[1] = 0.0;
      PG->AssembleBEMMatrix(Omega, BlochP, M);
      M->LUFactorize();

     /*******************************************************************/
     /*******************************************************************/
     /*******************************************************************/
     if (WriteCache)
      { StoreCache( WriteCache );
        WriteCache=0;       
      };

      // set plane wave direction 
      nHat[0] = +SinTheta;
      nHat[1] = 0.0;
      nHat[2] = -CosTheta;
      PW.SetnHat(nHat);

      // solve with E-field perpendicular to plane of incidence
      E0[0]=0.0;
      E0[1]=1.0;
      E0[2]=0.0;
      PW.SetE0(E0);
      PG->AssembleRHSVector(Omega, &PW, KN);
      M->LUSolve(KN);

GetFieldObject=-1;
     PG->GetFields(0, KN, Omega, BlochP, XAbove, EHAbovePerp);
GetFieldObject=-1;
     PG->GetFields(0, KN, Omega, BlochP, XBelow, EHBelowPerp);

     // solve with E-field parallel to plane of incidence
     E0[0]=CosTheta;
     E0[1]=0.0;
     E0[2]=SinTheta;
     PW.SetE0(E0);
     PG->AssembleRHSVector(Omega, &PW, KN);
     M->LUSolve(KN);
GetFieldObject=-1;
     PG->GetFields(0, KN, Omega, BlochP, XAbove, EHAbovePar);
GetFieldObject=-1;
     PG->GetFields(0, KN, Omega, BlochP, XBelow, EHBelowPar);

     // extract r- and t- coefficients 
     ExpFacAbove = exp( II*kExterior*(SinTheta*XAbove[0] + CosTheta*XAbove[2]) );
     rPerp = - EHAbovePerp[1] / ExpFacAbove;
     rPar  = - (CosTheta*EHAbovePar[0] - SinTheta*EHAbovePar[2]) / ExpFacAbove;

     ExpFacBelow= exp( II*kExterior*(SinTheta*XBelow[0] - CosTheta*XBelow[2]) );
     tPerp = EHBelowPerp[1] / ExpFacBelow;
     tPar  = (CosThetaP*EHBelowPar[0] + SinThetaP*EHBelowPar[2]) / ExpFacBelow;

     fprintf(f,"%s %e %e %e %e %e %e %e %e %e\n",z2s(Omega),Theta,
                abs(rPerp), arg(rPerp), abs(rPar), arg(rPar), 
                abs(tPerp), arg(tPerp), abs(tPar), arg(tPar));

   }; 

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Thank you for your support.\n");

}
