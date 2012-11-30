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
 *                       -- through thin films and metamaterial arrays
 *
 * homer reid            -- 9/2012
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libIncField.h>

using namespace scuff;

#define II cdouble (0.0, 1.0)

#define MAXFREQ 10
#define MAXCACHE 10    // max number of cache files for preload

/*******************************************************************/
/* this is a 9th-order, 17-point cubature rule for the unit square */
/* with corners {(0,0) (1,0) (1,1) (1,0)}.                         */
/* array entries:                                                  */
/*  x_0, y_0, w_0,                                                 */
/*  x_1, y_1, w_1,                                                 */
/*  ...                                                            */
/*  x_16, y_16, w_16                                               */
/* where (x_n, y_n) and w_n are the nth cubature point and weight. */
/*******************************************************************/
double SCR9[]={
  +5.0000000000000000e-01, +5.0000000000000000e-01, +1.3168724279835392e-01,
  +9.8442498318098881e-01, +8.1534005986583447e-01, +2.2219844542549678e-02,
  +9.8442498318098881e-01, +1.8465994013416559e-01, +2.2219844542549678e-02,
  +1.5575016819011134e-02, +8.1534005986583447e-01, +2.2219844542549678e-02,
  +1.5575016819011134e-02, +1.8465994013416559e-01, +2.2219844542549678e-02,
  +8.7513854998945029e-01, +9.6398082297978482e-01, +2.8024900532399120e-02,
  +8.7513854998945029e-01, +3.6019177020215176e-02, +2.8024900532399120e-02,
  +1.2486145001054971e-01, +9.6398082297978482e-01, +2.8024900532399120e-02,
  +1.2486145001054971e-01, +3.6019177020215176e-02, +2.8024900532399120e-02,
  +7.6186791010721466e-01, +7.2666991056782360e-01, +9.9570609815517519e-02,
  +7.6186791010721466e-01, +2.7333008943217640e-01, +9.9570609815517519e-02,
  +2.3813208989278534e-01, +7.2666991056782360e-01, +9.9570609815517519e-02,
  +2.3813208989278534e-01, +2.7333008943217640e-01, +9.9570609815517519e-02,
  +5.3810416409630857e-01, +9.2630786466683113e-01, +6.7262834409945196e-02,
  +5.3810416409630857e-01, +7.3692135333168873e-02, +6.7262834409945196e-02,
  +4.6189583590369143e-01, +9.2630786466683113e-01, +6.7262834409945196e-02,
  +4.6189583590369143e-01, +7.3692135333168873e-02, +6.7262834409945196e-02 
};

/*******************************************************************/
/* get the transmitted and reflected flux by integrating the       */
/* scattered poynting vector over the area of the unit cell.       */
/*                                                                 */
/* more specifically, the transmitted power is the integral of     */
/* the upward-directed poynting vector at ZAbove, while the        */
/* reflected power is the integral of the downward-directed        */
/* poynting vector at ZBelow; the flux is the power divided by     */
/* the area of the unit cell.                                      */
/*                                                                 */
/* Return values: TRFlux[0,1] = transmitted, reflected flux        */
/*******************************************************************/
void GetTRFlux(RWGGeometry *G, IncField *IF, HVector *KN, cdouble Omega, 
               double *kBloch, double ZAbove, double ZBelow, double *TRFlux)
{
  int NCP=17; // number of cubature points
  double *SCR=SCR9;

  /***************************************************************/
  /* on the first invocation we allocate space for the matrices  */
  /* of evaluation points and fields.                            */
  /***************************************************************/
  static HMatrix *XMatrixAbove = 0, *XMatrixBelow = 0;
  static HMatrix *FMatrixAbove = 0, *FMatrixBelow = 0;
  if (XMatrixAbove==0)
   { XMatrixAbove = new HMatrix(NCP, 3 ); 
     XMatrixBelow = new HMatrix(NCP, 3 ); 
     FMatrixAbove = new HMatrix(NCP, 6, LHM_COMPLEX);
     FMatrixBelow = new HMatrix(NCP, 6, LHM_COMPLEX);
   };

  /***************************************************************/ 
  /* fill in coordinates of evaluation points.                   */ 
  /* the first NCP points are for the upper surface; the next    */ 
  /* NCP points are for the lower surface.                       */ 
  /***************************************************************/ 
  double x, y, *LBV[2];
  LBV[0]=G->LatticeBasisVectors[0];
  LBV[1]=G->LatticeBasisVectors[1];
  for(int ncp=0; ncp<NCP; ncp++)
   { 
     x=SCR[3*ncp+0];
     y=SCR[3*ncp+1];

     XMatrixAbove->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
     XMatrixAbove->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
     XMatrixAbove->SetEntry(ncp, 2, ZAbove);

     XMatrixBelow->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
     XMatrixBelow->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
     XMatrixBelow->SetEntry(ncp, 2, ZBelow);

   };

  /***************************************************************/ 
  /* get scattered fields at all cubature points                 */ 
  /***************************************************************/ 
  G->GetFields(0, KN, Omega, kBloch, XMatrixAbove, FMatrixAbove);
  G->GetFields(0, KN, Omega, kBloch, XMatrixBelow, FMatrixBelow);

  /***************************************************************/
  /* integrate poynting vector over upper and lower surfaces.    */
  /* Note: The jacobian in this cubature should be the area of   */
  /*       the unit cell, so omitting that factor is equivalent  */
  /*       to dividing the integrated power by the unit-cell     */
  /*       area, which is what we want to do anyway.             */
  /***************************************************************/
  double w, PTransmitted=0.0, PReflected=0.0;
  cdouble E[3], H[3];
  for(int ncp=0; ncp<NCP; ncp++)
   {
     w=SCR[3*ncp+2];  // cubature weight

     E[0]=FMatrixAbove->GetEntry(ncp, 0);
     E[1]=FMatrixAbove->GetEntry(ncp, 1);
     H[0]=FMatrixAbove->GetEntry(ncp, 3);
     H[1]=FMatrixAbove->GetEntry(ncp, 4);
     PTransmitted += 0.5*w*real( E[0]*conj(H[1]) - E[1]*conj(H[0]) );

     E[0]=FMatrixBelow->GetEntry(ncp, 0);
     E[1]=FMatrixBelow->GetEntry(ncp, 1);
     H[0]=FMatrixBelow->GetEntry(ncp, 3);
     H[1]=FMatrixBelow->GetEntry(ncp, 4);
     PReflected -= 0.5*w*real( E[0]*conj(H[1]) - E[1]*conj(H[0]) );

   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  TRFlux[0]=PTransmitted;
  TRFlux[1]=PReflected;

}

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
  cdouble OmegaVals[MAXFREQ];	int nOmegaVals;
  char *OmegaFile=0;
  double Theta=0.0;
  double ThetaMin=0.0;
  double ThetaMax=0.0;
  int ThetaPoints=25;
  double ZAbove=2.0;
  double ZBelow=-1.0;
  char *OutFileName=0;
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",    PA_STRING,  1, 1,       (void *)&GeoFileName,  0,       ".scuffgeo file"},
/**/
     {"Omega",       PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,     &nOmegaVals, "(angular) frequency"},
     {"OmegaFile",   PA_STRING,  1, 1,       (void *)&OmegaFile,    0,       "list of (angular) frequencies"},
/**/
     {"Theta",       PA_DOUBLE,  1, 1,       (void *)&Theta,        0,       "incident angle in degrees"},
     {"ThetaMin",    PA_DOUBLE,  1, 1,       (void *)&ThetaMin,     0,       "minimum incident angle in degrees"},
     {"ThetaMax",    PA_DOUBLE,  1, 1,       (void *)&ThetaMax,     0,       "maximum incident angle in degrees"},
     {"ThetaPoints", PA_INT,     1, 1,       (void *)&ThetaPoints,  0,       "number of incident angles"},
/**/
     {"ZAbove",      PA_DOUBLE,  1, 1,       (void *)&ZAbove,       0,       "Z-coordinate of upper integration plane"},
     {"ZBelow",      PA_DOUBLE,  1, 1,       (void *)&ZBelow,       0,       "Z-coordinate of lower integration plane"},
/**/
     {"OutFile",     PA_STRING,  1, 1,       (void *)&OutFileName,  0,       "output file name"},
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
  /*- create the RWGGeometry                                        -*/
  /*******************************************************************/
  SetLogFileName("scuff-transmission.log");
  RWGGeometry::AssignBasisFunctionsToExteriorEdges=false;
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  if (G->NumLatticeBasisVectors!=2)
   ErrExit("%s: geometry must have two-dimensional lattice periodicity",GeoFileName);
  HMatrix *M  = G->AllocateBEMMatrix();
  HVector *KN = G->AllocateRHSVector();

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run calculations                        */
  /*******************************************************************/
  HVector *OmegaVector=0;
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
     HVector *OmegaVector0=OmegaVector;
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
  /* incident angles at which to run calculations.                   */
  /* Note: The --ThetaMin and --ThetaMax arguments are interpreted   */
  /*       in degrees (i.e. they should be numbers between 0 and 90),*/
  /*       but internally the entries of the ThetaVector vector are  */
  /*       in radians (values between 0 and Pi/2).                   */
  /*******************************************************************/
  HVector *ThetaVector;
  if ( ThetaMax!=0.0 )
   { if ( ThetaMax<ThetaMin )
      OSUsage(argv[0], OSArray, "--ThetaMin must not be greater than --ThetaMax");
     if ( ThetaMax==ThetaMin )  
      ThetaPoints=1;
     ThetaVector=LinSpace(ThetaMin*DEG2RAD, ThetaMax*DEG2RAD, ThetaPoints);
   }
  else 
   { ThetaVector=new HVector(1);
     ThetaVector->SetEntry(0,Theta*DEG2RAD);
   };
  if (ThetaPoints==1)
   Log("Calculating at the single incident angle Theta=%e degrees",ThetaVector->GetEntryD(0)*RAD2DEG);
  else 
   Log("Calculating at %i incident angles in range Theta=(%e,%e) degrees",
       ThetaVector->N, ThetaVector->GetEntryD(0)*RAD2DEG,
       ThetaVector->GetEntryD(ThetaVector->N-1)*RAD2DEG);

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

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f;
  if (OutFileName)
   { f=fopen(OutFileName,"w");
     if (!f) ErrExit("could not open file %s",f);
   }
  else
   { char buffer[1000];
     snprintf(buffer,1000,"%s.transmission",GetFileBase(GeoFileName));
     f=CreateUniqueFile(buffer,1,buffer);
     OutFileName=strdup(buffer);
   };
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1: omega \n");
  fprintf(f,"# 2: theta (incident angle) (theta=0 --> normal incidence)\n");
  fprintf(f,"# 3: transmitted flux / incident flux (perpendicular polarization)\n");
  fprintf(f,"# 4: reflected flux   / incident flux (perpendicular polarization)\n");
  fprintf(f,"# 5: transmitted flux / incident flux (parallel polarization)\n");
  fprintf(f,"# 6: reflected flux   / incident flux (parallel polarization)\n");
  fflush(f);

  cdouble EpsExterior, MuExterior, kExterior;

  /*--------------------------------------------------------------*/
  /*- loop over frequencies and incident angles ------------------*/
  /*--------------------------------------------------------------*/
  double kBloch[2];
  double SinTheta, CosTheta;
  cdouble Omega;
  double TRFluxPerp[2], TRFluxPar[2], IncFlux;
  for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)
   for(int nTheta=0; nTheta<ThetaVector->N; nTheta++)
    { 
      Omega = OmegaVector->GetEntry(nOmega);
      Theta = ThetaVector->GetEntryD(nTheta);
      SinTheta=sin(Theta);
      CosTheta=cos(Theta);
      Log("Solving the scattering problem at (Omega,Theta)=(%g,%g)",real(Omega),Theta*RAD2DEG);

      // set bloch wavevector and assemble BEM matrix 
      G->RegionMPs[0]->GetEpsMu(Omega, &EpsExterior, &MuExterior);
      kExterior = csqrt2(EpsExterior*MuExterior)*Omega;
      kBloch[0] = real(kExterior)*SinTheta;
      kBloch[1] = 0.0;
      G->AssembleBEMMatrix(Omega, kBloch, M);
      if (WriteCache)
       { StoreCache( WriteCache );
         WriteCache=0;       
       };
      M->LUFactorize();

      // set plane wave direction 
      nHat[0] = SinTheta;
      nHat[1] = 0.0;
      nHat[2] = CosTheta;
      PW.SetnHat(nHat);

      // solve with E-field perpendicular to plane of incidence 
      E0[0]=0.0;
      E0[1]=1.0;
      E0[2]=0.0;
      PW.SetE0(E0);
      G->AssembleRHSVector(Omega, kBloch, &PW, KN);
      M->LUSolve(KN);
      GetTRFlux(G, &PW, KN, Omega, kBloch, ZAbove, ZBelow, TRFluxPerp);

      // solve with E-field parallel to plane of incidence
      E0[0]=CosTheta;
      E0[1]=0.0;
      E0[2]=-SinTheta;
      PW.SetE0(E0);
      G->AssembleRHSVector(Omega, kBloch, &PW, KN);
      M->LUSolve(KN);
      GetTRFlux(G, &PW, KN, Omega, kBloch, ZAbove, ZBelow, TRFluxPar);
   
      IncFlux = CosTheta/(2.0*ZVAC);

      fprintf(f,"%s %e %e %e %e %e\n",z2s(Omega),Theta*RAD2DEG,
                 TRFluxPerp[0]/IncFlux, TRFluxPerp[1]/IncFlux,
                  TRFluxPar[0]/IncFlux,  TRFluxPar[1]/IncFlux);
      fflush(f);

   }; 
  fclose(f);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Transmission/reflection data written to %s.\n",OutFileName);
  printf("Thank you for your support.\n");

}
