/* Copyright (C) 2005-2011 M. T. Homer ReidElectricOnly=true;
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <libIncField.h>
#include "scuff-transmission.h"

using namespace scuff;

#define II cdouble (0.0, 1.0)

#define MAXFREQ 10

/***************************************************************/
/* Figure out which of the various regions in the geometry are */
/* the ultimate origin of the incident wave and the ultimate   */
/* destination of the transmitted wave.                        */
/***************************************************************/
bool GetSourceDestRegions(RWGGeometry *G, bool FromAbove,
                          int SourceDestRegions[2])
{
  double XSource[3]={0.0, 0.0, -1.0e6};
  double XDest[3]  ={0.0, 0.0, +1.0e6};
  if (FromAbove)
   { XSource[2]*=-1.0;
     XDest  [2]*=-1.0;
   };

  int SourceRegionIndex=G->GetRegionIndex(XSource);
  int   DestRegionIndex=G->GetRegionIndex(XDest);

  if (    SourceRegionIndex==-1
       || G->RegionMPs[SourceRegionIndex]==0
       || G->RegionMPs[SourceRegionIndex]->IsPEC()
     ) 
   ErrExit("wave cannot emanate from a PEC region");

  Log("Identified source region as region # %i (%s).",
       SourceRegionIndex, G->RegionLabels[SourceRegionIndex]);

  bool DestIsPEC=false;
  if (    DestRegionIndex==-1
       || G->RegionMPs[DestRegionIndex]==0 
       || G->RegionMPs[DestRegionIndex]->IsPEC() 
     ) 
   DestIsPEC=true;

  if (DestIsPEC)
   Log("Destination region is PEC (transmission coefficients identically zero");
  else
   Log("Identified dest region as region   # %i (%s).",
        DestRegionIndex, G->RegionLabels[DestRegionIndex]);

  SourceDestRegions[0]=SourceRegionIndex;
  SourceDestRegions[1]=DestRegionIndex;
  return DestIsPEC;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void WriteFilePreamble(FILE *f)
{
  fprintf(f,"# scuff-transmission " VERSION " run on %s (%s)\n",
             GetHostName(),GetTimeString());
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1:      omega \n");
  fprintf(f,"# 2:      theta (incident angle) (theta=0 --> normal incidence)\n");
  fprintf(f,"# 3:      upward   flux / vacuum plane-wave flux (TE)\n");
  fprintf(f,"# 4:      downward flux / vacuum plane-wave flux (TE)\n");
  fprintf(f,"# 5:      upward   flux / vacuum plane-wave flux (TM)\n");
  fprintf(f,"# 6:      downward flux / vacuum plane-wave flux (TM)\n");
  fprintf(f,"# \n");
  fprintf(f,"# 7,8:    mag, phase a_Upper(TE -> TE)\n");
  fprintf(f,"# 9,10:   mag, phase a_Upper(TE -> TM)\n");
  fprintf(f,"# 11,12   mag, phase a_Upper(TM -> TE)\n");
  fprintf(f,"# 13,14   mag, phase a_Upper(TM -> TM)\n");
  fprintf(f,"# 15,16   mag, phase a_Lower(TE -> TE)\n");
  fprintf(f,"# 17,18   mag, phase a_Lower(TE -> TM)\n");
  fprintf(f,"# 19,20   mag, phase a_Lower(TM -> TE)\n");
  fprintf(f,"# 21,22   mag, phase a_Lower(TM -> TM)\n");
  fprintf(f,"# \n");
  fprintf(f,"# 23,24   mag, phase tTETE\n");
  fprintf(f,"# 25,26   mag, phase tTETM\n");
  fprintf(f,"# 27,28   mag, phase tTMTE\n");
  fprintf(f,"# 29,30   mag, phase tTMTM\n");
  fprintf(f,"# 31,32   mag, phase rTETE\n");
  fprintf(f,"# 33,34   mag, phase rTETM\n");
  fprintf(f,"# 35,36   mag, phase rTMTE\n");
  fprintf(f,"# 37,38   mag, phase rTMTM\n");
  fflush(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void UpdateStandardPWs(cdouble Omega, double Theta, bool FromAbove,
                       MatProp *SourceMP, MatProp *DestMP,
                       PlaneWave *IncidentPW[2],
                       PlaneWave *ReflectedPW[2],
                       PlaneWave *TransmittedPW[2])
{
  cdouble IndexRatio
   = SourceMP->GetRefractiveIndex(Omega) / DestMP->GetRefractiveIndex(Omega);

  double SinTheta  = sin(Theta);
  double CosTheta  = cos(Theta);
  double CosThetaT = real(sqrt(1.0 - IndexRatio*IndexRatio*SinTheta*SinTheta));

  double nHat[3];
  nHat[0] = SinTheta;
  nHat[1] = 0.0;

  cdouble E0TE[3]={0.0, 1.0, 0.0};
  cdouble E0TM[3]={0.0, 0.0, 0.0};
  E0TM[2] = SinTheta;

  // incident
  nHat[2] = (FromAbove ? -1.0 : 1.0) * CosTheta;
  E0TM[0] = -nHat[2]; 

  IncidentPW[POL_TE]->SetE0(E0TE);
  IncidentPW[POL_TE]->SetnHat(nHat);
  IncidentPW[POL_TE]->SetFrequency(Omega);

  IncidentPW[POL_TM]->SetE0(E0TM);
  IncidentPW[POL_TM]->SetnHat(nHat);
  IncidentPW[POL_TM]->SetFrequency(Omega);

  // reflected
  nHat[2] = (FromAbove ? +1.0 : -1.0) * CosTheta;
  E0TM[0] = -nHat[2];

  ReflectedPW[POL_TE]->SetE0(E0TE);
  ReflectedPW[POL_TE]->SetnHat(nHat);
  ReflectedPW[POL_TE]->SetFrequency(Omega);

  ReflectedPW[POL_TM]->SetE0(E0TM);
  ReflectedPW[POL_TM]->SetnHat(nHat);
  ReflectedPW[POL_TM]->SetFrequency(Omega);

  // transmitted
  nHat[2] = (FromAbove ? -1.0 : 1.0) * CosThetaT;
  E0TM[0] = -nHat[2];

  TransmittedPW[POL_TE]->SetE0(E0TE);
  TransmittedPW[POL_TE]->SetnHat(nHat);
  TransmittedPW[POL_TE]->SetFrequency(Omega);

  TransmittedPW[POL_TM]->SetE0(E0TM);
  TransmittedPW[POL_TM]->SetnHat(nHat);
  TransmittedPW[POL_TM]->SetFrequency(Omega);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  InstallHRSignalHandler();
  InitializeLog(argv[0]);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;
  cdouble OmegaVals[MAXFREQ];	int nOmegaVals;
  cdouble LambdaVals[MAXFREQ];	int nLambdaVals;
  char *OmegaFile=0;
  char *LambdaFile=0;
  double Theta=0.0;
  double ThetaMin=0.0;
  double ThetaMax=0.0;
  int ThetaPoints=25;
  double ZAbove=2.0;
  double ZBelow=-1.0;
  int NQPoints=0;
  char *FileBase=0;
  bool FromAbove=false;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",    PA_STRING,  1, 1,       (void *)&GeoFileName,  0,       ".scuffgeo file"},
/**/
     {"Omega",       PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,     &nOmegaVals, "(angular) frequency"},
     {"Lambda",      PA_CDOUBLE, 1, MAXFREQ, (void *)LambdaVals,     &nLambdaVals, "(free-space) wavelength"},
     {"OmegaFile",   PA_STRING,  1, 1,       (void *)&OmegaFile,    0,       "list of (angular) frequencies"},
     {"LambdaFile",  PA_STRING,  1, 1,       (void *)&LambdaFile,   0,       "list of (free-space) wavelengths"},
/**/
     {"Theta",       PA_DOUBLE,  1, 1,       (void *)&Theta,        0,       "incident angle in degrees"},
     {"ThetaMin",    PA_DOUBLE,  1, 1,       (void *)&ThetaMin,     0,       "minimum incident angle in degrees"},
     {"ThetaMax",    PA_DOUBLE,  1, 1,       (void *)&ThetaMax,     0,       "maximum incident angle in degrees"},
     {"ThetaPoints", PA_INT,     1, 1,       (void *)&ThetaPoints,  0,       "number of incident angles"},
/**/
     {"ZAbove",      PA_DOUBLE,  1, 1,       (void *)&ZAbove,       0,       "Z-coordinate of upper integration plane"},
     {"ZBelow",      PA_DOUBLE,  1, 1,       (void *)&ZBelow,       0,       "Z-coordinate of lower integration plane"},
/**/
/**/
     {"NQPoints",    PA_INT,     1, 1,       (void *)&NQPoints,     0,       "number of quadrature points per dimension"},
/**/
     {"FileBase",    PA_STRING,  1, 1,       (void *)&FileBase,     0,       "base file name for output files"},
/**/
     {"FromAbove",   PA_BOOL,    0, 1,       (void *)&FromAbove,    0,       "plane wave impinges from above"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  if (GeoFileName==0)
   OSUsage(argv[0], VERSION, OSArray, "--geometry option is mandatory");
  if (!FileBase)
   FileBase=vstrdup(GetFileBase(GeoFileName));
  SetLogFileName("%s.log",FileBase);
  Log("scuff-transmission" VERSION " running on %s",GetHostName());
  
  /*******************************************************************/
  /*- create the RWGGeometry                                        -*/
  /*******************************************************************/
  RWGGeometry::UseHRWGFunctions=false;
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  if (G->LDim!=2)
   ErrExit("%s: geometry must have two-dimensional lattice periodicity",GeoFileName);

  /*******************************************************************/
  /* determine the indices of the regions from which the plane wave  */
  /* emanates and into which it eventually propagates                */
  /*******************************************************************/
  int SDIndex[2]; // "source, dest index"
  bool DestIsPEC        = GetSourceDestRegions(G, FromAbove, SDIndex);
  int SourceRegionIndex = SDIndex[0];
  int   DestRegionIndex = SDIndex[1];
  int  UpperRegionIndex = (FromAbove) ? SourceRegionIndex : DestRegionIndex;
  int  LowerRegionIndex = (FromAbove) ?   DestRegionIndex : SourceRegionIndex;
  MatProp     *SourceMP = G->RegionMPs[ SourceRegionIndex ];
  MatProp       *DestMP = (DestIsPEC ? 0 : G->RegionMPs[ DestRegionIndex ]);

  /*******************************************************************/
  /* process frequency/wavelength options to construct a list of     */
  /* frequencies at which to run calculations.                       */
  /*******************************************************************/
  HVector *OmegaVector=GetOmegaList(OmegaFile, OmegaVals, nOmegaVals, LambdaFile, LambdaVals, nLambdaVals);
  if ( !OmegaVector || OmegaVector->N==0)
   OSUsage(argv[0], VERSION, OSArray, "you must specify at least one frequency");

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
      OSUsage(argv[0], VERSION, OSArray, "--ThetaMin must not be greater than --ThetaMax");
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
  /*- create the incident field and allocate matrices and vectors    */
  /*******************************************************************/
  cdouble E0[3]={1.0, 0.0, 0.0};
  double nHat[3]={0.0, 0.0, 1.0};
  PlaneWave PW(E0, nHat, G->RegionLabels[SourceRegionIndex]);

  HMatrix *M   = G->AllocateBEMMatrix();
  HVector *KN  = G->AllocateRHSVector();

  PlaneWave *IncidentPW[2];
  PlaneWave *ReflectedPW[2];
  PlaneWave *TransmittedPW[2];
  for(int np=0; np<NUMPOLS; np++)
   { IncidentPW[np] = new PlaneWave(E0, nHat, G->RegionLabels[SourceRegionIndex]);
     ReflectedPW[np] = new PlaneWave(E0, nHat, G->RegionLabels[SourceRegionIndex]);
     TransmittedPW[np] = new PlaneWave(E0, nHat, G->RegionLabels[DestRegionIndex]);
   };

  /*******************************************************************/
  /* set up output files *********************************************/
  /*******************************************************************/
  FILE *f;
  if (!FileBase)
   FileBase=GetFileBase(GeoFileName);
  f=vfopen("%s.transmission","a",FileBase);
  if (!f) ErrExit("could not open file %s",f);
  WriteFilePreamble(f);

  /*******************************************************************/
  /* loop over frequencies and incident angles   *********************/
  /*******************************************************************/
  for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)
   for(int nTheta=0; nTheta<ThetaVector->N; nTheta++)
    { 
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      cdouble Omega = OmegaVector->GetEntry(nOmega);
      Theta = ThetaVector->GetEntryD(nTheta);
      double SinTheta=sin(Theta);
      double CosTheta=cos(Theta);
      Log("Solving the scattering problem at (Omega,Theta)=(%g,%g)",real(Omega),Theta*RAD2DEG);

      UpdateStandardPWs(Omega, Theta, FromAbove, SourceMP, DestMP,
                        IncidentPW, ReflectedPW, TransmittedPW);

      /*--------------------------------------------------------------*/
      /* set bloch wavevector and assemble BEM matrix                 */
      /*--------------------------------------------------------------*/
      cdouble kSource  = SourceMP->GetRefractiveIndex(Omega) * Omega;
      if ( imag(kSource)!=0.0 )
       Warn("complex wavenumber in source region (behavior undefined)");
      double kBloch[2] = {0.0, 0.0};
      kBloch[0] = real(kSource)*SinTheta;
      G->AssembleBEMMatrix(Omega, kBloch, M);
      M->LUFactorize();

      /*--------------------------------------------------------------*/
      /* set plane wave direction and compute polarization vectors    */
      /*--------------------------------------------------------------*/
      nHat[0] = SinTheta;
      nHat[1] = 0.0;
      nHat[2] = FromAbove ? -CosTheta : CosTheta;
      PW.SetnHat(nHat);
      double EpsTE[3]={0.0, 1.0, 0.0}, EpsTM[3], *EpsVectors[2]={EpsTE, EpsTM};
      VecCross(EpsTE, nHat, EpsTM);

      /*--------------------------------------------------------------*/
      /*- loop over the two polarizations of the incident field       */
      /*--------------------------------------------------------------*/
      double UpperFluxRatio[NUMPOLS], LowerFluxRatio[NUMPOLS];
      cdouble UpperAmplitude[NUMPOLS][NUMPOLS], LowerAmplitude[NUMPOLS][NUMPOLS];
      cdouble UpperAmplitude2[NUMPOLS][NUMPOLS], LowerAmplitude2[NUMPOLS][NUMPOLS];
      cdouble UpperAmplitude3[NUMPOLS][NUMPOLS], LowerAmplitude3[NUMPOLS][NUMPOLS];
      cdouble tIntegral[NUMPOLS][NUMPOLS], rIntegral[NUMPOLS][NUMPOLS];
      
      for(int IncPol = POL_TE; IncPol<=POL_TM; IncPol++)
       { 
         E0[0]=EpsVectors[IncPol][0];
         E0[1]=EpsVectors[IncPol][1];
         E0[2]=EpsVectors[IncPol][2];
         PW.SetE0(E0);
         G->AssembleRHSVector(Omega, kBloch, &PW, KN);
         M->LUSolve(KN);

         double Flux[NUMREGIONS];
         GetFlux(G, &PW, KN, Omega, kBloch, NQPoints, 
                 ZAbove, ZBelow, FromAbove,
                 ReflectedPW, TransmittedPW, 
                 Flux, tIntegral[IncPol], rIntegral[IncPol]);
         UpperFluxRatio[IncPol] = Flux[REGION_UPPER];
         LowerFluxRatio[IncPol] = Flux[REGION_LOWER];

         cdouble aTETM[2];
         GetPlaneWaveAmplitudes(G, KN, Omega, kBloch, UpperRegionIndex, true, aTETM, true);
         UpperAmplitude[IncPol][POL_TE] = aTETM[POL_TE];
         UpperAmplitude[IncPol][POL_TM] = aTETM[POL_TM];

         GetPlaneWaveAmplitudes(G, KN, Omega, kBloch, LowerRegionIndex, false, aTETM, true);
         LowerAmplitude[IncPol][POL_TE] = aTETM[POL_TE];
         LowerAmplitude[IncPol][POL_TM] = aTETM[POL_TM];

         for(int DestPol = POL_TE; DestPol<=POL_TM; DestPol++)
          { 
            PlaneWave *UpperPW = FromAbove ? ReflectedPW[DestPol]   : TransmittedPW[DestPol];
            PlaneWave *LowerPW = FromAbove ? TransmittedPW[DestPol] : ReflectedPW[DestPol];
            cdouble Terms[2];
            UpperAmplitude2[IncPol][DestPol] = G->GetPlaneWaveAmplitude(Omega, kBloch, KN, UpperPW, Terms);
            UpperAmplitude3[IncPol][DestPol] = Terms[0];
            LowerAmplitude2[IncPol][DestPol] = G->GetPlaneWaveAmplitude(Omega, kBloch, KN, LowerPW, Terms);
            LowerAmplitude3[IncPol][DestPol] = Terms[0];
          }
       }

      // write results to file
      fprintf(f,"%s %e ", z2s(Omega), Theta*RAD2DEG);

      fprintf(f,"%e %e ", UpperFluxRatio[POL_TE], LowerFluxRatio[POL_TE]);
      fprintf(f,"%e %e ", UpperFluxRatio[POL_TM], LowerFluxRatio[POL_TM]);

      fprintf(f,"%e %e ", abs(UpperAmplitude[POL_TE][POL_TE]), arg(UpperAmplitude[POL_TE][POL_TE]));
      fprintf(f,"%e %e ", abs(UpperAmplitude[POL_TE][POL_TM]), arg(UpperAmplitude[POL_TE][POL_TM]));
      fprintf(f,"%e %e ", abs(UpperAmplitude[POL_TM][POL_TE]), arg(UpperAmplitude[POL_TM][POL_TE]));
      fprintf(f,"%e %e ", abs(UpperAmplitude[POL_TM][POL_TM]), arg(UpperAmplitude[POL_TM][POL_TM]));
      fprintf(f,"%e %e ", abs(LowerAmplitude[POL_TE][POL_TE]), arg(LowerAmplitude[POL_TE][POL_TE]));
      fprintf(f,"%e %e ", abs(LowerAmplitude[POL_TE][POL_TM]), arg(LowerAmplitude[POL_TE][POL_TM]));
      fprintf(f,"%e %e ", abs(LowerAmplitude[POL_TM][POL_TE]), arg(LowerAmplitude[POL_TM][POL_TE]));
      fprintf(f,"%e %e ", abs(LowerAmplitude[POL_TM][POL_TM]), arg(LowerAmplitude[POL_TM][POL_TM]));

      fprintf(f,"%e %e ", abs(tIntegral[POL_TE][POL_TE]), arg(tIntegral[POL_TE][POL_TE]));
      fprintf(f,"%e %e ", abs(tIntegral[POL_TE][POL_TM]), arg(tIntegral[POL_TE][POL_TM]));
      fprintf(f,"%e %e ", abs(tIntegral[POL_TM][POL_TE]), arg(tIntegral[POL_TM][POL_TE]));
      fprintf(f,"%e %e ", abs(tIntegral[POL_TM][POL_TM]), arg(tIntegral[POL_TM][POL_TM]));
      fprintf(f,"%e %e ", abs(rIntegral[POL_TE][POL_TE]), arg(rIntegral[POL_TE][POL_TE]));
      fprintf(f,"%e %e ", abs(rIntegral[POL_TE][POL_TM]), arg(rIntegral[POL_TE][POL_TM]));
      fprintf(f,"%e %e ", abs(rIntegral[POL_TM][POL_TE]), arg(rIntegral[POL_TM][POL_TE]));
      fprintf(f,"%e %e ", abs(rIntegral[POL_TM][POL_TM]), arg(rIntegral[POL_TM][POL_TM]));

      fprintf(f,"%e %e ", abs(UpperAmplitude2[POL_TE][POL_TE]), arg(UpperAmplitude2[POL_TE][POL_TE]));
      fprintf(f,"%e %e ", abs(UpperAmplitude2[POL_TE][POL_TM]), arg(UpperAmplitude2[POL_TE][POL_TM]));
      fprintf(f,"%e %e ", abs(UpperAmplitude2[POL_TM][POL_TE]), arg(UpperAmplitude2[POL_TM][POL_TE]));
      fprintf(f,"%e %e ", abs(UpperAmplitude2[POL_TM][POL_TM]), arg(UpperAmplitude2[POL_TM][POL_TM]));
      fprintf(f,"%e %e ", abs(LowerAmplitude2[POL_TE][POL_TE]), arg(LowerAmplitude2[POL_TE][POL_TE]));
      fprintf(f,"%e %e ", abs(LowerAmplitude2[POL_TE][POL_TM]), arg(LowerAmplitude2[POL_TE][POL_TM]));
      fprintf(f,"%e %e ", abs(LowerAmplitude2[POL_TM][POL_TE]), arg(LowerAmplitude2[POL_TM][POL_TE]));
      fprintf(f,"%e %e ", abs(LowerAmplitude2[POL_TM][POL_TM]), arg(LowerAmplitude2[POL_TM][POL_TM]));

      fprintf(f,"%e %e ", abs(UpperAmplitude3[POL_TE][POL_TE]), arg(UpperAmplitude3[POL_TE][POL_TE]));
      fprintf(f,"%e %e ", abs(UpperAmplitude3[POL_TE][POL_TM]), arg(UpperAmplitude3[POL_TE][POL_TM]));
      fprintf(f,"%e %e ", abs(UpperAmplitude3[POL_TM][POL_TE]), arg(UpperAmplitude3[POL_TM][POL_TE]));
      fprintf(f,"%e %e ", abs(UpperAmplitude3[POL_TM][POL_TM]), arg(UpperAmplitude3[POL_TM][POL_TM]));
      fprintf(f,"%e %e ", abs(LowerAmplitude3[POL_TE][POL_TE]), arg(LowerAmplitude3[POL_TE][POL_TE]));
      fprintf(f,"%e %e ", abs(LowerAmplitude3[POL_TE][POL_TM]), arg(LowerAmplitude3[POL_TE][POL_TM]));
      fprintf(f,"%e %e ", abs(LowerAmplitude3[POL_TM][POL_TE]), arg(LowerAmplitude3[POL_TM][POL_TE]));
      fprintf(f,"%e %e ", abs(LowerAmplitude3[POL_TM][POL_TM]), arg(LowerAmplitude3[POL_TM][POL_TM]));

      fprintf(f,"\n");
      fflush(f);

   }; 
  fclose(f);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Transmission/reflection data written to %s.transmission\n",FileBase);
  printf("Thank you for your support.\n");

}
