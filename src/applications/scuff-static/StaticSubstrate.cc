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
 * StaticSubstrateGF -- compute the electrostatic green's function
 *                   -- above a single-layer dielectric substrate with
 *                   -- ground plane
 *
 * homer reid        -- 3/2017
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <fenv.h>

#include "libhrutil.h"
#include "libSGJC.h"
#include "libscuff.h"
#include "libSpherical.h"
#include "SSSolver.h"
#include "PanelCubature.h"
#include "StaticSubstrate.h"

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef void (*PPC2Function)(double XA[3], double XB[3], void *UserData,
                             double Weight, double *Integral);

void GetPanelPanelCubature2(RWGSurface *SA, int npA,
                            RWGSurface *SB, int npB,
                            PPC2Function Integrand,
                            void *UserData, int IDim,
                            int Order, double *Integral);

/***************************************************************/
/***************************************************************/
/***************************************************************/
SubstrateData *CreateSubstrateData(const char *FileName, char **pErrMsg)
{
  MatProp *MPMedium=0;
  MatProp **MPLayer=0;
  double *zLayer=0;
  int NumLayers=0;
  double zGP=HUGE_VAL;

  FILE *f=fopen(FileName,"r");
  if (!f)
   { *pErrMsg=vstrdup("could not open file %s",FileName);
     return 0;
   };
  Log("Reading substrate description from file %s.",FileName);

#define MAXSTR 1000
  char Line[MAXSTR];
  int LineNum=0;
  while( fgets(Line,MAXSTR,f) )
   { 
     /*--------------------------------------------------------------*/
     /*- skip blank lines and constants -----------------------------*/
     /*--------------------------------------------------------------*/
     LineNum++;
     int NumTokens;
     char *Tokens[2];
     int Length=strlen(Line);
     if (Length==0) continue;
     Line[(Length-1)]=0; // remove trailing carriage return
     NumTokens=Tokenize(Line, Tokens, 2);
     if ( NumTokens==0 || Tokens[0][0]=='#' )
      continue; 

     /*--------------------------------------------------------------*/
     /*- all lines must be of the form   ----------------------------*/
     /*-   zValue  MaterialName          ----------------------------*/
     /*- or                              ----------------------------*/
     /*-   MEDIUM  MaterialName          ----------------------------*/
     /*- or                              ----------------------------*/
     /*-   zValue  GROUNDPLANE           ----------------------------*/
     /*--------------------------------------------------------------*/
     if ( NumTokens!=2 )
      { *pErrMsg=vstrdup("%s:%i syntax error",FileName,LineNum);
        return 0;
      };

     if ( !strcasecmp(Tokens[0],"MEDIUM") )
      { MPMedium = new MatProp(Tokens[1]);
        if (MPMedium->ErrMsg)
          { *pErrMsg=vstrdup("%s:%i: %s",FileName,LineNum,MPMedium->ErrMsg);
            return 0;
          }
        Log("Setting upper half-space medium to %s.",MPMedium->Name);
        continue;
      };

     double z;
     if ( 1!=sscanf(Tokens[0],"%le",&z) )
      { *pErrMsg=vstrdup("%s:%i bad z-value %s",FileName,LineNum,Tokens[0]);
        return 0;
      };

     if ( !strcasecmp(Tokens[1],"GROUNDPLANE") )
      { zGP = z;
        Log(" Ground plane at z=%e.",zGP);
      }
     else
      { MatProp *MP = new MatProp(Tokens[1]);
        if (MP->ErrMsg)
         { *pErrMsg=vstrdup("%s:%i: %s",FileName,LineNum,MP->ErrMsg);
            return 0;
         };
        if (NumLayers==MAXLAYER)
         { *pErrMsg=vstrdup("%s:%i: too many layers",FileName,LineNum);
            return 0;
         };
        NumLayers++;
        MPLayer=(MatProp **)reallocEC(MPLayer,NumLayers*sizeof(MatProp *));
         zLayer=(double  *)reallocEC(zLayer, NumLayers*sizeof(double));
        MPLayer[NumLayers-1]=MP;
         zLayer[NumLayers-1]=z;
        Log(" Layer #%i: %s at z=%e.",NumLayers,MP->Name,z);
      };
   };
  fclose(f);

  /*--------------------------------------------------------------*/
  /*- sanity check that ground plane lies below all layers       -*/
  /*--------------------------------------------------------------*/
  if (zGP!=HUGE_VAL)
   for(int n=0; n<NumLayers; n++)
    if ( zLayer[n] < zGP )
     { *pErrMsg=vstrdup("%s: ground plane must lie below all dielectric layers",FileName);
       return 0;
     };

  /*--------------------------------------------------------------*/
  /*- prefill buffer of DC permittivities for each layer ---------*/
  /*--------------------------------------------------------------*/
  cdouble EpsMedium = (MPMedium) ? MPMedium->GetEps(0.0) : 1.0;
  cdouble *EpsLayer = (cdouble *)mallocEC(NumLayers*sizeof(cdouble));
  for(int n=0; n<NumLayers; n++)
   EpsLayer[n]=MPLayer[n]->GetEps(0.0);

  int qMaxEval     = 10000;
  double qAbsTol   = 1.0e-12;
  double qRelTol   = 1.0e-8;
  int PPIOrder     = 9;
  int PhiEOrder    = 9;
  //double PPIAbsTol = 1.0e-12;
  //double PPIRelTol = 1.0e-8;
  char *s;
  if ((s=getenv("SCUFF_SUBSTRATE_QMAXEVAL")))
   sscanf(s,"%i",&qMaxEval);
  if ((s=getenv("SCUFF_SUBSTRATE_QABSTOL")))
   sscanf(s,"%le",&qAbsTol);
  if ((s=getenv("SCUFF_SUBSTRATE_QRELTOL")))
   sscanf(s,"%le",&qRelTol);
  if ((s=getenv("SCUFF_SUBSTRATE_PPIORDER")))
   sscanf(s,"%i",&PPIOrder);
  if ((s=getenv("SCUFF_SUBSTRATE_PHIEORDER")))
   sscanf(s,"%i",&PhiEOrder);
  // if ((s=getenv("SCUFF_SUBSTRATE_PPIABSTOL")))
  // sscanf(s,"%le",&PPIAbsTol);
  // if ((s=getenv("SCUFF_SUBSTRATE_PPIRELTOL")))
  // sscanf(s,"%le",&PPIRelTol);

  /*--------------------------------------------------------------*/
  /*- create and return data structure ---------------------------*/
  /*--------------------------------------------------------------*/
  SubstrateData *SD = (SubstrateData *)mallocEC(sizeof(SubstrateData));
  SD->NumLayers     = NumLayers;
  SD->MPMedium      = MPMedium;
  SD->EpsMedium     = EpsMedium;
  SD->MPLayer       = MPLayer;
  SD->EpsLayer      = EpsLayer;
  SD->zLayer        = zLayer;
  SD->zGP           = zGP;
  SD->qMaxEval      = qMaxEval;
  SD->qAbsTol       = qAbsTol;
  SD->qRelTol       = qRelTol;
  SD->PPIOrder      = PPIOrder;
  SD->PhiEOrder     = PhiEOrder;
  //SD->PPIAbsTol     = PPIAbsTol;
  //SD->PPIRelTol     = PPIRelTol;

  *pErrMsg=0;
  return SD;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DestroySubstrateData(SubstrateData *SD)
{
  if (SD->MPMedium)
   delete SD->MPMedium;
  for(int n=0; n<SD->NumLayers; n++)
   delete SD->MPLayer[n];
  free(SD->MPLayer);
  free(SD->EpsLayer);
  free(SD->zLayer);
  free(SD);
}

/*****************************************************************/
/* add the contributions of a charge of strength Q at (xs,ys,zs) */
/* to the potential and E-field of a charge of strength Q at XDest*/
/*****************************************************************/
void AddPhiE0(double XDest[3], double xs, double ys, double zs, double Q, double PhiE[4])
{
  double R[3];
  R[0]=XDest[0]-xs;
  R[1]=XDest[1]-ys;
  R[2]=XDest[2]-zs;
  double r2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double Term = Q/(4.0*M_PI*sqrt(r2));
  PhiE[0] += Term;
  Term/=r2;
  PhiE[1] += R[0]*Term;
  PhiE[2] += R[1]*Term;
  PhiE[3] += R[2]*Term;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetZetaXi(SubstrateData *SD, double q, double zD, double zS, double ZetaXi[2], double Sign=0.0)
{
  double zGP  = SD->zGP;
  bool HaveGP = (zGP!=HUGE_VAL);

  double Term1 = exp(-q*fabs(zD-zS));
  double Term2 = HaveGP ? exp(-q*(zD + zS - 2.0*zGP)) : 0.0;
  if (Sign==0.0)
   Sign  = (zD >= zS) ? 1.0 : -1.0;

  ZetaXi[0] =      Term1 - Term2;
  ZetaXi[1] = Sign*Term1 - Term2;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetSigmaTwiddle(SubstrateData *SD, double zS, double q, double *SigmaTwiddle)
{
  int NumLayers     = SD->NumLayers;
  cdouble *EpsLayer = SD->EpsLayer;
  double *zLayer    = SD->zLayer;

  /*--------------------------------------------------------------*/
  /* assemble RHS vector                                          */
  /*--------------------------------------------------------------*/
  double RHS[MAXLAYER];
  for(int m=0; m<NumLayers; m++)
   { double ZetaXi[2];
     GetZetaXi(SD, q, zLayer[m], zS, ZetaXi);
     RHS[m] = -1.0*ZetaXi[1];
   };

  /*--------------------------------------------------------------*/
  /* assemble M matrix                                           -*/
  /*--------------------------------------------------------------*/
  double M[MAXLAYER*MAXLAYER];
  bool Degenerate=false;
  for(int m=0; m<NumLayers; m++)
   for(int n=0; n<NumLayers; n++)
    { 
      if (m==n)
       { double Epsmm1   = (m==0) ? 1.0 : real(EpsLayer[m-1]);
         double Epsm     = real(EpsLayer[m]);
         double DeltaEps = Epsmm1 - Epsm;
         if (DeltaEps==0.0)
          { Degenerate=true;
            continue; // this is handled below
          };
         double ZetaXiP[2], ZetaXiM[2];
         GetZetaXi(SD, q, zLayer[m], zLayer[m], ZetaXiP, +1.0);
         GetZetaXi(SD, q, zLayer[m], zLayer[m], ZetaXiM, -1.0);
         M[m+m*NumLayers] = (Epsmm1*ZetaXiP[1] - Epsm*ZetaXiM[1])/DeltaEps;
       }
      else
       { double ZetaXi[2];
         GetZetaXi(SD, q, zLayer[m], zLayer[n], ZetaXi);
         M[m+n*NumLayers]=ZetaXi[1];
       };
    };

  /*--------------------------------------------------------------*/
  /*- handle degenerate cases in which Eps_m = Eps_{m-1} by       */
  /*- replacing equation #m with 1*Sigma_m = 0                    */
  /*--------------------------------------------------------------*/
  if (Degenerate)
   for(int m=0; m<NumLayers; m++)
    { double Epsmm1 = (m==0) ? 1.0 : real(EpsLayer[m-1]);
      double Epsm   = real(EpsLayer[m]);
      if (EqualFloat(Epsm, Epsmm1))
       { RHS[m]=0.0;
         for(int n=0; n<NumLayers; n++)
          M[m + n*NumLayers]=(m==n) ? 1.0 : 0.0;
      };
   };

  /*--------------------------------------------------------------*/
  /*- solve the system--------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (NumLayers==1)
   { 
     SigmaTwiddle[0] = RHS[0] / M[0];
   }
  else if (NumLayers==2)
   { double M11=M[0], M21=M[1], M12=M[2], M22=M[3];
     double R1=RHS[0], R2=RHS[1];
     double Denom = M11*M22 - M21*M12;
     SigmaTwiddle[0] = ( M22*R1 - M12*R2 )/Denom;
     SigmaTwiddle[1] = (-M21*R1 + M11*R2 )/Denom;
   }
  else
   { HMatrix MMatrix(NumLayers, NumLayers, LHM_REAL, M);
     HVector RVector(NumLayers, LHM_REAL, RHS);
     MMatrix.LUFactorize();
     MMatrix.LUSolve(&RVector);
     memcpy(SigmaTwiddle, RHS, NumLayers*sizeof(double));
   };
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct qIntegrandData
 { double Rho[2];   // evaluation point - source point (transverse)
   double zD;       // evaluation point z coordinate 
   double zS;       // source point z coordinate 
   SubstrateData *SD;
   FILE *LogFile;
   int NCalls;
 } qIntegrandData;

int qIntegrand(unsigned ndim, const double *u, void *UserData,
               unsigned fdim, double *PhiE)
{
  (void) fdim; // unused 
  (void) ndim; // unused 

  PhiE[0]=PhiE[1]=PhiE[2]=PhiE[3]=0.0;

  /*--------------------------------------------------------------*/
  /* integrand vanishes at q->infinity                            */
  /*--------------------------------------------------------------*/
  if ( EqualFloat(u[0],1.0) )
   return 0;

  double Denom   = 1.0-u[0];
  double q       = u[0] / Denom;
  double Jac     = 1.0/(Denom*Denom);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  qIntegrandData *qID= (qIntegrandData *)UserData;
  qID->NCalls++;

  double *Rho       = qID->Rho;
  double zD         = qID->zD;
  double zS         = qID->zS;
  SubstrateData *SD = qID->SD;
  FILE *LogFile     = qID->LogFile;
  int NumLayers     = SD->NumLayers;
  double *zLayer    = SD->zLayer;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double SigmaTwiddle[MAXLAYER];
  GetSigmaTwiddle(SD, zS, q, SigmaTwiddle);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double RhoMag = sqrt(Rho[0]*Rho[0] + Rho[1]*Rho[1]);
  double qRho = q*RhoMag;
  double J0, J1OverRho;
  if (qRho<1.0e-4) // Bessel-function series expansions for small arguments
   { double qRho2=qRho*qRho;
     J0 = 1.0 - qRho2/4.0;
     J1OverRho = 0.5*q*(1.0 - qRho2/8.0);
   }
  else if (qRho>1.0e2) // asymptotic forms for large argument
   {
     double Factor = sqrt(2.0/(M_PI*qRho));
     J0        = Factor * cos(qRho - 0.25*M_PI);
     J1OverRho = Factor * cos(qRho - 0.75*M_PI) / RhoMag;
   }
  else
  { cdouble J0J1[2];
    double Workspace[16];
    AmosBessel('J', qRho, 0.0, 2, false, J0J1, Workspace);
    J0=real(J0J1[0]);
    J1OverRho= real(J0J1[1])/RhoMag;
  };
  J0*=Jac/(4.0*M_PI);
  J1OverRho*=Jac/(4.0*M_PI);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int n=0; n<NumLayers; n++)
   { 
     double ZetaXi[2];
     GetZetaXi(SD, q, zD, zLayer[n], ZetaXi);

     PhiE[0] +=                 J0 * ZetaXi[0] * SigmaTwiddle[n];
     PhiE[1] += q*Rho[0]*J1OverRho * ZetaXi[0] * SigmaTwiddle[n];
     PhiE[2] += q*Rho[1]*J1OverRho * ZetaXi[0] * SigmaTwiddle[n];
     PhiE[3] +=               q*J0 * ZetaXi[1] * SigmaTwiddle[n];
     if ( EqualFloat(zD,zLayer[n]) && EqualFloat(zS,zLayer[n]) )
      { cdouble *EpsLayer = SD->EpsLayer;
        double Epsnm1  = (n==0) ? 1.0 : real(EpsLayer[n-1]);
        double Epsn    = real(EpsLayer[n]);
        double Sign    = 1.0; // assumes zS > zLayer[n]
        double Factor  = Sign*(Epsn - Epsnm1) / (Epsn + Epsnm1);
        PhiE[0] -= J0*Factor;
        PhiE[1] -= q*Rho[0]*J1OverRho*Factor;
        PhiE[2] -= q*Rho[1]*J1OverRho*Factor;
        PhiE[3] -= Sign*q*J0*Factor;
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (LogFile)
   { fprintf(LogFile,"%e %e %e %e %e ", q, Rho[0], Rho[1], zD, zS);
     fprintf(LogFile,"%e %e ", J0, J1OverRho);
     fprintf(LogFile,"%e %e %e %e ",PhiE[0],PhiE[1],PhiE[2],PhiE[3]);
     for(int n=0; n<NumLayers; n++)
      { double ZetaXi[2];
        GetZetaXi(SD, q, zD, zLayer[n], ZetaXi);
        fprintf(LogFile,"%e %e %e ",SigmaTwiddle[n],ZetaXi[0],ZetaXi[1]);  
      };
     fprintf(LogFile,"\n");
   };

  return 0;
}

/***************************************************************/
/* Compute the corrections to the potential and E-field at     */
/* XDest due to a point charge at XSource in the presence of a */
/* dielectric subtrate of thickness T and relative dielectric  */
/* constant Eps above an infinite ground plane at z=0.         */
/* 'Corrections' means everything but the direct contributions */
/* of the point charge in vacuum.                              */
/***************************************************************/
void GetStaticSubstrateGFCorrection(SubstrateData *SD,
                                    double XD[3], double XS[3],
                                    double PhiE[4],
                                    bool RetainSameLayerContributions)
{
  qIntegrandData MyqID, *qID=&MyqID;
  qID->Rho[0]  = XD[0]-XS[0];
  qID->Rho[1]  = XD[1]-XS[1];
  qID->zD      = XD[2];
  qID->zS      = XS[2];
  qID->SD      = SD;
  qID->NCalls  = 0;
  qID->LogFile = 0;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *s=getenv("SSGF_LOGFILE");
  if (s) 
   { qID->LogFile=fopen(s,"a");
     if (qID->LogFile==0) ErrExit("could not open file %s",s);
     fprintf(qID->LogFile,"#1 q \n");
     fprintf(qID->LogFile,"#2 3 4 5     Rhox Rhoy zDest zSource\n");
     fprintf(qID->LogFile,"#6 7         J0 J1/Rho\n");
     fprintf(qID->LogFile,"#8 9 10 11   PhiE[0..3]\n");
     fprintf(qID->LogFile,"#12 13 14    Sigma,Zeta,Xi (layer 0) \n");
     fprintf(qID->LogFile,"#15 16 17    Sigma,Zeta,Xi (layer 1) \n");
     fprintf(qID->LogFile,"# ...\n");
   };

  /*--------------------------------------------------------------*/
  /*- get contributions of charges at dielectric interfaces      -*/
  /*--------------------------------------------------------------*/
  double uMin=0.0, uMax=1.0;
  int FDim=4;
  int NDim=1;
  double Error[4];
  int MaxEval   = SD->qMaxEval;
  double AbsTol = SD->qAbsTol;
  double RelTol = SD->qRelTol;
  hcubature(FDim, qIntegrand, (void *)qID, NDim, &uMin, &uMax,
	    MaxEval, AbsTol, RelTol, ERROR_INDIVIDUAL, PhiE, Error);
  if (qID->LogFile)
   fclose(qID->LogFile);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (RetainSameLayerContributions && EqualFloat(XD[2],XS[2]) )
   for(int n=0; n<SD->NumLayers; n++)
    if ( EqualFloat(XD[2],SD->zLayer[n]) )
     { 
       double *Rho = qID->Rho;
       double Rho2 = Rho[0]*Rho[0] + Rho[1]*Rho[1];
       if (Rho2==0.0) continue;
       double RhoMag = sqrt(Rho2);
       double CosTheta = Rho[0]/RhoMag, SinTheta=Rho[1]/RhoMag;

       cdouble *EpsLayer = SD->EpsLayer;
       double Epsnm1     = (n==0) ? 1.0 : real(EpsLayer[n-1]);
       double Epsn       = real(EpsLayer[n]);
       double Sign       = 1.0;
       double Factor     = Sign*(Epsn - Epsnm1) / (Epsn + Epsnm1);

       PhiE[0] += Factor/(4.0*M_PI*RhoMag);
       PhiE[1] += CosTheta*Factor/(4.0*M_PI*Rho2);
       PhiE[2] += SinTheta*Factor/(4.0*M_PI*Rho2);
     };

  // contribution of image charge if present
  double zGP = SD->zGP;
  if (zGP != HUGE_VAL)
   AddPhiE0(XD, XS[0], XS[1], 2.0*zGP-XS[2], -1.0, PhiE);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SubstratePPIntegrand(double *xA, double *xB, void *UserData, 
                          double Weight, double *Result)
{
  SubstrateData *SD = (SubstrateData *)UserData;

  double PhiE[4];
  GetStaticSubstrateGFCorrection(SD, xA, xB, PhiE, false);
  double Integrand = (SD->WhichIntegral==0) ? PhiE[0] : PhiE[3];
  Result[0] += Weight * Integrand;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetSubstratePPI(RWGSurface *SA, int npA,
                       RWGSurface *SB, int npB,
                       SubstrateData *Substrate,
                       double *pFactor)
{
  /***************************************************************/
  /* test for source and destination panel on dielectric intrface*/
  /***************************************************************/
  int NumLayers     = Substrate->NumLayers;
  double *zLayer    = Substrate->zLayer;
  cdouble *EpsLayer = Substrate->EpsLayer;
  double EpsMedium  = real(Substrate->EpsMedium);
  RWGPanel *PA = SA->Panels[npA], *PB = SB->Panels[npB];
  double zA = PA->Centroid[2], *zHatA = PA->ZHat;
  double zB = PB->Centroid[2], *zHatB = PB->ZHat;
  if (    EqualFloat(zA, zB)
       && EqualFloat(fabs(zHatA[2]), 1.0)
       && EqualFloat(fabs(zHatB[2]), 1.0)
     )
   {
     for(int nl=0; nl<NumLayers; nl++)
      if ( EqualFloat(zA, zLayer[nl]) )
       { double EpsAbove = (nl==0) ? EpsMedium : real(EpsLayer[nl-1]);
         double EpsBelow = real(EpsLayer[nl]);
         double Numerator = (zA>=zLayer[nl]) ? EpsAbove : EpsBelow;
         *pFactor = 2.0*Numerator/(EpsAbove+EpsBelow);
       };
   };
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int IDim=1;
  int PPIOrder = Substrate->PPIOrder;
  //double RelTol = Substrate->PPIRelTol;
  //double AbsTol = Substrate->PPIAbsTol;
  double Result;
  GetPanelPanelCubature2(SA, npA, SB, npB, SubstratePPIntegrand,
                         (void *)Substrate, IDim, PPIOrder, &Result);

  return Result;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SSSolver::AddSubstrateContributionsToBEMMatrixBlock(int nsa, int nsb, HMatrix *M, int RowOffset, int ColOffset)
{
  Log("Adding substrate to M(%i,%i)...",nsa,nsb);

  RWGSurface *Sa = G->Surfaces[nsa];
  RWGSurface *Sb = G->Surfaces[nsb];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  SurfType SurfaceType;
  double Delta=0.0;
  if (Sa->IsPEC)
   SurfaceType = PEC;
  else
   { SurfaceType = DIELECTRIC;
     double EpsR  = real( G->RegionMPs[ Sa->RegionIndices[0] ] -> GetEps(0.0) );
     double EpsRP = real( G->RegionMPs[ Sa->RegionIndices[1] ] -> GetEps(0.0) );
     Delta = 2.0*(EpsR - EpsRP) / (EpsR + EpsRP);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
#ifdef USE_OPENMP
  int NumThreads = GetNumThreads();
  Log("OpenMP multithreading (%i threads)",NumThreads);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int npa=0; npa<Sa->NumPanels; npa++)
   for(int npb=0; npb<Sb->NumPanels; npb++)
    { 
      if (npb==0) LogPercent(npa, Sa->NumPanels);

      double MEPreFactor=1.0, CorrectionFactor=1.0;
      if (SurfaceType==PEC)
       {
         Substrate->WhichIntegral = PHIINTEGRAL;
         MEPreFactor = 1.0;
       }
      else // SurfaceType==DIELECTRIC
       { Substrate->WhichIntegral = ENORMALINTEGRAL;
         MEPreFactor = Delta;
       };
      double MatrixEntry = MEPreFactor*GetSubstratePPI(Sa,npa,Sb,npb,Substrate,&CorrectionFactor);
      if (CorrectionFactor!=1.0)
       M->ScaleEntry(RowOffset + npa, ColOffset + npb, CorrectionFactor);
      M->AddEntry(RowOffset + npa, ColOffset + npb, MatrixEntry);
    };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/ 
typedef void (*PC2Function)(double X[3], void *UserData, double Weight, double *Integral);

void GetPanelCubature2(RWGSurface *S, int np, PC2Function Integrand,
                       void *UserData, int IDim, int Order,
                       double *Integral);

typedef struct SPEIntegrandData
 {
   SubstrateData *Substrate;
   double *XDest;
 } SPEIntegrandData;

void SubstratePhiEIntegrand(double XS[3], void *UserData,
                            double Weight, double *Integral)
{
  SPEIntegrandData *SPEIData = (SPEIntegrandData *) UserData;
  SubstrateData *SD          = SPEIData->Substrate;
  double *XD                 = SPEIData->XDest;

  double DeltaPhiE[4];
  GetStaticSubstrateGFCorrection(SD, XD, XS, DeltaPhiE, true);
  VecPlusEquals(Integral, Weight, DeltaPhiE, 4);
}

void SSSolver::AddSubstratePhiE(int ns, int np, double *XD, double PhiE[4])
{
  SPEIntegrandData MyData, *SPEIData = &MyData;
  SPEIData->Substrate = Substrate;
  SPEIData->XDest     = XD;

  int IDim=4;
  int Order = Substrate->PhiEOrder;
  double DeltaPhiE[4];
  GetPanelCubature2(G->Surfaces[ns], np, SubstratePhiEIntegrand,
                    (void *)SPEIData, IDim, Order, DeltaPhiE);
  VecPlusEquals(PhiE, 4.0*M_PI, DeltaPhiE, 4);
}
