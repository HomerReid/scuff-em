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
#include "libMDInterp.h"
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

  SD->I1D = 0;
  SD->I1DRhoMin=HUGE_VAL;
  SD->I1DRhoMax=0;

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

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
double GetG0Correction(SubstrateData *SD, double z)
{
  for(int n=0; n<SD->NumLayers; n++)
   if ( EqualFloat(z,SD->zLayer[n]) )
    { 
      // break ties by assuming source and observation points
      // lie in whichever region has the lower permittivity
      double EpsA = real( (n==0) ? SD->EpsMedium : SD->EpsLayer[n-1]);
      double EpsB = real(SD->EpsLayer[n]);
      return 2.0*fmin(EpsA, EpsB)/(EpsA + EpsB);
    };
 return 1.0;
}

double GetG0Correction(SubstrateData *SD, double zD, double zS)
{ 
  if (!EqualFloat(zD, zS))
   return 1.0;
  return GetG0Correction(SD, zD);
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
  double EpsMedium  = real(SD->EpsMedium);
  double *zLayer    = SD->zLayer;

  /*--------------------------------------------------------------*/
  /* assemble RHS vector                                          */
  /*--------------------------------------------------------------*/
  double RHS[MAXLAYER];
  for(int m=0; m<NumLayers; m++)
   { 
     double Sign=0.0;
     if (EqualFloat(zLayer[m], zS))
      { double EpsA = real( (m==0) ? SD->EpsMedium : SD->EpsLayer[m-1]);
        double EpsB = real( SD->EpsLayer[m]);
        Sign = (EpsA < EpsB) ? -1.0 : 1.0;
      };

     double ZetaXi[2];
     GetZetaXi(SD, q, zLayer[m], zS, ZetaXi, Sign);
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
       { double Epsmm1   = (m==0) ? EpsMedium : real(EpsLayer[m-1]);
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
 { double RhoMag;   // transverse distance dest-source
   double zD;       // evaluation point z coordinate 
   double zS;       // source point z coordinate 
   SubstrateData *SD;
   FILE *LogFile;
   int NCalls;
 } qIntegrandData;

int qIntegrand(unsigned ndim, const double *u, void *UserData,
               unsigned fdim, double *qIntegral)
{
  (void) fdim; // unused 
  (void) ndim; // unused 

  qIntegral[0]=qIntegral[1]=qIntegral[2]=0.0;

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

  double RhoMag     = qID->RhoMag;
  double zD         = qID->zD;
  double zS         = qID->zS;
  SubstrateData *SD = qID->SD;
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
  double qRho = q*RhoMag;
  double J0, J1;
  if (qRho<1.0e-4) // Bessel-function series expansions for small arguments
   { double qRho2=qRho*qRho;
     J0 = 1.0 - qRho2/4.0;
     J1 = 0.5*qRho*(1.0 - qRho2/8.0);
   }
  else if (qRho>1.0e2) // asymptotic forms for large argument
   {
     double Factor = sqrt(2.0/(M_PI*qRho));
     J0 = Factor * cos(qRho - 0.25*M_PI);
     J1 = Factor * cos(qRho - 0.75*M_PI);
   }
  else
  { cdouble J0J1[2];
    double Workspace[16];
    AmosBessel('J', qRho, 0.0, 2, false, J0J1, Workspace);
    J0=real(J0J1[0]);
    J1=real(J0J1[1]);
  };
  J0*=Jac/(4.0*M_PI);
  J1*=Jac/(4.0*M_PI);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int n=0; n<NumLayers; n++)
   { 
     double ZetaXi[2];
     GetZetaXi(SD, q, zD, zLayer[n], ZetaXi);

     qIntegral[0] +=     J0 * ZetaXi[0] * SigmaTwiddle[n];
     qIntegral[1] += q * J1 * ZetaXi[0] * SigmaTwiddle[n];
     qIntegral[2] += q * J0 * ZetaXi[1] * SigmaTwiddle[n];

     if ( EqualFloat(zD,zLayer[n]) && EqualFloat(zS,zLayer[n]) )
      { cdouble *EpsLayer = SD->EpsLayer;
        cdouble EpsMedium = SD->EpsMedium;
        double EpsA   = real( (n==0) ? EpsMedium: EpsLayer[n-1]);
        double EpsB   = real(EpsLayer[n]);
        double Sign   = (EpsA < EpsB) ? 1.0 : -1.0;
        double Factor = Sign*(EpsA- EpsB) / (EpsA + EpsB);
        qIntegral[0] -= J0*Factor;
        qIntegral[1] -= q*J1*Factor;
        qIntegral[2] -= q*J0*Sign*Factor;
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (qID->LogFile)
   { FILE *LogFile = qID->LogFile;
     static bool init=true;
     if (init)
      { init=false;
        fprintf(LogFile,"\n\n");
        fprintf(LogFile,"#1 2 3 4    q  Rho zD zS\n");
        fprintf(LogFile,"#6 7 8 9 10 J0 J1 qI[0] qI[1] qI[2]\n");
        fprintf(LogFile,"#11 12 13  Sigma,Zeta,Xi (layer 0) \n");
        fprintf(LogFile,"#14 15 16  Sigma,Zeta,Xi (layer 1) \n");
      };
     fprintf(LogFile,"%e %e %e %e ", q, RhoMag, zD, zS);
     fprintf(LogFile,"%e %e ", J0, J1);
     fprintf(LogFile,"%e %e %e ",qIntegral[0],qIntegral[1],qIntegral[2]);
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
/***************************************************************/
/***************************************************************/
void GetqIntegral(SubstrateData *SD, double RhoMag, double zD, double zS, 
                  double qIntegral[3])
{ 
  char *LogFileName = getenv("SSGF_LOGFILE");
  FILE *LogFile     = LogFileName ? fopen(LogFileName, "a") : 0;
  qIntegrandData MyqID, *qID=&MyqID;
  qID->RhoMag  = RhoMag;
  qID->zD      = zD;
  qID->zS      = zS;
  qID->SD      = SD;
  qID->LogFile = LogFile;
  qID->NCalls  = 0;
  double uMin=0.0, uMax=1.0;
  int FDim=3;
  int NDim=1;
  double Error[3];
  int MaxEval   = SD->qMaxEval;
  double AbsTol = SD->qAbsTol;
  double RelTol = SD->qRelTol;
  hcubature(FDim, qIntegrand, (void *)qID, NDim, &uMin, &uMax,
	    MaxEval, AbsTol, RelTol, ERROR_INDIVIDUAL, qIntegral, Error);
  if (LogFile)
   fclose(LogFile);

}

/***************************************************************/
/* Compute the extra contributions to the potential and E-field*/
/* at XDest due to a point charge at XSource in the presence of*/
/* a substrate.                                                */
/*                                                             */
/* 'Extra contributions' include everything but the direct     */
/* contributions of the point charge in vacuum.                */
/***************************************************************/
void GetDeltaPhiESubstrate(SubstrateData *SD,
                           double XD[3], double XS[3],
                           double PhiE[4], double *pG0Correction)
{  
  double Rho[2], ZD=XD[2], ZS=XS[2];
  Rho[0] = XD[0]-XS[0];
  Rho[1] = XD[1]-XS[1];
  double RhoMag = sqrt(Rho[0]*Rho[0] + Rho[1]*Rho[1]);
  double CosTheta = (RhoMag==0.0) ? 1.0 : Rho[0]/RhoMag;
  double SinTheta = (RhoMag==0.0) ? 0.0 : Rho[1]/RhoMag;

  /***************************************************************/
  /* try to get values of q integral using lookup table          */
  /***************************************************************/
  double qIntegral[3];
  bool GotqIntegral=false;
  if (     SD->I1D 
        && EqualFloat(ZD, ZS) && EqualFloat(ZD, SD->I1DZ)
        && SD->I1DRhoMin<=RhoMag && RhoMag<=SD->I1DRhoMax
     )
   {
     GotqIntegral = SD->I1D->Evaluate(RhoMag, qIntegral);
   };
 
  /***************************************************************/
  /*- evaluate q integral to get contributions of surface        */
  /*- charges at dielectric interface                            */
  /***************************************************************/
  if (!GotqIntegral) 
   GetqIntegral(SD, RhoMag, ZD, ZS, qIntegral);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PhiE[0] = qIntegral[0];
  PhiE[1] = CosTheta*qIntegral[1];
  PhiE[2] = SinTheta*qIntegral[1];
  PhiE[3] = qIntegral[2];
  if (pG0Correction)
   *pG0Correction=GetG0Correction(SD,ZD,ZS);

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
  GetDeltaPhiESubstrate(SD, xA, xB, PhiE);
  double Integrand = (SD->WhichIntegral==0) ? PhiE[0] : PhiE[3];
  Result[0] += Weight*Integrand;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double GetSubstratePPI(RWGSurface *SA, int npA,
                       RWGSurface *SB, int npB,
                       SubstrateData *Substrate,
                       double *pG0Correction)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int IDim=1;
  int PPIOrder = Substrate->PPIOrder;
  double Result;
  GetPanelPanelCubature2(SA, npA, SB, npB, SubstratePPIntegrand,
                         (void *)Substrate, IDim, PPIOrder, &Result);

  /***************************************************************/
  /* test for source and destination panel on dielectric intrface*/
  /***************************************************************/
  if (pG0Correction)
   {   
     double G0Correction=1.0;

     RWGPanel *PA = SA->Panels[npA], *PB = SB->Panels[npB];
     double zA = PA->Centroid[2], *zHatA = PA->ZHat;
     double zB = PB->Centroid[2], *zHatB = PB->ZHat;
     if (    EqualFloat(zA, zB)
          && EqualFloat(fabs(zHatA[2]), 1.0)
          && EqualFloat(fabs(zHatB[2]), 1.0)
        ) G0Correction=GetG0Correction(Substrate, zA);

     *pG0Correction=G0Correction;
   };

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
  double Rho2Min=HUGE_VAL, Rho2Max=0.0;
  double DeltazMin=HUGE_VAL, DeltazMax=0.0;
  double z=HUGE_VAL;
  for(int npa=0; npa<Sa->NumPanels; npa++)
   for(int via=0; via<3; via++)
    for(int npb=0; npb<Sb->NumPanels; npb++)
     for(int vib=0; vib<3; vib++)
      { double *VA = Sa->Vertices + 3*(Sa->Panels[npa]->VI[via]);
        double *VB = Sb->Vertices + 3*(Sb->Panels[npb]->VI[vib]);
        double Rho2 = (VA[0]-VB[0])*(VA[0]-VB[0]) 
                     +(VA[1]-VB[1])*(VA[1]-VB[1]); 
        Rho2Min = fmin(Rho2, Rho2Min);
        Rho2Max = fmax(Rho2, Rho2Max);
        double Deltaz = fabs(VA[2]-VB[2]);
        DeltazMin=fmin(Deltaz, DeltazMin);
        DeltazMax=fmax(Deltaz, DeltazMax);
        if (z==HUGE_VAL)
         z=VA[2];
      };
  double RhoMin = sqrt(Rho2Min), RhoMax=sqrt(Rho2Max);
  if (DeltazMax==0.0)
   InitSubstrateAccelerator1D(Substrate, RhoMin, RhoMax,
                              Sa->Panels[0]->Centroid[2]);

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

      double MEPreFactor=1.0;
      if (SurfaceType==PEC)
       {
         Substrate->WhichIntegral = PHIINTEGRAL;
         MEPreFactor = 1.0;
       }
      else // SurfaceType==DIELECTRIC
       { Substrate->WhichIntegral = ENORMALINTEGRAL;
         MEPreFactor = Delta;
       };
      double G0Correction=1.0;
      double MatrixEntry = MEPreFactor*GetSubstratePPI(Sa,npa,Sb,npb,Substrate,&G0Correction);
      if (G0Correction!=1.0)
       M->ScaleEntry(RowOffset + npa, ColOffset + npb, G0Correction);
      M->AddEntry(RowOffset + npa, ColOffset + npb, MatrixEntry);
    };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
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
  GetDeltaPhiESubstrate(SD, XD, XS, DeltaPhiE);
  VecPlusEquals(Integral, Weight, DeltaPhiE, 4);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SSSolver::AddSubstrateContributionToPanelPhiE(int ns, int np, double *XD, double PhiE[4])
{
  SPEIntegrandData MyData, *SPEIData = &MyData;
  SPEIData->Substrate = Substrate;
  SPEIData->XDest     = XD;

  int IDim=4;
  int Order = Substrate->PhiEOrder;
  double DeltaPhiE[4];
  GetPanelCubature2(G->Surfaces[ns], np, SubstratePhiEIntegrand,
                    (void *)SPEIData, IDim, Order, DeltaPhiE);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGPanel *P = G->Surfaces[ns]->Panels[np];
  if ( EqualFloat(P->ZHat[2], 1.0) && EqualFloat(P->Centroid[2], XD[2]) )
   { double G0Correction = GetG0Correction(Substrate, XD[2]);
     if (G0Correction!=1.0)
      VecScale(PhiE, G0Correction, 4);
   };

  VecPlusEquals(PhiE, 1.0, DeltaPhiE, 4);
}

/***************************************************************/
/* entry point for that has the proper prototype for           */
/* passage to the Interp1D() initialization routine.           */
/***************************************************************/
typedef struct fInterpData
 {
   SubstrateData *SD;
   double zD, zS;
 } fInterpData;

void fInterp1D(double Rho, void *UserData, double *fInterp)
{
  fInterpData *fID   = (fInterpData *)UserData;
  SubstrateData *SD  = fID->SD;
  double zD          = fID->zD;
  double zS          = fID->zS;

  double qI[3], dqI[3];
  if (Rho==0.0)
   { double DeltaRho=1.0e-5;
     GetqIntegral(SD, Rho,          zD, zS, qI);
     GetqIntegral(SD, Rho+DeltaRho, zD, zS, dqI);
     VecPlusEquals(dqI, -1.0, qI, 3);
     VecScale(dqI, 1.0/DeltaRho, 3);
   }
  else
   { double DeltaRho=1.0e-5 * fabs(Rho);
     GetqIntegral(SD, Rho+DeltaRho, zD, zS, dqI);
     GetqIntegral(SD, Rho-DeltaRho, zD, zS, qI);
     VecPlusEquals(dqI, -1.0, qI, 3);
     VecScale(dqI, 1.0/(2.0*DeltaRho), 3);
     GetqIntegral(SD, Rho,          zD, zS, qI);
   };
 
  fInterp[0] =  qI[0];
  fInterp[1] = dqI[0];
  fInterp[2] =  qI[1];
  fInterp[3] = dqI[1];
  fInterp[4] =  qI[2];
  fInterp[5] = dqI[2];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitSubstrateAccelerator1D(SubstrateData *SD,
                                double RhoMin, double RhoMax,
                                double z)
{
  if (      SD->I1D 
       &&  (SD->I1DRhoMin <= RhoMin )
       &&  (SD->I1DRhoMax >= RhoMax )
       &&  (SD->I1DZ      == z      )
     )
   { Log(" reusing existing substrate interpolation table");
     return;
   };

  SD->I1DRhoMin = fmin(RhoMin, SD->I1DRhoMin);
  SD->I1DRhoMax = fmax(RhoMax, SD->I1DRhoMax);
  
  if (SD->I1D) delete SD->I1D;

  Log(" (re)allocating substrate I1D(%g,%g,%g)...",SD->I1DRhoMin,SD->I1DRhoMax,z);

  int NGrid=100;
  char *s=getenv("SCUFF_SUBSTRATE_NGRID");
  if (s)
  sscanf(s,"%i",&NGrid);
   
  struct fInterpData MyfID, *fID=&MyfID;
  fID->SD=SD;
  fID->zD=fID->zS=z;
  SD->I1D= new Interp1D(RhoMin, RhoMax, NGrid, 3, fInterp1D,
                        (void *)fID);
  SD->I1DZ      = z;

  Log(" ...done with substrate I1D");

}
