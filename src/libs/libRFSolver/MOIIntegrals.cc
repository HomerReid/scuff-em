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
 * homer reid    -- 11/2005 -- 3/2018
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "PanelCubature.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "RFSolver.h"

namespace scuff {
void GetReducedPotentials_Nearby(RWGSurface *S, const int ne,
                                 const double X0[3],  const cdouble k,
                                 cdouble *p, cdouble a[3],
                                 cdouble dp[3], cdouble da[3][3],
                                 cdouble ddp[3][3], cdouble dcurla[3][3],
                                 bool *IncludeTerm=0);

using namespace scuff;

#define II cdouble(0.0,1.0)

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct MOIIntegrandData
 {
   LayeredSubstrate *S;
   cdouble Omega;
   cdouble Eps;
   bool Subtract;

   // used for GetMOI1BFFields
   double *XDest;

   // used for GetMOIMatrixElement and GetMOIMatrixElement_NC
   bool SeparateTerms;

   // used for GetMOIMatrixElement_NC
   int nsInner, neInner, InnerOrder;
   RWGGeometry *G;

   int CubaturePoints, qPoints;
   
 } MOIIntegrandData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void MOIPFIntegrand(double xS[3], double b[3], double Divb,
                    void *UserData, double Weight, double *Integral)
{
   MOIIntegrandData *Data = (MOIIntegrandData *)UserData;
   LayeredSubstrate *S    = Data->S;
   cdouble Omega          = Data->Omega;
   bool Subtract          = Data->Subtract;
   double *xD             = Data->XDest;

   if ( fabs(xS[2]) > 1.0e-12 )
    ErrExit("%s:%i: internal error",__FILE__,__LINE__);

   double Rhox  = xD[0]-xS[0], Rhoy=xD[1]-xS[1];
   double Rho   = sqrt(Rhox*Rhox + Rhoy*Rhoy);
   double CT    = (Rho==0.0) ? 0.0 : Rhox/Rho;
   double ST    = (Rho==0.0) ? 0.0 : Rhoy/Rho;
   double zDest = xD[2];

   ScalarGFOptions Options;
   InitScalarGFOptions(&Options);
   Options.PPIsOnly=false;
   Options.Subtract=Subtract;
   Options.RetainSingularTerms=!Subtract;

   cdouble V[NUMSGFS_MOI];
   int qPoints = S->GetScalarGFs_MOI(Omega, Rho, zDest, V, &Options);

   cdouble *PFCorrections=(cdouble *)Integral;
   cdouble iw=II*Omega;
   Weight*=ZVAC;
   PFCorrections[_PF_AX]    += Weight*V[_SGF_APAR ]*b[0];
   PFCorrections[_PF_AY]    += Weight*V[_SGF_APAR ]*b[1];
   PFCorrections[_PF_AZ]    += Weight*V[_SGF_AZ   ]*(b[0]*CT + b[1]*ST);
   PFCorrections[_PF_PHI]   += Weight*V[_SGF_PHI  ]*Divb/iw;
   PFCorrections[_PF_DXPHI] += Weight*V[_SGF_DRPHI]*Divb*CT/iw;
   PFCorrections[_PF_DYPHI] += Weight*V[_SGF_DRPHI]*Divb*ST/iw;
   PFCorrections[_PF_DZPHI] += Weight*V[_SGF_DZPHI]*Divb/iw;

   Data->CubaturePoints++;
   Data->qPoints+=qPoints;
}

/*--------------------------------------------------------------*/
/* this routine computes the *bare* (substrate-independent)     */
/* contributions to the potential and fields, in order to add   */
/* them back in after the corresponding terms have been         */
/* subtracted from the Sommerfeld integral for the full         */
/* substrate GF                                                 */
/*--------------------------------------------------------------*/
typedef struct RPIData
 { cdouble k;
   double *XDest;
 } RPIData;

void RPIntegrand(double XSource[3], double b[3], double Divb,
                 void *UserData, double Weight, double *Integral)
{
  RPIData *Data = (RPIData *)UserData;
  cdouble k     = Data->k;
  double *XDest = Data->XDest;

  double R[3]; VecSub(XDest, XSource, R);
  double r2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  if (r2==0.0) return;
  double r=sqrt(r2);
  cdouble ikr=II*k*r;
  cdouble G0 = Weight*exp(ikr)/(4.0*M_PI*r), dG0 = G0*(ikr-1.0)/r2;

  cdouble *apdp = (cdouble *)Integral;
  apdp[_PF_AX]    += G0       * b[0];
  apdp[_PF_AY]    += G0       * b[1];
  apdp[_PF_AX]    += G0       * b[2];
  apdp[_PF_PHI]   += G0       * Divb;
  apdp[_PF_DXPHI] += R[0]*dG0 * Divb;
  apdp[_PF_DYPHI] += R[1]*dG0 * Divb;
  apdp[_PF_DZPHI] += R[2]*dG0 * Divb;
}
  

void GetReducedPotentials(RWGGeometry *G, int ns, int ne, cdouble k,
                          double *XDest,
                          cdouble *p, cdouble a[3], cdouble dp[3])
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  static double rRelOuterThreshold=4.0;
  static double rRelInnerThreshold=1.0;
  static int LowOrder=7;
  static int HighOrder=20;

  static bool Init=true;
  if (Init)
   { Init=false;
     char *s1=getenv("SCUFF_RREL_OUTER_THRESHOLD");
     char *s2=getenv("SCUFF_RREL_INNER_THRESHOLD");
     char *s3=getenv("SCUFF_LOWORDER");
     char *s4=getenv("SCUFF_HIGHORDER");
     if (s1) sscanf(s1,"%le",&rRelOuterThreshold);
     if (s2) sscanf(s2,"%le",&rRelInnerThreshold);
     if (s3) sscanf(s3,"%i",&LowOrder);
     if (s4) sscanf(s4,"%i",&HighOrder);
     if (s1||s2||s3||s4)
      Log("({O,I}rRelThreshold | LowOrder | HighOrder)=(%e,%e,%i,%i)",
          rRelOuterThreshold,rRelInnerThreshold,LowOrder,HighOrder);
   }
   
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGSurface *S = G->Surfaces[ns];
  RWGEdge *E = S->GetEdgeByIndex(ne);
  double rRel = VecDistance(XDest, E->Centroid) / E->Radius;
  if (rRel>=rRelInnerThreshold)
   { 
     int Order = rRel >= rRelOuterThreshold ? LowOrder : HighOrder;
     RPIData Data;
     Data.k     = k;
     Data.XDest = XDest;
     cdouble apdp[7];
     int zfdim=7, fdim=2*zfdim;
     GetBFCubature2(G, ns, ne, RPIntegrand, (void *)&Data,
                    fdim, Order, (double *)apdp);
     *p=apdp[_PF_PHI];
     for(int i=0; i<3; i++)
      {  a[i] = apdp[_PF_AX+i];
        dp[i] = apdp[_PF_DXPHI+i];
      }
   }
  else
   { cdouble da[3][3], ddp[3][3], dcurla[3][3];
     GetReducedPotentials_Nearby(S, ne, XDest, k, p, a, dp,
                                 da, ddp, dcurla);
   }

}

/***************************************************************/
/* Compute the MOI potential and fields due to a single basis  */
/* function populated with unit strength.                      */
/* PFVector = 'potential and field vector'                     */
/*                                                             */
/* PFVector[0,1,2] = Ax, Ay, Az                                */
/* PFVector[3]     = Phi                                       */
/* PFVector[4,5,6] = dxPhi, dyPhi, dzPhi                       */
/*                                                             */
/* if SingularTerms != NULL, it should be a buffer of length 7 */
/* which is filled in on return with the contributions of just */
/* the most singular terms in the GFs (if Subtract==true).     */
/***************************************************************/
int Get1BFMOIFields(RWGGeometry *G, LayeredSubstrate *Substrate,
                    int ns, int ne, cdouble Omega, double *XDest,
                    cdouble *PFVector,
                    int Order, bool Subtract, cdouble *SingularTerms)
{
  /***************************************************************/
  /* get substrate contributions                                 */
  /*  (minus the most singular terms if Subtract==true)          */
  /***************************************************************/
  Substrate->UpdateCachedEpsMu(Omega);
  cdouble Eps = Substrate->EpsLayer[1];

  MOIIntegrandData MOIData;
  MOIData.S              = Substrate;
  MOIData.Omega          = Omega;
  MOIData.Eps            = Eps;
  MOIData.Subtract       = Subtract;
  MOIData.XDest          = XDest;
  MOIData.qPoints        = 0;
  MOIData.CubaturePoints = 0;

  if (Order==-1) Order = 9;

  int zfdim=NPFC, fdim=2*zfdim;
  GetBFCubature2(G, ns, ne, MOIPFIntegrand, (void *)&MOIData, fdim, Order, (double *)PFVector);

  /***************************************************************/
  /* add contributions of most singular terms in GF              */
  /***************************************************************/
  if (SingularTerms) memset(SingularTerms, 0, NPFC*sizeof(cdouble));
  if ( Subtract && (XDest[2] > -1.0e-12) )
   { 
     cdouble p=0.0, a[3]={0.0,0.0,0.0}, dp[3]={0.0,0.0,0.0};
     GetReducedPotentials(G, ns, ne, Omega, XDest, &p, a, dp);

     cdouble STBuffer[NPFC];
     if (!SingularTerms) SingularTerms=STBuffer;
     cdouble iw=II*Omega;
     cdouble Eta = (Eps-1.0)/(Eps+1.0);
     SingularTerms[_PF_AX]    = ZVAC*a[0];
     SingularTerms[_PF_AY]    = ZVAC*a[1];
     SingularTerms[_PF_AZ]    = 0.0; // we don't subtract the q integral for Az
     SingularTerms[_PF_PHI]   = (1.0-Eta)*ZVAC*p/iw;
     SingularTerms[_PF_DXPHI] = (1.0-Eta)*ZVAC*dp[0]/iw;
     SingularTerms[_PF_DYPHI] = (1.0-Eta)*ZVAC*dp[1]/iw;
     SingularTerms[_PF_DZPHI] = (1.0-Eta)*ZVAC*dp[2]/iw;

     VecPlusEquals(PFVector, 1.0, SingularTerms, NPFC);
   }

  return MOIData.qPoints;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void MOIMEIntegrand(double xA[3], double bA[3], double DivbA,
                    double xB[3], double bB[3], double DivbB,
                    void *UserData, double Weight, double *Integral)
{
   MOIIntegrandData *Data = (MOIIntegrandData *)UserData;
   LayeredSubstrate *S    = Data->S;
   cdouble Omega          = Data->Omega;
   cdouble Eps            = Data->Eps;
   bool Subtract          = Data->Subtract;
   bool SeparateTerms     = Data->SeparateTerms;

   double Rhox = xA[0]-xB[0], Rhoy=xA[1]-xB[1];
   double Rho = sqrt(Rhox*Rhox + Rhoy*Rhoy);
   if ( fabs(xA[2])>1.0e-12 || fabs(xB[2])>1.0e-12 )
    ErrExit("%s:%i: internal error",__FILE__,__LINE__);

   ScalarGFOptions Options;
   InitScalarGFOptions(&Options);
   Options.PPIsOnly = true;
   Options.Subtract = Subtract;
   Options.RetainSingularTerms = !Subtract;
   Options.MaxTerms = S->qMaxEval;
   Options.RelTol   = S->qRelTol;
   Options.AbsTol   = S->qAbsTol;

   cdouble V[2];
   int qPoints = S->GetScalarGFs_MOI(Omega, Rho, 0.0, V, &Options);

   Data->CubaturePoints += 1;
   Data->qPoints += qPoints;

   cdouble iw = II*Omega, ooiw=1.0/iw;
   double bDb = Weight*VecDot(bA, bB), DivbDivb = Weight*DivbA*DivbB;
   cdouble ATerm   =   iw*bDb      * V[_SGF_APAR];
   cdouble PhiTerm = ooiw*DivbDivb * V[_SGF_PHI];
   
   cdouble *zIntegral = (cdouble *)Integral;
   if (!SeparateTerms)
    { 
      zIntegral[0] += ATerm + PhiTerm;
    }
   else
    { zIntegral[0] += ATerm;
      zIntegral[1] += PhiTerm;
      if (!Subtract) return;

      cdouble DeltaVS[2], DeltaVNS[2];
      double h=S->zInterface[0] - S->zGP;
      Options.RetainSingularTerms=true;
      Options.MaxTerms = 1;
      GetSGFCorrection_MOI(Omega, Rho, 0.0, Eps, h, DeltaVS, &Options);
      Options.MaxTerms = S->qMaxEval;
      GetSGFCorrection_MOI(Omega, Rho, 0.0, Eps, h, DeltaVNS, &Options);
      VecPlusEquals(DeltaVNS, -1.0, DeltaVS, 2);

      zIntegral[2] += iw*bDb        * DeltaVNS[_SGF_APAR];
      zIntegral[3] += ooiw*DivbDivb * DeltaVNS[_SGF_PHI];
      zIntegral[4] += iw*bDb        * DeltaVS[_SGF_APAR];
      zIntegral[5] += ooiw*DivbDivb * DeltaVS[_SGF_PHI];
    }
}

/***************************************************************/
/* if Terms is nonzero, its entries are filled in as follows:  */
/*  Terms[0] = ATerm   (subtracted q integral)                 */
/*  Terms[2] = ATerm   (correction minus singular part)        */
/*  Terms[4] = ATerm   (singular part of correction, cubature) */
/*  Terms[6] = ATerm   (singular part of correction, exact)    */
/*  Terms[1] = PhiTerm (subtracted q integral)                 */
/*  Terms[3] = PhiTerm (correction minus singular part)        */
/*  Terms[5] = PhiTerm (singular part of correction, cubature) */
/*  Terms[7] = PhiTerm (singular part of correction, exact)    */
/***************************************************************/
int GetMOIMatrixElement(RWGGeometry *G, LayeredSubstrate *Substrate,
                        int nsa, int nea, int nsb, int neb,
                        cdouble Omega, cdouble *ME,
                        int Order, bool Subtract, cdouble *Terms)
{
  /***************************************************************/
  /* if cubature order unspecified, autodetermine it based on    */
  /* distance between panels                                     */
  /***************************************************************/
  if (Order==-1) 
   { double rRel;
     int ncv=AssessBFPair(G->Surfaces[nsa], nea, G->Surfaces[nsb], neb, &rRel);
     Order = ( (ncv>0 || rRel <=1.0) ? 9 : 4);
   }

  /***************************************************************/
  /* get substrate contributions *********************************/
  /***************************************************************/
  Substrate->UpdateCachedEpsMu(Omega);
  cdouble Eps = Substrate->EpsLayer[1];
  
  bool SeparateTerms = (Terms!=0);

  MOIIntegrandData MOIData;
  MOIData.S              = Substrate;
  MOIData.Omega          = Omega;
  MOIData.Eps            = Eps;
  MOIData.Subtract       = Subtract;
  MOIData.SeparateTerms  = SeparateTerms;
  MOIData.qPoints        = 0;
  MOIData.CubaturePoints = 0;

  int zfdim = SeparateTerms ? 6 : 1;
  int fdim = 2*zfdim;
  GetBFBFCubature2(G, nsa, nea, nsb, neb,
                   MOIMEIntegrand, (void *)&MOIData,
                   fdim, Order, (double *)(SeparateTerms ? Terms : ME));
  if (SeparateTerms)
   ME[0]=Terms[0]+Terms[1];

  /***************************************************************/
  /* add 'bare' (substrate-independent) EEIs                     */
  /***************************************************************/
  if (Subtract)
   { 
     GetGCMEArgStruct GCMEArgs, *Args=&GCMEArgs;
     InitGetGCMEArgs(Args);
     Args->nsa = nsa;
     Args->nsb = nsb;
     Args->NumRegions = 1;
     Args->k[0] = Omega;
     Args->NeedGC=true;
     Args->FIBBICache = (nsa==nsb) ? G->FIBBICaches[nsa] : 0;
     cdouble GabArray[2][NUMGCMES];
     cdouble ikCabArray[2][NUMGCMES];
     cdouble GPhiTerm[2];
     GetGCMatrixElements(G, Args, nea, neb, GabArray, ikCabArray, GPhiTerm);
     GabArray[0][GCME_GC]*=II*Omega;
     GPhiTerm[0]*=II*Omega;
     cdouble Eta = (Eps-1.0)/(Eps+1.0);
     if (Terms)
      { Terms[6] = GabArray[0][GCME_GC] - GPhiTerm[0];
        Terms[7] = (1.0-Eta)*GPhiTerm[0];
      }
     ME[0] += GabArray[0][GCME_GC] - Eta*GPhiTerm[0];
   }

  return MOIData.qPoints;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void MOIMEIntegrand_NC(double XDest[3], double b[3], double Divb,
                       void *UserData, double Weight, double *Integral)
{
   MOIIntegrandData *OuterData = (MOIIntegrandData *)UserData;
   LayeredSubstrate *S          = OuterData->S;
   cdouble Omega                = OuterData->Omega;
   cdouble Eps                  = OuterData->Eps;
   bool Subtract                = OuterData->Subtract;
   bool SeparateTerms           = OuterData->SeparateTerms;

   int nsInner                  = OuterData->nsInner;
   int neInner                  = OuterData->neInner;
   int InnerOrder               = OuterData->InnerOrder;
   RWGGeometry *G               = OuterData->G;

   cdouble PFVector[NPFC], STBuffer[NPFC], *SingularTerms = (SeparateTerms ? STBuffer : 0);
   int qPoints = Get1BFMOIFields(G, S, nsInner, neInner, Omega, XDest, PFVector, InnerOrder, Subtract, SingularTerms);

   cdouble *METerms = (cdouble *)Integral;
   cdouble iw=II*Omega;
   Weight/=ZVAC;
   if (SeparateTerms)
    { cdouble NSTerms[NPFC];
      VecSub(PFVector, SingularTerms, NSTerms, NPFC);
      METerms[0] += Weight*iw*(b[0]*NSTerms[_PF_AX] + b[1]*NSTerms[_PF_AY]);
      METerms[1] += Weight*Divb*NSTerms[_PF_PHI];
      METerms[2] += Weight*iw*(b[0]*SingularTerms[_PF_AX] + b[1]*SingularTerms[_PF_AY]);
      METerms[3] += Weight*Divb*SingularTerms[_PF_PHI];
    }
   else
    {
      METerms[0]+= Weight * ( iw*(b[0]*PFVector[_PF_AX] + b[1]*PFVector[_PF_AY])
                              +Divb*PFVector[_PF_PHI]
                            );
    }

   OuterData->qPoints        += qPoints;
   OuterData->CubaturePoints += 1;
}

/***************************************************************/
/* Get substrate EEIs via nested cubature, mostly just as a    */
/* sanity check of GetMOIMatrixElement                         */
/* if Terms is nonzero, its entries are filled in as follows:  */
/*  Terms[0] = ATerm   (subtracted q integral)                 */
/*  Terms[2] = ATerm   (singular part of correction, exact)    */
/*  Terms[1] = PhiTerm (subtracted q integral)                 */
/*  Terms[3] = PhiTerm (singular part of correction, exact)    */
/***************************************************************/
int GetMOIMatrixElement_NC(RWGGeometry *G, LayeredSubstrate *Substrate,
                           int nsa, int nea, int nsb, int neb,
                           cdouble Omega, cdouble *MEE,
                           int InnerOrder, int OuterOrder, bool Subtract, cdouble *Terms)
{
  /***************************************************************/
  /* get substrate contributions *********************************/
  /***************************************************************/
  Substrate->UpdateCachedEpsMu(Omega);
  cdouble Eps = Substrate->EpsLayer[1];

  if (InnerOrder==-1) InnerOrder=9;
  if (OuterOrder==-1) OuterOrder=9;
  bool SeparateTerms = (Terms!=0);

  MOIIntegrandData MOIData;
  MOIData.S              = Substrate;
  MOIData.Omega          = Omega;
  MOIData.Eps            = Eps;
  MOIData.Subtract       = Subtract;
  MOIData.SeparateTerms  = SeparateTerms;
  MOIData.nsInner        = nsb;
  MOIData.neInner        = neb;
  MOIData.InnerOrder     = InnerOrder;
  MOIData.G              = G;
  MOIData.qPoints        = 0;
  MOIData.CubaturePoints = 0;
 
  int zfdim = (SeparateTerms ? 4 : 1);
  int fdim = 2*zfdim;
  cdouble TermBuffer[4];
  GetBFCubature2(G, nsa, nea, MOIMEIntegrand_NC, (void *)&MOIData,
                 fdim, OuterOrder, (double *)(Terms ? Terms : TermBuffer));
  if (Terms)
   MEE[0] = Terms[0]+Terms[1]+Terms[2]+Terms[3];

  return MOIData.qPoints;

}

/***************************************************************/
/* shortest, longest distances between points on given surfaces*/
/***************************************************************/
void GetRhoMinMax(RWGGeometry *G, int nsa, int nsb, double RhoMinMax[2])
{ 
  double RMax[3], RMin[3];
  RWGSurface *Sa = G->Surfaces[nsa];
  RWGSurface *Sb = G->Surfaces[nsb];
  VecSub(Sa->RMax, Sb->RMin, RMax);
  VecSub(Sa->RMin, Sb->RMax, RMin);
  RhoMinMax[0] = sqrt(RMin[0]*RMin[0] + RMin[1]*RMin[1]);
  RhoMinMax[1] = sqrt(RMax[0]*RMax[0] + RMax[1]*RMax[1]);
}

/***************************************************************/
/* shortest, longest distances between points on any surfaces  */
/***************************************************************/
void GetRhoMinMax(RWGGeometry *G, double RhoMinMax[2])
{  
  RhoMinMax[0]=1.0e89;
  RhoMinMax[1]=-1.0e89;
  for(int nsa=0; nsa<G->NumSurfaces; nsa++)
   for(int nsb=nsa; nsb<G->NumSurfaces; nsb++)
    { double ThisRhoMinMax[2];
      GetRhoMinMax(G, nsa, nsb, ThisRhoMinMax);
      RhoMinMax[0]=0.0; //fmin(RhoMinMax[0],ThisRhoMinMax[0]);
      RhoMinMax[1]=fmax(RhoMinMax[1],ThisRhoMinMax[1]);
    }
}

/***************************************************************/
/* shortest, longest distances between points in XMatrix and   */
/* points on any surface. Note XMatrix is 3 x NX.              */
/***************************************************************/
void GetRzMinMax(RWGGeometry *G, HMatrix *XMatrix, double RhoMinMax[2], double zMinMax[2])
{
  double PointXMax[3]={-1.0e89, -1.0e89, -1.0e89};
  double PointXMin[3]={+1.0e89, +1.0e89, +1.0e89};
  for(int nx=0; nx<XMatrix->NC; nx++)
   for(int i=0; i<3; i++)
    { double Xi = XMatrix->GetEntryD(i,nx);
      PointXMax[i] = fmax(PointXMax[i], Xi);
      PointXMin[i] = fmin(PointXMin[i], Xi);
    }
  double RMax[3]={-1.0e89, -1.0e89, -1.0e89};
  double RMin[3]={+1.0e89, +1.0e89, +1.0e89};
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
    double *SurfaceXMax = G->Surfaces[ns]->RMax;
    double *SurfaceXMin = G->Surfaces[ns]->RMin;
    for(int i=0; i<3; i++)
      { RMax[i] = fmax(RMax[i], PointXMax[i] - SurfaceXMin[i]);
        if (SurfaceXMin[i]>PointXMax[i]) 
         RMin[i] = fmin(RMin[i], fabs(SurfaceXMin[i]-PointXMax[i]));
        else if (PointXMin[i]>SurfaceXMax[i])
         RMin[i] = fmin(RMin[i], fabs(PointXMin[i]-SurfaceXMax[i]));
        else
         RMin[i] = 0.0;
      }
   }
  RhoMinMax[0] = sqrt(RMin[0]*RMin[0] + RMin[1]*RMin[1]);
  RhoMinMax[1] = sqrt(RMax[0]*RMax[0] + RMax[1]*RMax[1]);
  zMinMax[0]   = RMin[2];
  zMinMax[1]   = RMax[2];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleMOIMatrixBlock(RWGGeometry *G, LayeredSubstrate *Substrate,
                            int nsa, int nsb, cdouble Omega, HMatrix *M)
{
  int Order=9;
  char *s=getenv("SCUFF_MOIME_ORDER");
  if (s && 1==sscanf(s,"%i",&Order))
   Log("Setting panel-cubature order=%i for MOI matrix elements",Order);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGSurface *Sa = G->Surfaces[nsa], *Sb=G->Surfaces[nsb];
  int OffsetA = G->BFIndexOffset[nsa], OffsetB=G->BFIndexOffset[nsb];
  int NEA = Sa->NumEdges, NEB=Sb->NumEdges;
  int nPair=0, NPair = (nsa==nsb) ? NEA*(NEA+1)/2 : NEA*NEB;
#ifndef USE_OPENMP
  Log("Computing MOI matrix block (%i,%i) (%i entries)",nsa,nsb,NEA*NEB);
#else
  int NumThreads=GetNumThreads();
  Log("Computing MOI matrix block (%i,%i) (%i entries) (%i threads)",nsa,nsb,NEA*NEB,NumThreads);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int neaneb=0; neaneb<NEA*NEB; neaneb++)
   { 
     int nea = neaneb / NEB;
     int neb = neaneb % NEB;
     if (nsa==nsb && neb<nea) continue;
     LogPercent(nPair++, NPair);

     cdouble ME;
     GetMOIMatrixElement(G, Substrate, nsa, nea, nsb, neb, Omega, &ME, Order, true);
     M->SetEntry(OffsetA + nea, OffsetB + neb, ME);
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleMOIMatrix(RWGGeometry *G, LayeredSubstrate *Substrate,
                       cdouble Omega, HMatrix *M)
{
  if ( !M || (M->NR != G->TotalBFs) || (M->NC != G->TotalBFs) )
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  /***************************************************************/
  /* pre-allocate interpolator for subtrate green's functions    */
  /***************************************************************/
  double RhoMinMax[2];
  GetRhoMinMax(G, RhoMinMax);
  bool PPIsOnly = true;
  bool Subtract = true;
  bool RetainSingularTerms = false;
  double DeltaRho=0.0, DeltaZ=0.0;
  bool Verbose = (G->LogLevel == SCUFF_VERBOSE2);
  Log("Initializing ScalarGF interpolator for Rho range (%e,%e)",RhoMinMax[0],RhoMinMax[1]);
  Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], 0.0, 0.0,
                                      PPIsOnly, Subtract, RetainSingularTerms,
                                      DeltaRho, DeltaZ, Verbose);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nsa=0; nsa<G->NumSurfaces; nsa++)
   for(int nsb=nsa; nsb<G->NumSurfaces; nsb++)
    { 
      // attempt to reuse the diagonal block of an identical previous object
      int nsm; // 'number of surface mate'
      if (nsa==nsb && (nsm=G->Mate[nsa])!=-1)
       { int ThisOffset = G->BFIndexOffset[nsa];
         int MateOffset = G->BFIndexOffset[nsm];
         int Dim = G->Surfaces[nsa]->NumBFs;
         Log("Block(%i,%i) is identical to block (%i,%i) (reusing)",nsa,nsa,nsm,nsm);
         M->InsertBlock(M, ThisOffset, ThisOffset, Dim, Dim, MateOffset, MateOffset);
       }
      else
       AssembleMOIMatrixBlock(G, Substrate, nsa, nsb, Omega, M);
    }

  // fill in lower triangle
  for(int nr=1; nr<M->NR; nr++)
   for(int nc=0; nc<nr; nc++)
    M->SetEntry(nr, nc, M->GetEntry(nc, nr) );
}

/***************************************************************/
/* Get 'Reduced Potential/Field' (RPF) matrices tabulating the */
/* contributions of individual basis functions to potentials   */
/* and fields at individual points.                            */
/*                                                             */
/* XMatrix = 3xNX matrix of evaluation point coordinates       */
/*                                                             */
/* RPFMatrix = 4*NX x NBFP matrix                              */
/*  where NBFP = NBF + NumPorts (if PortList!=0)               */
/*  or    NBFP = NBF            (if PortList==0)               */
/*                                                             */
/* RPFMatrix[4*nx + {0,1,2,3}, nbfp]                           */
/*  = contribution to {Phi,Ex,Ey,Ez} at eval point #nx due to  */
/*     (a) unit-strength basis function #nbfp (if nbfp<NBF)    */
/*     (b) RWG port #np (if nbfp = NBF + np)                   */
/*                                                             */
/* Thus, if KNP is a vector of length NBFP storing the         */
/* surface-current coefficients of the basis functions followed*/
/* by the port currents, then the (4*nx + Mu) the column       */
/* of RPFMatrix is a vector that may be dotted into KNP to     */
/* yield the Muth potential/field component at eval point #nx. */
/***************************************************************/
void GetMOIRPFMatrices(RWGGeometry *G, LayeredSubstrate *Substrate, RWGPortList *PortList,
                       cdouble Omega, HMatrix *XMatrix,
                       HMatrix **pBFRPFMatrix, HMatrix **pPortRPFMatrix)
{
  int NX           = XMatrix->NC;
  int NPFCX        = NPFC*NX; // total number of potential/field components at all points

  int NBFEdges     = G->TotalBFs;

  int NPorts       = PortList ? PortList->Ports.size() : 0;
  int NPortEdges   = PortList ? PortList->PortEdges.size() : 0;

  int NEdges       = NBFEdges + NPortEdges;
  
  HMatrix *BFRPFMatrix = *pBFRPFMatrix;
  if (BFRPFMatrix && (BFRPFMatrix->NR!=NPFCX || BFRPFMatrix->NC!=NBFEdges || BFRPFMatrix->RealComplex!=LHM_COMPLEX))
   { Warn("wrong-size BFRPFMatrix in GetMOIRPFMatrices");
     delete BFRPFMatrix;
     BFRPFMatrix=0;
   }
  if (BFRPFMatrix==0)
   BFRPFMatrix = *pBFRPFMatrix = new HMatrix(NPFCX, NBFEdges, LHM_COMPLEX);

  HMatrix *PortRPFMatrix = *pPortRPFMatrix;
  if (PortRPFMatrix && (PortRPFMatrix->NR!=NPFCX || PortRPFMatrix->NC!=NPorts || PortRPFMatrix->RealComplex!=LHM_COMPLEX) )
   { Warn("wrong-size PortRPFMatrix in GetMOIRPFMatrices");
     delete PortRPFMatrix;
     PortRPFMatrix=0;
   }
  if (PortRPFMatrix==0)
   PortRPFMatrix = *pPortRPFMatrix = new HMatrix(NPFCX, NPorts, LHM_COMPLEX);

  HMatrix *PERPFMatrix=0;
  if (NPorts>0)
   PERPFMatrix = new HMatrix(NPFCX, NPortEdges, LHM_COMPLEX);
 
  int Order=9; 

  /***************************************************************/
  /* pre-allocate interpolator for subtrate green's functions    */
  /***************************************************************/
  double RhoMinMax[2], zMinMax[2];
  GetRzMinMax(G, XMatrix, RhoMinMax, zMinMax);
  bool PPIsOnly = false;
  bool Subtract = true;
  bool RetainSingularTerms = false;
  double DeltaRho=0.0, DeltaZ=0.0;
  bool Verbose = (G->LogLevel==SCUFF_VERBOSE2);
  Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], zMinMax[0], zMinMax[1],
                                      PPIsOnly, Subtract, RetainSingularTerms,
                                      DeltaRho, DeltaZ, Verbose);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NXNE = NX*NEdges;
#ifndef USE_OPENMP
  Log("Computing MOI fields at %i points",NX);
#else
  int NumThreads=GetNumThreads();
  Log("Computing MOI fields at %i points (%i threads)",NX,NumThreads);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nxne=0; nxne<NXNE; nxne++)
   {
     int nx    = nxne / NEdges;
     int nEdge = nxne % NEdges;

     LogPercent(nxne, NXNE);

     int ns, ne;
     cdouble *RPFVector;
     if (nEdge < NBFEdges)
      { G->ResolveEdge(nEdge, &ns, &ne);
        RPFVector = (cdouble *)BFRPFMatrix->GetColumnPointer(nEdge);
      }
     else
      { int nPE = nEdge - NBFEdges;
        ns = PortList->PortEdges[nPE]->ns;
        ne = PortList->PortEdges[nPE]->ne;
        RPFVector = (cdouble *)PERPFMatrix->GetColumnPointer(nPE);
      }
    
     double *XDest = (double *)XMatrix->GetColumnPointer(nx);
     Get1BFMOIFields(G, Substrate, ns, ne, Omega, XDest, RPFVector + NPFC*nx, Order, Subtract);
   }
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (PERPFMatrix)
   { PortRPFMatrix->Zero();
     for(int nPE=0; nPE<NPortEdges; nPE++)
      { RWGPortEdge *PE = PortList->PortEdges[nPE];
        int nPort       = PE->nPort;
        double Sign     = (PE->Pol ==_PLUS ? +1.0 : -1.0);
        for(int npfcx=0; npfcx<NPFCX; npfcx++)
         PortRPFMatrix->AddEntry(npfcx,nPort,
                                 -1.0*Sign*PERPFMatrix->GetEntry(npfcx,nPE));
      }
     delete PERPFMatrix;
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetMOIFields(RWGGeometry *G, LayeredSubstrate *S, cdouble Omega,
                      HMatrix *XMatrix, HVector *KN, RWGPortList *PortList,
                      cdouble *PortCurrents, HMatrix **PFContributions)
{
  if (XMatrix->NR!=3)
   ErrExit("wrong-size (%ix%i) XMatrix in GetMOIFields",XMatrix->NR,XMatrix->NC);
  int NX = XMatrix->NC;
  int NumPorts = PortList->Ports.size();

  HMatrix *BFPFMatrix=0, *PortPFMatrix=0;
  if (PFContributions)
   { BFPFMatrix   = PFContributions[0];
     PortPFMatrix = PFContributions[1];
     if ( BFPFMatrix->NR != NPFC || BFPFMatrix->NC != NX || BFPFMatrix->RealComplex != LHM_COMPLEX)
      ErrExit("wrong-size BFPFMatrix(%i,%i) in GetMOIFields",BFPFMatrix->NR,BFPFMatrix->NC);
     if ( PortPFMatrix->NR != NPFC || PortPFMatrix->NC != NX || BFPFMatrix->RealComplex!=LHM_COMPLEX)
      ErrExit("wrong-size PortPFMatrix(%i,%i) in GetMOIFields",PortPFMatrix->NR,PortPFMatrix->NC);
   }

  HMatrix *BFRPFMatrix=0, *PortRPFMatrix=0;
  GetMOIRPFMatrices(G, S, PortList, Omega, XMatrix, &BFRPFMatrix, &PortRPFMatrix);

  cdouble *BFPFBuffer   = BFPFMatrix   ? BFPFMatrix->ZM   : new cdouble[NPFC*NX];
  HVector BFPFVector(NPFC*NX, BFPFBuffer);
  BFRPFMatrix->Apply(KN, &BFPFVector);

  cdouble *PortPFBuffer = PortPFMatrix ? PortPFMatrix->ZM : new cdouble[NPFC*NX];
  HVector PortPFVector(NPFC*NX, PortPFBuffer);
  HVector PCVector(NumPorts, PortCurrents);
  PortRPFMatrix->Apply(&PCVector, &PortPFVector);

  HMatrix *PFMatrix = new HMatrix(NPFC, NX, LHM_COMPLEX);
  VecAdd(BFPFBuffer, PortPFBuffer, PFMatrix->ZM, NPFC*NX);
  
  delete BFRPFMatrix;
  delete PortRPFMatrix;
  if (!BFPFMatrix) delete[] BFPFBuffer;
  if (!PortPFMatrix) delete[] PortPFBuffer;
  return PFMatrix;
}

} // namespace scuff

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
HMatrix *GetMOIFields(RWGGeometry *G, LayeredSubstrate *S, cdouble Omega,
                      HMatrix *XMatrix, HVector *KN, RWGPortList *PortList,
                      cdouble *PortCurrents, HMatrix *PFMatrix)
{
  int NX           = XMatrix->NC;

  int NBFEdges     = KN ? G->TotalBFs : 0;

  int NPorts       = PortList ? PortList->Ports.size() : 0;
  int NPortEdges   = PortList ? PortList->PortEdges.size() : 0;
  int NEdges       = NBFEdges + NPortEdges;
  
  int Order=9; 

  /***************************************************************/
  /* pre-allocate interpolator for subtrate green's functions    */
  /***************************************************************/
  double RhoMinMax[2], zMinMax[2];
  GetRzMinMax(G, XMatrix, RhoMinMax, zMinMax);
  bool PPIsOnly = false;
  bool Subtract = true;
  bool RetainSingularTerms = false;
  Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], zMinMax[0], zMinMax[1],
                                      PPIsOnly, Subtract, RetainSingularTerms);

  if ( PFMatrix && (PFMatrix->NR!=NPFC || PFMatrix->NC!=NX) )
   { Warn("wrong-size PFMatrix in GetMOIFields (reallocating)");
     delete PFMatrix;
     PFMatrix=0;
   }
  if (!PFMatrix)
   PFMatrix=new HMatrix(NPFC, NX, LHM_COMPLEX);
  PFMatrix->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NXNE = NX*NEdges;
#ifndef USE_OPENMP
  Log("Computing MOI fields at %i points",NX);
#else
  int NumThreads=GetNumThreads();
  HMatrix **PFByThread=0;
  bool AllocateThreadBuffers=true;
  if (AllocateThreadBuffers)
   { PFByThread = new HMatrix *[NumThreads];
     PFByThread[0] = PFMatrix;
     for(int nt=1; nt<NumThreads; nt++)
      PFByThread=new HMatrix(NPFC, NX, LHM_COMPLEX);
   }
  Log("Computing MOI fields at %i points (%i threads)",NX,NumThreads);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nxne=0; nxne<NXNE; nxne++)
   {
     int nx    = nxne / NEdges;
     int nEdge = nxne % NEdges;

     LogPercent(nxne, NXNE);

     int ns, ne;
     cdouble Weight;
     if (KN && nEdge < NBFEdges)
      { G->ResolveEdge(nEdge, &ns, &ne);
        Weight=KN->GetEntry(nEdge);
      }
     else
      { int nPE = nEdge - (KN==0 ? 0 : NBFEdges);
        RWGPortEdge *PE = PortList->PortEdges[nPE];
        ns = PE->ns;
        ne = PE->ne;
        int nPort = PE->nPort;
        double Sign = (PE->Pol==_PLUS) ? 1.0 : -1.0;
        double Perimeter = PortList->Ports[nPort]->Perimeter[
        Weight = -1.0*
      }
    
     double *XDest = (double *)XMatrix->GetColumnPointer(nx);
     Get1BFMOIFields(G, Substrate, ns, ne, Omega, XDest, RPFVector + NPFC*nx, Order, Subtract);
   }
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (PERPFMatrix)
   { PortRPFMatrix->Zero();
     for(int nPE=0; nPE<NPortEdges; nPE++)
      { RWGPortEdge *PE = PortList->PortEdges[nPE];
        int nPort       = PE->nPort;
        double Sign     = (PE->Pol ==_PLUS ? +1.0 : -1.0);
        for(int npfcx=0; npfcx<NPFCX; npfcx++)
         PortRPFMatrix->AddEntry(npfcx,nPort,
                                 -1.0*Sign*PERPFMatrix->GetEntry(npfcx,nPE));
      }
     delete PERPFMatrix;
   }
}
#endif
