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
   RWGGeometry *G;
   int nsInner, neInner, InnerOrder;

   // used for diagnostics
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
 
   PFCorrections[_PF_AX]    += Weight*V[_SGF_APAR ]*b[0];
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
     cdouble apdp[NPFC];
     int zfdim=NPFC, fdim=2*zfdim;
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
int Get1BFMOIFields(RWGGeometry *G, int ns, int ne,
                    cdouble Omega, double *XDest, cdouble *PFVector,
                    int Order, bool Subtract, cdouble *SingularTerms)
{
  //FIXME
  if (XDest[2] < -1.0e-10)
   Subtract=false;

  LayeredSubstrate *Substrate = G->Substrate;
  cdouble Eps                 = Substrate ? Substrate->EpsLayer[1] : 1.0;
  bool NeedSITerms            = (Substrate==0 || Subtract);
  bool NeedSubstrateTerms     = Substrate && (Eps!=1.0 || !isinf(Substrate->zGP));

  memset(PFVector, 0, NPFC*sizeof(cdouble));

  /***************************************************************/
  /* substrate-independent (singular) terms                      */
  /***************************************************************/
  cdouble SITermBuffer[NPFC], *SITerms = SingularTerms ? SingularTerms : SITermBuffer;
  memset(SITerms, 0, NPFC*sizeof(cdouble));
  if ( NeedSITerms )
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

  if (!NeedSubstrateTerms) return 0;

  /***************************************************************/
  /* get substrate contributions                                 */
  /*  (minus the most singular terms if Subtract==true)          */
  /***************************************************************/
  if (Substrate) Substrate->UpdateCachedEpsMu(Omega);
  MOIIntegrandData MOIData;
  MOIData.S              = Substrate;
  MOIData.Omega          = Omega;
  MOIData.Eps            = Eps;
  MOIData.Subtract       = Subtract;
  MOIData.XDest          = XDest;
  MOIData.qPoints        = 0;
  MOIData.CubaturePoints = 0;

  if (Order==-1)
   { RWGEdge *E = G->Surfaces[ns]->GetEdgeByIndex(ne);
     Order = ( VecDistance(XDest, E->Centroid) < 2.0*E->Radius ) ? 9 : 4;
   }
  int zfdim=NPFC, fdim=2*zfdim;
  cdouble SubstrateTerms[NPFC];
  GetBFCubature2(G, ns, ne, MOIPFIntegrand, (void *)&MOIData, fdim, Order, (double *)SubstrateTerms);
  VecPlusEquals(PFVector, 1.0, SubstrateTerms, NPFC);

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
int GetMOIMatrixElement(RWGGeometry *G,
                        int nsa, int nea, int nsb, int neb,
                        cdouble Omega, cdouble *ME,
                        int Order, bool Subtract, cdouble *Terms)
{
  ME[0]=0.0;
  if (Terms) 
   memset(Terms, 0, 8*sizeof(cdouble));

  LayeredSubstrate *Substrate = G->Substrate;
  cdouble Eps = (Substrate && Substrate->NumLayers>0 ? Substrate->EpsLayer[1] : 1.0);

  bool NeedSITerms        = (Substrate==0 || Subtract);
  bool NeedSubstrateTerms = (Substrate && (Eps!=1.0 || !isinf(Substrate->zGP)));
  int qPoints=0;

  /***************************************************************/
  /* get substrate contributions *********************************/
  /***************************************************************/ 
  if (Substrate)
   { 
     Substrate->UpdateCachedEpsMu(Omega);
     Eps = Substrate->EpsLayer[1];

     bool SeparateTerms = (Terms!=0);

     MOIIntegrandData MOIData;
     MOIData.S              = Substrate;
     MOIData.Omega          = Omega;
     MOIData.Eps            = Eps;
     MOIData.Subtract       = Subtract;
     MOIData.SeparateTerms  = SeparateTerms;
     MOIData.qPoints        = 0;
     MOIData.CubaturePoints = 0;

     if (Order==-1) 
      { double rRel;
        int ncv=AssessBFPair(G->Surfaces[nsa], nea, G->Surfaces[nsb], neb, &rRel);
        Order = ( (ncv>0 || rRel <=1.0) ? 9 : 4);
      }

     int zfdim = SeparateTerms ? 6 : 1, fdim=2*zfdim;
     GetBFBFCubature2(G, nsa, nea, nsb, neb, MOIMEIntegrand, (void *)&MOIData,
                      fdim, Order, (double *)(SeparateTerms ? Terms : ME));

     qPoints += MOIData.qPoints;

     if (SeparateTerms)
      ME[0] = Terms[0] + Terms[1];
   }

  /***************************************************************/
  /* add 'bare' (substrate-independent) EEIs                     */
  /***************************************************************/
  if (NeedSITerms)
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

  return qPoints;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void MOIMEIntegrand_NC(double XDest[3], double b[3], double Divb,
                       void *UserData, double Weight, double *Integral)
{
   MOIIntegrandData *OuterData = (MOIIntegrandData *)UserData;
   cdouble Omega                = OuterData->Omega;
   bool Subtract                = OuterData->Subtract;
   bool SeparateTerms           = OuterData->SeparateTerms;

   RWGGeometry *G               = OuterData->G;
   int nsInner                  = OuterData->nsInner;
   int neInner                  = OuterData->neInner;
   int InnerOrder               = OuterData->InnerOrder;

   cdouble PFVector[NPFC], STBuffer[NPFC], *SingularTerms = (SeparateTerms ? STBuffer : 0);
   int qPoints = Get1BFMOIFields(G, nsInner, neInner, Omega, XDest, PFVector, InnerOrder, Subtract, SingularTerms);

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
      METerms[0]+= Weight * ( iw*(b[0]*PFVector[_PF_AX] + b[1]*PFVector[_PF_AY]) + Divb*PFVector[_PF_PHI] );
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
int GetMOIMatrixElement_NC(RWGGeometry *G,
                           int nsa, int nea, int nsb, int neb,
                           cdouble Omega, cdouble *ME,
                           int InnerOrder, int OuterOrder, bool Subtract, cdouble *Terms)
{
  /***************************************************************/
  /* get substrate contributions *********************************/
  /***************************************************************/
  LayeredSubstrate *Substrate = G->Substrate;
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
 
  int zfdim = (SeparateTerms ? 4 : 1), fdim=2*zfdim;
  GetBFCubature2(G, nsa, nea, MOIMEIntegrand_NC, (void *)&MOIData,
                 fdim, OuterOrder, (double *)(Terms ? Terms : ME));
  if (Terms)
   ME[0] = Terms[0]+Terms[1]+Terms[2]+Terms[3];

  return MOIData.qPoints;

}

/***************************************************************/
/* shortest, longest distances between points on given surfaces*/
/***************************************************************/
void GetRhoMinMax(RWGGeometry *G, int nsa, int nsb, double RhoMinMax[2])
{ 
  double RMax[3], RMin[3];
  double *RaMax = G->Surfaces[nsa]->RMax, *RaMin = G->Surfaces[nsa]->RMin;
  double *RbMax = G->Surfaces[nsb]->RMax, *RbMin = G->Surfaces[nsb]->RMin;
  for(int Mu=0; Mu<3; Mu++)
   { RMax[Mu] = fmax( fabs(RbMax[Mu]-RaMin[Mu]), fabs(RaMax[Mu]-RbMin[Mu]) );
     RMin[Mu] = (   ( RaMax[Mu] < RbMin[Mu] ) ? RbMin[Mu]-RaMax[Mu]
                  : ( RbMax[Mu] < RaMin[Mu] ) ? RaMin[Mu]-RbMax[Mu]
                  : 0.0
                );
   }
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
      { RMax[i] = fmax(RMax[i], fabs(PointXMax[i] - SurfaceXMin[i]) );
        RMax[i] = fmax(RMax[i], fabs(PointXMin[i] - SurfaceXMax[i]) );
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
void AssembleMOIMatrixBlock(RWGGeometry *G, int nsa, int nsb,
                            cdouble Omega, HMatrix *Block, int OffsetA, int OffsetB,
                            EquivalentEdgePairTable *EEPTable)
{
  int Order=-1;
  char *s=getenv("SCUFF_MOIME_ORDER");
  if (s && 1==sscanf(s,"%i",&Order))
   Log("Setting panel-cubature order=%i for MOI matrix elements",Order);

  /***************************************************************/
  /* pre-allocate interpolator for subtrate green's functions if */
  /* it is not already there                                     */
  /***************************************************************/
  if (G->Substrate)
   { 
     double RhoMinMax[2];
     GetRhoMinMax(G, nsa, nsb, RhoMinMax);
     bool PPIsOnly = true;
     bool Subtract = true;
     bool RetainSingularTerms = false;
     bool Verbose = (G->LogLevel == SCUFF_VERBOSE2);
     Log("Initializing ScalarGF interpolator for Rho range (%e,%e)",RhoMinMax[0],RhoMinMax[1]);
     G->Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], 0.0, 0.0,
                                            PPIsOnly, Subtract, RetainSingularTerms, Verbose);
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGSurface *Sa = G->Surfaces[nsa], *Sb=G->Surfaces[nsb];
  int NEA=Sa->NumEdges, NEB=Sb->NumEdges;

  IrreducibleEdgePairList *IEPList = EEPTable ? EEPTable->GetIrreducibleEdgePairList(nsa, nsb) : 0;
  int NumPairs = IEPList ? IEPList->size() : NEA*NEB;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumThreads=GetNumThreads();
  int FullNumPairs = (nsa==nsb ? NEA*(NEA+1)/2 : NEA*NEB);
  Log("Computing MOI matrix block (%i,%i) (%i/%i entries) (%i threads)",nsa,nsb,NumPairs,FullNumPairs,NumThreads);
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nPair=0; nPair<NumPairs; nPair++)
   { 
     LogPercent(nPair, NumPairs);

     int nea, neb;
     if (IEPList)
      { int nfeaParent, nfebParent, nsaParent, nsbParent;
        EEPTable->ResolveEdgePairIndex( (*IEPList)[nPair].ParentPair, &nfeaParent, &nfebParent);
        G->ResolveEdge(nfeaParent, &nsaParent, &nea);
        G->ResolveEdge(nfebParent, &nsbParent, &neb);
        if (nsaParent!=nsa || nsbParent!=nsb) continue; // this should never happen
      }
     else
      { nea = nPair/NEB;
        neb = nPair%NEB;
        if (nsa==nsb && neb<nea) continue;
      }

     cdouble ME;
     GetMOIMatrixElement(G, nsa, nea, nsb, neb, Omega, &ME, Order, true);
     Block->SetEntry(OffsetA + nea, OffsetB + neb, ME);
   }

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (IEPList)
   { Log("Filling in matrix entries for child edge pairs ...");
     for(int nPair=0; nPair<NumPairs; nPair++)
      { 
        int nfeaParent, nsaParent, neaParent, nfebParent, nsbParent, nebParent;
        EEPTable->ResolveEdgePairIndex( (*IEPList)[nPair].ParentPair, &nfeaParent, &nfebParent);
        G->ResolveEdge(nfeaParent, &nsaParent, &neaParent);
        G->ResolveEdge(nfebParent, &nsbParent, &nebParent);
        if (nsaParent!=nsa || nsbParent!=nsb) continue; // this should never happen

        cdouble ME = Block->GetEntry(OffsetA + neaParent, OffsetB + nebParent);

        iVec *ChildPairs = (*IEPList)[nPair].ChildPairs;
        if (!ChildPairs) continue;

        for(size_t nc=0; nc<ChildPairs->size(); nc++)
         { int ChildPair = (*ChildPairs)[nc];
           double Sign = (ChildPair < 0 ? -1.0 : 1.0);
           if (ChildPair<0) ChildPair*=-1;
           int nfeaChild, nfebChild;
           EEPTable->ResolveEdgePairIndex(ChildPair,&nfeaChild, &nfebChild);
           int nsaChild, neaChild; G->ResolveEdge(nfeaChild, &nsaChild, &neaChild);
           int nsbChild, nebChild; G->ResolveEdge(nfebChild, &nsbChild, &nebChild);
           if (nsaChild!=nsa || nsbChild!=nsb) continue;
           Block->SetEntry(OffsetA + neaChild, OffsetB + nebChild, Sign*ME);
         }
      }
     Log(" ...done with equivalent pairs");
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AssembleMOIMatrix(RWGGeometry *G, cdouble Omega, HMatrix *M, EquivalentEdgePairTable *EEPTable)
{
  if ( !M || (M->NR != G->TotalBFs) || (M->NC != G->TotalBFs) )
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  /***************************************************************/
  /* pre-allocate interpolator for subtrate green's functions    */
  /***************************************************************/
  if (G->Substrate)
   { 
     double RhoMinMax[2];
     GetRhoMinMax(G, RhoMinMax);
     bool PPIsOnly = true;
     bool Subtract = true;
     bool RetainSingularTerms = false;
     bool Verbose = (G->LogLevel == SCUFF_VERBOSE2);
     Log("Initializing ScalarGF interpolator for Rho range (%e,%e)",RhoMinMax[0],RhoMinMax[1]);
     G->Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], 0.0, 0.0,
                                            PPIsOnly, Subtract, RetainSingularTerms, Verbose);
   }

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
       AssembleMOIMatrixBlock(G, nsa, nsb, Omega, M, G->BFIndexOffset[nsa], G->BFIndexOffset[nsb], EEPTable);
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
/* RPFMatrix[0][7*nx + npfc, nbfp]                             */
/*  = contribution of basis function #nbfp, weighted with unit */
/*    strength, to potential/field component #npfc @ point #nx */
/*                                                             */
/* RPFMatrix[1][7*nx + npfc, np]                               */
/*  = contribution of port #np, driven with unit input current,*/
/*    to potential/field component #npfc @ point #nx           */
/*                                                             */
/* Thus, if KN is a vector of length NBF storing the           */
/* surface-current coefficients of the basis functions,        */
/* then RPFMatrix[0] * KN is the payload of a 7xNX matrix      */
/* whose (npfc,nx) entry is the potential/field component npfc */
/* at point #nx due to RWG functions weighted by KN.           */
/*                                                             */
/* Similarly, if PC is a vector of length NPorts storing       */
/* port currents, then RPFMatrix[1] * KN is the payload of a   */
/* 7xNX matrix giving the contributions of driven ports.       */
/***************************************************************/
void GetMOIRPFMatrices(RWGGeometry *G, RWGPortList *PortList,
                       cdouble Omega, HMatrix *XMatrix,
                       HMatrix **pBFRPFMatrix, HMatrix **pPortRPFMatrix)
{
  int NX               = XMatrix->NC;
  int NPFCX            = NPFC*NX; // total number of potential/field components at all points

  int NBFEdges         = 0;
  HMatrix *BFRPFMatrix = 0;
  bool NeedBFRPFMatrix = (pBFRPFMatrix!=0); 
  if (NeedBFRPFMatrix)
   { NBFEdges = G->TotalBFs;
     BFRPFMatrix = *pBFRPFMatrix = CheckHMatrix(*pBFRPFMatrix, NPFCX, NBFEdges, LHM_COMPLEX);
   }

  // PERPFMatrix ("PortEdgeRPFMatrix") is a port-edge-resolved version of PortRPFMatrix
  // that is reduced after the loop to obtain PortRPFMatrix[1].
  int NPortEdges         = 0;
  HMatrix *PERPFMatrix   = 0;
  bool NeedPortRPFMatrix = (pPortRPFMatrix!=0);
  if (NeedPortRPFMatrix)
   { NPortEdges   = PortList->PortEdges.size();
     PERPFMatrix  = new HMatrix(NPFCX, NPortEdges, LHM_COMPLEX);
   }

  /***************************************************************/
  /* pre-allocate interpolator for subtrate green's functions    */
  /***************************************************************/
  bool Subtract = (XMatrix->GetEntryD(2,0) >= 0.0);
  if (G->Substrate)
   { 
     double RhoMinMax[2], zMinMax[2];
     GetRzMinMax(G, XMatrix, RhoMinMax, zMinMax);
     bool PPIsOnly = false;
     bool RetainSingularTerms = false;
     bool Verbose = (G->LogLevel==SCUFF_VERBOSE2);
     G->Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], zMinMax[0], zMinMax[1],
                                            PPIsOnly, Subtract, RetainSingularTerms, Verbose);
   }

  /***************************************************************/
  /* loop over all (edge, eval point) pairs to get contributions */
  /* of individual basis functions (and optionally individual    */
  /* port edges) to potentials and fields at eval points         */
  /***************************************************************/
  int NEdges = NBFEdges + NPortEdges;
  int NXNE   = NX*NEdges;
  int NumThreads=GetNumThreads();
  Log("Computing (%i,%i) MOIRPF matrix (%i threads)",BFRPFMatrix->NR,BFRPFMatrix->NC,NumThreads);
#ifdef USE_OPENMP
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
     int Order=-1;
     Get1BFMOIFields(G, ns, ne, Omega, XDest, RPFVector + NPFC*nx, Order, Subtract);
   }
  
  /***************************************************************/
  /* collapse the port-edge resolved RPF matrix, which gives the */
  /* potentials and fields due to individual port edges, into the*/
  /* port RPF matrix (NumPorts columns), which gives the full PFs*/
  /* due to the port.                                            */
  /***************************************************************/
  if (NeedPortRPFMatrix)
   { int NPorts = PortList->Ports.size();
     HMatrix *PortRPFMatrix = *pPortRPFMatrix = CheckHMatrix(*pPortRPFMatrix, NPFCX, NPorts, LHM_COMPLEX);
     PortRPFMatrix->Zero();
     for(int nPE=0; nPE<NPortEdges; nPE++)
      { RWGPortEdge *PE        = PortList->PortEdges[nPE];
        int nPort              = PE->nPort;
        int Pol                = PE->Pol;
        double PolSign         = (Pol ==_PLUS ? +1.0 : -1.0);
        double Perimeter       = PortList->Ports[nPort]->Perimeter[Pol];
        double Weight          = PE->Sign*PolSign / Perimeter;
        cdouble *PortRPFVector = (cdouble *)PortRPFMatrix->GetColumnPointer(nPort);
        cdouble *PERPFVector   = (cdouble *)PERPFMatrix->GetColumnPointer(nPE);
        VecPlusEquals(PortRPFVector, Weight, PERPFVector, NPFCX);
      }
     delete PERPFMatrix;
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RFSolver::GetFieldsViaRPFMatrices(HMatrix *XMatrix)
{
  if (XMatrix->NR!=3)
   ErrExit("wrong-size (%ix%i) XMatrix in GetMOIFields",XMatrix->NR,XMatrix->NC);
  int NX = XMatrix->NC;

  // compute reduced potential-field matrices
  Log("Getting RPF matrices for %i points...",NX);
  HMatrix *BFRPFMatrix = 0, *PortRPFMatrix = 0;
  GetMOIRPFMatrices( G, PortList, Omega, XMatrix,
                     (RetainContributions & CONTRIBUTION_BF)   ? &BFRPFMatrix   : 0,
                     (RetainContributions & CONTRIBUTION_PORT) ? &PortRPFMatrix : 0
                   );
                    
  // create output matrix
  HMatrix *PFMatrix = new HMatrix(NPFC, NX, LHM_COMPLEX);

  if (BFRPFMatrix)
   { HVector Payload(NPFC*NX, PFMatrix->ZM);
     BFRPFMatrix->Apply(KN, &Payload);
     delete BFRPFMatrix;
   }

  // add port contribution
  if (PortRPFMatrix)
   { HVector Payload(NPFC*NX, LHM_COMPLEX);
     HVector PCVector(NumPorts, PortCurrents);
     PortRPFMatrix->Apply(&PCVector, &Payload);
     HMatrix PortPFContribution(NPFC, NX, Payload.ZV);
     PFMatrix->Add(&PortPFContribution);
     delete PortRPFMatrix;
   }
  
  Log("...done with MOI fields via RPF  matrices");
  return PFMatrix;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *RFSolver::GetFields(HMatrix *XMatrix, HMatrix *PFMatrix)
{
  int NX           = XMatrix->NC;

  int NBFEdges     = ( (RetainContributions & CONTRIBUTION_BF)   ? G->TotalBFs : 0 ); 
  int NPortEdges   = ( (RetainContributions & CONTRIBUTION_PORT) ? PortList->PortEdges.size() : 0);
  int NEdges       = NBFEdges + NPortEdges;
  
  int Order=-1; CheckEnv("SCUFF_MOIFIELDS_ORDER",&Order);

  /***************************************************************/
  /* pre-allocate interpolator for subtrate green's functions    */
  /***************************************************************/
  bool Subtract = (XMatrix->GetEntryD(2,0)>=0.0);
  if (G->Substrate)
   {
     double RhoMinMax[2], zMinMax[2];
     GetRzMinMax(G, XMatrix, RhoMinMax, zMinMax);
     bool PPIsOnly = false;
     bool RetainSingularTerms = false;
     G->Substrate->InitScalarGFInterpolator(Omega, RhoMinMax[0], RhoMinMax[1], zMinMax[0], zMinMax[1],
                                            PPIsOnly, Subtract, RetainSingularTerms);
   }

  PFMatrix=CheckHMatrix(PFMatrix,NPFC,NX,LHM_COMPLEX,"GetMOIFields");
  PFMatrix->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NXNE = NX*NEdges;
  int NumThreads=GetNumThreads();
  bool AllocateThreadBuffers=false, NonCritical=false;
  HMatrix **PFByThread=0;
  Log("Computing MOI fields at %ix%i = %i edge-point pairs (%i threads)",NEdges,NX,NXNE,NumThreads);
#ifdef USE_OPENMP
  NonCritical = CheckEnv("SCUFF_NONCRITICAL");
  AllocateThreadBuffers = !CheckEnv("SCUFF_NO_THREAD_BUFFERS");
  Log("%s thread buffers",AllocateThreadBuffers ? "Allocating" : "Not allocating");
  if (AllocateThreadBuffers)
   { PFByThread = new HMatrix *[NumThreads];
     PFByThread[0] = PFMatrix;
     for(int nt=1; nt<NumThreads; nt++)
      PFByThread[nt] = new HMatrix(NPFC, NX, LHM_COMPLEX);
   }
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int nxne=0; nxne<NXNE; nxne++)
   {
     int nx    = nxne / NEdges;
     int nEdge = nxne % NEdges;
     LogPercent(nxne, NXNE);

     int ns, ne;
     cdouble Weight=0.0;
     if (nEdge<NBFEdges)
      { G->ResolveEdge(nEdge, &ns, &ne);
        Weight = KN->GetEntry(nEdge);
      }
     else
      { int nPE          = nEdge - NBFEdges;
        RWGPortEdge *PE  = PortList->PortEdges[nPE];
        int nPort        = PE->nPort;
        double PolSign   = (PE->Pol==_PLUS ? 1.0 : -1.0);
        double Perimeter = PortList->Ports[nPort]->Perimeter[PE->Pol];
        ns               = PE->ns;
        ne               = PE->ne;
        Weight           = PE->Sign*PolSign*PortCurrents[nPort] / Perimeter;
      }
     if (Weight==0.0) continue;
    
     double *XDest = (double *)XMatrix->GetColumnPointer(nx);
     cdouble DeltaPF[NPFC];
     Get1BFMOIFields(G, ns, ne, Omega, XDest, DeltaPF, Order, Subtract);

     if (AllocateThreadBuffers)
      { cdouble *PF = (cdouble *)PFByThread[GetThreadNum()]->GetColumnPointer(nx);
        VecPlusEquals(PF, Weight, DeltaPF, NPFC);
      }
     else
      { cdouble *PF = (cdouble *)PFMatrix->GetColumnPointer(nx);
        if (NonCritical)
         VecPlusEquals(PF, Weight, DeltaPF, NPFC);
        else
         {
#pragma omp critical
          {
            VecPlusEquals(PF, Weight, DeltaPF, NPFC);
          }
         }
      }
   }

  /***************************************************************/
  /* accumulate per-thread contributions                         */
  /***************************************************************/
  if (AllocateThreadBuffers)
   { for(int nt=1; nt<NumThreads; nt++)
      { PFMatrix->Add(PFByThread[nt]);
        delete PFByThread[nt];
      }
     delete[] PFByThread;
   }
  return PFMatrix;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RFSolver::GetFields(double X[3], cdouble PF[NPFC])
{ 
  HMatrix XMatrix(3,1,X);
  HMatrix PFMatrix(NPFC,1,PF);
  GetFields(&XMatrix, &PFMatrix);
}

} // namespace scuff
