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
 * GCMatrixElements.cc      -- routine for evaluating matrix
 *                          -- elements of the the G and C dyadics
 *                          -- between RWG basis functions
 *                          
 *                          -- (before 1/2016 this routine was 
 *                          -- known as "GetEdgeEdgeInteractions")
 * 
 * homer reid 11/2005 -- 1/2016
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "PanelCubature.h"
#include "TaylorDuffy.h"

namespace scuff {

#define II cdouble(0,1)

#define FIBBI_PEFIE1_RM1      0
#define FIBBI_PEFIE1_R0       1
#define FIBBI_PEFIE1_R1       2
#define FIBBI_PEFIE1_R2       3

#define FIBBI_PEFIE2_RM1      4
#define FIBBI_PEFIE2_R0       5
#define FIBBI_PEFIE2_R1       6
#define FIBBI_PEFIE2_R2       7

#define FIBBI_PMFIE_RM3       8
#define FIBBI_PMFIE_RM1       9
#define FIBBI_PMFIE_R0       10
#define FIBBI_PMFIE_R1       11

#define FIBBI_PEFIE1_RX_RM3  12
#define FIBBI_PEFIE1_RX_RM1  15
#define FIBBI_PEFIE1_RX_R0   18
#define FIBBI_PEFIE1_RX_R1   21

#define FIBBI_PEFIE2_RX_RM3  24
#define FIBBI_PEFIE2_RX_RM1  27
#define FIBBI_PEFIE2_RX_R0   30
#define FIBBI_PEFIE2_RX_R1   33

#define FIBBI_BXBX_RM1       36
#define FIBBI_BXBX_R0        39
#define FIBBI_BXBX_R1        42

#define FIBBI_PMFIE_RX_RM3   45
#define FIBBI_PMFIE_RX_R0    48

#define NUMFIBBIS            51

cdouble ExpRel(cdouble x, int n);

/***************************************************************/
/* integrand routine passed to GetBFBFCubature2 for computing  */
/* matrix elements of the G, C kernels                         */
/***************************************************************/
typedef struct GetGCMEData
 {
   GetGCMEArgStruct *Args;
   bool GetFIBBIs;
   bool DeSingularize;
   double *ZHatAP, *ZHatAM, *ZHatBP, *ZHatBM;
 } GetGCMEData;

void GCMEIntegrand(double xA[3], double bA[3], double DivbA,
                   double xB[3], double bB[3], double DivbB,
                   void *UserData, double Weight, double *Integral)
{
  GetGCMEData *Data = (GetGCMEData *)UserData;

  double *Displacement        = Data->Args->Displacement;
  int NumRegions              = Data->Args->NumRegions;
  cdouble *kVector            = Data->Args->k;
  GBarAccelerator **GBA       = Data->Args->GBA;

  bool NeedGC                 = Data->Args->NeedGC;
  bool NeedDW                 = Data->Args->NeedDW;
  bool NeedForce              = Data->Args->NeedForce;
  bool NeedTorque             = Data->Args->NeedTorque;
  bool NeedSpatialDerivatives = (NeedForce || NeedTorque);
  if (NeedDW) NeedGC=true;

  bool DeSingularize          = Data->DeSingularize;
  
  /***************************************************************/
  /* compute polynomial factors for this pair of RWG functions   */
  /***************************************************************/
  double R[3];
  VecSub(xA, xB, R);
  if (Displacement)
   VecPlusEquals(R, -1.0, Displacement);
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double r=sqrt(r2);

  double bdb  = VecDot(bA, bB);
  double DbDb = DivbA*DivbB;
  double bxb[3]; VecCross(bA, bB, bxb);
  double PMFIE = VecDot(bxb, R);

  double RPA[3], RPB[3];
  if (NeedSpatialDerivatives)
   { double *ZHatA = (DivbA > 0.0 ? Data->ZHatAP : Data->ZHatAM);
     double *ZHatB = (DivbB > 0.0 ? Data->ZHatBP : Data->ZHatBM);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if (DivbA<0.0) return;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     memcpy(RPA, R, 3*sizeof(double));
     memcpy(RPB, R, 3*sizeof(double));
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { RPA[Mu] -= ZHatA[Mu]*ZHatA[Nu]*R[Nu];
         RPB[Mu] -= ZHatB[Mu]*ZHatB[Nu]*R[Nu];
       };
   };

  /***************************************************************/
  /* compute frequency-independent contributions if that's what  */
  /* was requested                                               */
  /***************************************************************/
  if (Data->GetFIBBIs)
   {
     double rm1 = (r==0.0 ? 0.0 : 1.0/r);
     double rm3 = (r==0.0 ? 0.0 : 1.0/(r*r2));
     Integral[ FIBBI_PEFIE1_RM1   ] += Weight*bdb*rm1;
     Integral[ FIBBI_PEFIE1_R0    ] += Weight*bdb;
     Integral[ FIBBI_PEFIE1_R1    ] += Weight*bdb*r;
     Integral[ FIBBI_PEFIE1_R2    ] += Weight*bdb*r2;

     Integral[ FIBBI_PEFIE2_RM1   ] += Weight*DbDb*rm1;
     Integral[ FIBBI_PEFIE2_R0    ] += Weight*DbDb;
     Integral[ FIBBI_PEFIE2_R1    ] += Weight*DbDb*r;
     Integral[ FIBBI_PEFIE2_R2    ] += Weight*DbDb*r2;

     Integral[ FIBBI_PMFIE_RM3 ] += Weight*PMFIE*rm3;
     Integral[ FIBBI_PMFIE_RM1 ] += Weight*PMFIE*rm1;
     Integral[ FIBBI_PMFIE_R0  ] += Weight*PMFIE;
     Integral[ FIBBI_PMFIE_R1  ] += Weight*PMFIE*r;

     if (NeedForce)
      for(int Mu=0; Mu<3; Mu++)
       { Integral[ FIBBI_PEFIE1_RX_RM3 + Mu ] += Weight*bdb*R[Mu]*rm3;
         Integral[ FIBBI_PEFIE1_RX_RM1 + Mu ] += Weight*bdb*R[Mu]*rm1;
         Integral[ FIBBI_PEFIE1_RX_R0  + Mu ] += Weight*bdb*R[Mu];
         Integral[ FIBBI_PEFIE1_RX_R1  + Mu ] += Weight*bdb*R[Mu]*r;
     
         Integral[ FIBBI_PEFIE2_RX_RM3 + Mu ] += Weight*DbDb*RPA[Mu]*rm3;
  //       Integral[ FIBBI_PEFIE2_RX_RM3 + Mu ] += Weight*DbDb*R[Mu]*rm3;
         Integral[ FIBBI_PEFIE2_RX_RM1 + Mu ] += Weight*DbDb*R[Mu]*rm1;
         Integral[ FIBBI_PEFIE2_RX_R0  + Mu ] += Weight*DbDb*R[Mu];
         Integral[ FIBBI_PEFIE2_RX_R1  + Mu ] += Weight*DbDb*R[Mu]*r;

         Integral[ FIBBI_BXBX_RM1      + Mu ] += Weight*bxb[Mu]*rm1;
         Integral[ FIBBI_BXBX_R0       + Mu ] += Weight*bxb[Mu];
         Integral[ FIBBI_BXBX_R1       + Mu ] += Weight*bxb[Mu]*r;

         Integral[ FIBBI_PMFIE_RX_RM3  + Mu ] += Weight*PMFIE*R[Mu]*rm3;
         Integral[ FIBBI_PMFIE_RX_R0   + Mu ] += Weight*PMFIE*R[Mu];
       };
     return;

   };

  /***************************************************************/
  /* compute wavenumber-dependent integrand quantities for all   */
  /* requested regions                                           */
  /***************************************************************/
  int nzi=0;
  cdouble *zIntegral = (cdouble *)Integral;
  for(int nr=0; nr<NumRegions; nr++)
   { 
     /***************************************************************/
     /* compute kernel factors for each region **********************/
     /***************************************************************/
     cdouble G0, dG[3], dGDS[3], PsiS, ddGBuffer[9], *ddG=(NeedSpatialDerivatives ? ddGBuffer : 0);
     cdouble k  = kVector[nr];
     cdouble ik = II*k, ikr=ik*r, k2=k*k, k3=k2*k;
   
     if (GBA[nr])
      { 
        G0=GetGBar(R, GBA[nr], dG, ddG);
        memset(dGDS,0,3*sizeof(cdouble));
      }
     else 
      { 
        cdouble Psi;
        if (r==0.0)
         { G0   = DeSingularize ? 0.0 : ik/(4.0*M_PI);
           Psi  = DeSingularize ? 0.0 : -II*k3/(12.0*M_PI);
           PsiS = 0.0;
         }
        else
         { G0 = DeSingularize ? ExpRel(ikr,4) : exp(ikr);
           G0 /= (4.0*M_PI*r);
           Psi = G0*(ikr-1.0)/r2;
           PsiS = DeSingularize ? 0.0 : -1.0/(4.0*M_PI*r*r*r);
         };
        dG[0]   = R[0]*Psi;
        dG[1]   = R[1]*Psi;
        dG[2]   = R[2]*Psi;
        dGDS[0] = R[0]*(Psi-PsiS);
        dGDS[1] = R[1]*(Psi-PsiS);
        dGDS[2] = R[2]*(Psi-PsiS);

        if (ddG)
         { 
           double ImPsi = imag(Psi);
           double ImZeta;
           double r4=r2*r2, kr5=pow( real(k), 5.0 );
           if (r==0.0)
            ImZeta = DeSingularize ? 0.0 : kr5 / (4.0*M_PI);
           else
            ImZeta = imag( G0*(ikr*ikr - 3.0*ikr + 3.0)/r4 );

           ddG[3*0+0] = II*(ImPsi + R[0]*R[0]*ImZeta);
           ddG[3*1+1] = II*(ImPsi + R[1]*R[1]*ImZeta);
           ddG[3*2+2] = II*(ImPsi + R[2]*R[2]*ImZeta);
           ddG[3*0+1] = ddG[3*1+0] = II*R[0]*R[1]*ImZeta;
           ddG[3*0+2] = ddG[3*2+0] = II*R[0]*R[2]*ImZeta;
           ddG[3*1+2] = ddG[3*2+1] = II*R[1]*R[2]*ImZeta;
         };
      };

     /***************************************************************/
     /* assemble integrands *****************************************/
     /***************************************************************/

     /*- G, C integrals ---------------------------------------------*/
     double PEFIE1  = bdb;
     cdouble PEFIE2 = DbDb/(k*k);
     cdouble PEFIE  = PEFIE1 - PEFIE2;
     if (NeedGC)
      { zIntegral[nzi++] += Weight*bdb*G0;
        zIntegral[nzi++] += Weight*DbDb*G0 / k2;
        zIntegral[nzi++] += Weight*(bxb[0]*dG[0] + bxb[1]*dG[1] + bxb[2]*dG[2]);
      };

     /*- \partial_\omega G,C ----------------------------------------*/
     if (NeedDW)
      { 
        cdouble G0Prime;
        if (GBA[nr])
         G0Prime = -II*r*(R[0]*dG[0] + R[1]*dG[1] + R[2]*dG[2]) + k*r2*G0;
        else
         G0Prime = II*exp(II*k*r)/(4.0*M_PI);

        zIntegral[nzi++] += Weight*PEFIE*G0Prime;
        if (DeSingularize && r>0.0)
         zIntegral[nzi++] += Weight*(-1.0)*k*PMFIE*exp(II*k*r)/(4.0*M_PI*r);
        else
         zIntegral[nzi++] += Weight*(-1.0)*k*PMFIE*G0;
      };

     /*- \partial_i G,C ---------------------------------------------*/
     if (NeedForce)
      { 
        zIntegral[nzi++] += Weight*(PEFIE1*dG[0] - PEFIE2*(dGDS[0] + RPA[0]*PsiS));
        zIntegral[nzi++] += Weight*(PEFIE1*dG[1] - PEFIE2*(dGDS[1] + RPA[1]*PsiS));
        zIntegral[nzi++] += Weight*(PEFIE1*dG[2] - PEFIE2*(dGDS[2] + RPA[2]*PsiS));
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//zIntegral[nzi++] += Weight*PEFIE*dG[0];
//zIntegral[nzi++] += Weight*PEFIE*dG[1];
//zIntegral[nzi++] += Weight*PEFIE*dG[2];
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
        zIntegral[nzi++] += Weight*(bxb[0]*ddG[0*3+0] + bxb[1]*ddG[0*3 + 1] + bxb[2]*ddG[0*3+2]);
        zIntegral[nzi++] += Weight*(bxb[0]*ddG[1*3+0] + bxb[1]*ddG[1*3 + 1] + bxb[2]*ddG[1*3+2]);
        zIntegral[nzi++] += Weight*(bxb[0]*ddG[2*3+0] + bxb[1]*ddG[2*3 + 1] + bxb[2]*ddG[2*3+2]);
      };

     /*- \tilde_i G, C ----------------------------------------------*/
     if (NeedTorque)
      { 
        double XT[3]; VecSub(xA, Data->Args->XTorque, XT);

        cdouble XTxdG[3], bxbDotddG[3], bADotdG=0.0;
        for(int Mu=0; Mu<3; Mu++)
         { int MP1=(Mu+1)%3, MP2=(Mu+2)%3;
           XTxdG[Mu]     = XT[MP1]*dG[MP2] - XT[MP2]*dG[MP1];
           bxbDotddG[Mu] = bxb[0]*ddG[Mu*3 + 0]
                          +bxb[1]*ddG[Mu*3 + 1]
                          +bxb[2]*ddG[Mu*3 + 2];
           bADotdG += bA[Mu]*dG[Mu];
         };
    
        cdouble CT3[3];
        CT3[0] = XT[1]*bxbDotddG[2] - XT[2]*bxbDotddG[1];
        CT3[1] = XT[2]*bxbDotddG[0] - XT[0]*bxbDotddG[2];
        CT3[2] = XT[0]*bxbDotddG[1] - XT[1]*bxbDotddG[0];

        zIntegral[nzi++] += Weight*( bxb[0]*G0 + PEFIE*XTxdG[0] );
        zIntegral[nzi++] += Weight*( bxb[1]*G0 + PEFIE*XTxdG[1] );
        zIntegral[nzi++] += Weight*( bxb[2]*G0 + PEFIE*XTxdG[2] );
        zIntegral[nzi++] += Weight*( II*imag( bADotdG*bB[0] - bdb*dG[0] + CT3[0]));
        zIntegral[nzi++] += Weight*( II*imag( bADotdG*bB[1] - bdb*dG[1] + CT3[1]));
        zIntegral[nzi++] += Weight*( II*imag( bADotdG*bB[2] - bdb*dG[2] + CT3[2]));

      };

   }; // for(int nr=0; nr<NumRegions; nr++)
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void AddFIBBIContributions(double *FIBBIs, cdouble k,
                           bool NeedGC, bool NeedForce,
                           cdouble Gab[NUMGCMES],
                           cdouble ikCab[NUMGCMES])
{
  cdouble ik=II*k, ik2=ik*ik, ik3=ik2*ik, ik4=ik3*ik;
  cdouble OOK2=1.0/(k*k);
  cdouble C[4], D[4];

  C[0]= 1.0/(4.0*M_PI);
  C[1]=  ik/(4.0*M_PI);
  C[2]= ik2/(8.0*M_PI);
  C[3]= ik3/(24.0*M_PI);

  D[0]= -1.0/(4.0*M_PI);
  D[1]=  ik2/(8.0*M_PI);
  D[2]=  ik3/(12.0*M_PI);
  D[3]=  ik4/(24.0*M_PI);
    
  if (NeedGC)
   { 
     Gab[GCME_GC] +=  C[0]*FIBBIs[ FIBBI_PEFIE1_RM1]
                     +C[1]*FIBBIs[ FIBBI_PEFIE1_R0 ]
                     +C[2]*FIBBIs[ FIBBI_PEFIE1_R1 ]
                     +C[3]*FIBBIs[ FIBBI_PEFIE1_R2 ];

     Gab[GCME_GC] -= OOK2* ( C[0]*FIBBIs[ FIBBI_PEFIE2_RM1]
                            +C[1]*FIBBIs[ FIBBI_PEFIE2_R0 ]
                            +C[2]*FIBBIs[ FIBBI_PEFIE2_R1 ]
                            +C[3]*FIBBIs[ FIBBI_PEFIE2_R2 ]);

   ikCab[GCME_GC] +=  D[0]*FIBBIs[ FIBBI_PMFIE_RM3 ]
                     +D[1]*FIBBIs[ FIBBI_PMFIE_RM1 ]
                     +D[2]*FIBBIs[ FIBBI_PMFIE_R0  ]
                     +D[3]*FIBBIs[ FIBBI_PMFIE_R1  ];
   };

  if (NeedForce)
   { 
     // combinations of real(k), imag(k) needed to
     // assemble singular contributions to imag( \partial_i ikC ) 
     double kr=real(k), kr2=kr*kr, kr4=kr2*kr2;
     double ki=imag(k), ki2=ki*ki, ki4=ki2*ki2;
     double krki=kr*ki, Term2=kr*(kr4 - 10.0*ki2*kr2 + 5.0*ki4);

     for(int Mu=0; Mu<3; Mu++)
      {
         Gab[GCME_FX + Mu]
           +=  D[0]*(        FIBBIs[FIBBI_PEFIE1_RX_RM3  + Mu] 
                      - OOK2*FIBBIs[FIBBI_PEFIE2_RX_RM3  + Mu] )
              +D[1]*(        FIBBIs[FIBBI_PEFIE1_RX_RM1  + Mu]
                      - OOK2*FIBBIs[FIBBI_PEFIE2_RX_RM1  + Mu] )
              +D[2]*(        FIBBIs[FIBBI_PEFIE1_RX_R0   + Mu]
                      - OOK2*FIBBIs[FIBBI_PEFIE2_RX_R0   + Mu] )
              +D[3]*(        FIBBIs[FIBBI_PEFIE1_RX_R1   + Mu]
                      - OOK2*FIBBIs[FIBBI_PEFIE2_RX_R1   + Mu] );

       ikCab[GCME_FX + Mu]
           +=II*( imag( D[1]*FIBBIs[FIBBI_BXBX_RM1 + Mu]
                       +D[2]*FIBBIs[FIBBI_BXBX_R0  + Mu]
                       +D[3]*FIBBIs[FIBBI_BXBX_R1  + Mu]
                      )
                  + krki*FIBBIs[FIBBI_PMFIE_RX_RM3 + Mu]/(4.0*M_PI)
                  + Term2*FIBBIs[FIBBI_PMFIE_RX_R0 + Mu]/(24.0*M_PI)
                );
      };
   };

}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
static int TDMaxEval=1000;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeFIBBIData(RWGSurface *Sa, int nea,
                      RWGSurface *Sb, int neb,
                      double *FIBBIs)
{
// FIXME
  bool NeedForce=true;

  int NumFIBBIs = NeedForce ? NUMFIBBIS : 12;
  memset(FIBBIs, 0, NumFIBBIs*sizeof(double));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  GetGCMEArgStruct MyGetGCMEArgs, *GetGCMEArgs=&MyGetGCMEArgs;
  InitGetGCMEArgs(GetGCMEArgs);
  GetGCMEData MyGetGCMEData, *Data=&MyGetGCMEData;
  Data->Args=GetGCMEArgs;
  Data->GetFIBBIs=true;
  Data->Args->NeedGC    = true;
  Data->Args->NeedForce = NeedForce;
  Data->ZHatAP = Sa->Panels[ Sa->Edges[nea]->iPPanel ]->ZHat;
  Data->ZHatAM = Sa->Panels[ Sa->Edges[nea]->iMPanel ]->ZHat;
  Data->ZHatBP = Sb->Panels[ Sb->Edges[neb]->iPPanel ]->ZHat;
  Data->ZHatBM = Sb->Panels[ Sb->Edges[neb]->iMPanel ]->ZHat;

  int ncv = AssessBFPair(Sa, nea, Sb, neb);

  memset(FIBBIs, 0, NUMFIBBIS*sizeof(double));
  if (ncv>=1)
   { 
     double Error[NUMFIBBIS];
     GetBFBFCubatureTD(Sa, nea, Sb, neb,
                       GCMEIntegrand, (void *)Data, NUMFIBBIS,
                       FIBBIs, Error, TDMaxEval);
   }
  else
   { int Order=9;
     GetBFBFCubature2(Sa, nea, Sb, neb,
                      GCMEIntegrand, (void *)Data, NUMFIBBIS,
                      Order, FIBBIs);
   };
  memset(FIBBIs, 0, 12*sizeof(double));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double *Va[3], *Vb[3];
  TaylorDuffyArgStruct TDArgStruct, *TDArgs=&TDArgStruct;
  InitTaylorDuffyArgs(TDArgs);
  int PIndex[12] = { TD_PMCHWG1, TD_PMCHWG1, TD_PMCHWG1, TD_PMCHWG1,
                     TD_UNITY, TD_UNITY, TD_UNITY, TD_UNITY,
                     TD_PMCHWC,  TD_PMCHWC,  TD_PMCHWC,  TD_PMCHWC};
  int KIndex[12] = { TD_RP, TD_RP, TD_RP, TD_RP, 
                     TD_RP, TD_RP, TD_RP, TD_RP, 
                     TD_RP, TD_RP, TD_RP, TD_RP};
  cdouble KParam[12]={ -1.0,  0.0, 1.0, 2.0,
                       -1.0,  0.0, 1.0, 2.0,
                       -3.0, -1.0, 0.0, 1.0};
  TDArgs->PIndex=PIndex;
  TDArgs->KIndex=KIndex;
  TDArgs->KParam=KParam;
  TDArgs->AbsTol=0.0;
  TDArgs->RelTol=1.0e-4;
  TDArgs->MaxEval=TDMaxEval;
  RWGEdge *Ea = Sa->Edges[nea];
  RWGEdge *Eb = Sb->Edges[neb];
  double LL = Ea->Length * Eb->Length;
  for(int npa=0; npa<2; npa++)
   for(int npb=0; npb<2; npb++)
    { 
      int iQa = (npa==0 ? Ea->iQP : Ea->iQM);
      int iQb = (npb==0 ? Eb->iQP : Eb->iQM);
      if (iQa==-1 || iQb==-1) continue;
      double *Qa = Sa->Vertices + 3*iQa;
      double *Qb = Sb->Vertices + 3*iQb;
      int ipa = (npa==0 ? Ea->iPPanel : Ea->iMPanel);
      int ipb = (npb==0 ? Eb->iPPanel : Eb->iMPanel);
      double rRel;
      ncv=AssessPanelPair(Sa,ipa,Sb,ipb,&rRel,Va,Vb);

      if (ncv==0)
       { 
         double PPContributions[NUMFIBBIS];
         int IDim=12;
         int Order=20;
         GetBFBFCubature2(Sa, nea, Sb, neb,
                          GCMEIntegrand, (void *)Data, IDim,
                          Order, PPContributions, npa, npb);
         for(int n=0; n<IDim; n++)
          FIBBIs[n] += PPContributions[n];
       }
      else 
       { TDArgs->WhichCase=ncv;
         TDArgs->NumPKs = (ncv==3) ? 8 : 12;
         TDArgs->V1=Va[0];
         TDArgs->V2=Va[1];
         TDArgs->V3=Va[2];
         TDArgs->V2P=Vb[1];
         TDArgs->V3P=Vb[2];
         TDArgs->Q  = Qa;
         TDArgs->QP = Qb;

         cdouble PPContributions[12], Error[12];
         memset(PPContributions, 0, 12*sizeof(cdouble));
         TDArgs->Result = PPContributions;
         TDArgs->Error  = Error;
         TaylorDuffy(TDArgs);

         PPContributions[4]*=4.0;
         PPContributions[5]*=4.0;
         PPContributions[6]*=4.0;
         PPContributions[7]*=4.0;
         double Sign = (npa==npb) ? 1.0 : -1.0;
         for(int n=0; n<12; n++)
          FIBBIs[n] += Sign*4.0*M_PI*LL*real(PPContributions[n]);
       };
    };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void dGSabIntegrand(double xA[3], double bA[3], double DivbA,
                    double xB[3], double bB[3], double DivbB,
                    void *UserData, double Weight, double *Integral)
{
  GetGCMEData *Data = (GetGCMEData *)UserData;

  double R[3];
  VecSub(xA, xB, R);
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  double r=sqrt(r2);
  if (r==0.0) return;

  double RPA[3], RPB[3];
  double *ZHatA = (DivbA > 0.0 ? Data->ZHatAP : Data->ZHatAM);
  double *ZHatB = (DivbB > 0.0 ? Data->ZHatBP : Data->ZHatBM);
  RPA[0]=RPA[1]=RPA[2]=RPB[0]=RPB[1]=RPB[2]=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { RPA[Mu] += (1.0 - ZHatA[Mu]*ZHatA[Nu])*R[Nu];
      RPB[Mu] += (1.0 - ZHatB[Mu]*ZHatB[Nu])*R[Nu];
    };

  double ScalarFactor = DivbA*DivbB/(4.0*M_PI*r*r2);

  Integral[0] += Weight*RPA[0]*ScalarFactor;
  Integral[1] += Weight*RPA[1]*ScalarFactor;
  Integral[2] += Weight*RPA[2]*ScalarFactor;
  Integral[3] += Weight*RPB[0]*ScalarFactor;
  Integral[4] += Weight*RPB[1]*ScalarFactor;
  Integral[5] += Weight*RPB[2]*ScalarFactor;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetGCMatrixElements(RWGGeometry *G, GetGCMEArgStruct *Args,
                         int nea, int neb,
                         cdouble Gab[2][NUMGCMES],
                         cdouble ikCab[2][NUMGCMES])
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int nsa        = Args->nsa;
  int nsb        = Args->nsb;
  int NumRegions = Args->NumRegions;

  RWGSurface  *Sa = G->Surfaces[nsa];
  RWGSurface  *Sb = G->Surfaces[nsb];

  /***************************************************************/
  /* figure out the dimension of the integrand vector            */
  /***************************************************************/
  bool NeedGC     = (Args->NeedGC || Args->NeedDW);
  bool NeedDW     = Args->NeedDW;
  bool NeedForce  = Args->NeedForce;
  bool NeedTorque = Args->NeedTorque;

  int NumRawMEs = 0;
  if (NeedGC)        NumRawMEs+=3;
  if (NeedDW)        NumRawMEs+=2;
  if (NeedForce)     NumRawMEs+=6;
  if (NeedTorque)    NumRawMEs+=6;

  int IDim = NumRegions * NumRawMEs;
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double rRel;
  int ncv=AssessBFPair(Sa, nea, Sb, neb, &rRel);
  bool DeSingularize = (ncv > 0);
  if (Args->ForceDeSingularize) DeSingularize=true;
  if (Args->DoNotDeSingularize) DeSingularize=false;
  int Order = (DeSingularize || rRel>4.0) ? 4 : (rRel>1.0 ? 9 : 25);
  if (Args->ForceOrder!=-1)
   Order=Args->ForceOrder;

  GetGCMEData MyData, *Data = &MyData;
  Data->Args          = Args;
  Data->GetFIBBIs     = false;
  Data->DeSingularize = DeSingularize;
  Data->ZHatAP        = Sa->Panels[ Sa->Edges[nea]->iPPanel ]->ZHat;
  Data->ZHatAM        = Sa->Panels[ Sa->Edges[nea]->iMPanel ]->ZHat;
  Data->ZHatBP        = Sb->Panels[ Sb->Edges[neb]->iPPanel ]->ZHat;
  Data->ZHatBM        = Sb->Panels[ Sb->Edges[neb]->iMPanel ]->ZHat;
  cdouble RawMEs[34]; 
  double Error[2*34];
  if (Order==0)
   GetBFBFCubatureTD(G, nsa, nea, nsb, neb, GCMEIntegrand,
                     (void *)Data, 2*IDim, (double *)RawMEs, Error);
  else
   GetBFBFCubature2(G, nsa, nea, nsb, neb, GCMEIntegrand,
                    (void *)Data, 2*IDim, Order, (double *)RawMEs);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double FIBBIs[NUMFIBBIS];
  if (DeSingularize)
   GetFIBBIData(/*Args->FIBBICache*/0, Sa, nea, Sb, neb, FIBBIs);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int nr=0, nrme=0; nr<NumRegions; nr++)
   {
     cdouble k=Args->k[nr], PEFIE1, PEFIE2;

     if (NeedGC)
      { PEFIE1             = RawMEs[nrme+0];
        PEFIE2             = RawMEs[nrme+1];
          Gab[nr][GCME_GC] = PEFIE1 - PEFIE2;
        ikCab[nr][GCME_GC] = RawMEs[nrme+2];
        nrme+=3;
      };

     if (NeedDW)
      {   Gab[nr][GCME_DW] = RawMEs[nrme+0] + 2.0*PEFIE2/k;
        ikCab[nr][GCME_DW] = RawMEs[nrme+1];
        nrme+=2;
      };

     if (NeedForce)
      { for(int Mu=0; Mu<3; Mu++)
         {   Gab[nr][GCME_FX + Mu] = RawMEs[nrme + 0 + Mu];
           ikCab[nr][GCME_FX + Mu] = RawMEs[nrme + 3 + Mu];
         };
        nrme+=6;
      };

     if (NeedTorque)
      { for(int Mu=0; Mu<3; Mu++)
         {   Gab[nr][GCME_TX + Mu] = RawMEs[nrme + 0 + Mu];
           ikCab[nr][GCME_TX + Mu] = RawMEs[nrme + 3 + Mu];
         };
        nrme+=6;
      };

     if (DeSingularize)
      AddFIBBIContributions(FIBBIs, k, NeedGC, NeedForce, Gab[nr], ikCab[nr]);

   }; // for(int nr=0, nrme=0; nr<NumRegions; nr++)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void InitGetGCMEArgs(GetGCMEArgStruct *Args)
{
  Args->NeedGC=false;
  Args->NeedDW=false;
  Args->NeedForce=false;
  Args->NeedTorque=false;
  Args->GBA[0]=Args->GBA[1]=0;
  Args->XTorque[0]=Args->XTorque[1]=Args->XTorque[2]=0.0;
  Args->Displacement=0;
  Args->FIBBICache=0;
  Args->ForceDeSingularize=false;
  Args->DoNotDeSingularize=false;
  Args->ForceOrder=-1;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
static int init=0;
if (init==0)
 { init=1;
   char *s=getenv("SCUFF_TDMAXEVAL");
   if (s)
    { sscanf(s,"%i",&TDMaxEval);
      Log("Setting TDMaxEval=%i.");
    };
 };
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

}

} // namespace scuff 
