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
 * GetPPIs_BruteForce.cc -- compute panel-panel integrals by brute force
 * 
 * homer reid -- 11/2005 -- 12/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libSGJC.h>
#include <libPolyFit.h>

#include "libscuff.h"
#include "libscuffInternals.h"

using namespace scuff;

#define ABSTOL 1.0e-8    // absolute tolerance
#define RELTOL 1.0e-4    // relative tolerance

#define II cdouble(0.0,1.0)

/***************************************************************/
/* data structure used to pass data to integrand routines      */
/***************************************************************/
typedef struct PPIBFData
 { 
   double *V0, A[3], B[3], *Q;
   double V0P[3], AP[3], BP[3], QP[3];
   cdouble k;
   int NeedGradient;
 } PPIBFData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void PPIBFIntegrand(unsigned ndim, const double *x, void *parms,
                           unsigned nfun, double *fval)
{
  PPIBFData *PPIBFD=(PPIBFData *)parms;
  int NeedGradient=PPIBFD->NeedGradient;
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double u=x[0];
  double v=u*x[1];
  double up=x[2];
  double vp=up*x[3];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double X[3], XmQ[3], XP[3], XPmQP[3], R[3], gXh[3];
  double r, r2, r4, kr;

  memcpy(X,PPIBFD->V0,3*sizeof(double));
  VecPlusEquals(X,u,PPIBFD->A);
  VecPlusEquals(X,v,PPIBFD->B);
  VecSub(X, PPIBFD->Q, XmQ);

  memcpy(XP,PPIBFD->V0P,3*sizeof(double));
  VecPlusEquals(XP,up,PPIBFD->AP);
  VecPlusEquals(XP,vp,PPIBFD->BP);
  VecSub(XP, PPIBFD->QP, XPmQP);

  VecSub(X, XP, R);
  r=VecNorm(R);
  r2=r*r;
  r4=r2*r2;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble ik=II*PPIBFD->k;
  cdouble ik2=ik*ik;
  cdouble ikr=ik*r;
  cdouble Phi=exp(ikr)/(4.0*M_PI*r);
  cdouble Psi=Phi*(ikr - 1.0) / r2;
  cdouble Zeta=Phi*(ikr*ikr - 3.0*ikr + 3.0) / r4;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double hDot=u*up*VecDot(XmQ, XPmQP); 
  double hNabla=u*up*4.0;
  cdouble hPlus=hDot+hNabla/ik2;
  double hTimes=u*up*VecDot( VecCross(XmQ, XPmQP, gXh), R );

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble *zf=(cdouble *)fval;

  zf[0] = hPlus*Phi;
  zf[1] = hTimes*Psi;

  if (NeedGradient)
   { int Mu; 
     for(Mu=0; Mu<3; Mu++)
      { zf[2 + 2*Mu + 0] = R[Mu]*hPlus*Psi;
        zf[2 + 2*Mu + 1] = u*up*gXh[Mu]*Psi + R[Mu]*hTimes*Zeta;
      };
   };

} 

/***************************************************************/
/* compute panel-panel integrals using brute-force technique   */
/* (adaptive cubature over both panels).                       */
/***************************************************************/
void GetPPIs_BruteForce(GetPPIArgStruct *Args, int PlotFits)
{ 
  /***************************************************************/
  /* extract fields from argument structure **********************/
  /***************************************************************/
  RWGSurface *Sa = Args->Sa;
  RWGSurface *Sb = Args->Sb;
  int npa = Args->npa;
  int npb = Args->npb;
  int iQa = Args->iQa;
  int iQb = Args->iQb;
  double *Qa = Sa->Vertices + 3*Sa->Panels[npa]->VI[iQa];
  double *Qb = Sb->Vertices + 3*Sb->Panels[npb]->VI[iQb];
  int NeedGradient= (Args->NumGradientComponents > 0 );
  cdouble k = Args->k;

  double rRel;
  double *Va[3], *Vb[3];
  int ncv=AssessPanelPair(Sa, npa, Sb, npb, &rRel, Va, Vb);

  /***************************************************************/
  /* setup for call to cubature routine    ***********************/
  /***************************************************************/
  PPIBFData MyPPIBFData, *PPIBFD=&MyPPIBFData;
 
  PPIBFD->V0 = Va[0];
  VecSub(Va[1], Va[0], PPIBFD->A);
  VecSub(Va[2], Va[1], PPIBFD->B);
  PPIBFD->Q  = Qa;

  // note that for Vb[0] and Qb we make copies of the 
  // entries (not just the pointers) because we may need
  // to displace them, see below.
  memcpy(PPIBFD->V0P,Vb[0],3*sizeof(double));
  VecSub(Vb[1], Vb[0], PPIBFD->AP);
  VecSub(Vb[2], Vb[1], PPIBFD->BP);
  memcpy(PPIBFD->QP,Qb,3*sizeof(double));
   
  PPIBFD->k = k;
  PPIBFD->NeedGradient = NeedGradient;

  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};

  int fDim=(PPIBFD->NeedGradient) ? 16 : 4;
  double Result[fDim], Error[fDim];

  /***************************************************************/
  /* switch off based on whether or not there are any common     */
  /* vertices                                                    */
  /***************************************************************/
  if (ncv==0)
   {
     /*--------------------------------------------------------------*/
     /* if there are no common vertices then we can just use naive   */
     /* cubature                                                     */
     /*--------------------------------------------------------------*/
 //    adapt_integrate(fDim, PPIBFIntegrand, (void *)PPIBFD, 4, Lower, Upper,
 //                    0, ABSTOL, RELTOL, Result, Error);
       adapt_integrate_log(fDim, PPIBFIntegrand, (void *)PPIBFD, 4, Lower, Upper,
                           0, ABSTOL, RELTOL, Result, Error, "SGJC.log", 15);

     Args->H[0] = cdouble(Result[0], Result[1]);
     Args->H[1] = cdouble(Result[2], Result[3]);
     if (NeedGradient)
      { Args->GradH[0] = cdouble(Result[ 4], Result[ 5]);
        Args->GradH[1] = cdouble(Result[ 6], Result[ 7]);
        Args->GradH[2] = cdouble(Result[ 8], Result[ 9]);
        Args->GradH[3] = cdouble(Result[10], Result[11]);
        Args->GradH[4] = cdouble(Result[12], Result[13]);
        Args->GradH[5] = cdouble(Result[14], Result[15]);
      };

   }
  else
   {
     /*--------------------------------------------------------------*/
     /* if there are common vertices then we estimate the panel-panel*/
     /* integrals using a limiting process in which we displace the  */
     /* second of the two panels through a distance Z in the         */
     /* direction of the panel normal and try to fit to Z==0         */
     /*--------------------------------------------------------------*/
     int nz, NZ=10;
     double Z[NZ], GR[NZ], GI[NZ], CR[NZ], CI[NZ];
     double DeltaZ=Sb->Panels[npb]->Radius/100.0;
     double *ZHat=Sb->Panels[npb]->ZHat;

     for(nz=0; nz<NZ; nz++)
      { 
        Z[nz]=((double)(nz+1))*DeltaZ;
        VecScaleAdd(Vb[0], Z[nz], ZHat, PPIBFD->V0P);
        VecScaleAdd(Qb,    Z[nz], ZHat, PPIBFD->QP);

 //       adapt_integrate(fDim, PPIBFIntegrand, (void *)PPIBFD, 4, Lower, Upper,
 //                       0, ABSTOL, RELTOL, Result, Error);
       adapt_integrate_log(fDim, PPIBFIntegrand, (void *)PPIBFD, 4, Lower, Upper,
                           100000, ABSTOL, RELTOL, Result, Error, "SGJC.log", 15);

        GR[nz]=Result[0];
        GI[nz]=Result[1];
        CR[nz]=Result[2];
        CI[nz]=Result[3];
      };
 
     PolyFit *PF=new PolyFit(Z, GR, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, GR, NZ, 0.0, Z[NZ-1], "real(<fa|G|fb>)");
     real(Args->H[0])=PF->f(0.0);
     delete PF;

     PF=new PolyFit(Z, GI, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, GI, NZ, 0.0, Z[NZ-1], "imag(<fa|G|fb>)");
     imag(Args->H[0])=PF->f(0.0);
     delete PF;

     PF=new PolyFit(Z, CR, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, CR, NZ, 0.0, Z[NZ-1], "imag(<fa|C|fb>)");
     real(Args->H[1])=PF->f(0.0);
     delete PF;

     PF=new PolyFit(Z, CI, NZ, 4);
     if (PlotFits) 
      PF->PlotFit(Z, CI, NZ, 0.0, Z[NZ-1], "imag(<fa|C|fb>)");
     imag(Args->H[1])=PF->f(0.0);
     delete PF;
     
   }; // if (ncv==0 ... else)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetEEIs_BruteForce(GetEEIArgStruct *EEIArgs) 
{
  GetPPIArgStruct MyPPIArgs, *PPIArgs=&MyPPIArgs;
  InitGetPPIArgs(PPIArgs);

  PPIArgs->Sa = EEIArgs->Sa;
  PPIArgs->Sb = EEIArgs->Sb;
  PPIArgs->k  = EEIArgs->k;
  PPIArgs->NumGradientComponents=3;

  RWGEdge *Ea = EEIArgs->Sa->Edges[EEIArgs->nea];
  RWGEdge *Eb = EEIArgs->Sb->Edges[EEIArgs->neb];

  PPIArgs->npa = Ea->iPPanel;
  PPIArgs->iQa = Ea->PIndex;
  PPIArgs->npb = Eb->iPPanel;
  PPIArgs->iQb = Eb->PIndex;
  GetPPIs_BruteForce(PPIArgs, 0);
  EEIArgs->GC[0]     = PPIArgs->H[0];
  EEIArgs->GC[1]     = PPIArgs->H[1];
  EEIArgs->GradGC[0] = PPIArgs->GradH[0];
  EEIArgs->GradGC[1] = PPIArgs->GradH[1];
  EEIArgs->GradGC[2] = PPIArgs->GradH[2];
  EEIArgs->GradGC[3] = PPIArgs->GradH[3];
  EEIArgs->GradGC[4] = PPIArgs->GradH[4];
  EEIArgs->GradGC[5] = PPIArgs->GradH[5];

  PPIArgs->npa = Ea->iPPanel;
  PPIArgs->iQa = Ea->PIndex;
  PPIArgs->npb = Eb->iMPanel;
  PPIArgs->iQb = Eb->MIndex;
  GetPPIs_BruteForce(PPIArgs, 0);
  EEIArgs->GC[0]     -= PPIArgs->H[0];
  EEIArgs->GC[1]     -= PPIArgs->H[1];
  EEIArgs->GradGC[0] -= PPIArgs->GradH[0];
  EEIArgs->GradGC[1] -= PPIArgs->GradH[1];
  EEIArgs->GradGC[2] -= PPIArgs->GradH[2];
  EEIArgs->GradGC[3] -= PPIArgs->GradH[3];
  EEIArgs->GradGC[4] -= PPIArgs->GradH[4];
  EEIArgs->GradGC[5] -= PPIArgs->GradH[5];

  PPIArgs->npa = Ea->iMPanel;
  PPIArgs->iQa = Ea->MIndex;
  PPIArgs->npb = Eb->iPPanel;
  PPIArgs->iQb = Eb->PIndex;
  GetPPIs_BruteForce(PPIArgs, 0);
  EEIArgs->GC[0]     -= PPIArgs->H[0];
  EEIArgs->GC[1]     -= PPIArgs->H[1];
  EEIArgs->GradGC[0] -= PPIArgs->GradH[0];
  EEIArgs->GradGC[1] -= PPIArgs->GradH[1];
  EEIArgs->GradGC[2] -= PPIArgs->GradH[2];
  EEIArgs->GradGC[3] -= PPIArgs->GradH[3];
  EEIArgs->GradGC[4] -= PPIArgs->GradH[4];
  EEIArgs->GradGC[5] -= PPIArgs->GradH[5];

  PPIArgs->npa = Ea->iMPanel;
  PPIArgs->iQa = Ea->MIndex;
  PPIArgs->npb = Eb->iMPanel;
  PPIArgs->iQb = Eb->MIndex;
  GetPPIs_BruteForce(PPIArgs, 0);
  EEIArgs->GC[0]     += PPIArgs->H[0];
  EEIArgs->GC[1]     += PPIArgs->H[1];
  EEIArgs->GradGC[0] += PPIArgs->GradH[0];
  EEIArgs->GradGC[1] += PPIArgs->GradH[1];
  EEIArgs->GradGC[2] += PPIArgs->GradH[2];
  EEIArgs->GradGC[3] += PPIArgs->GradH[3];
  EEIArgs->GradGC[4] += PPIArgs->GradH[4];
  EEIArgs->GradGC[5] += PPIArgs->GradH[5];

  double GPreFac=Ea->Length * Eb->Length;
  cdouble CPreFac=GPreFac / (II*PPIArgs->k);
  EEIArgs->GC[0] *= GPreFac;
  EEIArgs->GC[1] *= CPreFac;
  EEIArgs->GradGC[0] *= GPreFac;
  EEIArgs->GradGC[1] *= CPreFac;
  EEIArgs->GradGC[2] *= GPreFac;
  EEIArgs->GradGC[3] *= CPreFac;
  EEIArgs->GradGC[4] *= GPreFac;
  EEIArgs->GradGC[5] *= CPreFac;
}
