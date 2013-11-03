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
 * GetPhiE.cc  -- routine for computing electrostatic potential and field
 *                due to a constant charge density on a single triangle,
 *                following these references:
 *                (1) Graglia, IEEE Trans. Ant. Prop. *41* 1448 (1993)
 *                (2) Wilton et al, IEEE Trans. Ant. Prop. *32* 276 (1984)
 *
 * homer reid  -- 5/2013 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libTriInt.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "TaylorDuffy.h"
#include "SSSolver.h"

int ForcePPIMethod=0;
namespace scuff {

/***********************************************************************/
/* compute panel-panel integrals using cubature over both panels.      */
/* note: 'CC' stands for 'cubature-cubature.'                          */
/***********************************************************************/
double GetPPI_CC(double *Va[3], double *Vb[3], double *nHat, int Order, int WhichIntegral)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double *V0, A[3], B[3]; 
  double *V0P, AP[3], BP[3]; 

  V0=Va[0];
  VecSub(Va[1], Va[0], A);
  VecSub(Va[2], Va[1], B);

  V0P=Vb[0];
  VecSub(Vb[1], Vb[0], AP);
  VecSub(Vb[2], Vb[1], BP);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double *TCR;
  int NumPts;
  if (Order==20)
   TCR=GetTCR(20, &NumPts);
  else
   TCR=GetTCR(4, &NumPts);
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Result=0.0;
  for(int np=0, ncp=0; np<NumPts; np++) 
   { 
     double u=TCR[ncp++]; 
     double v=TCR[ncp++]; 
     double w=TCR[ncp++];
     double X[3];
     for(int Mu=0; Mu<3; Mu++)
      X[Mu] = V0[Mu] + u*A[Mu] + v*B[Mu];

     for(int npp=0, ncpp=0; npp<NumPts; npp++)
      { double up=TCR[ncpp++]; 
        double vp=TCR[ncpp++]; 
        double wp=TCR[ncpp++];
        double XP[3], R[3], r2=0.0, r;
        for(int Mu=0; Mu<3; Mu++)
         { XP[Mu] = V0P[Mu] + up*AP[Mu] + vp*BP[Mu];
           R[Mu] = X[Mu] - XP[Mu];
           r2 += R[Mu]*R[Mu];
         };
        r=sqrt(r2);

        if (WhichIntegral==0)
         Result += w*wp/r;
        else
         Result += w*wp*(nHat[0]*R[0] + nHat[1]*R[1] + nHat[2]*R[2])/(r*r2);
      };
   };

  return Result/(4.0*M_PI);

}

/***********************************************************************/
/* compute panel-panel integrals using cubature over the first panel   */
/* but the exact potential/field calculation for the second panel.     */
/* note: 'ce' stands for 'cubature/exact'                              */
/***********************************************************************/
double GetPPI_CE(SSSolver *SSS, double *Va[3], double *nHat, int ns2, int np2, 
                 int Order, int WhichIntegral)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double *V0, A[3], B[3]; 

  V0=Va[0];
  VecSub(Va[1], Va[0], A);
  VecSub(Va[2], Va[1], B);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double *TCR;
  int NumPts;
  if (Order==20)
   TCR=GetTCR(20, &NumPts);
  else
   TCR=GetTCR(4, &NumPts);
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Result=0.0;
  for(int np=0, ncp=0; np<NumPts; np++) 
   { 
     double u=TCR[ncp++]; 
     double v=TCR[ncp++]; 
     double w=TCR[ncp++];
     double X[3];
     for(int Mu=0; Mu<3; Mu++)
      X[Mu] = V0[Mu] + u*A[Mu] + v*B[Mu];

     double PhiE[4];
     SSS->GetPhiE(ns2, np2, X, PhiE);

     if (WhichIntegral==0)
      Result += w*PhiE[0]/(4.0*M_PI);
     else
      Result += w*(nHat[0]*PhiE[1+0] + nHat[1]*PhiE[1+1] + nHat[2]*PhiE[1+2])/(4.0*M_PI);
   };

  return Result;

}

/***********************************************************************/
/* WhichIntegral:                                                      */
/*  0 --> average over Panel A of potential due to panel B             */
/*  1 --> average over Panel A of E-field due to panel B dotted into   */
/*        normal to Panel A                                            */
/***********************************************************************/
#define PPIMETHOD_CE20 1
#define PPIMETHOD_CE4  2
#define PPIMETHOD_CC20 3
#define PPIMETHOD_CC4  4
double SSSolver::GetPPI(RWGSurface *Sa, int npa, 
                        RWGSurface *Sb, int npb, 
                        int WhichIntegral)
{
  double rRel, *Va[3], *Vb[3];

  int ncv=AssessPanelPair(Sa,npa,Sb,npb,&rRel,Va,Vb);
  double *nHat = Sa->Panels[npa]->ZHat;
  double Jacobian = 4.0 * Sa->Panels[npa]->Area * Sb->Panels[npb]->Area;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (ForcePPIMethod==PPIMETHOD_CE20)
   { Jacobian=2.0*Sa->Panels[npa]->Area;
     return Jacobian*GetPPI_CE(this, Va, nHat, Sb->Index, npb, 20, WhichIntegral);
   }
  else if (ForcePPIMethod==PPIMETHOD_CE4)
   { Jacobian=2.0*Sa->Panels[npa]->Area;
     return Jacobian*GetPPI_CE(this ,Va, nHat, Sb->Index, npb, 4, WhichIntegral);
   }
  else if (ForcePPIMethod==PPIMETHOD_CC20)
   return Jacobian*GetPPI_CC(Va, Vb, nHat, 20, WhichIntegral);
  else if (ForcePPIMethod==PPIMETHOD_CC4)
   return Jacobian*GetPPI_CC(Va, Vb, nHat, 4, WhichIntegral);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (rRel>10.0)
   { // distant panels 
     return Jacobian*GetPPI_CC(Va, Vb, nHat, 4, WhichIntegral);
   }
  else if (ncv==0)
   { // nearby panels, no common vertices 
     return Jacobian*GetPPI_CC(Va, Vb, nHat, 20, WhichIntegral);
   }
  else
   { 
     // common vertices are present; use Taylor-Duffy
     int PIndex, KIndex;
     cdouble KParam;
     cdouble Result, Error;
     TaylorDuffyArgStruct TDArgs, *Args=&TDArgs;
     InitTaylorDuffyArgs(Args);
     Args->WhichCase = ncv;
     Args->NumPKs    = 1;
     Args->PIndex    = &PIndex;
     Args->KIndex    = &KIndex;
     Args->KParam    = &KParam;
     Args->V1        = Va[0];
     Args->V2        = Va[1];
     Args->V3        = Va[2];
     Args->V2P       = Vb[1];
     Args->V3P       = Vb[2];
     Args->nHat      = nHat;
     Args->Result    = &Result;
     Args->Error     = &Error; 

     KIndex = TD_RP;
     if (WhichIntegral==0)
      { PIndex = TD_UNITY;
        KParam = -1.0;
      }
     else // WhichIntegral==1
      { PIndex = TD_RNORMAL;
        KParam = -3.0;
      };

     TaylorDuffy(Args);
     return Jacobian*real(Result);
   };
  
}

} // namespace scuff
