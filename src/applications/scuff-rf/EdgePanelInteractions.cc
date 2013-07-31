/* Copyright (C) 2005-2011 M. T. Homer Reid *
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
 * EdgePanelInteractions.cc -- routines for evaluating the interactions
 *                          -- between a single RWG panel and the line charge
 *                          -- associated with a half-RWG basis function     
 * 
 * homer reid               -- 6/2012
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libhrutil.h"
#include "libTriInt.h"
#include "libSGJC.h"

#include "libscuff.h"
//#include "libscuffInternals.h"

namespace scuff{

#define ABSTOL 1.0e-12
#define RELTOL 1.0e-8
#define MAXEVALS 100000

int nCalls;

cdouble ExpRelV3P0(int n, cdouble Z);

#define II cdouble(0.0,1.0)

/*--------------------------------------------------------------*/
/* 5- and 10-pt legendre quadrature rules for the interval [0,1]*/
/* format: QR[2*n+0, 2*n+1] = nth quadrature point, weight      */
/*--------------------------------------------------------------*/
double QR5[]=
 { 4.691007703066802e-02, 1.184634425280909e-01, 
   2.307653449471584e-01, 2.393143352496832e-01, 
   5.000000000000000e-01, 2.844444444444444e-01,
   7.692346550528415e-01, 2.393143352496832e-01, 
   9.530899229693319e-01, 1.184634425280909e-01
 };

double QR10[]=
 { 1.304673574141418e-02, 3.333567215434143e-02, 
   6.746831665550773e-02, 7.472567457529027e-02, 
   1.602952158504878e-01, 1.095431812579910e-01,
   2.833023029353764e-01, 1.346333596549959e-01, 
   4.255628305091844e-01, 1.477621123573765e-01, 
   5.744371694908156e-01, 7.166976970646236e-01, 
   8.397047841495122e-01, 9.325316833444923e-01,
   9.869532642585859e-01, 1.477621123573765e-01,
   1.346333596549959e-01, 1.095431812579910e-01, 
   7.472567457529027e-02, 3.333567215434143e-02
 };

double QR20[]=
 { 3.435700407452558e-03, 8.807003569575289e-03,
   1.801403636104310e-02, 2.030071490019352e-02,
   4.388278587433703e-02, 3.133602416705452e-02,
   8.044151408889055e-02, 4.163837078835237e-02,
   1.268340467699246e-01, 5.096505990861661e-02,
   1.819731596367425e-01, 5.909726598075887e-02,
   2.445664990245864e-01, 6.584431922458830e-02,
   3.131469556422902e-01, 7.104805465919108e-02,
   3.861070744291775e-01, 7.458649323630188e-02,
   4.617367394332513e-01, 7.637669356536299e-02,
   5.382632605667487e-01, 7.637669356536299e-02,
   6.138929255708225e-01, 7.458649323630188e-02,
   6.868530443577098e-01, 7.104805465919108e-02,
   7.554335009754136e-01, 6.584431922458830e-02,
   8.180268403632576e-01, 5.909726598075887e-02,
   8.731659532300754e-01, 5.096505990861661e-02,
   9.195584859111094e-01, 4.163837078835237e-02,
   9.561172141256630e-01, 3.133602416705452e-02,
   9.819859636389570e-01, 2.030071490019352e-02,
   9.965642995925474e-01, 8.807003569575289e-03
};

int QROrder=-1, TCROrder=-1;

/*--------------------------------------------------------------*/
/*- PART 1: Routine to compute edge-panel interaction using     */
/*-         fixed-order numerical cubature.                     */
/*-                                                             */
/*- PV[0][0..2] = cartesian coordinates of panel vertex 1       */
/*- PV[1][0..2] = cartesian coordinates of panel vertex 2       */
/*- PV[2][0..2] = cartesian coordinates of panel vertex 3       */
/*- EV[0][0..2] = cartesian coordinates of edge vertex 1        */
/*- EV[1][0..2] = cartesian coordinates of edge vertex 2        */
/*-                                                             */
/*--------------------------------------------------------------*/
cdouble GetEPI_Cubature(double **PV, double **EV, cdouble K)
{ 

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double *V0, A[3], B[3], *V0P, L[3];
  V0=PV[0];
  VecSub(PV[1], PV[0], A);
  VecSub(PV[2], PV[0], B);
  V0P=EV[0];
  VecSub(EV[1],EV[0],L);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double *TCR;  // triangle cubature rule 
  int NumTCPts;
  TCR=GetTCR(4, &NumTCPts);

if (TCROrder!=-1)
 TCR=GetTCR(TCROrder, &NumTCPts);
if (TCR==0) ErrExit("unknown TCR order"); 

  double *QR=QR5;   // quadrature rule 
  int NumQPts=5;

if (QROrder==5)
 QR=QR5, NumQPts=5;
else if (QROrder==10)
 QR=QR10, NumQPts=10;
else if (QROrder==20)
 QR=QR20, NumQPts=20;

  /***************************************************************/
  /* outer loop (cubature over triangle) *************************/
  /***************************************************************/
  int np, ncp, npp, ncpp;
  int Mu;
  double u, v, w, up, wp;
  double X[3], XP[3], R[3];
  double r;
  cdouble Phi;
  cdouble GOuter, GInner;
  GOuter=0.0;
  for(np=ncp=0; np<NumTCPts; np++) 
   { 
     u=TCR[ncp++]; v=TCR[ncp++]; w=TCR[ncp++];

     /***************************************************************/
     /* set X and F=X-Q *********************************************/
     /***************************************************************/
     for(Mu=0; Mu<3; Mu++)
      X[Mu] = V0[Mu] + u*A[Mu] + v*B[Mu];

     /***************************************************************/
     /* inner loop (quadrature over edge)                           */
     /***************************************************************/
     GInner=0.0;
     for(npp=ncpp=0; npp<NumQPts; npp++)
      { 
        up=QR[ncpp++]; wp=QR[ncpp++];

        /***************************************************************/ 
        /* set XP and FP=XP-QP *****************************************/
        /***************************************************************/
        for(Mu=0; Mu<3; Mu++)
         { XP[Mu] = V0P[Mu] + up*L[Mu];
           R[Mu] = X[Mu] - XP[Mu];
         };
      
        /***************************************************************/
        /* inner integrand  ********************************************/
        /***************************************************************/
        r=VecNorm(R);
        Phi = exp(II*K*r)/(4.0*M_PI*r);
        GInner += wp*Phi;

      }; /* for(npp=ncpp=0; npp<NumPts; npp++) */

     /*--------------------------------------------------------------*/
     /*- accumulate contributions to outer integral                  */
     /*--------------------------------------------------------------*/
     GOuter+=w*GInner;

   }; // for(np=ncp=0; np<nPts; np++) 

  return GOuter;

}

/***************************************************************/
/* 'edge-panel integrand, taylor-duffy data'                   */
/***************************************************************/
typedef struct EPITDData
{ 
  double A[3], B[3], L[3];
  cdouble K;
  int WhichCase; // = 1 for common-vertex, 2 for common-edge 

} EPITDData;

/***************************************************************/
/* 'edge-panel integral, taylor-duffy integrand'               */
/***************************************************************/
void EPITDIntegrand(unsigned ndim, const double *xy, void *params,
                    unsigned fdim, double *fval)
{

nCalls++;

  (void) ndim;
  (void) fdim;

  EPITDData *Data= (EPITDData *)params;

  double *A     = Data->A;
  double *B     = Data->B;
  double *L     = Data->L;
  cdouble IK    = II*Data->K;
  int WhichCase = Data->WhichCase;

  double RR1[3], RR2[3], RR3[3], R1, R2, R3;
  cdouble *zfval=(cdouble *)fval;
  int p;

  if (WhichCase==1) // common-vertex case 
   { 
     double x=xy[0];
     double y=xy[1];

     for(R1=R2=R3=0.0, p=0; p<3; p++)
      { RR1[p] = x*A[p] + x*y*B[p]       - L[p];
        RR2[p] =   A[p] + (1.0-x)*y*B[p] - x*L[p];
        RR3[p] =   A[p] + (1.0-x*y)*B[p] - x*L[p];
        R1 += RR1[p]*RR1[p];
        R2 += RR2[p]*RR2[p];
        R3 += RR3[p]*RR3[p];
      };
     R1=sqrt(R1);
     R2=sqrt(R2);
     R3=sqrt(R3);
  
     zfval[0] =        x*exp(IK*R1)*ExpRelV3P0(2,-IK*R1) / (8.0*M_PI*R1)
                +(1.0-x)*exp(IK*R2)*ExpRelV3P0(2,-IK*R2) / (8.0*M_PI*R2)
                      +x*exp(IK*R3)*ExpRelV3P0(2,-IK*R3) / (8.0*M_PI*R3);
   }
  else
   { 
     double y=xy[0];

     for(R1=R2=R3=0.0, p=0; p<3; p++)
      { RR1[p] = (y-1.0)*A[p] + y*B[p];
        RR2[p] =       y*A[p] +   B[p];
        RR3[p] =         A[p] + y*B[p];
        R1 += RR1[p]*RR1[p];
        R2 += RR2[p]*RR2[p];
        R3 += RR3[p]*RR3[p];
      };
     R1=sqrt(R1);
     R2=sqrt(R2);
     R3=sqrt(R3);
  
     zfval[0] = exp(IK*R1)*(2.0*ExpRelV3P0(1,-IK*R1) - 0.5*(y+1.0)*ExpRelV3P0(2,-IK*R1)) / (8.0*M_PI*R1)
               +exp(IK*R2)*(2.0*ExpRelV3P0(1,-IK*R2) -             ExpRelV3P0(2,-IK*R2)) / (8.0*M_PI*R2)
               +exp(IK*R3)*(2.0*ExpRelV3P0(1,-IK*R3) -             ExpRelV3P0(2,-IK*R3)) / (8.0*M_PI*R3);
   };
}

 
/***************************************************************/
/* PV[0][0..2] = cartesian coordinates of panel vertex 1       */
/* PV[1][0..2] = cartesian coordinates of panel vertex 2       */
/* PV[2][0..2] = cartesian coordinates of panel vertex 3       */
/*                                                             */
/* EVL[0][0..2] = cartesian coordinates of edge vertex 1       */
/* EV[1][0..2] = cartesian coordinates of edge vertex 2        */
/*                                                             */
/* Note:  What this routine computes is simply the unadorned   */
/*        integral                                             */
/*                                                             */
/*        \int_{panel} dx \int_{edge} dy G(x-y)                */
/*                                                             */
/*       where G(r)=e^{ikr}/(4*pi*r).                          */
/*                                                             */
/*       This means that the calling routine must multiply the */
/*       return value by a prefactor of                        */
/*                                                             */
/*         +2*l_\alpha l_\beta/k^2                             */
/*                                                             */
/*       to get the value of the <f|G|h> inner product.        */
/***************************************************************/
cdouble GetEdgePanelInteraction(double **PV, double **EV, cdouble K)
{ 
  /*--------------------------------------------------------------*/
  /*- look for common vertices -----------------------------------*/
  /*--------------------------------------------------------------*/
  int EICV[2]; // 'edge index of common vertex' 
  int PICV[2]; // 'panel index of common vertex'
  int ncv=0;   // 'number of common vertices'
  for (int i=0; i<2; i++)
   for (int j=0; j<3; j++)
    if ( VecEqualFloat(EV[i], PV[j]) )
     { if (ncv==2) ErrExit("%s:%i: internal error \n",__FILE__,__LINE__);
       EICV[ncv]=i;
       PICV[ncv]=j;
       ncv++;
     };

  /*--------------------------------------------------------------*/
  /*- use simple cubature if there are no common vertices --------*/
  /*--------------------------------------------------------------*/
  if ( ncv == 0 )
   return GetEPI_Cubature(PV, EV, K);

  /*--------------------------------------------------------------*/
  /*- use taylor-duffy method if there are common vertices -------*/
  /*--------------------------------------------------------------*/
  EPITDData MyData, *Data=&MyData;

  int i, ip1, ip2, j, jp1;
  // define the A and B vectors from the memo
  if (ncv==1)
   { i=PICV[0], ip1=(i+1)%3, ip2=(i+2)%3;
     j=EICV[0], jp1=(j+1)%2;
     VecSub(PV[ip1], PV[i],   Data->A);
     VecSub(PV[ip2], PV[ip1], Data->B);
     VecSub(EV[jp1], EV[j],   Data->L);
   }
  else // ncv==2
   { i=PICV[0], ip1=PICV[1], ip2=3-i-ip1;
     VecSub(PV[ip1], PV[i],   Data->A);
     VecSub(PV[ip2], PV[ip1], Data->B);
   };
   
  Data->K=K;
  Data->WhichCase = ncv;
 
  double Lower[2] = { 0.0, 0.0 };
  double Upper[2] = { 1.0, 1.0 };
  cdouble G, Error;
  // note: the dimension of the integral is 2 for ncv==1, 1 for ncv==2
  adapt_integrate(2, EPITDIntegrand, (void *)Data, 3-ncv, Lower, Upper,
		  MAXEVALS, ABSTOL, RELTOL, (double *)&G, (double *)&Error);

  return G;

}

} // namespace scuff;
