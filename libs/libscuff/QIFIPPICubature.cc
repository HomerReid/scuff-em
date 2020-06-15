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
 * QIFIPPICubature.cc  -- libscuff routines for computing Q-independent 
 *                     -- frequency-independent panel-panel integrals 
 *                     -- (QIFIPPIs) using simple cubature 
 * 
 * homer reid   -- 11/2005 -- 11/2011
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libSGJC.h>

#include "libscuff.h"
#include "libscuffInternals.h"

namespace scuff {

#define ABSTOL   1.0e-10
#define RELTOL   1.0e-6
#define MAXEVALS 1000000

#define _ONE 0
#define _UP  1
#define _VP  2
#define _U   3
#define _UUP 4
#define _UVP 5
#define _V   6
#define _VUP 7
#define _VVP 8

#define uvupvpR0_ONE (1.0/4.0)
#define uvupvpR0_UP  (1.0/6.0)
#define uvupvpR0_VP  (1.0/12.0)
#define uvupvpR0_U   (1.0/6.0)
#define uvupvpR0_UUP (1.0/9.0)
#define uvupvpR0_UVP (1.0/18.0)
#define uvupvpR0_V   (1.0/12.0)
#define uvupvpR0_VUP (1.0/18.0)
#define uvupvpR0_VVP (1.0/36.0)

void ComputeQIFIPPIData_TaylorDuffy(double *V1, double *V2, double *V3, 
                                    double *V2P, double *V3P, 
                                    QIFIPPIData *QIFD);

void ComputeQIFIPPIData_TaylorDuffyV2P0(double *V1, double *V2, double *V3, 
                                        double *V2P, double *V3P, 
                                        QIFIPPIData *QIFD);


/*--------------------------------------------------------------*/
/*- integrand routine used for evaluating FIPPIs by cubature.  -*/
/*- note: CFD = 'compute FIPPI data.'                          -*/
/*--------------------------------------------------------------*/
typedef struct CFDData
 {
   double *V0, A[3], B[3];
   double *V0P, AP[3], BP[3];
   int nCalls;
 } CFDData;

#if 0
 // this is the integrand routine for evaluating the FIPPIs
 // using 4-dimensional cubature; i have since replaced
 // it with the much faster 3-dimensional cubature scheme
 // described in the memo
void CFDIntegrand4D(unsigned ndim, const double *x, void *params,
                    unsigned fdim, double *fval)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  CFDData *CFDD=(CFDData *)params;
  CFDD->nCalls++;

  double *V0=CFDD->V0;
  double *A=CFDD->A;
  double *B=CFDD->B;

  double *V0P=CFDD->V0P;
  double *AP=CFDD->AP;
  double *BP=CFDD->BP;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double u, v, up, vp, uup;
  u=x[0];
  v=u*x[1];
  up=x[2];
  vp=up*x[3];
  uup=u*up; // jacobian

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double X[3], XP[3], R[3], XxXP[3];
  double r, r2=0.0;
  int Mu;
  for (Mu=0; Mu<3; Mu++)
   { X[Mu]  = V0[Mu]  + u*A[Mu]   + v*B[Mu];
     XP[Mu] = V0P[Mu] + up*AP[Mu] + vp*BP[Mu];
     R[Mu]  = X[Mu] - XP[Mu];
     r2    += R[Mu]*R[Mu];
   };
  r=sqrt(r2);
  VecCross(X,XP,XxXP);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nf=0;
  double oor, oor3;

  // we put the jacobian factors into these quantities for convenience
  oor=u*up/r;
  oor3=oor/r2;
  r*=u*up;
  r2*=u*up;

  fval[nf++] = R[0] * oor3;
  fval[nf++] = R[1] * oor3;
  fval[nf++] = R[2] * oor3;
  fval[nf++] = XxXP[0] * oor3;
  fval[nf++] = XxXP[1] * oor3;
  fval[nf++] = XxXP[2] * oor3;

  fval[nf++] = oor;
  fval[nf++] = up*oor;
  fval[nf++] = vp*oor;
  fval[nf++] = u*oor;
  fval[nf++] = u*up*oor;
  fval[nf++] = u*vp*oor;
  fval[nf++] = v*oor;
  fval[nf++] = v*up*oor;
  fval[nf++] = v*vp*oor;

  fval[nf++] = r;
  fval[nf++] = up*r;
  fval[nf++] = vp*r;
  fval[nf++] = u*r;
  fval[nf++] = u*up*r;
  fval[nf++] = u*vp*r;
  fval[nf++] = v*r;
  fval[nf++] = v*up*r;
  fval[nf++] = v*vp*r;

  fval[nf++] = r2;
  fval[nf++] = up*r2;
  fval[nf++] = vp*r2;
  fval[nf++] = u*r2;
  fval[nf++] = u*up*r2;
  fval[nf++] = u*vp*r2;
  fval[nf++] = v*r2;
  fval[nf++] = v*up*r2;
  fval[nf++] = v*vp*r2;

} 
#endif

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void CFDIntegrand3D(unsigned ndim, const double *x, void *params,
                    unsigned fdim, double *fval)
{
  (void) ndim; (void) fdim; // unused;
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  CFDData *CFDD=(CFDData *)params;
  CFDD->nCalls++;

  double *V0=CFDD->V0;
  double *A=CFDD->A;
  double *B=CFDD->B;

  double *V0P=CFDD->V0P;
  double *AP=CFDD->AP;
  double *BP=CFDD->BP;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double u, v, up;
  u=x[0];
  v=u*x[1];
  up=x[2];
  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int Mu;
  double X[3], Y[3], Z[3], BPdY, Y2, a2, a, a3, vp0, vp02, vp03, vp04, b2; 
  a2=BPdY=Y2=0.0;
  for (Mu=0; Mu<3; Mu++)
   { 
     //Y[Mu]  = V0[Mu] - V0P[Mu] + u*A[Mu] + v*B[Mu] - up*AP[Mu];
     X[Mu]  = V0[Mu] + u*A[Mu] + v*B[Mu];
     Z[Mu]  = V0P[Mu] + up*AP[Mu];
     Y[Mu]  = X[Mu] - Z[Mu];
     a2    += BP[Mu]*BP[Mu];
     BPdY  += BP[Mu]*Y[Mu];
     Y2    += Y[Mu]*Y[Mu];
   };
  a=sqrt(a2);
  a3=a*a2;
  vp0=-BPdY/a2;
  vp02=vp0*vp0;
  vp03=vp02*vp0;
  vp04=vp03*vp0;
  b2=Y2/a2 - vp02;

  // the following quantities are the integrals evaluated in Section 10 of the memo.
  // note we put the jacobian factor (u) into these quantities for convenience.
  double OneRM3Int, vpRM3Int;
  double OneRM1Int, vpRM1Int;
  double OneR1Int,  vpR1Int;
  double OneR2Int,  vpR2Int;
  if ( b2 > 1.0e-10)
   { 
     double S1, S2, LogFac, Sum, Sum3, Sum4;
     S1=sqrt( b2 + vp02 );
     S2=sqrt( b2 + (vp0+up)*(vp0+up) );
     LogFac=log( (S2 + (up+vp0)) / (S1 + vp0) );
     Sum=vp0+up;
     Sum3=Sum*Sum*Sum;
     Sum4=Sum3*Sum;
   
     OneRM3Int = u*( (up+vp0)/S2 - vp0/S1 ) / (a3*b2);
     vpRM3Int  = -vp0*OneRM3Int + u*( 1.0/S1 - 1.0/S2 ) / a3;
     OneRM1Int = u*LogFac/a;
     vpRM1Int  = -vp0*OneRM1Int + u*(S2-S1)/a;
     OneR1Int  = u*a*0.5*(b2*LogFac + (up+vp0)*S2 - vp0*S1);
     vpR1Int   = -vp0*OneR1Int + u*a*(S2*S2*S2 - S1*S1*S1) / 3.0;
     OneR2Int  = u*a2*( (Sum3-vp03)/3.0 + up*b2 );
     vpR2Int   = -vp0*OneR2Int + u*a2*( (Sum4-vp04)/4.0 + up*(vp0 + 0.5*up)*b2 );
   }
  else // b is close to zero
   {
     double vp0pup = vp0 + up;
     double vp0pup2 = vp0pup*vp0pup;
     double vp0pup3 = vp0pup2*vp0pup;
     double vp0pup4 = vp0pup3*vp0pup;
     OneRM3Int = u*fabs( 1.0/vp02 - 1.0/vp0pup2 ) / (2.0*a3);
     vpRM3Int  = -vp0*OneRM3Int + u*fabs(1.0/vp0 - 1.0/vp0pup) / a3;
     OneRM1Int = u*fabs(log(vp0pup/vp0)) / a;
     vpRM1Int  = -vp0*OneRM1Int + u*up/a;
     OneR1Int  =  u*a*fabs(vp0pup2 - vp02) / 2.0;
     OneR2Int  = u*a2*fabs(vp0pup3 - vp03) / 3.0;
     vpR1Int   = -vp0*OneR1Int + OneR2Int/a;
     vpR2Int   = -vp0*OneR2Int + u*a2*fabs(vp0pup4 - vp04) / 4.0;

   };
  
  double CPCT[3], CPLT[3]; // 'cross product constant term' and 'cross product linear term'
  VecCross(X,Z,CPCT);
  VecCross(X,BP,CPLT);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nf=0;

  fval[nf++] = Y[0]*OneRM3Int - BP[0]*vpRM3Int;
  fval[nf++] = Y[1]*OneRM3Int - BP[1]*vpRM3Int;
  fval[nf++] = Y[2]*OneRM3Int - BP[2]*vpRM3Int;
  fval[nf++] = CPCT[0]*OneRM3Int + CPLT[0]*vpRM3Int;
  fval[nf++] = CPCT[1]*OneRM3Int + CPLT[1]*vpRM3Int;
  fval[nf++] = CPCT[2]*OneRM3Int + CPLT[2]*vpRM3Int;

  fval[nf++] = OneRM1Int;
  fval[nf++] = up*OneRM1Int;
  fval[nf++] = vpRM1Int;
  fval[nf++] = u*OneRM1Int;
  fval[nf++] = u*up*OneRM1Int;
  fval[nf++] = u*vpRM1Int;
  fval[nf++] = v*OneRM1Int;
  fval[nf++] = v*up*OneRM1Int;
  fval[nf++] = v*vpRM1Int;

  fval[nf++] = OneR1Int;
  fval[nf++] = up*OneR1Int;
  fval[nf++] = vpR1Int;
  fval[nf++] = u*OneR1Int;
  fval[nf++] = u*up*OneR1Int;
  fval[nf++] = u*vpR1Int;
  fval[nf++] = v*OneR1Int;
  fval[nf++] = v*up*OneR1Int;
  fval[nf++] = v*vpR1Int;

  fval[nf++] = OneR2Int;
  fval[nf++] = up*OneR2Int;
  fval[nf++] = vpR2Int;
  fval[nf++] = u*OneR2Int;
  fval[nf++] = u*up*OneR2Int;
  fval[nf++] = u*vpR2Int;
  fval[nf++] = v*OneR2Int;
  fval[nf++] = v*up*OneR2Int;
  fval[nf++] = v*vpR2Int;

} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void ComputeQIFIPPIData_Cubature(double **Va, double **Vb, QIFIPPIData *QIFD)
{
  /*--------------------------------------------------------------*/
  /* fill in the data structure used to pass parameters to the    */
  /* integrand routine                                            */
  /*--------------------------------------------------------------*/
  CFDData MyCFDData, *CFDD=&MyCFDData;

  CFDD->V0 = Va[0];
  VecSub(Va[1], Va[0], CFDD->A);
  VecSub(Va[2], Va[1], CFDD->B);

  CFDD->V0P = Vb[0];
  VecSub(Vb[1], Vb[0], CFDD->AP);
  VecSub(Vb[2], Vb[1], CFDD->BP);

  /*--------------------------------------------------------------*/
  /*- evaluate the adaptive cubature over the pair of triangles  -*/
  /*--------------------------------------------------------------*/
  const int fdim = 33;
  double F[fdim], E[fdim];

  CFDD->nCalls=0;
//#define CUBATURE4D
#ifdef CUBATURE4D
  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};
  adapt_integrate(fdim, CFDIntegrand4D, CFDD, 4, Lower, Upper,
                  MAXEVALS, ABSTOL, RELTOL, F, E);
#else
  double Lower[3]={0.0, 0.0, 0.0};
  double Upper[3]={1.0, 1.0, 1.0};
  CFDD->nCalls=0;
  adapt_integrate(fdim, CFDIntegrand3D, CFDD, 3, Lower, Upper,
                  MAXEVALS, ABSTOL, RELTOL, F, E);
#endif
  
  // error logging 
  if ( CFDD->nCalls >= MAXEVALS ) 
   { 
     int IsBad[fdim], NumBad=0, nn;
     memset(IsBad, 0, fdim*sizeof(int));
     for (nn=0; nn<fdim; nn++)
      { 
        if ( fabs(F[nn])>1.0e-6 && E[nn]>1.0e-3*fabs(F[nn]) )
         { IsBad[nn]=1;
           NumBad++;
         };

#ifndef _WIN32
        if (NumBad>0)
         { 
           char NewFileName[20];
           sprintf(NewFileName,"/tmp/FIPPI.XXXXXX");
           FILE *f=fdopen( mkstemp(NewFileName), "w");
           if (!f)
            { fprintf(stderr,"WARNING: panel integration irregularities detected\n");
            }
           else
            { fprintf(f,"%.12e %.12e %.12e \n", Va[0][0], Va[0][1], Va[0][2]);
              fprintf(f,"%.12e %.12e %.12e \n", Va[1][0], Va[1][1], Va[1][2]);
              fprintf(f,"%.12e %.12e %.12e \n", Va[2][0], Va[2][1], Va[2][2]);
              fprintf(f,"%.12e %.12e %.12e \n", Vb[0][0], Vb[0][1], Vb[0][2]);
              fprintf(f,"%.12e %.12e %.12e \n", Vb[1][0], Vb[1][1], Vb[1][2]);
              fprintf(f,"%.12e %.12e %.12e \n", Vb[2][0], Vb[2][1], Vb[2][2]);
              for(nn=0; nn<fdim; nn++)
               if (IsBad[nn])
                fprintf(f,"%2i: %+10.e %+10.e\n",nn,F[nn],E[nn]);
              fclose(f);
              fprintf(stderr,"WARNING: panel integration irregularities detected (%s)\n",NewFileName);
            };

         }; // if (NumBad>0) 
#endif // _WIN32

      }; // for(nn=0) ... 

   }; // if (nCalls > MAXEVALS)


  /*--------------------------------------------------------------*/
  /*- unpack the results into the output data record -------------*/
  /*--------------------------------------------------------------*/
  int nf=0;

  QIFD->xMxpRM3[0]   = F[nf++];
  QIFD->xMxpRM3[1]   = F[nf++];
  QIFD->xMxpRM3[2]   = F[nf++];
  QIFD->xXxpRM3[0]   = F[nf++];
  QIFD->xXxpRM3[1]   = F[nf++];
  QIFD->xXxpRM3[2]   = F[nf++];

  QIFD->uvupvpRM1[0] = F[nf++];
  QIFD->uvupvpRM1[1] = F[nf++];
  QIFD->uvupvpRM1[2] = F[nf++];
  QIFD->uvupvpRM1[3] = F[nf++];
  QIFD->uvupvpRM1[4] = F[nf++];
  QIFD->uvupvpRM1[5] = F[nf++];
  QIFD->uvupvpRM1[6] = F[nf++];
  QIFD->uvupvpRM1[7] = F[nf++];
  QIFD->uvupvpRM1[8] = F[nf++];

  QIFD->uvupvpR1[0] = F[nf++];
  QIFD->uvupvpR1[1] = F[nf++];
  QIFD->uvupvpR1[2] = F[nf++];
  QIFD->uvupvpR1[3] = F[nf++];
  QIFD->uvupvpR1[4] = F[nf++];
  QIFD->uvupvpR1[5] = F[nf++];
  QIFD->uvupvpR1[6] = F[nf++];
  QIFD->uvupvpR1[7] = F[nf++];
  QIFD->uvupvpR1[8] = F[nf++];

  QIFD->uvupvpR2[0] = F[nf++];
  QIFD->uvupvpR2[1] = F[nf++];
  QIFD->uvupvpR2[2] = F[nf++];
  QIFD->uvupvpR2[3] = F[nf++];
  QIFD->uvupvpR2[4] = F[nf++];
  QIFD->uvupvpR2[5] = F[nf++];
  QIFD->uvupvpR2[6] = F[nf++];
  QIFD->uvupvpR2[7] = F[nf++];
  QIFD->uvupvpR2[8] = F[nf++];

}

/*--------------------------------------------------------------*/
/*- routine for computing Q-independent FIPPIs. this is a       */
/*- switchboard routine that calls one of two actual            */
/*- computational routines.                                     */
/*-                                                             */
/*- inputs:                                                     */
/*-                                                             */
/*-  Va[i][j] = jth cartesian coord of ith vertex of panel A    */
/*-  Vb       = similarly for panel B                           */
/*-  FD       = must point to a preallocated QIFIPPIData        */
/*--------------------------------------------------------------*/
void ComputeQIFIPPIData(double **Va, double **Vb, int ncv, QIFIPPIData *QIFD)
{ 
  /*--------------------------------------------------------------*/
  /*- if there are no common vertices, then use 4-dimensional    -*/
  /*- adaptive cubature over both triangles to compute the FIPPIs-*/
  /*--------------------------------------------------------------*/
  if (ncv==0)
   ComputeQIFIPPIData_Cubature(Va, Vb, QIFD);
  else
   { if (RWGGeometry::UseTaylorDuffyV2P0)
      ComputeQIFIPPIData_TaylorDuffyV2P0(Va[0], Va[1], Va[2], Vb[1], Vb[2], QIFD);
     else
      ComputeQIFIPPIData_TaylorDuffy(Va[0], Va[1], Va[2], Vb[1], Vb[2], QIFD);
   };

#define FIPPIFIX_20150306
#ifdef FIPPIFIX_20150306
  QIFD->x0[0]=(Va[0][0]+Va[1][0]+Va[2][0])/3.0;
  QIFD->x0[1]=(Va[0][1]+Va[1][1]+Va[2][1])/3.0;
  QIFD->x0[2]=(Va[0][2]+Va[1][2]+Va[2][2])/3.0;
#endif

#define FIPPIFIX_20200606
#ifdef FIPPIFIX_20200606
  /* Check NAN for RM1 and RM3 in QIFD, Override to 0 if there is NAN*/
  int j=0;
  for (j=0; j<3; j++){
    if (std::isnan(QIFD->xMxpRM3[j]))
      QIFD->xMxpRM3[j] = 0;
    if (std::isnan(QIFD->xXxpRM3[j]))
      QIFD->xXxpRM3[j] = 0;
  }
  for (j=0; j<9; j++){
    if (std::isnan(QIFD->uvupvpRM1[j]))
      QIFD->uvupvpRM1[j] = 0;
  }
#endif

}

} // namespace scuff
