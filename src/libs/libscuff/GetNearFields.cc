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
 * GetNearFields.cc -- libscuff routines for computing the E and H
 *                  -- fields due to RWG basis functions at evaluation 
 *                  -- points lying on or near the source panels
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libscuffInternals.h>
#include <libTriInt.h>

#define II cdouble(0.0,1.0)

namespace scuff{

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
#define PINDEX_M5 0
#define PINDEX_M3 1
#define PINDEX_M1 2
#define PINDEX_P1 3
#define NUMPS 4
static double PValues[NUMPS]={-5.0,-3.0,-1.0,1.0};

#define NINDEX_M1 0
#define NINDEX_P1 1
#define NUMNS 2

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
cdouble ExpRel(cdouble x, int n);

/*--------------------------------------------------------------*/
/* This routine computes the I^p, J^p, and dJ^p / d\ell         */
/* integrals for a single panel edge.                           */
/*                                                              */
/* I[0,1,2,3] = I_{-5}, I_{-3}, I_{-1}, I_{1}                   */
/* J[0,1,2,3] = J_{-5}, J_{-3}, J_{-1}, J_{1}                   */
/*--------------------------------------------------------------*/
void GetEdgeIntegrals(const double sM, const double sP, 
                      const double t, const double w,
                      double I[NUMPS], double J[NUMPS], double dJdl[NUMPS])
{
#if 0
20140915 Why is this not necessary?
  if ( (sP>0.0 && sM<0.0) || (sP<0.0 && sM>0.0) )
   { double I2[NUMPS], J2[NUMPS], dJdl2[NUMPS];
     GetEdgeIntegrals(0.0, sP, t, w, I, J, dJdl);
     GetEdgeIntegrals(0.0, sM, t, w, I2, J2, dJdl2);
     for(int np=0; np<NUMPS; np++)
      { I[np]    += I2[np];
        J[np]    += J2[np];
        dJdl[np] += dJdl2[np];
      };
     return; 
   };
#endif

  double LengthScale = fabs(sP - sM);
  double sP2 = sP*sP, sP3=sP2*sP;
  double sM2 = sM*sM, sM3=sM2*sM;
  double t2  = t*t;
  double w2  = w*w;
  double w3  = w2*w;
  double X2  = t2 + w2;
  double X4  = X2*X2;
  double X   = sqrt(X2);
  double RP  = sqrt(sP2 + X2);
  double RM  = sqrt(sM2 + X2);

  double LogFac;
  if ( X > 1.0e-4*LengthScale )
   LogFac = log( (RP+sP) / (RM+sM) );
  else if ( (fabs(sM)>1.0e-4*LengthScale) && (fabs(sP)>1.0e-4*LengthScale) )
   LogFac = log( fabs(sP) / fabs(sM) );
  else
   LogFac = 0.0;

  double atanFac;
  if ( fabs(w)<1.0e-8*LengthScale  || fabs(t)<1.0e-8*LengthScale )
   atanFac = 0.0;
  else 
   atanFac = atan( (sP*w)/(RP*t) ) - atan( (sM*w)/(RM*t) );

  if ( fabs(t) < 1.0e-8*LengthScale )
   { 
     memset(I, 0, NUMPS*sizeof(double));
   }
  else
   { if ( fabs(w)<1.0e-8*LengthScale )
      { I[0]=  sP*(2.0*sP*sP + 3.0*t*t)/(3.0*t*t*t*RP*RP*RP)
              -sM*(2.0*sM*sM + 3.0*t*t)/(3.0*t*t*t*RM*RM*RM);
        I[1]= (sP/RP - sM/RM) / t;
      }
     else
      { I[0]= t*(sM/RM - sP/RP)/(w2*X2) + atanFac/w3;
        I[1]= atanFac / w;
      };
     I[2]= w*atanFac + t*LogFac;
     I[3]= t*(RP*sP - RM*sM + LogFac*(t2 + 3.0*w2))/2.0 
             +w3*atanFac;
   };

  if ( X < 1.0e-4*LengthScale )
   J[0] = 0.0;
  else
   J[0]= (sP/RP - sM/RM) / X2;
  J[1]= LogFac;
  J[2]= (RP*sP - RM*sM + LogFac*X2)/2.0;
  J[3]= (2.0*RP*sP3 - 2.0*RM*sM3 + 5.0*X2*(RP*sP - RM*sM)
        + 3.0*LogFac*X4)/8.0;

  memset(dJdl,0,NUMPS*sizeof(double));
  if ( fabs(RM) > 1.0e-8*LengthScale )
   { for(int np=1; np<NUMPS; np++)
      dJdl[np] += pow(RM, PValues[np]+2.0);
   };
  if ( fabs(RP) > 1.0e-8*LengthScale )
   { for(int np=1; np<NUMPS; np++)
      dJdl[np] -= pow(RP, PValues[np]+2.0);
   };

}

/*--------------------------------------------------------------*/
/* This routine computes the \mathcal I^p, \mathcal J^p         */
/* integrals for a single panel.                                */
/*                                                              */
/* The panel lives in the (u,v) plane with vertices u_i, v_i    */
/* (i=0, 1, 2). The coordinates of the observation point in the */
/* (u,v,w) coordinate system are (0,0,w0).                      */
/*                                                              */
/* LengthScale is used as a reference to determine when         */
/* quantities with length dimensions are to be considered       */
/* 'small.'                                                     */
/*                                                              */
/* Note that the components of vector-valued quantities are     */
/* returned *in the panel coordinate system* (in which the      */
/* panel normal defines the Z axis and one panel edge defines   */
/* the X axis).                                                 */
/*                                                              */
/* For labeling, I use a Penrose-style tensor notation in       */
/* which vector-valued quantities (whose components are indexed */
/* by an index Mu=0,1,2) have 'Mu' built into their symbol name,*/
/* while tensor-valued quantities (whose components are indexed */
/* by a pair of indices Mu, Nu=0,1,2) have 'MuNu' built into the*/
/* symbol name, etc.                                            */
/*                                                              */
/* ScriptI[p]              = I^p                                */
/*                         = \int_T r^p dA / 4\pi               */
/*                                                              */
/* ScriptJMu[p][Mu]        = J^p_\Mu                            */
/*                         = \int_T \rho_\Mu r^p dA / 4\pi      */
/*                                                              */
/* (here \int_T ... dA denotes an integral over the area of the */
/* triangle, r is the distance from the integration point to    */
/* the observation point, and \rho_\Mu (\Mu=0,1) are the        */
/* in-plane components of the r vector.                         */
/*                                                              */
/* dScriptIMu[p][Mu]       = d_\Mu I^p                          */
/* dScriptJMuNu[p][Mu][Nu] = d_\Mu J^p_\Nu                      */
/*--------------------------------------------------------------*/
void GetScriptIJ(const double u[3], const double v[3], const double w0, const double LengthScale,
                 double ScriptI[NUMPS], double dScriptIMu[NUMPS][3],
                 double ScriptJMu[NUMPS][3], double dScriptJMuNu[NUMPS][3][3])
{ 
  for(int np=0; np<NUMPS; np++)
   { ScriptI[np]=0.0;
     for(int Mu=0; Mu<3; Mu++)
      { dScriptIMu[np][Mu]=ScriptJMu[np][Mu]=0.0;
        for(int Nu=0; Nu<3; Nu++)
         dScriptJMuNu[np][Mu][Nu]=0.0;
      };
   };

  int NumEdges=3; // this could be changed to handle higher polygons
  double Beta=0.0; // angle for interior point correction
  bool OnEdgeOrVertex=false;
  for(int i=0; i<NumEdges; i++)
   {
     double ui=u[i], uip1=u[ (i+1) % 3 ];
     double vi=v[i], vip1=v[ (i+1) % 3 ];

     double lHat[2], mHat[2], Length;
     lHat[0] = (uip1-ui);
     lHat[1] = (vip1-vi);
     Length = sqrt(lHat[0]*lHat[0] + lHat[1]*lHat[1]);
     lHat[0] /= Length;
     lHat[1] /= Length;
     mHat[0] = lHat[1];
     mHat[1] = -lHat[0];
    
     double t  = -(ui*mHat[0] + vi*mHat[1]);
     if ( fabs(t) < 1.0e-8*LengthScale )
      OnEdgeOrVertex=true;
     else
      Beta +=   acos( (ui*uip1 + vi*vip1)
                  / ( sqrt(ui*ui+vi*vi)*sqrt(uip1*uip1 + vip1*vip1) ) );

     double sM = ui*lHat[0] + vi*lHat[1];
     double sP = uip1*lHat[0] + vip1*lHat[1];

     double I[NUMPS], J[NUMPS], dJdl[NUMPS];
     GetEdgeIntegrals(sM, sP, t, w0, I, J, dJdl);

     for(int np=0; np<NUMPS; np++)
      { 
        double P = PValues[np];
        double Denom = (P+2.0);

        ScriptI[np]      -= I[np] / Denom;

        ScriptJMu[np][0] += mHat[0] * J[np] / Denom;
        ScriptJMu[np][1] += mHat[1] * J[np] / Denom;

        if (np>0)
         { 
           double dJdm = (P+2.0)*t*J[np-1];
           for(int Mu=0; Mu<2; Mu++)
            for(int Nu=0; Nu<2; Nu++)
             dScriptJMuNu[np][Mu][Nu] 
              += ( dJdl[np]*lHat[Mu] + dJdm*mHat[Mu] ) * mHat[Nu] / Denom;
         };

      };

   };

  /***************************************************************/
  /* handle derivatives defined in terms of other I, J functions */
  /***************************************************************/
  for(int np=1; np<NUMPS; np++)
   {
     double P = PValues[np];

     dScriptIMu[np][0] = -P*ScriptJMu[np-1][0];
     dScriptIMu[np][1] = -P*ScriptJMu[np-1][1];
     dScriptIMu[np][2] =  P*w0*ScriptI[np-1];

     dScriptJMuNu[np][2][0] = P*w0*ScriptJMu[np-1][0];
     dScriptJMuNu[np][2][1] = P*w0*ScriptJMu[np-1][1];
     dScriptJMuNu[np][2][2] = P*w0*ScriptJMu[np-1][2];

   };

  /***************************************************************/
  /* interior point correction ***********************************/
  /***************************************************************/
  if ( EqualFloat(Beta, M_PI) || EqualFloat(Beta, 2.0*M_PI) )
   { 
     double Sign  = w0 > 0.0 ? 1.0 : -1.0;
     double w0Abs = fabs(w0);

     if (w0Abs > 1.0e-8*LengthScale)
      {
        double w02  = w0*w0, w0Powers[NUMPS];
        w0Powers[0] = pow( fabs(w0), PValues[0]+2.0 );
        for(int np=1; np<NUMPS; np++)
         w0Powers[np] = w0Powers[np-1]*w02;

        for(int np=0; np<NUMPS; np++)
         { double PP2=PValues[np] + 2.0;
           ScriptI[np] -= Beta * w0Powers[np] / PP2;
           dScriptIMu[np][2] -= Sign * Beta * w0Powers[np] / w0Abs;
         };
      };
   };

}

/***************************************************************/
/* Compute the first few frequency-independent terms in the    */
/* small-k expansion of the reduced scalar and vector          */
/* potentials at X0 due to RWG currents on a single triangular */
/* panel. These are the quantities named p^{(n)} and a^{(n)}_i */
/* in the writeup. We compute them for n=-1, 1. We also        */
/* compute their derivatives with respect to displacements of  */
/* the evaluation point.                                       */
/***************************************************************/
void GetStaticPanelIntegrals(RWGSurface *S, const int nPanel, const int iQ,
                             const double X0[3],
                             double pn[NUMNS], double an[NUMNS][3],
                             double dpn[NUMNS][3], double dan[NUMNS][3][3],
                             double ddpn[NUMNS][3][3], double dcurlan[NUMNS][3][3])
{
  /***************************************************************/
  /* unpack panel vertices ***************************************/
  /***************************************************************/
  RWGPanel *Pan = S->Panels[nPanel];
  double *V[3], *Q;
  V[0] = S->Vertices + 3*(Pan->VI[0]);
  V[1] = S->Vertices + 3*(Pan->VI[1]);
  V[2] = S->Vertices + 3*(Pan->VI[2]);
  Q    = V[iQ];
  double *XC = Pan->Centroid;

  /***************************************************************/
  /* compute the uHat, vHat, wHat vectors that define the        */
  /* modified coordinate system (in which the panel lives in     */
  /* the uv plane with one edge in the u direction).             */
  /* note we can't use P->nHat as our wHat vector because it     */
  /* doesn't necessarily satisfy a right-hand rule with respect  */
  /* to the vertex ordering.                                     */
  /***************************************************************/
  double uHat[3], vHat[3], wHat[3];
  VecSub(V[1], V[0], uHat);
  VecNormalize(uHat);
  VecSub(V[2], V[0], vHat); // this is just temporary
  VecCross(uHat, vHat, wHat);
  VecNormalize(wHat);
  VecCross(wHat, uHat, vHat);

  /***************************************************************/
  /* compute coordinates of evaluation point in modified system, */
  /* x0 = (Rho0, w0)                                             */
  /***************************************************************/
  double Rho0[3];
  VecSub(X0, XC, Rho0);
  double w0 = VecDot(Rho0, wHat);
  VecScaleAdd(X0, -w0, wHat, Rho0);

  /***************************************************************/
  /* QBar is the vector Q-Rho0 where Q is the RWG source/sink    */
  /* vertex. This vector only has components in the plane of the */
  /* triangle. QBar[0,1,2] are its components in the default     */
  /* coordinate system, while QBarMu[0,1,2] are its components in*/
  /* the (u,v,w) system.                                         */
  /***************************************************************/
  double QBar[3], QBarMu[3];
  VecScaleAdd(Q, -1.0, Rho0, QBar);
  QBarMu[0] = VecDot(QBar, uHat);
  QBarMu[1] = VecDot(QBar, vHat);
  QBarMu[2] = 0.0;

  /***************************************************************/
  /* compute coordinates of panel vertices in modified system,   */
  /* V_i = (u_i, v_i, 0)                                         */
  /***************************************************************/
  double u[3], v[3];
  for(int i=0; i<3; i++)
   { double VmRho[3];
     VecSub(V[i], Rho0, VmRho);
     u[i] = VecDot(VmRho, uHat);
     v[i] = VecDot(VmRho, vHat);
   };

  /***************************************************************/
  /* get the scriptI, scriptJ integrals for this panel ***********/
  /***************************************************************/
  double ScriptI[NUMPS];
  double dScriptIMu[NUMPS][3];
  double ScriptJMu[NUMPS][3];
  double dScriptJMuNu[NUMPS][3][3];

  GetScriptIJ(u, v, w0, Pan->Radius,
              ScriptI, dScriptIMu, ScriptJMu, dScriptJMuNu);

  /***************************************************************/
  /* assemble p^p and a^p_\mu functions and their derivatives    */
  /* with components of vectors still referred to the panel      */
  /* coordinate system (hence the Mu, MuNu suffixes). This is    */
  /* equation (C...) in the memo.                                */
  /***************************************************************/
  double anMu[NUMNS][3], dpnMu[NUMNS][3], danMuNu[NUMNS][3][3];
  for(int NIndex=NINDEX_M1; NIndex<=NINDEX_P1; NIndex++)
   { 
     int PIndex = (NIndex==NINDEX_M1) ? PINDEX_M1 : PINDEX_P1;
     double P = PValues[PIndex];

     pn[NIndex] = ScriptI[PIndex];

     for(int Mu=0; Mu<3; Mu++)
      { 
        dpnMu[NIndex][Mu] = dScriptIMu[PIndex][Mu];

        anMu[NIndex][Mu] 
         = ScriptJMu[PIndex][Mu] - QBarMu[Mu]*ScriptI[PIndex];

        danMuNu[NIndex][Mu][2]=0.0;

        danMuNu[NIndex][2][Mu]
         = P*w0*(  ScriptJMu[PIndex-1][Mu]
                  -QBarMu[Mu]*ScriptI[PIndex-1]
                );
      };

     for(int Mu=0; Mu<2; Mu++)
      for(int Nu=0; Nu<2; Nu++)
       danMuNu[NIndex][Mu][Nu]
        = dScriptJMuNu[PIndex][Mu][Nu]
          +( (Mu==Nu) ? ScriptI[PIndex] : 0.0 )
          +P*QBarMu[Nu]*ScriptJMu[PIndex-1][Mu];

   };

  /***************************************************************/
  /* second derivatives ******************************************/
  /***************************************************************/
  double ddpnMuNu[2][3][3], ddanMuNuRho[2][3][3][3];
  for(int NIndex=NINDEX_M1; NIndex<=NINDEX_P1; NIndex++)
   { 
     int PIndex = (NIndex==NINDEX_M1) ? PINDEX_M1 : PINDEX_P1;
     double P=PValues[PIndex];

     for(int Mu=0; Mu<2; Mu++)
      for(int Nu=0; Nu<2; Nu++)
       ddpnMuNu[NIndex][Mu][Nu] 
        = -P*dScriptJMuNu[PIndex-1][Mu][Nu];

     for(int Mu=0; Mu<3; Mu++)
      ddpnMuNu[NIndex][2][Mu]
       = ddpnMuNu[NIndex][Mu][2]
       = -P*(P-2.0)*w0*ScriptJMu[PIndex-2][Mu];

     ddpnMuNu[NIndex][2][2]
      = P*ScriptI[PIndex-1] + P*(P-2.0)*w0*w0*ScriptI[PIndex-2];

     // equation C...
     for(int A=0; A<2; A++)
      {
        ddanMuNuRho[NIndex][2][2][A]
         =  P*( ScriptJMu[PIndex-1][A] - QBarMu[A]*ScriptI[PIndex-1])
           +P*(P-2.0)*w0*w0*( ScriptJMu[PIndex-2][A] - QBarMu[A]*ScriptI[PIndex-2]);

        for(int B=0; B<3; B++)
          ddanMuNuRho[NIndex][2][A][B]
           = ddanMuNuRho[NIndex][A][2][B]
           = P*w0*(  dScriptJMuNu[PIndex-1][A][B]
                    + ( (A==B) ? ScriptI[PIndex-1] : 0.0)
                    + (P-2.0)*QBarMu[B]*ScriptJMu[PIndex-2][A]
                  );
      };
     ddanMuNuRho[NIndex][2][2][2]=0.0;
      
   };

  double dcurlanMuNu[NUMNS][3][3];
  for(int NIndex=NINDEX_M1; NIndex<=NINDEX_P1; NIndex++)
   { 
     int PIndex = (NIndex==NINDEX_M1) ? PINDEX_M1 : PINDEX_P1;
     double P=PValues[PIndex];

     dcurlanMuNu[NIndex][0][0] = -ddanMuNuRho[NIndex][0][2][1];
     dcurlanMuNu[NIndex][0][1] = +ddanMuNuRho[NIndex][0][2][0];
     dcurlanMuNu[NIndex][0][2] = P*( ScriptJMu[PIndex-1][1]
                                     +QBarMu[1]*dScriptJMuNu[PIndex-1][0][0]
                                     -QBarMu[0]*dScriptJMuNu[PIndex-1][0][1]
                                   );

     dcurlanMuNu[NIndex][1][0] = -ddanMuNuRho[NIndex][1][2][1];
     dcurlanMuNu[NIndex][1][1] = +ddanMuNuRho[NIndex][1][2][0];
     dcurlanMuNu[NIndex][1][2] = P*( -ScriptJMu[PIndex-1][0]
                                     +QBarMu[1]*dScriptJMuNu[PIndex-1][1][0]
                                     -QBarMu[0]*dScriptJMuNu[PIndex-1][1][1]
                                   );

 
     dcurlanMuNu[NIndex][2][0] = -ddanMuNuRho[NIndex][2][2][1];
     dcurlanMuNu[NIndex][2][1] = +ddanMuNuRho[NIndex][2][2][0];
     dcurlanMuNu[NIndex][2][2] =  ddanMuNuRho[NIndex][2][0][1] 
                                  - ddanMuNuRho[NIndex][2][1][0];

   };

  /***************************************************************/
  /* finally, convert vector and tensor components back to       */
  /* original component system.                                  */
  /* note: the transformation matrix M is defined such that      */
  /* vector components transform like                            */
  /*  v_i = \sum_\mu M_{\mu i} v_\mu                             */
  /* which means derivatives should transform like               */
  /*  d_i f = \sum_\mu M_{i\mu} d_mu f                           */
  /* i.e. with the inverse (transpose) of M, but that doesn't    */
  /* seem to be the case! what am i doing wrong?                 */
  /***************************************************************/
  double M[3][3]; 
  memcpy(M[0], uHat, 3*sizeof(double));
  memcpy(M[1], vHat, 3*sizeof(double));
  memcpy(M[2], wHat, 3*sizeof(double));

  for(int nn=0; nn<NUMNS; nn++)
   { 
     // one-index quantities
     for(int i=0; i<3; i++)
      { an[nn][i] = dpn[nn][i] = 0.0;
        for(int Mu=0; Mu<3; Mu++)
         { an[nn][i] += M[Mu][i] * anMu[nn][Mu];
           dpn[nn][i] += M[Mu][i] * dpnMu[nn][Mu];
         };
      };

     // two-index quantities
     for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
       { dan[nn][i][j] = 0.0;
         ddpn[nn][i][j] = 0.0;
         dcurlan[nn][i][j] = 0.0;
         for(int Mu=0; Mu<3; Mu++)
          for(int Nu=0; Nu<3; Nu++)
           { dan[nn][i][j] += M[Mu][i]*M[Nu][j]*danMuNu[nn][Mu][Nu];
             ddpn[nn][i][j] += M[Mu][i]*M[Nu][j]*ddpnMuNu[nn][Mu][Nu];
             dcurlan[nn][i][j] += M[Mu][i]*M[Nu][j]*dcurlanMuNu[nn][Mu][Nu];
           };
       };
    };

}

/*--------------------------------------------------------------*/
/*- get the desingularized scalar green's function and its first*/
/*- and second derivatives. (The 'desingularized scalar green's */
/*- function) is the usual Helmholtz kernel minus the first     */
/*- three terms of its small-k expansion.)                      */
/*--------------------------------------------------------------*/
#define EXPRELTOL  1.0e-8
#define EXPRELTOL2 EXPRELTOL*EXPRELTOL
cdouble GetGDS(double R[3], cdouble k, cdouble dG[3],
               cdouble ddG[3][3], bool Need_ddG=true)
{
  double r = VecNorm(R);
  if (r==0.0)
   { memset(dG, 0, 3*sizeof(cdouble));
     if (Need_ddG)
      { memset(ddG[0], 0, 3*sizeof(cdouble));
        memset(ddG[1], 0, 3*sizeof(cdouble));
        memset(ddG[2], 0, 3*sizeof(cdouble));
      };
     return 0.0; 
   };

  cdouble IK = II*k;

  /*--------------------------------------------------------------*/
  /* computation of ExpRelBar(ikr,1...3)                          */
  /*--------------------------------------------------------------*/
  cdouble x = IK*r;
  cdouble Term1 = x;
  cdouble Term2 = x*Term1/2.0;
  cdouble Term  = x*Term2/3.0;
  cdouble ERB3;
  if ( abs(x) > 0.1 )
   { 
     ERB3 = exp(x) - 1.0 - Term1 - Term2;
   }
  else // use series expansion
   { 
     cdouble Sum  = Term;
     for(int m=4; m<100; m++)
      { Term*=x/((double)m);
        Sum+=Term;
        if ( norm(Term) < EXPRELTOL2*norm(Sum) )
         break;
      };
     ERB3 = Sum;
   };

  cdouble ERB2 = ERB3 + Term2;
  cdouble ERB1 = ERB2 + Term1;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Denom=4.0*M_PI*r, Denom2=Denom*r, Denom3=Denom2*r;

  cdouble Factor1 = IK*ERB2/Denom2 - ERB3/Denom3;

  dG[0]     = R[0]*Factor1;
  dG[1]     = R[1]*Factor1;
  dG[2]     = R[2]*Factor1;

  if (Need_ddG)
   { 
     double Denom4=Denom3*r, Denom5=Denom4*r;
     cdouble Factor2 = IK*IK*ERB1/Denom3 - 3.0*IK*ERB2/Denom4 + 3.0*ERB3/Denom5;
     ddG[0][0] = Factor1 + R[0]*R[0]*Factor2;
     ddG[1][1] = Factor1 + R[1]*R[1]*Factor2;
     ddG[2][2] = Factor1 + R[2]*R[2]*Factor2;
     ddG[0][1] = ddG[1][0] = R[0]*R[1]*Factor2;
     ddG[0][2] = ddG[2][0] = R[0]*R[2]*Factor2;
     ddG[1][2] = ddG[2][1] = R[1]*R[2]*Factor2;
   };

  return ERB3 / Denom;
} 


/*--------------------------------------------------------------*/
/* get desingularized contributions to reduced potentials       */
/*--------------------------------------------------------------*/
void GetDSReducedPotentials(RWGSurface *S, const int np, const int iQ,
                            const double X0[3], const cdouble k,
                            cdouble *pDS, cdouble aDS[3],
                            cdouble dpDS[3], cdouble daDS[3][3],
                            cdouble ddpDS[3][3], cdouble dcurlaDS[3][3])
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGPanel *P=S->Panels[np];
  double *Q  = S->Vertices + 3*P->VI[iQ];
  double *V1 = S->Vertices + 3*P->VI[(iQ+1)%3];
  double *V2 = S->Vertices + 3*P->VI[(iQ+2)%3];
  double PreFactor = 2.0*P->Area;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double A[3], B[3];
  VecSub(V1, Q, A);
  VecSub(V2, V1, B);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  *pDS=0;
  for(int Mu=0; Mu<3; Mu++)
   { aDS[Mu]=dpDS[Mu]=0.0;
     for(int Nu=0; Nu<3; Nu++)
      daDS[Mu][Nu]=ddpDS[Mu][Nu]=dcurlaDS[Mu][Nu]=0.0;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumPts;
  int Order=9;
  double *TCR=GetTCR(Order, &NumPts);
  for(int nPt=0, ncp=0; nPt<NumPts; nPt++)
   { 
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     double u=TCR[ncp++];
     double v=TCR[ncp++];
     double w=TCR[ncp++] * PreFactor;
     u+=v;

     double X[3], F[3], R[3];
     for(int Mu=0; Mu<3; Mu++)
      { 
        F[Mu] = u*A[Mu] + v*B[Mu];
        X[Mu] = F[Mu] + Q[Mu];
        R[Mu] = X0[Mu] - X[Mu];
      };

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     cdouble dG[3], ddG[3][3];
     cdouble GG=GetGDS(R, k, dG, ddG);
     *pDS += w*GG;
     for(int Mu=0; Mu<3; Mu++)
      { dpDS[Mu] += w*dG[Mu];
        aDS[Mu]  += w*GG*F[Mu];
        for(int Nu=0; Nu<3; Nu++)
         { daDS[Mu][Nu] += w*dG[Mu]*F[Nu];
           ddpDS[Mu][Nu]+= w*ddG[Mu][Nu];
           int NP1=(Nu+1)%3, NP2=(Nu+2)%3;
           dcurlaDS[Mu][Nu] 
            += w*(ddG[Mu][NP1]*F[NP2] - ddG[Mu][NP2]*F[NP1]);
         };
      };
     
   };

}

/*--------------------------------------------------------------*/
/*- fvrp = 'four-vector reduced potential' with components      */
/*-  fvrp[0]    = p       (reduced scalar potential)            */
/*-  fvrp[1..3] = a[0..2] (reduced vector potential)            */
/*- dfvrp = derivative of fvp with components                   */
/*-  dfvrp[Mu][Nu] = d_\mu fvrp[Nu]                             */
/*--------------------------------------------------------------*/
void GetReducedPotentials_Nearby(RWGSurface *S, const int np, const int iQ,
                                 const double X0[3], const cdouble k,
                                 cdouble *p, cdouble a[3],
                                 cdouble dp[3], cdouble da[3][3],
                                 cdouble ddp[3][3], cdouble dcurla[3][3])
{
   /*--------------------------------------------------------------*/
   /*- get singular contributions ---------------------------------*/
   /*--------------------------------------------------------------*/
   double pn[NUMNS], an[NUMNS][3];
   double dpn[NUMNS][3], dan[NUMNS][3][3];
   double ddpn[NUMNS][3][3], dcurlan[NUMNS][3][3];
   GetStaticPanelIntegrals(S, np, iQ, X0,
                           pn, an, dpn, dan, ddpn, dcurlan);
   double pM1=pn[0],             pP1=pn[1];
   double *aM1=an[0],            *aP1=an[1];
   double *dpM1=dpn[0],          *dpP1=dpn[1];

   double *daM1[3], *daP1[3], *ddpM1[3], *ddpP1[3], *dcurlaM1[3], *dcurlaP1[3];
   for(int Mu=0; Mu<3; Mu++)
    { daM1[Mu]=dan[0][Mu];           daP1[Mu]=dan[1][Mu];
      ddpM1[Mu]=ddpn[0][Mu];         ddpP1[Mu]=ddpn[1][Mu];
      dcurlaM1[Mu]=dcurlan[0][Mu];   dcurlaP1[Mu]=dcurlan[1][Mu];
    };
   
   /*--------------------------------------------------------------*/
   /* get a^{(0)} and p^{(0)} contributions                        */
   /*--------------------------------------------------------------*/
   RWGPanel *P   = S->Panels[np];
   double *Q     = S->Vertices + 3*P->VI[iQ];
   double *XC    = P->Centroid;
   double p0, a0[3];
   p0 = P->Area;
   a0[0] = P->Area*(XC[0]-Q[0]);
   a0[1] = P->Area*(XC[1]-Q[1]);
   a0[2] = P->Area*(XC[2]-Q[2]);

   /*--------------------------------------------------------------*/
   /* get desingularized contributions                             */
   /*--------------------------------------------------------------*/
   cdouble pDS, aDS[3], dpDS[3], daDS[3][3], ddpDS[3][3], dcurlaDS[3][3];
   GetDSReducedPotentials(S, np, iQ, X0, k,
                          &pDS, aDS, dpDS, daDS, ddpDS, dcurlaDS);
   
   /*--------------------------------------------------------------*/ 
   /* assemble the results                                         */
   /*--------------------------------------------------------------*/
   cdouble IK  = II*k;
   cdouble CM1 = 1.0   / (4.0*M_PI);
   cdouble C0  = IK    / (4.0*M_PI);
   cdouble CP1 = IK*IK / (8.0*M_PI);
   *p = CM1*pM1 + C0*p0 + CP1*pP1 + pDS;
   for(int Mu=0; Mu<3; Mu++)
    { 
      a[Mu] = CM1*aM1[Mu] + C0*a0[Mu] + CP1*aP1[Mu] + aDS[Mu];

      dp[Mu] = CM1*dpM1[Mu] + CP1*dpP1[Mu] + dpDS[Mu];

      for(int Nu=0; Nu<3; Nu++)
       { 
         da[Mu][Nu] = CM1*daM1[Mu][Nu] + CP1*daP1[Mu][Nu] + daDS[Mu][Nu];
         ddp[Mu][Nu] = CM1*ddpM1[Mu][Nu] + CP1*ddpP1[Mu][Nu] + ddpDS[Mu][Nu];
         dcurla[Mu][Nu] = CM1*dcurlaM1[Mu][Nu] + CP1*dcurlaP1[Mu][Nu] + dcurlaDS[Mu][Nu];
       };
    };

   /*--------------------------------------------------------------*/
   /*- put in the RWG basis function prefactors -------------------*/
   /*--------------------------------------------------------------*/
   double aPrefac = S->Edges[P->EI[iQ]]->Length / (2.0*P->Area);
   double pPrefac = 2.0*aPrefac;
   (*p)*=pPrefac;
   for(int Mu=0; Mu<3; Mu++)
    { dp[Mu] *= pPrefac;
      a[Mu]  *= aPrefac;
      for(int Nu=0; Nu<3; Nu++)
       { da[Mu][Nu]     *= aPrefac;
         ddp[Mu][Nu]    *= pPrefac;
         dcurla[Mu][Nu] *= aPrefac;
       };
    };
}

/***************************************************************/
/* get the reduced potentials of a full RWG basis function by  */
/* summing contributions from the positive and negative panels */
/***************************************************************/
void GetReducedPotentials_Nearby(RWGSurface *S, const int ne,
                                 const double X0[3],  const cdouble k,
                                 cdouble *p, cdouble a[3],
                                 cdouble dp[3], cdouble da[3][3],
                                 cdouble ddp[3][3], cdouble dcurla[3][3])
{
  RWGEdge *E = S->Edges[ne];

  GetReducedPotentials_Nearby(S, E->iPPanel, E->PIndex, X0, k,
                              p, a, dp, da, ddp, dcurla);

  cdouble pM, aM[3], dpM[3], daM[3][3], ddpM[3][3], dcurlaM[3][3];
  GetReducedPotentials_Nearby(S, E->iMPanel, E->MIndex, X0, k,
                              &pM, aM, dpM, daM, ddpM, dcurlaM);

  *p -= pM;
  for(int Mu=0; Mu<3; Mu++)
   { a[Mu]  -= aM[Mu];
     dp[Mu] -= dpM[Mu];
     for(int Nu=0; Nu<3; Nu++)
      { da[Mu][Nu]     -= daM[Mu][Nu];
        ddp[Mu][Nu]    -= ddpM[Mu][Nu];
        dcurla[Mu][Nu] -= dcurlaM[Mu][Nu];
      };
   };

}

/***************************************************************/
/* E = ikZ * kAlpha*e - nAlpha*h                               */
/* H =      +kAlpha*h + (ik/Z)*nAlpha*e                        */
/***************************************************************/
void GetReducedFields_Nearby(RWGSurface *S, const int ne,
                             const double X0[3], const cdouble k,
                             cdouble e[3], cdouble h[3],
                             cdouble de[3][3], cdouble dh[3][3])
{
  cdouble p[1], a[3], dp[3], da[3][3], ddp[3][3], dcurla[3][3];
  GetReducedPotentials_Nearby(S, ne, X0, k, p, a, dp, da, ddp, dcurla);

  cdouble k2=k*k;
  for(int Mu=0; Mu<3; Mu++)
   { 
     int MP1 = (Mu+1)%3, MP2=(Mu+2)%3;

     e[Mu] = a[Mu] + dp[Mu]/k2;
     h[Mu] = da[MP1][MP2] - da[MP2][MP1];

     for(int Nu=0; Nu<3; Nu++)
      { de[Mu][Nu] = da[Mu][Nu] + ddp[Mu][Nu]/k2;
        dh[Mu][Nu] = dcurla[Mu][Nu];
      };
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetReducedFields_Nearby(RWGSurface *S, const int ne,
                             const double X0[3], const cdouble k,
                             cdouble e[3], cdouble h[3])
{
  cdouble de[3][3], dh[3][3];
  GetReducedFields_Nearby(S, ne, X0, k, e, h, de, dh);
}

/***************************************************************/
/* 20141005 get desingularized contributions to far fields.    */
/* TODO: merge me into one of the other very similar routines  */
/* floating around libscuff.                                   */
/***************************************************************/
void GetDSFarFields(RWGSurface *S, const int np, const int iQ,
                    const double X0[3],  const cdouble k, const bool Desingularize,
                    cdouble eDS[3], cdouble hDS[3])
{ 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGPanel *P=S->Panels[np];
  double *Q  = S->Vertices + 3*P->VI[iQ];
  double *V1 = S->Vertices + 3*P->VI[(iQ+1)%3];
  double *V2 = S->Vertices + 3*P->VI[(iQ+2)%3];
  double PreFactor = 2.0*P->Area;

  double A[3], B[3];
  VecSub(V1, Q, A);
  VecSub(V2, V1, B);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumPts;
  int Order=9;
  double *TCR=GetTCR(Order, &NumPts);
  eDS[0]=eDS[1]=eDS[2]=hDS[0]=hDS[1]=hDS[2]=0.0;
  for(int nPt=0, ncp=0; nPt<NumPts; nPt++)
   { 
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     double u=TCR[ncp++];
     double v=TCR[ncp++];
     double w=TCR[ncp++] * PreFactor;
     u+=v;

     double X[3], F[3], R[3], r2=0.0;
     for(int Mu=0; Mu<3; Mu++)
      { 
        F[Mu] = u*A[Mu] + v*B[Mu];
        X[Mu] = F[Mu] + Q[Mu];
        R[Mu] = X0[Mu] - X[Mu];
        r2   += R[Mu]*R[Mu];
      };

     double r=sqrt(r2);
     cdouble ExpFac=exp(II*k*r) / (4.0*M_PI*r);
     cdouble DSExpFac=ExpFac;
     if (Desingularize)
      DSExpFac = ExpRel(II*k*r, 1) / (4.0*M_PI*r);

     cdouble GMuNu[3][3], CMuNu[3][3];
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { 
         GMuNu[Mu][Nu] = ( ( (Mu==Nu) ? 1.0 : 0.0) - R[Mu]*R[Nu]/r2
                         ) * DSExpFac;
       };

     CMuNu[0][0]=CMuNu[1][1]=CMuNu[2][2] = 0.0;
     CMuNu[0][1]= R[2]*ExpFac/r;
     CMuNu[0][2]=-R[1]*ExpFac/r;
     CMuNu[1][2]= R[0]*ExpFac/r;
     CMuNu[1][0]=-R[2]*ExpFac/r;
     CMuNu[2][0]= R[1]*ExpFac/r;
     CMuNu[2][1]=-R[0]*ExpFac/r;
     
     cdouble MIK=-II*k;
     for(int Mu=0; Mu<3; Mu++)
      for(int Nu=0; Nu<3; Nu++)
       { eDS[Mu] += w*GMuNu[Mu][Nu]*F[Nu];
         hDS[Mu] += w*MIK*CMuNu[Mu][Nu]*F[Nu];
       };
   };

}

/***************************************************************/
/* 20141005 get the reduced (e and h) fields, retaining only   */
/* far-field contributions, but using desingularization to     */
/* ensure that the calculation is accurate even when the       */
/* evalution point lies on or near the panel.                  */
/***************************************************************/
void GetReducedFarFields(RWGSurface *S, const int np, const int iQ,
                         const double X0[3],  const cdouble k,
                         cdouble e[3], cdouble h[3])
{
   /*--------------------------------------------------------------*/
   /*- get singular contributions ---------------------------------*/
   /*--------------------------------------------------------------*/
   RWGPanel *P = S->Panels[np];
   bool DeSingularize = (VecDistance(X0,P->Centroid) < 5.0*P->Radius);

   GetDSFarFields(S, np, iQ, X0, k, DeSingularize, e, h);

   if (DeSingularize)
    { 
      double pn[NUMNS], an[NUMNS][3];
      double dpn[NUMNS][3], dan[NUMNS][3][3];
      double ddpn[NUMNS][3][3], dcurlan[NUMNS][3][3];
      GetStaticPanelIntegrals(S, np, iQ, X0,
                              pn, an, dpn, dan, ddpn, dcurlan);

      for(int Mu=0; Mu<3; Mu++)
       e[Mu] += 2.0*dpn[1][Mu]/(4.0*M_PI);
    };
   
   /*--------------------------------------------------------------*/
   /*- put in the RWG basis function prefactors -------------------*/
   /*--------------------------------------------------------------*/
   double Prefac = S->Edges[P->EI[iQ]]->Length / (2.0*P->Area);
   for(int Mu=0; Mu<3; Mu++)
    { e[Mu] *= Prefac;
      h[Mu] *= Prefac;
    };
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetReducedFarFields(RWGSurface *S, const int ne,
                         const double X0[3], const cdouble k,
                         cdouble e[3], cdouble h[3])
{
  RWGEdge *E = S->Edges[ne];
  GetReducedFarFields(S, E->iPPanel, E->PIndex, X0, k, e, h);

  cdouble eM[3], hM[3];
  GetReducedFarFields(S, E->iMPanel, E->MIndex, X0, k, eM, hM);

  for(int Mu=0; Mu<3; Mu++)
   { e[Mu] -= eM[Mu];
     h[Mu] -= hM[Mu];
   };

}

} // namespace scuff
