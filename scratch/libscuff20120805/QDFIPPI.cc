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
 * QDFIPPI.cc  -- libscuff routines for computing Q-dependent 
 *             -- frequency-independent panel-panel integrals 
 *             -- (QDFIPPIs)
 * 
 * homer reid  -- 11/2005 -- 1/2012
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "libscuffInternals.h"

namespace scuff {

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

/*--------------------------------------------------------------*/
/*- routine for computing Q-dependent FIPPIs.                   */
/*-                                                             */
/*- note: on entry, ncv is the number of common vertices, and   */
/*-       if ncv>0 then the Va and Vb arrays are assumed to be  */
/*-       in the order computed by AssessPanelPair().           */
/*--------------------------------------------------------------*/
void GetQDFIPPIData(double **Va, double *Qa,
                    double **Vb, double *Qb,
                    int ncv, void *opFC, QDFIPPIData *QDFD)
{
  QIFIPPIData MyQIFD, *QIFD;
  double *OVa[3], *OQa, *OVb[3], *OQb;
  int Flipped;

  /*--------------------------------------------------------------*/
  /*- get the Q-independent FIPPIs by looking them up in a table  */
  /*- if we have one or by computing them if we don't             */
  /*--------------------------------------------------------------*/
  if (opFC)
   { 
     Flipped=CanonicallyOrderVertices(Va, Vb, ncv, OVa, OVb);
     if (Flipped)
      { OQa=Qb; OQb=Qa; }
     else
      { OQa=Qa; OQb=Qb; };

     QIFD=((FIPPICache *)opFC)->GetQIFIPPIData(OVa, OVb, ncv);
   }
  else
   { 
     Flipped=0;
     OVa[0]=Va[0]; OVa[1]=Va[1]; OVa[2]=Va[2]; OQa=Qa;
     OVb[0]=Vb[0]; OVb[1]=Vb[1]; OVb[2]=Vb[2]; OQb=Qb;

     ComputeQIFIPPIData(OVa, OVb, ncv, &MyQIFD);
     QIFD=&MyQIFD;
   };
  
  /*--------------------------------------------------------------*/
  /*- now assemble the Q-dependent FIPPIs.                        */
  /*--------------------------------------------------------------*/
  double A[3],  B[3],  V0mQ[3];
  double AP[3], BP[3], V0PmQP[3]; 

  VecSub(OVa[1], OVa[0], A);
  VecSub(OVa[2], OVa[1], B);
  VecSub(OVa[0], OQa, V0mQ);

  VecSub(OVb[1], OVb[0], AP);
  VecSub(OVb[2], OVb[1], BP);
  VecSub(OVb[0], OQb, V0PmQP);

  double VdVP = VecDot(V0mQ, V0PmQP);
  double VdAP = VecDot(V0mQ, AP);
  double VdBP = VecDot(V0mQ, BP);
  double AdVP = VecDot(A, V0PmQP);
  double AdAP = VecDot(A, AP);
  double AdBP = VecDot(A, BP);
  double BdVP = VecDot(B, V0PmQP);
  double BdAP = VecDot(B, AP);
  double BdBP = VecDot(B, BP);

  double QamQb[3], QaxQb[3], V0mV0P[3], Scratch[3];
  double *V0=OVa[0], *V0P=OVb[0];
  VecSub(V0, V0P, V0mV0P);
  VecSub(Qa, Qb, QamQb);
  VecCross(Qa, Qb, QaxQb);
  double hTimes0000 = VecDot(QaxQb, V0mV0P)  + VecDot(QamQb,  VecCross(V0,V0P, Scratch));
  double hTimes1000 = VecDot(QaxQb, A     )  + VecDot(QamQb,  VecCross(A,V0P,Scratch));
  double hTimes0100 = VecDot(QaxQb, B     )  + VecDot(QamQb,  VecCross(B,V0P,Scratch));
  double hTimes0010 = -VecDot(QaxQb, AP    ) + VecDot(QamQb,  VecCross(V0,AP,Scratch));
  double hTimes0001 = -VecDot(QaxQb, BP    ) + VecDot(QamQb,  VecCross(V0,BP,Scratch));
  double hTimes1010 = VecDot(QamQb,  VecCross(A,AP,Scratch));
  double hTimes1001 = VecDot(QamQb,  VecCross(A,BP,Scratch));
  double hTimes0110 = VecDot(QamQb,  VecCross(B,AP,Scratch));
  double hTimes0101 = VecDot(QamQb,  VecCross(B,BP,Scratch));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hTimesRM3 = VecDot(QaxQb, QIFD->xMxpRM3) + VecDot(QamQb, QIFD->xXxpRM3);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hDotRM1 =   QIFD->uvupvpRM1[_ONE] * VdVP
                  + QIFD->uvupvpRM1[_UP ] * VdAP
                  + QIFD->uvupvpRM1[_VP ] * VdBP
                  + QIFD->uvupvpRM1[_U  ] * AdVP
                  + QIFD->uvupvpRM1[_UUP] * AdAP
                  + QIFD->uvupvpRM1[_UVP] * AdBP
                  + QIFD->uvupvpRM1[_V  ] * BdVP
                  + QIFD->uvupvpRM1[_VUP] * BdAP
                  + QIFD->uvupvpRM1[_VVP] * BdBP;

  QDFD->hNablaRM1 =  4.0*QIFD->uvupvpRM1[_ONE];

  QDFD->hTimesRM1 =   QIFD->uvupvpRM1[_ONE] * hTimes0000
                    + QIFD->uvupvpRM1[_U  ] * hTimes1000
                    + QIFD->uvupvpRM1[_V  ] * hTimes0100
                    + QIFD->uvupvpRM1[_UP ] * hTimes0010
                    + QIFD->uvupvpRM1[_VP ] * hTimes0001
                    + QIFD->uvupvpRM1[_UUP] * hTimes1010
                    + QIFD->uvupvpRM1[_UVP] * hTimes1001
                    + QIFD->uvupvpRM1[_VUP] * hTimes0110
                    + QIFD->uvupvpRM1[_VVP] * hTimes0101;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hDotR0  =   uvupvpR0_ONE * VdVP
                  + uvupvpR0_UP  * VdAP
                  + uvupvpR0_VP  * VdBP
                  + uvupvpR0_U   * AdVP
                  + uvupvpR0_UUP * AdAP
                  + uvupvpR0_UVP * AdBP
                  + uvupvpR0_V   * BdVP
                  + uvupvpR0_VUP * BdAP
                  + uvupvpR0_VVP * BdBP;

  QDFD->hNablaR0 = 1.0;

  QDFD->hTimesR0  =   uvupvpR0_ONE * hTimes0000
                    + uvupvpR0_U   * hTimes1000
                    + uvupvpR0_V   * hTimes0100
                    + uvupvpR0_UP  * hTimes0010
                    + uvupvpR0_VP  * hTimes0001
                    + uvupvpR0_UUP * hTimes1010
                    + uvupvpR0_UVP * hTimes1001
                    + uvupvpR0_VUP * hTimes0110
                    + uvupvpR0_VVP * hTimes0101;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  QDFD->hDotR1 =   QIFD->uvupvpR1[_ONE] * VdVP
                 + QIFD->uvupvpR1[_UP ] * VdAP
                 + QIFD->uvupvpR1[_VP ] * VdBP
                 + QIFD->uvupvpR1[_U  ] * AdVP
                 + QIFD->uvupvpR1[_UUP] * AdAP
                 + QIFD->uvupvpR1[_UVP] * AdBP
                 + QIFD->uvupvpR1[_V  ] * BdVP
                 + QIFD->uvupvpR1[_VUP] * BdAP
                 + QIFD->uvupvpR1[_VVP] * BdBP;

  QDFD->hNablaR1 =  4.0*QIFD->uvupvpR1[_ONE];

  QDFD->hTimesR1 =   QIFD->uvupvpR1[_ONE] * hTimes0000
                   + QIFD->uvupvpR1[_U  ] * hTimes1000
                   + QIFD->uvupvpR1[_V  ] * hTimes0100
                   + QIFD->uvupvpR1[_UP ] * hTimes0010
                   + QIFD->uvupvpR1[_VP ] * hTimes0001
                   + QIFD->uvupvpR1[_UUP] * hTimes1010
                   + QIFD->uvupvpR1[_UVP] * hTimes1001
                   + QIFD->uvupvpR1[_VUP] * hTimes0110
                   + QIFD->uvupvpR1[_VVP] * hTimes0101;

  /*--------------------------------------------------------------*/
  /*-- NOTE: the R2 integrals may actually be done analytically,  */
  /*--       but the calculation is so cumbersome that i do them  */
  /*--       numerically for now. in the future, maybe explore    */
  /*--       the costs and benefits of doing them analytically.   */
  /*--------------------------------------------------------------*/
  QDFD->hDotR2 =   QIFD->uvupvpR2[_ONE] * VdVP
                 + QIFD->uvupvpR2[_UP ] * VdAP
                 + QIFD->uvupvpR2[_VP ] * VdBP
                 + QIFD->uvupvpR2[_U  ] * AdVP
                 + QIFD->uvupvpR2[_UUP] * AdAP
                 + QIFD->uvupvpR2[_UVP] * AdBP
                 + QIFD->uvupvpR2[_V  ] * BdVP
                 + QIFD->uvupvpR2[_VUP] * BdAP
                 + QIFD->uvupvpR2[_VVP] * BdBP;

  QDFD->hNablaR2 =  4.0*QIFD->uvupvpR2[_ONE];

  if (Flipped)
   { QDFD->hTimesRM3 *= -1.0; 
     QDFD->hTimesRM1 *= -1.0; 
     QDFD->hTimesR0  *= -1.0; 
     QDFD->hTimesR1  *= -1.0; 
   };

}

} // namespace scuff
