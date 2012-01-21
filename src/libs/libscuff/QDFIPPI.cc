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

/***************************************************************/
/* 'vertex less than.' returns 1 if V1<=V2, 0 otherwise.       */
/* vertices are sorted using a fairly obvious sorting scheme.  */
/***************************************************************/
static int VLT(double *V1, double *V2)
{
  double DV;

  DV=V1[0]-V2[0];
  if ( fabs(DV) > 1.0e-6*fabs(V1[0]) )
   return DV<0.0 ? 1 : 0;

  DV=V1[1]-V2[1];
  if ( fabs(DV) > 1.0e-6*fabs(V1[1]) )
   return DV<0.0 ? 1 : 0;
 
  DV=V1[2]-V2[2];
  if ( fabs(DV) > 1.0e-6*fabs(V1[2]) )
   return DV<0.0 ? 1 : 0;

  return 1;
}

/***************************************************************/
/* given two sets of panel vertices, this routine puts them    */
/* in a canonical ordering.                                    */
/*                                                             */
/*                                                             */
/* my 'canonical ordering' for panel vertices obeys the        */
/* following properties.                                       */
/*                                                             */
/*  (a) if the panels have no common vertices, then the        */
/*      vertices of each panel are sorted in ascending order   */
/*      and the smallest vertex of OVa is less than OVb, i.e.  */
/*      we have                                                */ 
/*        OVa[0] < OVa[1] < OVa[2]                             */ 
/*        OVb[0] < OVb[1] < OVb[2]                             */
/*      and                                                    */
/*        OVa[0] < OVb[0].                                     */
/*                                                             */
/*  (b) if the panels have any common vertices, then those     */
/*      common vertices appear first in both lists; within     */
/*      each panel, the subset of vertices that are common     */
/*      with the other panel are sorted in ascending order,    */
/*      the subset of vertices that are not common are         */
/*      separately sorted in ascending order, and the first    */
/*      noncommon vertex in OVa is less than the first         */
/*      noncommon vertex in OVb.                               */
/*                                                             */
/*      more specifically, what this means is:                 */
/*                                                             */
/*       (b1) if there are three common vertices, we have      */
/*                                                             */
/*             OVa[0] < OVa[1] < OVa[2]                        */
/*                                                             */
/*            and                                              */
/*                                                             */
/*             OVb[i] = OVa[i]  for i=1,2,3.                   */
/*                                                             */
/*       (b2) if there are two common vertices, we have        */
/*                                                             */
/*             OVa[0] < OVa[1],                                */
/*                                                             */
/*             OVb[i] = OVa[i]  for i=1,2                      */
/*                                                             */
/*             and OVa[2] < OVb[2].                            */
/*                                                             */
/*       (b3) if there is one common vertex, we have           */
/*                                                             */
/*             OVa[0] = OVb[0],                                */
/*             OVa[1] < OVa[2],                                */
/*             OVb[1] < OVb[2],                                */
/*            and                                              */
/*             OVa[1] < OVb[1].                                */
/***************************************************************/
int CanonicallyOrderVertices(double **Va, double **Vb,
                             double **OVa, double **OVb, int *pncv)
{
  double *TV, *TVa[3], *TVb[3];
  TVa[0]=Va[0]; TVa[1]=Va[1]; TVa[2]=Va[2];
  TVb[0]=Vb[0]; TVb[1]=Vb[1]; TVb[2]=Vb[2];

  int ncv=AssessPanelPair(TVa, TVb);
  if (pncv) *pncv=ncv;

  int iMina, iMeda, iMaxa;
  int iMinb, iMedb, iMaxb;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( ncv == 3 )
   { 
     if ( VLT(TVa[0], TVa[1]) )
      { iMina=0; iMaxa=1; }
     else
      { iMina=1; iMaxa=0; }

     if ( VLT(TVa[2], TVa[iMina]) )
      { iMeda=iMina; iMina=2; }
     else if ( VLT(TVa[2], TVa[iMaxa]) )
      { iMeda=2; }
     else
      { iMeda=iMaxa; iMaxa=2; }

     OVa[0]=TVa[iMina]; OVa[1]=TVa[iMeda]; OVa[2]=TVa[iMaxa];
     OVb[0]=TVb[iMina]; OVb[1]=TVb[iMeda]; OVb[2]=TVb[iMaxa];
     return 0;
   }
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  else if ( ncv == 2 )
   {
     if ( VLT(TVa[0], TVa[1]) )
      { OVa[0]=TVa[0]; OVa[1]=TVa[1];
        OVb[0]=TVb[0]; OVb[1]=TVb[1];
      }
     else
      { OVa[0]=TVa[1]; OVa[1]=TVa[0];
        OVb[0]=TVb[1]; OVb[1]=TVb[0];
      };

     if ( VLT(TVa[2], TVb[2]) )
      { OVa[2]=TVa[2];
        OVb[2]=TVb[2];
        return 0;
      }
     else
      { TV=OVa[0]; OVa[0]=OVb[0]; OVb[0]=TV;
        TV=OVa[1]; OVa[1]=OVb[1]; OVb[1]=TV;
        OVa[2]=TVb[2];
        OVb[2]=TVa[2];
        return 1;
      };  
   }
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  else if ( ncv == 1 )
   {
     if ( VLT(TVa[1], TVa[2]) )
      { TV=TVa[1]; TVa[1]=TVa[2]; TVa[2]=TV; };

     if ( VLT(TVb[1], TVb[2]) )
      { TV=TVb[1]; TVb[1]=TVb[2]; TVb[2]=TV; };

     if ( VLT(TVa[0], TVb[0]) )
      { OVa[0]=TVa[0]; OVa[1]=TVa[1]; OVa[2]=TVa[2];
        OVb[0]=TVb[0]; OVb[1]=TVb[1]; OVb[2]=TVb[2];
        return 0;
      }
     else
      { OVa[0]=TVb[0]; OVa[1]=TVb[1]; OVa[2]=TVb[2];
        OVb[0]=TVa[0]; OVb[1]=TVa[1]; OVb[2]=TVa[2];
        return 1;
      };
   }
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  else // ( ncv == 0 )
   { 
     if ( VLT(TVa[0], TVa[1]) )
      { iMina=0; iMaxa=1; }
     else
      { iMina=1; iMaxa=0; }

     if ( VLT(TVa[2], TVa[iMina]) )
      { iMeda=iMina; iMina=2; }
     else if ( VLT(TVa[2], TVa[iMaxa]) )
      { iMeda=2; }
     else
      { iMeda=iMaxa; iMaxa=2; }

     if ( VLT(TVb[0], TVb[1]) )
      { iMinb=0; iMaxb=1; }
     else
      { iMinb=1; iMaxb=0; }

     if ( VLT(TVb[2], TVb[iMinb]) )
      { iMedb=iMinb; iMinb=2; }
     else if ( VLT(TVb[2], TVb[iMaxb]) )
      { iMedb=2; }
     else
      { iMedb=iMaxb; iMaxb=2; }

     if ( VLT(TVa[iMina], TVa[iMinb]) )
      { OVa[0] = TVa[iMina]; OVa[1] = TVa[iMeda]; OVa[2] = TVa[iMaxa]; 
        OVb[0] = TVb[iMina]; OVb[1] = TVb[iMeda]; OVb[2] = TVb[iMaxa]; 
        return 0;
      }
     else
      { OVb[0] = TVa[iMina]; OVb[1] = TVa[iMeda]; OVb[2] = TVa[iMaxa]; 
        OVa[0] = TVb[iMina]; OVa[1] = TVb[iMeda]; OVa[2] = TVb[iMaxa]; 
        return 1;
      };

   }; // if ncv ... else
  
}

/*--------------------------------------------------------------*/
/*- routine for computing Q-dependent FIPPIs.                   */
/*--------------------------------------------------------------*/
void GetQDFIPPIData(double **Va, double *Qa, double **Vb, double *Qb, 
                    void *opFDT, QDFIPPIData *QDFD)
{
  int Flipped;
  double *OVa[3], *OQa, *OVb[3], *OQb;  // 'ordered vertices'
  QIFIPPIData MyQIFD, *QIFD;

  Flipped=CanonicallyOrderVertices(Va, Vb, OVa, OVb, &ncv);

  if (Flipped) 
   { Oqa=Qb;
     Oqb=Qa;
   }
  else
   { Oqa=Qa;
     Oqb=Qb;
   };

  /*--------------------------------------------------------------*/
  /*- get the Q-independent FIPPIs by looking them up in a table  */
  /*- if we have one or by computing them if we don't             */
  /*--------------------------------------------------------------*/
  if (opFDT)
   { 
     QIFD=((FIPPITable *)opFT)->GetQIFIPPIData(OVa, OVb, ncv);
   }
  else
   { 
     ComputeQIFIPPIData(OVa, OVb, ncv, &MyQIFD);
     QIFD=&MyQIFD;
   };
  
  /*--------------------------------------------------------------*/
  /*- now assemble the Q-dependent FIPPIs.                        */
  /*--------------------------------------------------------------*/
  double A[3],  B[3],  V0mQ[3];
  double AP[3], BP[3], V0PmQP[3]; 
  double Delta;

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

  double V0QaxQb[3], QamQb[3], QaxQb[3], V0mV0P[3], Scratch[3];
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
