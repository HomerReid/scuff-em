/*
 * Overlap.cc  -- computation of overlap integrals between RWG basis functions
 *
 * homer reid  -- 5/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"

namespace scuff {

/***************************************************************/
/* Compute the overlap integral between the RWG basis functions*/
/* associated with two edges in an RWG object.                 */
/***************************************************************/
#if 0
double RWGObject::GetOverlapOld(int neAlpha, int neBeta, double *pOTimes)
{ 
  RWGEdge *EAlpha=Edges[neAlpha], *EBeta=Edges[neBeta];

  /*--------------------------------------------------------------*/
  /*- handle the diagonal case -----------------------------------*/
  /*--------------------------------------------------------------*/
  if ( EAlpha==EBeta )
   { 
     double *QP = Vertices + 3*(EAlpha->iQP);
     double *V1 = Vertices + 3*(EAlpha->iV1);
     double *V2 = Vertices + 3*(EAlpha->iV2);
     double *QM = Vertices + 3*(EAlpha->iQM);

     double PArea = Panels[EAlpha->iPPanel]->Area;
     double MArea = Panels[EAlpha->iMPanel]->Area;

     double lA2 = (EAlpha->Length) * (EAlpha->Length);

     double LA[3], LBP[3], LBM[3]; 
     double LAdLBP=0.0, LAdLBM=0.0, lBP2=0.0, lBM2=0.0;
     int i;
     for(i=0; i<3; i++)
      { LA[i]  = V2[i] - V1[i];
        LBP[i] = V1[i] - QP[i];
        LBM[i] = V1[i] - QM[i];

        LAdLBP += LA[i] * LBP[i];
        lBP2   += LBP[i] * LBP[i];
        LAdLBM += LA[i] * LBM[i];
        lBM2   += LBM[i] * LBM[i];
      };
    
     if (pOTimes) 
      *pOTimes=0.0;
     return lA2 * (   ( lA2 + 3.0*lBP2 + 3.0*LAdLBP ) / PArea
                    + ( lA2 + 3.0*lBM2 + 3.0*LAdLBM ) / MArea
                  ) / 24.0;
   };

  /*--------------------------------------------------------------*/
  /*- figure out if there is nonzero overlap ---------------------*/
  /*--------------------------------------------------------------*/
  double Sign, Area, *QA, *QB; 
  int IndexA, IndexB;
  if ( EAlpha->iPPanel == EBeta->iPPanel )
   { Sign = 1.0;
     Area = Panels[EAlpha->iPPanel]->Area;
     QA = Vertices + 3*(EAlpha->iQP);
     QB = Vertices + 3*(EBeta ->iQP);
     IndexA = EAlpha->PIndex;
     IndexB = EBeta->PIndex;
   }
  else if ( EAlpha->iPPanel == EBeta->iMPanel )
   { Sign = -1.0;
     Area = Panels[EAlpha->iPPanel]->Area;
     QA = Vertices + 3*(EAlpha->iQP);
     QB = Vertices + 3*(EBeta ->iQM);
     IndexA = EAlpha->PIndex;
     IndexB = EBeta->MIndex;
   }
  else if ( EAlpha->iMPanel == EBeta->iPPanel )
   { Sign = -1.0;
     Area = Panels[EAlpha->iMPanel]->Area;
     QA = Vertices + 3*(EAlpha->iQM);
     QB = Vertices + 3*(EBeta ->iQP);
     IndexA = EAlpha->MIndex;
     IndexB = EBeta->PIndex;
   }
  else if ( EAlpha->iMPanel == EBeta->iMPanel )
   { Sign = +1.0;
     Area = Panels[EAlpha->iMPanel]->Area;
     QA = Vertices + 3*(EAlpha->iQM);
     QB = Vertices + 3*(EBeta ->iQM);
     IndexA = EAlpha->MIndex;
     IndexB = EBeta->MIndex;
   }
  else
   {  if (pOTimes)
       *pOTimes=0.0; 
      return 0.0;
   };
 /*--------------------------------------------------------------*/
  /*- do the computation -----------------------------------------*/
  /*--------------------------------------------------------------*/
  double *V1 = Vertices + 3*(EAlpha->iV1);
  double *V2 = Vertices + 3*(EAlpha->iV2);
  double *QI=0; // 'QIntermediate' is the common vertex of L_\alpha, L_\beta
  if ( QB == V1 ) 
   QI = V2;
  else if ( QB == V2 ) 
   QI = V1;
  else
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  double lA = EAlpha->Length;
  double lB = EBeta->Length;
  double DotProduct =  (QI[0]-QA[0])*(QB[0]-QI[0])
                      +(QI[1]-QA[1])*(QB[1]-QI[1])
                      +(QI[2]-QA[2])*(QB[2]-QI[2]);

  if (pOTimes)
   { 
     double SignPrime = ( ((IndexB-IndexA+3)%3) == 2) ? 1.0 : -1.0;
     *pOTimes = Sign*SignPrime*lA*lB/6.0;
   };

  return -1.0*Sign*lA*lB*( lA*lA + lB*lB + 3.0*DotProduct ) / (24.0*Area);

}
#endif

/***************************************************************/
/* this is a helper function for GetOverlaps that computes the */
/* contributions of a single panel to the overlap integrals    */
/***************************************************************/
void AddOverlapContributions(RWGObject *O, RWGPanel *P, int iQa, int iQb, 
                             double Sign, double LL, double Overlaps[11])
{
  double *Qa   = O->Vertices + 3*P->VI[ iQa ];
  double *QaP1 = O->Vertices + 3*P->VI[ (iQa+1)%3 ];
  double *QaP2 = O->Vertices + 3*P->VI[ (iQa+2)%3 ];
  double *Qb   = O->Vertices + 3*P->VI[ iQb ];
  double *ZHat = P->ZHat;

  double A[3], B[3], DQ[3];
  VecSub(QaP1, Qa, A);
  VecSub(QaP2, QaP1, B);
  VecSub(Qa, Qb, DQ);

  double ZxA[3], ZxB[3], ZxDQ[3];
  VecCross(ZHat, A, ZxA);
  VecCross(ZHat, B, ZxB);
  VecCross(ZHat, DQ, ZxDQ);

  double PreFac = Sign * LL / (2.0*P->Area);

  double Bullet=0.0, Times=0.0;
  for(int Mu=0; Mu<3; Mu++)
   { Bullet += A[Mu]*(A[Mu]/4.0 + B[Mu]/4.0 + DQ[Mu]/3.0) + B[Mu]*(B[Mu]/12.0 + DQ[Mu]/6.0);
     Times  += (A[Mu]+0.5*B[Mu])*ZxDQ[Mu]/3.0;
   };

  Overlaps[0]  += PreFac*Bullet;
  Overlaps[1]  += PreFac*Times;

  Overlaps[2]  += PreFac * ZHat[0] * Bullet;
  Overlaps[3]  += PreFac * ZHat[0] * 2.0;
  Overlaps[4]  += PreFac * (ZxA[0]/3.0 + ZxB[0]/6.0) * 2.0;

  Overlaps[5]  += PreFac * ZHat[1] * Bullet;
  Overlaps[6]  += PreFac * ZHat[1] * 2.0;
  Overlaps[7]  += PreFac * (ZxA[1]/3.0 + ZxB[1]/6.0) * 2.0;

  Overlaps[8]  += PreFac * ZHat[2] * Bullet;
  Overlaps[9]  += PreFac * ZHat[2] * 2.0;
  Overlaps[10] += PreFac * (ZxA[2]/3.0 + ZxB[2]/6.0) * 2.0;

}


/***************************************************************/
/* entries of output array:                                    */
/*                                                             */
/*  Overlaps [0] = O^{\bullet}_{\alpha\beta}                   */
/*   = \int f_a \cdot f_b                                      */
/*                                                             */
/*  Overlaps [1] = O^{\times}_{\alpha\beta}                    */
/*   = \int f_a \cdot (nHat \times f_b)                        */
/*                                                             */
/*  Overlaps [2] = O^{x,\bullet}_{\alpha\beta}                 */
/*   = \int nHat_x f_a \cdot f_b                               */
/*                                                             */
/*  Overlaps [3] = O^{x,\nabla\nabla}_{\alpha\beta}            */
/*   = \int nHat_x (\nabla \cdot f_a) (\nabla \cdot f_b)       */
/*                                                             */
/*  Overlaps [4] = O^{x,\times\nabla}_{\alpha\beta}            */
/*   = \int (nHat \times \cdot f_a)_x (\nabla \cdot f_b)       */
/*                                                             */
/*  [5,6,7]  = like [2,3,4] but with x-->y                     */
/*  [8,9,10] = like [2,3,4] but with x-->z                     */
/***************************************************************/
void RWGObject::GetOverlaps(int neAlpha, int neBeta, double Overlaps[11])
{
  RWGEdge *EAlpha = Edges[neAlpha];
  RWGEdge *EBeta  = Edges[neBeta];

  RWGPanel *PAlphaP=Panels[EAlpha->iPPanel];
  RWGPanel *PAlphaM=Panels[EAlpha->iMPanel];

  int iQPAlpha = EAlpha->PIndex;
  int iQMAlpha = EAlpha->MIndex;
  int iQPBeta  = EBeta->PIndex;
  int iQMBeta  = EBeta->MIndex;

  double LL = EAlpha->Length * EBeta->Length;

  memset(Overlaps,0,11*sizeof(double));

  if ( EAlpha->iPPanel == EBeta->iPPanel )
   AddOverlapContributions(this, PAlphaP, iQPAlpha, iQPBeta,  1.0, LL, Overlaps);
  if ( EAlpha->iPPanel == EBeta->iMPanel )
   AddOverlapContributions(this, PAlphaP, iQPAlpha, iQMBeta, -1.0, LL, Overlaps);
  if ( EAlpha->iMPanel == EBeta->iPPanel )
   AddOverlapContributions(this, PAlphaM, iQMAlpha, iQPBeta, -1.0, LL, Overlaps);
  if ( EAlpha->iMPanel == EBeta->iMPanel )
   AddOverlapContributions(this, PAlphaM, iQMAlpha, iQMBeta,  1.0, LL, Overlaps);
}

/*****************************************************************/
/* this is a simpler interface to the above routine that returns */
/* the simple overlap integral and sets *pOTimes = crossed       */
/* overlap integral if it is non-NULL                            */
/*****************************************************************/
double RWGObject::GetOverlap(int neAlpha, int neBeta, double *pOTimes)
{
  double Overlaps[11];
  GetOverlaps(neAlpha, neBeta, Overlaps);
  if (pOTimes) *pOTimes=Overlaps[1];
  return Overlaps[0];
}


} // namespace scuff
