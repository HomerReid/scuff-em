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
double RWGObject::GetOverlap(int neAlpha, int neBeta, double *pOTimes)
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

/***************************************************************/
/*                                                             */
/* entries of output array:                                    */
/*                                                             */
/***************************************************************/
void AddOverlapContributions(RWGEdge *E, RWGEdge 
{
  PreFac = Sign*EAlpha->Length*EBeta->Length / P->Area;
  OTimes       += PreFac * ( VecDot(ApB02,ApB02)/8.0 + VecDot(ApB02,DeltaQ)/6.0 )
  OiBullet     += PreFac * ( VecDot(ApB02,ApB02)/8.0 + VecDot(ApB02,DeltaQ)/6.0 )
  OiNablaNabla += PreFac * ZHat[i] / 4.0;
  OiTimesNabla += PreFac * ( B[i]/3.0 - A[i]/6.0 );
}


void RWGObject::GetOverlaps(int neAlpha, int neBeta, int i, double Overlaps[5])
{
  RWGEdge *EAlpha = Edges[neAlpha];
  RWGEdge *EBeta  = Edges[neBeta];

  double LL = EAlpha->L

  memset(Overlaps,0,5*sizeof(double));
  if ( EAlpha->iPPanel == EBeta->iPPanel )
   { AddOverlapContributions(LL, Panels[EAlpha->iPPanel], 

  O            += PreFac * ( VecDot(ApB02,ApB02)/8.0 + VecDot(ApB02,DeltaQ)/6.0 )
}

} // namespace scuff
