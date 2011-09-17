/*
 * GetMatrixElement.cc -- routine that computes the contributions to the 
 *                        BEM matrix from a single pair of basis functions
 * 
 * homer reid          -- 05/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <libhmat.h>
#include <libhrutil.h>

#include "libscuff.h"

#define II cdouble(0,1)

/***************************************************************/
/* NOTE: because the T matrix is symmetric, we only fill in    */
/* its upper triangle, and it is stored in packed format       */ 
/* inside the HMatrix class. also, in this routine, we compute */
/* the matrix entries in 2x2 blocks at a time (corresponding   */
/* to a single pair of basis functions.) when we go to stamp   */
/* a given 2x2 block of entries into the matrix, there is a    */
/* subtlety that affects blocks on the diagonal.               */
/* consider for example what happens when we go to stamp in    */
/* the 2x2 block corresponding to the interactions of basis    */
/* function #7 with itself. this block corresponds to elements */
/*  (14,14) (14,15)                                            */
/*  (15,14) (15,15)                                            */
/* of the larger matrix. but notice that (15,14) lies in the   */
/* lower triangle of the matrix. hence, if we try to set this  */
/* entry to a value, that value will override the value that   */
/* is already set for (14,15) (because of the way HMatrix is   */
/* implemented.) hence, whenever we go to stamp in a 2x2 block */
/* corresponding to ne==nep in the loop below, we need NOT to  */
/* stamp in the (2*ne+1, 2*nep) element.                       */
/***************************************************************/

/***************************************************************/
/* note: iwe = i * omega * epsilon *****************************/
/* note: iwu = i * omega * mu       ****************************/
/***************************************************************/
#define STAMP_SET 0
#define STAMP_ADD 1
static void Stamp(HMatrix *M, 
                  int OffsetA, int NBFA, int OffsetB, int NBFB,
                  cdouble iwe, cdouble iwu, double Sign,
                  cdouble *L, int SetOrAdd)
{
  cdouble LL[3];

  LL[0] = Sign*L[0];
  LL[1] = Sign*L[1];
  LL[2] = Sign*L[2];

  if ( SetOrAdd == STAMP_ADD )
   { M->AddEntry(OffsetA, OffsetB, iwu*LL[0] + LL[1]/iwe );
     if (NBFB==2)
      M->AddEntry(OffsetA, OffsetB+1, -LL[2] );
     if (NBFA==2 && (M->StorageType==LHM_NORMAL || OffsetA!=OffsetB) )
      M->AddEntry(OffsetA+1, OffsetB, -LL[2] );
     if (NBFA==2 && NBFB==2)
      M->AddEntry(OffsetA+1, OffsetB+1, -(iwe*LL[0] + LL[1]/iwu ) );
   }
  else
   { M->SetEntry(OffsetA, OffsetB, iwu*LL[0] + LL[1]/iwe );
     if (NBFB==2)
      M->SetEntry(OffsetA, OffsetB+1, -LL[2] );
     if (NBFA==2)
      M->SetEntry(OffsetA+1, OffsetB, -LL[2] );
     if (NBFA==2 && NBFB==2)
      M->SetEntry(OffsetA+1, OffsetB+1, -( iwe*LL[0] + LL[1]/iwu ) );
   }

}

/***************************************************************/
/* important note: this routine uses the values of EpsThisFreq */
/* and MuThisFreq that are cached inside Oa, Ob, and (this);   */
/* you need to make sure those values are correct before       */
/* calling this routine.                                       */
/***************************************************************/
void RWGGeometry::StampMatrixElements(void *pLFW, 
                                      RWGObject *Oa, int nea, int OffsetA, 
                                      RWGObject *Ob, int neb, int OffsetB, 
                                      double Frequency, int RealFreq,
                                      int NumTorqueAxes, double *GammaMatrix,
                                      HMatrix *M, 
                                      HMatrix *dMdX, HMatrix *dMdY, HMatrix *dMdZ, 
                                      HMatrix *dMdT1, HMatrix *dMdT2, HMatrix *dMdT3)
{

  /***************************************************************/
  /* first switch off based on the containership relation between*/
  /* objects Oa and Ob to determine the medium through which     */
  /* the objects interact                                        */
  /***************************************************************/
  cdouble Eps;
  double Mu;
  double Sign;
  if ( Oa==Ob->ContainingObject )
   {
     /*--------------------------------------------------------------*/
     /*- object b is contained in object a --------------------------*/
     /*--------------------------------------------------------------*/
     //Oa->MP->GetEpsMu(Frequency, RealFreq, &Eps, &Mu);
     Eps=Oa->EpsThisFreq;
     Mu=Oa->MuThisFreq;
     Sign = -1.0;
   }
  else if ( Ob==Oa->ContainingObject )
   {
     /*--------------------------------------------------------------*/
     /*- object a is contained in object b --------------------------*/
     /*--------------------------------------------------------------*/
     //Ob->MP->GetEpsMu(Frequency, RealFreq, &Eps, &Mu);
     Eps=Ob->EpsThisFreq;
     Mu=Ob->MuThisFreq;
     Sign = -1.0;
   }
  else if ( Oa->ContainingObject==Ob->ContainingObject )
   {   
     /*--------------------------------------------------------------*/
     /*- the objects are both contained in the same object (or in the*/
     /*- external medium). note that this includes the case Oa==Ob.  */
     /*--------------------------------------------------------------*/
     if ( Oa->ContainingObject==0 )
      { //MP->GetEpsMu(Frequency, RealFreq, &Eps, &Mu);
        Eps=EpsThisFreq;
        Mu=MuThisFreq;
      }
     else
      { //Oa->ContainingObject->MP->GetEpsMu(Frequency, RealFreq, &Eps, &Mu);
        Eps=Oa->EpsThisFreq;
        Mu=Oa->MuThisFreq;
      }
     Sign=+1.0;
   }
  else 
   {
     /*--------------------------------------------------------------*/
     /*- if none of the above are true, the objects do not interact. */
     /*--------------------------------------------------------------*/
     return;
   }

  /***************************************************************/
  /* figure out whether we need derivatives of matrix elements   */
  /* (if so we have to pass nonzero buffers to GetLFunctions)    */
  /***************************************************************/
  cdouble L[3], GradLBuffer[9], dLdTBuffer[9];
  cdouble *GradL = (dMdX || dMdY || dMdZ) ? GradLBuffer : 0;
  cdouble *dLdT  = NumTorqueAxes>0 ? dLdTBuffer : 0;

  int NeedCross; 
  if ( Oa->MP->IsPEC() && Ob->MP->IsPEC() )
   NeedCross=0;
  else
   NeedCross=1;

  int NBFA, NBFB;
  NBFA = Oa->MP->IsPEC() ? 1 : 2;
  NBFB = Ob->MP->IsPEC() ? 1 : 2;

  /***************************************************************/
  /* compute the edge--edge interaction through the medium       */
  /* whose properties (Eps,Mu) we identified above               */
  /***************************************************************/
  cdouble K;
  K=RealFreq ? csqrt2(Eps*Mu)*Frequency : II*sqrt(real(Eps)*Mu)*Frequency;
  GetLFunctions(pLFW, Oa, nea, Ob, neb, K, NeedCross, 
                NumTorqueAxes, GammaMatrix, L, GradL, dLdT);

  /***************************************************************/
  /* stamp in those contributions ********************************/
  /***************************************************************/
  cdouble iwe, iwu;
  if(RealFreq)
   { iwe=II*Frequency*Eps;
     iwu=II*Frequency*Mu;
   }
  else
   { iwe=-1.0*Frequency*Eps;
     iwu=-1.0*Frequency*Mu;
   };
  Stamp(M, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, L, 0);
  if (dMdX) 
   Stamp(dMdX, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, GradL+0, 0);
  if (dMdY)
   Stamp(dMdY, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, GradL+3, 0);
  if (dMdZ)
   Stamp(dMdZ, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, GradL+6, 0);
  if (dMdT1)
   Stamp(dMdT1, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, dLdT+0, 0);
  if (dMdT2)
   Stamp(dMdT2, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, dLdT+3, 0);
  if (dMdT3)
   Stamp(dMdT3, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, dLdT+6, 0);

  /***************************************************************/
  /* if the two BFs lie on the surface of the same object and    */
  /* the object is non-PEC, then there is an additional          */
  /* contribution coming from the interior medium                */
  /***************************************************************/
  if (Oa==Ob && !(Oa->MP->IsPEC()) )
   {
     //Oa->MP->GetEpsMu(Frequency, RealFreq, &Eps, &Mu);
     Eps=Oa->EpsThisFreq;
     Mu=Oa->MuThisFreq;
     K=RealFreq ? csqrt2(Eps*Mu)*Frequency : II*sqrt(real(Eps)*Mu)*Frequency;
     GetLFunctions(pLFW, Oa, nea, Ob, neb, K, NeedCross,
                   NumTorqueAxes, GammaMatrix, L, GradL, dLdT);

     if(RealFreq)
      { iwe=II*Frequency*Eps;
        iwu=II*Frequency*Mu;
      }
     else
      { iwe=-1.0*Frequency*Eps;
        iwu=-1.0*Frequency*Mu;
      };

     Stamp(M, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, L, 1);
     if (dMdX) 
      Stamp(dMdX, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, GradL+0, 1);
     if (dMdY)
      Stamp(dMdY, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, GradL+3, 1);
     if (dMdZ)
      Stamp(dMdZ, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, GradL+6, 1);
     if (dMdT1)
      Stamp(dMdT1, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, dLdT+0, 1);
     if (dMdT2)
      Stamp(dMdT2, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, dLdT+3, 1);
     if (dMdT3)
      Stamp(dMdT3, OffsetA, NBFA, OffsetB, NBFB, iwe, iwu, Sign, dLdT+6, 1);
   
   };

}
