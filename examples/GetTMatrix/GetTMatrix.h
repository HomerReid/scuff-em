/*
 * GetTMatrix.h    -- header file for the GetTMatrix program
 *
 * Homer Reid      -- 7/2012
 */

#ifndef GETTMATRIX_H 
#define GETTMATRIX_H

#include <libhrutil.h>
#include <libscuff.h>
#include <libIncField.h>
#include <libhmat.h>

// values for the Type field in the SphericalWave class 
#define SW_MAGNETIC 0
#define SW_ELECTRIC 1

using namespace scuff;

/**********************************************************************/
/* SphericalWave is an implementation of the IncField class that      */
/* describes an incoming spherical wave with given spherical wave     */
/* indices L and M.                                                   */
/**********************************************************************/
class SphericalWave : public IncField
 { 
 public:
   int L, M;        // spherical wave indices 
   int Type;        // either SW_MAGNETIC or SW_ELECTRIC

   SphericalWave(const int L=1, const int M=0, const int Type=SW_MAGNETIC);

   void SetL(int NewL);
   void SetM(int NewM);
   void SetType(int NewType);

   void GetFields(const double X[3], cdouble EH[6]);

 };

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
HVector *GetSphericalMoments(RWGObject *O, cdouble k, int lMax,
                             HVector *KNVector, int BFIndexOffset, 
                             HVector *AVector=0, int NumThreads=0);

#endif // #ifndef GETTMATRIX_H
