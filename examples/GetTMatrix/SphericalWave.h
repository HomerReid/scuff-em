/*
 * SphericalWave.h -- definition for the SphericalWave implementation
 *                 -- of the IncField class 
 *
 * Homer Reid      -- 7/2012
 */

#ifndef SPHERICALWAVE_H
#define SPHERICALWAVE_H 

#include <libhrutil.h>
#include <libIncField.h>

#define SW_MAGNETIC 0
#define SW_ELECTRIC 1

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
class SphericalWave : public IncField
 { 
 public:
   int L, M;        // spherical wave indices 
   int Type;        // either SW_ELECTRIC or SW_MAGNETIC

   SphericalWave(const cdouble E0[3], const double nHat[3], const char *Label = 0);

   void SetL(int NewL);
   void SetM(int NewM);
   void SetType(int NewType);

   void GetFields(const double X[3], cdouble EH[6]);

 };

#endif // #ifndef SPHERICALWAVE_H
