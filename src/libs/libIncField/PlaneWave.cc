/*
 * PlaneWave.cc   -- plane wave implementation of IncField
 *
 * homer reid     -- 11/2009 -- 2/2012
 */

#include <string.h>
#include <math.h>

#include "libIncField.h"
#define II cdouble(0.0,1.0)

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
PlaneWave::PlaneWave(const cdouble pE0[3], const double pnHat[3],
		     const char *Label)
{
  memcpy(E0, pE0, 3*sizeof(cdouble));
  memcpy(nHat, pnHat, 3*sizeof(double));
  SetObjectLabel(Label);
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void PlaneWave::GetFields(const double X[3], cdouble EH[6])
{

  cdouble K=sqrt(Eps*Mu) * Omega;
  cdouble Z=ZVAC*sqrt(Mu/Eps);
  cdouble ExpFac=exp(II*K*(nHat[0]*X[0] + nHat[1]*X[1] + nHat[2]*X[2]));

  EH[0] = E0[0] * ExpFac;
  EH[1] = E0[1] * ExpFac;
  EH[2] = E0[2] * ExpFac;

  /* H = (nHat \cross E) / Z */
  EH[3] = (nHat[1]*EH[2] - nHat[2]*EH[1]) / Z;
  EH[4] = (nHat[2]*EH[0] - nHat[0]*EH[2]) / Z ;
  EH[5] = (nHat[0]*EH[1] - nHat[1]*EH[0]) / Z;
  
} 
