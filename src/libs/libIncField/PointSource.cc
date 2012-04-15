/*
 * PointSource.cc -- point source implementation of IncField
 *
 * homer reid     -- 11/2009 -- 2/2012
 */

#include "libIncField.h"
#include <math.h>
#include <string.h>

#define II cdouble(0.0,1.0)

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
PointSource::PointSource(const double pX0[3], const cdouble pP[3], int pType,
			 const char *Label)
{
  memcpy(X0, pX0, 3*sizeof(double));
  memcpy(P,  pP, 3*sizeof(cdouble));
  Type=pType;
  SetObjectLabel(Label);
}

/**********************************************************************/
/* fields of a complex point source.                                  */
/*                                                                    */
/* NOTE: for the default case of an electric dipole, the quantity P   */
/* in the PointSource structure is assumed to be the dipole           */
/* moment divided by \epsilon_0, which means that P has units of      */
/* voltage*length^2.                                                  */
/*                                                                    */
/**********************************************************************/
void PointSource::GetFields(const double X[3], cdouble EH[6])
{
  /* construct R, RHat, etc. */
  double RHat[3], R;
  RHat[0]=X[0] - X0[0];
  RHat[1]=X[1] - X0[1];
  RHat[2]=X[2] - X0[2];
  R=sqrt(  RHat[0]*RHat[0] + RHat[1]*RHat[1] + RHat[2]*RHat[2] );
  RHat[0]/=R;
  RHat[1]/=R;
  RHat[2]/=R;

  cdouble PDotR, RCrossP[3];
  PDotR=P[0]*RHat[0] + P[1]*RHat[1] + P[2]*RHat[2];
  RCrossP[0]= RHat[1]*P[2] - RHat[2]*P[1];
  RCrossP[1]= RHat[2]*P[0] - RHat[0]*P[2];
  RCrossP[2]= RHat[0]*P[1] - RHat[1]*P[0];
  
  cdouble k      = Omega*sqrt(Eps*Mu);
  cdouble ikr    = II*k*R;
  cdouble ikr2   = ikr*ikr;
  cdouble ExpFac = k*k*exp(ikr) / (4.0*M_PI*R);

  cdouble Z      = ZVAC*sqrt(Mu/Eps);

  /* compute the various scalar quantities in the point source formulae */
  cdouble Term1=  1.0 - 1.0/ikr + 1.0/ikr2; 
  cdouble Term2= (-1.0 + 3.0/ikr - 3.0/ikr2) * PDotR; 
  cdouble Term3= (1.0 - 1.0/ikr);

  /* now assemble everything based on source type */
  if ( Type == LIF_ELECTRIC_DIPOLE )
   { 
     EH[0]=ExpFac*( Term1*P[0] + Term2*RHat[0] );
     EH[1]=ExpFac*( Term1*P[1] + Term2*RHat[1] );
     EH[2]=ExpFac*( Term1*P[2] + Term2*RHat[2] );

     EH[3]=ExpFac*Term3*RCrossP[0] / Z;
     EH[4]=ExpFac*Term3*RCrossP[1] / Z;
     EH[5]=ExpFac*Term3*RCrossP[2] / Z;
   }
  else // ( Type == LIF_TYPE_PSMC )
   { 
     EH[0]=-1.0*ExpFac*Term3*RCrossP[0] / Z;
     EH[1]=-1.0*ExpFac*Term3*RCrossP[1] / Z;
     EH[2]=-1.0*ExpFac*Term3*RCrossP[2] / Z;

     EH[3]=ExpFac*( Term1*P[0] + Term2*RHat[0] ) / (Z*Z);
     EH[4]=ExpFac*( Term1*P[1] + Term2*RHat[1] ) / (Z*Z);
     EH[5]=ExpFac*( Term1*P[2] + Term2*RHat[2] ) / (Z*Z);
   };

}
