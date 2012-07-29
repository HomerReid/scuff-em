/*
 * SphericalWave .cc   -- spherical wave implementation of IncField
 *
 * homer reid          -- 11/2009 -- 2/2012
 */

#include <string.h>
#include <math.h>

#include <libIncField.h>
#include "SphericalWave.h"

#define II cdouble(0.0,1.0)

/**********************************************************************/
/* forward definitions of a couple of utility routines defined below  */
/**********************************************************************/

// vector spherical harmonics 
void GetLYlm(int L, int M, double Theta, double Phi, cdouble LYlm[3]);

// spherical hankel functions 
void Gethl(int l, cdouble z, cdouble *hl, cdouble *hlPrime);

/**********************************************************************/
/* constructor and field-setting routines *****************************/
/**********************************************************************/
SphericalWave::SphericalWave(int NewL, int NewM, int NewType)
{
  L=NewL;
  M=NewM;
  Type=NewType;
}

void SphericalWave::SetL(int NewL) { L=NewL; }
void SphericalWave::SetM(int NewM) { M=NewM; }
void SphericalWave::SetType(int NewType) { Type=NewType; }

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void SphericalWave::GetFields(const double X[3], cdouble EH[6])
{

  cdouble K=sqrt(Eps*Mu) * Omega;
  cdouble Z=ZVAC*sqrt(Mu/Eps);

  /**********************************************************************/
  /* convert cartesian coordinates to spherical *************************/
  /**********************************************************************/
  double Rho2 = X[0]*X[0] + X[1]*X[1], Rho=sqrt(Rho2);
  double r = sqrt(Rho2 + X[2]*X[2]);
  double CT = X[2] / r, ST = sqrt(1.0-CT*CT); // cos(Theta), sin(Theta)
  double Phi = atan2(X[1], X[0]);
  double CP = cos(Phi), SP = sin(Phi);

  /**********************************************************************/
  /* get spherical components of E and H fields                         */
  /* (note: ES, HS = 'E, spherical' and 'H, spherical')                 */
  /**********************************************************************/
  cdouble hl, hlPrime;
  cdouble LY[3], rxLY[3];
  VectorSphericalHarmonic(L, K*r, &hl, &hlPrime);
  SphericalHankelH1(L, K*r, &hl, &hlPrime);

  cdouble ES[3], HS[3];
  if ( Type == SW_ELECTRIC ) // jackson equation 9.116 
   { 
     ES[0] = 0.0; 
     ES[1] = ZVAC*hl*LY[1];
     ES[2] = ZVAC*hl*LY[2];

     HS[0] = hlPrime * Y;
     HS[1] = hl*rxLY[1] / ZVAC;
     HS[2] = hl*rxLY[2] / ZVAC;
     
   }
  else  // Type == SW_MAGNETIC, use jackson equation 9.118
   {
     HS[0] = 0.0; 
     HS[1] = hl*LY[1];
     HS[2] = hl*LY[2];

     ES[0] = hlPrime * Y;
     ES[1] = hl*rxLY[1] / ZVAC;
     ES[2] = hl*rxLY[2] / ZVAC;

   };

  /**********************************************************************/
  /* convert spherical vector components to cartesian                   */
  /* (note: EC, HC = 'E, cartesian' and 'H, cartesian')                 */
  /**********************************************************************/
  cdouble EC[3], HC[3]; 
  memset(EC, 0, 3*sizeof(cdouble));
  memset(HC, 0, 3*sizeof(cdouble));

  M[0][0]=ST*CP;    M[0][1]=CT*CP;    M[0][2]=-SP;
  M[1][0]=ST*SP;    M[1][1]=CT*SP;    M[1][2]=CP;
  M[2][0]=CT;       M[2][1]=-ST;      M[2][2]=0.0;
  
  for(int mu=0; mu<3; mu++)
   for(int nu=0; nu<3; nu++)
    { EC[mu] += M[mu][nu]*ES[nu];
      HC[mu] += M[mu][nu]*HS[nu];
    };

  /**********************************************************************/
  /* pack cartesian components of E and H fields into EH output vector **/
  /**********************************************************************/
  EH[0]=EC[0];
  EH[1]=EC[1];
  EH[2]=EC[2];
  EH[3]=HC[0];
  EH[4]=HC[1];
  EH[5]=HC[2];


}

/**********************************************************************/
/* ********************************************************************/
/**********************************************************************/
void GetLYlm(int L, int M, double Theta, double Phi, cdouble LYlm[3])
{
}

/**********************************************************************/
/* sets LY[0..2]   = r, theta, phi components of L x Y_{lm}           */
/* sets rxLY[0..2] = \hat{r} \times LY                                */
/**********************************************************************/
void VectorSphericalHarmonic(int L, int M, double CosTheta, double SinTheta,
                             double Phi, cdouble LY[3], cdouble rxLY[3])
{
}


/**********************************************************************/
/* this and the following routine provide a simple implementation of  */
/* spherical hankel functions, not intended to be highly accurate or  */
/* efficient, based on these recurrences:                             */
/*  h_{l}(z)        = (2l-1)*h_{l-1}(z)/z - hl(l-2, z)                */
/*  h_{l}^\prime(z) = h_{h-1}(z) - (l+1)*h_l(z) /z                    */
/**********************************************************************/
cdouble Gethl(int l, cdouble z)
{ 
  if (l<0) return 0.0;

  cdouble ooz = 1.0/z;
  cdouble ExpFac = ooz*exp(i*z);

  switch(l)
   {
     case 0:  return  -II*ExpFac;
     case 1:  return -1.0*ExpFac*(1.0 + II*ooz);
     case 2:  return   II*ExpFac*(1.0 + ooz*(3.0*II - 3.0*ooz));
     case 3:  return      ExpFac*(1.0 + ooz*(6.0*II + ooz*(-15.0 -15.0*II*ooz)));
     default: return (2.0*((double)l)-1.0)*Gethl(l-1,z)/z - Gethl(l-2,z);
   };
}

void SphericalHankelH1(int l, cdouble z, cdouble *phl, cdouble *phlPrime)
{
  cdouble hl   = SphericalHankelH1(l, z);
  cdouble hlM1 = SphericalHankelH1(l-1, z);
  *phl = hl;
  *phlPrime = hlM1 - ( double(l) + 1 )*hl/z;
}
