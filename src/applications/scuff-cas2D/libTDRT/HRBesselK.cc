/* 
 * HRBesselK.cc -- an accelerated routine for computing BesselK
 *
 * homer reid   -- 6/2011
 */

#include <math.h>
#include <string.h>
#include <libAmosBessel.h>
#include <libMDInterp.h>

/*--------------------------------------------------------------*/
// constants for the small-z series expansions, defined as follows:
// BesselK[0,z] 
//   = A00 + A02*z^2 + A04*z^2 + log(z/2)*(B00 + B02*z^2 + B04*z^4)
// BesselK[1,z]/z - 1/z^2  
//   = A10 + A12*z^2 + A14*z^2 + log(z/2)*(B10 + B12*z^2 + B14*z^4)
// BesselK[2,z]/z^2 - (2/z^4 - 1/2z^2)
//   = A20 + A22*z^2 + A24*z^2 + log(z/2)*(B20 + B22*z^2 + B24*z^4)
/*--------------------------------------------------------------*/
#define EULERGAMMA 0.57721566490153286061
#define A00 -EULERGAMMA
#define A02 (1.0-EULERGAMMA)/4.0
#define A04 (3.0-2.0*EULERGAMMA)/128.0
#define B00 -1.0
#define B02 -1.0/4.0
#define B04 -1.0/64.0
#define A10 (2.0*EULERGAMMA-1.0)/4.0
#define A12 (2.0*EULERGAMMA-2.5)/32.0
#define A14 (2.0*EULERGAMMA-10.0/3.0)/768.0
#define B10 1.0/2.0
#define B12 1.0/16.0
#define B14 1.0/384.0
#define A20 -(2.0*EULERGAMMA-3.0/2.0)/16.0
#define A22 -(2.0*EULERGAMMA-17.0/6.0)/192.0
#define A24 -(2.0*EULERGAMMA-43.0/12.0)/6144.0
#define B20 -1.0/8.0
#define B22 -1.0/96.0
#define B24 -1.0/3072.0

/*--------------------------------------------------------------*/
// constants for the large-z series expansions, defined as follows:
// BesselK[0,z] 
//   = sqrt(pi/2z)* Exp[-z] * ( C00 + C01/z + C02/z^2 + \cdots )
// BesselK[1,z]
//   = sqrt(pi/2z)* Exp[-z] * ( C10 + C11/z + C12/z^2 + \cdots )
// BesselK[2,z]
//   = sqrt(pi/2z)* Exp[-z] * ( C20 + C21/z + C22/z^2 + \cdots )
/*--------------------------------------------------------------*/
#define C00 1.0
#define C01 -1.0/8.0
#define C02 9.0/128.0
#define C03 -75.0/1024.0
#define C10 1.0
#define C11 3.0/8.0
#define C12 -15.0/128.0
#define C13 105.0/1024.0
#define C20 1.0
#define C21 15.0/8.0
#define C22 105.0/128.0
#define C23 -315.0/1024.0

#define RTPI2 1.25331413731550

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define LOWZTHRESH  0.2
#define HIGHZTHRESH 25.0
#define NINTERP     5000
static Interp1D *BesselKInterp=0;

/***************************************************************/
/* accelerated routine for computing BesselK                   */
/*                                                             */
/* return values:                                              */
/*                                                             */
/*  KArray [0] = BesselK[ 0, z ]                               */
/*  KArray [1] = BesselK[ 1, z ] / z                           */
/*  KArray [2] = BesselK[ 2, z ] / z^2                         */
/*                                                             */
/***************************************************************/
void HRBesselK(double z, int NeedK2, double *KArray)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (BesselKInterp==0)
   { cdouble ZKArray[3];
     if (NeedK2)
      { AmosBessel('K', z, 0, 3, 0, ZKArray);
        KArray[0]=real(ZKArray[0]);
        KArray[1]=real(ZKArray[1]) / z;
        KArray[2]=real(ZKArray[2]) / (z*z);
      }
     else
      { AmosBessel('K', z, 0, 2, 0, ZKArray);
        KArray[0]=real(ZKArray[0]);
        KArray[1]=real(ZKArray[1]) / z;
      };
     return;
   };

  /*--------------------------------------------------------------*/
  /*- low-z expansion --------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( z<=LOWZTHRESH )
   { double lzo2=log(0.5*z);
     double z2=z*z;
     double ooz2=1.0/z2;

     KArray[0] =        A00 + lzo2*B00 + z2*((A02 + lzo2*B02) + z2*(A04 + lzo2*B04));
     KArray[1] = ooz2 + A10 + lzo2*B10 + z2*((A12 + lzo2*B12) + z2*(A14 + lzo2*B14));

     if (NeedK2)
      KArray[2]=(2.0*ooz2-0.5)*ooz2 
                 + A20 + lzo2*B20 + z2*((A22 + lzo2*B22) + z2*(A24 + lzo2*B24));
   }
  /*--------------------------------------------------------------*/
  /*- high-z expansion  ------------------------------------------*/
  /*--------------------------------------------------------------*/
  else if ( z>=HIGHZTHRESH )
   {
     double ExpFac=RTPI2 * exp(-z) / sqrt(z) ;
     double ooz=1.0/z;
     KArray[0]=ExpFac*(C00 + ooz*(C01 + ooz*(C02 + ooz*C03)));
     KArray[1]=ExpFac*(C10 + ooz*(C11 + ooz*(C12 + ooz*C13)))/z;
     if (NeedK2)
      KArray[2]=ExpFac*(C20 + ooz*(C21 + ooz*(C22 + ooz*C23)))/(z*z);
   }
  /*--------------------------------------------------------------*/
  /*- in the intermediate region we use an interpolating table.   */
  /*- note that the values cached in the table are not exactly    */
  /*- the bessel functions, but rather the following quantities:  */
  /*-  K[0]                                                       */
  /*-  K[1]/z   - 1/z^2                                           */
  /*-  K[2]/z^2 - 2/z^4 + 0.5/z^2                                 */
  /*--------------------------------------------------------------*/
  else
   { 
     double ooz2=1.0/(z*z);
     if (NeedK2)
      { BesselKInterp->Evaluate(z, KArray);
        KArray[1]+=ooz2;
        KArray[2]+=ooz2*(-0.5 + 2.0*ooz2);
      }
     else
      { double KKArray[3];
        BesselKInterp->Evaluate(z, KKArray);
        KArray[0]=KKArray[0];
        KArray[1]=KKArray[1] + ooz2;
      };
   };

} 

/***************************************************************/
/***************************************************************/
/***************************************************************/ 
void MyPhi1D(double z, void *UserData, double *PhiVD)
{ 

  cdouble ZKArray[4];
  double KArray[4];
  double z2=z*z;
  double z3=z2*z;
  double z4=z3*z;
  double z5=z4*z;

  AmosBessel('K', z, 0, 4, 0, ZKArray);
  
  KArray[0]=real(ZKArray[0]);
  KArray[1]=real(ZKArray[1]);
  KArray[2]=real(ZKArray[2]);
  KArray[3]=real(ZKArray[3]);

  PhiVD[0] = KArray[0];
  PhiVD[1] = -KArray[1];

  PhiVD[2] = KArray[1]/z - 1.0/z2;
  PhiVD[3] = -KArray[2]/z + 2.0/z3;

  PhiVD[4] = KArray[2]/z2 - 2.0/z4 + 0.5/z2;
  PhiVD[5] = -KArray[3]/z2 + 8.0/z5 - 1.0/z3;
} 

int InitHRBesselK()
{
  BesselKInterp=new Interp1D(LOWZTHRESH, HIGHZTHRESH, NINTERP, 3, MyPhi1D, 0);
}
