/*
 * uupIntegralCubature.cc -- code for evaluating the (u,uprime) integral
 *                           using numerical cubature 
 * 
 * homer reid             -- 11/2008 -- 10/2010
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <gsl/gsl_sf.h>

#include "libTDRT.h"
#include "libhrutil.h"

// ln(2) - euler's constant 
#define LN2MGAMMA 0.115931515658412 
#define EULERGAMMA 0.57721566490153286061

/***************************************************************/
/* Use fixed-order numerical cubature to evaluate the integrals*/
/* I^{pq}_r over the line segments (Xs,Xe) and (Xsp,Xep).      */
/*                                                             */
/* I^{pq}_r                                                    */
/*  = l*lp*\int_0^1 du \int_0^1 dup u^p up^q KK_r(R(u,up))     */
/*                                                             */
/* where                                                       */
/*         R(u,up) = X(u) - Xp(up),                            */
/*            X(u) = Xs  +  u*(Xe-Xs)                          */
/*           Xp(u) = Xsp + up*(Xep-Xsp)                        */
/*           Xp(u) = Xsp + up*(Xep-Xsp)                        */
/*           KK_0  = K_0(Alpha*R)                              */
/*           KK_1  = Alpha*K_1(Alpha*R)/R                      */
/*           KK_2  = Alpha^2*K_2(Alpha*R)/(R^2)                */
/*                                                             */
/* and K_{0-2} are modified bessel functions of the second     */
/* kind.                                                       */
/*                                                             */
/* Order=(3,5,7,9) is the order of the numerical cubature      */
/* scheme to use.                                              */
/*                                                             */
/* If SSSIDR is non-null, it must point to a                   */
/* StaticSSIDataRecord for the given pair of line segments.    */
/*                                                             */
/* note: this is intended to be a fairly low-level,            */
/* unintelligent workhorse routine. in particular, this routine*/
/* does not attempt to assess how accurate the cubature is --  */
/* it just applies the cubature rule of the requested order.   */
/* the assumption is that the calling routine will have done   */
/* some work to figure out how close the segments are to one   */
/* another and will have chosen an appropriate cubature order  */
/* accordingly.                                                */
/* the same goes for the SSSIDR -- the assumption of this      */
/* routine is that the caller will have already figured out    */
/* whether or not we need to desingularize for this pair of    */
/* segments, and will have provided an appropriate SSSIDR if   */
/* so. if the SSSIDR is NULL even though the line segments are */
/* nearby one another, this routine will happily compute the   */
/* non-desingularized cubature and return the highly inaccurate*/
/* result.                                                     */
/***************************************************************/
void uupIntegralCubature(double *Xs, double *Xe, double *Xsp, double *Xep,
                         double Alpha, int Order, int NeedDerivatives,
                         StaticSSIDataRecord *SSSIDR, IPQRData *I)
{    
  int i, j, np, nqr, NumPts;
  const double *SCR;
  double u, up, w;
  double D[2], L[2], LP[2], XmXP[2];
  double l, lp, R2, R, Arg, Arg2o4, Factor;
  double K0, K1oR, K2oR2;
  double A2, LogArg, LogAlpha, Gamma1, Gamma2, Gamma3;
  int DeSingularize;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double MaxD;
  MaxD=VecDistance(Xs,Xsp);
  MaxD=fmax(MaxD,VecDistance(Xs,Xep));
  MaxD=fmax(MaxD,VecDistance(Xe,Xsp));
  MaxD=fmax(MaxD,VecDistance(Xe,Xep));
  // do not bother to compute the integrals if 
  // Alpha*MaxD < 45 (note e^{-45} ~ 2e-20)
  if ( Alpha*MaxD > 45.0 ) 
   { memset(I,0,sizeof(*I)); 
     return;
   };

  /***************************************************************/
  /* preliminary setup *******************************************/
  /***************************************************************/
  memset(I,0,sizeof(*I));
  DeSingularize = (SSSIDR==0) ? 0 : 1;

  /***************************************************************/
  /* choose a cubature scheme ************************************/
  /***************************************************************/
  switch(Order)
   { case 3:   SCR=SCR3; NumPts=5; break;
     case 5:   SCR=SCR5; NumPts=8; break;
     case 7:   SCR=SCR7; NumPts=12; break;
     case 9:   SCR=SCR9; NumPts=17; break;
     default:  ErrExit("unsupported cubature order %i",Order);
   };
  
  /***************************************************************/
  /* in the main loop we take X=Xs + u*L, Xp=Xsp+up*Lp (where    */
  /* L=(Xe-Xs), Lp=(Xep-Xsp)), so X-Xp = D + u*L - up*Lp, where  */
  /* D=Xs-Xsp.                                                   */
  /***************************************************************/
  for(i=0; i<2; i++)
   { L[i]=Xe[i]-Xs[i];
     LP[i]=Xep[i]-Xsp[i];
     D[i]=Xs[i]-Xsp[i];
   };
  l=sqrt(L[0]*L[0] + L[1]*L[1]);
  lp=sqrt(LP[0]*LP[0] + LP[1]*LP[1]);

  /***************************************************************/
  /* loop over all cubature points in the cubature scheme        */
  /***************************************************************/
  for(nqr=np=0; np<NumPts; np++)
   { 
     /* unpack next cubature point and weight */
     u=SCR[nqr++]; up=SCR[nqr++]; w=SCR[nqr++];
 
     /* compute R(u,up)=|X(u)-Xp(up)| */
     XmXP[0]=D[0] + u*L[0] - up*LP[0];
     XmXP[1]=D[1] + u*L[1] - up*LP[1];
     R2=XmXP[0]*XmXP[0] + XmXP[1]*XmXP[1];
     R=sqrt( XmXP[0]*XmXP[0] + XmXP[1]*XmXP[1] );
     Arg=Alpha*R;
    
     /* compute K0-1 factors. these are the modified bessel        */
     /* functions, possibly with the first several terms removed   */
     /* for desingularization purposes.                            */
     /* note: both desingularized functions are <1e-12 for         */
     /* Arg<=1e-4, so in this case we just set them to 0; even if  */
     /* they were worth computing, we would need to compute them   */
     /* another way, since the method we use here yields only      */ 
     /* numerical noise for Arg<=1e-4.                             */
     if ( Arg<=0.0 || Arg>50.0 )
      K0=K1oR=0.0;
     else if (DeSingularize)
      { if ( Arg<=1.0e-4 )
         K0=K1oR=0.0;
        else
         { 
           LogArg=log(Arg);

           K0 = gsl_sf_bessel_K0(Arg); 
           K0 += LogArg - LN2MGAMMA;

           K1oR = gsl_sf_bessel_K1(Arg)/R; 
           K1oR -= 1.0/(Arg*R) + 0.5*Alpha*(LogArg - LN2MGAMMA - 0.5);
         };
      }
     else
      { K0=gsl_sf_bessel_K0(Arg); 
        K1oR=gsl_sf_bessel_K1(Arg)/R;
      };
  
     /* compute contributions to I_i integrals from this cubature point */
     I->I_00_0 += w*K0;
     I->I_10_0 += w*u*K0;
     I->I_01_0 += w*up*K0;
     I->I_11_0 += w*u*up*K0;
     I->I_X_1  += w*u*up*XmXP[0]*K1oR;
     I->I_Y_1  += w*u*up*XmXP[1]*K1oR;

     if (NeedDerivatives)
      { 
        if ( Arg<=0.0 || Arg>50.0 )
         K2oR2=0.0;
        else if (DeSingularize)
         { if ( Arg<=0.05 )
            K2oR2=0.0;
           else
            { 
              K2oR2 = gsl_sf_bessel_Kn(2,Arg)/R2
                      -2.0/(Arg*Arg*R2) + 0.5/R2
                      +Alpha*Alpha*(LogArg-LN2MGAMMA-0.75)/8.0;
            };
         }
        else
         K2oR2 = gsl_sf_bessel_Kn(2,Arg)/R2;

        I->I_00_1 += w*K1oR;
        I->I_10_1 += w*u*K1oR;
        I->I_01_1 += w*up*K1oR;
        I->I_11_1 += w*u*up*K1oR;
        I->I_20_1 += w*u*u*K1oR;
        I->I_21_1 += w*u*u*up*K1oR;
        I->I_12_1 += w*u*up*up*K1oR;
        I->I_02_1 += w*up*up*K1oR;

        I->I_11_2 += w*u*up*K2oR2;
        I->I_21_2 += w*u*u*up*K2oR2;
        I->I_22_2 += w*u*u*up*up*K2oR2;
        I->I_12_2 += w*u*up*up*K2oR2;
        I->I_31_2 += w*u*u*u*up*K2oR2;
        I->I_13_2 += w*u*up*up*up*K2oR2;

     };

   }; // for(nqr=np=0; np<NumPts; np++)

  /***************************************************************/
  /*- if we desingularized, add back in the contributions of    -*/
  /*- the singular terms                                        -*/
  /***************************************************************/
  if ( DeSingularize && Arg<50.0 )
   { 
     A2=Alpha*Alpha;
     Gamma1=log(0.5*Alpha) + EULERGAMMA;
     Gamma2=0.5*Alpha*(Gamma1 - 0.5);
     Gamma3=0.125*A2*(Gamma1 - 0.75);
  
     I->I_00_0 -= SSSIDR->J_00_2 + Gamma1*SSSIDR->J_00_1;
     I->I_10_0 -= SSSIDR->J_10_2 + Gamma1*SSSIDR->J_10_1;
     I->I_01_0 -= SSSIDR->J_01_2 + Gamma1*SSSIDR->J_01_1;
     I->I_11_0 -= SSSIDR->J_11_2 + Gamma1*SSSIDR->J_11_1;
     I->I_X_1  += SSSIDR->J_X_3/Alpha + 0.5*Alpha*SSSIDR->J_X_2 + Gamma2*SSSIDR->J_X_1;
     I->I_Y_1  += SSSIDR->J_Y_3/Alpha + 0.5*Alpha*SSSIDR->J_Y_2 + Gamma2*SSSIDR->J_Y_1;

     if (NeedDerivatives)
      { 
        I->I_00_1 += SSSIDR->J_00_3/Alpha + 0.5*Alpha*SSSIDR->J_00_2 + Gamma2*SSSIDR->J_00_1;
        I->I_10_1 += SSSIDR->J_10_3/Alpha + 0.5*Alpha*SSSIDR->J_10_2 + Gamma2*SSSIDR->J_10_1;
        I->I_01_1 += SSSIDR->J_01_3/Alpha + 0.5*Alpha*SSSIDR->J_01_2 + Gamma2*SSSIDR->J_01_1;
        I->I_11_1 += SSSIDR->J_11_3/Alpha + 0.5*Alpha*SSSIDR->J_11_2 + Gamma2*SSSIDR->J_11_1;
        I->I_20_1 += SSSIDR->J_20_3/Alpha + 0.5*Alpha*SSSIDR->J_20_2 + Gamma2*SSSIDR->J_20_1;
        I->I_21_1 += SSSIDR->J_21_3/Alpha + 0.5*Alpha*SSSIDR->J_21_2 + Gamma2*SSSIDR->J_21_1;
        I->I_12_1 += SSSIDR->J_12_3/Alpha + 0.5*Alpha*SSSIDR->J_12_2 + Gamma2*SSSIDR->J_12_1;
        I->I_02_1 += SSSIDR->J_02_3/Alpha + 0.5*Alpha*SSSIDR->J_02_2 + Gamma2*SSSIDR->J_02_1;

        I->I_11_2 += 2.0*SSSIDR->J_11_4/A2 - 0.5*SSSIDR->J_11_3 - A2*SSSIDR->J_11_2/8.0 - Gamma3*SSSIDR->J_11_1;
        I->I_21_2 += 2.0*SSSIDR->J_21_4/A2 - 0.5*SSSIDR->J_21_3 - A2*SSSIDR->J_21_2/8.0 - Gamma3*SSSIDR->J_21_1;
        I->I_22_2 += 2.0*SSSIDR->J_22_4/A2 - 0.5*SSSIDR->J_22_3 - A2*SSSIDR->J_22_2/8.0 - Gamma3*SSSIDR->J_22_1;
        I->I_12_2 += 2.0*SSSIDR->J_12_4/A2 - 0.5*SSSIDR->J_12_3 - A2*SSSIDR->J_12_2/8.0 - Gamma3*SSSIDR->J_12_1;
        I->I_31_2 += 2.0*SSSIDR->J_31_4/A2 - 0.5*SSSIDR->J_31_3 - A2*SSSIDR->J_31_2/8.0 - Gamma3*SSSIDR->J_31_1;
        I->I_13_2 += 2.0*SSSIDR->J_13_4/A2 - 0.5*SSSIDR->J_13_3 - A2*SSSIDR->J_13_2/8.0 - Gamma3*SSSIDR->J_13_1;

      };

   };

  /***************************************************************/
  /* finally, put in prefactors **********************************/
  /***************************************************************/
  Factor=l*lp;
  I->I_00_0 *= Factor;
  I->I_10_0 *= Factor;
  I->I_01_0 *= Factor;
  I->I_11_0 *= Factor;

  Factor=l*lp*Alpha;
  I->I_00_1 *= Factor;
  I->I_10_1 *= Factor;
  I->I_01_1 *= Factor;
  I->I_11_1 *= Factor;
  I->I_20_1 *= Factor;
  I->I_21_1 *= Factor;
  I->I_12_1 *= Factor;
  I->I_02_1 *= Factor;
  I->I_X_1  *= Factor;
  I->I_Y_1  *= Factor;

  Factor=l*lp*Alpha*Alpha;
  I->I_11_2 *= Factor;
  I->I_21_2 *= Factor;
  I->I_22_2 *= Factor;
  I->I_12_2 *= Factor;
  I->I_31_2 *= Factor;
  I->I_13_2 *= Factor;
    
}
