/*
 * LFunctionsMultipole.cc -- multipole approximations to L-functions
 *                           between RWG basis functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>

#include "libscuff.h"

#define II cdouble (0,1)

/***************************************************************/
/***************************************************************/
/***************************************************************/
static double LeviCivita[3][3][3]=
{ { { 0.0,  0.0,  0.0 }, {  0.0, 0.0, +1.0 }, {  0.0, -1.0, +0.0 }  },
  { { 0.0,  0.0, -1.0 }, {  0.0, 0.0,  0.0 }, { +1.0, +0.0, +0.0 }  },
  { { 0.0, +1.0,  0.0 }, { -1.0, 0.0,  0.0 }, {  0.0,  0.0,  0.0 }  }
};

/***************************************************************/
/* Use multipole expansions to compute L-functions and their   */
/* derivatives. (Angular derivatives not currently supported.) */
/* Set Order = 1 or 2 to retain the first nonvanishing 1 or 2  */
/*             terms in the multipole expansions.              */
/***************************************************************/
void GetLFunctions_Multipole(int Order,
                             RWGObject *O1, int ne1, RWGObject *O2, int ne2,
                             double Wavevector, int RealFreq, int NeedCross,
                             cdouble L[3], cdouble *GradL)
{ 
  RWGEdge *E1, *E2;
  double *QPAlpha, *V1Alpha, *V2Alpha, *QMAlpha;
  double *QPBeta, *V1Beta, *V2Beta, *QMBeta;
  double APAlpha[3], AMAlpha[3], BPAlpha[3], BMAlpha[3];
  double APBeta[3], AMBeta[3], BPBeta[3], BMBeta[3];
  double LAlpha, LBeta;

  double MMuAlpha[3], MMuBeta[3];
  double MMuNuAlpha[3][3], MMuNuBeta[3][3];
  double SMuNuAlpha[3][3], SMuNuBeta[3][3];
  double MMuNuRhoAlpha[3][3][3], MMuNuRhoBeta[3][3][3];
  
  int Mu, Nu, Rho, Sigma, Tau;

  double R0[3], r0, r03, r05, r07;
  double Term, dL1, dL2, LC;
  cdouble ikr, ExpFac, Phi0, Psi0, Zeta0, Upsilon0, Theta0;

  /***************************************************************/
  /* unpack data on edges and vertices ***************************/
  /***************************************************************/
  E1=O1->Edges[ne1];
  QPAlpha=O1->Vertices + 3*(E1->iQP);
  V1Alpha=O1->Vertices + 3*(E1->iV1);
  V2Alpha=O1->Vertices + 3*(E1->iV2);
  QMAlpha=O1->Vertices + 3*(E1->iQM);
  VecSub(V1Alpha, QPAlpha, APAlpha);
  VecSub(V1Alpha, QMAlpha, AMAlpha);
  VecSub(V2Alpha, QPAlpha, BPAlpha);
  VecSub(V2Alpha, QMAlpha, BMAlpha);
  LAlpha=E1->Length;

  E2=O2->Edges[ne2];
  QPBeta=O2->Vertices + 3*(E2->iQP);
  V1Beta=O2->Vertices + 3*(E2->iV1);
  V2Beta=O2->Vertices + 3*(E2->iV2);
  QMBeta=O2->Vertices + 3*(E2->iQM);
  VecSub(V1Beta, QPBeta, APBeta);
  VecSub(V1Beta, QMBeta, AMBeta);
  VecSub(V2Beta, QPBeta, BPBeta);
  VecSub(V2Beta, QMBeta, BMBeta);
  LBeta=E2->Length;

  /***************************************************************/
  /* compute multipole moments needed for 1st nonvanishing terms */
  /***************************************************************/
  for(Mu=0; Mu<3; Mu++)
   { 
     MMuAlpha[Mu] = LAlpha*(QMAlpha[Mu]-QPAlpha[Mu])/3.0;
     MMuBeta[Mu] = LBeta*(QMBeta[Mu]-QPBeta[Mu])/3.0;

     for(Nu=0; Nu<3; Nu++)
      { 
        MMuNuAlpha[Mu][Nu] = LAlpha*(   AMAlpha[Mu]*BMAlpha[Nu]
                                      + BMAlpha[Mu]*AMAlpha[Nu]
                                      - APAlpha[Mu]*BPAlpha[Nu]
                                      - BPAlpha[Mu]*APAlpha[Nu] ) / 24.0;

        MMuNuBeta[Mu][Nu] = LBeta*(     AMBeta[Mu]*BMBeta[Nu]
                                      + BMBeta[Mu]*AMBeta[Nu]
                                      - APBeta[Mu]*BPBeta[Nu]
                                      - BPBeta[Mu]*APBeta[Nu] ) / 24.0;

      };
    };
  
  /***************************************************************/
  /* compute the scalar kernel functions                         */
  /***************************************************************/
  VecSub(E1->Centroid, E2->Centroid, R0);
  r0=VecNorm(R0);
  r03=r0*r0*r0;
  r05=r03*r0*r0;
  r07=r05*r0*r0;
  if (RealFreq)
   { ikr=II*Wavevector*r0;
     ExpFac=expi(Wavevector*r0) / (4.0*M_PI);
   }
  else
   { ikr=-1.0*Wavevector*r0;
     ExpFac=exp(-Wavevector*r0) / (4.0*M_PI);
   };
  Phi0=ExpFac/r0;
  Psi0=(-1.0 + ikr)*ExpFac/r03;
  Zeta0=(3.0 + ikr*(-3.0 + ikr))*ExpFac/r05;
  Upsilon0=(-15.0 + ikr*(15.0 + ikr*(-6.0 + ikr)))*ExpFac/r07;

  /***************************************************************/
  /* assemble first nonvanishing terms in multipole expansions   */
  /***************************************************************/

  /*==================================================*/
  /*= L_\bullet ======================================*/
  /*==================================================*/
  Term=MMuAlpha[0]*MMuBeta[0]+MMuAlpha[1]*MMuBeta[1]+MMuAlpha[2]*MMuBeta[2];
  L[0] = Phi0*Term;

  if (GradL)
   { GradL[3*0 + 0]=Term*R0[0]*Psi0;
     GradL[3*1 + 0]=Term*R0[1]*Psi0;
     GradL[3*2 + 0]=Term*R0[2]*Psi0;
   };

  /*==================================================*/
  /*= L_\nabla =======================================*/
  /*==================================================*/
  dL1=dL2=0.0;
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    { Term = ( MMuAlpha[Mu]*MMuBeta[Nu] + MMuAlpha[Nu]*MMuBeta[Mu] );
      if (Mu==Nu) dL1+=Term;
      dL2 += R0[Mu]*R0[Nu]*Term;
    };
  L[1] = -0.5*(Psi0*dL1 + Zeta0*dL2);

  if (GradL)
   { 
     GradL[3*0 + 1] = -0.5*R0[0]*(Zeta0*dL1 + Upsilon0*dL2);
     GradL[3*1 + 1] = -0.5*R0[1]*(Zeta0*dL1 + Upsilon0*dL2);
     GradL[3*2 + 1] = -0.5*R0[2]*(Zeta0*dL1 + Upsilon0*dL2);

     for(Tau=0;Tau<3; Tau++)
      { Term=  R0[0]*(MMuAlpha[Tau]*MMuBeta[0] + MMuAlpha[0]*MMuBeta[Tau])
              +R0[1]*(MMuAlpha[Tau]*MMuBeta[1] + MMuAlpha[1]*MMuBeta[Tau])
              +R0[2]*(MMuAlpha[Tau]*MMuBeta[2] + MMuAlpha[2]*MMuBeta[Tau]);
        GradL[3*Tau + 1] += -1.0*Zeta0*Term;
      };
   };

  /*==================================================*/
  /*= L_\times =======================================*/
  /*==================================================*/
  if (NeedCross)
   {
     Term=0.0;
     for(Mu=0; Mu<3; Mu++)
      for(Nu=0; Nu<3; Nu++)
       for(Rho=0; Rho<3; Rho++)
        { if ((LC=LeviCivita[Mu][Nu][Rho])==0.0) 
           continue;
          Term+=LC*(   R0[Mu]*MMuAlpha[Nu]*MMuBeta[Rho]
                     + MMuNuAlpha[Mu][Nu]*MMuBeta[Rho]
                     - MMuAlpha[Nu]*MMuNuBeta[Mu][Rho]);
        };
     L[2] = Psi0*Term;
  
     if (GradL)
      { 
        GradL[3*0 + 2] = R0[0]*Zeta0*Term;
        GradL[3*1 + 2] = R0[1]*Zeta0*Term;
        GradL[3*2 + 2] = R0[2]*Zeta0*Term;
  
       /* extra term that comes from differentiating R_0^\mu in */
       /* the integrand                                         */
       for(Tau=0; Tau<3; Tau++)
        { 
          dL1=dL2=0.0;
          for(Nu=0; Nu<3; Nu++)
           for(Rho=0; Rho<3; Rho++)
            { if ((LC=LeviCivita[Tau][Nu][Rho])==0.0) 
               continue;
              dL1+=LC*MMuAlpha[Nu]*MMuBeta[Rho];
  
              for(Sigma=0; Sigma<3; Sigma++)
               dL2+=R0[Sigma]*(   MMuNuAlpha[Sigma][Nu]*MMuBeta[Rho]
                                - MMuAlpha[Nu]*MMuNuBeta[Sigma][Rho]);
            };
          GradL[3*Tau+2] += dL1*Psi0 + dL2*Zeta0;
        };
      };
 
   }; // if (NeedCross)

  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  if (Order==1) 
   return;
  
  /****************************************************************/
  /* compute multipole moments needed for 2nd nonvanishing term   */
  /****************************************************************/
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    {
      SMuNuAlpha[Mu][Nu] = LAlpha*(   APAlpha[Mu]*APAlpha[Nu]
                                    + BPAlpha[Mu]*BPAlpha[Nu]
                                    - AMAlpha[Mu]*AMAlpha[Nu]
                                    - BMAlpha[Mu]*BMAlpha[Nu] ) / 12.0;

      SMuNuBeta[Mu][Nu] =   LBeta*(   APBeta[Mu]*APBeta[Nu]
                                    + BPBeta[Mu]*BPBeta[Nu]
                                    - AMBeta[Mu]*AMBeta[Nu]
                                    - BMBeta[Mu]*BMBeta[Nu] ) / 12.0;
    };

  if (NeedCross)
   {
     for(Mu=0; Mu<3; Mu++)
      for(Nu=0; Nu<3; Nu++)
       for(Rho=0; Rho<3; Rho++)
        { 
          MMuNuRhoAlpha[Mu][Nu][Rho] 
            = LAlpha*(  2.0*APAlpha[Mu]*APAlpha[Nu]*APAlpha[Rho]
                       +4.0*APAlpha[Mu]*APAlpha[Nu]*BPAlpha[Rho]
                           -APAlpha[Mu]*BPAlpha[Nu]*APAlpha[Rho]
                           -APAlpha[Mu]*BPAlpha[Nu]*BPAlpha[Rho]
                           -BPAlpha[Mu]*APAlpha[Nu]*APAlpha[Rho]
                           -BPAlpha[Mu]*APAlpha[Nu]*BPAlpha[Rho]
                       +4.0*BPAlpha[Mu]*BPAlpha[Nu]*APAlpha[Rho]
                       +2.0*BPAlpha[Mu]*BPAlpha[Nu]*BPAlpha[Rho]
                       -2.0*AMAlpha[Mu]*AMAlpha[Nu]*AMAlpha[Rho]
                       -4.0*AMAlpha[Mu]*AMAlpha[Nu]*BMAlpha[Rho]
                           +AMAlpha[Mu]*BMAlpha[Nu]*AMAlpha[Rho]
                           +AMAlpha[Mu]*BMAlpha[Nu]*BMAlpha[Rho]
                           +BMAlpha[Mu]*AMAlpha[Nu]*AMAlpha[Rho]
                           +BMAlpha[Mu]*AMAlpha[Nu]*BMAlpha[Rho]
                       -4.0*BMAlpha[Mu]*BMAlpha[Nu]*AMAlpha[Rho]
                       -2.0*BMAlpha[Mu]*BMAlpha[Nu]*BMAlpha[Rho] 
                     )/240.0;
 
          MMuNuRhoBeta[Mu][Nu][Rho] 
            = LBeta*(   2.0*APBeta[Mu]*APBeta[Nu]*APBeta[Rho]
                       +4.0*APBeta[Mu]*APBeta[Nu]*BPBeta[Rho]
                           -APBeta[Mu]*BPBeta[Nu]*APBeta[Rho]
                           -APBeta[Mu]*BPBeta[Nu]*BPBeta[Rho]
                           -BPBeta[Mu]*APBeta[Nu]*APBeta[Rho]
                           -BPBeta[Mu]*APBeta[Nu]*BPBeta[Rho]
                       +4.0*BPBeta[Mu]*BPBeta[Nu]*APBeta[Rho]
                       +2.0*BPBeta[Mu]*BPBeta[Nu]*BPBeta[Rho]
                       -2.0*AMBeta[Mu]*AMBeta[Nu]*AMBeta[Rho]
                       -4.0*AMBeta[Mu]*AMBeta[Nu]*BMBeta[Rho]
                           +AMBeta[Mu]*BMBeta[Nu]*AMBeta[Rho]
                           +AMBeta[Mu]*BMBeta[Nu]*BMBeta[Rho]
                           +BMBeta[Mu]*AMBeta[Nu]*AMBeta[Rho]
                           +BMBeta[Mu]*AMBeta[Nu]*BMBeta[Rho]
                       -4.0*BMBeta[Mu]*BMBeta[Nu]*AMBeta[Rho]
                       -2.0*BMBeta[Mu]*BMBeta[Nu]*BMBeta[Rho] 
                    )/240.0;

         }; /* for (Rho=0; Rho<3; Rho++) */

   }; // if (NeedCross)

  /***************************************************************/
  /* add second nonvanishing terms in multipole expansions       */
  /***************************************************************/
  Theta0=(105.0 + ikr*(-105.0 + ikr*(45.0 + ikr*(-10.0 + ikr))))*ExpFac/(r07*r0*r0);

  /*==================================================*/
  /*= L_\bullet=======================================*/
  /*==================================================*/
  Term=0.0;
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    Term+=R0[Mu]*(MMuNuAlpha[Mu][Nu]*MMuBeta[Nu] - MMuAlpha[Nu]*MMuNuBeta[Mu][Nu]);
  L[0] += Psi0*Term;

  if (GradL)
   { 
     for(Tau=0; Tau<3; Tau++)
      { GradL[3*Tau + 0] += R0[Tau]*Zeta0*Term;
        GradL[3*Tau + 0] += Psi0*(  MMuNuAlpha[Tau][0]*MMuBeta[0]-MMuAlpha[0]*MMuNuBeta[Tau][0]
                                   +MMuNuAlpha[Tau][1]*MMuBeta[1]-MMuAlpha[1]*MMuNuBeta[Tau][1]
                                   +MMuNuAlpha[Tau][2]*MMuBeta[2]-MMuAlpha[2]*MMuNuBeta[Tau][2]);
      };
   };

  /*==================================================*/
  /*= L_\nabla =======================================*/
  /*==================================================*/
  dL1=dL2=0.0;
  for(Mu=0; Mu<3; Mu++)
   for(Nu=0; Nu<3; Nu++)
    for(Rho=0; Rho<3; Rho++)
     {
       Term=(  SMuNuAlpha[Mu][Nu ]*MMuBeta[Rho]
              +SMuNuAlpha[Mu][Rho]*MMuBeta[Nu ]
              +SMuNuAlpha[Nu][Rho]*MMuBeta[Mu ]
              -MMuAlpha[Mu]*SMuNuBeta[Nu][Rho ]
              -MMuAlpha[Nu]*SMuNuBeta[Mu][Rho ]
              -MMuAlpha[Rho]*SMuNuBeta[Mu][Nu ] ) / 6.0;

       dL1 += Term*(   Mu==Nu  ? R0[Rho] : 0.0 
                     + Mu==Rho ? R0[Nu ] : 0.0
                     + Nu==Rho ? R0[Mu ] : 0.0 );

       dL2 += Term*R0[Mu]*R0[Nu]*R0[Rho];
     };
  L[1] += Zeta0*dL1 + Upsilon0*dL2;

  if (GradL)
   { 
     for(Tau=0; Tau<3; Tau++)
      { 
        GradL[3*Tau + 1] += R0[Tau]*Upsilon0*dL1 + R0[Tau]*Theta0*dL2;
        
        for(Mu=0; Mu<3; Mu++)
         for(Nu=0; Nu<3; Nu++)
          { Term= ( SMuNuAlpha[Mu][Nu ]*MMuBeta[Tau]
                   +SMuNuAlpha[Mu][Tau]*MMuBeta[Nu ]
                   +SMuNuAlpha[Nu][Tau]*MMuBeta[Mu ]
                   -MMuAlpha[Mu]*SMuNuBeta[Nu][Tau ]
                   -MMuAlpha[Nu]*SMuNuBeta[Mu][Tau ]
                   -MMuAlpha[Tau]*SMuNuBeta[Mu][Nu ] ) / 2.0;

            GradL[3*Tau + 1] += Term*(  Mu==Nu ? Zeta0 : 0.0
                                       +R0[Mu]*R0[Nu]*Upsilon0 );
          };
      };
   };

  /*==================================================*/
  /*= L_\times =======================================*/
  /*==================================================*/
  if (NeedCross)
   {
     Term=0.0;
     for(Mu=0; Mu<3; Mu++)
      for(Nu=0; Nu<3; Nu++)
       for(Rho=0; Rho<3; Rho++)
        { 
          if ( (LC=LeviCivita[Mu][Nu][Rho]) == 0.0 ) 
           continue;
   
          for(Sigma=0; Sigma<3; Sigma++)
           { Term+=LC*( R0[Mu]*R0[Sigma]*( MMuNuAlpha[Sigma][Nu]*MMuBeta[Rho]
                                          -MMuAlpha[Nu]*MMuNuBeta[Sigma][Rho] 
                                         )
                       +R0[Sigma]*( MMuNuRhoAlpha[Mu][Sigma][Nu]*MMuBeta[Rho]
                                   -MMuNuAlpha[Mu][Nu]*MMuNuBeta[Sigma][Rho]
                                   -MMuNuAlpha[Sigma][Nu]*MMuNuBeta[Mu][Rho]
                                   +MMuAlpha[Nu]*MMuNuRhoBeta[Mu][Sigma][Rho]
                                  )
                      );
           };
        };
     L[2] += Zeta0*Term;
     if (GradL)
      { 
      };

   }; // if (NeedCross)

}
