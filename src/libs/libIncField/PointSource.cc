/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * PointSource.cc -- point source implementation of IncField
 *
 *           Note:   Unlike the other IncField implementations,
 *                   this one depends on libscuff because it
 *                   requires Ewald summation in the periodic
 *                   case.
 *
 * homer reid     -- 11/2009 -- 2/2012
 */

#include "libIncField.h"
#include <math.h>
#include <string.h>
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#include <stdio.h>
#include <libhrutil.h>
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#define II cdouble(0.0,1.0)

namespace scuff {

void GBarVDEwald(double *R, cdouble k, double *kBloch,
                 double **LBV, int LDim,
                 double E, bool ExcludeInnerCells, cdouble *GBarVD);

                }

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
PointSource::PointSource(const double pX0[3], const cdouble pP[3], int pType,
			 const char *Label)
{
  memcpy(X0, pX0, 3*sizeof(double));
  memcpy(P,  pP, 3*sizeof(cdouble));
  Type=pType;
  SetRegionLabel(Label);
}

PointSource::~PointSource()
{ 
  // no malloc'ed data to free
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void PointSource::SetX0(double pX0[3]) { memcpy(X0, pX0, 3*sizeof(double)); }
void PointSource::SetP(cdouble pP[3])  { memcpy(P,  pP, 3*sizeof(cdouble)); }
void PointSource::SetType(int pType)   { Type=pType; }

/**********************************************************************/
/* For this implementation of IncField we want to override the default*/
/* implementation of GetSourcePoint so that libscuff can use the      */
/* coordinates of the source point to determine automatically which   */
/* (if any) object contains the source point.                         */
/**********************************************************************/
bool PointSource::GetSourcePoint(double X[3]) const
{ 
  memcpy(X, X0, 3*sizeof(double)); 
  return true;
}

/**********************************************************************/
/* fields of a point source.                                          */
/*                                                                    */
/* NOTE: for the default case of an electric dipole, the quantity P   */
/* in the PointSource structure is assumed to be the dipole           */
/* moment divided by \epsilon_0, which means that P has units of      */
/* voltage*length^2.                                                  */
/*                                                                    */
/**********************************************************************/
void PointSource::GetFields(const double X[3], cdouble EH[6])
{
  if (LDim>0)
   { GetFields_Periodic(X, EH);
     return; 
   };

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
     ExpFac /= Eps;

     EH[0]=ExpFac*( Term1*P[0] + Term2*RHat[0] );
     EH[1]=ExpFac*( Term1*P[1] + Term2*RHat[1] );
     EH[2]=ExpFac*( Term1*P[2] + Term2*RHat[2] );

     EH[3]=ExpFac*Term3*RCrossP[0] / Z;
     EH[4]=ExpFac*Term3*RCrossP[1] / Z;
     EH[5]=ExpFac*Term3*RCrossP[2] / Z;
   }
  else // ( Type == LIF_TYPE_PSMC )
   { 
     ExpFac /= Mu;

     EH[0]=-1.0*ExpFac*Term3*RCrossP[0] / Z;
     EH[1]=-1.0*ExpFac*Term3*RCrossP[1] / Z;
     EH[2]=-1.0*ExpFac*Term3*RCrossP[2] / Z;

     EH[3]=ExpFac*( Term1*P[0] + Term2*RHat[0] ) / (Z*Z);
     EH[4]=ExpFac*( Term1*P[1] + Term2*RHat[1] ) / (Z*Z);
     EH[5]=ExpFac*( Term1*P[2] + Term2*RHat[2] ) / (Z*Z);
   };

}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void PointSource::GetFields_Periodic(const double X[3], cdouble EH[6])
{
  cdouble k    = sqrt(Eps*Mu) * Omega;
  cdouble k2   = k*k;
  cdouble ZRel = sqrt(Mu/Eps);

  double R[3];
  R[0]= X[0]-X0[0];
  R[1]= X[1]-X0[1];
  R[2]= X[2]-X0[2];

  /***************************************************************/
  /* get the scalar green's function, its first derivatives, and */
  /* its mixed second partials by Ewald summation.               */
  /***************************************************************/
  cdouble GArray[8];
  scuff::GBarVDEwald(R, k, kBloch, LBV, LDim, -1.0, false, GArray);

  cdouble G = GArray[0];

  cdouble dG[3];
  dG[0] = GArray[1];
  dG[1] = GArray[2];
  dG[2] = GArray[3];

  cdouble ddG[3][3];
  ddG[0][1] = ddG[1][0] = GArray[4];
  ddG[0][2] = ddG[2][0] = GArray[5];
  ddG[1][2] = ddG[2][1] = GArray[6];

  /***************************************************************/
  /* The Ewald routine only computes the mixed second partials,  */
  /* so we do finite differencing to get unmixed second partials.*/
  /***************************************************************/
  for(int i=0; i<3; i++)
   { 
     double DeltaR = (R[i]==0.0 ? 1.0e-4 : 1.0e-4*fabs(R[i]) );
     double RTweaked[3];

     memcpy(RTweaked, R, 3*sizeof(double));
     RTweaked[i]+=DeltaR;
     cdouble GPArray[8];
     scuff::GBarVDEwald(RTweaked, k, kBloch, LBV, LDim, -1.0, false, GPArray);
     cdouble GP=GPArray[0];

     memcpy(RTweaked, R, 3*sizeof(double));
     RTweaked[i]-=DeltaR;
     cdouble GMArray[8];
     scuff::GBarVDEwald(RTweaked, k, kBloch, LBV, LDim, -1.0, false, GMArray);
     cdouble GM=GMArray[0];

     ddG[i][i] = (GP + GM - 2.0*G) / (DeltaR*DeltaR);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
printf("ddG_{%i%i} (method 1) = %s \n",i,i,CD2S(ddG[i][i])); 
cdouble ddG2 = (GPArray[1+i] - GMArray[1+i]) / (2.0*DeltaR);
printf("ddG_{%i%i} (method 2) = %s \n",i,i,CD2S(ddG2));
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   };

  /***************************************************************/
  /* now assemble the derivatives of G0 appropriately to form the*/
  /* dyadic GFs and read off the E and H fields due to the source*/
  /* note the scuff convention that dipole moment is measured    */
  /* in units of volts/um^2 instead of coulomb*um; what this     */
  /* means is that the numerical value of the dipole moment you  */
  /* specify to scuff is the dipole moment in coulombs*microns   */
  /* divided by 377 (the impedance of free space).               */
  /***************************************************************/
  //cdouble PreFac1 = Omega*k*ZVAC*ZRel;
  //cdouble PreFac2 = -II*Omega;
  cdouble PreFac1 = Omega*k*ZRel;
  cdouble PreFac2 = II*Omega/ZVAC;

  EH[0*3 + 0 ]
   = PreFac1 * (G*P[0] + (ddG[0][0]*P[0]+ddG[0][1]*P[1]+ddG[0][2]*P[2])/k2 );
  EH[0*3 + 1 ] 
   = PreFac1 * (G*P[1] + (ddG[1][0]*P[0]+ddG[1][1]*P[1]+ddG[1][2]*P[2])/k2 );
  EH[0*3 + 2 ] 
   = PreFac1 * (G*P[2] + (ddG[2][0]*P[0]+ddG[2][1]*P[1]+ddG[2][2]*P[2])/k2 );

  EH[1*3 + 0] = PreFac2 * (P[1]*dG[2] - P[2]*dG[1]);
  EH[1*3 + 1] = PreFac2 * (P[2]*dG[0] - P[0]*dG[2]);
  EH[1*3 + 2] = PreFac2 * (P[0]*dG[1] - P[1]*dG[0]);

}
