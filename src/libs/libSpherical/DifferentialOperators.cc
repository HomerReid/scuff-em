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
 * DifferentialOperators.cc -- libSpherical routines for finite-difference
 *                          -- evaluation of vector-calculus operators
 *                          -- in spherical coordinates
 *
 * homer reid               -- 4/2005 -- 2/2010
 */
#include "libSpherical.h"
 
/***************************************************************/
/* compute the curl of a vector-valued function ****************/
/***************************************************************/
void GetCurl(double r, double Theta, double Phi,
             void (*F)(double r, double Theta, double Phi, cdouble *V),
             cdouble *CurlF, cdouble *CurlF2)
{
 
  /***************************************************************/
  /* METHOD 1: use the expression for the curl in spherical      */
  /*           coordinates                                       */
  /***************************************************************/
#define METHOD1
#ifdef METHOD1
{
  cdouble F0[3], FP[3], FM[3], dFdr[3], dFdTheta[3], dFdPhi[3];
  double Delta;
  double CT, ST;
  int Mu;

  CT=cos(Theta);
  ST=sin(Theta);

  F(r, Theta, Phi, F0);
 
  Delta = (r==0.0) ? 1.0e-4 : 1.0e-4*fabs(r);
  F(r+Delta, Theta, Phi, FP);
  F(r-Delta, Theta, Phi, FM); 
  for(Mu=0; Mu<3; Mu++)
   dFdr[Mu] = (FP[Mu]-FM[Mu]) / (2.0*Delta);
 
  Delta = (Theta==0.0) ? 1.0e-4 : 1.0e-4*fabs(Theta);
  F(r, Theta+Delta, Phi, FP);
  F(r, Theta-Delta, Phi, FM);
  for(Mu=0; Mu<3; Mu++)
   dFdTheta[Mu] = (FP[Mu]-FM[Mu]) / (2.0*Delta);
 
  Delta = (Phi==0.0) ? 1.0e-4 : 1.0e-4*fabs(Phi);
  F(r, Theta, Phi+Delta, FP);
  F(r, Theta, Phi-Delta, FM);
  for(Mu=0; Mu<3; Mu++)
   dFdPhi[Mu] = (FP[Mu]-FM[Mu]) / (2.0*Delta);

  /*
  F(r+DeltaR, Theta, Phi, dFdTheta); 
  F(r, Theta+DeltaT, Phi, dFdTheta); 
  F(r, Theta, Phi+DeltaP, dFdPhi); 
  for(Mu=0; Mu<3; Mu++)
   { dFdr[Mu]     = (dFdr[Mu]     - F0[Mu]) / DeltaR;
     dFdTheta[Mu] = (dFdTheta[Mu] - F0[Mu]) / DeltaT;
     dFdPhi[Mu]   = (dFdPhi[Mu]   - F0[Mu]) / DeltaP;
   };
  */

  CurlF[0]=(ST*dFdTheta[2] + CT*F0[2] - dFdPhi[1]) / (r*ST);
  CurlF[1]=dFdPhi[0]/(r*ST) - F0[2]/r - dFdr[2];
  CurlF[2]=F0[1]/r + dFdr[1] - dFdTheta[0]/r;

 };
#endif 

  /***************************************************************/
  /* METHOD 2: convert back and forth to cartesian components    */
  /***************************************************************/
#define METHOD2
#ifdef METHOD2
{
  double X[3], rr, tt, pp;
  cdouble F0[3], F0C[3], dFdX[3], dFdY[3], dFdZ[3], CurlFC[3];
  double Delta;
  int mu;

  CoordinateS2C(r, Theta, Phi, X);

  F(r,Theta,Phi,F0);
  VectorS2C(Theta, Phi, F0, F0C);

  Delta = (X[0]==0.0) ? 1.0e-4 : 1.0e-4*fabs(X[0]);
  X[0]+=Delta;
  CoordinateC2S(X, &rr, &tt, &pp);
  F(rr, tt, pp, F0);
  VectorS2C(tt, pp, F0, dFdX);
  for(mu=0; mu<3; mu++)
   dFdX[mu] = (dFdX[mu] - F0C[mu]) / Delta;
  X[0]-=Delta;

  Delta = (X[1]==0.0) ? 1.0e-4 : 1.0e-4*fabs(X[1]);
  X[1]+=Delta;
  CoordinateC2S(X, &rr, &tt, &pp);
  F(rr, tt, pp, F0);
  VectorS2C(tt, pp, F0, dFdY);
  for(mu=0; mu<3; mu++)
   dFdY[mu] = (dFdY[mu] - F0C[mu]) / Delta;
  X[1]-=Delta;

  Delta = (X[2]==0.0) ? 1.0e-4 : 1.0e-4*fabs(X[2]);
  X[2]+=Delta;
  CoordinateC2S(X, &rr, &tt, &pp);
  F(rr, tt, pp, F0);
  VectorS2C(tt, pp, F0, dFdZ);
  for(mu=0; mu<3; mu++)
   dFdZ[mu] = (dFdZ[mu] - F0C[mu]) / Delta;
  X[2]-=Delta;

  CurlFC[0]=dFdY[2] - dFdZ[1];
  CurlFC[1]=dFdZ[0] - dFdX[2];
  CurlFC[2]=dFdX[1] - dFdY[0];
  VectorC2S(Theta, Phi, CurlFC, CurlF2);
}
#endif

}
 
/***************************************************************/
/* compute the divergence of a vector-valued function.         */
/* Terms[0..2] are the three separate terms in the expression  */
/* for the divergence in spherical coordinates. Use the second */
/* of the two calling conventions if you don't need these.     */
/***************************************************************/
void GetDiv(double r, double Theta, double Phi,
            void (*F)(double r, double Theta, double Phi, cdouble *V),
            cdouble Terms[3], cdouble *DivF)
{
  cdouble F0[3], dFdr[3], dFdTheta[3], dFdPhi[3];
  double CT, ST;
  int mu;

  CT=cos(Theta);
  ST=sin(Theta);
  F(r,Theta,Phi,F0);
  F(r+1.0e-6, Theta, Phi, dFdr); 
  F(r, Theta+1.0e-6, Phi, dFdTheta); 
  F(r, Theta, Phi+1.0e-6, dFdPhi); 
  for(mu=0; mu<3; mu++)
   { dFdr[mu]     = (dFdr[mu]     - F0[mu]) / 1.0e-6;
     dFdTheta[mu] = (dFdTheta[mu] - F0[mu]) / 1.0e-6;
     dFdPhi[mu]   = (dFdPhi[mu]   - F0[mu]) / 1.0e-6;
   };
  Terms[0]=dFdr[0] + 2.0*F0[0]/r;
  Terms[1]=(CT*F0[1] + ST*dFdTheta[1])/(r*ST);
  Terms[2]=dFdPhi[2]/(r*ST);
  DivF[0]=Terms[0]+Terms[1]+Terms[2];
}



cdouble GetDiv(double r, double Theta, double Phi,
               void (*F)(double r, double Theta, double Phi, cdouble *V))
{ 
  cdouble Terms[3], DivF;
  GetDiv(r, Theta, Phi, F, Terms, &DivF);
  return DivF;
}
