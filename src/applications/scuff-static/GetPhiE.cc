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
 * GetPhiE.cc  -- routine for computing electrostatic potential and field
 *                due to a constant charge density on a single triangle,
 *                following these references:
 *                (1) Graglia, IEEE Trans. Ant. Prop. *41* 1448 (1993)
 *                (2) Wilton et al, IEEE Trans. Ant. Prop. *32* 276 (1984)
 *
 * homer reid  -- 5/2013 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "SSSolver.h"

namespace scuff {

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
double EdgeInt(double x, double ya, double yb, double z, double smin)
{ 
  if ( x<=1.0e-8 )
   return  (0.0<smin && smin<1.0)  ? M_PI*z : 0.0;

  double x2=x*x;
  double r0=sqrt(x2+z*z), ra=sqrt(ya*ya+r0*r0), rb=sqrt(yb*yb+r0*r0);
  double Ia=z*atan( ya*z/(x*ra) ) + x*log( (ya+ra)/r0 );
  double Ib=z*atan( yb*z/(x*rb) ) + x*log( (yb+rb)/r0 ); 
  if ( 0.0<=smin && smin<=1.00 )
   return Ia+Ib;
  else
   return fabs(Ia-Ib);
}

/***********************************************************************/
/* get potential and field due to a single panel.                     */
/*                                                                    */
/* note: the underlying algorithm here is based on a line-integral    */
/* technique that is applicable to polygons of an arbitrary number    */
/* of vertices, not just triangles. in this particular case we are    */
/* only applying it to the case of triangles, but for the purpose     */
/* of future generality we leave the number of polygon vertices as a  */
/* variable 'N' (equal to 3 in this case.)                            */
/* loop indices i,j run over cartesian coordinates (and thus are      */
/* delimited by 0..2) while loop index n runs over polygon edges      */
/* (and thus is delimited by 0..N-1).                                 */
/***********************************************************************/
void SSSolver::GetPhiE(int ns, int np, double *X, double PhiE[4])
{ 
  // int N=3; // hardcoded here but could change in the future 

  /***************************************************************/
  /* unpack panel vertices ***************************************/
  /***************************************************************/
  RWGSurface *S = G->Surfaces[ns];
  RWGPanel *P = S->Panels[np];
  double *PV[3]; // panel vertices
  PV[0] = S->Vertices + 3*(P->VI[0]);
  PV[1] = S->Vertices + 3*(P->VI[1]);
  PV[2] = S->Vertices + 3*(P->VI[2]);

  /***************************************************************/
  /* compute geometric quantities, graglia equations (3)--(5) and*/
  /* (53--59). Note: this code could be significantly optimized. */
  /***************************************************************/
  double uHat[3], vHat[3], nHat[3], rmP1[3], Temp[3], mHat[3][3];
  double u0, v0, w0, w0Sign, l1, l2, l3, u3, v3, nMag;
  VecSub(PV[1],PV[0],uHat);
  VecSub(PV[2],PV[1],Temp);
  VecCross(uHat, Temp, nHat);
  l3 = VecNormalize(uHat);
  l1 = VecNormalize(Temp);
  nMag = VecNormalize(nHat);
  VecCross(nHat, uHat, vHat);
  VecSub(PV[2], PV[0], Temp );
  l2 = VecNorm(Temp);
  u3 = VecDot(Temp, uHat);
  v3 = nMag / l3;
  VecSub(X,PV[0],rmP1);
  u0 = VecDot(uHat, rmP1);
  v0 = VecDot(vHat, rmP1);
  w0 = VecDot(nHat, rmP1);
  w0Sign = w0 >= 0.0 ? 1.0 : -1.0;
  VecCross( VecSub(PV[2], PV[1], Temp), nHat, mHat[0]); VecNormalize(mHat[0]);
  VecCross( VecSub(PV[0], PV[2], Temp), nHat, mHat[1]); VecNormalize(mHat[1]);
  VecCross( uHat, nHat, mHat[2]); VecNormalize(mHat[2]);
  
  double sPlus[3], sMinus[3], tZero[3], tPlus[3], tMinus[3]; 
  
  sMinus[0] = -( (l3-u0)*(l3-u3) + v0*v3 ) / l1;
  sPlus[0]  =  ( (u3-u0)*(u3-l3) + v3*(v3-v0) ) / l1;
  sMinus[1] = -( u3*(u3-u0) + v3*(v3-v0) ) / l2;
  sPlus[1]  =  ( u0*u3 + v0*v3 ) / l2;
  sMinus[2] = -u0;
  sPlus[2]  =  l3 - u0;
  tZero[0]  = (v0*(u3-l3) + v3*(l3-u0)) / l1;
  tZero[1]  = (u0*v3 - v0*u3) / l2;
  tZero[2]  = v0;
  tPlus[0]  = sqrt( (u3-u0)*(u3-u0) + (v3-v0)*(v3-v0) );
  tPlus[1]  = sqrt( u0*u0 + v0*v0 );
  tPlus[2]  = sqrt( (l3-u0)*(l3-u0) + v0*v0 );
  tMinus[0] = tPlus[2];
  tMinus[1] = tPlus[0];
  tMinus[2] = tPlus[1];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double f2i, Betai; 
  PhiE[0]=PhiE[1]=PhiE[2]=PhiE[3]=0.0;
  for(int n=0; n<3; n++)
   { 
     if ( fabs(tZero[n]) < 1.0e-6 ) // on panel edge 
      { Betai = f2i = 0.0;
      } 
     else if ( fabs(w0) < 1.0e-6 )
      { Betai = atan( sPlus[n] / tZero[n] ) - atan( sMinus[n] / tZero[n] );
        f2i   = log ( (tPlus[n] + sPlus[n]) / (tMinus[n] + sMinus[n]) ); 
      }
     else
      { double RPlus  = sqrt( tPlus[n]*tPlus[n]   + w0*w0 );
        double RMinus = sqrt( tMinus[n]*tMinus[n] + w0*w0 );
        double RZero  = sqrt( tZero[n]*tZero[n]   + w0*w0 );
        f2i    = log( (RPlus + sPlus[n]) / (RMinus + sMinus[n]) );
        Betai  = atan( (tZero[n]*sPlus[n])  / ( RZero*RZero + fabs(w0)*RPlus) )
                -atan( (tZero[n]*sMinus[n]) / ( RZero*RZero + fabs(w0)*RMinus) );
      };

     PhiE[0] += tZero[n]*f2i - fabs(w0)*Betai;
     PhiE[1] += w0Sign*Betai*nHat[0] + f2i*mHat[n][0];
     PhiE[2] += w0Sign*Betai*nHat[1] + f2i*mHat[n][1];
     PhiE[3] += w0Sign*Betai*nHat[2] + f2i*mHat[n][2];

   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (Substrate)
   AddSubstratePhiE(ns, np, X, PhiE);

}

} // namespace scuff
