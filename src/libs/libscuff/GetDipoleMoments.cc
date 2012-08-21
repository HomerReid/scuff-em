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
 * GetDipoleMoments.cc  -- libscuff class methods for computing electric 
 *                      -- and magnetic dipole moments induced on 
 *                      -- scattering geometry by incident field
 *
 * homer reid    -- 3/2007  -- 11/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libTriInt.h>

#include "libscuff.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

namespace scuff {

#define II cdouble(0,1)


/***************************************************************/
/* get electric and magnetic dipole moments of a current       */
/* distribution described by a set of RWG basis functions      */
/* populated with a vector of basis-function weights.          */
/*                                                             */
/* The return value is an HVector of length 6N, where N is the */
/* number of surfaces in the geometry. The entries of this     */
/* vector are:                                                 */
/*                                                             */
/*  PM[ 6*n + 0..2 ] = x,y,z components of electric dipole     */
/*                     moment induced on nth surface           */
/*  PM[ 6*n + 3..5 ] = x,y,z components of magnetic dipole     */
/*                     moment induced on nth surface           */
/*                                                             */
/* If the input parameter PM is NULL on entry (or if PM points */
/* to an HVector of the wrong size), then the returned HVector */
/* is newly allocated.                                         */
/* Otherwise, the contents of PM are overwritten, and the      */
/* return value is PM.                                         */
/***************************************************************/
HVector *RWGGeometry::GetDipoleMoments(cdouble Omega, HVector *KN, HVector *PM)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( PM==0 || PM->N!=6*NumSurfaces)
   PM=new HVector(6*NumSurfaces, LHM_COMPLEX);
  PM->Zero(); 
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGSurface *S;
  RWGEdge *E;
  double QP[3], V1[3], V2[3], QM[3];
  double QPmQM[3], V1pV2[3], QQxVV[3];
  double PreFac;
  cdouble KAlpha, NAlpha;
  cdouble pRWG[3], mRWG[3];
  cdouble IK=II*Omega; 
  int nbf, ns, ne, Mu;
  for(nbf=ns=0; ns<NumSurfaces; ns++)
   for(S=Surfaces[ns], ne=0; ne<S->NumEdges; ne++)
    { 
      E=S->Edges[ne];

      memcpy(QP,S->Vertices + 3*E->iQP, 3*sizeof(double));
      memcpy(V1,S->Vertices + 3*E->iV1, 3*sizeof(double));
      memcpy(V2,S->Vertices + 3*E->iV2, 3*sizeof(double));
      memcpy(QM,S->Vertices + 3*E->iQM, 3*sizeof(double));
      if (S->GT)
       { S->GT->UnApply(QP);
         S->GT->UnApply(V1);
         S->GT->UnApply(V2);
         S->GT->UnApply(QM);
       };
      VecSub(QP, QM, QPmQM);
      VecAdd(V1, V2, V1pV2);
      VecCross(QPmQM, V1pV2, QQxVV);
      PreFac=E->Length / 3.0;

      pRWG[0] = PreFac * QPmQM[0] / (IK);
      pRWG[1] = PreFac * QPmQM[1] / (IK);
      pRWG[2] = PreFac * QPmQM[2] / (IK);
      mRWG[0] = PreFac * QQxVV[0] / 4.0;
      mRWG[1] = PreFac * QQxVV[1] / 4.0;
      mRWG[2] = PreFac * QQxVV[2] / 4.0;

      KAlpha=KN->GetEntry( nbf++ );
      if ( S->IsPEC )
       NAlpha = 0.0;
      else 
       NAlpha = -ZVAC*KN->GetEntry( nbf++ );

      for(Mu=0; Mu<3; Mu++)
       { PM->AddEntry( 6*ns + Mu + 0, KAlpha*pRWG[Mu] - NAlpha*mRWG[Mu]/ZVAC);
         PM->AddEntry( 6*ns + Mu + 3, KAlpha*mRWG[Mu] + NAlpha*pRWG[Mu]/ZVAC);
       };
 
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
PM->Scale(ZVAC); /* ? */ 
  return PM;

}

} // namespace scuff
