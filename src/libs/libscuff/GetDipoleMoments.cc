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
/* number of objects in the geometry. The entries of this      */
/* vector are:                                                 */
/*                                                             */
/*  PM[ 6*n + 0..2 ] = x,y,z components of electric dipole     */
/*                     moment induced on nth object            */
/*  PM[ 6*n + 3..5 ] = x,y,z components of magnetic dipole     */
/*                     moment induced on nth object            */
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
  if ( PM==0 || PM->N!=6*NumObjects )
   PM=new HVector(6*NumObjects, LHM_COMPLEX);
  PM->Zero(); 
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGObject *O;
  RWGEdge *E;
  double QP[3], V1[3], V2[3], QM[3];
  double QPmQM[3], V1pV2[3], QQxVV[3];
  double PreFac;
  cdouble KAlpha, NAlpha;
  cdouble pRWG[3], mRWG[3];
  cdouble IK=II*Omega; 
  int nbf, no, ne, Mu;
  for(nbf=no=0; no<NumObjects; no++)
   for(O=Objects[no], ne=0; ne<O->NumEdges; ne++)
    { 
      E=O->Edges[ne];

      memcpy(QP,O->Vertices + 3*E->iQP, 3*sizeof(double));
      memcpy(V1,O->Vertices + 3*E->iV1, 3*sizeof(double));
      memcpy(V2,O->Vertices + 3*E->iV2, 3*sizeof(double));
      memcpy(QM,O->Vertices + 3*E->iQM, 3*sizeof(double));
      if (O->GT)
       { O->GT->UnApply(QP);
         O->GT->UnApply(V1);
         O->GT->UnApply(V2);
         O->GT->UnApply(QM);
       };
      VecSub(QP, QM, QPmQM);
      VecAdd(V1, V2, V1pV2);
      VecCross(QPmQM, V1pV2, QQxVV);
      PreFac=E->Length / 3.0;

      pRWG[0] = PreFac * QPmQM[0] / (IK);
      pRWG[1] = PreFac * QPmQM[1] / (IK);
      pRWG[2] = PreFac * QPmQM[2] / (IK);
      mRWG[0] = PreFac * QQxVV[0];
      mRWG[1] = PreFac * QQxVV[1];
      mRWG[2] = PreFac * QQxVV[2];

      KAlpha=KN->GetEntry( nbf++ );
      if ( O->MP->Type==MP_PEC )
       NAlpha=0.0;
      else 
       NAlpha=KN->GetEntry( nbf++ ) * (-1.0*ZVAC);

      for(Mu=0; Mu<3; Mu++)
       { PM->AddEntry( 6*no + Mu + 0, KAlpha*pRWG[Mu] + NAlpha*mRWG[Mu]/ZVAC);
         PM->AddEntry( 6*no + Mu + 3, KAlpha*mRWG[Mu] - NAlpha*pRWG[Mu]/ZVAC);
       };
 
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  return PM;

}

} // namespace scuff
