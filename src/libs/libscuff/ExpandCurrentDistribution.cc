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
 * ExpandCurrentDistribution.cc 
 *
 * homer reid  -- 12/2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <libhrutil.h>
#include <libTriInt.h>
#include <libhmat.h>

#include "libscuff.h"

namespace scuff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::ExpandCurrentDistribution(IncField *IF, HVector *KNVec, int nThread)
{ 
  int n, ne, nep;
  RWGObject *O;
  HMatrix *M;
  double OVLP;

  /* FIXME */
  if (NumObjects>1)
   ErrExit("%s:%i: ExpandCurrentDistribution only implemented for single-object geometries");
  O=Objects[0];

  Log("ExpandCD: Assembling RHS");
  AssembleRHSVector(1.0, IF, KNVec, nThread);
  if (KNVec->RealComplex==LHM_COMPLEX)
   { for(n=0; n<(KNVec->N)/2; n++)
      { KNVec->ZV[2*n] *= -1.0*ZVAC;
        KNVec->ZV[2*n+1] /= ZVAC;
      };
   }
  else
   { for(n=0; n<(KNVec->N)/2; n++)
      { KNVec->DV[2*n] *= -1.0*ZVAC;
        KNVec->DV[2*n+1] /= ZVAC;
      };
   };

  M=new HMatrix(2*O->NumEdges, 2*O->NumEdges, KNVec->RealComplex);
  M->Zero();

  Log("ExpandCD: Assembling M");
  for(ne=0; ne<O->NumEdges; ne++)
   for(nep=ne; nep<O->NumEdges; nep++)
    { 
      OVLP=O->GetOverlap(ne, nep);
      M->SetEntry( 2*ne, 2*nep, OVLP);
      M->SetEntry( 2*nep, 2*ne, OVLP);

      M->SetEntry( 2*ne+1, 2*nep+1, OVLP);
      M->SetEntry( 2*nep+1, 2*ne+1, OVLP);
    };

  Log("ExpandCD: LU factorizing");
  M->LUFactorize();

  Log("ExpandCD: LU solving");
  M->LUSolve(KNVec);
  
  Log("ExpandCD: done ");
  delete M;

}

/* return 1 if X lies inside the triangle with vertices V1, V2, V3. */
/* return 0 otherwise.                                              */
/* X is assumed to lie in the plane of the triangle.                */
int InsideTriangle(const double *X, const double *V1, const double *V2, const double *V3)
{
  double V1mX[3], V2mX[3], V3mX[3];
  double Length1, Length2, Length3;
  double Angle1, Angle2, Angle3;

  VecSub(V1,X,V1mX);
  VecSub(V2,X,V2mX);
  VecSub(V3,X,V3mX);

  Length1=VecNorm(V1mX);
  Length2=VecNorm(V2mX);
  Length3=VecNorm(V3mX);

  Angle1=acos( VecDot(V1mX, V2mX) / (Length1*Length2) );
  Angle2=acos( VecDot(V2mX, V3mX) / (Length2*Length3) );
  Angle3=acos( VecDot(V3mX, V1mX) / (Length3*Length1) );

  if ( fabs(Angle1+Angle2+Angle3 - 2.0*M_PI) < 1.0e-6 )
   return 1;
  return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::EvalCurrentDistribution(const double X[3], HVector *KNVec, cdouble KN[6])
{ 
  int no, ne, Offset;
  RWGObject *O;
  RWGEdge *E;
  double fRWG[3];
  double *QP, *V1, *V2, *QM, *Q;
  double Area, Sign;
  cdouble KAlpha, NAlpha;

  memset(KN,0,6*sizeof(cdouble));
  for(no=0; no<NumObjects; no++)
   for(O=Objects[no], Offset=BFIndexOffset[no], ne=0; ne<O->NumEdges; ne++)
    { 
      E=O->Edges[ne];
      QP=O->Vertices + 3*E->iQP;
      V1=O->Vertices + 3*E->iV1;
      V2=O->Vertices + 3*E->iV2;
      QM=O->Vertices + 3*E->iQM;

      if ( InsideTriangle(X,QP,V1,V2) )
       Q=QP;
      else if ( InsideTriangle(X,QM,V1,V2) )
       Q=QM;
      else
       continue;

      Sign = (Q==QP) ? +1.0 : -1.0;
      Area = (Q==QP) ? O->Panels[E->iPPanel]->Area : O->Panels[E->iMPanel]->Area;
   
      VecSub(X,Q,fRWG);
      VecScale(fRWG, E->Length / (2.0*Area) );

      if ( O->IsPEC )
       { KAlpha = KNVec->GetEntry( Offset + ne );
         NAlpha = 0.0;
       }
      else
       { KAlpha = KNVec->GetEntry( Offset + 2*ne);
         NAlpha = -1.0*ZVAC*KNVec->GetEntry( Offset + 2*ne+1);
       };

      KN[0] += Sign * KAlpha * fRWG[0]; 
      KN[1] += Sign * KAlpha * fRWG[1]; 
      KN[2] += Sign * KAlpha * fRWG[2]; 
      KN[3] += Sign * NAlpha * fRWG[0]; 
      KN[4] += Sign * NAlpha * fRWG[1]; 
      KN[5] += Sign * NAlpha * fRWG[2]; 

   }; // for(no=0; ... for(ne=0 ... 
}  

} // namespace scuff
