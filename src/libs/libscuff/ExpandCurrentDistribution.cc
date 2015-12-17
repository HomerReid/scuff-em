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

#include <libhrutil.h>
#include <libTriInt.h>
#include <libhmat.h>

#include "libscuff.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#define II cdouble(0.0,1.0)

namespace scuff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::ExpandCurrentDistribution(IncField *IF, HVector *KNVec, cdouble Omega)
{ 
  int n, ne, nep;
  RWGSurface *S;
  HMatrix *M;
  double OVLP;

  /* FIXME */
  if (NumSurfaces>1)
   ErrExit("%s:%i: ExpandCurrentDistribution only implemented for single-object geometries");
  S=Surfaces[0];

  Log("ExpandCD: Assembling RHS");
  AssembleRHSVector(Omega, IF, KNVec);
  // undo the factors of -1/ZVAC and ZVAC that AssembleRHSVector
  // automatically puts into the the KNVector
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

  M=new HMatrix(2*S->NumEdges, 2*S->NumEdges, KNVec->RealComplex);
  M->Zero();

  Log("ExpandCD: Assembling M");
  for(ne=0; ne<S->NumEdges; ne++)
   for(nep=ne; nep<S->NumEdges; nep++)
    { 
      OVLP=S->GetOverlap(ne, nep);
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

/********************************************************************/
/* return 0 if X lies outside the triangle with the given vertices, */
/* or a positive integer otherwise.                                 */
/*                                                                  */
/* Points on edges or vertices are considered to lie inside the     */
/* triangle.                                                        */
/*                                                                  */
/* X is assumed to lie in the plane of the triangle.                */
/*                                                                  */
/* If L is nonnull, then the triangle is translated through +L      */
/* (actually X is translated through -L).                           */
/********************************************************************/
#define IT_EXTERIOR 0
#define IT_ONVERTEX 1
#define IT_ONEDGE   2
#define IT_INTERIOR 3
int InsideTriangle(const double *X,
                   const double *V1, const double *V2, const double *V3,
                   const double *L=0)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  const double *XX; 
  double XXBuffer[3];
  if (L==0)
   XX=X;
  else
   { XXBuffer[0] = X[0]-L[0];
     XXBuffer[1] = X[1]-L[1];
     XXBuffer[2] = X[2]-L[2];
     XX=XXBuffer;
   };
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double V1mX[3], V2mX[3], V3mX[3];
  VecSub(V1,XX,V1mX);
  VecSub(V2,XX,V2mX);
  VecSub(V3,XX,V3mX);

  double Length1=VecNorm(V1mX);
  double Length2=VecNorm(V2mX);
  double Length3=VecNorm(V3mX);

  /***************************************************************/
  /* detect point on vertex. FIXME: the comparison here should be*/
  /* relative to an appropriate lengthscale, probably the minimum*/
  /* edge length.                                                */
  /***************************************************************/
  if ( Length1<=1.0e-6 || Length2<=1.0e-6 || Length3<=1.0e-6 )
   return IT_ONVERTEX;

  /***************************************************************/
  /* compute angles subtended at vertex pairs ********************/
  /***************************************************************/
  double Angle1=acos( ((float)VecDot(V1mX, V2mX)) / ((float)(Length1*Length2)) );
  double Angle2=acos( ((float)VecDot(V2mX, V3mX)) / ((float)(Length2*Length3)) );
  double Angle3=acos( ((float)VecDot(V3mX, V1mX)) / ((float)(Length3*Length1)) );

  /***************************************************************/
  /* detect point on edge  ***************************************/
  /***************************************************************/
  if ( EqualFloat(Angle1, M_PI ) ) return IT_ONEDGE;
  if ( EqualFloat(Angle2, M_PI ) ) return IT_ONEDGE;
  if ( EqualFloat(Angle3, M_PI ) ) return IT_ONEDGE;

  /***************************************************************/
  /* detect point in interior ************************************/
  /***************************************************************/
  if ( fabs(Angle1+Angle2+Angle3 - 2.0*M_PI) < 1.0e-6 )
   return IT_INTERIOR;

  return IT_EXTERIOR;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::EvalCurrentDistribution(const double X[3],
                                          HVector *KNVec,
                                          double *kBloch,
                                          cdouble KN[6])
{ 
  if (LBasis && !kBloch)
   ErrExit("%s:%i: null KBloch in PBC geometry",__FILE__,__LINE__);
  if (!LBasis && kBloch)
   ErrExit("%s:%i: non-null KBloch in non-PBC geometry",__FILE__,__LINE__);

  /***************************************************************/
  /* If a lattice is present, translate the evaluation point into*/
  /* the unit cell and compute the corresponding Bloch factor.   */
  /* (We find the lattice vector L such that X = X_0 + L with    */
  /* X_0 living in the unit cell; then currents and fields at X  */
  /* are related to those at X_0 according to                    */
  /* F(X) = F(X_0 + L ) = e^{ i k\cdot L } F(X_0).               */
  /* double BPArg=kBloch[0]*LVector[0];                          */
  /***************************************************************/
  double EvalPoint[3];
  cdouble BlochPhase=1.0;
  int LDim=0;
  double LBV[3][3];
  if (LBasis)
   { 
     LDim = LBasis->NC;
     for(int nd=0; nd<LDim; nd++)
      for(int j=0; j<3; j++)
       LBV[nd][j] = LBasis->GetEntryD(j,nd);

     double LVector[3];
     GetUnitCellRepresentative(X, EvalPoint, LVector);
     double BPArg = kBloch[0]*LVector[0];
     if (LDim>1) BPArg+=kBloch[1]*LVector[1];
     if (LDim>2) BPArg+=kBloch[2]*LVector[2];
     BlochPhase = exp( II * BPArg );
   }
  else
   { LDim=0;
     EvalPoint[0] = X[0];
     EvalPoint[1] = X[1];
     EvalPoint[2] = X[2];
   };
 
  /***************************************************************/
  /* Note: If the evaluation point lies inside the translate of  */
  /* the negative panel of a straddler edge, we need to include  */
  /* a Bloch phase factor and translate the QM vertex.           */
  /***************************************************************/
  memset(KN,0,6*sizeof(cdouble));
  RWGSurface *S; 
  int Offset, ne;
  cdouble StraddlerPhase=1.0;
  for(int ns=0; ns<NumSurfaces; ns++)
   for(S=Surfaces[ns], Offset=BFIndexOffset[ns], ne=0; ne<S->NumEdges; ne++)
    { 
      RWGEdge *E=S->Edges[ne];
      double *QP=S->Vertices + 3*E->iQP;
      double *V1=S->Vertices + 3*E->iV1;
      double *V2=S->Vertices + 3*E->iV2;
      double *QM= (E->iQM==-1) ? 0 : S->Vertices + 3*E->iQM;

      /***************************************************************/
      /* Note: If the point lies on an RWG edge, it will be deemed   */
      /* to lie 'inside' both the positive and negative triangles    */
      /* for that edge. To avoid double-counting, for on-edge points */
      /* we only add the contribution of the positive panel.         */
      /***************************************************************/
      double *Q, QBuffer[3];
      if ( InsideTriangle(EvalPoint,QP,V1,V2) )
       Q=QP;
      else if ( QM && InsideTriangle(EvalPoint,QM,V1,V2)==IT_INTERIOR )
       Q=QM;
      else if ( QM && LDim>=1 && InsideTriangle(EvalPoint,QM,V1,V2,LBV[0])==IT_INTERIOR )
       { 
         QBuffer[0] = QM[0] + LBV[0][0];
         QBuffer[1] = QM[1] + LBV[0][1];
         QBuffer[2] = QM[2] + LBV[0][2];
         Q = QBuffer;
         StraddlerPhase = exp( II * ( kBloch[0]*LBV[0][0] + kBloch[1]*LBV[0][1]) );
       }
      else if ( QM && LDim>=2 && InsideTriangle(EvalPoint,QM,V1,V2,LBV[1])==IT_INTERIOR )
       { 
         QBuffer[0] = QM[0] + LBV[1][0];
         QBuffer[1] = QM[1] + LBV[1][1];
         QBuffer[2] = QM[2] + LBV[1][2];
         Q = QBuffer;
         StraddlerPhase = exp( II * ( kBloch[0]*LBV[1][0] + kBloch[1]*LBV[1][1]) );
       }
      else
       continue;

      double Sign = (Q==QP) ? +1.0 : -1.0;
      double Area = (Q==QP) ? S->Panels[E->iPPanel]->Area : S->Panels[E->iMPanel]->Area;
   
      double fRWG[3];
      VecSub(EvalPoint,Q,fRWG);
      VecScale(fRWG, E->Length / (2.0*Area) );

      cdouble KAlpha, NAlpha;
      if ( S->IsPEC )
       { KAlpha = KNVec->GetEntry( Offset + ne );
         NAlpha = 0.0;
       }
      else
       { KAlpha = KNVec->GetEntry( Offset + 2*ne);
         NAlpha = -1.0*ZVAC*KNVec->GetEntry( Offset + 2*ne+1);
       };
      KAlpha *= BlochPhase * StraddlerPhase;
      NAlpha *= BlochPhase * StraddlerPhase;

      KN[0] += Sign * KAlpha * fRWG[0]; 
      KN[1] += Sign * KAlpha * fRWG[1]; 
      KN[2] += Sign * KAlpha * fRWG[2]; 
      KN[3] += Sign * NAlpha * fRWG[0]; 
      KN[4] += Sign * NAlpha * fRWG[1]; 
      KN[5] += Sign * NAlpha * fRWG[2]; 

   }; // for(ns=0; ... for(ne=0 ... 
}  

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::EvalCurrentDistribution(const double X[3], 
                                          HVector *KNVec, 
                                          cdouble KN[6])
{
  EvalCurrentDistribution(X, KNVec, 0, KN);
}

} // namespace scuff
