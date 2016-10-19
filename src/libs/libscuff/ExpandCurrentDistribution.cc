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
#include "PanelCubature.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif

#define II cdouble(0.0,1.0)

namespace scuff {

void KNProjectionIntegrand(double x[3], double b[3], double Divb,
                           void *UserData, double Weight, double *Integral)
{
  (void )Divb;

  IncField *IF = (IncField *)UserData;
  cdouble EH[6];
  IF->GetFields(x,EH);
  cdouble *zIntegral=(cdouble *)Integral;
  zIntegral[0] += Weight*(b[0]*EH[0] + b[1]*EH[1] + b[2]*EH[2]);
  zIntegral[1] += Weight*(b[0]*EH[3] + b[1]*EH[4] + b[2]*EH[5]);
}

void EHProjectionIntegrand(double *x, PCData *PCD,
                           void *UserData, double *Integral)
{
  IncField *IF = (IncField *)UserData;

  cdouble EH[6], *E=EH+0, *H=EH+3;
  IF->GetFields(x,EH);

  cdouble *b   = PCD->K;
  double *nHat = PCD->nHat;

  cdouble K[3], N[3];
  K[0] = nHat[1]*H[2] - nHat[2]*H[1];
  K[1] = nHat[2]*H[0] - nHat[0]*H[2];
  K[2] = nHat[0]*H[1] - nHat[1]*H[0];

  N[0] = -(nHat[1]*E[2] - nHat[2]*E[1]);
  N[1] = -(nHat[2]*E[0] - nHat[0]*E[2]);
  N[2] = -(nHat[0]*E[1] - nHat[1]*E[0]);

  N[0] *= (-1.0/ZVAC);
  N[1] *= (-1.0/ZVAC);
  N[2] *= (-1.0/ZVAC);

  cdouble *zIntegral=(cdouble *)Integral;
  zIntegral[0] = b[0]*K[0] + b[1]*K[1] + b[2]*K[2];
  zIntegral[1] = b[0]*N[0] + b[1]*N[1] + b[2]*N[2];
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::ExpandCurrentDistribution(IncField *IF,
                                            HVector *KNVector,
                                            cdouble Omega,
                                            bool IsEHField)
{ 
  if (KNVector->ZV==0)
   ErrExit("%s:%i: internal error");

  /***************************************************************/
  /* project user's current distribution onto the RWG basis      */
  /* using AssembleRHSVector. Note we then need to               */
  /* undo the factors of -1/ZVAC and ZVAC that AssembleRHSVector */
  /* automatically puts into the KNVector                        */
  /***************************************************************/
  Log("ExpandCD: Assembling RHS");
#if 0
  AssembleRHSVector(Omega, IF, KNVector);
  for(int ns=0; ns<NumSurfaces; ns++)
   { RWGSurface *S = Surfaces[ns];
     for(int ne=0, nbf=0; ne<S->NumEdges; ne++)
      { KNVector->ZV[nbf++] *= -1.0*ZVAC;
        if(!S->IsPEC) KNVector->ZV[nbf++] /= ZVAC;
      };
   };
#endif

  Log("Computing projection of current onto RWG basis");
#ifdef USE_OPENMP
  int NumThreads=GetNumThreads();
  LogC(" (%i threads)",NumThreads);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
 for(int neFull=0; neFull<TotalEdges; neFull++)
  {
    int ns, ne, KNIndex;
    ResolveEdge(neFull, &ns, &ne, &KNIndex);
    cdouble KNProjections[2];
    int IDim=2*2;
    int Order=9;
    int MaxEvals=21;
    if (IsEHField)
     GetBFCubature(this, ns, ne,
                   EHProjectionIntegrand, (void *)IF, IDim,
                   MaxEvals, 0.0, 0.0, Omega,
                   0, (double *)KNProjections);
                   
    else
     GetBFCubature2(this, ns, ne,
                    KNProjectionIntegrand, (void *)IF, IDim,
                    Order, (double *)KNProjections);

    KNVector->SetEntry(KNIndex,KNProjections[0]);
    if ( !(Surfaces[ns]->IsPEC) )
     KNVector->SetEntry(KNIndex+1,KNProjections[1]);
  };
     

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int MaxBFs = Surfaces[0]->NumBFs;
  for(int ns=1; ns<NumSurfaces; ns++)
   if (Surfaces[ns]->NumBFs > MaxBFs)
    MaxBFs = Surfaces[ns]->NumBFs;
  cdouble *MBuffer = (cdouble *)mallocEC(MaxBFs*MaxBFs*sizeof(cdouble));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int ns=0; ns<NumSurfaces; ns++)
   { 
     RWGSurface *S = Surfaces[ns];
     int Offset    = BFIndexOffset[ns];
     int NE        = S->NumEdges;
     int NBF       = S->NumBFs;

     Log("ExpandCD: Assembling M%i",ns);
     HMatrix M(NBF, NBF, LHM_COMPLEX, LHM_NORMAL, (void *)MBuffer);
     for(int ne=0; ne<NE; ne++)
      for(int nep=ne; nep<NE; nep++)
       { 
         double OVLP=S->GetOverlap(ne, nep);
         if (S->IsPEC)
          { M.SetEntry(ne, nep, OVLP);
            M.SetEntry(nep, ne, OVLP);
          }
         else
          { M.SetEntry( 2*ne, 2*nep, OVLP);
            M.SetEntry( 2*nep, 2*ne, OVLP);
            M.SetEntry( 2*ne+1, 2*nep+1, OVLP);
            M.SetEntry( 2*nep+1, 2*ne+1, OVLP);
          };
       };

      Log("ExpandCD: LU factorizing");
      M.LUFactorize();

      HVector PartialKN(NBF, LHM_COMPLEX, (void *)(KNVector->ZV + Offset) );
      Log("ExpandCD: LU solving");
      M.LUSolve(&PartialKN);
   };
  Log("ExpandCD: done ");
  free(MBuffer);

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
                                          HVector *KNVector,
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
       { KAlpha = KNVector->GetEntry( Offset + ne );
         NAlpha = 0.0;
       }
      else
       { KAlpha = KNVector->GetEntry( Offset + 2*ne);
         NAlpha = -1.0*ZVAC*KNVector->GetEntry( Offset + 2*ne+1);
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
                                          HVector *KNVector, 
                                          cdouble KN[6])
{
  EvalCurrentDistribution(X, KNVector, 0, KN);
}

} // namespace scuff
