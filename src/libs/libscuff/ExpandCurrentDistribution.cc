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
void RWGGeometry::ExpandCurrentDistribution(EHFuncType KNFunc, void *KNFuncUD,
                                            int nThread, HVector *KNVec)
{ 
  int n, no, nop, ne, nep;
  RWGObject *O;
  int Offset, OffsetP;
  HMatrix *M;
  double OVLP;

  /* FIXME */
  if (NumObjects>1)
   ErrExit("%s:%i: ExpandCurrentDistribution only implemented for single-object geometries");
  O=Objects[0];

  Log("ExpandCD: Assembling RHS");
  AssembleRHSVector(KNFunc, KNFuncUD, nThread, KNVec);
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
int InsideTriangle(double *X, double *V1, double *V2, double *V3)
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
void RWGGeometry::EvalCurrentDistribution(double *X, HVector *KNVec, cdouble *KN)
{ 
  int ne; 
  RWGObject *O;
  RWGEdge *E;
  RWGPanel *P;
  double fRWG[3];
  double *QP, *V1, *V2, *QM;
  double PArea, MArea;
  cdouble KAlpha, NAlpha;

  /* FIXME */
  if (NumObjects>1)
   ErrExit("%s:%i: EvalCurrentDistribution only implemented for single-object geometries");
  O=Objects[0];

  memset(KN,0,6*sizeof(cdouble));
  for(ne=0; ne<O->NumEdges; ne++)
   { 
     E=O->Edges[ne];
     QP=O->Vertices + 3*E->iQP;
     V1=O->Vertices + 3*E->iV1;
     V2=O->Vertices + 3*E->iV2;
     QM=O->Vertices + 3*E->iQM;
     PArea=O->Panels[E->iPPanel]->Area;
     MArea=O->Panels[E->iMPanel]->Area;

     /* if X lies within positive panel */
     if ( InsideTriangle(X,QP,V1,V2) )
      { VecSub(X,QP,fRWG);
        VecScale(fRWG, E->Length / (2.0*PArea) );
        KAlpha=KNVec->GetEntry(2*ne);
        NAlpha=-1.0*ZVAC*KNVec->GetEntry(2*ne+1);
        KN[0] += KAlpha * fRWG[0]; 
        KN[1] += KAlpha * fRWG[1]; 
        KN[2] += KAlpha * fRWG[2]; 
        KN[3] += NAlpha * fRWG[0]; 
        KN[4] += NAlpha * fRWG[1]; 
        KN[5] += NAlpha * fRWG[2]; 
      }
     else if ( InsideTriangle(X,QM,V1,V2) )
      { VecSub(X,QM,fRWG);
        VecScale(fRWG, E->Length / (2.0*MArea) );
        KAlpha=KNVec->GetEntry(2*ne);
        NAlpha=-1.0*ZVAC*KNVec->GetEntry(2*ne+1);
        KN[0] -= KAlpha * fRWG[0]; 
        KN[1] -= KAlpha * fRWG[1]; 
        KN[2] -= KAlpha * fRWG[2]; 
        KN[3] -= NAlpha * fRWG[0]; 
        KN[4] -= NAlpha * fRWG[1]; 
        KN[5] -= NAlpha * fRWG[2]; 
      };

   };
}  

} // namespace scuff
