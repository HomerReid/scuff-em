/* Copyright (C) 2005-2011 M. T. Homer ReidElectricOnly=true;
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
 * GetAmplitudes.cc  --
 *                   --
 *
 * homer reid        -- 5/2015
 */
#include <stdio.h>
#include <stdlib.h>

#include "scuff-transmission.h"

using namespace scuff;

#define II cdouble (0.0, 1.0)

/***************************************************************/
/* f1 = \int_0^1 \int_0^u u*e^{-i(uX+vY)} dv du                */
/* f2 = \int_0^1 \int_0^u v*e^{-i(uX+vY)} dv du                */
/***************************************************************/
namespace scuff{ cdouble ExpRel(cdouble x, int n); }

void f1f2(double X, double Y, cdouble *f1, cdouble *f2)
{
  if (X==0.0 && Y==0.0)
   { *f1=1.0/3.0;
     *f2=1.0/6.0;
   }
  else if (Y==0.0)
   { 
     *f1 = 2.0*II*exp(-II*X)*ExpRel(II*X, 3) / (X*X*X);
     *f2 = (*f1)/2.0;
   }
  else if (X==0.0)
   { double Y3=Y*Y*Y;
     cdouble ER1=ExpRel(II*Y,1);
     cdouble ER2=ExpRel(II*Y,2);
     *f1= -II*exp(-II*Y)*ER2/Y3 - 0.5*II/Y;
     *f2= -2.0*II*exp(-II*Y)*( ER2 - 0.5*II*Y*ER1) / Y3;
   }
  else
   { 
     cdouble IX   = II*X;
     cdouble IY   = II*Y;
     double XPY   = X+Y;
     cdouble IXPY = II*XPY;

     cdouble ExpMIX=exp(-IX);
     cdouble ExpMIY=exp(-IY);
     cdouble ExpMIXPY=ExpMIX * ExpMIY;

     cdouble Term1 = X==0.0     ? -0.5 : ExpMIX*ExpRel(IX, 2) / (X*X);
     cdouble Term2 = Y==0.0     ? -0.5 : ExpMIY*ExpRel(IY, 2) / (Y*Y);
     cdouble Term3 = (XPY)==0.0 ? -0.5 : ExpMIXPY*ExpRel(IXPY, 2) / (XPY*XPY);

     *f1 = (Term1 - Term3) * II / Y;
  
     cdouble fFull = II*ExpRel(IY,2)*ExpRel(IX,1)*ExpMIXPY/(X*Y*Y);

     *f2 = fFull - (Term2 - Term3)*II/X;
   }
  
}

/***************************************************************/
/* compute the vector-valued integral                          */
/*  \int e^{-i*(q \cdot X)} b(x) dx                            */
/* where b is the RWG basis function associated with edge #ne  */
/* of the given RWGSurface.                                    */
/***************************************************************/
void GetbTwiddle(RWGSurface *S, int ne, double q[3], cdouble bTwiddle[3])
{
  RWGEdge *E    = S->Edges[ne];
  double *QP    = S->Vertices + 3*E->iQP;
  double *V1    = S->Vertices + 3*E->iV1;
  double *V2    = S->Vertices + 3*E->iV2;
  double *QM    = (E->iQM==-1) ? 0 : S->Vertices + 3*E->iQM;
  double Length = E->Length;

  double AP[3], AM[3], B[3];
  double qQP=0.0, qAP=0.0, qQM=0.0, qAM=0.0, qB=0.0;
  for(int Mu=0; Mu<3; Mu++)
   { 
     AP[Mu] = V1[Mu] - QP[Mu];
     B[Mu]  = V2[Mu] - V1[Mu];
      qQP  += q[Mu]*QP[Mu];
      qAP  += q[Mu]*AP[Mu];
       qB  += q[Mu]*B[Mu];

     if (!QM) continue;

     AM[Mu] = V1[Mu] - QM[Mu];
     qQM += q[Mu]*QM[Mu];
     qAM += q[Mu]*AM[Mu];
   };

  cdouble ExpFac, f1, f2;

  f1f2(qAP, qB, &f1, &f2);
  ExpFac = exp(-II*qQP);
  bTwiddle[0] = Length*ExpFac*( f1*AP[0] + f2*B[0] );
  bTwiddle[1] = Length*ExpFac*( f1*AP[1] + f2*B[1] );
  bTwiddle[2] = Length*ExpFac*( f1*AP[2] + f2*B[2] );

  if (!QM) return;

  f1f2(qAM, qB, &f1, &f2);
  ExpFac = exp(-II*qQM);
  bTwiddle[0] -= Length*ExpFac*( f1*AM[0] + f2*B[0] );
  bTwiddle[1] -= Length*ExpFac*( f1*AM[1] + f2*B[1] );
  bTwiddle[2] -= Length*ExpFac*( f1*AM[2] + f2*B[2] );
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPlaneWaveAmplitudes(RWGGeometry *G, HVector *KN,
                            cdouble Omega, double *kBloch,
                            int WhichRegion, bool IsUpper,
                            cdouble TETM[NUMPOLS], bool WriteByKNFile)
{
  if (G->LDim!=2)
   ErrExit("lattice must have 2D periodicity");

  if (    G->LBasis->GetEntryD(0,1)!=0.0
       || G->LBasis->GetEntryD(1,0)!=0.0
     )
   ErrExit("non-square lattices are not currently supported");

  double VUnitCell = G->LVolume;

  cdouble ZRel;
  cdouble nn=G->RegionMPs[WhichRegion]->GetRefractiveIndex(Omega, &ZRel);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double k0=real(nn*Omega);
  double k02=k0*k0;
  double kz2 = k02 - kBloch[0]*kBloch[0] - kBloch[1]*kBloch[1];
  if ( kz2 <= 0.0 ) // return zero amplitude for the evanescent case
   return;
 
  double kz = sqrt( kz2 );

  double q3D[3];
  q3D[0] = kBloch[0];
  q3D[1] = kBloch[1];
  q3D[2] = (IsUpper ? 1.0 : -1.0) * kz;

  /*--------------------------------------------------------------*/
  /*- loop over all edges on all surfaces that bound the region  -*/
  /*- in question to compute KTilde and NTilde (2D Fourier       -*/
  /*- transforms of surface currents).                           -*/
  /*--------------------------------------------------------------*/
  cdouble KTilde[3]={0.0,0.0,0.0}, NTilde[3]={0.0,0.0,0.0};
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     RWGSurface *S=G->Surfaces[ns];
     double Sign;
     if (S->RegionIndices[0]==WhichRegion)
      Sign=1.0;
     else if (S->RegionIndices[1]==WhichRegion)
      Sign=-1.0;
     else
      continue; // surface does not bound region

     for(int ne=0; ne<S->NumEdges; ne++)
      { 
        cdouble KAlpha, NAlpha;
        G->GetKNCoefficients(KN, ns, ne, &KAlpha, &NAlpha);

        cdouble bTwiddle[3]; 
        GetbTwiddle(S, ne, q3D, bTwiddle);

        VecPlusEquals(KTilde, Sign*KAlpha, bTwiddle);
        VecPlusEquals(NTilde, Sign*NAlpha, bTwiddle);
      };
   } 

  /*--------------------------------------------------------------*/
  /*- compute vector-matrix-vector products                       */
  /*-  eps * G * KTilde                                           */
  /*-  eps * C * NTilde                                           */
  /*- where eps = \eps^{TE}, \eps^{TM} are polarization vectors.  */
  /*--------------------------------------------------------------*/
  double EpsTE[3] = {0.0, 1.0, 0.0}, EpsTM[3];
  VecCross(EpsTE, q3D, EpsTM);
  VecScale(EpsTM, -1.0/k0);

  cdouble EGKTE=0.0, EGKTM=0.0;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    { cdouble GMuNu = -q3D[Mu]*q3D[Nu]/k02;
      if (Mu==Nu) GMuNu += 1.0;
      EGKTE += EpsTE[Mu]*GMuNu*KTilde[Nu];
      EGKTM += EpsTM[Mu]*GMuNu*KTilde[Nu];
    };
  EGKTE *= II/(2.0*kz);
  EGKTM *= II/(2.0*kz);

  cdouble ECNTE=0.0, ECNTM=0.0;
  for(int Mu=0; Mu<3; Mu++)
   { int MP1 = (Mu+1)%3;
     int MP2 = (Mu+2)%3;
     ECNTE += EpsTE[Mu]*(NTilde[MP1]*q3D[MP2] - NTilde[MP2]*q3D[MP1]);
     ECNTM += EpsTM[Mu]*(NTilde[MP1]*q3D[MP2] - NTilde[MP2]*q3D[MP1]);
   };
  ECNTE *= II/(2.0*k0*kz);
  ECNTM *= II/(2.0*k0*kz);
 
  TETM[POL_TE] = II*k0*(ZVAC*ZRel*EGKTE + ECNTE) / VUnitCell;
  TETM[POL_TM] = II*k0*(ZVAC*ZRel*EGKTM + ECNTM) / VUnitCell;

  if (WriteByKNFile)
   {
     FILE *f=vfopen("%s.byKN","r",GetFileBase(G->GeoFileName));
     if (f)
      { fprintf(f,"#  1,2,3  omega sin(theta) upper/lower\n");
        fprintf(f,"#  4,5    K contribution to aTE\n");
        fprintf(f,"#  6,7    N contribution to aTE\n");
        fprintf(f,"#  8,9    K contribution to aTM\n");
        fprintf(f,"# 10,11   N contribution to aTM\n");
        fclose(f);
        SetDefaultCD2SFormat("{%+.4e %+.4e}");
      };

     SetDefaultCD2SFormat("%e %e");
     f=vfopen("%s.byKN","a",GetFileBase(G->GeoFileName));
     fprintf(f,"%e %e %i ",real(Omega),(180.0/M_PI)*asin(kBloch[0]/k0),IsUpper);
     fprintf(f,"%s ",CD2S(II*k0*(ZVAC*ZRel*EGKTE)/VUnitCell));
     fprintf(f,"%s ",CD2S(II*k0*(ECNTE)/VUnitCell));
     fprintf(f,"%s ",CD2S(II*k0*(ZVAC*ZRel*EGKTM/VUnitCell)));
     fprintf(f,"%s ",CD2S(II*k0*(ECNTM)/VUnitCell));
     fprintf(f,"\n");
     fclose(f);
   };

}
