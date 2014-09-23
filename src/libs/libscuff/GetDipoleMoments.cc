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

/***************************************************************/
/* This routine computes the induced electric and magnetic     */
/* surface charge and current densities at the centroid of     */
/* each panel in an RWGGeometry. As an added bonus, it         */
/* computes the normally-directed poynting flux into the       */
/* surface. (This is the flux, i.e. the power per unit area;   */
/* to get the total power delivered to the object through this */
/* panel you multiply the flux by the panel area).             */
/*                                                             */
/* PSD is a complex-valued matrix with G->TotalPanels rows and */
/* 13 columns.                                                 */
/*                                                             */
/* (If PSD is NULL on entry, or if it points to an HMatrix of  */
/*  the wrong size, it is reallocated).                        */
/*                                                             */
/* On output, the columns of the nth row of PSD are filled in  */
/* with data on the values of the induced source distributions */
/* at the centroid of the nth panel in the geometry, as follows*/
/*                                                             */
/*  columns 0, 1, 2 : cartesian coords of nth panel centroid   */
/*  column  3       : area of nth panel                        */
/*  column  4       : electric charge density                  */
/*  columns 5,6,7   : components of electric current density   */
/*  column  8       : magnetic charge density                  */
/*  columns 9,10,11 : components of magnetic current density   */
/*  column  12      : normally-directed poynting flux INTO the */
/*                    surface                                  */
/***************************************************************/
#define PSD_MATRIX_WIDTH 13
HMatrix *RWGGeometry::GetPanelSourceDensities2(cdouble Omega,
                                               HVector *KN,
                                               HMatrix *PSD)

{ 
  cdouble iw = II*Omega;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (    PSD==0 
       || PSD->NR!=TotalPanels 
       || PSD->NC!=PSD_MATRIX_WIDTH 
       || PSD->RealComplex!=LHM_COMPLEX
     )
   { if (PSD) 
      Warn("invalid PSD matrix passed to GetPanelSourceDensities (reallocating)");
     PSD=new HMatrix(TotalPanels, 13, LHM_COMPLEX);
   };
  PSD->Zero();

  /***************************************************************/
  /* loop over entries in the KN vector to add the               */
  /* contributions of each basis function to the charge and      */
  /* current densities at the centroids of its two panels        */
  /***************************************************************/
  RWGSurface *S;
  RWGEdge *E;
  RWGPanel *PPanel, *MPanel;
  double *QP, *QM;
  int ne, Offset;
  cdouble KWeight, NWeight, PreFac;
  double XmQ[3];
  int nbf=0;
  for(int ns=0; ns<NumSurfaces; ns++)
   for(S=Surfaces[ns], Offset=PanelIndexOffset[ns], ne=0; ne<S->NumEdges; ne++)
    { 
      E=S->Edges[ne];
      PPanel = S->Panels[E->iPPanel];
      MPanel = S->Panels[E->iMPanel];
      QP=S->Vertices + 3*(E->iQP);
      QM=S->Vertices + 3*(E->iQM);

      /*--------------------------------------------------------------*/
      /*- this code ensures that the panel centroid coordinates are   */
      /*- filled in redundantly several times over, but i think this  */
      /*- is faster than having a separate loop over panels.          */
      /*--------------------------------------------------------------*/
      PSD->SetEntry(PPanel->Index + Offset, 0, PPanel->Centroid[0]);
      PSD->SetEntry(PPanel->Index + Offset, 1, PPanel->Centroid[1]);
      PSD->SetEntry(PPanel->Index + Offset, 2, PPanel->Centroid[2]);
      PSD->SetEntry(PPanel->Index + Offset, 3, PPanel->Area);
      PSD->SetEntry(MPanel->Index + Offset, 0, MPanel->Centroid[0]);
      PSD->SetEntry(MPanel->Index + Offset, 1, MPanel->Centroid[1]);
      PSD->SetEntry(MPanel->Index + Offset, 2, MPanel->Centroid[2]);
      PSD->SetEntry(MPanel->Index + Offset, 3, MPanel->Area);

      /*--------------------------------------------------------------*/
      /*- contributions of electric surface current associated with  -*/
      /*- this edge                                                  -*/
      /*--------------------------------------------------------------*/
      KWeight = KN->GetEntry(nbf++);

      PreFac  = E->Length * KWeight / (2.0*PPanel->Area);
      VecSub(PPanel->Centroid, QP, XmQ);
      PSD->AddEntry( PPanel->Index + Offset, 4, 2.0*PreFac/ iw);  // rho 
      PSD->AddEntry( PPanel->Index + Offset, 5, PreFac*XmQ[0] );  // K_x
      PSD->AddEntry( PPanel->Index + Offset, 6, PreFac*XmQ[1] );  // K_y
      PSD->AddEntry( PPanel->Index + Offset, 7, PreFac*XmQ[2] );  // K_z

      PreFac  = -1.0*E->Length * KWeight / (2.0*MPanel->Area);
      VecSub(MPanel->Centroid, QM, XmQ);
      PSD->AddEntry( MPanel->Index + Offset, 4, 2.0*PreFac/ iw);  // rho 
      PSD->AddEntry( MPanel->Index + Offset, 5, PreFac*XmQ[0] );  // K_x
      PSD->AddEntry( MPanel->Index + Offset, 6, PreFac*XmQ[1] );  // K_y
      PSD->AddEntry( MPanel->Index + Offset, 7, PreFac*XmQ[2] );  // K_z

      if (S->IsPEC) continue;

      /*--------------------------------------------------------------*/
      /*- contributions of magnetic surface current associated with  -*/
      /*- this edge                                                  -*/
      /*--------------------------------------------------------------*/
      NWeight = -ZVAC*KN->GetEntry(nbf++);

      PreFac  = E->Length * NWeight / (2.0*PPanel->Area);
      VecSub(PPanel->Centroid, QP, XmQ);
      PSD->AddEntry( PPanel->Index + Offset, 8,  2.0*PreFac/ iw);  // eta
      PSD->AddEntry( PPanel->Index + Offset, 9,  PreFac*XmQ[0] );  // N_x
      PSD->AddEntry( PPanel->Index + Offset, 10, PreFac*XmQ[1] );  // N_y
      PSD->AddEntry( PPanel->Index + Offset, 11, PreFac*XmQ[2] );  // N_z

      PreFac  = -1.0*E->Length * NWeight / (2.0*MPanel->Area);
      VecSub(MPanel->Centroid, QM, XmQ);
      PSD->AddEntry( MPanel->Index + Offset,  8, 2.0*PreFac/ iw);  // eta
      PSD->AddEntry( MPanel->Index + Offset,  9, PreFac*XmQ[0] );  // N_x
      PSD->AddEntry( MPanel->Index + Offset, 10, PreFac*XmQ[1] );  // N_y
      PSD->AddEntry( MPanel->Index + Offset, 11, PreFac*XmQ[2] );  // N_z

    }; // for (ns=0; ns<NumSurfaces; ns++) ... for(ne=0; ne<S->NumEdges; ne++)

  return PSD;
}

/***************************************************************/
/* This routine computes the induced electric and magnetic     */
/* surface charge and current densities at the centroid of     */
/* each panel in an RWGGeometry. As an added bonus, it         */
/* computes the normally-directed poynting flux into the       */
/* surface. (This is the flux, i.e. the power per unit area;   */
/* to get the total power delivered to the object through this */
/* panel you multiply the flux by the panel area).             */
/*                                                             */
/* PSD is a complex-valued matrix with G->TotalPanels rows and */
/* 13 columns.                                                 */
/*                                                             */
/* (If PSD is NULL on entry, or if it points to an HMatrix of  */
/*  the wrong size, it is reallocated).                        */
/*                                                             */
/* On output, the columns of the nth row of PSD are filled in  */
/* with data on the values of the induced source distributions */
/* at the centroid of the nth panel in the geometry, as follows*/
/*                                                             */
/*  columns 0, 1, 2 : cartesian coords of nth panel centroid   */
/*  column  3       : area of nth panel                        */
/*  column  4       : electric charge density                  */
/*  columns 5,6,7   : components of electric current density   */
/*  column  8       : magnetic charge density                  */
/*  columns 9,10,11 : components of magnetic current density   */
/*  column  12      : normally-directed poynting flux INTO the */
/*                    surface                                  */
/***************************************************************/
HMatrix *RWGGeometry::GetPBCPanelSourceDensities(cdouble Omega,
                                                 double *kBloch,
                                                 HVector *KN,
                                                 HMatrix *PSD)
{ 
  cdouble iw = II*Omega;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (    PSD==0 
       || PSD->NR!=TotalPanels 
       || PSD->NC!=PSD_MATRIX_WIDTH 
       || PSD->RealComplex!=LHM_COMPLEX
     )
   { if (PSD) 
      Warn("invalid PSD matrix passed to GetPanelSourceDensities (reallocating)");
     PSD=new HMatrix(TotalPanels, PSD_MATRIX_WIDTH, LHM_COMPLEX);
   };
  PSD->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGSurface *S;
  int ns, np, POffset, BFOffset;
  for(ns=0; ns<NumSurfaces; ns++)
   for(S=Surfaces[ns], POffset=PanelIndexOffset[ns], BFOffset=BFIndexOffset[ns], np=0; np<S->NumPanels; np++)
    { 
      RWGPanel *P  = S->Panels[np];

      /*--------------------------------------------------------------*/
      /* add the contributions of all panel edges (i.e. RWG basis     */
      /* functions associated with all panel edges) to source         */
      /* densities at panel centroid                                  */
      /*--------------------------------------------------------------*/
      cdouble K[3], Sigma, N[3], Eta;
      K[0]=K[1]=K[2]=Sigma=N[0]=N[1]=N[2]=Eta=0.0;
      for(int nce=0; nce<3; nce++) // 'number of contributing edges'
       { 
         int ne = P->EI[nce];
         if (ne < 0) continue; // panel edge #nce is an exterior edge

         // get the value of the RWG basis function associated with panel edge #nce
         // at the panel centroid
         RWGEdge *E    = S->Edges[ne];
         double *Q     = S->Vertices + 3*(P->VI[nce]);
         double Sign   = ( (np == E->iMPanel) ? -1.0 : 1.0);
         double PreFac = Sign * E->Length / (2.0*P->Area);

         double fRWG[3];
         fRWG[0] = PreFac * (P->Centroid[0] - Q[0]);
         fRWG[1] = PreFac * (P->Centroid[1] - Q[1]);
         fRWG[2] = PreFac * (P->Centroid[2] - Q[2]);
         
         // look up the coefficients of this RWG basis function in the 
         // expansions of the electric and magnetic surface currents
         cdouble KAlpha, NAlpha;
         if (S->IsPEC)
          { KAlpha = KN->GetEntry(BFOffset + ne);
            NAlpha = 0.0;
          }
         else
          { KAlpha =       KN->GetEntry(BFOffset + 2*ne + 0);
            NAlpha = -ZVAC*KN->GetEntry(BFOffset + 2*ne + 1);
          };

         // add the contributions of this RWG basis function to 
         // the source densities at the panel centroid
         K[0]  += KAlpha * fRWG[0];
         K[1]  += KAlpha * fRWG[1];
         K[2]  += KAlpha * fRWG[2];
         Sigma += 2.0*KAlpha*PreFac / iw;
         N[0]  += NAlpha * fRWG[0];
         N[1]  += NAlpha * fRWG[1];
         N[2]  += NAlpha * fRWG[2];
         Eta   += 2.0*NAlpha*PreFac / iw;
       };

      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      /*--------------------------------------------------------------*/
      if (kBloch)
       { 
         for(int nst=0; nst<S->TotalStraddlers; nst++)
          if ( S->PhasedBFCs[3*nst + 0] == np )
           { int ne      = S->PhasedBFCs[3*nst + 1];
             int WhichBV = S->PhasedBFCs[3*nst + 2];
   
             RWGEdge *E  = S->Edges[ne];
             double *L   = LBasis[WhichBV];
  
             double Q[3];
             Q[0] = S->Vertices[3*(E->iQM)+0] + L[0];
             Q[1] = S->Vertices[3*(E->iQM)+1] + L[1];
             Q[2] = S->Vertices[3*(E->iQM)+2] + L[2];
  
             double Sign = -1.0;
             double PreFac = Sign * E->Length / (2.0*P->Area);
  
             double fRWG[3];
             fRWG[0] = PreFac * (P->Centroid[0] - Q[0]);
             fRWG[1] = PreFac * (P->Centroid[1] - Q[1]);
             fRWG[2] = PreFac * (P->Centroid[2] - Q[2]);
  
             cdouble BlochPhase = exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );
  
             cdouble KAlpha, NAlpha;
             if (S->IsPEC)
              { KAlpha = KN->GetEntry(BFOffset + ne);
                NAlpha = 0.0;
              }
             else
              { KAlpha =       KN->GetEntry(BFOffset + 2*ne + 0);
                NAlpha = -ZVAC*KN->GetEntry(BFOffset + 2*ne + 1);
              };
             KAlpha*=BlochPhase;
             NAlpha*=BlochPhase;

             // add the contributions of this RWG basis function to 
             // the source densities at the panel centroid
             K[0]  += KAlpha * fRWG[0];
             K[1]  += KAlpha * fRWG[1];
             K[2]  += KAlpha * fRWG[2];
             Sigma += 2.0*KAlpha*PreFac / iw;
             N[0]  += NAlpha * fRWG[0];
             N[1]  += NAlpha * fRWG[1];
             N[2]  += NAlpha * fRWG[2];
             Eta   += 2.0*NAlpha*PreFac / iw;

         }; // for(int ns=0; ns<S->NumStraddlers; ns++)

       };

      /*--------------------------------------------------------------*/
      /* compute the inward-directed normal poynting flux             */
      /*  = (1/2) re K^* \cdot (nHat \times N)                        */
      /*--------------------------------------------------------------*/
      double *nHat = P->ZHat;
      double IDNPF = 0.5*real( conj(K[0]) * (nHat[1]*N[2] - nHat[2]*N[1])  
                              +conj(K[1]) * (nHat[2]*N[0] - nHat[0]*N[2])  
                              +conj(K[2]) * (nHat[0]*N[1] - nHat[1]*N[0]));

      /*--------------------------------------------------------------*/
      /*- fill in the rows of the matrix corresponding to this panel  */
      /*--------------------------------------------------------------*/
      PSD->SetEntry(POffset + np,  0, P->Centroid[0]);
      PSD->SetEntry(POffset + np,  1, P->Centroid[1]);
      PSD->SetEntry(POffset + np,  2, P->Centroid[2]);
      PSD->SetEntry(POffset + np,  3, P->Area);
      PSD->SetEntry(POffset + np,  4, Sigma);
      PSD->SetEntry(POffset + np,  5, K[0]);
      PSD->SetEntry(POffset + np,  6, K[1]);
      PSD->SetEntry(POffset + np,  7, K[2]);
      PSD->SetEntry(POffset + np,  8, Eta);
      PSD->SetEntry(POffset + np,  9, N[0]);
      PSD->SetEntry(POffset + np, 10, N[1]);
      PSD->SetEntry(POffset + np, 11, N[2]);
      PSD->SetEntry(POffset + np, 12, IDNPF);

    }; // for (ns=0; ns<NumSurfaces; ns++) ... for(np=0; np<S->NumPanels; np++)

  return PSD;
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
HMatrix *RWGGeometry::GetPanelSourceDensities(cdouble Omega, HVector *KN, HMatrix *PSD)
{ return GetPBCPanelSourceDensities(Omega, 0, KN, PSD); }

} // namespace scuff
