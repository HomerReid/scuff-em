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
  cdouble iw=II*Omega; 
  RWGSurface *S;
  int ne;
  for(int ns=0; ns<NumSurfaces; ns++)
   for(S=Surfaces[ns], ne=0; ne<S->NumEdges; ne++)
    { 
      RWGEdge *E=S->Edges[ne];

      double QP[3], V1[3], V2[3], QM[3]={0.0,0.0,0.0};
      memcpy(QP,S->Vertices + 3*E->iQP, 3*sizeof(double));
      memcpy(V1,S->Vertices + 3*E->iV1, 3*sizeof(double));
      memcpy(V2,S->Vertices + 3*E->iV2, 3*sizeof(double));
      if ( E->iQM!= -1 )
       memcpy(QM,S->Vertices + 3*E->iQM, 3*sizeof(double));
      if (S->GT)
       { S->GT->UnApply(QP);
         S->GT->UnApply(V1);
         S->GT->UnApply(V2);
         S->GT->UnApply(QM);
       };

      cdouble pRWG[3]={0.0, 0.0, 0.0}, mRWG[3]={0.0, 0.0, 0.0};
      double QPmQM[3], QQxX0[3];
      double PreFac=E->Length / 3.0;
      if ( E->iQM == -1 )
       VecSub(QP, E->Centroid, QPmQM);
      else
       VecSub(QP, QM, QPmQM);
      VecCross(QPmQM, E->Centroid, QQxX0);
      PlusEqualsVec(pRWG, PreFac/iw,  QPmQM);
      PlusEqualsVec(mRWG, PreFac/2.0, QQxX0);

      cdouble KAlpha, NAlpha;
      GetKNCoefficients(KN, ns, ne, &KAlpha, &NAlpha);

      for(int Mu=0; Mu<3; Mu++)
       { PM->AddEntry(6*ns + Mu + 0, KAlpha*pRWG[Mu] - NAlpha*mRWG[Mu]/ZVAC);
         PM->AddEntry(6*ns + Mu + 3, KAlpha*mRWG[Mu] + NAlpha*pRWG[Mu]/ZVAC);
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
/*                                                             */
/* Note: Because G->TotalPanels includes straddler panels,     */
/*       the number of rows in the PSD matrix is the total     */
/*       number of panels including straddler panels. However, */
/*       the contributions of straddler basis functions to the */
/*       source densities is already properly included in the  */
/*       values computed for the non-straddler panels that are */
/*       the periodic images of stradder panels. For this      */
/*       reason, all rows of the PSD matrix corresponding to   */
/*       straddler panels are returned as zeros.               */
/***************************************************************/
#define PSD_MATRIX_COLUMNS 13
HMatrix *RWGGeometry::GetPanelSourceDensities(cdouble Omega,
                                              double *kBloch,
                                              HVector *KN,
                                              HMatrix *PSD)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (    PSD==0 
       || PSD->NR!=TotalPanels 
       || PSD->NC!=PSD_MATRIX_COLUMNS
       || PSD->RealComplex!=LHM_COMPLEX
     )
   { if (PSD) 
      Warn("invalid PSD matrix passed to GetPanelSourceDensities (reallocating)");
     PSD=new HMatrix(TotalPanels, PSD_MATRIX_COLUMNS, LHM_COMPLEX);
   };
  PSD->Zero();

  /***************************************************************/
  /* loop over all panels on all surfaces ************************/
  /***************************************************************/
  cdouble iw = II*Omega;
  RWGSurface *S;
  int np;
  for(int ns=0; ns<NumSurfaces; ns++)
   for(S=Surfaces[ns], np=0; np<S->NumPanels; np++)
    { 
      /*--------------------------------------------------------------*/
      /* add the contributions of all three panel edges (i.e. RWG     */
      /* functions associated with all panel edges) to source         */
      /* densities at panel centroid                                  */
      /*--------------------------------------------------------------*/
      cdouble K[3]={0.0,0.0,0.0}, Sigma=0.0; 
      cdouble N[3]={0.0,0.0,0.0}, Eta=0.0;
      RWGPanel *P  = S->Panels[np];
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
         GetKNCoefficients(KN, ns, ne, &KAlpha, &NAlpha);

         // add the contributions of this RWG basis function to 
         // the source densities at the panel centroid
         PlusEqualsVec(K, KAlpha, fRWG);
         Sigma += 2.0*KAlpha*PreFac / iw;
         PlusEqualsVec(N, NAlpha, fRWG);
         Eta   += 2.0*NAlpha*PreFac / iw;
       }; // for(int nce=0; nce<3; nce++)

      /*--------------------------------------------------------------*/
      /*- store quantities associated with panel in the PSD matrix  --*/
      /*--------------------------------------------------------------*/
      int POffset=PanelIndexOffset[ns];
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

    }; // for(int ns=... for (int np=...

  /***************************************************************/
  /* handle straddler panels on periodic geometries              */
  /***************************************************************/
  if (kBloch)
   { 
     for(int ns=0; ns<NumSurfaces; ns++)
      { 
        S=Surfaces[ns];
        int Offset = PanelIndexOffset[ns];
        for(int nst=0; nst<S->TotalStraddlers; nst++)
         { 
           int npDest   = S->PhasedBFCs[3*nst + 0] + Offset;
           int ne       = S->PhasedBFCs[3*nst + 1];
           int WhichBV  = S->PhasedBFCs[3*nst + 2];

           int npSource = Offset + S->Edges[ne]->iMPanel;

           double L[3]; 
           L[0] = LBasis->GetEntryD(0,WhichBV);
           L[1] = LBasis->GetEntryD(1,WhichBV);
           L[2] = LBasis->GetEntryD(2,WhichBV);
           // folowing assumes that kBloch[2]==; FIXME for LDim=3
           cdouble BlochPhase = exp( II*(kBloch[0]*L[0] + kBloch[1]*L[1]) );

           for(int nc=4; nc<=11; nc++)
            { PSD->AddEntry(npDest,   nc, BlochPhase*PSD->GetEntry(npSource, nc) );
              PSD->SetEntry(npSource, nc, 0.0);
            }; 
         }; 
      };
   }; //  if (kBloch)...for(int ns=0; ns<NumSurfaces; ns++)...

  /***************************************************************/
  /* go back and compute the inward-directed normal pointing flux*/
  /*  = (1/2) re K^* \cdot (nHat \times N)                       */
  /***************************************************************/
  for(int ns=0; ns<NumSurfaces; ns++)
   { S=Surfaces[ns];
     int Offset = PanelIndexOffset[ns];
     for(np=0; np<S->NumPanels; np++)
      { double *nHat = S->Panels[np]->ZHat;
        cdouble K[3];
        cdouble N[3];
        PSD->GetEntries(Offset+np, "5:7",  K);
        PSD->GetEntries(Offset+np, "9:11", N);
        double IDNPF = 0.5*real( conj(K[0]) * (nHat[1]*N[2] - nHat[2]*N[1])
                                +conj(K[1]) * (nHat[2]*N[0] - nHat[0]*N[2])
                                +conj(K[2]) * (nHat[0]*N[1] - nHat[1]*N[0]));
        PSD->SetEntry(Offset + np, 12, IDNPF);
      };
   };

  return PSD;
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
HMatrix *RWGGeometry::GetPanelSourceDensities(cdouble Omega, HVector *KN, HMatrix *PSD)
{ return GetPanelSourceDensities(Omega, 0, KN, PSD); }

} // namespace scuff
