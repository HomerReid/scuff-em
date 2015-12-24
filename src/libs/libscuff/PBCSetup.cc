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
 * PBCSupport.cc -- a collection of various support routines used to 
 *               -- implement periodic boundary conditions in scuff-EM
 *
 * homer reid    -- 8/2011 -- 8/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libhrutil.h>
#include <libMDInterp.h> 
#include <libscuff.h>

namespace scuff{

/***************************************************************/
/* return 1 if the point with cartesian coordinates X lies on  */
/* the line connecting the origin to point L.                  */
/* Note that this routine is only looking at the first two     */
/* cartesian components of X; the z component is arbitrary.    */
/***************************************************************/
#if 0
static int PointOnLine(double *X, double *L)
{ 
   return ((float)(L[1]*X[0])) == ((float)(L[0]*X[1]));
} 
#endif

/***************************************************************/
/* given a single exterior edge on an RWGSurface, look for     */
/* a partner of this edge -- that is, another exterior edge    */
/* that is a translate of the original edge through a lattice  */
/* basis vector (LBV).                                         */
/*                                                             */
/* if a partner edge is found, its index within S's            */
/* ExteriorEdges array is returned. otherwise, -1 is returned. */
/*                                                             */
/* if a partner edge is found, then NumStraddlers[d] is        */
/* incremented, where d is the index of the LBV, and V[0..2]   */
/* is filled in with the cartesian coordinates of the new      */
/* vertex that must be added to turn the half-RWG basis        */
/* function associated with edge #nei into a full RWG basis    */
/* function.                                                   */
/*                                                             */
/* if a partner edge is found, then on return *pWhichBV is set */
/* to d.                                                       */
/***************************************************************/
static int FindPartnerEdge(RWGSurface *S, int nei, HMatrix *LBasis,
			   int NumStraddlers[MAXLDIM],
			   int *pWhichBV, double *V)
{
  RWGEdge *E = S->ExteriorEdges[nei];
  double *V1 = S->Vertices + 3*(E->iV1);
  double *V2 = S->Vertices + 3*(E->iV2);
  const double tolvc = S->tolVecClose;
  
  /*--------------------------------------------------------------*/
  /* Look for an exterior edge that is a translate through a      */
  /* lattice basis vector of the given edge. This could be more   */
  /* efficient.                                                   */
  /*--------------------------------------------------------------*/
  int LDim=LBasis->NC;
  for(int nd=0; nd<LDim; nd++)
   { 
     double V1T[3], V2T[3]; // 'V12, translated'
     for(int nc=0; nc<3; nc++)
      { V1T[nc] = V1[nc] + LBasis->GetEntryD(nc,nd);
        V2T[nc] = V2[nc] + LBasis->GetEntryD(nc,nd);
      };

     for(int neip=0; neip<S->NumExteriorEdges; neip++)
      { 
        if (S->ExteriorEdges[neip]==0) 
         continue;

        double *V1P = S->Vertices + 3*(S->ExteriorEdges[neip]->iV1);
        double *V2P = S->Vertices + 3*(S->ExteriorEdges[neip]->iV2);
        if (   (VecClose(V1T, V1P, tolvc) && VecClose(V2T, V2P, tolvc))
            || (VecClose(V1T, V2P, tolvc) && VecClose(V2T, V1P, tolvc))
           )
         { 
           /*--------------------------------------------------------------*/
           /*- found a translate of the edge in question.                  */
           /*--------------------------------------------------------------*/
           memcpy(V, S->Vertices + 3*(S->ExteriorEdges[neip]->iQP), 3*sizeof(double));
           for(int nc=0; nc<3; nc++)
            V[nc] -= LBasis->GetEntryD(nc,nd);
           if (NumStraddlers) NumStraddlers[nd]++;
           if (pWhichBV) *pWhichBV=nd;
           return neip;
         };
      };
   };

  return -1; // no partner found

}

/*--------------------------------------------------------------*/
/* A 'straddler' is an external edge that has a partner ('image')*/
/* edge translated by a lattice vector.                         */
/* This routine goes through the list of exterior edges         */
/* in the given RWGSurface and identifies all straddlers. For   */
/* each straddler, the routine adds one new panel and one new   */
/* vertex to the Surface. The new panel is the translate of the */
/* panel to which the partner edge is attached, while the new   */
/* vertex is the vertex opposite the partner edge within that   */
/* panel. Once we have added this new vertex and panel, we can  */
/* promote the exterior edge to an interior edge.               */
/*                                                              */
/* Parameters:                                                  */
/*  S points to the RWGSurface in question; its contents are    */
/*  modified if any straddlers were detected.                   */
/*                                                              */
/*  On return, NumStraddlers[i] is the number of straddlers     */
/*  detected on the unit-cell boundary normal to the ith        */
/*  lattice basis vector.                                       */
/*--------------------------------------------------------------*/
#define CHUNK 100
void RWGSurface::AddStraddlers(HMatrix *LBasis, HMatrix *RLBasis,
                               int NumStraddlers[MAXLDIM])
{ 
  int NumNew=0, NumAllocated=0;
  double V[3], *NewVertices=0;
  RWGPanel *P, **NewPanels=0;
  RWGEdge *E, **NewEdges=0;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  memset(NumStraddlers, 0, MAXLDIM*sizeof(int));
  for(int nei=0; nei<NumExteriorEdges; nei++)
   { 
      if ( ExteriorEdges[nei]==0 )
       continue;

      // see if this edge is a straddler, i.e. it lies on a face
      // of the unit cell and it has an image (a translated
      // version of itself) on the opposite side of the unit cell
      int WhichBV=0;
      int neip=FindPartnerEdge(this, nei, LBasis,
                               NumStraddlers, &WhichBV, V);

      // if so, add a new vertex, panel, and interior edge to the RWGSurface.
      if (neip!=-1)
       { 
          // if necessary, expand local arrays of new vertices, panels, and edges 
          if( NumAllocated == NumNew )
           { NumAllocated+=CHUNK;
             NewVertices = (double *)reallocEC(NewVertices, 3*NumAllocated*sizeof(double));
             NewPanels   = (RWGPanel **)reallocEC(NewPanels, NumAllocated*sizeof(RWGPanel *));
             NewEdges    = (RWGEdge **)reallocEC(NewEdges, NumAllocated*sizeof(RWGEdge *));
             PhasedBFCs  = (int *)reallocEC(PhasedBFCs, 3*NumAllocated*sizeof(int ));
           };

          // add a new vertex
          memcpy( NewVertices + 3*NumNew, V, 3*sizeof(double));
          RMax[0] = fmax(RMax[0], V[0] );
          RMax[1] = fmax(RMax[1], V[1] );
          RMax[2] = fmax(RMax[2], V[2] );
          RMin[0] = fmin(RMin[0], V[0] );
          RMin[1] = fmin(RMin[1], V[1] );
          RMin[2] = fmin(RMin[2], V[2] );
 
          // add a new edge. actually, we simply appropriate the 
          // existing RWGEdge structure for edge #nei, since we 
          // will be removing it from object S's set of exterior 
          // edges anyway. note that the iQP, iV1, iV2, iPPanel, PIndex, 
          // Centroid, and Length fields of this structure will be already 
          // correctly initialized. 
          E=ExteriorEdges[nei];
          E->iQM     = NumVertices + NumNew;
          E->iMPanel = NumPanels + NumNew;
          E->Index   = NumEdges + NumNew;
          E->MIndex  = 2; // because below we hard-code P->VI[2] = E->iQM;
          E->Radius=VecDistance(E->Centroid, Vertices+3*E->iQP);
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,V));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,Vertices+3*E->iV1));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,Vertices+3*E->iV2));
          NewEdges[NumNew] = E;

          // remove edge #nei from the list of exterior edges 
          ExteriorEdges[nei]=0;

          // add a new panel. Note that we can't call InitRWGPanel yet 
          // because the new vertex has not yet been added to the Vertices 
          // array; this happens later, below.
          P=(RWGPanel *)mallocEC(sizeof(RWGPanel));
          P->VI[0] = E->iV1;
          P->VI[1] = E->iV2;
          P->VI[2] = E->iQM;
          P->EI[0] = P->EI[1] = P->EI[2] = -1;
          P->Index = NumPanels + NumNew;
          NewPanels[NumNew] = P;

          // update our list of 'phased basis-function contributions'
          // to make a note of the fact that the newly-added edge
          // will make a phased contribution to a panel on the
          // opposite side of the unit cell
          //PhasedBFCs[NumNew].WhichPanel = Panels[neip]->iPPanel;
          //PhasedBFCs[NumNew].WhichEdge  = NumEdges + NumNew;
          //PhasedBFCs[NumNew].WhichBV;
          PhasedBFCs[3*NumNew + 0] = ExteriorEdges[neip]->iPPanel;
          PhasedBFCs[3*NumNew + 1] = NumEdges + NumNew;
          PhasedBFCs[3*NumNew + 2] = WhichBV;
  
          NumNew++;

       };
     
   }; // for(nei=0; nei<NumExteriorEdges; nei++)

  TotalStraddlers=NumNew;

  if (NumNew==0)
   return;

  /*--------------------------------------------------------------*/
  /*- replace internal arrays of Vertices, Panels, Edges, and     */
  /*- ExteriorEdges with new arrays that include the straddlers   */
  /*--------------------------------------------------------------*/
  Vertices = (double *)reallocEC( Vertices, 3*(NumVertices+NumNew) * sizeof(double));
  memcpy( &(Vertices[3*NumVertices]) , NewVertices, 3*NumNew*sizeof(double));
  NumVertices+=NumNew;

  Edges    = (RWGEdge **)reallocEC( Edges, (NumEdges+NumNew) * sizeof(RWGEdge *));
  memcpy( &(Edges[NumEdges]), NewEdges, NumNew*sizeof(RWGEdge *));
  NumEdges+=NumNew;

  Panels    = (RWGPanel **)reallocEC( Panels, (NumPanels+NumNew) * sizeof(RWGPanel *));
  memcpy( &(Panels[NumPanels]), NewPanels, NumNew*sizeof(RWGPanel *));
  NumPanels+=NumNew;
  for(int np=NumPanels-NumNew; np<NumPanels; np++)
   InitRWGPanel(Panels[np], Vertices);

  /*--------------------------------------------------------------*/
  /*- defragment the ExteriorEdges array -------------------------*/
  /*--------------------------------------------------------------*/
  int NewNumExteriorEdges=NumExteriorEdges-NumNew;
  RWGEdge **NewExteriorEdges=(RWGEdge **)mallocEC(NewNumExteriorEdges*sizeof(RWGEdge *));
  for(int nei=0, neip=0; nei<NumExteriorEdges; nei++)
   if (ExteriorEdges[nei]!=0)
    { NewExteriorEdges[neip]=ExteriorEdges[nei];
      NewExteriorEdges[neip]->Index=neip;
      neip++;
    };

  free(ExteriorEdges);
  ExteriorEdges=NewExteriorEdges;
  NumExteriorEdges=NewNumExteriorEdges;
 
  NumTotalEdges=NumEdges + NumExteriorEdges;
  NumBFs = ( IsPEC ? NumEdges : 2*NumEdges );

  /*------------------------------------------------------------*/
  /*- repeat the calculation to initialize the EI arrays in each*/
  /*- RWGPanel structure (this was originally done back in      */
  /*- the RWGSurface class constructor).                        */
  /*- Note: EI[i] = index of panel edge #i within the Edges[]   */
  /*- list for the parent RWGSurface. Panel edge #i is the      */
  /*- edge opposite vertex #i.                                  */
  /*------------------------------------------------------------*/
  for(int ne=0; ne<NumEdges; ne++)
   { E = Edges[ne];
     Panels[ E->iPPanel ] -> EI[ E->PIndex ] = ne;
     if ( E->iMPanel>=0 )
      Panels[ E->iMPanel ] -> EI[ E->MIndex ] = ne;
   };
  for(int ne=0; ne<NumExteriorEdges; ne++)
   { E=ExteriorEdges[ne];
     Panels[E->iPPanel]->EI[E->PIndex] = -(ne+1);
   };

  Log(" Surface %s: \n",Label);
  int LDim=LBasis->NC;
  for(int nd=0; nd<LDim; nd++)
   Log("  %i straddlers normal to lattice vector #%i",NumStraddlers[nd],nd);
  Log("  %i total straddlers ",TotalStraddlers);

}

/***************************************************************/
/* This is a helper function in the RWGGeometry class that     */
/* initializes internal data fields needed for working with    */
/* periodic boundary conditions.                               */
/***************************************************************/
void RWGGeometry::InitPBCData()
{
  /*--------------------------------------------------------------*/
  /* Step 1: Addition of 'straddlers.'                            */
  /* For all surfaces, we detect exterior triangle edges that lie */
  /* on the unit-cell boundary ('straddlers'), and we add new     */  
  /* panels and vertices to the surface to allow those edges to   */
  /* be promoted from exterior to interior edges.                 */
  /*                                                              */
  /* Note: NumStraddlers is an array stored within RWGGeometry    */
  /* that tabulates how many straddlers each surface has in each  */
  /* possible direction.                                          */
  /* NumStraddlers[MAXLDIM*ns + nd] = number of stradders for     */
  /*                                  surface #ns that straddle   */
  /*                                  the unit cell boundary      */
  /*                                  normal to basis vector #nd  */
  /*                                                              */
  /* Also, RegionIsExtended[MAXLDIM*nr + i] = true if region nr   */
  /* is extended (as opposed to compact) in the direction of      */
  /* lattice basis vector i; =false otherwise.                    */
  /*--------------------------------------------------------------*/
  TotalBFs=TotalPanels=0;
  for(int nd=0; nd<LDim; nd++)
   { NumStraddlers[nd]       = (int *)mallocEC(NumSurfaces*sizeof(int));
     RegionIsExtended[nd]    = (bool *)mallocEC(NumRegions*sizeof(bool));
     RegionIsExtended[nd][0] = true; // exterior medium is always extended
   };
  #define ONEMEG 1048576
  Log(" Mem before straddlers: %lu",GetMemoryUsage()/ONEMEG);
  for(int ns=0; ns<NumSurfaces; ns++)
   { 
     RWGSurface *S=Surfaces[ns];
     int NumStraddlersThisSurface[MAXLDIM];
     S->AddStraddlers(LBasis, RLBasis, NumStraddlersThisSurface);

     int nr1=S->RegionIndices[0];
     int nr2=S->RegionIndices[1];
     for(int nd=0; nd<LDim; nd++)
      { 
        NumStraddlers[nd][ns]=NumStraddlersThisSurface[nd];
        if ( NumStraddlers[nd][ns] > 0 )
         { RegionIsExtended[nd][nr1]=true;
           if (nr2!=-1) RegionIsExtended[nd][nr2]=true;
         };
      };
     TotalBFs+=S->NumBFs;
     TotalPanels+=S->NumPanels;

     int TotalStraddlers=0;
     for(int nd=0; nd<LDim; nd++)
      TotalStraddlers += NumStraddlersThisSurface[nd];

     S->IsClosed = (TotalStraddlers == S->NumExteriorEdges);
     Log("Surface %s: {Straddlers, ExteriorEdges}={%i,%i} {%s}",
          S->Label,TotalStraddlers,S->NumExteriorEdges,
         (S->IsClosed ? "closed" : "open"));

   };

}

/***************************************************************/
/* for periodic geometries, decompose X = L + XBar where       */
/* L is a lattice vector and XBar lies within the lattice      */
/* unit cell.                                                  */
/*                                                             */
/* On return, nVector={n1,n2,n3} where L=n1*L1 + n2*L2 + n3*L3 */
/* where {L1,L2,L3} are the lattice vectors.                   */
/*                                                             */
/* if WignerSeitz =true, translate X into the Wigner-Seitz cell*/
/* instead of the unit cell.                                   */
/*                                                             */
/* algorithm: suppose X = \sum_i (n_i + \nu_i) L_i where       */
/*            n_i = integer and 0 <= \nu_i < 1.                */
/*                                                             */
/*            then since \Gamma_i \cdot L_j = 2\pi \delta_{ij} */
/*            (where \Gamma_i = reciprocal lattice vectors)    */
/*            I have simply                                    */
/*            (n_i + \nu_i) = (\Gamma_i \cdot X) / (2\pi)      */
/*                                                             */
/* and then L=\sum n_i L_i, XBar = \sum \nu_i L_i.             */
/***************************************************************/
void RWGGeometry::GetUnitCellRepresentative(const double X[3],
                                            double XBar[3],
                                            double LVector[3],
                                            int nVector[MAXLDIM],
                                            bool WignerSeitz)
{
  XBar[0]=XBar[1]=XBar[2]=LVector[0]=LVector[1]=LVector[2]=0.0;
  for(int nd=0; nd<LDim; nd++)
   {
     double npNu=0.0;
     for(int nc=0; nc<3; nc++)
      npNu += RLBasis->GetEntryD(nc,nd)*X[nc];
     npNu /= 2.0*M_PI;

     double n  = WignerSeitz ? lround(npNu) : floor(npNu);
     double Nu = npNu - n;

     nVector[nd] = (int)n;
     for(int nc=0; nc<3; nc++)
      { XBar[nc]    += Nu*LBasis->GetEntryD(nc,nd);
        LVector[nc] +=  n*LBasis->GetEntryD(nc,nd);
      };
   };
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/* delete me after checking that the previous routine works     */
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void OldGetUnitCellRepresentative(RWGGeometry *G, const double X[3],
                                  double XBar[3], double LVector[3],
                                  int nVector[MAXLDIM],
                                  bool WignerSeitz)
{
  XBar[0] = X[0];
  XBar[1] = X[1];
  XBar[2] = X[2];
  
  HMatrix *LBasis=G->LBasis;
  int LDim=LBasis->NC;
  if (LDim>=1)
   { double L = LBasis->GetEntryD(0,0);
     double N = WignerSeitz ? lround(X[0]/L) : floor(X[0]/L);
     XBar[0] -= N*L;
     nVector[0] = N;
     LVector[0] = N*L;
   };
  if (LDim>=2)
   { double L = LBasis->GetEntryD(1,1);
     double N = WignerSeitz ? lround(X[1]/L) : floor(X[1]/L);
     XBar[1] -= N*L;
     nVector[1] = N;
     LVector[1] = N*L;
   };

}

void RWGGeometry::GetUnitCellRepresentative(const double X[3],
                                            double XBar[3],
                                            double LVector[3],
                                            bool WignerSeitz)
{ int NVector[3];
  GetUnitCellRepresentative(X, XBar, LVector, NVector, WignerSeitz);
}

void RWGGeometry::GetUnitCellRepresentative(const double X[3],
                                            double XBar[3],
                                            bool WignerSeitz)
{ int NVector[3];
  double LVector[3];
  GetUnitCellRepresentative(X, XBar, LVector, NVector, WignerSeitz);
}

} // namespace scuff
