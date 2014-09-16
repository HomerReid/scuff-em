#define ONEMEG 1048576
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
static int PointOnLine(double *X, double *L)
{ 
   return ((float)(L[1]*X[0])) == ((float)(L[0]*X[1]));
} 

/***************************************************************/
/* given a single exterior edge on an RWGSurface, look for     */
/* a partner of this edge -- that is, another exterior edge    */
/* that is a translate of the original edge through either     */
/* LBV[0] or LBV[1].                                           */
/*                                                             */
/* (note: LBV = 'lattice basis vector.')                       */
/*                                                             */
/* if a partner edge is found, its index within S's            */
/* ExteriorEdges array is returned. otherwise, -1 is returned. */
/*                                                             */
/* if a partner edge is found, then NumStraddlers[i] is        */
/* incremented (where i=0 or 1 depending on whether the        */
/* translation happened through LBV[0] or LBV[1]) and V[0..2]  */
/* is filled in with the cartesian coordinates of the new      */
/* vertex that must be added to turn the half-RWG basis        */
/* function associated with edge #nei into a full RWG basis    */
/* function.                                                   */
/*                                                             */
/* if a partner edge is found, then on return *pWhichBV is set */
/* to 0 or 1 depending on whether the translation happened     */
/* through LBV[0] or LBV[1].                                   */
/*                                                             */
/* Update 20140721: For 1D-periodic geometries I don't think   */
/* it's necessary to establish whether or not the edge lives   */
/* on the unit-cell boundary; indeed, without the second       */
/* lattice vector I'm not even sure the unit-cell boundary     */
/* is well-defined.                                            */
/***************************************************************/
static int FindPartnerEdge(RWGSurface *S, int nei,
			   double LBV[MAXLDIM][2],
                           double LBVi[MAXLDIM][2],
                           int LDim,
			   int NumStraddlers[MAXLDIM],
			   int *pWhichBV,
                           double *V)
{
  RWGEdge *E = S->ExteriorEdges[nei];
  double *V1 = S->Vertices + 3*(E->iV1);
  double *V2 = S->Vertices + 3*(E->iV2);
  const double tolvc = S->tolVecClose;
  
  /*--------------------------------------------------------------*/
  /*- determine whether or not the edge lies on the unit cell     */
  /*- boundary, and if so which face of that boundary it lies on. */
  /*-                                                             */
  /*- After this code snippet,                                    */
  /*-  LTranslate = the basis vector through which we expect to   */
  /*-               translate the given edge to find an image     */
  /*-  *pWhich    = the index of LTranslate within LBV (0 or 1)   */
  /*--------------------------------------------------------------*/
  double *LTranslate=LBV[0];
  int WhichBV=0;
  if (LDim==1)
   {
     WhichBV=*pWhichBV=0;
     LTranslate = LBV[0];
   }
  else if (LDim==2)
   { 
     // convert V1 and V2 to lattice basis
     double V1L[2], V2L[2];
     V1L[0] = LBVi[0][0]*V1[0] + LBVi[0][1]*V1[1];
     V1L[1] = LBVi[1][0]*V1[0] + LBVi[1][1]*V1[1];
     V2L[0] = LBVi[0][0]*V2[0] + LBVi[0][1]*V2[1];
     V2L[1] = LBVi[1][0]*V2[0] + LBVi[1][1]*V2[1];
     if ( fabs(V1L[0]) < tolvc && fabs(V2L[0]) < tolvc )
      { WhichBV=*pWhichBV = 0;
        LTranslate=LBV[0];
      }
     else if ( fabs(V1L[1]) < tolvc && fabs(V2L[1]) < tolvc )
      { WhichBV=*pWhichBV = 1;
        LTranslate=LBV[1];
      }
     else
      return -1; // edge does not lie on unit cell boundary
   }
  else
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);

  /*--------------------------------------------------------------*/
  /* Look for an exterior edge that is the translate through      */
  /* LTranslate of the present exterior edge.                     */
  /*--------------------------------------------------------------*/
  double V1T[3], V2T[3]; // 'V12, translated'
  V1T[0] = V1[0] + LTranslate[0]; V1T[1] = V1[1] + LTranslate[1]; V1T[2] = V1[2];
  V2T[0] = V2[0] + LTranslate[0]; V2T[1] = V2[1] + LTranslate[1]; V2T[2] = V2[2];
  double *V1P, *V2P; // 'V1,V2, primed'
  for(int neip=0; neip<S->NumExteriorEdges; neip++)
   { 
     if (S->ExteriorEdges[neip]==0) 
      continue;

     V1P = S->Vertices + 3*(S->ExteriorEdges[neip]->iV1);
     V2P = S->Vertices + 3*(S->ExteriorEdges[neip]->iV2);
     if (   (VecClose(V1T, V1P, tolvc) && VecClose(V2T, V2P, tolvc))
         || (VecClose(V1T, V2P, tolvc) && VecClose(V2T, V1P, tolvc))
        )
      { 
        /*--------------------------------------------------------------*/
        /*- found a translate of the edge in question.                  */
        /*--------------------------------------------------------------*/
        memcpy(V, S->Vertices + 3*(S->ExteriorEdges[neip]->iQP), 3*sizeof(double));
        V[0] -= LTranslate[0]; 
        V[1] -= LTranslate[1]; 
        if (NumStraddlers) NumStraddlers[WhichBV]++;
        return neip;
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (LDim==2)
   ErrExit("exterior edge %i of object %s has no image "
           "on the opposite side of the unit cell", nei, S->Label);

  return -1; // i am indecisive about how to handle this error

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
/*  LBV[i][j] is the jth cartesian component of the ith lattice */
/*  basis vector.                                               */
/*  LDim is the length of the first dimension of                */
/*  the LBV array.                                              */
/*                                                              */
/*  On return, NumStraddlers[i] is the number of straddlers     */
/*  detected on the unit-cell boundary normal to LBV[i].        */
/*--------------------------------------------------------------*/
#define CHUNK 100
void RWGSurface::AddStraddlers(double LBV[MAXLDIM][2],
                               int LDim,
                               int NumStraddlers[MAXLDIM])
{ 
  int NumNew=0, NumAllocated=0;
  double V[3], *NewVertices=0;
  RWGPanel *P, **NewPanels=0;
  RWGEdge *E, **NewEdges=0;

  memset(NumStraddlers, 0, MAXLDIM*sizeof(int));

  /***************************************************************/
  /* find the LBVi vectors, which satisfy dot(LBVi[k],LBV[j])    */
  /* = delta(j,k), and are used to convert points to the lattice */
  /* basis.                                                      */
  /***************************************************************/
  double LBVi[2][2] = {{0.0,0.0}, {0.0,0.0}};
  if (LDim==1) 
   {
     double L1Norm2 = LBV[0][0]*LBV[0][0] + LBV[0][1]*LBV[0][1];
     LBVi[0][0] = LBV[0][0] / L1Norm2;
     LBVi[0][1] = LBV[0][1] / L1Norm2;
   }
  else if (LDim==2) 
   {
     // LBVi[0] = (LBV[0] x LBV[1]) x LBV[1] / (LBV[0] * (... x ... x ...))
     // LBVi[1] = (LBV[0] x LBV[1]) x LBV[0] / (LBV[1] * (... x ... x ...))
     //double perp[3];
     //VecCross(LBV[0], LBV[1], perp);
     double L1Norm2 = LBV[0][0]*LBV[0][0] + LBV[0][1]*LBV[0][1];
     double L2Norm2 = LBV[1][0]*LBV[1][0] + LBV[1][1]*LBV[1][1];
     double zPerp = LBV[0][0]*LBV[1][1] - LBV[0][1]*LBV[1][0];
     if ( fabs(zPerp) < 1e-8 * sqrt(L1Norm2*L2Norm2) )
      ErrExit("Lattice vectors close to parallel.");
     for (int i=0; i<2; ++i)
      { 
        //VecCross(perp, LBV[1-i], LBVi[i]);
        LBVi[i][0] = -zPerp * LBV[1-i][1];
        LBVi[i][1] = +zPerp * LBV[1-i][0];

        //VecScale(LBVi[i], 1 / VecDot(LBV[i], LBVi[i]));
        double DotProd = LBV[i][0]*LBVi[i][0] + LBV[i][1]*LBVi[i][1];
        LBVi[i][0] /= DotProd;
        LBVi[i][1] /= DotProd;
      }
   }
  else // (LDim> 2)
   ErrExit("%d lattice vectors unsupported\n", LDim);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nei=0; nei<NumExteriorEdges; nei++)
   { 
      if ( ExteriorEdges[nei]==0 )
       continue;

      // see if this edge is a straddler, i.e. it lies on a face
      // of the unit cell and it has an image (a translated 
      // version of itself) on the opposite side of the unit cell
      int WhichBV=0;
      int neip=FindPartnerEdge(this, nei, LBV, LBVi, 
                               LDim, NumStraddlers, &WhichBV, V);

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
          E->iQM  = NumVertices + NumNew;
          E->iMPanel = NumPanels + NumNew;
          E->Index = NumEdges + NumNew;
          E->MIndex = 2; // because below we hard-code P->VI[2] = E->iQM;
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

  if (NumStraddlers)
   Log(" Detected (%i,%i) straddlers for surface %s", NumStraddlers[0],NumStraddlers[1],Label);

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
  /* NumStraddlers[MAXLDIM*ns + i] = number of stradders for   */
  /*                                    surface #ns that straddle */
  /*                                    the unit cell boundary    */
  /*                                    normal to basis vector #i */
  /*                                                              */
  /* Also, RegionIsExtended[MAXLDIM*nr + i] = true if region nr*/
  /* is extended (as opposed to compact) in the direction of      */
  /* lattice basis vector i; =false otherwise.                    */
  /*--------------------------------------------------------------*/
  TotalBFs=TotalPanels=0;
  for(int nd=0; nd<LDim; nd++)
   { NumStraddlers[nd]       = (int *)mallocEC(NumSurfaces*sizeof(int));
     RegionIsExtended[nd]    = (bool *)mallocEC(NumRegions*sizeof(bool));
     RegionIsExtended[nd][0] = true; // exterior medium is always extended
   };
  Log(" Mem before straddlers: %lu",GetMemoryUsage()/ONEMEG);
  for(int ns=0; ns<NumSurfaces; ns++)
   { 
     RWGSurface *S=Surfaces[ns];
     int NumStraddlersThisSurface[2];
     S->AddStraddlers(LBasis, LDim, NumStraddlersThisSurface);

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

   };

}

} // namespace scuff
