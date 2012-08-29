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
/* do something about me please ********************************/
/***************************************************************/
int RWGGeometry::PBCCubatureOrder=4;
double RWGGeometry::DeltaInterp=0.05;

/***************************************************************/
/* get the maximum and minimum cartesian coordinates obtained  */
/* by points on an RWGSurface. (we do this by looking at the    */
/* panel vertices, which suffices because the panels are       */
/* flat and thus convex)                                       */
/***************************************************************/
void GetSurfaceExtents(RWGSurface *S, double RMax[3], double RMin[3])
{ 
  memcpy(RMax, S->Vertices + 0, 3*sizeof(double));
  memcpy(RMin, S->Vertices + 0, 3*sizeof(double));
  for(int nv=1; nv<S->NumVertices; nv++)
   { 
     RMax[0] = fmax(RMax[0], S->Vertices[3*nv + 0]);
     RMin[0] = fmin(RMin[0], S->Vertices[3*nv + 0]);
     RMax[1] = fmax(RMax[1], S->Vertices[3*nv + 1]);
     RMin[1] = fmin(RMin[1], S->Vertices[3*nv + 1]);
     RMax[2] = fmax(RMax[2], S->Vertices[3*nv + 2]);
     RMin[2] = fmin(RMin[2], S->Vertices[3*nv + 2]);
   };
}

/***************************************************************/
/* return 1 if the point with cartesian coordinates X lies on  */
/* the line connecting the origin to point L.                  */
/* Note that this routine is only looking at the first two     */
/* cartesian components of X; the z component is arbitrary.    */
/***************************************************************/
static int PointOnLine(double *X, double *L)
{ 
   return ((float)(L[1]*X[0])) == ((float)(L[0]*X[1])) ? 1 : 0;
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
/* incremented (where i=0 or 1 depending on which unit cell    */
/* boundary the original edge lay on) and V[0..2] is filled in */
/* with the cartesian coordinates of the new vertex that must  */
/* be added to turn the half-RWG basis function associated with*/
/* edge #nei into a full RWG basis function.                   */
/***************************************************************/
static int FindPartnerEdge(RWGSurface *S, int nei, double *LBV[2], 
                           int NumLatticeVectors, int NumStraddlers[2], double *V)
{
  if (NumLatticeVectors!=2)
   ErrExit("%s:%i: NumLatticeVectors != 2 not yet supported",__FILE__,__LINE__);

  RWGEdge *E = S->ExteriorEdges[nei];
  double *V1 = S->Vertices + 3*(E->iV1);
  double *V2 = S->Vertices + 3*(E->iV2);
  
  /*--------------------------------------------------------------*/
  /*- determine whether or not the edge lies on the unit cell     */
  /*- boundary, and if so which face of that boundary it lies on. */
  /*--------------------------------------------------------------*/
  int WhichBV;
  double *ThisBV, *OtherBV;
  if ( PointOnLine(V1, LBV[0]) && PointOnLine(V2, LBV[0]) )
   { WhichBV=0;
     ThisBV=LBV[0]; 
     OtherBV=LBV[1];
   }
  else if ( PointOnLine(V1, LBV[1]) && PointOnLine(V2, LBV[1]) )
   { WhichBV=1;
     ThisBV=LBV[1]; 
     OtherBV=LBV[0];
   }
  else
   return -1; // edge does not lie on unit cell boundary 

  /*--------------------------------------------------------------*/
  /* Look for an exterior edge that is the translate through      */
  /* OtherBV of the present exterior edge.                        */
  /*--------------------------------------------------------------*/
  double V1T[3], V2T[3]; // 'V12, translated'
  V1T[0] = V1[0] + OtherBV[0]; V1T[1] = V1[1] + OtherBV[1]; V1T[2] = V1[2];
  V2T[0] = V2[0] + OtherBV[0]; V2T[1] = V2[1] + OtherBV[1]; V2T[2] = V2[2];
  double *QP, *V1P, *V2P; // 'Q,V1,V2, primed'
  for(int neip=0; neip<S->NumExteriorEdges; neip++)
   { 
     if (S->ExteriorEdges[neip]==0) 
      continue;

     V1P = S->Vertices + 3*(S->ExteriorEdges[neip]->iV1);
     V2P = S->Vertices + 3*(S->ExteriorEdges[neip]->iV2);
     if (   (VecEqualFloat(V1T, V1P) && VecEqualFloat(V2T, V2P))
         || (VecEqualFloat(V1T, V2P) && VecEqualFloat(V2T, V1P))
        )
      { 
        /*--------------------------------------------------------------*/
        /*- found a translate of the edge in question.                  */
        /*--------------------------------------------------------------*/
        memcpy(V, S->Vertices + 3*(S->ExteriorEdges[neip]->iQP), 3*sizeof(double));
        V[0] -= OtherBV[0]; 
        V[1] -= OtherBV[1]; 
        if (NumStraddlers) NumStraddlers[WhichBV]++;
        return neip;
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  ErrExit("exterior edge %i of object %s has no image "
          "on the opposite side of the unit cell",nei, S->Label);

  return -1; // i am indecisive about how to handle this error

}

/*--------------------------------------------------------------*/
/* A 'straddler' is a panel edge that lies on a face of the unit*/
/* cell and has a partner edge on the opposite face of the unit */
/* cell. This routine goes through the list of exterior edges   */
/* in the given RWGSurface and identifies all straddlers. For   */
/* each straddler, the routine adds one new panel and one new   */
/* vertex to the Surface. The new panel is the translate of the  */
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
/*  NumLatticeVectors is the length of the first dimension of   */
/*  the LBV array.                                              */
/*                                                              */
/*  If NumStraddlers is non-NULL, then on return                */ 
/*  NumStraddlers[i] is the number of straddlers detected on    */
/*  the unit-cell boundary normal to LBV[i].                    */
/*--------------------------------------------------------------*/
#define CHUNK 100
void AddStraddlers(RWGSurface *S, double **LBV, int NumLatticeVectors, int *NumStraddlers)
{ 
  int NumNew=0, NumAllocated=0;
  double V[3], *NewVertices=0;
  RWGPanel *P, **NewPanels=0;
  RWGEdge *E, **NewEdges=0;

  if (NumStraddlers)
   memset(NumStraddlers, 0, 2*sizeof(int));

  int nei, neip;
  for(nei=0; nei<S->NumExteriorEdges; nei++)
   { 
      if ( S->ExteriorEdges[nei]==0 )
       continue;

      // see if this edge is a straddler, i.e. it lies on a face
      // of the unit cell and it has a partner (a translated 
      // version of itself) on the opposite side of the unit cell
      neip=FindPartnerEdge(S, nei, LBV, NumLatticeVectors, NumStraddlers, V);

      // if so, add a new vertex, panel, and interior edge to the RWGSurface.
      if (neip!=-1)
       { 
          // if necessary, expand local arrays of new vertices, panels, and edges 
          if( NumAllocated == NumNew )
           { NumAllocated+=CHUNK;
             NewVertices = (double *)reallocEC(NewVertices, 3*NumAllocated*sizeof(double));
             NewPanels   = (RWGPanel **)reallocEC(NewPanels, NumAllocated*sizeof(RWGPanel *));
             NewEdges    = (RWGEdge **)reallocEC(NewEdges, NumAllocated*sizeof(RWGEdge *));
           };

          // add a new vertex
          memcpy( NewVertices + 3*NumNew, V, 3*sizeof(double));
 
          // add a new edge. actually, we simply appropriate the 
          // existing RWGEdge structure for edge #nei, since we 
          // will be removing it from object S's set of exterior 
          // edges anyway. note that the iQP, iV1, iV2, iPPanel, PIndex, 
          // Centroid, and Length fields of this structure will be already 
          // correctly initialized. 
          E=S->ExteriorEdges[nei];
          E->iQM  = S->NumVertices + NumNew;
          E->iMPanel = S->NumPanels + NumNew;
          E->Index = S->NumEdges + NumNew;
          E->MIndex = 2; // because below we hard-code P->VI[2] = E->iQM;
          E->Radius=VecDistance(E->Centroid, S->Vertices+3*E->iQP);
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,S->Vertices+3*E->iQM));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,S->Vertices+3*E->iV1));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,S->Vertices+3*E->iV2));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,S->Vertices+3*E->iV2));
          NewEdges[NumNew] = E;

          // what we just did was to combine two exterior edges 
          // into a single interior edge, so we now remove those
          // two exterior edges from S's list of exterior edges
          free(S->ExteriorEdges[neip]);
          S->ExteriorEdges[nei]=S->ExteriorEdges[neip]=0;

          // add a new panel. Note that we can't call InitRWGPanel yet 
          // because the new vertex has not yet been added to the Vertices 
          // array; this happens later, below.
          P=(RWGPanel *)mallocEC(sizeof(RWGPanel));
          P->VI[0] = E->iV1;
          P->VI[1] = E->iV2;
          P->VI[2] = E->iQM;
          P->Index = S->NumPanels + NumNew;
          NewPanels[NumNew] = P;
  
          NumNew++;

       };
     
   }; // for(nei=0; nei<S->NumExteriorEdges; nei++)

  if (NumNew==0)
   return;

  /*--------------------------------------------------------------*/
  /*- replace S's arrays of Vertices, Panels, Edges, and          */
  /*- ExteriorEdges with new arrays that include the straddlers   */
  /*--------------------------------------------------------------*/
  double *Vertices        = S->Vertices;
  RWGPanel **Panels       = S->Panels;
  RWGEdge **Edges         = S->Edges;
  RWGEdge **ExteriorEdges = S->ExteriorEdges;
  int NumVertices         = S->NumVertices;
  int NumPanels           = S->NumPanels;
  int NumEdges            = S->NumEdges;
  int NumExteriorEdges    = S->NumExteriorEdges;

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
  int NewNumExteriorEdges=NumExteriorEdges-2*NumNew;
  RWGEdge **NewExteriorEdges=(RWGEdge **)mallocEC(NewNumExteriorEdges*sizeof(RWGEdge *));
  for(nei=neip=0; nei<NumExteriorEdges; nei++)
   if (ExteriorEdges[nei]!=0)
    { NewExteriorEdges[neip]=ExteriorEdges[nei];
      NewExteriorEdges[neip]->Index=neip;
      neip++;
    };
      
  free(ExteriorEdges);
  ExteriorEdges=NewExteriorEdges;
  NumExteriorEdges=NewNumExteriorEdges;
 
  S->Vertices         = Vertices;
  S->Panels           = Panels;
  S->Edges            = Edges;
  S->ExteriorEdges    = ExteriorEdges;
  S->NumVertices      = NumVertices;
  S->NumPanels        = NumPanels;
  S->NumEdges         = NumEdges;
  S->NumExteriorEdges = NumExteriorEdges;

  S->NumTotalEdges=NumEdges + NumExteriorEdges;
  S->NumBFs = ( S->IsPEC ? NumEdges : 2*NumEdges );

  if (NumStraddlers)
   Log(" Detected (%i,%i) straddlers for object %s", NumStraddlers[MAXLATTICE*ns+0],NumStraddlers[MAXLATTICE*ns+1],S->Label);

}

/***************************************************************/
/* This is a helper function in the RWGGeometry class that     */
/* initializes internal data fields needed for working with    */
/* periodic boundary conditions.                               */
/***************************************************************/
RWGGeometry::InitPBCData()
{
  /*--------------------------------------------------------------*/
  /* Step 1: Addition of 'straddlers.'                            */
  /* For all surfaces, we detect exterior triangle edges that lie */
  /* on the unit-cell boundary ('straddlers'), and we add new     */
  /* panels and vertices to the surface to allow those edges to   */
  /* be promoted from exterior to interior edges.                 */
  /*--------------------------------------------------------------*/
  TotalBFs=TotalPanels=0;
  for(int ns=0; ns<NumSurfaces; ns++)
   { AddStraddlers(Surfaces[ns], LBV, NumLatticeVectors, 0);
     TotalBFs+=S->NumBFs;
     TotalPanels+=S->NumPanels;
     if ( ns+1 < G->NumSurfaces )
      { G->BFIndexOffset[ns+1]=G->BFIndexOffset[ns] + S->NumBFs;
        G->PanelIndexOffset[ns+1]=G->PanelIndexOffset[ns] + S->NumPanels;
      };
   };

  /*--------------------------------------------------------------*/
  /*- allocate memory for the contributions of the innermost     -*/
  /*- lattice cells to the BEM matrix.                           -*/
  /*- note: P, M, Z stand for 'plus 1, minus 1, zero.'           -*/
  /*- Mab is the BEM interaction matrix between the unit-cell    -*/
  /*- geometry and a copy of itself translated through vector    -*/
  /*- a*LBV[0] + b*LBV[1].                                       -*/
  /*--------------------------------------------------------------*/
  MPP=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MPM=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MPZ=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MZP=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MZZ=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  // this one could be symmetric ...

  /*--------------------------------------------------------------*/
  /*- allocate interpolators for each region in the geometry.     */
  /*--------------------------------------------------------------*/
#if 0
  GBarAB9_Interior=(Interp3D **)mallocEC(G->NumSurfaces * sizeof(Interp3D *));
  double RMaxTS[3], RMinTS[3]; // max/min coordinates values for this surface
  //double RMax[3], RMin[3];   // max/min coord values for entire geometry (now class fields)
  double DeltaR[3];
  int NPoints[3];
  RWGSurface *S;
  RMax[0] = RMax[1] = RMax[2] = -1.0e+9;
  RMin[0] = RMin[1] = RMin[2] = +1.0e+9;
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     S=G->Surfaces[ns];
     GetSurfaceExtents(S, RMaxTS, RMinTS);
     for(int i=0; i<3; i++)
      {  DeltaR[i]   = fmax( RMaxTS[i] - RMinTS[i], PBCGeometry::DeltaInterp );
         NPoints[i]  = 1 + (2.0*DeltaR[i] / PBCGeometry::DeltaInterp );
         if (NPoints[i] < 2)
          NPoints[i]=2;
         RMax[i]     = fmax(RMax[i], RMaxTS[i]);
         RMin[i]     = fmin(RMin[i], RMinTS[i]);
      };

     // FIXME to handle 1D periodicity
     if ( S->IsPEC )
      { 
        Log(" Surface %s is PEC; no interpolation table needed",S->Label);
        GBarAB9_Interior[ns]=0;
      }
     else if ( NumStraddlers[MAXLATTICE*ns+0]==0 && NumStraddlers[MAXLATTICE*ns+1]==0 )
      { 
        Log(" Surface %s straddles no unit cell boundaries; no interpolation table needed",S->Label);
        GBarAB9_Interior[ns]=0;
      }
     else
      { 
        Log("Creating %ix%ix%i i-table for object %s",NPoints[0],NPoints[1],NPoints[2],S->Label);
        GBarAB9_Interior[ns]=new Interp3D( -DeltaR[0], DeltaR[0], NPoints[0],
                                           -DeltaR[1], DeltaR[1], NPoints[1],
                                                  0.0, DeltaR[2], 1 + NPoints[2]/2,
                                           2, 0, 0, 0);
      };

   };

  /*--------------------------------------------------------------*/
  /*- allocate interpolator for exterior medium ------------------*/
  /*--------------------------------------------------------------*/
  for(int i=0; i<3; i++)
   { DeltaR[i] = fmax( RMax[i] - RMin[i], PBCGeometry::DeltaInterp );
     NPoints[i]  = 1 + (2.0*DeltaR[i] / PBCGeometry::DeltaInterp );
     if (NPoints[i] < 2)
      NPoints[i]=2;
   };
  GBarAB9_Exterior=new Interp3D( -DeltaR[0], DeltaR[0], NPoints[0],
                                 -DeltaR[1], DeltaR[1], NPoints[1],
                                        0.0, DeltaR[2], 1 + NPoints[2]/2,
                                 2, 0, 0, 0);
#endif

}

} // namespace scuff
