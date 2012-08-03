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
 * AddStraddlers.cc -- a routine to analyze a unit-cell geometry
 *                     and add new vertices, edges, and panels for 
 *                     RWG basis functions that straddle the unit 
 *                     cell boundary
 *
 * homer reid       -- 7/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include <libhrutil.h>
#include <libscuff.h>
#include "PBCGeometry.h"

namespace scuff {

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
/* given a single exterior edge in a meshed geometry, look for */
/* a partner of this edge -- that is, another exterior edge    */
/* that is a translate of the original edge through either     */
/* LBV[0] or LBV[1].                                           */
/*                                                             */
/* (note: LBV = 'lattice basis vector.')                       */
/*                                                             */
/* if a partner edge is found, its index within O's            */
/* ExteriorEdges array is returned. otherwise, -1 is returned. */
/*                                                             */
/* if a partner edge is found, then NumStraddlers[i] is        */
/* incremented (where i=0 or 1 depending on which unit cell    */
/* boundary the original edge lay on) and V[0..2] is filled in */
/* with the cartesian coordinates of the new vertex that must  */
/* be added to turn the half-RWG basis function associated with*/
/* edge #nei into a full RWG basis function.                   */
/***************************************************************/
static int FindPartnerEdge(RWGObject *O, int nei, double *LBV[2], 
                           int NumStraddlers[2], double *V)
{
  RWGEdge *E = O->ExteriorEdges[nei];
  double *V1 = O->Vertices + 3*(E->iV1);
  double *V2 = O->Vertices + 3*(E->iV2);
  
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
  for(int neip=0; neip<O->NumExteriorEdges; neip++)
   { 
     if (O->ExteriorEdges[neip]==0) 
      continue;

     V1P = O->Vertices + 3*(O->ExteriorEdges[neip]->iV1);
     V2P = O->Vertices + 3*(O->ExteriorEdges[neip]->iV2);
     if (   (VecEqualFloat(V1T, V1P) && VecEqualFloat(V2T, V2P))
         || (VecEqualFloat(V1T, V2P) && VecEqualFloat(V2T, V1P))
        )
      { 
        /*--------------------------------------------------------------*/
        /*- found a translate of the edge in question.                  */
        /*--------------------------------------------------------------*/
        memcpy(V, O->Vertices + 3*(O->ExteriorEdges[neip]->iQP), 3*sizeof(double));
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
          "on the opposite side of the unit cell",nei, O->Label);

  return -1; // i am indecisive about how to handle this error

}


/*--------------------------------------------------------------*/
/* A 'straddler' is a panel edge that lies on a face of the unit*/
/* cell and has a partner edge on the opposite face of the unit */
/* cell. This routine goes through the list of exterior edges   */
/* in the given RWGObject and identifies all straddlers. For    */
/* each straddler, the routine adds one new panel and one new   */
/* vertex to the Object. The new panel is the translate of the  */
/* panel to which the partner edge is attached, while the new   */
/* vertex is the vertex opposite the partner edge within that   */
/* panel. Once we have added this new vertex and panel, we can  */
/* promote the exterior edge to an interior edge.               */
/*                                                              */
/* Parameters:                                                  */
/*  O points to the RWGObject in question; its contents are     */
/*  modified if any straddlers were detected.                   */
/*                                                              */
/*  LBV[i][j] is the jth cartesian component of the ith lattice */
/*  basis vector (i,j=0,1.)                                     */
/*                                                              */
/*  On return, NumStraddlers[i] is the number of straddlers     */
/*  that were detected on unit cell boundary #i (i=0,1.)        */
/*--------------------------------------------------------------*/
#define CHUNK 100
void AddStraddlers(RWGObject *O, double **LBV, int *NumStraddlers)
{ 
  int NumNew=0, NumAllocated=0;
  double V[3], *NewVertices=0;
  RWGPanel *P, **NewPanels=0;
  RWGEdge *E, **NewEdges=0;

  if (NumStraddlers)
   memset(NumStraddlers, 0, 2*sizeof(int));

  int nei, neip;
  for(nei=0; nei<O->NumExteriorEdges; nei++)
   { 
      if ( O->ExteriorEdges[nei]==0 )
       continue;

      // see if this edge is a straddler, i.e. it lies on a face
      // of the unit cell and it has a partner (a translated 
      // version of itself) on the opposite side of the unit cell
      neip=FindPartnerEdge(O, nei, LBV, NumStraddlers, V);

      // if so, add a new vertex, panel, and edge.
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
          // will be removing it from object O's set of exterior 
          // edges anyway. note that the iQP, iV1, iV2, iPPanel, PIndex, 
          // Centroid, and Length fields of this structure will be already 
          // correctly initialized. 
          E=O->ExteriorEdges[nei];
          E->iQM  = O->NumVertices + NumNew;
          E->iMPanel = O->NumPanels + NumNew;
          E->Index = O->NumEdges + NumNew;
          E->MIndex = 2; // because below we hard-code P->VI[2] = E->iQM;
          E->Radius=VecDistance(E->Centroid, O->Vertices+3*E->iQP);
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,O->Vertices+3*E->iQM));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,O->Vertices+3*E->iV1));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,O->Vertices+3*E->iV2));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,O->Vertices+3*E->iV2));
          NewEdges[NumNew] = E;

          // what we just did was to combine two exterior edges 
          // into a single interior edge, so we now remove those
          // two exterior edges from O's list of exterior edges
          free(O->ExteriorEdges[neip]);
          O->ExteriorEdges[nei]=O->ExteriorEdges[neip]=0;

          // add a new panel. Note that we can't call InitRWGPanel yet 
          // because the new vertex has not yet been added to the Vertices 
          // array; this happens later, below.
          P=(RWGPanel *)mallocEC(sizeof(RWGPanel));
          P->VI[0] = E->iV1;
          P->VI[1] = E->iV2;
          P->VI[2] = E->iQM;
          P->Index = O->NumPanels + NumNew;
          NewPanels[NumNew] = P;
  
          NumNew++;

       };
     
   }; // for(nei=0; nei<O->NumExteriorEdges; nei++)

  if (NumNew==0)
   return;

  /*--------------------------------------------------------------*/
  /*- replace O's arrays of Vertices, Panels, Edges, and          */
  /*- ExteriorEdges with new arrays that include the straddlers   */
  /*--------------------------------------------------------------*/
  double *Vertices        = O->Vertices;
  RWGPanel **Panels       = O->Panels;
  RWGEdge **Edges         = O->Edges;
  RWGEdge **ExteriorEdges = O->ExteriorEdges;
  int NumVertices         = O->NumVertices;
  int NumPanels           = O->NumPanels;
  int NumEdges            = O->NumEdges;
  int NumExteriorEdges    = O->NumExteriorEdges;

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
 
  O->Vertices         = Vertices;
  O->Panels           = Panels;
  O->Edges            = Edges;
  O->ExteriorEdges    = ExteriorEdges;
  O->NumVertices      = NumVertices;
  O->NumPanels        = NumPanels;
  O->NumEdges         = NumEdges;
  O->NumExteriorEdges = NumExteriorEdges;

  O->NumTotalEdges=NumEdges + NumExteriorEdges;
  O->NumBFs = ( O->IsPEC ? NumEdges : 2*NumEdges );

  /*--------------------------------------------------------------*/
  /*- implement the kdtri thing here ... -------------------------*/
  /*--------------------------------------------------------------*/
}

} // namespace scuff
