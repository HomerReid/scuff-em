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
 * AddStraddlers.cc 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include <libhrutil.h>

#include "libscuff.h"
#include "cmatheval.h"

using namespace scuff;

/***************************************************************/
/* return 1 if the point with cartesian coordinates X lies on  */
/* the line connecting the origin to point L.                  */
/* Note that this routine is only looking at the first two     */
/* cartesian components of X; the z component is arbitrary.    */
/***************************************************************/
int PointOnLine(double *X, double *L)
{ 
   return ((float)(L[1]*X[0])) == ((float)(L[0]*X[1])) ? 1 : 0;
} 

/***************************************************************/
/* LBV = 'lattice basis vectors'                               */
/***************************************************************/
int FindPartnerEdge(RWGObject *O, int nei, double *LBV[2], double *V)
{
  RWGEdge *E = O->ExteriorEdges[nei];
  double *V1 = O->Vertices + 3*(E->iV1);
  double *V2 = O->Vertices + 3*(E->iV2);
  
  /*--------------------------------------------------------------*/
  /*- determine whether or not the edge lies on the unit cell     */
  /*- boundary, and if so which face of that boundary it lies on. */
  /*--------------------------------------------------------------*/
  double *ThisBV, *OtherBV;
  if ( PointOnLine(V1, LBV[0]) && PointOnLine(V2, LBV[0]) )
   { ThisBV=LBV[0]; 
     OtherBV=LBV[1];
   }
  else if ( PointOnLine(V1, LBV[1]) && PointOnLine(V2, LBV[1]) )
   { ThisBV=LBV[1]; 
     OtherBV=LBV[0];
   }
  else
   return -1; // edge does not lie on unit cell boundary 

  /*--------------------------------------------------------------*/
  /* Look for an exterior edge that is the translate through      */
  /* OtherBV of the present exterior edge.                        */
  /*--------------------------------------------------------------*/
  double V1T[3], V2T[3]; // 'V12, translated'
  VecScaleAdd(V1, 1.0, OtherBV, V1T);
  VecScaleAdd(V2, 1.0, OtherBV, V2T);
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
        QP = O->Vertices + 3*(O->ExteriorEdges[neip]->iQP);
        VecScaleAdd(QP, -1.0, OtherBV, V);
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
/* LBV[0][0..1] and LBV[1][0..1] are two 2-dimensional vectors  */
/* defining a lattice                                           */
/*                                                              */
/*--------------------------------------------------------------*/
#define CHUNK 100
void AddStraddlers(RWGObject *O, double **LBV)
{ 
  int NumNew=0, NumAllocated=0;
  double V[3], *NewVertices=0;
  RWGPanel *P, **NewPanels=0;
  RWGEdge *E, **NewEdges=0;

  int nei, neip;
  for(nei=0; nei<O->NumExteriorEdges; nei++)
   { 
      if ( O->ExteriorEdges[nei]==0 )
       continue;

      // see if this edge is a straddler, i.e. it lies on a face
      // of the unit cell and it has a partner (a translated 
      // version of itself) on the opposite side of the unit cell
      neip=FindPartnerEdge(O, nei, LBV, V);

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
          // edges anyway. 
          E=O->ExteriorEdges[nei];
          E->iQM  = O->NumVertices + NumNew;
          E->iMPanel = O->NumPanels + NumNew;
          E->Index = O->NumEdges + NumNew;
          E->Radius=VecDistance(E->Centroid, O->Vertices+3*E->iQP);
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,O->Vertices+3*E->iQM));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,O->Vertices+3*E->iV1));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,O->Vertices+3*E->iV2));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,O->Vertices+3*E->iV2));
          NewEdges[NumNew] = E;

          // what we just did was to combine two exterior edges 
          // into a single interior edge, so we now remove the 
          // exterior edges 
          free(O->ExteriorEdges[neip]);
          O->ExteriorEdges[nei]=O->ExteriorEdges[neip]=0;

          // add a new panel. Note that we can't call InitRWGPanel yet 
          // because the new vertices have yet been added to the Vertices 
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
  //NumStraddlers=NumNew;
  if (NumNew==0)
   return;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
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

  //NumStraddlers=NumNew;
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
  O->NumBFs = ( O->MP->IsPEC() ? NumEdges : 2*NumEdges );

  /*--------------------------------------------------------------*/
  /*- implement the kdtri thing here ... -------------------------*/
  /*--------------------------------------------------------------*/
}
