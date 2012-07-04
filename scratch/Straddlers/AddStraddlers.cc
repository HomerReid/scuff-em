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

namespace scuff {

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
typedef struct NewVPEs
 {
  
 } 


/*--------------------------------------------------------------*/
/*- Result= 0: the edge is not a straddler.                     */
/*-         1: the edge is a proper straddler. in this case,    */
/*-            on return Q[0..2] are the coordinates of the     */
/*-            new vertex.                                      */
/*-            coordinates of the dler and Q[0..2] are the      */
/*-         2: the edge is an improper straddler.               */
/*--------------------------------------------------------------*/
int StraddleCheck(int nei, double Q[3])

/*--------------------------------------------------------------*/
/*                                                              */
/* L1 and L2 are two-dimensional lattice vectors                */
/*                                                              */
/* algorithm:                                                   */
/*                                                              */
/*--------------------------------------------------------------*/
#define CHUNK 100
void AddStraddlers(RWGObject *O, double *L1, double *L2)
{ 
  int NumNew=0, NumAllocated=0;
  double V[3], *NewVertices=0;
  RWGPanel *P, **NewPanels=0;
  RWGEdge *E, **NewEdges=0;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int Result;
  for(int nei=0; nei<NumExteriorEdges; nei++)
   { 
      Result=StraddleCheck(O, nei, L1, L2);
                           
      if ( Result==1 )
       { 
          // add a new vertex, panel, and edge
          if( NumAllocated == NumNew )
           { NumAllocated+=CHUNK;
             NewVertices = (double *)reallocEC(NewVertices, NumAllocated*sizeof(double));
             NewPanels   = (RWGPanel **)reallocEC(NewPanels, NumAllocated*sizeof(NewPanels *));
             NewEdges    = (RWGEdge **)reallocEC(NewEdges, NumAllocated*sizeof(NewEdges *));
           };
      
          memcpy( NewVertices + 3*NumNew, V, 3*sizeof(double));
 
          E=(RWGEdge *)mallocEC(sizeof(RWGEdge));
          memcpy(E, O->ExteriorEdges[nei], sizeof(RWGEdge) );
          E->iQM  = O->NumVertices + NumNew;
          E->iMPanel = O->Panels + NumNew;
          E->Index = O->NumEdges + NumNew;
          E->Radius=VecDistance(E->Centroid, Vertices+3*E->iQP);
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,Vertices+3*E->iQM));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,Vertices+3*E->iV1));
          E->Radius=fmax(E->Radius, VecDistance(E->Centroid,Vertices+3*E->iV2));
          NewEdges[NumNew] = E;

          (RWGPanel *)mallocEC(sizeof(RWGPanel));
          P->VI[0] = E->iV1;
          P->VI[1] = E->iV2;
          P->VI[2] = E->iQM;
          NewPanels[NumNew] = P;
          // can't call InitRWGPanel yet because not all the panel vertices
          // have been added to the Vertices array; this happens later, below.

       };
     
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  

}

  // namespace scuff
