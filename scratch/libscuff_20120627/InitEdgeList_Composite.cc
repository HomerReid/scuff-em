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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPPSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * InitEdgeList.cc -- RWGComposite class method for gathering 
 *                 -- information on the interior and exterior
 *                 -- edges of an PartialSurface.
 *                 
 *                 -- This is really just a part of the RWGComposite
 *                 -- class constructor, but we put it in a separate
 *                 -- file for clarity.
 *
 * homer reid      -- 6/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libscuff.h"
#include "RWGComposite.h"

namespace scuff {

/***************************************************************/
/* InitEdgeList: After identifying all the PANELS on an open   */
/* surface, we call this routine to classify all the EDGES on  */
/* that surface, and in particular to construct arrays of      */
/* full and half RWG basis functions.                          */
/***************************************************************/
void RWGComposite::InitEdgeList(PartialSurface *PS)
{ 
  RWGPanel *P; 
  RWGEdge *E;
  int i, np, ne, nv, nvp, iVLesser, iVGreater;
  double *VLesser, *VGreater;
  char *MFN=MeshFileName;
  RWGEdge **EdgeLists;

  /***************************************************************/
  /***************************************************************/
  /*- allocate temporary storage arrays  *************************/
  /***************************************************************/
  /***************************************************************/

  /*--------------------------------------------------------------*/
  /*- EdgeList[nv] is a linked list of all Edge structures for   -*/
  /*- which the lesser-numbered of the two vertex indices is nv. -*/
  /*--------------------------------------------------------------*/
  EdgeLists=(RWGEdge **)mallocEC(NumVertices*sizeof(RWGEdge *));
  memset(EdgeLists,0,NumVertices*sizeof(RWGEdge*));

  /****************************************************************/
  /* construct a list of RWGEdge structures by going over every  */
  /* edge of every panel:                                         */ 
  /*                                                              */
  /*  a. if the edge does not exist in the list (we have not      */
  /*     visited this edge before) then add it to the list        */
  /*  b. if the edge does exist in the list, then we are visiting */
  /*     it for the second time, so it is an interior edge, and   */
  /*     we mark it as such.                                      */
  /*  c. if we find ourselves visiting an edge more than two      */
  /*     times then the topology of the panel set is defective.   */
  /*                                                              */
  /* note that when we first create a new RWGEdge structure for  */
  /* edge, we assign the panel it came from as its 'positive'     */
  /* panel, but we don't yet assign it an 'Index' in the overall  */
  /* problem, because it may turn out to be an exterior edge. we  */
  /* also tentatively set the 'EdgeSign' and 'EdgeIndex' fields   */
  /* for this edge in the corresponding panel structure to 0.0    */
  /* and -1 in view of this possibility.                          */
  /*                                                              */
  /* when we come across the edge for the second time (because it */
  /* turned out to be attached to a second triangle) we call the  */
  /* second triangle the 'negative' panel associated with the     */
  /* edge, and now we do assign an Index to this edge.            */
  /*                                                              */
  /* thus, at the conclusion of this loop, any Edge structure that*/
  /* still has Index==-1 is an exterior edge. such a structure    */
  /* will have valid data stored for its iV1, iV2, iQP, PPanel,   */
  /* iPPanel, and PIndex fields, but not its iQM, MPanel, iMPanel,*/
  /* or MIndex fields.                                            */
  /****************************************************************/
  int NumInteriorEdges=0, NumTotalEdges=0;
  for(np=0; np<PS->NumPanels; np++)
   for(ne=0; ne<3; ne++)   /* loop over panel edges */
    { 
      P=PS->Panels[np];

      /***************************************************************/
      /* get lesser & greater of the two vertex indices for this edge*/
      /***************************************************************/
      iVLesser=P->VI[ne];
      iVGreater=P->VI[(ne+1)%3];
      if ( iVLesser > iVGreater ) 
       { iVLesser=P->VI[(ne+1)%3];
         iVGreater=P->VI[ne];
       };

      /**********************************************************************/
      /* look for this edge in list of edges connected to vertex # iVLesser */
      /**********************************************************************/
      for(E=EdgeLists[iVLesser]; E; E=E->Next)
       if (E->iV2==iVGreater) 
        break;
      
      if ( E )
       { 
         /***************************************************************/
         /* we have encountered this edge twice before ******************/
         /***************************************************************/
         if ( E->iMPanel != -1 )
          ErrExit("%s: invalid mesh topology: edge %i of panel %i also belongs to panels %i and %i ",
                      MFN,ne,np,E->iPPanel,E->iMPanel);

         /***************************************************************/
         /* we have encountered this edge once before                   */
         /***************************************************************/
         E->iQM=P->VI[(ne+2)%3];
         E->iMPanel=P->Index;
         E->MIndex=(ne+2)%3;

         E->Index=NumInteriorEdges++;

         /* bounding radius is max distance from centroid to any vertex */
         E->Radius=VecDistance(E->Centroid, Vertices+3*E->iQP);
         E->Radius=fmax(E->Radius, VecDistance(E->Centroid,Vertices+3*E->iQM));
         E->Radius=fmax(E->Radius, VecDistance(E->Centroid,Vertices+3*E->iV1));
         E->Radius=fmax(E->Radius, VecDistance(E->Centroid,Vertices+3*E->iV2));

       }
      else
       { 
         /******************************************************************/
         /* we are encountering this edge for the first time. create a     */ 
         /* new RWGEdge structure for this edge and chain it in to the     */
         /* linked list of RWGEdge structures connected to vrtx #iVLesser  */
         /******************************************************************/
         NumTotalEdges++;
         E=(RWGEdge *)mallocEC(sizeof *E);
         E->Next=EdgeLists[iVLesser];
         EdgeLists[iVLesser]=E;

         E->iV1=iVLesser;
         E->iV2=iVGreater;
         E->iQP=P->VI[ (ne+2)%3 ];

	 VLesser=Vertices + 3*iVLesser;
         VGreater=Vertices + 3*iVGreater;
         for(i=0; i<3; i++)
          E->Centroid[i]=(VLesser[i] + VGreater[i]) / 2.0;
         E->Length=VecDistance(VLesser, VGreater);

         E->iPPanel=P->Index;
         E->PIndex=(ne+2)%3;

         E->iQM=-1;
         E->MIndex=-1;
         E->iMPanel=-1;
         E->Index=-1;

       };
    };

  /*--------------------------------------------------------------*/
  /*- now go back through our list of all edges:                  */
  /*-  a. put each RWGEdge structure corresponding to an interior */ 
  /*      edge into the Edges array.                              */
  /*-  b. put each RWGEdge structure corresponding to an exterior */
  /*      edge into the HEdges array, and assign the Index field  */
  /*      of that structure to its index witin the HEdges array.  */
  /*      Note that the Edge structures in the HEdges entry have  */
  /*      the value -1 for their iQM, iMPanel, and MIndex fields. */
  /*--------------------------------------------------------------*/
  int NumExteriorEdges=NumTotalEdges-NumInteriorEdges; 
  PS->Edges=(RWGEdge **)mallocEC(NumInteriorEdges*sizeof(RWGEdge *));
  PS->HEdges=(RWGEdge **)mallocEC(NumExteriorEdges*sizeof(RWGEdge *));
  int nee2=0; // at the end of the operation, we should have nee2=NumExteriorEdges
  for(nv=0; nv<NumVertices; nv++)
   for(E=EdgeLists[nv]; E; E=E->Next)
    { 
      if (E->Index==-1)  /* E is an exterior edge */
       { 
         E->Index=nee2;
         PS->HEdges[nee2]=E;
         nee2++;
       }
      else               /* E is an interior edge */
       PS->Edges[E->Index]=E;
    };
  if ( (nee2!=NumExteriorEdges ) )
   ErrExit("%s:%i: internal error (%i!=%i)",__FILE__,__LINE__,nee2,NumExteriorEdges);

  PS->NumEdges=NumInteriorEdges;
  PS->NumHEdges=NumExteriorEdges;
  PS->NumTotalEdges=NumTotalEdges;

  /*--------------------------------------------------------------*/
  /*- deallocate temporary storage -------------------------------*/
  /*--------------------------------------------------------------*/
  free(EdgeLists);

} // InitEdgeList

} // namespace scuff
