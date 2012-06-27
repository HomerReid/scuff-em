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
 * InitEdgeList.cc -- RWGComposite class method for gathering 
 *                 -- information on the interior and exterior
 *                 -- edges of an OpenSurface.
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

namespace scuff {

/***************************************************************/
/* InitEdgeList: After reading in a set of panels, we call     */
/* this function to extract all necessary information          */
/* regarding interior edges, exterior edges, boundary contours,*/
/* etc.                                                        */
/***************************************************************/
void RWGComposite::InitEdgeList(OpenSurface *OS)
{ 
  RWGPanel *P; 
  RWGEdge *E, ***EVEdges, *BCEdgeList;
  int i, np, ne, nv, nvp, iVLesser, iVGreater;
  int NumExteriorVertices, NumUnusedVertices;
  int *VertexUsed;
  double *VLesser, *VGreater;
  int *EVNumEdges;
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

  /*--------------------------------------------------------------*/
  /*- VertexUsed[nv] = 1 if vertex # nv is a vertex of any panel  */
  /*- on the object. (Used below in the determination of the      */
  /*- number of "interior" vertices.)                             */
  /*--------------------------------------------------------------*/
  VertexUsed=(int *)mallocEC(NumVertices*sizeof(int));
  memset(VertexUsed,0,NumVertices*sizeof(int));

  /*--------------------------------------------------------------*/
  /*- EVNumEdges[nv] is the number of edges connected to vertex   */
  /*- nv (only filled in for nv=exterior vertex). This number     */
  /*- should be 2 for all exterior vertices.                      */
  /*- EVEdges[nv][0] and EVEdges[nv][1] are pointers to the two   */
  /*- RWGEdge structures connected to vertex nv (again only      */
  /*- used for nv=exterior vertex).                               */
  /*--------------------------------------------------------------*/
  EVNumEdges=(int *)mallocEC(NumVertices*sizeof(RWGEdge));
  memset(EVNumEdges,0,NumVertices*sizeof(int));
  EVEdges=(RWGEdge ***)mallocEC(NumVertices*sizeof(RWGEdge **));
  EVEdges[0]=(RWGEdge **)mallocEC(2*NumVertices*sizeof(RWGEdge *)); 
  for(nv=1; nv<NumVertices; nv++)
   EVEdges[nv]=EVEdges[nv-1]+2;

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
  NumEdges=NumTotalEdges=0;
  for(np=0; np<NumPanels; np++)
   for(ne=0; ne<3; ne++)   /* loop over panel edges */
    { 
      P=Panels[np];

      /***************************************************************/
      /* get lesser & greater of the two vertex indices for this edge*/
      /***************************************************************/
      iVLesser=P->VI[ne];
      iVGreater=P->VI[(ne+1)%3];
      if ( iVLesser > iVGreater ) 
       { iVLesser=P->VI[(ne+1)%3];
         iVGreater=P->VI[ne];
       };

      VertexUsed[iVLesser]=VertexUsed[iVGreater]=1;

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

         E->Index=NumEdges++;

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
         E->iQM=-1;

	 VLesser=Vertices + 3*iVLesser;
         VGreater=Vertices + 3*iVGreater;
         for(i=0; i<3; i++)
          E->Centroid[i]=(VLesser[i] + VGreater[i]) / 2.0;
         E->Length=VecDistance(VLesser, VGreater);

         E->iPPanel=P->Index;
         E->PIndex=(ne+2)%3;

         E->iMPanel=-1;
         E->Index=-1;

       };
    };

  /*--------------------------------------------------------------*/
  /*- now go back through our list of all edges:                  */
  /*-  a. put each RWGEdge structure corresponding to an interior */ 
  /*      edge into the Edges array.                              */
  /*-  b. put each RWGEdge structure corresponding to an exterior */ 
  /*      edge into ExteriorEdges array and also into the EVEdges */
  /*      lists for the edge's 2 vertices.                        */
  /*      note that the Index field of the RWGEdge struct         */
  /*      for an exterior edge is set to -(i+1), where i is the   */
  /*      index of the edge in the ExteriorEdges array. (thus the */
  /*      first exterior edge has Index=-1, the second has        */
  /*      Index=-2, etc.)                                         */
  /*--------------------------------------------------------------*/
  Edges=(RWGEdge **)mallocEC(NumEdges*sizeof(Edges[0]));
  ExteriorEdges=(RWGEdge **)mallocEC((NumTotalEdges-NumEdges)*sizeof(Edges[0]));
  NumExteriorEdges=0;
  NumExteriorVertices=0;
  for(nv=0; nv<NumVertices; nv++)
   for(E=EdgeLists[nv]; E; E=E->Next)
    { 
      if (E->Index==-1)  /* E is an exterior edge */
       { 
         ExteriorEdges[NumExteriorEdges]=E;
         E->Index=-(NumExteriorEdges+1);
         NumExteriorEdges++;

         if (EVNumEdges[E->iV1]==2) 
          ErrExit("%s: invalid mesh topology: vertex %i",MFN,E->iV1);
         EVEdges[E->iV1][ EVNumEdges[E->iV1]++ ] = E;

         if (EVNumEdges[E->iV2]==2) 
          ErrExit("%s: invalid mesh topology: vertex %i",MFN,E->iV2);
         EVEdges[E->iV2][ EVNumEdges[E->iV2]++ ] = E;

         NumExteriorVertices++;
       }
      else               /* E is an interior edge */
       Edges[E->Index]=E;
    };
  if ( (NumExteriorEdges+NumEdges) != NumTotalEdges )
   ErrExit("%s:%i: internal error (%i!=%i)",__FILE__,__LINE__,NumExteriorEdges,NumTotalEdges-NumEdges);

  /*--------------------------------------------------------------*/
  /*- now go through and classify exterior boundary contours     -*/
  /*- using the following algorithm:                             -*/
  /*-  1. pick, at random, an exterior vertex (call it vertex    -*/
  /*-     #nv).                                                  -*/ 
  /*   2. the EVEdges array for vertex nv contains two edges.    -*/
  /*-     take the first of these edges (call it edge E) and     -*/
  /*-     follow it to its other vertex (call it vertex #nvp).   -*/
  /*-  3. the EVEdges array for vertex nvp contains two edges,   -*/
  /*-     one of which is E. call the other one Ep. follow Ep    -*/
  /*-     to its other vertex (call it vertex #nvpp).            -*/
  /*-  4. continue in this way to traverse the edges on a single -*/
  /*-     boundary contour until we come back to vertex #nv. at  -*/
  /*-     this point we have just classified a single exterior   -*/
  /*-     boundary contour. we collect all the edges in this     -*/
  /*-     countour and identify them as belonging to exterior    -*/
  /*-     boundary contour #1. also, we mark all the exterior    -*/
  /*-     vertices in this contour as having been already        -*/
  /*-     dealt with.                                            -*/
  /*-  5. now repeat from step 1: choose another exterior vertex -*/
  /*-     (one that is not contained in the list of exterior     -*/
  /*-      vertices comprising the boundary contour we just      -*/
  /*-      traversed) and traverse a loop of edges until we come -*/
  /*-      back to the original vertex, then call this whole     -*/
  /*-      boundary contour #2.                                  -*/
  /*-  6. continue in this way until all exterior boundary       -*/
  /*-     contours have been identified.                         -*/
  /*--------------------------------------------------------------*/
  WhichBC=(int *)mallocEC(NumVertices*sizeof(int));
  memset(WhichBC,0,NumVertices*sizeof(int));
  NumBCs=0;
  NumExteriorVertices=0;
  NumBCEdges=0;
  BCEdges=0;
  for(;;)
   {   
     /***************************************************************/
     /* step 1: find an exterior vertex                             */
     /***************************************************************/
     for(nv=0; nv<NumVertices; nv++)
      if (EVNumEdges[nv]>0) 
       break;

     if (nv==NumVertices) break; /* no more exterior vertices left */

     NumBCEdges=(int *)realloc(NumBCEdges, (NumBCs+1)*sizeof(int));
     BCEdges=(RWGEdge ***)realloc(BCEdges, (NumBCs+1)*sizeof(RWGEdge **));

     /*****************************************************************/
     /* step 2--4: traverse the boundary contour containing vertex nv */
     /*****************************************************************/
     BCEdgeList=0;
     E=EVEdges[nv][0];
     nvp=nv;
     NumBCEdges[NumBCs]=0;
     do
      { 
        /* verify that this vertex is connected to exactly 2 exterior edges */
        if ( EVNumEdges[nvp]!=2 )
         ErrExit("%s: invalid mesh topology: vertex %i",MFN,nvp);

        /* mark this vertex as having been visited */
        EVNumEdges[nvp]=0; 
        WhichBC[nvp]=NumBCs;
        NumExteriorVertices++;

        /* add E to list of edges for this boundary contour. */
        E->Next=BCEdgeList;
        BCEdgeList=E;
        NumBCEdges[NumBCs]++;
       
        /* set nvp equal to next vertex in boundary contour */
        if ( nvp==E->iV1 )
         nvp=E->iV2;
        else if ( nvp==E->iV2 )
         nvp=E->iV1;
        else
         ErrExit("%s:%i: internal error",__FILE__,__LINE__);

        /* set E equal to next edge in boundary contour */
        if ( E==EVEdges[nvp][0] )  
         E=EVEdges[nvp][1];
        else if ( E==EVEdges[nvp][1] )  
         E=EVEdges[nvp][0];
        else
         ErrExit("%s:%i: internal error",__FILE__,__LINE__);

      } while (nvp!=nv);

     /*****************************************************************/
     /* step 5: create an array containing all the edges we just      */
     /* visited in this boundary contour.                             */
     /*****************************************************************/
     BCEdges[NumBCs]=(RWGEdge **)mallocEC(NumBCEdges[NumBCs]*sizeof(RWGEdge *));
     for(ne=0, E=BCEdgeList; E; E=E->Next)
      BCEdges[NumBCs][ne++]=E;

     NumBCs++;
 
   }; // for(;;)

  /*--------------------------------------------------------------*/
  /*- count unused vertices --------------------------------------*/
  /*--------------------------------------------------------------*/
  NumUnusedVertices=0;
  for(nv=0; nv<NumVertices; nv++)
   if (VertexUsed[nv]==0) 
    NumUnusedVertices++;

  NumInteriorVertices=NumVertices-NumExteriorVertices-NumUnusedVertices;

  /*--------------------------------------------------------------*/
  /*- deallocate temporary storage -------------------------------*/
  /*--------------------------------------------------------------*/
  free(EdgeLists);
  free(VertexUsed);
  free(EVNumEdges);
  free(EVEdges[0]);
  free(EVEdges);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
#if 0
  FILE *f=vfopen("%s.EdgeInfo","w",GetFileBase(MeshFileName)); 
  fprintf(f,"%i total edges\n",NumTotalEdges);
  fprintf(f,"%i interior edges\n",NumEdges);
  fprintf(f,"%i interior edges, take 2\n",NumTotalEdges-NumExteriorEdges);
  fprintf(f,"%i panels\n",NumPanels);
  fprintf(f,"%i total vertices\n",NumVertices);
  fprintf(f,"%i interior vertices\n",NumInteriorVertices);
  fprintf(f,"%i boundary contours\n",NumBCs);
  for(int nbc=1; nbc<=NumBCs; nbc++)
   { fprintf(f,"\n** boundary contour %i: \n",NumBCs);

     fprintf(f,"    vertices:");
     for(nv=0; nv<NumVertices; nv++)
      if(WhichBC[nv]==nbc)
       fprintf(f," %i",nv);
     fprintf(f,"\n");

     fprintf(f,"    edges:");
     for(ne=0; ne<NumBCEdges[nbc]; ne++)
      fprintf(f," %i ",-(BCEdges[nbc][ne]->Index));
     fprintf(f,"\n");

   };
  fclose(f);
#endif

}

} // namespace scuff
