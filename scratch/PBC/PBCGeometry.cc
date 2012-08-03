/*
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libhrutil.h>
#include <libMDInterp.h>
#include <libscuff.h>
#include "PBCGeometry.h"

namespace scuff{

/***************************************************************/
/* do something about me please ********************************/
/***************************************************************/
int PBCGeometry::TriangleCubatureOrder=4;
double PBCGeometry::DeltaInterp=0.05;

/***************************************************************/
/* get the maximum and minimum cartesian coordinates obtained  */
/* by points on an RWGObject. (we do this by looking at the    */
/* panel vertices, which suffices because the panels are       */
/* flat and thus convex)                                       */
/***************************************************************/
void GetObjectExtents(RWGObject *O, double RMax[3], double RMin[3])
{ 
  memcpy(RMax, O->Vertices + 0, 3*sizeof(double));
  memcpy(RMin, O->Vertices + 0, 3*sizeof(double));
  for(int nv=1; nv<O->NumVertices; nv++)
   { 
     RMax[0] = fmax(RMax[0], O->Vertices[3*nv + 0]);
     RMin[0] = fmin(RMin[0], O->Vertices[3*nv + 0]);
     RMax[1] = fmax(RMax[1], O->Vertices[3*nv + 1]);
     RMin[1] = fmin(RMin[1], O->Vertices[3*nv + 1]);
     RMax[2] = fmax(RMax[2], O->Vertices[3*nv + 2]);
     RMin[2] = fmin(RMin[2], O->Vertices[3*nv + 2]);
   };
}

/***************************************************************/
/* PBCGeometry class constructor *******************************/
/***************************************************************/
PBCGeometry::PBCGeometry(RWGGeometry *pG, double **pLBV)
{
  /*--------------------------------------------------------------*/
  /*- initialize simple class fields -----------------------------*/
  /*--------------------------------------------------------------*/
  G=pG;
  LBV[0][0]=pLBV[0][0];
  LBV[0][1]=pLBV[0][1];
  LBV[1][0]=pLBV[1][0];
  LBV[1][1]=pLBV[1][1];
  CurrentOmega=-1.0;
  EpsTF = (cdouble *)mallocEC( (G->NumObjects+1)*sizeof(cdouble));
  MuTF  = (cdouble *)mallocEC( (G->NumObjects+1)*sizeof(cdouble));

  Log("Creating PBC geometry: unit cell geometry %s, lattice (%g,%g) x (%g,%g)",
       G->GeoFileName,LBV[0][0],LBV[0][1],LBV[1][0],LBV[1][1]); 

  /*--------------------------------------------------------------*/
  /*- add straddlers to all objects ... note that this modifies  -*/
  /*- G and its constituent objects. it would be better program  -*/
  /*- design to make a copy of G?                                -*/
  /*--------------------------------------------------------------*/
  NumStraddlers=(int *)mallocEC(2*G->NumObjects*sizeof(int));
  G->TotalBFs=G->TotalPanels=0;
  double *MyLBVs[2];
  MyLBVs[0]=LBV[0];
  MyLBVs[1]=LBV[1];
  for(int no=0; no<G->NumObjects; no++)
   {
     RWGObject *O = G->Objects[no];
     AddStraddlers(O, MyLBVs, NumStraddlers + 2*no);
 
     // FIXME
     if ( (NumStraddlers[2*no+0]==0) != (NumStraddlers[2*no+1]==0) )
      ErrExit("object %s: 1D straddling periodicity is not yet supported",O->Label); 

     G->TotalBFs    += O->NumBFs;
     G->TotalPanels += O->NumPanels;
     if ( no+1 < G->NumObjects )
      { G->BFIndexOffset[no+1]=G->BFIndexOffset[no] + O->NumBFs;
        G->PanelIndexOffset[no+1]=G->PanelIndexOffset[no] + O->NumPanels;
      };

     Log(" Detected (%i,%i) straddlers for object %s", NumStraddlers[2*no+0],NumStraddlers[2*no+1],O->Label);

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
  MZZ=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX); // this one could be symmetric ...

  /*--------------------------------------------------------------*/
  /*- allocate interpolators for each object interior             */
  /*--------------------------------------------------------------*/
  GBarAB9_Interior=(Interp3D **)mallocEC(G->NumObjects * sizeof(Interp3D *));
  double RMaxTO[3], RMinTO[3]; // max/min coordinates values for this object
  //double RMax[3], RMin[3];   // max/min coord values for entire geometry (now class fields)
  double DeltaR[3];
  int NPoints[3];
  RWGObject *O;
  RMax[0] = RMax[1] = RMax[2] = -1.0e+9;
  RMin[0] = RMin[1] = RMin[2] = +1.0e+9;
  for(int no=0; no<G->NumObjects; no++)
   { 
     O=G->Objects[no];
     GetObjectExtents(O, RMaxTO, RMinTO);
     for(int i=0; i<3; i++)
      {  DeltaR[i]   = fmax( RMaxTO[i] - RMinTO[i], PBCGeometry::DeltaInterp );
         NPoints[i]  = 1 + (2.0*DeltaR[i] / PBCGeometry::DeltaInterp );
         if (NPoints[i] < 2)
          NPoints[i]=2;
         RMax[i]     = fmax(RMax[i], RMaxTO[i]);
         RMin[i]     = fmin(RMin[i], RMinTO[i]);
      };

     // FIXME to handle 1D periodicity
     if ( O->IsPEC )
      { 
        Log(" Object %s is PEC; no interpolation table needed",O->Label);
        GBarAB9_Interior[no]=0;
      }
     else if ( NumStraddlers[2*no+0]==0 && NumStraddlers[2*no+1]==0 )
      { 
        Log(" Object %s straddles no unit cell boundaries; no interpolation table needed",O->Label);
        GBarAB9_Interior[no]=0;
      }
     else
      { 
        Log("Creating %ix%ix%i i-table for object %s",NPoints[0],NPoints[1],NPoints[2],O->Label);
        GBarAB9_Interior[no]=new Interp3D( -DeltaR[0], DeltaR[0], NPoints[0],
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
/* note: this routine is a class method within PBCGeometry, but */
/* it actually doesn't depend on anything within the PBCGeometry*/
/* class and could just as easily be standalone.                */
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
void PBCGeometry::AddStraddlers(RWGObject *O, double **LBV, int *NumStraddlers)
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
