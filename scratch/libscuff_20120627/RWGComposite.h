/*
 *
 */

#ifndef RWGCOMPOSITE_H
#define RWGCOMPOSITE_H

#include <libscuff.h>

/***************************************************************/
/* An PartialSurface is a basically just a collection of RWG   */
/* and half-RWG basis functions.                               */
/*                                                             */
/* Note: (slightly tricky) All vertex indices in the RWGPanel, */
/* RWGEdge, and RWGHEdge structures contained within an        */
/* PartialSurface refer to the Vertex[] array in the parent    */
/* RWGComposite structure. However, the *panel* indices in     */
/* the RWGEdge and RWGHEdge structures refer to the Panels[]   */
/* array within the PartialSurface structure, NOT the Panels[] */
/* array within the parent RWGComposite structure.             */
/***************************************************************/
typedef struct PartialSurface
{  
  int NumPanels;
  RWGPanel *Panels;

  int NumEdges;
  RWGEdge *Edges;   // internal edges to which we assign a full RWG function

  int NumHEdges;
  RWGEdge *HEdges;  // external edges to which we assign a half-RWG function

  int NumTotalEdges;

} PartialSurface;

/***************************************************************/
/* RWGComposite is a generalization of RWGObject.              */
/***************************************************************/
class RWGComposite
 {
   /*--------------------------------------------------------------*/
   /*- class methods ----------------------------------------------*/
   /*--------------------------------------------------------------*/
   RWGComposite::RWGComposite(FILE *f, const char *pLabel, int *LineNum);
   ~RWGComposite::RWGComposite();

   /*--------------------------------------------------------------*/
   /*- class data -------------------------------------------------*/
   /*--------------------------------------------------------------*/
   // the x,y,z coordinate of the nth vertex are Vertices[3*n + 0,1,2]
   double *Vertices;  
   int NumVertices;

   // this is a full list of *all* panels on all PartialSurfaces 
   // in the composite. 
   RWGPanel *Panels;
   int NumPanels;

   // for nps = 0 , 1, ..., NumPartialSurfaces-1, 
   // SubRegions[2*nps+0] and Regions[2*nps+1] are the indices of 
   // the two SubRegions bounded by PartialSurface #nps.
   int *SubRegions;
   int NumSubRegions;

   // SubRegionMPs[0]     = material properties of exterior medium
   // SubRegionMPs[nsr+1] = material properties of subregion #nsr (nsr=0,...,NumSubRegions-1)
   MatProp **SubRegionMPs;

   // PartialSurface structures for each section of the composite surface
   PartialSurface *PartialSurfaces;
   char *PartialSurfaceLabels;
   int NumPartialSurfaces;
   int *NumPanelsPerPartialSurface;

   // BFIndexOffset[npsa] is the index of the first basis function
   // on partial surface #npsa within the overall vector of basis 
   // functions forthis RWGComposite.
   int *BFIndexOffset;

   int TotalBFs;

   // EpsTF[nr] = Epsilon for subregion #nsr at the present frequency
   //             (nr==0 for exterior medium)                   
   //             (note: TF stands for 'this frequency')
   cdouble *EpsTF;
   cdouble *MuTF;

   char *Label; 
   char *MeshFileName;
   char *ErrMsg;        // used to indicate when something goes wrong

   /* GT encodes any transformation that has been carried out since */
   /* the composite was read from its mesh file (not including a    */
   /* possible one-time GTransformation that may have been specified*/
   /* in the .scuffgeo file when the object was first created.)     */
   GTransformation *GT;

   /*--------------------------------------------------------------*/
   /*- private class methods --------------------------------------*/
   /*--------------------------------------------------------------*/
   void InitRWGComposite(const char *pMeshFileName, const GTransformation *OTGT=0);
   void ReadGMSHFile(FILE *MeshFile, char *FileName, const GTransformation *GT);
   void InitEdgeList();

 };
