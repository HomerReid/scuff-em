/*
 *
 */

#ifndef RWGCOMPOSITE_H
#define RWGCOMPOSITE_H

#include <libscuff.h>

/***************************************************************/
/* An OpenSurface is a basically just a collection of RWG and  */
/* half-RWG basis functions.                                   */
/*                                                             */
/* Note: (slightly tricky) All vertex indices in the RWGPanel, */
/* RWGEdge, and RWGHEdge structures contained within an        */
/* OpenSurface refer to the Vertex[] array in the parent       */
/* RWGComposite structure. However, the *panel* indices in     */
/* the RWGEdge and RWGHEdge structures refer to the Panels[]   */
/* array within the OpenSurface structure, NOT the Panels[]    */
/* array within the parent RWGComposite structure.             */
/***************************************************************/
typedef struct OpenSurface
{  
  int NumPanels;
  RWGPanel *Panels;

  int NumEdges;
  RWGEdge *Edges;   // internal edges to which we assign a full RWG function

  int NumHEdges;
  RWGHEdge *HEdges; // external edges to which we assign a half-RWG function

  int NumTotalEdges;

} OpenSurface;

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
   /*--------------------------------------------------------------*/
   /*--------------------------------------------------------------*/

   /*--------------------------------------------------------------*/
   /*- class data -------------------------------------------------*/
   /*--------------------------------------------------------------*/
   // the x,y,z coordinate of the nth vertex are Vertices[3*n + 0,1,2]
   double *Vertices;  
   int NumVertices;

   // this is a full list of *all* panels on all OpenSurfaces 
   // in the composite. 
   RWGPanel *Panels;
   int NumPanels;

   // for nos = 0 , 1, ..., NumOpenSurfaces-1, 
   // Regions[2*nos+0] and Regions[2*nos+1] are the indices of 
   // the two regions bounded by OpenSurface #nos
   int *Regions;
   int NumRegions;

   // RegionMPs[0] = material properties of EXTERIOR medium
   // RegionMPs[r] = material properties of internal region r for r=1,...,NumRegions
   MatProp **RegionMPs;

   // OpenSurface structures for each open boundary subsurface
   OpenSurface *OpenSurfaces;
   char *OpenSurfaceLabels;
   int NumOpenSurfaces;
   int *NumPanelsPerOpenSurface;

   int TotalBFs;

   // BFIndexOffset[nosa] is the index of the first basis function
   // on open surface #nosa within the vector of all basis functions 
   // on this RWGComposite.
   int *BFIndexOffset;

   // EpsTF[nr] = Epsilon for region #nr at the present frequency
   //             (nr==0 for exterior medium)                   
   //             (note: TF stands for 'this frequency')
   cdouble *EpsTF;
   cdouble *MuTF;

   char *Label;
   char *ErrMsg; // used to indicate when something goes wrong
   char *MeshFileName;

   /* GT encodes any transformation that has been carried out since */
   /* the object was read from its mesh file (not including a       */
   /* possible one-time GTransformation that may have been specified*/
   /* in the .scuffgeo file when the object was first created.)     */
   GTransformation *GT;

   /*--------------------------------------------------------------*/
   /*- private class methods --------------------------------------*/
   /*--------------------------------------------------------------*/
   void InitRWGObject(const char *pMeshFileName, const GTransformation *OTGT=0);
   void ReadGMSHFile(FILE *MeshFile, char *FileName, const GTransformation *GT);
   void InitEdgeList();

 };
