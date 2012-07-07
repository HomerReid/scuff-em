/*
 *
 */

#ifndef RWGCOMPOSITE_H
#define RWGCOMPOSITE_H

#include "libscuff.h"

namespace scuff {

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
  RWGPanel **Panels;

  int NumEdges;
  RWGEdge **Edges;   // internal edges to which we assign a full RWG function

  int NumHEdges;
  RWGEdge **HEdges;  // external edges to which we assign a half-RWG function

  int NumTotalEdges;

} PartialSurface;

/***************************************************************/
/* RWGComposite is a generalization of RWGObject.              */
/***************************************************************/
class RWGComposite
 {
public:
   /*--------------------------------------------------------------*/
   /*- class methods ----------------------------------------------*/
   /*--------------------------------------------------------------*/
   RWGComposite(FILE *f, const char *pLabel, int *LineNum);
   ~RWGComposite();

   /*--------------------------------------------------------------*/
   /*- class data -------------------------------------------------*/
   /*--------------------------------------------------------------*/
   // the x,y,z coordinate of the nth vertex are Vertices[3*n + 0,1,2]
   double *Vertices;  
   int NumVertices;

   // this is a full list of *all* panels on all PartialSurfaces 
   // in the composite. 
   RWGPanel **Panels;
   int NumPanels;

   // for nps = 0 , 1, ..., NumPartialSurfaces-1, 
   // PSSubRegions[2*nps+0] and PSSubRegions[2*nps+1] are the indices of 
   // the two SubRegions bounded by PartialSurface #nps.
   int *PSSubRegions;
   int NumSubRegions;

   // SubRegionMPs[0]   = material properties of exterior medium
   // SubRegionMPs[nsr] = material properties of subregion #nsr for nsr=1,...,NumSubRegions
   MatProp **SubRegionMPs;

   // PartialSurface structures for each section of the composite surface
   PartialSurface **PartialSurfaces;
   char **PartialSurfaceLabels;
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
   void InitEdgeList(PartialSurface *PS);

 };

typedef struct ACCMBArgStruct 
{
  //RWGGeometry *G;
  RWGComposite *CA, *CB;
  cdouble Omega;
  HMatrix *B;
  int RowOffset, ColOffset;
  int nThread;

} ACCMBArgStruct;

void AssembleCCMatrixBlock(ACCMBArgStruct *Args, int nThread=0);
void AddEdgePanelContributions(ACCMBArgStruct *Args, int nThread=0);

HMatrix *GetFields(RWGComposite *C, int SubRegion,
                   IncField *IF, HVector *KN,
                   cdouble Omega, HMatrix *XMatrix,
                   HMatrix *FMatrix=0, char *FuncString=0,
                   int nThread=0);

HVector AssembleRHSVector_Composite(RWGComposite *C,
                                                  cdouble Omega,
                                                  IncField *IF,
                                                  HVector *RHS,
                                                  int nThread=0);

cdouble GetEdgePanelInteraction(double **PV, double **EV, cdouble K);

void PlotSurfaceCurrents(RWGComposite *C, HVector *KN, cdouble Omega,
                         const char *format, ...);

} // namespace scuff

#endif // ifdef RWGCOMPOSITE
