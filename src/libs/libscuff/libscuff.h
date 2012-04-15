/*
 * libscuff.h    -- header file for libscuff
 *
 * this file is long, but basically is divided into several digestible
 * sections:
 *
 *  1. class definitions for RWGObject and supporting classes
 *  2. class definition for RWGGeometry
 *  3. non-class methods that operate on RWGPanels or RWGObjects 
 *  4. other lower-level non-class methods 
 *
 * homer reid  -- 3/2007 -- 9/2011
 */

#ifndef LIBSCUFF_H
#define LIBSCUFF_H

#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <complex>
#include <cmath>

#include <libhrutil.h>
#include <libhmat.h>
#include <libMatProp.h>
#include <libIncField.h>

#include "GTransformation.h"
#include "FieldGrid.h"

namespace scuff {

/*--------------------------------------------------------------*/
/*- some constants used to pass values to libscuff functions   -*/
/*--------------------------------------------------------------*/
#define SCUFF_GENERALFREQ  0
#define SCUFF_PUREIMAGFREQ 1

#define SCUFF_NOLOGGING      0
#define SCUFF_TERSELOGGING   1
#define SCUFF_VERBOSELOGGING 2

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* impedance of free space */
#ifndef ZVAC
#define ZVAC 376.73031346177
#endif

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- 1. Class definitions for RWGObject and supporting classes  -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/ 

/***************************************************************/
/* RWGPanel is a structure containing data on a single        */
/* triangular panel in the meshed geometry.                    */
/*                                                             */
/* note: after the following code snippet                      */
/*  double *R= O->Vertices + 3*O->Panels[np]->VI[i]            */
/* we have that R[0..2] are the cartesian coordinates of the   */
/* ith vertex (i=0,1,2) of the npth panel in object O.         */ 
/*                                                             */
/***************************************************************/
typedef struct RWGPanel
 { 
   int VI[3];                   /* indices of vertices in Vertices array */
   double Centroid[3];          /* panel centroid */
   double ZHat[3];              /* normal vector */
   double Radius;               /* radius of enclosing sphere */
   double Area;                 /* panel area */

   int Index;                   /* index of this panel within object(0..NP-1)*/

} RWGPanel;

/***************************************************************/
/* RWGEdge is a structure containing data on a single          */
/* internal edge in the meshed geometry.                       */ 
/*                                                             */
/* note: after the following code snippet                      */
/*  double *R= O->Vertices + 3*E->IQP                          */
/* we have that R[0..2] are the cartesian coordinates of the   */
/* QP vertex in the basis function corresponding to edge E.    */ 
/***************************************************************/
typedef struct RWGEdge 
 { 
   int iV1, iV2, iQP, iQM;	/* indices of panel vertices (iV1<iV2) */
   double Centroid[3];          /* edge centroid */
   double Length;               /* length of edge */
   double Radius;               /* radius of enclosing sphere */

   int iPPanel;                 /* index of PPanel within object (0..NP-1)*/
   int iMPanel;                 /* index of MPanel within object (0..NP-1)*/
   int PIndex;                  /* index of this edge within PPanel (0..2)*/
   int MIndex;                  /* index of this edge within MPanel (0..2)*/
   int Index;                   /* index of this edge within object(0..NE-1)*/

   RWGEdge *Next;               /* pointer to next edge in linked list */

} RWGEdge;

/***************************************************************/
/* fast kd-tree based point-in-object calculations             */
/***************************************************************/
typedef struct kdtri_s *kdtri;
void kdtri_destroy(kdtri t); // destructor

// some tree statistics for informational purposes:
unsigned kdtri_maxdepth(kdtri t);
size_t kdtri_maxleaf(kdtri t);
double kdtri_meandepth(kdtri t);
double kdtri_meanleaf(kdtri t);

/***************************************************************/
/* RWGObject is a class describing a single physical object    */
/* with a surface mesh read in from a mesh file.               */
/***************************************************************/
class RWGObject
 { 
   /*--------------------------------------------------------------*/ 
   /*- class methods ----------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
  public:  

   /* constructor entry point 1: construct from an 'OBJECT...ENDOBJECT' */
   /* section in a .scuffgeo file                                       */ 
   RWGObject(FILE *f, const char *Label, int *LineNum);

   /* constructor entry points 2 and 3: construct from a given mesh file */
   RWGObject(const char *pMeshFileName);
   RWGObject(const char *pMeshFileName, const char *pLabel,
             const char *Material, const GTransformation *OTGT);

   /* constructor entry point 3: construct from a list of vertices */
   RWGObject(double *pVertices, int pNumVertices, 
             int **PanelVertexIndices, int pNumPanels);

   /* destructor */
   ~RWGObject();

   /* get overlap integral between two basis functions */
   double GetOverlap(int neAlpha, int neBeta, double *pOTimes = NULL);

   /* calculate the inner product of a single basis function with  */
   /* given incident electric and magnetic fields                  */
   void GetInnerProducts(int nbf, EHFuncType2 EHFunc, void *EHFuncUD,
                         int PureImagFreq, cdouble *EProd, cdouble *HProd,
			 int exterior_index, int interior_index);

   /* apply a general transformation (rotation+displacement) to the object */
   void Transform(const GTransformation *GT);
   void Transform(char *format, ...);
   void UnTransform();

   /* fast inclusion tests */
   bool Contains(const double X[3]);
   bool Contains(const RWGObject *O);

   /* visualization */
   // void Visualize(double *KVec, double Kappa, char *format, ...);
   void WriteGPMesh(const char *format, ...);
   void WritePPMesh(const char *FileName, const char *Tag, int PlotNormals);
   void WritePPMesh(const char *FileName, const char *Tag);
   void WritePPMeshLabels(const char *FileName, const char *Tag, int WhichLabels);
   void WritePPMeshLabels(const char *FileName, const char *Tag);

   /* calculate spherical multipole moments due to a single basis function */
#if 0
   void Get1BFSphericalMoments(int ne, double *X0, int lMax,
                               double Wavevector, int RealFreq, 
                               cdouble *aE, cdouble *aM);
#endif

//  private:

   /*--------------------------------------------------------------*/
   /*- private data fields  ---------------------------------------*/
   /*--------------------------------------------------------------*/
   MatProp *MP;                   /* material properties */
   cdouble EpsThisFreq;           // permittivity and permeability
   double MuThisFreq;             //   at the current frequency

   double *Vertices;              /* Vertices[3*n,3*n+1,3*n+2]=nth vertex coords */
   RWGPanel **Panels;             /* array of pointers to panels */
   RWGEdge **Edges;               /* array of pointers to edges */
   RWGEdge **ExteriorEdges;       /* array of pointers to exterior edges */

   int NumVertices;                /* number of vertices in mesh */
   int NumInteriorVertices;        /* number of interior vertices*/
   int NumRefPts;                  /* number of vertices used as reference points */

   int NumTotalEdges;              /* total number of edges */
   int NumEdges;                   /* number of interior edges */
   int NumExteriorEdges;           /* number of exterior edges */

   int NumBFs;                     /* number of basis functions */
   int NumPanels;                  /* number of panels */
   int NumBCs;                     /* number of boundary countours */

   int Index;                      /* index of this object in geometry  */

   int *WhichBC;                   /* WhichBC[nv] = index of boundary contour */
                                   /* on which vertex #nv lies (=0 if vertex  */
                                   /* #nv is an internal vertex)              */
   int *NumBCEdges;                /* NumBCEdges[2] is the number of edges in */
                                   /* boundary contour #2                     */
   RWGEdge ***BCEdges;             /* BCEdges[2][3] is a pointer to the 3rd   */
                                   /* edge in boundary contour #2             */

   int NumRedundantVertices;

   char *MeshFileName;             /* saved name of mesh file */
   char *Label;                    /* unique label identifying object */

   RWGObject *ContainingObject;    /* pointer to object containing this object, if any */

   char *ContainingObjectLabel;    /* these fields are only used by */
   char *MPName;                   /* the class constructor         */
   char *ErrMsg;

   kdtri kdPanels; /* kd-tree of panels */
   void InitkdPanels(bool reinit = false, int LogLevel = SCUFF_NOLOGGING);
  
   /* GT encodes any transformation that has been carried out since */
   /* the object was read from its mesh file (not including a       */
   /* possible one-time GTransformation that may have been specified*/
   /* in the .scuffgeo file when the object was first created.)     */
   GTransformation *GT;

   /*--------------------------------------------------------------*/ 
   /*- private class methods --------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
   /* the actual body of the class constructor */
   void InitRWGObject(const char *pMeshFileName, const char *pLabel, 
                      const char *Material, const GTransformation *GT);

   /* constructor subroutines */
   void InitEdgeList();
   void ReadGMSHFile(FILE *MeshFile, char *FileName, const GTransformation *GT);
   void ReadComsolFile(FILE *MeshFile, char *FileName, const GTransformation *GT);

   /* calculate reduced potentials due to a single basis function */
   /* (this is a helper function used to implement the            */
   /*  GetInnerProducts() class method)                           */
   void GetReducedPotentials(int ne, const double *X, cdouble K,
                             cdouble *a, cdouble *Curla, cdouble *Gradp);
 
 };

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- 3. Class definition for RWGGeometry                        -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/***************************************************************/
/* an RWGGeometry is a collection of RWGObjects.               */
/***************************************************************/
class RWGGeometry 
 { 
   /*--------------------------------------------------------------*/ 
   /*- class methods ----------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
  public:  

   /* constructor / destructor */
   RWGGeometry(const char *GeoFileName, int pLogLevel = SCUFF_NOLOGGING);
   ~RWGGeometry();

   /* geometrical transformations */
   void Transform(GTComplex *GTC);
   void UnTransform();
   char *CheckGTCList(GTComplex **GTCList, int NumGTCs);

   /* visualization */
   void WritePPMesh(const char *FileName, const char *Tag, int PlotNormals);
   void WritePPMesh(const char *FileName, const char *Tag);
   void WriteGPMesh(const char *format, ...);
   void WriteGPMeshPlus(const char *format, ...);
   void PlotSurfaceCurrents(HVector *KN, double Frequency, int RealFreq, 
                            const char *format, ...);

   /* routines for allocating, and then filling in, the BEM matrix */
   HMatrix *AllocateBEMMatrix(bool PureImagFreq = false, bool Packed = false);

   HMatrix *AssembleBEMMatrix(cdouble Omega, HMatrix *M = NULL, int nThread = 0);

#if 0 // delete me 20120408
   void AssembleBEMMatrix(cdouble Frequency, HMatrix *M, int nThread = 0);
   HMatrix *AssembleBEMMatrix(cdouble Frequency, int nThread = 0) {
	HMatrix *M = AllocateBEMMatrix(real(Frequency) == 0.0);
	AssembleBEMMatrix(Frequency, M, nThread);
	return M;
   }
#endif
   
#if 0
   /* routines for allocating, and then filling in, the derivative */
   /* of the bem matrix w.r.t. the coordinates of a mesh vertex    */
   HMatrix *AllocateDMDVMatrix(int RealFreq);
   void AssembleDMDVMatrix(int ObjectIndex, int VertexIndex, int Mu, 
                           double Frequency, int RealFreq,
                           HMatrix *DMDV, int nThread = 0);
   void GetMEVertexDerivative(void *pLFW, 
                              RWGObject *O, int ne, 
                              RWGObject *OP, int nep,
                              int VertexIndex, int Mu, 
                              double AverageRadius,
                              double Frequency, int RealFreq,
                              cdouble *dmdv);
#endif

   /* routines for allocating, and then filling in, the RHS vector */
   HVector *AllocateRHSVector(int PureImagFreq);
   HVector *AllocateRHSVector() { return AllocateRHSVector(0); }

   HVector *AssembleRHSVector(cdouble Omega, IncField *IF,
                              HVector *RHS = NULL, int nThread = 0);

   // update ObjectIndex, and optionally Omega/Eps/Mu, of IF:
   int UpdateIncFields(IncField *IF, cdouble Omega);

   // get the index of the object containing point X
   int GetObjectIndex(const double X[3]);
   RWGObject *GetObject(const double X[3]);

   /* routine for evaluating scattered fields. */
   /* in the first two entry points, the caller already knows     */
   /* which object the evaluation point lies inside (possibly the */
   /* exterior medium), which saves time.                         */
   /* in the third entry point, the code automatically determines */
   /* which object the evaluation point lies inside.              */
   void GetFields(const double X[3], int ObjectIndex,
                  cdouble Omega, HVector *KN, cdouble EH[6], int nThread = 0);
   void GetFields(const double X[3], const char *ObjectLabel,
                  cdouble Omega, HVector *KN, cdouble EH[6], int nThread = 0);
   void GetFields(const double X[3], 
                  cdouble Omega, HVector *KN, cdouble EH[6], int nThread = 0);

   /****************************************************/

   /* Routine for evaluating arbitrary functions of the fields on a 2d
      surface grid; see FieldGrid.h/cc.  Returns a NULL-terminated
      (malloc'ed) array of HMatrix pointers.  If a non-NULL incident
      field function is supplied, then the total of scattered plus
      incident fields is used.  If KN == NULL, then the scattered
      fields are set to zero. */
   HMatrix **GetFieldsGrids(SurfaceGrid &grid, int nfuncs, FieldFunc **funcs,
			    cdouble Omega, HVector *KN,
			    EHFuncType2 EHFunc, void *EHFuncUD,
			    int nThread = 0);
   HMatrix **GetFieldsGrids(SurfaceGrid &grid, int nfuncs, FieldFunc **funcs,
			    cdouble Omega, HVector *KN,
			    EHFuncType EHFunc, void *EHFuncUD, // in exterior
			    int nThread = 0);
   HMatrix **GetFieldsGrids(SurfaceGrid &grid, int nfuncs, FieldFunc **funcs,
			    cdouble Omega, HVector *KN = NULL,
			    IncField *inc = NULL, int nThread = 0);

   // exprs is a string of COMMA-SEPARATED expressions
   HMatrix **GetFieldsGrids(SurfaceGrid &grid, const char *exprs,
                            cdouble Omega, HVector *KN=NULL, IncField *inc=NULL,
			    int nThread = 0);

   // variants that only compute one function and return one matrix
   HMatrix *GetFieldsGrid(SurfaceGrid &grid, FieldFunc &func,
			  cdouble Omega, HVector *KN=NULL, IncField *inc=NULL,
			  int nThread = 0);
   HMatrix *GetFieldsGrid(SurfaceGrid &grid, const char *expr,
			  cdouble Omega, HVector *KN=NULL, IncField *inc=NULL,
			  int nThread = 0);

   /****************************************************/

   /* routine for calculating electric and magnetic dipole moments */
   void GetDipoleMoments(double Frequency, int RealFreq, HVector *KN, 
                         cdouble (*PM)[6], int nThread = 0);


   /* routine for calculating spherical multipole moments */
#if 0
   void GetSphericalMoments(int WhichObject, const double X0[3], 
                            double Frequency, int RealFreq,
                            HVector *KN, 
                            cdouble *aE, cdouble *aM, int lMax,
			    int nThread = 0);
#endif

   /* routine for computing the expansion coefficients in the RWG basis */
   /* of an arbitrary user-supplied surface-tangential vector field     */
   void ExpandCurrentDistribution(EHFuncType EHFunc, void *EHFuncUD, 
                                  HVector *KNVec, int nThread = 0);
   void ExpandCurrentDistribution(IncField *inc, HVector *KNv, int nT=0) {
	ExpandCurrentDistribution(EHIncField, (void*)inc, KNv, nT);
   }

   /* evaluate the surface currents at a given point X on an object */
   /* surface, given a vector of RWG expansion coefficients         */
   void EvalCurrentDistribution(const double X[3], HVector *KNVec, cdouble KN[6]);

   /* routine for setting logging verbosity */
   void SetLogLevel(int LogLevel);

   /* some simple utility functions */
   int GetDimension();
   RWGObject *GetObjectByLabel(char *Label, int *WhichObject = NULL);

   /*--------------------------------------------------------------*/ 
   /*- private data fields  ---------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
//private:
   RWGObject **Objects;             /* array of pointers to objects */
   int NumObjects;
   int TotalBFs;
   int TotalPanels;
   double AveragePanelArea;

   MatProp *ExteriorMP;             /* material properties of exterior medium*/
   cdouble EpsThisFreq;             /* permittivity and permeability of      */
   double MuThisFreq;               /* exterior medium at this frequency     */
   int AllPEC;                      /* = 1 if all objects are PEC bodies     */

   int Verbose;

   char *GeoFileName;

   /* BFIndexOffset[2] is the index within the overall BEM        */
   /* system vector of the first basis function on object 2. thus */
   /*  BFIndexOffset[0]=0                                         */
   /*  BFIndexOffset[1]=# BFs for object 0                        */
   /*  BFIndexOffset[2]=# BFs for object 0 + # BFs for object 1   */
   /* etc.                                                        */
   int *BFIndexOffset;
   int *PanelIndexOffset;

   /* Mate[i]=j if j<i and object j is identical to object i. */
   /* Mate[i]=-1 if there is no object identical to object i, */
   /* OR if there is an object identical to object i but its  */
   /* index is greater than i.                                */
   /* for example if objects 2 and 3 are identical then       */
   /*  Mate[2]=-1, Mate[3]=2.                                 */
   /* note: two objects are identical if  	              */
   /*  (1) they have the same mesh file, and                  */
   /*  (2) they have the same material properties.            */
   int *Mate;

   /* ObjectMoved[i] = 1 if object #i was moved on the most   */
   /* recent call to Transform(). otherwise ObjectMoved[i]=0. */
   int *ObjectMoved;
  
   // int GetObjectAndEdgeIndex(int ei, RWGObject **pO);

   int LogLevel; 

   // short-wavelength panel-panel-integral tolerance
   static double SWPPITol; 

 };

/***************************************************************/
/* non-class methods that operate on RWGPanels                 */
/***************************************************************/
RWGPanel *NewRWGPanel(double *Vertices, int iV1, int iV2, int iV3);
void InitRWGPanel(RWGPanel *P, double *Vertices);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- 4. some other lower-level non-class methods                -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
  
/* 3D vector manipulations */
void VecZero(double v[3]);
double *VecCopy(const double v1[3], double v2[3]);
double *VecScale(const double v1[3], double alpha, double v2[3]);
double *VecScale(double v[3], double alpha);
double *VecScaleAdd(const double v1[3], double alpha, const double v2[3], double v3[3]);
double *VecAdd(const double v1[3], const double v2[3], double v3[3]);
double *VecSub(const double v1[3], const double v2[3], double v3[3]);
double *VecPlusEquals(double v1[3], double alpha, const double v2[3]);
double *VecCross(const double v1[3], const double v2[3], double v3[3]);
double *VecLinComb(double alpha, const double v1[3], double beta, const double v2[3], 
                   double v3[3]);
double VecDot(const double v1[3], const double v2[3]);
double VecDistance(const double v1[3], const double v2[3]);
double VecDistance2(const double v1[3], const double v2[3]);
double VecNorm(const double v[3]);
double VecNorm2(const double v[3]);
double VecNormalize(double v[3]);

/* routines for creating the 'Gamma Matrix' used for torque calculations */
void CreateGammaMatrix(double *TorqueAxis, double *GammaMatrix);
void CreateGammaMatrix(double TorqueAxisX, double TorqueAxisY, 
                       double TorqueAxisZ, double *GammaMatrix);
void CreateGammaMatrix(double Theta, double Phi, double *GammaMatrix);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void PreloadGlobalFIPPICache(const char *FileName);
void StoreGlobalFIPPICache(const char *FileName);

} // namespace scuff

#endif // #ifndef LIBSCUFF_H
