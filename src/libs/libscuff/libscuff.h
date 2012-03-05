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

#include "GTransformation.h"

namespace scuff{

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
#define ZVAC 376.73031346177

/* prototype for incident field routine passed to AssembleRHS */
typedef void (*EHFuncType)(double *R, void *UserData, cdouble *EH); 

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
             const char *Material, GTransformation *OTGT);

   /* constructor entry point 3: construct from a list of vertices */
   RWGObject(double *pVertices, int pNumVertices, 
             int **PanelVertexIndices, int pNumPanels);

   /* destructor */
   ~RWGObject();

   /* get overlap integral between two basis functions */
   double GetOverlap(int neAlpha, int neBeta);

   /* calculate the inner product of a single basis function with  */
   /* given incident electric and magnetic fields                  */
   void GetInnerProducts(int nbf, EHFuncType EHFunc, void *EHFuncUD,
                         int PureImagFreq, cdouble *EProd, cdouble *HProd);

   /* apply a general transformation (rotation+displacement) to the object */
   void Transform(GTransformation *GT);
   void UnTransform();

   /* visualization */
   void Visualize(double *KVec, double Kappa, char *format, ...);
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
                      const char *Material, GTransformation *GT);

   /* constructor subroutines */
   void InitEdgeList();
   void ReadGMSHFile(FILE *MeshFile, char *FileName, GTransformation *GT);
   void ReadComsolFile(FILE *MeshFile, char *FileName, GTransformation *GT);

   /* calculate reduced potentials due to a single basis function */
   /* (this is a helper function used to implement the            */
   /*  GetInnerProducts() class method)                           */
   void GetReducedPotentials(int ne, double *X, cdouble K,
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
   RWGGeometry(const char *GeoFileName);
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
   HMatrix *AllocateBEMMatrix(int PureImagFreq);
   HMatrix *AllocateBEMMatrix() { return AllocateBEMMatrix(0); };

   void AssembleBEMMatrix(cdouble Frequency, int nThread, HMatrix *M);

   /* routines for allocating, and then filling in, the derivative */
   /* of the bem matrix w.r.t. the coordinates of a mesh vertex    */
   HMatrix *AllocateDMDVMatrix(int RealFreq);
   void AssembleDMDVMatrix(int ObjectIndex, int VertexIndex, int Mu, 
                           double Frequency, int RealFreq,
                           int nThread, HMatrix *DMDV);
   void GetMEVertexDerivative(void *pLFW, 
                              RWGObject *O, int ne, 
                              RWGObject *OP, int nep,
                              int VertexIndex, int Mu, 
                              double AverageRadius,
                              double Frequency, int RealFreq,
                              cdouble *dmdv);

   /* routines for allocating, and then filling in, the RHS vector */
   HVector *AllocateRHSVector(int PureImagFreq);
   HVector *AllocateRHSVector() { return AllocateRHSVector(0); }

   void AssembleRHSVector(EHFuncType EHFunc, void *UserDataD, 
                          int nThread, HVector *B);

   /* routine for evaluating scattered fields. */
   /* in the first two entry points, the caller already knows     */
   /* which object the evaluation point lies inside (possibly the */
   /* exterior medium), which saves time.                         */
   /* in the third entry point, the code automatically determines */
   /* which object the evaluation point lies inside.              */
   void GetFields(double *X, int ObjectIndex,
                  cdouble Omega, HVector *KN, int nThread, cdouble *EH);
   void GetFields(double *X, const char *ObjectLabel,
                  cdouble Omega, HVector *KN, int nThread, cdouble *EH);
   void GetFields(double *X, 
                  cdouble Omega, HVector *KN, int nThread, cdouble *EH);

   /* routine for calculating electric and magnetic dipole moments */
   void GetDipoleMoments(double Frequency, int RealFreq, HVector *KN, 
                         int nThread, cdouble (*PM)[6]);


   /* routine for calculating spherical multipole moments */
#if 0
   void GetSphericalMoments(int WhichObject, double *X0, int lMax,
                            double Frequency, int RealFreq,
                            HVector *KN, int nThread,
                            cdouble *aE, cdouble *aM);
#endif

   /* routine for computing the expansion coefficients in the RWG basis */
   /* of an arbitrary user-supplied surface-tangential vector field     */
   void ExpandCurrentDistribution(EHFuncType KNFunc, void *KNFuncUD, 
                                  int nThread, HVector *KNVec);

   /* evaluate the surface currents at a given point X on an object */
   /* surface, given a vector of RWG expansion coefficients         */
   void EvalCurrentDistribution(double *X, HVector *KNVec, cdouble *KN);

   /* routines for evaluating the scattering portions of the electric and magnetic */
   /* dyadic green's functions and the VEV of the maxwell stress tensor            */
   void GetGij(double *R, HMatrix *M, int Cholesky, HVector *KN, 
               double Frequency, int RealFreq, int nThread,
               cdouble GE[3][3], cdouble GM[3][3]);

   void GetTijNj(double R[3], double nHat[3], HMatrix *M, int Cholesky, 
                 HVector *KN, double Frequency, int RealFreq, 
                 int nThread, double TE[3], double TM[3]);

   /* alternate entry points to the above two routines in which the Cholesky */
   /* parameter is taken equal to its default value of 0                     */
   void GetGij(double *R, HMatrix *M, HVector *KN, 
               double Frequency, int RealFreq, int nThread,
               cdouble GE[3][3], cdouble GM[3][3]);

   void GetTijNj(double R[3], double nHat[3], HMatrix *M,
                 HVector *KN, double Frequency, int RealFreq, 
                 int nThread, double TE[3], double TM[3]);

   /* optimized and accelerated array-based routines for evaluating the */
   /* stress tensor at several spatial points all at once               */
#if 0
   void AssembleRHSVectorArray(double *RArray, int NumPts, double Xi, 
                               int nThread, HMatrix *KNArray);

   void GetFieldsArray(double *RArray, int NumPts, double Xi, 
                       HMatrix *KNArray, int nThread, 
                       double EArray[][3][3], double HArray[][3][3]);

   void GetGijArray(double *RArray, int NumPts, HMatrix *M, int Cholesky, 
                    HMatrix *KNArray, double Xi, int nThread,
                    double GEArray[][3][3], double GMArray[][3][3]);

   void GetTijNjArray(double *RArray, double *nHatArray, int NumPts,
                      HMatrix *M, int Cholesky, HMatrix *KNArray, 
                      double Xi, int nThread, 
                      double TEArray[][3], double TMArray[][3]);
#endif

   /* routine for setting logging verbosity */
   void SetLogLevel(int LogLevel);

   /* some simple utility functions */
   int GetDimension();
   RWGObject *GetObjectByLabel(char *Label);
   RWGObject *GetObjectByLabel(char *Label, int *WhichObject);

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
  
   int GetObjectAndEdgeIndex(int ei, RWGObject **pO);

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
void VecZero(double *v);
double *VecScale(double *v, double alpha);
double *VecScaleAdd(double *v1, double alpha, double *v2, double *v3);
double *VecAdd(double *v1, double *v2, double *v3);
double *VecSub(double *v1, double *v2, double *v3);
double *VecPlusEquals(double *v1, double alpha, double *v2);
double *VecCross(double *v1, double *v2, double *v3);
double *VecLinComb(double alpha, double *v1, double beta, double *v2, 
                   double *v3);
double VecDot(double *v1, double *v2);
double VecDistance(double *v1, double *v2);
double VecDistance2(double *v1, double *v2);
double VecNorm(double *v);
double VecNorm2(double *v);
double VecNormalize(double *v);

/* routines for creating the 'Gamma Matrix' used for torque calculations */
void CreateGammaMatrix(double *TorqueAxis, double *GammaMatrix);
void CreateGammaMatrix(double TorqueAxisX, double TorqueAxisY, 
                       double TorqueAxisZ, double *GammaMatrix);
void CreateGammaMatrix(double Theta, double Phi, double *GammaMatrix);

/* miscellaneous miscellany */
void RWGErrExit(const char *format, ...);
void *RWGMalloc(int size);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void PreloadGlobalFIPPICache(char *FileName);
void StoreGlobalFIPPICache(char *FileName);

} // namespace scuff

#endif // #ifndef LIBSCUFF_H
