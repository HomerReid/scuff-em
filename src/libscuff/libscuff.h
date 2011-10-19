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

#include <libhmat.h>
#include <libMatProp.h>

#include "StaticPPI.h"

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- 0 couple of quick things before we begin                   -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

#ifndef cdouble
  typedef std::complex<double> cdouble;
#endif 

/* values for the RealFreq parameter to functions */
#define IMAG_FREQ 0
#define REAL_FREQ 1

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
#define MAXBCS 10  // maximum number of external boundary contours

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
/* an RWGObject is an object describing a single contiguous    */
/* object read in from a mesh file.                            */
/***************************************************************/
class RWGObject
 { 
   /*--------------------------------------------------------------*/ 
   /*- class methods ----------------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
  public:  

   /* constructor entry points, type 1: construct from mesh file   */
   RWGObject(const char *pMeshFileName);
   RWGObject(const char *pMeshFileName, const char *pLabel, 
             const char *Material, const char *RotFileName, double *DX);

   /* constructor entry points, type 2: construct from list of vertices */
   RWGObject(double *pVertices, int pNumVertices, 
             int **PanelVertexIndices, int pNumPanels);

   /* destructor */
   ~RWGObject();

   /* move the object */
   void Displace(double DX[3]); /* displace object through DX */
   void Displace(double dx, double dy, double dz);
   void UnDisplace();

   /* rotate the object around an axis */
   void Rotate(double *ZHat, double Theta);

   /* apply a general transformation (rotation+displacement) to the object */
   int Transform(const char *TransLine, ... );
   void UnTransform();

   void TransformPoint(double *X);
   void UnTransformPoint(double *X);

   /* visualization */
   void Visualize(double *KVec, double Kappa, char *format, ...);
   void WriteGPMesh(const char *format, ...);
   void WritePPMesh(const char *FileName, const char *Tag, int PlotNormals);
   void WritePPMesh(const char *FileName, const char *Tag);
       
   /* calculate reduced potentials due to a single basis function */
   void GetReducedPotentials(int ne, double *X, cdouble K,
                             cdouble *a, cdouble *Curla, cdouble *Gradp);

   /* calculate spherical multipole moments due to a single basis function */
   void Get1BFSphericalMoments(int ne, double *X0, int lMax,
                               double Wavevector, int RealFreq, 
                               cdouble *aE, cdouble *aM);

   /* calculate the inner product of a single basis function with 
      given incident electric and magnetic fields 
   */
   void GetInnerProducts(int nbf, EHFuncType EHFunc, void *EHFuncUD, int RealFreq,
                         cdouble *EProd, cdouble *HProd);

   /* get overlap between two basis functions */
   double GetOverlap(int neAlpha, int neBeta);

//  private:

   /*--------------------------------------------------------------*/
   /*- private data fields  ---------------------------------------*/
   /*--------------------------------------------------------------*/
   MatProp *MP;                   /* material properties */
   cdouble EpsThisFreq;           // permittivity and permeability
   double MuThisFreq;             // at the current frequency

   double *Vertices;              /* Vertices[3*n,3*n+1,3*n+2]=nth vertex coords */
   RWGPanel **Panels;             /* array of pointers to panels */
   RWGEdge **Edges;               /* array of pointers to edges */

   int NumVertices;                /* number of vertices in mesh */
   int NumInteriorVertices;        /* number of interior vertices*/
   int NumRefPts;                  /* number of vertices used as reference points */

   int NumTotalEdges;              /* total number of edges */
   int NumEdges;                   /* number of interior edges */

   int NumBFs;                     /* number of basis functions */
   int NumPanels;                  /* number of panels */
   int NumBCs;                     /* number of boundary countours */

   int *WhichBC;                   /* WhichBC[nv] = index of boundary contour */
                                   /* on which vertex #nv lies (=0 if vertex  */
                                   /* #nv is an internal vertex)              */
   int NumBCEdges[MAXBCS];         /* NumBCEdges[2] is the number of edges in */
                                   /* boundary contour #2                     */
   RWGEdge **BCEdges[MAXBCS];      /* BCEdges[2][3] is a pointer to the 3rd   */
                                   /* edge in boundary contour #2             */

   int NumRedundantVertices;

   char *MeshFileName;             /* saved name of mesh file */
   char *Label;                    /* unique label identifying object */

   RWGObject *ContainingObject;    /* pointer to object containing this object, if any */
  
   /* MT and VT encode any transformation that has been carried out */
   /* since the object was read from its mesh file. if X=[X1 X2 X3] */
   /* are the coordinates of a vertex in the original mesh, then    */
   /* XT=MT*X + VT are the transformed coordinates.                 */
   double MT[3][3];
   double VT[3];

   int Index;                      /* index of this object in geometry  */
 
   StaticPPIDataTable *SPPIDTable;

   /*--------------------------------------------------------------*/ 
   /*- private class methods --------------------------------------*/ 
   /*--------------------------------------------------------------*/ 

   /* actual body of constructor */
   void InitRWGObject(const char *pMeshFileName, const char *pLabel, 
                      const char *Material, const char *RotFileName, double *DX);

   /* constructor subroutines */
   void InitEdgeList();
   void InitSPPIDTable();
   void ReadGMSHFile(FILE *MeshFile, char *FileName, double *RotMat, double *DX);
   void ReadComsolFile(FILE *MeshFile, char *FileName, double *RotMat, double *DX);

   /* utility routines */
   int CountCommonVertices(int np1, int np2, int *Index1, int *Index2);
   int CountCommonVertices(int np1, int np2) 
    { return CountCommonVertices(np1, np2, 0, 0); }
 
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

   /* get dimension of linear system */
   int GetDimension();

   /* geometrical transformations */
   int Transform(char *TransLine, char *Tag, char *ErrMsg);
   int Transform(char *TransLine) { return Transform(TransLine, 0, 0); }
   void UnTransform();

   /* visualization */
   void WritePPMesh(const char *FileName, const char *Tag, int PlotNormals);
   void WritePPMesh(const char *FileName, const char *Tag);
   void WriteGPMesh(const char *format, ...);
   void WriteGPMeshPlus(const char *format, ...);
   void WriteMLMesh(const char *format, ...);
   void PlotSurfaceCurrents(HVector *KN, double Frequency, int RealFreq, 
                            const char *format, ...);

   /* initialize tables for accelerating computation of matrix elements */
   void PreCompute(int nThread);

   /* routines for allocating, and then filling in, the BEM matrix */
   HMatrix *AllocateBEMMatrix(int RealFreq);
   void AssembleBEMMatrix(double Frequency, int RealFreq, 
                          int nThread, HMatrix *M);

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

   /* routines for filling in subblocks of the BEM matrix */
   void AssembleT(int no, double Frequency, int RealFreq, 
                  int nThread, HMatrix *T);

   void AssembleU(int noa, int nob, double Frequency, int RealFreq,
                  int NumTorqueAxes, double *GammaMatrix,
                  int nThread, HMatrix *Uab,
                  HMatrix *dUdX, HMatrix *dUdY, HMatrix *dUdZ,
                  HMatrix *dUdTheta1, HMatrix *dUdTheta2, HMatrix *dUdTheta3);

   /* lowest-level matrix-element routine that stamps a 1x1, 1x2, 2x1, or 2x2 */
   /* block of matrix elements (corresponding to the interactions of a single */
   /* pair of basis functions) into the BEM matrix or a subblock thereof      */
   void StampMatrixElements(void *pLFW, 
                            RWGObject *Oa, int nea, int OffsetA, 
                            RWGObject *Ob, int neb, int OffsetB, 
                            double Frequency, int RealFreq,
                            int NumTorqueAxes, double *GammaMatrix,
                            HMatrix *M, 
                            HMatrix *dMdX, HMatrix *dMdY, HMatrix *dMdZ, 
                            HMatrix *dMdT1, HMatrix *dMdT2, HMatrix *dMdT3);

   /* routines for allocating, and then filling in, the RHS vector */
   HVector *AllocateRHSVector(int RealFreq);
   void AssembleRHSVector(EHFuncType EHFunc, void *UserDataD, int RealFreq,
                          int nThread, HVector *B);

   /* routine for evaluating scattered fields */
   void GetFields(double *X, double Frequency, int RealFreq, 
                  HVector *KN, int nThread, cdouble *EH);
   void GetFields(double *X, int WhichObject, 
                  double Frequency, int RealFreq,
                  HVector *KN, int nThread, cdouble *EH);

   /* routine for calculating electric and magnetic dipole moments */
   void GetDipoleMoments(double Frequency, int RealFreq, HVector *KN, 
                         int nThread, cdouble (*PM)[6]);

   /* routine for calculating spherical multipole moments */
   void GetSphericalMoments(int WhichObject, double *X0, int lMax,
                            double Frequency, int RealFreq,
                            HVector *KN, int nThread,
                            cdouble *aE, cdouble *aM);

   /* justify or delete me ****************************************/
   void ExpandCurrentDistribution(EHFuncType KNFunc, void *KNFuncUD, int RealFreq, 
                                  int nThread, HVector *KNVec);

   void EvalCurrentDistribution(double *X, int RealFreq, HVector *KNVec, cdouble *KN);

   /* routines for evaluating the scattering portions of the electric and magnetic */
   /* dyadic green's functions and the VEV of the maxwell stress tensor            */
   /* note: the integer parameter after the M m is assumed tress tensor            */
   /* these are alternate versions of the above two routines in which the */
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

   /*--------------------------------------------------------------*/ 
   /*- private data fields  ---------------------------------------*/ 
   /*--------------------------------------------------------------*/ 
//private:
   RWGObject **Objects;             /* array of pointers to objects */
   int NumObjects;
   int TotalBFs;
   int TotalPanels;
   double AveragePanelArea;

   MatProp *MP;                     /* material properties of exterior medium*/
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

 };

/***************************************************************/
/* non-class methods that operate on RWGPanels                 */
/***************************************************************/
RWGPanel *NewRWGPanel(double *Vertices, int iV1, int iV2, int iV3);
void InitRWGPanel(RWGPanel *P, double *Vertices);

/***************************************************************/
/* non-class methods that operate on RWGObjects ***************/
/***************************************************************/

// argument structure for GetPanelPanelInts() routine
typedef struct GPPIArgStruct
 { 
   // inputs 
   RWGObject *O1;
   int ne1; 
   RWGObject *O2; 
   int ne2;
   cdouble K; 
   int NeedCross;
   int NumTorqueAxes; 
   double *GammaMatrix;

   // outputs
   cdouble L[3]; 
   double *GradL; 
   cdouble *dLdT;

 } GPPIArgStruct;

void InitGLFArgs(GLFArgs *AS);

// argument structure for GetLFunctions() routine
typedef struct GLFArgStruct
 { 
   // inputs 
   RWGObject *O1;
   int ne1; 
   RWGObject *O2; 
   int ne2;
   cdouble K; 
   int NeedCross;
   int NumTorqueAxes; 
   double *GammaMatrix;

   // outputs
   cdouble L[3]; 
   double *GradL; 
   cdouble *dLdT;

 } GLFArgStruct;

void InitGLFArgs(GLFArgStruct *AS);
void GetLFunctions(GLFArgStruct *AS);

// argument structure for AssembleBEMMatrixBlock() routine
typedef struct ABMBArgStruct
 {
   // input fields to be filled in by caller
   RWGGeometry *G;
   RWGObject *Oa, *Ob;
   cdouble Frequency;
   int nThread;

   int NumTorqueAxes;
   double *GammaMatrix;
  
   int RowOffset, ColOffset;

   int ForceSymmetric;

   // output fields filled in by routine
   HMatrix *B;
   HMatrix **GradB;
   HMatrix **dBdTheta;

   // additional fields used internally
   double Sign;
   cdouble EpsA, EpsB; 
   double MuA, MuA;
   int OaIsPEC, ObIsPEC;

 } ABMBArgStruct;

void InitABMBArgs(ABMBArgStruct *Args);
void AssembleBEMMatrixBlock(ABMBArgStruct *Args);


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
double VecNorm(double *v);
double VecNorm2(double *v);
double VecNormalize(double *v);
void VecZero(double *v);

/* routines for creating the 'Gamma Matrix' used for torque calculations */
void CreateGammaMatrix(double *TorqueAxis, double *GammaMatrix);
void CreateGammaMatrix(double TorqueAxisX, double TorqueAxisY, 
                       double TorqueAxisZ, double *GammaMatrix);
void CreateGammaMatrix(double Theta, double Phi, double *GammaMatrix);

/* miscellaneous miscellany */
void RWGErrExit(const char *format, ...);
void *RWGMalloc(int size);

/*--------------------------------------------------------------*/
/*- prototypes for taylor master routines                      -*/
/*--------------------------------------------------------------*/
void *CreateTMWorkspace();  
void FreeTMWorkspace(void *pTMW);  
cdouble TaylorMaster(void *pTMW, int WhichCase, int WhichG, int WhichH,
                     cdouble GParam, double HParam, 
                     double *V1, double *V2, double *V3, 
                     double *V2P, double *V3P, double *Q, double *QP,
                     double RefVal);


#endif // #ifndef LIBSCUFF_H
