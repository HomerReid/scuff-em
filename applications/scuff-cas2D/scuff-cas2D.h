/*
 * scuff-cas2D.h -- header file for scuff-cas2D
 * 
 * homer reid -- 2/2007 -- 3/2009
 *
 */

#ifndef CASIMIR2D_H
#define CASIMIR2D_H

#include <stdio.h>
#include <stdlib.h>

#include <libTDRT.h>

/***************************************************************/
/* constants                                                   */
/***************************************************************/
#define QUANTITY_ENERGY 1
#define QUANTITY_XFORCE 2
#define QUANTITY_YFORCE 4
#define QUANTITY_ANYFORCE 6

/***************************************************************/
/* C2DWorkspace is the primary workspace structure passed      */
/* around among the various routines in Casimir2D.             */
/***************************************************************/
typedef struct C2DWorkspace
 { 
   /* information on the geometry we are simulating */
   TDRTGeometry *G;
   int N; /* total dimension of BEM matrix */
   char *TransListName;

   /* images for each object if a ground plane is present */
   TDRTObject **ImageObjects;

   /* information on the transformations the user requested */
   char **TransLines, **Tags;
   int WhichQuantities, TETM, GroundPlane;
   int NumTransforms, NumQuantities, NTNQ;

   /* these fields are used in case the user asks for the */
   /* contributions to the force/energy from basis        */
   /* functions associated with a subset of the set of    */
   /* all internal vertices (IVs).                        */
   /* (If the user did not request this then we set       */
   /* NumContributingIVs to 0).                           */
   double Rectangle[4];
   int NumContributingIVs;
   int *ContributingIVIndices;

   /* matrices and vectors for eigenvalue computations */
   HMatrix **T, **TI, *dT0dY;
   HMatrix ***Uab, **dU0bdX, **dU0bdY;
   HMatrix *M, *dM;
   int *ipiv;
   double *DRMInf;

   /* tables of static segment-segment integral data */
   StaticSSIDataTable **TSSSIDataTables;
   StaticSSIDataTable ***USSSIDataTables;

   /* integration tolerances */
   double AbsTol, RelTol, XQMin;

   /* miscellaneous stuff */
   int NumThreads;
   int *Converged;
   char *ByXQFileName;
   char *ByXiFileName;
   char *ProfileFileName;
   char *CurrentTag;
   int WriteHDF5;
   int IntCache;
   double FixedXi, FixedQ;

} C2DWorkspace;

/***************************************************************/
/* function prototypes *****************************************/
/***************************************************************/
C2DWorkspace *CreateC2DWorkspace(TDRTGeometry *G, char *TransListName,
                                 int WhichQuantities, double *Rectangle,
                                 int nThread, int TETM, 
                                 int GroundPlane, int WriteHDF5, int IntCache,
                                 int VisualizeOnly);
 
void XQIntegrand(C2DWorkspace *EFW, double Xi, double q, double *EF);
void PrintConsoleOutput(C2DWorkspace *W, double *I);

/*--------------------------------------------------------------*/
/* routines for dealing with image objects in the case of a     */
/* metallic ground plane                                        */
/*--------------------------------------------------------------*/
TDRTObject *CreateImageObject(TDRTObject *O);
void WriteImageObjectPPMeshes(C2DWorkspace *W, const char *FileName, const char *Tag);
void AddImageContributionToU(TDRTObject *Oa, TDRTObject *ObImage, double Xi, double q,
                             double EpsOut, double MuOut, int nThread, 
                             StaticSSIDataTable *SSSIDT,
                             HMatrix *Uab, HMatrix *dUabdX, HMatrix *dUabdY);

void AssembleTI(TDRTObject *O, TDRTObject *OImage, double Xi, double q,
                double EpsOut, double MuOut, int nThread, 
                StaticSSIDataTable *SSSIDT, HMatrix *TI, HMatrix *dTIdY);

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define CMETHOD_DC      1

void EvaluateXQIntegral(C2DWorkspace *W, double *I, double *E);
void EvaluateQIntegral(C2DWorkspace *W, double Xi, double *I, double *E);
void EvaluateXiIntegral(C2DWorkspace *W, double Q, double *I, double *E);
void EvaluateMatsubaraSum(C2DWorkspace *W, double T, double *I, double *E);
void ProcessXQList(C2DWorkspace *W, char *FileName);
int CacheRead(C2DWorkspace *W, double Xi, double q, double *EF);

#endif
