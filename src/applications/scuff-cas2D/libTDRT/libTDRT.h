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
 * libTDRT.h -- header file for libTDRT library for working with
 *           -- two-dimensional rooftop functions
 *
 * homer reid -- 11/2008
 */

#ifndef LIBTDRT_H
#define LIBTDRT_H

#include <stdio.h>
#include <stdarg.h>
#include <libhmat.h>
#include <libMatProp.h>
#include <libMDInterp.h>
#include <libhrutil.h>
#if HAVE_CXX11
#include <unordered_map>
#elif HAVE_TR1
#include <tr1/unordered_map>
#endif

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- libTDRT.h section 1: constants ----------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

#define ZVAC 376.7     /* impedance of the vacuum */

#define DESINGULARIZATION_RADIUS 10.0

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- libTDRT.h section 2: type and structure definitions  ------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/***************************************************************/
/* a StaticSSIDataRecord contains all static segment-segment   */
/* integrals between a given pair of line segments             */
/***************************************************************/
typedef struct StaticSSIDataRecord 
 { 
   double J_00_1, J_00_2, J_00_3;
   double J_10_1, J_10_2, J_10_3;
   double J_01_1, J_01_2, J_01_3;
   double J_11_1, J_11_2, J_11_3, J_11_4;
   double J_20_1, J_20_2, J_20_3;
   double J_21_1, J_21_2, J_21_3, J_21_4;
   double J_12_1, J_12_2, J_12_3, J_12_4;
   double J_02_1, J_02_2, J_02_3;
   double J_31_1, J_31_2, J_31_3, J_31_4;
   double J_22_1, J_22_2, J_22_3, J_22_4;
   double J_13_1, J_13_2, J_13_3, J_13_4;
   double J_X_1, J_X_2, J_X_3;
   double J_Y_1, J_Y_2, J_Y_3;

 } StaticSSIDataRecord;

/***************************************************************/
/*- a StaticSSIDataTable is an array of StaticSSIDataRecords,  */
/*- one for each nearby pair of line segments within a single  */
/*- object or on a pair of objects, together with a hash table */
/*- that stores pointers into the array in a way that allows   */
/*- for efficient retrieval.                                   */
/***************************************************************/
//typedef google::dense_hash_map <unsigned long, StaticSSIDataRecord*> StaticSSIDataMap;
#if HAVE_CXX11
typedef std::unordered_map <unsigned long, StaticSSIDataRecord*> StaticSSIDataMap;
#elif HAVE_TR1
typedef std::tr1::unordered_map <unsigned long, StaticSSIDataRecord*> StaticSSIDataMap;
#endif
typedef struct StaticSSIDataTable
 { 
    StaticSSIDataMap *Map;
    StaticSSIDataRecord *Buffer;

 } StaticSSIDataTable;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- libTDRT.h section 3: class definitions ---------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/***************************************************************/
/* a 'TDRTObject' describes a discretized 2D object.           */
/***************************************************************/
class TDRTObject
 { 
   public:

    /***************************************************************/
    /* class methods ***********************************************/
    /***************************************************************/

    /* constructor entry points */
    TDRTObject(const char *MeshFileName);
    TDRTObject(const char *MeshFileName, const char *Label, 
               const char *Material);

    /* actual body of constructor */
    void InitTDRTObject(const char *pMeshFileName, const char *pLabel, 
                        const char *Material);
 
    /* constructor helper functions */
    void ReadComsolFile(FILE *f);
    void ReadGMSHFile(FILE *f);

    /* class destructor */
    ~TDRTObject();

    /* apply rigid displacements */
    void Displace(double *DX);
    void Displace(double x, double y);
    void UnDisplace();

    /* rotate the object around the Z axis */
    void Rotate(double Theta);
    void UnRotate();

    /* apply a general transformation (rotation+displacement) to the object */
    int Transform(char *TransLine);
    void UnTransform();

    /* visualization */
    void WriteGPMesh(char *FileName);
  
// private:
    /***************************************************************/
    /* private data fields                                         */
    /* (i.e. fields that would be private if we were conscientious */
    /* about hononoring standard C++ abstraction conventions, which*/ 
    /* we are not)                                                 */ 
    /***************************************************************/
    // Vertices[2*n], Vertices[2*n+1] are the x and y cartesian components
    // of the nth vertex. 
    double *Vertices;
    int NumVertices;

    // Segments[2*n], Segments[2*n+1] are the indices (within the Vertices
    // table) of the endpoints of the nth line segment.
    int *Segments;
    int NumSegments;
  
    // IVs[n]=index of nth interior vertex.
    // Neighbors[2*n], Neighbors[2*n+1] = indices of neighbors of IVs[n].
    // note: 'interior vertices' are what we called 'control points'
    // in the paper.
    int *IVs, *Neighbors;
    int NumIVs;                   /* number of interior vertices */ 
    int NumBFs;                   /* number of basis functions   */
                                  /*  =2*NumIVs for PEC objects  */
                                  /*  =4*NumIVs otherwise        */

    double Displacement[2];
    double Rotation;

    /* ObjectWasRotated=1 if object was rotated on most recent call to Transform() */
    int ObjectWasRotated; 

    char *MeshFileName;
    char *Label;
    char UnTransLine[1000];

    MatProp *MP;                   /* material properties */

 };

/***************************************************************/
/* a 'TDRTGeometry' is a collection of discretized 2D objects.   */
/***************************************************************/
class TDRTGeometry
 { 
   public:

    /***************************************************************/
    /* class methods ***********************************************/
    /***************************************************************/
    TDRTGeometry(const char *pGeoFileName);
    ~TDRTGeometry();

    void WriteGPMesh(const char *format, ...);
    void WriteGPMeshPlus(const char *format, ...);
    void WritePPMesh(char *FileName, char *Tag);

    int Transform(char *TransLine, char *Tag, char *ErrMsg);
    int Transform(char *TransLine) { return Transform(TransLine, 0, 0); }
    void UnTransform();

//  private:
    /***************************************************************/
    /* private data fields ******************************************/
    /***************************************************************/
    TDRTObject **Objects;            /* array of pointers to objects */
    int NumObjects;
    int TotalBFs;                  /* total # of basis functions (dimension of BEM matrix) */

    int *Mate;
    char *GeoFileName;
  
    MatProp *MP;                   /* material properties of external medium*/
    int AllPEC;                    /* =0 if all objects are PEC */
    
    /* BFIndexOffset[n] is the index in the overall system vector of */
    /* the first basis function corresponding to object N.           */
    /* thus BFIndexOffset[0]=0                                       */
    /*      BFIndexOffset[1]=Objects[0]->NumBFs                      */
    /*      BFIndexOffset[2]=Objects[0]->NumBFs + Objects[1]->NumBFs */
    /* ...                                                           */
    int *BFIndexOffset;          

    /* ObjectMoved[n]=1 if the nth object was moved on the most recent */
    /* call to Transform().                                            */
    int *ObjectMoved;

    /* ObjectRotated=1 if ANY object was rotated on the most recent */
    /* call to Transform().                                         */
    int ObjectRotated;
/* */
    static int LogLevel;

    static Interp1D *G1G2Interp;
    void InitG1G2Interp();

 };

/*--------------------------------------------------------------*/
/*- 'square cubature rules' (rules for 2D numerical cubature   -*/
/*- over the square with vertices at (0,0) and (1,1))          -*/
/*--------------------------------------------------------------*/
extern const double SCR3[], SCR5[], SCR7[], SCR9[];

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- libTDRT.h section 4: prototypes for routines that assemble -*/
/*- the T and U blocks of the BEM matrix                       -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void AssembleT(TDRTObject *O, double Xi, double q,
               double EpsOut, double MuOut, int nThread,
               StaticSSIDataTable *SSSIDT, HMatrix *T);

void AssembleU(TDRTObject *Oa, TDRTObject *Ob, double Xi, double q,
               double EpsOut, double MuOut, int nThread, 
               StaticSSIDataTable *SSSIDT,
               HMatrix *Uab, HMatrix *dUabdX, HMatrix *dUabdY);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- libTDRT.h section 5: routines for computing the             */
/*- 'L-functions' that go into the BEM matrix elements          */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/*-----------------------------------------------------------------*/
/*- data structure that stores all L-functions between pairs of   -*/
/*- basis functions associated with a given pair of control points.*/
/*- note that LPZNabla, LZPNabla have additional factors of i     -*/
/*- (i.e. the quantity stored here as 'LPZNabla' is really the    -*/
/*- imaginary part of LPZNabla).                                  -*/
/*-----------------------------------------------------------------*/
typedef struct LFBuffer
 { 
   double LPPBullet, LPPNabla,  LPPTimes;
   double LPZBullet, LPZNabla,  LPZTimes;
   double LZPBullet, LZPNabla,  LZPTimes;
   double LZZBullet, LZZNabla,  LZZTimes;

 } LFBuffer;

/*-----------------------------------------------------------------*/
/*- routine for computing L-functions between basis functions      */
/*-----------------------------------------------------------------*/
void ComputeLFunctions(TDRTObject *Oa, int niva, TDRTObject *Ob, int nivb,
                       double Kappa, double q, StaticSSIDataTable *SSSIDT,
                       LFBuffer *L, LFBuffer *dLdX, LFBuffer *dLdY);

/*-----------------------------------------------------------------*/
/*- data structures and routines for evaluating the (u,up)         */
/*- integrals that go into the L-functions                         */
/*-----------------------------------------------------------------*/

typedef struct IPQRData
 { 
    double I_00_0, I_10_0, I_01_0, I_11_0;
    double I_X_1, I_Y_1;

    double I_00_1, I_10_1, I_01_1, I_11_1;
    double I_20_1, I_21_1, I_12_1, I_02_1;
    double I_11_2, I_21_2, I_22_2, I_12_2; 
    double I_31_2, I_13_2;
     
 } IPQRData;

void uupIntegralCubature(double *Xs, double *Xe, double *Xsp, double *Xep,
                         double Alpha, int Order, int NeedDerivatives,
                         StaticSSIDataRecord *SSSIDR, IPQRData *I);

void uupIntegralDuffy(double *Xs, double *Xe, double *Xsp, double *Xep,
                      double Alpha, IPQRData *IPQRD);

void uupIntegralSameSegment(double *Xs, double *Xe, double Alpha, int Flip, 
                            IPQRData *IPQRD);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- libTDRT.h section 6: routines for computing static          */
/*- segment-segment integral data for a given pair of segments  */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
StaticSSIDataRecord *ComputeStaticSSIData(double *Xs, double *Xe,
                                          double *Xsp, double *Xep,
                                          int NeedDerivatives,
                                          StaticSSIDataRecord *SSIDR);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- libTDRT.h section 7: routines for working with tables of    */
/*- segment-segment integral data for pairs of objects          */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
int Nearby(double *Xs, double *Xe, double *Xsp, double *Xep);
StaticSSIDataTable *CreateStaticSSIDataTable(TDRTObject *O, int nThread);
StaticSSIDataTable *CreateStaticSSIDataTable(TDRTObject *Oa, TDRTObject *Ob, int nThread);
void DestroyStaticSSIDataTable(StaticSSIDataTable *SSSIDT);
StaticSSIDataRecord *GetStaticSSIData(StaticSSIDataTable *SSSIDT,
                                      TDRTObject *Oa, int iXs, int iXe, 
                                      TDRTObject *Ob, int iXsp, int iXep, 
                                      StaticSSIDataRecord *OutputBuffer);
void ComputeStaticSSIData_CommonVertex(double *Xs, double *Xe, double *Xsp, double *Xep,
                                       StaticSSIDataRecord *SSSIDR);
void ComputeStaticSSIData_SameSegment(double *Xs, double *Xe, int Flip,
                                      StaticSSIDataRecord *SSSIDR);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void InitHRBesselK();
void HRBesselK(double z, int NeedK2, double *KArray);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- libTDRT.h section 8: 2D vector manipulations  --------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
double VecDot(double *v1, double *v2);
double VecNorm(double *v);
double VecNormalize(double *v);
void VecSub(double *v1, double *v2, double *v3);
void VecScale(double *v, double alpha);
void VecScaleAdd(double *v1, double alpha, double *v2, double *v3);
void VecPlusEquals(double *v1, double alpha, double *v2);
double VecDistance(double *v1, double *v2);
double VecD2(double *v1, double *v2);

#endif  
